#!/usr/bin/env python3
"""
compare_interfaces.py
---------------------
Automated Comparative Structural Interface Analyzer for Protein-DNA Complexes.
1. Isolates user-specified protein chains to prevent dimer/complex mismatches.
2. Executes Foldseek to superpose structures in a self-destructing directory.
3. Identifies reference interface residues (<= 4.5A discovery radius).
4. Spatially maps query equivalent residues via CA Euclidean distances (<= 6.0A).
5. Breaks down contact mechanics (Major, Minor, Backbone, H-Bonds).
6. Packages clean, pre-aligned subset PDBs, a TSV, and a PyMOL script into a dedicated folder.
"""

import os
import sys
import argparse
import tempfile
import subprocess
import warnings
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select, NeighborSearch
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.simplefilter('ignore', PDBConstructionWarning)

# --- CONFIG ---
STANDARD_AA = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
DNA_RESIDUES = {'DA','DT','DG','DC','DI','A','T','G','C','I'}
DNA_BACKBONE_ATOMS = {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"}
MAJOR_GROOVE_ATOMS = {'N7','O6','N4','O4','C8'}
MINOR_GROOVE_ATOMS = {'N3','O2','N2',"O4'","C1'","C2'"}
POLAR_ELEMENTS = {'N', 'O'}

COLORS = {
    'H': "0.0, 0.48, 1.0",   # Blue
    'B': "0.86, 0.2, 0.27",  # Red
    'M': "0.15, 0.65, 0.27", # Green
    'm': "1.0, 0.75, 0.0"    # Yellow
}

# --- PDB EXTRACTION CLASSES ---
class SingleChainSelect(Select):
    """Filter to save strictly a specific Protein chain for the Foldseek calculation."""
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return 1 if chain.get_id() == self.chain_id else 0
        
    def accept_residue(self, residue):
        return 1 if residue.get_resname().strip() in STANDARD_AA else 0

class ComplexSelect(Select):
    """Filter to save ONLY the specified Protein and DNA chains for the final clean PDB."""
    def __init__(self, target_chains):
        self.target_chains = set(target_chains)

    def accept_chain(self, chain):
        return 1 if chain.get_id() in self.target_chains else 0
        
    def accept_residue(self, residue):
        resname = residue.get_resname().strip()
        return 1 if (resname in STANDARD_AA or resname in DNA_RESIDUES) else 0

def prep_temp_pdb(in_pdb, out_pdb, chain_id, parser):
    struct = parser.get_structure("tmp", in_pdb)
    io = PDBIO()
    io.set_structure(struct)
    io.save(out_pdb, select=SingleChainSelect(chain_id))

# --- CONTACT LOGIC ---
def classify_dna_atom(atom_name):
    atom_name = atom_name.strip()
    if atom_name in DNA_BACKBONE_ATOMS: return 'B'
    if atom_name in MAJOR_GROOVE_ATOMS: return 'M'
    if atom_name in MINOR_GROOVE_ATOMS: return 'm'
    return 'M'

def infer_hbond(prot_atom, dna_atom, cutoff=3.5):
    if (prot_atom - dna_atom) <= cutoff:
        if prot_atom.element in POLAR_ELEMENTS and dna_atom.element in POLAR_ELEMENTS:
            return True
    return False

def get_atomic_contacts(residue, dna_atoms, radius=4.5):
    ns = NeighborSearch(dna_atoms)
    counts = {'M': 0, 'm': 0, 'B': 0, 'H': 0}
    dna_details = []
    atom_pairs = []

    for p_atom in residue.get_atoms():
        nearby_dna = ns.search(p_atom.coord, radius)
        for d_atom in nearby_dna:
            dtype = classify_dna_atom(d_atom.get_name())
            d_res = d_atom.get_parent()
            d_res_name = f"{d_res.get_resname().strip()}{d_res.id[1]}"
            
            counts[dtype] += 1
            dna_details.append(f"{d_res_name}({dtype})")
            atom_pairs.append((p_atom, d_atom, dtype))
            
            if infer_hbond(p_atom, d_atom):
                counts['H'] += 1
                dna_details.append(f"{d_res_name}(H)")
                atom_pairs.append((p_atom, d_atom, 'H'))
                
    dna_str = ",".join(sorted(list(set(dna_details)))) if dna_details else "-"
    return counts, dna_str, atom_pairs

def find_equivalent_residue(ref_res, query_chain, max_cutoff=6.0):
    if 'CA' not in ref_res: return None, None
    ref_ca = ref_res['CA'].coord
    
    best_res, best_dist = None, float('inf')
    for q_res in query_chain:
        if q_res.get_resname().strip() not in STANDARD_AA or 'CA' not in q_res: continue
        dist = np.linalg.norm(ref_ca - q_res['CA'].coord)
        if dist < best_dist:
            best_dist, best_res = dist, q_res
            
    if best_res and best_dist <= max_cutoff:
        return best_res, best_dist
    return None, None

def generate_cgo_lines(pairs):
    lines = []
    for p_atom, d_atom, ctype in pairs:
        p = p_atom.coord
        d = d_atom.coord
        c = COLORS[ctype]
        lines.append(f"CYLINDER, {p[0]:.3f}, {p[1]:.3f}, {p[2]:.3f}, {d[0]:.3f}, {d[1]:.3f}, {d[2]:.3f}, 0.12, {c}, {c},")
        lines.append(f"COLOR, 1.0, 0.55, 0.0, SPHERE, {p[0]:.3f}, {p[1]:.3f}, {p[2]:.3f}, 0.20,")
        lines.append(f"COLOR, 1.0, 0.55, 0.0, SPHERE, {d[0]:.3f}, {d[1]:.3f}, {d[2]:.3f}, 0.20,")
    return lines

def align_query_with_foldseek(ref_pdb, ref_chain, query_pdb, query_chain, parser):
    with tempfile.TemporaryDirectory() as tmpdir:
        ref_tmp = os.path.join(tmpdir, "ref_isolated.pdb")
        query_tmp = os.path.join(tmpdir, "query_isolated.pdb")
        
        prep_temp_pdb(ref_pdb, ref_tmp, ref_chain, parser)
        prep_temp_pdb(query_pdb, query_tmp, query_chain, parser)
        
        aln_out = os.path.join(tmpdir, "aln.m8")
        tmp_cache = os.path.join(tmpdir, "tmp")
        
        cmd = [
            "foldseek", "easy-search", ref_tmp, query_tmp, 
            aln_out, tmp_cache, "--format-output", "query,target,alntmscore,u,t"
        ]
        
        print(f"[*] Running Foldseek structural alignment on Chains {ref_chain} and {query_chain}...")
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as e:
            print("\n[!] FOLDSEEK CRASHED. ERROR LOG:")
            print(e.stderr)
            sys.exit(1)

        if not os.path.exists(aln_out) or os.stat(aln_out).st_size == 0:
            sys.exit("[!] Foldseek found no structural homology hits between these chains.")

        with open(aln_out, 'r') as f:
            top_hit = f.readline().strip().split('\t')
            
        u_mat = np.array([float(x) for x in top_hit[3].split(',')]).reshape(3, 3)
        t_vec = np.array([float(x) for x in top_hit[4].split(',')])
        
        return u_mat, t_vec

def main():
    parser = argparse.ArgumentParser(description="Comparative Structural Interface Analyzer.")
    parser.add_argument("--ref_pdb", required=True)
    parser.add_argument("--ref_prot", required=True)
    parser.add_argument("--ref_dna1", required=True)
    parser.add_argument("--ref_dna2", required=True)
    parser.add_argument("--query_pdb", required=True)
    parser.add_argument("--query_prot", required=True)
    parser.add_argument("--query_dna1", required=True)
    parser.add_argument("--query_dna2", required=True)
    parser.add_argument("--tsv", default="interface_comparison.tsv")
    parser.add_argument("--pymol", default="render_interfaces.py")
    args = parser.parse_args()

    # --- PROJECT FOLDER SETUP ---
    ref_basename = os.path.splitext(os.path.basename(args.ref_pdb))[0]
    query_basename = os.path.splitext(os.path.basename(args.query_pdb))[0]
    out_dir = f"{ref_basename}_vs_{query_basename}"
    os.makedirs(out_dir, exist_ok=True)

    out_tsv = os.path.join(out_dir, args.tsv)
    out_pymol = os.path.join(out_dir, args.pymol)
    out_ref_pdb = os.path.join(out_dir, f"{ref_basename}_clean.pdb")
    out_query_pdb = os.path.join(out_dir, f"{query_basename}_aligned.pdb")

    pdb_parser = PDBParser(QUIET=True)
    
    # 1. Obtain Alignment Matrix using ISOLATED chains
    u_mat, t_vec = align_query_with_foldseek(args.ref_pdb, args.ref_prot, args.query_pdb, args.query_prot, pdb_parser)

    # 2. Parse FULL Structures
    ref_struct = pdb_parser.get_structure('ref', args.ref_pdb)
    query_struct = pdb_parser.get_structure('query', args.query_pdb)
    ref_model = ref_struct[0]
    query_model = query_struct[0]
    
    # 3. Superpose FULL Query Structure in-memory to shared reference space
    for atom in query_model.get_atoms():
        atom.coord = np.dot(u_mat, atom.coord) + t_vec

    # 4. Save Cleaned, Pre-aligned PDBs to the Project Folder
    io = PDBIO()
    io.set_structure(ref_struct)
    io.save(out_ref_pdb, select=ComplexSelect([args.ref_prot, args.ref_dna1, args.ref_dna2]))
    
    io.set_structure(query_struct)
    io.save(out_query_pdb, select=ComplexSelect([args.query_prot, args.query_dna1, args.query_dna2]))

    ref_dna_atoms = [a for c_id in [args.ref_dna1, args.ref_dna2] if c_id in ref_model for r in ref_model[c_id] if r.get_resname().strip() in DNA_RESIDUES for a in r.get_atoms()]
    query_dna_atoms = [a for c_id in [args.query_dna1, args.query_dna2] if c_id in query_model for r in query_model[c_id] if r.get_resname().strip() in DNA_RESIDUES for a in r.get_atoms()]
    
    if not ref_dna_atoms or args.ref_prot not in ref_model:
        sys.exit("[!] Error: Specified reference chains not found or contain no valid molecular data.")

    ref_interface = []
    for res in ref_model[args.ref_prot]:
        if res.get_resname().strip() not in STANDARD_AA: continue
        counts, dna_str, pairs = get_atomic_contacts(res, ref_dna_atoms)
        if (counts['M'] + counts['m'] + counts['B']) > 0:
            ref_interface.append((res, counts, dna_str, pairs))
            
    ref_interface.sort(key=lambda x: x[0].id[1])

    tsv_headers = [
        "Ref_Anchor", "Ref_M", "Ref_m", "Ref_B", "Ref_H", "Ref_DNA_Targets",
        "Query_Equivalent", "CA_Distance", "Query_M", "Query_m", "Query_B", "Query_H", "Query_DNA_Targets"
    ]

    totals = {'r_M': 0, 'r_m': 0, 'r_B': 0, 'r_H': 0, 'q_M': 0, 'q_m': 0, 'q_B': 0, 'q_H': 0}
    pymol_sites = []

    with open(out_tsv, 'w') as f_tsv:
        f_tsv.write("\t".join(tsv_headers) + "\n")
        
        for ref_res, r_counts, r_dna, r_pairs in ref_interface:
            ref_id = f"{ref_res.get_resname().strip()}{ref_res.id[1]}"
            q_chain = query_model[args.query_prot] if args.query_prot in query_model else []
            q_res, dist = find_equivalent_residue(ref_res, q_chain)
            
            if q_res:
                q_id = f"{q_res.get_resname().strip()}{q_res.id[1]}"
                dist_str = f"{dist:.2f}"
                q_counts, q_dna, q_pairs = get_atomic_contacts(q_res, query_dna_atoms)
            else:
                q_id, dist_str, q_dna = "-", "-", "-"
                q_counts = {'M': 0, 'm': 0, 'B': 0, 'H': 0}
                q_pairs = []

            # Math Tracking
            totals['r_M'] += r_counts['M']; totals['r_m'] += r_counts['m']
            totals['r_B'] += r_counts['B']; totals['r_H'] += r_counts['H']
            totals['q_M'] += q_counts['M']; totals['q_m'] += q_counts['m']
            totals['q_B'] += q_counts['B']; totals['q_H'] += q_counts['H']

            row = [
                ref_id, str(r_counts['M']), str(r_counts['m']), str(r_counts['B']), str(r_counts['H']), r_dna,
                q_id, dist_str, str(q_counts['M']), str(q_counts['m']), str(q_counts['B']), str(q_counts['H']), q_dna
            ]
            f_tsv.write("\t".join(row) + "\n")

            # --- DYNAMIC PYMOL OBJECT CREATION ---
            stick_sel = []
            stick_sel.append(f"(reference_complex and chain {args.ref_prot} and resi {ref_res.id[1]})")
            if q_res:
                stick_sel.append(f"(query_complex and chain {args.query_prot} and resi {q_res.id[1]})")
            
            for _, d_atom, _ in r_pairs:
                d_res = d_atom.get_parent()
                stick_sel.append(f"(reference_complex and chain {d_res.get_parent().id} and resi {d_res.id[1]})")
            for _, d_atom, _ in q_pairs:
                d_res = d_atom.get_parent()
                stick_sel.append(f"(query_complex and chain {d_res.get_parent().id} and resi {d_res.id[1]})")

            pymol_sites.append({
                'name': ref_id,
                'sele_string': " or ".join(set(stick_sel)),
                'cgo_lines': generate_cgo_lines(r_pairs) + generate_cgo_lines(q_pairs)
            })

        total_row = [
            "TOTAL", str(totals['r_M']), str(totals['r_m']), str(totals['r_B']), str(totals['r_H']), "-",
            "TOTAL", "-", str(totals['q_M']), str(totals['q_m']), str(totals['q_B']), str(totals['q_H']), "-"
        ]
        f_tsv.write("\t".join(total_row) + "\n")

    # Generate the Geometric PyMOL Script with Toggleable Groups
    with open(out_pymol, 'w') as f_py:
        f_py.write("import os\nfrom pymol import cmd\nfrom pymol.cgo import *\n\ndef render_scene():\n")
        
        f_py.write("    cmd.bg_color('white')\n    cmd.hide('everything')\n")
        f_py.write(f"    cmd.load('{ref_basename}_clean.pdb', 'reference_complex')\n")
        f_py.write(f"    cmd.load('{query_basename}_aligned.pdb', 'query_complex')\n\n")

        f_py.write("    cmd.show('cartoon', 'reference_complex or query_complex')\n")
        f_py.write("    cmd.color('gray70', 'reference_complex')\n")
        f_py.write("    cmd.color('gray40', 'query_complex')\n")
        f_py.write("    cmd.util.cnc('reference_complex')\n")
        f_py.write("    cmd.util.cnc('query_complex')\n\n")

        for site in pymol_sites:
            f_py.write(f"    # --- SITE: {site['name']} ---\n")
            f_py.write(f"    cmd.create('Sticks_{site['name']}', '{site['sele_string']}')\n")
            f_py.write(f"    cmd.show('sticks', 'Sticks_{site['name']}')\n")
            f_py.write(f"    cmd.hide('cartoon', 'Sticks_{site['name']}')\n")
            
            f_py.write(f"    cgo_{site['name']} = [\n")
            for line in site['cgo_lines']:
                f_py.write(f"        {line}\n")
            f_py.write("    ]\n")
            
            if site['cgo_lines']:
                f_py.write(f"    cmd.load_cgo(cgo_{site['name']}, 'Mech_{site['name']}')\n")
                f_py.write(f"    cmd.group('Site_{site['name']}', 'Sticks_{site['name']} Mech_{site['name']}')\n")
            else:
                f_py.write(f"    cmd.group('Site_{site['name']}', 'Sticks_{site['name']}')\n")
                
            f_py.write(f"    cmd.disable('Site_{site['name']}')\n\n")

        f_py.write("    print('[+] Clean subset PDBs loaded perfectly superposed.')\n")
        f_py.write("    print('[*] USE THE RIGHT PANEL: Click on a Site (e.g., Site_GLN73) to toggle that specific interface.')\n\n")
        f_py.write("render_scene()\n")

    print(f"\n[+] SUCCESS! Created standalone analysis package in: ./{out_dir}/")
    print(f"    1. Cleaned Reference:   {ref_basename}_clean.pdb")
    print(f"    2. Aligned Query:       {query_basename}_aligned.pdb")
    print(f"    3. Publication Data:    {args.tsv}")
    print(f"    4. PyMOL Scene File:    {args.pymol}")
    print("\n[*] To view results, run:")
    print(f"    cd {out_dir}")
    print(f"    pymol {args.pymol}")

if __name__ == "__main__":
    main()
