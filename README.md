# Comparative Structural Interface Analyzer for Protein–DNA Complexes

> **Manuscript Context & Disclaimer**
> This repository contains the standalone code submitted in support of the manuscript:
> *"Central role of the transcriptional regulator MprA in coordinating antibiotic resistance, sialic acid metabolism and capsule production"* (currently under review).
>
> The codebase was designed specifically to generate the comparative contact mechanics data hosted on the interactive web resource: [https://biosig.lab.uq.edu.au/mpra/](https://biosig.lab.uq.edu.au/mpra/). While this pipeline is fully functional and can be adapted for other protein–DNA complexes, including monomeric and homodimeric proteins bound to single- or double-stranded DNA, it has not been exhaustively tested outside of the parameters of this study. Researchers applying this tool to their own bespoke structural datasets should do so with standard computational care.

---

## Table of Contents

- [System Requirements \& Installation](#1-system-requirements--installation)
- [Running the Example Dataset](#2-running-the-example-dataset)
- [Interpreting the Output Data](#3-interpreting-the-output-data-tsv)
- [Visualizing in PyMOL](#4-visualizing-the-mechanisms-in-pymol)
- [Working With Your Own Structures](#5-working-with-your-own-structures)
- [Correspondence \& Contact](#6-correspondence--contact)

---

## 1. System Requirements & Installation

This pipeline physically extracts interface geometries, calculates Euclidean distances, and executes [Foldseek](https://github.com/steineggerlab/foldseek) to rigorously align structural homologs.

To ensure strict cross-platform compatibility across modern macOS (Intel and Apple Silicon/M-series) and Linux systems without requiring manual compilation, the setup is handled entirely via Conda.

### Setup Instructions

Ensure you have [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/) installed.

Clone this repository and navigate to the root directory:

```bash
git clone https://github.com/proteinmechanic/Structural-Comparative-Interface-Analyzer.git
cd Structural-Comparative-Interface-Analyzer
```

Run the automated environment setup script:

```bash
# Make the setup script executable and run it
chmod +x setup.sh
./setup.sh

# Activate the dedicated SCIA environment
conda activate scia_env
```

---

## 2. Running the Example Dataset

An example dataset is provided in the `/example` folder to verify your installation. This command compares the *E. coli* OhrR reference complex (PDB: 1Z9C) against the *B. subtilis* OhrR query homolog (PDB: 4LLN).

Run the following command from the repository root:

```bash
python compare_interfaces.py \
  --ref_pdb example/1Z9C.pdb --ref_prot A --ref_dna1 H --ref_dna2 G \
  --query_pdb example/4LLN.pdb --query_prot A --query_dna1 H --query_dna2 G \
  --tsv output_table.tsv \
  --pymol render_interfaces.py
```

### Expected Output

The script prevents workspace clutter by dynamically generating a self-contained project folder (`./1Z9C_vs_4LLN/`) containing:

- Clean subset PDBs
- The TSV contact data
- A PyMOL visualization state script

---

## 3. Interpreting the Output Data (`.tsv`)

The generated tab-separated matrix breaks down the precise atomic contact mechanisms using an explicit **≤ 4.5 Å** discovery radius and maps them across evolutionary homologs.

| Metric | Description |
|---|---|
| `Ref_Anchor` | The specific amino acid in the reference protein mediating the contact. |
| `M` / `m` | Counts of structural contacts localized inside the DNA **M**ajor (M) or **M**inor (m) grooves. The sum of these represents the **Specificity Score**. |
| `B` | Counts of contacts strictly targeting the sugar-phosphate **B**ackbone. This represents the **Anchoring Score**. |
| `H` | Stringent structural **H**ydrogen bonds filtered to polar N/O atom pairs (≤ 3.5 Å). |
| `Query_Equivalent` | The exact functional equivalent residue in the query homolog, mapped via Cα Euclidean distance (≤ 6.0 Å) after Foldseek rigid-body affine transformation. |
| `TOTAL` (bottom row) | High-level sums across all interface contacts. |

---

## 4. Visualizing the Mechanisms in PyMOL

The pipeline automatically generates a highly specific, toggleable PyMOL script for visual validation of the numeric TSV data.

Navigate into the newly created output directory:

```bash
cd 1Z9C_vs_4LLN
```

Launch the PyMOL validation engine:

```bash
pymol render_interfaces.py
```

### How to Use the Viewer

PyMOL will open with the pre-aligned structures displayed as cartoons only. Use the right-hand control panel to click and enable specific interface sites (e.g., `Site_<RESNAME><RESNUM>`). The script will instantly isolate and render only the sticks and geometric cylinders for that specific evolutionary anchor.

**Color Key:**

| Color | Contact Type |
|---|---|
| Blue cylinders | Hydrogen bonds |
| Red | Backbone anchors |
| Green | Major groove contacts |
| Yellow | Minor groove contacts |

---

## 5. Working With Your Own Structures

This pipeline can be applied to X-ray crystallography, Cryo-EM, or AlphaFold3 models. Please note the following guidelines.

### Identifying Chain IDs (CLI Flags)

PDB files often contain biological assemblies with multiple chains. The script deliberately isolates specific chains to prevent Foldseek from mismatching homodimer subunits.

1. Open your custom `.pdb` file in a visualizer such as [PyMOL](https://pymol.org/) or [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/).
2. Click on your primary protein monomer to identify its Chain ID (e.g., `Chain A`).
3. Click on the two complementary strands of your target DNA double helix (e.g., `Chain H` and `Chain G`).
4. Pass these exact IDs to `--ref_prot`, `--ref_dna1`, and `--ref_dna2` (and equivalently for the query structure).

### Watch Out for Missing Regions

When working with experimental PDBs, highly flexible loops or DNA-binding heads are occasionally unresolved due to missing electron density.

If the script reports a `-` for a `Query_Equivalent`, it means Foldseek successfully aligned the core structural folds, but the specific homologous Cα atom was either:

- Entirely absent from the query `.pdb`, or
- Subject to a large conformational hinge-shift placing it beyond the 6.0 Å discovery cutoff.

---

## 6. Correspondence & Contact

For questions regarding the manuscript, the structural biology findings, or the associated interactive web resource, please contact the corresponding author:

**Professor Mark Schembri**
[m.schembri@uq.edu.au](mailto:m.schembri@uq.edu.au)
School of Chemistry and Molecular Biosciences
The University of Queensland, Brisbane, Australia
