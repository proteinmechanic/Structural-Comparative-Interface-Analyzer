#!/bin/bash
# setup_env.sh - Environment Configurator for SCIA

echo "[*] Initializing SCIA Environment Setup..."

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "[!] Error: 'conda' is not installed or not in your PATH."
    echo "    Please install Miniconda or Anaconda first: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "[*] Creating conda environment 'scia_env'..."
conda create -n scia_env -c conda-forge -c bioconda python=3.10 biopython numpy foldseek -y

echo ""
echo "=========================================================================="
echo "[+] Setup Complete!"
echo ""
echo "[+] STEP 1: Activate the environment by running:"
echo "    conda activate scia_env"
echo ""
echo "[+] STEP 2: Run the test example using 1Z9C as Reference and 4LLN as Query:"
echo "    python compare_interfaces.py \\"
echo "      --ref_pdb example/1Z9C.pdb --ref_prot A --ref_dna1 H --ref_dna2 G \\"
echo "      --query_pdb example/4LLN.pdb --query_prot A --query_dna1 H --query_dna2 G \\"
echo "      --tsv output_table.tsv \\"
echo "      --pymol render_interfaces.py"
echo ""
echo "[+] STEP 3: Visually verify the mechanism alignment:"
echo "    pymol render_interfaces.py"
echo "=========================================================================="
