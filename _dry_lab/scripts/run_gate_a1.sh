#!/bin/bash
#SBATCH --job-name=gate_a1
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16GB
#SBATCH --partition=Centaurus
#SBATCH --output=~/vamp/binf6350_001/logs/gate_a1_%j.out

set -euo pipefail
P=~/vamp/binf6350_001

module load anaconda3/2023.09
python --version
python -c "import pandas; print('pandas', pandas.__version__)"

python $P/scripts/gate_a1.py \
  $P/T23_66_Diyo_RNA_spleen \
  $P/T23_66_Diyo_RNA_spleen_2ndRun \
  --out $P/results/gate_a1_results.tsv
