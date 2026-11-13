#!/bin/bash
#SBATCH --job-name=qc_a1
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16GB
#SBATCH --partition=Centaurus
#SBATCH --output~/vamp/binf6350_001/logs/qc_a1_%j.out

set -euo pipefail
P=~/vamp/binf6350_001
OUT=$P/results/qc
mkdir -p $OUT

module load anaconda3/2023.09
python --version
python -c "import pandas, numpy, matplotlib; print('pandas', pandas.__version__, '| numpy', numpy.__version__, '| matplotlib', matplotlib.__version__)"

R1=$(find $P/T23_66_Diyo_RNA_spleen -name "sequencing_summary_*.txt")
R2=$(find $P/T23_66_Diyo_RNA_spleen_2ndRun -name "sequencing_summary_*.txt")
echo "Run 1 summary: $R1"
echo "Run 2 summary: $R2"

echo "########## RUN 1 ##########"
python $P/bat_spleens/qc_nanopore.py $R1 $OUT/run1

echo "########## RUN 2 ##########"
python $P/bat_spleens/qc_nanopore.py $R2 $OUT/run2