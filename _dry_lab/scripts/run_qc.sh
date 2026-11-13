#!/bin/bash
#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16GB
#SBATCH --partition=Centaurus
#SBATCH --job-name=nanopore_qc
#SBATCH --output=nanopore_qc_%j.out
#SBATCH --error=nanopore_qc_%j.err


set -euo pipefail   # stop on the first error

module load anaconda3/2023.09

cd "$1"

SUMMARY=/path/to/sequencing_summary.txt   # <-- change this
PREFIX=bat_spleen

if [[ ! -f "$SUMMARY" ]]; then
    echo "ERROR: $SUMMARY not found" >&2
    exit 1
fi

python qc_nanopore.py "$SUMMARY" "$PREFIX"