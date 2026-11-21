#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16GB
#SBATCH --partition=Centaurus
#SBATCH --job-name=diyo_bin_run
#SBATCH --output=diyo_bin_run_%j.out
#SBATCH --error=diyo_bin_run_%j.err


# ==========================================================
# DIRECTORY PATHS
# ==========================================================

SCRIPT_DIR="path/to/working/dir" # change this for your directory

FILE1="path/to/run1/data" # change this for your directory
FILE2="path/to/run2/data" # change this for your directory

# ==========================================================
# CHECK FILES
# ==========================================================

for FILE in "$FILE1" "$FILE2"; do
    if [ ! -f "$FILE" ]; then
        echo "Error: File not found: $FILE"
        exit 1
    fi
done

if [ ! -f "$SCRIPT_DIR/binning.py" ]; then # change binning.py to actual script name
    echo "Error: Python script not found: $SCRIPT_DIR/binning.py" # change binning.py to actual script name
    exit 1
fi

# ==========================================================
# RUN PYTHON SCRIPT
# ==========================================================

echo "=========================================================="
echo "        SEQUENCING RUN COMPARISON"
echo "=========================================================="

python3 "$SCRIPT_DIR/binning.py" "$FILE1" "$FILE2"

if [ $? -ne 0 ]; then
    echo "Error: Python script failed."
    exit 1
fi

echo ""
echo "QC comparison completed successfully."
