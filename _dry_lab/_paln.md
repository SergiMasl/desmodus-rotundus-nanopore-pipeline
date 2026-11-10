# HPC:

## Step 1:

Folder directory is located in folder **bat_spleens**, for copy into work space run comand:

```
$ cp -r /projects/class/binf6350_001/bat_spleens /users/student_id/
```

- change student_id to student/worker workspace

**Output:** directory was copied into workspace and now we can work without changing main project

## Step 2:

- Because, _qc_nanopore.py_ is a quality-control (QC) report for a Nanopore direct RNA sequencing run. It reads the sequencing_summary.txt file that MinKNOW writes during a run (one row per read) and produces summary statistics and a 6-panel figure. It only needs pandas, numpy, and matplotlib, with no bioinformatics tools.

- NEVER run .py files outsite of .sh, next step is create a .sh file for safy run qc_nanopore.py.

On cluster run next:

**Option 1:**

```
nano run_qc.sh
```

In nano write:

```
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


# Load Python environment
module load anaconda3/2023.09

# Navigate to working directory
cd "$1"

SUMMARY=/path/to/sequencing_summary.txt   # <-- change this
PREFIX=bat_spleen

if [[ ! -f "$SUMMARY" ]]; then
    echo "ERROR: $SUMMARY not found" >&2
    exit 1
fi

python qc_nanopore.py "$SUMMARY" "$PREFIX"

```

$ run run_qc.sh **your work directory**

**Option 2:**

Copy **run_qc.sh** from _\_dry_lab/bash_files/run_qc.sh_ and run on cluster

**_Output:_**

-nanopore_qc_25084.err
-nanopore_qc_25084.out
