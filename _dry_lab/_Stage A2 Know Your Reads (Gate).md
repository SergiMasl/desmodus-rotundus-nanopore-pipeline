# Lab Notebook Entry — 2026-11-21 — Stage A2: Know Your Reads

## 1. Goal

## Send in as Deliverable #2:

1. Bin the reads and submit as a .csv
2. Create an annotated length distribution plot
3. Get the .fastq files of the bins and send in one example bin
4. Examine the median and N50 values in the read data

## 2. Environment Add-Ons

The following module was included for Stage A2's plotting requirement:

| Item                 | Value                 |
| -------------------- | --------------------- |
| matplotlib           | 3. 9. 4               |



###Created a virtual environment to run matplotlib in:

*Example Environment: 

source ~/venvs/sequencing_qc/bin/activate  

## 3. Script Creation:

Python and bash scripts were made to bin the reads. They can be found in the Scripts reference section:

## 4. Scripts reference

The scripts are in [`scripts/`](scripts/):

- [`run_gate_a2_binning.py`](scripts/run_gate_a2_binning.py) utilizes modules such as pandas and matplotlib to bin the reads, export a .png plot, and organize the bins into respective .fastq files
- [`run_gate_a2_binning.sh`](run_gate_a2_binning.sh) is the Slurm/Bash wrapper for running the binning Python script


## 5. Issues and their Fixes 

### 5.1 Matplotlib Plotting problems: 

The initial .png was unreadable because the descriptive text for each bin count was too close together. Fixing the scrunched text required additional code to the run_gate_a2_binning.py script. 


### 5.2 Fix the Scrunched Text by Adding Additional Padding: 


        counts = df["length_bin"].value_counts(sort=False)

        # Add headroom above the tallest bar so labels sit clear of the data
        ax.set_ylim(0, ax.get_ylim()[1] * 1.35)
        ymax = ax.get_ylim()[1]
        levels = [0.97, 0.88, 0.79]   # stagger label heights so neighbors don't collide

        for i, ((lo, hi), label) in enumerate(zip(spans, BIN_LABELS)):
            lo_c, hi_c = max(lo, xmin), min(hi, xmax * 1.05)
            if lo_c >= hi_c:
                continue
            if i % 2 == 0:
                ax.axvspan(lo_c, hi_c, color="grey", alpha=0.08, lw=0)
            if lo > xmin:
                ax.axvline(lo, color="black", ls="--", lw=0.8)
            x_text = np.sqrt(lo_c * hi_c)  # geometric centre suits a log axis
            ax.text(
                x_text, ymax * levels[i % 3],
                f"{label}\n{int(counts.get(label, 0)):,} reads",
                ha="center", va="top", fontsize=8,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5),
            )

## 6. What Happened/Observations:

Some observations include:
- A large visual spike found in the read length distribution plot at the 500-650 range for Runs 1 and 2
- Overall very similar read shares between Run 1 and Run 2 Bins (Around 10.7 - 23.9 in total range)
- N50 is far larger than the median value (Found in Gate A1 and revisited) 


## 7. Interpretations:

- N50 being higher than the median value means the assembly quality may be good, as reads pass the 50% threshold
- Many read shares are in the 1000-2000 bin range, explaining the high N50 value in comparison to the median value. 
 

## 8. Next

- Create the Stage 2 Deliverable Folder with four files 
- Move to Stage 3
