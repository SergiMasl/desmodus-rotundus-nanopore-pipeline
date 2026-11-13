# Lab Notebook — 2026-11-13 — Stage A1: First Contact (Gate)

## 1. Goal

## Report, for each run:

1. Passed reads
2. Passed bases (Gb, 3 significant figures)
3. % of reads with `end_reason == unblock_mux_change`

## 2. Environment

| Item                 | Value                 |
| -------------------- | --------------------- |
| Cluster node (login) | `gal-i9`              |
| Slurm partition      | `Centaurus`           |
| Module               | `anaconda3/2023.09`   |
| Python               | 3.11.7                |
| pandas               | 2.2.0                 |
| Project folder       | `~/vamp/binf6350_001` |

## 3. Input files

The summary files are **nested inside MinKNOW run subfolders**, not directly in the run folder as the handout implies.

| Run          | Path (relative to project)                                                                                                                                                        | Size   | Lines     | md5                                |
| ------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------ | --------- | ---------------------------------- |
| Run 1 (19 h) | `T23_66_Diyo_RNA_spleen/20261102_1734_MD-101579_FBC09298_cb4c8893/sequencing_summary_FBC09298_cb4c8893_9b6d8c39.txt`                                                              | 323 MB | 1,089,201 | `28707ada50f79640304a6a8364d7a7f5` |
| Run 2 (69 h) | `T23_66_Diyo_RNA_spleen_2ndRun/BINF6350_Fall2026_RNA01/T23_66_Diyo_RNA_spleen_2ndrun/20261103_1357_MD-101579_FBC09298_013e8950/sequencing_summary_FBC09298_013e8950_973790c0.txt` | 30 MB  | 105,294   | `4f8bf12e48de922b29e2c2329847965d` |

Both runs: same flow cell (FBC09298), same device (MD-101579), 22 columns, `end_reason` is the last column.

## 4. Commands

### 4.1 Inspect the files

```bash
cd ~/vamp/binf6350_001
find T23_66_Diyo_RNA_spleen T23_66_Diyo_RNA_spleen_2ndRun -name "sequencing_summary_*.txt" -exec ls -lh {} \;

for f in $(find T23_66_Diyo_RNA_spleen T23_66_Diyo_RNA_spleen_2ndRun -name "sequencing_summary_*.txt"); do
  echo "== $f"; file "$f"; wc -l "$f"
  head -1 "$f" | tr '\t' '\n' | cat -n | tail -3
done
```

`file` output for both: `ASCII text, with CRLF line terminators`.

### 4.2 Reproduce the CRLF bug

```bash
F1=T23_66_Diyo_RNA_spleen/20261102_1734_MD-101579_FBC09298_cb4c8893/sequencing_summary_FBC09298_cb4c8893_9b6d8c39.txt
F2=T23_66_Diyo_RNA_spleen_2ndRun/BINF6350_Fall2026_RNA01/T23_66_Diyo_RNA_spleen_2ndrun/20261103_1357_MD-101579_FBC09298_013e8950/sequencing_summary_FBC09298_013e8950_973790c0.txt

head -2 $F1 | tail -1 | cat -A | tail -c 60                        # ...signal_positive^M$
awk -F'\t' '$22=="unblock_mux_change"' $F1 | wc -l                 # 0      <- BUG
awk -F'\t' '{sub(/\r$/,"")} $22=="unblock_mux_change"' $F1 | wc -l # 63756  <- FIXED
md5sum $F1 $F2
```

## 5. The CRLF bug

|                  |                                                                                                                                                                   |
| ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Symptom**      | `awk` filter `$22=="unblock_mux_change"` returned **0 rows, no error**                                                                                            |
| **Cause**        | MinKNOW ran on Windows, so lines end in `\r\n`. Tools that split on `\n` leave an invisible `\r` on the last field: `unblock_mux_change\r` ≠ `unblock_mux_change` |
| **Diagnosis**    | `file` reports "CRLF line terminators"; `cat -A` shows `^M$`; pandas with `lineterminator="\n"` shows the column name as `'end_reason\r'`                         |
| **Fix**          | Strip `\r` before comparing: `sub(/\r$/,"")` in awk, `.str.strip()` in pandas, or `sed 's/\r$//'`                                                                 |
| **Confirmation** | Run 1: 0 → 63,756. Run 2: 0 → 15,013. awk and pandas give the same Run 1 count (63,756)                                                                           |

Note: pandas `read_csv` with default settings handles CRLF itself (the provided `qc_nanopore.py` says so too). The bug shows up with awk, cut, grep, and pandas with `lineterminator="\n"`. Always check the last column.

## 6. Definitions used

- **Passed read:** `passes_filtering == TRUE` (MinKNOW Q9 filter)
- **Passed bases:** sum of `sequence_length_template` over passed reads
- **unblock_mux_change %:** count ÷ **all** reads in the summary × 100

## 7. Results

Source: `results/gate_a1_results.tsv` (not typed from screen).

| Metric                                  | Run 1       | Run 2      |
| --------------------------------------- | ----------- | ---------- |
| Total reads                             | 1,089,200   | 105,293    |
| **Passed reads**                        | **799,535** | **50,482** |
| Pass rate                               | 73.4%       | 47.9%      |
| **Passed bases (Gb, 3 s.f.)**           | **0.796**   | **0.0528** |
| Passed bases (exact)                    | 795,999,099 | 52,845,928 |
| unblock_mux_change reads                | 63,756      | 15,013     |
| **unblock_mux_change (% of all reads)** | **5.85%**   | **14.26%** |
| unblock_mux_change (% of passed reads)  | 5.34%       | 6.22%      |
| Duplicate read_ids                      | 0           | 0          |

### end_reason breakdown

| end_reason             | Run 1     | Run 2  |
| ---------------------- | --------- | ------ |
| signal_positive        | 1,015,056 | 89,405 |
| unblock_mux_change     | 63,756    | 15,013 |
| device_data_error      | 8,341     | 467    |
| signal_negative        | 1,533     | 289    |
| mux_change             | 250       | 13     |
| analysis_config_change | 220       | 106    |
| paused                 | 44        | 0      |

### Cross-check

The handout's extra-credit table lists 850,017 passing reads and 0.85 Gb.
Run 1 + Run 2 = 799,535 + 50,482 = **850,017** and 0.796 + 0.0528 = **0.849 Gb**.
So that table is both runs combined, and my numbers reconcile exactly.

## 8. Observations (for the Stage 1 interpretation)

- Run 2 (the restart on the same flow cell) produced about a tenth of Run 1's reads despite running 69 h vs 19 h.
- Run 2 pass rate is much lower (47.9% vs 73.4%).
- Run 2 unblock rate is much higher (14.3% vs 5.85% of all reads), consistent with more pore blocking after the restart. Relevant to Pair F's question.
- Watch-out: `qc_nanopore.py` prints passed bases with `:.3f` (3 decimal places), giving `0.053` for Run 2, which is only 2 significant figures. Gate asks for 3 s.f. → use **0.0528**.

## 9. What broke

| #   | Problem                                                               | Fix                                                                                                                    |
| --- | --------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| 1   | `ls $D/...` failed: `$D` never set                                    | Set the variable / use relative paths                                                                                  |
| 2   | Glob `T23_66_Diyo_RNA_spleen*/sequencing_summary_*.txt` found nothing | Files are nested in MinKNOW subfolders; used `find` and a recursive glob                                               |
| 3   | Job 25340: `syntax error near unexpected token 'newline'`             | `module load <python-module>` placeholder: `<` `>` parsed as redirection. Replaced with `anaconda3/2023.09`            |
| 4   | Job 25341: `can't open file ... gate_a1.py`                           | Python script was never copied to the cluster                                                                          |
| 5   | Job 25342 and next attempt: `SyntaxError ... set -euo pipefail`       | Bash wrapper content pasted into `gate_a1.py`                                                                          |
| 6   | `sbatch gate_a1.py`: "Time limit specification required"              | Submitted the `.py` instead of the `.sh`                                                                               |
| —   | Clean restart                                                         | Deleted scripts/logs/results, recreated both files one at a time, checked `head -1` and `py_compile` before submitting |

Lessons: check the first line of a file to see what it is (`#!/usr/bin/env python3` vs `#!/bin/bash`); run `python -m py_compile` before `sbatch`; **`sbatch` the `.sh`, and the `.sh` runs the `.py`.**

## 10. Gate submission

- Submitted: `<date / time>`
- Result: `<opened / wrong — what was fixed>`

## 11. Next

- Run `bat_spleens/qc_nanopore.py` on both runs → six-panel figure and summary table (Stage 1 deliverable)
- Write the one-page QC report: is this a good run, and the single biggest problem
- Stage 2: length bins, the 500–650 bp gate
