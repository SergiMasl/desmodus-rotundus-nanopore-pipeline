### Deliverable:

**(Figures for Run 1 and Run 2 submitted separately)**

### Environment

| Item                 | Value                 |
| -------------------- | --------------------- |
| Cluster node (login) | `gal-i9`              |
| Slurm partition      | `Centaurus`           |
| Module               | `anaconda3/2023.09`   |
| Python               | 3.11.7                |
| pandas               | 2.2.0                 |
| Project folder       | `~/vamp/binf6350_001` |


### Summary Table:

| Metric                                  | Run 1       | Run 2      |
| --------------------------------------- | ----------- | ---------- |
| Total reads                             | 1,089,200   | 105,293    |
| **Passed reads**                        | **799,535** | **50,482** |
| Pass rate                               | 73.4%       | 47.9%      |
| Failed reads                            | 289,665     | 54,811     |
| **Passed bases (Gb, 3 s.f.)**           | **0.796**   | **0.0528** |
| Passed bases (exact)                    | 795,999,099 | 52,845,928 |
| Mean read length                        | 996         | 1,047      |
| Median read length                      | 624         | 662        |
| Read length N50                         | 1,568       | 1,632      |
| Longest read                            | 27,953      | 14,362     |
| Median qscore (passed)                  | 16.3        | 16.1       |
| Run duration (hours)                    | 18.9        | 68.7       |
| Channels producing reads                | 441         | 423        |
| Split reads                             | 118,311     | 7,556      |
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

## Catching the Carriage Return Bug

|                  |                                                                                                                                                                   |
| ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Symptom**      | `awk` filter `$22=="unblock_mux_change"` returned **0 rows, no error**                                                                                            |
| **Cause**        | MinKNOW ran on Windows, so lines end in `\r\n`. Tools that split on `\n` leave an invisible `\r` on the last field: `unblock_mux_change\r` ≠ `unblock_mux_change` |
| **Diagnosis**    | `file` reports "CRLF line terminators"; `cat -A` shows `^M$`; pandas with `lineterminator="\n"` shows the column name as `'end_reason\r'`                         |
| **Fix**          | Strip `\r` before comparing: `sub(/\r$/,"")` in awk, `.str.strip()` in pandas, or `sed 's/\r$//'`                                                                 |
| **Confirmation** | Run 1: 0 → 63,756. Run 2: 0 → 15,013. awk and pandas give the same Run 1 count (63,756)                                                                           |


### Analysis:

The first run is significantly better than the second run due to the far greater percentage of passed reads. The MinION’s hardware issues were evident during the first run due to a 
partial reading failure in a sizable portion of the device, likely due to pore blockages as later evidenced by the increase of mux_change readings. The first run is decent, with a 
passed read percentage of 73.4%, but could be improved. The second run does not yield an adequate passed read percentage to be considered decent with a value at only 47.9%. 

The most significant issue is the substantial sequencing yield drop, and the end_reason breakdown alongside the pore inactivity witnessed when running the MinION support the conclusion 
that malfunctioning or clogged pores were a larger factor than deteriorated RNA quality. Only 441 pores produced reads in the first run, while the second one dropped further in active 
pore count and had a far greater run time. 
