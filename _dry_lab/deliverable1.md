### Deliverable:

**(Figures submitted separately)**

### Summary Table (Need to add a few more summary values):

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

### Analysis:

The first run is significantly better than the second run due to the far greater percentage of passed reads. The MinION’s hardware issues were evident during the first run due to a 
partial reading failure in a sizable portion of the device, likely due to pore blockages as later evidenced by the increased mux_change readings. The first run is decent, with a 
passed read percentage of 73.4%, but could be improved. The second run does not yield an adequate passed read percentage to be considered decent with a value at only 47.9%. 

The end_reason breakdowns support the conclusion that clogged pores were the most significant issue during the sequencing runs, accumulating a far greater percentage in the second run 
than in the first. A flow-cell problem is the likely culprit, explaining the significant loss of proper sequencing capacity in the second run. 