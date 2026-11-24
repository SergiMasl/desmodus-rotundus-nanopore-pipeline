# Results — Gate B1

| File | What it is | Source |
|---|---|---|
| `mhc_survey.tsv` | Keyword survey of the annotation (41 hits, incl. 6 false hits) | `codes/run_b1.sh` part 1 |
| `mhc_symbol_check.tsv` | Symbol check: TAP1, TAP2, PSMB8, PSMB9 (+ IDH1/2 false hits from the first pattern) | `codes/run_b1.sh` part 2 |
| `gene_list_pairD.tsv` | Final gene list, 39 entries, with a note per gene | decided in Step 4 of the plan |
| `b1_gene_table.tsv` | **Deliverable table**: coordinates, transcript IDs, min/max transcript length | `b1_find_genes.py` (copy from cluster `results/b1/`) |
| `b1_transcripts.tsv` | Every transcript: ID, length (nt), exon count | `b1_find_genes.py` (copy from cluster) |
| `b1_not_found.tsv` | List genes not in the annotation (expected: DPA1, DPB1, DOA, DOB) | `b1_find_genes.py` (copy from cluster) |

Gate B1: **39 genes on the list, 35 located.**
