# Plan: Gate B1 — Define the family

**Pair D · Family: MHC class I and class II**
**Due:** Wed 25 Nov 2026 (B1 gate cleared)
**Reference:** `GCF_022682495.2_HLdesRot8A.1_genomic.gtf` (_Desmodus rotundus_; our sample is _Diaemus youngi_)

---

## What the stage asks for (from the handout)

- Read at least five papers about the gene family (some on bats).
- Build a gene list: start from the provided family, then justify additions and removals.
- Find each gene in the reference annotation: coordinates, transcript IDs, expected transcript lengths.
- Document what could not be found, and keep two claims apart:
  - **not annotated in the reference** (a B1 finding)
  - **not expressed in this sample** (a B2 question, not answered here)

**Gate:** how many genes are on the final list, and how many of those were located in the reference annotation.

---

## Steps

### Step 0 — Set up

```bash
cd /users/smaslini/vamp/binf6350_001
mkdir -p results/b1
G=bat_spleens/Desmodus_rotundus_references/GCF_022682495.2_HLdesRot8A.1_genomic.gtf
```

> Comment: the handout's accession (`GCF_022628495.2`) has transposed digits. The real file is `GCF_022682495.2`.
> Logged as a bonus-point candidate in the A3 notebook.

### Step 1 — Check how genes are named in the annotation

```bash
awk -F'\t' '$3=="gene"' $G | head -2
awk -F'\t' '$3=="gene"' $G | wc -l                       # 30,007 genes
awk -F'\t' '$3=="gene"' $G | grep -c 'gene_id "LOC'       # 13,505 are LOC-only
```

> Comment: 45% of genes have no symbol, only a `LOC` number. We already met this in Stage A3:
> hemoglobin alpha/beta were `LOC112298250/52` and `LOC112317579/80`.
> So a search by gene symbol alone would undercount. We search the `description` field too.

### Step 2 — Keyword survey of the annotation (`codes/run_b1.sh`, part 1)

Search every gene line for MHC-related words: `histocompatib`, `HLA`, `beta-2-microglobulin`,
`TAP binding`, `CD74`, `class II transactivator`.
Output: `results/mhc_survey.tsv` (41 hits).

> Comment: a keyword hit is a lead, not a result. Six hits were false:
> `SLFN11`, `SLFN5`, `SLFNL1` (sc**hla**fen contains "hla"), `HHLA1`, `HHLA2` (contain "HLA"),
> `HM13` (_minor_ histocompatibility antigen, not MHC).

### Step 3 — Symbol check for genes whose descriptions lack the keywords (`codes/run_b1.sh`, part 2)

Search by symbol for `TAP1`, `TAP2`, `PSMB8`, `PSMB9`, `HLA-D*`, `MR1`, and by description for DP/DO chains.
Output: `results/mhc_symbol_check.tsv`.

> Comment: TAP1/TAP2 descriptions read "transporter 1/2, ATP binding cassette subfamily B member", so Step 2 missed them.
> Found: TAP1, TAP2, PSMB8, PSMB9 (all on NC_071397.1, the MHC region).
> DP and DO chains: **not found** by either search.
> Our first version of this search also hit IDH1/IDH2 because the pattern `DP(` matched "NA**DP(**+)". That pattern was removed.

### Step 4 — Decide the gene list (the defended part)

Rule: **genes needed for MHC antigen presentation** = class I molecules + class II molecules +
obligate accessory genes. Written to `results/gene_list_pairD.tsv` (one row per gene, with a note).

| Part                 | Genes                                            | On list |
| -------------------- | ------------------------------------------------ | ------- |
| Class I heavy chains | 23 annotated `LOC` loci (each locus = one gene)  | 23      |
| Class II             | DRA, DRB, DQA, DQB, DMA, DMB, DPA, DPB, DOA, DOB | 10      |
| Accessory            | B2M, CD74, CIITA, TAPBP, TAP1, TAP2              | 6       |
| **Total**            |                                                  | **39**  |

Excluded, with reasons:

- MR1-like (`LOC112296868`): MHC-like, but presents microbial metabolites, not peptides.
- TAPBPL: optional class I peptide editor, not obligate.
- PSMB8/PSMB9: immunoproteasome; defensible to add (then the list is 41/37). We keep them out to keep the list to presentation molecules + obligate partners.
- The false keyword hits from Steps 2–3.

> Comment: each class I locus counts as its own gene because Stage B2 counts reads per locus.
> The names (BoLA = cow, Patr = chimp, Gogo = gorilla, Saoe = monkey, DLA = dog) are NCBI's
> closest-hit labels from other species, **not** orthology.

### Step 5 — Locate every list gene and build the deliverable table (`codes/run_b1.sh`, part 3)

```bash
cp codes/b1_find_genes.py scripts/          # if not already there
cp results/gene_list_pairD.tsv results/b1/
bash codes/run_b1.sh     # or: sbatch codes/run_b1.sh
```

`b1_find_genes.py` matches each list entry by gene ID/symbol, or by description regex, then reports
coordinates, strand, biotype, all transcript IDs, and transcript lengths (sum of exons).

Outputs (in `results/b1/` on the cluster; copy to `results/` here):

- `b1_gene_table.tsv`: one row per located gene → **the B1 deliverable table**
- `b1_transcripts.tsv`: every transcript with length (nt) and exon count
- `b1_not_found.tsv`: list genes with no match (expected: HLA-DPA1, HLA-DPB1, HLA-DOA, HLA-DOB)

### Step 6 — Gate

| Question                        | Answer |
| ------------------------------- | ------ |
| Genes on final list             | **39** |
| Located in reference annotation | **35** |

Check: the script's printed `===== GATE B1 =====` block must show the same two numbers.

### Step 7 — Write the inclusion and missing-gene note

See `explanation.md`. Cite the ≥5 papers there.

---

## Environment

| Item       | Value                                                         |
| ---------- | ------------------------------------------------------------- |
| Cluster    | gal-i9, partition Centaurus                                   |
| Module     | anaconda3/2023.09 (Python 3.11.7)                             |
| Tools      | awk, grep, sort (GNU coreutils); Python standard library only |
| Annotation | GCF_022682495.2_HLdesRot8A.1_genomic.gtf (505 MB)             |

## Checklist

- [x] Annotation naming checked (LOC problem)
- [x] Keyword survey + false hits removed
- [x] Symbol check (TAP1/2 found; DP/DO not annotated)
- [x] Gene list decided and written
- [ ] `run_b1.sh` run on the cluster; outputs copied to `results/`
- [ ] Gate submitted (39 / 35); result recorded here
- [ ] ≥5 papers read and cited in `explanation.md`
