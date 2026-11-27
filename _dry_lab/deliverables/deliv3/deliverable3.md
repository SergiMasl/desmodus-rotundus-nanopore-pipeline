# Stage 3 Deliverable: The Anomaly


## 1. Molecule Identity: 

We found the molecule to be hemoglobin subunit beta, or HBB. Lesser similarities were found among different sequences of hemoglobin subunit alpha, albeit still at a near
90% coverage. The runs for Run 1 and Run 2 had the following results, matching with HBB on NCBI:

### Run 1 Top 5 Results:

| Species Name         | Accession             | Percent Identity(%)  | Coverage(%)          | E-Value              | Max Score            | Molecule Name                |
| -------------------- | --------------------- |--------------------- |--------------------- |--------------------- |--------------------- |------------------------------|
| Desmodus rotundus    | XM_024574662.3        | 95.32                | 99                   | 0.0                  | 941                  | Hemoglobin Subunit Beta      |
| Phyllostomus hastatus| XM_045839476.1        | 94.15                | 99                   | 0.0                  | 902                  | Hemoglobin Subunit Beta      |
| Artibeus jamaicensis | XM_037154031.2        | 93.96                | 99                   | 0.0                  | 893                  | Hemoglobin Subunit Beta      |
| Phyllostomus discolor| XM_028515665.2        | 93.65                | 99                   | 0.0                  | 885                  | Hemoglobin Subunit Beta      |
| Sturnira hondurensis | XM_037040040.1        | 93.74                | 98                   | 0.0                  | 878                  | Hemoglobin Subunit Beta      |


### Run 2 Top 5 Results:

| Species Name         | Accession             | Percent Identity(%)  | Coverage(%)          | E-Value              | Max Score            | Molecule Name                |
| -------------------- | --------------------- |--------------------- |--------------------- |--------------------- |--------------------- |------------------------------|
| Desmodus rotundus    | XM_024574662.3        | 96.19                | 100                  | 0.0                  | 981                  | Hemoglobin Subunit Beta      |
| Phyllostomus hastatus| XM_045839476.1        | 93.38                | 100                  | 0.0                  | 887                  | Hemoglobin Subunit Beta-like |
| Sturnira hondurensis | XM_037154031.2        | 93.32                | 99                   | 0.0                  | 876                  | Hemoglobin Subunit Beta      |
| Artibeus jamaicensis | XM_028515665.2        | 93.11                | 98                   | 0.0                  | 865                  | Hemoglobin Subunit Beta      |
| Phyllostomus discolor| XM_037040040.1        | 93.07                | 98                   | 0.0                  | 859                  | Hemoglobin Subunit Beta-like |


The best match was from the closely related and well-studied cousin of the white-winged vampire bat (Diaemus youngi), the common vampire bat (Desmodus rotundus). 
Based on phylogenetics and evolutionary relationships, matching with the closest relative is a good indication that the molecular transcript identity went accordingly.
All of the matches mapped back to hemoglobin subunit beta transcripts of species belonging to the New World leaf-nosed bat family, Phyllostomidae. In Run 2's results, 
the Hemoglobin Subunit Beta-like results in the top 2nd and top 5th results may be paralogs or loosely assigned from sequence similarity instead of experimental evidence. 
When running MegaBLAST for high sequence similarities, the core nucleotide database and reference RNA sequences were considered. Both had similar results, with the 
core nucleotide databases returning slightly higher values. Both databases had all predicted results, so the core nucleotide database was chosen due to its relevance and 
frequent updating. 


## 2. Quantification: Share of Reads and Share of Bases:

### Run 1 Share of Reads and Bases:

| Run    | Species Name           | Accession        | Mapped Reads   | Mapped Bases   | % Passed Reads   | % Passed Bases   | Molecule Name                            |
| ------ | ---------------------- | ---------------- | -------------- | -------------- | ---------------- | ---------------- | ---------------------------------------- |
| Run 1  | Desmodus rotundus      | XM_024574662.3   | 43,330         | 21,505,376     | 5.419            | 2.702            | Hemoglobin Subunit Beta (LOC112317579)   |
| Run 1  | Desmodus rotundus      | XM_024574664.3   | 49,849         | 24,790,054     | 6.235            | 3.114            | Hemoglobin Subunit Beta (LOC112317580)   |
| Run 1  | Artibeus jamaicensis   | XM_037157363.2   | 58,924         | 28,585,381     | 7.370            | 3.591            | Hemoglobin Subunit Alpha (LOC119058257)  |
| Run 1  | All globin (Total)     | -                | 152,103        | 74,880,811     | 19.024           | 9.407            | Beta + Alpha                             |


### Run 2 Share of Reads and Bases:

| Run    | Species Name           | Accession        | Mapped Reads   | Mapped Bases   | % Passed Reads   | % Passed Bases   | Molecule Name                            |
| ------ | ---------------------- | ---------------- | -------------- | -------------- | ---------------- | ---------------- | ---------------------------------------- |
| Run 2  | Desmodus rotundus      | XM_024574662.3   | 2,735          | 1,377,914      | 5.418            | 2.607            | Hemoglobin Subunit Beta (LOC112317579)   |
| Run 2  | Desmodus rotundus      | XM_024574664.3   | 3,185          | 1,595,650      | 6.309            | 3.019            | Hemoglobin Subunit Beta (LOC112317580)   |
| Run 2  | Artibeus jamaicensis   | XM_037157363.2   | 3,711          | 1,837,036      | 7.351            | 3.476            | Hemoglobin Subunit Alpha (LOC119058257)  |
| Run 2  | All globin (Total)     | -                | 9,631          | 4,810,600      | 19.078           | 9.103            | Beta + Alpha                             |


Hemoglobin Subunit B has a passed read total of ~11.7% in both Run 1 and Run 2. Hemoglobin Subunit Alpha appears to map to a different genus of Phyllostomid bat, Artibeus 
jamaicensis. The reads have respective HBA totals of 7.4%, totaling in both runs to around 19%. Greater variation can be seen in the percentage of passed bases, with the
total globin passed bases measuring to 9.4% in Run 1, and 9.1% in Run 2. When taking the Gate 3 quiz online, it appears that the correct answer of 19% accounts for the sum
of HBA and HBB, rather than just the ~11.7% measure from HBB. Although all results are predicted, the mapping to HBB across the board from closely-related species is a 
proper indication that the highly abundant molecule in the 500-650 bin is indeed HBB. 


## 3. Spleen Abundance, Why?

HBB, as well as HBA, are abundant in the spleen due to its blood-filtering functions. The spleen is predominantly composed of red pulp and white pulp, allowing for the removal
of damaged or abnormal erythrocytes and the production of lymphocytes. Since erythrocytes contain hemoglobin for oxygen transportation, lysed erythrocytes in the spleen allow 
for high amounts of debris containing hemoglobin subunits A and B to be located in the tissue. 
 

## 4. Alphabet Problem and RNA vs DNA:

The issue in the calling is with direct RNA reads being used instead of DNA reads, which are expected for BLAST, minimap2, and many reference databases. When a uracil
is called instead of an expected thymine from DNA, the uracil is treated as an unknown base or mismatch. Although I did not encounter the issue on my run, I would
expect a significant identity score drop in the case that the following fix was not applied: 

seqkit seq --rna2dna

*Module load seqkit must be applied for the conversion to be applied.*


## 5. Protocol Changes for More Sequencing: 

To sequence 29 more spleens, a lot of the pipeline processing done for the Diaemus youngi spleen sample, such as the QC and quantification portions, 
would need to be automated. Two runs were completed for one spleen sample, and in the case of a MinION technical error resulting in a needed rerun, 
a fraction of the 29 runs may have to be processed again, totaling around 30-40 runs in total for processing. Therefore, automating the process by 
flagging abnormalities, such as the large spike in the 500-650 bp bin, or potentially subdividing the bins for easier diagnostics may warrant a 
faster and more streamlined approach. Flagging the finer subdivisions if they were well above other bins, alongside automating a U to T conversion
would save time and reduce errors. 




