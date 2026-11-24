#!/bin/bash
#SBATCH --job-name=gate_b1
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8GB
#SBATCH --partition=Centaurus
#SBATCH --output=/users/smaslini/vamp/binf6350_001/logs/gate_b1_%j.out

# Gate B1: Define the family (Pair D: MHC class I and class II)
# Fast (< 1 min): can run as  bash run_b1.sh  or  sbatch run_b1.sh
set -euo pipefail

P=/users/smaslini/vamp/binf6350_001
R=$P/bat_spleens/Desmodus_rotundus_references
G=$R/GCF_022682495.2_HLdesRot8A.1_genomic.gtf
S=$P/scripts
O=$P/results/b1
mkdir -p $O

module load anaconda3/2023.09
python --version

# pull gene_id, location, description from GTF gene lines
gene_fields () {
  awk -F'\t' '{a=$9; g=a; sub(/.*gene_id "/,"",g); sub(/".*/,"",g);
    d=a; sub(/.*description "/,"",d); sub(/".*/,"",d);
    print g"\t"$1":"$4"-"$5"\t"d}'
}

echo "===== 1. KEYWORD SURVEY of the annotation ====="
K='histocompatib|HLA|beta-2-microglobulin|antigen peptide transporter'
K="$K|tapasin|TAP binding|CD74|class II transactivator"
awk -F'\t' '$3=="gene"' $G | grep -i -E "$K" | gene_fields \
 | sort -t$'\t' -k3,3 > $O/mhc_survey.tsv
echo "hits: $(wc -l < $O/mhc_survey.tsv)"
cut -f1,3 $O/mhc_survey.tsv

echo "===== 2. SYMBOL CHECK (genes whose description lacks the keywords) ====="
SYM='gene_id "(TAP1|TAP2|PSMB8|PSMB9|HLA-D[A-Z0-9]+|MR1)"'
DSC='description "[^"]*(DO alpha|DO beta|DP alpha|DP beta|transporter [12], ATP binding)'
awk -F'\t' '$3=="gene"' $G | grep -E "$SYM|$DSC" | gene_fields \
 > $O/mhc_symbol_check.tsv
cat $O/mhc_symbol_check.tsv

echo "===== 3. LOCATE THE FINAL GENE LIST ====="
python $S/b1_find_genes.py $G $O/gene_list_pairD.tsv $O/b1

echo "Done."
