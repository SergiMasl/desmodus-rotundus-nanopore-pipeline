#!/bin/bash
#SBATCH --job-name=stage_a3
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32GB
#SBATCH --partition=Centaurus
#SBATCH --output=/users/smaslini/vamp/binf6350_001/logs/stage_a3_%j.out

# Stage A3: identify the 500-650 bp spike and quantify its share.
set -euo pipefail
P=/users/smaslini/vamp/binf6350_001
R=$P/bat_spleens/Desmodus_rotundus_references
REF=$R/GCF_022628495.2_HLdesRot8A.1_rna.fna
SPIKE=$P/bat_spleens/spike_reads_500_650bp.fasta
OWN=$P/results/a2/run1_bin_500-650.fastq.gz
O=$P/results/a3
T=8
mkdir -p $O/db

module load blast/2.15.0+
module load minimap2/2.28
module load samtools/1.19
module load seqkit/2.9.0
echo "===== VERSIONS ====="
blastn -version | head -1
echo "minimap2 $(minimap2 --version)"
samtools --version | head -1
echo "seqkit $(seqkit version)"

echo "===== 0. TRANSCRIPT -> GENE TABLE ====="
grep '^>' $REF \
 | sed -E 's/^>([^ ]+) .*\(([^()]+)\)[^()]*$/\1\t\2/' \
 | awk -F'\t' 'NF==2' > $O/tx2gene.tsv
echo "transcripts with a gene symbol: $(wc -l < $O/tx2gene.tsv)"

echo "===== 1. ALPHABET FIX (U -> T) ====="
sed '/^>/!y/Uu/Tt/' $SPIKE > $O/spike_dna.fa
seqkit head -n 200 $OWN | seqkit fq2fa \
 | sed '/^>/!y/Uu/Tt/' > $O/own200_dna.fa
for f in $SPIKE $O/spike_dna.fa $O/own200_dna.fa; do
  nU=$(grep -v '^>' $f | tr -cd 'Uu' | wc -c)
  nT=$(grep -v '^>' $f | tr -cd 'Tt' | wc -c)
  echo "$(basename $f): U=$nU T=$nT"
done

echo "===== 2. BLAST DB ====="
makeblastdb -in $REF -dbtype nucl -parse_seqids \
 -out $O/db/desrot_rna > $O/db/makeblastdb.log

FMT="6 qseqid sseqid pident length qlen slen qcovs evalue bitscore stitle"
blast_top () {
  blastn -query $1 -db $O/db/desrot_rna -outfmt "$FMT" \
   -max_target_seqs 5 -max_hsps 1 -evalue 1e-10 \
   -num_threads $T > $2
  awk -F'\t' '!s[$1]++' $2 > ${2%.tsv}_top.tsv
}

echo "===== 3. BLAST: RAW U reads (the bug) ====="
blast_top $SPIKE $O/blast_spike_RAW.tsv \
 || echo "blastn FAILED on raw U reads (that is also evidence)"
touch $O/blast_spike_RAW_top.tsv
echo "queries with a hit: $(wc -l < $O/blast_spike_RAW_top.tsv) / 30"

echo "===== 4. BLAST: fixed spike reads (30) ====="
blast_top $O/spike_dna.fa $O/blast_spike.tsv
echo "queries with a hit: $(wc -l < $O/blast_spike_top.tsv) / 30"
echo "top-hit gene per query (count gene):"
cut -f10 $O/blast_spike_top.tsv \
 | sed -E 's/.*\(([^()]+)\)[^()]*$/\1/' | sort | uniq -c | sort -rn
echo "best 5 hits (qseqid sseqid pident alnlen qlen slen qcovs evalue bits title):"
sort -t$'\t' -k9,9gr $O/blast_spike_top.tsv | head -5

echo "===== 5. BLAST: own 200 reads from run1 500-650 bin ====="
blast_top $O/own200_dna.fa $O/blast_own200.tsv
echo "queries with a hit: $(wc -l < $O/blast_own200_top.tsv) / 200"
cut -f10 $O/blast_own200_top.tsv \
 | sed -E 's/.*\(([^()]+)\)[^()]*$/\1/' | sort | uniq -c | sort -rn | head -10

GENE=$(cut -f10 $O/blast_spike_top.tsv \
 | sed -E 's/.*\(([^()]+)\)[^()]*$/\1/' | sort | uniq -c \
 | sort -rn | head -1 | awk '{print $2}')
echo ">>> most common spike gene: $GENE"

echo "===== 6. MAP ALL PASSED READS (minimap2 -> transcripts) ====="
for RUN in run1 run2; do
  if [ $RUN = run1 ]; then D=$P/T23_66_Diyo_RNA_spleen
  else D=$P/T23_66_Diyo_RNA_spleen_2ndRun; fi
  FQ=$(find $D -path '*fastq_pass*' -name '*.fastq.gz' | sort)

  zcat $FQ | awk 'NR%4==2{n++; b+=length($0)} END{print n"\t"b}' \
   > $O/${RUN}_passed_totals.tsv
  read NR NB < $O/${RUN}_passed_totals.tsv

  zcat $FQ | awk 'NR%4==2{gsub(/[Uu]/,"T")}1' \
   | minimap2 -ax map-ont -t $T --secondary=no $REF - \
     2> $O/${RUN}_minimap2.log \
   | samtools view -F 0x904 - \
   | awk -F'\t' -v OFS='\t' '{print $1,$3,length($10)}' \
   > $O/${RUN}_primary.tsv

  awk -F'\t' -v OFS='\t' 'NR==FNR{g[$1]=$2; next}
   {k=($2 in g)?g[$2]:$2; n[k]++; b[k]+=$3}
   END{for(k in n) print k,n[k],b[k]}' \
   $O/tx2gene.tsv $O/${RUN}_primary.tsv \
   | sort -t$'\t' -k2,2nr > $O/${RUN}_gene_counts.tsv

  NM=$(wc -l < $O/${RUN}_primary.tsv)
  echo "----- $RUN -----"
  echo "passed reads (fastq_pass): $NR   passed bases: $NB"
  echo "mapped (primary): $NM  = $(awk -v a=$NM -v b=$NR 'BEGIN{printf "%.2f", 100*a/b}')%"
  echo "top 10 genes: gene  reads  %passed_reads  bases  %passed_bases"
  head -10 $O/${RUN}_gene_counts.tsv | awk -F'\t' -v N=$NR -v B=$NB \
   '{printf "  %-14s %9d %6.2f%% %12d %6.2f%%\n",$1,$2,100*$2/N,$3,100*$3/B}'
  echo ">>> $GENE:"
  awk -F'\t' -v G="$GENE" -v N=$NR -v B=$NB '$1==G{
   printf "  reads %d = %.2f%% of passed reads (rounded %d%%)\n", $2,100*$2/N,int(100*$2/N+0.5);
   printf "  bases %d = %.2f%% of passed bases\n", $3,100*$3/B}' \
   $O/${RUN}_gene_counts.tsv
done

echo "Done."