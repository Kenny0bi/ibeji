#!/usr/bin/env bash
# Build the genotype files for ibeji from 1000 Genomes phase 3 (GRCh37).
#
# 1. Allele frequencies separately in EUR358 and YRI87 (biallelic autosomal SNVs).
# 2. Union of SNVs with MAF >= 0.01 in either population, so the decomposition
#    has frequencies and LD in both populations for every SNP either model can use.
# 3. One PLINK 2 fileset for all 445 matched individuals on that union.
# 4. Genotype PCs within each population from LD-pruned common SNPs.
set -euo pipefail
cd "$(dirname "$0")/.."

P=tools/plink2
RAW=data/raw
S=data/processed/samples
G=data/processed/geno
mkdir -p "$G"

BASE=(--pgen "$RAW/all_phase3.pgen" --pvar "$RAW/all_phase3_noannot.pvar.zst"
      --psam "$RAW/phase3_corrected.psam" --autosome --snps-only just-acgt
      --max-alleles 2 --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 10 truncate
      --rm-dup exclude-all --threads 7 --memory 9000)

for pop in EUR358 YRI87; do
  [ -s "$G/freq_$pop.afreq" ] || "$P" "${BASE[@]}" --keep "$S/$pop.txt" --freq --out "$G/freq_$pop"
done

# Streamed union: the .afreq files hold ~78M rows each, too large to load in memory.
# Columns: #CHROM ID REF ALT PROVISIONAL_REF? ALT_FREQS OBS_CT (ALT_FREQS located by header name).
for pop in EUR358 YRI87; do
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="ALT_FREQS") c=i; next}
              {f=$c; m=(f>0.5)?1-f:f; if(m>=0.01) print $2}' "$G/freq_$pop.afreq" > "$G/maf01_$pop.txt"
done
LC_ALL=C sort "$G/maf01_EUR358.txt" > "$G/maf01_EUR358.sorted"
LC_ALL=C sort "$G/maf01_YRI87.txt" > "$G/maf01_YRI87.sorted"
LC_ALL=C comm -12 "$G/maf01_EUR358.sorted" "$G/maf01_YRI87.sorted" | wc -l | awk '{print "both populations MAF>=0.01: " $1}'
LC_ALL=C sort -m -u "$G/maf01_EUR358.sorted" "$G/maf01_YRI87.sorted" > "$G/union_maf01.txt"
wc -l "$G/maf01_EUR358.txt" "$G/maf01_YRI87.txt" "$G/union_maf01.txt"
rm -f "$G/maf01_EUR358.txt" "$G/maf01_YRI87.txt"

"$P" "${BASE[@]}" --keep "$S/ALL.txt" --extract "$G/union_maf01.txt" \
  --make-pgen --out "$G/geuv445"

for pop in EUR358 YRI87; do
  "$P" --pfile "$G/geuv445" --keep "$S/$pop.txt" --maf 0.05 \
    --indep-pairwise 1000kb 1 0.2 --threads 7 --out "$G/prune_$pop"
  "$P" --pfile "$G/geuv445" --keep "$S/$pop.txt" --extract "$G/prune_$pop.prune.in" \
    --pca 10 --threads 7 --out "$G/pca_$pop"
done

ls -l "$G"
