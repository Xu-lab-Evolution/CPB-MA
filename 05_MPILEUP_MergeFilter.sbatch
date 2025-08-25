#!/usr/bin/env bash
#SBATCH -J mergeChr
#SBATCH -p smp
#SBATCH -A m2_jgu-evolpest
#SBATCH -c 4
#SBATCH --mem 40G
#SBATCH -t 8:00:00
#SBATCH -o mergeChr.%j.out
#SBATCH -e mergeChr.%j.err

set -euo pipefail

OUTDIR=${OUTDIR:-jointcall_chr2}
CHRLIST=${CHRLIST:-"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18"}
THREADS=${THREADS:-$SLURM_CPUS_PER_TASK}

module purge
ml bio/BCFtools/1.14-GCC-11.2.0
ml bio/tabixpp/1.1.0-GCC-11.2.0

TMP=$(mktemp)
for C in $CHRLIST; do
  f="$OUTDIR/${C}.flt.vcf.gz"
  [[ -s $f ]] || { echo "[merge] missing $f" >&2; exit 1; }
  echo "$f" >> "$TMP"
done
# merge vcf files from each chr
MERGED=$OUTDIR/joint.chr1-18.flt.vcf.gz
bcftools concat -f "$TMP" -Oz -o "$MERGED" --threads "$THREADS" -a
bcftools index -t --threads "$THREADS" "$MERGED"
rm -f "$TMP"
#keep variants that have only one alternative allele count
bcftools view --min-ac=1 --max-ac=1 --threads $THREADS --no-update -Oz $MERGED > $OUTDIR/joint.chr1-18.flt.AC1.vcf.gz
bcftools index -t --threads $THREADS -f $OUTDIR/joint.chr1-18.flt.AC1.vcf.gz
echo "[merge] merged VCF → $MERGED"
