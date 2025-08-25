#!/usr/bin/env bash
#SBATCH -J jointChr
#SBATCH -p smp
#SBATCH -A m2_jgu-evolpest
#SBATCH -c 24
#SBATCH --mem 80G
#SBATCH -t 48:00:00
#SBATCH -o jointChr.%A_%a.out
#SBATCH -e jointChr.%A_%a.err
#SBATCH --array=1-18           # 18 chromosomes → 18 tasks

set -euo pipefail

# ---------- user-tunable paths ---------------------------------------------
REF=${REF:-/lustre/project/m2_jgu-evolpest/somaia/full_mapping/Leptinotarsa_decemlineata.fa}
PED=${PED:-pedigree.sub.ped}
BAMDIR=${BAMDIR:-/lustre/project/m2_jgu-evolpest/somaia/full_mapping}
OUTDIR=${OUTDIR:-jointcall_chr2}
THREADS=${THREADS:-$SLURM_CPUS_PER_TASK}
CHRLIST=${CHRLIST:-"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18"}
mkdir -p "$OUTDIR"

# ---------- modules ---------------------------------------------------------
module purge
ml bio/SAMtools/1.17-GCC-12.2.0
ml bio/BCFtools/1.14-GCC-11.2.0
ml bio/tabixpp/1.1.0-GCC-11.2.0

# ---------- build BAM list (once, guarded by flock) -------------------------
LOCK="$OUTDIR/lists.lock"
BAMLIST="$OUTDIR/bams.final.list"

(
 flock -w 300 9 || { echo "[lock] timeout" >&2; exit 1; }

 if [[ ! -s $BAMLIST || $BAMLIST -ot $PED ]]; then
   echo "[build] regenerating BAM list"
   awk 'NF>=5 && !/^#/ {print $2"\n"$3"\n"$4}' "$PED" | sort -u > "$OUTDIR/samples.tmp"

   > "$BAMLIST.tmp"
   while read -r sid; do
     bam="$BAMDIR/${sid}_sorted_RG_noDup.bam"
     [[ -s $bam && ( -s $bam.bai || -s ${bam%.bam}.bai ) ]] \
        || { echo "Missing BAM $sid" >&2; exit 1; }
     echo "$bam" >> "$BAMLIST.tmp"
   done < "$OUTDIR/samples.tmp"

   mv "$BAMLIST.tmp" "$BAMLIST"
   rm -f "$OUTDIR/samples.tmp"
 fi
) 9>"$LOCK"

# ---------- map array index -> chromosome -----------------------------------
read -ra CHR_ARR <<< "$CHRLIST"
CHR=${CHR_ARR[$((SLURM_ARRAY_TASK_ID-1))]}
[[ -n $CHR ]] || { echo "No chromosome for task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }

RAW_BCF=$OUTDIR/${CHR}.raw.bcf
FLT_VCF=$OUTDIR/${CHR}.flt.vcf.gz

echo "[$(date)] task=$SLURM_ARRAY_TASK_ID  chr=$CHR"

# ---------- mpileup → call → fill-tags → filter -----------------------------
bcftools mpileup -Ou --threads 4 -a "FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR" \
   -f "$REF" \
   -b "$BAMLIST" \
   -r "$CHR" \
   -Ou --threads "$THREADS" |
bcftools call -mv --threads "$THREADS" -Ob -o "$RAW_BCF"


# Convert BCF to VCF
bcftools view --threads 4 "$OUTPUT.bcf" > "$OUTPUT.vcf"

bcftools +fill-tags "$RAW_BCF" --threads "$THREADS" \
  -Oz -o "$FLT_VCF" -- -t AC,AN

bcftools index -t --threads "$THREADS" "$FLT_VCF"

# per-chr filter
bcftools view --threads "$THREADS" \
   -i 'QUAL>30 && INFO/DP>60 && INFO/MQ>30' \
   "$FLT_VCF" -Oz -o "${FLT_VCF%.vcf.gz}.vcf.tmp.gz" && \
mv "${FLT_VCF%.vcf.gz}.vcf.tmp.gz" "$FLT_VCF" && \
bcftools index -t --threads "$THREADS" "$FLT_VCF"

echo "[$(date)] finished $CHR"
