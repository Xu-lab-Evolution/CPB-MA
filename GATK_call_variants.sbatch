#!/bin/bash

#SBATCH -J variant_processing        # Job name
#SBATCH -o variant_processing.%A_%a.out  # Output file
#SBATCH -p smp                       # Partition
#SBATCH -A m2_jgu-evolpest           # Account
#SBATCH -n 8			   
#SBATCH -c 1                          # Cores per task
#SBATCH --mem 40G                     # Memory
#SBATCH -t 48:00:00                   # Time
#SBATCH --array=1-19                 # Array for chromosomes (1-18 + 1 job for all contigs)
#SBATCH --mail-type=END
#SBATCH --mail-user=shuqing.xu@uni-mainz.de

# Load GATK module
module load bio/GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8

# Paths and file definitions
REFERENCE="/lustre/project/m2_jgu-evolpest/somaia/full_mapping/Leptinotarsa_decemlineata.fa"
GVCF_DIR="/lustre/project/m2_jgu-evolpest/xus/CheckCPB/01_Variantscalling/gvcfs"
GVCF_LIST="gvcf_list.txt"
DICT_FILE="/lustre/project/m2_jgu-evolpest/somaia/full_mapping/Leptinotarsa_decemlineata.dict"
CHROMOSOME_INTERVALS="chromosome_intervals.list"
CONTIG_INTERVALS="contig_intervals.list"

# Process .dict file and split into chromosome and contig intervals
awk '/@SQ/ { 
    split($2, sn, ":");  
    split($3, ln, ":");   
    if (sn[2] ~ /^chr/) { 
        # Write chromosome intervals
        print sn[2] ":1-" ln[2] >> "'$CHROMOSOME_INTERVALS'"
    } else if (sn[2] ~ /^utg/) { 
        # Write contig intervals
        print sn[2] ":1-" ln[2] >> "'$CONTIG_INTERVALS'"
    }
}' "$DICT_FILE"



# Step 0: Create sample map file (SAMPLE_ID and path to each GVCF file)
if [ ! -f "$GVCF_LIST" ]; then
  for file in "$GVCF_DIR"/*.g.vcf.gz; do
    sample_id=$(basename "$file" | sed 's/.g.vcf.gz//')
    echo -e "${sample_id}\t${file}" >> "$GVCF_LIST"
  done
fi



if [ "$SLURM_ARRAY_TASK_ID" -le 18 ]; then
    INTERVAL="chr$SLURM_ARRAY_TASK_ID"
else
    INTERVAL="contig_intervals.list"  # Combined entry for all contigs in INTERVAL_FILE
fi

# Step 1: GenomicsDBImport
gatk GenomicsDBImport \
    --genomicsdb-workspace-path GenomicsDB_${SLURM_ARRAY_TASK_ID} \
    --intervals "$INTERVAL" \
    --reader-threads 8 \
    --sample-name-map $GVCF_LIST

# Step 2: GenotypeGVCFs
gatk GenotypeGVCFs \
    -R "$REFERENCE" \
    -V gendb://GenomicsDB_${SLURM_ARRAY_TASK_ID} \
    -L "$INTERVAL" \
    -O "${INTERVAL}_output.vcf.gz"

# Step 3: VariantFiltration
gatk VariantFiltration \
    -R "$REFERENCE" \
    -V "${INTERVAL}_output.vcf.gz" \
    -L "$INTERVAL" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O "${INTERVAL}.flt.vcf.gz"

# Step 4: SelectVariants
gatk SelectVariants \
    -R "$REFERENCE" \
    -V "${INTERVAL}.flt.vcf.gz" \
    -O "${INTERVAL}.flt.pass.vcf.gz" \
    --exclude-filtered

# Final MergeVcfs for last task in array
if [ "$SLURM_ARRAY_TASK_ID" -eq 19 ]; then
    echo "Merging all filtered VCFs..."
    gatk MergeVcfs \
        $(for chr in $(seq 1 18); do echo "-I chr${chr}.flt.pass.vcf.gz"; done) \
        -I contigs.flt.pass.vcf.gz \
        -O All_filtered_variants.vcf.gz
fi

echo "Workflow complete for array task $SLURM_ARRAY_TASK_ID."
