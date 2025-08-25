#!/bin/bash

#SBATCH --job-name=gatk_variant_merge       # Job name
#SBATCH --output=gatk_variant_filter_%A_%a.out  # Standard output file
#SBATCH --error=gatk_variant_filter_%A_%a.err   # Standard error file
#SBATCH --time=4:00:00                         # Max runtime (adjust as needed)
#SBATCH --mem=40G                               # Memory per job
#SBATCH --cpus-per-task=2                       # Number of CPUs per job
#SBATCH -A m2_jgu-evolpest                      # Allocation account
#SBATCH -p smp                                  # Partition/queue (adjust for your HPC cluster)
#SBATCH --dependency=afterok:16406482		# Run this after the filtering is completed

# Load necessary modules
module load bio/GATK/4.6.0.0-GCCcore-13.3.0-Java-17.0.1
module load bio/tabixpp/1.1.0-GCC-11.2.0   # To use bgzip and tabix for compression and indexing

# Define the list of chromosomes you want to process in parallel
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18")   # Adjust according to your data


# Perform merging
echo "Merging all filtered VCF files..."
gatk MergeVcfs \
  $(for chr in "${chromosomes[@]}"; do echo "-I ${chr}.flt.pass.vcf.gz"; done) \
  -O Chr_all.flt.pass.vcf.gz

echo "Merging complete."
