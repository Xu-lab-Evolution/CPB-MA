#!/bin/bash

#SBATCH -J gvcfRepeat                  # Job name
#SBATCH -o gvcfRepeat.%j.out           # Specify stdout output file (%j expands to jobId)
#SBATCH -p smp                   # Queue name 'smp' or 'parallel' on Mogon II
#SBATCH -n 12                     # Total number of tasks, here explicitly 12
#SBATCH -c 1                    # Total number of cores for the single task
#SBATCH --mem 50G                # The default is 300M memory per job. You'll likely have to adapt this to your needs
#SBATCH -t 36:00:00              # Run time (hh:mm:ss)

#SBATCH -A  m2_jgu-evolpest     # Specify allocation to charge against

#SBATCH --mail-type=END
#SBATCH --mail-user=shuqinxu@uni-mainz.de

#SBATCH --array=1-30

## Load all necessary modules if needed
module load bio/GATK/4.6.0.0-GCCcore-13.3.0-Java-17.0.1 

#param_store=chrom_BAM_list.txt #list with read files (change --array above according to the number of lines in the file passed here)
param_store=bam.list                 # list with read files
BAMfile=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $param_store)  # Get the BAM file for the current task


reference=/lustre/project/m2_jgu-evolpest/somaia/full_mapping/Leptinotarsa_decemlineata.fa
ID=$(basename "$BAMfile" | sed 's/_sorted_RG_noDup.bam//')  # Get the base name for output

## Launch the executable
gatk HaplotypeCaller -I $BAMfile -R $reference -O ${ID}.g.vcf.gz --linked-de-bruijn-graph true --native-pair-hmm-threads 12 --standard-min-confidence-threshold-for-calling 40  -ERC BP_RESOLUTION -G AS_StandardAnnotation -G StandardAnnotation
