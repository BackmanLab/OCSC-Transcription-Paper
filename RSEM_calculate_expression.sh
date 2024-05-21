#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --mail-user=janefrederick2021@u.northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH -t 24:00:00
#SBATCH --job-name="RSEM expression"

# Load modules
module load rsem/1.3.0

# Set working directory
cd /projects/b1042/BackmanLab/Jane/Chemo_OCSC/seq_raw_output

# Commands to execute
genomefold=/projects/b1042/BackmanLab/Jane/ref_hg38ensembl
for file in *.toTranscriptome.out.bam;
do rsem-calculate-expression -p 10 --bam $file "$genomefold/hg38_ensembl" ${file%.*.*.*};
done
