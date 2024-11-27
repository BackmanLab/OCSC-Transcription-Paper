#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --mail-user=janefrederick2021@u.northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH -t 24:00:00
#SBATCH --job-name="RSEM expression"

# Load the RSEM module version 1.3.0
module load rsem/1.3.0

# Change to the working directory where the raw sequence output files are located
cd /projects/b1042/BackmanLab/Jane/Chemo_OCSC/seq_raw_output

# Define the path to the genome folder
genomefold=/projects/b1042/BackmanLab/Jane/ref_hg38ensembl

# Loop through each BAM file and calculate expression using RSEM
for file in *.toTranscriptome.out.bam;
do 
    # Run RSEM calculate expression with 10 threads, specifying the BAM file and genome reference
    rsem-calculate-expression -p 10 --bam $file "$genomefold/hg38_ensembl" ${file%.*.*.*}
done
