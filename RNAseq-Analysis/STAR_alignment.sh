#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --mail-user=janefrederick2021@u.northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH -t 24:00:00
#SBATCH --job-name="STAR Alignment Test"

# Load the STAR module version 2.5.2
module load STAR/2.5.2

# Change to the working directory containing raw sequence files
cd /projects/b1042/BackmanLab/Jane/Chemo_OCSC/seq_raw

# Define the path to the genome directory
genomefold=/projects/b1042/BackmanLab/Jane/ref_hg38ensembl

# Define the output directory based on the current directory name
outfold="$(dirname "$(pwd)")/${PWD##*/}_output"

# Create the output directory if it doesn't exist
[ ! -d $outfold ] && mkdir $outfold

# Loop through each FASTQ file in the directory
for file in *.fastq.gz;
do
    # Print the file name
    echo "$file"
    
    # Run STAR alignment for each file
    STAR --runThreadN 16 --quantMode TranscriptomeSAM --genomeDir $genomefold --readFilesIn <(gunzip -c "${file}") --outFileNamePrefix "$outfold/${file%.*.*}"
done
