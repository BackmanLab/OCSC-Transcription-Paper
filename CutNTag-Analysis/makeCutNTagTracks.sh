#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH --mem=50GB
#SBATCH --chdir=/projects/b1025/etb/Matei.CSCs.GEO/NFcoreOut.CutNTag/02_alignment/bowtie2/target/markdup/
#SBATCH -o "%x.o%j.log"
#SBATCH --nodes=1
#SBATCH -n 12
#SBATCH --job-name=makeTracks
#SBATCH --mail-user=ebartom@northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module load deeptools/3.5.1

# Make CPM bigwigs with deepTools
for bam in *.markdup.sorted.bam
do
    sample=${bam%.target.markdup.sorted.bam}
    echo $sample 
    bamCoverage -b $bam -o $sample.bw -of bigwig --normalizeUsing CPM
done

