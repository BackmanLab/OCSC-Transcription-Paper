#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH --mem=3GB
#SBATCH --chdir=/projects/b1025/etb/Matei.CSCs.GEO
#SBATCH -o "%x.o%j.log"
#SBATCH --nodes=1
#SBATCH -n 4
#SBATCH --job-name=RNAseq.Matei
#SBATCH --mail-user=ebartom@northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL

#If there are any modules loaded, remove them.

module purge
module load singularity/latest
module load graphviz/2.40.1
module load nextflow/22.10.5

#use local genomics at /projects/genomicsshare/AWS_iGenomes, simply adding the genomeID below will do this

nextflow run nf-core/rnaseq \
	 --input /projects/b1025/etb/Matei.CSCs.GEO/sampleSheet.RNA.csv \
	 --outdir /projects/b1025/etb/Matei.CSCs.GEO/NFcoreOut \
	 --genome hg38 -profile nu_genomics \
	 -c nextflow.config \
	 -resume

