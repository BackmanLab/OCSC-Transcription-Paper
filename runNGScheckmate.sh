#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH --mem=50000
#SBATCH --chdir=/projects/b1025/etb/Matei.CSCs.GEO
#SBATCH -o "%x.o%j"
#SBATCH --job-name=Matei.CSCs.NCM
#SBATCH --nodes=1
#SBATCH -n 24

#If there are any modules loaded, remove them.
module purge

pattern="/projects/p20742/tools/bin/NGSCheckMate/SNP/SNP.pt"

# Setup paths for the analysis
export PATH=$PATH:/projects/p20742//tools/bin/
export NCM_HOME=/projects/p20742//tools/bin/NGSCheckMate/

module load python/anaconda
module load R/3.2.2

# Run NGSCheckmate on the fastq files listed in fastqList.txt
python $NCM_HOME/ncm_fastq.py -p 24 -l /projects/b1025/etb/Matei.CSCs.GEO/fastqList.dataset3.txt -pt $pattern -O /projects/b1025/etb/Matei.CSCs.GEO/NGSCheckMate/ >& fastqs.ncm.err.log

module load R/4.3.0

# Compare NCM profile to those from CCLE to validate OVCAR5 identity.
Rscript /projects/p20742/celID.testData/testAgainstCCLE_cm.R --queryDir=/projects/b1025/etb/Matei.CSCs.GEO/NGSCheckMate/
