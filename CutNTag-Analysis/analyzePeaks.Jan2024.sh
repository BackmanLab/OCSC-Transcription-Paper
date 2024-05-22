#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 8:00:00
#SBATCH --mem=96GB
#SBATCH --chdir=/projects/b1025/etb/MateiCutNTag2/NFcoreOut/
#SBATCH -o "%x.o%j.log"
#SBATCH --nodes=1
#SBATCH -n 12
#SBATCH --job-name=CutAndTag.Matei.Plots
#SBATCH --mail-user=ebartom@northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL


module load R/3.2.2
module load python/anaconda3
module load deeptools

# Collect a list of canonical transcripts from the UCSC genome browser.
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/knownCanonical.txt.gz
gunzip knownCanonical.txt.gz
sort -k 5 knownCanonical.txt  >  knownCanonical.geneSorted.txt
join -1 5 -2 1 knownCanonical.geneSorted.txt knownGene.geneSorted.txt | awk '{printf "%s\t%d\t%d\t%s\t%d\t%s\n",$2,$3,$4,$1,1,$8}' > knownCanonical.withStrand.bed

# Make tracks for CutNTag data
for bam in 02_alignment/bowtie2/target/markdup/*markdup.sorted.bam
do
     sample=${bam%%.target.markdup.sorted.bam}
     sample=$(basename $sample)
     echo $bam
     echo $sample
     bamCoverage -bl hg38.nonCanonicalChr.genome.blacklist.bed -b $bam -of bigwig --normalizeUsing CPM -o $sample.bw
done

# Get Canonical gene coordinates for hg38 and make bed files with strand.
sort -k 5 knownCanonical.txt  >  knownCanonical.enstSorted.txt
join -1 5 -2 1 knownCanonical.enstSorted.txt knownGene.geneSorted.txt | awk '{printf "%s\t%d\t%d\t%s.%s\t%d\t%s\n",$2,$3,$4,$6,$1,1,$8}' > knownCanonical.withStrand.bed

# # Pull out up-regulated and down-regulated genes.
ln -s /projects/b1025/etb/Matei.RNAseq2/cetoOut/newFastq/analysis/*geneList.txt .
grep -f ALDH.PlusVsMinus.fixed.htseq.01.dn.geneList.txt knownCanonical.withStrand.bed > ALDH.PlusVsMinus.htseq.01.dn.genes.bed
grep -f ALDH.PlusVsMinus.fixed.htseq.01.up.geneList.txt knownCanonical.withStrand.bed > ALDH.PlusVsMinus.htseq.01.up.genes.bed

# Set up log2 ratio bigwigs comparing epigenetic signature in ALDH plus vs ALDH minus cells.
bigwigCompare --bigwig2 K27ac.A_minus_R1.bw --bigwig1 K27ac.A_plus_R1.bw \
 	      -o K27ac.log2PlusOverMinus.R1.bw
bigwigCompare --bigwig2 K27ac.A_minus_R2.bw --bigwig1 K27ac.A_plus_R2.bw \
 	      -o K27ac.log2PlusOverMinus.R2.bw
bigwigCompare --bigwig2 K27me3.A_minus_R1.bw --bigwig1 K27me3.A_plus_R1.bw \
 	      -o K27me3.log2PlusOverMinus.R1.bw
bigwigCompare --bigwig2 K27me3.A_minus_R2.bw --bigwig1 K27me3.A_plus_R2.bw \
 	      -o K27me3.log2PlusOverMinus.R2.bw
bigwigCompare --bigwig2 K4me3.A_minus_R1.bw --bigwig1 K4me3.A_plus_R1.bw \
 	      -o K4me3.log2PlusOverMinus.R1.bw
bigwigCompare --bigwig2 K4me3.A_minus_R2.bw --bigwig1 K4me3.A_plus_R2.bw \
 	      -o K4me3.log2PlusOverMinus.R2.bw


# Plot the log2 ratio for each set of genes using deeptools computeMatrix.
computeMatrix reference-point --referencePoint TSS \
  	      -S K27me3.log2PlusOverMinus.R1.bw \
  	      K4me3.log2PlusOverMinus.R1.bw \
  	      K27ac.log2PlusOverMinus.R1.bw \
  	      -R ALDH.PlusVsMinus.htseq.01.up.genes.bed ALDH.PlusVsMinus.htseq.01.dn.genes.bed \
  	      -b 2000 -a 2000 -o log2K27me3K4me3K27ac.R1.deGenes.tss.2kb.R1.computeMatrix

# First plot the epigenetic marks for up genes and down genes separately.
 plotHeatmap --matrixFile log2K27me3K4me3K27ac.R1.deGenes.tss.2kb.R1.computeMatrix \
  	    --outFileName log2K27me3K4me3K27ac.deGenes.tss.2kb.R1.computeMatrix.pdf \
  	    --outFileSortedRegions log2K27me3K4me3K27ac.deGenes.tss.2kb.R1.computeMatrix.sorted.txt \
  	    --samplesLabel log2.K27me3.R1 log2.K4me3.R1 log2.K27ac.R1  \
  	    --regionsLabel up.ALDH.PlusVsMinus down.ALDH.PlusVsMinus \
  	    >& plotHeatmap.logK27me3K4me3K27ac.deGenes.tss.2kb.R1.log

# Then combine all of the DE genes and use kmeans to divide the TSS's into clusters of similarly regulated genes.
 plotHeatmap --matrixFile log2K27me3K4me3K27ac.R1.deGenes.tss.2kb.R1.computeMatrix \
  	    --outFileName log2K27me3K4me3K27ac.deGenes.tss.2kb.R1.km6b.computeMatrix.pdf \
  	    --outFileSortedRegions log2K27me3K4me3K27ac.deGenes.tss.2kb.R1.km6b.computeMatrix.sorted.txt \
  	    --samplesLabel log2.K27me3.R1 log2.K4me3.R1 log2.K27ac.R1  \
 	    --kmeans 6 \
  	    >& plotHeatmap.logK27me3K4me3K27ac.deGenes.tss.2kb.R1.km6b.log

# Plot the RNA seq data for each cluster of similarly regulated genes.
 for cluster in cluster_1 cluster_2 cluster_3 cluster_4 cluster_5 cluster_6
 do
     echo $cluster
     grep $cluster log2K27me3K4me3K27ac.deGenes.tss.2kb.R1.km6.computeMatrix.sorted.txt | awk '{print $4}' | perl -pe "s/\./\t/g" | awk '{print $1}' > log2K27me3K4me3K27ac.deGenes.tss.2kb.R1.km6.$cluster.geneList.txt
     Rscript makeBigHeatmapFromGeneList.R --geneListFile=log2K27me3K4me3K27ac.deGenes.tss.2kb.R1.km6.$cluster.geneList.txt --countFile=/projects/b1025/etb/Matei.RNAseq2/cetoOut/newFastq/analysis/htseq.normCounts.txt --clusterSamples=0 --clusterGenes=0 --assembly=hg38
 done
