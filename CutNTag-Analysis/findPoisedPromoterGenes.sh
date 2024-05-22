# Add genomic context to consensus peaks from NFcore output, using homer annotatePeaks.pl
# awk.bed files are in 03_peak_calling/05_consensus_peaks
for bed in *.awk.bed  
do
     echo $bed
     annotatePeaks.pl $bed hg38 > $bed.anno.txt   
done

module load bedtools/2.30.0

# Calculate the number of reads for each of the CutNTag samples for K4me3 and K27me3.
for mark in "K4me3" "K27me3"
do
      for cell in "A_plus" "A_minus"
      do
# 	 awk -F "\t" '$10 < 2000 && $10 > -2000 {print $16}' $mark.$cell.consensus.peaks.awk.bed.anno.txt | sort | uniq > $mark.$cell.consensus.peaks.awk.bed.anno.geneList.txt
 	 awk -F "\t" '$10 < 2000 && $10 > -2000 && length($16) > 0 {printf "%s\t%s\t%s\t%s\n",$2,$3,$4,$16}' $mark.$cell.consensus.peaks.awk.bed.anno.txt | sort | uniq > $mark.$cell.consensus.peaks.promoter.bed
 	 printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Chr" "Start" "Stop" "K27me3.Minus.R1" "K27me3.Minus.R2" "K27me3.Plus.R1" "K27me3.Plus.R2" "K4me3.Minus.R1" "K4me3.Minus.R2" "K4me3.Plus.R1" "K4me3.Plus.R2" > $mark.$cell.consensus.peaks.promoter.readCounts.table.txt
     bedtools multicov -bed $mark.$cell.consensus.peaks.promoter.bed -bams ../../02_alignment/bowtie2/target/K27me3.A_minus_R1.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K27me3.A_minus_R2.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K27me3.A_plus_R1.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K27me3.A_plus_R2.target.dedup.sorted.bam  ../../02_alignment/bowtie2/target/K4me3.A_minus_R1.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K4me3.A_minus_R2.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K4me3.A_plus_R1.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K4me3.A_plus_R2.target.dedup.sorted.bam >> $mark.$cell.consensus.peaks.promoter.readCounts.table.txt
      done
done

# The problem with the .peaks.promoter.bed files is that they are peak centered rather than gene centered, so different genomic regions are analyzed depending on the peak set that we started with.  For a more unbiased view, I think it makes sense to start with the gene lists and go from there to promoter regions to gene counts.

# Instead, pull out the lists of genes near the consensus peaks for each mark in each of the samples.
for cell in "A_plus" "A_minus"
do
     cat K4me3.$cell.consensus.peaks.awk.bed.anno.geneList.txt K27me3.$cell.consensus.peaks.awk.bed.anno.geneList.txt | sort | uniq -c | awk '$1 > 1' > poisedPromoter.$cell.geneList.txt
      awk '$2 !~ "" {print $2}' poisedPromoter.$cell.geneList.txt > poisedPromoter.$cell.geneList.geneNames.txt 
done

 wc -l poisedPromoter.*.geneList.txt

# # Used https://bioinformatics.psb.ugent.be/cgi-bin/liste/Venn/calculate_venn.htpl to find the venn diagram comparing the gene lists.

printf "%s:%s-%s\t%s\n" "Chr:Start-Stop" "Gene" > allGenes.coord.txt
cat *.anno.txt | sort | uniq | awk '{printf "%s:%s-%s\t%s\n",$2,$3,$4,$16}' | sort > allGenes.coord.txt


module load bedtools/2.30.0

for bed in *.awk.bed 
do
     echo $bed
     cut -f 1,2,3 $bed > $bed.simplified.bed
     printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Chr" "Start" "Stop" "K27me3.Minus.R1" "K27me3.Minus.R2" "K27me3.Plus.R1" "K27me3.Plus.R2" "K4me3.Minus.R1" "K4me3.Minus.R2" "K4me3.Plus.R1" "K4me3.Plus.R2" > $bed.readCounts.table.txt
     bedtools multicov -bed $bed.simplified.bed -bams ../../02_alignment/bowtie2/target/K27me3.A_minus_R1.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K27me3.A_minus_R2.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K27me3.A_plus_R1.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K27me3.A_plus_R2.target.dedup.sorted.bam  ../../02_alignment/bowtie2/target/K4me3.A_minus_R1.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K4me3.A_minus_R2.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K4me3.A_plus_R1.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K4me3.A_plus_R2.target.dedup.sorted.bam >> $bed.readCounts.table.txt

     awk '{printf "%s:%s-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' $bed.readCounts.table.txt | sort > $bed.readCounts.table.id.txt
     join --header -j 1 allGenes.coord.txt $bed.readCounts.table.id.txt > $bed.readCounts.labeled.table.txt
done

ls poisedPromoter.A_*.geneList.geneNames.txt

# Use reference annotations to build a bed file for a canonical TSS for all genes
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/knownCanonical.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/knownToRefSeq.txt.gz
gunzip knownCanonical.txt.gz
gunzip knownToRefSeq.txt.gz
sort -k 2 knownToRefSeq.txt > knownToRefSeq.RSsorted.txt
sort -k 1 knownToRefSeq.txt > knownToRefSeq.ENSTsorted.txt
sort -k 5 knownCanonical.txt  | awk -F "\t" '$1 !~ "_"' > knownCanonical.ENSTsorted.txt
#awk '{print $5}' knownCanonical.ENSTsorted.txt > knownCanonical.ENST.txt
join -1 5 -2 1 knownCanonical.ENSTsorted.txt knownToRefSeq.ENSTsorted.txt | sort -k 7 | perl -pe "s/ /\t/g" > knownCanonical.ENST.ENSG.RS.txt
sort -k 2 refGene.txt > refGene.RSsorted.txt
join -1 7 -2 2 knownCanonical.ENST.ENSG.RS.txt refGene.RSsorted.txt | awk -F " " '$10 ~ "+" {printf "%s\t%d\t%d\t%s.%s.%s\t%d\t%s\n",$3,$4-2000,$4+2000,$19,$2,$1,100,$10}' | sort | uniq > knownCanonicalTSS.posStrand.bed
join -1 7 -2 2 knownCanonical.ENST.ENSG.RS.txt refGene.RSsorted.txt | awk -F " " '$10 ~ "-" {printf "%s\t%d\t%d\t%s.%s.%s\t%d\t%s\n",$3,$5-2000,$5+2000,$19,$2,$1,100,$10}' | sort | uniq > knownCanonicalTSS.negStrand.bed
cat knownCanonicalTSS.posStrand.bed knownCanonicalTSS.negStrand.bed | bedtools sort > knownCanonicalTSS.bed

# Count the number of reads mapping to each canonical TSS from the different epigenetic markets
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Chr" "Start" "Stop" "Name" "Score" "Strand" "K27me3.Minus.R1" "K27me3.Minus.R2" "K27me3.Plus.R1" "K27me3.Plus.R2" "K4me3.Minus.R1" "K4me3.Minus.R2" "K4me3.Plus.R1" "K4me3.Plus.R2" > knownCanonicalTSS.readCounts.table.txt
bedtools multicov -bed knownCanonicalTSS.bed -bams ../../02_alignment/bowtie2/target/K27me3.A_minus_R1.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K27me3.A_minus_R2.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K27me3.A_plus_R1.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K27me3.A_plus_R2.target.dedup.sorted.bam  ../../02_alignment/bowtie2/target/K4me3.A_minus_R1.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K4me3.A_minus_R2.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K4me3.A_plus_R1.target.dedup.sorted.bam ../../02_alignment/bowtie2/target/K4me3.A_plus_R2.target.dedup.sorted.bam >> knownCanonicalTSS.readCounts.table.txt

# Calculate the mean for the replicates.
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Chr" "Start" "Stop" "Name" "Score" "Strand" "K27me3.Minus.RepMean" "K27me3.Plus.RepMean" "K4me3.Minus.RepMean" "K4me3.Plus.RepMean" > knownCanonicalTSS.readCounts.table.repMeans.txt
awk '$1 !~ "Chr" {printf "%s\t%s\t%s\t%s\t%s\t%s\t%8.2f\t%8.2f\t%8.2f\t%8.2f\n",$1,$2,$3,$4,$5,$6,($7+$8)/2,($9+$10)/2,($11+$12)/2,($13+$14)/2}' knownCanonicalTSS.readCounts.table.txt >> knownCanonicalTSS.readCounts.table.repMeans.txt

# grep is misbehaving, and reading in all lines, not just matching ones, so not using grep -f as I usually would.
while IFS= read -r line; do     echo $line; grep -w $line refGene.txt; done < poisedPromoter.A_minus.geneList.geneNames.txt > poisedPromoter.A_minus.refGene.txt
while IFS= read -r line; do     echo $line; grep -w $line refGene.txt; done < poisedPromoter.A_plus.geneList.geneNames.txt > poisedPromoter.A_plus.refGene.txt

for l in A_minus A_plus
do
     echo $l
     # # Put together the right list, in bed file format
     awk -F "\t" '$3 !~ "_" {printf "%s\t%d\t%d\t%s\t1\t%s\t%s\n",$3,$5,$6,$13,$4,$2}' poisedPromoter.$l.refGene.txt | sort | uniq | sort -k 7 > poisedPromoter.$l.almostBed
     join -1 7 -2 2 poisedPromoter.$l.almostBed knownToRefSeq.RSsorted.txt | awk '{printf "%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\n",$2,$3,$4,$5,$6,$7,$8,$1}' | sort -k 7 > poisedPromoter.$l.ENST.RS.almostBed
     grep -w -f knownCanonical.ENST.txt poisedPromoter.$l.ENST.RS.almostBed | sort -k 4 | awk '{printf "%s\t%d\t%d\t%s.%s.%s\t%d\t%s\n",$1,$2,$3,$4,$7,$8,$5,$6}' > poisedPromoter.$l.bed
done
