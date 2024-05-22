# Baseline_GO_Analysis.R --------------------------------------------------
# Gene ontology analysis of OVCAR5 ALDH+ vs. ALDH- cells
# Created on 12-02-2023 by Jane Frederick, Paola Carrillo Gonzalez, and Rachel 
# Ye
# Last updated 5-20-2024 by Jane Frederick


# Import Libraries --------------------------------------------------------

# if (!requireNamespace("BiocManager", quietly = TRUE))     
#   install.packages("BiocManager")

# BiocManager::install("clusterProfiler")
library("clusterProfiler")
# BiocManager::install("enrichplot")
library("enrichplot")
# BiocManager::install("msigdbr")
library("msigdbr")
# BiocManager::install('EnhancedVolcano')
library("EnhancedVolcano")
# BiocManager::install("pheatmap")
library("pheatmap")
# BiocManager::install('apeglm')
library("apeglm")
# BiocManager::install("DESeq2")
library("DESeq2")
# BiocManager::install("tximport")
library("tximport")
# BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")
library("ggplot2")
library("RColorBrewer")


# Set Directories ---------------------------------------------------------

# Directory with all files
parentdir <- "~/Ovarian CSCs/Baseline"
# Path to folder with the .genes.results files output by RSEM
RSEMdir <- "~/Ovarian CSCs/Baseline/TPM Files"
# Path to folder to save results
savedir <- "~/Ovarian CSCs/Baseline/Plots"

# Check if folder to save files exists 
if (!file.exists(savedir)){
  # Create the folder if it does not exist
  dir.create(savedir)
}


# Set Global Variables ----------------------------------------------------

# Set limit for p-value to determine significance
pval_lim <- 0.05
# Set limit for adjusted p-value to determine significance
padj_lim <- 0.05


# Read In Data ------------------------------------------------------------

# Read in samples file containing information about data files
samples <- read.csv(paste(parentdir,"samples.csv", sep = "/"), header = TRUE)
# Set the CSC column as a factor
samples$CSC <- factor(samples$CSC)
# Change the row names of the samples list to be CSC and replicate
rownames(samples) <- paste(samples$CSC, samples$Replicate, sep="_")
# Create an array with the file paths to read in data
files <- paste0(paste(RSEMdir,samples$ID,sep="/"), "Aligned.genes.results")
# Add names to each path to match the row names of the samples list
names(files) <- paste(samples$CSC, samples$Replicate, sep="_")

# Initialize the RNA TPM data frame
read_TPM <- data.frame(gene_id=NA)
# Read data into a single data frame
for (data in names(files)) {
  # Read in data from .genes.results file
  sample_vals <- read.table(files[data],
                            col.names=c("gene_id","transcript_id(s)","length",
                                        "effective_length",
                                        paste0(data,"_count"),
                                        paste0(data,"_TPM"),"FPKM"))
  # Extract the TPM column and add to the TPM data frame
  read_TPM <- merge(read_TPM, 
                    sample_vals[ , c("gene_id", paste0(data,"_TPM"))], 
                    by='gene_id', all=T)
}

# Remove the TPM substring from column names
colnames(read_TPM) <- sub("_T.*", "", colnames(read_TPM))
# Get the list of gene names
genes <- read_TPM$gene_id
# Convert the ENSEMBL gene IDs to HGNC symbols
annots <- AnnotationDbi::select(org.Hs.eg.db, keys=genes, 
                 columns="SYMBOL", keytype="ENSEMBL")
# Add the HGNC symbols to the TPM dataframe
TPMresult <- merge(read_TPM, annots, by.x="gene_id", by.y="ENSEMBL")
# Export the TPM dataframe as a CSV file
write.csv(TPMresult, file=paste(savedir,"All_Samples_TPM.csv",sep="/"))


# Use DESeq for DGE Analysis ----------------------------------------------

# Use the tximport package to read in RSEM files
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
# Find genes with zero length or that are not expressed
zero_length_and_unexpressed <- (apply(txi.rsem$abundance, 1, max) == 0) &
  (apply(txi.rsem$length, 1, min) == 0)
# Remove the zero length and unexpressed genes from the dataset
txi.rsem$length <- txi.rsem$length[!zero_length_and_unexpressed,]
txi.rsem$abundance <- txi.rsem$abundance[!zero_length_and_unexpressed,]
txi.rsem$counts <- txi.rsem$counts[!zero_length_and_unexpressed,]
# Create a DESeq dataset from the tximport object
ddsTxi <- DESeqDataSetFromTximport(txi.rsem,colData = samples,design = ~CSC)
# Run DESeq for differential gene expression analysis
dds <- DESeq(ddsTxi)
# Define the correct contrast (comparison) for ALDH+ differential expression
baseline_contrast <- c("CSC","ALDHplus","ALDHminus")
# Extract the DESeq results using the correct contrast
res <- results(dds, contrast=baseline_contrast)
# Display the results
head(res)
# Sort the results by p-value
resOrdered <- res[order(res$pvalue),]
# Filter the results using the adjusted p-value
resfilt <- subset(resOrdered, padj < padj_lim)
# Convert the filtered results to a dataframe
res_out <- as.data.frame(resfilt)
# Create a dataframe to convert ENSEMBL gene IDs to HGNC symbols
sym_conv <- annots[match(rownames(res), annots$ENSEMBL),]
# Add the HGNC symbols to the results dataframe
res_out <- merge(res_out, sym_conv, by.x=0, by.y="ENSEMBL")
# Export the filtered results dataframe
write.csv(res_out, 
          file=paste0(savedir,"/",baseline_contrast[2],"_vs_",
                      baseline_contrast[3],"_DGE_adjpval",
                      substr(padj_lim,nchar(padj_lim[1])-1,nchar(padj_lim[1])),
                      ".csv"))
# Get a list of HGNC gene symbols
symbols <- sym_conv$SYMBOL
# Extract genes with HGNC symbols from the DESeq results
ressym <- res[!is.na(symbols),]
# Remove genes that do not have HGNC symbols from the symbols list
symbols <- symbols[!is.na(symbols)]
# Add the HGNC symbols as row names for the results
rownames(ressym) <- symbols
# Perform shifted logarithm transformation
ntd <- normTransform(dds)
# Perform variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)
# Create a dataframe with the name of the samples
df <- as.data.frame(colData(dds)[,"Group"])
# Add the column names from the DESeq object as row names
rownames(df) <- colnames(ntd)
# Rename column
colnames(df) <- "Group"


# PCA and Distance Plots --------------------------------------------------

# Perform principal components analysis on VST data
pcaData <- plotPCA(vsd, intgroup=c("CSC"), returnData=TRUE)
# Calculate the variance between groups
percentVar <- round(100 * attr(pcaData, "percentVar"))
# Save PCA plot as a PNG file
png(file=paste0(savedir,"/PCA.png"))
# Plot PCA data as a scatter plot
ggplot(pcaData, aes(PC1, PC2, color=CSC, shape=CSC)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

# Calculate Euclidean distance between samples using VST data
sampleDists <- dist(t(assay(vsd)))
# Convert distances into a matrix
sampleDistMatrix <- as.matrix(sampleDists)
# Remove column names
colnames(sampleDistMatrix) <- NULL
# Define color palette for distance heatmap
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# Save Euclidean distance heatmap as a PNG file
png(file=paste0(savedir,"/Distance_Heatmap.png"))
# Use pheatmap to plot Euclidean distances on both axes
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()


# Create Volcano Plot -----------------------------------------------------

# Find number of significantly upregulated genes
sigup = sum(ressym$pvalue<10e-6 & ressym$log2FoldChange>1 , na.rm=TRUE)
# Find number of upregulated genes that are not significant
nonsigup = sum(ressym$log2FoldChange>1 , na.rm=TRUE)-sigup
# Find number of significantly downregulated genes
sigdown = sum(ressym$pvalue<10e-6 & ressym$log2FoldChange<(-1) , na.rm=TRUE)
# Find number of downregulated genes that are not significant
nonsigdown = sum(ressym$log2FoldChange<(-1) , na.rm=TRUE)-sigdown
# Display the number of genes in each category of the volcano plot
print("Number of Genes:")
print(paste0("P-value<10e-6 & log2FC>1: ",sigup," genes"))
print(paste0("P-value>10e-6 & log2FC>1: ",nonsigup," genes"))
print(paste0("P-value<10e-6 & log2FC<-1: ",sigdown," genes"))
print(paste0("P-value>10e-6 & log2FC<-1: ",nonsigdown," genes"))

# Save volcano plot as a TIFF file
tiff(filename = paste0(savedir,"/Volcano_Plot.tiff"), width = 12, height = 10,
     units="in", res = 1000, compression = "lzw")
par(mar = c(3, 3, 0.5, 0.5))
# Use EnhancedVolcano package to create volcano plot
EnhancedVolcano(ressym,lab = rownames(ressym),x = 'log2FoldChange',y = 'pvalue',
                pCutoff = 10e-6,FCcutoff = 1,boxedLabels = TRUE,
                drawConnectors = TRUE,title = "ALDH+ vs. ALDH-",
                subtitle = "Differential Gene Expression",
                caption = paste0(nrow(ressym)," Total Genes"))
dev.off()


# Enrichment Gene Ontology Analysis ---------------------------------------

# Extract fold change values from the filtered results object for genes with 
# and adjusted p value <0.05
geneList_filt <- resfilt$log2FoldChange
# Add the ENSEMBL gene IDs to the fold change list
names(geneList_filt) <- rownames(resfilt)
# Sort the list by fold change
geneList_filt <- sort(geneList_filt, decreasing = TRUE)
# Perform enrichment analysis on list of genes for biological process GO terms
ego <- enrichGO(gene          = names(geneList_filt),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                keyType       = "ENSEMBL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
# Find semantically similar GO terms
ego <- pairwise_termsim(ego)
# Remove redundant or parent/child terms
ego <- simplify(ego, cutoff=0.5, by="p.adjust", select_fun=min)

# Save a bar plot showing the Q score of top enriched terms as a PNG file
png(filename = paste0(savedir,"/enrichGO_Qscore_Bar_Plot.png"))
# Create bar plot for Q scores of top terms
mutate(ego, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore",xlab=bquote(-Log[10]~"P"))
dev.off()

# Convert ENSEMBL gene ID to HGNC symbol in the GO enrichment object
egox <- setReadable(ego, 'org.Hs.eg.db', 'ENSEMBL')
# Save a circular network plot from the GO enrichment as a TIFF file
tiff(filename = paste0(savedir,"/enrichGO_Circular_Network_Plot.tiff"), 
     width = 9, height = 9, units="in", res = 1000, compression = "lzw")
par(mar = c(3, 3, 0.5, 0.5))
# Create circular network plot
cnetplot(egox, color.params = list(foldChange = geneList_filt), circular = TRUE,
         colorEdge = TRUE) 
dev.off()

# Filter the gene list to find upregulated genes
geneList_upreg <- geneList_filt[geneList_filt>0]
# Perform enrichment analysis on upregulated genes for BP GO terms
ego_upreg <- enrichGO(gene          = names(geneList_upreg),
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      keyType       = "ENSEMBL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
# Save a dot plot showing GO terms for upregulated genes as a TIFF file
tiff(filename = paste0(savedir,"/enrichGO_Upregulated_Dot_Plot.tiff"), 
     width = 8, height = 11, units="in", res = 1000, compression = "lzw")
par(mar = c(3, 3, 0.5, 0.5))
# Create dot plot for top upregulated enriched terms
dotplot(ego_upreg, showCategory=20)
dev.off()

# Filter the gene list to find downregulated genes
geneList_downreg <- geneList_filt[geneList_filt<0]
# Perform enrichment analysis on downregulated genes for BP GO terms
ego_downreg <- enrichGO(gene          = names(geneList_downreg),
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      keyType       = "ENSEMBL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
# Save a dot plot showing GO terms for downregulated genes as a TIFF file
tiff(filename = paste0(savedir,"/enrichGO_Downregulated_Dot_Plot.tiff"), 
     width = 8, height = 11, units="in", res = 1000, compression = "lzw")
par(mar = c(3, 3, 0.5, 0.5))
# Create dot plot for top downregulated enriched terms
dotplot(ego_downreg, showCategory=20)
dev.off()


# GO Gene Set Enrichment Analysis -----------------------------------------

# Extract fold change values from the results object
geneList_baseline <- res$log2FoldChange
# Add the ENSEMBL gene IDs to the fold change list
names(geneList_baseline) <- rownames(res)
# Sort the list by fold change
geneList_baseline <- sort(geneList_baseline, decreasing = TRUE)

# Gene Set Enrichment Analysis for biological process GO terms
gse <- gseGO(geneList=geneList_baseline, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")
# Use DOSE package
require(DOSE)
# Filter terms to only keep those with a positive enrichment score
posenrich <-  filter(gse, NES>0)
# Save a ridge plot showing positive GO enrichment as a TIFF file
tiff(filename = paste0(savedir,"/gseGO_Ridge_Plot_Positive_NES.tiff"), 
     width = 8, height = 10, units="in", res = 1000, compression = "lzw")
par(mar = c(3, 3, 0.5, 0.5))
# Create ridge plot for top enriched terms with positive enrichment
ridgeplot(posenrich, showCategory=15) + labs(x = "Enrichment Score")
dev.off()


# Gene Set Enrichment Analysis --------------------------------------------

# Get gene sets from the molecular signatures database
m_df <- msigdbr(species = "Homo sapiens") %>% 
  dplyr::select(gs_name, ensembl_gene)
# Runs GSEA on all gene sets
em <- GSEA(geneList_baseline, TERM2GENE = m_df)
# Export GSEA results for all gene sets
write.csv(as.data.frame(em), file=paste0(savedir,"/GSEA_List.csv"))

# Save GSEA plot for one gene set as a TIFF file
tiff(filename = paste0(savedir,"/BOQUEST_STEM_CELL_UP.tiff"), width = 7, 
     height = 7, units="in", res = 1000, compression = "lzw")
par(mar = c(3, 3, 0.5, 0.5))
# Create GSEA plot for specified gene set
gseaplot2(em, geneSetID = "BOQUEST_STEM_CELL_UP", title = "BOQUEST_STEM_CELL_UP")
dev.off()
# Save GSEA plot for one gene set as a TIFF file
tiff(filename = paste0(savedir,"/LIM_MAMMARY_STEM_CELL_UP.tiff"), width = 7, 
     height = 7, units="in", res = 1000, compression = "lzw")
par(mar = c(3, 3, 0.5, 0.5))
# Create GSEA plot for specified gene set
gseaplot2(em, geneSetID = "LIM_MAMMARY_STEM_CELL_UP", 
          title = "LIM_MAMMARY_STEM_CELL_UP")
dev.off()
