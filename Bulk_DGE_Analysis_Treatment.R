# Combo_Treatment_GO_Analysis.R --------------------------------------------
# Gene ontology analysis of OVCAR5 ALDH+ vs. ALDH- cells treated with DMSO,
# cisplatin, EPZ-5676, or both cisplatin and EPZ-5676
# Created on 8-17-2023 by Jane Frederick
# Last updated 5-21-2024 by Jane Frederick


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
parentdir <- "~/Ovarian CSCs/Combo_Treatment"
# Path to folder with the .genes.results files output by RSEM
RSEMdir <- "~/Ovarian CSCs/Combo_Treatment/TPM Files"
# Path to folder to save results
savedir <- "~/Ovarian CSCs/Combo_Treatment/Plots"

# Check if folder to save files exists 
if (!file.exists(savedir)){
  # Create the folder if it does not exist
  dir.create(savedir)
}


# Set Global Variables ----------------------------------------------------

# Group to use as the reference for initial DGE analyses
deseq_ref <- "NonCSC_DMSO"
# Set limit for p-value to determine significance
pval_lim <- 0.1
# Set limit for adjusted p-value to determine significance
padj_lim <- 0.05


# Read In Data ------------------------------------------------------------

# Read in samples file containing information about data files
samples <- read.table(paste(parentdir,"samples.txt", sep = "/"), header = TRUE)
# Set the Group, CSC, and Treatment columns as factor columns
samples$Group <- factor(samples$Group)
samples$CSC <- factor(samples$CSC)
samples$Treatment <- factor(samples$Treatment)
# Change the row names of the samples list to be group and replicate
rownames(samples) <- paste(samples$Group, samples$Replicate, sep="_")
# Create an array with the file paths to read in data
files <- paste0(RSEMdir,"/",samples$ID, "Aligned.genes.results")
# Add names to each path to match the row names of the samples list
names(files) <- paste(samples$Group, samples$Replicate, sep="_")

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
ddsTXi <- DESeqDataSetFromTximport(txi.rsem,colData = samples,design = ~Group)
# Relevel the contrast to compare to the reference group
ddsTXi$Group <- relevel(ddsTXi$Group, ref = deseq_ref)
# Run DESeq for differential gene expression analysis
dds <- DESeq(ddsTXi)
# Extract the DESeq results
res <- results(dds)
# Remove genes without values
res <- na.omit(res)
# Display the results
head(res)
# Sort the results by p-value
resOrdered <- res[order(res$pvalue),]
# Filter the results using the adjusted p-value
resfilt <- subset(resOrdered, padj < padj_lim)
# Perform shifted logarithm transformation
ntd <- normTransform(dds)
# Perform variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)
# Create a dataframe with the name of the samples
df <- as.data.frame(colData(dds)[,c("CSC","Treatment")])
# Add the column names from the DESeq object as row names
rownames(df) <- colnames(ntd)


# PCA and Distance Plots --------------------------------------------------

# Perform principal components analysis on VST data
pcaData <- plotPCA(vsd, intgroup=c("Treatment", "CSC"), returnData=TRUE)
# Calculate the variance between groups
percentVar <- round(100 * attr(pcaData, "percentVar"))
# Save PCA plot as a PNG file
png(file=paste0(savedir,"/PCA.png"))
# Plot PCA data as a scatter plot
ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=CSC)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

# Calculate Euclidean distance between samples using VST data
sampleDists <- dist(t(assay(vsd)))
# Convert distances into a matrix
sampleDistMatrix <- as.matrix(sampleDists)
# Add groups to the row names
rownames(sampleDistMatrix) <- paste(vsd$Group, vsd$type, sep="-")
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


# Enrichment Gene Ontology Analysis ---------------------------------------

# List of comparisons for DGE analysis
contrast_list <- list(c("Group","NonCSC_Chemo","NonCSC_DMSO"),
                      c("Group","NonCSC_DOT1Li","NonCSC_DMSO"),
                      c("Group","NonCSC_Combo","NonCSC_DMSO"),
                      c("Group","CSC_Chemo","CSC_DMSO"),
                      c("Group","CSC_DOT1Li","CSC_DMSO"),
                      c("Group","CSC_Combo","CSC_DMSO"))
# Initialize empty dataframe to contain results for all comparisons
go_df <- data.frame(Ensembl=character(),
                    FC=double(),
                    reg=character(),
                    cond=character(),
                    stringsAsFactors=FALSE)
# Loop through list of comparisons
for (cont in contrast_list) {
  # Extract the result for the specific comparison
  res_cont <- results(dds, contrast=cont)
  # Sort the results by p-value
  resOrdered_cont <- res_cont[order(res_cont$pvalue),]
  # Filter out genes using prespecified p-value threshold
  resfilt_cont <- subset(resOrdered_cont, pvalue < pval_lim)
  # Create dataframe using fold change and gene IDs
  mydf <- data.frame(Ensembl=rownames(resfilt_cont),
                     FC=resfilt_cont$log2FoldChange)
  # Use comparison name to set cell type and treatment
  if (strsplit(cont[2], split = "_")[[1]][1]=="CSC") {
    mydf$type <- "ALDH+"
  }
  else {
    mydf$type <- "ALDH-"
  }
  if (strsplit(cont[2], split = "_")[[1]][2]=="Chemo") {
    mydf$cond <- "Cisplatin"
  }
  else if (strsplit(cont[2], split = "_")[[1]][2]=="DOT1Li") {
    mydf$cond <- "EPZ-5676"
  }
  else {
    mydf$cond <- "Combo"
  } 
  # Add values into the full dataframe for all comparisons
  go_df <- rbind(go_df, mydf)
}
# Sort the fold change for all comparisons dataframe
go_sort <- go_df[order(go_df$FC, decreasing = TRUE),]
# Perform BP GO enrichment analysis using fold change dataframe
# This function takes several minutes to run
formula_res <- compareCluster(Ensembl~type+cond,
                              data=go_sort,
                              fun="enrichGO",
                              ont ="BP", 
                              keyType = "ENSEMBL", 
                              OrgDb = org.Hs.eg.db)
# View GO results
head(formula_res)
# Save a dot plot showing GO terms as a TIFF file
tiff(filename = paste0(savedir,"/BP_GO_Dot_Plot.tiff"), width = 9, height = 9,
     units="in", res = 1000, compression = "lzw")
par(mar = c(3, 3, 0.5, 0.5))
# Create split dot plot for cell types with treatments on the x-axis
dotplot(formula_res, x="cond", showCategory=3, split=".type") + 
  facet_grid(~type) + 
  scale_x_discrete(limits = c("Cisplatin", "EPZ-5676", "Combo")) + 
  theme(axis.title.x = element_blank())
dev.off()

# Get GO results for chemotherapy treatment
chemo_res <- filter(formula_res,cond=="Cisplatin")
# Filter results for non-CSC ALDH- cells treated with chemotherapy
chemo_noncsc_res <- filter(chemo_res,type=="ALDH-")
# Save a dot plot showing GO terms as a TIFF file
tiff(filename = paste0(savedir,"/NonCSC_Chemo_BP_GO_Dot_Plot.tiff"), width = 12,
     height = 6, units="in", res = 1000, compression = "lzw")
par(mar = c(3, 3, 0.5, 0.5))
# Create split dot plot for cell types with treatments on the x-axis
dotplot(formula_res, x="cond", showCategory=10, split=".type") + 
  facet_grid(~type) + 
  scale_y_discrete(limits = chemo_noncsc_res[1:10,c('Description')]) + 
  scale_x_discrete(limits = c("Cisplatin", "EPZ-5676", "Combo")) + 
  theme(axis.title.x = element_blank())
dev.off()

# Filter results for CSC ALDH+ cells treated with chemotherapy
chemo_csc_res <- filter(chemo_res,type=="ALDH+")
# Save a dot plot showing GO terms as a TIFF file
tiff(filename = paste0(savedir,"/CSC_Chemo_BP_GO_Dot_Plot.tiff"), width = 12,
     height = 6, units="in", res = 1000, compression = "lzw")
par(mar = c(3, 3, 0.5, 0.5))
# Create split dot plot for cell types with treatments on the x-axis
dotplot(formula_res, x="cond", showCategory=10, split=".type") + 
  facet_grid(~type) + 
  scale_y_discrete(limits = chemo_csc_res[1:10,c('Description')]) + 
  scale_x_discrete(limits = c("Cisplatin", "EPZ-5676", "Combo")) + 
  theme(axis.title.x = element_blank())
dev.off()