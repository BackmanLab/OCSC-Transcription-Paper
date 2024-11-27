# -------------------------------------------------------------------------
# Script Name: Bulk_DGE_Analysis_Treatment.R
# Description: Gene ontology and differential gene expression (DGE) analysis of
# OVCAR5 ALDH+ vs. ALDH- cells treated with DMSO, cisplatin, EPZ-5676, or both.
# 
# Author: Jane Frederick
# Created on: 8-17-2023
# Last updated: 11-20-2024
# 
# Dependencies: BiocManager, clusterProfiler, enrichplot, msigdbr, pheatmap,
# apeglm, DESeq2, tximport, org.Hs.eg.db, ggplot2, svglite, RColorBrewer
#
# -------------------------------------------------------------------------

# Load Libraries ----------------------------------------------------------

# List of required libraries
required_packages <- c("BiocManager", "clusterProfiler", "enrichplot", 
             "msigdbr", "pheatmap", "apeglm", 
             "DESeq2", "tximport", "org.Hs.eg.db", "ggplot2",
             "svglite", "RColorBrewer")

# Check and install missing packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
  if (pkg == "BiocManager") {
    install.packages(pkg)
  } else {
    BiocManager::install(pkg)
  }
  }
}

# Load all libraries
lapply(required_packages, library, character.only = TRUE)

# Set Directories ---------------------------------------------------------

# Define directories
parentdir <- "C:/Users/janef/OneDrive - Northwestern University/Documents - Backman Lab - Shared Folders/Lab Paper Drafts/OCSC Chromatin Transcription Paper 2024"
seqdir <- file.path(parentdir, "Data/Bulk RNA-seq/Treated") # Folder with all files
RSEMdir <- file.path(seqdir, "TPM Files") # Folder with .genes.results files
savedir <- file.path(seqdir, "DGE Analysis") # Folder to save results

# Set working directory to parent directory
setwd(parentdir)

# Ensure save directory exists
if (!dir.exists(savedir)) dir.create(savedir, recursive = TRUE)

# Set Global Variables ----------------------------------------------------

# Analysis parameters
deseq_ref <- "NonCSC_DMSO"  # Reference group
pval_lim <- 0.1             # p-value significance threshold
padj_lim <- 0.05            # Adjusted p-value significance threshold

# Read and Process Sample Metadata ----------------------------------------

# Read in samples file containing information about data files
samples_file <- file.path(seqdir, "samples.txt")
samples <- read.table(samples_file, header = TRUE, sep = "\t")
# Set the Group, CSC, and Treatment columns as factor columns
samples$Group <- factor(samples$Group)
samples$CSC <- factor(samples$CSC)
samples$Treatment <- factor(samples$Treatment)
# Change the row names of the samples list to be group and replicate
rownames(samples) <- paste(samples$Group, samples$Replicate, sep = "_")

# Create file paths for RSEM data
files <- file.path(RSEMdir, paste0(samples$ID, "Aligned.genes.results"))
names(files) <- rownames(samples)

# Read and Annotate RNA TPM Data ------------------------------------------

# Initialize an empty data frame for TPM data
read_TPM <- data.frame(gene_id = NA)

# Loop through each sample and merge TPM data into a single data frame
for (sample in names(files)) {
  sample_data <- read.table(
  files[sample],
  col.names = c("gene_id", "transcript_id(s)", "length", "effective_length",
          paste0(sample, "_count"), paste0(sample, "_TPM"), "FPKM")
  )
  # Extract the TPM column and add to the TPM data frame
  read_TPM <- merge(
  read_TPM,
  sample_data[, c("gene_id", paste0(sample, "_TPM"))],
  by = "gene_id",
  all = TRUE
  )
}

# Remove TPM suffix from column names
colnames(read_TPM) <- sub("_TPM", "", colnames(read_TPM))

# Annotate with HGNC gene symbols
genes <- read_TPM$gene_id
annots <- AnnotationDbi::select(org.Hs.eg.db, keys = genes, columns = "SYMBOL",
                keytype = "ENSEMBL")
TPMresult <- merge(read_TPM, annots, by.x = "gene_id", by.y = "ENSEMBL")

# Save TPM results to a CSV file
write.csv(TPMresult, file = file.path(savedir, "All_Samples_TPM.csv"))

# Differential Gene Expression Analysis -----------------------------------

# Use the tximport package to read in RSEM files
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

# Remove zero-length or unexpressed genes
valid_genes <- !(apply(txi.rsem$abundance, 1, max) == 0 & 
           apply(txi.rsem$length, 1, min) == 0)
# Subset elements with two dimensions
txi.rsem <- lapply(txi.rsem, function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
  x[valid_genes, , drop = FALSE] # Ensure it retains two dimensions
  } else {
  x
  }
})

# Create a DESeq dataset from the tximport object
dds <- DESeqDataSetFromTximport(txi.rsem, colData = samples, design = ~Group)
# Re-level the contrast to compare to the reference group
dds$Group <- relevel(dds$Group, ref = deseq_ref)
dds <- DESeq(dds)
res <- na.omit(results(dds))
# Display the results
head(res)

# Sort the results by p-value
resOrdered <- res[order(res$pvalue), ]
# Filter the results using the adjusted p-value
resFiltered <- subset(resOrdered, padj < padj_lim)

# PCA and Distance Plots --------------------------------------------------

# Perform variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)
# Perform principal components analysis on VST data
pcaData <- plotPCA(vsd, intgroup = c("Treatment", "CSC"), returnData = TRUE)
# Calculate the variance between groups
percentVar <- round(100 * attr(pcaData, "percentVar"))
# Save PCA plot
p <- ggplot(pcaData, aes(PC1, PC2, color = Treatment, shape = CSC)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
ggsave(file.path(savedir, "PCA.svg"), p, device = "svg")
ggsave(file.path(savedir, "PCA.tiff"), p, device = "tiff")

# Calculate Euclidean distance between samples using VST data
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(samples$Group, samples$Replicate, sep = "-")
colnames(sampleDistMatrix) <- NULL
# Define color palette for distance heatmap
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
# Save Euclidean distance heatmap
p <- pheatmap(sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors)
ggsave(file.path(savedir, "Distance_Heatmap.svg"), p, device = "svg")
ggsave(file.path(savedir, "Distance_Heatmap.tiff"), p, device = "tiff")

# Gene Ontology Enrichment Analysis ---------------------------------------

# List of comparisons for DGE analysis
contrast_list <- list(
  c("Group", "NonCSC_Chemo", "NonCSC_DMSO"),
  c("Group", "NonCSC_DOT1Li", "NonCSC_DMSO"),
  c("Group", "NonCSC_Combo", "NonCSC_DMSO"),
  c("Group", "CSC_Chemo", "CSC_DMSO"),
  c("Group", "CSC_DOT1Li", "CSC_DMSO"),
  c("Group", "CSC_Combo", "CSC_DMSO")
)

# Initialize an empty data frame for GO analysis results
go_df <- data.frame()

# Loop through each contrast and perform GO enrichment analysis
for (cont in contrast_list) {
  # Extract the result for the specific comparison
  res_cont <- results(dds, contrast = cont)
  resFiltered_cont <- subset(res_cont, pvalue < pval_lim)
  mydf <- data.frame(
  Ensembl = rownames(resFiltered_cont),
  FC = resFiltered_cont$log2FoldChange,
  type = ifelse(grepl("^CSC_", cont[2]), "ALDH+", "ALDH-"),
  cond = dplyr::recode(
    sub("NonCSC_|CSC_", "", cont[2]),
    "Chemo" = "Cisplatin",
    "DOT1Li" = "EPZ-5676",
    "Combo" = "Combo"
  )
  )
  go_df <- rbind(go_df, mydf)
}

# Perform BP GO enrichment analysis using fold change dataframe
formula_res <- compareCluster(
  Ensembl ~ type + cond,
  data = go_df,
  fun = "enrichGO",
  ont = "BP",
  keyType = "ENSEMBL",
  OrgDb = org.Hs.eg.db
)

# View GO results
head(formula_res)

# Save a dot plot showing GO terms
p <- dotplot(formula_res, x = "cond", showCategory = 3, split = ".type") + 
  facet_grid(~type) + 
  scale_x_discrete(limits = c("Cisplatin", "EPZ-5676", "Combo")) + 
  theme(axis.title.x = element_blank())
ggsave(
  filename = file.path(savedir, "BP_GO_Dot_Plot.svg"),
  plot = p,
  width = 9,
  height = 9,
  units = "in",
  device = "svg"
)
ggsave(
  filename = file.path(savedir, "BP_GO_Dot_Plot.tiff"),
  plot = p,
  width = 9,
  height = 9,
  units = "in",
  device = "tiff"
)

# Gene Set Enrichment Analysis --------------------------------------------

# Get gene sets from the molecular signatures database
m_df <- msigdbr(species = "Homo sapiens") %>% 
  dplyr::select(gs_name, ensembl_gene)

# Define the correct contrast (comparison) for ALDH+ differential expression
csc_dot1li_contrast <- c("Group", "CSC_DOT1Li", "CSC_DMSO")
# Extract the DESeq results using the correct contrast
res_csc_dot1li <- results(dds, contrast=csc_dot1li_contrast)
# Display the results
head(res_csc_dot1li)

# Extract fold change values from the results object
geneList_csc_dot1li <- res_csc_dot1li$log2FoldChange
# Add the ENSEMBL gene IDs to the fold change list
names(geneList_csc_dot1li) <- rownames(res_csc_dot1li)
# Sort the list by fold change
geneList_csc_dot1li <- sort(geneList_csc_dot1li, decreasing = TRUE)

# Runs GSEA on all gene sets
em <- GSEA(geneList_csc_dot1li, TERM2GENE = m_df)
# Convert the gseaResult object to a data frame
gsea_df <- as.data.frame(em)
# Export GSEA results for all gene sets
write.csv(gsea_df, file=paste0(savedir,"/All_GSEA_List.csv"))

# Extract the adjusted p-value
p_value <- gsea_df["BHATTACHARYA_EMBRYONIC_STEM_CELL","p.adjust"]
# Generate the GSEA plot for the gene set
p <- gseaplot2(em, geneSetID = "BHATTACHARYA_EMBRYONIC_STEM_CELL",
         title = paste0("BHATTACHARYA_EMBRYONIC_STEM_CELL",
                "\nAdjusted p-value: ", signif(p_value, 3)))
# Save the plot
ggsave(
  filename = file.path(savedir, "BHATTACHARYA_EMBRYONIC_STEM_CELL.svg"),
  plot = p,
  width = 7,
  height = 7,
  units = "in",
  device = "svg"
)
ggsave(
  filename = file.path(savedir, "BHATTACHARYA_EMBRYONIC_STEM_CELL.tiff"),
  plot = p,
  width = 7,
  height = 7,
  units = "in",
  device = "tiff"
)

# Extract the adjusted p-value
p_value <- gsea_df["BROWN_MYELOID_CELL_DEVELOPMENT_DN","p.adjust"]
# Generate the GSEA plot for the gene set
p <- gseaplot2(em, geneSetID = "BROWN_MYELOID_CELL_DEVELOPMENT_DN",
         title = paste0("BROWN_MYELOID_CELL_DEVELOPMENT_DN",
                "\nAdjusted p-value: ", signif(p_value, 3)))
# Save the plot
ggsave(
  filename = file.path(savedir, "BROWN_MYELOID_CELL_DEVELOPMENT_DN.svg"),
  plot = p,
  width = 7,
  height = 7,
  units = "in",
  device = "svg"
)
ggsave(
  filename = file.path(savedir, "BROWN_MYELOID_CELL_DEVELOPMENT_DN.tiff"),
  plot = p,
  width = 7,
  height = 7,
  units = "in",
  device = "tiff"
)