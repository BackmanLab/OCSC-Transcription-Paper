# 
# This script performs single-cell RNA sequencing (scRNA-seq) analysis using the Seurat package.
# The analysis includes data loading, quality control, normalization, clustering, and differential expression analysis.
# 
# Author: Jane Frederick
# Date: 2024-5-20
# 
# Libraries Used:
# - Matrix: For handling sparse matrices.
# - Seurat: For scRNA-seq data analysis.
# - tidyverse: For data manipulation and visualization.
# - magrittr: For pipe operations.
# - data.table: For efficient data manipulation.
# - clusterProfiler: For gene ontology enrichment analysis.
# - org.Hs.eg.db: For human gene annotation.
# 
# Steps:
# 1. List sample directories and check their existence.
# 2. Create Seurat objects for each sample and merge them.
# 3. Calculate quality control (QC) metrics and plot violin plots for these metrics.
# 4. Subset the Seurat object based on QC metrics and plot violin plots again.
# 5. Set identities and plot UMAP and t-SNE for visualization.
# 6. Perform differential expression analysis and gene ontology enrichment analysis.
# 7. Load data for individual samples, create Seurat objects, and combine them.
# 8. Normalize the data and perform differential expression analysis again.
# 9. Pseudobulk the counts based on donor-condition-celltype.
#

# Load Data

# install.packages('Seurat')
library(Matrix)
library(Seurat)

library(tidyverse)
library(magrittr)

# Step 1: List sample directories
dir.ls <- list.dirs(path = '/home/jfi977/Ovarian CSCs/scRNAseq/Matrix',
      full.names = T,
      recursive = F)

dir.ls %<>% map( ~ paste0(.x, "/filtered_feature_bc_matrix"))
names(dir.ls) <- c("carb_dn", "carb_dp", "pbs_dn", "pbs_dp")

# Step 2: Check whether directories exist
dir.ls %>% map( ~ dir.exists(.x))

# Step 3: Create Seurat objects per sample
obj.ls <- dir.ls %>% map( ~ Read10X(.x)) %>% map( ~ CreateSeuratObject(.x, min.cells = 3))

# Merge Seurat objects
combined <- merge(x = obj.ls[[1]],
      y = obj.ls[2:4],
      add.cell.ids = names(dir.ls))

# Calculate mitochondrial percentage
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

# Calculate log10 genes per UMI
combined$log10GenesPerUMI <- log10(combined$nFeature_RNA) / log10(combined$nCount_RNA)

# Add sample identifiers to metadata
sample_names <- sub("(.*?_.*?)_.*", "\\1", colnames(combined))
names(sample_names) <- colnames(combined)
combined <- AddMetaData(
  object = combined,
  metadata = sample_names,
  col.name = 'sample.idents'
)
head(combined[[]])

# Convert metadata to data.table
library(data.table)
df <- as.data.table(combined@meta.data)
sel <- c("sample.idents", "nCount_RNA", "nFeature_RNA", "percent.mt", "log10GenesPerUMI")
df <- df[, sel, with = FALSE]
df

# Plot violin plots for QC metrics
fontsize <- 10
linesize <- 0.35

gp.ls <- df[, 2:5] %>% imap( ~ {
  
  # Define label function
  give.n <- function(x) {
  return(c(y = median(x) + max(x) / 10, label = round(median(x), 2)))
  }
  
  # Assign colors
  col.ls <-
  setNames(
  c('lightpink2', 'lightblue2', 'lightgreen', 'coral1'),
  c("nCount_RNA", "nFeature_RNA", "percent.mt", "log10GenesPerUMI")
  )
  
  ggplot(data = df, aes(x = sample.idents, y = .x)) +
  geom_violin(trim = FALSE, fill = col.ls[.y]) +
  ggtitle(label = .y) + ylab(label = .y) +
  theme_bw() +
  theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_blank()
  ) +
  theme(
  axis.text = element_text(size = fontsize),
  axis.line = element_line(colour = "black", size = linesize),
  axis.ticks = element_line(size = linesize),
  axis.title.x = element_blank(),
  axis.ticks.length = unit(.05, "cm"),
  plot.title = element_text(size = fontsize + 2, hjust = 0.5),
  legend.position = 'none'
  ) +
  stat_summary(fun = median, geom = "point", col = "black") +  # Add points to plot
  stat_summary(fun.data = give.n,
     geom = "text",
     col = "black")
})

library(gridExtra)
grid.arrange(gp.ls[[1]], gp.ls[[2]], gp.ls[[3]], gp.ls[[4]], ncol = 2)

# Subset Seurat object based on QC metrics
combined <- subset(combined, subset = nFeature_RNA > 200 & nCount_RNA < 1e05 & nCount_RNA > 20000 & percent.mt < 15)

# Convert metadata to data.table
df <- as.data.table(combined@meta.data)
sel <- c("sample.idents", "nCount_RNA", "nFeature_RNA", "percent.mt", "log10GenesPerUMI")
df <- df[, sel, with = FALSE]
df

# Plot violin plots for QC metrics after filtering
fontsize <- 10
linesize <- 0.35

gp.ls <- df[, 2:5] %>% imap( ~ {
  
  # Define label function
  give.n <- function(x) {
  return(c(y = median(x) + max(x) / 10, label = round(median(x), 2)))
  }
  
  # Assign colors
  col.ls <-
  setNames(
  c('lightpink2', 'lightblue2', 'lightgreen', 'coral1'),
  c("nCount_RNA", "nFeature_RNA", "percent.mt", "log10GenesPerUMI")
  )
  
  ggplot(data = df, aes(x = sample.idents, y = .x)) +
  geom_violin(trim = FALSE, fill = col.ls[.y]) +
  ggtitle(label = .y) + ylab(label = .y) +
  theme_bw() +
  theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_blank()
  ) +
  theme(
  axis.text = element_text(size = fontsize),
  axis.line = element_line(colour = "black", size = linesize),
  axis.ticks = element_line(size = linesize),
  axis.title.x = element_blank(),
  axis.ticks.length = unit(.05, "cm"),
  plot.title = element_text(size = fontsize + 2, hjust = 0.5),
  legend.position = 'none'
  ) +
  stat_summary(fun = median, geom = "point", col = "black") +  # Add points to plot
  stat_summary(fun.data = give.n,
     geom = "text",
     col = "black")
})

grid.arrange(gp.ls[[1]], gp.ls[[2]], gp.ls[[3]], gp.ls[[4]], ncol = 2)

# Set identities and plot UMAP and t-SNE
Idents(object = combined) <- 'Individual'
DimPlot(object = combined, reduction = "umap", label = T) 
DimPlot(object = combined, reduction = "tsne", label = T) 

# Convert metadata to data.table
df <- data.table(combined@meta.data)
sel.meta <- c("Individual", str_c('SCT_snn_res.', c(1, 1.5, 1.8)))
df <- df[, sel.meta, with = FALSE]

# Plot bar plots for cluster proportions
df[, 2:4] %>% imap(~ {
  freq1 <- df[, .N, keyby = .(.x, Individual)]
  freq1[, total := sum(N), by = .(.x)]
  freq1[, ratio := N / total]
  
  linesize = .35
  fontsize = 8
  
  ggplot(freq1, aes(fill = Individual, y = ratio, x = .x)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_brewer(palette = "Dark2") +
  xlab('Clsuter') +
  ggtitle(.y) +
  scale_y_continuous(breaks = seq(0, 1, 0.1),
       expand = c(0, 0),
       name = 'Percentage') +
  theme_bw() +
  theme(
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(size = linesize),
  axis.ticks = element_blank(),
  axis.text.x = element_text(size = 5)
  ) +
  coord_flip()
})

# Normalize data
combined <- NormalizeData(combined)
Idents(combined) <- sample_names

# Differential expression analysis
deg.ls <- split(rownames(combined), f = combined$sample.idents)

library(clusterProfiler)
library(org.Hs.eg.db)

# Perform GO enrichment analysis
compGO <- compareCluster(geneCluster = deg.ls,
       fun = "enrichGO",
       keyType = "SYMBOL",
       OrgDb = org.Hs.eg.db, 
       ont = 'BP')

# Dot plot for GO enrichment results
dotplot(compGO, showCategory = 5) 

# Join RNA layers
combined[["RNA"]] <- JoinLayers(combined[["RNA"]])

# Find DE features between PBS-DP and PBS-DN
csc.de.markers <- FindMarkers(combined, ident.1 = "pbs_dp", ident.2 = "pbs_dn")
head(csc.de.markers)

# Find DE features between Carb-DN and PBS-DN
noncsccarb.de.markers <- FindMarkers(all_seurat, ident.1 = "carb_dn", ident.2 = "pbs_dn")
head(noncsccarb.de.markers)

# Find DE features between Carb-DP and PBS-DP
csccarb.de.markers <- FindMarkers(all_seurat, ident.1 = "carb_dp", ident.2 = "pbs_dp")
head(csccarb.de.markers)

# Plot DE genes
de.genes <- rownames(csccarb.de.markers)[grep('MT',rownames(csccarb.de.markers),invert=TRUE)]
RidgePlot(combined, features = de.genes[1])
DotPlot(combined, features = de.genes[1:10]) + RotatedAxis()

# Enrichment analysis using enrichR
DEenrichRPlot(combined, ident.1 = "pbs_dp", ident.2 = "pbs_dn", enrich.database = "GO_Biological_Process_2023", max.genes = 100)

# SCTransform and clustering
combined %<>% 
  SCTransform(return.only.var.genes = FALSE) %>%
  RunPCA(features = VariableFeatures(object = .)) %>%
  FindNeighbors(dims = 1:40) %>%
  FindClusters(resolution = c(0.5, 0.6, 0.8, 1, 1.2, 1.5, 1.8,2,2.5,3)) %>%
  RunUMAP(dims = 1:40) %>%
  RunTSNE(dims = 1:40)

# Load data for individual samples
pbs_dn_data <- Read10X(data.dir = "/home/jfi977/Ovarian CSCs/scRNAseq/Matrix/PBS-DN/filtered_feature_bc_matrix/")
pbs_dp_data <- Read10X(data.dir = "/home/jfi977/Ovarian CSCs/scRNAseq/Matrix/PBS-DP/filtered_feature_bc_matrix/")
carb_dn_data <- Read10X(data.dir = "/home/jfi977/Ovarian CSCs/scRNAseq/Matrix/Carb-DN/filtered_feature_bc_matrix/")
carb_dp_data <- Read10X(data.dir = "/home/jfi977/Ovarian CSCs/scRNAseq/Matrix/Carb-DP/filtered_feature_bc_matrix/")

# Create Seurat objects for individual samples
pbs_dn <- CreateSeuratObject(counts = pbs_dn_data, project = "pbs_dn", assay = "RNA", min.cells = 3, min.features = 200)
pbs_dp <- CreateSeuratObject(counts = pbs_dp_data, project = "pbs_dp", assay = "RNA", min.cells = 3, min.features = 200)
carb_dn <- CreateSeuratObject(counts = carb_dn_data, project = "carb_dn", assay = "RNA", min.cells = 3, min.features = 200)
carb_dp <- CreateSeuratObject(counts = carb_dp_data, project = "carb_dp", assay = "RNA", min.cells = 3, min.features = 200)

# Combine Seurat objects
all_seurat <- merge(pbs_dn, y = c(pbs_dp, carb_dn, carb_dp), add.cell.ids = c("pbs_dn", "pbs_dp", "carb_dn","carb_dp"), project = "CSC", merge.data = TRUE)

# Add mitochondrial read percentage metadata
all_seurat[["percent.mt"]] <- PercentageFeatureSet(all_seurat, pattern = "^MT-")

# Subset Seurat objects based on QC metrics
all_seurat <- subset(all_seurat, subset = nFeature_RNA > 200 & nCount_RNA < 1e05 & nCount_RNA > 20000 & percent.mt < 15)

# Normalize the data
all_seurat <- NormalizeData(all_seurat)

# Join RNA layers
all_seurat[["RNA"]] <- JoinLayers(all_seurat[["RNA"]])

# Find DE features between PBS-DP and PBS-DN
csc.de.markers <- FindMarkers(all_seurat, ident.1 = "pbs_dp", ident.2 = "pbs_dn")
head(csc.de.markers)

# Enrichment analysis using enrichR
DEenrichRPlot(all_seurat, ident.1 = "pbs_dp", ident.2 = "pbs_dn", enrich.database = "GO_Biological_Process_2023", max.genes = 100)

# Find DE features between Carb-DN and PBS-DN
noncsccarb.de.markers <- FindMarkers(all_seurat, ident.1 = "carb_dn", ident.2 = "pbs_dn")
head(noncsccarb.de.markers)

# Find DE features between Carb-DP and PBS-DP
csccarb.de.markers <- FindMarkers(all_seurat, ident.1 = "carb_dp", ident.2 = "pbs_dp")
head(csccarb.de.markers)

# Pseudobulk the counts based on donor-condition-celltype
pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, group.by = c("stim", "donor_id", "seurat_annotations"))

# Each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_ifnb))
