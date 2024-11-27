# scRNAseq Analysis Script
# Authors: Rachel Ye and Jane Frederick
# Date: 2024-05-20
# Description: This script performs comprehensive scRNAseq analysis including data loading, normalization, PCA, t-SNE, and visualization.
# The script workflow includes:
# 1. Reading data from filtered feature barcode matrices.
# 2. Creating Seurat objects for each dataset.
# 3. Calculating mitochondrial gene percentages.
# 4. Filtering cells based on quality control metrics.
# 5. Normalizing the data.
# 6. Merging datasets and finding variable features.
# 7. Scaling the data and performing PCA.
# 8. Visualizing PCA loadings and generating an elbow plot.
# 9. Running t-SNE and extracting t-SNE embeddings.
# 10. Saving metadata to CSV files and plotting 3D t-SNE.
# 11. Performing weighted t-SNE analysis and exporting the weighted t-SNE embeddings as CSV files.

# Install necessary packages
install.packages(c('dplyr','ggplot2','ggpubr','EnvStats','glmGamPoi'))
install.packages('Seurat')
devtools::install_github("hhoeflin/hdf5r")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("glmGamPoi")
install.packages('scatterplot3d')
install.packages('rgl')
install.packages('rmarkdown')
install.packages('magick')
install.packages("data.table")

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(EnvStats)
library(glmGamPoi)
library(hdf5r)
library(magick)
library(rmarkdown)
library(data.table)

# Set working directory
setwd("/home/jfi977/Ovarian CSCs/scRNAseq/")

# Read data from filtered feature barcode matrices
pbs_dn_data_Alex <- Read10X(data.dir = "Matrix/PBS-DN/filtered_feature_bc_matrix/")
pbs_dp_data_Alex <- Read10X(data.dir = "Matrix/PBS-DP/filtered_feature_bc_matrix/")
carb_dn_data_Alex <- Read10X(data.dir = "Matrix/Carb-DN/filtered_feature_bc_matrix/")
carb_dp_data_Alex <- Read10X(data.dir = "Matrix/Carb-DP/filtered_feature_bc_matrix/")

# Create Seurat objects for each dataset
pbs_dn <- CreateSeuratObject(counts = pbs_dn_data_Alex, project = "pbs_dn", assay = "RNA", min.cells = 3, min.features = 200)
pbs_dp <- CreateSeuratObject(counts = pbs_dp_data_Alex, project = "pbs_dp", assay = "RNA", min.cells = 3, min.features = 200)
carb_dn <- CreateSeuratObject(counts = carb_dn_data_Alex, project = "carb_dn", assay = "RNA", min.cells = 3, min.features = 200)
carb_dp <- CreateSeuratObject(counts = carb_dp_data_Alex, project = "carb_dp", assay = "RNA", min.cells = 3, min.features = 200)

# Calculate the percentage of mitochondrial genes
pbs_dn[["percent.mt"]] <- PercentageFeatureSet(pbs_dn, pattern = "^MT-")
pbs_dp[["percent.mt"]] <- PercentageFeatureSet(pbs_dp, pattern = "^MT-")
carb_dn[["percent.mt"]] <- PercentageFeatureSet(carb_dn, pattern = "^MT-")
carb_dp[["percent.mt"]] <- PercentageFeatureSet(carb_dp, pattern = "^MT-")

# Filter cells based on quality control metrics
pbs_dn <- subset(pbs_dn, subset = nFeature_RNA > 200 & nCount_RNA < 1e05 & nCount_RNA > 20000 & percent.mt < 15)
pbs_dp <- subset(pbs_dp, subset = nFeature_RNA > 200 & nCount_RNA < 1e05 & nCount_RNA > 20000 & percent.mt < 15)
carb_dn <- subset(carb_dn, subset = nFeature_RNA > 200 & nCount_RNA < 1e05 & nCount_RNA > 20000 & percent.mt < 15)
carb_dp <- subset(carb_dp, subset = nFeature_RNA > 200 & nCount_RNA < 1e05 & nCount_RNA > 20000 & percent.mt < 15)

# Normalize the data
pbs_dn_norm <- NormalizeData(pbs_dn, normalization.method = "LogNormalize", scale.factor = 10000)
pbs_dp_norm <- NormalizeData(pbs_dp, normalization.method = "LogNormalize", scale.factor = 10000)
carb_dn_norm <- NormalizeData(carb_dn, normalization.method = "LogNormalize", scale.factor = 10000)
carb_dp_norm <- NormalizeData(carb_dp, normalization.method = "LogNormalize", scale.factor = 10000)

# Merge normalized datasets
seurat_ind_norm_merge <- merge(pbs_dn_norm, y = c(pbs_dp_norm, carb_dn_norm, carb_dp_norm), project = "CSC")

# Find variable features
seurat_ind_norm_merge <- FindVariableFeatures(seurat_ind_norm_merge, selection.method = "vst", nfeatures = 2000)

# Scale the data
all.genes <- rownames(seurat_ind_norm_merge)
merged_seurat_scale <- ScaleData(seurat_ind_norm_merge, features = all.genes)

# Perform PCA
merged_seurat_pca <- RunPCA(merged_seurat_scale, features = VariableFeatures(seurat_ind_norm_merge))

# Visualize PCA loadings and elbow plot
VizDimLoadings(merged_seurat_pca, dims = 1:3, reduction = "pca")
ElbowPlot(merged_seurat_pca)

# Run t-SNE
merged_seurat_tsne <- RunTSNE(object = merged_seurat_pca, dims = 1:20, dim.embed = 3)

# Extract t-SNE embeddings
merged_seurat_tsne_1 <- merged_seurat_tsne[["tsne"]]@cell.embeddings[,1]
merged_seurat_tsne_2 <- merged_seurat_tsne[["tsne"]]@cell.embeddings[,2]
merged_seurat_tsne_3 <- merged_seurat_tsne[["tsne"]]@cell.embeddings[,3]

# Save metadata to CSV
data_to_write_out <- as.data.frame(as.matrix(merged_seurat_tsne@meta.data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "tsne.csv")

# Plot 3D t-SNE
library(rgl) # interactive 3D plotting
plot3d(x = merged_seurat_tsne_1, y = merged_seurat_tsne_2, z = merged_seurat_tsne_3, col = as.numeric(1:4)[merged_seurat_pca@active.ident], type="p", radius=0.3)
cell.num <- table(merged_seurat_pca@active.ident)
seurat_label <- paste(names(cell.num), paste0("(n=", cell.num, ")"))
legend3d("topright", legend = seurat_label, pch = 16, cex=1, inset=c(0.02))
rgl.snapshot('3dplot.png', fmt = 'png')
# Individual t-SNE plot

# Run t-SNE for each dataset and save the Seurat objects
pbs_dn_tsne <- RunTSNE(object = pbs_dn_pca, dims = 1:20, dim.embed = 3)
SaveSeuratRds(pbs_dn_tsne, file = "PBS_DN_Seurat_Object.rds")
pbs_dp_tsne <- RunTSNE(object = pbs_dp_pca, dims = 1:20, dim.embed = 3)
SaveSeuratRds(pbs_dp_tsne, file = "PBS_DP_Seurat_Object.rds")
carb_dn_tsne <- RunTSNE(object = carb_dn_pca, dims = 1:20, dim.embed = 3)
SaveSeuratRds(carb_dn_tsne, file = "Cisplatin_DN_Seurat_Object.rds")
carb_dp_tsne <- RunTSNE(object = carb_dp_pca, dims = 1:20, dim.embed = 3)
SaveSeuratRds(carb_dp_tsne, file = "Cisplatin_DP_Seurat_Object.rds")

# Save metadata to CSV files
data_to_write_out <- as.data.frame(as.matrix(pbs_dn_tsne@meta.data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "pbs_dn_tsne.csv")
data_to_write_out <- as.data.frame(as.matrix(pbs_dp_tsne@meta.data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "pbs_dp_tsne.csv")
data_to_write_out <- as.data.frame(as.matrix(carb_dn_tsne@meta.data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "carb_dn_tsne.csv")
data_to_write_out <- as.data.frame(as.matrix(carb_dp_tsne@meta.data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "carb_dp_tsne.csv")

# Extract t-SNE embeddings for each dataset
pbs_dn_tsne_1 <- pbs_dn_tsne[["tsne"]]@cell.embeddings[,1]
pbs_dn_tsne_2 <- pbs_dn_tsne[["tsne"]]@cell.embeddings[,2]
pbs_dn_tsne_3 <- pbs_dn_tsne[["tsne"]]@cell.embeddings[,3]

pbs_dp_tsne_1 <- pbs_dp_tsne[["tsne"]]@cell.embeddings[,1]
pbs_dp_tsne_2 <- pbs_dp_tsne[["tsne"]]@cell.embeddings[,2]
pbs_dp_tsne_3 <- pbs_dp_tsne[["tsne"]]@cell.embeddings[,3]

carb_dn_tsne_1 <- carb_dn_tsne[["tsne"]]@cell.embeddings[,1]
carb_dn_tsne_2 <- carb_dn_tsne[["tsne"]]@cell.embeddings[,2]
carb_dn_tsne_3 <- carb_dn_tsne[["tsne"]]@cell.embeddings[,3]

carb_dp_tsne_1 <- carb_dp_tsne[["tsne"]]@cell.embeddings[,1]
carb_dp_tsne_2 <- carb_dp_tsne[["tsne"]]@cell.embeddings[,2]
carb_dp_tsne_3 <- carb_dp_tsne[["tsne"]]@cell.embeddings[,3]

# Plot 3D t-SNE for each dataset
open3d()
plot3d(x = pbs_dn_tsne_1, y = pbs_dn_tsne_2, z = pbs_dn_tsne_3, col = as.numeric(1:4)[pbs_dn_pca@active.ident], type="p", radius=0.3)
cell.num_pbs_dn <- table(pbs_dn_pca@active.ident)
seurat_label_pbs_dn <- paste(names(cell.num_pbs_dn), paste0("(n=", cell.num_pbs_dn, ")"))
legend3d("topright", legend = seurat_label_pbs_dn, pch = 16, cex=1, inset=c(0.02))

open3d()
plot3d(x = pbs_dp_tsne_1, y = pbs_dp_tsne_2, z = pbs_dp_tsne_3, col = as.numeric(1:4)[pbs_dp_pca@active.ident], type="p", radius=0.3)
cell.num_pbs_dp <- table(pbs_dp_pca@active.ident)
seurat_label_pbs_dp <- paste(names(cell.num_pbs_dp), paste0("(n=", cell.num_pbs_dp, ")"))
legend3d("topright", legend = seurat_label_pbs_dp, pch = 16, cex=1, inset=c(0.02))

open3d()
plot3d(x = carb_dn_tsne_1, y = carb_dn_tsne_2, z = carb_dn_tsne_3, col = as.numeric(1:4)[carb_dn_pca@active.ident], type="p", radius=0.3)
cell.num_carb_dn <- table(carb_dn_pca@active.ident)
seurat_label_carb_dn <- paste(names(cell.num_carb_dn), paste0("(n=", cell.num_carb_dn, ")"))
legend3d("topright", legend = seurat_label_carb_dn, pch = 16, cex=1, inset=c(0.02))

open3d()
plot3d(x = carb_dp_tsne_1, y = carb_dp_tsne_2, z = carb_dp_tsne_3, col = as.numeric(1:4)[carb_dp_pca@active.ident], type="p", radius=0.3)
cell.num_carb_dp <- table(carb_dp_pca@active.ident)
seurat_label_carb_dp <- paste(names(cell.num_carb_dp), paste0("(n=", cell.num_carb_dp, ")"))
legend3d("topright", legend = seurat_label_carb_dp, pch = 16, cex=1, inset=c(0.02))

rgl.snapshot('3dplot.png', fmt = 'png')

#################################################################################################

# Weighted t-SNE

# Step 1: Calculate weights for each dataset
group_size_pbs_dn <- table(pbs_dn_pca$orig.ident)
weight_pbs_dn <- 1/group_size_pbs_dn[as.character(pbs_dn_pca$orig.ident)]

group_size_pbs_dp <- table(pbs_dp_pca$orig.ident)
weight_pbs_dp <- 1/group_size_pbs_dp[as.character(pbs_dp_pca$orig.ident)]

group_size_carb_dn <- table(carb_dn_pca$orig.ident)
weight_carb_dn <- 1/group_size_carb_dn[as.character(carb_dn_pca$orig.ident)]

group_size_carb_dp <- table(carb_dp_pca$orig.ident)
weight_carb_dp <- 1/group_size_carb_dp[as.character(carb_dp_pca$orig.ident)]

# Apply weights to t-SNE embeddings
weighted_tsne_pbs_dn_1 <- pbs_dn_tsne_1 * sqrt(weight_pbs_dn)
weighted_tsne_pbs_dn_2 <- pbs_dn_tsne_2 * sqrt(weight_pbs_dn)
weighted_tsne_pbs_dn_3 <- pbs_dn_tsne_3 * sqrt(weight_pbs_dn)

weighted_tsne_pbs_dp_1 <- pbs_dp_tsne_1 * sqrt(weight_pbs_dp)
weighted_tsne_pbs_dp_2 <- pbs_dp_tsne_2 * sqrt(weight_pbs_dp)
weighted_tsne_pbs_dp_3 <- pbs_dp_tsne_3 * sqrt(weight_pbs_dp)

weighted_tsne_carb_dn_1 <- carb_dn_tsne_1 * sqrt(weight_carb_dn)
weighted_tsne_carb_dn_2 <- carb_dn_tsne_2 * sqrt(weight_carb_dn)
weighted_tsne_carb_dn_3 <- carb_dn_tsne_3 * sqrt(weight_carb_dn)

weighted_tsne_carb_dp_1 <- carb_dp_tsne_1 * sqrt(weight_carb_dp)
weighted_tsne_carb_dp_2 <- carb_dp_tsne_2 * sqrt(weight_carb_dp)
weighted_tsne_carb_dp_3 <- carb_dp_tsne_3 * sqrt(weight_carb_dp)

# Plot weighted 3D t-SNE for each dataset
open3d()
plot3d(x = weighted_tsne_pbs_dn_1, y = weighted_tsne_pbs_dn_2, z = weighted_tsne_pbs_dn_3, col = as.numeric(1:4)[pbs_dn_pca@active.ident], type="p", radius=0.3)
cell.num_pbs_dn <- table(pbs_dp_pca@active.ident)
seurat_label_pbs_dn <- names(cell.num_pbs_dn)
legend3d("topright", legend = seurat_label_pbs_dn, pch = 16, cex=1, inset=c(0.02))

open3d()
plot3d(x = weighted_tsne_pbs_dp_1, y = weighted_tsne_pbs_dp_2, z = weighted_tsne_pbs_dp_3, col = as.numeric(1:4)[pbs_dp_pca@active.ident], type="p", radius=0.3)
cell.num_pbs_dp <- table(pbs_dp_pca@active.ident)
seurat_label_pbs_dp <- names(cell.num_pbs_dp)
legend3d("topright", legend = seurat_label_pbs_dp, pch = 16, cex=1, inset=c(0.02))

open3d()
plot3d(x = weighted_tsne_carb_dn_1, y = weighted_tsne_carb_dn_2, z = weighted_tsne_carb_dn_3, col = as.numeric(1:4)[carb_dn_pca@active.ident], type="p", radius=0.3)
cell.num_carb_dn <- table(carb_dn_pca@active.ident)
seurat_label_carb_dn <- names(cell.num_carb_dn)
legend3d("topright", legend = seurat_label_carb_dn, pch = 16, cex=1, inset=c(0.02))

open3d()
plot3d(x = weighted_tsne_carb_dp_1, y = weighted_tsne_carb_dp_2, z = weighted_tsne_carb_dp_3, col = as.numeric(1:4)[carb_dp_pca@active.ident], type="p", radius=0.3)
cell.num_carb_dp <- table(carb_dp_pca@active.ident)
seurat_label_carb_dp <- names(cell.num_carb_dp)
legend3d("topright", legend = seurat_label_carb_dp, pch = 16, cex=1, inset=c(0.02))

rgl.snapshot('3dplot.png', fmt = 'png')

# Export weighted t-SNE embeddings as CSV files
write.csv(weighted_tsne_pbs_dn_1, row.names = FALSE, file = "weighted_pbs_dn_tsne_1.csv")
write.csv(weighted_tsne_pbs_dn_2, row.names = FALSE, file = "weighted_pbs_dn_tsne_2.csv")
write.csv(weighted_tsne_pbs_dn_3, row.names = FALSE, file = "weighted_pbs_dn_tsne_3.csv")

write.csv(weighted_tsne_pbs_dp_1, row.names = FALSE, file = "weighted_pbs_dp_tsne_1.csv")
write.csv(weighted_tsne_pbs_dp_2, row.names = FALSE, file = "weighted_pbs_dp_tsne_2.csv")
write.csv(weighted_tsne_pbs_dp_3, row.names = FALSE, file = "weighted_pbs_dp_tsne_3.csv")

write.csv(weighted_tsne_carb_dn_1, row.names = FALSE, file = "weighted_carb_dn_tsne_1.csv")
write.csv(weighted_tsne_carb_dn_2, row.names = FALSE, file = "weighted_carb_dn_tsne_2.csv")
write.csv(weighted_tsne_carb_dn_3, row.names = FALSE, file = "weighted_carb_dn_tsne_3.csv")

write.csv(weighted_tsne_carb_dp_1, row.names = FALSE, file = "weighted_carb_dp_tsne_1.csv")
write.csv(weighted_tsne_carb_dp_2, row.names = FALSE, file = "weighted_carb_dp_tsne_2.csv")
write.csv(weighted_tsne_carb_dp_3, row.names = FALSE, file = "weighted_carb_dp_tsne_3.csv")
