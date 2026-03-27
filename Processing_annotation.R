library(dplyr)
library(Seurat)
library(SeuratObject)
library(sctransform)
library(ggplot2)
library(ggpubr)
library(NMF)
library(remotes)
library(SeuratDisk)
library(Matrix)
library(openxlsx)
library(Nebulosa)
library(HGNChelper)
library(scCustomize)
library(AnnotationDbi)
library(CellChat)
library(scRNAtoolVis)
library(cowplot)
library(pheatmap)
library(dittoSeq)
library(genekitr)
library(enrichplot)
library(SCP)


# Function to process each object
process_sample <- function(sample_id) {
  path <- sprintf("%s.raw_feature_bc_matrix.h5", sample_id)
  obj <- CreateSeuratObject(Read10X_h5(path), min.features = 200, project = sample_id)
  obj <- RenameCells(obj, add.cell.id = sample_id)
  
  # Calculate percentages
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HB")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  
  # Gene annotations
  mito_genes <- grep("^MT-", rownames(obj), value = TRUE)
  ribo_genes <- grep("^RP[SL]", rownames(obj), value = TRUE)
  hemo_genes <- grep("^HB", rownames(obj), value = TRUE)
  
  # Apply QC filters
  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 35 & percent.hb < 0.5)
  
  # SCTransform
  obj <- SCTransform(obj, vars.to.regress = c("percent.mt", "percent.ribo", "percent.hb"), verbose = TRUE, variable.features.n = NULL)
  new_var_genes <- setdiff(VariableFeatures(obj), c(mito_genes, ribo_genes, hemo_genes, "XIST"))
  VariableFeatures(obj) <- new_var_genes
  
  return(obj)
}

# Sample IDs
sample_ids <- c("P1004", "P1005", "P1006", "P1007", "P1008", "P1010", "P1011", "P1016", "AT1010D","AT1011D", "AT1013D", "AT1016D")
samples_list <- lapply(sample_ids, process_sample)

# Assign original sample ID and status
statuses <- c("Non-permissive", "Non-permissive", "Permissive", "Non-permissive", "Permissive", "Non-permissive", "Permissive", "Permissive", "Non-permissive_0", "Permissive_0", "Non-permissive_0", "Permissive_0")
names(samples_list) <- sample_ids
for (i in seq_along(samples_list)) {
  samples_list[[i]]$status <- statuses[i]
}

VlnPlot(yearstatus, group.by = "orig.ident", features = "percent.hb")

################################Analysis#########################################
# Integration features selection
features <- SelectIntegrationFeatures(object.list = samples_list, nfeatures = 3000)

# Preparing for SCT integration
yearstatus.list <- PrepSCTIntegration(object.list = samples_list, anchor.features = features)

# Finding integration anchors
immune.anchors <- FindIntegrationAnchors(object.list = yearstatus.list, normalization.method = "SCT", anchor.features = features)

# Integrating data
yearstatus <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

# Principal Component Analysis
yearstatus <- RunPCA(yearstatus, features = VariableFeatures(object = yearstatus), npcs = 50, verbose = TRUE)

# Elbow Plot to determine the number of PCs to use
ElbowPlot(object = yearstatus, ndims = 50)

# UMAP for visualization
yearstatus <- RunUMAP(object = yearstatus, dims = 1:30, verbose = TRUE, n.neighbors = 20)

# Finding neighbors
yearstatus <- FindNeighbors(object = yearstatus, dims = 1:30, verbose = TRUE, k.param = 20)

# Clustering
yearstatus <- FindClusters(object = yearstatus, verbose = TRUE, resolution = 1.0)


######################################################Cell typist annotation##################################################### 

#Labeling with celltypes; Access assay for celltypist
raw_counts_yearstatus <- GetAssayData(object = yearstatus, assay = "RNA", layer = "counts")
dense_counts <- as.data.frame(as.matrix(raw_counts_yearstatus))
write.csv(dense_counts, file = "raw_counts_yearstatus.csv", row.names = TRUE)

#Added labels
predicted_labels <- read.csv("predicted_labels.csv", row.names = 1)
yearstatus <- AddMetaData(yearstatus, metadata = predicted_labels$predicted_label, col.name = "predicted_cell_type")
yearstatus <- AddMetaData(yearstatus, metadata = predicted_labels$majority_voting, col.name = "majority_voting")
Idents(object = yearstatus) <- "majority_voting"
