library(dplyr)
library(Seurat)
library(SeuratObject)
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
library(scRNAtoolVis)
library(cowplot)
library(pheatmap)
library(dittoSeq)
library(genekitr)
library(enrichplot)
library(SCP)

######################################Subsetting data by types#############################

# Split into permissive, non-permissive groups
yearstatus <- NormalizeData(yearstatus, assay = "RNA")
permissive <- subset(yearstatus, subCSDCDset = status == "Permissive")
nonpermissive <- subset(yearstatus, subset = status == "Non-permissive")
permissive0 <- subset(yearstatus, subset = status == "Permissive_0")
nonpermissive0 <- subset(yearstatus, subset = status == "Non-permissive_0")
PversusNP <- subset(yearstatus, status %in% c("Permissive", "Non-permissive"))
PversusNP0 <- subset(yearstatus, status %in% c("Permissive_0", "Non-permissive_0"))

############################################DE testing_12M#############################################
#Subseting CD4 T cells
CD4 <- subset(PversusNP, subset = majority_voting == "CD4 T cells")
FeatureDimPlot(CD4, features = "CD4", reduction = "UMAP")
CD4 <- subset(x = CD4, subset = CD4 > 0.5)
CD4 <- SCTransform(CD4, vars.to.regress = c("percent.mt", "percent.ribo", "percent.hb"), verbose = TRUE)
CD4 <- RunPCA(CD4, verbose = FALSE)
CD4 <- RunUMAP(CD4, dims = 1:50, verbose = FALSE)
CD4 <- FindNeighbors(CD4, dims = 1:50, verbose = FALSE)
CD4 <- FindClusters(CD4, verbose = FALSE)
CD4markers <- FindAllMarkers(object = CD4, min.pct = 0.25, logfc.threshold = 0.25)

CD4markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::filter(p_val < 0.0001) %>%
  slice_head(n = 3) %>%
  ungroup() -> top3

# final heatmap 
get_unique_genes <- function(data, gene_column) {
  genes <- data[[gene_column]]
  unique_genes <- unique(genes)
  return(unique_genes)
}

unique_genes_list <- get_unique_genes(top3, "gene")

DotPlot(CD4, features = unique_genes_list, col.min = 0, col.max = 5, cols = c("white", "#8b0000"), dot.scale = 3) + 
  rotate() + 
  RotatedAxis() + 
  theme(text = element_text(size = 9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9))

write.csv(top10, file = "top10.csv", row.names = TRUE)

#celltypist labeling
raw_counts_CD4 <- GetAssayData(object = CD4, assay = "RNA", layer = "counts")
dense_counts <- as.data.frame(as.matrix(raw_counts_CD4))
write.csv(dense_counts, file = "raw_counts_CD4.csv", row.names = TRUE)

predicted_labels <- read.csv("predicted_labels.csv", row.names = 1)
CD4 <- AddMetaData(CD4, metadata = predicted_labels$predicted_label, col.name = "predicted_cell_type_sub")
CD4 <- AddMetaData(CD4, metadata = predicted_labels$majority_voting, col.name = "majority_voting_sub")

#view results from celltypis
CellDimPlot(CD4, group.by = c("majority_voting_sub", "seurat_clusters"), reduction = "UMAP", theme_use = "theme_blank")
CellStatPlot(CD4, stat.by = "status", group.by = "majority_voting_sub", stat_type = "count", position = "dodge", label = TRUE, palcolor = c("#D2042D", "#0000FF"))
CellDimPlot(CD4, group.by = "majority_voting_sub", split.by = "status", reduction = "UMAP", cells.highlight = TRUE) 
CD4@meta.data$majority_voting_sub <- factor(CD4@meta.data$majority_voting_sub, levels = c("Tcm/Naive helper T cells", "Tem/Effector helper T cells", "Tem/Trm cytotoxic T cells", "Type 17 helper T cells", "Regulatory T cells", "MAIT cells"))

FeatureDimPlot(CD4, features = c("TCF7", "CXCR6"), split.by = "status", keep_scale = "feature")   
FeatureDimPlot(yearstatus, features = "IFNG", split.by = "status")
FeatureDimPlot(CD4, features = c("TCF7", "CXCR6"), split.by = "status", reduction = "UMAP", cells.highlight = TRUE) 

#DE 
CD4$celltype.stim <- paste(CD4$majority_voting, CD4$status, sep = "_")
DefaultAssay(object = CD4) <- "SCT"
Idents(CD4) <- "celltype.stim"
unique_titles <- unique(CD4$celltype.stim)
print(unique_titles)

CD4_T_cells_DE <- FindMarkers(CD4, ident.1 = "CD4 T cells_Permissive", ident.2 = "CD4 T cells_Non-permissive", verbose = TRUE, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
CD4_T_cells_DE$Significant <- with(CD4_T_cells_DE, ifelse(p_val_adj < 0.05, "Yes", "No"))
CD4_T_cells_filtered <- subset(CD4_T_cells_DE, Significant == "Yes")
Idents(CD4) <- "majority_voting"
genes1 <- c("IKZF2", "ID3", "HLA-B", "IRF1", "PTPN6", "ISG15", 
            "GBP4", "ITGB7", "STAT1", "CASP1", "SELL", "KLRG1", "CXCR3", "KLF2", 
            "LTB", "CD27", "HMGB1", "CD69", "PTGER4", "CXCR4", "ITGAE", 
            "TSC22D3", "IKZF3", "VSIR", "NME2", "CXCR6", "TCF7")
CD4@meta.data$status <- factor(CD4@meta.data$status, levels = c('Permissive', 'Non-permissive'))

GroupHeatmap(CD4, features = genes1,  group.by = "status", exp_method = "log2fc", 
             slot = "counts", nlabel = 35, 
             heatmap_palette = "material-blue", 
             group_palcolor = list(c("#0000FF", "#D2042D")), cluster_rows = TRUE)

CD4 <- RunDEtest(srt = CD4, group_by = "status", fc.threshold = 1)
VolcanoPlot(srt = CD4, group_by = "status", nlabel = 20)
CD4 <- RunEnrichment(
  srt = CD4, group_by = "status", db = c("GO_BP", "KEGG","MSigDB"), species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05")
EnrichmentPlot(srt = CD4, group_by = "status", plot_type = "comparison", db = c("GO_BP", "MSigDB"), 
               id_use = c("M5936", "M5890", "M5932", "M5935", "M5902", "M5913", "M5950", "M5911", "M5890", "M5902", "GO:0045577",
                          "GO:1902150", "GO:0045619", "GO:0009206", "GO:0009145", "GO:0002181", "GO:0042274", "GO:0071346", 
                          "GO:0034341", "GO:0042254", "GO:0032615"))

#Subsetting CD8 T cells
CD8 <- subset(PversusNP, subset = majority_voting == "CD8 T cells")
FeatureDimPlot(CD8, features = "CD8A", reduction = "UMAP")
CD8 <- subset(x = CD8, subset = CD8A > 0.5)
CD8 <- SCTransform(CD8, vars.to.regress = c("percent.mt", "percent.ribo", "percent.hb"), verbose = TRUE)
CD8 <- RunPCA(CD8, verbose = FALSE)
CD8 <- RunUMAP(CD8, dims = 1:50, verbose = FALSE)
CD8 <- FindNeighbors(CD8, dims = 1:50, verbose = FALSE)
CD8 <- FindClusters(CD8, verbose = FALSE)
CD8markers <- FindAllMarkers(object = CD8, min.pct = 0.25, logfc.threshold = 0.25)
Idents(object = CD8) <- "majority_voting_sub"
CellDimPlot(CD8, group.by = "majority_voting_sub", split.by = "status", reduction = "UMAP", cells.highlight = TRUE)
CellDimPlot(CD8, group.by = c("majority_voting_sub","seurat_clusters"), reduction = "UMAP")
CellStatPlot(CD8, stat.by = "status", group.by = "majority_voting_sub", stat_type = "count", position = "dodge", label = TRUE, palcolor = c("#D2042D", "#0000FF"))
CellDimPlot(CD8, group.by = "majority_voting_sub", split.by = "status", reduction = "UMAP", cells.highlight = TRUE)
FeatureDimPlot(CD8, features = "ITGAE")

CD8markers <- FindAllMarkers(object = CD8, min.pct = 0.25, logfc.threshold = 0.25)

CD8markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::filter(p_val < 0.00001) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

# final heatmap 
get_unique_genes <- function(data, gene_column) {
  genes <- data[[gene_column]]
  unique_genes <- unique(genes)
  return(unique_genes)
}

unique_genes_list <- get_unique_genes(top10, "gene")

DotPlot(CD8, features = unique_genes_list, col.min = 0, col.max = 5, cols = c("white", "#8b0000"), dot.scale = 3) + 
  rotate() + 
  RotatedAxis() + 
  theme(text = element_text(size = 9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9))

write_csv(top5, file = "cd8markers_top10.csv")

raw_counts_CD8 <- GetAssayData(object = CD8, assay = "RNA", layer = "counts")
dense_counts <- as.data.frame(as.matrix(raw_counts_CD8))
write.csv(dense_counts, file = "raw_counts_CD8.csv", row.names = TRUE)

predicted_labels <- read.csv("predicted_labels.csv", row.names = 1)
CD8 <- AddMetaData(CD8, metadata = predicted_labels$predicted_label, col.name = "predicted_cell_type_sub")
CD8 <- AddMetaData(CD8, metadata = predicted_labels$majority_voting, col.name = "majority_voting_sub")

CellDimPlot(CD8, group.by = "majority_voting_sub", split.by = "status", reduction = "UMAP", cells.highlight = TRUE)
CellDimPlot(CD8, group.by = c("majority_voting_sub","seurat_clusters"), reduction = "UMAP")
CellStatPlot(CD8, stat.by = "status", group.by = "seurat_clusters", stat_type = "count", position = "dodge", label = TRUE, palcolor = c("#0000FF", "#D2042D"))
CellDimPlot(CD8, group.by = "majority_voting_sub", split.by = "status", reduction = "UMAP", cells.highlight = TRUE)

# Unique titles
CD8$majority_voting_sub <- paste0(CD8$majority_voting_sub, "-", CD8$seurat_clusters)
unique_titles <- unique(CD8$majority_voting_sub)
print(unique_titles)

# Individual scattered cell relabeling
target_cell_types <- c("Tem/Temra cytotoxic T cells-7", "Tem/Trm cytotoxic T cells-7")
new_label <- "PDCD1+ CD8 T cell"
CD8$majority_voting_sub [CD8$majority_voting_sub %in% target_cell_types] <- new_label
CD8$majority_voting_sub <- sub("-.*", "", CD8$majority_voting_sub)

#####DE expression
CD8$celltype.stim <- paste(CD8$majority_voting, CD8$status, sep = "_")
DefaultAssay(object = CD8) <- "SCT"
Idents(CD8) <- "celltype.stim"
unique_titles <- unique(CD8$celltype.stim)
print(unique_titles)
CD8_cells_DE <- FindMarkers(CD8, ident.1 = "CD8 T cells_Permissive", ident.2 = "CD8 T cells_Non-permissive", verbose = TRUE, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
CD8_cells_DE$Significant <- with(CD8_cells_DE, ifelse(p_val_adj < 0.05, "Yes", "No"))
CD8_cells_filtered <- subset(CD8_cells_DE, Significant == "Yes")
Idents(CD8) <- "majority_voting"
genes1 <- c("HLA-B", "CD74", "HLA-DRB1", "HLA-DQA1", 
            "STAT1"," IRF1", "ITGB7", "IFI35", "CCL3L3", 
            "CXCL13", "PECAM1", "ADA2", "ID3", "CD70", "PDCD1", "IL32", "PDCD4", "ISG20", "CD2", 
            "CD69", "CASP8", "SOD2", "ITGB2", "CST7", "TNFRSF9", "CYBA", "TUBB", "ANXA5", "SELL", "FOXP1", "KIR2DL2", "TSC22D3", "STAT1", "SOCS1")
write.csv(CD8_cells_filtered, file = "CD8_cells_filtered.csv")
CD8@meta.data$status <- factor(CD8@meta.data$status, levels = c('Permissive', 'Non-permissive'))

GroupHeatmap(CD8, features = genes1,  group.by = "status", exp_method = "log2fc", 
             slot = "counts", nlabel = 35, features_fontsize = 4,
             heatmap_palette = "material-blue", 
             group_palcolor = list(c("#0000FF", "#D2042D")), cluster_rows = TRUE)

CD8 <- RunDEtest(srt = CD8, group_by = "status", fc.threshold = 1)
VolcanoPlot(srt = CD8, group_by = "status", nlabel = 20)
CD8 <- RunEnrichment(
  srt = CD8, group_by = "status", db = c("GO_BP","MSigDB"), species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05")

EnrichmentPlot(srt = CD8, group_by = "status", plot_type = "comparison", db = c("GO_BP", "MSigDB"), 
               id_use = c("GO:0002181", "GO:0002396", "GO:0002501", "GO:0048002", "GO:0019882", "GO:0034341",
                          "GO:0019646", "GO:0022904", "GO:0041773", "GO:0042775", "M5890","M5950", "M5902","M5943","M5938",
                          "M5913", "M5911", "M5936", "M5915"))

#subsetting CD16 positive NK cells
CD16positiveNK <- subset(PversusNP, subset = majority_voting == "CD16+ NK cells")
FeatureDimPlot(CD16positiveNK, features = "FCGR3A", reduction = "UMAP")
CD16positiveNK <- subset(x = CD16positiveNK, subset = FCGR3A > 0.5)
CD16positiveNK <- SCTransform(CD16positiveNK, vars.to.regress = c("percent.mt", "percent.ribo", "percent.hb"), verbose = TRUE)
CD16positiveNK <- RunPCA(CD16positiveNK, verbose = FALSE)
CD16positiveNK <- RunUMAP(CD16positiveNK, dims = 1:50, verbose = FALSE)
CD16positiveNK <- FindNeighbors(CD16positiveNK, dims = 1:50, verbose = FALSE)
CD16positiveNK <- FindClusters(CD16positiveNK, verbose = FALSE)
CD16positiveNKmarkers <- FindAllMarkers(object = CD16positiveNK, min.pct = 0.25, logfc.threshold = 0.25)

##DE expression
CD16positiveNK$celltype.stim <- paste(CD16positiveNK$majority_voting, CD16positiveNK$status, sep = "_")
DefaultAssay(object = CD16positiveNK) <- "SCT"
Idents(CD16positiveNK) <- "celltype.stim"
unique_titles <- unique(CD16positiveNK$celltype.stim)
print(unique_titles)
CD16positiveNK_cells_DE <- FindMarkers(CD16positiveNK, ident.1 = "CD16+ NK cells_Permissive", ident.2 = "CD16+ NK cells_Non-permissive", verbose = TRUE, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
CD16positiveNK_cells_DE$Significant <- with(CD16positiveNK_cells_DE, ifelse(p_val_adj < 0.05, "Yes", "No"))
CD16positiveNK_cells_filtered <- subset(CD16positiveNK_cells_DE, Significant == "Yes")
Idents(CD16positiveNK) <- "majority_voting"
CD16positiveNK@meta.data$status <- factor(CD16positiveNK@meta.data$status, levels = c('Permissive', 'Non-permissive'))
genes1 <- c("CCL4L2", "KLRC2", "FKBP5", "CXCR4", "ADGRG1", "ADGRE5", "DTHD1", "TSC22D3", 
            "CCL4", "CCL3", "IFITM1", "KLRD1", "PRF1", "CCL5", "NKG7", "CST7", "ITGB2", "FKBP8", 
            "CAPZB", "S100A4", "HLA-B", "BST2", "B2M", "GZMM", "CD81", "CARD16", "IRF1", "GZMA", 
            "HLA-DQB1", "HLA-DRB1", "LTB", "GZMK", "CXCR3", "SPG20", "IL7R", "KLRG1", "HLA-DRA", 
            "XCL2", "TIGIT", "ITM2A", "CD169", "CASP1", "CSTB", "ISG15", "GSDMD", "CD97", "PDCD5", 
            "CD74", "PDCD6", "CASP4")

GroupHeatmap(CD16positiveNK, features = genes1,  group.by = "status", exp_method = "log2fc", 
             slot = "counts", nlabel = 50, features_fontsize = 4,
             heatmap_palette = "material-blue", 
             group_palcolor = list(c("#0000FF", "#D2042D")), cluster_rows = TRUE)

CD16positiveNK <- RunDEtest(srt = CD16positiveNK, group_by = "status", fc.threshold = 1)
VolcanoPlot(srt = CD16positiveNK, group_by = "status", nlabel = 20)
CD16positiveNK <- RunEnrichment(
  srt = CD16positiveNK, group_by = "status", db = c("GO_BP","MSigDB", "KEGG"), species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05")

EnrichmentPlot(srt = CD16positiveNK, group_by = "status", plot_type = "comparison", db = c("GO_BP", "MSigDB"), 
               id_use = c("GO:0006119", "GO:0009060", "GO:0045333", "GO:0042776", "GO:1902600", "GO:0002181",
                          "GO:0022613", "GO:0022618", "GO:0071826", "GO:0042254", "M5926", "M5913", "M5936", "M5950", "M5911",
                          "M5939", "M5950", "M5911", "M5946", "M5902"))

###############subsetting CD16 negative NK cells
CD16negNK <- subset(PversusNP, subset = majority_voting == "CD16- NK cells")
FeatureDimPlot(CD16negNK, features = "NCAM1", reduction = "UMAP")
CD16negNK <- subset(x = CD16negNK, subset = FCGR3A < 4.0)
CD16negNK <- subset(x = CD16negNK, subset = NCAM1 > 0.5)
CD16negNK <- SCTransform(CD16negNK, vars.to.regress = c("percent.mt", "percent.ribo", "percent.hb"), verbose = TRUE)
CD16negNK <- RunPCA(CD16negNK, verbose = FALSE)
CD16negNK <- RunUMAP(CD16negNK, dims = 1:50, verbose = FALSE)
CD16negNK <- FindNeighbors(CD16negNK, dims = 1:50, verbose = FALSE)
CD16negNK <- FindClusters(CD16negNK, verbose = FALSE)
CD16negNKmarkers <- FindAllMarkers(object = CD16negNK, min.pct = 0.25, logfc.threshold = 0.25)

#####DE expression
CD16negNK@meta.data$status <- factor(CD16negNK@meta.data$status, levels = c('Permissive', 'Non-permissive'))
FeatureDimPlot(yearstatus, features = "FCGR3A", reduction = "UMAP")
CD16negNK$celltype.stim <- paste(CD16negNK$majority_voting, CD16negNK$status, sep = "_")
DefaultAssay(object = CD16negNK) <- "SCT"
Idents(CD16negNK) <- "celltype.stim"
unique_titles <- unique(CD16negNK$celltype.stim)
print(unique_titles)
CD16negNK_cells_DE <- FindMarkers(CD16negNK, ident.1 = "CD16- NK cells_Permissive", ident.2 = "CD16- NK cells_Non-permissive", verbose = TRUE, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
CD16negNK_cells_DE$Significant <- with(CD16negNK_cells_DE, ifelse(p_val_adj < 0.05, "Yes", "No"))
CD16negNK_cells_filtered <- subset(CD16negNK_cells_DE, Significant == "Yes")
Idents(CD16negNK) <- "majority_voting"
genes1 <- c("CCL3L3", "KDM5D", "DKK3", "CSF1", "ID3", "CXCR3", "LTB", "CD52", "HLA-DRB1", 
            "HLA-DRA", "S100A6", "ISG15", "GZMM", "CCL5", "ITGB2", "KLRD1", 
            "FCER1G", "CCL4", "CTSW", "CCL3", "IFITM2", "IKZF3", "IFITM1", 
            "KLRC2", "FKBP5", "CCL3L1", "ADGRE5", "CCL4L2", "XCL2", "TNFRSF18")

GroupHeatmap(CD16negNK, features = genes1,  group.by = "status", exp_method = "log2fc", 
             slot = "counts", nlabel = 45, features_fontsize = 4,
             heatmap_palette = "material-blue", 
             group_palcolor = list(c("#0000FF", "#D2042D")), cluster_rows = TRUE)

CD16negNK@meta.data$status <- factor(CD16negNK@meta.data$status, levels = c('Permissive', 'Non-permissive'))

CD16negNK <- RunDEtest(srt = CD16negNK, group_by = "status", fc.threshold = 1)
VolcanoPlot(srt = CD16negNK, group_by = "status", nlabel = 20)
CD16negNK <- RunEnrichment(
  srt = CD16negNK, group_by = "status", db = c("GO_BP","MSigDB"), species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05")

EnrichmentPlot(srt = CD16negNK, group_by = "status", plot_type = "comparison", db = c("GO_BP", "MSigDB"), 
               id_use = c("GO:0009205", "GO:0009206", "GO:0009145", "GO:0009144", "GO:0009199", "GO:0002483", 
                          "GO:0048002", "GO:0001913", "GO:0019883", "GO:0001912", "M5950", "M5913", "M5911", "M5947", "M5924", "M5936", "M5891", "M5950", "M5919", "M5890"))

#subsetting CD16+ MONOCYTES
CD16monocytes <- subset(PversusNP, subset = majority_voting == "CD16+ Monocytes")
FeatureDimPlot(CD16monocytes, features = "FCGR3A", reduction = "UMAP")
CD16monocytes <- SCTransform(CD16monocytes, vars.to.regress = c("percent.mt", "percent.ribo", "percent.hb"), verbose = TRUE)
CD16monocytes <- RunPCA(CD16monocytes, verbose = FALSE)
CD16monocytes <- RunUMAP(CD16monocytes, dims = 1:50, verbose = FALSE)
CD16monocytes <- FindNeighbors(CD16monocytes, dims = 1:50, verbose = FALSE)
CD16monocytes <- FindClusters(CD16monocytes, verbose = FALSE)
CD16monocytesNKmarkers <- FindAllMarkers(object = CD16monocytes, min.pct = 0.25, logfc.threshold = 0.25)

#DE expression
CD16monocytes$celltype.stim <- paste(CD16monocytes$majority_voting, CD16monocytes$status, sep = "_")
DefaultAssay(object = CD16monocytes) <- "SCT"
Idents(CD16monocytes) <- "celltype.stim"
unique_titles <- unique(CD16monocytes$celltype.stim)
print(unique_titles)
CD16monocytes_cells_DE <- FindMarkers(CD16monocytes, ident.1 = "CD16+ Monocytes_Permissive", ident.2 = "CD16+ Monocytes_Non-permissive", verbose = TRUE, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
CD16monocytes_cells_DE$Significant <- with(CD16monocytes_cells_DE, ifelse(p_val_adj < 0.05, "Yes", "No"))
CD16monocytes_cells_filtered <- subset(CD16monocytes_cells_DE, Significant == "Yes")
Idents(CD16monocytes) <- "majority_voting"
CD16monocytes@meta.data$status <- factor(CD16monocytes@meta.data$status, levels = c('Permissive', 'Non-permissive'))
genes1 <- row.names(CD16monocytes_cells_filtered)
genes1 <- c("PECAM1", "TAP2.1", "LILRA3", "VIMP", "CECR1", "CD97", "LSMD1", "CASP3", "CTSL", "TNFAIP2", "CYBA", "HLA-A", "ADGRE2", "SARAF", "JPT1", "VSIR", "ADGRE5", "STMP1", "JAML", "TMEM176A", "FYB1", "RACK1", 
            "TYMP", "CASP3")

GroupHeatmap(CD16monocytes, features = genes1,  group.by = "status", exp_method = "log2fc", 
             slot = "counts", nlabel = 40, features_fontsize = 4,
             heatmap_palette = "material-blue", 
             group_palcolor = list(c("#0000FF", "#D2042D")), cluster_rows = TRUE)

CD16monocytes <- RunDEtest(srt = CD16monocytes, group_by = "status", fc.threshold = 1.0, assay = "SCT")
VolcanoPlot(srt = CD16monocytes, group_by = "status", nlabel = 20)
CD16monocytes <- RunEnrichment(
  srt = CD16monocytes, group_by = "status", db = c("GO_BP", "GO_CC", "GO_MF","MSigDB"), species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05")

EnrichmentPlot(srt = CD16monocytes, group_by = "status", plot_type = "comparison", db = c("GO_BP", "MSigDB"), 
               id_use = c("GO:0015986", "GO:0009206", "GO:0009145", "GO:0009201", "GO:0006754", "GO:2001198", 
                          "GO:0002478",  "GO:0097067", "GO:0019884", "GO:0042542", "M5921", "M5890", "M5913", "M5905", "M5950", "M5936", "M5926", "M5939", "M5915", "M5905"))

###############subsetting CD14+ MONOCYTES
CD14monocytes <- subset(PversusNP, subset = majority_voting == "CD14+ Monocytes")
FeatureDimPlot(CD14monocytes, features = "CD14", reduction = "UMAP")
CD14monocytes <- SCTransform(CD14monocytes, vars.to.regress = c("percent.mt", "percent.ribo", "percent.hb"), verbose = TRUE)
CD14monocytes <- RunPCA(CD14monocytes, verbose = FALSE)
CD14monocytes <- RunUMAP(CD14monocytes, dims = 1:50, verbose = FALSE)
CD14monocytes <- FindNeighbors(CD14monocytes, dims = 1:50, verbose = FALSE)
CD14monocytes <- FindClusters(CD14monocytes, verbose = FALSE)
CD14monocytesmarkers <- FindAllMarkers(object = CD14monocytes, min.pct = 0.25, logfc.threshold = 0.25)
CD14monocytes@meta.data$status <- factor(CD14monocytes@meta.data$status, levels = c('Permissive', 'Non-permissive'))

#DE expression
CD14monocytes$celltype.stim <- paste(CD14monocytes$majority_voting, CD14monocytes$status, sep = "_")
DefaultAssay(object = CD14monocytes) <- "SCT"
Idents(CD14monocytes) <- "celltype.stim"
unique_titles <- unique(CD14monocytes$celltype.stim)
print(unique_titles)
CD14monocytes_DE <- FindMarkers(CD14monocytes, ident.1 = "CD14+ Monocytes_Permissive", ident.2 = "CD14+ Monocytes_Non-permissive", verbose = TRUE, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
CD14monocytes_DE$Significant <- with(CD14monocytes_DE, ifelse(p_val_adj < 0.05, "Yes", "No"))
CD14monocytes_filtered <- subset(CD14monocytes_DE, Significant == "Yes")
Idents(CD14monocytes) <- "majority_voting"
genes1 <- row.names(CD14monocytes_filtered)
gene1 <- c("IRF1", "IRF7", "FYB", "IFI30", "TNFAIP2", "CD81", "FCGRT", "CD9", "VSIR", 
           "TCEB2", "TMBIM4", "RACK1", "TIMP1", "BAX", "SQRDL", "SARAF", "HLA-A", 
           "HLA-B", "HLA-DQB1", "SEPTIN7", "JAML", "CYBA", "CD52", "CD97", "FTL", 
           "ADA2", "GRK3", "TSC22D4", "LILRA3", "S100A9", "TMEM176B")
CD14monocytes@meta.data$status <- factor(CD14monocytes@meta.data$status, levels = c('Permissive', 'Non-permissive'))

GroupHeatmap(CD14monocytes, features = gene1,  group.by = "status", exp_method = "log2fc", 
             slot = "counts", nlabel = 39, features_fontsize = 4,
             heatmap_palette = "material-blue", 
             group_palcolor = list(c("#0000FF", "#D2042D")), cluster_rows = TRUE)

CD14monocytes <- RunDEtest(srt = CD14monocytes, group_by = "status", fc.threshold = 1)

CD14monocytes <- RunEnrichment(
  srt = CD14monocytes, group_by = "status", db = c("GO_BP", "MSigDB"), species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05")

EnrichmentPlot(srt = CD14monocytes, group_by = "status", plot_type = "comparison", db = c("GO_BP", "MSigDB"), 
               id_use = c("GO:0015986", "GO:0009206", "GO:0009145",	"GO:0009201", "GO:0006754", "GO:0031349", "GO:0045088", 
                          "GO:0002757",	"GO:0019882", "GO:0002764", "M5913", "M5921", "M5902", "M5911", "M5950", "M5936", "M5913", "M5915", "M5926", "M5939"))
