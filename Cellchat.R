library(dplyr)
library(Seurat)
library(SeuratObject)
library(sctransform)
library(ggplot2)
library(Matrix)
library(CellChat)
library(readr)
library(ComplexHeatmap)
library(patchwork)

###############################################12M results#################################################################
#Create separate CellChat objects for each condition from merged Seurat object
Permissive <- subset(yearstatus, subset = status == "Permissive")
Nonpermissive <- subset(yearstatus, subset = status == "Non-permissive")
DefaultAssay(object = Permissive) <- "SCT"
DefaultAssay(object = Nonpermissive) <- "SCT"

# Preparing data input and metadata for permissive
data.input_Permissive <- GetAssayData(Permissive, assay = "SCT", slot = "data")
Idents(Permissive) <- Permissive@meta.data[["majority_voting"]]
labels_Permissive <- Idents(Permissive)
meta_Permissive <- data.frame(labels = labels_Permissive, row.names = names(labels_Permissive))

# Preparing data input and metadata for nonpermissive
data.input_Nonpermissive <- GetAssayData(Nonpermissive, assay = "SCT", slot = "data")
Idents(Nonpermissive) <- Nonpermissive@meta.data[["majority_voting"]]
labels_Nonpermissive <- Idents(Nonpermissive)
meta_Nonpermissive <- data.frame(labels = labels_Nonpermissive, row.names = names(labels_Nonpermissive))

# Create cellchat
cellChat_Permissive <- createCellChat(object = data.input_Permissive, meta = meta_Permissive, group.by = "labels")
cellChat_Nonpermissive <- createCellChat(object = data.input_Nonpermissive, meta = meta_Nonpermissive, group.by = "labels")

cellChat_Permissive <- createCellChat(object = Permissive, group.by = "majority_voting", assay = "SCT")
cellChat_Nonpermissive <- createCellChat(object = Nonpermissive, group.by = "majority_voting", assay = "SCT")

# use human info
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB)
cellChat_Permissive@DB <- CellChatDB.use
cellChat_Nonpermissive@DB <- CellChatDB.use

# Analysis for permissive and nonpermissive samples
cellChat_Permissive <- subsetData(cellChat_Permissive)
future::plan("multisession", workers = 4)
cellChat_Permissive <- identifyOverExpressedGenes(cellChat_Permissive)
cellChat_Permissive <- identifyOverExpressedInteractions(cellChat_Permissive)
cellChat_Permissive <- computeCommunProb(cellChat_Permissive)
cellChat_Permissive <- filterCommunication(cellChat_Permissive, min.cells = 10)
cellChat_Permissive <- computeCommunProbPathway(cellChat_Permissive)
cellChat_Permissive <- netAnalysis_computeCentrality(cellChat_Permissive, slot.name = "netP")
cellChat_Permissive <- aggregateNet(cellChat_Permissive)

cellChat_Nonpermissive <- subsetData(cellChat_Nonpermissive)
cellChat_Nonpermissive <- identifyOverExpressedGenes(cellChat_Nonpermissive)
cellChat_Nonpermissive <- identifyOverExpressedInteractions(cellChat_Nonpermissive)
cellChat_Nonpermissive <- computeCommunProb(cellChat_Nonpermissive)
cellChat_Nonpermissive <- filterCommunication(cellChat_Nonpermissive, min.cells = 10)
cellChat_Nonpermissive <- computeCommunProbPathway(cellChat_Nonpermissive)
cellChat_Nonpermissive <- netAnalysis_computeCentrality(cellChat_Nonpermissive, slot.name = "netP")
cellChat_Nonpermissive <- aggregateNet(cellChat_Nonpermissive)

object.list <- list(Nonpermissive = cellChat_Nonpermissive, Permissive = cellChat_Permissive)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(object.list, file = "object.list.rds")
saveRDS(cellchat, file = "cellchat_merged.rds")
cellchat <- read_rds("cellchat_merged.rds")
object.list <- read_rds("object.list.rds")

#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2)) 
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#Circle plot showing differential number of interactions or interaction strength among different cell populations across two datasets
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, 
                          arrow.width = 0, arrow.size = 0, 
                          remove.isolate = TRUE,
                          color.edge = c("#0000FF", "#D2042D"), 
                          alpha.edge = 0.5, vertex.label.cex = 0.7)

netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", arrow.width = 0, arrow.size = 0, remove.isolate = TRUE,color.edge = c("#0000FF", "#D2042D"), 
                          alpha.edge = 0.5, vertex.label.cex = 0.7)


#Identify cell populations with significant changes in sending or receiving signals
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax) + scale_x_continuous(limits = c(0,7)) + scale_y_continuous(limits = c(0,12))
}

patchwork::wrap_plots(gg) + coord_fixed()


#Compare the overall information flow of each signaling pathway or ligand-receptor pair
rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE, cutoff.pvalue = 0.00001)

###############################################PRB results#################################################################
permissive0 <- subset(yearstatus, subset = status == "Permissive_0")
nonpermissive0 <- subset(yearstatus, subset = status == "Non-permissive_0")
DefaultAssay(object = permissive0) <- "SCT"
DefaultAssay(object = nonpermissive0) <- "SCT"
CellDimPlot(nonpermissive0, group.by = "majority_voting")
CellDimPlot(permissive0, group.by = "majority_voting")

# Create cellchat
cellChat_Permissive0 <- createCellChat(object = permissive0, group.by = "majority_voting", assay = "SCT")
cellChat_Nonpermissive0 <- createCellChat(object = nonpermissive0, group.by = "majority_voting", assay = "SCT")

# use human info
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB)
cellChat_Permissive0@DB <- CellChatDB.use
cellChat_Nonpermissive0@DB <- CellChatDB.use

# Analysis for permissive and nonpermissive samples
cellChat_Permissive0 <- subsetData(cellChat_Permissive0)
future::plan("multisession", workers = 4)
cellChat_Permissive0 <- identifyOverExpressedGenes(cellChat_Permissive0)
cellChat_Permissive0 <- identifyOverExpressedInteractions(cellChat_Permissive0)
cellChat_Permissive0 <- computeCommunProb(cellChat_Permissive0)
cellChat_Permissive0 <- filterCommunication(cellChat_Permissive0, min.cells = 10)
cellChat_Permissive0 <- computeCommunProbPathway(cellChat_Permissive0)
cellChat_Permissive0 <- netAnalysis_computeCentrality(cellChat_Permissive0, slot.name = "netP")
cellChat_Permissive0 <- aggregateNet(cellChat_Permissive0)

cellChat_Nonpermissive0 <- subsetData(cellChat_Nonpermissive0)
cellChat_Nonpermissive0 <- identifyOverExpressedGenes(cellChat_Nonpermissive0)
cellChat_Nonpermissive0 <- identifyOverExpressedInteractions(cellChat_Nonpermissive0)
cellChat_Nonpermissive0 <- computeCommunProb(cellChat_Nonpermissive0)
cellChat_Nonpermissive0 <- filterCommunication(cellChat_Nonpermissive0, min.cells = 10)
cellChat_Nonpermissive0 <- computeCommunProbPathway(cellChat_Nonpermissive0)
cellChat_Nonpermissive0 <- netAnalysis_computeCentrality(cellChat_Nonpermissive0, slot.name = "netP")
cellChat_Nonpermissive0 <- aggregateNet(cellChat_Nonpermissive0)

object.list0 <- list(Nonpermissive0 = cellChat_Nonpermissive0, Permissive0 = cellChat_Permissive0)
cellchat0 <- mergeCellChat(object.list0, add.names = names(object.list0))
saveRDS(object.list0, file = "object.list0.rds")
saveRDS(cellchat0, file = "cellchat_merged0.rds")
object.list0 < - read_rds("object.list0.rds")
cellchat0 <- read_rds("cellchat_merged0.rds")


#interaction number
gg3 <- compareInteractions(cellchat0, show.legend = F, group = c(1,2), color.use = c("#D2042D", "#0000FF")) 
gg4 <- compareInteractions(cellchat0, show.legend = F, group = c(1,2), color.use = c("#D2042D", "#0000FF"), measure = "weight")
gg3 + gg4


#Circle plot showing the number of interactions or interaction strength among different cell populations across multiple datasets
weight.max <- getMaxWeight(object.list0, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list0)) {
  netVisual_circle(object.list0[[i]]@net$count, arrow.size = 0, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list0)[i]))
}

#Identify cell populations with significant changes in sending or receiving signals
num.link <- sapply(object.list0, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg1 <- list()
for (i in 1:length(object.list0)) {
  gg1[[i]] <- netAnalysis_signalingRole_scatter(object.list0[[i]], title = names(object.list0)[i], weight.MinMax = weight.MinMax) + scale_x_continuous(limits = c(0,3.5)) + scale_y_continuous(limits = c(0,4.5))
}

patchwork::wrap_plots(plots = gg1) + coord_fixed()


#Compare the overall information flow of each signaling pathway or ligand-receptor pair
rankNet(cellchat0, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE, cutoff.pvalue = 0.001)
