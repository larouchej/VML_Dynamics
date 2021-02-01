# Read in loom files from velocyto analysis of 10x sequencing results using Seurat package.
# QC, combine and perform batch correction using CCA, filter, calculate marker gene expression, annotate using scCATCH, save object, and make plotting dataframe.
# Perform regeneration/degeneration differential expression using MAST on each cell type.
# December 12, 2020
# Jacqueline Larouche

library(Seurat) #v3.2.2
library(SeuratDisk) #v0.0.0.9013
library(SeuratWrappers) #v0.3.0
library(dplyr) #v1.0.2
library(scCATCH) #v2.1
library(ggplot2) # v3.2.1
library(patchwork) # v1.0.0
library(tidyverse) # v1.3.0
library(tibble) # v3.0.1
library(RColorBrewer) # v1.1-2
library(viridis) # v0.5.1
library(scales) # v1.1.1
library(NCmisc) #v1.1.5
sessionInfo()

setwd("~/Documents/UMichigan/Aguilar_Lab/ResearchProjects/VMLTimecourse/Scanpy")
cells.combined <- readRDS('vml_cca_dec11.rds')
#====================================================================
# Read in loom files and merge
#====================================================================
ldat_d0_R1 <- ReadVelocity(file = "LoomFile/UninjuredR1.loom")
ldat_d0_R2 <- ReadVelocity(file = "LoomFile/UninjuredR2.loom")
ldat_d7_3mm_M <- ReadVelocity(file = "LoomFile/D73mmM.loom")
ldat_d7_3mm_F <- ReadVelocity(file = "LoomFile/D73mmF.loom")
ldat_d14_3mm <- ReadVelocity(file = "LoomFile/D143mm.loom")
ldat_d28_3mm <- ReadVelocity(file = "LoomFile/D283mm.loom")
ldat_d42_3mm <- ReadVelocity(file = "LoomFile/D423mmM.loom")
bm_d0_R2 <- as.Seurat(x = ldat_d0_R2)
bm_d0_R1 <- as.Seurat(x = ldat_d0_R1)
bm_d7_3mm_M <- as.Seurat(x = ldat_d7_3mm_M)
bm_d7_3mm_F <- as.Seurat(x = ldat_d7_3mm_F)
bm_d14_3mm <- as.Seurat(x = ldat_d14_3mm)
bm_d28_3mm <- as.Seurat(x = ldat_d28_3mm)
bm_d42_3mm <- as.Seurat(x = ldat_d42_3mm)

ldat_d7_2mm <- ReadVelocity(file = "LoomFile/D72mm.loom")
ldat_d14_2mm <- ReadVelocity(file = "LoomFile/D142mm.loom")
ldat_d28_2mm <- ReadVelocity(file = "LoomFile/D282mm.loom")
bm_d7_2mm <- as.Seurat(x = ldat_d7_2mm)
bm_d14_2mm <- as.Seurat(x = ldat_d14_2mm)
bm_d28_2mm <- as.Seurat(x = ldat_d28_2mm)

ds.list <- list(bm_d0_R1, bm_d0_R2, bm_d7_3mm_F, bm_d7_3mm_M, bm_d7_2mm, bm_d14_3mm, bm_d14_2mm, bm_d28_3mm, bm_d28_2mm, bm_d42_3mm)

bm_d0_R1[["orig.ident"]] <- "UninjuredR1"
bm_d0_R2[["orig.ident"]] <- "UninjuredR2"
bm_d7_3mm_F[["orig.ident"]] <- "D7.3mm.F"
bm_d7_3mm_M[["orig.ident"]] <- "D7.3mm.M"
bm_d7_2mm[["orig.ident"]] <- "D7.2mm"
bm_d14_3mm[["orig.ident"]] <- "D14.3mm"
bm_d14_2mm[["orig.ident"]] <- "D14.2mm"
bm_d28_3mm[["orig.ident"]] <- "D28.3mm"
bm_d28_2mm[["orig.ident"]] <- "D28.2mm"
bm_d42_3mm[["orig.ident"]] <- "D42.3mm"
#merge into one seurat object
cells.combined <- merge(bm_d0_R1, y= c(bm_d0_R2, bm_d7_3mm_F, bm_d7_3mm_M, bm_d7_2mm, bm_d14_3mm, bm_d14_2mm, bm_d28_3mm, bm_d28_2mm, bm_d42_3mm), add.cell.ids = c('d0', 'd0', 'd7.F.3mm', 'd7.M.3mm', 'd7.2mm', 'd14.3mm', 'd14.2mm', 'd28.3mm', 'd28.2mm', 'd42.3mm'), project = 'VMLTimecourse')
table(cells.combined$orig.ident)
# D14.2mm     D14.3mm     D28.2mm     D28.3mm     D42.3mm      D7.2mm    D7.3mm.F    D7.3mm.M 
# 6938        3224        5678        6824        7543        7343        5276        2764 
# UninjuredR1 UninjuredR2 
# 6485        6330

# Add other metadata for plotting later
Idents(cells.combined) <- 'orig.ident'
current.cluster.ids <- c("UninjuredR1","UninjuredR2", "D7.3mm.F", "D7.3mm.M", "D7.2mm", "D14.3mm", "D14.2mm", "D28.3mm", "D28.2mm", "D42.3mm")
defect.cluster.ids <- c("Uninjured", "Uninjured", "3mm", "3mm", "2mm", "3mm", "2mm", "3mm", "2mm", "3mm")
timepoint.cluster.ids <- c('0', '0', '7', '7', '7', '14', '14', '28', '28', '42')
condition.cluster.ids <- c("Uninjured","Uninjured", "D7.3mm", "D7.3mm", "D7.2mm", "D14.3mm", "D14.2mm", "D28.3mm", "D28.2mm", "D42.3mm")
cells.combined$condition <- plyr::mapvalues(x = cells.combined$orig.ident, from = current.cluster.ids, to = condition.cluster.ids)
cells.combined$defect <- plyr::mapvalues(x = cells.combined$orig.ident, from = current.cluster.ids, to = defect.cluster.ids)
cells.combined$timepoint <- plyr::mapvalues(x = cells.combined$orig.ident, from = current.cluster.ids, to = timepoint.cluster.ids)
condition.cluster.ids <- c("Uninjured","Uninjured", "D7.3mm", "D7.3mm", "D7.2mm", "D14.3mm", "D14.2mm", "D28.3mm", "D28.2mm", "D42.3mm")
cells.combined$condition <- plyr::mapvalues(x = cells.combined$orig.ident, from = current.cluster.ids, to = condition.cluster.ids)

cells.combined[["RNA"]] <- cells.combined[["spliced"]]
saveRDS(cells.combined, file = "VML_raw_Dec07.rds")
cells.combined <- readRDS('VML_raw_Dec07.rds')

#====================================================================
# CCA batch correction
#====================================================================
merged.list <- SplitObject(cells.combined, split.by = "orig.ident")
merged.list <- lapply(X = merged.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 10000)
})
merged.anchors <- FindIntegrationAnchors(object.list = merged.list, dims = 1:30, k.filter = 200, k.score = 30, anchor.features = 2000)
cells.combined <- IntegrateData(anchorset = merged.anchors, dims = 1:30)
saveRDS(cells.combined, 'vml_cca_dec11.rds')

DefaultAssay(cells.combined) <- "integrated"
cells.combined <- ScaleData(cells.combined, verbose = TRUE)
cells.combined <- RunPCA(cells.combined, verbose = TRUE)
cells.combined <- RunUMAP(cells.combined, reduction = "pca", dims = 1:30)
cells.combined <- FindNeighbors(cells.combined, reduction = "pca", dims = 1:30)
cells.combined <- FindClusters(cells.combined, resolution = 0.5)
p1 <- DimPlot(cells.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(cells.combined, reduction = "umap", label = TRUE)
p1+p2

#====================================================================
# Annotate using scCATCH and markergene expression
#====================================================================

clu_markers <- findmarkergenes(cells.combined,
                               species = 'Mouse',
                               cluster = 'All',
                               match_CellMatch = TRUE,
                               cancer = NULL,
                               tissue = 'Muscle',
                               cell_min_pct = 0.25,
                               logfc = 0.25,
                               pvalue = 0.05)
clu_ann <- scCATCH(clu_markers$clu_markers,
                   species = 'Mouse',
                   cancer = NULL,
                   tissue = 'Muscle')
write.table(clu_ann, file = 'scCATCH_celltype_annotations.csv', sep = ',')

cells.combined.markers <- FindAllMarkers(cells.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
top10 <- cells.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.table(cells.combined.markers, file = 'SeuratOutput/markergenes_celltype.csv', sep = ',')

# label celltypes
Idents(cells.combined) <- 'seurat_clusters'
current.cluster.ids <- 0:26
celltype.ids <- c('Endothelial_Stem', 'Macrophage', 'FAP_Pro-Remodeling', 'Tenocyte', 'Erythroblast', 'Endothelial_Stem', 'Satellite_Cell', 'Smooth_Muscle', 'Endothelial_Vein', 'Endothelial_Capillary', 'Monocyte', 'T/NK_Cell', 'FAP_Matrix', 'Dendritic', 'Pericyte', 'Endothelial_Artery', 'Neutrophil', 'B_Cell', 'Erythroblast', 'FAP_Stem', 'Neutrophil', 'FAP_Adipogenic', 'Myonuclei', 'Lymph_Vessel', 'Neural_Progenitor', 'Neural', 'Smooth_Muscle')
cells.combined$Celltype <- plyr::mapvalues(x = cells.combined$seurat_clusters, from = current.cluster.ids, to = celltype.ids)
DimPlot(cells.combined, reduction = "umap", group.by = "Celltype", label = FALSE)

#====================================================================
# Defect-associated differential expression for each cell type
#====================================================================
# Adapted from Dulken et al. 2019
## Append defect to celltype label and add to metadata
cellDefect <- paste0(cells.combined@meta.data$Celltype, "_", cells.combined@meta.data$defect)
m <- data.frame("Celltype_Defect" = cellDefect)
rownames(m) <- rownames(cells.combined@meta.data)
cells.combined <- AddMetaData(object = cells.combined, metadata = m)
DimPlot(cells.combined, reduction = "umap", group.by = "Celltype_Defect", label = FALSE)

# Assign as main identity
Idents(cells.combined) <- "Celltype_Defect"

CELLTYPES <- unique(cells.combined@meta.data$Celltype)

mast_list <- vector(mode="list", length = 15)
names(mast_list) = CELLTYPES

# Compare 2mm/3mm for each celltype, save each matrix in a list. 
# This step takes a few hours on 1 core
for (CELLTYPE in CELLTYPES) {
  print(CELLTYPE)
  # Find  cluster marker genes
  cells.mast.de <- FindMarkers(object = cells.combined,
                               ident.1 = paste0(CELLTYPE, "_3mm"),
                               ident.2 = paste0(CELLTYPE, "_2mm"),
                               only.pos = FALSE, 
                               min.pct = 0,
                               logfc.threshold = 0,
                               test.use = "MAST")
  cells.mast.de <- as.data.frame(cells.mast.de)
  cells.mast.de <- rownames_to_column(cells.mast.de, var = "gene")
  cells.mast.de$celltype <- rep(CELLTYPE, dim(cells.mast.de)[1])
  mast_list[[CELLTYPE]] <- cells.mast.de
}

# Combine all matrices into one dataframe
mast_df <- data.frame()
for (CELLTYPE in CELLTYPES) {
  print(CELLTYPE)
  mast_df  <- rbind(mast_df, mast_list[[CELLTYPE]])
}
dim(mast_df); head(mast_df) # 159516 by 7

# Save differential expression results
write.csv(mast_df, "mast_df.csv")

