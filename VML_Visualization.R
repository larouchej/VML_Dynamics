# Read in Seurat objects created in VML_analysis.R and generate figure plots.
# Jacqueline Larouche
# December 12, 2020

library(Seurat) #v3.1.2
library(dplyr) # v0.8.3
library(cowplot) # v1.0.0
library(ggplot2) # v3.2.1
library(patchwork) # v1.0.0
library(tidyverse) # v1.3.0
library(tibble) # v3.0.1
library(RColorBrewer) # v1.1-2
library(viridis) # v0.5.1
library(scales) # v1.1.1
library(NCmisc) #v1.1.5
library(nichenetr) #v1.0.0
library(SeuratWrappers) #v0.3.0
library(SeuratDisk) #v0.0.0.9013
library(dittoSeq) #v1.0.2
sessionInfo()

setwd('~/Documents/UMichigan/Aguilar_Lab/ResearchProjects/VMLTimecourse/Scanpy')

# Read in Seurat objects
cells.combined <- readRDS('vml_cca_dec11.rds')
muscs <- readRDS('vml_musc_integrated.rds')
lymphocytes <- readRDS("vml_lymphocytes.rds")

#=====================================================================
# Figure 1
#=====================================================================
Endothelial_Stem <- '#5F9C04'
Endothelial_Vein <- '#139C04'
Endothelial_Capillary <- '#82CEA0'
Endothelial_Artery <- '#049C8D'
Lymph_Vessel <- '#04139C'

Smooth_Muscle <- '#C1DA26'
Pericyte <- '#26C1DA'

Macrophage <- '#F92B39'
Monocyte <- '#EB2BF9'
T_Cell <- '#F9842B'
Dendritic <- '#9C1B64'
Neutrophil <- '#842BF9'
B_Cell <- '#F9EB2B'

FAP_Remodeling <- '#5993FF'
FAP_Stem <- '#7259FF'
FAP_Adipogenic <- '#FF59E6'
FAP_Matrix <- '#C559FF'
Tenocyte <- '#FF7259'
Neural_Progenitor <- '#FF5993'
Neural <- '#BF436E'

Satellite_Cell <- '#00FFAF'
Myonuclei <- '#0050FF'

Erythroblast <- '#BABABA'

colors <- c(Monocyte, Dendritic, Macrophage, Neutrophil, T_Cell, B_Cell,  Satellite_Cell, Myonuclei, Endothelial_Stem, Endothelial_Vein, Endothelial_Capillary, Endothelial_Artery, Pericyte,  Smooth_Muscle, Lymph_Vessel, FAP_Stem, FAP_Remodeling, FAP_Matrix, FAP_Adipogenic, Tenocyte, Neural_Progenitor, Neural, Erythroblast)

# Figure 1b: UMAP colored by celltype
Idents(cells.combined) <- 'Celltype'
my_levels <- c('Monocyte', 'Dendritic', 'Macrophage', 'Neutrophil',  'B_Cell', 'T/NK_Cell', 'Satellite_Cell', 'Myonuclei', 'Endothelial_Stem', 'Endothelial_Vein', 'Endothelial_Capillary', 'Endothelial_Artery', 'Pericyte',  'Smooth_Muscle', 'Lymph_Vessel', 'FAP_Stem', 'FAP_Pro-Remodeling', 'FAP_Matrix', 'FAP_Adipogenic', 'Tenocyte', 'Neural_Progenitor', 'Neural',  'Erythroblast')
Idents(cells.combined) <- factor(Idents(cells.combined), levels= my_levels)
tiff("SeuratOutput/umap.celltype.colored.tiff", width = 500, height = 500)
DimPlot(cells.combined, reduction = "umap", label = FALSE, cols = colors) + NoLegend()
dev.off()

# Figure 1c-e: Stacked bar graphs of cell abundance
write.table(table(cells.combined$Celltype, cells.combined$timepoint, cells.combined$defect), file = 'SeuratOutput/abundance.csv', sep = ',')
# Graphs made in Excel

# Figure 1f: differential gene expression within celltypes (adapted from Dulken et al. 2019)
df <- read_csv("mast_df.csv")
df$z <- p.to.Z(df$p_val) * sign(df$avg_logFC)
df$z.adj <- p.to.Z(df$p_val_adj) * sign(df$avg_logFC)
df <- df[sample(nrow(df)), ]
df$celltype_factor <- factor(df$celltype,  levels=my_levels, ordered=T)

# Make custom color column to facilitate grey coloring by threshold.
col <- colors[df$celltype_factor]
col[df$p_val_adj > 0.05] <- "#D3D3D3" # grey
df$col <- as.factor(col)

q <- ggplot(df, aes(x = celltype_factor, y = z, color = col)) +
  geom_jitter(width = 0.40, alpha = .55, size = 1) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(axis.title.y = element_text(size = 20, face = "plain")) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(plot.title = element_text(size=20, face = "plain")) +
  theme(legend.position="none") +
  labs(y = "Z-score") +
  scale_color_manual(values = levels(df$col)) +
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed")
tiff("SeuratOutput/mast.degs.tiff", width = 800, height = 400)
q
dev.off()

#=====================================================================
# Figure 2
#=====================================================================
# Figure 2b: Nichenet results for Neutrophils (bonafide ligand/receptor interaction predictions)
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))

Idents(cells.combined) <- 'timepoint'
DefaultAssay(cells.combined) <- 'spliced'
seuratObj <- subset(cells.combined, idents = c('7' ,'14'))
seuratObj <- ScaleData(seuratObj)
seuratObj[['integrated']] <- NULL
Idents(seuratObj) <- 'Celltype'

nichenet_output <- nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver = 'Neutrophil', 
  condition_colname = 'defect', condition_oi = '3mm', condition_reference = '2mm', 
  sender = 'all', filter_top_ligands = TRUE, top_n_ligands = 50, expression_pct = 0.05,
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = 'mouse')

p <- nichenet_output$ligand_receptor_heatmap_bonafide
nichenet_output$top_ligands_filtered <- c('Bmp2', 'Apoe', 'Gzmb', 'Ifng', 'Il6', 'Pgf', 'Plau', 'Tnf', 'Tgfb2', 'Tgfb1', 'Il1b', 'Il1a', 'Dll4', 'Dll1', 'Ccl3', 'Ccl4')
nichenet_output$top_receptors_filtered <- c("Bmpr1a", 'Bmpr2', 'Ldlr', 'Sorl1', 'Igf2r', 'Ifngr1', 'Ifngr2', 'Il6ra', 'Nrp1', 'Plaur', 'Tnfrsf21', 'Tnfrsf1b', 'Ltbr', 'Tnfrsf1a', 'Tgfbr1', 'Tgfbr2', 'Acvrl1', 'Il1r2', 'Notch1', 'Notch2', 'Ccr1')
p2 <- nichenet_output$ligand_receptor_heatmap_bonafide
b <- p2$data$y %in% nichenet_output$top_ligands_filtered
c <- p2$data$x %in% nichenet_output$top_receptors_filtered
p2$data <- p2$data[(b&c),]
p3 <- nichenet_output$ligand_activity_target_heatmap

# Figure 2c: Violin plot of CCL5 expression by celltype
Idents(cells.combined) <- 'Celltype'
my_levels <- c('Monocyte', 'Dendritic', 'Macrophage', 'Neutrophil',  'B_Cell', 'T/NK_Cell', 'Satellite_Cell', 'Myonuclei', 'Endothelial_Stem', 'Endothelial_Vein', 'Endothelial_Capillary', 'Endothelial_Artery', 'Pericyte',  'Smooth_Muscle', 'Lymph_Vessel', 'FAP_Stem', 'FAP_Pro-Remodeling', 'FAP_Matrix', 'FAP_Adipogenic', 'Tenocyte', 'Neural_Progenitor', 'Neural',  'Erythroblast')
Idents(cells.combined) <- factor(Idents(cells.combined), levels= my_levels)

tiff("SeuratOutput/vlnplot.ccl5.tiff", width = 800, height = 250)
VlnPlot(cells.combined, features = c("Ccl5"), pt.size = 0)
dev.off()

# Figure 2f: Violin plot of NK cell markers by defect size
DimPlot(lymphocytes, reduction = "umap", label = FALSE) 
# Isolate NK cells, re-cluster
nk.cells <- subset(lymphocytes, idents = c('1', '2'))
nk.cells <- NormalizeData(nk.cells)
nk.cells <- FindVariableFeatures(nk.cells)
all.genes <- rownames(nk.cells)
nk.cells <- ScaleData(nk.cells, features = all.genes, split.by = 'batch')
nk.cells <- RunPCA(nk.cells, features = VariableFeatures(object = nk.cells))
ElbowPlot(nk.cells, ndims = 40)
nk.cells <- FindNeighbors(nk.cells, dims = 1:40)
nk.cells <- FindClusters(nk.cells, resolution = 1)
nk.cells <- RunUMAP(nk.cells, dims = 1:40)
DimPlot(nk.cells, reduction = "umap", pt.size = 5)
Idents(nk.cells) <- 'seurat_clusters'
VlnPlot(nk.cells, features = c('Ncr1','Nkg7'), pt.size = 0)
nk.cells <- subset(nk.cells, idents = c('5', '6'), invert = TRUE) #First round filtering
nk.cells <- subset(nk.cells, idents = c('4'), invert = TRUE) #First round filtering

Idents(nk.cells) <- 'defect'
VlnPlot(nk.cells, features = c('Ncr1', 'Klrk1', 'Ccl5', 'Klrc1', 'Ifng', 'Ltb'), pt.size = 0, group.by = 'defect', idents = c('2mm', '3mm'))

#=====================================================================
# Figure 5
#=====================================================================
# Figure 5a: NicheNet results for satellite cells
peach <- '#f7766d'
teal <- '#088992'

ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))

seuratObj <- subset(seuratObj, subset = timepoint == "7")
seuratObj <- ScaleData(seuratObj)
seuratObj[['integrated']] <- NULL
Idents(seuratObj) <- 'Celltype'

# Perform NicheNet Analysis
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver = 'Satellite_Cell', 
  condition_colname = 'defect', condition_oi = '3mm', condition_reference = '2mm', 
  sender = 'all', 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = 'mouse')

nichenet_output$top_ligands_filtered <- c('Tgfb1', 'Edn1', 'Il1b', 'Ptgs2', 'Il1a', 'Hmgb1', 'Tgfb2', 'Tgfb3', 'Dll1', 'Il10')
p2 <- nichenet_output$ligand_target_heatmap
b <- p2$data$y %in% nichenet_output$top_ligands_filtered
p2$data <- p2$data[b,]
p3 <- nichenet_output$ligand_activity_target_heatmap

ECM <- c("Col1a1", "Col1a2", "Col3a1", "Fn1", "Postn", "Vim", "Vcam1", "Lcn2", "Sptbn1", "Dcn", "Zeb2", "Fosb", "Adamts5")
Cytoskeleton <- c("Des", "Capn2", "Dst", "Utrn")
Stress <- c("C1qb", "Hsp90aa1")
Inflammation <- c("Itgb1", 'S100a8', 'S100a9', 'Il1r2', 'Socs3')
Cellcycle <- c("Klf6", "Egr1",  "Zfp36l1", "Fgfr1", "Igfbp5", "Bmpr2", "Ddx5", "Prox1", "Acvr2a", "Ash1l", "Nedd4", "Zbtb20")
other <- c("Ar","Lmna","Hp","Tcf4", "Pgyrp1", "Mbnl2", "Gnas", "Med13", "Tubb2b")
ordered_response <- c(ECM, Cytoskeleton, Stress, Inflammation, Cellcycle, other)
c <- p2$data$x %in% ordered_response
p2$data <- p2$data[c,]
p2$data$x <- factor(p2$data$x, levels = ordered_response)

tiff("SeuratOutput/nichenet.musc.d7.response.tiff", width = 680, height = 350)
p2 + scale_fill_gradient2(low = "whitesmoke",  high = teal) + xlab("") + ylab("") + theme(legend.position = "bottom")
dev.off()

pearson_heatmap <- data.frame(nichenet_output$top_ligands_filtered, nichenet_output$ligand_activities$pearson[nichenet_output$ligand_activities$test_ligand %in% nichenet_output$top_ligands_filtered])
colnames(pearson_heatmap) <- c('Ligand', 'Pearson')
pearson_heatmap$Ligand <- factor(pearson_heatmap$Ligand, levels = pearson_heatmap$Ligand[order(pearson_heatmap$Pearson)])

p3 <- ggplot(pearson_heatmap, aes(x = 1, y = Ligand, fill = Pearson)) + geom_tile(size = 1, color= "white") + scale_fill_gradient2(low = "whitesmoke",  high = peach, limits = c(0.05, 0.12), midpoint = 0.06)

tiff("SeuratOutput/nichenet.musc.d7.pearson.tiff", width = 275, height = 350)
p3 + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background = element_blank(), axis.ticks = element_blank()) +
  ggtitle("Ligand Activity") +
  theme(axis.title.y = element_text(size = 10, face = "plain", family = 'Arial')) +
  theme(axis.text.y = element_text(size = 10, family = 'Arial')) +
  theme(plot.title = element_text(size=15, face = "plain", family = 'Arial')) + 
  theme(legend.title = element_text(size=10, face = "plain", family = 'Arial'), legend.position = "bottom") +
  labs(y = "Ligand")
dev.off()

#=====================================================================
# Supplemental Figure 2 (RNA-Seq QC, Markergene Expression)
#=====================================================================
# Supplemental Figure 2a-b: UMAPs with nCounts and nFeatures
tiff("SeuratOutput/umap.qc.tiff", width = 300, height = 500)
p1 <- FeaturePlot(cells.combined, features = c("nFeature_RNA"), min.cutoff = "q10", max.cutoff = "q90")
p2 <- FeaturePlot(cells.combined, features = c("nCount_RNA"), min.cutoff = "q10", max.cutoff = "q90")
p1+p2 + plot_layout(ncol = 1)
dev.off()

# Supplemental Figure 2c: UMAP colored by dataset 
tiff("SeuratOutput/umap.dataset.tiff", width = 500, height = 400)
DimPlot(cells.combined, reduction = "umap", group.by = "orig.ident", label = FALSE)
dev.off()

# Figure 2d: Dot Plots with marker gene overlays
Idents(cells.combined) <- 'Celltype'
p_endothelial <- DotPlot(cells.combined, 
                         features = c('Cd36', 'Kdr', 'Myh11', 'Myl9', 'Vwf', 'Lrg1', 'Rgs5', 'Pdgfrb',
                                      'Alpl', 'Hey1', 'Pecam1', 'Lyve1'), 
                         idents = c('Endothelial_Stem', 'Endothelial_Vein', 'Endothelial_Artery',
                                    'Endothelial_Capillary', 'Pericyte', 'Smooth_Muscle', 'Lymph_Vessel'),
                         cols = c("lightgrey", "red")) + NoLegend()
p_mesenchymal <- DotPlot(cells.combined, 
                         features = c('Col3a1', 'Pdgfra', 'Mmp14', 'Col1a1', 'Scx', 'Tnmd', 'Gdf10', 'Dpp4',
                                      'Thy1', 'Cd34', 'Apod', 'Dcn'),
                         idents = c('FAP_Pro-Remodeling', 'Tenocyte', 'FAP_Matrix', 'FAP_Stem', 'FAP_Adipogenic',
                                    'Neural_Progenitor', 'Neural'),
                         cols = c("lightgrey", "red")) + NoLegend()
p_immune <- DotPlot(cells.combined, 
                    features = c('Ptprc', 'Adgre1', 'Cd68','Ms4a3', 'Cd3d', 'Nkg7', 'Cd74', 'H2-Aa', 'S100a8',
                                 'Pglyrp1', 'Cd79a', 'Ighm'), 
                    idents = c('Macrophage', 'Myeloid_Progenitor', 'T/NK_Cell', 'Dendritic', 'Neutrophil',
                               'B_Cell'),
                    cols = c("lightgrey", "red")) + NoLegend()

p_myogenic <- DotPlot(cells.combined, 
                      features = c('Pax7', 'Myf5', 'Cd34', 'Cdh15', 'Sdc4', 'Vcam1', 'Myod1', 'Cdkn1a', 'Myog', 
                                   'Myh4', 'Myh1', 'Tnnt3'), 
                      idents = c('Satellite_Cell', 'Myonuclei'),
                      cols = c("lightgrey", "red")) + NoLegend()

tiff("SeuratOutput/Dotplot.markergenes.tiff", width = 850, height = 750)
p_endothelial + p_mesenchymal + p_immune + p_myogenic +  plot_layout(ncol = 1)
dev.off()

#=====================================================================
# Supplemental Figure 3 (NK cells)
#=====================================================================
# Supplemental Figure 3a: Feature plot of lymphocyte marker genes
Idents(cells.combined) <- 'Celltype'
lymphocytes <- subset(cells.combined, idents = c("T/NK_Cell"), invert = FALSE)
DefaultAssay(lymphocytes) <- 'spliced'

merged.list <- SplitObject(lymphocytes, split.by = "orig.ident")
merged.list <- lapply(X = merged.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})
merged.anchors <- FindIntegrationAnchors(object.list = merged.list, dims = 1:10, k.filter = 20, k.score = 30, anchor.features = 2000)
lymphocytes <- IntegrateData(anchorset = merged.anchors, dims = 1:30)
DefaultAssay(lymphocytes) <- 'integrated'

lymphocytes <- subset(lymphocytes, idents = c('6', '10', '11', '12', '13', '14'), invert = TRUE) #round 1 filtering
lymphocytes <- ScaleData(lymphocytes, verbose = TRUE, features = all.genes)#, split.by = 'orig.ident')
lymphocytes <- RunPCA(lymphocytes, verbose = TRUE)
lymphocytes <- RunUMAP(lymphocytes, reduction = "pca", dims = 1:30)
lymphocytes <- FindNeighbors(lymphocytes)
lymphocytes <- FindClusters(lymphocytes, resolution = 0.2)
DimPlot(lymphocytes, reduction = "umap", label = TRUE)
FeaturePlot(lymphocytes, features = c('Ncr1','Nkg7', 'Cd3d', 'Cd4', 'Cd8a', 'Foxp3'), min.cutoff = 'q10', max.cutoff = 'q90')
VlnPlot(lymphocytes, features = c('Ncr1','Gzma', 'Cd3d', 'Cd4', 'Cd8a', 'Foxp3'), group.by = 'seurat_clusters', pt.size = 0)

current.cluster.ids <- c("0", "1", "2", "3", "4", "5")
new.cluster.ids <- c("Cytotoxic_Tcell", "NKcell", "Helper_Tcell", "Cytotoxic_Tcell", "Tcell", "NKTcell")
lymphocytes$Celltype <- plyr::mapvalues(x = lymphocytes$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
Idents(lymphocytes) <- 'Celltype'
VlnPlot(lymphocytes, features = c("Ccl5"), pt.size = 0, group.by = 'Celltype') + NoLegend() + theme(text = element_text(family = 'Arial')) + stat_summary(fun = median, geom='point', size = 25, colour = "black", shape = 95)

# Supplemental Figure 3b: Stacked bar graph of cell abundances
t1 <- table(lymphocytes$Celltype, lymphocytes$condition)
t2 <- table(cells.combined$condition)
t3 <- sweep(t1, 2, t2, FUN = '/')
df <- data.frame(t3)
names(df) <- c("cluster", "condition", "frequency")
df$condition <- factor(df$condition,levels = c("Uninjured", "D7.2mm", "D7.3mm", "D14.2mm", "D14.3mm", "D28.2mm", "D28.3mm", "D42.3mm"))
df2 <- df
df$defect <- c(2,2,2,2,2,3,3,3,3,3,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3)
df2$defect <- c(2,2,2,2,2,3,3,3,3,3,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,3,3,3,3,3,2,2,2,2,2)
df_2 <- df2[df2$defect == 2,]
df_2$timepoint <- c(rep(14, 5),rep(28, 5), rep(7,5), rep(0,5))
df_3 <- df[df$defect == 3,]
df_3$timepoint <- c(rep(14, 5),rep(28, 5), rep(42,5), rep(7,5), rep(0,5))
ratios <- NA
for (a in c(14, 28, 7, 0)) {
  b <- c(df_3[df_3$timepoint == a,]$frequency/df_2[df_2$timepoint == a,]$frequency)
  ratios <- c(ratios, b)
}
df_2$ratios <- ratios
df_2$timepoint <- factor(df_2$timepoint,levels = c("0", "7", "14", "28"))
df_2$cluster <- factor(df_2$cluster,levels = c("Cytotoxic_Tcell", "Helper_Tcell", "Tcell", "NKTcell", "NKcell"))

p <- ggplot(df_2, aes(x = timepoint, y = cluster, fill = ratios)) + 
  geom_tile(colour="white",size=0.5) + scale_fill_viridis(na.value = "grey90") +   theme(text = element_text(family = 'Arial')) + theme(panel.background = element_blank(), plot.background=element_blank())

tiff("SeuratOutput/heatmap_lymphocyte_ratios.tiff", height = 400, width = 300)
p 
dev.off()

bp3<- ggplot(df[df$defect == c(3),], aes(x=condition, y=frequency, fill=cluster)) +
  geom_bar(width = 1, stat = "identity") + 
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) + ggtitle("Lymphocyte Populations 3mm") +
  theme(axis.title.y = element_text(size = 10, face = "plain")) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 7.5)) +
  theme(plot.title = element_text(size=10, face = "plain")) +
  scale_y_continuous(limits = c(0,0.12)) + NoLegend() + theme(text = element_text(family = 'Arial'))
bp2<- ggplot(df2[df2$defect == c(2),], aes(x=condition, y=frequency, fill=cluster)) +
  geom_bar(width = 1, stat = "identity") + 
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) + ggtitle("Lymphocyte Populations 2mm") +
  theme(axis.title.y = element_text(size = 10, face = "plain")) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(plot.title = element_text(size=10, face = "plain")) +
  scale_y_continuous(limits = c(0,1)) + theme(text = element_text(family = 'Arial')) + NoLegend()
bp2+bp3

#Supplemental Figure 3c: Lymphocyte markergene heatmap
Idents(lymphocytes) <- "Celltype"
seuratObj.markers.all <- FindAllMarkers(lymphocytes, verbose = TRUE, logfc.threshold = 1, only.pos = TRUE)
heatmap_features <- seuratObj.markers.all$gene
tiff("SeuratOutput/heatmap_lymphocyte_markergenes.tiff", width = 550, height = 800)
DoHeatmap(lymphocytes, features = heatmap_features, slot = 'scale.data') + scale_fill_gradient2(low = "blue", mid = "white", high = "red") + theme(text = element_text(family = 'Arial'))
dev.off()

#=====================================================================
# Supplemental Figures 5 (Satellite Cells)
#=====================================================================
# Supplemental Figure 5b: Violin plots of MuSC differentially expressed genes
Idents(muscs) <- 'timepoint'
my_levels <- c('0','7', '14', '28')
Idents(muscs) <- factor(Idents(muscs), levels= my_levels)
VlnPlot(muscs, features = c('S100a8', 'S100a9', 'Gamt','Cebpb'), pt.size = 0, split.by = 'defect', idents = c('D7', 'D14', 'D28'), ncol = 2) + theme(text = element_text(family = 'Arial')) + NoLegend()

# Supplemental Figure 5c: Fusion Gene Expression
Idents(muscs)<- "condition"
my_levels <- c("Uninjured", "D7.2mm", "D7.3mm", "D14.2mm", "D14.3mm", "D28.2mm", "D28.3mm", "D42.3mm")
Idents(muscs) <- factor(Idents(muscs), levels= my_levels)

fusion.filtered <- c('Adgrb3', 'Cacna1s', 'Capn2','Cav3', 'Ccl8', 'Cd9', 'Cd53', 'Cd81', 'Cdon', 'Cflar', 'Cxcl9', 'Cxcl10', 'Cxcl12', 'Dyrk1b', 'Ehd1', 'Ehd2', 'Fer1l5','Gsk3b', 'Il4', 'Il4ra', 'Ins2', 'Mapk14', 'Myh9', 'Mymk', 'Mymx', 'Myod1', 'Myof', 'Myog', 'Neo1', 'Nfatc2', 'Nphs1', 'Pitx2', 'Plekho1', 'Ptgfrn','Tanc1')
d7 <- FindMarkers(muscs, features = fusion.filtered, assay = 'spliced', ident.1 = 'D7.3mm', ident.2 = 'D7.2mm', logfc.threshold = 0, min.pct = 0, only.pos = F, return.thresh = 0, slot = "data")
d14 <- FindMarkers(muscs, features = fusion.filtered, assay = 'spliced', ident.1 = 'D14.3mm', ident.2 = 'D14.2mm', logfc.threshold = 0, min.pct = 0, only.pos = F, return.thresh = 0, slot = "data")
d28 <- FindMarkers(muscs, features = fusion.filtered, assay = 'spliced', ident.1 = 'D28.3mm', ident.2 = 'D28.2mm', logfc.threshold = 0, min.pct = 0, only.pos = F, return.thresh = 0, slot = "data")

timepoint <- c(rep('D7', 30), rep('D14', 30), rep('D28', 30))
genes <- c(rownames(d7), rownames(d14), rownames(d28))
avg_logFC <- c(d7$avg_logFC, d14$avg_logFC, d28$avg_logFC)
df <- data.frame(timepoint, genes, logfc)
df$timepoint <- factor(df$timepoint, levels=c('D7', 'D14', 'D28'))

p <- ggplot(df, aes(x=timepoint, y=genes, fill=avg_logFC)) +
  geom_tile(colour="white", size=0.5) + 
  scale_fill_gradient2(low="blue", mid = 'white', high="red") + 
  labs(x = "", y = "") +
  theme(text = element_text(family = 'Arial', size = 15)) +
  theme(panel.background = element_blank(), plot.background=element_blank()) +
  theme(axis.text.x = element_text(color = 'black')) +
  theme(axis.text.y = element_text(color = 'black')) +
  theme(legend.position="bottom")

tiff('SeuratOutput/heatmap_fusiongenes.tiff', width = 800, height = 400)
p
dev.off()
#=====================================================================
# Supplemental Figures 6 (TGFb)
#=====================================================================
# Supplemental Figure 6a: NicheNet Target Gene Expression
receiver = "Satellite_Cell"
cells.combined.d7 <- subset(cells.combined, idents = '7')
Idents(cells.combined.d7)<-'Celltype'
seurat_obj_receiver7 = subset(cells.combined.d7, idents = receiver)
Idents(seurat_obj_receiver7) <- 'defect'
DE_table_receiver = FindAllMarkers(object = seurat_obj_receiver7, return.thresh = 0.05, logfc.threshold = 0.1)
DE_subset <- filter(DE_table_receiver, gene %in% ordered_response)

p <- ggplot(DE_subset, aes(x = cluster, y = gene)) + 
  geom_tile(aes(fill = avg_logFC), color = 'white') + 
  scale_fill_gradient2(low="blue", mid = 'white', high="red") + 
  labs(x = "", y = "") +
  theme(text = element_text(family = 'Arial', size = 15)) +
  theme(panel.background = element_blank(), plot.background=element_blank()) +
  theme(axis.text.x = element_text(color = 'black')) +
  theme(axis.text.y = element_text(color = 'black')) +
  theme(legend.position="bottom")
p

tiff("SeuratOutput/heatmap.logFC.musc.nichenet.targets.tiff", width = 300, height = 800)
p
dev.off()

# Supplemental Figure 6b: Violin plots of TGFb1 expression
Idents(cells.combined) <- 'timepoint'
cells.combined.d7 <- subset(cells.combined, idents = '7')
cells.combined.d14 <- subset(cells.combined, idents = '14')
cells.combined.d28 <- subset(cells.combined, idents = '28')
Idents(cells.combined.d7) <- 'Celltype'
Idents(cells.combined.d14) <- 'Celltype'
Idents(cells.combined.d28) <- 'Celltype'
my_levels <- c('Monocyte', 'Dendritic', 'Macrophage', 'Neutrophil',  'B_Cell', 'T/NK_Cell', 'Satellite_Cell', 'Myonuclei', 'Endothelial_Stem', 'Endothelial_Vein', 'Endothelial_Capillary', 'Endothelial_Artery', 'Pericyte',  'Smooth_Muscle', 'Lymph_Vessel', 'FAP_Stem', 'FAP_Pro-Remodeling', 'FAP_Matrix', 'FAP_Adipogenic', 'Tenocyte', 'Neural_Progenitor', 'Neural',  'Erythroblast')
Idents(cells.combined.d7) <- factor(Idents(cells.combined.d7), levels= my_levels)
Idents(cells.combined.d14) <- factor(Idents(cells.combined.d14), levels= my_levels)
Idents(cells.combined.d28) <- factor(Idents(cells.combined.d28), levels= my_levels)

p1 <- VlnPlot(cells.combined.d7, features = 'Tgfb1', pt.size = 0, split.by = 'defect', split.plot = TRUE)
p2 <- VlnPlot(cells.combined.d14, features = 'Tgfb1', pt.size = 0, split.by = 'defect', split.plot = TRUE)
p3 <- VlnPlot(cells.combined.d28, features = 'Tgfb1', pt.size = 0, split.by = 'defect', split.plot = TRUE)
p1+p2+p3
