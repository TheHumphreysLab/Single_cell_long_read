# Load libraries
library(ggplot2)
library(dplyr)
library(Seurat)
library(presto)
library(future)
library(harmony)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize) # Needed for ComplexHeatmap()
library(scCustomize)
library(DTUrtle)

# Settings
options(future.globals.maxSize = 1e+11)

# Load data
srat1 <- readRDS("scWT1.rds")
srat2 <- readRDS("scWT2.rds")
srat3 <- readRDS("scWT3.rds")
srat4 <- readRDS("snWT1.rds")
srat5 <- readRDS("snWT2.rds")
srat6 <- readRDS("snWT3.rds")

# Add metadata
srat1 <- RenameCells(srat1, add.cell.id = "scWT1")
srat2 <- RenameCells(srat2, add.cell.id = "scWT2")
srat3 <- RenameCells(srat3, add.cell.id = "scWT3")
srat4 <- RenameCells(srat4, add.cell.id = "snWT1")
srat5 <- RenameCells(srat5, add.cell.id = "snWT2")
srat6 <- RenameCells(srat6, add.cell.id = "snWT3")
srat1@meta.data$intervention <- "none"
srat2@meta.data$intervention <- "none"
srat3@meta.data$intervention <- "none"
srat4@meta.data$intervention <- "none"
srat5@meta.data$intervention <- "none"
srat6@meta.data$intervention <- "none"

# Find variable features
srat1 <- FindVariableFeatures(srat1, selection.method = "vst", nfeatures = 1000)
srat2 <- FindVariableFeatures(srat2, selection.method = "vst", nfeatures = 1000)
srat3 <- FindVariableFeatures(srat3, selection.method = "vst", nfeatures = 1000)
srat4 <- FindVariableFeatures(srat4, selection.method = "vst", nfeatures = 1000)
srat5 <- FindVariableFeatures(srat5, selection.method = "vst", nfeatures = 1000)
srat6 <- FindVariableFeatures(srat6, selection.method = "vst", nfeatures = 1000)

# Find variable features that are common to both single-cell and single-nucleus datasets
genes_sc <- union(union(VariableFeatures(srat1), VariableFeatures(srat2)), VariableFeatures(srat3))
genes_sn <- union(union(VariableFeatures(srat4), VariableFeatures(srat5)), VariableFeatures(srat6))
genes <- intersect(genes_sc, genes_sn)

# Merge data
srat <- merge(srat1, y = c(srat2, srat3, srat4, srat5, srat6), project = "NULL")
srat$sample <- factor(srat$sample, levels = c("scWT1", "scWT2", "scWT3", "snWT1", "snWT2", "snWT3"))

# Remove unnecessary objects
remove(srat1, srat2, srat3, srat4, srat5, srat6)
gc()

# Remove unnecessary metadata
DefaultAssay(srat) <- "RNA"
srat@assays["SCT"] <- NULL
Idents(srat) <- "NULL"
srat@meta.data <- srat@meta.data[, !(colnames(srat@meta.data) %in% c("nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.8", "seurat_clusters"))]

# Clustering before integration
srat <- SCTransform(srat)
srat <- RunPCA(srat, features = genes)  # Running PCA based on variable genes that are common to both methods

# Integrate samples - Harmony integration
srat_harmony <- IntegrateLayers(srat, method = HarmonyIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "integrated_harmony", verbose = T)
srat_harmony <- RunUMAP(srat_harmony, dims = 1:30, reduction = "integrated_harmony", reduction.name = "umap_harmony")
srat_harmony <- FindNeighbors(srat_harmony, reduction = "integrated_harmony", dims = 1:30)
srat_harmony <- FindClusters(srat_harmony, resolution = 0.4)

# Check culstering after integration
DimPlot(srat_harmony, reduction = "umap_harmony", group.by = "sample", pt.size = .1)
DimPlot(srat_harmony, reduction = "umap_harmony", group.by = "method", pt.size = .1)

# Quality check by clusters - UMAP plot
FeaturePlot(srat_harmony, reduction = "umap_harmony", features = "nFeature_SCT", cols = c("grey", "red"))
FeaturePlot(srat_harmony, reduction = "umap_harmony", features = "nCount_SCT", cols = c("grey", "red"))
FeaturePlot(srat_harmony, reduction = "umap_harmony", features = "percent.mt", cols = c("grey", "red"))

# Quality check by clusters - violin plot
VlnPlot(srat_harmony, feature = "nFeature_SCT", group.by = "SCT_snn_res.0.4", pt.size = 0.3)
VlnPlot(srat_harmony, feature = "nCount_SCT", group.by = "SCT_snn_res.0.4", pt.size = 0.3)
VlnPlot(srat_harmony, feature = "percent.mt", group.by = "SCT_snn_res.0.4", pt.size = 0.3)

# Profile of clusters
table(srat_harmony@meta.data$SCT_snn_res.0.4)
features <- c("Slc5a2", "Slc5a12", "Slc22a6", "Slc22a13", "Havcr1", "Vcam1", "Top2a", "Cryab", "Slc4a11", "Psca", "Slc12a1", "Umod", "Slc12a3", "Slc8a1", "Aqp2", "Slc26a7", "Slc26a4", "Ncam1", "Nphs1", "Nphs2", "Flt1", "Emcn", "Col1a1", "Col1a2", "Pdgfra", "Dgkb", "Upk1b", "Ptprc")
Idents(srat_harmony) = srat_harmony$SCT_snn_res.0.4
levels(srat_harmony) <- rev(levels(srat_harmony))  # Code for visualization purposes
DotPlot(srat_harmony, features = features, cols = c("lightyellow","royalblue")) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
levels(srat_harmony) <- rev(levels(srat_harmony))  # Code for visualization purposes

# Checking DoubletFinder annotations
table(srat_harmony@meta.data$DoubletFinder)
table(srat_harmony@meta.data$sample, srat_harmony@meta.data$DoubletFinder)
table(srat_harmony@meta.data$SCT_snn_res.0.4, srat_harmony@meta.data$DoubletFinder)
DimPlot(srat_harmony, reduction = "umap_harmony", group.by = "DoubletFinder", shuffle = T)

# Clusters 11, 18, 20 and 22 are likely doublets because two or more marker genes are present. More than 30% of the cells in each cluster are annotated as doublets by DoubletFinder. Therefore, all cells in these clusters are assumed to be doublets and will be removed from the analysis. All cells annotated as doublets by DoubletFinder will be removed from the analysis as well.

# Removing doublets
srat_harmony_c <- subset(srat_harmony, srat_harmony$SCT_snn_res.0.4 %in% c("11", "18", "20", "22"), invert = T)
srat_harmony_c <- subset(srat_harmony_c, subset = srat_harmony_c$DoubletFinder == "Singlet")
srat_harmony_c@meta.data$oldGroup <- srat_harmony_c@meta.data$SCT_snn_res.0.4

# Tidy up the data
DefaultAssay(srat_harmony_c) <- "RNA"
srat_harmony_c@assays["SCT"] <- NULL
srat_harmony_c@reductions[c("pca", "umap_unintegrated", "integrated_harmony", "umap_harmony")] <- NULL
srat_harmony_c@graphs[c("SCT_nn", "SCT_snn")] <- NULL
Idents(srat_harmony_c) <- "NULL"
colnames_meta <- colnames(srat_harmony_c@meta.data)
cols_to_remove <- colnames_meta[grepl("^SCT_snn_res_", colnames_meta)]
srat_harmony_c@meta.data <- srat_harmony_c@meta.data[, !(colnames(srat_harmony_c@meta.data) %in% c("nCount_SCT", "nFeature_SCT", "seurat_clusters", cols_to_remove))]

# Run PCA before reintegration
srat_harmony_c <- SCTransform(srat_harmony_c)
srat_harmony_c <- RunPCA(srat_harmony_c, features = genes)  # Running PCA based on variable genes that are common to both methods

# Reintegrate samples - Harmony integration
srat_harmony_re <- IntegrateLayers(srat_harmony_c, method = HarmonyIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "integrated_harmony", verbose = T)
srat_harmony_re <- RunUMAP(srat_harmony_re, dims = 1:20, reduction = "integrated_harmony", reduction.name = "umap_harmony")
srat_harmony_re <- FindNeighbors(srat_harmony_re, reduction = "integrated_harmony", dims = 1:20)
srat_harmony_re <- FindClusters(srat_harmony_re, resolution = 0.3)

# Check culstering after reintegration
DimPlot(srat_harmony_re, reduction = "umap_harmony", group.by = "sample", pt.size = .1)
DimPlot(srat_harmony_re, reduction = "umap_harmony", group.by = "method", pt.size = .1)
Idents(srat_harmony_re) = srat_harmony_re$SCT_snn_res.0.3
table(Idents(srat_harmony_re))
DimPlot(srat_harmony_re, reduction = "umap_harmony", label = TRUE,  pt.size = .1)

# Quality check by clusters - UMAP plot
FeaturePlot(srat_harmony_re, reduction = "umap_harmony", features = "nFeature_SCT", cols = c("grey", "red"))
FeaturePlot(srat_harmony_re, reduction = "umap_harmony", features = "nCount_SCT", cols = c("grey", "red"))
FeaturePlot(srat_harmony_re, reduction = "umap_harmony", features = "percent.mt", cols = c("grey", "red"))

# Quality check by clusters - violin plot
VlnPlot(srat_harmony_re, feature = "nFeature_SCT", group.by = "SCT_snn_res.0.3", pt.size = 0.3)
VlnPlot(srat_harmony_re, feature = "nCount_SCT", group.by = "SCT_snn_res.0.3", pt.size = 0.3)
VlnPlot(srat_harmony_re, feature = "percent.mt", group.by = "SCT_snn_res.0.3", pt.size = 0.3)

# Profile of clusters
table(srat_harmony_re@meta.data$SCT_snn_res.0.3)
features <- c("Slc5a2", "Slc5a12", "Slc22a6", "Slc22a13", "Havcr1", "Vcam1", "Top2a", "Cryab", "Slc4a11", "Psca", "Slc12a1", "Umod", "Slc12a3", "Slc8a1", "Aqp2", "Slc26a7", "Slc26a4", "Ncam1", "Nphs1", "Nphs2", "Flt1", "Emcn", "Col1a1", "Col1a2", "Pdgfra", "Dgkb", "Upk1b", "Ptprc")
Idents(srat_harmony_re) = srat_harmony_re$SCT_snn_res.0.3
levels(srat_harmony_re) <- rev(levels(srat_harmony_re))  # Code for visualization purposes
DotPlot(srat_harmony_re, features = features, cols = c("lightyellow","royalblue")) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
levels(srat_harmony_re) <- rev(levels(srat_harmony_re))  # Code for visualization purposes

# UMAP with cell cluster annotation
Idents(srat_harmony_re) = srat_harmony_re$SCT_snn_res.0.3
new_cluster_ids_re <- c("PCT", "PST1", "TAL", "Endo", "PST1", "DCT", "PST2", "Fibro", "PCT", "CNT", "CD PC", "Immune", "DTL & ATL", "DCT/CNT", "CD ICB", "CD ICA", "Podo", "Low quality")
names(new_cluster_ids_re) <- levels(srat_harmony_re)
srat_harmony_re <- RenameIdents(srat_harmony_re, new_cluster_ids_re)
DimPlot(srat_harmony_re, reduction = "umap_harmony", label = TRUE) + NoLegend() + labs(x = "UMAP1", y = "UMAP2")
levels(srat_harmony_re) <- c("PCT", "PST1", "PST2", "DTL & ATL", "TAL", "DCT", "DCT/CNT", "CNT", "CD PC", "CD ICA", "CD ICB", "Fibro", "Podo", "Endo", "Immune", "Low quality")

# DotPlot with cell cluster annotation
features_re <- c("Slc5a12", "Slc22a6", "Slc22a13", "Cryab", "Slc12a1", "Umod", "Slc12a3", "Calb1", "Aqp2", "Slc26a7", "Slc26a4", "Pdgfra", "Nphs1", "Flt1", "Ptprc")
levels(srat_harmony_re) <- rev(levels(srat_harmony_re))  # Code for visualization purposes
DotPlot(srat_harmony_re, features = features_re, cols = c("lightyellow","royalblue")) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "Gene", y = "Cell type")
levels(srat_harmony_re) <- rev(levels(srat_harmony_re))  # Code for visualization purposes

# Violin plot with cell cluster annotation
VlnPlot(srat_harmony_re, feature = "nFeature_SCT", pt.size = 0.01) + labs(title = "Number of expressed genes", x = "Cell type") + NoLegend()

# Cluster 17 is an unknown cluster. FindMarkers function will be used to find out the identity of this cluster.

# Join layers
srat_harmony_pre <- JoinLayers(srat_harmony_re, assay = "RNA")
srat_harmony_pre <- JoinLayers(srat_harmony_pre, assay = "ISOFORM")

# Find markers
DefaultAssay(srat_harmony_pre) <- "SCT"
Idents(srat_harmony_pre) = srat_harmony_pre$SCT_snn_res.0.3
srat_harmony_pre <- PrepSCTFindMarkers(srat_harmony_pre, assay = "SCT", verbose = TRUE)
cluster_markers_pre <- list()
for (i in 0:17) {
  cluster_markers_pre[[paste0("cluster", i, ".pre.markers")]] <- FindMarkers(srat_harmony_pre, ident.1 = i)
  write.csv(cluster_markers_pre[[paste0("cluster", i, ".pre.markers")]], paste0("cluster", i, "premarkers.csv"))
}
cluster17vs5premarkers <- FindMarkers(srat_harmony_pre, ident.1 = 17, ident.2 = 5)
write.csv(cluster17vs5premarkers, "cluster17vs5premarkers.csv")
cluster17vs13premarkers <- FindMarkers(srat_harmony_pre, ident.1 = 17, ident.2 = 13)
write.csv(cluster17vs13premarkers, "cluster17vs13premarkers.csv")

# Cluster 17 is a cluster of low quality cells. Most of the differentially expressed genes are those expresed in DCT, CNT, TAL and other PTs. nFeatures values are mostly lower than 1,000. Therefore, cluster 17 is removed from the analysis.

# Removing low quality cells
srat_harmony_re <- subset(srat_harmony_re, srat_harmony_re$SCT_snn_res.0.3 %in% "17", invert = T)

# Join layers
srat_harmony_re <- JoinLayers(srat_harmony_re, assay = "RNA")
srat_harmony_re <- JoinLayers(srat_harmony_re, assay = "ISOFORM")

# Find markers
DefaultAssay(srat_harmony_re) <- "SCT"
Idents(srat_harmony_re) = srat_harmony_re$SCT_snn_res.0.3
srat_harmony_re <- PrepSCTFindMarkers(srat_harmony_re, assay = "SCT", verbose = TRUE)
cluster_markers_re <- list()
for (i in 0:16) {
  cluster_markers_re[[paste0("cluster", i, ".markers")]] <- FindMarkers(srat_harmony_re, ident.1 = i)
  write.csv(cluster_markers_re[[paste0("cluster", i, ".markers")]], paste0("cluster", i, "markers.csv"))
}

# Annotate clusters
new_cluster_ids_re <- c("PCT", "PST1", "TAL", "Endo", "PST1", "DCT", "PST2", "Fibro", "PCT", "CNT", "CD PC", "Immune", "DTL & ATL", "DCT/CNT", "CD ICB", "CD ICA", "Podo")
names(new_cluster_ids_re) <- levels(srat_harmony_re)
srat_harmony_re <- RenameIdents(srat_harmony_re, new_cluster_ids_re)
levels(srat_harmony_re) <- c("PCT", "PST1", "PST2", "DTL & ATL", "TAL", "DCT", "DCT/CNT", "CNT", "CD PC", "CD ICA", "CD ICB", "Fibro", "Podo", "Endo", "Immune")
srat_harmony_re@meta.data$celltype <- Idents(srat_harmony_re)

# UMAP with cell cluster annotation
levels(srat_harmony_re) <- c("PCT", "PST1", "TAL", "Endo", "DCT", "PST2", "Fibro", "CNT", "CD PC", "Immune", "DTL & ATL", "DCT/CNT", "CD ICB", "CD ICA", "Podo")  # Code for visualization purposes
DimPlot(srat_harmony_re, reduction = "umap_harmony", label = TRUE) + NoLegend() + labs(x = "UMAP1", y = "UMAP2")
DimPlot(srat_harmony_re, reduction = "umap_harmony", label = TRUE, split.by = "method") + NoLegend() + labs(x = "UMAP1", y = "UMAP2") + theme(strip.text = element_blank())
levels(srat_harmony_re) <- c("PCT", "PST1", "PST2", "DTL & ATL", "TAL", "DCT", "DCT/CNT", "CNT", "CD PC", "CD ICA", "CD ICB", "Fibro", "Podo", "Endo", "Immune")  # Code for visualization purposes

# DotPlot with cell cluster annotation
features_re <- c("Slc5a12", "Slc22a6", "Slc22a13", "Cryab", "Slc12a1", "Slc12a3", "Calb1", "Aqp2", "Slc26a7", "Slc26a4", "Pdgfra", "Nphs1", "Flt1", "Ptprc")
levels(srat_harmony_re) <- rev(levels(srat_harmony_re))  # Code for visualization purposes
DotPlot(srat_harmony_re, features = features_re, cols = c("lightyellow","royalblue")) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "Gene", y = "Cell type") + guides(color = guide_colorbar(title = "Average\nExpression"), size  = guide_legend(title = "Percent\nExpressed"))
levels(srat_harmony_re) <- rev(levels(srat_harmony_re))  # Code for visualization purposes

# Save data
saveRDS(srat_harmony_re, file = "single_cell_and_nucleus_integrated.rds")

# Copy data
srat <- srat_harmony_re

# Save separate data after running SCTransform on the ISOFORM assay
DefaultAssay(srat_harmony_re) <- "ISOFORM"
srat_harmony_re[["ISOFORM"]] <- split(srat_harmony_re[["ISOFORM"]], f = srat_harmony_re$sample)
srat_harmony_re <- SCTransform(srat_harmony_re, assay = "ISOFORM", new.assay.name = "SCT_i")
saveRDS(srat_harmony_re, "single_cell_and_nucleus_integrated_SCTi.rds")

# UMAP for Cldn10-201 and Cldn10-203
max(GetAssayData(srat_harmony_re, assay = "SCT_i", layer = "data")["ENSMUST00000047761", ])  # Returns 1.94591014905531
max(GetAssayData(srat_harmony_re, assay = "SCT_i", layer = "data")["ENSMUST00000100314", ])  # Returns 2.99573227355399
FeaturePlot(srat_harmony_re, reduction = "umap_harmony", features = "ENSMUST00000047761", cols = c("grey", "red"), order = TRUE) + labs(x = "UMAP1", y = "UMAP2") + theme(plot.title = element_blank()) + scale_color_gradientn(colors = c("grey", "red"), limits = c(0, 3))  # Cldn10-201
FeaturePlot(srat_harmony_re, reduction = "umap_harmony", features = "ENSMUST00000100314", cols = c("grey", "red"), order = TRUE) + labs(x = "UMAP1", y = "UMAP2") + theme(plot.title = element_blank()) + scale_color_gradientn(colors = c("grey", "red"), limits = c(0, 3))  # Cldn10-203

# Removing unnecessary data for visualizing isoform expression of cell type marker genes
DefaultAssay(srat) <- "ISOFORM"
srat@assays["SCT"] <- NULL

# Import GTF file
tx2gene <- import_gtf(gtf_file = "gencode.vM33.primary_assembly.annotation.gtf.gz")

# Remove period and following numbers in transcript_id
tx2gene$transcript_id <- gsub("\\..*", "", tx2gene$transcript_id)
head(tx2gene)

# Isoforms of conventional cell type marker genes with more than 100 counts
markers <- c("ENSMUST00000111026", "ENSMUST00000045972", "ENSMUST00000127822", "ENSMUST00000146664", "ENSMUST00000010250", "ENSMUST00000084797", "ENSMUST00000216755", "ENSMUST00000217475", "ENSMUST00000034562", "ENSMUST00000110494", "ENSMUST00000028630", "ENSMUST00000110495", "ENSMUST00000147095", "ENSMUST00000029876", "ENSMUST00000136266", "ENSMUST00000141336", "ENSMUST00000023752", "ENSMUST00000042221", "ENSMUST00000001253", "ENSMUST00000202186", "ENSMUST00000202681", "ENSMUST00000168162", "ENSMUST00000000476", "ENSMUST00000149086", "ENSMUST00000006825", "ENSMUST00000123880", "ENSMUST00000031653", "ENSMUST00000110529", "ENSMUST00000200895", "ENSMUST00000182283", "ENSMUST00000182138", "ENSMUST00000112036", "ENSMUST00000187586")

# Running SCTransform for genes including the cell type marker genes
srat <- SCTransform(srat, assay = "ISOFORM", residual.features = unique(c(VariableFeatures(srat), markers)))

# Rename cell clusters
srat$celltype2 <- paste(srat$method, srat$celltype, sep = "_")
srat$celltype2 <- factor(srat$celltype2, levels = c(
  "sc_PCT", "sn_PCT",
  "sc_PST1", "sn_PST1",
  "sc_PST2", "sn_PST2",
  "sc_DTL & ATL", "sn_DTL & ATL",
  "sc_TAL", "sn_TAL",
  "sc_DCT", "sn_DCT",
  "sc_DCT/CNT", "sn_DCT/CNT",
  "sc_CNT", "sn_CNT",
  "sc_CD PC", "sn_CD PC",
  "sc_CD ICA", "sn_CD ICA",
  "sc_CD ICB", "sn_CD ICB",
  "sc_Fibro", "sn_Fibro",
  "sc_Podo", "sn_Podo",
  "sc_Endo", "sn_Endo",
  "sc_Immune", "sn_Immune"
))
Idents(srat) <- srat$celltype2

# Average expression by groups
expression <- AverageExpression(srat, assays = "SCT", slot = "scale.data", features = markers, group.by = "celltype2")$SCT

# Convert the object to matrix
expression <- as.matrix(expression)
rownames(expression) <- tx2gene$transcript_name[match(rownames(expression), tx2gene$transcript_id)]

# Define condition
condition <- ifelse(grepl("sn", colnames(expression)), "sn", "sc")
celltype <- gsub("sn-|sc-", "", colnames(expression))

# Define color scale using Spectral palette
col_fun <- colorRamp2(seq(-2.5, 2.5, length.out = 11), rev(brewer.pal(11, "Spectral")))

# PCT marker genes
expressionsub <- expression[c("Slc5a12-201", "Slc5a12-202", "Slc5a12-203", "Slc5a12-205"), c("sc-PCT", "sn-PCT")]
conditionsub <- ifelse(grepl("sn", colnames(expressionsub)), "sn", "sc")
celltypesub <- gsub("sn-|sc-", "", colnames(expressionsub))
top_ann <- HeatmapAnnotation(Condition = conditionsub,
col = list(Condition = setNames(brewer.pal(8, "Set2")[1:2], c("sc", "sn")),
Celltype = setNames(colorRampPalette(brewer.pal(10, "Set3"))(length(unique(celltypesub))),unique(celltypesub))),
annotation_name_side = "left", show_annotation_name = FALSE)

Heatmap(expressionsub, col = col_fun, top_annotation = top_ann,
column_split = factor(celltypesub, levels = unique(celltypesub)), cluster_rows = FALSE, cluster_columns = FALSE,
show_column_names = FALSE, show_row_names = TRUE, row_names_side = "right", border = TRUE,
heatmap_legend_param = list(title = "Gene\nexpression", at = c(-4, 4), labels = c("Low", "High")))

# TAL marker genes
expressionsub <- expression[c("Slc12a1-201", "Slc12a1-202", "Slc12a1-203", "Slc12a1-204"), c("sc-TAL", "sn-TAL")]
conditionsub <- ifelse(grepl("sn", colnames(expressionsub)), "sn", "sc")
celltypesub <- gsub("sn-|sc-", "", colnames(expressionsub))
top_ann <- HeatmapAnnotation(Condition = conditionsub,
col = list(Condition = setNames(brewer.pal(8, "Set2")[1:2], c("sc", "sn")),
Celltype = setNames(colorRampPalette(brewer.pal(10, "Set3"))(length(unique(celltypesub))),unique(celltypesub))),
annotation_name_side = "left", show_annotation_name = FALSE)

Heatmap(expressionsub, col = col_fun, top_annotation = top_ann,
column_split = factor(celltypesub, levels = unique(celltypesub)), cluster_rows = FALSE, cluster_columns = FALSE,
show_column_names = FALSE, show_row_names = TRUE, row_names_side = "right", border = TRUE,
heatmap_legend_param = list(title = "Gene\nexpression", at = c(-4, 4), labels = c("Low", "High")))

# CNT marker genes
expressionsub <- expression[c("Calb1-201", "Calb1-202", "Calb1-203"), c("sc-CNT", "sn-CNT")]
conditionsub <- ifelse(grepl("sn", colnames(expressionsub)), "sn", "sc")
celltypesub <- gsub("sn-|sc-", "", colnames(expressionsub))
top_ann <- HeatmapAnnotation(Condition = conditionsub,
col = list(Condition = setNames(brewer.pal(8, "Set2")[1:2], c("sc", "sn")),
Celltype = setNames(colorRampPalette(brewer.pal(10, "Set3"))(length(unique(celltypesub))),unique(celltypesub))),
annotation_name_side = "left", show_annotation_name = FALSE)

Heatmap(expressionsub, col = col_fun, top_annotation = top_ann,
column_split = factor(celltypesub, levels = unique(celltypesub)), cluster_rows = FALSE, cluster_columns = FALSE,
show_column_names = FALSE, show_row_names = TRUE, row_names_side = "right", border = TRUE,
heatmap_legend_param = list(title = "Gene\nexpression", at = c(-4, 4), labels = c("Low", "High")))

# Fibroblasts marker genes
expressionsub <- expression[c("Pdgfra-201", "Pdgfra-202", "Pdgfra-207", "Pdgfra-208"), c("sc-Fibro", "sn-Fibro")]
conditionsub <- ifelse(grepl("sn", colnames(expressionsub)), "sn", "sc")
celltypesub <- gsub("sn-|sc-", "", colnames(expressionsub))
top_ann <- HeatmapAnnotation(Condition = conditionsub,
col = list(Condition = setNames(brewer.pal(8, "Set2")[1:2], c("sc", "sn")),
Celltype = setNames(colorRampPalette(brewer.pal(10, "Set3"))(length(unique(celltypesub))),unique(celltypesub))),
annotation_name_side = "left", show_annotation_name = FALSE)

Heatmap(expressionsub, col = col_fun, top_annotation = top_ann,
column_split = factor(celltypesub, levels = unique(celltypesub)), cluster_rows = FALSE, cluster_columns = FALSE,
show_column_names = FALSE, show_row_names = TRUE, row_names_side = "right", border = TRUE,
heatmap_legend_param = list(title = "Gene\nexpression", at = c(-4, 4), labels = c("Low", "High")))
