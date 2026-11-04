# Load libraries
library(ggplot2)
library(dplyr)
library(Seurat)
library(presto)
library(future)
library(harmony)
library(ggalluvial)
library(DTUrtle)

# Settings
options(future.globals.maxSize = 1e+11)

# Load data
srat1 <- readRDS("snWT1.rds")
srat2 <- readRDS("snWT2.rds")
srat3 <- readRDS("snWT3.rds")
srat4 <- readRDS("snIRI1.rds")
srat5 <- readRDS("snIRI2.rds")
srat6 <- readRDS("snIRI3.rds")

# Add metadata
srat1 <- RenameCells(srat1, add.cell.id = "snWT1")
srat2 <- RenameCells(srat2, add.cell.id = "snWT2")
srat3 <- RenameCells(srat3, add.cell.id = "snWT3")
srat4 <- RenameCells(srat4, add.cell.id = "snIRI1")
srat5 <- RenameCells(srat5, add.cell.id = "snIRI2")
srat6 <- RenameCells(srat6, add.cell.id = "snIRI3")
srat1@meta.data$intervention <- "none"
srat2@meta.data$intervention <- "none"
srat3@meta.data$intervention <- "none"
srat4@meta.data$intervention <- "IRI"
srat5@meta.data$intervention <- "IRI"
srat6@meta.data$intervention <- "IRI"

# Merge data and remove unnecessary data
srat <- merge(srat1, y = c(srat2, srat3, srat4, srat5, srat6), project = "NULL")
srat$intervention <- factor(srat$intervention, levels = c("none", "IRI"))
srat$sample <- factor(srat$sample, levels = c("snWT1", "snWT2", "snWT3", "snIRI1", "snIRI2", "snIRI3"))
DefaultAssay(srat) <- "ISOFORM"
srat@assays["SCT"] <- NULL
Idents(srat) <- "NULL"
srat@meta.data <- srat@meta.data[, !(colnames(srat@meta.data) %in% c("nCount_SCT", "nFeature_SCT", "seurat_clusters"))]
remove(srat1, srat2, srat3, srat4, srat5, srat6)
gc()

# Clustering before integration
srat <- SCTransform(srat, assay = "ISOFORM", new.assay.name = "SCT_i")
srat <- RunPCA(srat, assay = "SCT_i", reduction.name = "pca_i")

# Integrate samples - Harmony integration
srat_harmony <- IntegrateLayers(srat, method = HarmonyIntegration, normalization.method = "SCT", orig.reduction = "pca_i", new.reduction = "integrated_harmony_i", verbose = T)
srat_harmony <- RunUMAP(srat_harmony, dims = 1:30, reduction = "integrated_harmony_i", reduction.name = "umap_harmony_i")
srat_harmony <- FindNeighbors(srat_harmony, reduction = "integrated_harmony_i", dims = 1:30)
srat_harmony <- FindClusters(srat_harmony, resolution = 0.8)

# Check culstering after integration
DimPlot(srat_harmony, reduction = "umap_harmony_i", group.by = "sample", pt.size = .1)
DimPlot(srat_harmony, reduction = "umap_harmony_i", group.by = "intervention", pt.size = .1)

# Quality check by clusters - UMAP plot
FeaturePlot(srat_harmony, reduction = "umap_harmony_i", features = "nFeature_SCT_i", cols = c("grey", "red"))
FeaturePlot(srat_harmony, reduction = "umap_harmony_i", features = "nCount_SCT_i", cols = c("grey", "red"))
FeaturePlot(srat_harmony, reduction = "umap_harmony_i", features = "percent.mt", cols = c("grey", "red"))

# Quality check by clusters - violin plot
VlnPlot(srat_harmony, feature = "nFeature_SCT_i", group.by = "SCT_i_snn_res.0.8", pt.size = 0.3)
VlnPlot(srat_harmony, feature = "nCount_SCT_i", group.by = "SCT_i_snn_res.0.8", pt.size = 0.3)
VlnPlot(srat_harmony, feature = "percent.mt", group.by = "SCT_i_snn_res.0.8", pt.size = 0.3)

# Profile of clusters
table(srat_harmony@meta.data$SCT_i_snn_res.0.8)
features <- c("ENSMUST00000111026", "ENSMUST00000010250", "ENSMUST00000084797", "ENSMUST00000109223", "ENSMUST00000029574", "ENSMUST00000068031", "ENSMUST00000099362", "ENSMUST00000023265", "ENSMUST00000110494", "ENSMUST00000034218", "ENSMUST00000029876", "ENSMUST00000023752", "ENSMUST00000042221", "ENSMUST00000001253", "ENSMUST00000202681", "ENSMUST00000221176", "ENSMUST00000006825", "ENSMUST00000031653", "ENSMUST00000057767", "ENSMUST00000183301")  # Corresponds to "Slc5a12-202", "Slc22a6-201", "Slc22a13-201", "Havcr1-203", "Vcam1-201", "Top2a-201", "Slc4a11-201", "Psca-201", "Slc12a1-202", "Slc12a3-201", "Calb1-201", "Aqp2-201", "Slc26a7-201", "Slc26a4-201", "Pdgfra-208", "Dgkb-205", "Nphs1-201", "Flt1-202", "Upk1b-201", "Ptprc-207"
Idents(srat_harmony) = srat_harmony$SCT_i_snn_res.0.8
levels(srat_harmony) <- rev(levels(srat_harmony))
DotPlot(srat_harmony, features = features, cols = c("lightyellow","royalblue")) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
levels(srat_harmony) <- rev(levels(srat_harmony))

# Checking DoubletFinder annotations
table(srat_harmony@meta.data$DoubletFinder)
table(srat_harmony@meta.data$sample, srat_harmony@meta.data$DoubletFinder)
table(srat_harmony@meta.data$SCT_i_snn_res.0.8, srat_harmony@meta.data$DoubletFinder)
DimPlot(srat_harmony, reduction = "umap_harmony_i", group.by = "DoubletFinder", shuffle = T)

# Clusters 13, 17 and 22 are likely doublets because two or more marker genes are present. More than 70% of the cells in each cluster are annotated as doublets by DoubletFinder. Therefore, all cells in these clusters are assumed to be doublets and will be removed from the analysis. All cells annotated as doublets by DoubletFinder will be removed from the analysis as well.

# Removing doubletes
srat_harmony_c <- subset(srat_harmony, srat_harmony$SCT_i_snn_res.0.8 %in% c("13", "17", "22"), invert = T)
srat_harmony_c <- subset(srat_harmony_c, subset = srat_harmony_c$DoubletFinder == "Singlet")

# Tidy up the data
DefaultAssay(srat_harmony_c) <- "ISOFORM"
srat_harmony_c@assays["SCT_i"] <- NULL
srat_harmony_c@reductions[c("pca_i", "umap_unintegrated", "integrated_harmony_i", "umap_harmony_i")] <- NULL
srat_harmony_c@graphs[c("SCT_i_nn", "SCT_i_snn")] <- NULL
Idents(srat_harmony_c) <- "NULL"
colnames_meta <- colnames(srat_harmony_c@meta.data)
cols_to_remove <- colnames_meta[grepl("^SCT_i_snn_res_", colnames_meta)]
srat_harmony_c@meta.data <- srat_harmony_c@meta.data[, !(colnames(srat_harmony_c@meta.data) %in% c("nCount_SCT_i", "nFeature_SCT_i", "seurat_clusters", cols_to_remove))]

# Run PCA before reintegration
srat_harmony_c <- SCTransform(srat_harmony_c, assay = "ISOFORM", new.assay.name = "SCT_i")
srat_harmony_c <- RunPCA(srat_harmony_c, assay = "SCT_i", reduction.name = "pca_i")

# Reintegrate samples - Harmony integration
srat_harmony_re <- IntegrateLayers(srat_harmony_c, method = HarmonyIntegration, normalization.method = "SCT_i", orig.reduction = "pca_i", new.reduction = "integrated_harmony_i", verbose = T)
srat_harmony_re <- RunUMAP(srat_harmony_re, dims = 1:30, reduction = "integrated_harmony_i", reduction.name = "umap_harmony_i")
srat_harmony_re <- FindNeighbors(srat_harmony_re, reduction = "integrated_harmony_i", dims = 1:30)
srat_harmony_re <- FindClusters(srat_harmony_re, resolution = 0.4)

# Check culstering after reintegration
DimPlot(srat_harmony_re, reduction = "umap_harmony_i", group.by = "sample", pt.size = .1)
DimPlot(srat_harmony_re, reduction = "umap_harmony_i", group.by = "intervention", pt.size = .1)
Idents(srat_harmony_re) = srat_harmony_re$SCT_i_snn_res.0.4
table(Idents(srat_harmony_re))
DimPlot(srat_harmony_re, reduction = "umap_harmony_i", label = TRUE,  pt.size = .1)

# Quality check by clusters - UMAP plot
FeaturePlot(srat_harmony_re, reduction = "umap_harmony_i", features = "nFeature_SCT_i", cols = c("grey", "red"))
FeaturePlot(srat_harmony_re, reduction = "umap_harmony_i", features = "nCount_SCT_i", cols = c("grey", "red"))
FeaturePlot(srat_harmony_re, reduction = "umap_harmony_i", features = "percent.mt", cols = c("grey", "red"))

# Quality check by clusters - violin plot
VlnPlot(srat_harmony_re, feature = "nFeature_SCT_i", group.by = "SCT_i_snn_res.0.4", pt.size = 0.3)
VlnPlot(srat_harmony_re, feature = "nCount_SCT_i", group.by = "SCT_i_snn_res.0.4", pt.size = 0.3)
VlnPlot(srat_harmony_re, feature = "percent.mt", group.by = "SCT_i_snn_res.0.4", pt.size = 0.3)

# Profile of clusters
table(srat_harmony_re@meta.data$SCT_i_snn_res.0.4)
features <- c("ENSMUST00000111026", "ENSMUST00000010250", "ENSMUST00000084797", "ENSMUST00000109223", "ENSMUST00000029574", "ENSMUST00000068031", "ENSMUST00000099362", "ENSMUST00000023265", "ENSMUST00000110494", "ENSMUST00000034218", "ENSMUST00000029876", "ENSMUST00000023752", "ENSMUST00000042221", "ENSMUST00000001253", "ENSMUST00000202681", "ENSMUST00000221176", "ENSMUST00000006825", "ENSMUST00000031653", "ENSMUST00000057767", "ENSMUST00000183301")  # Corresponds to "Slc5a12-202", "Slc22a6-201", "Slc22a13-201", "Havcr1-203", "Vcam1-201", "Top2a-201", "Slc4a11-201", "Psca-201", "Slc12a1-202", "Slc12a3-201", "Calb1-201", "Aqp2-201", "Slc26a7-201", "Slc26a4-201", "Pdgfra-208", "Dgkb-205", "Nphs1-201", "Flt1-202", "Upk1b-201", "Ptprc-207"
levels(srat_harmony_re) <- rev(levels(srat_harmony_re))
DotPlot(srat_harmony_re, features = features, cols = c("lightyellow","royalblue")) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
levels(srat_harmony_re) <- rev(levels(srat_harmony_re))

# Join layers
srat_harmony_re <- JoinLayers(srat_harmony_re, assay = "RNA") %>% JoinLayers(assay = "ISOFORM")

# Find markers
Idents(srat_harmony_re) = srat_harmony_re$SCT_i_snn_res.0.4
for (i in 1:6) {
    slot(object = srat_harmony_re@assays$SCT_i@SCTModel.list[[i]], name = "umi.assay") <- "ISOFORM"
}
srat_harmony_re <- PrepSCTFindMarkers(srat_harmony_re, assay = "SCT_i", verbose = TRUE)
cluster_markers <- list()
for (i in 0:16) {
  cluster_markers[[paste0("cluster", i, ".markers")]] <- FindMarkers(srat_harmony_re, ident.1 = i)
  write.csv(cluster_markers[[paste0("cluster", i, ".markers")]], paste0("cluster", i, "markers.csv"))
}

# Annotate clusters
new.cluster.ids <- c("PCT", "PST1", "TAL", "Endo", "CD PC", "DCT", "Injured PT", "PST2", "Fibro", "Proliferating PT", "Immune", "DTL & ATL", "CD ICA", "CD ICB", "CNT1", "CNT2", "Podo")
names(new.cluster.ids) <- levels(srat_harmony_re)
srat_harmony_re <- RenameIdents(srat_harmony_re, new.cluster.ids)
levels(srat_harmony_re) <- c("PCT", "PST1", "PST2", "Injured PT", "Proliferating PT", "DTL & ATL", "TAL", "DCT", "CNT1", "CNT2", "CD PC", "CD ICA", "CD ICB", "Fibro", "Podo", "Endo", "Immune")
srat_harmony_re@meta.data$celltype_isoform <- Idents(srat_harmony_re)

# UMAP with cell cluster annotation
levels(srat) <- c("PCT", "PST1", "TAL", "Endo", "CD PC", "DCT", "Injured PT", "PST2", "Fibro", "Proliferating PT", "Immune", "DTL & ATL", "CD ICA", "CD ICB", "CNT1", "CNT2", "Podo")
DimPlot(srat, reduction = "umap_harmony_i", label = TRUE) + NoLegend() + labs(x = "UMAP1", y = "UMAP2")
levels(srat) <- c("PCT", "PST1", "PST2", "Injured PT", "Proliferating PT", "DTL & ATL", "TAL", "DCT", "CNT1", "CNT2", "CD PC", "CD ICA", "CD ICB", "Fibro", "Podo", "Endo", "Immune")

# DotPlot with cell cluster annotation
features <- c("ENSMUST00000111026", "ENSMUST00000010250", "ENSMUST00000084797", "ENSMUST00000109223", "ENSMUST00000029574", "ENSMUST00000068031", "ENSMUST00000099362", "ENSMUST00000023265", "ENSMUST00000110494", "ENSMUST00000034218", "ENSMUST00000029876", "ENSMUST00000023752", "ENSMUST00000042221", "ENSMUST00000001253", "ENSMUST00000202681", "ENSMUST00000006825", "ENSMUST00000031653", "ENSMUST00000183301")
labels <- c("Slc5a12-202", "Slc22a6-201", "Slc22a13-201", "Havcr1-203", "Vcam1-201", "Top2a-201", "Slc4a11-201", "Psca-201", "Slc12a1-202", "Slc12a3-201", "Calb1-201", "Aqp2-201", "Slc26a7-201", "Slc26a4-201", "Pdgfra-208", "Nphs1-201", "Flt1-202", "Ptprc-207")
levels(srat) <- rev(levels(srat))  # Code for visualization purposes
DotPlot(srat, features = features, cols = c("lightyellow","royalblue")) + scale_x_discrete(labels = setNames(labels, features)) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "Isoform", y = "Cell type") + guides(color = guide_colorbar(title = "Average\nExpression"), size = guide_legend(title = "Percent\nExpressed"))
levels(srat) <- rev(levels(srat))  # Code for visualization purposes

# Rename object
srat_i <- srat
remove(srat)
gc()

# Load data
srat <- readRDS("healthy_and_post-IRI_integrated.rds")

# Import GTF file
tx2gene <- import_gtf(gtf_file = "gencode.vM33.primary_assembly.annotation.gtf.gz")

# Move gene and transcript identifier columns to the beginning
tx2gene <- move_columns_to_front(df = tx2gene, columns = c("transcript_id", "gene_name", "transcript_name", "gene_id"))

# Remove period and following numbers in transcript_id
tx2gene$transcript_id <- gsub("\\..*", "", tx2gene$transcript_id)

# Make sure all gene_id and transcript_id has a unique gene_name and transcript_name
tx2gene$gene_name <- one_to_one_mapping(name = tx2gene$gene_name, id = tx2gene$gene_id, ext = "_overlap")
tx2gene$transcript_name <- one_to_one_mapping(name = tx2gene$transcript_name, id = tx2gene$transcript_id, ext = "_overlap")

# Cell type
table1 <- merge(srat@meta.data["celltype"], srat_i@meta.data["celltype_isoform"], by = "row.names", all = TRUE, sort = FALSE)

# Remove cells with N/A
table1 <- dplyr::filter(table1, !is.na(celltype), !is.na(celltype_isoform))

# Create a summary
table_summary <- dplyr::count(table1, celltype, celltype_isoform)
table_summary$celltype <- factor(table_summary$celltype, levels = c("PCT", "PST1", "PST2", "Injured PT", "Proliferating PT", "DTL & ATL", "TAL", "DCT", "CNT", "CNT1", "CNT2", "CD PC", "CD ICA", "CD ICB", "Fibro", "JGA", "Podo", "Endo", "Uro", "Immune"))
table_summary$celltype_isoform <- factor(table_summary$celltype_isoform, levels = c("PCT", "PST1", "PST2", "Injured PT", "Proliferating PT", "DTL & ATL", "TAL", "DCT", "CNT", "CNT1", "CNT2", "CD PC", "CD ICA", "CD ICB", "Fibro", "JGA", "Podo", "Endo", "Uro", "Immune"))

# Create Sankey diagram
celltypes <- c("PCT", "PST1", "PST2", "Injured PT", "Proliferating PT", "DTL & ATL", "TAL", "DCT", "CNT", "CNT1", "CNT2", "CD PC", "CD ICA", "CD ICB", "Fibro", "JGA", "Podo", "Endo", "Uro", "Immune")
type_colors <- setNames(scales::hue_pal()(length(celltypes)), celltypes)

ggplot(table_summary, aes(axis1 = celltype, axis2 = celltype_isoform, y = n)) +
geom_alluvium(aes(fill = celltype), width = 0.25) +
geom_stratum(aes(fill = after_stat(stratum)), width = 0.25, color = "transparent") +
geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
scale_x_discrete(limits = c("Gene-based clustering", "Isoform-based clustering"), expand = c(.1, .1)) +
scale_fill_manual(values = type_colors) +
labs(y = "Number of Cells") +
theme_minimal() +
theme(legend.position = "none", axis.text.x = element_text(color = "black"), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
