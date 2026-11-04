# Load libraries
library(ggplot2)
library(dplyr)
library(Seurat)
library(presto)
library(future)
library(harmony)
library(Matrix)
library(clusterProfiler)
library(fgsea)
library(AnnotationHub)
library(AnnotationDbi)

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
DefaultAssay(srat) <- "RNA"
srat@assays["SCT"] <- NULL
Idents(srat) <- "NULL"
srat@meta.data <- srat@meta.data[, !(colnames(srat@meta.data) %in% c("nCount_SCT", "nFeature_SCT", "seurat_clusters"))]
remove(srat1, srat2, srat3, srat4, srat5, srat6)
gc()

# Clustering before integration
srat <- SCTransform(srat)
srat <- RunPCA(srat)

# Integrate samples - Harmony integration
srat_harmony <- IntegrateLayers(srat, method = HarmonyIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "integrated_harmony", verbose = T)
srat_harmony <- RunUMAP(srat_harmony, dims = 1:30, reduction = "integrated_harmony", reduction.name = "umap_harmony")
srat_harmony <- FindNeighbors(srat_harmony, reduction = "integrated_harmony", dims = 1:30)
srat_harmony <- FindClusters(srat_harmony, resolution = 0.4)

# Check culstering after integration
DimPlot(srat_harmony, reduction = "umap_harmony", group.by = "sample", pt.size = .1)
DimPlot(srat_harmony, reduction = "umap_harmony", group.by = "intervention", pt.size = .1)

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
levels(srat_harmony) <- rev(levels(srat_harmony))
DotPlot(srat_harmony, features = features, cols = c("lightyellow","royalblue")) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
levels(srat_harmony) <- rev(levels(srat_harmony))

# Checking DoubletFinder annotations
table(srat_harmony@meta.data$DoubletFinder)
table(srat_harmony@meta.data$sample, srat_harmony@meta.data$DoubletFinder)
table(srat_harmony@meta.data$SCT_snn_res.0.4, srat_harmony@meta.data$DoubletFinder)
DimPlot(srat_harmony, reduction = "umap_harmony", group.by = "DoubletFinder", shuffle = T)

# Clusters 13, 15, 18, 20 and 21 are likely doublets because two or more marker genes are present. More than 45% of the cells in each cluster are annotated as doublets by DoubletFinder. Therefore, all cells in these clusters are assumed to be doublets and will be removed from the analysis. All cells annotated as doublets by DoubletFinder will be removed from the analysis as well.

# Removing doubletes
srat_harmony_c <- subset(srat_harmony, srat_harmony$SCT_snn_res.0.4 %in% c("13", "15", "18", "20", "21"), invert = T)
srat_harmony_c <- subset(srat_harmony_c, subset = srat_harmony_c$DoubletFinder == "Singlet")

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
srat_harmony_c <- RunPCA(srat_harmony_c)

# Reintegrate samples - Harmony integration
srat_harmony_re <- IntegrateLayers(srat_harmony_c, method = HarmonyIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "integrated_harmony", verbose = T)
srat_harmony_re <- RunUMAP(srat_harmony_re, dims = 1:30, reduction = "integrated_harmony", reduction.name = "umap_harmony")
srat_harmony_re <- FindNeighbors(srat_harmony_re, reduction = "integrated_harmony", dims = 1:30)
srat_harmony_re <- FindClusters(srat_harmony_re, resolution = 0.3)

# Check culstering after reintegration
DimPlot(srat_harmony_re, reduction = "umap_harmony", group.by = "sample", pt.size = .1)
DimPlot(srat_harmony_re, reduction = "umap_harmony", group.by = "intervention", pt.size = .1)
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
levels(srat_harmony_re) <- rev(levels(srat_harmony_re))
DotPlot(srat_harmony_re, features = features, cols = c("lightyellow","royalblue")) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
levels(srat_harmony_re) <- rev(levels(srat_harmony_re))

# Join layers
srat_harmony_re <- JoinLayers(srat_harmony_re, assay = "RNA") %>% JoinLayers(assay = "ISOFORM")

# Find markers
Idents(srat_harmony_re) = srat_harmony_re$SCT_snn_res.0.3
srat_harmony_re <- PrepSCTFindMarkers(srat_harmony_re, assay = "SCT", verbose = TRUE)
cluster_markers <- list()
for (i in 0:17) {
  cluster_markers[[paste0("cluster", i, ".markers")]] <- FindMarkers(srat_harmony_re, ident.1 = i)
  write.csv(cluster_markers[[paste0("cluster", i, ".markers")]], paste0("cluster", i, "markers.csv"))
}

# Annotate clusters
new.cluster.ids <- c("PCT", "TAL", "PST1", "Endo", "DCT", "Injured PT", "PST2", "Fibro", "CD PC", "CNT", "Proliferating PT", "DTL & ATL", "Immune", "CD ICA", "CD ICB", "Podo", "Uro", "JGA")
names(new.cluster.ids) <- levels(srat_harmony_re)
srat_harmony_re <- RenameIdents(srat_harmony_re, new.cluster.ids)
levels(srat_harmony_re) <- c("PCT", "PST1", "PST2", "Injured PT", "Proliferating PT", "DTL & ATL", "TAL", "DCT", "CNT", "CD PC", "CD ICA", "CD ICB", "Fibro", "JGA", "Podo", "Endo", "Uro", "Immune")
srat_harmony_re@meta.data$celltype <- Idents(srat_harmony_re)

# UMAP with cell cluster annotation
levels(srat) <- c("PCT", "TAL", "PST1", "Endo", "DCT", "Injured PT", "PST2", "Fibro", "CD PC", "CNT", "Proliferating PT", "DTL & ATL", "Immune", "CD ICA", "CD ICB", "Podo", "Uro", "JGA")  # Code for visualization purposes
DimPlot(srat, reduction = "umap_harmony", label = TRUE) + NoLegend() + labs(x = "UMAP1", y = "UMAP2")
DimPlot(srat, reduction = "umap_harmony", label = FALSE, repel = TRUE, split.by = "intervention") + NoLegend() + labs(x = "UMAP1", y = "UMAP2") + theme(strip.text = element_blank(), strip.background = element_blank())
levels(srat) <- c("PCT", "PST1", "PST2", "Injured PT", "Proliferating PT", "DTL & ATL", "TAL", "DCT", "CNT", "CD PC", "CD ICA", "CD ICB", "Fibro", "JGA", "Podo", "Endo", "Uro", "Immune")  # Code for visualization purposes

# DotPlot with cell cluster annotation
features <- c("Slc5a12", "Slc22a6", "Slc22a13", "Havcr1", "Vcam1", "Top2a", "Cryab", "Slc12a1", "Slc12a3", "Calb1", "Aqp2", "Slc26a7", "Slc26a4", "Pdgfra", "Dgkb", "Nphs1", "Flt1", "Upk1b", "Ptprc")
levels(srat_harmony_re) <- rev(levels(srat_harmony_re))  # Code for visualization purposes
DotPlot(srat_harmony_re, features = features, cols = c("lightyellow","royalblue")) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
levels(srat_harmony_re) <- rev(levels(srat_harmony_re))  # Code for visualization purposes

# Save data
saveRDS(srat_harmony_re, file = "healthy_and_post-IRI_integrated.rds")

# Rename object
srat <- srat_harmony_re

# Save data post-SCTransform on isoform expression
DefaultAssay(srat_harmony_re) <- "ISOFORM"
srat_harmony_re[["ISOFORM"]] <- split(srat_harmony_re[["ISOFORM"]], f = srat_harmony_re$sample)
srat_harmony_re <- SCTransform(srat_harmony_re, assay = "ISOFORM", new.assay.name = "SCT_i")
saveRDS(srat_harmony_re, "healthy_and_post-IRI_integrated_SCTi.rds")

# Obtain count matrix for calculating Shannon entropy
isoformCounts <- GetAssayData(srat, assay = "ISOFORM", layer = "counts")
isoformCounts <- t(isoformCounts)

# Total counts per cell (row sums of cells x isoforms)
cell_totals <- rowSums(isoformCounts)

# Divide each row by its total to get proportions — stays sparse
p_matrix <- isoformCounts / cell_totals

# Drop zero to avoid NA in the following step
p_matrix <- drop0(p_matrix)

# Compute p * log(p) only for non-zero entries (zero entries remain 0)
p_matrix@x <- p_matrix@x * log(p_matrix@x)

# Shannon entropy: -sum(p * log(p)) per cell
shannon_entropy <- -rowSums(p_matrix)

# Store in Seurat object
srat <- AddMetaData(srat, metadata = shannon_entropy, col.name = "diversity")

# Violin plot
VlnPlot(srat, feature = "diversity", group.by = "celltype", pt.size = 0, split.by = "intervention") + NoLegend() +
scale_y_continuous(limits = c(0, 8), expand = expansion(mult = c(0, 0.05))) +
labs(title = NULL, x = NULL, y = "Diversity")

### Differentially expressed gene analysis and gene set enrichment analysis (proximal tubule cells from healthy kidneys vs post-IRI kidneys)
# Load data and rename groups
srat <- readRDS("healthy_and_post-IRI_integrated.rds")
srat$celltype_group <- paste(srat$celltype, srat$intervention, sep = "_")
Idents(srat) <- srat$celltype_group
srat <- RenameIdents(srat, "PCT_none" = "PT_none") %>% RenameIdents("PST1_none" = "PT_none") %>% RenameIdents("PST2_none" = "PT_none") %>% RenameIdents("Injured PT_none" = "PT_none") %>% RenameIdents("Proliferating PT_none" = "PT_none") %>% RenameIdents("PCT_IRI" = "PT_IRI") %>% RenameIdents("PST1_IRI" = "PT_IRI") %>% RenameIdents("PST2_IRI" = "PT_IRI") %>% RenameIdents("Injured PT_IRI" = "PT_IRI") %>% RenameIdents("Proliferating PT_IRI" = "PT_IRI")
srat@meta.data$celltype_group <- Idents(srat)

# Find markers (differentially expressed gene analysis)
srat <- PrepSCTFindMarkers(srat, assay = "SCT", verbose = TRUE)
markers.IriPTvsNonePT.gene <- FindMarkers(srat, ident.1 = "PT_IRI", ident.2 = "PT_none")

# Prepare dataframe and reference
dat <- data.frame(gene = rownames(markers.IriPTvsNonePT.gene), value = markers.IriPTvsNonePT.gene$avg_log2FC)
ah <- AnnotationHub()
query(ah, c("Mus musculus", "OrgDb"))
orgdb <- ah[["AH119530"]]
dat <- mutate(dat, entrez_id = mapIds(orgdb, keys = gene, column  = "ENTREZID", keytype = "SYMBOL", multiVals = "first"))  # keep first match if multiple
dat <- filter(dat, !is.na(entrez_id))
geneList <- dat$value
geneList <- as.numeric(geneList)
names(geneList) <- dat$entrez_id
geneList <- geneList[!is.na(geneList)]
geneList <- geneList[!is.na(names(geneList))]  # remove NA names
geneList <- geneList[names(geneList) != ""]    # remove empty names
geneList <- sort(geneList, decreasing = TRUE)

# Gene set enrichment analysis
gsea <- gseKEGG(geneList = geneList, organism = "mmu", minGSSize = 15, maxGSSize = 500, eps = 0, nPermSimple = 10000, seed = TRUE, by = "fgsea", pvalueCutoff = 1, pAdjustMethod = "BH")

# Visualizastion
dat2 <- gsea@result %>%
arrange(p.adjust) %>%
slice_head(n = 19) %>%
mutate(gene_count  = sign(NES) * setSize)

ggplot(dat2, aes(x = NES, y = reorder(Description, NES), fill = p.adjust)) +
geom_bar(stat = "identity") +
scale_fill_gradient2(mid = "red", high = "#FDDBC7", name = "FDR", limits = c(0, 0.05), breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05)) +
geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
labs(x = "Normalized Enrichment Score", y = "") +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.y = element_text(size = 9, color = "black"), axis.text.x = element_text(size = 9), legend.position = "right")
