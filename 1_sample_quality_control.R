# Load libraries
library(ggplot2)
library(dplyr)
library(Seurat)
library(sctransform)
library(SoupX)
library(Matrix)
library(plotly)
library(DoubletFinder)

# Load data - Sample code for sample scWT1 (single-cell WT, duplicate #1). Replace "sc" with "sn" (single-nucleus) and/or "WT1" with "WT2", "WT3", "IRI1", "IRI2" or "IRI3" for other samples.
tod = ReadMtx(mtx = "data/sc_WT1_gene_matrix.mtx.gz", features = "data/sc_WT1_gene_features.tsv.gz", cells = "data/sc_WT1_gene_barcodes.tsv.gz")
srat <- CreateSeuratObject(counts = tod)
tod_i = ReadMtx(mtx = "data/sc_WT1_isoform_matrix.mtx.gz", features = "data/sc_WT1_isoform_features.tsv.gz", cells = "data/sc_WT1_isoform_barcodes.tsv.gz")
isoform <- CreateAssayObject(counts = tod_i)
srat[["ISOFORM"]] <- isoform
srat@meta.data$sample <-"scWT1"
srat@meta.data$method <- "sc"

# Ambient RNA removal with SoupX (gene) - clustering
srat <- SCTransform(srat, verbose = F)
srat <- RunPCA(srat, verbose = F)
srat <- RunUMAP(srat, dims = 1:30, verbose = F)
srat <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat <- FindClusters(srat, verbose = T)

# Ambient RNA removal with SoupX (gene) - estimating contamination fraction
sn = SoupChannel(tod, tod, calcSoupProfile = FALSE)
soupProf = data.frame(row.names = rownames(tod), est = rowSums(tod)/sum(tod), counts = rowSums(tod))
sn = setSoupProfile(sn, soupProf)
meta <- srat@meta.data
umap <- srat@reductions$umap@cell.embeddings
sn <- setClusters(sn, setNames(meta$seurat_clusters, rownames(meta)))
sn <- setDR(sn, umap)
sn = autoEstCont(sn)

# Ambient RNA removal with SoupX (gene) - adjust counts and merge to the original Seurat object
out = adjustCounts(sn, roundToInt = T)
srat1 <- CreateSeuratObject(counts = out)
srat[["RNA"]] <- srat1[["RNA"]]  # This code was run because srat1 does not contain isoform data.

# Ambient RNA removal with SoupX (isoform) - clustering
DefaultAssay(srat) <- "ISOFORM"
srat <- SCTransform(srat, new.assay.name = "SCTiso", verbose = F)
srat <- RunPCA(srat, reduction.name = "pcaiso", verbose = F)
srat <- RunUMAP(srat, reduction = "pcaiso", reduction.name = "umapiso", dims = 1:30, verbose = F)
srat <- FindNeighbors(srat, reduction = "pcaiso", dims = 1:30, verbose = F)
srat <- FindClusters(srat, graph.name = "SCTiso_snn", verbose = T)

# Ambient RNA removal with SoupX (isoform) - estimating contamination fraction
sn_i = SoupChannel(tod_i, tod_i, calcSoupProfile = FALSE)
soupProf_i = data.frame(row.names = rownames(tod_i), est = rowSums(tod_i)/sum(tod_i), counts = rowSums(tod_i))
sn_i = setSoupProfile(sn_i, soupProf_i)
meta_i <- srat@meta.data  # seurat_clusters is now SCTiso_snn_res.0.8 and not SCT_snn_res.0.8
umap_i <- srat@reductions$umapiso@cell.embeddings
sn_i <- setClusters(sn_i, setNames(meta_i$seurat_clusters, rownames(meta_i)))
sn_i <- setDR(sn_i, umap_i, reductName = c("umapiso_1", "umapiso_2"))
sn_i = autoEstCont(sn_i)

# Ambient RNA removal with SoupX (isoform) - adjust counts and merge to the original Seurat object
out_i = adjustCounts(sn_i, roundToInt = T)
srat2 <- CreateSeuratObject(counts = out_i, assay = "ISOFORM")
srat[["ISOFORM"]] <- srat2[["ISOFORM"]]

# Tidy up the data
DefaultAssay(srat) <- "RNA"
srat@assays[c("SCT", "SCTiso")] <- NULL
srat@reductions[c("pca", "umap", "pcaiso", "umapiso")] <- NULL
srat@graphs[c("SCT_nn", "SCT_snn", "SCTiso_nn", "SCTiso_snn")] <- NULL
Idents(srat) <- "NULL"
srat@meta.data <- srat@meta.data[, !(colnames(srat@meta.data) %in% c("orig.ident", "nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.8", "seurat_clusters", "nCount_SCTiso", "nFeature_SCTiso", "SCTiso_snn_res.0.8"))]

# Knee plot
counts = data.frame(counts = srat@meta.data$nCount_RNA)
counts = counts[order(counts$counts, decreasing = TRUE), , drop = FALSE]
counts$ID <- seq.int(nrow(counts))
ggplot(counts, aes(x = ID, y = counts)) +
geom_point(color = "blue3", size = 0.8) +
scale_x_continuous(trans = "log10", limits = c(1, NA), minor_breaks = NULL) +
scale_y_continuous(trans = "log10", limits = c(1, NA), minor_breaks = NULL) +
theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "#B5B5B5", linewidth = 0.2), axis.line.x.top = element_blank(), axis.line.x.bottom = element_line(color = "black", linewidth = 0.3), axis.line.y.left = element_line(color = "black", linewidth = 0.2), axis.line.y.right = element_blank(), axis.ticks = element_line(color = "black", linewidth = 0.2)) +
labs(x = "Cell rank", y = "UMI counts")

# Cell quality control - precheck
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^mt-")
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt")

# Cell quality control - nCount_RNA threshold was set based on the knee plot.
srat <- subset(srat, subset = nFeature_RNA > 200 & nCount_RNA > 281)  # nCount_RNA threshold was 281 for scWT1, 193 for scWT2, 231 for scWT3, 514 for snWT1, 538 for snWT2, 367 for snWT3, 631 for snIRI1, 678 for snIRI2 and 699 for snIRI3. Also, "percent.mt < 1" was added for single-nucleus samples.
nFeatureLimit = quantile(srat$nFeature_RNA, 0.99)
nCountLimit = quantile(srat$nCount_RNA, 0.99)
srat <- subset(srat, subset = nFeature_RNA < nFeatureLimit & nCount_RNA < nCountLimit)
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Normalization and clustering
srat <- SCTransform(srat, verbose = F)
srat <- RunPCA(srat, verbose = F)
srat <- RunUMAP(srat, dims = 1:30, verbose = F)
srat <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat <- FindClusters(srat, resolution = 0.8, verbose = T)

# Quality check by clusters - UMAP plot
DimPlot(srat, reduction = "umap", label = T)
FeaturePlot(srat, reduction = "umap", features = "nFeature_RNA", cols = c("grey", "red"))
FeaturePlot(srat, reduction = "umap", features = "nCount_RNA", cols = c("grey", "red"))
FeaturePlot(srat, reduction = "umap", features = "percent.mt", cols = c("grey", "red"))

# Quality check by clusters - violin plot
VlnPlot(srat, feature = "nFeature_RNA", pt.size = 0.3)
VlnPlot(srat, feature = "nCount_RNA", pt.size = 0.3)
VlnPlot(srat, feature = "percent.mt", pt.size = 0.3)

# Cluster profile
table(srat@meta.data$SCT_snn_res.0.8)
features <- c("Slc5a2", "Slc5a12", "Slc22a6", "Slc22a13", "Havcr1", "Vcam1", "Cryab", "Slc4a11", "Psca", "Slc12a1", "Umod", "Slc12a3", "Slc8a1", "Aqp2", "Slc26a7", "Slc26a4", "Ncam1", "Nphs1", "Nphs2", "Flt1", "Emcn", "Col1a1", "Col1a2", "Pdgfra", "Upk1b", "Ptprc")
levels(srat) <- rev(levels(srat))  # Code for visualization purposes
DotPlot(srat, features = features, cols = c("lightyellow","royalblue")) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
levels(srat) <- rev(levels(srat))  # Code for visualization purposes

# Identify low quality cell clusters from the outputs above
# Single-cell samples: The following clusters were identified as low quality cell clusters because of no cell-type marker gene expression, high percent.mt, low nFeature_RNA and low nCount_RNA. scWT1, clusters 0 and 2 (1,694 cells); scWT2, clusters 0, 2 and 13 (1,763 cells); scWT3, clusters 0, 2 and 14 (2,601 cells).
# Single-nucleus samples: The following clusters were identified as low quality cell clusters because of expression of multiple marker genes, no marker genes that are specific to these clusters, high percent.mt, low nFeature_RNA and low nCount_RNA. snWT1, clusters 1 and 12 (1,372 cells); snWT2, clusters 7 and 9 (1,237 cells); snWT3, clusters 7 and 12 (933 cells).
# The numbers of recovered cells after removing low quality cell clusters are as follows. scWT1, 8,421 - 1,694 = 6,727 cells; scWT2, 7,489 - 1,763 = 5,726 cells; scWT3, 10,347 - 2,601 = 7,746 cells; snWT1, 11,130 - 1,372 = 9,758 cells; snWT2, 12,682 - 1,237 = 11,445 cells; snWT3, 12,853 - 933 = 11,920 cells; snIRI1, 10,009 - 2,003 = 8,006 cells; snIRI2, 9,119 - 788 = 8,331 cells; snIRI3, 9,823 - 2,834 = 6,989 cells.
# According to the multiple rate table from 10x, the estimated multiplet rates are as follows. scWT1, 5.4%; scWT2, 4.6%; scWT3, 6.2%; snWT1, 8.0%; snWT2, 8.8%; snWT3, 9.6%; snIRI1, 6.4%; snIRI2, 6.4%; snIRI3, 5.6%.

# Removal of low quality cells
srat <- subset(srat, idents = c(0, 2), invert = T)  # Replace c(0, 2) with c(0, 2, 13) for scWT2, c(0, 2, 14) for scWT3, c(1, 12) for snWT1, c(7, 9) for snWT2, c(7, 12) for snWT3, c(0, 9, 10) for snIRI1, c(8, 11, 22) for snIRI2 and c(1, 3, 7) for snIRI3.

# Doublet estimation with DoubletFinder - step 1
sweep.res.list_kidney <- paramSweep(srat, PCs = 1:20, sct = TRUE)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)
ggplot(bcmvn_kidney, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
pK <- bcmvn_kidney %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))
annotations <- srat@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*nrow(srat@meta.data))  # Estimated multiplet rate (see above): replace 0.054 with 0.046 for scWT2, 0.062 for scWT3, 0.08 for snWT1, 0.088 for snWT2, 0.096 for snWT3, 0.064 for snIRI1, 0.064 for snIRI2 and 0.056 for snIRI3.
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
srat_df <- doubletFinder(srat, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi, sct = TRUE)

# Doublet estimation with DoubletFinder - step 2
DimPlot(srat_df, reduction = 'umap', group.by = "DF.classifications_0.25_0.05_363")  # For other samples, change group.by value according to the previous output.
srat_df2 <- doubletFinder(srat_df, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.05_363", sct = TRUE)  # For other samples, change reuse.pANN value according to the previous output.
head(srat_df2@meta.data, 3)

# Doublet estimation with DoubletFinder - step 3
DimPlot(srat_df2, reduction = 'umap', group.by = "DF.classifications_0.25_0.05_339")  # For other samples, change group.by value according to the previous output.
table(srat_df2@meta.data$DF.classifications_0.25_0.05_339)  # For other samples, change DF.classifications_0.25_0.05_339 according to the previous output.

# Tidy up the data
srat_df2$DoubletFinder <- srat_df2$DF.classifications_0.25_0.05_339  # For other samples, change DF.classifications_0.25_0.05_339 according to the previous output.
srat_df2@meta.data <- srat_df2@meta.data[, !(colnames(srat_df2@meta.data) %in% c("pANN_0.25_0.05_363", "DF.classifications_0.25_0.05_363", "DF.classifications_0.25_0.05_339"))]  # For other samples, change the column names according to the previous output.

# Save data - Sample code for scWT1. Replace scWT1 with scWT2, scWT3, snWT1, snWT2, snWT3, snIRI1, snIRI2 or snIRI3 for other samples.
saveRDS(srat_df2, file = "scWT1.rds")
