# Load libraries
library(ggplot2)
library(dplyr)
library(Seurat)
library(presto)
library(future)
library(Matrix)
library(RColorBrewer)
library(tidyr)
library(hexbin)
library(ggthemes)

# Load data
srat <- readRDS("single_cell_and_nucleus_integrated.rds")

# Subsetting the object
sratc <- subset(srat, subset = method == "sc")
sratn <- subset(srat, subset = method == "sn")

# Get the counts matrix (raw counts)
countc <- GetAssayData(sratc, assay = "RNA", layer = "counts")
countn <- GetAssayData(sratn, assay = "RNA", layer = "counts")

# Count number of cells expressing each gene (non-zero counts)
num_cells_expressing_c <- rowSums(countc > 0)
num_cells_expressing_n <- rowSums(countn > 0)

# Convert to data frames
df_c <- data.frame(gene = names(num_cells_expressing_c), value_c = num_cells_expressing_c)
df_n <- data.frame(gene = names(num_cells_expressing_n), value_n = num_cells_expressing_n)

# Merge by gene name, keeping all genes from both (union)
df <- merge(df_c, df_n, by = "gene", all = TRUE)

# Replace NA with 0
df$value_c[is.na(df$value_c)] <- 0
df$value_n[is.na(df$value_n)] <- 0

# Count the number of cells and nuclei
noc <- sum(srat$method == "sc")
non <- sum(srat$method == "sn")

# Calculate the proportion of cells and nuclei expressing each gene
df$value_c <- df$value_c / noc
df$value_n <- df$value_n / non

# Set color palette
pal.spectral <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(100)

# Binned scatter plot showing the proportion of cells or nuclei that express each gene
ggplot(df, aes(x = value_c, y = value_n)) +
stat_bin_hex(bins = 50) +
scale_fill_gradientn(colours = pal.spectral, trans = "log10", name = "Number\nof genes") +
geom_abline(intercept = 0, slope = 1, color = "black") +
labs(x = "Proportion of cells expressing", y = "Proportion of nuclei expressing") +
theme_bw() + 
theme(panel.grid.minor = element_blank(), axis.title = element_text(size=12))

# Find markers
DefaultAssay(srat) <- "SCT"
Idents(srat) <- "method"
srat <- PrepSCTFindMarkers(srat, verbose = TRUE)
markers.sn.sc <- FindMarkers(srat, ident.1 = "sn", ident.2 = "sc")

# Replace 0 with 10^(-320) in p_val_adj
markers.sn.sc$p_val_adj[markers.sn.sc$p_val_adj == 0] <- 10^(-320)

# Binned scatter plot of differentially detected genes
ggplot(markers.sn.sc, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
stat_bin_hex(bins = 50) +
geom_hline(yintercept = -log10(0.05), color = "black", size = 0.5) +
geom_hline(yintercept = 320, linetype = "dashed", color = "black", size = 0.5) +
geom_vline(xintercept = c(-log2(1.5), log2(1.5)), color = "black", size = 0.5) +
xlim(c(-15, 15)) +
scale_fill_gradientn(colours = pal.spectral, trans = "log10", name = "Number\nof genes") +
labs(x = "Nuclei vs Cells (log2 fold change)", y = "Significance (-log10 adjusted P value)") +
theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Violin plots
genes <- c("mt-Nd1", "mt-Cytb", "mt-Co1")  # The following sets were also used for visualization: c("Rplp1", "Rpl32", "Rps3"), c("Hspa8", "Hsp90ab1", "Hspe1"), c("Lncpint", "Mir99ahg", "Gm15713"), c("Slc12a6", "Slc25a21", "Slc47a1"), c("Zfand3", "Zcchc7", "Zfr")
dfgene <- FetchData(srat, vars = c(genes, "method"), layer = "data")
dfgene <- pivot_longer(dfgene, cols = all_of(genes), names_to = "gene", values_to = "value")
dfgene$gene <- factor(dfgene$gene, levels = genes) 
ggplot(dfgene, aes(x = method, y = value, fill = method)) +
facet_wrap(~ gene, nrow = 3, dir = "v") +
geom_violin(scale = "width", adjust = 2, show.legend = FALSE) +
scale_fill_manual(values = brewer.pal(3, "Set1")[2:1]) +
labs(x = "", y = "Normalized UMI count") +
theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text = element_text(face = "italic"))

# Obtain count matrix
isoformCountsC <- GetAssayData(sratc, assay = "ISOFORM", layer = "counts")
isoformCountsN <- GetAssayData(sratn, assay = "ISOFORM", layer = "counts")

# Keep isoforms with non-zero expression levels
isoformCountsC <- isoformCountsC[Matrix::rowSums(isoformCountsC) > 0, ]
isoformCountsN <- isoformCountsN[Matrix::rowSums(isoformCountsN) > 0, ]
dfc <- data.frame(isoform_id = rownames(isoformCountsC))
dfn <- data.frame(isoform_id = rownames(isoformCountsN))

# Add meta data
dfc$isoform_name <- tx2gene$transcript_name[match(dfc$isoform_id, tx2gene$transcript_id)]
dfn$isoform_name <- tx2gene$transcript_name[match(dfn$isoform_id, tx2gene$transcript_id)]
dfc$gene_name <- tx2gene$gene_name[match(dfc$isoform_id, tx2gene$transcript_id)]
dfn$gene_name <- tx2gene$gene_name[match(dfn$isoform_id, tx2gene$transcript_id)]
dfc$gene_id <- tx2gene$gene_id[match(dfc$isoform_id, tx2gene$transcript_id)]
dfn$gene_id <- tx2gene$gene_id[match(dfn$isoform_id, tx2gene$transcript_id)]
dfc$method <- "sc"
dfn$method <- "sn"

# Number of detected isoforms for each gene
nIsoformTypesC <- table(dfc$gene_name)
nIsoformTypesN <- table(dfn$gene_name)
max(nIsoformTypesC)  # Returns 92
max(nIsoformTypesN)  # Returns 111
nIsoformTypesC <- sort(nIsoformTypesC, decreasing = TRUE)
nIsoformTypesN <- sort(nIsoformTypesN, decreasing = TRUE)
head(nIsoformTypesC, 1)  # Returns Gas5, 92
head(nIsoformTypesN, 1)  # Returns Gas5, 111
nIsoformTypes <- data.frame(
  count = c(as.numeric(nIsoformTypesC), as.numeric(nIsoformTypesN)),
  method = rep(c("sc", "sn"), c(length(nIsoformTypesC), length(nIsoformTypesN)))
)

# Bar chart
ggplot(nIsoformTypes, aes(x = count, fill = method)) +
geom_bar(width = 0.8, color = "white", position = "dodge") +
labs(x = "Number of isoforms per gene", y = "Number of genes") +
theme_few() +
theme(text = element_text(size = 12)) +
scale_x_continuous(expand = c(0,0), limits = c(0,15), breaks=seq(1,14,length=14)) +  # Note that genes outside this range are not plotted.
scale_y_continuous(expand = c(0,0), limits = c(0,10000)) +
NoLegend()

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
VlnPlot(srat, feature = "diversity", group.by = "celltype", pt.size = 0, split.by = "method") + NoLegend() +
scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
labs(title = NULL, x = NULL, y = "Diversity")

# UMAP plot
FeaturePlot(srat, reduction = "umap_harmony", features = "diversity", cols = c("grey", "red"), split.by = "method") & labs(title = NULL, x = "UMAP1", y = "UMAP2") & theme(legend.position = "right")
