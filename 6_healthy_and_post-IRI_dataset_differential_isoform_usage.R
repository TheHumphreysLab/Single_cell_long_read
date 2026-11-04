# Load libraries
library(ggplot2)
library(dplyr)
library(Seurat)
library(DTUrtle)
library(clusterProfiler)
library(fgsea)
library(AnnotationHub)
library(AnnotationDbi)

# Parallelize the computations
biocpar <- BiocParallel::MulticoreParam(16)

# Import GTF file
tx2gene <- import_gtf(gtf_file = "gencode.vM33.primary_assembly.annotation.gtf.gz")

# Move gene and transcript identifier columns to the beginning
tx2gene <- move_columns_to_front(df = tx2gene, columns = c("transcript_id", "gene_name", "transcript_name", "gene_id"))

# Remove period and following numbers in transcript_id
tx2gene$transcript_id <- gsub("\\..*", "", tx2gene$transcript_id)

# Make sure all gene_id and transcript_id has a unique gene_name and transcript_name
tx2gene$gene_name <- one_to_one_mapping(name = tx2gene$gene_name, id = tx2gene$gene_id, ext = "_overlap")
tx2gene$transcript_name <- one_to_one_mapping(name = tx2gene$transcript_name, id = tx2gene$transcript_id, ext = "_overlap")

# Load data
srat <- readRDS("healthy_and_post-IRI_integrated.rds")

# Obtain count matrix
isoformCounts <- GetAssayData(srat, assay = "ISOFORM", layer = "counts")

# Obtain metadata
metadata <- srat@meta.data
metadata <- metadata[, c("sample", "intervention", "celltype")]

# Transcript-level count matrix
cts <- combine_to_matrix(tx_list = isoformCounts)

# Estimate influence of priming-bias
unfilt_gtf <- import_gtf(gtf_file = "gencode.vM33.primary_assembly.annotation.gtf.gz", feature_type = NULL, out_df = FALSE)
unfilt_gtf$transcript_id <- gsub("\\..*", "", unfilt_gtf$transcript_id)
tx2gene2 <- priming_bias_detection_probability(counts = cts, tx2gene = tx2gene, gtf = unfilt_gtf, one_to_one = "_overlap", priming_enrichment = "3", add_to_table = tx2gene, BPPARAM = biocpar)

# Subsetting all proximal tubule cells
metadata_allPT <- subset(metadata, metadata$celltype == "PCT" | metadata$celltype == "PST1" | metadata$celltype == "PST" | metadata$celltype == "Proliferating PT" | metadata$celltype == "Injured PT")
cts_allPT <- cts[, colnames(cts) %in% rownames(metadata_allPT)]

# Differential isoform usage analysis
dtu_allPT <- run_drimseq(counts = cts_allPT, tx2gene = tx2gene2, pd = metadata_allPT, cond_col = "intervention", cond_levels = c("none", "IRI"), filtering_strategy = "sc", BPPARAM = biocpar)
dtu_allPT <- posthoc_and_stager(dturtle = dtu_allPT, ofdr = 0.05, posthoc = 0.1)
dtu_allPT <- create_dtu_table(dturtle = dtu_allPT, add_gene_metadata = list("chromosome" = "seqnames"), add_tx_metadata = list("max_tx_expr_in_none" = c("exp_in_none", max), "max_tx_expr_in_IRI" = c("exp_in_IRI", max)))
dtu_allPT$dtu_table$abs_diff_expr_in <- abs(dtu_allPT$dtu_table$max_tx_expr_in_none - dtu_allPT$dtu_table$max_tx_expr_in_IRI)
dtu_allPT <- plot_proportion_barplot(dturtle = dtu_allPT, meta_gene_id = "gene", savepath = "images_allPT", add_to_table = "barplot", BPPARAM = biocpar)
dtu_allPT <- plot_proportion_pheatmap(dturtle = dtu_allPT, include_expression = TRUE, treeheight_col = 20, sample_meta_table_columns = c("sample_id", "sample", "condition"), savepath = "images_allPT", add_to_table = "pheatmap", BPPARAM = biocpar)
dtu_allPT <- plot_transcripts_view(dturtle = dtu_allPT, gtf = unfilt_gtf, genome = "mm39", one_to_one = "_overlap", savepath = "images_allPT", add_to_table = "transcript_view", BPPARAM = biocpar)
column_formatter_list <- list(
      "gene_qvalue" = table_pval_tile("white", "orange", digits = 3),
      "minimal_tx_qvalue" = table_pval_tile("white", "orange", digits = 3),
      "number_tx" = formattable::color_tile('white', "lightblue"),
      "number_significant_tx" = formattable::color_tile('white', "lightblue"),
      "max(IRI)" = table_percentage_bar('lightgreen', "#FF9999", digits=2),
      "max_tx_expr_in_none" = table_percentage_bar('white', "lightblue", color_break = 0, digits=2),
      "max_tx_expr_in_IRI" = table_percentage_bar('white', "lightblue", color_break = 0, digits=2),
      "abs_diff_expr_in" = table_percentage_bar('white', "lightblue", color_break = 0, digits=2))
plot_dtu_table(dturtle = dtu_allPT, savepath = "dtu_allPT.html", column_formatters = column_formatter_list)

# Genes with two isoforms detected
dtu_allPT_meta_table_tx_2 <- dtu_allPT$meta_table_tx[dtu_allPT$meta_table_tx$gene %in% rownames(subset(dtu_allPT$dtu_table, number_tx == 2)), ]
write.csv(dtu_allPT_meta_table_tx_2, "dtu_allPT_meta_table_tx_2.csv")

# Classify genes in dtu_allPT_meta_table_tx_2.csv to four categories depending on their shifts in isoform expression: protein-coding to protein-coding (PP), protein-coding to non-protein-coding (PN), non-protein-coding to protein-coding (NP), non-protein-coding to non-protein-coding (NN). Create a CSV file with the gene names in a column named gene and their corresponding categories in a column named type. Save the file as ORA_data.csv.

# Load data for over-representation analysis
dat <- read.csv("data_ORA.csv")
ah <- AnnotationHub()
query(ah, c("Mus musculus", "OrgDb"))
orgdb <- ah[["AH119530"]]
dat <- mutate(dat, entrez_id = mapIds(orgdb, keys = gene, column  = "ENTREZID", keytype = "SYMBOL", multiVals = "first"))
dat <- filter(dat, !is.na(entrez_id))

# Over-representation analysis: protein-coding to non-protein-coding
genes <- as.vector(subset(dat, type == "PN")$entrez_id)
ORApn <- enrichKEGG(gene = genes, organism = "mmu", keyType = "ncbi-geneid", pvalueCutoff = 0.05)
dotplot(ORApn, x = "FoldEnrichment", showCategory = 10) + xlim(1, 7) + scale_fill_gradient2(mid = "red", high = "#FDDBC7", name = "FDR")

# Over-representation analysis: protein-coding to protein-coding
genes <- as.vector(subset(dat, type == "PP")$entrez_id)
ORApn <- enrichKEGG(gene = genes, organism = "mmu", keyType = "ncbi-geneid", pvalueCutoff = 0.05)
dotplot(ORApn, x = "FoldEnrichment", showCategory = 10) + xlim(1, 7) + scale_fill_gradient2(mid = "red", high = "#FDDBC7", name = "FDR")

### Differentially expressed isoform analysis (proximal tubule cells from healthy kidneys vs post-IRI kidneys)
# Load data and rename groups
srat <- readRDS("healthy_and_post-IRI_integrated_SCTi.rds")
srat$celltype_group <- paste(srat$celltype, srat$intervention, sep = "_")
Idents(srat) <- srat$celltype_group
srat <- RenameIdents(srat, "PCT_none" = "PT_none") %>% RenameIdents("PST1_none" = "PT_none") %>% RenameIdents("PST2_none" = "PT_none") %>% RenameIdents("Injured PT_none" = "PT_none") %>% RenameIdents("Proliferating PT_none" = "PT_none") %>% RenameIdents("PCT_IRI" = "PT_IRI") %>% RenameIdents("PST1_IRI" = "PT_IRI") %>% RenameIdents("PST2_IRI" = "PT_IRI") %>% RenameIdents("Injured PT_IRI" = "PT_IRI") %>% RenameIdents("Proliferating PT_IRI" = "PT_IRI")
srat@meta.data$celltype_group <- Idents(srat)

# Find markers
for (i in 1:6) {
    slot(object = srat@assays$SCT_i@SCTModel.list[[i]], name = "umi.assay") <- "ISOFORM"
}
srat <- PrepSCTFindMarkers(srat, assay = "SCT_i", verbose = TRUE)
markers.IriPTvsNonePT <- FindMarkers(srat, assay = "SCT_i", ident.1 = "PT_IRI", ident.2 = "PT_none")

# Add metadata
markers.IriPTvsNonePT$transcript_name <- tx2gene$transcript_name[match(markers.IriPTvsNonePT$X, tx2gene$transcript_id)]
markers.IriPTvsNonePT$transcript_type <- tx2gene$transcript_type[match(markers.IriPTvsNonePT$X, tx2gene$transcript_id)]
markers.IriPTvsNonePT$gene_id <- tx2gene$gene_id[match(markers.IriPTvsNonePT$X, tx2gene$transcript_id)]
markers.IriPTvsNonePT$gene_name <- tx2gene$gene_name[match(markers.IriPTvsNonePT$X, tx2gene$transcript_id)]
markers.IriPTvsNonePT$protein_id <- tx2gene$protein_id[match(markers.IriPTvsNonePT$X, tx2gene$transcript_id)]
write.csv(markers.IriPTvsNonePT, "markers_IriPTvsNonePT_detail.csv")

# Among oxidative phosphorylation-associated genes defined in KEGG (mmu00190) and detected in our datasets, the isoform expression of genes with isoform-level regulation, defined as concordant up- or down-regulation of all protein-coding isoforms accompanied by opposite regulation of at least one non-protein-coding isoform, was plotted in a heatmap. The average log2 fold-change (column avg_log2FC in markers_IriPTvsNonePT_detail.csv) was plotted.

### Differentially expressed isoform analysis (injured proximal tubule cells vs healthy proximal tubule cells)
# Reload data and rename groups
srat <- readRDS("healthy_and_post-IRI_integrated_SCTi.rds")
srat <- RenameIdents(srat, "PCT" = "Healthy PT") %>% RenameIdents("PST1" = "Healthy PT") %>% RenameIdents("PST2" = "Healthy PT")
srat@meta.data$celltype_injPT <- Idents(srat)

# Find markers
for (i in 1:6) {
    slot(object = srat@assays$SCT_i@SCTModel.list[[i]], name = "umi.assay") <- "ISOFORM"
}
srat <- PrepSCTFindMarkers(srat, assay = "SCT_i", verbose = TRUE)
markers.InjuredPTvsHealthyPT <- FindMarkers(srat, assay = "SCT_i", ident.1 = "Injured PT", ident.2 = "Healthy PT")

# Add metadata
markers.InjuredPTvsHealthyPT$transcript_name <- tx2gene$transcript_name[match(markers.InjuredPTvsHealthyPT$X, tx2gene$transcript_id)]
markers.InjuredPTvsHealthyPT$transcript_type <- tx2gene$transcript_type[match(markers.InjuredPTvsHealthyPT$X, tx2gene$transcript_id)]
markers.InjuredPTvsHealthyPT$gene_id <- tx2gene$gene_id[match(markers.InjuredPTvsHealthyPT$X, tx2gene$transcript_id)]
markers.InjuredPTvsHealthyPT$gene_name <- tx2gene$gene_name[match(markers.InjuredPTvsHealthyPT$X, tx2gene$transcript_id)]
markers.InjuredPTvsHealthyPT$protein_id <- tx2gene$protein_id[match(markers.InjuredPTvsHealthyPT$X, tx2gene$transcript_id)]

# markers.InjuredPTvsHealthyPT was used for creating a volcano plot.

### Plotting Fut9 and Prrc2b expression
# Reload data
srat <- readRDS("healthy_and_post-IRI_integrated.rds")
srat_i <- readRDS("healthy_and_post-IRI_integrated_SCTi.rds")

# Check maximum expression levels of Fut9 isoforms
max(GetAssayData(srat_i, assay = "SCT_i", layer = "data")["ENSMUST00000084770", ])  # Returns 3.36729582998647
max(GetAssayData(srat_i, assay = "SCT_i", layer = "data")["ENSMUST00000108199", ])  # Returns 1.09861228866811
max(GetAssayData(srat_i, assay = "SCT_i", layer = "data")["ENSMUST00000149875", ])  # Returns 3.68887945411394

# Fut9 isoforms, UMAP
FeaturePlot(srat_i, reduction = "umap_harmony", features = "ENSMUST00000084770", cols = c("grey", "red"), split.by = "intervention", order = FALSE) & theme(plot.title = element_blank()) & scale_color_gradientn(colors = c("grey", "red"), limits = c(0, 3.7))
FeaturePlot(srat_i, reduction = "umap_harmony", features = "ENSMUST00000108199", cols = c("grey", "red"), split.by = "intervention", order = FALSE) & theme(plot.title = element_blank()) & scale_color_gradientn(colors = c("grey", "red"), limits = c(0, 3.7))
FeaturePlot(srat_i, reduction = "umap_harmony", features = "ENSMUST00000149875", cols = c("grey", "red"), split.by = "intervention", order = FALSE) & theme(plot.title = element_blank()) & scale_color_gradientn(colors = c("grey", "red"), limits = c(0, 3.7))

# Fut9 isoforms, violin plot
VlnPlot(srat_i, features = "ENSMUST00000084770", sort = FALSE, flip = TRUE, add.noise = FALSE, pt.size = 0, split.by = "intervention") + theme(legend.position = "none", axis.title.x = element_blank(), plot.title = element_blank())
VlnPlot(srat_i, features = "ENSMUST00000108199", sort = FALSE, flip = TRUE, add.noise = FALSE, pt.size = 0, split.by = "intervention") + theme(legend.position = "none", axis.title.x = element_blank(), plot.title = element_blank())
VlnPlot(srat_i, features = "ENSMUST00000149875", sort = FALSE, flip = TRUE, add.noise = FALSE, pt.size = 0, split.by = "intervention") + theme(legend.position = "none", axis.title.x = element_blank(), plot.title = element_blank())

# Fut9 gene, UMAP
FeaturePlot(srat, reduction = "umap_harmony", features = "Fut9", cols = c("grey", "red"), split.by = "intervention") & theme(plot.title = element_blank())

# Fut9 gene, violin plot
VlnPlot(srat, features = "Fut9", sort = FALSE, flip = TRUE, add.noise = FALSE, pt.size = 0, split.by = "intervention") + theme(legend.position = "none", axis.title.x = element_blank(), plot.title = element_blank())

# Check maximum expression levels of Prrc2b isoforms
max(GetAssayData(srat_i, assay = "SCT_i", layer = "data")["ENSMUST00000036691", ])  # Returns 1.79175946922805
max(GetAssayData(srat_i, assay = "SCT_i", layer = "data")["ENSMUST00000069817", ])  # Returns 1.6094379124341

# Prrc2b isoforms, UMAP
FeaturePlot(srat_i, reduction = "umap_harmony", features = "ENSMUST00000036691", cols = c("grey", "red"), split.by = "intervention", order = TRUE) & theme(plot.title = element_blank()) & scale_color_gradientn(colors = c("grey", "red"), limits = c(0, 1.8))
FeaturePlot(srat_i, reduction = "umap_harmony", features = "ENSMUST00000069817", cols = c("grey", "red"), split.by = "intervention", order = TRUE) & theme(plot.title = element_blank()) & scale_color_gradientn(colors = c("grey", "red"), limits = c(0, 1.8))
