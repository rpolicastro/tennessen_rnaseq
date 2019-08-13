#!/usr/bin/env Rscript

library("tidyverse")
library("DESeq2")

###########################
## Making Sample PCA Plots
###########################

## Loading Counts
## ----------

## Load Count Files.

raw_count_files <- list.files(
	file.path("..", "results", "diff_expression"),
	pattern=".*raw.*\\.tsv", full.names = TRUE
)

sample_names <- raw_count_files %>%
	basename %>%
	str_replace("_raw_counts.tsv", "")

counts <- map(
	raw_count_files,
	~read.delim(., sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
		as_tibble(.name_repair = "unique")
)

## Convert counts to single matrix.

counts <- counts %>%
	plyr::join_all(by = "gene_id", type = "left") %>%
	as_tibble(.name_repair = "unique") %>%
	column_to_rownames("gene_id") %>%
	as.matrix

## Creating PCA Plot
## ----------

## Sample info.

sample_sheet <- data.frame(
	condition = c(
		rep("Gpdh.het_Ldh.het", 3),
		rep("Gpdh.mut_Ldh.mut", 3),
		rep("Gpdh.het", 3),
		rep("Gpdh.mut", 3),
		rep("Ldh.mut", 4),
		rep("Ldh.ctrl", 4)
	)
)

rownames(sample_sheet) <- colnames(counts)

## DESeq2 analysis in preparation of PCA plot.

deseq.obj <- counts %>%
	DESeqDataSetFromMatrix(
		countData = .,
		colData = sample_sheet,
		design = ~ condition
	) %>%
	DESeq

## PCA analysis.

# Variance stabalizing transformation.
vsd <- vst(deseq.obj, blind=FALSE)

# Get PCA results.
pca.results <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pca.var <- round(100 * attr(pca.results, "percentVar"))

pca.results$condition <- factor(
	pca.results$condition,
	levels = c(
		"Gpdh.het", "Gpdh.mut",
		"Ldh.ctrl", "Ldh.mut",
		"Gpdh.het_Ldh.het", "Gpdh.mut_Ldh.mut"
	)
)

## Plot PCA.

p <- ggplot(pca.results, aes(x = PC1, y = PC2, color = condition)) +
	geom_point(size = 3) +
	labs(
		x = paste0("PC1: ",pca.var[1],"% variance"),
		y = paste0("PC2: ",pca.var[2],"% variance")
	) +
	theme_bw() +
	scale_color_viridis_d()

dir.create(file.path("..", "results", "pca"), showWarnings = FALSE)

ggsave(
	file.path("..", "results", "pca", "sample_pca.pdf"),
	plot = p, device = cairo_pdf, height = 4, width = 5
)
