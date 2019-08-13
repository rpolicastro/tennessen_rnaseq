#!/usr/bin/env Rscript

library("clusterProfiler")
library("tidyverse")
library("org.Dm.eg.db")

#######################
## Enrichment Analysis
#######################

## Load DEGs
## ----------

## Get file info.

deg_files <- list.files(
	file.path("..", "results", "diff_expression"),
	pattern = ".*_DEGs_only\\.tsv", full.names = TRUE
)

deg_file_names <- deg_files %>%
	basename %>%
	str_replace("_DEGs_only.tsv", "")

## Load files.

degs <- map(
	deg_files,
	~read.delim(., sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
		as_tibble(.name_repair = "unique") %>%
		select(gene_id, logFC, Change)
)

names(degs) <- deg_file_names

## Ontology Analysis
## ----------

## Convert to entrezid.

entrez <- map(degs,
	~ pull(., gene_id) %>%
		bitr(
			geneID = .,
			fromType = "FLYBASE",
			toType = "ENTREZID",
			OrgDb = "org.Dm.eg.db",
			drop = TRUE
		) %>%
		as_tibble(.name_repair = "unique") %>%
		left_join(., .x, by = c("FLYBASE" = "gene_id"))
) %>% bind_rows(.id = "sample")

## Compare clusters.

cluster_comparison <- compareCluster(
	ENTREZID ~ sample + Change,
	data = entrez,
	fun = "enrichGO",
	OrgDb = "org.Dm.eg.db",
	ont = "BP",
	readable = TRUE,
	pAdjustMethod = "fdr"
)

## Plot clusters.

dir.create(file.path("..", "results", "term_enrichment"), showWarnings = FALSE)

p <- dotplot(cluster_comparison, showCategory = 20, x = ~sample) +
	facet_grid(~ Change) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_color_viridis_c(name = "FDR")

ggsave(
	file.path("..", "results", "term_enrichment", "go_cluster_comparison.pdf"),
	plot = p, device = cairo_pdf, height = 15, width = 10
)
