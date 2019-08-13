#!/usr/bin/env Rscript

library("clusterProfiler")
library("tidyverse")
library("org.Dm.eg.db")
library("ReactomePA")

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
		dplyr::select(gene_id, logFC, Change)
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

go_categories <- c("BP", "MF", "CC")

go_enrichment <- map(
	go_categories,
	~compareCluster(
		ENTREZID ~ sample + Change,
		data = entrez,
		fun = "enrichGO",
		OrgDb = "org.Dm.eg.db",
		ont = .,
		readable = TRUE,
		pAdjustMethod = "fdr"
	)
)

names(go_enrichment) <- go_categories

## Export tibble as table.

map2(
	go_enrichment, names(go_enrichment),
	~ as_tibble(.x, .name_repair = "universal") %>%
		dplyr::rename("FDR" = p.adjust) %>%
		write.table(
			., file.path("..", "results", "term_enrichment", paste0("GO_", .y, "_enrichment_table.tsv")),
			sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
		)
)

## Plot clusters.

dir.create(file.path("..", "results", "term_enrichment"), showWarnings = FALSE)

enrichment_dotplot <- function(x, y, cats) {
	p <- dotplot(x, showCategory = cats, x = ~ sample) +
		facet_grid(~ Change) +
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		scale_color_viridis_c(name = "FDR")

	ggsave(
		file.path("..", "results", "term_enrichment", paste0("GO_", y, "_enrichment_dotplot.pdf")),
		plot = p, device = cairo_pdf, height = 15, width = 10
	)
}

map2(go_enrichment, names(go_enrichment), ~enrichment_dotplot(.x, .y, 20))

## ReactomeDB Analysis
## ----------

## ReactomeDB pathway enrichment.

pathway_enrichment <- compareCluster(
	ENTREZID ~ sample + Change,
	data = entrez,
	fun = "enrichPathway",
	organism = "fly",
	pAdjustMethod = "fdr",
	readable = TRUE
)

## Export table of results.

pathway_enrichment %>%
	as_tibble(.name_repair = "universal") %>%
	dplyr::rename("FDR" = p.adjust) %>%
	write.table(
		., file.path("..", "results", "term_enrichment", "REACTOME_enrichment_table.tsv"),
		sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
	)

## Plot reactome pathway dotplot.

p <- dotplot(pathway_enrichment, showCategory = 25, x = ~ sample) +
	facet_grid(~ Change) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_color_viridis_c(name = "FDR")

ggsave(
	file.path("..", "results", "term_enrichment", "REACTOME_enrichment_dotplot.pdf"),
	plot = p, device = cairo_pdf, height = 15, width = 15
)

## KEGG Analysis.
## ----------

## Convert IDs to uniprot.

uniprot <- map(degs,
	~ pull(., gene_id) %>%
		bitr(
			geneID = .,
			fromType = "FLYBASE",
			toType = "UNIPROT",
			OrgDb = "org.Dm.eg.db",
			drop = TRUE
		) %>%
		as_tibble(.name_repair = "unique") %>%
		left_join(., .x, by = c("FLYBASE" = "gene_id"))
) %>% bind_rows(.id = "sample")

## KEGG pathway enrichment.

kegg_enrichment <- compareCluster(
	UNIPROT ~ sample + Change,
	data = uniprot,
	fun = "enrichKEGG",
	organism = "dme",
	keyType = "uniprot",
	pAdjustMethod = "fdr"
)

## Export table of results.

kegg_enrichment %>%
	as_tibble(.name_repair = "universal") %>%
	dplyr::rename("FDR" = p.adjust) %>%
	write.table(
		., file.path("..", "results", "term_enrichment", "KEGG_enrichment_table.tsv"),
		sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
	)

## Plot KEGG pathway enrichment.

p <- dotplot(kegg_enrichment, showCategory = 20, x = ~ sample) +
	facet_grid(~ Change) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_color_viridis_c(name = "FDR")

ggsave(
	file.path("..", "results", "term_enrichment", "KEGG_enrichment_dotplot.pdf"),
	plot = p, device = cairo_pdf, height = 10, width = 10
)
