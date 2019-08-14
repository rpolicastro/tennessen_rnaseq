#!/usr/bin/env Rscript

library("tidyverse")

#####################
## DEG Volcano Plots
#####################

## Load DEGs
## ----------

## Get deg files names and paths.

deg_files <- list.files(
	file.path("..", "results", "diff_expression"),
	pattern = ".*all_genes\\.tsv$",
	full.names = TRUE
)

## Extract sample names from DEG files.

deg_samples <- deg_files %>%
	basename %>%
	str_replace("_all_genes.tsv", "")

## Load DEG files.

degs <- map(
	deg_files,
	~read.delim(., sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
		as_tibble(.name_repair = "unique")
)

names(degs) <- deg_samples

## Volcano Plot
## ----------

dir.create(file.path("..", "results", "volcano_plots"), showWarnings = FALSE)

plot_volcano <- function(x, y) {
	p <- ggplot(x, aes(x = logFC, y = -log10(FDR))) +
		geom_point(aes(color = Change), size = 0.75) +
		scale_color_viridis_d() +
		theme_bw() +
		geom_vline(xintercept = -1, lty = 2) +
		geom_vline(xintercept = 1, lty = 2) +
		geom_hline(yintercept = -log10(0.05), lty = 2)

	ggsave(
		file.path("..", "results", "volcano_plots", paste0(y, "_volcano_plot.pdf")),
		plot = p, device = cairo_pdf, height = 4, width = 8
	)
}

map2(degs, names(degs), ~plot_volcano(.x, .y))
