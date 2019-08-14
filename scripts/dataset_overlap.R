#!/usr/bin/env Rscript

library("tidyverse")


#######################
## Overlap of Datasets
#######################

## Loading DEGs
## ----------

## Get DEG files paths.

deg_files <- list.files(
	file.path("..", "results", "diff_expression"),
	pattern = ".*all_genes\\.tsv$",
	full.names = TRUE
)

## Extract sample names.

deg_samples <- deg_files %>%
	basename %>%
	str_replace("_all_genes.tsv", "")

## Import DEGs.

degs <- map(
	deg_files,
	~ read.delim(., sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
		as_tibble(.name_repair = "unique")
)

names(degs) <- deg_samples

## Merge DEG files into one.

degs_renamed <- map2(
	degs, names(degs),
	~ rename_at(.x, vars(logFC:Change), function(x) paste0(x, "_", .y))
)

degs_merged <- degs_renamed %>%
	reduce(
		full_join,
		by = c(
			"gene_id", "gene_name", "chr", "start",
			"end", "strand", "gene_biotype"
		)
	) %>%
	replace_na(list(
		Change_Gpdh = "unchanged",
		Change_Ldh = "unchanged",
		Change_Gpdh_Ldh = "unchanged"
	))

## Export merged deg files.

dir.create(file.path("..", "results", "dataset_overlap"), showWarnings = FALSE)

write.table(
	degs_merged, file.path("..", "results", "dataset_overlap", "merged_dataset.tsv"),
	sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Overlap Information
## ----------

## Get overlap summary.

overlap_summary <- degs_merged %>%
	count(Change_Gpdh, Change_Ldh, Change_Gpdh_Ldh, name = "count") %>%
	mutate(percent = ((count / sum(count)) * 100) %>% round(3))

## Export overlap summary.

dir.create(file.path("..", "results", "dataset_overlap"), showWarnings = FALSE)

write.table(
	overlap_summary, file.path("..", "results", "dataset_overlap", "overlap_summary.tsv"),
	sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)
