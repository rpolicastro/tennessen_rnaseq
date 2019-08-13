#!/usr/env/bin Rscript

library("tidyverse")
library("rtracklayer")
library("GenomicRanges")
library("edgeR")

###########################
## Differential Expression
###########################

## Biomart Gene Info
## ----------

biomart <- read.delim(
	"../files/Drosophila_melanogaster.BDGP6.22.97_biomart.tsv",
	sep = "\t", header = TRUE, stringsAsFactors = FALSE
) %>%
	as_tibble(.name_repair = "universal") %>%
	dplyr::rename(
		gene_id = Gene.stable.ID,
		gene_name = Gene.name,
		chr = Chromosome.scaffold.name,
		start = Gene.start..bp.,
		end = Gene.end..bp.,
		strand = Strand,
		gene_biotype = Gene.type
	) %>%
	mutate(strand = case_when(
		strand == 1 ~ "+",
		strand == -1 ~ "-"
	))

## Genes to Filter
## ----------

## Loading annotations file.

annotation <- import("../genome/Drosophila_melanogaster.BDGP6.22.97.gtf", "gtf")

## Getting list of pseudogenes, tRNA, and rRNA.

to_filter <- annotation %>%
	as_tibble(.name_repair = "unique") %>%
	dplyr::select(gene_id, gene_biotype) %>%
	filter(gene_biotype %in% c("pseudogene", "tRNA", "rRNA")) %>%
	distinct(gene_id, gene_biotype)

## Preparing Counts File
## ----------

## Loading counts file.

# Load counts as tibble.
counts <- read.delim("../files/counts_merged.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
	as_tibble(.name_repair = "universal")

## Cleaning counts file.

# Remove the directory part of the column names.
colnames(counts)[2:ncol(counts)] <- str_replace(
	colnames(counts)[2:ncol(counts)],
	"X.N.dc2.scratch.rpolicas.jason_collab.results.aligned.",
	""
)

# Remove the bam file part of the column names.
colnames(counts)[2:ncol(counts)] <- str_replace(
        colnames(counts)[2:ncol(counts)],
        "_Aligned.sortedByCoord.out.bam",
        ""
)

# Properly order the columns.
counts <- counts %>%
	dplyr::select(tidyselect::peek_vars() %>% sort) %>%
	dplyr::select(Geneid, everything()) %>%
	dplyr::rename(gene_id = Geneid)

## Splitting samples.

counts <- counts %>%
	gather(key = "sample", value = "counts", -gene_id) %>%
	mutate(group = case_when(
		grepl(sample, pattern = "_Gpdh.(het|mut)_\\d") ~ "Gpdh",
		grepl(sample, pattern = "_Gpdh.(het|mut)_Ldh.(het|mut)_") ~ "Gpdh_Ldh",
		grepl(sample, pattern = "_(Pre|Delta)_") ~ "Ldh"
	))

group_names <- counts %>% group_keys(group) %>% pull(group)

counts <- counts %>%
	group_split(group) %>%
	setNames(group_names)

## Making count matrix.

counts <- map(
	counts,
	~dplyr::select(., -group) %>%
		spread(key = sample, value = counts) %>%
		anti_join(., to_filter, by = "gene_id") %>%
		column_to_rownames("gene_id") %>%
		as.matrix
)

## Export raw count matrices.

dir.create("../results", showWarnings = FALSE)
dir.create("../results/diff_expression", showWarnings = FALSE)

map2(
	counts, names(counts),
	~as_tibble(.x, .name_repair = "unique", rownames = "gene_id") %>%
		write.table(
			., file.path("..", "results", "diff_expression", paste0(.y, "_raw_counts.tsv")),
			sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
		)
)

## Differential Expression
## ----------

## Establishing sample groups and designs.

samp_groups = list(
	"Gpdh" = factor(c(rep(1, 3), rep(2, 3)), levels = c(1, 2)),
	"Gpdh_Ldh" = factor(c(rep(1, 3), rep(2, 3)), levels = c(1, 2)),
	"Ldh" = factor(c(rep(2, 4), rep(1, 4)), levels = c(1, 2))
)

samp_design <- map(samp_groups, function(group) model.matrix(~ group))

## Getting Gene TMM expression values.

TMM <- map(
	counts,
	~DGEList(.) %>%
		calcNormFactors %>%
		cpm %>%
		as_tibble(.name_repair = "unique", rownames = "gene_id")
)

## Differential Expression of genes.

DEGs <- map2(
	counts, names(counts),
	~DGEList(.x, group = samp_groups[[.y]]) %>%
		.[filterByExpr(.), , keep.lib.sizes = FALSE] %>%
		calcNormFactors %>%
		estimateDisp(design = samp_design[[.y]]) %>%
		glmQLFit(design = samp_design[[.y]]) %>%
		glmQLFTest(coef=2)
)

## Prepare DEGs for export.

DEGs.df <- map2(
	DEGs, TMM,
	~.x$table %>%
		as_tibble(.name_repair = "unique", rownames = "gene_id") %>%
		mutate(
			"FDR" = p.adjust(PValue, method = "fdr"),
			"Change" = case_when(
				FDR <= 0.05 & logFC >= 1 ~ "upregulated",
				FDR <= 0.05 & logFC <= -1 ~ "downregulated",
				TRUE ~ "unchanged"
			)
		) %>%
		left_join(., .y, by = "gene_id") %>%
		inner_join(biomart, ., by = "gene_id") %>%
		arrange(desc(logFC))
)

## Export DEGs.

# Number of DEGs.
n_DEGs <- map(DEGs.df, ~count(., Change, name = "count", sort = TRUE))

map2(
	n_DEGs, names(n_DEGs),
	~write.table(
		.x, file.path("..", "results", "diff_expression", paste0(.y, "_DEG_numbers.tsv")),
		sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
	)
)

# Status of all genes.
map2(
	DEGs.df, names(DEGs.df),
	~write.table(
		.x, file.path("..", "results", "diff_expression", paste0(.y, "_all_genes.tsv")),
		sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
	)
)

# Only DEGs.
map2(
	DEGs.df, names(DEGs.df),
	~filter(.x, Change == "downregulated" | Change == "upregulated") %>%
		write.table(
			., file.path("..", "results", "diff_expression", paste0(.y, "_DEGs_only.tsv")),
			sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
		)
)

