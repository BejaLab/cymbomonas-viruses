
library(dplyr)
library(tidyr)
library(igraph)
library(ggplot2)
# library(gggenes)

input <- snakemake@input
output <- snakemake@output
params <- snakemake@params

hmmsearch_files <- unlist(input)

# figure_file <- unlist(output["figure"])
bed_file <- unlist(output["bed"])

E.value.threshold <- unlist(params["e_value_threshold"]) # 0.001
distance.threshold <- unlist(params["distance_threshold"]) # 200000
genes.threshold <- unlist(params["genes_threshold"]) # 2

read.hmmer.tblout <- function(fname) {
	col.names <- c("target.name", "target.accession", "query.name", "query.accession", "full.sequence.E.value", "full.sequence.score", "full.sequence.bias", "best.1.domain.E.value", "best.1.domain.score", "best.1.domain.bias", "domain.number.estimation.exp", "domain.number.estimation.reg", "domain.number.estimation.clu", "domain.number.estimation.ov", "domain.number.estimation.env", "domain.number.estimation.dom", "domain.number.estimation.rep", "domain.number.estimation.inc", "description.of.target")
	numeric.cols <- which(col.names == "full.sequence.E.value")        : which(col.names == "domain.number.estimation.exp")
	integer.cols <- which(col.names == "domain.number.estimation.reg") : which(col.names == "domain.number.estimation.inc")
	readLines(fname) %>%
		data.frame(line = .) %>%
		filter(!grepl("^#", line)) %>%
		separate(line, into = col.names, sep = " +", extra = "merge", convert = F) %>%
		mutate_at(numeric.cols, as.numeric) %>%
		mutate_at(integer.cols, as.integer)
}

data <- lapply(hmmsearch_files, read.hmmer.tblout) %>%
	bind_rows %>%
	extract(description.of.target, into = c("start", "end"), regex = "(\\d+) - (\\d+)", convert = T) %>%
	extract(target.name, into = "scaffold", regex = "^(.+)_\\d+$", remove = F) %>%
	extract(query.name, into = "gene", regex = "^([^_]+)", remove = F) %>%
	arrange(best.1.domain.E.value) %>%
	select(target.name, scaffold, query.name, gene, best.1.domain.E.value, best.1.domain.score, start, end) %>%
	filter(best.1.domain.E.value < E.value.threshold) %>%
	distinct(target.name, .keep_all = T)

relations <- left_join(data, data, by = "scaffold") %>%
	filter(start.x < start.y) %>%
	mutate(distance = pmin(abs(start.x - start.y), abs(start.x - end.y), abs(end.x - start.y), abs(end.x - end.y))) %>%
	filter(distance < distance.threshold) %>%
	arrange(start.x, start.y) %>%
	distinct(target.name.x, .keep_all = T) %>%
	select(target.name.x, target.name.y, scaffold, start.x, start.y, distance)

g <- graph_from_data_frame(relations, vertices = data)
V(g)$membership <- components(g)$membership

genes <- as_data_frame(g, "vertices") %>%
	mutate(molecule = paste(scaffold, membership)) %>%
	group_by(molecule) %>%
	mutate(n_genes = n_distinct(gene)) %>%
	filter(n_genes > 1)

group_by(genes, scaffold, membership) %>%
	summarize(start = min(start, end) - 1, end = max(start, end), .groups = "drop") %>%
	select(scaffold, start, end) %>%
	write.table(bed_file, row.names = F, col.names = F, quote = F, sep = "\t")

#p <- ggplot(genes, aes(xmin = start, xmax = end, x = (start + end) / 2, y = molecule, fill = gene, label = query.name)) +
# 	geom_gene_arrow() +
#	geom_text(vjust = -2, size = 2) +
#	facet_wrap(~ molecule, scales = "free", ncol = 1) +
#	scale_fill_brewer(palette = "Set3") +
#	xlab("ORF coordinate") +
#	scale_x_continuous(expand = c(.2, .2)) +
#	theme_genes()
#ggsave(figure_file, p, height = length(unique(genes$molecule)))
