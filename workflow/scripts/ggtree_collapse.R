
library(dplyr)
library(tidyr)
library(ape)
library(ggtree)
library(treeio)
library(phangorn)
library(stringr)
library(ggplot2)

#setClass("snake", slots = list(input = "list", output = "list"))
#snakemake <- new("snake", input  = list(
#        tree = "analysis/phylogeny/MCP_PLV.treefile",
#        fasta = Sys.glob("analysis/blast_algae/*-MCP_PLV.faa"),
#        tblout = Sys.glob("analysis/hmm_algae/*.tblout"),
#        synonyms = "annotations/organisms.txt",
#        hmm = Sys.glob("hmm_algae/*.hmm")
#), output = list("test.svg"))

with(snakemake@input, {
	tree_file     <<- tree
	tblout_files  <<- tblout
	fasta_files   <<- fasta
	synonyms_file <<- synonyms
	hmm_files     <<- hmm
})
output_file <- unlist(snakemake@output)

synonyms <- read.table(synonyms_file, header = T, sep = "\t") %>%
	mutate(Len = nchar(Match))

read.hmmer.tblout <- function(fname) {
	col.names <- c(
		"target.name",
		"target.accession",
		"query.name",
		"query.accession",
		"full.sequence.E.value",
		"full.sequence.score",
		"full.sequence.bias",
		"best.1.domain.E.value",
		"best.1.domain.score",
		"best.1.domain.bias",
		"domain.number.estimation.exp",
		"domain.number.estimation.reg",
		"domain.number.estimation.clu",
		"domain.number.estimation.ov",
		"domain.number.estimation.env",
		"domain.number.estimation.dom",
		"domain.number.estimation.rep",
		"domain.number.estimation.inc",
		"description.of.target"
	)
	numeric.cols <- which(col.names == "full.sequence.E.value")        : which(col.names == "domain.number.estimation.exp")
	integer.cols <- which(col.names == "domain.number.estimation.reg") : which(col.names == "domain.number.estimation.inc")
	readLines(fname) %>%
		data.frame(line = .) %>%
		filter(!grepl("^#", line)) %>%
		separate(line, into = col.names, sep = " +", extra = "merge", convert = F) %>%
		mutate_at(numeric.cols, as.numeric) %>%
		mutate_at(integer.cols, as.integer)
}

headers <- file.info(fasta_files) %>%
	filter(size > 0) %>%
	rownames %>%
	lapply(treeio::read.fasta) %>%
	lapply(names) %>%
	unlist %>%
	data.frame(full.label = .) %>%
	extract(full.label, into = c("label", "start", "end", "desc"), regex = "^([^ ]+) \\[(\\d+) - (\\d+)\\](.*)", convert = T, remove = F) %>%
	crossing(synonyms) %>%
	mutate(Start = str_locate(full.label, Match)[,1]) %>%
	arrange(is.na(Start), -Len) %>%
	mutate(label = ifelse(is.na(Start), full.label, label)) %>%
	distinct(label, .keep_all = T)

tree <- read.tree(tree_file) %>%
	midpoint %>%
	as_tibble %>%
	mutate(support = ifelse(node %in% parent & label != "", label, NA)) %>%
	separate(support, into = c("SH_aLRT", "UFboot"), sep = "/", convert = T) %>%
	mutate(contig = sub("_\\d+$", "", label)) %>%
	left_join(headers, by = "label") %>%
	mutate(label.show = Name) %>% # sprintf("%s (%s)", Name, contig)) %>%
	mutate(isInternal = node %in% parent) %>%
	`class<-`(c("tbl_tree", "tbl_df", "tbl", "data.frame"))

hmm <- lapply(hmm_files, readLines) %>%
	lapply(as.data.frame) %>%
	bind_rows %>%
	separate(1, c("key", "value"), sep = "  ", fill = "right", extra = "drop") %>%
	filter(key %in% c("NAME","DESC")) %>%
	mutate(NAME = ifelse(key == "NAME", value, NA)) %>%
	fill(NAME, .direction = "down") %>%
	filter(key == "DESC") %>%
	select(query.name = NAME, DESC = value) %>%
	mutate(DESC = substr(DESC, 1, 50)) %>%
	mutate(category = case_when(
		grepl("Major Capsid Protein", DESC) ~ "MCP",
		grepl("polymerase",           DESC) ~ "Pol",
		grepl("packaging ATPase",     DESC) ~ "A32",
		T ~ "Other")
	)

#genes <- lapply(tblout_files, read.hmmer.tblout) %>%
#	setNames(tblout_files) %>%
#	bind_rows(.id = "fname") %>%
#	mutate(contig = sub("_\\d+$", "", target.name)) %>%
#	extract(description.of.target, into = c("start", "end"), regex = "\\[(\\d+) - (\\d+)\\]", convert = T) %>%
#	filter(full.sequence.E.value < 1e-5) %>%
#	arrange(-full.sequence.score) %>%
#	distinct(target.name, .keep_all = T) %>%
#	left_join(tree, by = "contig") %>%
#	left_join(hmm, by = "query.name") %>%
#	filter(!is.na(label)) %>%
#	filter(pmin(abs(start.x - start.y), abs(end.x - end.y)) < 200000) %>%
#	select(label, DESC, query.name, category)

ntaxa <- filter(tree, ! node %in% parent) %>% nrow

colors <- list(
	Haptophyta = "orange",
	Chlorophyta = "green",
	Streptophyta = "darkgreen",
	PLV = "purple",
	Ochrophyta = "brown",
	Cryptophyta = "red"
)

tree_data <- as.treedata(tree)

multi_scale <- function(p, multi) {
	with(multi, Reduce(function(.x, .y) {
		scaleClade(.x, .y, 1 / (num_tips[node == .y] - 1))
	}, node, p, accumulate = T)) %>%
	`[[`(length(.))
}
multi_collapse <- function(p, multi) {
	with(multi, Reduce(function(.x, .y) {
		fill <- unlist(colors[Group[node == .y]])
		.x$data[.x$data$node == .y, "label.show"] <- label.show[node == .y]
		.x <- collapse(.x, .y, "mixed", fill = fill)
		return(.x)
	}, node, p, accumulate = T)) %>%
	`[[`(length(.))
}

multi_species <- allDescendants(tree_data@phylo) %>%
	lapply(function(x) filter(tree, node %in% x)) %>%
	bind_rows(.id = "ancestor") %>%
	group_by(ancestor) %>%
	mutate(Group = first(na.omit(Group))) %>%
	filter(n_distinct(Name, na.rm = T) == 1, sum(!isInternal) > 2, !any(Group == "Haptophyta")) %>%
	ungroup %>%
	filter(! ancestor %in% node) %>%
	group_by(ancestor, Group) %>%
	summarize(num_tips = sum(!isInternal), Name = first(na.omit(Name)), label.show = sprintf("%s (%d)", Name, num_tips)) %>%
	mutate(node = as.numeric(ancestor))

p <- ggtree(as.treedata(tree)) +
	geom_nodepoint(aes(x = branch, subset = !is.na(UFboot) & UFboot >= 90, size = UFboot)) +
	geom_tiplab(aes(label = label.show), size = 3, align = T, linesize = 0) +
	geom_text2(aes(subset = node %in% multi_species$node, x = max(x, na.rm = T), label = label.show), nudge_x = 0.01, size = 3, hjust = 0) +
	geom_tippoint(aes(color = Group), size = 3) +
	geom_treescale(width = 0.5) +
	scale_size_continuous(range = c(1, 3)) +
	scale_shape_manual(values = seq(0,15)) +
	scale_color_manual(values = colors)
	# geom_hilight()

p <- multi_scale(p, multi_species)
p <- multi_collapse(p, multi_species)

# p <- facet_plot(p, mapping = aes(x = as.numeric(as.factor(query.name)), shape = DESC), data = genes, geom = geom_point, panel = 'Genes')

ggsave(output_file, p, height = ntaxa * 0.1, width = 7, limitsize = F)

