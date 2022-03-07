library(dplyr)
library(tidyr)
library(ape)
library(ggtree)
library(treeio)
library(phangorn)
library(stringr)
library(ggplot2)

input  <- snakemake@input
output <- snakemake@output

tree_file     <- input$tree
tblout_files  <- input$tblout
fasta_files   <- input$fasta
synonyms_file <- input$synonyms
output_file   <- unlist(output)

hmm_files     <- input$hmm

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
	mutate(label.show = sprintf("%s (%s)", Name, contig)) %>%
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

genes <- lapply(tblout_files, read.hmmer.tblout) %>%
	setNames(tblout_files) %>%
	bind_rows(.id = "fname") %>%
	mutate(contig = sub("_\\d+$", "", target.name)) %>%
	extract(description.of.target, into = c("start", "end"), regex = "\\[(\\d+) - (\\d+)\\]", convert = T) %>%
	filter(full.sequence.E.value < 1e-5) %>%
	arrange(-full.sequence.score) %>%
	distinct(target.name, .keep_all = T) %>%
	left_join(tree, by = "contig") %>%
	left_join(hmm, by = "query.name") %>%
	filter(!is.na(label)) %>%
	filter(pmin(abs(start.x - start.y), abs(end.x - end.y)) < 200000) %>%
	select(label, DESC, query.name, category)

ntaxa <- filter(tree, ! node %in% parent) %>% nrow

colors <- list(
	Haptophyta = "orange",
	Chlorophyta = "green",
	Streptophyta = "darkgreen",
	PLV = "purple",
	Ochrophyta = "brown",
	Cryptophyta = "red"
)

p <- ggtree(as.treedata(tree)) +
	geom_nodepoint(aes(x = branch, subset = !is.na(UFboot) & UFboot >= 90, size = UFboot)) +
	geom_tiplab(aes(label = label.show), size = 5, align = T) +
	geom_tippoint(aes(color = Group), size = 3) +
	geom_treescale(width = 0.5) +
	scale_size_continuous(range = c(1, 3)) +
	scale_shape_manual(values = seq(0,15)) +
	scale_color_manual(values = colors)
p <- facet_plot(p, mapping = aes(x = as.numeric(as.factor(query.name)), shape = DESC), data = genes, geom = geom_point, panel = 'Genes')

ggsave(output_file, p, height = ntaxa * 0.2, width = 10, limitsize = F)
