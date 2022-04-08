
library(dplyr)
library(tidyr)
library(ape)
library(ggtree)
library(treeio)
library(phangorn)
library(stringr)
library(ggplot2)
library(phytools)

#setClass("snake", slots = list(input = "list", output = "list"))
#snakemake <- new("snake", input  = list(
#        tree = "analysis/phylogeny/MCP_PLV.treefile",
#        fasta_algae  = Sys.glob("analysis/phylogeny/algae/hmm/*-MCP_PLV.faa"),
#        fasta_PLVs   = Sys.glob("analysis/phylogeny/PLVs/hmm/*/*-MCP_PLV.faa"),
#        blast_algae  = Sys.glob("analysis/phylogeny/algae/blast/*-MCP_PLV.blast"),
#        blast_PLVs   = Sys.glob("analysis/phylogeny/PLVs/blast/*/*-MCP_PLV.blast"),
#        cp_queries   = "analysis/CPs/queries/MCP_PLV.faa",
#        root_file    = "annotations/MCP_PLV_root.txt",
#        neighbors_file = ""
#        synonyms = "annotations/organisms.txt",
#        hmm = Sys.glob("hmm_algae/*.hmm")
#), output = list("test.svg"))

with(snakemake@input, {
	tree_file     <<- tree
	fasta_files   <<- c(fasta_PLVs, fasta_algae)
	synonyms_file <<- synonyms
	hmm_files     <<- hmm
	blast_files   <<- c(blast_algae, blast_PLVs)
	cp_file       <<- cp_queries
	root_file     <<- root_file
})
output_file <- unlist(snakemake@output)

outgroups <- readLines(root_file)

synonyms <- read.table(synonyms_file, header = T, sep = "\t", fill = T, na.strings = "") %>%
	mutate(Len = nchar(Match)) %>%
	mutate(Collapse = ifelse(is.na(Collapse), Name, Collapse))
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

read.outfmt6 <- function(fname, col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"), extra = c()) {
	tbl <- read.table(fname, sep = "\t", col.names = c(col.names, extra)) %>%
	`if`(nrow(.) > 0, ., NULL)
}

cp_queries <- treeio::read.fasta(cp_file) %>%
	{data.frame(title = names(.))} %>%
	extract(title, into = c("qseqid", "clade.color"), regex = "([^ ]+) .*group=([^ ]+)") %>%
	filter(!is.na(clade.color))

headers <- file.info(fasta_files) %>%
	filter(size > 0) %>%
	rownames %>%
	lapply(treeio::read.fasta) %>%
	lapply(names) %>%
	unlist %>%
	data.frame(full.label = .) %>%
	extract(full.label, into = c("label", "desc"), regex = "^([^ ]+) ?(.+)?", remove = F) %>%
	extract(desc, into = c("start", "end"), regex = "\\[(\\d+) - (\\d+)\\]", convert = T, remove = F) %>%
	crossing(synonyms) %>%
	mutate(Start = str_locate(full.label, Match)[,1], End = str_locate(full.label, Match)[,2]) %>%
	arrange(is.na(Start), Start, Start - End) %>%
	distinct(full.label, .keep_all = T)

no_name <- filter(headers, is.na(Name)) %>%
	pull(label) %>%
	paste(collapse = ", ")
if (no_name > 0) {
	write(paste("No aliases found for: ", no_name))
	quit(status = 1)
}

tree <- read.tree(tree_file)
if (length(outgroups) > 0) {
	tree <- reroot(tree, MRCA(tree, outgroups), position = 0.1, resolve.root = T)
} else {
	tree <- midpoint(tree)
}
tree <- as_tibble(tree) %>%
	mutate(support = ifelse(node %in% parent & label != "", label, NA)) %>%
	separate(support, into = c("SH_aLRT", "UFboot"), sep = "/", convert = T) %>%
	mutate(contig = sub("_\\d+$", "", label)) %>%
	left_join(headers, by = "label") %>%
	mutate(label.show = ifelse(is.na(Name), desc, Name)) %>% # sprintf("%s (%s)", Name, contig)) %>%
	mutate(isInternal = node %in% parent) %>%
	`class<-`(c("tbl_tree", "tbl_df", "tbl", "data.frame"))
tree_data <- as.treedata(tree)

protein_blast <- lapply(blast_files, read.outfmt6, extra = "stitle") %>%
	bind_rows %>%
	filter(evalue < 1e-40, pident > 52) %>%
	left_join(cp_queries, by = "qseqid") %>%
	filter(!is.na(clade.color), sseqid %in% tree_data@phylo$tip.label) %>%
        select(label = sseqid, clade.color) %>%
        distinct(label, .keep_all = T) %>%
	group_by(clade.color) %>%
	group_map(function(x, y) {
		if (length(x$label)) data.frame(node = MRCA(tree_data@phylo, x$label), clade.color = y)
	}) %>%
	bind_rows

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
#	filter(pmin(abs(start.p - start.node), abs(end.p - end.node)) < 200000) %>%
#	select(label, DESC, query.name, category)

ntaxa <- filter(tree, ! node %in% parent) %>% nrow

colors <- list(
	Haptophyta = "orange",
	Chlorophyta = "green",
	Streptophyta = "darkgreen",
	MAG = "purple",
	Stramenopiles = "brown",
	Cryptophyta = "red",
	Amoebozoa = "gold4",
	Euglenozoa = "yellow",
	Choanoflagellata = "darkslateblue",
	Glaucophyta = "cyan",
	Animals = "blue",
	Dinoflagellata = "gray",
	Rhizaria = "gray"
)

scaleClades <- function(p, df) {
	with(df, Reduce(function(.p, .node) {
		offs <- offspring(.p$data, .node)
		scale <- 1 / (nrow(offs) - 1)
		scaleClade(.p, .node, scale)
	}, node, p))
}
collapseClades <- function(p, df) {
	with(df, Reduce(function(.p, .node) {
		fill <- unlist(colors[Host[node == .node]])
		.p$data[.p$data$node == .node, "label.show"] <- label.show[node == .node]
		collapse(.p, .node, "mixed", height = 1, fill = fill)
	}, node, p))
}
labelClades <- function(p) {
	with(df, Reduce(function(.p, .node) {
		.p + geom_cladelab(node = .node, label = label[node == .node], align = T, offset = .2, textcolor = 'blue')
	}, node, p))
}

multi_species <- allDescendants(tree_data@phylo) %>%
	lapply(function(x) filter(tree, node %in% x)) %>%
	bind_rows(.id = "ancestor") %>%
	group_by(ancestor) %>%
	filter(n_distinct(Collapse, na.rm = T) == 1, sum(!isInternal) > 2) %>% # , !any(Group == "Haptophyta")) %>%
	ungroup %>%
	mutate(ancestor = as.numeric(ancestor)) %>%
	filter(! ancestor %in% node) %>%
	filter(!is.na(Collapse)) %>%
	group_by(ancestor, Collapse) %>%
	summarize(num_tips = sum(!isInternal), Host = first(na.omit(Host))) %>%
	mutate(label.show = sprintf("%s (%d)", Collapse, num_tips)) %>%
	rename(node = ancestor)
p <- ggtree(tree_data) +
	geom_nodepoint(aes(x = branch, subset = !is.na(UFboot) & UFboot >= 90, size = UFboot)) +
	geom_tiplab(aes(label = label.show), size = 4, align = T, linesize = 0) +
	geom_text2(aes(subset = node %in% multi_species$node, x = max(x, na.rm = T), label = label.show), nudge_x = 0.01, size = 4.5, hjust = 0) +
	geom_tippoint(aes(color = Host), size = 3) +
	geom_treescale(width = 0.5) +
	scale_size_continuous(range = c(1, 3)) +
	scale_shape_manual(values = seq(0,15)) +
	scale_color_manual(values = colors)

p <- scaleClades(p, multi_species)
p <- collapseClades(p, multi_species)
if (length(protein_blast) > 0) {
	p <- p + geom_hilight(data = protein_blast, aes(node = node, fill = clade.color))
}

# p <- facet_plot(p, mapping = aes(x = as.numeric(as.factor(query.name)), shape = DESC), data = genes, geom = geom_point, panel = 'Genes')

ggsave(output_file, p, height = ntaxa * 0.1, width = 7, limitsize = F)

