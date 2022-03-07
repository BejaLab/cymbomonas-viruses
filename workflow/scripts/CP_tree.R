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

tree1_file    <- unlist(input["tree_MCP"])
tree2_file    <- unlist(input["tree_mCP"])
fasta_files   <- unlist(input["fasta"])
synonyms_file <- unlist(input["synonyms"])
output_file   <- unlist(output)

synonyms <- read.table(synonyms_file, header = T, sep = "\t") %>%
	mutate(Len = nchar(Match))

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
	scale_color_manual(values = colors)
ggsave(output_file, p, height = ntaxa * 0.2, width = 10)
