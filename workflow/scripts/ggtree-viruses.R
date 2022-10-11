
library(dplyr)
library(tidyr)
library(ape)
library(ggtree)
library(treeio)
library(phangorn)
library(stringr)
library(ggplot2)
library(phytools)

if (interactive()) {
    setClass("snake", slots = list(input = "list", output = "list"))
    snakemake <- new("snake", input  = list(
            tree = "analysis/phylogeny_viruses/MCP-small.treefile",
            tsv  = "analysis/phylogeny_viruses/MCP-small.tsv",
            virus_metadata = "metadata/viruses.tsv",
            segment_metadata = "metadata/viral_segments.tsv",
            colors = "metadata/subclade_colors.txt"
    ), output = list(
        image = "test.svg",
        jtree = "test.jtree"
    ))
}

with(snakemake@input, {
    tree_file <<- tree
    tsv_file  <<- tsv
    virus_metadata_file <<- virus_metadata
    segment_metadata_file <<- segment_metadata
    colors_file <<- colors
})
with(snakemake@output, {
    out_image_file <<- image
    out_jtree_file <<- jtree
})

colors <- read.table(colors_file, col.names = c("subclade", "color"), comment.char = "") %>%
    with(setNames(color, subclade))
viruses  <- read.table(virus_metadata_file, header = T, sep = "\t", na.strings = "")
segments <- read.table(segment_metadata_file, header = T, sep = "\t", na.strings = "")
metadata <- bind_rows(viruses, segments) %>%
    mutate(subclade = ifelse(is.na(subclade), clade, subclade))

orgs <- read.table(tsv_file, col.names = c("label", "Organism"), sep = "\t")

tree <- read.tree(tree_file) %>%
    phangorn::midpoint(node.labels = "support") %>%
    as_tibble %>%
    mutate(support = ifelse(node %in% parent & label != "", label, NA)) %>%
    separate(support, into = c("SH_aLRT", "UFboot"), sep = "/", convert = T) %>%
    left_join(orgs, by = "label") %>%
    left_join(metadata, by = c(Organism = "short")) %>%
    mutate(isInternal = node %in% parent) %>%
    `class<-`(c("tbl_tree", "tbl_df", "tbl", "data.frame"))
ntaxa <- filter(tree, ! node %in% parent) %>% nrow
tree_data <- as.treedata(tree)
write.jtree(tree_data, file = out_jtree_file)

p <- ggtree(tree_data) +
    geom_nodepoint(aes(x = branch, subset = !is.na(UFboot) & UFboot >= 90, size = UFboot)) +
    geom_tiplab(aes(label = Organism, color = subclade), size = 2, align = F, linesize = 0) +
    scale_color_manual(values = colors) +
    geom_treescale(width = 0.5) +
    scale_size_continuous(limits = c(90, 100), range = c(1, 2))

ggsave(out_image_file, p, height = ntaxa * 0.1, width = 4, limitsize = F)
