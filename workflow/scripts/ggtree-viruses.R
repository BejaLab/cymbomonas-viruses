
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
            tsv  = "analysis/phylogeny_viruses/MCP-small.tsv"
    ), output = list(
        image = "test.svg",
        jtree = "test.jtree"
    ))
}

with(snakemake@input, {
    tree_file <<- tree
    tsv_file  <<- tsv
})
with(snakemake@output, {
    out_image_file <<- image
    out_jtree_file <<- jtree
})
with(snakemake@params, {
    outgroup_rooting <<- outgroup_rooting
})

orgs <- read.table(tsv_file, col.names = c("label", "Organism"), sep = "\t")

tree <- read.tree(tree_file)
tree <- phangorn::midpoint(tree, node.labels = "support")
if (outgroup_rooting) {
    outgroup_df <- read.fasta.headers(outgroup_file)
    outgroups <- with(outgroup_df, sub(" .*", "", title))
    tree <- ape::root(tree, node = MRCA(tree, outgroups), edgelabel = T, resolve.root = T)
}
tree <- as_tibble(tree) %>%
    mutate(support = ifelse(node %in% parent & label != "", label, NA)) %>%
    separate(support, into = c("SH_aLRT", "UFboot"), sep = "/", convert = T) %>%
    left_join(orgs, by = "label") %>%
    mutate(isInternal = node %in% parent) %>%
    `class<-`(c("tbl_tree", "tbl_df", "tbl", "data.frame"))
ntaxa <- filter(tree, ! node %in% parent) %>% nrow
tree_data <- as.treedata(tree)
write.jtree(tree_data, file = out_jtree_file)

p <- ggtree(tree_data) +
    geom_nodepoint(aes(x = branch, subset = !is.na(UFboot) & UFboot >= 90, size = UFboot)) +
    geom_tiplab(aes(label = Organism), size = 2, align = F, linesize = 0) +
    geom_treescale(width = 0.5) +
    scale_size_continuous(limits = c(90, 100), range = c(1, 2))

ggsave(out_image_file, p, height = ntaxa * 0.1, width = 4, limitsize = F)
