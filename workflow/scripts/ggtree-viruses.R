
library(dplyr)
library(tidyr)
library(ape)
library(ggtree)
library(treeio)
library(phangorn)
library(stringr)
library(ggplot2)
library(phytools)
library(tibble)

read.cdhit.clstr <- function(fname) {
    data.fields <- c("E.Value", "Aln", "Identity")
    read.table(fname, sep = "\t", comment.char = "", quote = "", fill = T, stringsAsFactors = F, col.names = c("Col1", "Col2")) %>%
        separate(Col1, into = c("Seq.Num", "Cluster"), sep = " ", fill = "right") %>%
        fill(Cluster) %>%
        filter(!grepl(">", Seq.Num)) %>%
        separate(Col2, into = c("Seq.Len", "Col2"), sep = "aa, >") %>%
        extract(Col2, into = c("Seq.Name", "Is.Representative", "Col2"), regex = "(.*?)[.]{3} ([*]|at) ?(.*)") %>%
        mutate(Is.Representative = Is.Representative == "*", Col2 = ifelse(Is.Representative, "100%", Col2)) %>%
        group_by(Cluster) %>%
        mutate(Representative = Seq.Name[which(Is.Representative)]) %>%
        separate_rows(Col2, sep = ",") %>%
        separate(Col2, into = data.fields, sep = "/", fill = "left", convert = T) %>%
        mutate(Identity = sub("%", "", Identity) %>% as.numeric)
}

with(snakemake@input, {
    tree_file <<- tree
    protein_files <<- proteins
    viruses_file <<- viruses
    img_file     <<- img
    clstr_files  <<- clstr
    lens_files   <<- lens
    aln_files    <<- alns
    markers_file <<- markers
})
with(snakemake@output, {
    out_image_file <<- image
    out_jtree_file <<- jtree
})
with(snakemake@params, {
    output_width <<- width
})

markers <- read.table(markers_file, col.names = c("marker", "desc", "score"))
lens <- lapply(lens_files, read.table, col.names = c("short", "contig", "len")) %>%
    bind_rows %>%
    group_by(short) %>%
    summarize(total_len = sum(len))
alns <- lapply(aln_files, read.fasta) %>%
    lapply(names) %>%
    setNames(aln_files) %>%
    lapply(as.data.frame) %>%
    bind_rows(.id = "marker") %>%
    setNames(c("marker", "short")) %>%
    mutate(marker = basename(tools::file_path_sans_ext(marker)), short = sub(" .+", "", short)) %>%
    left_join(markers, by = "marker")
img_locations <- read.table(img_file, sep = "\t", col.names = c("Analysis", "Location", "Status"))
clstr <- lapply(clstr_files, read.cdhit.clstr) %>%
    bind_rows %>%
    group_by(Representative) %>%
    summarize(n_clstr = n())
viruses  <- read.table(viruses_file, header = T, sep = "\t", na.strings = "")
metadata <- protein_files %>%
    `[`(file.size(.)>0) %>%
    lapply(read.table, sep = "\t") %>%
    bind_rows %>%
    setNames(c("short", "protein")) %>%
    left_join(viruses, by = "short") %>%
    left_join(lens, by = "short") %>%
    distinct(protein, .keep_all = T)

tree <- read.tree(tree_file) %>%
    phangorn::midpoint(node.labels = "support") %>%
    as_tibble %>%
    mutate(support = ifelse(node %in% parent & label != "", label, NA)) %>%
    separate(support, into = c("SH_aLRT", "UFboot"), sep = "/", convert = T) %>%
    left_join(metadata, by = c(label = "protein")) %>%
    left_join(clstr, by = c(label = "Representative")) %>%
    extract(label, into = "IMGVR", regex = "IMGVR_UViG_[\\d_]+\\|\\d+\\|(.+)_\\d+$", remove = F) %>%
    extract(IMGVR, into = "Analysis", regex = "^([^_]+)", remove = F) %>%
    left_join(img_locations, by = "Analysis") %>%
    mutate(Label = case_when(
	!is.na(IMGVR) ~ paste(IMGVR, ifelse(is.na(Location), "", sprintf("(%s)", Location))),
        # !is.na(name)  ~ name,
        !is.na(short) ~ short,
        !is.na(n_clstr) ~ sprintf('(%d)', n_clstr),
        T ~ label)
    ) %>%
    mutate(isInternal = node %in% parent) %>%
    `class<-`(c("tbl_tree", "tbl_df", "tbl", "data.frame"))
ntaxa <- filter(tree, ! node %in% parent) %>% nrow
tree_data <- as.treedata(tree)
write.jtree(tree_data, file = out_jtree_file)

lens_df <- select(metadata, id = short, x = total_len)
mark_df <- select(alns, short, desc) %>%
    mutate(present = T) %>%
    spread(desc, present) %>%
    column_to_rownames("short")
p <- ggtree(tree_data) +
    geom_nodepoint(aes(x = branch, subset = !is.na(UFboot) & UFboot >= 90, size = UFboot)) +
    geom_tiplab(aes(label = Label), size = 2, align = F, linesize = 0) +
    # scale_color_manual(values = colors) +
    geom_treescale(width = 0.5) +
    scale_size_continuous(limits = c(95, 100), range = c(1, 2)) +
    geom_facet(panel = "Assembly length", data = lens_df, geom = geom_col, orientation = 'y') +
    # geom_facet(panel = "Markers", data = mark_df, geom = geom_point) +
    theme_tree2() # theme(legend.position = "bottom")
p <- gheatmap(p, data = mark_df)
ggsave(out_image_file, p, height = ntaxa * 0.1, width = output_width, limitsize = F)
