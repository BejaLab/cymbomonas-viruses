suppressPackageStartupMessages(library(dplyr))
library(XML)
library(ggplot2)
library(ggseqlogo)
library(stringr)
library(tidyr)
suppressPackageStartupMessages(library(cowplot))

input_file <- unlist(snakemake@input)
output_file <- unlist(snakemake@output)
with(snakemake@params, {
	motif_res <<- data.frame(motif_name = motif_res, regex = motif_res)
	meme_name <<- meme_name
})
data <- xmlToList(xmlParse(input_file))
n_seq <- as.numeric(data$training_set$.attrs["primary_count"])

motif.data <- data$motifs %>%
	`[`(names(.) == "motif") %>%
	setNames(sapply(., function(x) x$.attrs["id"]))

motifs <- motif.data %>%
	lapply(function(x) as.data.frame(t(x$.attrs))) %>%
	bind_rows %>%
	mutate_at(vars(-id, -name, -alt), as.numeric) %>%
	mutate(sites_r = sites / n_seq)
motifs.filtered <- left_join(motifs, motif_res, by = character()) %>%
	filter(str_match(name, regex) > 0) %>%
	distinct(motif_name, .keep_all = T)

site.attrs <- function(site) {
	attrs <- as.data.frame(t(site$.attrs))
	attrs$left_flank  <- site$left_flank
	attrs$site        <- site$site[names(site$site) == "letter_ref"] %>% paste(collapse = "")
	attrs$right_flank <- site$right_flank
	return(attrs)
}
fetch.sites <- function(sites) {
	sites <- sites$contributing_sites
	sites[names(sites) == "contributing_site"] %>% lapply(site.attrs) %>% bind_rows
}

sites <- lapply(motif.data, fetch.sites) %>%
	bind_rows(.id = "id") %>%
	filter(id %in% motifs.filtered$id) %>%
	mutate_at(vars(position, pvalue), as.numeric)

plot.motif <- function(motif, .y) {
	motif.sites <- sites %>% filter(id == .y$id)
	label <- with(motif, paste(
		gsub("_", " ", meme_name),
		name,
		paste("Sites =", sites, paste0("(", round(sites_r * 100, 1), "%)")),
		# paste("IC =", round(ic, 3)),
		paste("E-value =", e_value),
		sep = "\n"
	))
	g1 <- ggplot() + annotate("text", label = label, x = 1, y = 1) +
		theme_bw() + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	g2 <- motif.sites %>% pull(site) %>% ggseqlogo(seq_type = "dna") +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	g3 <- ggplot(motif.sites, aes(x = position)) + geom_histogram(binwidth=10) + xlim(0, 150) +
		theme_bw() +
		theme(axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
	g <- plot_grid(g1, g2, g3, nrow = 1, align = 'h')
	return(g)
}

if (nrow(motifs.filtered) == 0) {
	g <- ggplot() + theme_void()
	ggsave(output_file, g, width = 1, height = 1)
} else {
	g.list <- motifs.filtered %>%
		group_by(id) %>%
		group_map(plot.motif)
	ggsave(output_file, plot_grid(plotlist = g.list, ncol = 1), width = 7, height = length(g.list) * 1.4)
}
