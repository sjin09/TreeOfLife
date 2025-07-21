
#!/usr/bin/env Rscript

library(plyr)
library(ggplot2)
library(magrittr)
library(data.table)
library(patchwork)

## options
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
sbs96_mutsig_colours = c("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#83A603", "#F5ABCC")

if (length(args) == 0) {
  cat("No arguments provided.\n")
} else {
    # import
    somatic_sbs96_df = fread(args[1])
    somatic_normalised_sbs96_df = fread(args[2])
    species_name = args[3]
    species_name = paste(unlist(strsplit(species_name, "_")), collapse = " ")
    sample_name = args[4]
    sample_name = sprintf(" (%s)", sample_name)
    output = args[5]

    # ggplot2
    sbs96_plot = ggplot(somatic_sbs96_df, aes(x=tri, y=normcounts, fill=sub)) +
        geom_bar(stat="identity")  +
        theme_bw(24) +
        facet_grid(. ~ sub, scales="free") +
        scale_fill_manual(values = sbs96_mutsig_colours) +
        labs(
             x = "\nTrinucleotide\n",
             y = "\nCounts\n",
             title = bquote(italic(.(species_name)) * .(sample_name) * " somatic mutational spectrum")
        ) +
        theme(
          axis.text.x = element_text(family = "mono", angle = 90, hjust = 1),
          axis.line = element_line(color = "black", linewidth = 1, linetype = "solid"),
          legend.position = "none",
          legend.title = element_blank()
        )

    # ggplot2
    normalised_sbs96_plot = ggplot(somatic_normalised_sbs96_df, aes(x=TRI, y=WEIGHTED_COUNT, fill=SUB)) +
        geom_bar(stat="identity")  +
        theme_bw(24) +
        facet_grid(. ~ SUB, scales="free") +
        scale_fill_manual(values = sbs96_mutsig_colours) +
        labs(
             x = "\nTrinucleotide\n", 
             y = "\nNormalised Counts\n",
             title = bquote(italic(.(species_name)) * .(sample_name) * " normalised somatic mutational spectrum")
        ) +
        theme(
          axis.text.x = element_text(family = "mono", angle = 90, hjust = 1),
          axis.line = element_line(color = "black", linewidth = 1, linetype = "solid"),
          legend.position = "none",
          legend.title = element_blank()
        )

    # combine plots
    plot = sbs96_plot / normalised_sbs96_plot

    # ggsave
    ggsave(output, plot, width = 24, height = 28)
}
