#!/usr/bin/env Rscript

library(plyr)
library(ggplot2)
library(magrittr)
library(data.table)

## options
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
sbs52_mutsig_colours = c("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#F5ABCC")

if (length(args) == 0) {
  cat("No arguments provided.\n")
} else {
    # import
    input = fread(args[1])
    sample = args[2]
    output = args[3]
    # ggplot2
    plot = ggplot(input, aes(x=TRI, y=COUNT, fill=SUB)) +
    geom_bar(stat="identity")  +
    theme_bw(24) + # maybe
    facet_grid(. ~ SUB, scales="free") + 
    scale_fill_manual(values = sbs52_mutsig_colours) +
    theme(legend.title = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "\nTrinucleotide\n", y = "\nCounts\n") +
    theme(
      axis.text.x = element_text(family = "mono", angle = 90, hjust = 1),
      axis.line = element_line(color = "black", linewidth = 1, linetype = "solid"),
    )
    # ggsave
    ggsave(output, plot, width = 26, height=15)
}
