#!/usr/bin/env Rscript

library(lsa)
library(plyr)
library(ggplot2)
library(magrittr)
library(data.table)

## options
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
sbs96_mutsig_colours = c("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#83A603", "#F5ABCC")

if (length(args) == 0) {
  cat("No arguments provided.\n")
} else {
    # import
    input = fread(args[1])
    sample = args[2]
    output = args[3]
    # ggplot2
    plot = ggplot(input, aes(x=TRI, y=PROBABILITY, fill=SUB)) +
    geom_bar(stat="identity")  +
    theme_classic() +
    facet_grid(. ~ SUB, scales = "free") +
    scale_fill_manual(values = sbs96_mutsig_colours) +
    theme(legend.title = element_blank()) +
    # theme(plot.title = element_text(hjust = 0.5)) +
    labs(title = sample, x = "\nTrinucleotide\n", y = "\nProportion of Mutations (%)\n") +
    theme(
      axis.text.x = element_text(family = "mono", size = 15, angle = 90, hjust = 1),
      axis.line = element_line(color = "black", linewidth = 1, linetype = "solid"),
      text = element_text(size=25)
    )
    # ggsave
    ggsave(output, plot, width = 26, height=15)
}
