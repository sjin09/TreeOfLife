#!/usr/bin/env Rscript

library(lsa)
library(plyr)
library(ggplot2)
library(magrittr)
library(data.table)

## options
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
sbs96_mutsig_colours = c("#98D7EC","#212121","#FF003A","#A6A6A6","#83A603","#F5ABCC")

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
    facet_wrap(~SUB, nrow=1, scales="free_x") + 
    scale_fill_manual(values = sbs96_mutsig_colours) +
    theme(
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(color = "black", linewidth = 1, linetype = "solid"),
        legend.position = "none",
        text = element_text(size=25)
    ) +
    labs(title = sample, x = "", y = "")
    # ggsave
    ggsave(output, plot, width = 26, height=15)
}
