#!/usr/bin/env Rscript

library(plyr)
library(ggplot2)
library(magrittr)
library(data.table)
library(patchwork)

## options
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
sbs52_mutsig_colours = c("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#F5ABCC")
sbs96_mutsig_colours = c("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#83A603", "#F5ABCC")

if (length(args) == 0) {
  cat("No arguments provided.\n")
} else {
    # import
    sbs52_df = fread(args[1])
    sbs96_df = fread(args[2])
    sample = args[3]
    output = args[4]
    
    # manipulate 
    sbs52_df$TRI <- sapply(strsplit(sbs52_df$TRI, " "), function(x) x[2])
    
    # ggplot2: sbs52 mutational spectrum
    sbs52_plot = ggplot(sbs52_df, aes(x=TRI, y=COUNT, fill=SUB)) +
      geom_bar(stat="identity")  +
      theme_bw(24) + # maybe
      facet_grid(. ~ SUB, scales="free") + 
      scale_fill_manual(values = sbs52_mutsig_colours) +
      theme(legend.title = element_blank()) +
      labs(
        x = "\nTrinucleotide\n", 
        y = "\nCounts\n", 
        title = paste(sample, "germline mutational spectrum")
      ) +
      theme(
        axis.text.x = element_text(family = "mono", angle = 90, hjust = 1),
        axis.line = element_line(color = "black", linewidth = 1, linetype = "solid"),
        legend.position = "none"
      )
    
    # ggplot2: sbs96 mutational spectrum
    sbs96_plot = ggplot(sbs96_df, aes(x=tri, y=normcounts, fill=sub)) +
      geom_bar(stat="identity")  +
      theme_bw(24) + 
      scale_fill_manual(values = sbs96_mutsig_colours) +
      theme(legend.title = element_blank()) +
      facet_grid(. ~ sub, scales="free") + 
      # theme(plot.title = element_text(hjust = 0.5)) +
      labs(
        x = "\nTrinucleotide\n", 
        y = "\nCounts\n",
        title = paste(sample, "somatic mutational spectrum")
      ) +
      theme(
        axis.text.x = element_text(family = "mono", angle = 90, hjust = 1),
        axis.line = element_line(color = "black", linewidth = 1, linetype = "solid"),
        legend.position = "none"
      )
   
    # combine plots  
    plot = sbs52_plot / sbs96_plot
    
    # ggsave
    ggsave(output, plot, width = 24, height = 28)
}


 