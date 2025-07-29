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

if (length(args) == 0) {
  cat("No arguments provided.\n")
} else {
    # import
    germline_signature_df = fread(args[1])
    germline_attribution_df = fread(args[2])
    sig_name = args[3]
    output = args[4]
   
    # manipulate signature
    germline_signature_df$TRI <- sapply(strsplit(germline_signature_df$TRI, " "), function(x) x[2])
     
    # manipulate attributions
    germline_attribution_df = germline_attribution_df[order(ATTRIBUTION, decreasing=TRUE)]
    germline_attribution_df[, x_label_expr := mapply(function(id, sp) {
      paste0("paste('", id, "', ' (', italic('", sp, "'), ')')")
    }, SAMPLE, SPECIES)]
    germline_attribution_df[, x_label_expr := factor(
      x_label_expr, 
      levels = x_label_expr[order(ATTRIBUTION, decreasing = TRUE)]
    )]
     
    # ggplot2
    sig_plot = ggplot(germline_signature_df, aes(x=TRI, y=PROBABILITY, fill=SUB)) +
        geom_bar(stat="identity")  +
        theme_bw(24) +
        facet_grid(. ~ SUB, scales="free") +
        scale_fill_manual(values = sbs52_mutsig_colours) +
        labs(
             x = "\nTrinucleotide\n",
             y = "\nProportion of Mutations (%)\n",
             title = bquote(.(sig_name) * " mutational signature")
        ) +
        theme(
          axis.text.x = element_text(family = "mono", angle = 90, hjust = 1),
          axis.line = element_line(color = "black", linewidth = 1, linetype = "solid"),
          legend.position = "none",
          legend.title = element_blank()
        )

    # ggplot2
    sig_attribution_plot = ggplot(germline_attribution_df[1:100,], aes(x=x_label_expr, y=ATTRIBUTION)) +
        geom_bar(stat="identity")  +
        theme_bw(24) +
        labs(
             x = "\nSample\n", 
             y = "\nMutational Signature Attribution (%)\n",
        ) + 
        theme(
          axis.text.x = element_text(size = 16, angle = 90, hjust = 1),
          axis.line = element_line(color = "black", linewidth = 1, linetype = "solid"),
          legend.position = "none",
          legend.title = element_blank()
        ) + 
        scale_x_discrete(labels = function(x) sapply(x, function(lbl) parse(text = lbl)))

    # combine plots
    plot = sig_plot / sig_attribution_plot 

    # ggsave
    ggsave(output, plot, width = 24, height = 28)
}
