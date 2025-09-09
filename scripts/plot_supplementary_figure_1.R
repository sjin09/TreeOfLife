#!/usr/bin/env Rscript

# Load packages
library(ape)
library(argparse)
library(ggtree)
library(ggtreeExtra)
library(tidyverse)
library(stringr)
library(dplyr)

## setwd
# setwd("/Users/sl666/Manuscript/My Drive/ToL_Nature/SI/SF1-3")

# Load data
# metadata = read.csv("dtol_all_samples.taxonomic_classification.csv")
# txt = readLines("tree.nwk")

# ---- CLI ----
parser <- ArgumentParser(description = "Plot Supplementary Figure 1")
parser$add_argument(
  "--newick", 
  required = TRUE,
  help = "Path to Newick tree"
)
parser$add_argument(
  "--taxonomic-classification", 
  required = TRUE,
  help = "Path to dtol_all_samples.taxonomic_classification.csv"
)
args <- parser$parse_args()

# get parsed labels
get_expr_labels <- function(labels){
  esc <- function(s) gsub("'", "\\\\'", s, perl = TRUE)  # escape single quotes

  has_paren = grepl("\\(", labels)
  
  # get species name
  species   <- sub("\\s*\\(.*$", "", labels) 
  
  # get sample name
  paren     <- sub("^[^\\(]*", "", labels)
  
  species_e <- esc(species)
  paren_e   <- esc(paren)
  
  ifelse(
    has_paren,
    paste0("italic('", species_e, "')~'", paren_e, "'"),  # italic species + plain "(...)"
    paste0("italic('", species_e, "')")                   # all italic when no "(...)"
  )
}

# Load data
metadata = read.csv(args$taxonomic_classification)
txt = readLines(args$newick)

# Manipulate data
txt = gsub("_", "^" , txt, fixed = TRUE)
txt2 <- gsub(" +", "_", txt)   # turn spaces inside labels into underscores
tr = read.tree(text=txt2)
org_labels =  tr$tip.label
org_labels <- str_replace_all(org_labels, "_", " ")
org_labels <- gsub("\\^(.*?)\\^", "(\\1)", org_labels)

# build annotation
anno <- tibble(label = org_labels) %>%
  mutate(Species = str_trim(str_replace(label, "\\s*\\(.*\\)$", "")))  # drop the (...) part
anno$Kingdom <- sapply(anno$Species,function(x){
  unique(metadata$Kingdom[metadata$Species==x])
})
anno <- anno %>% mutate(label_expr = get_expr_labels(org_labels))

# assign expression label
tr$tip.label = anno$label_expr

# Palette
pal <- c(
  Fungi = "#7F2704",
  Metazoa = "#001F54",
  Viridiplantae = "#115923"
)

# annotate by shade
king_groups <- split(anno$label, anno$Kingdom)
# king_groups <- split(anno$label_expr, anno$Kingdom)
tr_g <- groupOTU(tr, king_groups, group_name = "kingdom")
hilight_df <- lapply(names(king_groups), function(k) {
  idx <- which(anno$label %in% king_groups[[k]])
  if (length(idx) >= 2) {
    node_id <- getMRCA(tr, idx)
    if (!is.na(node_id)) {
      return(data.frame(node = node_id, Kingdom = k))
    }
  }
  NULL
}) %>% bind_rows()
# 
p <- ggtree(tr, layout = "circular")
p <- p + geom_hilight(data = hilight_df,
                      aes(node = node, fill = Kingdom),
                      alpha = 0.6) +
  scale_fill_manual(values = pal, name = "Kingdom") +
  geom_tiplab(
    aes(label = label),
    parse = TRUE,
    family = "Helvetica", 
    size = 0.8
  ) +
  theme(
    legend.position = "right",
    text = element_text(family = "Helvetica")
    )
ggsave("SF1.pdf", plot = p, width = 10, height = 10, units = "in")


# annotate by rectangle
# p <- ggtree(tr, layout = "circular", size = 0.5) +
#   geom_fruit(
#     data = anno,
#     geom = geom_tile,
#     mapping = aes(y = label_expr, fill = Kingdom),
#     width = 0.2, 
#     offset = 0.02, 
#     color = NA
#   ) + 
#   geom_tiplab(
#     aes(label = label),
#     parse = TRUE,
#     offset = 0.37,
#     family = "Helvetica",
#     size = 0.8
#   ) +
#   theme(legend.position = "right",
#         text = element_text(family = "Helvetica"))+
#   scale_fill_manual(values = pal, name = "Kingdom")
# ggsave("SF1_rectangle.pdf", plot = p, width = 10, height = 10, units = "in") 
# 
# # annotate by circle
# p <- ggtree(tr, layout = "circular", size = 0.5) +
#   geom_fruit(
#     data    = anno,
#     geom    = geom_point, 
#     mapping = aes(y = label_expr, color = Kingdom),  # map to color, not fill
#     size    = 0.3,                  # adjust circle size
#     offset  = 0.02
#   ) + 
#   geom_tiplab(
#     aes(label = label),
#     parse = TRUE,
#     offset = 0.3,    # smaller offset since no wide tiles
#     family = "Helvetica",
#     size = 0.8
#   ) +
#   theme(legend.position = "right",
#         text = element_text(family = "Helvetica")) +
#   scale_color_manual(values = pal, name = "Kingdom")
# ggsave("SF1_point.pdf", plot = p, width = 10, height = 10, units = "in") 
