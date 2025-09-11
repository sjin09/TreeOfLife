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
parser <- ArgumentParser(description = "Plot Supplementary Figure 2")
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
tr$tip.label <- str_replace_all(tr$tip.label, "_", " ")
tr$tip.label <- gsub("\\^(.*?)\\^", "(\\1)", tr$tip.label)

# Define palette
pal <- c(
  Viridiplantae = "#004300",
  Fungi = "#8C564B",
  Arthropoda = "#B80080",
  Chordata = "#30A2DA",
  Mollusca = "#000097",
  Annelida = "#73168C",
  Bryozoa = "#0B7A75",
  Cnidaria = "#A9D2D5"
)

# build annotation
anno <- tibble(label = tr$tip.label) %>%
  mutate(Species = str_trim(str_replace(label, "\\s*\\(.*\\)$", "")))  # drop the (...) part
anno$Kingdom <- sapply(anno$Species,function(x){
  unique(metadata$Kingdom[metadata$Species==x])
})
anno$Phylum <- sapply(anno$Species,function(x){
  unique(metadata$Phylum[metadata$Species==x])
})
anno <- anno %>% mutate(label_expr = get_expr_labels(anno$label))
anno <- anno %>%
  mutate(
    ColorGroup = case_when(
      Phylum %in% names(pal)        ~ Phylum,                # use phylum if in palette
      Kingdom == "Viridiplantae"    ~ "Viridiplantae",
      Kingdom == "Fungi"            ~ "Fungi",
      TRUE ~ NA_character_ 
    )
  ) %>%
  mutate(
    ColorGroup = factor(ColorGroup, levels = names(pal))  # enforce order
  ) %>% 
  select(label, label_expr, ColorGroup)

# build data frame
grp_list <- split(anno$label, anno$ColorGroup)
hilight_df <- lapply(names(grp_list), function(g) {
  tips <- grp_list[[g]]
  idx  <- which(tr$tip.label %in% tips)
  if (length(idx) >= 2) {
    node_id <- ape::getMRCA(tr, idx)
    if (!is.na(node_id)) data.frame(node = node_id, ColorGroup = g)
  }
}) %>% bind_rows()

# assign expression label
tr$tip.label = anno$label_expr

# plot
p <- ggtree(tr, layout = "circular") +
  geom_hilight(
    data = hilight_df,
    aes(node = node, fill = ColorGroup),
    alpha = 0.6
  ) +
  geom_tiplab(
    aes(label = label),
    parse = TRUE,
    family = "Helvetica", 
    size = 0.8
  ) + 
  theme(
    legend.position = "right",
    text = element_text(family = "Helvetica")
  ) +
  scale_fill_manual(values = pal, name = "Phylum")
ggsave("SF2.pdf", plot = p, width = 10, height = 10, units = "in") 


