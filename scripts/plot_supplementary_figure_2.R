#!/usr/bin/env Rscript

# Load packages
library(ape)
library(argparse)
library(ggtree)
library(ggtreeExtra)
library(tidyverse)
library(stringr)
library(dplyr)

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
metadata = read.csv("dtol_all_samples.taxonomic_classification.csv")
txt = readLines("tree.nwk")

# Manipulate data
txt = gsub("_", "^" , txt, fixed = TRUE)
txt2 <- gsub(" +", "_", txt)   # turn spaces inside labels into underscores
tr = read.tree(text=txt2)
org_labels =  tr$tip.label
org_labels <- str_replace_all(org_labels, "_", " ")
org_labels <- gsub("\\^(.*?)\\^", "(\\1)", org_labels)

# Define palette

# ----------------By Phylum-------------------
pal <- c(
  Arthropoda = "#CE6C47",
  Chordata = "#EABE63",
  Mollusca = "#C04BF2",
  Annelida = "#73168C",
  Bryozoa = "#0B7A75",
  Cnidaria = "#A9D2D5",
  Viridiplantae = "#115923",
  Fungi = "#A67458",
  `Other Metazoa` = "#CECBD0"
)

# build annotation
anno <- tibble(label = org_labels) %>%
  mutate(Species = str_trim(str_replace(label, "\\s*\\(.*\\)$", "")))  # drop the (...) part
anno$Kingdom <- sapply(anno$Species,function(x){
  unique(metadata$Kingdom[metadata$Species==x])
})
anno$Phylum <- sapply(anno$Species,function(x){
  unique(metadata$Phylum[metadata$Species==x])
})
anno <- anno %>% mutate(label_expr = get_expr_labels(org_labels))
anno <- anno %>%
  mutate(
    ColorGroup = case_when(
      Phylum %in% names(pal)        ~ Phylum,                # use phylum if in palette
      Kingdom == "Viridiplantae"    ~ "Viridiplantae",
      Kingdom == "Fungi"            ~ "Fungi",
      TRUE                          ~ "Other Metazoa"          # anything else -> uncoloured
    )
  ) %>%
  mutate(
    ColorGroup = factor(ColorGroup, levels = names(pal))  # enforce order
  )

# assign expression label
tr$tip.label = anno$label_expr

# plot
p <- ggtree(tr, layout = "circular", size = 0.5) +
  geom_fruit(
    data    = anno,
    geom    = geom_point,
    mapping = aes(y = label_expr, color = ColorGroup),
    size    = 0.3,
    offset  = 0.02
  ) +
  geom_tiplab(
    aes(label = label),
    parse = TRUE,
    offset = 0.3,
    family = "Helvetica", 
    size = 0.8
  ) +
  theme(
    legend.position = "right",
    text = element_text(family = "Helvetica")
  ) +
  scale_color_manual(values = pal, name = "Taxonomic rank")
ggsave("SF2.pdf", plot = p, width = 10, height = 10, units = "in") 


