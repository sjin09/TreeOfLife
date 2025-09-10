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
setwd("/Users/sl666/Manuscript/My Drive/ToL_Nature/SI/SF1-3")

# Load data
# metadata = read.csv("dtol_all_samples.taxonomic_classification.csv")
# txt = readLines("tree.nwk")

# ---- CLI ----
parser <- ArgumentParser(description = "Plot Supplementary Figure 3")
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

# ----------------By Class and Order-------------------
pal1 <- c(
  # Kingdom
  Fungi = "#7F2704", # warm brown
  Viridiplantae = "#115923", # deep green
  # Annelida
  Clitellata = "#3F007D",
  Polychaeta = "#6A51A3",
  # Chordata
  Aves = "#031CA6", # navy blue
  Actinopteri = "#0468BF", # blue
  Mammalia = "#04B2D9", # cyan
  # Mollusca
  Bivalvia = "#B39DD8",
  Gastropoda = "#BCBDDC",
  # Bryozoa
  Gymnolaemata = "#807DBA",
  # Insecta by Order
  Lepidoptera = "#8C510A",
  Coleoptera = "#EC7014",
  Diptera = "#A63603",
  Hemiptera = "#FFD92F",
  Trichoptera = "#FEC44F",
  Plecoptera = "#FE9929",
  Hymenoptera = "#CC4C02",
  # Arthropoda
  Arachnida = "#FEE725"
)

pal2 <- c(
  Fungi = "#8b4513",  # saddlebrown
  Viridiplantae = "#556b2f",  # darkolivegreen
  # Clitellata = "#228b22",  # forestgreen
  # Polychaeta = "#483d8b",  # darkslateblue
  Clitellata = "#ffe4b5",  # forestgreen
  Polychaeta = "#ffb6c1",  # darkslateblue
  Aves = "#ff1493",  # darkcyan
  Actinopteri = "#da70d6",  # steelblue
  # Aves = "#008b8b",  # darkcyan
  # Actinopteri = "#4682b4",  # steelblue
  Mammalia = "#000080",  # navy
  Bivalvia = "#9acd32",  # yellowgreen
  Gastropoda = "#7f007f",  # purple2
  Gymnolaemata = "#8fbc8f",  # darkseagreen
  Lepidoptera = "#ff8c00",  # darkorange
  Coleoptera = "#ffd700",  # gold
  Diptera = "#7fff00",  # chartreuse
  Hemiptera = "#00ff7f",  # springgreen
  Trichoptera = "#00ffff",  # aqua
  Plecoptera = "#f4a460",  # sandybrown
  Hymenoptera = "#0000ff",  # blue
  Arachnida = "#7b68ee"  # purple3
)
#   "#da70d6",  # orchid
#   "#1e90ff",  # dodgerblue
#   "#90ee90",  # lightgreen
#   "#add8e6",  # lightblue
#   "#ff1493",  # deeppink
#   "#7b68ee",  # mediumslateblue
#   "#ffe4b5",  # moccasin
#   "#ffb6c1"   # lightpink

# Build data e annotation
anno <- tibble(label = tr$tip.label) %>%
  mutate(Species = str_trim(str_replace(label, "\\s*\\(.*\\)$", "")))  # drop the (...) part
anno$Kingdom <- sapply(anno$Species,function(x){
  unique(metadata$Kingdom[metadata$Species==x])
})
anno$Phylum <- sapply(anno$Species,function(x){
  unique(metadata$Phylum[metadata$Species==x])
})
anno$Class <- sapply(anno$Species,function(x){
  unique(metadata$Class[metadata$Species==x])
})
anno$Order <- sapply(anno$Species,function(x){
  unique(metadata$Order[metadata$Species==x])
})
anno <- anno %>% mutate(label_expr = get_expr_labels(anno$label))
# palette 1
# anno <- anno %>%
#   mutate(
#     ColorGroup = case_when(
#       Phylum %in% names(pal1) ~ Phylum,                # use phylum if in pal1ette
#       Order %in% names(pal1) ~ Order,                # use phylum if in pal1ette
#       Class %in% names(pal1) ~ Class,                # use phylum if in pal1ette
#       Kingdom == "Viridiplantae" ~ "Viridiplantae",
#       Kingdom == "Fungi" ~ "Fungi",
#       TRUE ~ NA_character_ 
#     )
#   ) %>%
#   mutate(
#     ColorGroup = factor(ColorGroup, levels = names(pal1))  # enforce order
#   ) %>% 
#   select(label, label_expr, ColorGroup)
# palette 2
anno <- anno %>%
  mutate(
    ColorGroup = case_when(
      Phylum %in% names(pal2) ~ Phylum,                # use phylum if in pal2ette
      Order %in% names(pal2) ~ Order,                # use phylum if in pal2ette
      Class %in% names(pal2) ~ Class,                # use phylum if in pal2ette
      Kingdom == "Viridiplantae" ~ "Viridiplantae",
      Kingdom == "Fungi" ~ "Fungi",
      TRUE ~ NA_character_ 
    )
  ) %>%
  mutate(
    ColorGroup = factor(ColorGroup, levels = names(pal2))  # enforce order
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

# dot
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
  # scale_fill_manual(values = pal1, name = "Class and Order")
  scale_fill_manual(values = pal2, name = "Class and Order")
ggsave("SF3_beta.pdf", plot = p, width = 10, height = 10, units = "in") 


