#!/usr/bin/env Rscript

# Load packages
library(ape)
library(argparse)
library(ggtree)
library(ggtreeExtra)
library(tidyverse)
library(stringr)
library(dplyr)
library(Polychrome)

## setwd
setwd("/Users/sl666/Manuscript/My Drive/ToL_Nature/SI/SF1-3")

# Load data
metadata = read.csv("dtol_all_samples.taxonomic_classification.csv")
txt = readLines("tree.nwk")

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

pal <- c(
  # Yellows / Oranges
  "#FFDC00", "#E5AE38", "#FF8002", 
  
  # Blues
  "#30A2DA", "#2B55B9", "#1F77B4", "#4D6EFF", "#6126FF",
  "#000097", "#0100F8", "#00457B", "#6E7CBB", "#888FFF",
  
  # Cyans / Teals / Aquas
  "#04E3C8", "#17BECF", "#01BE8A", "#00BF00", "#00E83B", 
  "#83D371", "#91FF00", "#B4FF92", 
  
  # Purples / Violets
  "#9467BD", "#7D00A0", "#9400F5", "#895DFF", "#520066",
  "#3A0183", "#A737AE", "#7E7CBB", "#AD8BB1", "#CE85FF",
 
  # Pinks / Magentas (strong, not pastel)
  "#FF55FF", "#DD00FF", "#FF1F83", "#F500C6", "#E377C2",
  "#B80080", "#84206F", "#FF798F", "#A56089", "#D796AB",
   
  # Blue-greens / Blue-purples (bridges)
  "#95D3FF", "#7CB2FF"
)

# taxanomic ranks
taxaonomic_ranks <- c(
  "Viridiplantae",
  "Fungi",
  "Lepidoptera",
  "Diptera",
  "Hymenoptera",
  "Coleoptera",
  "Hemiptera",
  "Plecoptera",
  "Trichoptera",
  "Arachnida",
  "Actinopteri",
  "Mammalia",
  "Aves",
  "Clitellata",
  "Polychaeta",
  "Bivalvia",
  "Gastropoda",
  "Gymnolaemata"
)

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
set.seed(42)
# rank_colours = c()
# rank_colours = c("#004300", "#8C564B")
# rank_colours = c(rank_colours, sample(pal, length(taxaonomic_ranks) - 2, replace=FALSE))
# names(rank_colours) = taxaonomic_ranks
rank_colours <- c(
  Viridiplantae = "#004300",
  Fungi        = "#8C564B",
  Lepidoptera  = "#B80080",
  Diptera      = "#FFDC00",
  Hymenoptera  = "#895DFF",
  Coleoptera   = "#0100F8",
  Hemiptera    = "#E377C2",
  Plecoptera   = "#00E83B",
  Trichoptera  = "#9400F5",
  Arachnida    = "#04E3C8",
  Actinopteri  = "#30A2DA",
  Mammalia     = "#AB7200",
  Aves         = "#520066",
  Clitellata   = "#84206F",
  Polychaeta   = "#17BECF",
  Bivalvia     = "#FF8002",
  Gastropoda   = "#000097",
  Gymnolaemata = "#F500C6"
)

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
anno <- anno %>%
  mutate(
    ColorGroup = case_when(
      Phylum %in% names(rank_colours) ~ Phylum,                # use phylum if in rank_coloursette
      Order %in% names(rank_colours) ~ Order,                # use phylum if in rank_coloursette
      Class %in% names(rank_colours) ~ Class,                # use phylum if in rank_coloursette
      Kingdom == "Viridiplantae" ~ "Viridiplantae",
      Kingdom == "Fungi" ~ "Fungi",
      TRUE ~ NA_character_ 
    )
  ) %>%
  mutate(
    ColorGroup = factor(ColorGroup, levels = names(rank_colours))  # enforce order
  ) %>% 
  select(label, label_expr, ColorGroup)
# build data frame
tr$tip.label = anno$label
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
  theme(
    legend.position = "right",
    text = element_text(family = "Helvetica")
  )
ggsave("EXF3_alpha.pdf", plot = p, width = 10, height = 10, units = "in")

# plot
p <- ggtree(tr, layout = "circular") +
  geom_hilight(
    data = hilight_df,
    aes(node = node, fill = ColorGroup),
    alpha = 0.6
  ) +
  theme(
    legend.position = "right",
    text = element_text(family = "Helvetica")
  ) +
  scale_fill_manual(values = rank_colours, name = "Taxanomic rank")
ggsave("EXF3_beta.pdf", plot = p, width = 10, height = 10, units = "in")

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
  # scale_fill_manual(values = pal, name = "Class and Order")
  scale_fill_manual(values = rank_colours, name = "Taxanomic rank")
ggsave("SF3.pdf", plot = p, width = 10, height = 10, units = "in")


