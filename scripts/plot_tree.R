library(ape)
library(ggtree)
library(tidyverse)
library(ggtreeExtra)
library(stringr)
library(dplyr)


metadata <- read.csv("../data/dtol_all_samples.taxonomic_classification.csv")

txt  <- readLines("../data/ToL_new.newick")
txt <- gsub('_', '^', txt, fixed = TRUE)
txt2 <- gsub(" +", "_", txt)   # turn spaces inside labels into underscores
tr   <- read.tree(text = txt2)
tr$tip.label <- str_replace_all(tr$tip.label, "_", " ")
tr$tip.label <- gsub("\\^(.*?)\\^", "(\\1)", tr$tip.label)

# -------By kingdom---------
# Palette
pal <- c(
  Fungi = "#A67458",
  Metazoa = "#001F54",
  Viridiplantae = "#115923"
)


anno <- tibble(label = tr$tip.label) %>%
  mutate(Species = str_trim(str_replace(label, "\\s*\\(.*\\)$", "")))  # drop the (...) part
anno$Kingdom <- sapply(anno$Species,function(x){
  unique(metadata$Kingdom[metadata$Species==x])
})
  

# colour by shade
# king_groups <- split(anno$label, anno$Kingdom)
# tr_g <- groupOTU(tr, king_groups, group_name = "kingdom")
# hilight_df <- lapply(names(king_groups), function(k) {
#   idx <- which(tr$tip.label %in% king_groups[[k]])
#   if (length(idx) >= 2) {
#     node_id <- getMRCA(tr, idx)
#     if (!is.na(node_id)) {
#       return(data.frame(node = node_id, Kingdom = k))
#     }
#   }
#   NULL
# }) %>% bind_rows()
# 
# p <- ggtree(tr, layout = "circular")
# p <- p + geom_hilight(data = hilight_df,
#                       aes(node = node, fill = Kingdom),
#                       alpha = 0.6) +
#   scale_fill_manual(values = pal, name = "Kingdom") +
#   geom_tiplab(aes(label = label),
#               fontface = "italic", family = "Helvetica", size = 1) +
#   theme(legend.position = "right",
#         text = element_text(family = "Helvetica"))
# 
# p

# colour by ring
p<-ggtree(tr, layout = "circular", size = 0.5) +
  geom_fruit(data = anno,
  geom = geom_tile,
  mapping = aes(y = label, fill = Kingdom),
  width = 0.2, offset = 0.02, color = NA) +
  geom_tiplab(aes(label = label),
              offset = 0.37,
              fontface = "italic",
              family = "Helvetica",
              size = 1) +
  theme(legend.position = "right",
        text = element_text(family = "Helvetica"))+
  scale_fill_manual(values = pal, name = "Kingdom")
p
ggsave("~/Documents/Postdoc/ToL/figs/ToL_kingdom_ring.pdf", plot = p,
       width = 10, height = 10, units = "in") 


# colour by dot
p<-ggtree(tr, layout = "circular", size = 0.25) +
  geom_fruit(
    data    = anno,
    geom    = geom_point,        # << use points instead of tiles
    mapping = aes(y = label, color = Kingdom),  # map to color, not fill
    size    = 0.2,                  # adjust circle size
    offset  = 0.02
  ) +
  geom_tiplab(aes(label = label),
              offset   = 0.3,    # smaller offset since no wide tiles
              fontface = "italic",
              family   = "Helvetica",
              size     = 1) +
  theme(legend.position = "right",
        text = element_text(family = "Helvetica")) +
  scale_color_manual(values = pal, name = "Kingdom")
p

ggsave("~/Documents/Postdoc/ToL/figs/ToL_kingdom_ring.pdf", plot = p,
       width = 10, height = 10, units = "in") 

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
  



anno <- tibble(label = tr$tip.label) %>%
  mutate(Species = str_trim(str_replace(label, "\\s*\\(.*\\)$", "")))  # drop the (...) part
anno$Kingdom <- sapply(anno$Species,function(x){
  unique(metadata$Kingdom[metadata$Species==x])
})
anno$Phylum <- sapply(anno$Species,function(x){
  unique(metadata$Phylum[metadata$Species==x])
})


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

# anno <- metadata %>%   
#   filter(Species %in% tr$tip.label) %>%
#   rename(label = Species)%>%
#   mutate(ColorGroup = ifelse(Phylum %in% names(pal), Phylum, Kingdom)) %>%
#   mutate(ColorGroup = ifelse(ColorGroup == "Viridiplantae", "Other Viridiplantae", ColorGroup)) %>%
#   mutate(ColorGroup = ifelse(ColorGroup == "Metazoa", "Other Metazoa", ColorGroup)) %>%
#   mutate(ColorGroup = ifelse(ColorGroup == "Fungi", "Other Fungi", ColorGroup)) %>%
#   select(label, ColorGroup)

# ring
p<-ggtree(tr, layout = "circular", size = 0.25) +
  geom_fruit(data = anno,
             geom = geom_tile,
             mapping = aes(y = label, fill = ColorGroup),
             width = 0.2, offset = 0.02, color = NA) +
  geom_tiplab(aes(label = label),
              offset = 0.37,
              fontface = "italic",
              family = "Helvetica",
              size = 1) +
  theme(legend.position = "right",
        text = element_text(family = "Helvetica"))+
  scale_fill_manual(values = pal, , name = "Taxonomic groups")
p
ggsave("~/Documents/Postdoc/ToL/figs/ToL_phylum_ring.pdf", plot = p,
       width = 10, height = 10, units = "in") 

# dot
p <- ggtree(tr, layout = "circular", size = 0.25) +
  geom_fruit(
    data    = anno,
    geom    = geom_point,
    mapping = aes(y = label, color = ColorGroup),
    size    = 2,
    offset  = 0.02
  ) +
  geom_tiplab(aes(label = label),
              offset = 0.3, fontface = "italic",
              family = "Helvetica", size = 1) +
  theme(legend.position = "right",
        text = element_text(family = "Helvetica")) +
  scale_color_manual(values = pal, name = "Taxonomic groups")
p
ggsave("~/Documents/Postdoc/ToL/figs/ToL_phylum_circular.pdf", plot = p,
       width = 10, height = 10, units = "in") 

# ----------------By Class-------------------
pal <- c(
  Aves = "#031CA6",
  Actinopteri = "#0468BF",
  Lepidoptera = "#F25C05",
  Coleoptera = "#324001",
  Diptera = "#B2BF4B",
  Hymenoptera = "#F2E205",
  Mammalia = "#04B2D9",
  Annelida = "#73168C",
  Mollusca= "#C04BF2",
  Fungi = "#A67458",
  Viridiplantae = "#115923",
  `Other Metazoa` = "#CECBD0"
)

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

anno <- anno %>%
  mutate(
    ColorGroup = case_when(
      Phylum %in% names(pal)        ~ Phylum,                # use phylum if in palette
      Order %in% names(pal)        ~ Order,                # use phylum if in palette
      Kingdom == "Viridiplantae"    ~ "Viridiplantae",
      Kingdom == "Fungi"            ~ "Fungi",
      TRUE                          ~ "Other Metazoa"          # anything else -> uncoloured
    )
  ) %>%
  mutate(
    ColorGroup = factor(ColorGroup, levels = names(pal))  # enforce order
  )

# anno <- metadata %>%   
#   filter(Species %in% tr$tip.label) %>%
#   rename(label = Species)%>%
#   mutate(ColorGroup = ifelse(Class == "Insecta", Order, Class)) %>% 
#   mutate(ColorGroup = ifelse(Kingdom == "Viridiplantae", Kingdom, ColorGroup)) %>% 
#   mutate(ColorGroup = ifelse(Kingdom == "Fungi", Kingdom, ColorGroup)) %>% 
#   mutate(ColorGroup = ifelse(Phylum == "Mollusca", Phylum, ColorGroup)) %>% 
#   mutate(ColorGroup = ifelse(Phylum == "Annelida", Phylum, ColorGroup)) %>% 
#   mutate(ColorGroup = ifelse(ColorGroup %in% names(pal), ColorGroup, "Other Metazoa")) %>% 
#   select(label, ColorGroup)

# ring
p<-ggtree(tr, layout = "circular", size = 0.25) +
  geom_fruit(data = anno,
             geom = geom_tile,
             mapping = aes(y = label, fill = ColorGroup),
             width = 0.2, offset = 0.02, color = NA) +
  geom_tiplab(aes(label = label),
              offset = 0.37,
              fontface = "italic",
              family = "Helvetica",
              size = 1) +
  theme(legend.position = "right",
        text = element_text(family = "Helvetica"))+
  scale_fill_manual(values = pal, , name = "Taxonomic groups")
p
ggsave("~/Documents/Postdoc/ToL/figs/ToL_class_ring.pdf", plot = p,
       width = 10, height = 10, units = "in") 
# dot
p <- ggtree(tr, layout = "circular", size = 0.25) +
  geom_fruit(
    data    = anno,
    geom    = geom_point,
    mapping = aes(y = label, color = ColorGroup),
    size    = 2,
    offset  = 0.02
  ) +
  geom_tiplab(aes(label = label),
              offset   = 0.3,
              fontface = "italic",
              family   = "Helvetica",
              size     = 1) +
  theme(legend.position = "right",
        text = element_text(family = "Helvetica")) +
  scale_color_manual(values = pal, name = "Taxonomic group")
p

ggsave("~/Documents/Postdoc/ToL/figs/ToL_class_circular.pdf", plot = p,
       width = 10, height = 10, units = "in") 


