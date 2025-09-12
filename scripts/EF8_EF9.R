library(ape)
library(ggtree)
library(tidyverse)
library(ggtreeExtra)
library(ggnewscale)
library(stringr)
library(dplyr)


metadata <- read.csv("../data/dtol_all_samples.taxonomic_classification.csv")

somatic_signatures <- read.csv("../data/somatic_mutational_signature_attributions.x0_excluded.rtol_filtered.csv",check.names = F)
colnames(somatic_signatures) = paste0('sToL',colnames(somatic_signatures) )
somatic_signatures[,1] <- ifelse(grepl("\\.", somatic_signatures[,1]),
                                 sub(".*\\.", "", somatic_signatures[,1]),
                                 somatic_signatures[,1])
somatic_signatures$label=NA
somatic_signatures$label=sapply(somatic_signatures[,1], function(x){
  Species = metadata$Species[metadata$Sample==x]
  if (length(metadata$Species[metadata$Species==Species]) == 1){
    label = Species
  }
  else{
    label = paste0(Species," (", x, ")")
  } 
  label
})
somatic_signatures$label[somatic_signatures$label=="Hemaris fuciformis"]="Hemaris fuciformis (iHemFuc2)"


txt  <- readLines("../data/ToL_new.newick")
txt <- gsub('_', '^', txt, fixed = TRUE)
txt2 <- gsub(" +", "_", txt)   # turn spaces inside labels into underscores
tr   <- read.tree(text = txt2)
tr$tip.label <- str_replace_all(tr$tip.label, "_", " ")
tr$tip.label <- gsub("\\^(.*?)\\^", "(\\1)", tr$tip.label)
# org_labels =  tr$tip.label

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

# tr$tip.label = get_expr_labels(org_labels)
# ------------------Extended Figure 8--------------
pal <- c(
  Fungi = "#8C564B",
  Metazoa = "#001F54",
  Viridiplantae = "#004300"
)


anno <- tibble(label = org_labels) %>%
  mutate(Species = str_trim(str_replace(label, "\\s*\\(.*\\)$", "")))  # drop the (...) part
anno$Kingdom <- sapply(anno$Species,function(x){
  unique(metadata$Kingdom[metadata$Species==x])
})

grp_list <- split(anno$label, anno$Kingdom)

hilight_df <- lapply(names(grp_list), function(g) {
  tips <- grp_list[[g]]
  idx  <- which(tr$tip.label %in% tips)
  if (length(idx) >= 2) {
    node_id <- ape::getMRCA(tr, idx)
    if (!is.na(node_id)) data.frame(node = node_id, ColorGroup = g)
  }
}) %>% bind_rows()

p <- ggtree(tr, layout = "circular", size = 0.25,color = "grey60") %<+% anno +
  geom_hilight(data = hilight_df,
               aes(node = node, fill = ColorGroup),
               alpha = 0.5) +
  # geom_tiplab(aes(label = label),
  #             offset = 0.3,
  #             fontface = "italic",
  #             family = "Helvetica",
  #             size = 1) +
  scale_fill_manual(values = pal, name = "Taxonomic group", 
                    guide = guide_legend(override.aes = list(alpha = 1))) #+
# theme(legend.position = "right",
#       text = element_text(family = "Helvetica"))

p

ggsave("../figs/ToL_EF8_noheatmap.pdf", plot = p,
       width = 7, height = 7, units = "in") 



hm <- somatic_signatures %>%
  select(label, sToL1) %>%             
  distinct(label, .keep_all = TRUE) %>% 
  filter(label %in% org_labels) %>%   
  column_to_rownames("label")


# --- add a NEW fill scale, then the heatmap ring ---
# custom palettes for each signature
pal_1  <- c("white", "#E41A1B", "#A51212")


# layout parameters
base_offset <- 0.05
band_width  <- 0.05
gap         <- 0.3

grid_col  <- "grey85" 


p <- p + guides(fill = "none", colour = "none")
p1 <- p + new_scale_fill()
p1 <- gheatmap(p1, hm[,"sToL1", drop = F],
               offset = base_offset,
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_1, name = " sToL1",na.value = "white",
                       limits  = c(0, 1),
                       breaks  = c(0, 0.25, 0.5, 0.75),
                       labels  = c("0", "0.25", "0.50", "0.75"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 1))
p1
ggsave("../figs/ToL_EF8_heatmap_scale_right.pdf", plot = p1,
       width = 7, height = 7, units = "in") 

p1 <- p1+theme(legend.position = "bottom", legend.box = "horizontal",legend.text = element_text(size = 6))
p1

ggsave("../figs/ToL_EF8_heatmap_scale_bottom.pdf", plot = p1,
       width = 7, height = 7, units = "in") 


# ------------------Extended Figure 9--------------
pal <- c(
  Viridiplantae = "#004300",
  Gastropoda = "#000097",
  Platycheirus = "#ED9B40"
)


anno <- tibble(label = org_labels) %>%
  mutate(Species = str_trim(str_replace(label, "\\s*\\(.*\\)$", "")))  # drop the (...) part
anno$Kingdom <- sapply(anno$Species,function(x){
  unique(metadata$Kingdom[metadata$Species==x])
})

anno$Class <- sapply(anno$Species,function(x){
  unique(metadata$Class[metadata$Species==x])
})
anno$Genus <- sapply(anno$Species,function(x){
  unique(metadata$Genus[metadata$Species==x])
})


anno <- anno %>%
  mutate(ColorGroup = case_when(
    Kingdom %in% names(pal) ~ Kingdom,
    Class  %in% names(pal) ~ Class,
    Genus  %in% names(pal) ~ Genus,
    TRUE ~ NA_character_   # << no group = no colour
  )) %>% filter(!is.na(ColorGroup)) 


grp_list <- split(anno$label, anno$ColorGroup)

hilight_df <- lapply(names(grp_list), function(g) {
  tips <- grp_list[[g]]
  idx  <- which(tr$tip.label %in% tips)
  if (length(idx) >= 2) {
    node_id <- ape::getMRCA(tr, idx)
    if (!is.na(node_id)) data.frame(node = node_id, ColorGroup = g)
  }
}) %>% bind_rows()

p <- ggtree(tr, layout = "circular", size = 0.25,color = "grey60") %<+% anno +
  geom_hilight(data = hilight_df,
               aes(node = node, fill = ColorGroup),
               alpha = 0.5) +
  # geom_tiplab(aes(label = label),
  #             offset = 0.3,
  #             fontface = "italic",
  #             family = "Helvetica",
  #             size = 1) +
  scale_fill_manual(values = pal, name = "Taxonomic group", 
                    guide = guide_legend(override.aes = list(alpha = 1))) #+
# theme(legend.position = "right",
#       text = element_text(family = "Helvetica"))

p

ggsave("../figs/ToL_EF9_noheatmap.pdf", plot = p,
       width = 7, height = 7, units = "in") 



hm <- somatic_signatures %>%
  select(label, sToL8, sToL10, sToL41) %>%             
  distinct(label, .keep_all = TRUE) %>% 
  filter(label %in% org_labels) %>%   
  column_to_rownames("label")


# --- add a NEW fill scale, then the heatmap ring ---
# custom palettes for each signature
pal_8  <- c("white", "#FEBC41", "#EB7822")
pal_10  <- c("white", "#AB739B", "#754668")
pal_41 <- c("white", "#E41A1B", "#A51212")


# layout parameters
base_offset <- 0.05
band_width  <- 0.05
gap         <- 0.3

grid_col  <- "grey85" 


p <- p + guides(fill = "none", colour = "none")
p1 <- p + new_scale_fill()
p1 <- gheatmap(p1, hm[,"sToL8", drop = F],
               offset = base_offset,
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_8, name = "sToL8",na.value = "white",
                       limits  = c(0, 0.7),
                       breaks  = c(0, 0.2, 0.4, 0.6),
                       labels  = c("0", "0.2", "0.4", "0.6"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 1))
p1

p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, hm[,"sToL10", drop = F],
               offset = base_offset + band_width + gap,
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_10, name = "sToL10",na.value = "white",
                       limits  = c(0, 0.45),
                       breaks  = c(0, 0.1, 0.2, 0.3, 0.4),
                       labels  = c("0", "0.1", "0.2", "0.3", "0.4"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 2))
p2

p3 <- p2 + ggnewscale::new_scale_fill()
p3 <- gheatmap(p3, hm["sToL41", drop = FALSE],
               offset = base_offset + 2*(band_width + gap),
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_41, name = "sToL41",na.value = "white",
                       limits  = c(0, 0.4),
                       breaks  = c(0, 0.1, 0.2, 0.3),
                       labels  = c("0", "0.1", "0.2", "0.3"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 4))
p3
ggsave("../figs/ToL_EF9_heatmap_scale_right.pdf", plot = p3,
       width = 7, height = 7, units = "in") 

p3 <- p3+theme(legend.position = "bottom", legend.box = "horizontal",legend.text = element_text(size = 6))
p3

ggsave("../figs/ToL_EF9_heatmap_scale_bottom.pdf", plot = p3,
       width = 7, height = 7, units = "in") 

