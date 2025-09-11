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



germline_signatures <- read.csv("../data/germline_mutational_signature_attributions.x0_excluded.csv",check.names = F)
colnames(germline_signatures) = paste0('gToL',colnames(germline_signatures))
germline_signatures[,1] <- ifelse(grepl("\\.", germline_signatures[,1]),
                                 sub(".*\\.", "", germline_signatures[,1]),
                                 germline_signatures[,1])
germline_signatures$label=NA
germline_signatures$label=sapply(germline_signatures[,1], function(x){
  Species = metadata$Species[metadata$Sample==x]
  if (length(metadata$Species[metadata$Species==Species]) == 1){
    label = Species
  }
  else{
    label = paste0(Species," (", x, ")")
  } 
  label
})

germline_signatures$label[germline_signatures$label=="Hemaris fuciformis"]="Hemaris fuciformis (iHemFuc2)"



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
# ------------------Fig3 taxa enriched signatures--------------
pal <- c(
  Coleoptera = "#0100F8",
  Chordata = "#30A2DA",
  Viridiplantae = "#004300",
  Vespidae = "#EE9480",
  Tenthredinidae = "#EFD2CB"
)


anno <- tibble(label = org_labels) %>%
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
anno$Family <- sapply(anno$Species,function(x){
  unique(metadata$Family[metadata$Species==x])
})

anno <- anno %>%
  mutate(ColorGroup = case_when(
    Kingdom %in% names(pal) ~ Kingdom,
    Phylum %in% names(pal) ~ Phylum,
    Class  %in% names(pal) ~ Class,
    Order  %in% names(pal) ~ Order,
    Family  %in% names(pal) ~ Family,
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

ggsave("../figs/ToL_Fig3_noheatmap.pdf", plot = p,
       width = 7, height = 7, units = "in") 



hm <- somatic_signatures %>%
  select(label, sToL4, sToL5, sToL15, sToL24) %>%             
  distinct(label, .keep_all = TRUE) %>% 
  filter(label %in% org_labels) %>%   
  column_to_rownames("label")


# --- add a NEW fill scale, then the heatmap ring ---
# custom palettes for each signature
pal_4  <- c("white", "skyblue", "navy")
pal_5  <- c("white", "#E41A1B", "#A51212")
pal_15 <- c("white", "#AB739B", "#754668")
pal_24 <- c("white", "#FEBC41", "#EB7822")

# layout parameters
base_offset <- 0.05
band_width  <- 0.05
gap         <- 0.3

grid_col  <- "grey85" 


p <- p + guides(fill = "none", colour = "none")
p1 <- p + new_scale_fill()
p1 <- gheatmap(p1, hm[,"sToL4", drop = F],
               offset = base_offset,
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_4, name = "sToL4",na.value = "white",
                       limits  = c(0, 0.7),
                       breaks  = c(0, 0.2, 0.4, 0.6),
                       labels  = c("0", "0.2", "0.4", "0.6"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 1))
p1

p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, hm[,"sToL5", drop = F],
               offset = base_offset + band_width + gap,
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_5, name = "sToL5",na.value = "white",
                       limits  = c(0, 0.8),
                       breaks  = c(0, 0.25, 0.5, 0.75),
                       labels  = c("0", "0.25", "0.50", "0.75"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 2))
p2

p3 <- p2 + ggnewscale::new_scale_fill()
p3 <- gheatmap(p3, hm["sToL15", drop = FALSE],
               offset = base_offset + 2*(band_width + gap),
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_15, name = "sToL15",na.value = "white",
                       limits  = c(0, 0.4),
                       breaks  = c(0, 0.1, 0.2, 0.3, 0.4),
                       labels  = c("0", "0.1", "0.2", "0.3","0.4"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 3))
p3

p4 <- p3 + ggnewscale::new_scale_fill()
p4 <- gheatmap(p4, hm["sToL24", drop = FALSE],
               offset = base_offset + 3*(band_width + gap),
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_24, name = "sToL24",na.value = "white",
                       limits  = c(0, 0.7),
                       breaks  = c(0, 0.2, 0.4, 0.6),
                       labels  = c("0", "0.2", "0.4", "0.6"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 4))
p4
ggsave("~/Documents/Postdoc/ToL/figs/ToL_Fig3_heatmap_scale_right.pdf", plot = p4,
       width = 7, height = 7, units = "in") 

p4 <- p4+theme(legend.position = "bottom", legend.box = "horizontal",legend.text = element_text(size = 6))
p4

ggsave("~/Documents/Postdoc/ToL/figs/ToL_Fig3_heatmap_scale_bottom.pdf", plot = p4,
       width = 7, height = 7, units = "in") 


# ------------------Fig4--------------
pal <- c(
  Actinopteri  = "#D796AB",
  Mollusca = "#000097",
  Actiniaria = "#ED9B40",
  Echinodermata = "#17BEBB"
)


anno <- tibble(label = org_labels) %>%
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
anno$Family <- sapply(anno$Species,function(x){
  unique(metadata$Family[metadata$Species==x])
})

anno <- anno %>%
  mutate(ColorGroup = case_when(
    Kingdom %in% names(pal) ~ Kingdom,
    Phylum %in% names(pal) ~ Phylum,
    Class  %in% names(pal) ~ Class,
    Order  %in% names(pal) ~ Order,
    Family  %in% names(pal) ~ Family,
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

ggsave("../figs/ToL_Fig4_noheatmap.pdf", plot = p,
       width = 7, height = 7, units = "in") 



hm <- somatic_signatures %>%
  select(label, sToL9, sToL13, sToL18) %>%             
  distinct(label, .keep_all = TRUE) %>% 
  filter(label %in% org_labels) %>%   
  column_to_rownames("label")


# --- add a NEW fill scale, then the heatmap ring ---
# custom palettes for each signature
pal_9  <- c("white", "skyblue", "navy")
pal_13  <- c("white", "#E41A1B", "#A51212")
pal_18 <- c("white", "#AB739B", "#754668")


# layout parameters
base_offset <- 0.05
band_width  <- 0.05
gap         <- 0.3

grid_col  <- "grey85" 


p <- p + guides(fill = "none", colour = "none")
p1 <- p + new_scale_fill()
p1 <- gheatmap(p1, hm[,"sToL9", drop = F],
               offset = base_offset,
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_9, name = "sToL9",na.value = "white",
                       limits  = c(0, 0.5),
                       breaks  = c(0, 0.1, 0.2, 0.3, 0.4),
                       labels  = c("0", "0.1", "0.2", "0.3", "0.4"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 1))
p1

p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, hm[,"sToL13", drop = F],
               offset = base_offset + band_width + gap,
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_13, name = "sToL13",na.value = "white",
                       limits  = c(0, 0.6),
                       breaks  = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                       labels  = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 2))
p2

p3 <- p2 + ggnewscale::new_scale_fill()
p3 <- gheatmap(p3, hm["sToL18", drop = FALSE],
               offset = base_offset + 2*(band_width + gap),
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_18, name = "sToL18",na.value = "white",
                       limits  = c(0, 0.5),
                       breaks  = c(0, 0.1, 0.2, 0.3, 0.4),
                       labels  = c("0", "0.1", "0.2", "0.3", "0.4"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 4))
p3
ggsave("~/Documents/Postdoc/ToL/figs/ToL_Fig4_heatmap_scale_right.pdf", plot = p3,
       width = 7, height = 7, units = "in") 

p3 <- p3+theme(legend.position = "bottom", legend.box = "horizontal",legend.text = element_text(size = 6))
p3

ggsave("~/Documents/Postdoc/ToL/figs/ToL_Fig4_heatmap_scale_bottom.pdf", plot = p3,
       width = 7, height = 7, units = "in") 

# ------------------Fig5--------------
pal <- c(
  Coleoptera = "#0100F8",
  Chordata = "#30A2DA",
  Viridiplantae = "#004300",
  Fungi        = "#8C564B"
)


anno <- tibble(label = org_labels) %>%
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
anno$Family <- sapply(anno$Species,function(x){
  unique(metadata$Family[metadata$Species==x])
})

anno <- anno %>%
  mutate(ColorGroup = case_when(
    Kingdom %in% names(pal) ~ Kingdom,
    Phylum %in% names(pal) ~ Phylum,
    Class  %in% names(pal) ~ Class,
    Order  %in% names(pal) ~ Order,
    Family  %in% names(pal) ~ Family,
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

ggsave("../figs/ToL_Fig5_noheatmap.pdf", plot = p,
       width = 7, height = 7, units = "in") 



hm <- germline_signatures %>%
  select(label, gToL1, gToL3, gToL4, gToL6) %>%             
  distinct(label, .keep_all = TRUE) %>% 
  filter(label %in% org_labels) %>%   
  column_to_rownames("label")


# --- add a NEW fill scale, then the heatmap ring ---
# custom palettes for each signature
pal_1  <- c("white", "skyblue", "navy")
pal_3  <- c("white", "#E41A1B", "#A51212")
pal_4 <- c("white", "#AB739B", "#754668")
pal_6 <- c("white", "#FEBC41", "#EB7822")

# layout parameters
base_offset <- 0.05
band_width  <- 0.05
gap         <- 0.3

grid_col  <- "grey85" 


p <- p + guides(fill = "none", colour = "none")
p1 <- p + new_scale_fill()
p1 <- gheatmap(p1, hm[,"gToL1", drop = F],
               offset = base_offset,
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_1, name = "gToL1",na.value = "white",
                       limits  = c(0, 1),
                       breaks  = c(0, 0.25, 0.5, 0.75),
                       labels  = c("0", "0.25", "0.50", "0.75"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 1))
p1

p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, hm[,"gToL3", drop = F],
               offset = base_offset + band_width + gap,
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_3, name = "gToL3",na.value = "white",
                       limits  = c(0, 0.8),
                       breaks  = c(0, 0.25, 0.5, 0.75),
                       labels  = c("0", "0.25", "0.50", "0.75"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 2))
p2

p3 <- p2 + ggnewscale::new_scale_fill()
p3 <- gheatmap(p3, hm["gToL4", drop = FALSE],
               offset = base_offset + 2*(band_width + gap),
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_4, name = "gToL4",na.value = "white",
                       limits  = c(0, 0.65),
                       breaks  = c(0,  0.2, 0.4, 0.6),
                       labels  = c("0", "0.2", "0.4", "0.6"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 3))
p3

p4 <- p3 + ggnewscale::new_scale_fill()
p4 <- gheatmap(p4, hm["gToL6", drop = FALSE],
               offset = base_offset + 3*(band_width + gap),
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_6, name = "gToL6",na.value = "white",
                       limits  = c(0, 0.4),
                       breaks  = c(0, 0.1, 0.2, 0.3),
                       labels  = c("0", "0.1", "0.2", "0.3"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 4))
p4
ggsave("~/Documents/Postdoc/ToL/figs/ToL_Fig5_heatmap_scale_right.pdf", plot = p4,
       width = 7, height = 7, units = "in") 

p4 <- p4+theme(legend.position = "bottom", legend.box = "horizontal",legend.text = element_text(size = 6))
p4

ggsave("~/Documents/Postdoc/ToL/figs/ToL_Fig5_heatmap_scale_bottom.pdf", plot = p4,
       width = 7, height = 7, units = "in") 

