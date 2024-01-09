### Integrating pySCENIC results in the nichenetr pipeline ----
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
sig_network = readRDS(url("https://zenodo.org/record/3260758/files/signaling_network.rds"))
gr_network = readRDS(url("https://zenodo.org/record/3260758/files/gr_network.rds"))
library(nichenetr)
library(readr)
regulons <- read_csv("~/HSC activation/Paper/data/regulons.csv")
niche_scenic <- regulons
colnames(niche_scenic) <- c("from","to")
niche_scenic$from <- substr(niche_scenic$from,1,nchar(niche_scenic$from)-3)
niche_scenic$source <- "scenic"
niche_scenic$database <- "scenic"
niche_scenic$from <- convert_mouse_to_human_symbols(niche_scenic$from )
niche_scenic$to <- convert_mouse_to_human_symbols(niche_scenic$to )

library(tidyr)
niche_scenic <- drop_na(niche_scenic)

gr_network2 <- rbind(gr_network, niche_scenic)

source_weights_df <- nichenetr::source_weights_df
source_weights_df[nrow(source_weights_df)+1,] <- list("scenic",1)


# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks = construct_weighted_networks(
  lr_network = lr_network,
  sig_network = sig_network,
  gr_network = gr_network2,
  source_weights_df = source_weights_df
)

# downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks = apply_hub_corrections(
  weighted_networks = weighted_networks,
  lr_sig_hub = hyperparameter_list[hyperparameter_list$parameter == "lr_sig_hub","avg_weight"],
  gr_hub = hyperparameter_list[hyperparameter_list$parameter == "gr_hub","avg_weight"]
)


x <- 1
ligandlist <- list()
for (i in lr_network$from) {
  ligandlist[[x]] <- i
  x <- x +1 
}

ligand_target_matrix = construct_ligand_target_matrix(
  weighted_networks = weighted_networks,
  ligands = ligandlist,
  algorithm = "PPR",
  damping_factor = hyperparameter_list$damping_factor,
  ltf_cutoff = hyperparameter_list$ltf_cutoff
)

save(ligand_target_matrix, file = "ligand_target_matrix.RData")
### Prediicting Top ligand activities ----
library(nichenetr)
library(Seurat) 
library(tidyverse)

#lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

Idents(hscs.filtered) <- "celltype"
expressed_genes_receiver = get_expressed_genes(ident = "HSC", seurat_obj = hscs.filtered,assay_oi = "RNA")
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

DefaultAssay(hscs.filtered) <- "RNA"
Idents(hscs.filtered) <- "subtype"
quiescent_to_mfb <- FindMarkers(hscs.filtered,
                                ident.1 = "myofibroblasts",
                                ident.2 = c("quiescent"),
                                group.by = "subtype",
                                logfc.threshold = log2(2),
                                min.pct = 0,
                                only.pos = T) %>% rownames_to_column("gene")
geneset_oi <- quiescent_to_mfb$gene[quiescent_to_mfb$gene %in% rownames(ligand_target_matrix)]
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = ligands
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) 

ligand_oi <- ligand_activities$test_ligand[1:10]


### Fig 6B: Top ligands inducing target genes  ----
best_upstream_ligands = ligand_activities %>% top_n(10, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 500) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network

save(ligand_activities, file = "ligand_activities_qui_to_mfb.RData")
save(lr_network_top_df, file = "lr_network_top_df.RData")


### Top ligands inducing activation per study ----
for (i in unique(hscs.filtered$study)) {
  Idents(hscs.filtered) <- "study"
  DefaultAssay(hscs.filtered) <- "RNA"
  .quiescent_to_myofibroblasts <- FindMarkers(hscs.filtered,
                                              ident.1 = "myofibroblasts",
                                              ident.2 = "quiescent",
                                              subset.ident = i, 
                                              group.by = "subtype",
                                              logfc.threshold = log2(2),
                                              min.pct = 0,
                                              only.pos = T) %>% rownames_to_column("gene")
  geneset_oi <- .quiescent_to_myofibroblasts$gene[.quiescent_to_myofibroblasts$gene %in% rownames(ligand_target_matrix)]
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()
  
  expressed_ligands = ligands
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  
  potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
  
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  
  ligand_activities = ligand_activities %>% arrange(-pearson) 
  ligand_activities$order <- rownames(ligand_activities)
  ligand_activities <- as.data.frame(ligand_activities)
  rownames(ligand_activities) <- ligand_activities$test_ligand
  assign(paste(i,"_quiescent_to_myofibroblasts", sep = ""),ligand_activities)
}


### Comparison ----
ligand_activity_matrix <- data.frame("APAP" = APAP_quiescent_to_myofibroblasts[APAP_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"pearson"],
"TAA" = TAA_quiescent_to_myofibroblasts[TAA_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"pearson"],
"CCl4_acute" = CCl4_acute_quiescent_to_myofibroblasts[CCl4_acute_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"pearson"],
"WD" = WD_quiescent_to_myofibroblasts[WD_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"pearson"],
"BDL" = BDL_quiescent_to_myofibroblasts[BDL_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"pearson"],
"CCl4_Dobie" = CCl4_Dobie_quiescent_to_myofibroblasts[CCl4_Dobie_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"pearson"],
"CCl4_Krenkel" = CCl4_Krenkel_quiescent_to_myofibroblasts[CCl4_Krenkel_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"pearson"],
"CCl4_Terkelsen" = CCl4_Terkelsen_quiescent_to_myofibroblasts[CCl4_Terkelsen_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"pearson"],
"CCl4_Yang" = CCl4_Yang_quiescent_to_myofibroblasts[CCl4_Yang_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"pearson"],
"Cancer" = Cancer_quiescent_to_myofibroblasts[Cancer_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"pearson"]
)

ligand_activity_matrix <- data.frame("APAP" = APAP_quiescent_to_myofibroblasts[APAP_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"order"],
                                     "TAA" = TAA_quiescent_to_myofibroblasts[TAA_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"order"],
                                     "CCl4_acute" = CCl4_acute_quiescent_to_myofibroblasts[CCl4_acute_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"order"],
                                     "WD" = WD_quiescent_to_myofibroblasts[WD_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"order"],
                                     "BDL" = BDL_quiescent_to_myofibroblasts[BDL_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"order"],
                                     "CCl4_Dobie" = CCl4_Dobie_quiescent_to_myofibroblasts[CCl4_Dobie_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"order"],
                                     "CCl4_Krenkel" = CCl4_Krenkel_quiescent_to_myofibroblasts[CCl4_Krenkel_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"order"],
                                     "CCl4_Terkelsen" = CCl4_Terkelsen_quiescent_to_myofibroblasts[CCl4_Terkelsen_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"order"],
                                     "CCl4_Yang" = CCl4_Yang_quiescent_to_myofibroblasts[CCl4_Yang_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"order"],
                                     "Cancer" = Cancer_quiescent_to_myofibroblasts[Cancer_quiescent_to_myofibroblasts$test_ligand %in% ligand_oi,"order"]
)
ligand_activity_matrix <- data.frame("APAP" = as.numeric(APAP_quiescent_to_myofibroblasts[ligand_oi,"order"]),
                                     "TAA" = as.numeric(TAA_quiescent_to_myofibroblasts[ligand_oi,"order"]),
                                     "CCl4_acute" = as.numeric(CCl4_acute_quiescent_to_myofibroblasts[ ligand_oi,"order"]),
                                     "WD" = as.numeric(WD_quiescent_to_myofibroblasts[ligand_oi,"order"]),
                                     "BDL" = as.numeric(BDL_quiescent_to_myofibroblasts[ ligand_oi,"order"]),
                                     "CCl4_Dobie" = as.numeric(CCl4_Dobie_quiescent_to_myofibroblasts[ligand_oi,"order"]),
                                     "CCl4_Krenkel" = as.numeric(CCl4_Krenkel_quiescent_to_myofibroblasts[ ligand_oi,"order"]),
                                     "CCl4_Terkelsen" = as.numeric(CCl4_Terkelsen_quiescent_to_myofibroblasts[ligand_oi,"order"]),
                                     "CCl4_Yang" = as.numeric(CCl4_Yang_quiescent_to_myofibroblasts[ligand_oi,"order"]),
                                     "Cancer" = as.numeric(Cancer_quiescent_to_myofibroblasts[ ligand_oi,"order"])
)

colnames(ligand_activity_matrix) <- c("APAP" ,
                                     "TAA" ,
                                     "CCl4_acute" ,
                                     "WD" ,
                                     "BDL",
                                     "CCl4_Dobie",
                                     "CCl4_Krenkel" ,
                                     "CCl4_Terkelsen" ,
                                     "CCl4_Yang",
                                     "Cancer")
rownames(ligand_activity_matrix) <- ligand_oi

ligand_activity_matrix$means <- unname(rowMeans(ligand_activity_matrix))
ligand_activity_matrix <- ligand_activity_matrix[order(ligand_activity_matrix$means),]
ligand_activity_matrix <- ligand_activity_matrix[,-11]
library(ComplexHeatmap)
library(circlize)

cond_annotation <- colnames(ligand_activity_matrix)
cond_annotation <- factor(cond_annotation, levels = unique(cond_annotation))

cols = colorRamp2(c(20, 10, 1), c("#fdfed5", "#72c8bd","#253997"))
annot <- HeatmapAnnotation("condition" = cond_annotation
                           #,col = list(subtype = c("myofibroblasts" = "#c87dff", "initiatory" = "#18bdb9","quiescent" = "#f8766d"))
)
## In vivo genesets
hmval <- Heatmap(as.matrix(ligand_activity_matrix), col = cols,
                 row_title = "",
                 column_split = factor(c(rep("Acute injury",3),
                                         rep("Chronic injury",7)),
                                       levels = c("Acute injury","Chronic injury")),
                 heatmap_legend_param = list(
                   title = "Ligand rank"),
                 # top_annotation = annot,
                 cluster_columns = T,
                 cluster_rows =F,
                 show_row_names = T,
                 show_column_dend = F)
pdf(file = "~/HSC activation/Paper/figures/Fig_6/ligands_per_study.pdf", width = 5, height = 4) 
hmval
dev.off() 

Idents(hscs.filtered) <- "subtype"
qui_cells <- WhichCells(hscs.filtered, idents = c("quiescent"))
myo_cells <- WhichCells(hscs.filtered, idents = c("myofibroblasts"))



pdf(file = "~/HSC activation/Paper/figures/Fig_6/dimplot.pdf", width = 5.5, height = 5)
hsc_nichenet <-DimPlot(hscs.filtered, cells.highlight= list(qui_cells,
                                                        myo_cells),
                   cols.highlight = c("#619cff","#f8786f"), cols= "lightgrey", 
                   sizes.highlight = 0.3)+ NoLegend()  +NoAxes()  +theme(plot.title = element_text(hjust = 0.5))
hsc_nichenet <- rasterize(hsc_nichenet, layers='Point', dpi=300)
plot(hsc_nichenet)
dev.off()

### Saving ligands per study ----
overview_studies <- data.frame("rank" = 1:20)
for (i in unique(hscs.filtered$study)) {
  a <- get(paste(i,"_quiescent_to_myofibroblasts", sep = ""))
  a <- as.data.frame(a)
  rownames(a) <- a$test_ligand
  assign(paste(i,"_quiescent_to_myofibroblasts", sep = ""),a)
  overview_studies[,i] <- a$test_ligand[1:20]
}
