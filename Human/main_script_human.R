library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(data.table)
library(stringr)
library(pheatmap)
library(readxl)
library(GSVA)
library(viridis)
library(scales)
library(ComplexHeatmap)
library(DESeq2)
library(circlize)
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)
library(RColorBrewer)
library(ggsci)
library(cowplot)
library(ggpubr)
library(clustree)
library(ggrastr)



### Creating the human mesenchymal atlas ----
load("~/HSC activation/Ramachandran/mesenchyme.RData")
load("~/HSC activation/Payen/Payenhscs_158.RData")
load("~/HSC activation/Fred/HSCs.RData")


ramachandran <- mesenchyme 
payen <- hscs_158
fred <- HSCs

datasets <- c("payen","fred","ramachandran")

for (i in datasets) {
  a <- get(i)
  a@meta.data$orig.ident <- i
  DefaultAssay(a) <- "RNA"
  mt.genes <- rownames(a)[grep("^mt-",rownames(a), ignore.case = TRUE)]
  rb.genes <- rownames(a)[grep("^Rp[sl][[:digit:]]",rownames(a))]
  C<-GetAssayData(object = a, slot = "counts")
  a@meta.data$percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
  a@meta.data$percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
  b <- which(datasets == i)
  assign(x = datasets[b],value = a)
  plot(VlnPlot(a, features = c("nCount_RNA", "nFeature_RNA","percent.mito"), pt.size = 0,group.by = "orig.ident"))
}

Idents(ramachandran) <- "orig.ident"
VlnPlot(ramachandran, "nCount_RNA", y.max = 20000) # 1500
ramachandran <- subset(ramachandran, subset = nCount_RNA < 20000 )

commongenes <- Reduce(intersect, list(rownames(payen),rownames(ramachandran),rownames(fred)))
payen <- subset(payen, features = commongenes)
ramachandran <- subset(ramachandran, features = commongenes)
fred <- subset(fred, features = commongenes)

human <- merge(payen, c(ramachandran,fred), add.cell.ids = c("Payen","Ramachandran","Fred"))
human$dataset <- substr(colnames(human),1,4)


Idents(human) <- "dataset"
human <- RenameIdents(human, "Paye" = "Payen", "Rama" = "Ramachandran", "Fred" = "Fred")
human$dataset <- Idents(human)

human.list <- SplitObject(human, split.by = "dataset")
human.list <- lapply(X = human.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = human.list, nfeatures = 3000)
human.list <- PrepSCTIntegration(object.list = human.list, anchor.features = features)
human.list <- lapply(X = human.list, FUN = RunPCA, features = features)


gc()

anchors <- FindIntegrationAnchors(object.list = human.list, normalization.method = "SCT",
                                  anchor.features = features, dims = 1:20,
                                  reduction = "rpca", k.anchor = 5)
human.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

DefaultAssay(human.sct) <- "RNA"
mt.genes <- rownames(human.sct)[grep("^mt-",rownames(human.sct), ignore.case = TRUE)]
C<-GetAssayData(object = human.sct, slot = "counts")
human.sct$percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100

rb.genes <- rownames(human.sct)[grep("^Rp[sl][[:digit:]]",rownames(human.sct))]
C<-GetAssayData(object = human.sct, slot = "counts")
human.sct$percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100


DefaultAssay(human.sct) <- "integrated"
human.sct <- RunPCA(human.sct)
human.sct <- FindNeighbors(human.sct, reduction = "pca", dims = 1:30)
DefaultAssay(human.sct) <- "integrated"
human.sct <- FindClusters(human.sct, resolution = 0.20)
human.sct <- RunUMAP(human.sct, reduction = "pca", dims = 1:30)
DimPlot(human.sct)
DimPlot(human.sct, group.by = "dataset")

print(human.sct[["pca"]], dims = 1:5, nfeatures = 20)

for (i in unique(human.sct$seurat_clusters)) {
  Idents(human.sct) <- "seurat_clusters"
  .markers <- FindMarkers(human.sct, ident.1 = i, min.pct = 0, only.pos = T)
  assign(paste0("cluster_",i), .markers)
}

Idents(human.sct) <- "seurat_clusters"
human.sct <- RenameIdents(human.sct, "0" = "VSMC", "1" = "B-cells","2" =  "HSC","3" = "Undefined",
                         "4" =  "HSC","5" = "Undefined","6" = "Undefined","7" = "Mesothelial cells", "8" = "Undefined")
DimPlot(human.sct)
human.sct$celltype <- Idents(human.sct)
human.sct$celltype <- factor(human.sct$celltype, levels = c("HSC","B-cells","VSMC","Undefined","Mesothelial cells"))
Idents(human.sct) <- "celltype"

### Fig. S7: Validation of celltypes ----

pdf(file = "~/HSC activation/Paper/Suppfigures/S10/datasetdim.pdf", width = 5, height = 4)
DimPlot(human.sct, label = F, group.by = "dataset") + NoAxes() + ggtitle("")
dev.off()

pdf(file = "~/HSC activation/Paper/Suppfigures/S10/celltypedim.pdf", width = 5, height = 4)
DimPlot(human.sct, label = F, group.by = "celltype") + NoAxes() + ggtitle("")
dev.off()

lijst <- list("HSC \n markers" = c("DCN", "COLEC11"),
              "B-cell \n markers" = c("IGKC", "IGHA1"), 
              "VSMC \n markers" = c("MYH11", "CNN1"), 
              "Mesothelial \n markers" = c("UPK3B","KRT8"))
Idents(human.sct) <- "celltype"
DotPlot(human.sct, features = lijst) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file = "~/HSC activation/Paper/Suppfigures/S10/upk.pdf", width = 4, height = 4)
FeaturePlot(human.sct, c("UPK3B")) + scale_color_viridis_c() + NoAxes()
dev.off()
pdf(file = "~/HSC activation/Paper/Suppfigures/S10/dcn.pdf", width = 4, height = 4)
FeaturePlot(human.sct, c("DCN")) + scale_color_viridis_c() + NoAxes()
dev.off()
pdf(file = "~/HSC activation/Paper/Suppfigures/S10/myh.pdf", width = 4, height = 4)
FeaturePlot(human.sct, c("MYH11")) + scale_color_viridis_c() + NoAxes()
dev.off()
pdf(file = "~/HSC activation/Paper/Suppfigures/S10/igha.pdf", width = 4, height = 4)
FeaturePlot(human.sct, c("IGHA1")) + scale_color_viridis_c() + NoAxes()
dev.off()



### Fig. 4A: Human HSC Atlas ----
human.hscs <- subset(human.sct, idents = c(2,4))
DimPlot(human.hscs)

DefaultAssay(human.hscs) <- "integrated"
human.hscs <- RunPCA(human.hscs)
human.hscs <- FindNeighbors(human.hscs, reduction = "pca", dims = 1:30)
DefaultAssay(human.hscs) <- "integrated"
human.hscs <- FindClusters(human.hscs, resolution = 0.20)
human.hscs <- RunUMAP(human.hscs, reduction = "pca", dims = 1:30)
DimPlot(human.hscs)
DimPlot(human.hscs, group.by = "dataset")


DefaultAssay(human.hscs) <- "RNA"
for (i in unique(human.hscs$seurat_clusters)) {
  Idents(human.hscs) <- "seurat_clusters"
  .markers <- FindMarkers(human.hscs, ident.1 = i, min.pct = 0, only.pos = T)
  assign(paste0("HSC_cluster_",i), .markers)
}

Idents(human.hscs) <- "seurat_clusters"
human.hscs <- RenameIdents(human.hscs, "0" = "initiatory", "1" = "quiescent", "2" = "myofibroblasts")
human.hscs$subtype <- Idents(human.hscs)
human.hscs$subtype <- factor(human.hscs$subtype, levels = c("quiescent","initiatory","myofibroblasts") )
Idents(human.hscs) <- human.hscs$subtype 


pdf(file = "~/HSC activation/Paper/figures/Fig_3/Dimplot.pdf", width = 3.5, height = 3.5)
dimplt_hsc <- DimPlot(human.hscs, group.by = "subtype", label = F, raster = F) + NoLegend() + NoAxes() +ggtitle("")
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off()


### Fig. 4B,D:  heatmap ----
for (i in unique(human.hscs$subtype)) {
  .generalmarkers <- FindMarkers(human.hscs,
                                 ident.1 = i,
                                 assay = "RNA",
                                 min.pct = 0,
                                 only.pos = TRUE,
                                 logfc.threshold = 1)
  tocompare <- unique(human.hscs$subtype)[!unique(human.hscs$subtype) %in% i ]
  .specific1 <- FindMarkers(human.hscs,
                            ident.1 = i,
                            ident.2 = tocompare[1],
                            min.pct = 0,
                            assay = "RNA",
                            only.pos = TRUE,
                            logfc.threshold = 1)
  .specific2 <- FindMarkers(human.hscs,
                            ident.1 = i,
                            ident.2 = tocompare[2],
                            min.pct = 0,
                            assay = "RNA",
                            only.pos = TRUE,
                            logfc.threshold = 1)
  .generalmarkers <- .generalmarkers[rownames(.generalmarkers) %in% intersect(rownames(.specific1),rownames(.specific2)),]

  assign(paste(i,"geneset_human", sep = "__"), .generalmarkers)
  
}


toplot <- c(rownames(quiescent__geneset_human), rownames(initiatory__geneset_human), rownames(myofibroblasts__geneset_human))

mat<- human.hscs[["RNA"]]@data[toplot, ] %>% as.matrix()
mat<- t(scale(t(mat)))
cluster_anno <- human.hscs@meta.data$subtype
quantile(mat, c(0.1, 0.95))

col_fun = circlize::colorRamp2(c(-1, 0, 2), rev(c("#d73027","#ffffbf","#4575b4")))


markers <- c("COLEC10","RELN","ECM1",
             "CCL2","THBS1","TIMP1",
             "COL1A1","S100A6","WT1")




ra <- rowAnnotation(foo = anno_mark(at = match(markers,rownames(mat)),
                                    labels = markers))

d <- Heatmap(mat, name = "Expression",  
             column_split = factor(cluster_anno),
             cluster_columns = F,
             show_column_dend = FALSE,
             right_annotation = ra,
             cluster_column_slices = F,
             column_title_gp = gpar(fontsize = 8),
             column_gap = unit(0.5, "mm"),
             cluster_rows = FALSE,
             show_row_dend = FALSE,
             col = col_fun,
             show_row_names = F,
             row_names_gp = gpar(fontsize = 7),
             column_title_rot = 0,
             top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#f8786f","#00ba38","#619cff")))),
             show_column_names = FALSE,
             use_raster = T,
             raster_quality = 10)

d

pdf(file = "~/HSC activation/Paper/figures/Fig_3/Heatmap.pdf", width = 4, height = 4)
plot(d)
dev.off()
DefaultAssay(human.hscs) <- "RNA"

pdf(file = "~/HSC activation/Paper/figures/Fig_3/hsc_rgs5.pdf", width = 3.5, height = 3)
hsc_rgs5 <- FeaturePlot(human.hscs, "RGS5") + scale_color_viridis() + NoAxes() + ggtitle("")
hsc_rgs5 <- rasterize(hsc_rgs5, layers='Point', dpi=300)
plot(hsc_rgs5)
dev.off()

pdf(file = "~/HSC activation/Paper/figures/Fig_3/hsc_ccl2.pdf", width = 3.5, height = 3)
hsc_ccl2 <- FeaturePlot(human.hscs, "CCL2") + scale_color_viridis() + NoAxes() + ggtitle("")
hsc_ccl2 <- rasterize(hsc_ccl2, layers='Point', dpi=300)
plot(hsc_ccl2)
dev.off()

pdf(file = "~/HSC activation/Paper/figures/Fig_3/hsc_col1a1.pdf", width = 3.5, height = 3)
hsc_col1a1 <- FeaturePlot(human.hscs, "COL1A1") + scale_color_viridis() + NoAxes() + ggtitle("")
hsc_col1a1 <- rasterize(hsc_col1a1, layers='Point', dpi=300)
plot(hsc_col1a1)
dev.off()


quiescent_df <- data.frame("Gene" = rownames(quiescent__geneset_human),"log2FC" = quiescent__geneset_human$avg_log2FC)
quiescent_df[,"HSC stage"] <- "quiescent"

initiatory_df <- data.frame("Gene" = rownames(initiatory__geneset_human),"log2FC" = initiatory__geneset_human$avg_log2FC)
initiatory_df[,"HSC stage"] <- "initiatory"

myofibroblasts_df <- data.frame("Gene" = rownames(myofibroblasts__geneset_human),"log2FC" = myofibroblasts__geneset_human$avg_log2FC)
myofibroblasts_df[,"HSC stage"] <- "myofibroblasts"
overview_df <- rbind(quiescent_df,initiatory_df, myofibroblasts_df)
overview_df <- overview_df[,c(3,1,2)]
write.csv(overview_df,file= "~/HSC activation/Paper/Excel_magic/Human_gene_overview.csv", row.names = F)


### Fig. 4C: Comparison between mouse and human HSC subtypes ----

load("~/HSC activation/Paper/RData_objects/quiescent_geneset.RData")
load("~/HSC activation/Paper/RData_objects/initiatory_geneset.RData")
load("~/HSC activation/Paper/RData_objects/myofibroblasts_geneset.RData")



library(nichenetr)

for (i in c("quiescent_geneset","initiatory_geneset","myofibroblasts_geneset")) {
  .module <- convert_mouse_to_human_symbols(rownames(get(i)))
  .module <- unname(.module)
  .module <- .module[!is.na(.module)]
  .module <- list(".module" = .module)
  
  human.hscs <- AddModuleScore(human.hscs, features = .module, name = i, assay = "RNA")
}

pdf(file = "~/HSC activation/Paper/figures/Fig_3/quiescent_enrichment.pdf", width = 3.5, height = 3)
VlnPlot(human.hscs, c("quiescent_geneset1"), pt.size = 0) +
  theme(axis.title.x = element_blank()) & labs(title = " ")  & theme(legend.position = "none", 
                                                                     axis.text.x = element_text(size = rel(1), angle = 0), 
                                                                     axis.title.y = element_text(size = rel(1), angle = 0), 
                                                                     axis.text.y = element_text(size = rel(1)) )

dev.off()

pdf(file = "~/HSC activation/Paper/figures/Fig_3/initiatory_enrichment.pdf", width = 3.5, height = 3)
VlnPlot(human.hscs, c("initiatory_geneset1"), pt.size = 0) +
  theme(axis.title.x = element_blank()) & labs(title = " ")  & theme(legend.position = "none", 
                                                                     axis.text.x = element_text(size = rel(1), angle = 0), 
                                                                     axis.title.y = element_text(size = rel(1), angle = 0), 
                                                                     axis.text.y = element_text(size = rel(1)) )

dev.off()

pdf(file = "~/HSC activation/Paper/figures/Fig_3/myofibroblasts_enrichment.pdf", width = 3.5, height = 3)
VlnPlot(human.hscs, c("myofibroblasts_geneset1"), pt.size = 0) +
  theme(axis.title.x = element_blank()) & labs(title = " ")  & theme(legend.position = "none", 
                                                                     axis.text.x = element_text(size = rel(1), angle = 0), 
                                                                     axis.title.y = element_text(size = rel(1), angle = 0), 
                                                                     axis.text.y = element_text(size = rel(1)) )

dev.off()

pdf(file = "~/HSC activation/Paper/figures/Fig_3/quiescent_fp.pdf", width = 4, height = 3.5)
dimplt_hsc <- FeaturePlot(human.hscs, c("quiescent_geneset1"), order = T) +
  scale_color_gradient(low = "#f2e2df",high = "#42419d") + NoAxes() + ggtitle("")
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off()

pdf(file = "~/HSC activation/Paper/figures/Fig_3/initiatory_fp.pdf", width = 4, height = 3.5)
dimplt_hsc <- FeaturePlot(human.hscs, c("initiatory_geneset1"), order = T) + 
  scale_color_gradient(low = "#f2e2df",high = "#42419d") + NoAxes()+ ggtitle("")
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off()


pdf(file = "~/HSC activation/Paper/figures/Fig_3/myofibroblast_fp.pdf", width = 4, height = 3.5)
dimplt_hsc <- FeaturePlot(human.hscs, c("myofibroblasts_geneset1"), order = T) +
  scale_color_gradient(low = "#f2e2df",high = "#42419d") + NoAxes()+ ggtitle("")
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off()


### SCENIC
library(data.table)
library(tidyverse)

data_to_write <- as.data.frame(as.matrix(human.hscs@assays$RNA@counts))
dim(data_to_write)
fwrite(x = data_to_write, file = "human.hscs.csv", row.names = T)

auc_mtx <- read_delim("~/HSC activation/Human_integrated/SCENIC/DATA/auc_mtx.tsv", 
                           delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
auc_mtx <- as.data.frame(auc_mtx)
rownames(auc_mtx) <- auc_mtx$Cell
auc_mtx <- auc_mtx[,-c(1)]

b <- t(auc_mtx)
auc_assay <- CreateAssayObject(counts = b)
human.hscs[["AUC"]] <- auc_assay

regulons <- read_csv("~/HSC activation/Human_integrated/SCENIC/DATA/regulons.csv")

alltfs <- list()
for (i in rownames(human.hscs@assays$AUC@counts)) {
  print(i)
  .a <- substring(i,1, nchar(i)-3)
  alltfs[[.a]] <- as.data.frame(regulons[regulons$name == i ,2])[,1]
}

lengtes <- lengths(alltfs)
names(lengtes) <- paste0(names(lengtes), "(+)")

DefaultAssay(human.hscs) <- "AUC"
Idents(human.hscs) <- "subtype"
for (i in unique(human.hscs$subtype)) {
  .generalmarkers <- FindMarkers(human.hscs,
                                 ident.1 = i,
                                 assay = "AUC",
                                 only.pos = TRUE,
                                 logfc.threshold = 0)
  tocompare <- unique(human.hscs$subtype)[!unique(human.hscs$subtype) %in% i ]
  .specific1 <- FindMarkers(human.hscs,
                            ident.1 = i,
                            ident.2 = tocompare[1],
                            assay = "AUC",
                            only.pos = TRUE,
                            logfc.threshold = 0)
  .specific2 <- FindMarkers(human.hscs,
                            ident.1 = i,
                            ident.2 = tocompare[2],
                            assay = "AUC",
                            only.pos = TRUE,
                            logfc.threshold = 0)
  .generalmarkers <- .generalmarkers[rownames(.generalmarkers) %in% intersect(rownames(.specific1),rownames(.specific2)),]
  .generalmarkers$regulonsize <-lengtes[match(rownames(.generalmarkers),names(lengtes))]
  
  assign(paste(i,"markers", sep = "__"), .generalmarkers)
  
}

### Fig. 5B,C: TRajectory inference ----

DefaultAssay(human.hscs) <- "AUC"
c <- DietSeurat(human.hscs, assays = "AUC",dimreducs = c("umap"))
sce <- as.SingleCellExperiment(c)
DimPlot(c, group.by = "subtype")

plot(reducedDim(sce,"UMAP"), col = brewer.pal(9,"Set1")[sce$subtype], pch=16, asp = 0.5)
# Calculate trajectory
sce <- getLineages(sce, clusterLabels = 'subtype', reducedDim = 'UMAP', start.clus = "quiescent")
sce <- getCurves(sce)

# plot the trajectory
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$UMAP, col = plotcol, pch=16)
lines(SlingshotDataSet(sce), lwd=2, col='black')

plot(reducedDims(sce)$UMAP, col = brewer.pal(9,'Set1')[sce$subtype], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

set.seed(3)
sce <- fitGAM(sce)
human.hscs$pseudotime <- sce$slingPseudotime_1

FeaturePlot(human.hscs, "pseudotime") + scale_color_viridis_c() + NoAxes()
FeaturePlot(wo_intermediate, "pseudotime") + scale_color_viridis_c() + NoAxes()

pdf(file = "~/HSC activation/Paper/figures/Fig_5/humanPT.pdf", width = 3, height = 2.5)
dimplt_hsc <- FeaturePlot(human.hscs, "pseudotime", cols = c("#f8786f","#619cff")) + NoAxes() + ggtitle("")
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off()

#### Make scatterplots

top3 <- c("RXRA(+)","FOXF1(+)","NR1H4(+)",
          "FOSL1(+)","EGR3(+)","NFKB2(+)",
          "WT1(+)", "PRRX1(+)")


DefaultAssay(human.hscs) <- "AUC"
hm_df <- FetchData(human.hscs, vars = top3 )
hm_df <- as.data.frame(scale(hm_df))
hm_df$pseudotime <- human.hscs$pseudotime


pdf(file = "~/HSC activation/Paper/figures/Fig_5/hQTFs.pdf", width = 5, height = 2.52)
sc_plot <- ggscatter(hm_df,"pseudotime", top3[1:3], size = 0, merge = T, point =F, 
                     add = "loess",conf.int = F) + labs(x = NULL, y = NULL) 
ggpar(sc_plot, legend = "right",
      font.legend = c(12.5, "plain", "black"))
dev.off()

pdf(file = "~/HSC activation/Paper/figures/Fig_5/hITFS.pdf", width = 5, height = 2.52)
sc_plot <- ggscatter(hm_df,"pseudotime", top3[4:6], size = 0, merge = T, point =F, 
                     add = "loess",conf.int = F) + labs(x = NULL, y = NULL) 
ggpar(sc_plot, legend = "right",
      font.legend = c(12.5, "plain", "black"))
dev.off()

pdf(file = "~/HSC activation/Paper/figures/Fig_5/hMTFS.pdf", width = 5, height = 2.52)
sc_plot <- ggscatter(hm_df,"pseudotime", top3[7:8], size = 0, merge = T, point =F, 
                     add = "loess",conf.int = F) + labs(x = NULL, y = NULL) 
ggpar(sc_plot, legend = "right",
      font.legend = c(12.5, "plain", "black"),
      palette = c("#ee756e","#3bab3c"))
dev.off()


### Fig. S8: Gene markers in the human dataset -----

markers2 <- c("Reln","Rgs5", "Ecm1","Pth1r","Vipr1","Calcrl","Agtr1a", "Ednrb",
              "Colec10","Colec11","Fcna","Ets1",
              "Thbs1","Ccl2","Timp1",
              "Col1a1","S100a6","Wt1","C3")

hmarkers <-
  unname(convert_mouse_to_human_symbols(markers2))
hmarkers[11] <- "FCN1"
hmarkers[14] <- "CCL2"
hmarkers <- hmarkers[!is.na(hmarkers)]


cols = circlize::colorRamp2(c(-1.2, 0, 1.2), c("#fdfed5", "#72c8bd","#253997"))

data_hm <- AverageExpression(human.hscs, assays = "RNA")
df <- as.data.frame(data_hm$RNA)


subtype_annotation <- c("quiescent","initiatory","myofibroblasts")
subtype_annotation <- factor(subtype_annotation, levels = c("quiescent","initiatory","myofibroblasts"))


data_hm <- df[hmarkers,]
data_hm <- t(scale(t(data_hm)))
data_ha <- HeatmapAnnotation("subtype" = subtype_annotation,
                             col = list(subtype = c("myofibroblasts" = "#619cff", "initiatory" = "#00ba38","quiescent" = "#f8766c")))
a <- Heatmap(data_hm, col = cols,
             row_title = "",
             cluster_columns = F,
             cluster_rows = F,
             top_annotation = data_ha)


ht_list = hm_APAP + hm_TAA + hm_CCl4_acute + hm_WD + hm_BDL + hm_CCl4_Dobie + hm_CCl4_Krenkel + hm_CCl4_Terkelsen + hm_CCl4_Yang + hm_Cancer
pdf(file = "~/HSC activation/REintegration/Cancer_inc/figures/qui_TFs_hm.pdf", width = 16, height = 8)
draw(a)

dev.off()
