library(Seurat)
library(dplyr)
library(readr)
library(Matrix)

annot_mouseStStAll <- read_csv("~/HSC activation/Paper/Guilliams/annot_mouseStStAll.csv")
annot_mouseStStAll <- as.data.frame(annot_mouseStStAll)
rownames(annot_mouseStStAll) <- annot_mouseStStAll$cell
annot_mouseStStAll <- annot_mouseStStAll[,-6]

setwd("~/HSC activation/Paper/Guilliams/rawData_mouseStSt/countTable_mouseStSt")

counts<- readMM("matrix.mtx.gz")
genes<- read_tsv("features.tsv.gz",col_names = F)
barcodes <- read_tsv("barcodes.tsv.gz", col_names = F)
rownames(counts) <- genes$X1
colnames(counts) <- barcodes$X1
SO_mouse <- CreateSeuratObject(counts = counts, meta.data = annot_mouseStStAll)

SO_mouse <- NormalizeData(SO_mouse)
SO_mouse <- FindVariableFeatures(SO_mouse)
SO_mouse <- ScaleData(SO_mouse)
SO_mouse <- RunPCA(SO_mouse)
SO_mouse <- FindNeighbors(SO_mouse,dims = 1:30, reduction = "pca")
SO_mouse <- FindClusters(SO_mouse, resolution = 0.05)
SO_mouse <- RunUMAP(SO_mouse, dims = 1:30, reduction = "pca")
Idents(SO_mouse) <-  "annot"


umap <- cbind(SO_mouse$UMAP_1, SO_mouse$UMAP_2)
colnames(umap) <- c("UMAP_1","UMAP_2")
SO_mouse@reductions$umap@cell.embeddings <- umap

unique(annot_mouseStStFibro$annot)

mfibroblasts <- as.vector(annot_mouseStStFibro[annot_mouseStStFibro$annot %in% "Capsule fibroblasts","cell"])
SO_mouse$annot[names(SO_mouse$annot) %in% mfibroblasts$cell] <- "Fibroblasts"

mfibroblasts <- as.vector(annot_mouseStStFibro[annot_mouseStStFibro$annot %in% "Fibroblast 2","cell"])
SO_mouse$annot[names(SO_mouse$annot) %in% mfibroblasts$cell] <- "Fibroblasts"

mfibroblasts <- as.vector(annot_mouseStStFibro[annot_mouseStStFibro$annot %in% "Fibroblast 1","cell"])
SO_mouse$annot[names(SO_mouse$annot) %in% mfibroblasts$cell] <- "Fibroblasts"

mVSMCs <- as.vector(annot_mouseStStFibro[annot_mouseStStFibro$annot %in% "VSMCs","cell"])
SO_mouse$annot[names(SO_mouse$annot) %in% mVSMCs$cell] <- "VSMCs"

mHSCs <- as.vector(annot_mouseStStFibro[annot_mouseStStFibro$annot %in% "Stellate cells","cell"])
SO_mouse$annot[names(SO_mouse$annot) %in% mHSCs$cell] <- "Stellate cells"

meso <- as.vector(annot_mouseStStFibro[annot_mouseStStFibro$annot %in% "Mesothelial cells","cell"])
SO_mouse$annot[names(SO_mouse$annot) %in% meso$cell] <- "Mesothelial cells"


SO_mouse2 <- subset(SO_mouse, cells = rownames(annot_mouseStStAll))
Idents(SO_mouse2) <- "annot"
DimPlot(SO_mouse2, raster = F, group.by = "annot")

unique(Idents(SO_mouse2))
SO_mouse2 <- RenameIdents(SO_mouse2, "Monocytes & Monocyte-derived cells" = "Monocytes",
                         "cDC1s" = "DCs","pDCs" = "DCs" , "cDC2s" = "DCs", "Kupffer cells" = "Macrophages",
                         "Mig. cDCs" = "DCs")

SO_mouse2$celltype <- Idents(SO_mouse2)
DimPlot(SO_mouse, group.by = "annot", raster = F, label = T)

SO_mouse2$celltype <- factor(SO_mouse2$celltype, levels = levels(SO_human2$celltype))

pdf(file = "~/HSC activation/Paper/figures/Fig_colec/Mviolin.pdf", width = 6, height = 4)
VlnPlot(SO_mouse2, "Colec10",pt.size = 0, group.by = "celltype", idents = a) + NoLegend() 
dev.off()

pdf(file = "~/HSC activation/Paper/figures/Fig_6/mPth1r.pdf", width = 3, height = 4)
VlnPlot(SO_mouse2, "Pth1r",pt.size = 0, group.by = "celltype", idents = c("Stellate cells","Hepatocytes")) &
  xlab("")  &
  theme(legend.position = "none", 
        axis.text.x = element_text(size = rel(1), angle = 0), 
        axis.title.y = element_text(size = rel(1), angle = 90), 
        axis.text.y = element_text(size = rel(1)) ) 
dev.off()

pdf(file = "~/HSC activation/Paper/Suppfigures/S13/mousedrivers.pdf", width = 18, height = 3.5)
VlnPlot(SO_mouse2, c("Colec10","Gfap","Lrat","Tcf21","Hand2"),
          pt.size = 0, group.by = "celltype", ncol = 5, idents = a)& NoLegend()&
  xlab("")  &
  theme(legend.position = "none", 
        axis.text.x = element_text(size = rel(1), angle = 45), 
        axis.title.y = element_text(size = rel(1), angle = 90), 
        axis.text.y = element_text(size = rel(1)) ) 
dev.off()

annot_humanAll <- read_csv("~/HSC activation/Paper/Guilliams/annot_humanAll.csv")
annot_humanAll <- as.data.frame(annot_humanAll)
rownames(annot_humanAll) <- annot_humanAll$cell
annot_humanAll <- annot_humanAll[,-7]

setwd("~/HSC activation/Paper/Guilliams/rawData_human/countTable_human")

counts<- readMM("matrix.mtx.gz")
genes<- read_tsv("features.tsv.gz",col_names = F)
barcodes <- read_tsv("barcodes.tsv.gz", col_names = F)
rownames(counts) <- genes$X1
colnames(counts) <- barcodes$X1
SO_human <- CreateSeuratObject(counts = counts, meta.data = annot_humanAll)

SO_human <- NormalizeData(SO_human)
SO_human <- FindVariableFeatures(SO_human)
SO_human <- ScaleData(SO_human)
SO_human <- RunPCA(SO_human)
SO_human <- FindNeighbors(SO_human,dims = 1:30, reduction = "pca")
SO_human <- FindClusters(SO_human, resolution = 0.05)
SO_human <- RunUMAP(SO_human, dims = 1:30, reduction = "pca")
DimPlot(SO_human, raster = F, group.by = "annot")

umap <- cbind(SO_human$UMAP_1, SO_human$UMAP_2)
colnames(umap) <- c("UMAP_1","UMAP_2")
SO_human@reductions$umap@cell.embeddings <- umap
unique(annot_humanCD45neg$annot)

hfibroblasts <- as.vector(annot_humanCD45neg[annot_humanCD45neg$annot == "Fibroblasts","cell"])
SO_human$annot[names(SO_human$annot) %in% hfibroblasts$cell] <- "Fibroblasts"

hVSMCs <- as.vector(annot_humanCD45neg[annot_humanCD45neg$annot == "VSMCs","cell"])
SO_human$annot[names(SO_human$annot)%in% hVSMCs$cell] <- "VSMCs"

hHSCs <- as.vector(annot_humanCD45neg[annot_humanCD45neg$annot == "Stellate cells","cell"])
SO_human$annot[names(SO_human$annot)%in% hHSCs$cell] <- "Stellate cells"

DimPlot(SO_human, raster = F, group.by = "annot")

Idents(SO_human) <-  "annot"
SO_human2 <- subset(SO_human, cells = rownames(annot_humanAll))
SO_human2 <- RenameIdents(SO_human2, "Mono+mono derived cells" = "Monocytes",
                          "Mig.cDCs" = "DCs","Circulating NK/NKT" = "NK cells","cDC1s" ="DCs",
                          "Resident NK" ="NK cells","cDC2s" = "DCs","pDCs" = "DCs")
SO_human2$celltype <- Idents(SO_human2) 
levels(SO_human2$celltype)
SO_human2$celltype <- factor(SO_human2$celltype, levels = rev(levels(SO_human2$celltype)))
b <- levels(SO_human2$celltype)
b[1] <-  "Stellate cells"
b[3] <-  "VSMCs"

SO_human2$celltype <- factor(SO_human2$celltype, levels = b)


DimPlot(SO_human2, raster = F, group.by = "celltype", label = T) 
VlnPlot(SO_human2, "COLEC10",pt.size = 0, group.by = "celltype")+ NoLegend()


pdf(file = "~/HSC activation/Paper/figures/Fig_colec/Hviolin.pdf", width = 6, height = 4)
VlnPlot(SO_human2, "COLEC10",pt.size = 0, group.by = "celltype", idents = a) + NoLegend() 
dev.off()
pdf(file = "~/HSC activation/Paper/figures/Fig_6/pth1rviolin.pdf", width = 3, height = 4)
VlnPlot(SO_human2, "PTH1R",pt.size = 0, group.by = "celltype", idents = c("Stellate cells","Hepatocytes"))&
  xlab("")  &
  theme(legend.position = "none", 
        axis.text.x = element_text(size = rel(1), angle = 0), 
        axis.title.y = element_text(size = rel(1), angle = 90), 
        axis.text.y = element_text(size = rel(1)) ) 
dev.off()

pdf(file = "~/HSC activation/Paper/Suppfigures/S13/humandrivers.pdf", width = 18, height = 3.5)
VlnPlot(SO_human2, c("COLEC10","GFAP","LRAT","TCF21","HAND2"),
        pt.size = 0, group.by = "celltype", ncol = 5, idents = a)& NoLegend()&
  xlab("")  &
  theme(legend.position = "none", 
        axis.text.x = element_text(size = rel(1), angle = 45), 
        axis.title.y = element_text(size = rel(1), angle = 90), 
        axis.text.y = element_text(size = rel(1)) ) 
dev.off()
FeaturePlot(SO_human2, "RELN", raster = F, order = F)

a <- intersect(unique(SO_mouse2$celltype), unique(SO_human2$celltype))


