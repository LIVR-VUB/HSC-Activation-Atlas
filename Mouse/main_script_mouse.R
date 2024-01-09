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


### Loading datasets and additional QC ----
load("~/HSC activation/Datasets/Dobie/dobieacute.RData")
load("~/HSC activation/Datasets/Dobie/dobiechronic.RData")
load("~/HSC activation/Datasets/Kolodziejczyk/kolo_apap.RData")
load("~/HSC activation/Datasets/Kolodziejczyk/kolo_aat.RData")
load("~/HSC activation/Datasets/Krenkel/krenkel.RData")
load("~/HSC activation/Datasets/Rosenthal/rosenthal.RData")
load("~/HSC activation/Datasets/Terkelsen/terkelsen.RData")
load("~/HSC activation/Datasets/Yang/yangbdl.RData")
load("~/HSC activation/Datasets/Yang/yangccl4.RData")
load("~/HSC activation/Datasets/Affo/affo_mes.RData")


# Quality control checks for each dataset
datasets <- c("aat","apap","dobieacute","dobiechronic","krenkel","rosenthal","terkelsen",
              "yangbdl","yangccl4","cancer")
cancer <- affo_mes

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
  plot(VlnPlot(a, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), pt.size = 0,group.by = "orig.ident"))
}

Idents(aat) <- "orig.ident"
VlnPlot(aat, "nCount_RNA", y.max = 2500) # 1750
aat <- subset(aat, subset = nCount_RNA > 1750 & nFeature_RNA > 950 )

Idents(apap) <- "orig.ident"
VlnPlot(apap, "nCount_RNA", y.max = 2500) # 1500
apap <- subset(apap, subset = nCount_RNA > 1500 & nFeature_RNA > 1000 )

Idents(apap) <- "orig.ident"
VlnPlot(apap, "nCount_RNA", y.max = 2500) # 1500
apap <- subset(apap, subset = nCount_RNA > 1500 & nFeature_RNA > 1000 )

Idents(dobieacute) <- "orig.ident"
VlnPlot(dobieacute, "nCount_RNA") 
dobieacute <- subset(dobieacute, subset = nCount_RNA < 50000)

Idents(dobiechronic) <- "orig.ident"
VlnPlot(dobiechronic, "nCount_RNA")
dobiechronic <- subset(dobiechronic, subset = nCount_RNA < 25000)

Idents(terkelsen) <- "orig.ident"
VlnPlot(terkelsen, "nCount_RNA")
terkelsen <- subset(terkelsen, subset = nCount_RNA < 12000)


Idents(rosenthal) <- "orig.ident"
VlnPlot(rosenthal, "nCount_RNA") 
rosenthal <- subset(rosenthal, subset = nCount_RNA < 48000)

cancer$condition <- "cancer"


### Creating murine mesenchymal liver atlas ----
# Only retain genes shared by each dataset
commongenes <- Reduce(intersect, list(rownames(yangccl4),rownames(krenkel),rownames(dobieacute),rownames(yangbdl),rownames(rosenthal),rownames(terkelsen), rownames(cancer), rownames(dobiechronic), rownames(aat), rownames(apap)))
yangccl4 <- subset(yangccl4, features = commongenes)
krenkel <- subset(krenkel, features = commongenes)
dobieacute <- subset(dobieacute, features = commongenes)
yangbdl <- subset(yangbdl, features = commongenes)
rosenthal <- subset(rosenthal, features = commongenes)
terkelsen <- subset(terkelsen, features = commongenes)
dobiechronic <- subset(dobiechronic, features = commongenes)
apap <- subset(apap, features = commongenes)
aat <- subset(aat, features = commongenes)
cancer <- subset(cancer, features = commongenes)


mes_atlas <- merge(krenkel, c(yangccl4, dobieacute,dobiechronic, yangbdl,rosenthal,terkelsen,aat, apap, cancer), add.cell.ids = c("Krenkel", "YangCCl4","Dobieacute","Chronicdobie","BdlYang","Rosenthal","Terkelsen","Kolotaa","Apapkolo","Cancer"))
mes_atlas$cond <- substr(colnames(mes_atlas),1,2)
mes_atlas$treatment <- paste(mes_atlas$cond, mes_atlas$condition, sep = "_")
mes_atlas$cond <- factor(mes_atlas$cond, levels = c("Do", "Ap", "Ko","Bd","Ro","Kr","Ya","Ch","Te","Ca"))

# rename dataset metadata
Idents(mes_atlas) <- "cond"
unique(Idents(mes_atlas))
new.cluster.ids <- c("CCl4_acute","APAP","TAA","BDL","WD","CCl4_Krenkel","CCl4_Yang", "CCl4_Dobie",
                     "CCl4_Terkelsen","Cancer")
names(new.cluster.ids) <- levels(mes_atlas)
mes_atlas <- RenameIdents(mes_atlas, new.cluster.ids)
mes_atlas$study<- Idents(mes_atlas)



# removal of batch effects + integration
mes_atlas.list <- SplitObject(mes_atlas, split.by = "study")
mes_atlas.list <- lapply(X = mes_atlas.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = mes_atlas.list, nfeatures = 3000)
mes_atlas.list <- PrepSCTIntegration(object.list = mes_atlas.list, anchor.features = features)
mes_atlas.list <- lapply(X = mes_atlas.list, FUN = RunPCA, features = features)

rm(krenkel)
rm(dobieacute)
rm(dobiechronic)
rm(rosenthal)
rm(yangbdl)
rm(yangccl4)
rm(terkelsen)
rm(aat)
rm(apap)
rm(mes_atlas)
rm(a)
rm(C)
rm(cancer)
rm(affo_mes)
gc()

anchors <- FindIntegrationAnchors(object.list = mes_atlas.list, normalization.method = "SCT",
                                  anchor.features = features, dims = 1:20,
                                  reduction = "rpca", k.anchor = 5)
rm(mes_atlas.list)
gc()
mes_atlas.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")



# Add percentage counts of mitochondrial and ribosomal genes as metadata
DefaultAssay(mes_atlas.sct) <- "RNA"
mt.genes <- rownames(mes_atlas.sct)[grep("^mt-",rownames(mes_atlas.sct), ignore.case = TRUE)]
C<-GetAssayData(object = mes_atlas.sct, slot = "counts")
mes_atlas.sct$percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100

rb.genes <- rownames(mes_atlas.sct)[grep("^Rp[sl][[:digit:]]",rownames(mes_atlas.sct))]
C<-GetAssayData(object = mes_atlas.sct, slot = "counts")
mes_atlas.sct$percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100


# Mesenchymal Atlas
DefaultAssay(mes_atlas.sct) <- "integrated"
mes_atlas.sct <- RunPCA(mes_atlas.sct)
mes_atlas.sct <- FindNeighbors(mes_atlas.sct, reduction = "pca", dims = 1:30)
DefaultAssay(mes_atlas.sct) <- "integrated"
mes_atlas.sct <- FindClusters(mes_atlas.sct, resolution = 0.20)
mes_atlas.sct <- RunUMAP(mes_atlas.sct, reduction = "pca", dims = 1:30)
DimPlot(mes_atlas.sct, raster = F)

# Fig.1A --> see script Fig1A_S1

# Reorder metadata
mes_atlas.sct$study <- factor(mes_atlas.sct$study, levels = c("CCl4_acute","APAP","TAA","BDL","WD","CCl4_Krenkel","CCl4_Yang", "CCl4_Dobie",
                                                    "CCl4_Terkelsen","Cancer"))


# Rename the Idents based on markers
mes_atlas.sct <- RenameIdents(mes_atlas.sct, "0" = "HSC","1" = "HSC","2" = "HSC", "3" = "HSC","4" = "Portal fibroblasts",
                         "5" = "VSMC","6" = "HSC","7" = "Proliferating cells", "8" = "Endothelial cells", "9" ="Mesothelial cells")
mes_atlas.sct$celltype<- Idents(mes_atlas.sct)
Idents(mes_atlas.sct) <- "celltype"



### Fig 1B,C: Dotplot & Violins ----
celltypemarkers <- list("HSC \n markers" = c("Lrat", "Reln", "Colec11"),
              "Portal fibroblast \n markers" = c("Eln", "Fbln2", "Entpd2", "Col15a1", "Cd34"), 
              "Vascular smooth muscle \n markers" = c("Myh11", "Cnn1"), 
              "Proliferation \n markers" = c("Stmn1", "Cdk1", "Mki67", "Mcm6"), 
              "LSEC \n markers" = c("Clec4g", "Lyve1", "Stab2"),
              "Mesothelial \n markers" = c("Upk3b","Krt7"))
DefaultAssay(mes_atlas.sct) <- "RNA"

DotPlot(mes_atlas.sct, features = celltypemarkers) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_viridis()

pdf(file = "~/HSC activation/Paper/figures/celltype_dotplot.pdf", width = 9, height = 4) 
print(DotPlot(mes_atlas.sct, features = celltypemarkers)+ theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off() 

pdf(file = "~/HSC activation/Paper/figures/violins_mes_atlas.pdf", width = 14, height = 4)
VlnPlot(mes_atlas.sct, c("nFeature_RNA","Lrat","Pecam1"), pt.size = 0, group.by = "celltype") 
dev.off()


### Fig. S2: Mesothelial cells per study ----
# Create vector with levels of object@ident
identities <- levels(mes_atlas.sct$celltype)

# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))

pdf(file = "~/HSC activation/Paper/Suppfigures/S4/upk3B.pdf", width = 7, height = 3.5)
VlnPlot(mes_atlas.sct, "Upk3b",group.by = "study", split.by = "celltype", pt.size = 0, cols = my_color_palette)&
  xlab("")  &
  theme(axis.text.x = element_text(size = rel(1), angle = 45), 
        axis.title.y = element_text(size = rel(1), angle = 90), 
        axis.text.y = element_text(size = rel(1)) ) 
dev.off()
VlnPlot(mes_atlas.sct, "Krt7",group.by = "study", split.by = "celltype", pt.size = 0, cols = my_color_palette)


### Fig 1D: Celltype composition per liver injury ----
Idents(mes_atlas.sct) <- "treatment"
unique(Idents(mes_atlas.sct))
new.cluster.ids <- c("Control","Chronic","Control","Chronic","Control","Acute","Control", "Chronic",
                     "Control","Chronic","Chronic","Chronic","Control","Control","Chronic","Chronic","Control", "Control",
                     "Control", "Acute", "Acute", "Acute","Acute","Control", "Control",
                     "Control", "Acute", "Acute","Acute","Acute","Chronic")
names(new.cluster.ids) <- levels(mes_atlas.sct)
mes_atlas.sct <- RenameIdents(mes_atlas.sct, new.cluster.ids)
mes_atlas.sct$acutevschronic<- Idents(mes_atlas.sct)

mes_atlas.sct <- RenameIdents(mes_atlas.sct,"Control" = "Control","Acute" = "Treatment","Chronic" = "Treatment")
mes_atlas.sct$treated <- Idents(mes_atlas.sct)


mes_atlas.treated <- subset(mes_atlas.sct, idents = "Treatment")
mes_atlas.control <- subset(mes_atlas.sct, idents = "Control")

# meta.data <- mes_atlas.sct[[]]
meta.data <- mes_atlas.treated[[]]
# meta.data <- mes_atlas.control[[]]


counts <- group_by(meta.data, study, celltype) %>% summarise(count = n())
totals <- counts %>% 
  group_by(study) %>% 
  summarise(Frequency = sum(count))

b <- as.data.frame(table(counts$study))


counts$total <- rep(totals$Frequency,times = b$Freq)
counts$percent <- (counts$count/counts$total) * 100

pdf(file = "~/HSC activation/Paper/figures/celltype_barplot2.pdf", width = 8, height = 2)
print(ggplot(counts, aes(study, counts$percent, fill = celltype)) +
        geom_bar(stat = 'identity') + theme_classic() + theme(axis.text.x=element_text(angle=0,hjust=0.7,vjust=0.85)) +
        theme(axis.text = element_text(size = 8))
)
dev.off()


### Creating HSC Atlas ----

# subsetting HSCs
Idents(mes_atlas.sct) <- "celltype"
sct.hscs <- subset(mes_atlas.sct, idents = c("HSC"))

# Removing ribosomal and mitochondrial genes from the top variable genes
mt.genes <- rownames(sct.hscs)[grep("^mt-",rownames(sct.hscs), ignore.case = TRUE)]
rb.genes <- rownames(sct.hscs)[grep("^Rp[sl][[:digit:]]",rownames(sct.hscs))]

DefaultAssay(sct.hscs) <- "integrated"
sct.hscs@assays$integrated@var.features
str(sct.hscs@assays$integrated@var.features)
sct.hscs@assays$integrated@var.features <- sct.hscs@assays$integrated@var.features[!sct.hscs@assays$integrated@var.features %in% c(mt.genes,rb.genes)]
str(sct.hscs@assays$integrated@var.features)


# Dimensional reduction
DefaultAssay(sct.hscs) <- "integrated"
sct.hscs <- RunPCA(sct.hscs)
ElbowPlot(sct.hscs, ndims = 50)
sct.hscs <- FindNeighbors(sct.hscs, dims = 1:20)
sct.hscs <- FindClusters(sct.hscs, resolution = 0.13)
sct.hscs <- RunUMAP(sct.hscs, dims = 1:20)

DimPlot(sct.hscs)
DefaultAssay(sct.hscs) <- "RNA"



### Determine HSC subtypes ----
Idents(sct.hscs) <- "seurat_clusters"
DimPlot(sct.hscs)
for (i in unique(sct.hscs$seurat_clusters)) {
  .subtypemarker <- FindMarkers(sct.hscs,
                                ident.1 = i, 
                                logfc.threshold = 1,
                                only.pos = TRUE,
                                min.pct = 0)
  assign(paste("markers_cl",i,sep = ""), .subtypemarker)
}

# Rename
Idents(sct.hscs) <- "seurat_clusters"
sct.hscs <- RenameIdents(sct.hscs, "0" = "inactivating","1" = "quiescent",
                         "2" = "myofibroblasts", "3" = "initiatory")
sct.hscs$subtype<- Idents(sct.hscs)
sct.hscs$subtype <- factor(sct.hscs$subtype, levels = c("quiescent","inactivating","initiatory","myofibroblasts"))
Idents(sct.hscs) <- "subtype"

DimPlot(sct.hscs, group.by = "subtype", label = T) + NoLegend()

# HSC atlas without inactivating cluster
Idents(sct.hscs) <- "subtype"
hscs.filtered <- subset(sct.hscs,
                          idents = c("quiescent","initiatory","myofibroblasts"))
Idents(hscs.filtered) <- "condition"
hscs.filtered <- subset(hscs.filtered,
                        idents = c("8w Regression"),
                        invert = TRUE)

# Fig 2: Dimensional reduction ----
DefaultAssay(hscs.filtered) <- "integrated"
str(hscs.filtered@assays$integrated@var.features)
hscs.filtered <- RunPCA(hscs.filtered)

ElbowPlot(hscs.filtered, ndims = 50)
hscs.filtered <- FindNeighbors(hscs.filtered, dims = 1:20)
hscs.filtered <- RunUMAP(hscs.filtered, dims = 1:20)

DimPlot(hscs.filtered, group.by = "subtype", label = TRUE) + NoLegend() + NoAxes() + ggtitle("HSC atlas")
DimPlot(hscs.filtered, group.by = "subtype", split.by = "study")
DimPlot(hscs.filtered, group.by = "subtype", reduction = "pca")

pdf(file = "~/HSC activation/Paper/figures/Fig_2/Dimplot.pdf", width = 5, height = 5)
dimplt_hsc <- DimPlot(hscs.filtered, group.by = "subtype", label = F, raster = F) + NoLegend() + NoAxes() +ggtitle("")
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off()


### Fig 2C: Location of Healthy, Acute or Chronic conditions on UMAP ----
Idents(hscs.filtered) <- "acutevschronic"
cellsofinterest <- WhichCells(hscs.filtered, idents = c("Control"))

for (i in unique(hscs.filtered$subtype)) {
  Idents(hscs.filtered) <- "subtype"
  .studycelltype <- WhichCells(hscs.filtered, cells = cellsofinterest, idents = i)
  assign(paste("Control",i, sep = "_"),.studycelltype)
}


DimPlot(hscs.filtered, cells.highlight= list(Control_quiescent,
                                               Control_initiatory,
                                               Control_myofibroblasts),
        cols.highlight = c("#619cff","#00ba38","#f8786f"), cols= "lightgrey", 
        sizes.highlight = 0.3)+ NoLegend()  +NoAxes() + ggtitle("Healthy") +theme(plot.title = element_text(hjust = 0.5))

Idents(hscs.filtered) <- "acutevschronic"
cellsofinterest <- WhichCells(hscs.filtered, idents = c("Acute"))

for (i in unique(hscs.filtered$subtype)) {
  Idents(hscs.filtered) <- "subtype"
  .studycelltype <- WhichCells(hscs.filtered, cells = cellsofinterest, idents = i)
  assign(paste("Acute",i, sep = "_"),.studycelltype)
}


DimPlot(hscs.filtered, cells.highlight= list(Acute_quiescent,
                                               Acute_initiatory,
                                               Acute_myofibroblasts),
        cols.highlight = c("#619cff","#00ba38","#f8786f"), cols= "lightgrey", 
        sizes.highlight = 0.3)+ NoLegend()  +NoAxes() + ggtitle("Acute injury") +theme(plot.title = element_text(hjust = 0.5))

Idents(hscs.filtered) <- "acutevschronic"
cellsofinterest <- WhichCells(hscs.filtered, idents = c("Chronic"))

for (i in unique(hscs.filtered$subtype)) {
  Idents(hscs.filtered) <- "subtype"
  .studycelltype <- WhichCells(hscs.filtered, cells = cellsofinterest, idents = i)
  assign(paste("Chronic",i, sep = "_"),.studycelltype)
}


DimPlot(hscs.filtered, cells.highlight= list(Chronic_quiescent,
                                               Chronic_initiatory,
                                               Chronic_myofibroblasts),
        cols.highlight = c("#619cff","#00ba38","#f8786f"), cols= "lightgrey", 
        sizes.highlight = 0.3)+ NoLegend()  +NoAxes() + ggtitle("Chronic injury") +theme(plot.title = element_text(hjust = 0.5))


### Define subtype markers ----
for (i in unique(hscs.filtered$subtype)) {
  Idents(hscs.filtered) <- "subtype"
  .generalmarkers <- FindMarkers(hscs.filtered,
                                 ident.1 = i,
                                 assay = "RNA",
                                 min.pct = 0,
                                 only.pos = TRUE,
                                 logfc.threshold = 1)
  tocompare <- unique(hscs.filtered$subtype)[!unique(hscs.filtered$subtype) %in% i ]
  .specific1 <- FindMarkers(hscs.filtered,
                            ident.1 = i,
                            ident.2 = tocompare[1],
                            min.pct = 0,
                            assay = "RNA",
                            only.pos = TRUE,)
  .specific2 <- FindMarkers(hscs.filtered,
                            ident.1 = i,
                            ident.2 = tocompare[2],
                            min.pct = 0,
                            assay = "RNA",
                            only.pos = TRUE,)
  .generalmarkers <- .generalmarkers[rownames(.generalmarkers) %in% intersect(rownames(.specific1),rownames(.specific2)),]
  assign(paste(i,"geneset", sep = "_"), .generalmarkers)
  
}


### Fig. 2B: Heatmap of subtype markers ----
toplot <- c(rownames(quiescent_geneset), rownames(initiatory_geneset), rownames(myofibroblasts_geneset))

mat<- hscs.filtered[["RNA"]]@data[toplot, ] %>% as.matrix()
mat<- t(scale(t(mat)))
cluster_anno <- hscs.filtered@meta.data$subtype
quantile(mat, c(0.1, 0.95))

col_fun = circlize::colorRamp2(c(-1.5, 0, 2), rev(c("#d73027","#ffffbf","#4575b4")))

markers <- c("Rgs5","Lrat","Reln","Ets1","Ntm","Ngfr",
             "Ccl2","Ankrd1","Timp1","Ccl7",
             "Col1a1","Loxl1","Mmp2","Spp1", "Col3a1","Lox")

ra <- rowAnnotation(foo = anno_mark(at = match(markers,rownames(mat)),
                                    labels = markers))

HSC_heatmap <- Heatmap(mat, name = "Expression",  
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

HSC_heatmap



pdf(file = "~/HSC activation/Paper/figures/Fig_2/Heatmap.pdf", width = 5.5, height = 4)
plot(HSC_heatmap)
dev.off()

### subtype markers table
quiescent_df <- data.frame("Gene" = rownames(quiescent_geneset),"log2FC" = quiescent_geneset$avg_log2FC)
quiescent_df[,"HSC stage"] <- "quiescent"

initiatory_df <- data.frame("Gene" = rownames(initiatory_geneset),"log2FC" = initiatory_geneset$avg_log2FC)
initiatory_df[,"HSC stage"] <- "initiatory"

myofibroblasts_df <- data.frame("Gene" = rownames(myofibroblasts_geneset),"log2FC" = myofibroblasts_geneset$avg_log2FC)
myofibroblasts_df[,"HSC stage"] <- "myofibroblasts"
overview_df <- rbind(quiescent_df,initiatory_df, myofibroblasts_df)
overview_df <- overview_df[,c(3,1,2)]
write.csv(overview_df,file= "~/HSC activation/Paper/Excel_magic/Gene_overview.csv", row.names = F)


### Fig. S4: Heatmaps per study of representative genes ----
data_hm <- AverageExpression(hscs.filtered, assays = "RNA", add.ident = "study")
df <- as.data.frame(data_hm$RNA)


subtype_annotation <- c("quiescent","initiatory","myofibroblasts")
subtype_annotation <- factor(subtype_annotation, levels = c("quiescent","initiatory","myofibroblasts"))
markers2 <- c("Lrat","Reln","Rgs5", "Ecm1","Pth1r","Vipr1","Calcrl","Agtr1a", "Adra2b", "Ednrb",
              "Colec10","Colec11","Fcna","Ets1",
              "Ankrd1","Thbs1","Ccl2","Ccl7","Cxcl10","Mt1","Mt2","Timp1",
              "Col1a1","Col1a2","Col3a1","Col5a2","Col6a2","Col15a1","Lox","Loxl1","Mmp2","S100a6","Gas6","Spp1")

cols = colorRamp2(c(-1.2, 0, 1.2), c("#fdfed5", "#72c8bd","#253997"))

for (i in unique(hscs.filtered$study)) {
  .study_hm <- dplyr::select(df, matches(i))
  .study_hm <- .study_hm[markers2,]
  .study_hm <- t(scale(t(.study_hm)))
  .study_ha <- HeatmapAnnotation("subtype" = subtype_annotation,
                                 col = list(subtype = c("myofibroblasts" = "#619cff", "initiatory" = "#00ba38","quiescent" = "#f8766c")))
  if (i == "Cancer") {
    assign(paste("hm",i,sep = "_"), Heatmap(.study_hm, col = cols,
                                            row_title = "",
                                            column_title = i,
                                            cluster_columns = F,
                                            cluster_rows = F,
                                            top_annotation = .study_ha,
                                            column_labels = c("","","")))
  }else {
    assign(paste("hm",i,sep = "_"), Heatmap(.study_hm, col = cols,
                                            row_title = "",
                                            column_title = i,
                                            cluster_columns = F,
                                            cluster_rows = F,
                                            top_annotation = .study_ha,
                                            column_labels = c("","","")))
  }
}
ht_list = hm_APAP + hm_TAA + hm_CCl4_acute + hm_WD + hm_BDL + hm_CCl4_Dobie + hm_CCl4_Krenkel + hm_CCl4_Terkelsen + hm_CCl4_Yang + hm_Cancer
pdf(file = "~/HSC activation/Paper/Suppfigures/S5/Heatmaps.pdf", width = 14, height = 8)
draw(ht_list)
dev.off()
### Fig. 2C: Barplot of subtypes ratio per injury type ----
meta.data <- hscs.filtered[[]]

meta_data <- meta.data[meta.data$acutevschronic == "Control",]
counts <- group_by(meta_data, subtype) %>% summarise(count = n())

total <- sum(counts$count)
counts$percents <- (counts$count/total)

ggplot(counts, aes(x= 1, counts$percents, fill = subtype)) +
  geom_bar(stat = 'identity') + theme_classic() + theme(axis.text.x=element_text(angle=0,hjust=0.7,vjust=0.85))

meta_data <- meta.data[meta.data$acutevschronic == "Acute",]
counts <- group_by(meta_data, subtype) %>% summarise(count = n())

total <- sum(counts$count)
counts$percents <- (counts$count/total)

ggplot(counts, aes(x= 1, counts$percents, fill = subtype)) +
  geom_bar(stat = 'identity') + theme_classic() + theme(axis.text.x=element_text(angle=0,hjust=0.7,vjust=0.85))

meta_data <- meta.data[meta.data$acutevschronic == "Chronic",]
counts <- group_by(meta_data, subtype) %>% summarise(count = n())

total <- sum(counts$count)
counts$percents <- (counts$count/total)

ggplot(counts, aes(x= 1, counts$percents, fill = subtype)) +
  geom_bar(stat = 'identity') + theme_classic() + theme(axis.text.x=element_text(angle=0,hjust=0.7,vjust=0.85))


### Fig. 2D: Featureplots of representative genes ----

pdf(file = "~/HSC activation/Paper/figures/Fig_2/hsc_rgs5.pdf", width = 5.5, height = 5)
hsc_rgs5 <- FeaturePlot(hscs.filtered, "Rgs5") + scale_color_viridis() + NoAxes() + ggtitle("")
hsc_rgs5 <- rasterize(hsc_rgs5, layers='Point', dpi=300)
plot(hsc_rgs5)
dev.off()

pdf(file = "~/HSC activation/Paper/figures/Fig_2/hsc_Col1a1.pdf", width = 5.5, height = 5)
hsc_Col1a1 <- FeaturePlot(hscs.filtered, "Col1a1") + scale_color_viridis() + NoAxes() + ggtitle("")
hsc_Col1a1 <- rasterize(hsc_Col1a1, layers='Point', dpi=300)
plot(hsc_Col1a1)
dev.off()

pdf(file = "~/HSC activation/Paper/figures/Fig_2/hsc_Ccl2.pdf", width = 5.5, height = 5)
hsc_Ccl2 <- FeaturePlot(hscs.filtered, "Ccl2") + scale_color_viridis() + NoAxes() + ggtitle("")
hsc_Ccl2 <- rasterize(hsc_Ccl2, layers='Point', dpi=300)
plot(hsc_Ccl2)
dev.off()

### Fig. S5: Subtypes ratio per specific injury type ----
DimPlot(hscs.filtered)
unique(hscs.filtered$treatment)
Idents(hscs.filtered) <- "treatment"
hscs.filtered <- RenameIdents(hscs.filtered, "Kr_MFB" = "3w_CCl4","Kr_Resting HSCs" = "Healthy","Ya_oil" = "Healthy",
                              "Ya_ccl4" = "3w_CCl4","Do_Uninjured" ="Healthy","Do_Acute CCl4" ="1x_CCl4",
                              "Ch_Uninjured" = "Healthy","Ch_Peak CCl4" ="6w_CCl4","Bd_oil" ="Healthy","Bd_bdl" = "BDL",
                              "Ro_12w Western diet" = "12w_WD","Ro_control" ="Healthy","Te_w0" = "Healthy", 
                              "Te_w2" = "2w_CCl4","Te_w4"="4w_CCl4","Ko_SPF_1" = "Healthy", "Ko_SPF_2" ="Healthy",
                              "Ko_SPF_3" = "Healthy","Ko_SPF_TAA_1" = "1x_TAA","Ko_SPF_TAA_2" = "1x_TAA",
                              "Ko_SPF_TAA_3" = "1x_TAA","Ko_SPF_TAA_4" = "1x_TAA","Ap_SPF_1" = "Healthy",
                              "Ap_SPF_2" = "Healthy","Ap_SPF_3" ="Healthy", "Ap_SPF_APAP_1" = "1x_APAP",
                              "Ap_SPF_APAP_2" = "1x_APAP", "Ap_SPF_APAP_3" = "1x_APAP", "Ap_SPF_APAP_4" = "1x_APAP",
                              "Ca_cancer" = "iCCA")
hscs.filtered$injury <- Idents(hscs.filtered)
a <- levels(hscs.filtered$injury )
a[1] <- "Healthy"
a[2] <- "3w_CCl4"
hscs.filtered$injury <- factor(hscs.filtered$injury, levels = a)

pdf(file = "~/HSC activation/Paper/Suppfigures/S6/ondimplots.pdf", width = 16, height = 3)
p <- DimPlot(hscs.filtered, split.by = "study", group.by = "injury", label = T, repel = F) + NoLegend() + NoAxes()
rasterize(p, layers='Point', dpi=300)
dev.off()

for (i in unique(hscs.filtered$study)) {
  meta.data <- hscs.filtered[[]]
  meta_data <- meta.data[meta.data$study == i,]
  
  counts <- group_by(meta_data, injury, subtype) %>% summarise(count = n())
  totals <- counts %>% 
    group_by(injury) %>% 
    summarise(Frequency = sum(count))
  
  b <- as.data.frame(table(counts$injury))
  b <- b[!b$Freq == 0,]
  
  counts$total <- rep(totals$Frequency,times = b$Freq)
  counts$percent <- (counts$count/counts$total) * 100
  .p <- ggplot(counts, aes(injury, percent, fill = subtype)) +
    geom_bar(stat = 'identity')  + theme_classic() + theme(axis.text.x=element_text(angle=0,hjust=0.7,vjust=0.85)) +
    theme(axis.text = element_text(size = 8)) + labs(y = "percent") + theme(legend.position = "none")
  assign(paste0("barplot_",i), .p)
}

pdf(file = "~/HSC activation/Paper/Suppfigures/S6/barplots.pdf", width = 16, height = 3)
plot_grid(barplot_CCl4_acute, barplot_APAP, barplot_TAA, barplot_BDL, barplot_WD, barplot_CCl4_Krenkel, barplot_CCl4_Yang,
          barplot_CCl4_Dobie, barplot_CCl4_Terkelsen, barplot_Cancer, ncol = 10)
dev.off()


pdf(file = "~/HSC activation/Paper/Suppfigures/S6/dimplots.pdf", width = 16, height = 3)
p <- DimPlot(hscs.filtered, split.by = "study", group.by = "subtype", label = F, repel = T) + NoLegend() + NoAxes()
rasterize(p, layers='Point', dpi=300)
dev.off()
### Fig. 2E,F: GSVA of subtype markers in time resolved CCl4/in vitro ----

# In vivo CCl4
load("~/HSC initiation/LIVR/LIVR_Seq_acute_chronic_CCl4_HSC.RData")
norm_counts <- counts(RNA_Seq_HSC, normalized=T)
rownames(norm_counts) <- mcols(RNA_Seq_HSC)$mgi_symbol
colnames(norm_counts) <- colData(RNA_Seq_HSC)$Condition

HSC_stages <- list("quiescent" = rownames(quiescent_geneset),
                   "initiatory" = rownames(initiatory_geneset)
                   ,"myofibroblasts" = rownames(myofibroblasts_geneset))


gsva_short_tp <- gsva(norm_counts, HSC_stages)
pheatmap(gsva_short_tp, cluster_rows = F, cluster_cols = F, annotation = as.factor(colnames(gsva_short_tp)))

df <- as.data.frame(gsva_short_tp)
cond_annotation <- colnames(gsva_short_tp)
cond_annotation <- factor(cond_annotation, levels = unique(cond_annotation))

quantile(as.matrix(df))
col_fun = circlize::colorRamp2(c(-0.78, 0, 0.76), rev(c("#d73027","#ffffbf","#4575b4")))
annot <- HeatmapAnnotation("condition" = cond_annotation
                           #,col = list(subtype = c("myofibroblasts" = "#c87dff", "initiatory" = "#18bdb9","quiescent" = "#f8766d"))
)

hmval <- Heatmap(as.matrix(df), col = col_fun,
                 row_title = "",
                 column_split = factor(c(rep("Healthy",4),
                                         rep("Acute injury",11),
                                         rep("Chronic injury",12)),
                                       levels = c("Healthy","Acute injury","Chronic injury")),
                 cluster_columns = F,
                 show_row_names = F,
                 heatmap_legend_param = list(title = "Enrichment"),
                 cluster_rows = F)
pdf(file = "~/HSC activation/Paper/figures/Fig_2/invivo_enrichment.pdf", width = 8, height = 4)
draw(hmval)
dev.off()

# In vitro genesets
load("~/HSC initiation/LIVR/In vitro Vincent/20201002 2D_TL RNAseq data.RData")

norm_counts <- counts(RNA_Seq_2DTL_mHSC, normalized=T)
rownames(norm_counts) <- mcols(RNA_Seq_2DTL_mHSC)$SYM
colnames(norm_counts) <- colData(RNA_Seq_2DTL_mHSC)$Time

gsva_short_tp <- gsva(norm_counts, HSC_stages)
pheatmap(gsva_short_tp, cluster_rows = F, cluster_cols = F, annotation = as.factor(colnames(gsva_short_tp)))

df <- as.data.frame(gsva_short_tp)
cond_annotation <- colnames(gsva_short_tp)
cond_annotation <- factor(cond_annotation, levels = unique(cond_annotation))

quantile(as.matrix(df))
col_fun = circlize::colorRamp2(c(-0.74, 0, 0.69), rev(c("#d73027","#ffffbf","#4575b4")))
annot <- HeatmapAnnotation("condition" = cond_annotation
                           #,col = list(subtype = c("myofibroblasts" = "#c87dff", "initiatory" = "#18bdb9","quiescent" = "#f8766d"))
)

hmval <- Heatmap(as.matrix(df), col = col_fun,
                 row_title = "",
                 
                 # top_annotation = annot,
                 # row_split = factor(rowsplitvector, levels = c("Quiescent","Initiatory","Myofibroblast")),
                 cluster_columns = F,
                 cluster_rows = F, 
                 heatmap_legend_param = list(title = "Enrichment"),
                 show_row_names =  F)

hmval
pdf(file = "~/HSC activation/Paper/figures/Fig_2/invitro_enrichment.pdf", width = 8, height = 2.7)
draw(hmval)
dev.off()

### Fig. 3A: COLEC10 violinplot ----
pdf(file = "~/HSC activation/Paper/figures/Fig_colec/vlns.pdf", width = 4, height = 9)
VlnPlot(hscs.filtered,c("Colec10","Colec11","Fcna","Masp1"), pt.size = 0, ncol = 1)&
  theme(legend.position = "none", 
        axis.text.x = element_text(size = rel(1), angle = 0), 
        axis.title.y = element_text(size = rel(1), angle = 90), 
        axis.text.y = element_text(size = rel(1)) )
dev.off()


FeaturePlot(hscs.filtered,c("Colec10","Colec11","Fcna","Masp1"), ncol = 1)
### pySCENIC: TF activity ----
hscs.filtered <- RenameCells(hscs.filtered,new.names = str_replace_all(colnames(hscs.filtered), fixed(" "), ""))
data_to_write <- as.data.frame(as.matrix(hscs.filtered@assays$RNA@counts))
dim(data_to_write)
fwrite(x = data_to_write, file = "hscs.filtered.csv", row.names = T)
# --> jupyter notebook

# jupyter notebook -->
auc_mtx <- read_table2("~/HSC activation/Paper/data/AUC_hsc.filtered.tsv",col_names = TRUE)

auc_mtx <- as.data.frame(auc_mtx)
rownames(auc_mtx) <- auc_mtx$Cell
auc_mtx <- auc_mtx[,-c(1)]

b <- t(auc_mtx)
auc_assay <- CreateAssayObject(counts = b)
hscs.filtered[["AUC"]] <- auc_assay


# Fig. 5A: differentially active regulons ----
regulons <- read_csv("~/HSC activation/Paper/data/regulons.csv")
alltfs <- list()
for (i in rownames(hscs.filtered@assays$AUC)) {
  print(i)
  .name <- substring(i,1, nchar(i)-3)
  alltfs[[.name]] <- as.data.frame(regulons[regulons$name == i ,2])[,1]
}

lengtes <- lengths(alltfs)
names(lengtes) <- paste0(names(lengtes), "(+)")

Idents(hscs.filtered) <- "subtype"
for (i in unique(hscs.filtered$subtype)) {
  .generalmarkers <- FindMarkers(hscs.filtered,
                                 ident.1 = i,
                                 assay = "AUC",
                                 only.pos = TRUE,
                                 logfc.threshold = 0,
                                 min.pct = 0)
  tocompare <- unique(hscs.filtered$subtype)[!unique(hscs.filtered$subtype) %in% i ]
  .specific1 <- FindMarkers(hscs.filtered,
                            ident.1 = i,
                            ident.2 = tocompare[1],
                            assay = "AUC",
                            min.pct = 0,
                            only.pos = TRUE,
                            logfc.threshold = 0)
  .specific2 <- FindMarkers(hscs.filtered,
                            ident.1 = i,
                            ident.2 = tocompare[2],
                            assay = "AUC",
                            min.pct = 0,
                            only.pos = TRUE,
                            logfc.threshold = 0)
  .generalmarkers <- .generalmarkers[rownames(.generalmarkers) %in% intersect(rownames(.specific1),rownames(.specific2)),]
  .generalmarkers$regulonsize <-lengtes[match(rownames(.generalmarkers),names(lengtes))]
  
  assign(paste(i,"markers", sep = "_"), .generalmarkers)
  
}
# Heatmap
quiescent_TFs <- quiescent_markers[quiescent_markers$regulonsize > 20,]
quiescent_TFs <- quiescent_TFs[order(-quiescent_TFs$avg_log2FC),]
quiescent_TFs <- rownames(quiescent_TFs)[1:15]

initiatory_TFs <- initiatory_markers[initiatory_markers$regulonsize > 20,]
initiatory_TFs <- initiatory_TFs[order(-initiatory_TFs$avg_log2FC),]
initiatory_TFs <- rownames(initiatory_TFs)[1:15]

myofibroblasts_TFs <- myofibroblasts_markers[myofibroblasts_markers$regulonsize > 20,]
myofibroblasts_TFs <- myofibroblasts_TFs[order(-myofibroblasts_TFs$avg_log2FC),]
myofibroblasts_TFs <- rownames(myofibroblasts_TFs)[1:15]

toplot <- c(quiescent_TFs[1:3], initiatory_TFs[1:3], myofibroblasts_TFs[1:3])
toplot <- c(quiescent_TFs, initiatory_TFs, myofibroblasts_TFs)


mat<- hscs.filtered[["AUC"]]@data[toplot, ] %>% as.matrix()
mat<- t(scale(t(mat)))
cluster_anno <- hscs.filtered@meta.data$subtype
quantile(mat, c(0.05,0.9))

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#0046ff", "black","#ffb900"))

tf_heatmap <- Heatmap(mat, name = "TF activity",  
             column_split = factor(cluster_anno),
             cluster_columns = F,
             show_column_dend = FALSE,
             cluster_column_slices = F,
             column_title_gp = gpar(fontsize = 8),
             column_gap = unit(0.5, "mm"),
             cluster_rows = F,
             show_row_dend = FALSE,
             col = col_fun,
             show_row_names = T,
             row_names_gp = gpar(fontsize = 7),
             column_title_rot = 0,
             top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#f8786f","#00ba38","#619cff")))),
             show_column_names = FALSE,
             use_raster = T,
             raster_quality = 10)


pdf(file = "~/HSC activation/Paper/figures/Fig_5/Heatmap.pdf", width = 5.5, height = 7)
plot(tf_heatmap)
dev.off()



### Fig. 5C: Trajectory inference ----
DefaultAssay(hscs.filtered) <- "AUC"
diethsc <- DietSeurat(hscs.filtered, assays = "AUC",dimreducs = c("umap"))
sce <- as.SingleCellExperiment(diethsc)
DimPlot(diethsc, group.by = "subtype")

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

hscs.filtered$pseudotime <- sce$slingPseudotime_1


pdf(file = "~/HSC activation/Paper/figures/Fig_5/mousePT.pdf", width = 3, height = 2.5)
dimplt_hsc <- FeaturePlot(hscs.filtered, "pseudotime", cols = c("#f8786f","#619cff")) + NoAxes() +ggtitle("")
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off()


### Determine number of knots
set.seed(3)
icMat <- evaluateK(sce,
                   conditions = factor(colData(sce)$study),
                   nGenes = 300,
                   k = 3:7, 
                   parallel=T)
icMat2 <- na.omit(icMat)
print(icMat)
aicMat <- icMat2
k <- 3:7
init_shape <- graphics::par()$mfrow
graphics::par(mfrow = c(1, 4))
# boxplots of AIC
devs <- matrix(NA, nrow = nrow(aicMat), ncol = length(k))
for (ii in seq_len(nrow(aicMat))) {
  devs[ii, ] <- aicMat[ii, ] - mean(aicMat[ii, ])
}
graphics::boxplot(devs, ylab = "Deviation from genewise average AIC",
                  xlab = "Number of knots", xaxt = "n")
graphics::axis(1, at = seq_len(length(k)), labels = k)

# scatterplot of average AIC
plot(x = k, y = colMeans(aicMat, na.rm = TRUE), type = "b",
     ylab = "Average AIC", xlab = "Number of knots")
# scatterplot of relative AIC
plot(x = k, y = colMeans(aicMat / aicMat[, 1], na.rm = TRUE), type = "b",
     ylab = "Relative AIC", xlab = "Number of knots")
# barplot of optimal AIC for genes with at least a difference of 2.
aicRange <- apply(apply(aicMat, 1, range), 2, diff)
varID <- which(aicRange > aicDiff)
if (length(varID) > 0) {
  aicMatSub <- aicMat[varID, ]
  tab <- table(k[apply(aicMatSub, 1, which.min)])
  graphics::barplot(tab, xlab = "Number of knots",
                    ylab = "# Genes with optimal k")
}
graphics::par(mfrow = init_shape)




top3 <- c("Rxra(+)","Foxf1(+)","Nr1h4(+)",
          "Fosl1(+)","Egr3(+)","Nfkb2(+)",
          "Wt1(+)","Prrx1(+)","Mef2c(+)")
sce <- fitGAM(sce, conditions = factor(colData(sce)$study), genes = top3,
              nknots = 5)





# association test
ATres <-  associationTest(sce)


# Do TFs vary per study?
condRes <- conditionTest(sce, l2fc = log2(2), global = TRUE, pairwise = F )

condRes <- condRes[rownames(condRes) %in% top3,]
condRes <- na.omit(condRes)
condRes$padj <- p.adjust(condRes$pvalue, "fdr")

for (i in top3) {
  col = pal_igv()
  .plotje <- plotSmoothers(sce, assays(sce)$counts,
                           size = 0,
                           gene = i,
                           sample = 0,
                           alpha = 0.3,
                           border = TRUE,
                           lwd = 0.75, 
                           curvesCols = col(n = 10)
                           ) + ggtitle(i) +
    theme(legend.position="none") +  scale_color_igv() + labs(x = NULL, y = NULL)
    
  assign(paste0(i,"_plot"),.plotje)
}
pdf(file = "~/HSC activation/Paper/figures/Fig_5/injuryTFs.pdf", width = 6, height = 4.5)
plot_grid(`Rxra(+)_plot`, `Foxf1(+)_plot`, `Nr1h4(+)_plot`,
          `Fosl1(+)_plot`, `Egr3(+)_plot`, `Nfkb2(+)_plot`,
          `Wt1(+)_plot`,`Prrx1(+)_plot`,`Mef2c(+)_plot`)
dev.off()


plot_grid(`Ets1(+)_plot`, `Ets1(+)_plot`, `Tcf21(+)_plot`,
          `Fosl2(+)_plot`, `Nfkb2(+)_plot`, `Nfkb1(+)_plot`)

example <- plotSmoothers(sce, assays(sce)$counts,
                         gene = "Foxf1(+)",
                         alpha = 0.3,
                         border = TRUE,
                         curvesCols = col(n = 10)
)  +  scale_color_igv()

e <- ggpubr::get_legend(example)

pdf(file = "~/HSC activation/Paper/figures/Fig_5/injuryleg end.pdf", width = 2, height = 8)
plot(e)
dev.off()


### Fig. 5B Scatterplot TFs ----
DefaultAssay(hscs.filtered) <- "AUC"
hm_df <- FetchData(hscs.filtered, vars = top3 )
hm_df <- as.data.frame(scale(hm_df))
hm_df$pseudotime <- hscs.filtered$pseudotime


pdf(file = "~/HSC activation/Paper/figures/Fig_5/QTFs.pdf", width = 5, height = 2.52)
sc_plot <- ggscatter(hm_df,"pseudotime", top3[1:3], size = 0, merge = T, point =F, 
          add = "loess",conf.int = F) + labs(x = NULL, y = NULL) 
ggpar(sc_plot, legend = "right",
                 font.legend = c(12.5, "plain", "black"))
dev.off()

pdf(file = "~/HSC activation/Paper/figures/Fig_5/ITFS.pdf", width = 5, height = 2.52)
sc_plot <- ggscatter(hm_df,"pseudotime", top3[4:6], size = 0, merge = T, point =F, 
                     add = "loess",conf.int = F) + labs(x = NULL, y = NULL) 
ggpar(sc_plot, legend = "right",
      font.legend = c(12.5, "plain", "black"))
dev.off()

pdf(file = "~/HSC activation/Paper/figures/Fig_5/MTFS.pdf", width = 5, height = 2.52)
sc_plot <- ggscatter(hm_df,"pseudotime", top3[7:9], size = 0, merge = T, point =F, 
                     add = "loess",conf.int = F) + labs(x = NULL, y = NULL) 
ggpar(sc_plot, legend = "right",
      font.legend = c(12.5, "plain", "black"))
dev.off()

### Fig. S9: top TFs across injuries ----

data_hm <- AverageExpression(hscs.filtered, assays = "AUC", add.ident = "study")
df <- as.data.frame(data_hm$AUC)


subtype_annotation <- c("quiescent","initiatory","myofibroblasts")
subtype_annotation <- factor(subtype_annotation, levels = c("quiescent","initiatory","myofibroblasts"))

tf_list <- list()
for (i in toplot) {
  print(i)
  tf_list[[i]] <- as.data.frame(regulons[regulons$name == i ,2])[,1]
}

cols = circlize::colorRamp2(c(-1.2, 0, 1.2), c("#0046ff", "black","#ffb900"))


for (i in unique(hscs.filtered$study)) {
  .study_hm <- dplyr::select(df, matches(i))
  .study_hm <- .study_hm[toplot,]
  .study_hm <- t(scale(t(.study_hm)))
  .study_ha <- HeatmapAnnotation("subtype" = subtype_annotation,
                                 col = list(subtype = c("myofibroblasts" = "#619cff", "initiatory" = "#00ba38","quiescent" = "#f8766c")))
  .row_ha = rowAnnotation(RegulonSize = anno_barplot(lengths(tf_list)))
  if (i == "Cancer") {
    assign(paste("hm",i,sep = "_"), Heatmap(.study_hm, col = cols,
                                            row_title = "",
                                            column_title = i,
                                            cluster_columns = F,
                                            cluster_rows = F,
                                            top_annotation = .study_ha,
                                            column_labels = c("","",""),
                                            right_annotation = .row_ha))
  }else {
    assign(paste("hm",i,sep = "_"), Heatmap(.study_hm, col = cols,
                                            row_title = "",
                                            column_title = i,
                                            cluster_columns = F,
                                            cluster_rows = F,
                                            top_annotation = .study_ha,
                                            column_labels = c("","","")))
  }
}
ht_list = hm_APAP + hm_TAA + hm_CCl4_acute + hm_WD + hm_BDL + hm_CCl4_Dobie + hm_CCl4_Krenkel + hm_CCl4_Terkelsen + hm_CCl4_Yang + hm_Cancer

pdf(file = "~/HSC activation/Paper/Suppfigures/S7/Heatmaps.pdf", width = 15, height = 8)
draw(ht_list)
dev.off()


save(toplot,file = "TFs.RData")

sink("sessionInfo.txt")
sessionInfo()
sink()

