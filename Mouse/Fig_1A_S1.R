### Highlighting cells per dataset ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)

meta.data <- mes_atlas.sct[[]]

counts <- group_by(meta.data, study, celltype) %>% summarise(count = n())
totals <- group_by(counts, count) %>% summarise(count = n())

b <- as.data.frame(table(counts$study))

totals <- counts %>% 
  group_by(study) %>% 
  summarise(Frequency = sum(count))

counts$total <- rep(totals$Frequency,times = b$Freq)
counts$percent <- (counts$count/counts$total) * 100
ggplot(counts, aes(study, counts$percent, fill = celltype)) +
  geom_bar(stat = 'identity') + theme_classic() + theme(axis.text.x=element_text(angle=0,hjust=0.7,vjust=0.85)) +
  theme(axis.text = element_text(size = 8))




load("~/HSC activation/Datasets/Affo/affo_yapakt.integrated.RData")
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("Cancer"))
affo_yapakt.integrated <- RenameCells(affo_yapakt.integrated, add.cell.id = "Cancer")
for (i in unique(mes_atlas.sct$celltype)) {
  Idents(mes_atlas.sct) <- "celltype"
  .studycelltype <- WhichCells(mes_atlas.sct, cells = instudy, idents = i)
  assign(paste("Cancer",i, sep = "_"),.studycelltype)
}


Cancer_orig <- DimPlot(affo_yapakt.integrated, cells.highlight= list(`Cancer_Endothelial cells`,
                                                      Cancer_HSC,
                                                      Cancer_VSMC,
                                                      `Cancer_Portal fibroblasts`,
                                                      `Cancer_Mesothelial cells`,
                                                      `Cancer_Proliferating cells`),
                       , raster = F,
        cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3)+ NoLegend() +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
Cancer_orig <- rasterize(Cancer_orig, layers='Point', dpi=300)


Cancer_pdgfrb <- FeaturePlot(affo_yapakt.integrated, "Pdgfrb", raster = F) + scale_color_viridis_c() + NoAxes()+ ggtitle("")
Cancer_pdgfrb <- rasterize(Cancer_pdgfrb, layers='Point', dpi=300)

FeaturePlot(affo_yapakt.integrated, "Msln",, raster = F) + scale_color_viridis_c() + NoAxes()

Cancer_integrated <- DimPlot(mes_atlas.sct,
                             cells.highlight= list(`Cancer_Endothelial cells`,
                                                      Cancer_HSC,
                                                      Cancer_VSMC,
                                                      `Cancer_Portal fibroblasts`,
                                                      `Cancer_Mesothelial cells`,
                                                      `Cancer_Proliferating cells`),
        cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3, raster = F)+ NoLegend() +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
Cancer_integrated <- rasterize(Cancer_integrated, layers='Point', dpi=300)

### KoloAPAP --------- 
load("~/HSC activation/Datasets/Kolodziejczyk/spf.RData")
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("APAP"))
colnames(spf)
spf <- RenameCells(spf, add.cell.id = "Apapkolo")


for (i in unique(mes_atlas.sct$celltype)[!unique(mes_atlas.sct$celltype) %in% "VSMC"]) {
  Idents(mes_atlas.sct) <- "celltype"
  .studycelltype <- WhichCells(mes_atlas.sct, cells = instudy, idents = i)
  assign(paste("APAP",i, sep = "_"),.studycelltype)
}


APAP_orig <- DimPlot(spf, cells.highlight= list(`APAP_Endothelial cells`,
                                                      APAP_HSC,
                                                      `APAP_Portal fibroblasts`,
                                                      `APAP_Mesothelial cells`,
                                                      `APAP_Proliferating cells`),
                     raster = F,
        cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3)+ NoLegend() + NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
APAP_orig <- rasterize(APAP_orig, layers='Point', dpi=300)

APAP_pdgfrb <- FeaturePlot(spf, "Pdgfrb", raster = F) + scale_color_viridis_c() + NoAxes()+ ggtitle("")
APAP_pdgfrb <- rasterize(APAP_pdgfrb, layers='Point', dpi=300)

FeaturePlot(spf, "Msln") + scale_color_viridis_c() + NoAxes()


APAP_integrated <- DimPlot(mes_atlas.sct, cells.highlight= list(`APAP_Endothelial cells`,
                                   APAP_HSC,
                                   `APAP_Portal fibroblasts`,
                                   `APAP_Mesothelial cells`,
                                   `APAP_Proliferating cells`),
        cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3, raster = F)+ NoLegend() + NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
APAP_integrated <- rasterize(APAP_integrated, layers='Point', dpi=300)

### KoloTAA ---------
load("~/HSC activation/Datasets/Kolodziejczyk/spf.RData")
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("TAA"))
colnames(spf)
spf <- RenameCells(spf, add.cell.id = "Kolotaa")


for (i in unique(mes_atlas.sct$celltype)[!unique(mes_atlas.sct$celltype) %in% c("VSMC", "Endothelial cells","Mesothelial cells")]) {
  Idents(mes_atlas.sct) <- "celltype"
  .studycelltype <- WhichCells(mes_atlas.sct, cells = instudy, idents = i)
  assign(paste("TAA",i, sep = "_"),.studycelltype)
}


TAA_orig <- DimPlot(spf ,cells.highlight= list(TAA_HSC,
                                             `TAA_Portal fibroblasts`,
                                             `TAA_Proliferating cells`),
        cols.highlight = c("#00bfc4","#b79f00","#f8766d"), cols= "lightgrey", 
        raster = F,
        sizes.highlight = 0.3)+ NoLegend() +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))

TAA_orig <- rasterize(TAA_orig, layers='Point', dpi=300)

TAA_pdgfrb <- FeaturePlot(spf, "Pdgfrb", raster  = T) + scale_color_viridis_c() + NoAxes()+ ggtitle("")
TAA_pdgfrb <- rasterize(TAA_pdgfrb, layers='Point', dpi=300)
FeaturePlot(spf, "Msln") + scale_color_viridis_c() + NoAxes()



TAA_integrated <- DimPlot(mes_atlas.sct ,cells.highlight= list(TAA_HSC,
                                   `TAA_Portal fibroblasts`,
                                   `TAA_Proliferating cells`),
        cols.highlight = c("#00bfc4","#b79f00","#f8766d"), cols= "lightgrey", 
        sizes.highlight = 0.3, raster = F)+ NoLegend() +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
TAA_integrated <- rasterize(TAA_integrated, layers='Point', dpi=300)

### Krenkel ---------
load("~/HSC activation/Datasets/Krenkel/krenkel_post.RData")
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("CCl4_Krenkel"))
colnames(krenkel)


for (i in unique(mes_atlas.sct$celltype)) {
  Idents(mes_atlas.sct) <- "celltype"
  .studycelltype <- WhichCells(mes_atlas.sct, cells = instudy, idents = i)
  assign(paste("Krenkel",i, sep = "_"),.studycelltype)
}


CCl4_Krenkel_orig <- DimPlot(krenkel, cells.highlight= list(`Krenkel_Endothelial cells`,
                                                      Krenkel_HSC,
                                                      Krenkel_VSMC,
                                                      `Krenkel_Portal fibroblasts`,
                                                      `Krenkel_Mesothelial cells`,
                                                      `Krenkel_Proliferating cells`),
                             raster = F,
        cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3)+ NoLegend() +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
CCl4_Krenkel_orig <- rasterize(CCl4_Krenkel_orig, layers='Point', dpi=300)

CCl4_Krenkel_pdgfrb <- FeaturePlot(krenkel, "Pdgfrb", raster = F) + scale_color_viridis_c() + NoAxes()+ ggtitle("")
CCl4_Krenkel_pdgfrb <- rasterize(CCl4_Krenkel_pdgfrb, layers='Point', dpi=300)

FeaturePlot(krenkel, "Msln") + scale_color_viridis_c() + NoAxes()


CCl4_Krenkel_integrated <- DimPlot(mes_atlas.sct, cells.highlight= list(`Krenkel_Endothelial cells`,
                                       Krenkel_HSC,
                                       Krenkel_VSMC,
                                       `Krenkel_Portal fibroblasts`,
                                       `Krenkel_Mesothelial cells`,
                                       `Krenkel_Proliferating cells`),
        cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3, raster = F)+ NoLegend() +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
CCl4_Krenkel_integrated <- rasterize(CCl4_Krenkel_integrated, layers='Point', dpi=300)

### Rosenthal ---------
load("~/HSC activation/Datasets/Rosenthal/pbmc.RData")
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("WD"))
colnames(pbmc)
pbmc <- RenameCells(pbmc, add.cell.id = "Rosenthal")


for (i in unique(mes_atlas.sct$celltype)) {
  Idents(mes_atlas.sct) <- "celltype"
  .studycelltype <- WhichCells(mes_atlas.sct, cells = instudy, idents = i)
  assign(paste("Rosenthal",i, sep = "_"),.studycelltype)
}


WD_orig <- DimPlot(pbmc, cells.highlight= list(`Rosenthal_Endothelial cells`,
                                       Rosenthal_HSC,
                                       Rosenthal_VSMC,
                                       `Rosenthal_Portal fibroblasts`,
                                       `Rosenthal_Mesothelial cells`,
                                       `Rosenthal_Proliferating cells`),
                   raster = F,
        cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3)+ NoLegend() +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
WD_orig <- rasterize(WD_orig, layers='Point', dpi=300)

WD_pdgfrb <- FeaturePlot(pbmc, "Pdgfrb", raster  = T) + scale_color_viridis_c() + NoAxes()+ ggtitle("")
WD_pdgfrb <- rasterize(WD_pdgfrb, layers='Point', dpi=300)

FeaturePlot(pbmc, "Stmn1") + scale_color_viridis_c() + NoAxes()
FeaturePlot(pbmc, "Stab2") + scale_color_viridis_c() + NoAxes()


WD_integrated <- DimPlot(mes_atlas.sct, cells.highlight= list(`Rosenthal_Endothelial cells`,
                                    Rosenthal_HSC,
                                    Rosenthal_VSMC,
                                    `Rosenthal_Portal fibroblasts`,
                                    `Rosenthal_Mesothelial cells`,
                                    `Rosenthal_Proliferating cells`),
        cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3, raster = F)+ NoLegend() +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
WD_integrated <- rasterize(WD_integrated, layers='Point', dpi=300)

### Terkelsen ---------
load("~/HSC activation/Datasets/Terkelsen/terkelsen0.RData")
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("CCl4_Terkelsen"))
terkelsen0 <- RenameCells(terkelsen0, add.cell.id = "Terkelsen")


for (i in unique(mes_atlas.sct$celltype)) {
  Idents(mes_atlas.sct) <- "celltype"
  .studycelltype <- WhichCells(mes_atlas.sct, cells = instudy, idents = i)
  assign(paste("Terkelsen",i, sep = "_"),.studycelltype)
}


CCl4_Terkelsen_orig <- DimPlot(terkelsen0, cells.highlight= list(`Terkelsen_Endothelial cells`,
                                    Terkelsen_HSC,
                                    Terkelsen_VSMC,
                                    `Terkelsen_Portal fibroblasts`,
                                    `Terkelsen_Mesothelial cells`,
                                    `Terkelsen_Proliferating cells`),
                               raster = F,
        cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3)+ NoLegend()  +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
CCl4_Terkelsen_orig <- rasterize(CCl4_Terkelsen_orig, layers='Point', dpi=300)

CCl4_Terkelsen_pdgfrb <- FeaturePlot(terkelsen0, "Pdgfrb", raster = F) + scale_color_viridis_c() + NoAxes()+ ggtitle("")
CCl4_Terkelsen_pdgfrb <- rasterize(CCl4_Terkelsen_pdgfrb, layers='Point', dpi=300)

FeaturePlot(terkelsen0, "Stmn1") + scale_color_viridis_c() + NoAxes()
FeaturePlot(terkelsen0, "Stab2") + scale_color_viridis_c() + NoAxes()


CCl4_Terkelsen_integrated <- DimPlot(mes_atlas.sct, cells.highlight= list(`Terkelsen_Endothelial cells`,
                                          Terkelsen_HSC,
                                          Terkelsen_VSMC,
                                          `Terkelsen_Portal fibroblasts`,
                                          `Terkelsen_Mesothelial cells`,
                                          `Terkelsen_Proliferating cells`),
        cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3, raster = F)+ NoLegend()  +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
CCl4_Terkelsen_integrated <- rasterize(CCl4_Terkelsen_integrated, layers='Point', dpi=300)


### Yangccl4  ---------
`GSE171904_2021-01-04_seurat` <- readRDS("~/HSC activation/Datasets/Yang/GSE171904_2021-01-04_seurat.rds")
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("CCl4_Yang"))
colnames(`GSE171904_2021-01-04_seurat`)
`GSE171904_2021-01-04_seurat` <- RenameCells(`GSE171904_2021-01-04_seurat`, add.cell.id = "YangCCl4")


for (i in unique(mes_atlas.sct$celltype)) {
  Idents(mes_atlas.sct) <- "celltype"
  .studycelltype <- WhichCells(mes_atlas.sct, cells = instudy, idents = i)
  assign(paste("Yang",i, sep = "_"),.studycelltype)
}


CCl4_Yang_orig <- DimPlot(`GSE171904_2021-01-04_seurat`, cells.highlight= list(`Yang_Endothelial cells`,
                                          Yang_HSC,
                                          Yang_VSMC,
                                          `Yang_Portal fibroblasts`,
                                          `Yang_Mesothelial cells`,
                                          `Yang_Proliferating cells`),
                          raster = F,
        cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3)+ NoLegend() +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
CCl4_Yang_orig <- rasterize(CCl4_Yang_orig, layers='Point', dpi=300)


CCl4_Yang_pdgfrb <- FeaturePlot(`GSE171904_2021-01-04_seurat`, "Pdgfrb", raster = F) + scale_color_viridis_c() + NoAxes() + ggtitle("")
CCl4_Yang_pdgfrb <- rasterize(CCl4_Yang_pdgfrb, layers='Point', dpi=300)

FeaturePlot(`GSE171904_2021-01-04_seurat`, "Stmn1") + scale_color_viridis_c() + NoAxes()
FeaturePlot(`GSE171904_2021-01-04_seurat`, "Stab2") + scale_color_viridis_c() + NoAxes()


CCl4_Yang_integrated <- DimPlot(mes_atlas.sct, cells.highlight= list(`Yang_Endothelial cells`,
                                                             Yang_HSC,
                                                             Yang_VSMC,
                                                             `Yang_Portal fibroblasts`,
                                                             `Yang_Mesothelial cells`,
                                                             `Yang_Proliferating cells`),
        cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3, raster = F)+ NoLegend() +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
CCl4_Yang_integrated <- rasterize(CCl4_Yang_integrated, layers='Point', dpi=300)

### YangBDL  ---------
`GSE171904_2021-01-04_seurat` <- readRDS("~/HSC activation/Datasets/Yang/GSE171904_2021-01-04_seurat.rds")
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("BDL"))
colnames(`GSE171904_2021-01-04_seurat`)
`GSE171904_2021-01-04_seurat` <- RenameCells(`GSE171904_2021-01-04_seurat`, add.cell.id = "BdlYang")


for (i in unique(mes_atlas.sct$celltype)) {
  Idents(mes_atlas.sct) <- "celltype"
  .studycelltype <- WhichCells(mes_atlas.sct, cells = instudy, idents = i)
  assign(paste("BDL",i, sep = "_"),.studycelltype)
}


BDL_orig <- DimPlot(`GSE171904_2021-01-04_seurat`, cells.highlight= list(`BDL_Endothelial cells`,
                                                             BDL_HSC,
                                                             BDL_VSMC,
                                                             `BDL_Portal fibroblasts`,
                                                             `BDL_Mesothelial cells`,
                                                             `BDL_Proliferating cells`),
                    raster = F,
        cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3)+ NoLegend()  +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
BDL_orig <- rasterize(BDL_orig, layers='Point', dpi=300)

BDL_pdgfrb <- FeaturePlot(`GSE171904_2021-01-04_seurat`, "Pdgfrb", raster = F) + scale_color_viridis_c() + NoAxes() + ggtitle("")
BDL_pdgfrb <- rasterize(BDL_pdgfrb, layers='Point', dpi=300)

FeaturePlot(`GSE171904_2021-01-04_seurat`, "Stmn1") + scale_color_viridis_c() + NoAxes()
FeaturePlot(`GSE171904_2021-01-04_seurat`, "Stab2") + scale_color_viridis_c() + NoAxes()



BDL_integrated <- DimPlot(mes_atlas.sct, cells.highlight= list(`BDL_Endothelial cells`,
                                                             BDL_HSC,
                                                             BDL_VSMC,
                                                             `BDL_Portal fibroblasts`,
                                                             `BDL_Mesothelial cells`,
                                                             `BDL_Proliferating cells`),
        cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3, raster = F)+ NoLegend()  +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
BDL_integrated <- rasterize(BDL_integrated, layers='Point', dpi=300)

### Dobieccl4  ---------
load("~/HSC activation/Datasets/Dobie/dobie.RData")
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("CCl4_Dobie"))
colnames(dobie)
dobie <- RenameCells(dobie, add.cell.id = "Chronicdobie")


for (i in unique(mes_atlas.sct$celltype)) {
  Idents(mes_atlas.sct) <- "celltype"
  .studycelltype <- WhichCells(mes_atlas.sct, cells = instudy, idents = i)
  assign(paste("Dob",i, sep = "_"),.studycelltype)
}


CCl4_Dobie_orig <- DimPlot(dobie, cells.highlight= list(`Dob_Endothelial cells`,
                                                             Dob_HSC,
                                                             Dob_VSMC,
                                                             `Dob_Portal fibroblasts`,
                                                             `Dob_Mesothelial cells`,
                                                             `Dob_Proliferating cells`),
                           raster = F,
        cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3)+ NoLegend() +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
CCl4_Dobie_orig <- rasterize(CCl4_Dobie_orig, layers='Point', dpi=300)

CCl4_Dobie_pdgfrb <- FeaturePlot(dobie, "Pdgfrb", raster = F) + scale_color_viridis_c() + NoAxes() + ggtitle("")
CCl4_Dobie_pdgfrb <- rasterize(CCl4_Dobie_pdgfrb, layers='Point', dpi=300)




CCl4_Dobie_integrated <- DimPlot(mes_atlas.sct, cells.highlight= list(`Dob_Endothelial cells`,
                                     Dob_HSC,
                                     Dob_VSMC,
                                     `Dob_Portal fibroblasts`,
                                     `Dob_Mesothelial cells`,
                                     `Dob_Proliferating cells`),
        cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3, raster = F)+ NoLegend() +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
CCl4_Dobie_integrated <- rasterize(CCl4_Dobie_integrated, layers='Point', dpi=300)

### Dobieacute  ---------
load("~/HSC activation/Datasets/Dobie/dobie.RData")
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("CCl4_acute"))
colnames(dobie)
dobie <- RenameCells(dobie, add.cell.id = "Dobieacute")



for (i in unique(mes_atlas.sct$celltype)[!unique(mes_atlas.sct$celltype) %in% c("Mesothelial cells")]) {
  Idents(mes_atlas.sct) <- "celltype"
  .studycelltype <- WhichCells(mes_atlas.sct, cells = instudy, idents = i)
  assign(paste("Acu",i, sep = "_"),.studycelltype)
}


CCl4_acute_orig <- DimPlot(dobie, cells.highlight= list(`Acu_Endothelial cells`,
                                     Acu_HSC,
                                     Acu_VSMC,
                                     `Acu_Portal fibroblasts`,
                                     `Acu_Proliferating cells`),
                           raster = F,
        cols.highlight = c("#00bfc4","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3)+ NoLegend() +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
CCl4_acute_orig <- rasterize(CCl4_acute_orig, layers='Point', dpi=300)

CCl4_acute_pdgfrb <- FeaturePlot(dobie, "Pdgfrb", raster = F) + scale_color_viridis_c() + NoAxes()+ ggtitle("")
CCl4_acute_pdgfrb <- rasterize(CCl4_acute_pdgfrb, layers='Point', dpi=300)




CCl4_acute_integrated <- DimPlot(mes_atlas.sct, cells.highlight= list(`Acu_Endothelial cells`,
                                     Acu_HSC,
                                     Acu_VSMC,
                                     `Acu_Portal fibroblasts`,
                                     `Acu_Proliferating cells`),
        cols.highlight = c("#00bfc4","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
        sizes.highlight = 0.3, raster = F)+ NoLegend() +NoAxes() + ggtitle("") +theme(plot.title = element_text(hjust = 0.5))
CCl4_acute_integrated <- rasterize(CCl4_acute_integrated, layers='Point', dpi=300)




for (i in unique(mes_atlas.sct$study)) {
  Idents(mes_atlas.sct) <- "study"
  counts_plot <- VlnPlot(
    mes_atlas.sct,
    idents = i,
    features =  c("nCount_RNA"),
    pt.size = 0
  ) &
    theme(axis.title.x = element_blank()) & labs(title = "")  &  theme(legend.position = 'none')&
    xlab("")  &
    theme(legend.position = "none", 
          axis.text.x = element_text(size = rel(1), angle = 0), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)) )
  
  features_plot <- VlnPlot(
    mes_atlas.sct,
    idents = i,
    features =  c("nFeature_RNA"),
    pt.size = 0
  ) &
    theme(axis.title.x = element_blank()) & labs(title = "") &  theme(legend.position = 'none')&
    xlab("")  &
    theme(legend.position = "none", 
          axis.text.x = element_text(size = rel(1), angle = 0), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)) )
  
  mito_plot <- VlnPlot(
    mes_atlas.sct,
    idents = i,
    features =  c("percent.mito"),
    pt.size = 0
  ) &
    theme(axis.title.x = element_blank()) & labs(title = "")&  theme(legend.position = 'none')&
    xlab("")  &
    theme(legend.position = "none", 
          axis.text.x = element_text(size = rel(1), angle = 0), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)) )
  
  assign(paste0("QC_counts_",i), counts_plot)
  assign(paste0("QC_features_",i), features_plot)
  assign(paste0("QC_mito_",i), mito_plot)
  
}



pdf(file = "~/HSC activation/Paper/Suppfigures/S1/Supp1.pdf", width = 9, height = 16)


for (i in unique(mes_atlas.sct$study)) {
  print(i)
  .qc <- ggarrange(get(paste0(i,"_orig")),get(paste0(i,"_pdgfrb")),
               get(paste0("QC_counts_",i)),get(paste0("QC_features_",i)), get(paste0("QC_mito_",i)),
               get(paste0(i,"_integrated")), nrow = 1, widths = c(2,2.75,1.75,1.75,1.75,2))
  assign(paste0("QC_",i), .qc)
}
grid.arrange(QC_APAP, QC_BDL, QC_Cancer, QC_CCl4_acute, QC_CCl4_Dobie, QC_CCl4_Krenkel, QC_CCl4_Terkelsen, QC_CCl4_Yang, QC_TAA, QC_WD, ncol = 1)


dev.off()


for (i in unique(mes_atlas.sct$study)) {
  .qc <- grid.arrange(get(paste0(i,"_pdgfrb")))
}

  
pdf(file = "~/HSC activation/Paper/figures/overview_origs.pdf", height = 4, width = 9)
grid.arrange(APAP_orig,   CCl4_acute_orig,      CCl4_Yang_orig,     CCl4_Krenkel_orig,     Cancer_orig,           
             TAA_orig,  CCl4_Dobie_orig, BDL_orig,            WD_orig, CCl4_Terkelsen_orig, nrow = 2)
dev.off()

pdf(file = "~/HSC activation/Paper/figures/mes_atlas_umap.pdf", height = 4, width = 5)
DimPlot(mes_atlas.sct, group.by = "celltype", ) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

pdf(file = "~/HSC activation/Paper/figures/mes_atlas_umap.pdf", height = 4, width = 5)
DimPlot(mes_atlas.sct, group.by = "celltype", split.by = "study", raster = F) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

