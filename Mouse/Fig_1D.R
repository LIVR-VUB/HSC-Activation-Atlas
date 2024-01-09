DimPlot(mes_atlas.sct, raster = F, group.by = "celltype")
mes_atlas.sct$treated_celltype <- paste(mes_atlas.sct$treated,  mes_atlas.sct$celltype,sep = "_")
unique(mes_atlas.sct$treated_celltype)
mes_atlas.sct$treated_celltype_study <- paste(mes_atlas.sct$study,  mes_atlas.sct$treated_celltype,sep = "_")
table(mes_atlas.sct$treated_celltype_study)



Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("Cancer"))
for (i in unique(mes_atlas.sct$celltype)) {
  Idents(mes_atlas.sct) <- "treated_celltype"
  .studycelltype <- WhichCells(mes_atlas.sct,cells = instudy,  idents = paste("Treatment",i, sep = "_"))
  assign(paste("Cancer",i, sep = "_"),.studycelltype)
}


Cancer_integrated <- DimPlot(mes_atlas.sct,
                             cells.highlight= list(`Cancer_Endothelial cells`,
                                                   Cancer_HSC,
                                                   Cancer_VSMC,
                                                   `Cancer_Portal fibroblasts`,
                                                   `Cancer_Mesothelial cells`,
                                                   `Cancer_Proliferating cells`),
                             cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
                             sizes.highlight = 0.1, raster = F)+ NoLegend() +NoAxes() + ggtitle("") 

### KoloAPAP --------- 
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("APAP"))
for (i in unique(mes_atlas.sct$celltype)[unique(mes_atlas.sct$celltype) %in% c("HSC","Proliferating cells")]) {
  Idents(mes_atlas.sct) <- "treated_celltype"
  .studycelltype <- WhichCells(mes_atlas.sct,cells = instudy,  idents = paste("Treatment",i, sep = "_"))
  assign(paste("APAP",i, sep = "_"),.studycelltype)
}


APAP_integrated <- DimPlot(mes_atlas.sct,
                             cells.highlight= list(APAP_HSC,
                                                   `APAP_Proliferating cells`),
                             cols.highlight = c("#00bfc4","#f8766d"), cols= "lightgrey", 
                             sizes.highlight = 0.1, raster = F)+ NoLegend() +NoAxes() + ggtitle("") 
### KoloTAA ---------
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("TAA"))
for (i in unique(mes_atlas.sct$celltype)[unique(mes_atlas.sct$celltype) %in% c("HSC")]) {
  Idents(mes_atlas.sct) <- "treated_celltype"
  .studycelltype <- WhichCells(mes_atlas.sct,cells = instudy,  idents = paste("Treatment",i, sep = "_"))
  assign(paste("TAA",i, sep = "_"),.studycelltype)
}


TAA_integrated <- DimPlot(mes_atlas.sct,
                           cells.highlight= list(TAA_HSC),
                           cols.highlight = c("#f8766d"), cols= "lightgrey", 
                           sizes.highlight = 0.1, raster = F)+ NoLegend() +NoAxes() + ggtitle("") 
### Krenkel ---------
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("CCl4_Krenkel"))

for (i in unique(mes_atlas.sct$celltype)) {
  Idents(mes_atlas.sct) <- "treated_celltype"
  .studycelltype <- WhichCells(mes_atlas.sct,cells = instudy,  idents = paste("Treatment",i, sep = "_"))
  assign(paste("Krenkel",i, sep = "_"),.studycelltype)
}

CCl4_Krenkel_integrated <- DimPlot(mes_atlas.sct, cells.highlight= list(`Krenkel_Endothelial cells`,
                                                                        Krenkel_HSC,
                                                                        Krenkel_VSMC,
                                                                        `Krenkel_Portal fibroblasts`,
                                                                        `Krenkel_Mesothelial cells`,
                                                                        `Krenkel_Proliferating cells`),
                                   cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
                                   sizes.highlight = 0.1, raster = F)+ NoLegend() +NoAxes() + ggtitle("")

### Rosenthal ---------
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("WD"))


for (i in unique(mes_atlas.sct$celltype)) {
  Idents(mes_atlas.sct) <- "treated_celltype"
  .studycelltype <- WhichCells(mes_atlas.sct,cells = instudy,  idents = paste("Treatment",i, sep = "_"))
  assign(paste("Rosenthal",i, sep = "_"),.studycelltype)
}


WD_integrated <- DimPlot(mes_atlas.sct, cells.highlight= list(`Rosenthal_Endothelial cells`,
                                                              Rosenthal_HSC,
                                                              Rosenthal_VSMC,
                                                              `Rosenthal_Portal fibroblasts`,
                                                              `Rosenthal_Mesothelial cells`,
                                                              `Rosenthal_Proliferating cells`),
                         cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
                         sizes.highlight = 0.1, raster = F)+ NoLegend() +NoAxes() + ggtitle("") 

### Terkelsen ---------
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("CCl4_Terkelsen"))

for (i in unique(mes_atlas.sct$celltype)) {
  Idents(mes_atlas.sct) <- "treated_celltype"
  .studycelltype <- WhichCells(mes_atlas.sct,cells = instudy,  idents = paste("Treatment",i, sep = "_"))
  assign(paste("Terkelsen",i, sep = "_"),.studycelltype)
}


CCl4_Terkelsen_integrated <- DimPlot(mes_atlas.sct, cells.highlight= list(`Terkelsen_Endothelial cells`,
                                                                          Terkelsen_HSC,
                                                                          Terkelsen_VSMC,
                                                                          `Terkelsen_Portal fibroblasts`,
                                                                          `Terkelsen_Mesothelial cells`,
                                                                          `Terkelsen_Proliferating cells`),
                                     cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
                                     sizes.highlight = 0.1, raster = F)+ NoLegend()  +NoAxes() + ggtitle("") 


### Yangccl4  ---------
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("CCl4_Yang"))

for (i in unique(mes_atlas.sct$celltype)) {
  Idents(mes_atlas.sct) <- "treated_celltype"
  .studycelltype <- WhichCells(mes_atlas.sct,cells = instudy,  idents = paste("Treatment",i, sep = "_"))
  assign(paste("Yang",i, sep = "_"),.studycelltype)
}


CCl4_Yang_integrated <- DimPlot(mes_atlas.sct, cells.highlight= list(`Yang_Endothelial cells`,
                                                                     Yang_HSC,
                                                                     Yang_VSMC,
                                                                     `Yang_Portal fibroblasts`,
                                                                     `Yang_Mesothelial cells`,
                                                                     `Yang_Proliferating cells`),
                                cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
                                sizes.highlight = 0.1, raster = F)+ NoLegend() +NoAxes() + ggtitle("") 

### YangBDL  ---------
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("BDL"))


for (i in unique(mes_atlas.sct$celltype)) {
  Idents(mes_atlas.sct) <- "treated_celltype"
  .studycelltype <- WhichCells(mes_atlas.sct,cells = instudy,  idents = paste("Treatment",i, sep = "_"))
  assign(paste("BDL",i, sep = "_"),.studycelltype)
}

BDL_integrated <- DimPlot(mes_atlas.sct, cells.highlight= list(`BDL_Endothelial cells`,
                                                               BDL_HSC,
                                                               BDL_VSMC,
                                                               `BDL_Portal fibroblasts`,
                                                               `BDL_Mesothelial cells`,
                                                               `BDL_Proliferating cells`),
                          cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
                          sizes.highlight = 0.1, raster = F)+ NoLegend()  +NoAxes() + ggtitle("") 

### Dobieccl4  ---------
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("CCl4_Dobie"))



for (i in unique(mes_atlas.sct$celltype)[!unique(mes_atlas.sct$celltype) %in% c("Endothelial cells")]) {
  Idents(mes_atlas.sct) <- "treated_celltype"
  .studycelltype <- WhichCells(mes_atlas.sct,cells = instudy,  idents = paste("Treatment",i, sep = "_"))
  assign(paste("Dob",i, sep = "_"),.studycelltype)
}

CCl4_Dobie_integrated <- DimPlot(mes_atlas.sct, cells.highlight= list(Dob_HSC,
                                                                      Dob_VSMC,
                                                                      `Dob_Portal fibroblasts`,
                                                                      `Dob_Mesothelial cells`,
                                                                      `Dob_Proliferating cells`),
                                 cols.highlight = c("#00bfc4", "#f564e3","#b79f00","#00ba38","#f8766d"), cols= "lightgrey", 
                                 sizes.highlight = 0.1, raster = F)+ NoLegend() +NoAxes() + ggtitle("") 

### Dobieacute  ---------
Idents(mes_atlas.sct) <- "study"
instudy <- WhichCells(mes_atlas.sct, idents = c("CCl4_acute"))




for (i in unique(mes_atlas.sct$celltype)[!unique(mes_atlas.sct$celltype) %in% c("Mesothelial cells")]) {
  Idents(mes_atlas.sct) <- "treated_celltype"
  .studycelltype <- WhichCells(mes_atlas.sct,cells = instudy,  idents = paste("Treatment",i, sep = "_"))
  assign(paste("Acu",i, sep = "_"),.studycelltype)
}


CCl4_acute_integrated <- DimPlot(mes_atlas.sct, cells.highlight= list(`Acu_Endothelial cells`,
                                                                      Acu_HSC,
                                                                      Acu_VSMC,
                                                                      `Acu_Portal fibroblasts`,
                                                                      `Acu_Proliferating cells`),
                                 cols.highlight = c("#00bfc4","#b79f00","#00ba38","#f8766d","#619cff"), cols= "lightgrey", 
                                 sizes.highlight = 0.1, raster = F)+ NoLegend() +NoAxes() + ggtitle("") 


library(gridExtra)

grid.arrange(
  CCl4_acute_integrated,
  APAP_integrated,
  TAA_integrated,
  BDL_integrated,
  WD_integrated,
  CCl4_Krenkel_integrated,
  CCl4_Yang_integrated,
  CCl4_Dobie_integrated,
  CCl4_Terkelsen_integrated,
  Cancer_integrated,
  nrow = 1
  )

library(ggrastr)
for (i in levels(mes_atlas.sct$study)) {
  .a <- rasterize(get(paste0(i,"_integrated")), layers='Point', dpi=300)
assign(paste0(i,"_integrated2"), .a)
}



grid.arrange(
  CCl4_acute_integrated2,
  APAP_integrated2,
  TAA_integrated2,
  BDL_integrated2,
  WD_integrated2,
  CCl4_Krenkel_integrated2,
  CCl4_Yang_integrated2,
  CCl4_Dobie_integrated2,
  CCl4_Terkelsen_integrated2,
  Cancer_integrated2,
  nrow = 1
)


pdf(file = "~/HSC activation/Paper/figures/Treated_umap_split.pdf", width = 14, height = 1.5)
print(grid.arrange(
  CCl4_acute_integrated2,
  APAP_integrated2,
  TAA_integrated2,
  BDL_integrated2,
  WD_integrated2,
  CCl4_Krenkel_integrated2,
  CCl4_Yang_integrated2,
  CCl4_Dobie_integrated2,
  CCl4_Terkelsen_integrated2,
  Cancer_integrated2,
  nrow = 1)
)
dev.off()

