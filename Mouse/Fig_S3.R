library(DESeq2)
library(Seurat)
library(ggplot2)
library(viridis)


colnames(RNA_Seq_HSC)
head(mcols(RNA_Seq_HSC))
head(colData(RNA_Seq_HSC))

dds <- DESeq(RNA_Seq_HSC, betaPrior = T)
res <- results(RNA_Seq_HSC, contrast=c("Condition","Healthy","Chronic Injury"), alpha = 0.05 )
res <- na.omit(res)

healthy_geneset <- rownames(res[res$padj<0.05 & res$log2FoldChange >1,])
healthy_geneset <- unique(Annotation[Annotation$ensembl_gene_id %in% healthy_geneset,"mgi_symbol"])
healthy_geneset <- na.omit(healthy_geneset)

RNA_Seq_HSC$condition <- "Healthy"
RNA_Seq_HSC$condition[RNA_Seq_HSC$Condition == "Acute Injury"] <- "Acute_injury"
RNA_Seq_HSC$condition[RNA_Seq_HSC$Condition == "Acute - 48h recovery"] <- "Acute_injury2d_recovery"
RNA_Seq_HSC$condition[RNA_Seq_HSC$Condition == "Acute - 1w recovery"] <- "Acute_injury1w_recovery"
RNA_Seq_HSC$condition[RNA_Seq_HSC$Condition == "Chronic Injury"] <- "Chronic_injury"
RNA_Seq_HSC$condition[RNA_Seq_HSC$Condition == "Chronic - 48h recovery"] <- "Chronic_injury2d_recovery"
RNA_Seq_HSC$condition[RNA_Seq_HSC$Condition == "Chronic - 2w recovery"] <- "Chronic_injury2w_recovery"
unique(RNA_Seq_HSC$condition)



for (i in unique(RNA_Seq_HSC$condition)) {
  colData(RNA_Seq_HSC)[ , ncol(colData(RNA_Seq_HSC)) + 1] <- "others"
  colnames(colData(RNA_Seq_HSC))[ncol(colData(RNA_Seq_HSC))] <- paste0( i,"_vsothers")
  colData(RNA_Seq_HSC)[RNA_Seq_HSC$condition == i, paste0( i,"_vsothers")] <- "group_of_interest"
}


countss <- counts(RNA_Seq_HSC)
metadata <- colData(RNA_Seq_HSC)
dds_Healthy <- DESeqDataSetFromMatrix(countData =  countss,colData =  metadata,
                                               design = ~Healthy_vsothers)
dds_Acute_injury <- DESeqDataSetFromMatrix(countData =  countss,colData =  metadata,
                                               design = ~Acute_injury_vsothers)

dds_Acute_injury2d_recovery <- DESeqDataSetFromMatrix(countData =  countss,colData =  metadata,
                                               design = ~Acute_injury2d_recovery_vsothers)

dds_Acute_injury1w_recovery <- DESeqDataSetFromMatrix(countData =  countss,colData =  metadata,
                                               design = ~Acute_injury1w_recovery_vsothers)

dds_Chronic_injury <- DESeqDataSetFromMatrix(countData =  countss,colData =  metadata,
                                               design = ~Chronic_injury_vsothers)

dds_Chronic_injury2d_recovery <- DESeqDataSetFromMatrix(countData =  countss,colData =  metadata,
                                               design = ~Chronic_injury2d_recovery_vsothers)

dds_Chronic_injury2w_recovery <- DESeqDataSetFromMatrix(countData =  countss,colData =  metadata,
                                               design = ~Chronic_injury2w_recovery_vsothers)


for (i in unique(RNA_Seq_HSC$condition)) {
  .dds <- get(paste0("dds_",i))
  .dds <- DESeq(.dds)
  res <- results(.dds,contrast= c(paste0( i,"_vsothers"),"group_of_interest", "others"),
                 listValues=c(1,-1/6), alpha = 0.05 )
  res <- na.omit(res)
  
  
  .temp_geneset <- rownames(res[res$padj<0.05 & res$log2FoldChange >1,])
  .temp_geneset <- unique(Annotation[Annotation$ensembl_gene_id %in% .temp_geneset,"mgi_symbol"])
  .temp_geneset <- na.omit(.temp_geneset)
  assign(paste0(i,"_geneset"), .temp_geneset)
}
colnames(colData(RNA_Seq_HSC))

for (i in unique(RNA_Seq_HSC$condition)) {
  conditions <- unique(RNA_Seq_HSC$Condition)
  conditions <- conditions[!conditions == i]  
  
  res <- results(RNA_Seq_HSC,contrast=list(c(as.character(i)), c(as.character(conditions))),listValues=c(1,-1/5), alpha = 0.05 )
  res <- na.omit(res)
  
  
  .temp_geneset <- rownames(res[res$padj<0.05 & res$log2FoldChange >1,])
  .temp_geneset <- unique(Annotation[Annotation$ensembl_gene_id %in% .temp_geneset,"mgi_symbol"])
  .temp_geneset <- na.omit(.temp_geneset)
  assign(paste0(i,"_geneset"), .temp_geneset)
}
res <- results(RNA_Seq_HSC, contrast=c("Condition","Healthy","Chronic Injury"), alpha = 0.05 )
res <- na.omit(res)

healthy_geneset <- rownames(res[res$padj<0.05 & res$log2FoldChange >1,])
healthy_geneset <- unique(Annotation[Annotation$ensembl_gene_id %in% healthy_geneset,"mgi_symbol"])
healthy_geneset <- na.omit(healthy_geneset)

genesets <- list("healthy" = Healthy_geneset, "acute_injury" = Acute_injury_geneset,
                 "2d_post_acute" = Acute_injury2d_recovery_geneset, "1w_post_acute" = Acute_injury1w_recovery_geneset,
                 "chronic_injury" = Chronic_injury_geneset, "2d_post_chronic" = Chronic_injury2d_recovery_geneset,
                 "2w_post_chronic" = Chronic_injury2w_recovery_geneset)


sink('bulk_genesets.csv')
print(genesets)
sink()

DefaultAssay(sct.hscs)
sct.hscs <- AddModuleScore(sct.hscs, features = genesets, name = c("Healthy", "Acute injury","2d_post_acute",
                                                                   "1w_post_acute","Chronic injury","2d Recovery","2 week recovery"))

FeaturePlot(sct.hscs, c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5","Cluster6", "Cluster7"),
            reduction = "pca", order = T)
FeaturePlot(sct.hscs, c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5","Cluster6", "Cluster7"),
            reduction = "umap")

sct.hscs$Healthy <- sct.hscs$Healthy1
sct.hscs$Acute_injury <- sct.hscs$Acute.injury2
sct.hscs$Chronic_injury <- sct.hscs$Chronic.injury5
sct.hscs$Recovery_2days <- sct.hscs$X2d.Recovery6
sct.hscs$Recovery_2weeks <- sct.hscs$X2.week.recovery7

FeaturePlot(sct.hscs, c("Healthy","Acute_injury","Chronic_injury","Recovery_2days","Recovery_2weeks"),
            reduction = "pca", order = T, ncol = 1) & scale_color_viridis() & NoAxes() & ggtitle("")
pdf(file = "~/HSC activation/Paper/Suppfigures/S2/PCA_enrichment1.pdf", width = 3, height = 2.5) 
dimplt_hsc <- FeaturePlot(sct.hscs, c("Healthy"),
            reduction = "pca", order = T, ncol = 1) & scale_color_viridis() & NoAxes() & ggtitle("")
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off()
pdf(file = "~/HSC activation/Paper/Suppfigures/S2/PCA_enrichment2.pdf", width = 3, height = 2.5) 
dimplt_hsc <- FeaturePlot(sct.hscs, c("Acute_injury"),
                          reduction = "pca", order = T, ncol = 1) & scale_color_viridis() & NoAxes() & ggtitle("")
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off()
pdf(file = "~/HSC activation/Paper/Suppfigures/S2/PCA_enrichment3.pdf", width = 3, height = 2.5) 
dimplt_hsc <- FeaturePlot(sct.hscs, c("Chronic_injury"),
                          reduction = "pca", order = T, ncol = 1) & scale_color_viridis() & NoAxes() & ggtitle("")
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off()
pdf(file = "~/HSC activation/Paper/Suppfigures/S2/PCA_enrichment4.pdf", width = 3, height = 2.5) 
dimplt_hsc <- FeaturePlot(sct.hscs, c("Recovery_2days"),
                          reduction = "pca", order = T, ncol = 1) & scale_color_viridis() & NoAxes() & ggtitle("")
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off()
pdf(file = "~/HSC activation/Paper/Suppfigures/S2/PCA_enrichment5.pdf", width = 3, height = 2.5) 
dimplt_hsc <- FeaturePlot(sct.hscs, c("Recovery_2weeks"),
                          reduction = "pca", order = T, ncol = 1) & scale_color_viridis() & NoAxes() & ggtitle("")
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off()



DefaultAssay(sct.hscs) <- "RNA"
Idents(sct.hscs) <- "subtype"
VlnPlot(sct.hscs, c("Rgs5","Ccl2","Col1a1"), pt.size = 0)
FeaturePlot(sct.hscs, c("Rgs5","Ccl2","Col1a1","Ctss"), reduction = "pca", order = T)

 
pdf(file = "~/HSC activation/Paper/Suppfigures/S2/Rgs5.pdf", width = 3, height = 3) 
dimplt_hsc <- FeaturePlot(sct.hscs, c("Rgs5"), reduction = "pca", order = F)  & NoAxes() & ggtitle("")
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off() 
pdf(file = "~/HSC activation/Paper/Suppfigures/S2/Ccl2.pdf", width = 3, height = 3) 
dimplt_hsc <- FeaturePlot(sct.hscs, c("Ccl2"), reduction = "pca", order = F)  & NoAxes() & ggtitle("")
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off() 
pdf(file = "~/HSC activation/Paper/Suppfigures/S2/Col1a1.pdf", width = 3, height = 3) 
dimplt_hsc <- FeaturePlot(sct.hscs, c("Col1a1"), reduction = "pca", order = F)  & NoAxes() & ggtitle("")
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off() 
pdf(file = "~/HSC activation/Paper/Suppfigures/S2/Ctss.pdf", width = 3, height = 3) 
dimplt_hsc <- FeaturePlot(sct.hscs, c("Ctss"), reduction = "pca", order = T) & NoAxes() & ggtitle("")
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off() 


pdf(file = "~/HSC activation/Paper/Suppfigures/S2/dimplot.pdf", width = 3, height = 3) 
dimplt_hsc <- DimPlot(sct.hscs, group.by = "subtype",reduction = "pca", cols = c("#f8786f","#c77cff","#00ba38","#619cff")) &NoAxes() & ggtitle("") &NoLegend()
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off() 

Idents(sct.hscs) <- "condition"
inactivating_cells <- WhichCells(sct.hscs,idents = "8w Regression" )

pdf(file = "~/HSC activation/Paper/Suppfigures/S2/8wregression.pdf", width = 6, height = 6) 
dimplt_hsc <- DimPlot(sct.hscs, group.by = "condition",reduction = "pca", cells.highlight = inactivating_cells, sizes.highlight = 0.1) &NoAxes() & ggtitle("") &NoLegend() 
dimplt_hsc <- rasterize(dimplt_hsc, layers='Point', dpi=300)
plot(dimplt_hsc)
dev.off() 

Idents(sct.hscs) <- "condition"
wd <- subset(sct.hscs, idents = "8w Regression")
meta.data <- wd[[]]

counts <- group_by(meta.data,  subtype) %>% summarise(count = n())
totals <- sum(counts$count)

counts$total <- rep(totals$Frequency,times = b$Freq)
counts$percent <- (counts$count/totals) * 100
pdf(file = "~/HSC activation/Paper/Suppfigures/S2/barplot.pdf", width = 3, height = 4) 
print(ggplot(counts, aes( x = "everything",y = counts$percent, fill = subtype)) + 
                      scale_fill_manual("legend", values = c("quiescent" = "#f8786f",
                                                             "intermediate" = "#c77cff",
                                                             "initiatory" = "#00ba38",
                                                             "myofibroblasts" = "#619cff")))+
  geom_bar(stat = 'identity') + theme_classic() + theme(axis.text.x=element_text(angle=0,hjust=0.7,vjust=0.85)) +
  theme(axis.text = element_text(size = 8))
dev.off() 


