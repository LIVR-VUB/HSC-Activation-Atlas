library(Seurat)
library(ggrastr)

unique(Idents(rs039_subset))

commongenes <- intersect(rownames(rs039_subset@assays$RNA), rownames(hscs.filtered@assays$RNA))
rs039_subset <- subset(rs039_subset, features = commongenes)
DefaultAssay(rs039_subset) <- "RNA"
rs039_subset <- SCTransform(rs039_subset)
DefaultAssay(hscs.filtered) <- "integrated"
DefaultAssay(rs039_subset) <- "SCT"


DefaultAssay(hscs.filtered) <- "integrated"
DefaultAssay(rs039_subset) <- "SCT"

Idents(hscs.filtered)


anchors <- FindTransferAnchors(reference = hscs.filtered, query = rs039_subset,
                               dims = 1:20, reference.reduction = "pca", 
                               normalization.method = "SCT")
predictions <- TransferData(anchorset = anchors, refdata = hscs.filtered$subtype,
                            dims = 1:20)
rs039_subset <- AddMetaData(rs039_subset, metadata = predictions)

pdf(file = "~/HSC activation/Paper/Suppfigures/S9/orig_dimplot_pred.pdf", width = 5, height = 4)
a <- DimPlot(rs039_subset, group.by = "predicted.id") + NoAxes() + ggtitle("")
rasterize(a, layers='Point', dpi=300)
dev.off()

DimPlot(rs039_subset, group.by = "predicted.id")
Idents(rs039_subset) <- "predicted.id"


hscs.filtered <- RunUMAP(hscs.filtered, dims = 1:20, reduction = "pca", return.model = TRUE)
rs039_subset <- MapQuery(anchorset = anchors, reference = hscs.filtered, query = rs039_subset,
                         refdata = list(celltype = "subtype"), reference.reduction = "pca",
                         reduction.model = "umap"
)


rs039_subset$predicted.id <- factor(rs039_subset$predicted.id ,
                                    levels = c("quiescent","initiatory","myofibroblasts"))
Idents(rs039_subset) <- "predicted.id"

p1 <- DimPlot(hscs.filtered, reduction = "umap", group.by = "subtype", label = F, label.size = 3,
              repel = TRUE, raster = F) + NoLegend() + ggtitle("") 
p2 <- DimPlot(rs039_subset, reduction = "ref.umap", group.by = "predicted.id", label = F,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("") 
a <- p1 + p2
pdf(file = "~/HSC activation/Paper/Suppfigures/S9/projection.pdf", width = 8, height = 4)
rasterize(a, layers='Point', dpi=300)
dev.off()

DefaultAssay(hscs.filtered)<- "RNA"
DefaultAssay(rs039_subset)<- "RNA"

p3 <- FeaturePlot(hscs.filtered, "Rgs5", raster = F) + scale_color_viridis() + NoLegend()
p4 <- FeaturePlot(rs039_subset, "Rgs5", raster = F, reduction = "ref.umap") + scale_color_viridis()

pdf(file = "~/HSC activation/Paper/Suppfigures/S9/rgs5.pdf", width = 8, height = 4)
a <- plot_grid(p3, p4, ncol=2, rel_widths = c(1,1.25))
rasterize(a, layers='Point', dpi=300)
dev.off()

p3 <- FeaturePlot(hscs.filtered, "Ccl2", raster = F) + scale_color_viridis() + NoLegend()
p4 <- FeaturePlot(rs039_subset, "Ccl2", raster = F, reduction = "ref.umap") + scale_color_viridis()
plot_grid(p3, p4, ncol=2, rel_widths = c(1,1.25))

pdf(file = "~/HSC activation/Paper/Suppfigures/S9/Ccl2.pdf", width = 8, height = 4)
a <- plot_grid(p3, p4, ncol=2, rel_widths = c(1,1.25))
rasterize(a, layers='Point', dpi=300)
dev.off()



p3 <- FeaturePlot(hscs.filtered, "Col1a1", raster = F) + scale_color_viridis() + NoLegend()
p4 <- FeaturePlot(rs039_subset, "Col1a1", raster = F, reduction = "ref.umap") + scale_color_viridis()
plot_grid(p3, p4, ncol=2, rel_widths = c(1,1.25))

pdf(file = "~/HSC activation/Paper/Suppfigures/S9/Col1a1.pdf", width = 8, height = 4)
a <- plot_grid(p3, p4, ncol=2, rel_widths = c(1,1.25))
rasterize(a, layers='Point', dpi=300)
dev.off()


DefaultAssay(rs039_subset) <- "RNA"

pdf(file = "~/HSC activation/Paper/Suppfigures/S9/qgenes.pdf", width = 6, height = 4)
VlnPlot(rs039_subset, c("Lrat","Reln","Rgs5"), pt.size = 0)&
  xlab("") 
dev.off()
pdf(file = "~/HSC activation/Paper/Suppfigures/S9/igenes.pdf", width = 6, height = 4)
VlnPlot(rs039_subset, c("Ccl2","Ccl7","Mt2"), pt.size = 0)&
  xlab("") 
dev.off()
pdf(file = "~/HSC activation/Paper/Suppfigures/S9/mgenes.pdf", width = 6, height = 4)
VlnPlot(rs039_subset, c("Col1a1","Loxl1","Mmp2"), pt.size = 0)&
  xlab("") 
dev.off()
