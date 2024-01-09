library(readr)
library(tidyr)
library(ggplot2)

Tissue_rna <- read.table("~/Nanobodies/rna_tissue_consensus.tsv",header = T,sep="\t")
Tissue_rna <- Tissue_rna[,-1]
Tissue_rna$duplicate <- paste0(Tissue_rna$Gene.name, Tissue_rna$Tissue)

length(Tissue_rna$duplicate)
length(unique(Tissue_rna$duplicate))
Tissue_rna <- Tissue_rna[!duplicated(Tissue_rna$duplicate), ]

Tissue_rna <- Tissue_rna[,-4]

Wide_rna <- reshape(data = Tissue_rna,                                 
                    idvar = "Gene.name",
                    timevar = "Tissue",
                    direction = "wide")


rownames(Wide_rna) <- Wide_rna$Gene.name
Wide_rna <- Wide_rna[,-1]
colnames(Wide_rna) <- substr(colnames(Wide_rna),6, nchar(colnames(Wide_rna)))

Wide_rna$sum <- rowSums(Wide_rna)
Wide_rna$Liver_min_rest <- (Wide_rna$liver)*2 - Wide_rna$sum 

gene <- "COLEC10"
a <- data.frame("Tissue" = colnames(Wide_rna[gene,1:54]), "nTPM" = t(Wide_rna[gene,1:54]))
a <- a[]
a <- a[order(a$COLEC10, decreasing = T),]
x = 
  p<-ggplot(data=a, aes(x=reorder(Tissue, +get(gene)), y=get(gene))) +
  geom_bar(stat="identity")
p + theme(axis.text.x=element_text(angle = 0, hjust = 0)) + coord_flip() + theme_classic()
p

pdf(file = "~/HSC activation/Paper/figures/Fig_colec/HPA.pdf", width = 4.5, height = 8)
p + theme(axis.text.x=element_text(angle = 0, hjust = 0)) + coord_flip() + theme_classic()

dev.off()
