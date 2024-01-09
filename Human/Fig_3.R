library(plotROC)
library(ggplot2)
library(readxl)
data <- read_excel("~/HSC activation/COLEC10/For_r.xlsx")

str(data)
unique(colnames(data))

longtest <- melt_roc(data, "F-stage_1vshigher", c("COLEC10(ng/ml)","APRI","FIB-4 score"))
head(longtest)
longtest$D <- longtest$D.F.stage_1vshigher
a <- ggplot(longtest, aes(d = D, m = M, color = name)) + style_roc() + 
  geom_roc(n.cuts = 0) + theme_classic() + theme(legend.position = "none")+ geom_abline(intercept = 0, slope = 1, colour="#f1f1f1")

a

calc_auc(a)



longtest <- melt_roc(data, "Fstage2vshigher", c("COLEC10(ng/ml)","APRI","FIB-4 score"))
head(longtest)
longtest$D <- longtest$D.Fstage2vshigher
b <- ggplot(longtest, aes(d = D, m = M, color = name)) + geom_roc(n.cuts = 0) +
  style_roc() + theme_classic() + theme(legend.position = "none")+ geom_abline(intercept = 0, slope = 1, colour="#f1f1f1")

b
calc_auc(b)


longtest <- melt_roc(data, "Fstage3vshigher", c( "COLEC10(ng/ml)","APRI","FIB-4 score"))
head(longtest)
longtest$D <- longtest$D.Fstage3vshigher
c <- ggplot(longtest, aes(d = D, m = M, color = name)) + geom_roc(n.cuts = 0) +
  style_roc() + theme_classic()  + theme(legend.position = "none") + geom_abline(intercept = 0, slope = 1, colour="#f1f1f1")

c
calc_auc(c)

library(cowplot)
plot_grid(a,b,c, ncol = 3, rel_widths	= c(1,1,1.38))
e <- ggpubr::get_legend(c)
plot(e)
pdf(file = "~/HSC activation/Paper/figures/Fig_colec/AUROC.pdf", width = 10, height = 3.25)
plot_grid(a,b,c, ncol = 3)
dev.off()
