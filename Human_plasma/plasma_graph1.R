aa<-read.table("plasma_graph1.txt",sep="\t",header=T)

library(ggplot2)
library(ggExtra)
library(RColorBrewer)
library(cowplot)
library(ggpubr)

##Graph theme 
my_theme <- theme(
   panel.grid = element_line(colour = "white"),
   panel.background = element_rect(fill = "white", color = "black"),
   panel.border = element_rect(fill = NA),
)



ggplot(aa, aes(x=log10(Ayuko_Ev_q), y=log10(Concentration))) + geom_point()+my_theme +
        geom_text_repel(size = 4.5) +
        scale_color_manual(values=c("blue","black", "red"))







aa<-read.table("Experiment_5top5bottom_genes_Ayuko_q_before.txt",sep="\t",header=T)
library(MASS) 
library(reshape2) 
library(reshape) 
library(RColorBrewer)
melt_aa <- melt(aa, id = c("cancer_type","Sample","type"))
ggplot(melt_aa, aes(x=variable, y=log10(value), fill=cancer_type)) + geom_boxplot()+my_theme +scale_fill_brewer(palette="Set3")





