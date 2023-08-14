library(ggrepel)
##Graph theme 
my_theme <- theme(
   panel.grid = element_line(colour = "white"),
   panel.background = element_rect(fill = "white", color = "black"),
   panel.border = element_rect(fill = NA),
)


de<-read.table("Ev_vesiclepedia_cancertype_protein_3250.txt",sep="\t",header=T)

# add a column of NAs
de$number <- "NO"
de$number[de$number.of.EXPERIMENT > 54] <- "UP"
de$number[de$GENE.SYMBOL=="EEF2"] <- "UP"
de$number[de$GENE.SYMBOL=="HSP90AB1"] <- "UP"

##########
de$number[de$GENE.SYMBOL=="UBAP1"] <- "down"
de$number[de$GENE.SYMBOL=="UHRF2"] <- "down"
de$number[de$GENE.SYMBOL=="VPS37A"] <- "down"
de$number[de$GENE.SYMBOL=="WDR3"] <- "down"
de$number[de$GENE.SYMBOL=="YTHDC1"] <- "down"
de$number[de$GENE.SYMBOL=="SERPINF2"] <- "down"
de$number[de$GENE.SYMBOL=="STAG1"] <- "down"
de$number[de$GENE.SYMBOL=="SUPT6H"] <- "down"
de$number[de$GENE.SYMBOL=="TP53BP1"] <- "down"
de$number[de$GENE.SYMBOL=="TXNDC9"] <- "down"



de$delabel <- NA
de$delabel[de$number  != "NO"] <- de$GENE.SYMBOL[de$number  != "NO"]

ggplot(de, aes(x=Seq, y=number.of.EXPERIMENT, col=number, label=delabel)) + geom_point()+my_theme +
        geom_text_repel(size=4, max.overlaps=Inf) +
        scale_color_manual(values=c("blue","black", "red"))






library(MASS) 
library(reshape2) 
library(reshape) 
library(RColorBrewer)
library(ggpubr)
library(patchwork)


aa<-read.table("Experiment_10top_genes_Ayuko_q.txt",sep="\t",header=T)
bb<-read.table("Experiment_10bottom_genes_Ayuko_q.txt",sep="\t",header=T)
melt_aa <- melt(aa, id = c("cancer_type","Sample","type"))
melt_bb <- melt(bb, id = c("cancer_type","Sample","type"))


p1<-ggplot(melt_aa, aes(x=variable, y=log2(value), fill=cancer_type)) + geom_boxplot(width=0.6,position=position_dodge(0.8),outlier.shape = NA)+ my_theme +scale_fill_brewer(palette="Set3")+ylim(20,40)+geom_dotplot(aes(fill = cancer_type),dotsize=0.6,binaxis='y', stackdir='center', position = position_dodge(0.8))
p2<-ggplot(melt_bb, aes(x=variable, y=log2(value), fill=cancer_type)) + geom_boxplot(width=0.6,position=position_dodge(0.8),outlier.shape = NA)+ my_theme +scale_fill_brewer(palette="Set3")+ylim(20,40)+geom_dotplot(aes(fill = cancer_type),dotsize=0.6,binaxis='y', stackdir='center', position = position_dodge(0.8))



p1+p2+plot_layout(ncol = 1)



#+my_theme 