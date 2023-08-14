library(ggplot2)
library(tidyverse)
library(officer)
library(rvg)


##Graph theme 
my_theme <- theme(
   panel.grid = element_line(colour = "white"),
   panel.background = element_rect(fill = "white", color = "black"),
   panel.border = element_rect(fill = NA),
)



aa<-read.table("surfaceome_SP_EV.txt",sep="\t",header=T)


# add a column of NAs
aa$grade.of.EXPERIMENT <- "NO"
aa$grade.of.EXPERIMENT[aa$number.of.EXPERIMENT<7] <- "-10"
aa$grade.of.EXPERIMENT[aa$number.of.EXPERIMENT>=10&aa$number.of.EXPERIMENT<20] <- "10-20"
aa$grade.of.EXPERIMENT[aa$number.of.EXPERIMENT>=20&aa$number.of.EXPERIMENT<30] <- "20-30"
aa$grade.of.EXPERIMENT[aa$number.of.EXPERIMENT>=30&aa$number.of.EXPERIMENT<40] <- "30-40"
aa$grade.of.EXPERIMENT[aa$number.of.EXPERIMENT>=40] <- "40-"

# add a column of NAs
aa$grade.of.Median <- "0"
aa$grade.of.Median[log10(aa$Median)<7] <- "0-7"
aa$grade.of.Median[log10(aa$Median)>=7&log10(aa$Median)<8] <- "7-8"
aa$grade.of.Median[log10(aa$Median)>=8] <- "8-"

p<-ggplot(data=aa, aes(x=number.of.EXPERIMENT, y=Final.GESP.Score, col=grade.of.Median,size = grade.of.Median)) +geom_point( ) + my_theme + scale_color_manual(values = c("#F4EDCA","#00AFBB", "#E7B800", "#FC4E07"))
plot1<-ggplot(data=aa, aes(x=log10(Median), y=Final.GESP.Score)) +geom_point(size=1 ) + my_theme +geom_smooth(method=lm) ## + scale_x_reverse()
plot2<-ggplot(data=aa, aes(x=number.of.EXPERIMENT, y=Final.GESP.Score)) +geom_point(size=1) + my_theme +geom_smooth(method=lm)

pp<-plot1+plot2

doc <- read_pptx()
doc <- add_slide(doc)
# Add the plot
doc <- ph_with(doc,dml(ggobj = pp), location = ph_location_type(type = "body"))  
# Write the document to a file
print(doc, target = 'SP_EVplot_reverse.pptx')



doc <- read_pptx()
doc <- add_slide(doc)
# Add the plot
doc <- ph_with(doc,dml(ggobj = pp), location = ph_location_type(type = "body"))  
# Write the document to a file
print(doc, target = 'SP_EVplot2.pptx')


##############


bb<-read.table("humanplasma_HPA_EV.csv",sep=",",header=T)


# add a column of NAs
bb$grade.of.EXPERIMENT <- "NO"
bb$grade.of.EXPERIMENT[aa$number.of.EXPERIMENT<10] <- "-10"
bb$grade.of.EXPERIMENT[aa$number.of.EXPERIMENT>=10&aa$number.of.EXPERIMENT<20] <- "10-20"
bb$grade.of.EXPERIMENT[aa$number.of.EXPERIMENT>=20&aa$number.of.EXPERIMENT<30] <- "20-30"
bb$grade.of.EXPERIMENT[aa$number.of.EXPERIMENT>=30&aa$number.of.EXPERIMENT<40] <- "30-40"
bb$grade.of.EXPERIMENT[aa$number.of.EXPERIMENT>=40] <- "40-"



# add a column of NAs
bb$grade.of.Median <- "0"
bb$grade.of.Median[log10(bb$Median)<7] <- "0-7"
bb$grade.of.Median[log10(bb$Median)>=7&log10(aa$Median)<8] <- "7-8"
bb$grade.of.Median[log10(bb$Median)>=8] <- "8-"

ggplot(data=bb, aes(x=number.of.EXPERIMENT, y=log10(Concentration.ng.L.), col=grade.of.Median,size = grade.of.Median)) +geom_point() + my_theme + scale_color_manual(values = c("#F4EDCA","#00AFBB", "#E7B800", "#FC4E07"))

plot1<-ggplot(data=bb, aes(x=log10(Median), y=log10(Concentration.ng.L.))) +geom_point(size=1 ) + my_theme +geom_smooth(method=lm)
plot2<-ggplot(data=bb, aes(x=number.of.EXPERIMENT, y=log10(Concentration.ng.L.))) +geom_point(size=1 ) + my_theme +geom_smooth(method=lm)

pp<-plot1+plot2


doc <- read_pptx()
doc <- add_slide(doc)
# Add the plot
doc <- ph_with(doc,dml(ggobj = pp), location = ph_location_type(type = "body"))  
# Write the document to a file
print(doc, target = 'HP_EVplot2_reverse.pptx')

