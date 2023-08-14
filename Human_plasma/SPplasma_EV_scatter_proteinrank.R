library(ggplot2)
library(tidyverse)
library(officer)
library(rvg)
library(patchwork)
library(purrr)


##Graph theme 
my_theme1 <- theme(
   panel.grid = element_line(colour = "white"),
   panel.background = element_rect(fill = "white", color = "black"),
   panel.border = element_rect(fill = NA),legend.position = "none"
)

my_theme2 <- theme(
   panel.grid = element_line(colour = "white"),
   panel.background = element_rect(fill = "white", color = "black"),
   panel.border = element_rect(fill = NA)
)


####################################
##############HP####################


bb<-read.table("humanplasma_HPA_EV_3250.csv",sep=",",header=T)


bb_sub <- bb[ which(bb$X3250_overlap !="NA"), ]
mid <- mean(-(bb_sub$protein_rank_3250))
plot1<-ggplot(data=bb_sub, aes(x=log2(Median), y=log2(Concentration.ng.L.),col=-(protein_rank_3250))) +geom_point(size=2 ) + my_theme1 +geom_smooth(method=lm)+ scale_color_gradient2(midpoint = mid, low = "blue", mid = "white", high = "red", space = "Lab" )

plot2<-ggplot(data=bb_sub, aes(x=number.of.EXPERIMENT, y=log2(Concentration.ng.L.),col=-(protein_rank_3250))) +geom_point(size=2 ) + my_theme1 +geom_smooth(method=lm)+ scale_color_gradient2(midpoint = mid, low = "blue", mid = "white", high = "red", space = "Lab" )
pp<-plot1+plot2



#############################
#############SP##############

aa<-read.table("C:\\Users\\kimj49\\Box\\Multicancer Biomarker\\Data_library\\Surfaceome\\surfaceome_SP_EV_3250.csv",sep=",",header=T)

aa_sub <- aa[ which(aa$X3250_overlap !="NA"), ]
mid <- mean(-(aa_sub$protein_rank_3250))
plot3<-ggplot(data=aa_sub, aes(x=log2(Median), y=Final.GESP.Score, col=-(protein_rank_3250))) +geom_point(size=2 ) + my_theme1 +geom_smooth(method=lm) + scale_color_gradient2(midpoint = mid, low = "blue", mid = "white", high = "red", space = "Lab" )
plot4<-ggplot(data=aa_sub, aes(x=number.of.EXPERIMENT, y=Final.GESP.Score,col=-(protein_rank_3250))) +geom_point(size=2) + my_theme1 +geom_smooth(method=lm)+ scale_color_gradient2(midpoint = mid, low = "blue", mid = "white", high = "red", space = "Lab" )


plot5<-ggplot(data=aa_sub, aes(x=number.of.EXPERIMENT, y=Final.GESP.Score,col=-(protein_rank_3250))) +geom_point(size=2) + my_theme2 +geom_smooth(method=lm)+ scale_color_gradient2(midpoint = mid, low = "blue", mid = "white", high = "red", space = "Lab" )






##############


doc <- read_pptx("HPSP_EVplot_commonprotein.pptx")
doc <- add_slide(doc, layout = "Four contents") |> 
  ph_with(dml(ggobj = plot1), location = ph_location_label(ph_label = "Content Placeholder 1")) |> 
  ph_with(dml(ggobj = plot2), location = ph_location_label(ph_label = "Content Placeholder 2")) |> 
  ph_with(dml(ggobj = plot3), location = ph_location_label(ph_label = "Content Placeholder 3")) |> 
  ph_with(dml(ggobj = plot4), location = ph_location_label(ph_label = "Content Placeholder 4"))


print(doc, target = 'HPSP_EVplot_commonprotein.pptx')


doc <- read_pptx()
doc <- add_slide(doc)
# Add the plot
doc <- ph_with(doc,dml(ggobj = plot5), location = ph_location_type(type = "body"))  
# Write the document to a file
print(doc, target = 'test_legend.pptx')

