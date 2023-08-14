####optimal cut off using cutoffr#######
library(survival)
library(Publish)
library(survminer)
library(reshape2)
library(reshape)
library(cutpointr)
library(prodlim)



###Example_data
#RCC<-read.table("Kidney_Cancer_summary_PSMs_grouping_upper70_tumor.txt",sep="\t",header=T)
#RCC_ID<-read.table("Kidney_sample_info.txt",sep="\t",header=T)


#####survival_function################

Input_data  <- function(a,b,pt) {        ##a=BRCA,b=BRCA_ID   # R function with return
a_new<- a[a$UNIPROT%in%pt,]
data<-a_new[,c(4:ncol(a_new))] 
rownames(data)<-a_new[,1]
targets<-b[,c(1,4,5)]

Input<-list()
Input$data<-data
Input$targets<-targets
return(Input)
}



HR_stat <- function(a,b) {  #a=data
HazardRatio <- sapply(1:nrow(a),function(i){

new<-merge(b,t(a[i,]),by.x = "Sample" , by.y =0,all.x=TRUE)

protein <- colnames(new)[4]
cp <-eval(parse(text=paste("cutpointr(new,", protein,",OS_Death,method = maximize_spline_metric, metric = sum_sens_spec, na.rm = TRUE)")))
cutopt_point <-cp$optimal_cutpoint
cutopt_point


new$cutopt <- factor(as.integer(new[,4]>= as.integer(cutopt_point)))

model<-coxph(Surv(OS_MONTHS,OS_Death)~ cutopt,data=new,na.action = na.omit)
x<-summary(model)
p.value<-signif(x$logtest["pvalue"], digits=2)
likelihood.test<-signif(x$logtest["test"], digits=2)
beta<-signif(x$coef[1], digits=2);#coeficient beta
HR <-signif(x$coef[2], digits=2);#exp(beta)
HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
HR_CI<-paste0(" (",HR.confint.lower, "-", HR.confint.upper, ")")
res<-c(beta, HR, HR_CI,likelihood.test, p.value)
names(res)<-c("beta", "HR","HR_CI(95% CI for HR)", "Likelihood.test","p.value")
return(res)
   })

HR_result<- t(as.data.frame(HazardRatio, check.names = FALSE))
rownames(HR_result)<-row.names(a)
HR_result <- cbind(rownames(HR_result), data.frame(HR_result, row.names=NULL))
return(HR_result)
}





##
create_save_survival_graphs <- function(a,b,output_dir) {
#Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
 
# Iterate 
create_survival_graphs  <- sapply(1:nrow(a),function(i){
gene_name <-rownames(a[i,])
new<-merge(b,t(a[i,]),by.x = "Sample" , by.y =0,all.x=TRUE)### uniprotiD 로 바꾸자


protein <- colnames(new)[4]
cp <-eval(parse(text=paste("cutpointr(new,", protein,",OS_Death,method = maximize_spline_metric, metric = sum_sens_spec, na.rm = TRUE)")))
cutopt_point <-cp$optimal_cutpoint
cutopt_point
new$cutopt <- factor(as.integer(new[,4]>= as.integer(cutopt_point)))
fit <- survfit(Surv(OS_MONTHS,OS_Death)~cutopt , data =new)

grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

plot_path <- file.path(output_dir, paste(gene_name, ".pdf", sep = ""))
pdf(plot_path)
p<-ggsurvplot(fit,data =new,palette =c("blue","red"),conf.int = TRUE,pval = TRUE, pval.method = TRUE,risk.table = TRUE)
print(p)
dev.off()

})

}


survival_anal<- function(a,b,pt,output_dir) {    

Input<-Input_data(a,b,pt)###cohort name만 바꾸면 됨. 
head(Input$data,5)
result<-HR_stat(Input$data,Input$targets)
create_save_survival_graphs(Input$data,Input$targets,output_dir)
return(result)

}




#####################################
################## Data input #######
RCC<-read.table("Kidney_Cancer_summary_PSMs_grouping_upper70_tumor.txt",sep="\t",header=T)
HCC<-read.table("Liver_Cancer_summary_PSMs_grouping_upper70_tumor.txt",sep="\t",header=T)
GBM<-read.table("GBM_Cancer_summary_PSMs_grouping_upper70_tumor.txt",sep="\t",header=T)
LUAD<-read.table("LUAD_Cancer_summary_PSMs_grouping_upper70_tumor.txt",sep="\t",header=T)
JHU<-read.table("PDCJHU_Cancer_summary_PSMs_grouping_upper70_tumor.txt",sep="\t",header=T)
LUSCC<-read.table("LUSCC_Cancer_summary_PSMs_grouping_upper70_tumor.txt",sep="\t",header=T)

BRCA<-read.table("Breast_Cancer_summary_PSMs_grouping_upper70_tumor.txt",sep="\t",header=T)
GA<-read.table("Gastric_Cancer_summary_PSMs_grouping_upper70_tumor.txt",sep="\t",header=T)
KU<-read.table("PDCKU_Cancer_summary_PSMs_grouping_upper70_tumor.txt",sep="\t",header=T)

OV<-read.table("PDCKU_Cancer_summary_PSMs_grouping_upper70_tumor.txt",sep="\t",header=T)


#######################################
###########sample info input /UCEC,HNSCC,Colon X
RCC_ID<-read.table("Kidney_sample_info.txt",sep="\t",header=T)
HCC_ID<-read.table("HCC_sample_info.txt",sep="\t",header=T)
GBM_ID<-read.table("GBM_sample_info.txt",sep="\t",header=T)
LUAD_ID<-read.table("LUAD_sample_info.txt",sep="\t",header=T)
JHU_ID<-read.table("JHU_sample_info.txt",sep="\t",header=T)
LUSCC_ID<-read.table("LUSCC_sample_info.txt",sep="\t",header=T)

BRCA_ID<-read.table("Breast_sample_info.txt",sep="\t",header=T)
GA_ID<-read.table("GA_sample_info.txt",sep="\t",header=T)
KU_ID<-read.table("KU_sample_info.txt",sep="\t",header=T)

OV_ID<-read.table("OV_sample_info.txt",sep="\t",header=T)



#####protein list to check ##############
pt <-c('P17948','P29323','P54753','Q04912','P21860','P21709','P54764','P35916','Q15303') ######cf. RCCdata; 'P17948','P21860','P35916' exist,'P21860'==not significant
pt <-c('P08581','P29317','P04626','Q13308') ###FDA
pt <-c('P08581','P29317','P04626','Q13308','P54760','Q08345','P17948','P00533','P29323','P54753','P36897','Q04912','Q04771','P21860','P21709','P06213','Q01974','Q12866','Q16832','P54764','P08069','P09619','P21802','P30530','P16234','P35916','Q15303','Q15375','P37173','P36896')
pt <-c('P50281')
#####Let's do######

RCC_result<-survival_anal(RCC,RCC_ID,pt,"RCC_graph")
write.table(RCC_result,"stat_RCC_HR_RKcandidate.txt",sep="\t",row.names=FALSE)

HCC_result<-survival_anal(HCC,HCC_ID,pt,"HCC_graph")
write.table(HCC_result,"stat_HCC_HR_RKcandidate.txt",sep="\t",row.names=FALSE)

GBM_result<-survival_anal(GBM,GBM_ID,pt,"GBM_graph")
write.table(GBM_result,"stat_GBM_HR_RKcandidate.txt",sep="\t",row.names=FALSE)

LUAD_result<-survival_anal(LUAD,LUAD_ID,pt,"LUAD_graph")
write.table(LUAD_result,"stat_LUAD_HR_RKcandidate.txt",sep="\t",row.names=FALSE)

JHU_result<-survival_anal(JHU,JHU_ID,pt,"JHU_graph")
write.table(JHU_result,"stat_JHU_HR_RKcandidate.txt",sep="\t",row.names=FALSE)

LUSCC_result<-survival_anal(LUSCC,LUSCC_ID,pt,"LUSCC_graph")
write.table(LUSCC_result,"stat_LUSCC_HR_RKcandidate.txt",sep="\t",row.names=FALSE)

BRCA_result<-survival_anal(BRCA,BRCA_ID,pt,"BRCA_graph")
write.table(BRCA_result,"stat_BRCA_HR_RKcandidate.txt",sep="\t",row.names=FALSE)

GA_result<-survival_anal(GA,GA_ID,pt,"GA_graph")
write.table(GA_result,"stat_GA_HR_RKcandidate.txt",sep="\t",row.names=FALSE)

KU_result<-survival_anal(KU,KU_ID,pt,"KU_graph")
write.table(KU_result,"stat_KU_HR_RKcandidate.txt",sep="\t",row.names=FALSE)


