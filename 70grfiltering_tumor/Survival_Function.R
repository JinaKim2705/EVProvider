####optimal cut off using cutoffr#######
library(survival)
library(Publish)
library(survminer)
library(reshape2)
library(reshape)

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
new$median <- factor(as.integer(new[,4]>= median(new[,4],na.rm=TRUE)))
model<-coxph(Surv(OS_MONTHS,OS_Death)~ median ,data=new,na.action = na.omit)
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
new$median <- factor(as.integer(new[,4]>= median(new[,4],na.rm=TRUE)))

fit <- survfit(Surv(OS_MONTHS,OS_Death)~median, data =new)

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
OV<-read.table("OV_Cancer_summary_PSMs_grouping_upper70_tumor.txt",sep="\t",header=T)
CRC<-read.table("Colon_Cancer_summary_PSMs_grouping_upper70_tumor.txt",sep="\t",header=T)


#######################################
###########sample info input /UCEC,HNSCC X 
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
CRC_ID<-read.table("CRC_sample_info.txt",sep="\t",header=T)



#####protein list to check ##############
pt <-c('P17948','P29323','P54753','Q04912','P21860','P21709','P54764','P35916','Q15303') ######cf. RCCdata; 'P17948','P21860','P35916' exist,'P21860'==not significant

pt <-c('P08581','P29317','P04626','Q13308') ###FDA


####RK70###
pt <-c('Q9UF33','P29322','Q5JZY3','P36888','P14616','P35968','P29376','P04629','P07949','P34925','Q6J9G0','Q8NER5','O00238','Q13705','Q16671','Q96Q04','Q8IWU2','Q12866','Q6ZMQ8','Q06418','P35590','P37173','P36897','Q02763','P08922','Q13308','P09619','P16234','Q16832','Q01974','Q01973','Q16288','Q16620','O15146','Q04912','P08581','P10721','P06213','P08069','P35916','P17948','P22455','P21802','P22607','P11362','Q15303','P21860','P04626','O15197','P54760','P54753','P29323','P54762','Q15375','P54756','P54764','P29320','P21709','P29317','P00533','P07333','Q08345','Q13873','P36894','P30530','Q9UM73','P37023','P27037','P36896','Q04771')

#####Let's do######

RCC_result<-survival_anal(RCC,RCC_ID,pt,"RCC_graph")
write.table(RCC_result,"stat_RCC_HR_70RKcandidate.txt",sep="\t",row.names=FALSE)

HCC_result<-survival_anal(HCC,HCC_ID,pt,"HCC_graph")
write.table(HCC_result,"stat_HCC_HR_70RKcandidate.txt",sep="\t",row.names=FALSE)

GBM_result<-survival_anal(GBM,GBM_ID,pt,"GBM_graph")
write.table(GBM_result,"stat_GBM_HR_70RKcandidate.txt",sep="\t",row.names=FALSE)

LUAD_result<-survival_anal(LUAD,LUAD_ID,pt,"LUAD_graph")
write.table(LUAD_result,"stat_LUAD_HR_70RKcandidate.txt",sep="\t",row.names=FALSE)

JHU_result<-survival_anal(JHU,JHU_ID,pt,"JHU_graph")
write.table(JHU_result,"stat_JHU_HR_70RKcandidate.txt",sep="\t",row.names=FALSE)

LUSCC_result<-survival_anal(LUSCC,LUSCC_ID,pt,"LUSCC_graph")
write.table(LUSCC_result,"stat_LUSCC_HR_70RKcandidate.txt",sep="\t",row.names=FALSE)

BRCA_result<-survival_anal(BRCA,BRCA_ID,pt,"BRCA_graph")
write.table(BRCA_result,"stat_BRCA_HR_70RKcandidate.txt",sep="\t",row.names=FALSE)

GA_result<-survival_anal(GA,GA_ID,pt,"GA_graph")
write.table(GA_result,"stat_GA_HR_70RKcandidate.txt",sep="\t",row.names=FALSE)

KU_result<-survival_anal(KU,KU_ID,pt,"KU_graph")
write.table(KU_result,"stat_KU_HR_70RKcandidate.txt",sep="\t",row.names=FALSE)

OV_result<-survival_anal(OV,OV_ID,pt,"OV_graph")
write.table(OV_result,"stat_OV_HR_70RKcandidate.txt",sep="\t",row.names=FALSE)

CRC_result<-survival_anal(CRC,CRC_ID,pt,"CRC_graph")
write.table(CRC_result,"stat_CRC_HR_70RKcandidate.txt",sep="\t",row.names=FALSE)


