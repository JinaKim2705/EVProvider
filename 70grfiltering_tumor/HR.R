BRCA<-read.table("Breast_Cancer_summary_PSMs_grouping_upper70_tumor.txt",sep="\t",header=T)
GA<-read.table("Gastric_Cancer_summary_PSMs_grouping_upper70_tumor.txt",sep="\t",header=T)
KU<-read.table("PDCKU_Cancer_summary_PSMs_grouping_upper70_tumor.txt",sep="\t",header=T)


BRCA_ID<-read.table("Breast_sample_info.txt",sep="\t",header=T)
GA_ID<-read.table("GA_sample_info.txt",sep="\t",header=T)
KU_ID<-read.table("KU_sample_info.txt",sep="\t",header=T)

pt_class<-read.table("C:\\Users\\kimj49\\Box\\Multicancer Biomarker\\receptor_kinase\\proteinclass_receptorkinase.txt",sep="\t",header=T)




####optimal cut off using cutoffr#######
library(survival)
library(Publish)
library(survminer)
library(reshape2)
library(reshape)


#####BRCA_HR#################
data<-BRCA[,c(4:ncol(BRCA))]
rownames(data)<-BRCA[,1]
targets<-BRCA_ID[,c(1,4,5)]


HazardRatio <- sapply(1:nrow(data),function(i){

new<-merge(targets,t(data[i,]),by.x = "Sample" , by.y =0,all.x=TRUE)
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
rownames(HR_result)<-BRCA[,1]
HR_result <- cbind(rownames(HR_result), data.frame(HR_result, row.names=NULL))

df<-merge(pt_class,HR_result,by.x="UNIPROT", by.y="rownames(HR_result)",all.y=TRUE)
new_df<-subset(df,df$class !="NA")


write.table(HR_result,"stat_BRCA_HR_total.txt",sep="\t",row.names=TRUE)
write.table(new_df,"stat_BRCA_HR_RK.txt",sep="\t",row.names=FALSE)


#####KU_HR#################
data<-KU[,c(4:ncol(KU))]
rownames(data)<-KU[,1]
targets<-KU_ID[,c(1,3,4)]


HazardRatio <- sapply(1:nrow(data),function(i){

new<-merge(targets,t(data[i,]),by.x = "Sample" , by.y =0,all.x=TRUE)
new$median <- factor(as.integer(new[,4]>= median(new[,4],na.rm=TRUE)))
model<-coxph(Surv(OS_MONTHS,OS_Death)~median,data=new,na.action = na.omit)
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
rownames(HR_result)<-KU[,1]
HR_result <- cbind(rownames(HR_result), data.frame(HR_result, row.names=NULL))

df<-merge(pt_class,HR_result,by.x="UNIPROT", by.y="rownames(HR_result)",all.y=TRUE)
new_df<-subset(df,df$class !="NA")


write.table(HR_result,"stat_KU_HR_total.txt",sep="\t",row.names=TRUE)
write.table(new_df,"stat_KU_HR_RK.txt",sep="\t",row.names=FALSE)



#####GA_HR#################
data<-GA[,c(4:ncol(GA))]
rownames(data)<-GA[,1]
targets<-GA_ID[,c(1,3,4)]



HazardRatio <- sapply(1:nrow(data),function(i){

new<-merge(targets,t(data[i,]),by.x = "Sample" , by.y =0,all.x=TRUE)
new$median <- factor(as.integer(new[,4]>= median(new[,4],na.rm=TRUE)))
model<-coxph(Surv(OS_MONTHS,OS_Death)~median,data=new,na.action = na.omit)
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
rownames(HR_result)<-GA[,1]
HR_result <- cbind(rownames(HR_result), data.frame(HR_result, row.names=NULL))

df<-merge(pt_class,HR_result,by.x="UNIPROT", by.y="rownames(HR_result)",all.y=TRUE)
new_df<-subset(df,df$class !="NA")


write.table(HR_result,"stat_GA_HR_total.txt",sep="\t",row.names=TRUE)
write.table(new_df,"stat_GA_HR_RK.txt",sep="\t",row.names=FALSE)



