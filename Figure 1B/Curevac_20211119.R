library(ggplot2)
library(plyr)
library('maxLik')
library(magic)

#### Set directory to same directory as the r-script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


### Import Tables Required for Estimating Mean and SD for each study
SummaryTable_Efficacy_NeutRatio_SD_SEM<-read.csv("SummaryTable_Efficacy_NeutRatio_SD_SEM.csv")
SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported=log10(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutMean/SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutConv)

VariantData<-read.csv("censoredChangeMeans_Var_20210812.csv")


####Reproduce fits from Nature Med (Khoury et al., 2021)
##########################################################################################
##################### The Models

####Logistic Model
ProbRemainUninfected=function(logTitre,logk,C50){1/(1+exp(-exp(logk)*(logTitre-C50)))}

LogisticModel_PercentUninfected=function(mu_titre,sig_titre,logk,C50){
  Output<-NULL

  if (length(C50)==1) {
    C50=rep(C50,length(mu_titre))
  }

  if (length(logk)==1) {
    logk=rep(logk,length(mu_titre))
  }

  for (i in 1:length(mu_titre)) {
    Step=sig_titre[i]*0.001
    IntegralVector=seq(mu_titre[i]-5*sig_titre[i],mu_titre[i]+5*sig_titre[i],by=Step)
    Output[i]=sum(ProbRemainUninfected(IntegralVector,logk[i],C50[i])*dnorm(IntegralVector,mu_titre[i],sig_titre[i]))*Step
  }
  Output
}


### Logistic model for Raw Efficacy Counts
FittingLogistic_Raw<-function(logRisk0,logk,C50,N_C,N_V,Inf_C,Inf_V,MeanVector,SDVector){

  Risk0=exp(logRisk0)

  if (length(SDVector)==1) {
    SDVector=rep(SDVector,length(N_C))
  }

  IndexNA=(is.na(N_C) | is.na(MeanVector) | is.na(SDVector))
  N_C=N_C[!IndexNA]
  N_V=N_V[!IndexNA]
  Inf_V=Inf_V[!IndexNA]
  Inf_C=Inf_C[!IndexNA]
  MeanVector=MeanVector[!IndexNA]
  SDVector=SDVector[!IndexNA]

  if (length(C50)==1) {
    C50=rep(C50,length(N_C))
  }

  if (length(logk)==1) {
    logk=rep(logk,length(N_C))
  }

  LL=0
  for (i in 1:length(N_C)) {

    LL=LL-log(dbinom(Inf_C[i],N_C[i],Risk0[i]))-log(dbinom(Inf_V[i],N_V[i],Risk0[i]*(1-LogisticModel_PercentUninfected(MeanVector[i],SDVector[i],logk[i],C50[i]))))
  }
  LL
}


#################################################################
####################################### Fit from Nature Med.
LogisticEstimate=c("logk"=log(2.7),"C50"=log10(0.5))
FittedLogistic_RawEfficacy_MeanRept_SDPool<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                                                                SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate),hessian=TRUE)

######Plotting
#Error bars
SummaryTable_Efficacy_NeutRatio_SD_SEM$RatioReported_LB=10^((SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported)-1.96*SummaryTable_Efficacy_NeutRatio_SD_SEM$SEM)
SummaryTable_Efficacy_NeutRatio_SD_SEM$RatioReported_UB=10^((SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported)+1.96*SummaryTable_Efficacy_NeutRatio_SD_SEM$SEM)

#Line from fitted model
NeutValue=seq(0.03,11,by=0.001)

Efficacy_Logistic_Raw<-NULL

for (i in 1:length(NeutValue)) {
  
  Efficacy_Logistic_Raw[i]=LogisticModel_PercentUninfected(log10(NeutValue[i]),SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1))
}

LogisticModel_withPoolSD=data.frame("NeutRatio_Reported"=log10(NeutValue),"Efficacy"=Efficacy_Logistic_Raw)
LogisticModel_withPoolSD$Study<-rep("LogisticModel",length(NeutValue))


##### Adding CureVac
CurevacTable<-SummaryTable_Efficacy_NeutRatio_SD_SEM[1,]
CurevacTable[1:nrow(VariantData),]<-NA
CurevacTable$Study<-"Curevac"
CurevacTable$TechnicalName<-"CVnCoV"
CurevacTable$Age<-"all"
CurevacTable$NeutMean<-92.93534
CurevacTable$NeutConv<-145.34387
CurevacTable$NeutRatio_Reported<-log10(CurevacTable$NeutMean/CurevacTable$NeutConv)
CurevacTable$RatioReported_LB<-CurevacTable$NeutRatio_Reported
CurevacTable$RatioReported_UB<-CurevacTable$NeutRatio_Reported
CurevacTable$Efficacy<-0.482#Paper
Curevac_PointEstimate=LogisticModel_PercentUninfected(CurevacTable$NeutRatio_Reported[1],SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1))  
CurevacTable$Variant<-factor(c("Ancestral","Alpha (B.1.1.7)","Beta (B.1.351)","Gamma (P.1)","Delta (B.1.617.2)"),levels=(c("Ancestral","Alpha (B.1.1.7)","Beta (B.1.351)","Gamma (P.1)","Delta (B.1.617.2)")))
CurevacTable$NeutRatio_Reported[CurevacTable$Variant=="Alpha (B.1.1.7)"]<-CurevacTable$NeutRatio_Reported[1]+VariantData$muL_changeFromWTL_Var[VariantData$Variant=="B.1.1.7"]
CurevacTable$NeutRatio_Reported[CurevacTable$Variant=="Beta (B.1.351)"]<-CurevacTable$NeutRatio_Reported[1]+VariantData$muL_changeFromWTL_Var[VariantData$Variant=="B.1.351"]
CurevacTable$NeutRatio_Reported[CurevacTable$Variant=="Gamma (P.1)"]<-CurevacTable$NeutRatio_Reported[1]+VariantData$muL_changeFromWTL_Var[VariantData$Variant=="P.1"]
CurevacTable$NeutRatio_Reported[CurevacTable$Variant=="Delta (B.1.617.2)"]<-CurevacTable$NeutRatio_Reported[1]+VariantData$muL_changeFromWTL_Var[VariantData$Variant=="B.1.617.2"]
CurevacTable$Lower<-31.0#Paper
CurevacTable$Upper<-61.4 #Paper

### Confidence bounds:
Cov<-solve(FittedLogistic_RawEfficacy_MeanRept_SDPool$hessian)[9:10,9:10]

#################################################### Adding 95% Prediction Intervals #########################################################
grad1<-NULL
grad2<-NULL
Lower_Pred<-NULL
Upper_Pred<-NULL
G<-NULL
Lower_Pred_Variant<-NULL
Upper_Pred_Variant<-NULL
Lower_Pred_ALL<-NULL
Upper_Pred_ALL<-NULL


for (i in 1:length(NeutValue))
{
  f_temp <- function(p_temp) LogisticModel_PercentUninfected(log10(NeutValue[i]),SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],p_temp[1],p_temp[2])
  grad1[i]<-numericGradient(f_temp, c(tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1)))[1]
  grad2[i]<-numericGradient(f_temp, c(tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1)))[2]
  G<-cbind(grad1[i],grad2[i])
  Lower_Pred[i]=Efficacy_Logistic_Raw[i]-1.96*sqrt(G%*%Cov%*%t(G))
  Upper_Pred[i]=Efficacy_Logistic_Raw[i]+1.96*sqrt(G%*%Cov%*%t(G))
  
}

LogisticModel_withPoolSD$Lower<-100*Lower_Pred
LogisticModel_withPoolSD$Upper<-100*Upper_Pred



group.colors=(c("gray27","navyblue","mediumpurple","dodgerblue","darkorange2")) #,"yellow","green","gray","black","purple")
group.colors.fill=(c("gray27","white","white","white","white"))


Figure1B<-ggplot(data=SummaryTable_Efficacy_NeutRatio_SD_SEM, aes(y=100*Efficacy,x=(10^NeutRatio_Reported),group=Study)) +
  #adding the bands
  annotate(geom="rect",xmin=10^CurevacTable$NeutRatio_Reported[CurevacTable$Variant=="Beta (B.1.351)"],xmax=10^CurevacTable$NeutRatio_Reported[CurevacTable$Variant=="Alpha (B.1.1.7)"],ymin=CurevacTable$Lower[1],ymax=CurevacTable$Upper[1],fill="yellow",alpha=0.2) +
  geom_ribbon(data=LogisticModel_withPoolSD,aes(ymin=Lower, ymax=Upper),fill="pink", alpha = 0.3)+
  geom_point(shape=21,fill="gray",color="gray") +
  geom_errorbar(aes(ymin=Lower,ymax=Upper),color="gray") +
  geom_errorbarh(aes(xmin=RatioReported_LB,xmax=RatioReported_UB),color="gray") +
  scale_x_log10(breaks=c(0.0625,0.125,0.25,0.5,1,2,4,8),labels=c(0.0625,0.125,0.25,0.5,1,2,4,8)) +
  scale_y_continuous() +
  coord_cartesian(ylim=c(20,100),xlim=c(0.0625/1.2,11)) +
  
  theme_linedraw() +
  geom_line(data=LogisticModel_withPoolSD,color="red") +

  # geom_point(data=CovaxinTable,aes(color=Study),size=4,shape=1) +
  geom_text(data=SummaryTable_Efficacy_NeutRatio_SD_SEM,aes(y=100*Efficacy-2*(0.5-1*(Efficacy<0.7))*3,x=(10^NeutRatio_Reported)*(1+0.1*2*(0.5-1*(Efficacy<0.7))),label=TechnicalName,hjust=1*(Efficacy<0.7),vjust=1*(Efficacy<0.7)),size=3,angle=-45,color="gray") +
  geom_hline(yintercept=100*CurevacTable$Efficacy[1],color="gray27",show.legend=FALSE,linetype=2) +
  geom_errorbar(data=CurevacTable,aes(ymin=Lower,ymax=Upper,color=Variant),width=0.05) +
  geom_point(data=CurevacTable,aes(color=Variant,fill=Variant),shape=21,size=4) +
  annotate("text",x=7, y=50,label=CurevacTable$TechnicalName[1]) +
  xlab("Neutralisation level") +
  ylab("Protective efficacy (%)") +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors.fill)



pdf("Figure1B.pdf",height=3.5,width=6)
Figure1B
dev.off()
Figure1B
