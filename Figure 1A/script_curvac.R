library(ggplot2)

rm(list=ls(all=TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

full_data <- read.csv("data_refinedConv.csv", fileEncoding="UTF-8-BOM") 

#negative log likelihood of normal distibrution model with censoring
Likelihood=function(p,log10_censT,data){-sum(log(dnorm(data[data>log10_censT],p[1],p[2])))-sum(log(pnorm(data[data<=log10_censT],p[1],p[2])))}

data_curevac<-full_data[full_data$Vaccine=='CureVac', ]
curevac_dose<-unique(data_curevac$Dose)

mean_curevac=NULL
sd_curevac=NULL
no_patients_curevac = NULL
mean_conv_adjusted=NULL
sem=NULL
upper=NULL
lower=NULL

#Last row is convalescent
for (i in 1:length(curevac_dose)) 
{
  temp_data<-data_curevac[data_curevac$Dose==curevac_dose[i],]
  Model_Fit<-nlm(function(p){Likelihood(p,log10(5),log10(temp_data$Neut))},c(mean(log10(temp_data$Neut)),sd(log10(temp_data$Neut))))    
  mean_curevac[i]<-Model_Fit$estimate[1]
  sd_curevac[i]<-Model_Fit$estimate[2]
  no_patients_curevac[i]<-length(temp_data$Neut)
}

for (i in 1:length(mean_curevac)-1) 
{
  mean_conv_adjusted[i]<-mean_curevac[i] - mean_curevac[length(mean_curevac)]
  sem[i]<-sqrt( (sd_curevac[i]^2)/no_patients_curevac[i] + (sd_curevac[length(mean_curevac)]^2)/no_patients_curevac[length(mean_curevac)] )
  upper[i]<-mean_conv_adjusted[i]+1.96*sem[i]
  lower[i]<-mean_conv_adjusted[i]-1.96*sem[i]
}

curevac_estimate<-data.frame("vaccine"="CVnCoV","dose"=curevac_dose[1:(length(curevac_dose)-1)],mean_conv_adjusted,upper,lower)
# write.csv(curevac_estimate,"curevac_estimate.csv") #To plot in Prism



################## PFIZER ##############################################
data_pfizer<-full_data[full_data$Vaccine=='Pfizer', ]
pfizer_dose<-unique(data_pfizer$Dose)

mean_pfizer=NULL
sd_pfizer=NULL
no_patients_pfizer = NULL
mean_conv_adjusted=NULL
sem=NULL
upper=NULL
lower=NULL

#Last row is convalescent
for (i in 1:length(pfizer_dose)) 
{
  temp_data<-data_pfizer[data_pfizer$Dose==pfizer_dose[i],]
  Model_Fit<-nlm(function(p){Likelihood(p,log10(10),log10(temp_data$Neut))},c(mean(log10(temp_data$Neut)),sd(log10(temp_data$Neut))))    
  mean_pfizer[i]<-Model_Fit$estimate[1]
  sd_pfizer[i]<-Model_Fit$estimate[2]
  no_patients_pfizer[i]<-length(temp_data$Neut)
}

for (i in 1:length(mean_pfizer)-1) 
{
  mean_conv_adjusted[i]<-mean_pfizer[i] - mean_pfizer[length(mean_pfizer)]
  sem[i]<-sqrt( (sd_pfizer[i]^2)/no_patients_pfizer[i] + (sd_pfizer[length(mean_pfizer)]^2)/no_patients_pfizer[length(mean_pfizer)] )
  upper[i]<-mean_conv_adjusted[i]+1.96*sem[i]
  lower[i]<-mean_conv_adjusted[i]-1.96*sem[i]
}

pfizer_estimate<-data.frame("vaccine"="BNT162b2","dose"=pfizer_dose[1:(length(pfizer_dose)-1)],mean_conv_adjusted,upper,lower)
# write.csv(pfizer_estimate,"pfizer_estimate.csv") #To plot in Prism

######################## MODERNA ##################################################

#negative log likelihood of normal distibrution model without censoring
Likelihood=function(p,data){-sum(log(dnorm(data,p[1],p[2])))}

data_moderna<-full_data[full_data$Vaccine=='Moderna', ]
moderna_dose<-unique(data_moderna$Dose)

mean_moderna=NULL
sd_moderna=NULL
no_patients_moderna = NULL
mean_conv_adjusted=NULL
sem=NULL
upper=NULL
lower=NULL

#Last row is convalescent
for (i in 1:length(moderna_dose)) 
{
  temp_data<-data_moderna[data_moderna$Dose==moderna_dose[i],]
  Model_Fit<-nlm(function(p){Likelihood(p,log10(temp_data$Neut))},c(mean(log10(temp_data$Neut)),sd(log10(temp_data$Neut))))    
  mean_moderna[i]<-Model_Fit$estimate[1]
  sd_moderna[i]<-Model_Fit$estimate[2]
  no_patients_moderna[i]<-length(temp_data$Neut)
}

for (i in 1:length(mean_moderna)-1) 
{
  mean_conv_adjusted[i]<-mean_moderna[i] - mean_moderna[length(mean_moderna)]
  sem[i]<-sqrt( (sd_moderna[i]^2)/no_patients_moderna[i] + (sd_moderna[length(mean_moderna)]^2)/no_patients_moderna[length(mean_moderna)] )
  upper[i]<-mean_conv_adjusted[i]+1.96*sem[i]
  lower[i]<-mean_conv_adjusted[i]-1.96*sem[i]
}

moderna_estimate<-data.frame("vaccine"="mRNA-1273","dose"=moderna_dose[1:(length(moderna_dose)-1)],mean_conv_adjusted,upper,lower)
# write.csv(moderna_estimate,"moderna_estimate.csv") #To plot in Prism

combined_estimate=rbind(curevac_estimate,pfizer_estimate,moderna_estimate)


combined_estimate$vaccine=factor(combined_estimate$vaccine,levels=c("CVnCoV","BNT162b2","mRNA-1273"))

colorlist=c("black","mediumpurple","red")
shapelist=c(19,15,18)
###Plot Figure 1A
Figure1A<-ggplot(data=combined_estimate,aes(x=dose,y=10^mean_conv_adjusted,color=vaccine))+
  geom_point(aes(shape=vaccine))+
  geom_line()+
  geom_errorbar(aes(ymin=10^lower,ymax=10^upper)) +
  geom_hline(yintercept = 1,linetype=3)+
  coord_cartesian(expand=FALSE,xlim=c(1,110)) +
  scale_x_log10(lim=c(1,120),breaks=c(2,4,6,8,10,12,20,25,30,100))+
  scale_y_log10(lim=c(0.0625,11),breaks=c(0.0625,0.125,0.25,0.5,1,2,4,8),labels=c(0.0625,0.125,0.25,0.5,1,2,4,8))+
  scale_color_manual(values=colorlist)+
  scale_shape_manual(values=shapelist)+
  theme_linedraw() +
  labs(y="Neutralisation level\n(fold of convalescencent plasma)",
       x="Dose") +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=45,vjust=1,hjust=1),
        axis.title.x =element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing.y = unit(-1.5, "lines"),
        panel.spacing.x = unit(-1, "lines"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key = element_rect(size = 1),
        legend.key.size = unit(0.45, "cm"),
        # legend.margin = margin(t=0,b=0.7,unit="cm"),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black"))


pdf("Figure1A.pdf",height=3.5,width=6)
Figure1A
dev.off()


################# T test for difference ######################################################
data_pfizer_dose10<-log10(data_pfizer[data_pfizer$Dose==10,]$Neut)
data_pfizer_conv<-log10(data_pfizer[data_pfizer$Dose==0,]$Neut)

data_curevac_dose12<-log10(data_curevac[data_curevac$Dose==12,]$Neut)
data_curevac_conv<-log10(data_curevac[data_curevac$Dose==0,]$Neut)

pfizer_diff_with_cov_ave<-data_pfizer_dose10-mean(data_pfizer_conv)
curevac_diff_with_cov_ave<-data_curevac_dose12-mean(data_curevac_conv)
test<-t.test(10^(pfizer_diff_with_cov_ave),10^(curevac_diff_with_cov_ave))
test$p.value #p=0.76

################### Fitting pfizer and curevac low dose ######################################

Likelihood_different_mean=function(p,data_conv_pfizer,data_vaccine_pfizer,data_conv_curevac,data_vaccine_curevac)
{
  LL_conv_pfizer = -sum(log(dnorm(data_conv_pfizer[data_conv_pfizer>log10(10)],p[1],p[2])))-sum(log(pnorm(data_conv_pfizer[data_conv_pfizer<=log10(10)],p[1],p[2])))
  LL_vaccine_pfizer = -sum(log(dnorm(data_vaccine_pfizer[data_vaccine_pfizer>log10(10)],p[1]*p[3],p[4])))-sum(log(pnorm(data_vaccine_pfizer[data_vaccine_pfizer<=log10(10)],p[1]*p[3],p[4])))
  LL_conv_curevac = -sum(log(dnorm(data_conv_curevac[data_conv_curevac>log10(5)],p[5],p[6])))-sum(log(pnorm(data_conv_curevac[data_conv_curevac<=log10(5)],p[5],p[6])))
  LL_vaccine_curevac = -sum(log(dnorm(data_vaccine_curevac[data_vaccine_curevac>log10(5)],p[5]*p[7],p[8])))-sum(log(pnorm(data_vaccine_curevac[data_vaccine_curevac<=log10(5)],p[5]*p[7],p[8])))
  LL_Total = LL_conv_pfizer+LL_vaccine_pfizer+LL_conv_curevac+LL_vaccine_curevac
}

Likelihood_same_mean=function(p,data_conv_pfizer,data_vaccine_pfizer,data_conv_curevac,data_vaccine_curevac)
{
  LL_conv_pfizer = -sum(log(dnorm(data_conv_pfizer[data_conv_pfizer>log10(10)],p[1],p[2])))-sum(log(pnorm(data_conv_pfizer[data_conv_pfizer<=log10(10)],p[1],p[2])))
  LL_vaccine_pfizer = -sum(log(dnorm(data_vaccine_pfizer[data_vaccine_pfizer>log10(10)],p[1]*p[3],p[4])))-sum(log(pnorm(data_vaccine_pfizer[data_vaccine_pfizer<=log10(10)],p[1]*p[3],p[4])))
  LL_conv_curevac = -sum(log(dnorm(data_conv_curevac[data_conv_curevac>log10(5)],p[5],p[6])))-sum(log(pnorm(data_conv_curevac[data_conv_curevac<=log10(5)],p[5],p[6])))
  LL_vaccine_curevac = -sum(log(dnorm(data_vaccine_curevac[data_vaccine_curevac>log10(5)],p[5]*p[3],p[7])))-sum(log(pnorm(data_vaccine_curevac[data_vaccine_curevac<=log10(5)],p[5]*p[3],p[7])))
  LL_Total = LL_conv_pfizer+LL_vaccine_pfizer+LL_conv_curevac+LL_vaccine_curevac
}

Fit_different_mean<-nlm(function(p){Likelihood_different_mean(p,data_pfizer_conv,data_pfizer_dose10,data_curevac_conv,data_curevac_dose12 )},c(mean(data_pfizer_conv),sd(data_pfizer_conv),1,0.2,mean(data_curevac_conv),sd(data_curevac_conv),1,0.2))    
Fit_same_mean<-nlm(function(p){Likelihood_same_mean(p,data_pfizer_conv,data_pfizer_dose10,data_curevac_conv,data_curevac_dose12 )},c(mean(data_pfizer_conv),sd(data_pfizer_conv),1,0.2,mean(data_curevac_conv),sd(data_curevac_conv),0.2))    

#Minus the minimum for the likelihood value
llr = 2*(-Fit_different_mean$minimum - -Fit_same_mean$minimum)
p_val_LLR_test<-pchisq(llr, df=1, lower.tail=FALSE) #p=0.18