library(grid)
library(ggpubr)
library(gridExtra)
source("Functions.R")
#import data
ChiF<-read.csv("bat_839.csv",header=T)

Chi_precip=cbind(rep(1,839),ChiF$precip_mean_mm) #x matrix
Chi_forearm=ChiF$adult_forearm_length_mm #response-y

#G from 2-8,fit GM, GMR, Homo-GMR
#Get LRT, BIC

#G=2
#alp_int<-c(0.6,0.4)
alp_int<-c(0.75,0.25)
#beta_int<-cbind(c(30,0),c(40,0.5)) try different sets of initial values
beta_int<-cbind(c(30,0),c(41,0.5))
sig_int<-10
mu_int<-c(30,100)
#mu_int<-c(20,80)
MG2<-EMGuasMix(Chi_forearm,alp_int,Chi_precip,beta_int,sig_int)
MG2N<-EMGuasMix_NoX(Chi_forearm,alp_int,mu_int,sig_int)
BIC2<-BICem(MG2)
BIC2N<-BICem(MG2N)
Lrt2<-LRT(MG2,MG2N)
beta_int<-cbind(c(30,0.1),c(100,0.1))
#beta_int<-cbind(c(50,0),c(150,0))
MG2S<-EMGuasMix_fix1(Chi_forearm,alp_int,Chi_precip,beta_int,sig_int)
BIC2S<-BICem_fix(MG2S)
Lrt2S<-LRT_fix(MG2,MG2S)

#G=3
#alp_int<-c(0.1,0.3,0.6)
alp_int<-c(0.2,0.25,0.55)
#beta_int<-cbind(c(50,0.1),c(20,0.3),c(10,0))
beta_int<-cbind(c(50,0.5),c(20,0.3),c(30,0))
sig_int<-20
#sig_int<-10
MG3<-EMGuasMix(Chi_forearm,alp_int,Chi_precip,beta_int,sig_int)
#mu_int<-c(100,50,30)
mu_int<-c(150,100,50)
MG3N<-EMGuasMix_NoX(Chi_forearm,alp_int,mu_int,sig_int)
BIC3<-BICem(MG3)
BIC3N<-BICem(MG3N)
Lrt3<-LRT(MG3,MG3N)
#beta_int<-cbind(c(1000,0),c(80,0),c(50,0))
beta_int<-cbind(c(150,0),c(100,0),c(50,0))
MG3S<-EMGuasMix_fix1(Chi_forearm,alp_int,Chi_precip,beta_int,sig_int)
BIC3S<-BICem_fix(MG3S)
Lrt3S<-LRT_fix(MG3,MG3S)

#G=4
alp_int<-c(0.2,0.2,0.45,0.15)
beta_int<-cbind(c(50,0.5),c(20,0.3),c(30,0),c(100,0))
#beta_int<-cbind(c(50,0.5),c(20,0.3),c(30,0),c(218,-0.5))
#beta_int<-cbind(c(50,0.5),c(20,0.3),c(60,0),c(168,-0.1))
#beta_int<-cbind(c(40,0.7),c(30,1.5),c(30,0),c(60,0))
sig_int<-20
MG4<-EMGuasMix(Chi_forearm,alp_int,Chi_precip,beta_int,sig_int)
alp_int<-c(0.1,0.2,0.25,0.45)
mu_int<-c(150,100,80,30)
#mu_int<-c(100,80,50,30)
MG4N<-EMGuasMix_NoX(Chi_forearm,alp_int,mu_int,sig_int)
BIC4<-BICem(MG4)
BIC4N<-BICem(MG4N)
Lrt4<-LRT(MG4,MG4N)
beta_int<-cbind(c(150,0),c(100,0),c(80,0),c(50,0))
#beta_int<-cbind(c(100,0),c(80,0),c(50,0),c(30,0))
MG4S<-EMGuasMix_fix1(Chi_forearm,alp_int,Chi_precip,beta_int,sig_int)
BIC4S<-BICem_fix(MG4S)
Lrt4S<-LRT_fix(MG4,MG4S)

#G=5
alp_int<-c(0.2,0.25,0.35,0.1,0.1)
#beta_int<-cbind(c(70,0.5),c(50,0.3),c(30,0),c(120,0.1),c(100,0))
beta_int<-cbind(c(40,0.7),c(30,1.5),c(30,0),c(120,0.1),c(60,0))
#beta_int<-cbind(c(40,0.2),c(60,-0.05),c(30,0),c(40,0.02),c(60,0.3))
#beta_int<-cbind(c(40,0.7),c(30,1.5),c(30,0),c(120,0.1),c(220,-0.5))
sig_int<-20
MG5<-EMGuasMix(Chi_forearm,alp_int,Chi_precip,beta_int,sig_int)
alp_int<-c(0.1,0.1,0.15,0.3,0.35)
#alp_int<-c(0.1,0.05,0.05,0.1,0.7)
#mu_int<-c(180,150,100,80,30)
mu_int<-c(180,150,120,80,30)
#mu_int<-c(180,150,80,60,30)
MG5N<-EMGuasMix_NoX(Chi_forearm,alp_int,mu_int,sig_int)
BIC5<-BICem(MG5)
BIC5N<-BICem(MG5N)
Lrt5<-LRT(MG5,MG5N)
#beta_int<-cbind(c(180,0),c(150,0),c(100,0),c(80,0),c(40,0))
beta_int<-cbind(c(180,0),c(150,0),c(120,0),c(80,0),c(30,0))
MG5S<-EMGuasMix_fix1(Chi_forearm,alp_int,Chi_precip,beta_int,sig_int)
BIC5S<-BICem_fix(MG5S)
Lrt5S<-LRT_fix(MG5,MG5S)


#G=6
alp_int<-c(0.2,0.15,0.35,0.1,0.1,0.1)
#alp_int<-c(0.2,0.2,0.35,0.1,0.1,0.05)
#beta_int<-cbind(c(40,0.7),c(210,-0.5),c(30,0),c(120,0.1),c(60,0),c(30,1.5))
beta_int<-cbind(c(40,0.7),c(30,1.5),c(30,0),c(120,0.1),c(60,0),c(106.5,0.16))
sig_int<-20
MG6<-EMGuasMix(Chi_forearm,alp_int,Chi_precip,beta_int,sig_int)
alp_int<-c(0.1,0.1,0.15,0.15,0.25,0.25)
mu_int<-c(180,150,120,100,80,30)
#beta_int<-cbind(c(200,0),c(100,0),c(60,0),c(80,0),c(40,0))
MG6N<-EMGuasMix_NoX(Chi_forearm,alp_int,mu_int,sig_int)
BIC6<-BICem(MG6)
BIC6N<-BICem(MG6N)
Lrt6<-LRT(MG6,MG6N)
beta_int<-cbind(c(180,0),c(150,0),c(120,0),c(100,0),c(80,0),c(30,0))
#beta_int<-cbind(c(180,0),c(150,0),c(120,0),c(100,0),c(60,0),c(30,0))
MG6S<-EMGuasMix_fix1(Chi_forearm,alp_int,Chi_precip,beta_int,sig_int)
BIC6S<-BICem_fix(MG6S)
Lrt6S<-LRT_fix(MG6,MG6S)

#G=7
alp_int<-c(0.05,0.05,0.1,0.1,0.15,0.25,0.3)
#mu_int<-c(200,150,130,120,100,50,30)
mu_int<-c(180,150,130,120,100,80,30)
MG7N<-EMGuasMix_NoX(Chi_forearm,alp_int,mu_int,sig_int)
alp_int<-c(0.2,0.15,0.35,0.1,0.1,0.05,0.05)
#beta_int<-cbind(c(40,0.7),c(30,1.5),c(30,0),c(120,0.1),c(60,0),c(106.5,0.16),c(210,-0.5))
beta_int<-cbind(c(40,0.7),c(30,1.5),c(30,0),c(120,0.1),c(60,0),c(106.5,0.16),c(80,0.1))
MG7<-EMGuasMix(Chi_forearm,alp_int,Chi_precip,beta_int,sig_int)
BIC7<-BICem(MG7)
BIC7N<-BICem(MG7N)
Lrt7<-LRT(MG7,MG7N)
alp_int<-c(0.05,0.05,0.1,0.1,0.15,0.25,0.3)
beta_int<-cbind(c(180,0),c(150,0),c(130,0),c(120,0),c(100,0),c(80,0),c(30,0))
#beta_int<-cbind(c(180,0),c(150,0),c(120,0),c(100,0),c(60,0),c(30,0))
MG7S<-EMGuasMix_fix1(Chi_forearm,alp_int,Chi_precip,beta_int,sig_int)
BIC7S<-BICem_fix(MG7S)
Lrt7S<-LRT_fix(MG7,MG7S)

#G=8
#alp_int<-c(0.2,0.15,0.35,0.1,0.08,0.05,0.05,0.02)
alp_int<-c(0.05,0.05,0.05,0.05,0.1,0.15,0.25,0.3)
mu_int<-c(180,150,130,120,100,90,80,30)
#mu_int<-c(200,150,130,120,100,80,50,30)
MG8N<-EMGuasMix_NoX(Chi_forearm,alp_int,mu_int,sig_int)
BIC8N<-BICem(MG8N)
alp_int<-c(0.2,0.15,0.35,0.1,0.08,0.05,0.05,0.02)
beta_int<-cbind(c(40,0.7),c(30,1.5),c(30,0),c(120,0.1),c(60,0),c(106.5,0.16),c(80,0.1),c(50,-0.5))
#beta_int<-cbind(c(40,0.7),c(30,1.5),c(30,0),c(120,0.1),c(60,0),c(106.5,0.16),c(80,0.1),c(210,-0.5))
#beta_int<-cbind(c(40,0.7),c(30,1.5),c(30,0),c(120,0.1),c(60,0),c(106.5,0.16),c(80,0.1),c(100,0))
MG8<-EMGuasMix(Chi_forearm,alp_int,Chi_precip,beta_int,sig_int)
BIC8<-BICem(MG8)
BIC8N<-BICem(MG8N)
Lrt8<-LRT(MG8,MG8N)
alp_int<-c(0.05,0.05,0.05,0.05,0.1,0.15,0.25,0.3)
beta_int<-cbind(c(180,0),c(150,0),c(130,0),c(120,0),c(100,0),c(90,0),c(80,0),c(30,0))
#beta_int<-cbind(c(180,0),c(150,0),c(120,0),c(100,0),c(60,0),c(30,0))
MG8S<-EMGuasMix_fix1(Chi_forearm,alp_int,Chi_precip,beta_int,sig_int)
BIC8S<-BICem_fix(MG8S)
Lrt8S<-LRT_fix(MG8,MG8S)


#Naive I-overall test
P_BIC=cbind(BIC2N,BIC2,BIC3N,BIC3,BIC4N,BIC4,BIC5N,BIC5,BIC6N,BIC6,BIC7N,BIC7,BIC8N,BIC8)
P_BICC=apply(P_BIC,1,which.min)
P_BICC 
P_BIC
Lrt5

#Naive II-overall test
P_BIC=cbind(BIC2N,BIC3N,BIC4N,BIC5N,BIC6N,BIC7N,BIC8N)
P_BICC=apply(P_BIC,1,which.min)
P_BICC 
P_BIC
Lrt7

#Naive III-overall test
cbind(Lrt2,Lrt3,Lrt4,Lrt5,Lrt6,Lrt7,Lrt8)
P_BIC=cbind(BIC2,BIC3,BIC4,BIC5,BIC6,BIC7,BIC8)
P_BICC=apply(P_BIC,1,which.min)
P_BICC 
Lrt5

#WEST - overall test
T_BIC=cbind(BIC2N,BIC3N,BIC4N,BIC5N,BIC6N,BIC7N,BIC8N)
colnames(T_BIC)=c("BIC2N","BIC3N","BIC4N","BIC5N","BIC6N","BIC7N","BIC8N")
nocol=ncol(T_BIC)
norow=nrow(T_BIC)
PI_all=matrix(0,norow,nocol) 
for(i in 1:norow ){
  for(j in 1:nocol ){
    PI_all[i,j]<-1/(sum(exp(-1/2*(T_BIC[i,]-T_BIC[i,j]))))
  }}
T_chipvalues=cbind(Lrt2[[4]],Lrt3[[4]],Lrt4[[4]],Lrt5[[4]],Lrt6[[4]],Lrt7[[4]],Lrt8[[4]])
colnames(T_chipvalues)=c("BIC2N","BIC3N","BIC4N","BIC5N","BIC6N","BIC7N","BIC8N")
T_Pva=apply(PI_all*T_chipvalues,1,sum)  

#Naive I-heterogeneous test
P_BIC=cbind(BIC2S,BIC2,BIC3S,BIC3,BIC4S,BIC4,BIC5S,BIC5,BIC6S,BIC6,BIC7S,BIC7,BIC8S,BIC8)
P_BICC=apply(P_BIC,1,which.min)
P_BICC 
P_BIC
Lrt5S

#Naive II-heterogeneous test
P_BIC=cbind(BIC2S,BIC3S,BIC4S,BIC5S,BIC6S,BIC7S,BIC8S)
P_BICC=apply(P_BIC,1,which.min)
P_BICC 
P_BIC
Lrt7S

#Naive III-heterogeneous test
cbind(Lrt2S,Lrt3S,Lrt4S,Lrt5S,Lrt6S,Lrt7S,Lrt8S)
P_BIC=cbind(BIC2,BIC3,BIC4,BIC5,BIC6,BIC7,BIC8)
P_BICC=apply(P_BIC,1,which.min)
P_BICC 
Lrt5S

#WEST-heterogeneous test
T_BIC=cbind(BIC2S,BIC3S,BIC4S,BIC5S,BIC6S,BIC7S,BIC8S)
colnames(T_BIC)=c("BIC2S","BIC3S","BIC4S","BIC5S","BIC6S","BIC7S","BIC8S")
nocol=ncol(T_BIC)
norow=nrow(T_BIC)
PI_all=matrix(0,norow,nocol) 
for(i in 1:norow ){
  for(j in 1:nocol ){
    PI_all[i,j]<-1/(sum(exp(-1/2*(T_BIC[i,]-T_BIC[i,j]))))
  }}
T_chipvalues=cbind(Lrt2S[[4]],Lrt3S[[4]],Lrt4S[[4]],Lrt5S[[4]],Lrt6S[[4]],Lrt7S[[4]],Lrt8S[[4]])
colnames(T_chipvalues)=c("BIC2S","BIC3S","BIC4S","BIC5S","BIC6S","BIC7S","BIC8S")
T_Pva=apply(PI_all*T_chipvalues,1,sum) 

#Creat plots
data_em1<-as.data.frame(cbind(ChiF$precip_mean_mm,ChiF$adult_forearm_length_mm,as.factor(MG7N[4][[1]])))
colnames(data_em1)<-c("precip","forearm","group")
data_em1$group=as.factor(data_em1$group)
p1<-ggplot(data=data_em1) +
  geom_point(aes(precip,forearm, colour=group))+
  #geom_text(aes(sspg,insulin,label=class,colour=group))+
  ggtitle("a) 7-GM")+
  labs(x= "Precip (mm)", y = "Forearm (mm)")+
  scale_colour_manual(values=c("#F0E442","#0072B2", "#E69F00","#999999","#CC79A7","#009E73", "#56B4E9"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") 

data_em1<-as.data.frame(cbind(ChiF$precip_mean_mm,ChiF$adult_forearm_length_mm,as.factor(MG7S[4][[1]])))
colnames(data_em1)<-c("precip","forearm","group")
data_em1$group=as.factor(data_em1$group)
p4<-ggplot(data=data_em1) +
  geom_point(aes(precip,forearm, colour=group))+
  #geom_text(aes(sspg,insulin,label=class,colour=group))+
  geom_abline(intercept = 186.46171720, slope = 0.01288193 , color="#F0E442", size=0.3)+
  geom_abline(intercept = 161.07538486, slope = 0.01288193  , color="#0072B2", size=0.3)+
  geom_abline(intercept = 134.26884637, slope = 0.01288193  , color="#E69F00", size=0.3)+
  geom_abline(intercept = 111.09349372, slope = 0.01288193 , color="#999999", size=0.3)+
  geom_abline(intercept = 79.72572516, slope = 0.01288193 , color="#CC79A7", size=0.3)+
  geom_abline(intercept = 57.83839984, slope = 0.01288193 , color="#009E73", size=0.3)+
  geom_abline(intercept = 38.36792375, slope = 0.01288193 , color="#56B4E9", size=0.3)+
  ggtitle("b) 7-Homo-GMR")+
  labs(x= "Precip (mm)", y = "Forearm (mm)")+
  scale_colour_manual(values=c("#F0E442","#0072B2", "#E69F00","#999999","#CC79A7","#009E73", "#56B4E9"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") 


data_em1<-as.data.frame(cbind(ChiF$precip_mean_mm,ChiF$adult_forearm_length_mm,as.factor(MG5[4][[1]])))
colnames(data_em1)<-c("precip","forearm","group")
data_em1$group=as.factor(data_em1$group)
p2<-ggplot(data=data_em1) +
  geom_point(aes(precip,forearm, colour=group))+
  #geom_text(aes(sspg,insulin,label=class,colour=group))+
  geom_abline(intercept = 146.9587056, slope = 0.1655061, color="#F0E442", size=0.3)+
  geom_abline(intercept = 106.5020481, slope = 0.1586673, color="#0072B2", size=0.3)+
  geom_abline(intercept = 60.0266181, slope = 0.2318351 , color="#E69F00", size=0.3)+
  geom_abline(intercept = 47.3913316, slope = 0.1232126, color="#999999", size=0.3)+
  geom_abline(intercept = 38.91420934, slope = 0.01602637, color="#CC79A7", size=0.3)+
  ggtitle("c) 5-GMR")+
  labs(x= "Precip (mm)", y = "Forearm (mm)")+
  scale_colour_manual(values=c("#0072B2","#F0E442","#CC79A7", "#E69F00","#999999"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") 


P2.1<-grid.arrange(p1, p4,p2,ncol=3)
ggsave("fig1_chip.eps",P2.1,width = 9, 
       height = 2.8,limitsize = FALSE )

