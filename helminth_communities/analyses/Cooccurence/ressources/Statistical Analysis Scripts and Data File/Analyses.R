########################################################
#
#       Analyses to accompany article
#
#   Pathogen community composition and co-infection patterns
#   in a wild community of rodents.
#
#   Jessica L. Abbate , Maxime Galan , Maria Razzauti , Tarja Sironen ,
#   Liina  Voutilainen , Heikki Henttonen , Patrick Gasqui ,
#   Jean-François Cosson , Nathalie Charbonnel
#
#   bioRxiv:  https://doi.org/10.1101/2020.02.09.940494
#   First published Feb 10, 2020
#   Last updated July 6, 2023
#
########################################################



########################################################
#       Load Packages & Source Codes           
########################################################


rm(list=ls(all=TRUE))

########################################################
#       Load Packages & Source Codes           
########################################################

#R version 4.0.2 (2020-06-22)

library(vegan)
library(agricolae)
library(glmulti)
library(entropart)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(boot)
library(MASS)
library(jtools)
library(flextable)
library(cowplot)
library(sjPlot)
library(rstudioapi)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("SCN.txt")
source("FctTestScreenENV.txt")

########################################################
#       Load Data           
########################################################

## MASTER data for pathobiome analysis
master<-read.csv("PA_DATA.csv",header=TRUE)
master$YEAR<-factor(master$YEAR)

# remove rare host species and VirusTBE
focal<- droplevels(master[which(master$HostSPP!="M_Ap_spp" & master$HostSPP!="C_On_zibe" & master$HostSPP!="C_Mi_subt"  & master$HostSPP!="Myo_My_coyp"),])
focal<- focal[,which(colnames(focal)!="VirusTBE")]

# remove individuals with missing data
dat<-subset(focal,focal$BactMiSeq==1)
dat<-subset(dat,dat$VirusHV!="NA")
dat<-subset(dat,dat$VirusCPXV!="NA")
dat<-subset(dat,dat$VirusLCMV!="NA")
dat<-subset(dat,dat$SEX!="NA")
dat<-subset(dat,dat$AGE!="NA")
dat<-subset(dat,dat$SITE!="NA")
dat<-subset(dat,dat$YEAR!="NA")
dat<-subset(dat,dat$HABITAT!="NA")

# remove R. norvegicus and Myco10 for extrinsic factor tests
datNORAT<-droplevels(dat[which(dat$HostSPP!="M_Ra_norv"),c(which(names(dat)!="Myco10"))])
# datNORAT<-droplevels(dat) # use to run analyses on Ra_norv

# identify columns with pathogen presence data (except for VirusTBE)
p1<-which(names(master)=="VirusHV") # first pathogen column
b1<-which(names(master)=="Bartonella")  # first bacteria column
pN<-which(names(master)=="Spiroplasma") # last pathogen column
pN_norat<-which(names(datNORAT)=="Spiroplasma") # last pathogen column w/o Myco10

##############################################################
#   3.1   Describe the Data        
##############################################################

# 3.1.1

# for CPXV
t<-table(master$VirusCPXV)
tab<-t(table(focal$HostSPP,focal$VirusCPXV))
prev<-tab[2,]/colSums(tab)
chisq.test(tab)

# for Hantavirus
table(master$VirusHV)
t(table(master$HostSPP,master$VirusHV))
tab<-t(table(focal$HostSPP,focal$VirusHV))
prev<-tab[2,]/colSums(tab)
chisq.test(tab)

# for LCMV
t(table(master$HostSPP,master$VirusLCMV))

# for TBE
t(table(master$HostSPP,master$VirusTBE))


# 3.1.2

# reduce data to only those with 16S data
masterbac<-subset(master,master$BactMiSeq==1)
focalbac<-subset(focal,focal$BactMiSeq==1)

#for Bact exposure 
t<-table(master$BactPA)
tab<-t(table(focalbac$HostSPP,focalbac$BactPA))
prev<-tab[2,]/colSums(tab)
chisq.test(tab)


##############################################################
#   3.2   Pathogen Diversity       
##############################################################

# 3.2.1

# Shannon Index
dat$divDUMMY<-paste(dat$YEAR,"-",dat$SITE,"-",dat$HABITAT,"-",dat$HostSPP,"-",dat$AGE,"-",dat$SEX)

results<-as.data.frame(levels(as.factor(dat$divDUMMY)))
names(results)<-c("divDUMMY")
for(i in p1:pN){
  results<-cbind(results,NA)
  colnames(results)[length(colnames(results))]<-names(dat)[i]
  if("1" %in% colnames(table(dat$divDUMMY,dat[,i]))){
    results[,dim(results)[2]]<-table(dat$divDUMMY,dat[,i])[,"1"]}else{
      results[,dim(results)[2]]<-c(0,0,0,0,0,0)
    }
}

# calculate Shannon Index for each group
require(vegan)
hist(diversity(results[,2:20],index="shannon"))
results$divGRP<-diversity(results[,2:20],index="shannon")

# check for outliers
OutVals=boxplot(results$divGRP)$out
which(results$divGRP %in% OutVals)

# break the DUMMY variable back into variable groups
results$YEAR<-NA
results$SITE<-NA
results$HABITAT<-NA
results$HostSPP<-NA
results$AGE<-NA
results$SEX<-NA

for(i in 1:dim(results)[1]){
  if(i %in% grep('2010',(results$divDUMMY))==TRUE){results$YEAR[i]<-"2010"}
  if(i %in% grep('2011',(results$divDUMMY))==TRUE){results$YEAR[i]<-"2011"}
  if(i %in% grep('BRIQUENAY',(results$divDUMMY))==TRUE){results$SITE[i]<-"BRIQUENAY"}
  if(i %in% grep('BOULT',(results$divDUMMY))==TRUE){results$SITE[i]<-"BOULT"}
  if(i %in% grep('Forest',(results$divDUMMY))==TRUE){results$HABITAT[i]<-"Forest"}
  if(i %in% grep('Hedge',(results$divDUMMY))==TRUE){results$HABITAT[i]<-"Hedge"}
  if(i %in% grep('Meadow',(results$divDUMMY))==TRUE){results$HABITAT[i]<-"Meadow"}
  if(i %in% grep('Farm',(results$divDUMMY))==TRUE){results$HABITAT[i]<-"Farm"}
  if(i %in% grep('C_Mi_arva',(results$divDUMMY))==TRUE){results$HostSPP[i]<-"C_Mi_arva"}
  if(i %in% grep('C_Ar_sche',(results$divDUMMY))==TRUE){results$HostSPP[i]<-"C_Ar_sche"}
  if(i %in% grep('C_Mi_agre',(results$divDUMMY))==TRUE){results$HostSPP[i]<-"C_Mi_agre"}
  if(i %in% grep('C_My_glar',(results$divDUMMY))==TRUE){results$HostSPP[i]<-"C_My_glar"}
  if(i %in% grep('M_Ap_flav',(results$divDUMMY))==TRUE){results$HostSPP[i]<-"M_Ap_flav"}
  if(i %in% grep('M_Ap_sylv',(results$divDUMMY))==TRUE){results$HostSPP[i]<-"M_Ap_sylv"}
  if(i %in% grep('M_Ra_norv',(results$divDUMMY))==TRUE){results$HostSPP[i]<-"M_Ra_norv"}
  if(i %in% grep('ADULT',(results$divDUMMY))==TRUE){results$AGE[i]<-"ADULT"}
  if(i %in% grep('JUVENILE',(results$divDUMMY))==TRUE){results$AGE[i]<-"JUVENILE"}
  if(i %in% grep('m',(results$divDUMMY))==TRUE){results$SEX[i]<-"MALE"}
  if(i %in% grep('f',(results$divDUMMY))==TRUE){results$SEX[i]<-"FEMALE"}
}

results$YEAR<-factor(results$YEAR)

# remove R. norvegicus
resultsNORAT<-results[which(results$HostSPP!="M_Ra_norv"),]
# resultsNORAT<-results  # use to run stats on Ra_norv

# Test for role of extrinsic factors on pathogen diversity
# full model without R. norvegicus:
mod<-lm(divGRP~SITE+YEAR+AGE+SEX+HABITAT+HostSPP,data=resultsNORAT)
summary(mod,test="F")
drop1(mod,test="F")
g<-glmulti(mod, family="binomial", crit="aic", level=1)
print(g)
weightable(g)[1:5,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[1]]
plot_model(best_model, type = "pred",title = "Diversity")
mod11<-stepAIC(mod)
summary(mod11)
plot_coefs(mod11, legend.title = "divgrp") 

# tukey tests for host species groups
mod<-lm(divGRP~HABITAT+SITE+YEAR+AGE+SEX,data=resultsNORAT)
resultsNORAT$resids<-residuals(mod)
mod<-lm(resids~HostSPP,data=resultsNORAT)
test<-aov(mod)
(posthoc <- TukeyHSD(x=test, 'HostSPP', conf.level=0.95))
(out<-HSD.test(test, 'HostSPP'))

# tukey tests for habitats
mod<-lm(divGRP~HostSPP+SITE+YEAR+AGE+SEX,data=resultsNORAT)
resultsNORAT$resids<-residuals(mod)
mod<-lm(resids~HABITAT,data=resultsNORAT)
test<-aov(mod)
(posthoc <- TukeyHSD(x=test, 'HABITAT', conf.level=0.95))
(out<-HSD.test(test, 'HABITAT'))

# tukey tests for host age groups
mod<-lm(divGRP~HostSPP+HABITAT+SITE+YEAR+SEX,data=resultsNORAT)
resultsNORAT$resids<-residuals(mod)
hist(resultsNORAT$resids)
mod<-lm(resids~AGE,data=resultsNORAT)
test<-aov(mod)
(posthoc <- TukeyHSD(x=test, 'AGE', conf.level=0.95))
(out<-HSD.test(test, 'AGE'))

# examine roles of R. norvegicus and Farm habitats

mod<-lm(divGRP~SITE+YEAR+AGE+SEX,data=results)
results$resids<-residuals(mod)
test<-aov(lm(resids~HostSPP,data=results))
(posthoc <- TukeyHSD(x=test, 'HostSPP', conf.level=0.95))
(out<-HSD.test(test, 'HostSPP'))

test<-aov(lm(resids~HABITAT,data=results))
(posthoc <- TukeyHSD(x=test, 'HABITAT', conf.level=0.95))
(out<-HSD.test(test, 'HABITAT'))

##############################################################
# test for host spp diversity correlation with pathogen species diversity 
# per community (year x site x habitat) (no R. norvegicus)
##############################################################
datNORAT$coinfDUMMY<-paste(datNORAT$YEAR,"-",datNORAT$SITE,"-",datNORAT$HABITAT)
test<-as.data.frame(table(datNORAT$coinfDUMMY))
names(test)<-c("coinfDUMMY","N_hosts")

# host species number and diversity per community

hostspp_per_comm<-table(datNORAT$HostSPP,datNORAT$coinfDUMMY)
for(i in 1:dim(hostspp_per_comm)[2]){test$hostChao[i]<-Shannon(NorP=hostspp_per_comm[,i],Correction = "ChaoShen")}
hostspp_per_comm[hostspp_per_comm > 0] <- 1 
(number_hostspp_per_comm<-colSums(hostspp_per_comm))
test$N_host_spp<-number_hostspp_per_comm

# pathogen species number and diversity per community

paths_per_comm<-as.data.frame(as.matrix(table(datNORAT[,p1],datNORAT$coinfDUMMY)[2,]))
colnames(paths_per_comm)[1]<-paste(colnames(datNORAT)[p1])
for(i in (p1+1):pN_norat)
{paths_per_comm2<-as.data.frame(as.matrix(table(datNORAT[,i],datNORAT$coinfDUMMY)[2,]))
colnames(paths_per_comm2)<-paste(names(datNORAT)[i])
paths_per_comm<-cbind(paths_per_comm,paths_per_comm2)}
paths_per_comm_Abundance<-paths_per_comm
test$PathoDivPerComm<-diversity(paths_per_comm_Abundance,index="shannon")

# test for host spp diversity correlation with pathogen species diversity 

cor.test(test$hostChao,test$PathoDivPerComm)

##############################################################
# co-infection/co-exposure frequency analyses per host species
##############################################################

# calculate the number of bacteria and pathogens per individual
dat$number_bacteria<-rowSums(dat[,b1:pN])
dat$number_pathogens<-rowSums(dat[,p1:pN])

#-------------------------------------------------------------
# test for correlation of bacterial co-infection frequency w/bacterial OTU diversity
coinf_table<-table(dat$number_bacteria, dat$HostSPP)
coinfs<-as.data.frame(colSums(coinf_table[c(3:6),])/colSums(coinf_table))
names(coinfs)<-c("coinfection_proportion")
coinfs$N_hosts<-colSums(coinf_table)

bacs_per_spp<-as.data.frame(as.matrix(table(dat[,b1],dat$HostSPP)[2,]))
colnames(bacs_per_spp)[1]<-paste(colnames(dat)[b1])
for(i in (b1+1):pN)
{bacs_per_spp2<-as.data.frame(as.matrix(table(dat[,i],dat$HostSPP)[2,]))
colnames(bacs_per_spp2)<-paste(names(dat)[i])
bacs_per_spp<-cbind(bacs_per_spp,bacs_per_spp2)}
bacs_per_spp_Abundance<-bacs_per_spp
coinfs$BactDivPerSPP<-diversity(bacs_per_spp_Abundance,index="shannon")

mod<-glm(coinfection_proportion~BactDivPerSPP,data=coinfs,family="binomial", weights=N_hosts)
summary(mod,test="Chisq")
preds<-predict(mod,se.fit=TRUE)
(rsq<-1-(mod$deviance/mod$null.deviance))  #Pseudo R2 for GLMs
hist(residuals(mod))

#outliers (Cook's distance)
cooksd<-cooks.distance(mod)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+0.75, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""))  # add labels
text(x=1:length(cooksd)+0.75, y=cooksd, labels=ifelse(cooksd>0.25*mean(cooksd, na.rm=T),names(cooksd),""))  # add labels

#-------------------------------------------------------------
# test for correlation of pathogen co-exposure frequency w/pathogen diversity
coexp_table<-table(dat$number_pathogens, dat$HostSPP)
coexps<-as.data.frame(colSums(coexp_table[c(3:6),])/colSums(coexp_table))
names(coexps)<-c("coexposure_proportion")
coexps$N_hosts<-colSums(coexp_table)

paths_per_spp<-as.data.frame(as.matrix(table(dat[,p1],dat$HostSPP)[2,]))
colnames(paths_per_spp)[1]<-paste(colnames(dat)[p1])
for(i in (p1+1):pN)
{paths_per_spp2<-as.data.frame(as.matrix(table(dat[,i],dat$HostSPP)[2,]))
colnames(paths_per_spp2)<-paste(names(dat)[i])
paths_per_spp<-cbind(paths_per_spp,paths_per_spp2)}
paths_per_spp_Abundance<-paths_per_spp
coexps$PathDivPerSPP<-diversity(paths_per_spp_Abundance,index="shannon")

mod<-glm(coexposure_proportion~PathDivPerSPP,data=coexps,family="binomial", weights=N_hosts)
summary(mod,test="Chisq")
preds<-predict(mod,se.fit=TRUE)
(rsq<-1-(mod$deviance/mod$null.deviance))  #Pseudo R2 for GLMs
hist(residuals(mod))

#outliers (Cook's distance)
cooksd<-cooks.distance(mod)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+0.75, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""))  # add labels
text(x=1:length(cooksd)+0.75, y=cooksd, labels=ifelse(cooksd>0.25*mean(cooksd, na.rm=T),names(cooksd),""))  # add labels

##############################################################
# co-infection/co-exposure frequency analyses (without R. norvegicus) per community
##############################################################
#--------------------------------------------------------------
# calculate the number of bacterial co-infections per individual (no R. norv)
datNORAT$number_bacteria<-rowSums(datNORAT[,b1:pN_norat])

# calculate the number of bacterial co-infections per community (year x site x habitat) (no R. norv)
coinf_table<-table(datNORAT$number_bacteria,datNORAT$coinfDUMMY)
coinfs<-as.data.frame(colSums(coinf_table[c(3:6),])/colSums(coinf_table))
names(coinfs)<-c("coinfection_proportion")
test<-cbind(test,coinfs)

# test against host diversity (no R. norv)
cor.test(asin(sqrt(test$coinfection_proportion)),test$hostChao) 

#--------------------------------------------------------------
# calculate the number of pathogen co-exposures per individual (no R. norv)
datNORAT$number_pathogens<-rowSums(datNORAT[,p1:pN_norat])

# calculate the number of pathogen co-exposures per community (year x site x habitat) (no R. norv)
coinf_table<-table(datNORAT$number_pathogens,datNORAT$coinfDUMMY)
coexp<-as.data.frame(colSums(coinf_table[c(3:6),])/colSums(coinf_table))
names(coexp)<-c("coexposure_proportion")
test<-cbind(test,coexp)

# test against host diversity (no R. norv)
cor.test(asin(sqrt(test$coexposure_proportion)),test$hostChao)


##############################################################
# Pathogen Community composition : PERMANOVA &  MCA
##############################################################

# 3.2.2

# define restricted dataset

pathocommNORAT<-c( "Bartonella","Myco1","Myco2","Myco3","Myco7","Myco9",
                   "Myco4","Myco6","Rickettsia","Neoehrlichia","Orientia",
                   "Spiroplasma","VirusCPXV","VirusHV")
extfacs<- c("SITE", "HABITAT","HostSPP","YEAR","SEX","AGE")

pcNORAT<-as.data.frame(unclass(datNORAT[,c(extfacs,pathocommNORAT)]),stringsAsFactors=TRUE)
pcNORAT$YEAR<-factor(pcNORAT$YEAR)

# Beta Diversity analysis ------------------------------------

# remove uninfected individuals
nozero <- pcNORAT[which(rowSums(pcNORAT[,c(pathocommNORAT)])>0),]

# set datasets (integer matrix and factor matrix)
spe<-nozero[,c(pathocommNORAT)]
env<-nozero[,c(extfacs)]

# Hellinger pre-transformation
spe.hell <- decostand (spe, 'hell')  

# PERMANOVA on Bray-Curtis Dissimilarity Matrix
set.seed<-1
adonis2(spe.hell~.,data=env,dist="bray",by="margin")

# MCA analysis ----------------------------------------------

pcNORAT[,c(pathocommNORAT)]<-lapply(pcNORAT[,c(pathocommNORAT)],factor)
compsev.active <- pcNORAT[,c(pathocommNORAT)]

# Eigenvalue scree plot
res.mca.compsev <- MCA(compsev.active, graph = FALSE,ncp=length(pathocommNORAT))
get_eigenvalue(res.mca.compsev)
p<-fviz_screeplot(res.mca.compsev,addlabels=TRUE,ncp=length(pathocommNORAT))
p + labs(title = "Variances - MCA",
         x = "MCA Dimensions", y = "% of explained variances")+geom_hline(yintercept=(100/length(pathocommNORAT)),linetype="dashed",col="red")
plot(res.mca.compsev,axes=c(1,2),choix="var",cex=1, ylim=c(0,0.35), xlim=c(0,0.45))

# Contribution of each variable to each dimension
var <- get_mca_var(res.mca.compsev)
corrplot(var$contrib,is.corr=FALSE,cl.cex=0.5)

habitats <- pcNORAT$HABITAT
ages<-pcNORAT$AGE
species<-pcNORAT$HostSPP

fviz_mca_biplot(res.mca.compsev, axes=c(1,2), invisible=c("var"),
                habillage = habitats, addEllipses = TRUE,
                label = "var", shape.var = 15) +theme_minimal()
fviz_mca_biplot(res.mca.compsev, axes=c(1,2), invisible=c("var"),
                habillage = ages, addEllipses = TRUE,
                label = "var", shape.var = 15) +theme_minimal()
fviz_mca_biplot(res.mca.compsev, axes=c(1,2), invisible=c("var"),
                habillage = species, addEllipses = TRUE,
                label = "var", shape.var = 15) +theme_minimal()

pcNORAT$Dim1<-res.mca.compsev$ind$coord[,1]
pcNORAT$Dim2<-res.mca.compsev$ind$coord[,2]
pcNORAT$Dim3<-res.mca.compsev$ind$coord[,3]
pcNORAT$Dim4<-res.mca.compsev$ind$coord[,4]
pcNORAT$Dim5<-res.mca.compsev$ind$coord[,5]
pcNORAT$Dim6<-res.mca.compsev$ind$coord[,6]
pcNORAT$Dim7<-res.mca.compsev$ind$coord[,7]

# Extrinsic Factors structuring MCA Dims ------------------

# MCA Dim 1
mod<-lm(pcNORAT$Dim1~pcNORAT$HostSPP+pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$YEAR+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
summary(mod)
drop1(mod,.~.,test="F")

# tukey tests
mod<-lm(pcNORAT$Dim1~pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$YEAR+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
pcNORAT$resids<-residuals(mod)
test<-aov(lm(pcNORAT$resids~pcNORAT$HostSPP,data=pcNORAT))
(out<-HSD.test(test, 'pcNORAT$HostSPP'))

mod<-lm(pcNORAT$Dim1~pcNORAT$HostSPP+pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$YEAR+pcNORAT$SEX,data=pcNORAT)
pcNORAT$resids<-residuals(mod)
test<-aov(lm(pcNORAT$resids~pcNORAT$AGE,data=pcNORAT))
(out<-HSD.test(test, 'pcNORAT$AGE'))

# MCA Dim 2
mod<-lm(pcNORAT$Dim2~pcNORAT$HostSPP+pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$YEAR+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
summary(mod)
drop1(mod,.~.,test="F")

mod_x<-stepAIC(mod)
summary(mod_x)
plot_coefs(mod_x)


# tukey tests
mod<-lm(pcNORAT$Dim2~pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$YEAR+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
pcNORAT$resids<-residuals(mod)
test<-aov(lm(pcNORAT$resids~pcNORAT$HostSPP,data=pcNORAT))
(out<-HSD.test(test, 'pcNORAT$HostSPP'))

mod<-lm(pcNORAT$Dim2~pcNORAT$HostSPP+pcNORAT$SITE+pcNORAT$YEAR+pcNORAT$SEX+pcNORAT$AGE,data=pcNORAT)
pcNORAT$resids<-residuals(mod)
test<-aov(lm(pcNORAT$resids~pcNORAT$HABITAT,data=pcNORAT))
(out<-HSD.test(test, 'pcNORAT$HABITAT'))

mod<-lm(pcNORAT$Dim2~pcNORAT$HostSPP+pcNORAT$HABITAT+pcNORAT$YEAR+pcNORAT$SEX+pcNORAT$AGE,data=pcNORAT)
pcNORAT$resids<-residuals(mod)
test<-aov(lm(pcNORAT$resids~pcNORAT$SITE,data=pcNORAT))
(out<-HSD.test(test, 'pcNORAT$SITE'))

# Dim 3
mod<-lm(pcNORAT$Dim3~pcNORAT$HostSPP+pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$YEAR+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
summary(mod)
drop1(mod,.~.,test="F")

mod_x<-stepAIC(mod)
summary(mod_x)
plot_coefs(mod_x)

#tukey tests
mod<-lm(pcNORAT$Dim3~pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$YEAR+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
pcNORAT$resids<-residuals(mod)
test<-aov(lm(pcNORAT$resids~pcNORAT$HostSPP,data=pcNORAT))
(out<-HSD.test(test, 'pcNORAT$HostSPP'))

mod<-lm(pcNORAT$Dim3~pcNORAT$HostSPP+pcNORAT$HABITAT+pcNORAT$YEAR+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
pcNORAT$resids<-residuals(mod)
test<-aov(lm(pcNORAT$resids~pcNORAT$SITE,data=pcNORAT))
(out<-HSD.test(test, 'pcNORAT$SITE'))

mod<-lm(pcNORAT$Dim3~pcNORAT$HostSPP+pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
pcNORAT$resids<-residuals(mod)
test<-aov(lm(pcNORAT$resids~pcNORAT$YEAR,data=pcNORAT))
(out<-HSD.test(test, 'pcNORAT$YEAR'))

# Dim 4
mod<-lm(pcNORAT$Dim4~pcNORAT$HostSPP+pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$YEAR+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
summary(mod)
drop1(mod,.~.,test="F")

mod_x<-stepAIC(mod, direction="both")
summary(mod_x)
plot_coefs(mod_x)
hist(residuals(mod_x))

#tukey tests
mod<-lm(pcNORAT$Dim4~pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$YEAR+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
pcNORAT$resids<-residuals(mod)
test<-aov(lm(pcNORAT$resids~pcNORAT$HostSPP,data=pcNORAT))
(out<-HSD.test(test, 'pcNORAT$HostSPP'))

mod<-lm(pcNORAT$Dim4~pcNORAT$HostSPP+pcNORAT$SITE+pcNORAT$YEAR+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
pcNORAT$resids<-residuals(mod)
test<-aov(lm(pcNORAT$resids~pcNORAT$HABITAT,data=pcNORAT))
(out<-HSD.test(test, 'pcNORAT$HABITAT'))

mod<-lm(pcNORAT$Dim4~pcNORAT$HostSPP+pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
pcNORAT$resids<-residuals(mod)
test<-aov(lm(pcNORAT$resids~pcNORAT$YEAR,data=pcNORAT))
(out<-HSD.test(test, 'pcNORAT$YEAR'))

mod<-lm(pcNORAT$Dim4~pcNORAT$HostSPP+pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$YEAR+pcNORAT$SEX,data=pcNORAT)
pcNORAT$resids<-residuals(mod)
test<-aov(lm(pcNORAT$resids~pcNORAT$AGE,data=pcNORAT))
(out<-HSD.test(test, 'pcNORAT$AGE'))

# Dim 5
mod<-lm(pcNORAT$Dim5~pcNORAT$HostSPP+pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$YEAR+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
summary(mod)
drop1(mod,.~.,test="F")

mod_x<-stepAIC(mod)
summary(mod_x)
# plot_coefs(mod_x) no factors remaining

# Dim 6
mod<-lm(pcNORAT$Dim6~pcNORAT$HostSPP+pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$YEAR+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
summary(mod)
drop1(mod,.~.,test="F")

mod_x<-stepAIC(mod)
summary(mod_x)
plot_coefs(mod_x)
hist(residuals(mod_x))

#tukey tests
mod<-lm(pcNORAT$Dim6~pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$YEAR+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
pcNORAT$resids<-residuals(mod)
test<-aov(lm(pcNORAT$resids~pcNORAT$HostSPP,data=pcNORAT))
(out<-HSD.test(test, 'pcNORAT$HostSPP'))

# Dim 7
mod<-lm(pcNORAT$Dim7~pcNORAT$HostSPP+pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$YEAR+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
summary(mod)
drop1(mod,.~.,test="F")

mod_x<-stepAIC(mod)
summary(mod_x)
plot_coefs(mod_x)

#tukey tests
mod<-lm(pcNORAT$Dim7~pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$YEAR+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
pcNORAT$resids<-residuals(mod)
test<-aov(lm(pcNORAT$resids~pcNORAT$HostSPP,data=pcNORAT))
(out<-HSD.test(test, 'pcNORAT$HostSPP'))

mod<-lm(pcNORAT$Dim7~pcNORAT$HostSPP+pcNORAT$HABITAT+pcNORAT$YEAR+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
pcNORAT$resids<-residuals(mod)
test<-aov(lm(pcNORAT$resids~pcNORAT$SITE,data=pcNORAT))
(out<-HSD.test(test, 'pcNORAT$SITE'))

mod<-lm(pcNORAT$Dim7~pcNORAT$HostSPP+pcNORAT$HABITAT+pcNORAT$SITE+pcNORAT$AGE+pcNORAT$SEX,data=pcNORAT)
pcNORAT$resids<-residuals(mod)
test<-aov(lm(pcNORAT$resids~pcNORAT$YEAR,data=pcNORAT))
(out<-HSD.test(test, 'pcNORAT$YEAR'))


##############################################################
#   3.3   Pathogen-Pathogen Associations  
##############################################################

# Myco1 vs. Myco3 vs. VirusHV --------------------------------

# define restricted dataset

pathos<-c( "Myco1","Myco3","VirusHV")
extfacs<- c("SITE", "HABITAT","HostSPP","YEAR","SEX","AGE")

# SCN on just Mi arvalis
test1<-datNORAT[which(datNORAT$AGE=="ADULT" & datNORAT$HostSPP=="C_Mi_arva"),c(extfacs,pathos)]
SCN(test1,pathos)

# SCN on just Mi arvalis
test2<-datNORAT[which(datNORAT$AGE=="ADULT" & datNORAT$HostSPP=="C_My_glar"),c(extfacs,pathos)]
SCN(test2,pathos)

# SCN on both Mi arvalis & My glareolus
test<-datNORAT[which(datNORAT$AGE=="ADULT" & (datNORAT$HostSPP=="C_Mi_arva" | datNORAT$HostSPP=="C_My_glar")),c(extfacs,pathos)]
SCN(test,pathos)

# global model Myco1
mod<-glm(Myco1~YEAR+SEX+HostSPP+SITE+HABITAT+VirusHV+Myco3,data=test,family="binomial")
summary(mod)
drop1(mod,.~.,test="Chisq")
confint.default(mod)
exp(cbind(OR = coef(mod), confint(mod)))

#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:4,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[1]]
plot_model(best_model, type = "pred",title = "Myco1")

#best model with pathogens of interest: "Myco1 ~ 1 + HostSPP + SITE + VirusHV + Myco3" 
mod<-glm(Myco1~HostSPP+SITE+VirusHV+Myco3,data=test,family="binomial")
drop1(mod,.~.,test="Chisq")
exp(cbind(OR = coef(mod), confint(mod)))

# global model Myco3
mod<-glm(Myco3~YEAR+HostSPP+HABITAT+SITE+SEX+VirusHV+Myco1,data=test,family="binomial")
summary(mod)
drop1(mod,.~.,test="Chisq")
confint.default(mod)
exp(cbind(OR = coef(mod), confint(mod)))
table(test$VirusHV,test$Myco3)

#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:5,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[1]]
plot_model(best_model, type = "pred",title = "Myco3")

#best model with pathogens of interest: "Myco3 ~ 1 + YEAR + VirusHV + Myco1"
mod<-glm(Myco3~YEAR + VirusHV + Myco1,data=test,family="binomial")
drop1(mod,.~.,test="Chisq")
exp(cbind(OR = coef(mod), confint(mod)))

# Hanta associated with Myco3, perfect positive association
test[which(test$VirusHV==1),c("AGE","HostSPP","Myco1","Myco3","YEAR","SITE","SEX","HABITAT")]

# global model anti-Hantavirus antibodies (Myco1 only)
mod<-glm(VirusHV~YEAR+SEX+HostSPP+SITE+HABITAT+Myco1,data=test,family="binomial")
drop1(mod,.~., test="Chisq")
coef(mod)
confint.default(mod)
# exp(cbind(OR = coef(mod), confint(mod)))

#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:2,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[1]]
plot_model(best_model, type = "pred",title = "Hanta Abs")

#best model: "VirusHV ~ 1 + HABITAT + Myco1"
mod<-glm(VirusHV~ HABITAT + Myco1,data=test,family="binomial")
drop1(mod,.~.,test="Chisq")
exp(cbind(OR = coef(mod), confint(mod)))


# Myco2 vs. Myco4 -------------------------------------------

# define restricted dataset

pathos<-c( "Myco2","Myco4")
extfacs<- c("SITE", "HABITAT","HostSPP","YEAR","SEX","AGE")

#only in Ap. sylvaticus
test<-datNORAT[which(datNORAT$HostSPP=="M_Ap_sylv"),c(extfacs,pathos)]

# SCN
SCN(test,pathos)


# global model Myco2
mod1<-glm(Myco2~YEAR+AGE+HABITAT+SITE+SEX+Myco4,data=test,family="binomial")
drop1(mod1,.~.,test="Chisq")

mod2<-glm(Myco4~YEAR+AGE+HABITAT+SITE+SEX+Myco2,data=test,family="binomial")
drop1(mod2,.~.,test="Chisq")

#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod1, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:4,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[1]]
plot_model(best_model, type = "pred",title = "Myco2")

#best model: "Myco2 ~ 1 + AGE + Myco4"
mod<-glm(Myco2 ~ AGE + Myco4,data=test,family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) # positive
confint.default(mod)
exp(cbind(OR = coef(mod), confint(mod)))

g<-glmulti(mod2, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:9,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[1]]
plot_model(best_model, type = "pred",title = "Myco4")

#best model: "Myco4 ~ 1 + Myco2"
mod<-glm(Myco4 ~ Myco2,data=test,family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) # positive
confint.default(mod)
exp(cbind(OR = coef(mod), confint(mod)))

# Bartonella vs. anti-VirusCPXV Abs --------------------------

# define restricted dataset

pathos<-c( "Bartonella","VirusCPXV")
extfacs<- c("SITE", "HABITAT","HostSPP","YEAR","SEX","AGE")

test<-dat[,c(extfacs,pathos)]

# test for interaction between host species and explanatory pathogen
mod1<-glm(Bartonella~YEAR+AGE+HABITAT+SITE+SEX+HostSPP*VirusCPXV,data=test,family="binomial")
drop1(mod1,.~.,test="Chisq")

mod2<-glm(VirusCPXV~YEAR+AGE+HABITAT+SITE+SEX+HostSPP*Bartonella,data=test,family="binomial")
drop1(mod2,.~.,test="Chisq")

# test for association in each host species separately

# --> only in C_Ar_sche  (only in Meadows)
SCN(test[which(test$HostSPP=="C_Ar_sche"),],pathos)

mod<-glm(Bartonella~YEAR+AGE+SITE+SEX+VirusCPXV,data=test[which(test$HostSPP=="C_Ar_sche"),],family="binomial")
drop1(mod,.~.,test="Chisq")
summary(mod)
coef(mod) # positive association
exp(cbind(OR = coef(mod), confint(mod)))

#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:2,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[1]]
plot_model(best_model, type = "pred",title = "Bartonella")
#best model: "Bartonella ~ 1 + YEAR + VirusCPXV"
mod<-glm(Bartonella ~ 1 + YEAR + VirusCPXV,data=test[which(test$HostSPP=="C_Ar_sche"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) # positive association
exp(cbind(OR = coef(mod), confint(mod)))

mod<-glm(VirusCPXV~YEAR+AGE+SITE+SEX+Bartonella,data=test[which(test$HostSPP=="C_Ar_sche"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) # positive association
exp(cbind(OR = coef(mod), confint(mod)))

#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:3,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[1]]
plot_model(best_model, type = "pred",title = "CPXV")
#best model: "VirusCPXV ~ 1 + YEAR + Bartonella"
mod<-glm(VirusCPXV ~ 1 + YEAR + Bartonella,data=test[which(test$HostSPP=="C_Ar_sche"),],family="binomial")
drop1(mod,.~.,test="Chisq")
exp(cbind(OR = coef(mod), confint(mod)))

# --> only in C_Mi_agre
SCN(test[which(test$HostSPP=="C_Mi_agre"),],pathos)

mod<-glm(Bartonella~AGE+VirusCPXV,data=test[which(test$HostSPP=="C_Mi_agre"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) #none
summary(mod)
#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:2,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[2]]
plot_model(best_model, type = "pred",title = "Bartonella")
#best model: "Bartonella ~ 1 + AGE + VirusCPXV"
mod<-glm(Bartonella~AGE+VirusCPXV,data=test[which(test$HostSPP=="C_Mi_agre"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) #none
exp(cbind(OR = coef(mod), confint(mod)))

mod<-glm(VirusCPXV~AGE+SITE+HABITAT+Bartonella,data=test[which(test$HostSPP=="C_Mi_agre"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) #none
#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:3,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[2]]
plot_model(best_model, type = "pred",title = "CPXV")
#best model: CPXV ~ AGE + Bart
mod<-glm(VirusCPXV~AGE+Bartonella,data=test[which(test$HostSPP=="C_Mi_agre"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) #none
exp(cbind(OR = coef(mod), confint(mod)))


# --> only in C_Mi_arva
SCN(test[which(test$HostSPP=="C_Mi_arva"),],pathos)

mod<-glm(Bartonella~YEAR+AGE+SITE+SEX+VirusCPXV,data=test[which(test$HostSPP=="C_Mi_arva"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) #none #w/o habitat (mostly meadow)
summary(mod)
exp(cbind(OR = coef(mod), confint(mod)))
#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:6,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[5]]
plot_model(best_model, type = "pred",title = "Bartonella")
#best model with pathogens of interest: "Bartonella ~ 1 + YEAR + AGE + SITE + VirusCPXV"
mod<-glm(Bartonella ~ YEAR + AGE + SITE + VirusCPXV,data=test[which(test$HostSPP=="C_Mi_arva"),],family="binomial")
drop1(mod,.~.,test="Chisq")

mod<-glm(VirusCPXV~YEAR+AGE+SITE+SEX+Bartonella,data=test[which(test$HostSPP=="C_Mi_arva"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) #none
summary(mod)
exp(cbind(OR = coef(mod), confint(mod)))
#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:3,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[3]]
plot_model(best_model, type = "pred",title = "CPXV")
#best model with pathogens of interest: "VirusCPXV ~ 1 + SITE+Bartonella"
mod<-glm(VirusCPXV ~ SITE + Bartonella,data=test[which(test$HostSPP=="C_Mi_arva"),],family="binomial")
drop1(mod,.~.,test="Chisq")

# --> only in C_My_glar
SCN(test[which(test$HostSPP=="C_My_glar"),],pathos)

mod<-glm(Bartonella~YEAR+AGE+SITE+SEX+HABITAT+VirusCPXV,data=test[which(test$HostSPP=="C_My_glar"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) #none
summary(mod)
exp(cbind(OR = coef(mod), confint(mod)))
#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:3,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[1]]
plot_model(best_model, type = "pred",title = "Bartonella")
#best model with pathogens of interest: "Bartonella ~ 1 + AGE + SITE + HABITAT + VirusCPXV"
mod<-glm(Bartonella ~ AGE + SITE + HABITAT + VirusCPXV,data=test[which(test$HostSPP=="C_My_glar"),],family="binomial")
drop1(mod,.~.,test="Chisq")

mod<-glm(VirusCPXV~YEAR+AGE+SITE+SEX+HABITAT+Bartonella,data=test[which(test$HostSPP=="C_My_glar"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) #none
summary(mod)
exp(cbind(OR = coef(mod), confint(mod)))
#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:4,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[1]]
plot_model(best_model, type = "pred",title = "CPXV")
#best model with pathogens of interest: "VirusCPXV ~ 1 + YEAR + Bartonella"
mod<-glm(VirusCPXV ~ YEAR+Bartonella,data=test[which(test$HostSPP=="C_My_glar"),],family="binomial")
drop1(mod,.~.,test="Chisq")

# --> only in M_Ap_flav
SCN(test[which(test$HostSPP=="M_Ap_flav"),],pathos)

mod<-glm(Bartonella~YEAR+AGE+SITE+SEX+HABITAT+VirusCPXV,data=test[which(test$HostSPP=="M_Ap_flav"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) #trend positive
summary(mod)
exp(cbind(OR = coef(mod), confint(mod)))
#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:9,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[8]]
plot_model(best_model, type = "pred",title = "Bartonella")
#best model with pathogens of interest: "Bartonella ~ 1 + SEX + YEAR+VirusCPXV"
mod<-glm(Bartonella ~ SEX+YEAR+VirusCPXV,data=test[which(test$HostSPP=="M_Ap_flav"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) #negative
exp(cbind(OR = coef(mod), confint(mod)))

mod<-glm(VirusCPXV~YEAR+AGE+SITE+SEX+HABITAT+Bartonella,data=test[which(test$HostSPP=="M_Ap_flav"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) #none
summary(mod)
exp(cbind(OR = coef(mod), confint(mod)))
#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:10,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[3]]
plot_model(best_model, type = "pred",title = "CPXV")
#best model with pathogens of interest: "VirusCPXV ~ 1 + AGE + SITE + SEX + Bartonella"
mod<-glm(VirusCPXV ~ AGE + SITE + SEX + Bartonella,data=test[which(test$HostSPP=="M_Ap_flav"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) #negative 
exp(cbind(OR = coef(mod), confint(mod)))

# --> only in M_Ap_sylv
SCN(test[which(test$HostSPP=="M_Ap_sylv"),],pathos)

mod<-glm(Bartonella~YEAR+AGE+SITE+SEX+HABITAT+VirusCPXV,data=test[which(test$HostSPP=="M_Ap_sylv"),],family="binomial")
summary(mod)
drop1(mod,.~.,test="Chisq")
coef(mod) 
exp(cbind(OR = coef(mod), confint(mod)))
#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:12,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[11]]
plot_model(best_model, type = "pred",title = "Bartonella")
#best model with pathogens of interest: "Bartonella ~ 1 + YEAR +AGE+SEX+ VirusCPXV"
mod<-glm(Bartonella ~ YEAR+ AGE+SEX+ VirusCPXV,data=test[which(test$HostSPP=="M_Ap_sylv"),],family="binomial")
drop1(mod,.~.,test="Chisq")

mod<-glm(VirusCPXV~YEAR+AGE+SITE+SEX+HABITAT+Bartonella,data=test[which(test$HostSPP=="M_Ap_sylv"),],family="binomial")
summary(mod)
drop1(mod,.~.,test="Chisq")
coef(mod) 
exp(cbind(OR = coef(mod), confint(mod)))
#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:7,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[2]]
plot_model(best_model, type = "pred",title = "CPXV")
#best model with pathogens of interest: "VirusCPXV ~ 1 + AGE + SEX + Bartonella"
mod<-glm(VirusCPXV ~ AGE + SEX + Bartonella,data=test[which(test$HostSPP=="M_Ap_sylv"),],family="binomial")
drop1(mod,.~.,test="Chisq")


# --> only in M_Ra_norv
SCN(test[which(test$HostSPP=="M_Ra_norv"),],pathos)

mod<-glm(Bartonella~YEAR+SITE+AGE+SEX+VirusCPXV,data=test[which(test$HostSPP=="M_Ra_norv"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod)
table(test$Bartonella[which(test$HostSPP=="M_Ra_norv")],test$VirusCPXV[which(test$HostSPP=="M_Ra_norv")])
#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:6,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[4]]
plot_model(best_model, type = "pred",title = "Bartonella")
#best model with pathogens of interest: "Bartonella ~ 1 + SITE + VirusCPXV"
mod<-glm(Bartonella ~ SITE + VirusCPXV,data=test[which(test$HostSPP=="M_Ra_norv"),],family="binomial")
drop1(mod,.~.,test="Chisq")

mod<-glm(VirusCPXV~YEAR+SITE+AGE+SEX+Bartonella,data=test[which(test$HostSPP=="M_Ra_norv"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod)
#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:5,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[5]]
plot_model(best_model, type = "pred",title = "CPXV")
#best model with pathogens of interest: "VirusCPXV ~ 1 + Bartonella"
mod<-glm(VirusCPXV ~ Bartonella,data=test[which(test$HostSPP=="M_Ra_norv"),],family="binomial")
drop1(mod,.~.,test="Chisq")

mod<-glm(VirusCPXV~AGE+Bartonella,data=test[which(test$HostSPP=="M_Ra_norv"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod)
summary(mod)
exp(cbind(OR = coef(mod), confint(mod)))


# Myco. haemomuris vs. Myco. coccoides  --------------------------

# define restricted dataset

dat$MH<-0
dat$MH[which(dat$Myco1==1|dat$Myco2==1|dat$Myco3==1|dat$Myco5==1|dat$Myco7==1|dat$Myco8==1|dat$Myco9==1)]<-1
dat$MC<-0
dat$MC[which(dat$Myco4==1|dat$Myco6==1)]<-1

pathos<-c( "MH","MC")
extfacs<- c("SITE", "HABITAT","HostSPP","YEAR","SEX","AGE")

# SCN on each host species
test1<-dat[which(dat$HostSPP=="C_My_glar"),c(extfacs,pathos)]
SCN(test1,pathos)
test2<-dat[which(dat$HostSPP=="M_Ap_flav"),c(extfacs,pathos)]
SCN(test2,pathos)
test3<-dat[which(dat$HostSPP=="M_Ap_sylv"),c(extfacs,pathos)]
SCN(test3,pathos)
test4<-dat[which(dat$HostSPP=="M_Ra_norv"),c(extfacs,pathos)]
SCN(test4,pathos)

# global logistic regression
test<-dat[which(dat$HostSPP=="C_My_glar" | dat$HostSPP=="M_Ap_flav" | dat$HostSPP=="M_Ap_sylv" | dat$HostSPP=="M_Ra_norv"),c(extfacs,pathos)]

mod<-glm(MH~YEAR+SITE+SEX+AGE+HABITAT+HostSPP*MC,data=test, family="binomial")
drop1(mod,.~.,test="Chisq") # no interaction

mod<-glm(MH~YEAR+SITE+SEX+AGE+HABITAT+HostSPP+MC,data=test, family="binomial")
summary(mod)
drop1(mod,.~.,test="Chisq")
coef(mod)
exp(cbind(OR = coef(mod), confint(mod)))

#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:4,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[2]]
plot_model(best_model, type = "pred",title = "MH")
#best model: "MH ~ 1 + YEAR + SITE + HABITAT + HostSPP + MC"
mod<-glm(MH ~ YEAR + SITE + HABITAT + HostSPP + MC,data=test,family="binomial")
drop1(mod,.~.,test="Chisq")
exp(cbind(OR = coef(mod), confint(mod)))

mod<-glm(MH ~ YEAR + SITE + HABITAT + HostSPP*MC,data=test,family="binomial")
drop1(mod,.~.,test="Chisq")
exp(cbind(OR = coef(mod), confint(mod)))

mod<-glm(MC~YEAR+SITE+SEX+AGE+HABITAT+HostSPP*MH,data=test, family="binomial")
drop1(mod,.~.,test="Chisq") # no interaction

mod<-glm(MC~YEAR+SITE+SEX+AGE+HABITAT+HostSPP+MH,data=test, family="binomial")
summary(mod)
drop1(mod,.~.,test="Chisq")
coef(mod)
exp(cbind(OR = coef(mod), confint(mod)))

#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:10,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[1]]
plot_model(best_model, type = "pred",title = "MC")
#best model: "MC ~ 1 + SEX+ HostSPP + MH"
mod<-glm(MC ~ SEX+HostSPP + MH,data=test,family="binomial")
drop1(mod,.~.,test="Chisq")
exp(cbind(OR = coef(mod), confint(mod)))

mod<-glm(MC ~ SEX+HostSPP*MH,data=test,family="binomial")
drop1(mod,.~.,test="Chisq")


# Bartonella vs. hemotropic Mycoplasma spp.  ------------------

# define restricted dataset

dat$HMyco<-0
dat$HMyco[which(dat$Myco1==1|dat$Myco2==1|dat$Myco3==1|dat$Myco5==1|dat$Myco7==1|dat$Myco8==1|dat$Myco9==1|dat$Myco4==1|dat$Myco6==1)]<-1

pathos<-c( "Bartonella","HMyco")
extfacs<- c("SITE", "HABITAT","HostSPP","YEAR","SEX","AGE")

test<-dat[,c(extfacs,pathos)]

# SCN on each host species
test1<-dat[which(dat$HostSPP=="C_Ar_sche"),c(extfacs,pathos)]
SCN(test1,pathos)
test2<-dat[which(dat$HostSPP=="C_Mi_agre"),c(extfacs,pathos)]
SCN(test2,pathos)
test3<-dat[which(dat$HostSPP=="C_Mi_arva"),c(extfacs,pathos)]
SCN(test3,pathos)
test4<-dat[which(dat$HostSPP=="C_My_glar"),c(extfacs,pathos)]
SCN(test4,pathos)
test5<-dat[which(dat$HostSPP=="M_Ap_flav"),c(extfacs,pathos)]
SCN(test5,pathos)
test6<-dat[which(dat$HostSPP=="M_Ap_sylv"),c(extfacs,pathos)]
SCN(test6,pathos)
test7<-dat[which(dat$HostSPP=="M_Ra_norv"),c(extfacs,pathos)]
SCN(test7,pathos)

# global logistic regressions
mod1<-glm(Bartonella~HABITAT+SITE+YEAR+AGE+SEX+HostSPP*HMyco,data=test,family="binomial")
drop1(mod1,.~.,test="Chisq")

mod2<-glm(HMyco~HABITAT+SITE+YEAR+AGE+SEX+HostSPP*Bartonella,data=test,family="binomial")
drop1(mod2,.~.,test="Chisq")


#model selection in response to reviewer comments (no qualitative impact on results)
g<-glmulti(mod1, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:4,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[3]]
plot_model(best_model, type = "pred",title = "Bartonella")
#best model: "Bartonella ~ 1 + HABITAT + YEAR + HMyco"
mod<-glm(Bartonella ~ HABITAT + YEAR + HostSPP*HMyco,data=test,family="binomial")
drop1(mod,.~.,test="Chisq")

g<-glmulti(mod2, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:10,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[5]]
plot_model(best_model, type = "pred",title = "Mycoplasma")
#best model with pathogens of interest: "HMyco ~ 1 + YEAR + SITE + HABITAT+ HostSPP + Bartonella"
mod<-glm(HMyco ~ YEAR + SITE + HABITAT + HostSPP*Bartonella,data=test,family="binomial")
drop1(mod,.~.,test="Chisq")


# excluding Ra_norv and Ar_sche
mod<-glm(Bartonella~HABITAT+SITE+YEAR+AGE+SEX+HostSPP*HMyco,data=test[which(test$HostSPP!="M_Ra_norv" & test$HostSPP!="C_Ar_sche"),],family="binomial")
summary(mod)
drop1(mod,.~.,test="Chisq")

g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:4,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[1]]
plot_model(best_model, type = "pred",title = "Bartonella")
#best model with pathogens of interest: "Bartonella ~ 1 + YEAR + HABITAT + HMyco"
mod<-glm(Bartonella ~ 1 + YEAR + HABITAT + HostSPP*HMyco,data=test[which(test$HostSPP!="M_Ra_norv" & test$HostSPP!="C_Ar_sche"),],family="binomial")
drop1(mod,.~.,test="Chisq")

mod<-glm(HMyco~HABITAT+SITE+YEAR+AGE+SEX+HostSPP*Bartonella,data=test[which(test$HostSPP!="M_Ra_norv"  & test$HostSPP!="C_Ar_sche"),],family="binomial")
summary(mod)
drop1(mod,.~.,test="Chisq")

g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:10,] %>% regulartable() %>% autofit()
plot(g, type="s")
best_model <- g@objects[[1]]
plot_model(best_model, type = "pred",title = "Mycoplasma")
#best model with pathogens of interest: "HMyco ~ 1 + YEAR + SITE + HABITAT + HostSPP + Bartonella"
mod<-glm(HMyco ~ HABITAT + SEX + YEAR + SITE + HostSPP*Bartonella,data=test[which(test$HostSPP!="M_Ra_norv" & test$HostSPP!="C_Ar_sche"),],family="binomial")
drop1(mod,.~.,test="Chisq")


# logistic regression for each host species

# --> only in C_Ar_sche  (too few)
table(test[which(test$HostSPP=="C_Ar_sche"),]$Bartonella,test[which(test$HostSPP=="C_Ar_sche"),]$HMyco)
# one myco (HC) in bart-infected individual

# --> only in C_Mi_agre (no HC)
mod<-glm(Bartonella~YEAR+HABITAT+SEX+SITE+AGE+HMyco,data=test[which(test$HostSPP=="C_Mi_agre"),],family="binomial")
drop1(mod,.~.,test="Chisq")
summary(mod)
coef(mod) #none
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:3,] %>% regulartable() %>% autofit()
plot(g, type="s")
#best model with pathogens of interest: "Bartonella ~ 1 + AGE + HMyco"
mod<-glm(Bartonella ~ AGE + HMyco,data=test[which(test$HostSPP=="C_Mi_agre"),],family="binomial")
drop1(mod,.~.,test="Chisq")

mod<-glm(HMyco~YEAR+HABITAT+SEX+SITE+AGE+Bartonella,data=test[which(test$HostSPP=="C_Mi_agre"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) #none
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:3,] %>% regulartable() %>% autofit()
plot(g, type="s")
#best model with pathogens of interest: "HMyco ~ 1 + HABITAT + SEX + Bartonella"
mod<-glm(HMyco ~ HABITAT + SEX + Bartonella,data=test[which(test$HostSPP=="C_Mi_agre"),],family="binomial")
drop1(mod,.~.,test="Chisq")
summary(mod)

# --> only in C_Mi_arva / Hedges excluded (no HC)
mod<-glm(Bartonella~SITE+YEAR+AGE+SEX+HMyco,data=test[which(test$HostSPP=="C_Mi_arva" & test$HABITAT!="Hedge"),],family="binomial")
drop1(mod,.~.,test="Chisq")
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:8,] %>% regulartable() %>% autofit()
plot(g, type="s")
#best model with pathogens of interest: "Bartonella ~ 1 + AGE + YEAR + SITE+ HMyco"
mod<-glm(Bartonella ~AGE + YEAR + SITE+ HMyco,data=test[which(test$HostSPP=="C_Mi_arva" & test$HABITAT!="Hedge"),],family="binomial")
drop1(mod,.~.,test="Chisq")
summary(mod)

mod<-glm(HMyco~SITE+YEAR+AGE+SEX+Bartonella,data=test[which(test$HostSPP=="C_Mi_arva" & test$HABITAT!="Hedge"),],family="binomial")
drop1(mod,.~.,test="Chisq")
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:5,] %>% regulartable() %>% autofit()
plot(g, type="s")
#best model with pathogens of interest: "HMyco ~ 1 + YEAR + SITE + Bartonella"
mod<-glm(HMyco ~ 1 + YEAR + SITE + Bartonella,data=test[which(test$HostSPP=="C_Mi_arva" & test$HABITAT!="Hedge"),],family="binomial")
drop1(mod,.~.,test="Chisq")
summary(mod)

table(test[which(test$HostSPP=="C_Mi_arva"),]$HMyco,test[which(test$HostSPP=="C_Mi_arva"),]$HABITAT)
table(test[which(test$HostSPP=="C_Mi_arva"),]$Bartonella,test[which(test$HostSPP=="C_Mi_arva"),]$HABITAT)
#all 3 hedge animals are infected with HMyco (HM)

mod<-glm(Bartonella~SITE+YEAR+AGE+SEX+HMyco,data=test[which(test$HostSPP=="C_Mi_arva"),],family="binomial")
drop1(mod,.~.,test="Chisq")
summary(mod)
mod<-glm(HMyco~SITE+YEAR+AGE+SEX+Bartonella,data=test[which(test$HostSPP=="C_Mi_arva"),],family="binomial")
drop1(mod,.~.,test="Chisq")

# interactions overfit the model

# --> only in C_My_glar
mod<-glm(Bartonella~HABITAT+SITE+YEAR+AGE+SEX+HMyco,data=test[which(test$HostSPP=="C_My_glar"),],family="binomial")
drop1(mod,.~.,test="Chisq") #NEG
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:6,] %>% regulartable() %>% autofit()
plot(g, type="s")
#best model with pathogens of interest: "Bartonella ~ 1 + YEAR + HABITAT + SITE + AGE + HMyco"
mod<-glm(Bartonella ~ YEAR + HABITAT + SITE + AGE + HMyco,data=test[which(test$HostSPP=="C_My_glar"),],family="binomial")
drop1(mod,.~.,test="Chisq")
summary(mod)
exp(cbind(OR = coef(mod), confint(mod)))

mod<-glm(HMyco~HABITAT+SITE+YEAR+AGE+SEX+Bartonella,data=test[which(test$HostSPP=="C_My_glar"),],family="binomial")
drop1(mod,.~.,test="Chisq") #NEG
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:3,] %>% regulartable() %>% autofit()
plot(g, type="s")
#best model with pathogens of interest: "HMyco ~ 1 + YEAR + SITE + Bartonella"
mod<-glm(HMyco ~ YEAR + SITE + Bartonella,data=test[which(test$HostSPP=="C_My_glar"),],family="binomial")
drop1(mod,.~.,test="Chisq")
exp(cbind(OR = coef(mod), confint(mod)))


# --> only in M_Ap_flav
table(test[which(test$HostSPP=="M_Ap_flav"),]$YEAR,test[which(test$HostSPP=="M_Ap_flav"),]$HMyco)
table(test[which(test$HostSPP=="M_Ap_flav"),]$Bartonella,test[which(test$HostSPP=="M_Ap_flav"),]$YEAR)
table(test[which(test$HostSPP=="M_Ap_flav"),]$HABITAT,test[which(test$HostSPP=="M_Ap_flav"),]$HMyco)
table(test[which(test$HostSPP=="M_Ap_flav"),]$Bartonella,test[which(test$HostSPP=="M_Ap_flav"),]$HABITAT)

mod<-glm(Bartonella~HABITAT+SITE+YEAR+AGE+SEX+HMyco,data=test[which(test$HostSPP=="M_Ap_flav"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) #none
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:5,] %>% regulartable() %>% autofit()
plot(g, type="s")
#best model with pathogens of interest: "Bartonella ~ 1 + YEAR + HMyco"
mod<-glm(Bartonella ~ YEAR + HMyco,data=test[which(test$HostSPP=="M_Ap_flav"),],family="binomial")
drop1(mod,.~.,test="Chisq")


mod<-glm(HMyco~HABITAT+SITE+YEAR+AGE+SEX+Bartonella,data=test[which(test$HostSPP=="M_Ap_flav"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) #none
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:7,] %>% regulartable() %>% autofit()
plot(g, type="s")
#best model with pathogens of interest: "HMyco ~ 1 + Bartonella"
mod<-glm(HMyco ~ Bartonella,data=test[which(test$HostSPP=="M_Ap_flav"),],family="binomial")
drop1(mod,.~.,test="Chisq")


# --> only in M_Ap_sylv
table(test[which(test$HostSPP=="M_Ap_sylv"),]$YEAR,test[which(test$HostSPP=="M_Ap_sylv"),]$HMyco)
table(test[which(test$HostSPP=="M_Ap_sylv"),]$Bartonella,test[which(test$HostSPP=="M_Ap_sylv"),]$YEAR)
table(test[which(test$HostSPP=="M_Ap_sylv"),]$HABITAT,test[which(test$HostSPP=="M_Ap_sylv"),]$HMyco)
table(test[which(test$HostSPP=="M_Ap_sylv"),]$Bartonella,test[which(test$HostSPP=="M_Ap_sylv"),]$HABITAT)

mod<-glm(Bartonella~YEAR+AGE+SITE+SEX+HABITAT*HMyco,data=test[which(test$HostSPP=="M_Ap_sylv"),],family="binomial")
drop1(mod,.~.,test="Chisq") 
coef(mod) # trend for interaction with habitat, only one Bart-uninfected in Forest
summary(mod)
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:10,] %>% regulartable() %>% autofit()
plot(g, type="s")
#best model with pathogens of interest: "Bartonella ~ 1 + YEAR + AGE+SEX+HMyco"
mod<-glm(Bartonella ~ YEAR +AGE + SEX + HMyco,data=test[which(test$HostSPP=="M_Ap_sylv"),],family="binomial")
drop1(mod,.~.,test="Chisq")


mod<-glm(HMyco~HABITAT+SITE+YEAR+AGE+SEX+Bartonella,data=test[which(test$HostSPP=="M_Ap_sylv"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) #none
g<-glmulti(mod, family="binomial", crit="aicc", level=1)
print(g)
weightable(g)[1:5,] %>% regulartable() %>% autofit()
plot(g, type="s")
#best model with pathogens of interest: "HMyco ~ 1 + Bartonella"
mod<-glm(HMyco ~ Bartonella,data=test[which(test$HostSPP=="M_Ap_flav"),],family="binomial")
drop1(mod,.~.,test="Chisq")


# --> only in M_Ra_norv
table(test[which(test$HostSPP=="M_Ra_norv"),]$Bartonella,test[which(test$HostSPP=="M_Ra_norv"),]$HMyco)

mod<-glm(Bartonella~SITE+YEAR+AGE+SEX+HMyco,data=test[which(test$HostSPP=="M_Ra_norv"),],family="binomial")
drop1(mod,.~.,test="Chisq")
coef(mod) #none
summary(mod)
weightable(glmulti(mod, family="binomial", level=1))
#best model with pathogens of interest: "Bartonella ~ 1 + SITE + HMyco"
mod<-glm(Bartonella ~ SITE + HMyco,data=test[which(test$HostSPP=="M_Ra_norv"),],family="binomial")
drop1(mod,.~.,test="Chisq")

mod<-glm(HMyco~SITE+YEAR+AGE+SEX+Bartonella,data=test[which(test$HostSPP=="M_Ra_norv"),],family="binomial")
drop1(mod,.~.,test="Chisq") # CAN"T BE DONE, ONLY 1 uninfected by HMyco
coef(mod) #none
weightable(glmulti(mod, family="binomial", level=1))
#best model with pathogens of interest: " HMyco ~ 1 + YEAR + Bartonella"
mod<-glm( HMyco ~ YEAR + Bartonella,data=test[which(test$HostSPP=="M_Ra_norv"),],family="binomial")
drop1(mod,.~.,test="Chisq")


##############################################################
#   3.4   False Discovery  
##############################################################
# Benjamini-Hochberg correction
pvals_reduced<-c(0.0001,
                 0.0096,
                 0.51,
                 0.24,
                 0.0001,
                 0.27,
                 0.032,
                 0.001,
                 0.024,
                 0.23,
                 0.15,
                 0.14,
                 0.67,
                 0.0000000000001,
                 0.00000001,
                 0.032,
                 0.92,
                 0.0001,
                 0.51,
                 0.21,
                 0.63,
                 0.013,
                 0.55,
                 0.002,
                 0.08,
                 0.1,
                 0.97,
                 0.52,
                 0.35,
                 0.0001,
                 0.49,
                 0.05,
                 0.032,
                 0.5,
                 0.97,
                 0.0001,
                 0.009,
                 0.18,
                 0.065,
                 0.053,
                 0.63,
                 0.21,
                 0.14,
                 0.57,
                 0.35,
                 0.37,
                 0.47,
                 0.0032,
                 0.44,
                 0.59,
                 0.17,
                 0.64,
                 0.94,
                 0.0045,
                 0.26,
                 0.068,
                 0.001,
                 0.45,
                 0.27,
                 0.8,
                 0.001,
                 0.017,
                 0.12,
                 0.021,
                 0.42,
                 0.39,
                 0.45,
                 0.053,
                 0.92,
                 0.36,
                 0.0021,
                 0.24,
                 0.4,
                 0.32,
                 0.017,
                 0.21,
                 0.67)
df<-as.data.frame(pvals_reduced)
df$pvals_BH<-p.adjust(pvals_reduced,method="BH")
max(df$pvals_reduced[which(df$pvals_BH<0.05)])
min(df$pvals_reduced[which(df$pvals_BH>=0.05)])

