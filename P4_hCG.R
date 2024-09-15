# Set theme----
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# Library packages----
library(readxl)
library(tidyverse)
library(mgcv)
library(geepack)
library(gtsummary)

# Import data----
df_antagonist<-read_xlsx("./df_antagonist.xlsx")
df_antagonist<-within(df_antagonist,{
  P10<-P_hCG_day*10#tranform of P4-hCG unit
})
df_antagonist_SB<-subset(df_antagonist,ET_protocol=="SB")
df_antagonist_SC<-subset(df_antagonist,ET_protocol=="SC")
df_antagonist_DC<-subset(df_antagonist,ET_protocol=="DC")

# Table 1
table1<-descrTable(ET_protocol~Female_age+BMI+Infertile_year+AMH+Gn_start_dose+Gn_total_dose+Gn_total_day+No_of_oocyte+Endometrial_thickness_hCG_day+E2_hCG_day+LH_hCG_day+P_hCG_day+No_of_ET+Embryo_phase+No_of_good_embryo+Biochemical_pregnancy_cat+Biochemical_pregnancy_loss_cat+Clinical_pregnancy_cat+Pregnancy_type+Multiple_pregnancy,data=df_antagonist,method="1",show.p.mul=T,digits=2,show.all=F, p.correcte=F)
table1
export2word(table1, file='./table1.docx')

# GAM For BP----
## all antagonist protocol----
### Smooth line between P-hCG level and BP
fml1<-'Biochemical_pregnancy~s(P_hCG_day)+Female_age+No_of_ET+Embryo_phase+No_of_good_embryo'
gam1<-mgcv::gam(formula(fml1),data=subset(df_antagonist,P_hCG_day<=1.5), family=binomial(link="logit"))
pdf(file='./gam_smoothline.antagonist.all.BP.P4.pdf',width=3.5,height=4)
plot(gam1)#Figure 1A
dev.off()

### Smooth line between P-hCG level and BP
fml2<-'Biochemical_pregnancy~P_hCG_day+s(Female_age)+No_of_ET+Embryo_phase+No_of_good_embryo'
gam2<-mgcv::gam(formula(fml2),data=subset(df_antagonist,P_hCG_day<=1.5),family=binomial(link="logit"))
pdf(file='./gam_smoothline.antagonist.all.BP.age.pdf',width=3.5,height=4)
plot(gam2)#Figure 1B
dev.off()

### Model1 treated P4-hCG as continuous variable
fml3<-'Biochemical_pregnancy~P10+Female_age+No_of_ET+Embryo_phase+No_of_good_embryo'
gam3<-mgcv::gam(formula(fml3),data=subset(df_antagonist,P_hCG_day<=1.5), family=binomial(link="logit"))
summary(gam3)
tbl_regression(gam3, exponentiate = TRUE)#Figure 1C

### Model2 treated P4-hCG as categorial variable
fml4<-'Biochemical_pregnancy~P.seg.1+Female_age_cat+No_of_ET+Embryo_phase+No_of_good_embryo'
gam4<-mgcv::gam(formula(fml4),data=subset(df_antagonist,P_hCG_day<=1.5), family=binomial(link="logit"))
tbl_regression(gam4, exponentiate = TRUE)#Figure 1D

## SB cycles----
### Model1 treated P4-hCG as continuous variable
fml5<-'Biochemical_pregnancy~P10+Female_age+No_of_good_embryo'
gam5 <-mgcv::gam(formula(fml5),data=subset(df_antagonist_SB,P_hCG_day<1.5), family=binomial(link="logit"))
tbl_regression(gam5, exponentiate = TRUE)#Figure 2A

### Model2 treated P4-hCG as categorial variable
fml6<-'Biochemical_pregnancy~P.seg.1+Female_age_cat+No_of_good_embryo'
gam6 <-mgcv::gam(formula(fml6),data=subset(df_antagonist_SB, P_hCG_day<=1.5), family=binomial(link="logit"))
tbl_regression(gam6, exponentiate = TRUE)#Figure 2B

## SC cycles----
### Model1 treated P4-hCG as continuous variable
fml7<-'Biochemical_pregnancy~P10+Female_age+No_of_good_embryo'
gam7 <-mgcv::gam(formula(fml7),data=subset(df_antagonist_SC,P_hCG_day<1.5), family=binomial(link="logit"))
tbl_regression(gam7, exponentiate = TRUE)#Figure 3A

### Model2 treated P4-hCG as categorial variable
fml8<-'Biochemical_pregnancy~P.seg.1+Female_age_cat+No_of_good_embryo'
gam8 <-mgcv::gam(formula(fml8),data=subset(df_antagonist_SC,P_hCG_day<1.5), family=binomial(link="logit"))
tbl_regression(gam8, exponentiate = TRUE)#Figure 3B

### Model3 treated P4-hCG as categorial variable, P4-hCG was split to <=0.7,0.7-1.0, and >1.0
fml9<-'Biochemical_pregnancy~P.seg.2+Female_age_cat+No_of_good_embryo'
gam9 <-mgcv::gam(formula(fml9),data=subset(df_antagonist_SC,P_hCG_day<1.5), family=binomial(link="logit"))
tbl_regression(gam9, exponentiate = TRUE)#Figure 3C

## DC cycles----
### Model1 treated P4-hCG as continuous variable
fml10<-'Biochemical_pregnancy~P10+Female_age+No_of_good_embryo'
gam10 <-mgcv::gam(formula(fml10),data=subset(df_antagonist_DC,P_hCG_day<=1.5), family=binomial(link="logit"))
tbl_regression(gam10, exponentiate = TRUE)#Figure 4A

### Model2 treated P4-hCG as categorial variable
fml11<-'Biochemical_pregnancy~P.seg.1+Female_age_cat+No_of_good_embryo'
gam11 <-mgcv::gam(formula(fml11),data=subset(df_antagonist_DC, P_hCG_day<1.5), family=binomial(link="logit"))
tbl_regression(gam11, exponentiate = TRUE)#Figure 4B

### Model3 treated P4-hCG as categorial variable, P4-hCG was split to <=0.7,0.7-1.0, and >1.0
fml12<-'Biochemical_pregnancy~P.seg.2+Female_age_cat+No_of_good_embryo'
gam12 <-mgcv::gam(formula(fml12),data=subset(df_antagonist_DC,P_hCG_day<1.5), family=binomial(link="logit"))
tbl_regression(gam12, exponentiate = TRUE)#Figure 4C


# GAM For CP----
## all antagonist protocol----
### Smooth line between P-hCG level and BP
fml13<-'Clinical_pregnancy~s(P_hCG_day)+Female_age+No_of_ET+Embryo_phase+No_of_good_embryo'
gam13<-mgcv::gam(formula(fml13),data=subset(df_antagonist,P_hCG_day<=1.5), family=binomial(link="logit"))
summary(gam13)
pdf(file='./gam_smoothline.antagonist.all.CP.P4.pdf',width=3.5,height=4)
plot(gam6)#Figure 1E
dev.off()

### Smooth line between age level and BP
fml14<-'Clinical_pregnancy~P_hCG_day+s(Female_age)+No_of_ET+Embryo_phase+No_of_good_embryo'
gam14<-mgcv::gam(formula(fml14),data=subset(df_antagonist,P_hCG_day<=1.5), family=binomial(link="logit"))
summary(gam14)
pdf(file='./gam_smoothline.antagonist.all.CP.age.pdf',width=3.5,height=4)
plot(gam14)#Figure 1F
dev.off()

### Model1 treated P4-hCG as continuous variable
fml15<-'Clinical_pregnancy~P10+Female_age+No_of_ET+Embryo_phase+No_of_good_embryo'
gam15<-mgcv::gam(formula(fml15),data=subset(df_antagonist,P_hCG_day<=1.5), family=binomial(link="logit"))
tbl_regression(gam15, exponentiate = TRUE)#Figure1G

### Model2 treated P4-hCG as categorial variable
fml16<-'Clinical_pregnancy~P.seg.1+Female_age_cat+No_of_ET+Embryo_phase+No_of_good_embryo'
gam16<-mgcv::gam(formula(fml16),data=subset(df_antagonist,P_hCG_day<=1.5), family=binomial(link="logit"))
tbl_regression(gam16, exponentiate = TRUE)#Figure1H

## SB cycles
### Model1
fml17<-'Clinical_pregnancy~P10+Female_age+No_of_good_embryo'
gam17 <-mgcv::gam(formula(fml17),data=subset(df_antagonist_SB,P_hCG_day<1.5), family=binomial(link="logit"))
tbl_regression(gam17, exponentiate = TRUE)#Figure 2C

### Model2
fml18<-'Clinical_pregnancy~P.seg.1+Female_age_cat+No_of_good_embryo'
gam18 <-mgcv::gam(formula(fml18),data=subset(df_antagonist_SB,P_hCG_day<1.5), family=binomial(link="logit"))
tbl_regression(gam18, exponentiate = TRUE)#Figure 2D

## SC cycles 
### Model1 
fml19<-'Clinical_pregnancy~P10+Female_age+No_of_good_embryo'
gam19 <-mgcv::gam(formula(fml19),data=subset(df_antagonist_SC,P_hCG_day<1.5), family=binomial(link="logit"))
tbl_regression(gam19, exponentiate = TRUE)#Figure 3D

### Model2
fml20<-'Clinical_pregnancy~P.seg.1+Female_age_cat+No_of_good_embryo'
gam20 <-mgcv::gam(formula(fml20),data=subset(df_antagonist_SC,P_hCG_day<1.5), family=binomial(link="logit"))
tbl_regression(gam20, exponentiate = TRUE)#Figure 3E

### Model3
fml21<-'Clinical_pregnancy~P.seg.2+Female_age_cat+No_of_good_embryo'
gam21 <-mgcv::gam(formula(fml21),data=subset(df_antagonist_SC,P_hCG_day<1.5), family=binomial(link="logit"))
tbl_regression(gam21, exponentiate = TRUE)#Figure 3F

## DC cycles
### Model1
fml22<-'Clinical_pregnancy~P10+Female_age+No_of_good_embryo'
gam22 <-mgcv::gam(formula(fml22),data=subset(df_antagonist_DC,P_hCG_day<=1.5), family=binomial(link="logit"))
tbl_regression(gam22, exponentiate = TRUE)#Figure 4D

### Model2
fml23<-'Clinical_pregnancy~P.seg.1+Female_age_cat+No_of_good_embryo'
gam23 <-mgcv::gam(formula(fml23),data=subset(df_antagonist_DC,P_hCG_day<=1.5), family=binomial(link="logit"))
tbl_regression(gam23, exponentiate = TRUE)#Figure 4E

### Model3
fml24<-'Clinical_pregnancy~P.seg.2+Female_age_cat+No_of_good_embryo'
gam24 <-mgcv::gam(formula(fml24),data=subset(df_antagonist_DC,P_hCG_day<=1.5), family=binomial(link="logit"))
tbl_regression(gam24, exponentiate = TRUE)#Figure 4F

# outcome comparisons bewteen P-hCG segmentations in SB, SC, and DC groups-----
## SB cycles----
## Table 2
Table.SB.all<-descrTable(P.seg.1~Biochemical_pregnancy_cat+Clinical_pregnancy_cat+Biochemical_pregnancy_loss_cat+Multiple_pregnancy+Endometrial_thickness_hCG_day+Female_age+AMH+Gn_start_dose+Gn_total_dose+Gn_total_day+No_of_oocyte,data=df_antagonist_SB,method="1",show.p.mul=T,digits=2,p.corrected = F)
Table.SB.all
export2word(Table.SB.all, file='./Table.SB.all.docx')

## SC cycles----
## Table 3
Table.SC.all<-descrTable(P.seg.2~Biochemical_pregnancy_cat+Clinical_pregnancy_cat+Biochemical_pregnancy_loss_cat+Multiple_pregnancy+Endometrial_thickness_hCG_day+Female_age+AMH+Gn_start_dose+Gn_total_dose+Gn_total_day+No_of_oocyte,data=df_antagonist_SC,method="1",show.p.mul=T,digits=2,p.corrected = F)
Table.SC.all
export2word(Table.SC.all, file='./Table.SC.all.docx')

## Table 4
Table.SC.good<-descrTable(P.seg.2~Biochemical_pregnancy_cat+Clinical_pregnancy_cat+Biochemical_pregnancy_loss_cat+Multiple_pregnancy+Endometrial_thickness_hCG_day+Female_age+AMH+Gn_start_dose+Gn_total_dose+Gn_total_day+No_of_oocyte,data=df_antagonist_SC,method="1",show.p.mul=T,digits=2,subset=No_of_good_embryo==1,p.corrected = F)
Table.SC.good
export2word(Table.SC.good, file='./Table.SC.good.docx')

## DC cycles----
## Table 5
Table.DC.all<-descrTable(P.seg.2~Biochemical_pregnancy_cat+Clinical_pregnancy_cat+Biochemical_pregnancy_loss_cat+Multiple_pregnancy+Endometrial_thickness_hCG_day+Female_age+AMH+Gn_start_dose+Gn_total_dose+Gn_total_day+No_of_oocyte,data=df_antagonist_DC,method="1",show.p.mul=T,digits=2,p.corrected = F)
Table.DC.all
export2word(Table.DC.all, file='./Table.DC.all.2.docx')

## Table 6
Table.DC.good<-descrTable(P.seg.2~Biochemical_pregnancy_cat+Clinical_pregnancy_cat+Biochemical_pregnancy_loss_cat+Multiple_pregnancy+Endometrial_thickness_hCG_day+Female_age+AMH+Gn_start_dose+Gn_total_dose+Gn_total_day+No_of_oocyte,data=df_antagonist_DC,method="1",show.p.mul=T,digits=2,subset=No_of_good_embryo==2,p.corrected = F)
Table.DC.good
export2word(Table.DC.good, file='./Table.DC.good.2.docx')


# Comparisons between ET strategies in P-hCG<=0.7ng/ml----
## Table 7
Table.P.0.7<-descrTable(ET_protocol~Biochemical_pregnancy_cat+Clinical_pregnancy_cat+Biochemical_pregnancy_loss_cat+Multiple_pregnancy,
                        data=df_antagonist, subset=P.seg.1=="<=0.7",method="1",show.p.mul=T,digits=2,p.corrected = F)
Table.P.0.7
export2word(Table.P.0.7, file='./Table7.docx')

# Comparisons between ET strategies in P-hCG>0.7ng/ml----
## Table 8
Table.P.1<-descrTable(ET_protocol~Biochemical_pregnancy_cat+Clinical_pregnancy_cat+Biochemical_pregnancy_loss_cat+Multiple_pregnancy,
                      data=df_antagonist, subset=P.seg.1==">0.7",method="1",show.p.mul=T,digits=2,p.corrected = F)
Table.P.1
export2word(Table.P.1, file='./Table8.docx')

# Statistical power estimation----
## power estimation for BP----
### SB cycles
library(pwr)
prob.SB.BP<-prop.table(table(df_antagonist_SB$P.seg.1,df_antagonist_SB$Biochemical_pregnancy_cat))
prob.SB.BP<-matrix(prob.SB.BP,byrow = F,nrow = 2)
prob.SB.BP
library(pwr)
sig<- c(0.05, 0.10, 0.15)
power.SB.BP <- NULL
for (i in 1:length(sig)){
  result <- pwr.chisq.test(w=0.20, N=nrow(df_antagonist_SB),df=(nrow(prob.SB.BP)-1)*(ncol(prob.SB.BP)-1), sig.level=sig[i])
  power.SB.BP[i] <- result$power
}
power.SB.BP

### SC cycles----
prob.SC.BP<-prop.table(table(df_antagonist_SC$P.seg.2,df_antagonist_SC$Biochemical_pregnancy_cat))
prob.SC.BP<-matrix(prob.SC.BP,byrow = F,nrow = 3)
prob.SC.BP
power.SC.BP <- NULL
for (i in 1:length_sig){
  result <- pwr.chisq.test(w=0.1, N=nrow(df_antagonist_SC),df=(nrow(prob.SC.BP)-1)*(ncol(prob.SC.BP)-1), sig.level=sig[i])
  power.SC.BP[i] <- result$power
}
power.SC.BP

### DC cycles----
prob.DC.BP<-prop.table(table(df_antagonist_DC$P.seg.2,df_antagonist_DC$Biochemical_pregnancy_cat))
prob.DC.BP<-matrix(prob.DC.BP,byrow = F,nrow = 3)
prob.DC.BP
power.DC.BP<- NULL
for (i in 1:length_sig){
  result <- pwr.chisq.test(w=0.1, N=nrow(df_antagonist_DC),df=(nrow(prob.DC.BP)-1)*(ncol(prob.DC.BP)-1), sig.level=sig[i])
  power.DC.BP[i] <- result$power
}
power.DC.BP

## Power curve for BP----
pdf(file="./Power.BP.pdf",width=5,heigh=5)
xrange <- c(0.05,0.15)                  
yrange <- c(0.70,1.0)
plot(xrange, yrange, type="n",
     xlab="Significance leve",
     ylab="Power" )
lines(sig, power.SB.BP, type = "l", lwd=2, col="pink2")
lines(sig, power.SC.BP, type = "l", lwd=2, col="skyblue")
lines(sig, power.DC.BP, type = "l", lwd=2, col="green3")
title("Statistical Power Estimation Curve\n for Each ET strategy on BP")                   
legend(x=0.12, y=0.85, title="ET strategy", legend = c("SB", "SC", "DC"),fill = c("pink2", "skyblue", "green3"), cex=1) 
dev.off()

## power estimation for CP----
### SB cycles----
prob.SB.CP<-prop.table(table(df_antagonist_SB$P.seg.1,df_antagonist_SB$Clinical_pregnancy_cat))
prob.SB.CP<-matrix(prob.SB.CP,byrow = F,nrow = 2)
prob.SB.CP
sig <- c(0.05, 0.10, 0.15)
power.SB.CP <- NULL
for (i in 1:length(sig)){
  result <- pwr.chisq.test(w=0.15, N=nrow(df_antagonist_SB),df=(nrow(prob.SB.CP)-1)*(ncol(prob.SB.CP)-1), sig.level=sig[i])
  power.SB.CP[i] <- result$power}
power.SB.CP

### SC cycles----
prob.SC.CP<-prop.table(table(df_antagonist_SC$P.seg.2,df_antagonist_SC$Clinical_pregnancy_cat))
prob.SC.CP<-matrix(prob.SC.CP,byrow = F,nrow = 3)
prob.SC.CP
sig <- c(0.05, 0.10, 0.15)
power.SC.CP <- NULL
for (i in 1:length(sig)){
  result <- pwr.chisq.test(w=0.08, N=nrow(df_antagonist_SC),df=(nrow(prob.SC.CP)-1)*(ncol(prob.SC.CP)-1), sig.level=sig[i])
  power.SC.CP[i] <- result$power}
power.SC.CP

### DC cycles----
prob.DC.CP<-prop.table(table(df_antagonist_DC$P.seg.2,df_antagonist_DC$Clinical_pregnancy_cat))
prob.DC.CP<-matrix(prob.DC.CP,byrow = F,nrow = 3)
prob.DC.CP
sig <- c(0.05, 0.10, 0.15)
power.DC.CP <- NULL
for (i in 1:length(sig)){
  result <- pwr.chisq.test(w=0.08, N=nrow(df_antagonist_DC),df=(nrow(prob.DC.CP)-1)*(ncol(prob.DC.CP)-1), sig.level=sig[i])
  power.DC.CP[i] <- result$power}
power.DC.CP

## Power curve for CP----
pdf(file="./Power.CP.pdf",width=5,heigh=5)
xrange <- c(0.05,0.15)                  
yrange <- c(0.50,1.0)
plot(xrange, yrange, type="n",
     xlab="Significance leve",
     ylab="Power" )
lines(sig, power.SB.CP, type = "l", lwd=2, col="pink2",lty=2)
lines(sig, power.SC.CP, type = "l", lwd=2, col="skyblue",lty=2)
lines(sig, power.DC.CP, type = "l", lwd=2, col="green3",lty=2)
title("Statistical Power Estimation Curve\n for Each ET strategy on CP")                   
legend(x=0.12, y=0.69, title="ET strategy", legend = c("SB", "SC", "DC"),fill = c("pink2", "skyblue", "green3"), cex=1,lty=2) 
dev.off()

