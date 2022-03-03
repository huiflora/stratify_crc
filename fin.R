rm(list=ls())
library("survminer")
require("survival")
library(survminer)
library("readxl")
library("My.stepwise") #tepwise Variable Selection Procedures for Regression Analysis
library("ggfortify") #Data Visualization Tools for Statistical Analysis Results
library("MASS")
library("qqman")
library("plyr")
library(reportRx) #Table One
library(mclust) #Model-based cluster
library(fpc) #Flexible Procedures for Clustering


####Table One####
setwd("Z:/New/GMM_730")
dt <- read.csv(file = "co20.csv", header=TRUE)
dt$ALLO1_CD<-as.factor(dt$ALLO1_CD)
dt$os_event<-as.factor(dt$os_event)
dt$sitenum<-as.factor(dt$sitenum)
dt$pcnum<-as.factor(dt$pcnum)
dt$livmeta<-as.factor(dt$livmeta)
dt$cetux_duration<-as.numeric(dt$cetux_duration)
dt$briv_duration<-as.numeric(dt$briv_duration)
dt$WEC_PERF_STA<-as.factor(dt$WEC_PERF_STA)
a<-covsum(data = dt,covs = c("os", "os_event", "response", "age", "SEX", 
                             "WEC_MALIG", "sitenum", "pcnum", "livmeta", "WEC_HEIGHT", "WEC_WEIGHT", 
                             "WEC_PERF_STA", "KRAS_RESULT", "cetux_duration", "briv_duration" ),maincov = "trt")

##
d1 <- read.csv(file = "co20_1.csv", header=TRUE)
d1$group<-as.factor(d1$group)
d1$os_event<-as.factor(d1$os_event)
d1$sitenum<-as.factor(d1$sitenum)
d1$pcnum<-as.factor(d1$pcnum)
d1$livmeta<-as.factor(d1$livmeta)
d1$cetux_duration<-as.numeric(d1$cetux_duration)
d1$briv_duration<-as.numeric(d1$briv_duration)
d1$WEC_PERF_STA<-as.factor(d1$WEC_PERF_STA)
b<-covsum(data = d1,covs = c("os", "os_event", "response", "age", "SEX", 
                             "WEC_MALIG", "sitenum", "pcnum", "livmeta", "WEC_HEIGHT", "WEC_WEIGHT", 
                             "WEC_PERF_STA", "KRAS_RESULT", "cetux_duration", "briv_duration" ),maincov = "group", testcat ="Chi-squared")

d2 <- read.csv(file = "co20_2.csv", header=TRUE)
d2$ALLO1_CD<-as.factor(d2$ALLO1_CD)
d2$os_event<-as.factor(d2$os_event)
d2$sitenum<-as.factor(d2$sitenum)
d2$pcnum<-as.factor(d2$pcnum)
d2$livmeta<-as.factor(d2$livmeta)
d2$cetux_duration<-as.numeric(d2$cetux_duration)
d2$briv_duration<-as.numeric(d2$briv_duration)
d2$WEC_PERF_STA<-as.factor(d2$WEC_PERF_STA)
a<-covsum(data = d2,covs = c("os", "os_event", "response", "age", "SEX", 
                             "WEC_MALIG", "sitenum", "pcnum", "livmeta", "WEC_HEIGHT", "WEC_WEIGHT", 
                             "WEC_PERF_STA", "KRAS_RESULT", "cetux_duration", "briv_duration" ),maincov = "trt")

####K-M plot for os by trt####
fit <- survfit(Surv(os, os_event) ~ trt, data = dt)
ggsurvplot(fit, data=dt, title = "Overall Survival",
           legend=c(0.9,0.9),
           xlab = "Time in months",
           legend.labs =c("Brivanib + Cet (n=376)","Placebo +Cet (n=374)"),
           risk.table=TRUE, risk.table.y.text =FALSE,
           pval=TRUE,size =1)

######Summary score of QLQ-C30########
library(PROscorer)
qol <- read.csv(file = "~/Desktop/R/Nco20_qol.csv",header = T)
sum(is.na(qol))
dat<-qol[,c(3:32)]
a<-qlq_c30(dat,iprefix = "QOL_")
dat_scored <- data.frame(qol, a)
str(dat_scored)
#write.csv(dat_scored, file = "qol_sum_2.csv")
#dat_scored $na_count <- apply(dat_scored[,c(4:33)] , 1, function(x) sum(is.na(x)))

###K-M curve for GMM group#########
dt$os_event<-as.numeric(dt$os_event)
d1<-dt[dt$ALLO1_CD==1,]
fit <- survfit(Surv(os, os_event) ~ group , data = d1)
summary(coxph(Surv(os, os_event) ~ group , data = d1))
ggsurv<-ggsurvplot(fit, data=d1, title = "Overall Survival",
           legend=c(0.9,0.9),
           xlab = "Time in months",
           legend.labs =c("Group 1","Group 2","Group 3"),
           risk.table=TRUE, risk.table.y.text =FALSE,
           palette =c("#f1948a","#a9cce3","#52be80"),
           size =1,
           conf.int = TRUE)
ggsurv
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 7, y = 0.05, 
                    label = "log-rank p = 0.00058,\n HR for group 2 vs group 1 = 0.61,\n HR for group 3 vs group 1 = 0.94", size = 3)
ggsurv

d2<-dt[dt$ALLO1_CD==2,]
fit1 <- survfit(Surv(os, os_event) ~ group , data = d2)
summary(coxph(Surv(os, os_event) ~ group , data = d2))
ggsurv<-ggsurvplot(fit1, data=d2, title = "Overall Survival",
                   legend=c(0.9,0.9),
                   xlab = "Time in months",
                   legend.labs =c("Group I","Group II","Group III"),
                   risk.table=TRUE, risk.table.y.text =FALSE,
                   palette =c("#f5b041","#4dbda7","#4d64bd"),
                   pval=F,size =1,
                   conf.int = TRUE)
ggsurv
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 5, y = 0.05, 
                    label = "log-rank p < 0.0001, \n HR for group I vs group III = 2.36,\n HR for group II vs group III = 0.65", size = 3)
ggsurv

####K-M plot for os by SNP####
rs7336411_0<-dt[dt$rs7336411_A==0,]
fit <- survfit(Surv(os, os_event) ~ base, data = rs7336411_0)
ggsurvplot(fit, data=dt, title = "Overall Survival for rs7336411_0",
           legend=c(0.9,0.9),
           xlab = "Time in months",
           legend.labs =c("High Baseline QoL","Low Baseline QoL"),
           risk.table=TRUE, risk.table.y.text =FALSE,
           pval=TRUE,size =1)

rs7336411_1<-dt[dt$rs7336411_A==1,]
fit <- survfit(Surv(os, os_event) ~ base, data = rs7336411_1)
ggsurvplot(fit, data=dt, title = "Overall Survival for rs7336411_1",
           legend=c(0.9,0.9),
           xlab = "Time in months",
           legend.labs =c("High Baseline QoL","Low Baseline QoL"),
           risk.table=TRUE, risk.table.y.text =FALSE,
           pval=TRUE,size =1)

rs7336411_2<-dt[dt$rs7336411_A==2,]
fit <- survfit(Surv(os, os_event) ~ base, data = rs7336411_2)
ggsurvplot(fit, data=dt, title = "Overall Survival for rs7336411_2",
           legend=c(0.9,0.9),
           xlab = "Time in months",
           legend.labs =c("High Baseline QoL","Low Baseline QoL"),
           risk.table=TRUE, risk.table.y.text =FALSE,
           pval=TRUE,size =1)

#########Trt 1
setwd("Z:/New/GMM_730/trt 1")
trt1<- read.table(file = "2.dat", header=F)
plot(tree[,1], tree[,2],type="n", xlab="Weeks",
     ylab="QoL", xlim=c(0, 9), ylim=c(12, 105) )
for (i in 5:366) {
  tree <- trt1[,c(1,i)]
  lines(tree[,1], tree[,2], type="b", lwd=1,
        lty=3, pch=1)
}
lines(trt1[,1], trt1[,2], type="b", lwd=2,
      lty=1, pch=0,col="#f1948a")
lines(trt1[,1], trt1[,3], type="b", lwd=2,
      lty=1, pch=0, col="#a9cce3")
lines(trt1[,1], trt1[,4], type="b", lwd=2,
      lty=1, pch=0, col="#52be80")
title("Estimated Means and Observed Individual Values for treatment 1")
legend(8, 105, legend=c("Group 1","Group 2","Group 3"),
       col=c("#f1948a","#a9cce3","#52be80"), lty=1, cex=0.5)

####
setwd("Z:/New/GMM_730/trt 1")
trt1_2<- read.table(file = "3.dat", header=F)
# set up the plot
plot(trt1_2[,1], trt1_2[,2],type="n", xlab="Weeks",
     ylab="QoL", xlim=c(0, 8), ylim=c(20, 90) )
# add lines
lines(trt1_2[,1], trt1_2[,2], type="b", lwd=2,
      lty=1, pch=0,col="#f1948a")
lines(trt1_2[,1], trt1_2[,3], type="b", lwd=2,
      lty=2, pch=0, col="#f1948a")
lines(trt1_2[,1], trt1_2[,4], type="b", lwd=2,
      lty=1, pch=0,col="#a9cce3")
lines(trt1_2[,1], trt1_2[,5], type="b", lwd=2,
      lty=2, pch=0, col="#a9cce3")
lines(trt1_2[,1], trt1_2[,6], type="b", lwd=2,
      lty=1, pch=0,col="#52be80")
lines(trt1_2[,1], trt1_2[,7], type="b", lwd=2,
      lty=2, pch=0, col="#52be80")
title("Estimated means and observed means for treatment 1")
legend(0,32, legend=c("Group 1 sample means","Group 1 estimated means","Group 2 sample means","Group 2 estimated means","Group 3 sample means","Group 3 estimated means"),
       col=c("#f1948a","#f1948a","#a9cce3","#a9cce3","#52be80","#52be80"), lty=1:2, cex=0.65)

#########Trt 2########
setwd("Z:/New/GMM_730/trt 2")
trt2<- read.table(file = "2.dat", header=F)
plot(tree[,1], tree[,2],type="n", xlab="Weeks",
     ylab="QoL", xlim=c(0, 9), ylim=c(12, 105) )
for (i in 5:371) {
  tree <- trt2[,c(1,i)]
  lines(tree[,1], tree[,2], type="b", lwd=1,
        lty=3, pch=1)
}
lines(trt2[,1], trt2[,2], type="b", lwd=2,
      lty=1, pch=0,col="#f5b041")
lines(trt2[,1], trt2[,3], type="b", lwd=2,
      lty=1, pch=0, col="#4dbda7")
lines(trt2[,1], trt2[,4], type="b", lwd=2,
      lty=1, pch=0, col="#4d64bd")
title("Estimated Means and Observed Individual Values for treatment 2")
legend(8, 105, legend=c("Group I","Group II","Group III"),
       col=c("#f5b041","#4dbda7","#4d64bd"), lty=1, cex=0.5)

####
trt2_2<- read.table(file = "3.dat", header=F)
# set up the plot
plot(trt2_2[,1], trt2_2[,2],type="n", xlab="Weeks",
     ylab="QoL", xlim=c(0, 9), ylim=c(20, 105) )
# add lines
lines(trt2_2[,1], trt2_2[,2], type="b", lwd=2,
      lty=1, pch=0,col="#f5b041")
lines(trt2_2[,1], trt2_2[,3], type="b", lwd=2,
      lty=2, pch=0, col="#f5b041")
lines(trt2_2[,1], trt2_2[,4], type="b", lwd=2,
      lty=1, pch=0,col="#4dbda7")
lines(trt2_2[,1], trt2_2[,5], type="b", lwd=2,
      lty=2, pch=0, col="#4dbda7")
lines(trt2_2[,1], trt2_2[,6], type="b", lwd=2,
      lty=1, pch=0,col="#4d64bd")
lines(trt2_2[,1], trt2_2[,7], type="b", lwd=2,
      lty=2, pch=0, col="#4d64bd")
title("Estimated means and observed means for treatment 2")
legend(0.5, 40, legend=c("Group I sample means","Group I estimated means","Group II sample means","Group II estimated means","Group III sample means","Group III estimated means"),
       col=c("#f5b041","#f5b041","#4dbda7","#4dbda7","#4d64bd","#4d64bd"), lty=1:2, cex=0.65)



#############t-test for two sample############
rm(list=ls())
sum_q<-read.csv(file = "~/Desktop/R/qol_sum.csv",header = TRUE)
a<-sum_q[sum_q$SCHED_TM==0,]
i<-t.test(Sum ~ trt, data = a)
b<-sum_q[sum_q$SCHED_TM==2,]
ii<-t.test(Sum ~ trt, data = b)
c<-sum_q[sum_q$SCHED_TM==16,]
iii<-t.test(Sum ~ trt, data = c)
d<-sum_q[sum_q$SCHED_TM==4,]
iv<-t.test(Sum ~ trt, data = d)
e<-sum_q[sum_q$SCHED_TM==6,]
v<-t.test(Sum ~ trt, data = e)
f<-sum_q[sum_q$SCHED_TM==8,]
vi<-t.test(Sum ~ trt, data = f)
g<-sum_q[sum_q$SCHED_TM==12,]
vi<-t.test(Sum ~ trt, data = g)
h<-sum_q[sum_q$SCHED_TM==12,]
vii<-t.test(Sum ~ trt, data = h)
cbind(i$p.value,ii$p.value,iii$p.value,iv$p.value,v$p.value,vi$p.value,vii$p.value)

##################Linear mixed effect model to choose covariates##########
rm(list=ls())
library(afex)
set.seed(0)
#####UVA####
sum_q<-read.csv(file = "~/Desktop/R/sum_1.csv",header = TRUE)
varlist<-names(sum_q)[7:16]
models<-lapply(varlist, function(x){
  fm3 <- summary(mixed(substitute(Sum ~ i + (SCHED_TM ||SID),list(i=as.name(x))), data=sum_q))})
#capture.output(models,file = "list2.txt")
####MVA###
fmf <- lmer(Sum ~ SEX+ WEC_PERF_STA + (SCHED_TM ||SID), data=sum_q)
# Backward elimination using terms with default alpha-levels:
step(fmf, direction = "backward", test = "F")
(step_res <- step(fmf))
final <- get_model(step_res)
anova(final)


#########Cox Regression for QoL######
QoL<-read.csv(file = "~/Desktop/R/Result/Nco20.csv",header = TRUE)
fit <- coxph(Surv(os, os_event) ~ Baseline+Slope+q, data = QoL)
summary(fit)


##############GWAS#######################
#UVA variable-> base QOL
rm(list=ls())
setwd("Z:/Thesis/Sum")
a<-read.csv(file="Z:/Thesis/Sum/sum.csv")
b<-a[,c(1:4,21)]
PR<-a[a$SCHED_TM==0,]
ia=c(4,7:9,11:15)
Table.uni=NULL
for (i in ia)
{linear <- lm(Sum~PR[,i], data = PR) 
Table.uni=rbind(Table.uni,
                cbind(index=i,
                      Gen=colnames(PR)[i],
                      Est=summary(linear)$coef[2,1],
                      P=summary(linear)$coef[2,4]))
}
write.csv(Table.uni,file="Linear.csv", row.names=F)

FullModel<-lm(Sum~age+ SEX+ WEC_HEIGHT+ WEC_WEIGHT+ WEC_PERF_STA, data = PR)

step(FullModel, direction = "backward", test = "F")
###SEX age ECOG
res1 <- assocRegression(genoData,
                        outcome="sum",
                        model.type="linear",
                        covar=c("sex", "age","ecog"))
b<-res1[order(res1$Wald.pval),]
c<-head(b,20)

pval=res1$Wald.pval
r<-gcontrol2(pval)
print(r$lambda)



###########
#UVA variable-> base QOL
a<-read.csv(file="Z:/Thesis/Sum/sum.csv")
ia=c(4,7:9,11:15)
Table.uni=NULL
for (i in ia)
{cox <- coxph(Surv(os, os_event) ~ a[,i], data =  a)
Table.uni=rbind(Table.uni,
                cbind(index=i,
                      Gen=colnames(a)[i],
                      Est=summary(cox)$coef[1,1],
                      P=summary(cox)$coef[1,5]))
}
write.csv(Table.uni,file="cox_slope.csv", row.names=F)

FullModel<-lm(slope~ trt+ age+ SEX+ WEC_MALIG +pcnum+ livmeta+ WEC_HEIGHT+ WEC_WEIGHT+ WEC_PERF_STA, data = a)
step(FullModel, direction = "backward", test = "F")


res1 <- assocCoxPH(genoData,
                   event="os_event", time.to.event="os",
                   covar=c("PC1","PC2","PC3","slope","sex","pcnum", "trt","age", "ecog"),ivar=c("slope"))


b<-res1[order(res1$GxE.pval),]
c<-head(b,20)

pval=res1$GxE.pval
r<-gcontrol2(pval)
print(r$lambda)
write.table(res1, file= "multivariable gen x slope.csv")




###########Gen->Baseline
base <- read.table("Z:/Thesis/plink/linearpca.base.glm.linear", header=F, stringsAsFactors=F)
# non-factor is important
base = base[complete.cases(base), ]
Add<-base[base$V7 == "ADD",]
b<-Add[order(Add$V12),]
c<-head(b,20)
# remove NA lines
Add$CHISQ <- qchisq(Add$V12,1,lower.tail=FALSE)
# qchisq(base$P,1,lower.tail=FALSE) can convert p-value (even < 5.5e-17) to chisq
# while qchisq(1-base$P,1) fails to convert small p-value (Inf in this case) 

lambda <- median(Add$CHISQ) / qchisq(0.5,1)

# qchisq(0.5,1) will give 0.4549364
# some people directly / 0.454 which may cause slightly different result
# chisq = z**2
# z <- qnorm(p-value/2)

#################### qqplot #####################
# there are various way to make qqplot, the quickest way is qqman package with less controls
qq(Add$V12)
# draw qqplot step by step
#m <- length(Add$V12)
#x <- rev(qchisq(1-(1:m)/(m), 1))
#y <- Add$CHISQ 
#qqplot(x, y, xlab="Expected chi-squared value", 
#ylab="Observed chi-squared value",        main = "qq-plot",
#sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))),
#cex.lab=1.2, cex.axis=1.2, font.lab=2,font.axis=2, lwd=1.3) 
# plot() will work also abline(0, 1, col='red') 

#Adjusted model baseline
adj_b <- read.table("Z:/Thesis/plink/adjustbase.base.glm.linear", header=F, stringsAsFactors=F)
adj_b = adj_b[complete.cases(adj_b), ]
adj_Add<-adj_b[adj_b$V7 == "ADD",]
b<-adj_Add[order(adj_Add$V12),]
c<-head(b,20)
adj_Add$CHISQ <- qchisq(adj_Add$V12,1,lower.tail=FALSE)
lambda <- median(adj_Add$CHISQ) / qchisq(0.5,1)
qq(adj_Add$V12)


###########Gen->slope
slope <- read.table("Z:/Thesis/plink/linearslope.slope.glm.linear", header=F, stringsAsFactors=F)
# non-factor is important
slope = slope[complete.cases(slope), ]
Add_s<-slope[slope$V7 == "ADD",]
b<-Add_s[order(Add_s$V12),]
c<-head(b,20)
Add_s$CHISQ <- qchisq(Add_s$V12,1,lower.tail=FALSE)
lambda <- median(Add_s$CHISQ) / qchisq(0.5,1)
qq(Add_s$V12)

#Adjusted model baseline
adj_s <- read.table("Z:/Thesis/plink/adjustlinearslope.slope.glm.linear", header=F, stringsAsFactors=F)
adj_s = adj_s[complete.cases(adj_s), ]
adj_Add_s<-adj_s[adj_s$V7 == "ADD",]
b<-adj_Add_s[order(adj_Add_s$V12),]
c<-head(b,20)
adj_Add_s$CHISQ <- qchisq(adj_Add_s$V12,1,lower.tail=FALSE)
lambda <- median(adj_Add_s$CHISQ) / qchisq(0.5,1)
qq(adj_Add_s$V12)



####Apendix###
##Mclust for 7 time points
#Agglomerative hierarchical clustering methods based on Gaussian probability models have recently shown promise in a variety of applications. 
#In this approach, a maximum-likelihood pair of clusters is chosen for merging at each stage. Unlike classical methods, model-based methods reduce to a recurrence relation only in the simplest case, 
#which corresponds to the classical sum of squares method. We show how the structure of the Gaussian model can be exploited to yield efficient algorithms for agglomerative hierarchical clustering.
#The optimal model according to BIC for EM initialized by hierarchical clustering for parameterized Gaussian mixture models.
rm(list=ls())
library(mclust)
sum_q<-read.csv(file = "~/Desktop/R/qol_sum.csv",header = TRUE)
a<-sum_q[sum_q$SCHED_TM==0,]
id_x<-a[,1]
X<-a[,3]
set.seed(0)
BIC<-mclustBIC(X)
summary(BIC)
mod1 <- Mclust(X, x = BIC)
summary(mod1, parameters = TRUE)
plot(mod1, what = "classification")
c<-mod1$classification
table(c)
i<-cbind(id_x,c)
zi<-as.data.frame(i)

a<-sum_q[sum_q$SCHED_TM==2,]
id_x<-a[,1]
X<-a[,3]
set.seed(0)
BIC<-mclustBIC(X)
summary(BIC)
mod1 <- Mclust(X, x = BIC)
summary(mod1, parameters = TRUE)
plot(mod1, what = "classification")
d<-mod1$classification
table(c)
ii<-cbind(id_x,d)
zii<-as.data.frame(ii)

a<-sum_q[sum_q$SCHED_TM==4,]
id_x<-a[,1]
X<-a[,3]
set.seed(0)
BIC<-mclustBIC(X)
summary(BIC)
mod1 <- Mclust(X, x = BIC)
summary(mod1, parameters = TRUE)
#plot(mod1, what = "classification")
e<-mod1$classification
iii<-cbind(id_x,e)
ziii<-as.data.frame(iii)

a<-sum_q[sum_q$SCHED_TM==6,]
id_x<-a[,1]
X<-a[,3]
set.seed(0)
BIC<-mclustBIC(X)
summary(BIC)
mod1 <- Mclust(X, x = BIC)
summary(mod1, parameters = TRUE)
#plot(mod1, what = "classification")
f<-mod1$classification
iv<-cbind(id_x,f)
ziv<-as.data.frame(iv)

a<-sum_q[sum_q$SCHED_TM==8,]
id_x<-a[,1]
X<-a[,3]
set.seed(0)
BIC<-mclustBIC(X)
summary(BIC)
mod1 <- Mclust(X, x = BIC)
summary(mod1, parameters = TRUE)
#plot(mod1, what = "classification")
g<-mod1$classification
v<-cbind(id_x,g)
zv<-as.data.frame(v)

a<-sum_q[sum_q$SCHED_TM==12,]
id_x<-a[,1]
X<-a[,3]
set.seed(0)
BIC<-mclustBIC(X)
summary(BIC)
mod1 <- Mclust(X, x = BIC)
summary(mod1, parameters = TRUE)
#plot(mod1, what = "classification")
h<-mod1$classification
vi<-cbind(id_x,h)
zvi<-as.data.frame(vi)

a<-sum_q[sum_q$SCHED_TM==16,]
id_x<-a[,1]
X<-a[,3]
set.seed(0)
BIC<-mclustBIC(X)
summary(BIC)
mod1 <- Mclust(X, x = BIC)
summary(mod1, parameters = TRUE)
#plot(mod1, what = "classification")
q<-mod1$classification
vi<-cbind(id_x,q)
zvii<-as.data.frame(vi)
#require(openxlsx)
#list_of_datasets <- list("PR" = zi, "W2" = zii, "W4" = ziii, "W6"=ziv, "W8" = zv, "W12" = zvi, "W16" = zvii)
#write.xlsx(list_of_datasets, file = "mclust_result.xlsx")

df<- read.xlsx(xlsxFile = "~/Desktop/R/mclust_result.xlsx", sheet = 1)
#count NAs per row
df$na_count <- apply(df, 1, function(x) sum(is.na(x)))
#delete rows with condition
#d<-df[!(df$na_count==6),]
#count of each value for each id
d$one<-rowSums(d[,c(2:8)] == 1, na.rm = T)
d$two<-rowSums(d[,c(2:8)] == 2, na.rm = T)
d$three<-rowSums(d[,c(2:8)] == 3, na.rm = T)
d$four<-rowSums(d[,c(2:8)] == 4, na.rm = T)
#get maximum value for each row
d[, "max"] <- apply(d[, 10:13], 1, max)
#get proportion
d$prop<-d$max/(7-d$na_count)
#
d$category <- cut(d$prop, 
                  breaks=c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9,1), 
                  labels=c("0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1"))

mytable<-table(d$category)
a<-data.frame(mytable)
ylim <- c(0, 1.1*max(a$Freq))
xx <- barplot(a$Freq, xaxt = 'n', width = 0.85, ylim = ylim,       
              main="Count for 6 time points n=678 ",
              xlab="range",
              ylab="Frequency")
## Add text at top of bars
text(x = xx, y = a$Freq, label = a$Freq, pos = 3, cex = 0.7, col = "blue")
## Add x-axis labels 
axis(1, at=xx, labels=a$Var1, tick=FALSE, las=2, line=-0.5, cex.axis=0.5)


#Appendix B#
####K-M plot for os by SNP####
rs12915196_0<-dt[dt$rs12915196_C==0,]
fit <- survfit(Surv(os, os_event) ~ s, data = rs12915196_0)
ggsurvplot(fit, data=dt, title = "Overall Survival for rs12915196_0",
           legend=c(0.9,0.9),
           xlab = "Time in months",
           legend.labs =c("Positive slope","Negative slope"),
           risk.table=TRUE, risk.table.y.text =FALSE,
           pval=TRUE,size =1)

rs12915196_1<-dt[dt$rs12915196_C==1,]
fit1 <- survfit(Surv(os, os_event) ~ s, data = rs12915196_1)
ggsurvplot(fit1, data=dt, title = "Overall Survival for rs12915196_1",
           legend=c(0.9,0.9),
           xlab = "Time in months",
           legend.labs =c("Positive slope","Negative slope"),
           risk.table=TRUE, risk.table.y.text =FALSE,
           pval=TRUE,size =1)

rs12915196_2<-dt[dt$rs12915196_C==2,]
fit2 <- survfit(Surv(os, os_event) ~ s, data = rs12915196_2)
ggsurvplot(fit2, data=dt, title = "Overall Survival for rs12915196_2",
           legend=c(0.9,0.9),
           xlab = "Time in months",
           legend.labs =c("Positive slope","Negative slope"),
           risk.table=TRUE, risk.table.y.text =FALSE,
           pval=TRUE,size =1)

###############?Appendix C#####
##Linear mixed effect model with quadratic term###
setwd("Z:/New/GMM_730")
a<-read.csv(file="sum.csv")
as.numeric(a$SCHED_TM)
#########group 1######
trt1<-a[a$trt==1,]
set.seed(0)
fm<-lmer(Sum ~ poly(SCHED_TM, 2, raw = TRUE) + gender + ecog + 
           (poly(SCHED_TM, 2, raw = TRUE) | SID),
         data=trt1)  
fm1<-lmer(Sum ~ poly(SCHED_TM, 2, raw = TRUE) + 
            (poly(SCHED_TM, 2, raw = TRUE) | SID),
          data=trt1)  
anova(fm,fm1)
trt1<-coef(fm)$SID
#write.csv(trt1, file = "Z:/New/trt1_quadratic.csv")

#########group 2######
trt2<-a[a$trt==2,]
set.seed(0)
fm<-lmer(Sum ~ poly(SCHED_TM, 2, raw = TRUE) + gender + ecog + 
           (poly(SCHED_TM, 2, raw = TRUE) | SID),
         data=trt2)  
fm1<-lmer(Sum ~ poly(SCHED_TM, 2, raw = TRUE) + 
            (poly(SCHED_TM, 2, raw = TRUE) | SID),
          data=trt2)  
anova(fm,fm1)
trt2<-coef(fm)$SID
#write.csv(trt2, file = "Z:/New/trt2_quadratic.csv")
