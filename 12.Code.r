setwd("")

a5 <- read.table("input1.txt",header = T,row.names = 1,sep = "\t", quote = "",fill = T)      #tcga    
head(a5)
a5$died <- a5$fustat==1
str(a5)
a5$futime <- as.numeric(a5$futime)
a5$fustat <- as.numeric(a5$fustat)
a5$age <- as.numeric(a5$age)
a5$grade <- as.numeric(a5$grade)
a5$stage <- as.numeric(a5$stage)
a5$riskscore <- as.numeric(a5$riskscore)
data6 <- a5
head(data6)
library(rms)

dd<-datadist(data6)
options(datadist="dd")
options(na.action="na.delete")
summary(data6$futime)

coxpbc<-cph(formula = Surv(futime,died) ~ age + gender + grade + stage + riskscore,data=data6,x=T,y=T,surv = T,na.action=na.delete)

print(coxpbc)

surv<-Survival(coxpbc) 
surv3<-function(x) surv(1095,x)
surv4<-function(x) surv(1825,x)

x<-nomogram(coxpbc,fun = list(surv3,surv4),lp=T,
            funlabel = c('3-year survival Probability','5-year survival Probability'),
            maxscale = 100,fun.at = c(0.95,0.85,0.7,0.5,0.3,0.1))

pdf("nomogram_classical.pdf",width = 12,height = 10)
plot(x, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()

###5y#####
f5<-cph(formula = Surv(futime,died) ~ age + gender + grade + stage + riskscore,data=data6,x=T,y=T,surv = T,na.action=na.delete,time.inc = 1095) 
cal5<-calibrate(f5, cmethod="KM", method="boot",u=1095,m=65,B=1000)

##
f8<-cph(formula = Surv(futime,died) ~ age + gender + grade + stage + riskscore,data=data6,x=T,y=T,surv = T,na.action=na.delete,time.inc = 1825) 
cal8<-calibrate(f8, cmethod="KM", method="boot",u=1825,m=65,B=1000)


######curve####
pdf("calibration_compare.pdf",width = 8,height = 8)
plot(cal5,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", 
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal8,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal8[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", 
       legend = c("3-year","5-year"), 
       col =c("#2166AC","#B2182B"), 
       lwd = 2,
       cex = 1.2,
       bty = "n")
dev.off()

####-------------------------####
library(timeROC)

ROC.DSST<-timeROC(T=data6$futime,
                  delta=data6$fustat,
                  marker=data6$riskscore,
                  other_markers=as.matrix(data6[,c("age", "gender", "grade", "stage")]),
                  weighting="cox",cause=1,
                  times=c(365,1095,1825),ROC=TRUE)

cols <- c("#FD0207", "#0122D4", "#26D401")
pdf("nomogram_ROC.pdf",width = 6,height = 6)
plot(ROC.DSST,time=365,col=cols[1],title = "")
plot(ROC.DSST,time=1095,add=TRUE,col=cols[2])
plot(ROC.DSST,time=1825,add=TRUE,col=cols[3])
legend("bottomright", legend=c(paste0("1-year AUC:",round(ROC.DSST$AUC["t=365"], 4)),
                               paste0("3-year AUC:",round(ROC.DSST$AUC["t=1095"], 4)),
                               paste0("5-year AUC:",round(ROC.DSST$AUC["t=1825"], 4))),
       col=cols, lwd=2)
dev.off()

##-------------##
library(ggplot2)
library(caret)
library(ggDCA)
library(rmda)
library(rms)

rt=data6
set.seed(123456)
lrm1 <- lrm(fustat ~ age, rt)
lrm2 <- lrm(fustat ~ gender, rt)
lrm3 <- lrm(fustat ~ grade, rt)
lrm4 <- lrm(fustat ~ stage, rt)
lrm8 <- lrm(fustat ~ riskscore, rt)

set.seed(123456)
dca_lrm <- dca(lrm1, lrm2, lrm3, lrm4, lrm8,
               model.names = c("age", "gender", "grade", "stage", "riskscore"))
p <- ggplot(dca_lrm)
ggsave(p, filename = "DCA plot.pdf", width = 6, height = 5)
