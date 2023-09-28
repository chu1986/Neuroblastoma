setwd("")
library(glmnet)
library(survival)
library(pheatmap)
library(gplots)
library(survcomp)
library(survivalROC)
# set colors
jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")

#---------------#
# training data #
load("./TCGA.expTime.rda")
load("./GSE.expTime.rda")

tpms <- TCGA.expTime
tpms <- tpms[which(tpms$futime > 0),]
colnames(tpms) <- gsub("-","_",colnames(tpms))
tum.sam <- rownames(tpms)

geo1.tpms <- GSE.expTime$`GSE62564-GPL11154-os`
colnames(geo1.tpms) <- gsub("-","_",colnames(geo1.tpms))
geo1.tpms <- geo1.tpms[which(geo1.tpms$futime > 0),]

geo2.tpms <- GSE.expTime$`GSE85047-GPL5175-os`
colnames(geo2.tpms) <- gsub("-","_",colnames(geo2.tpms))
geo2.tpms <- geo2.tpms[which(geo2.tpms$futime > 0),]

comgene <- intersect(intersect(colnames(tpms),colnames(geo1.tpms)),colnames(geo2.tpms));length(comgene)

tpms <- tpms[,comgene]##448-2
geo1.tpms <- geo1.tpms[,comgene]
geo2.tpms <- geo2.tpms[,comgene]

max(tpms[,3:ncol(tpms)]);max(geo1.tpms[,3:ncol(geo1.tpms)]);max(geo2.tpms[,3:ncol(geo2.tpms)])
tpms <- cbind(tpms[,comgene[1:2]], log2(tpms[,comgene[3:length(comgene)]] + 1))
#geo1.tpms <- cbind(geo1.tpms[,comgene[1:2]], log2(geo1.tpms[,comgene[3:length(comgene)]] + 1))
#geo2.tpms <- cbind(geo2.tpms[,comgene[1:2]], log2(geo2.tpms[,comgene[3:length(comgene)]] + 1))
max(tpms[,3:ncol(tpms)]);max(geo1.tpms[,3:ncol(geo1.tpms)]);max(geo2.tpms[,3:ncol(geo2.tpms)])

# univariate cox regression
Coxoutput <- NULL
for(i in colnames(tpms[,3:ncol(tpms)])){
  cox <- coxph(Surv(futime, fustat) ~ as.numeric(scale(tpms[,i])), data = tpms)
  coxSummary = summary(cox)
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene=i,
                                           HR=as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z=as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue=as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower=as.numeric(coxSummary$conf.int[,3][1]),
                                           upper=as.numeric(coxSummary$conf.int[,4][1]),
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}
Coxoutput <- Coxoutput[order(Coxoutput$pvalue),]
write.table(Coxoutput,"Coxoutput.tcga.txt",sep = "\t",row.names = F,quote = F)
rownames(Coxoutput) <- Coxoutput$gene

# survival-related genes
surv.genes <- Coxoutput[which(Coxoutput$pvalue < 0.05),]
table(surv.genes$HR < 1)#134+20
write.table(surv.genes,"surv.genes.txt",sep = "\t",row.names = F,quote = F)

#-----------------#
# validation data #
geo1.tpms <- geo1.tpms[,c("fustat","futime",rownames(surv.genes))]
geo2.tpms <- geo2.tpms[,c("fustat","futime",rownames(surv.genes))]


##
risk <- NULL;seed = 5214; i=1

createFolds <- function(strat_id, k) {
  set.seed(seed)
  if(k > length(strat_id)) {
    k <- length(strat_id)
  }	
  perm <- sample(length(strat_id))
  strat_id <- factor(strat_id, levels=sample(unique(strat_id)))
  
  strat_order <- order(strat_id[perm])
  
  num_item_per_fold_ceil <- ceiling(length(strat_id) / k)
  
  fold_ids <- rep(seq_len(k), times= num_item_per_fold_ceil)
  fold_ids <- fold_ids[seq_along(strat_id)]
  
  folds <- split(perm[strat_order], fold_ids)
  names(folds) <- paste0("Fold", seq_len(k))	
  return(folds)
}

fold <- createFolds(tum.sam,3)

train_sam <- tum.sam[-fold[[i]]]
test_sam <- setdiff(tum.sam,train_sam)

set.seed(seed = seed)
tmp <- tpms[train_sam,c("futime","fustat",rownames(surv.genes))]
colnames(tmp) <- make.names(colnames(tmp))
#--------------------------------------------------------------------------#
#             you can change other model to select features                #
cvfit = cv.glmnet(scale(as.matrix(tmp[,-c(1:2)])),
                  Surv(tmp$futime,tmp$fustat),
                  family = "cox",
                  alpha = 1,
                  nfold = 10) #
myCoefs <- coef(cvfit, s="lambda.min");
fea <- rownames(coef(cvfit, s = 'lambda.min'))[coef(cvfit, s = 'lambda.min')[,1]!= 0]
if(is.element("(Intercept)", fea)) {
  lasso_fea <- fea[-1] # 去掉截距项并排序
  lasso_coef <- myCoefs@x; names(lasso_coef) <- lasso_fea
} else {
  lasso_fea <- fea
  lasso_coef <- myCoefs@x; names(lasso_coef) <- lasso_fea
}

lasso_coef.hr <- data.frame(gene = names(lasso_coef),
                            coef = lasso_coef,
                            hr = Coxoutput[names(lasso_coef),"HR"],
                            low.ci = Coxoutput[names(lasso_coef),"lower"],
                            upp.ci = Coxoutput[names(lasso_coef),"upper"],
                            stringsAsFactors = F)
lasso_coef.hr <- lasso_coef.hr[order(lasso_coef.hr$coef,decreasing = F),]
#--------------------------------------------------------------------------#
# risk score in training
tmp <- scale(as.matrix(tpms[train_sam,lasso_fea]))
risk.score <- apply(tmp,1,function(x) {x %*% myCoefs@x})
risk.score.train <- risk.score

tmp <- tpms[names(risk.score),1:2]
tmp$futime <- tmp$futime/30.5
tmp$risk.score <- as.numeric(risk.score)
tmp$RiskGroup <- ifelse(tmp$risk.score > median(risk.score) ,"HRisk","LRisk")
risk <- rbind.data.frame(risk,
                         data.frame(samID = train_sam,
                                    riskscore = tmp$risk.score,
                                    riskgroup = tmp$RiskGroup,
                                    cohort = "TCGA training",
                                    stringsAsFactors = F),
                         stringsAsFactors = F)

fitd <- survdiff(Surv(futime, fustat) ~ RiskGroup, data=tmp, na.action=na.exclude)
p1 <- 1-pchisq(fitd$chisq, length(fitd$n)-1)
fit <- survfit(Surv(futime, fustat)~ RiskGroup, data=tmp, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
pdf("KM for training dataset.pdf",width = 4.5,height = 4)
par(bty="n", mgp = c(1.9,.33,0), mar=c(4.1,4.1,2.1,2.1)+.1, las=1, tcl=-.25)
plot(fit, col = jco[2:1], lwd = 1.4, xlab="Time (Months)", ylab="Overall survival",mark.time = T)
par(xpd=TRUE)
legend(x=100, y=1.05, bty="n", "Risk", cex=1, text.font=2)
legend(x=100, y=0.97, bty="n", text.col = jco[2:1], c("High","Low"), cex=0.9)
text(x=0, y=0.05, paste0("P ",ifelse(p1 < 0.001,"< 0.001",paste0("= ",round(p1,3)))), cex=1, pos=4)
invisible(dev.off())

fit <- coxph(Surv(futime,fustat) ~ risk.score,data = tmp)
cindex.train <- concordance.index(predict(fit),surv.time = tmp$futime,surv.event = tmp$fustat,method = "noether")

train.roc5 <- survivalROC(Stime=tmp$futime,
                          status=tmp$fustat,
                          marker = tmp$risk.score,
                          predict.time =5*12, 
                          method="KM")

train.roc3 <- survivalROC(Stime=tmp$futime,
                          status=tmp$fustat,
                          marker = tmp$risk.score,
                          predict.time =3*12, 
                          method="KM")

train.roc1 <- survivalROC(Stime=tmp$futime,
                          status=tmp$fustat,
                          marker = tmp$risk.score,
                          predict.time =1*12, 
                          method="KM")

pdf(file = "survivalROC for training dataset.pdf",width = 4,height = 4)
par(bty="o", mgp = c(2,0.5,0), mar = c(3.1,3.1,2.1,2.1),tcl=-.25)
plot(train.roc5$FP, train.roc5$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=jco[1],
     xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)",
     lwd = 2,main = "Training dataset")
lines(x=c(0,1),y=c(0,1),lwd=1.5,lty=2,col="grey40")
text(0.5,0.18,paste0("5-year AUC: ", round(train.roc5$AUC,2)),cex = 1,col = jco[1],adj = 0)
lines(train.roc3$FP, train.roc3$TP, col=jco[2],lwd = 2)
text(0.5,0.12,paste0("3-year AUC: ", round(train.roc3$AUC,2)),cex = 1,col = jco[2],adj = 0)
lines(train.roc1$FP, train.roc1$TP, col=jco[4],lwd = 2)
text(0.5,0.06,paste0("1-year AUC: ",round(train.roc1$AUC,2)),cex = 1,col = jco[4],adj = 0)
#text(0.5,0,paste0("C-index: ", round(cindex.train$c.index,2)),adj = 0)
invisible(dev.off())

# risk score in internal testing
tmp <- scale(as.matrix(tpms[test_sam,lasso_fea]))
risk.score <- apply(tmp,1,function(x) {x %*% myCoefs@x})
risk.score.test <- risk.score

tmp <- tpms[names(risk.score),1:2]
tmp$futime <- tmp$futime/30.5
tmp$risk.score <- as.numeric(risk.score)
tmp$RiskGroup <- ifelse(tmp$risk.score > median(risk.score) ,"HRisk","LRisk")

risk <- rbind.data.frame(risk,
                         data.frame(samID = test_sam,
                                    riskscore = tmp$risk.score,
                                    riskgroup = tmp$RiskGroup,
                                    cohort = "testing dataset",
                                    stringsAsFactors = F),
                         stringsAsFactors = F)

fitd <- survdiff(Surv(futime, fustat) ~ RiskGroup, data=tmp, na.action=na.exclude)
p2 <- 1-pchisq(fitd$chisq, length(fitd$n)-1)
fit <- survfit(Surv(futime, fustat)~ RiskGroup, data=tmp, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
pdf("KM for testing dataset.pdf",width = 4.5,height = 4)
par(bty="n", mgp = c(1.9,.33,0), mar=c(4.1,4.1,2.1,2.1)+.1, las=1, tcl=-.25)
plot(fit, col = jco[2:1], lwd = 1.4, xlab="Time (Months)", ylab="Overall survival",mark.time = T)
par(xpd=TRUE)
legend(x=100, y=1.05, bty="n", "Risk", cex=1, text.font=2)
legend(x=100, y=0.97, bty="n", text.col = jco[2:1], c("High","Low"), cex=0.9)
text(x=0, y=0.05, paste0("P ",ifelse(p2 < 0.001,"< 0.001",paste0("= ",round(p2,3)))), cex=1, pos=4)
invisible(dev.off())

fit <- coxph(Surv(futime,fustat) ~ risk.score,data = tmp)
cindex.test <- concordance.index(predict(fit),surv.time = tmp$futime,surv.event = tmp$fustat,method = "noether")

test.roc5 <- survivalROC(Stime=tmp$futime,
                         status=tmp$fustat,
                         marker = tmp$risk.score,
                         predict.time =5*12, 
                         method="KM")

test.roc3 <- survivalROC(Stime=tmp$futime,
                         status=tmp$fustat,
                         marker = tmp$risk.score,
                         predict.time =3*12, 
                         method="KM")

test.roc1 <- survivalROC(Stime=tmp$futime,
                         status=tmp$fustat,
                         marker = tmp$risk.score,
                         predict.time =1*12, 
                         method="KM")

pdf(file = "survivalROC for internal testing dataset.pdf",width = 4,height = 4)
par(bty="o", mgp = c(2,0.5,0), mar = c(3.1,3.1,2.1,2.1),tcl=-.25)
plot(test.roc5$FP, test.roc5$TP, type="l", xlim=c(0,1), ylim=c(0,1.035),col=jco[1],
     xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)",
     lwd = 2,main = "Testing dataset")
lines(x=c(0,1),y=c(0,1),lwd=1.5,lty=2,col="grey40")
text(0.5,0.18,paste0("5-year AUC: ", round(test.roc5$AUC,2)),cex = 1,col = jco[1],adj = 0)
lines(test.roc3$FP, test.roc3$TP, col=jco[2],lwd = 2)
text(0.5,0.12,paste0("3-year AUC: ", round(test.roc3$AUC,2)),cex = 1,col = jco[2],adj = 0)
lines(test.roc1$FP, test.roc1$TP, col=jco[4],lwd = 2)
text(0.5,0.06,paste0("1-year AUC: ",round(test.roc1$AUC,2)),cex = 1,col = jco[4],adj = 0)
#text(0.5,0,paste0("C-index: ", round(cindex.test$c.index,2)),adj = 0)
invisible(dev.off())

# risk score in external validation 1
tmp <- scale(as.matrix(geo1.tpms[,lasso_fea]))
risk.score <- apply(tmp,1,function(x) {x %*% myCoefs@x})
risk.score.validate <- risk.score

tmp <- geo1.tpms[names(risk.score),1:2]
tmp$futime <- tmp$futime/30.5
tmp$risk.score <- as.numeric(risk.score)
tmp$RiskGroup <- ifelse(tmp$risk.score > median(risk.score) ,"HRisk","LRisk")

risk <- rbind.data.frame(risk,
                         data.frame(samID = rownames(geo1.tpms),
                                    riskscore = tmp$risk.score,
                                    riskgroup = tmp$RiskGroup,
                                    cohort = "Validating set1",
                                    stringsAsFactors = F),
                         stringsAsFactors = F)

fitd <- survdiff(Surv(futime, fustat) ~ RiskGroup, data=tmp, na.action=na.exclude)
p3 <- 1-pchisq(fitd$chisq, length(fitd$n)-1)
fit <- survfit(Surv(futime, fustat)~ RiskGroup, data=tmp, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
pdf("KM for validation set1.pdf",width = 4.5,height = 4)
par(bty="n", mgp = c(1.9,.33,0), mar=c(4.1,4.1,2.1,2.1)+.1, las=1, tcl=-.25)
plot(fit, col = jco[2:1], lwd = 1.4, xlab="Time (Months)", ylab="Overall survival",mark.time = T)
par(xpd=TRUE)
legend(x=150, y=1.05, bty="n", "Risk", cex=1, text.font=2)
legend(x=150, y=0.97, bty="n", text.col = jco[2:1], c("High","Low"), cex=0.9)
text(x=0, y=0.05, paste0("P ",ifelse(p3 < 0.001,"< 0.001",paste0("= ",round(p3,3)))), cex=1, pos=4)
invisible(dev.off())

fit <- coxph(Surv(futime,fustat) ~ risk.score,data = tmp)
cindex.validate <- concordance.index(predict(fit),surv.time = tmp$futime,surv.event = tmp$fustat,method = "noether")

validate.roc5 <- survivalROC(Stime=tmp$futime,
                             status=tmp$fustat,
                             marker = tmp$risk.score,
                             predict.time =5*12, 
                             method="KM")

validate.roc3 <- survivalROC(Stime=tmp$futime,
                             status=tmp$fustat,
                             marker = tmp$risk.score,
                             predict.time =3*12, 
                             method="KM")

validate.roc1 <- survivalROC(Stime=tmp$futime,
                             status=tmp$fustat,
                             marker = tmp$risk.score,
                             predict.time =1*12, 
                             method="KM")

pdf(file = "survivalROC for external validation dataset1.pdf",width = 4,height = 4)
par(bty="o", mgp = c(2,0.5,0), mar = c(3.1,3.1,2.1,2.1),tcl=-.25)
plot(validate.roc5$FP, validate.roc5$TP, type="l", xlim=c(0,1), ylim=c(0,1.015),col=jco[1],
     xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)",
     lwd = 2,main = "Validation dataset1")
lines(x=c(0,1),y=c(0,1),lwd=1.5,lty=2,col="grey40")
text(0.5,0.18,paste0("5-year AUC: ", round(validate.roc5$AUC,2)),cex = 1,col = jco[1],adj = 0)
lines(validate.roc3$FP, validate.roc3$TP, col=jco[2],lwd = 2)
text(0.5,0.12,paste0("3-year AUC: ", round(validate.roc3$AUC,2)),cex = 1,col = jco[2],adj = 0)
lines(validate.roc1$FP, validate.roc1$TP, col=jco[4],lwd = 2)
text(0.5,0.06,paste0("1-year AUC: ",round(validate.roc1$AUC,2)),cex = 1,col = jco[4],adj = 0)
#text(0.5,0,paste0("C-index: ", round(cindex.validate$c.index,2)),adj = 0)
invisible(dev.off())

# risk score in external validation 2
tmp <- scale(as.matrix(geo2.tpms[,lasso_fea]))
risk.score <- apply(tmp,1,function(x) {x %*% myCoefs@x})
risk.score.validate2 <- risk.score

tmp <- geo2.tpms[names(risk.score),1:2]
tmp$futime <- tmp$futime/30.5
tmp$risk.score <- as.numeric(risk.score)
tmp$RiskGroup <- ifelse(tmp$risk.score > median(risk.score) ,"HRisk","LRisk")

risk <- rbind.data.frame(risk,
                         data.frame(samID = rownames(geo2.tpms),
                                    riskscore = tmp$risk.score,
                                    riskgroup = tmp$RiskGroup,
                                    cohort = "Validating set2",
                                    stringsAsFactors = F),
                         stringsAsFactors = F)

fitd <- survdiff(Surv(futime, fustat) ~ RiskGroup, data=tmp, na.action=na.exclude)
p4 <- 1-pchisq(fitd$chisq, length(fitd$n)-1)
fit <- survfit(Surv(futime, fustat)~ RiskGroup, data=tmp, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)
pdf("KM for validation set2.pdf",width = 4.5,height = 4)
par(bty="n", mgp = c(1.9,.33,0), mar=c(4.1,4.1,2.1,2.1)+.1, las=1, tcl=-.25)
plot(fit, col = jco[2:1], lwd = 1.4, xlab="Time (Months)", ylab="Overall survival",mark.time = T)
par(xpd=TRUE)
legend(x=150, y=1.10, bty="n", "Risk", cex=1, text.font=2)
legend(x=150, y=1.02, bty="n", text.col = jco[2:1], c("High","Low"), cex=0.9)
text(x=0, y=0.05, paste0("P ",ifelse(p4 < 0.001,"< 0.001",paste0("= ",round(p4,3)))), cex=1, pos=4)
invisible(dev.off())

fit <- coxph(Surv(futime,fustat) ~ risk.score,data = tmp)
cindex.validate2 <- concordance.index(predict(fit),surv.time = tmp$futime,surv.event = tmp$fustat,method = "noether")

validate2.roc5 <- survivalROC(Stime=tmp$futime,
                              status=tmp$fustat,
                              marker = tmp$risk.score,
                              predict.time =5*12,
                              method="KM")

validate2.roc3 <- survivalROC(Stime=tmp$futime,
                              status=tmp$fustat,
                              marker = tmp$risk.score,
                              predict.time =3*12,
                              method="KM")

validate2.roc1 <- survivalROC(Stime=tmp$futime,
                              status=tmp$fustat,
                              marker = tmp$risk.score,
                              predict.time =1*12,
                              method="KM")

pdf(file = "survivalROC for external validation dataset2.pdf",width = 4,height = 4)
par(bty="o", mgp = c(2,0.5,0), mar = c(3.1,3.1,2.1,2.1),tcl=-.25)
plot(validate2.roc5$FP, validate2.roc5$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=jco[1],
     xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)",
     lwd = 2,main = "Validation dataset2")
lines(x=c(0,1),y=c(0,1),lwd=1.5,lty=2,col="grey40")
text(0.5,0.18,paste0("5-year AUC: ", round(validate2.roc5$AUC,2)),cex = 1,col = jco[1],adj = 0)
lines(validate2.roc3$FP, validate2.roc3$TP, col=jco[2],lwd = 2)
text(0.5,0.12,paste0("3-year AUC: ", round(validate2.roc3$AUC,2)),cex = 1,col = jco[2],adj = 0)
lines(validate2.roc1$FP, validate2.roc1$TP, col=jco[4],lwd = 2)
text(0.5,0.06,paste0("1-year AUC: ",round(validate2.roc1$AUC,2)),cex = 1,col = jco[4],adj = 0)
#text(0.5,0,paste0("C-index: ", round(cindex.validate2$c.index,2)),adj = 0)
invisible(dev.off())

##
write.table(risk,"risk.txt",sep = "\t",row.names = F,quote = F)

# figures
# 1. coef_hr bar plot
darkred   <- "#F2042C"
darkblue   <- "#21498D"

cutoff <- 0
lasso_coef.hr$group <- as.character(cut(lasso_coef.hr$coef, breaks = c(-Inf, cutoff, Inf),labels = c("#EABF00","#21498D")))
pdf("lasso_coef_hr.pdf",width = 5,height = 2.5)
par(bty="n", mgp = c(1.7,.33,0),mar=c(2.5,2.7,1,1)+.1, las=1, tcl=-.25,xpd = T)
a <- barplot(lasso_coef.hr$coef,col = lasso_coef.hr$group,border = NA,
             horiz = T,xlim = c(-1,1),add=F,xaxt = "n")
axis(side = 1, at = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1),
     labels = c("-1","-.8","-.6","-.4","-.2","0",".2",".4",".6",".8","1"))
points(y = a[,1], x = log2(lasso_coef.hr$hr),
       col = lasso_coef.hr$group,
       pch = 19,cex = 1.5)
for (i in 1:nrow(lasso_coef.hr)) {
  text(y = a[,1][i],x = lasso_coef.hr$coef[i],labels = lasso_coef.hr$gene[i],adj = ifelse(lasso_coef.hr$coef[i]>0,0,1))
}
points(0.6,2,pch = 15, cex = 1.5)
points(0.6,1,pch = 19, cex = 1.5)
text(0.6,2,"Coefficient",pos = 4)
text(0.6,1,"log2(HR)",pos = 4)
invisible(dev.off())
write.table(lasso_coef.hr[,1:5], "lasso coefficient.txt",sep = "\t",row.names = F,col.names = T,quote = F)

# 2. lasso details
pdf("lasso.pdf",width = 4.5,height = 4)
par(bty="o", mgp = c(1.9,.33,0), mar=c(4.1,4.1,2.1,2.1)+.1, las=1, tcl=-.25,xpd = F)
plot(cvfit$glmnet.fit, "lambda", label=F)
abline(v=log(cvfit$lambda.min),lwd=2,col="grey60",lty=4)
invisible(dev.off())

pdf("cvfit.pdf",width = 4.5,height = 4)
par(bty="o", mgp = c(1.9,.33,0), mar=c(4.1,4.1,2.1,2.1)+.1, las=1, tcl=-.25)
plot(cvfit)
abline(h=min(cvfit$cvm),lwd=2,col="black",lty=4)
points(log(cvfit$lambda.min),min(cvfit$cvm),pch=18,cex=2,col="black")
points(log(cvfit$lambda.min),min(cvfit$cvm),pch=18,cex=1.5,col="#008B8A")
invisible(dev.off())
