###-----------------------####
setwd("")
library(reshape2)
library(ggplot2)
library(ggpubr)


Data <- read.table("CIBERSORT-Results.txt",row.names = 1,sep="\t",header=T,check.names=F,quote = "")
head(Data)
Data <- Data[,setdiff(colnames(Data), c("P-value", "Correlation", "RMSE"))]
Data$Acc <- rownames(Data)
dim(Data)
head(Data)

RiskGroup <- read.table("./risk.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,quote = "")
colnames(RiskGroup)[1] <- "Acc"
data_p <- dplyr::inner_join(Data, RiskGroup, by = c("Acc"))
data_p <- melt(data_p, id.vars = colnames(RiskGroup))
head(data_p)

p <- ggplot(data_p,
            aes(x=variable, y=value, 
                fill = riskgroup,
                color = riskgroup)) + 
  geom_boxplot(notch = F, alpha = 0.95, 
               outlier.shape = 16,
               outlier.size = 0.65) +
  xlab("") + ylab("Fraction") +
  scale_fill_manual(values= c("#D5EBFB","#FBEEB7","#B4FBCD","#F5B3FC")) +
  scale_color_manual(values= c("#0073C2","#EFC000","#00C244","#C501D7")) +
  ggtitle("") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "top") +
  stat_compare_means(method ="wilcox.test",hide.ns = TRUE,label = "p.signif")
ggsave(p, filename = "CIBERSORT ~ RiskGroup.pdf", width = 12, height = 6)

##
library(Hmisc)

##-----##
rownames(Data) <- Data$Acc
tmp <- Data[,setdiff(colnames(Data), colnames(RiskGroup))]
res <- rcorr(as.matrix(tmp))
##result_1 <- CorMatrix(res$r, res$P)
res_1 <- melt(res$r)
colnames(res_1) <- c("row", "column", "cor")
res_2 <- melt(res$P)
colnames(res_2) <- c("row", "column", "p")
result_1 <- dplyr::inner_join(res_1,res_2,by=c("row", "column"))
head(result_1)
data <- result_1
#data[is.na(data$cor),]$cor <- 0
data$pv <- ""
data[which(data$p < 0.05),]$pv <- "*"
data[which(data$p < 0.01),]$pv <- "**"
data[which(data$p < 0.001),]$pv <- "***"
#data[which(data$p < 0.0001),]$pv <- "****"

p <- ggplot(data,aes(x=row,y=column)) +
  geom_tile(aes(fill = cor))+
  geom_text(aes(label=pv),col ="black",size = 5) +
  labs(x="",y="") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, space = "Lab",#, limit = c(-1, 1)
                       name="Pearson\nCorrelation") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
  theme(axis.text=element_text(size = 15)) +
  theme(axis.text.x=element_text(colour = "black",angle=45,hjust=1,size = 20)) +
  theme(axis.text.y=element_text(colour = "black", vjust=0,size = 20)) +
  theme(axis.title =element_text(size = 25)) +
  theme(text = element_text(size = 20))
ggsave(p, filename = "CIBERSORT ~ cormap-R.pdf", width = 15, height = 12)

###########-------------------------------------------########
library(dplyr)

Expr <- read.table("./risk.txt",row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F,header = T,quote = "")
Expr <- Expr[,"riskscore",drop=FALSE]
ciber <- read.table("CIBERSORT-Results.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# extract common samples
comsam <- intersect(rownames(Expr),rownames(ciber))
expr <- as.data.frame(t(Expr[comsam,,drop=FALSE]))
ciber <- ciber[comsam,setdiff(colnames(ciber), c("P-value", "Correlation", "RMSE"))]

for (i in rownames(expr)) {
  message(paste0("analysis of ",i," starts..."))
  subexpr <- as.numeric(expr[i,])
  names(subexpr) <- colnames(expr)
  
  # immune correlation
  dat <- as.numeric(expr[i,]); names(dat) <- colnames(expr)
  comsam <- intersect(names(dat), rownames(ciber))
  tmp1 <- dat[comsam]
  tmp2 <- ciber[comsam,]
  
  var <- setdiff(colnames(ciber),"CancerType")
  data <- data.frame(var)
  for (j in 1:length(var)){
    test <- cor.test(as.numeric(tmp2[,j]),tmp1,method = "pearson") # you could change to pearson if you want
    data[j,2] <- test$estimate                                            
    data[j,3] <- test$p.value
  }
  names(data) <- c("symbol","correlation","pvalue")
  data <- as.data.frame(na.omit(data))
  data <- data[data$pvalue < 0.05,]
  if(nrow(data) > 0){
    data %>% 
      ggplot(aes(correlation,forcats::fct_reorder(symbol,correlation))) +
      geom_segment(aes(xend=0,yend=symbol)) +
      geom_point(aes(col=pvalue,size=abs(correlation))) +
      scale_colour_gradientn(colours=c("#7fc97f","#984ea3")) +
      scale_size_continuous(range =c(2,8))  +
      theme_minimal() +
      ylab(NULL)
    ggsave(paste0("correlation/correlation between cibersort and ", i,".pdf"),width = 8,height = 6)
  }
}

######################-------------------------------########
library(plyr)
load("./tcga.exp.rda")

RiskGroup <- read.table("./risk.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,quote = "")
colnames(RiskGroup)[1] <- "Acc"
Expr <- as.data.frame(t(tcga.exp[rowSums(tcga.exp) > 0,]))
Expr <- log2(Expr + 1)
Expr$Acc <- rownames(Expr)
head(Expr[,1:3])
dim(Expr)
DatGroup <- dplyr::inner_join(RiskGroup, Expr, by="Acc")
head(DatGroup[,1:3])
dim(DatGroup)
data_p <- melt(DatGroup, id.vars = colnames(RiskGroup))
head(data_p[,1:3])
dim(data_p)

c <- read.table("Immunomodulator_and_chemokines.txt",header = F,sep = "\t", quote = "",fill = T,check.names=F)
colnames(c) <- c("id","type","ID")
head(c)
dim(c)
#
for (i in unique(c$type)) {
  
  ##-----##
  plotgene <- as.character(c[c$type == i,]$ID)
  tmp <- data_p[data_p$variable %in% plotgene,]
  a <- intersect(data_p$variable, plotgene)
  if(nrow(tmp) > length(unique(data_p$Acc))){
    
    p <- ggplot(tmp,aes(x=variable,y=value,
                        fill = riskgroup,
                        color = riskgroup)) + 
      geom_boxplot(notch = F, alpha = 0.95, 
                   outlier.shape = 16,
                   outlier.size = 0.65) +
      scale_fill_manual(values= c("#D5EBFB","#FBEEB7","#B4FBCD","#F5B3FC")) +
      scale_color_manual(values= c("#0073C2","#EFC000","#00C244","#C501D7")) +
      ggtitle("") + 
      labs(x=i,y="")
    p <- p + theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      theme(axis.text=element_text(size = 15)) +
      theme(axis.text.x=element_text(colour = "black",angle=45,hjust=1,size = 15)) +
      theme(axis.text.y=element_text(colour = "black", vjust=0,size = 15)) +
      theme(axis.title =element_text(size = 20)) +
      theme(text = element_text(size = 15)) +
      stat_compare_means(aes(group=riskgroup),method ="wilcox.test",hide.ns = TRUE,label = "p.signif")
    ggsave(p, filename = paste0(i,".pdf"), width = ((length(a)/3)+5), height = 6)
  }
}

###########-------------------------------------------########
setwd("E:/case_scRNA_NBL_lasso/ciber")
library(reshape2)
library(ggplot2)

ciber <- read.table("CIBERSORT-Results.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ciber <- ciber[,setdiff(colnames(ciber), c("P-value", "Correlation", "RMSE"))]
ciber$Acc <- rownames(ciber)
RiskGroup <- read.table("./risk.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,quote = "")
colnames(RiskGroup)[1] <- "Acc"
ciber <- dplyr::inner_join(ciber, RiskGroup, by="Acc")
ciber <- ciber[order(ciber$riskgroup, ciber$Acc, decreasing = T),]
data_p <- melt(ciber, id.vars = colnames(RiskGroup))
head(data_p)
data_p$Acc <- factor(data_p$Acc, levels = ciber$Acc)

p <- ggplot(data_p, aes(Acc, value, fill=variable)) +
  geom_bar(stat="identity", position = "fill", width = 0.5) +
  geom_col(position = 'stack', width = 0.6) +
  guides(fill=guide_legend(title = NULL)) +
  ylab("Relative Percent") + xlab("") +
  theme_bw() +
  theme(axis.ticks.length=unit(0.5,'cm')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
  theme(axis.text=element_text(size = 15)) +
  theme(axis.text.x=element_text(colour = "black",angle=45,hjust=1,size = 15)) +
  theme(axis.text.y=element_text(colour = "black", vjust=0,size = 15)) +
  theme(axis.title =element_text(size = 20)) +
  theme(text = element_text(size = 15)) +
  scale_y_continuous(expand=c(0,0))
ggsave(p, filename = "Immune infiltration-R.pdf", width = 30, height = 5)

legendcol <- c(rep("#1CFA04",length(unique(data_p[data_p$riskgroup == "LRisk",]$Acc))),
               rep("#C705FF",length(unique(data_p[data_p$riskgroup != "LRisk",]$Acc))))
p <- ggplot(data_p, aes(Acc, value, fill=variable)) +
  geom_bar(stat="identity", position = "fill", width = 0.5) +
  geom_col(position = 'stack', width = 0.6) +
  guides(fill=guide_legend(title = NULL)) +
  ylab("Relative Percent") + xlab("") +
  theme_bw() +
  theme(axis.ticks.length=unit(0.5,'cm')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
  theme(axis.text=element_text(size = 15)) +
  theme(axis.text.x=element_text(colour = legendcol,angle=45,hjust=1,size = 15)) +
  theme(axis.text.y=element_text(colour = "black", vjust=0,size = 15)) +
  theme(axis.title =element_text(size = 20)) +
  theme(text = element_text(size = 15)) +
  scale_y_continuous(expand=c(0,0))
ggsave(p, filename = "Immune infiltration-R__.pdf", width = 30, height = 6)
