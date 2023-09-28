setwd("")
library(reshape2)
library(ggplot2)
library(ggpubr)
load("../data_TCGA/RNAseq/tcga.exp.rda")


geneList <- read.table("./lasso coefficient.txt",row.names = 1,sep="\t",header=T,check.names=F,quote = "")
geneList <- rownames(geneList)
geneCard <- read.table("GeneCards-SearchResults.csv",sep=",",header=T,check.names=F)

Expr <- tcga.exp[rowMeans(tcga.exp) > 0.3,];max(Expr)
Expr <- log2(Expr + 1);max(Expr)
Group <- data.frame(Acc=colnames(Expr),
                    Tissue=ifelse(substr(colnames(Expr),18,19)==11,"Control","NBL"))
geneCard <- geneCard[geneCard$`Gene Symbol` %in% rownames(Expr),]
geneCard <- geneCard[order(geneCard$`Relevance score`, decreasing = T),]
geneCard <- as.vector(geneCard$`Gene Symbol`[1:20])
save(geneCard,file = "geneCard.rda")#load("geneCard.rda")

Expr <- as.data.frame(t(Expr[unique(c(geneList,geneCard)),colnames(Expr) %in% Group$Acc]))
Expr$Acc <- rownames(Expr)
DatGroup <- dplyr::inner_join(Group, Expr, by="Acc")

#####-------------#####
library(Hmisc)
head(DatGroup)
rownames(DatGroup) <- DatGroup$Acc
DatGroup <- DatGroup[,setdiff(colnames(DatGroup),colnames(Group))]
head(DatGroup)

##-----##
res <- rcorr(as.matrix(DatGroup))
##result_1 <- CorMatrix(res$r, res$P)
res_1 <- melt(res$r)
colnames(res_1) <- c("row", "column", "cor")
res_2 <- melt(res$P)
colnames(res_2) <- c("row", "column", "p")
result_1 <- dplyr::inner_join(res_1,res_2,by=c("row", "column"))
head(result_1)
data <- result_1[result_1$row %in% geneList & result_1$column %in% geneCard,]
#data[is.na(data$cor),]$cor <- 0
data$pv <- ""
data[which(data$p < 0.05),]$pv <- "*"
data[which(data$p < 0.01),]$pv <- "**"
data[which(data$p < 0.001),]$pv <- "***"
#data[which(data$p < 0.0001),]$pv <- "****"


p <- ggplot(data,aes(x=row,y=column)) +
  geom_point(aes(colour = cor, size=pv)) +
  labs(x="",y="") +
  scale_colour_gradient2(low = "blue", high = "red", mid = "white",
                                midpoint = 0, space = "Lab",#, limit = c(-0.6, 0.6)
                                name="Pearson\nCorrelation") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
  theme(axis.text=element_text(size = 15)) +
  theme(axis.text.x=element_text(colour = "black",size = 15,angle=45,hjust=1)) +#
  theme(axis.text.y=element_text(colour = "black", vjust=0,size = 15)) +
  theme(axis.title =element_text(size = 20)) +
  theme(text = element_text(size = 15))
ggsave(p, filename = "DiseaseGene ~ cormap2.pdf", width = 10, height = 8)

##-----##
tmp <- DatGroup
for (i in geneList) {
  for (j in geneCard) {
    
    test <- cor.test(tmp[,i],tmp[,j],method = "pearson")
    #test$estimate;test$p.value
    data <- data.frame(moduleScore=tmp[,j],
                       PCAScore=tmp[,i])
    ggplot(data,aes(moduleScore, PCAScore)) +
      xlim(((min(tmp[,j])*1.05)-(max(tmp[,j])*0.05)),
           ((max(tmp[,j])*1.05)-(min(tmp[,j])*0.05))) +
      ylim(((min(tmp[,i])*1.05)-(max(tmp[,i])*0.05)),
           ((max(tmp[,i])*1.05)-(min(tmp[,i])*0.05))) +
      #xlim(round(min(tmp[,j]),2),round(max(tmp[,j]),2)) +
      #ylim(round(min(tmp[,i]),2),round(max(tmp[,i]),2)) +
      xlab(j) +
      ylab(i) +
      geom_point() +
      geom_smooth() +
      theme_bw() +
      annotate("text",x=(max(tmp[,j])+min(tmp[,j]))*0.5,y=((max(tmp[,i])*0.95)+(min(tmp[,i])*0.05)),
               label = paste0("cor= ",round(test$estimate,3),"\n","pvalue= ",format(test$p.value, digits = 2)))
    ggsave(paste0(i," ~ ",j,".pdf"),width = 4,height = 5)
    #,x=((max(tmp[,j])+min(tmp[,j]))*0.5),y=((max(tmp[,i])*0.9)+(min(tmp[,i])*0.1))
  }
}

data <- result_1[result_1$row %in% geneList & result_1$column %in% geneCard,]
tmp <- data[which(data$p < 0.05),]
tmp <- tmp[order(tmp$cor),]
head(tmp);tail(tmp)
