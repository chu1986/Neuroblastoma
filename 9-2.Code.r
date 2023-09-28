setwd("")
library(ComplexHeatmap)
library(pRRophetic)
library(SimDesign)
library(cowplot)
library(ggplot2)
library(stringr)
library(ggpubr)
library(dplyr)
library(limma)

# set color
jco <- c("#BDD5EA","#FFA5AB","#011627","#2874C5","#EABF00","#868686","#C6524A","#80A7DE")

# load data
load("./tcga.exp.rda")

RiskGroup <- read.table("./risk.txt",row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F,header = T,quote = "")
RiskGroup$riskgroup <- ifelse(RiskGroup$riskgroup == "HRisk",1,0)
comsam <- intersect(rownames(RiskGroup), colnames(tcga.exp));length(comsam)
expr <- tcga.exp[,comsam]
expr["riskgroup",] <- RiskGroup[comsam,"riskgroup"]
gene <- "riskgroup"

#################################
drug <- read.table("drug-new.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = F)

for (i in gene) {
  message(paste0("analysis of ",i," starts..."))
  subexpr <- as.numeric(expr[i,])
  names(subexpr) <- colnames(expr)
  
  # drug analysis
  predictedPtype <- predictedBoxdat <- list() 
  dat <- log2(expr + 1)
  hsam <- colnames(dat)[as.numeric(dat[i,]) == 1 ]
  lsam <- setdiff(colnames(dat),hsam)

  plotp <- list()
  for (d in drug$V1[1:6]) { # I choose six drugs here
    set.seed(20201013) 
    cat(paste0("-- drug of ",d," is calculating...\n")) 
    
    predictedPtype[[d]] <- quiet(pRRopheticPredict(testMatrix = as.matrix(dat[,c(hsam,lsam)]),
                                                   drug = d,
                                                   tissueType = "allSolidTumors",
                                                   selection = 1)) 
    
    predictedBoxdat[[d]] <- data.frame("est.ic50" = predictedPtype[[d]],
                                       "group" = rep(c("HRisk","LRisk"),c(length(hsam),length(lsam))), 
                                       row.names = names(predictedPtype[[d]])) 
    predictedBoxdat[[d]]$group <- factor(predictedBoxdat[[d]]$group,levels = c("HRisk","LRisk")) 
    
    p <- ggplot(data = predictedBoxdat[[d]],aes(x = group, y = est.ic50, fill = group))+
      scale_fill_manual(values = jco[1:3]) + 
      geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.8,color="black") +
      geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=0.8, alpha = 0.7)+
      geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
      theme_pubr()+
      ylab(bquote("Estimated IC"[50]~"of"~.(d))) +
      xlab(paste0("Group by riskscore")) +#,i
      rremove("legend.title")+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.2),
            axis.ticks = element_line(size=0.2,color="black"),
            axis.ticks.length=unit(0.2,"cm"),
            legend.position = "none")+
      font("xylab",size=15)+  
      font("xy",size=15)+ 
      font("xy.text", size = 15) +  
      font("legend.text",size = 15) + 
      stat_compare_means(method = "wilcox", # you could change to t.test if you want
                         label.x = 1.5,hjust = 0.5)
    plotp[[d]] <- p 
    ggsave(paste0("violin plot of IC50 for ", d, " between expression groups of ", i,".pdf"),width = 4,height = 4)
  }
  p2 <- plot_grid(plotlist = plotp, ncol = 3)
  ggsave(paste0("violin plot of IC50 between expression groups of ", i,".pdf"),width = 10,height = 8)
  
}
