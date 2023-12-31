setwd("./plot_go")
library(plyr)
library(ggplot2)
library(grid)
library(stringr)
library(gridExtra)

folder=grep("^go_analysis",dir("../"),value=T)
files=grep("^gsea_report_for_",grep(".tsv",dir(paste0("../",folder,"/")),value=T),value=T)
data = lapply(paste0("../",folder,"/",files),read.delim)
names(data) = str_split_fixed(files,"_",5)[,4]
dataSet = ldply(data, data.frame)
dataSetH = dataSet[dataSet$.id == "h",]
dataSetL = dataSet[dataSet$.id == "l",]
if(nrow(dataSetH[which(dataSetH$NOM.p.val < 0.05),]) > 4){
  dataSetH <- dataSetH[which(dataSetH$NOM.p.val < 0.05),]
}else{dataSetH <- dataSetH}
if(nrow(dataSetL[which(dataSetL$NOM.p.val < 0.05),]) > 2){
  dataSetL <- dataSetL[which(dataSetL$NOM.p.val < 0.05),]
}else{dataSetL <- dataSetL}

files=c(as.character(dataSetH[order(dataSetH$NES, decreasing = T),]$NAME[1:4]),
        as.character(dataSetL[order(dataSetL$NES),]$NAME[1:2]))
#files=paste0("../",folder,"/",files,".tsv")
for (i in files) {
  file.copy(from = paste0("../",folder,"/",files,".tsv"),
            to =paste0("./",files,".tsv"))
}

files=grep(".tsv",dir(),value=T)
data = lapply(files,read.delim)
names(data) = files

dataSet = ldply(data, data.frame)
dataSet$pathway = gsub(".tsv","",dataSet$.id)

gseaCol=c("#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
pGsea=ggplot(dataSet,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,colour=pathway,group=pathway))+
  geom_line(size = 1.5) + scale_color_manual(values = gseaCol[1:nrow(dataSet)]) +   
  labs(x = "", y = "Enrichment Score", title = "") + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0),limits = c(min(dataSet$RUNNING.ES - 0.02), max(dataSet$RUNNING.ES + 0.02))) +   
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black")) + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + 
  geom_hline(yintercept = 0) +   theme(legend.position = c(0,0),legend.justification = c(0,0)) + #legendע?͵?λֵ
  guides(colour = guide_legend(title = NULL)) + theme(legend.background = element_blank()) + theme(legend.key = element_blank())+theme(legend.key.size=unit(0.5,'cm'))
pGene=ggplot(dataSet,aes(RANK.IN.GENE.LIST,pathway,colour=pathway))+geom_tile()+
  scale_color_manual(values = gseaCol[1:nrow(dataSet)]) + 
  labs(x = "high riskscore<----------->low riskscore", y = "", title = "") + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +  
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black"))+
  theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+ guides(color=FALSE)

gGsea = ggplot_gtable(ggplot_build(pGsea))
gGene = ggplot_gtable(ggplot_build(pGene))
maxWidth = grid::unit.pmax(gGsea$widths, gGene$widths)
gGsea$widths = as.list(maxWidth)
gGene$widths = as.list(maxWidth)
dev.off()

#
pdf('1.GSEA GO.pdf',
     width=7,
     height=5.5)
par(mar=c(5,5,2,5))
grid.arrange(arrangeGrob(gGsea,gGene,nrow=2,heights=c(.8,.3)))
dev.off()

