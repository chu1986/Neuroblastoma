#!/usr/local/bin Rscript
setwd("")
library(reshape2)
library(celldex)
library(SingleR)
library(monocle)
library(ggplot2)
library(Seurat)
library(dplyr)

logFCfilter=0.585
adjPvalFilter=0.05
sample_color <- c("#6A2A0E","#BA3E1A","#EC7A00","#FCD46D","#FFF7C6")

###---------###
load("../GSE192906_RAW/pbmc.rda")
pbmc <- subset(x = pbmc, subset = Group == "NB")
save(pbmc,file = "pbmc_raw.rda")#load("pbmc_raw.rda")


##
pdf(file="01.featureViolin.pdf", width=8, height=5)
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf(file="01.featureCor.pdf",width=8,height=5)
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

pbmc=subset(x = pbmc, subset = nFeature_RNA > 100 & percent.mt < 15)
pbmc$orig.ident <- factor(pbmc$orig.ident)
save(pbmc,file = "pbmc.rda")#load("pbmc.rda")

pdf(file="01.featureViolin-f.pdf", width=8, height=5)
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf(file="01.featureCor-f.pdf",width=8,height=5)
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "disp", nfeatures = 5000)
#
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="01.featureVar.pdf",width=10,height=6)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()


###################################02.PCA###################################
pbmc=ScaleData(pbmc)
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))

##
pdf(file="02.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()

#
pdf(file="02.PCA.pdf",width=6.5,height=6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()

#
pdf(file="02.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()

#
pdf(file="02.pcaElbowPlot.pdf",width=10,height=8)
ElbowPlot(pbmc, ndims = 20)
dev.off()

#
pbmc <- JackStraw(object = pbmc, dims = 20, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf(file="02.pcaJackStraw.pdf",width=7,height=6)
JackStrawPlot(object = pbmc, dims = 1:20)
dev.off()


save(pbmc,file = "pbmc_Pca.rda")#load("pbmc_Pca.rda")
###################################03.TSNE###################################
##TSNE
pcSelect=17
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)
pbmc <- FindClusters(object = pbmc, resolution = 0.5)
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)
pdf(file="03.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)
dev.off()
write.table(pbmc$seurat_clusters,file="03.tsneCluster.txt",quote=F,sep="\t",col.names=F)

save(pbmc,file = "pbmc_Tsne.rda")#load("pbmc_Tsne.rda")

##
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(cbind(Symbol=rownames(sig.markers), sig.markers),file="03.clusterMarkers.txt",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#
pdf(file="03.tsneHeatmap.pdf",width=length(table(pbmc.markers$cluster)),height=nrow(top10)/5)
DoHeatmap(object = pbmc, features = top10$gene, hjust = 0, angle = 0) + NoLegend()
dev.off()

###################################04.SingleR##################################################
counts<-pbmc@assays$RNA@counts
clusters<-pbmc@meta.data$seurat_clusters
ref=celldex::BlueprintEncodeData()
labels=ref$label.main
save(ref,file = "ref.rda")#load("ref.rda")

singler=SingleR(test=counts,
                ref =ref,
                labels=labels,
                clusters = clusters)
clusterAnn=as.data.frame(singler)
clusterAnn=cbind(id=row.names(clusterAnn), clusterAnn)
clusterAnn=clusterAnn[,c("id", "labels")]
write.table(clusterAnn,file="04.clusterAnn.txt",quote=F,sep="\t", row.names=F)
singler2=SingleR(test=counts,
                 ref =ref,
                 labels=labels)
cellAnn=as.data.frame(singler2)
cellAnn=cbind(id=row.names(cellAnn), cellAnn)
cellAnn=cellAnn[,c("id", "labels")]
write.table(cellAnn, file="04.cellAnn.txt", quote=F, sep="\t", row.names=F)

#cluster
newLabels=singler$labels
names(newLabels)=levels(pbmc)
pbmc.plot2=RenameIdents(pbmc, newLabels)
pdf(file="04.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = pbmc.plot2, pt.size = 2, label = TRUE)
dev.off()
save(newLabels,file = "newLabels_t.rda")#load("newLabels_t.rda")
save(pbmc.plot2,file = "pbmc.plot2.rda")#load("pbmc.plot2.rda")

pbmc <- pbmc.plot2
##cluster
pbmc.markers=FindAllMarkers(object = pbmc,
                            only.pos = FALSE,
                            min.pct = 0.25,
                            logfc.threshold = logFCfilter)
sig.cellMarkers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(cbind(Symbol=rownames(sig.cellMarkers), sig.cellMarkers),file="04.cellMarkers.txt",sep="\t",row.names=F,quote=F)

###################################05.monocle R###################################
#
monocle.matrix=as.matrix(pbmc@assays$RNA@data)
monocle.sample=pbmc@meta.data
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.clusterAnn=clusterAnn
monocle.markers=sig.cellMarkers

monocle <- list()
monocle[["monocle.matrix"]] <- monocle.matrix
monocle[["monocle.sample"]] <- monocle.sample
monocle[["monocle.geneAnn"]] <- monocle.geneAnn
monocle[["monocle.clusterAnn"]] <- monocle.clusterAnn
monocle[["monocle.markers"]] <- monocle.markers
save(monocle,file = "monocle.rda")#load("monocle.rda")

#Seurat#monocle
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])

#
clusterAnn=as.character(monocle.clusterAnn[,2])
names(clusterAnn)=paste0("cluster",monocle.clusterAnn[,1])
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)

#
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, as.vector(sig.cellMarkers$gene))
#plot_ordering_genes(cds)
pdf(file="05.trajectory.pdf",width = 8,height = 8)
plot_ordering_genes(cds)
dev.off()
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree')
cds <- orderCells(cds)
#
pdf(file="05.trajectory.State.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "State")
dev.off()
#
pdf(file="05.trajectory.cellType.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "cell_type2")
dev.off()
#
save(cds,file = "cds.rda")#load("cds.rda")

quit();
