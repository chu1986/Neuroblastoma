setwd("")
library(RcisTarget.hg19.motifDBs.cisbpOnly.500bp)
library(RcisTarget)
library(visNetwork)
library(reshape2)
library(AUCell)
library(igraph)
library(doRNG)
library(doMC)
library(DT)

###########-------------################
gene <- read.table("./lasso coefficient.txt",row.names = 1,sep="\t",header=T,check.names=F,quote = "")
(geneList1 <- rownames(gene))
geneLists <- list(Module_Gene=geneList1)

## human:
data(motifAnnotations_hgnc)
data(hg19_500bpUpstream_motifRanking_cispbOnly)
motifRankings <- hg19_500bpUpstream_motifRanking_cispbOnly
motifRankings

##
motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                         motifAnnot=motifAnnotations_hgnc)
head(motifEnrichmentTable_wGenes)

##--plot--##
motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)

resultsSubset <- motifEnrichmentTable_wGenes_wLogo#[1:10,]

datatable(resultsSubset, #[,-c("enrichedGenes", "TF_lowConf"), with=FALSE]
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))

#####################################
motifs_AUC <- calcAUC(geneLists, motifRankings, nCores=1)
(motifs_AUC )

###########motif#annotation#############################

motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=3,
                                           motifAnnot=motifAnnotations_hgnc,
                                           highlightTFs=list(Key_Gene="AATF"))
head(motifEnrichmentTable[,-"TF_lowConf", with=FALSE])

motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                   rankings=motifRankings, 
                                                   geneSets=geneLists)

motifEnrichmentTable_wGenes[1:4,]

geneSetName <- names(geneLists)[1]
selectedMotifs <- c("cisbp__M2248","cisbp__M5576","cisbp__M4828")

pdf(file="motif enrichment best.pdf",width=8,height=8)
par(mfrow=c(2,2))
getSignificantGenes(geneLists[[geneSetName]], 
                    motifRankings,
                    signifRankingNames=selectedMotifs,
                    plotCurve=TRUE, maxRank=5000, genesFormat="none",
                    method="aprox")
dev.off()

