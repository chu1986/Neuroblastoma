#!/usr/local/bin Rscript
setwd("")
library(stringr)
library(Seurat)

(folders=list.files('./',pattern='^GSM576'))
sceList = lapply(folders,function(folder){ 
  proName=paste(unlist(strsplit(folder, "\\-|\\_"))[1:2], collapse = "-")
  rt <- read.table(folder,row.names = 1,sep="\t",header=T,check.names=F)
  CreateSeuratObject(counts = rt, 
                     project = proName, names.delim = "-")})

pbName <- gsub("\\-|\\_","-",gsub("\\_UMI_COUNTS_RAW.txt$","",folders))
pbmc <- merge(sceList[[1]], 
                 y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],sceList[[6]],
                       sceList[[7]],sceList[[8]],sceList[[9]],sceList[[10]]), 
                 add.cell.ids = pbName, 
                 project = "Seurat")

##
sample <- str_split_fixed(names(pbmc$orig.ident),"\\-|\\_",2)[,1]
Group <- gsub("1|2|3|4|5", "", str_split_fixed(names(pbmc$orig.ident),"\\-|\\_",3)[,2])

pbmc[["sample"]] <- data.frame(row.names = rownames(pbmc[["orig.ident"]]),sample)
pbmc[["Group"]] <- data.frame(row.names = rownames(pbmc[["orig.ident"]]),Group)
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
save(pbmc,file = "pbmc.rda")#load("pbmc.rda")

quit();

