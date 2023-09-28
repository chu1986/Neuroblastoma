setwd("")
library(ggstatsplot)


ann <- read.table("./risk.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,quote = "")
ann <- ann[substr(ann$samID,1,6) == "TARGET",]
ann <- data.frame(ID = ann$samID,
                  Group = ann$riskgroup,
                  stringsAsFactors = T)
head(ann)
TIDE.res <- read.csv("TIDE_output.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
TIDE.res$ID <- rownames(TIDE.res);dim(TIDE.res)
TIDE.res <- dplyr::inner_join(ann,TIDE.res,by="ID")
head(TIDE.res)
dim(TIDE.res)
##--Responder--##
data_p <- TIDE.res[,c("Group", "Responder")]
head(data_p)
colnames(data_p) <- c("Group", "status")

pdf("Responder.pdf",width = 8,height = 8)
ggstatsplot::ggbarstats(
  data = data_p,
  x = status,
  y = Group)
dev.off()
##--No benefits--##
data_p <- TIDE.res[,c("Group", "No benefits")]
head(data_p)
colnames(data_p) <- c("Group", "status")

pdf("No benefits.pdf",width = 8,height = 8)
ggstatsplot::ggbarstats(
  data = data_p,
  x = status,
  y = Group)
dev.off()

####------#Exclusion#------####
data_p <- TIDE.res[,c("Group", "Exclusion")]
head(data_p)
colnames(data_p) <- c("Group", "status")

data_p[,2][data_p[,2] > median(data_p$status)] <- c("HScore")
data_p[,2][data_p[,2] != c("HScore")] <- c("LScore")
table(data_p$status)

pdf("Exclusion.pdf",width = 8,height = 8)
ggstatsplot::ggbarstats(
  data = data_p,
  x = status,
  y = Group)
dev.off()
####-------#Dysfunction#------####
data_p <- TIDE.res[,c("Group", "Dysfunction")]
head(data_p)
colnames(data_p) <- c("Group", "status")

data_p[,2][data_p[,2] > median(data_p$status)] <- c("HScore")
data_p[,2][data_p[,2] != c("HScore")] <- c("LScore")
table(data_p$status)

pdf("Dysfunction.pdf",width = 8,height = 8)
ggstatsplot::ggbarstats(
  data = data_p,
  x = status,
  y = Group)
dev.off()
