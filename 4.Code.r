setwd("")
library(patchwork)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(cowplot)
library(stringr)
library(ggpubr)
library(dplyr)


all_pval = read.table("out/pvalues.txt", header=T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names=F)
all_means = read.table("out/means.txt", header=T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names=F)
all_pval <- as.matrix(all_pval)
all_means <- as.matrix(all_means)
rownames(all_pval) <- all_pval[,"interacting_pair"]
rownames(all_means) <- all_means[,"interacting_pair"]

all_pval = all_pval[,-c(1:11)]
all_means = all_means[,-c(1:11)]

pval_row_0.05 <- apply(all_pval,1,
                       function(x) sum(x < 0.05))
pval_col_0.05 <- apply(all_pval,2,
                         function(x) sum(x < 0.05))
sel_pval_row <- names(pval_row_0.05[order(pval_row_0.05, decreasing = T)])[1:60]
sel_pval_col <- names(pval_col_0.05[order(pval_col_0.05, decreasing = T)])[1:30]

means_row <- apply(all_means,1,
                       function(x) sum(abs(as.numeric(x)) > 0.5))
means_col <- apply(all_means,2,
                       function(x) sum(abs(as.numeric(x)) > 0.5))

sel_means_row <- names(means_row[order(means_row, decreasing = T)])[1:60]
sel_means_col <- names(means_col[order(means_col, decreasing = T)])[1:30]

selected_rows <- intersect(sel_means_row, sel_pval_row)[!is.na(intersect(sel_means_row, sel_pval_row))];length(selected_rows)
write.table(data.frame(Symbol=selected_rows),file="in/easy_input_rows.txt",sep="\t",quote=F,row.names=F,col.names = F)
selected_columns <- intersect(sel_means_col, sel_pval_col)[!is.na(intersect(sel_means_col, sel_pval_col))];length(selected_columns)
write.table(data.frame(Symbol=selected_columns),file="in/easy_input_columns.txt",sep="\t",quote=F,row.names=F,col.names = F)

sel_pval <- all_pval[selected_rows, selected_columns]
outTab1 <- data.frame(row.names = row.names(sel_pval))
for (i in colnames(sel_pval)) {
  outTab1[,i] <- as.numeric(sel_pval[,i])
}
sel_means <- all_means[selected_rows, selected_columns]
#sel_means <- as.data.frame(sel_means)
outTab2 <- data.frame(row.names = row.names(sel_means))
for (i in colnames(sel_means)) {
  outTab2[,i] <- as.numeric(sel_means[,i])
}

df_names = expand.grid(selected_rows, selected_columns)
#########################

pval = unlist(outTab1)
pval[pval==0] = 0.0009
plot.data = cbind(df_names,pval)
pr = unlist(as.data.frame(outTab2))
pr[pr==0] = 1
plot.data = cbind(plot.data,log2(as.numeric(pr)))
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
head(plot.data)
##
my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

pdf("receptorLigand.pdf", width = 12, height = 10)
ggplot(plot.data, aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', 
                        colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(size=12, angle=90, hjust=1, vjust = 0.2),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

dev.off()

###########################
library(psych)
library(qgraph)
library(igraph)
library(tidyverse)

mynet <- read.delim("out/count_network.txt", check.names = FALSE)
table(mynet$count)
mynet %>% filter(count>0) -> mynet
head(mynet)
net<- graph_from_data_frame(mynet)

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
            "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
            "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493",
            "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
            "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

karate_groups <- cluster_optimal(net)
coords <- layout_in_circle(net, order =
                             order(membership(karate_groups)))

E(net)$width  <- E(net)$count/10
pdf(file="Cell net.pdf",width=12,height=10)
par(mar=c(.3,.3,.3,.3))
plot(net, edge.arrow.size=.1, 
     edge.curved=0,loop.angle2=0.1,edge.curved=0.5,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",#allcolour
     layout = coords,
     vertex.label.cex=1.5)
dev.off()

########################
data = read.table("out/interaction_count.txt",row.names = 1 , header=T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names=F)
data$Cell_Type <- rownames(data)
data <- data[order(data$all_sum, decreasing = T),]
data$Cell_Type <- factor(data$Cell_Type, levels = data$Cell_Type)

sample_color <- c("#FB040B","#F6A717","#BA06FA","#172D7A")
sample_color1 <- ifelse(data$all_sum > 350,ifelse(data$all_sum > 450,ifelse(data$all_sum > 550,sample_color[1],sample_color[2]),sample_color[3]),sample_color[4])

pB2 <- ggplot(data = data, aes(x = Cell_Type, y = all_sum, fill = Cell_Type, colour = Cell_Type)) +
  geom_bar(stat = "identity", width=0.8)+
  geom_hline(yintercept = c(350,450,550), 
             color="white",
             linetype = 2, 
             size = 0.3) + 
  scale_fill_manual(values=sample_color1) +
  scale_colour_manual(values=sample_color1) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="Cell Type",y="Counts")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, angle = 60, hjust = 1, vjust = 1 , colour = "black"))
pdf(file="Interaction Count.pdf",width=8.5,height=8)
print(pB2)
dev.off()
