rm(list = ls())

load(file = "../data/HCC_ssgsea_cell.Rdata")

im_timer = as.data.frame(t(cibersort))

im_timer[1:4,1:4]
im_timer$sample = rownames(im_timer)
load(file = "../data/HCC_metabolic.Rdata")
metabolic[1:3,1:3]
metabolic = as.data.frame(t(metabolic))
metabolic$sample = rownames(metabolic)

metabolic[1:3,1:3]


tmp = merge(im_timer, metabolic, by = "sample")

library(dplyr)
library(data.table)
library(tidyr)
library(tibble)
library(GSVA)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

rownames(tmp) = tmp$sample
tmp$sample = NULL
results = ConsensusClusterPlus(as.matrix(t(tmp)),
                               maxK=6,
                               reps=1000,
                               pItem=0.8,
                               pFeature=1,
                               #tmyPal = c('navy','darkred'),
                               title='Immune_ConsensusCluster/',
                               clusterAlg="km",
                               distance="spearman",
                               innerLinkage = "ward.D2", # 内部链接函数，可修改
                               finalLinkage = "ward.D2", # 最终链接函数，可修改
                               seed="123456",
                               plot = "pdf")
icl <- calcICL(results,title = 'Immune_ConsensusCluster/',plot = 'pdf')
load(file = "../data/HCC_result.Rdata")

Kvec = 2:6
x1 = 0.1; x2 = 0.9 
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") 
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}
optK = Kvec[which.min(PAC)]
optK


PAC <- as.data.frame(PAC)
PAC$K <- 2:6
library(ggplot2)

pdf(file = "./figtures/HCC_PAC_curve.pdf")
ggplot(PAC,aes(factor(K),PAC,group=1))+
  geom_line()+
  theme_bw(base_rect_size = 1.5)+
  geom_point(size=4,shape=21,color='darkred',fill='red')+
  ggtitle('Proportion of ambiguous clustering')+
  xlab('Cluster number K')+ylab(NULL)+
  theme(axis.text = element_text(size=12),
        plot.title = element_text(hjust=0.5),
        axis.title = element_text(size=13))
dev.off()

clusterNum=2
cluster=results[[clusterNum]][["consensusClass"]]

sub <- data.frame(Sample=names(cluster),Cluster=cluster)
sub$Cluster <- paste0('C',sub$Cluster)
table(sub$Cluster)

head(sub)

my <- results[[2]][["ml"]]
library(pheatmap)
rownames(my) <- sub$Sample
colnames(my) <- sub$Sample
pheatmap(1-my,show_colnames = F,show_rownames = F,
         cluster_cols = T,
         treeheight_row = 20,treeheight_col = 20,
         clustering_method = 'complete',
         color = colorRampPalette(c("white","#C75D30"))(50),
         annotation_names_row = F,annotation_names_col = F,
         annotation_row = data.frame(Cluster=sub$Cluster,row.names = sub$Sample),
         annotation_col = data.frame(Cluster=sub$Cluster,row.names = sub$Sample),
         annotation_colors = list(Cluster=c('C2'='#B5739D','C1'='#4E8279')))



rm(list = ls())



load(file = "../data/HCC_ssgsea_cell.Rdata")

im_timer = as.data.frame(t(cibersort))

im_timer[1:4,1:4]
im_timer$sample = rownames(im_timer)
load(file = "../data/HCC_metabolic.Rdata")
metabolic[1:3,1:3]
metabolic = as.data.frame(t(metabolic))
metabolic$sample = rownames(metabolic)


metabolic[1:3,1:3]
tmp = merge(im_timer, metabolic, by = "sample")
rownames(tmp) = tmp$sample
tmp$sample = NULL
library(NbClust)
library(NbClust)
set.seed(1234)
nc <- NbClust(scale(tmp), min.nc = 2, max.nc = 6, method = "kmeans")
table(nc$Best.nc[1,])
pdf(file = "../revised/HCC_Nbclust.pdf")
barplot(table(nc$Best.n[1,]),xlab="Number of Clusters",
        ylab="Number of Criteria",main="Number of Clusters Chosen by 26 Criteria")
dev.off()
sub = as.data.frame(nc$Best.partition)
sub$sample = rownames(sub)
colnames(sub)[1] = "NbClust"
load(file = "../data/HCC_pdata.Rdata")
dim(pdata)
colnames(pdata)
pdata = pdata[,c('sample',"subtype")]
table(pdata$subtype)
pdata = merge(pdata, sub, by = "sample")
table(pdata$subtype,pdata$NbClust)
library(cluster)
set.seed(1234)
fit.pam <- pam(scale(tmp), k=2, stand=TRUE)
pam = as.data.frame(fit.pam$clustering)
head(pam)
pam$sample = rownames(pam)
colnames(pam)[1] = "pam"
pdata = merge(pdata, pam, by = "sample")
table(pdata$pam,pdata$subtype)

tmp[1:3,1:3]
km_model <- kmeans(tmp, centers = 2, nstart = 25)
km_model
km = as.data.frame(km_model$cluster)
table(km$`km_model$cluster`)
km = as.data.frame(km)
km$sample = rownames(km)
colnames(km)[1] = "km"
pdata = merge(pdata, km, by = "sample")
table(pdata$subtype, pdata$pam)

library(ggplot2)
library(ggalluvial)
head(pdata)
rownames(pdata) = pdata$sample
pdata$sample = NULL
pdata$pam = ifelse(pdata$pam == "Pam_C2", "C2","C1")
data = pdata
df <- to_lodes_form(data[,1:ncol(data)],
                    axes = 1:ncol(data),
                    id = "value")

col<- rep(c("#2E9FDF", "#FFD121",  '#f00a36', 
            '#4a8594', '#051736', '#dbe0e3'), 3)#

pdf("../revised/sankey.pdf",width = 8, height = 6)

ggplot(df, aes(x = x, fill=stratum, label=stratum,
               stratum = stratum, alluvium  = value))+#
  geom_flow(width = 0.3,#
            curve_type = "sine",#
            alpha = 0.5,
            color = 'white',#
            size = 0.1)+#
  geom_stratum(width = 0.28)+#
  geom_text(stat = 'stratum', size = 2, color = 'black')+
  scale_fill_manual(values = col)+#
  theme_void()+#
  theme(legend.position = 'none')#
dev.off()#
