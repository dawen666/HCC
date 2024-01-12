
rm(list = ls())
load(file = "../data/HCC_counts.Rdata")

library(tidyverse)
library(DESeq2)





exprSet[1:3,1:3]
IM_counts = as.data.frame(t(exprSet))
IM_counts[1:3,1:3]


load(file = "../data/HCC_pdata.Rdata")
dim(pdata)
pdata = pdata[,c('sample',"CTL_subtype")]
table(pdata$CTL_subtype)
pdata$CTL_subtype = ifelse(pdata$CTL_subtype == "S1","S1","S2")
pdata = na.omit(pdata)
table(pdata$CTL_subtype)
dim(pdata)
rownames(IM_counts) = substr(rownames(IM_counts),1,12)
IM_counts = IM_counts[rownames(IM_counts) %in% pdata$sample,]
dim(IM_counts)
dim(pdata)
km = arrange(pdata, desc(pdata$CTL_subtype))
colData = data.frame(row.names = km$sample,
                     subtype = km$CTL_subtype)


IM_counts = IM_counts[match(rownames(colData),rownames(IM_counts)),]
IM_counts = as.data.frame(t(IM_counts))
head(colData)
IM_counts[1:3,1:3]
table(colData$subtype)
expr_df = IM_counts

library(DESeq2)
dim(expr_df)
expr_df[1:3,1:3]
expr_df = ceiling(expr_df)
pick_row <- apply( expr_df, 1, function(x){
  sum(x == 0) < 366
})
expr_df <- expr_df[pick_row, ]
dds <-DESeqDataSetFromMatrix(countData=expr_df, 
                             colData=colData, 
                             design=~subtype,
                             tidy=FALSE)

dds <- DESeq(dds)
dds 
vsd <- vst(dds)
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))

contrast <- c("subtype", "S2", "S1")
dd1 <- results(dds, contrast=contrast, alpha = 0.05)
plotMA(dd1, ylim=c(-2,2))
resultsNames(dds)
dd2 <- lfcShrink(dds, coef="subtype_S2_vs_S1", type = "apeglm")
plotMA(dd2, ylim=c(-2,2))
summary(dd2, alpha = 0.05)

library(dplyr)
library(tibble)
res <- dd2 %>% 
  data.frame() %>% 
  rownames_to_column("gene_id")
res0rdered <- res[order(res$padj),]
head(res0rdered)

DEG=as.data.frame(res0rdered)
DEG=na.omit(DEG)

logFC_cutoff <- with( DEG, mean( abs( log2FoldChange ) ) + 2 * sd( abs( log2FoldChange ) ) )
logFC_cutoff = 1
diffSig <- DEG[with(DEG, (abs(log2FoldChange)>logFC_cutoff & padj < 0.05 )), ]
DEG$change = as.factor( ifelse( DEG$padj< 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
                                ifelse( DEG$log2FoldChange > logFC_cutoff , 'UP', 'DOWN' ), 'NOT' ) )
table(DEG$change)
library(ggplot2)
colnames(DEG)

ggplot(data = DEG, aes(x = log2FoldChange, y = -log10(padj), color = change)) +
  geom_point(size = 1) +  #绘制散点图
  scale_color_manual(values = c('red', 'gray', 'green'), limits = c('UP', 'NOT', 'DOWN')) +  #
  labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'C2 vs C1', color = '') +  #
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  #
  geom_hline(yintercept = 2, lty = 3, color = 'black') +
  xlim(-12, 12) + ylim(0, 35)  #



load(file = "../data/HCC_S1_DEG.Rdata")
load("../data/HCC_WGCNA_features.Rdata")

DEG=as.data.frame(res0rdered)
DEG=na.omit(DEG)


logFC_cutoff <- with( DEG, mean( abs( log2FoldChange ) ) + 2 * sd( abs( log2FoldChange ) ) )
logFC_cutoff = 1
diffSig <- DEG[with(DEG, (abs(log2FoldChange)>logFC_cutoff & padj < 0.05 )), ]
diffSig = diffSig[diffSig$log2FoldChange < 0 ,]
table(diffSig$gene_id %in% C_gene | diffSig$gene_id %in% D_gene)
diffgene = diffSig[diffSig$gene_id %in% C_gene | diffSig$gene_id %in% D_gene,]



load(file = "../data/HCC_tpm.Rdata")
exprSet[1:3,1:3]
table(rownames(exprSet) %in% diffgene$gene_id)
exprSet = exprSet[rownames(exprSet) %in% diffgene$gene_id,]
exprSet[1:3,1:3]
exprSet = as.data.frame(t(exprSet))
exprSet$sample = rownames(exprSet)

load(file = "../data/HCC_pdata.Rdata")
colnames(pdata)
pdata = pdata[,c("sample","CTL_subtype")]
head(pdata)
pdata = na.omit(pdata)
pdata$CTL_subtype = ifelse(pdata$CTL_subtype == "S1","S1","S2/S3")
tmp = merge(pdata, exprSet, by = "sample")
tmp = tmp[order(pdata$CTL_subtype),]
tmp[1:3,1:3]
rownames(tmp) = tmp$sample
tmp$sample = NULL

annotation_col = data.frame(row.names = rownames(tmp),
                            subtype = tmp$CTL_subtype)

ml = t(scale(tmp[,-1]))
ml[ml > 1] <- 1
ml[ml < -1] <- -1
ml[1:3,1:3]

p <- pheatmap::pheatmap(ml,cluster_cols = F,
                        show_rownames = T,
                        show_colnames = FALSE,
                        annotation_col = annotation_col,
                        color = colorRampPalette(c("#2E9FDF", "white", "#FFD121"))(50))
ggsave(p,filename = "HCC_S1_heatmap.pdf", width = 8,height = 6)






load(file = "./S1_Score.Rdata")
load(file = "../data/HCC_tpm.Rdata")
gsva_gs.up <- GSVA::gsva(data.matrix(exprSet), 
                         immunity, 
                         method="ssgsea"
)
cibersort = as.data.frame(t(gsva_gs.up))
head(cibersort)
cibersort$sample = rownames(cibersort)
head(cibersort)
colnames(cibersort) = c("Score","sample")
load(file = "../data/HCC_pdata.Rdata")
colnames(pdata)
surv = pdata[,c(1,3,2)]
head(surv)
surv = na.omit(surv)

colnames(surv) = c("sample","times","status")
head(surv)

tmp = merge(surv, cibersort, by = "sample")
rownames(tmp) = tmp$sample
tmp$sample = NULL






library(survminer)
library(survival)
bb=tmp
bb$times = bb$times/30
y<-Surv(bb$times,bb$status)
iscutoff<-surv_cutpoint(bb,time = "times",event = "status",variables ="Score")  
t<-iscutoff$cutpoint$cutpoint
palette = "npg"
plot(iscutoff,"Score",palette="npg")
summary(coxph(y~bb$Score,data =bb))
bb$Ba<-ifelse(bb$Score<t,"Low","High")
bb$Ba<-as.factor(bb$Ba)
summary(bb$Ba)
summary(coxph(y~bb$Ba,data =bb))
fit<-survfit(y~bb$Ba,data =bb)
pdf(file = "HCC_S1_survival.pdf")
ggsurvplot(fit, data = bb,pval = TRUE,
           palette = c("#E7B800", "#2E9FDF"),
           legend.title = "HCC",legend.labs = c("High-S1Score", "Low-S1Score"),
           conf.int = TRUE,
           surv.median.line = "hv",
           risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),
           ggtheme = theme_bw(),
           break.time.by = 30
)
dev.off()
