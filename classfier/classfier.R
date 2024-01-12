rm(list = ls())
load( file = "../code/Fig5/classfier_plot.Rdata")
data = df
rs <- data$PC1
names(rs) <- rownames(data)
rs_data <- data.frame(x=1:length(rs),rs=as.numeric(sort(rs)))
rs_data$Risk <- ifelse(rs_data$rs>=PC_t, "High-risk", "Low-risk")
head(rs_data)
table(rs_data$Risk)
library(ggplot2)
plotA = ggplot(rs_data, aes(x=x,y=rs))+
  geom_point(aes(col=Risk),size=0.5)+
  scale_color_manual(labels=c("S1/S2","S3"), 
                     #guide_legend(guide = NULL), 
                     name="PC1", values =c("grey", "#00A087FF")) + 
  
  # 画竖向虚线
  geom_segment(aes(x = sum(rs_data$Risk=="Low-risk"),
                   y = min(df$PC1), 
                   xend = sum(rs_data$Risk=="Low-risk"), 
                   yend = max(rs_data$rs)), linetype="dashed", size = 0.6)+

theme(axis.title.x=element_blank()) +
  scale_x_continuous(limits = c(0,NA),expand = c(0,0)) +
  labs(y="PC1",x="",fill="Risk") +
  #scale_colour_discrete(name="Risk scores") +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #
        axis.text.x=element_blank())


plotA

data = df[df$Ba == "S3",]
rs <- data$CTL
names(rs) <- rownames(data)
rs_data_1 <- data.frame(x=1:length(rs),rs=as.numeric(sort(rs,decreasing = T)))

rs_data_1$Risk <- "S3"
head(rs_data_1)
table(rs_data_1$Risk)


data = df[df$Ba != "S3",]
rs <- data$CTL
names(rs) <- rownames(data)
rs_data_2 <- data.frame(x=129:360,rs=as.numeric(sort(rs)))

rs_data_2$Risk <- ifelse(rs_data_2$rs>=ctl_t, "S1", "S2")
head(rs_data_2)
tmp = rbind(rs_data_1,rs_data_2)
table(tmp$Risk)

plotB = ggplot(tmp, aes(x=x,y=rs))+
  geom_point(aes(col=Risk),size=0.5)+
  scale_color_manual(labels=c("S1","S2","S3"), 
                     #guide_legend(guide = NULL), 
                     name="CTL", values =c("#FFD121","#2E9FDF",
                                           "#088247")) + 
  

  geom_segment(aes(x = 282,
                   y = min(df$CTL), 
                   xend = 282, 
                   yend = max(df$CTL)), linetype="dashed", size = 0.6)+

  
theme(axis.title.x=element_blank()) +
  scale_x_continuous(limits = c(0,NA),expand = c(0,0)) +
  labs(y="CTL",x="",fill="Risk") +
  #scale_colour_discrete(name="Risk scores") +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #
        axis.text.x=element_blank())
plotB