rm(list=ls())
library(RColorBrewer)
library(SingleCellExperiment)

setwd("/Users/yiyang/Desktop/2020-9-19/2020/spatialLIBD/realdata/CRC3_out/");

CRC3 = readRDS(paste0("CRC3_K_8_beta_2.5.rds"))
imagerow = CRC3@colData@listData$imagerow
imagecol = CRC3@colData@listData$imagecol

x_giotto = as.matrix(read.table("CRC3_giotto.txt"))
x_louvain = as.matrix(read.table("CRC3_louvain.txt"))

#################################################################
## adjust the order of clustering
clu1 = colData(CRC3)$x_gmrf
clu2 = colData(CRC3)$x_gmm
clu3 = colData(CRC3)$x_bs
clu4 = x_giotto
clu5 = x_louvain

for (k in 1:10){
  for (l in 1:7){
    for (m in (l+1):8){
      tmp = clu2
      tmp[tmp==l] = 0
      tmp[tmp==m] = l
      tmp[tmp==0] = m
      if (sum(clu1==clu2) < sum(clu1==tmp)){
        clu2 = tmp
      }
    }
  }
}

for (k in 1:10){
  for (l in 1:7){
    for (m in (l+1):8){
      tmp = clu3
      tmp[tmp==l] = 0
      tmp[tmp==m] = l
      tmp[tmp==0] = m
      if (sum(clu1==clu3) < sum(clu1==tmp)){
        clu3 = tmp
      }
    }
  }
}

for (k in 1:10){
  for (l in 1:7){
    for (m in (l+1):8){
      tmp = clu4
      tmp[tmp==l] = 0
      tmp[tmp==m] = l
      tmp[tmp==0] = m
      if (sum(clu1==clu4) < sum(clu1==tmp)){
        clu4 = tmp
      }
    }
  }
}

for (k in 1:10){
  for (l in 1:7){
    for (m in (l+1):8){
      tmp = clu5
      tmp[tmp==l] = 0
      tmp[tmp==m] = l
      tmp[tmp==0] = m
      if (sum(clu1==clu5) < sum(clu1==tmp)){
        clu5 = tmp
      }
    }
  }
}

####################################################
blue = brewer.pal(9,"Blues")[c(4,5,7)]
Green = brewer.pal(9,"Greens")[c(4,7)]
gray = brewer.pal(9,"Greys")[c(3,5)]
color = brewer.pal(8, "Dark2")
color[c(6,7,8)] = blue
color[c(4,5)] = Green
color[c(1,2)] = gray

dat = data.frame(imagerow,imagecol,factor(clu1, levels = c(1,5,2,7,6,3,8,4), labels = c("Stroma1","Stroma2","Muscle","Epithelial1","Epithelial2","Immune1","Immune2","Immune3")))
names(dat)= c("imagerow","imagecol","cluster")          

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagecol, y=600-imagerow, color=cluster)) +
  geom_point(size = 5) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size=20,face="bold"),
        legend.title = element_text(size=25,face="bold")) + 
  scale_color_manual(values=color)

ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2020/spatialLIBD/realdata/CRC3_selectK_out/heatmap_gmrf_CRC3_reorder_widthheight_imm3.png"), plot = p1, width = 26.5, height = 20, units = "cm",dpi = 600)



dat = data.frame(imagerow,imagecol,factor(clu2, levels = c(1,5,2,7,6,3,8,4), labels = c("Stroma1","Stroma2","Muscle","Epithelial1","Epithelial2","Immune1","Immune2","Immune3")))
names(dat)= c("imagerow","imagecol","cluster")

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagecol, y=600-imagerow, color=cluster)) +
  geom_point(size = 5) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size=20,face="bold"),
        legend.title = element_text(size=25,face="bold")) + 
  scale_color_manual(values=color)

ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2020/spatialLIBD/realdata/CRC3_selectK_out/heatmap_gmm_CRC3_reorder_widthheight_imm3.png"), plot = p1, width = 26.5, height = 20, units = "cm",dpi = 600)


dat = data.frame(imagerow,imagecol,factor(clu3, levels = c(1,5,2,7,6,3,8,4), labels = c("Stroma1","Stroma2","Muscle","Epithelial1","Epithelial2","Immune1","Immune2","Immune3")))
names(dat)= c("imagerow","imagecol","cluster")

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagecol, y=600-imagerow, color=cluster)) +
  geom_point(size = 5) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size=20,face="bold"),
        legend.title = element_text(size=25,face="bold")) + 
  scale_color_manual(values=color)

ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2020/spatialLIBD/realdata/CRC3_selectK_out/heatmap_bs_CRC3_reorder_widthheight_imm3.png"), plot = p1, width = 26.5, height = 20, units = "cm", dpi = 600)


dat = data.frame(imagerow,imagecol,factor(clu4, levels = c(1,5,2,7,6,3,8,4), labels = c("Stroma1","Stroma2","Muscle","Epithelial1","Epithelial2","Immune1","Immune2","Immune3")))
names(dat)= c("imagerow","imagecol","cluster")

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagecol, y=600-imagerow, color=cluster)) +
  geom_point(size = 5) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size=20,face="bold"),
        legend.title = element_text(size=25,face="bold")) + 
  scale_color_manual(values=color)

ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2020/spatialLIBD/realdata/CRC3_selectK_out/heatmap_giotto_CRC3_reorder_widthheight_imm3.png"), plot = p1, width = 26.5, height = 20, units = "cm",dpi = 600)


dat = data.frame(imagerow,imagecol,factor(clu5, levels = c(1,5,2,7,6,3,8,4), labels = c("Stroma1","Stroma2","Muscle","Epithelial1","Epithelial2","Immune1","Immune2","Immune3")))
names(dat)= c("imagerow","imagecol","cluster")

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagecol, y=600-imagerow, color=cluster)) +
  geom_point(size = 5) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size=20,face="bold"),
        legend.title = element_text(size=25,face="bold")) + 
  scale_color_manual(values=color)

ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2020/spatialLIBD/realdata/CRC3_selectK_out/heatmap_louvain_CRC3_reorder_widthheight_imm3.png"), plot = p1, width = 26.5, height = 20, units = "cm",dpi = 600)




dat = data.frame(imagerow,imagecol,factor(colData(CRC3)$x_bs, labels = c("Stroma1","Stroma2","Stroma3","Muscle","Epithelial1","Epithelial2","Immune1","Immune2")))
names(dat)= c("imagerow","imagecol","cluster")

library(ggplot2)
p1 <- ggplot(dat, aes(x=imagecol, y=600-imagerow, color=cluster)) +
  geom_point(size = 5) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank()) + 
  scale_color_manual(values=color)

ggsave(paste0("/Users/yiyang/Desktop/2020-9-19/2020/spatialLIBD/realdata/CRC3_markergene/heatmap_bs_CRC3_enrichment_reorder_widthheight_imm3.png"), plot = p1, width = 20, height = 20, units = "in")

