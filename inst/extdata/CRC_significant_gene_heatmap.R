setwd("/Users/yiyang/Desktop/2020-9-19/2020/spatialLIBD/realdata/CRC3_markergene/")

library(Seurat)
set.seed(100)
library(BPSC)
library(SingleCellExperiment)

sig_gene = c()
f = c(1,5,2,7,6,3,8,4)
for (i in 1:8){
  foldchange = read.csv(paste0("CRC3_gmrf_foldchange_cluster_", f[i], ".csv"))
  data = read.csv(paste0("CRC3_betapossion_gmrf_cluster_", f[i],  ".csv"))
  idx = match(data[,1], foldchange[,1])
  fd = foldchange[idx,3]
  print(length(data[fd>0.5,1]))
  sig_gene = c(sig_gene, data[fd>0.5,1])
}

idx_siggene = match(unique(sig_gene), sig_gene)
num_gene_each_cluster = c(19,10,21,6,3,10,8,1)
cum_gene_each_cluster = cumsum(c(19,10,21,6,3,10,8,1))
num_unique_gene_each_cluster = rep(0,8)
for (i in 1){
  num_unique_gene_each_cluster[i] = sum(idx_siggene<cum_gene_each_cluster[i]+1)
}
for (i in 2:8){
  num_unique_gene_each_cluster[i] = sum(idx_siggene<cum_gene_each_cluster[i]+1 &
                                          idx_siggene>cum_gene_each_cluster[i-1])
}
cum_unique_gene_each_cluster = cumsum(num_unique_gene_each_cluster)

data = read.table("CRC3_top2000_normalization.txt")
dim(data)
idx = match(unique(sig_gene), rownames(data))

data2 = as.matrix(data[idx,])
colnames(data2) = rep("", 2988)
library(gplots)

CRC3_sce = readRDS(paste0("/Users/yiyang/Desktop/2020-9-19/2020/spatialLIBD/realdata/CRC3_out/CRC3_K_8_beta_2.5.rds"))
clu1 = colData(CRC3_sce)$x_bs
clu2 = clu1
clu2[clu1==4] = 1
clu2[clu1==2] = 2
clu2[clu1==5] = 3
clu2[clu1==8] = 4
clu2[clu1==7] = 5
clu2[clu1==1] = 6
clu2[clu1==3] = 7
clu2[clu1==6] = 8
out = sort(clu2, index.return = T)
cumsum(table(out$x))

# creates a own color palette from red to green
colors = c(seq(-1.5,1.5,length=300))
my_palette <- colorRampPalette(c("magenta1", "black", "yellow"))(n = 299)

data3 = data2[,out$ix]

library(RColorBrewer)
blue = brewer.pal(9,"Blues")[c(4,5,7)]
Green = brewer.pal(9,"Greens")[c(4,7)]
gray = brewer.pal(9,"Greys")[c(3,5)]
color = brewer.pal(8, "Dark2")
color[c(6,7,8)] = blue
color[c(4,5)] = Green
color[c(1,2)] = gray

len_reorder = c(table(clu2)[1], table(clu2)[2], table(clu2)[3],
                table(clu2)[4], table(clu2)[5], table(clu2)[6], 
                table(clu2)[7], table(clu2)[8])

png(filename = "CRC3_gmrf_betapossion_heatmap_reorder_keysize_2_bar_imm3_SC_to_BS.png", width = 400, height = 360, units = "mm",res=240)
par(cex.main=3)
full<-heatmap.2(data3, breaks = colors, col=my_palette, scale="row",dendrogram="none", Rowv = F, Colv = F, key.par = list(cex.main=1.5, cex.axis = 2), 
                key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1.6,cexCol=1.6, keysize = 0.6,labRow=FALSE,
                margins=c(14,15), colsep = cumsum(len_reorder), lhei=c(0.9,9,0.9),
                ColSideColors = c(rep(color[1], len_reorder[1]), rep(color[2], len_reorder[2]),
                                  rep(color[3], len_reorder[3]), rep(color[4], len_reorder[4]),
                                  rep(color[5], len_reorder[5]), rep(color[6], len_reorder[6]),
                                  rep(color[7], len_reorder[7]), rep(color[8], len_reorder[8])))
text(c(0.13,0.35,0.58,0.69,0.75,0.81,0.88,0.92)*0.95, 0.98, srt = 0, labels = c("S1", "S2", "M", "E1","E2","I1","I2", "I3"), xpd = TRUE, cex = 2)
dev.off()


