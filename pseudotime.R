rm(list = ls())
#######load data###########
tf<-unique(tf_network[,1])
#####creat sce object##########
library(slingshot)
sim <- SingleCellExperiment(datExpr)
sim
length(intersect(rownames(assays(sim)[[1]]),tf))
colData(sim)$celltype = type
set.seed(111)
library(uwot)
rd1 <- umap(t(log1p(assays(sim)[[1]])))
colnames(rd1) <- c('UMAP1', 'UMAP2')

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

reducedDims(sim) <- SimpleList(UMAP = rd1)
sim


library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
table(cl1)

colData(sim)$GMM <- cl1
sim

cl2 <- kmeans(rd1, centers = 6)$cluster
table(cl2)
table(type)
colData(sim)$kmeans <- cl2
sim <- slingshot(sim, 
                 clusterLabels = 'kmeans',  
                 reducedDim = 'UMAP',  
                 start.clus= "00h",  
                 end.clus = NULL     
)     
sim
