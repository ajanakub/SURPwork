## Bukola Ajanaku
## July 12, 2021
## making code based on the loads from Kiran
## on screen 3 as testKiran

# --------- loading data form Kiran
load("/sc/arion/projects/epigenAD/Bukola/DNA_methylation/CRD_DNA_methylation.RDATA")

load("/sc/arion/projects/epigenAD/Bukola/DNA_methylation/CRD_DNA_methylation_plot.RDATA")

# ----------
# newFigure <- gatherData("50", "chr12", 0.17)


library(EnsDb.Hsapiens.v75)

ensdb = EnsDb.Hsapiens.v75

gatherData <- function(clusterSize, chrome, MACscore) {

  df <- clstScore[[clusterSize]][clstScore[[clusterSize]]$chrom == chrome,]
  MACscoreLL <- MACscore - 0.05
  MACscoreHL <- MACscore + 0.05

  newdf <- df[df$mean_abs_corr >= MACscoreLL & df$mean_abs_corr <= MACscoreHL,]

  otherdf <- as.data.frame(Clusters_filter_list[[clusterSize]][[chrome]])

  colnames(otherdf) <- "numberCluster"
  otherdf$ID <- rownames(otherdf)

  probe <- unique(newdf$cluster)[1:5]
  otherdf <- otherdf[otherdf$numberCluster %in% probe,]

  simLocation <- peaks_GR[peaks_GR$names %in% otherdf$ID]
  query <- range(simLocation)
  print(identical(otherdf$ID, names(simLocation)))

  fig1 = plotDecorate( ensdb, Clusters_genome_wide_list, Clusters_filter_list, peaks_GR, query, data=res_mat)
# simLocation$names = names(simLocation)
  fig1
}


##---
call_clusters_collapse=function(input_mat,peaksGR,meanCluster_Size){
treeList = runOrderedClusteringGenome(input_mat,peaksGR)
names(peaksGR)=peaksGR$names
treeListClusters = createClusters(treeList, method = "meanClusterSize", meanClusterSize=meanCluster_Size)
clstScore =scoreClusters(treeList, treeListClusters)
return(list(treeList,treeListClusters,clstScore))
}
