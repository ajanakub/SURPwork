## Bukola Ajanaku
## July 12, 2021
## finalPlay.R putting the code all together
## on screen 1 as DNAMethFin
## module load udunits proj gdal geos
## module load R/4.0.2
## R
## save(list = ls(), file = "/sc/arion/work/ajanab01/newFinalPlayData.RDATA")

.libPaths(c("/hpc/users/ajanab01/.Rlib", .libPaths()))

suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(decorate))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(dplyr))

options(xtable.type="html")

knitr::opts_chunk$set(
  echo=TRUE,
  warning=FALSE,
  message=TRUE,
  error = FALSE,
  tidy = FALSE,
  cache = TRUE,
  cache.lazy = FALSE,
  dev = c("png", "pdf"),
  fig.width=7, fig.height=7)

options(markdown.HTML.stylesheet = 'css/custom.css')

register(SnowParam(4, "SOCK", progressbar=TRUE))

library(EnsDb.Hsapiens.v75)

ensdb=EnsDb.Hsapiens.v75

# loading data
response = readRDS("/sc/arion/projects/epigenAD/Bukola/DNA_methylation/response.RDS" )
metadata = readRDS("/sc/arion/projects/epigenAD/Bukola/DNA_methylation/metadata.RDS")
featureLocation = readRDS("/sc/arion/projects/epigenAD/Bukola/DNA_methylation/featureLocation.RDS")

# creates factors for the target conditions
metadata$Dx = factor(metadata$Dx, c("Control", "Schizo"))

# making sure each feature has a width of at least two
end(featureLocation) = end(featureLocation) + 1

# Calculating residuals (become residValues):
design = model.matrix(~ Dx + Age + Race + negControl_PC1 + negControl_PC2 +
  negControl_PC3 + negControl_PC4, metadata)

# Getting cleaned (residualed) data
fit = lmFit(response, design)
fit = eBayes(fit)
residValues = residuals( fit, response)

# compact function used to produce treeList, lusters, and clstScore
# treeList = clusters peaks for methylation markers of CRDs for ONLY the features
#            that are adjacent to one another (adjacent hierarchical clustering)
# lusters = actually makes different clusters of genes based on the
#            mean clusters sizes we choose (for this script: 10, 25, 50, 100)
# clstScore = Will score clusters based on their correlation in the structure.
#             IMPORTANT because this gives the MAC and LEF scores.
call_clusters_collapse = function(input_mat, peaksGR, meanCluster_Size){
  treeList = runOrderedClusteringGenome(input_mat, peaksGR)
  names(peaksGR)=peaksGR$names

  treeListClusters = createClusters(treeList, method = "meanClusterSize",
    meanClusterSize = meanCluster_Size)

  clstScore = scoreClusters(treeList, treeListClusters)

  return(list(treeList, treeListClusters, clstScore))
}
 # treeList = runOrderedClusteringGenome(residValues, peaks_GR)
peaks = as.data.frame(featureLocation)
peaks$PeaksID = rownames(peaks)
peaks$names = peaks$PeakIDs
peaks_GR = as(peaks, "GRanges")

clstSizer = c(100, 125, 150, 200, 250, 300)
mainRunner <- call_clusters_collapse(residValues, peaks_GR, clstSizer)
# changed featureLocation to peaks_GR here

treeList = mainRunner[[1]]
treeListClusters = mainRunner[[2]]
clstScore = mainRunner[[3]]


 clstScore = as.data.frame(clstScore)
 clearing out all the clusters with only 1 feature
 cleanclst100 = clstScore[["100"]][clstScore[["100"]]$N != 1, ]
 cleanclst125 = clstScore[["125"]][clstScore[["125"]]$N != 1, ]
 cleanclst150 = clstScore[["150"]][clstScore[["150"]]$N != 1, ]
 cleanclst200 = clstScore[["200"]][clstScore[["200"]]$N != 1, ]
 cleanclst250 = clstScore[["250"]][clstScore[["250"]]$N != 1, ]
 cleanclst300 = clstScore[["300"]][clstScore[["300"]]$N != 1, ]

 clstScore = list(cleanclst100, cleanclst125, cleanclst150, cleanclst200, cleanclst250, cleanclst300)
 names(clstScore) = c("100", "125", "150", "200", "250", "300")

n_clusters = countClusters( treeListClusters )

clustInclude = retainClusters(clstScore, "LEF", 0.05 )
treeListClusters_filter = filterClusters( treeListClusters, clustInclude )

treeListClusters_collapse = collapseClusters( treeListClusters_filter, featureLocation )

# n_clusters = countClusters( treeListClusters_collapse )
# ecdBox = evalDiffCorr( residValues, metadata$Dx, featureLocation, treeListClusters_collapse, npermute, method = "Box.permute", method.corr="spearman")
# df = summary( ecdBox )
# df_results = combineResults( ecdBox, clstScore, treeListClusters, featureLocation)

# clstScore_filter = scoreClusters(treeList, treeListClusters_filter)
# clstScore <- clstScore_filter

# clusterSize = "50"
# MACscore = 0.17

# chrome = "chr1"
gatherData <- function(clusterSize, chrome, MACscore, LEFscore, rangel, rangeh) {
    # clusterSize can only be one of c("10", "25", "50", "100")
    # chrome can only be chr1 to chr 22
    # rangel: (lower) starting position for range of clusters (regarding number of clusters)
    # rangel: (higher) ending position for range of clusters (regarding number of clusters)

  df <- clstScore[[clusterSize]][clstScore[[clusterSize]]$chrom == chrome,]

  MACscoreHL <- MACscore + 0.05
  LEFscoreHL <- LEFscore + 0.05
  MACscoreLL <- MACscore - 0.05
  LEFscoreLL <- LEFscore - 0.05

  newdf <- df[df$mean_abs_corr >= MACscoreLL & df$mean_abs_corr <= MACscoreHL,]
  newdf <- newdf[newdf$LEF >= LEFscoreLL & newdf$LEF <= LEFscoreHL,]

  otherdf <- as.data.frame(treeListClusters_filter[[clusterSize]][[chrome]])

  colnames(otherdf) <- "numberCluster"
  otherdf$ID <- rownames(otherdf)

  probe <- unique(newdf$cluster)[rangel : rangeh]
  otherdf <- otherdf[otherdf$numberCluster %in% probe,]

  simLocation <- peaks_GR[peaks$PeaksID %in% otherdf$ID]
  query <- range(simLocation)
  print(identical(otherdf$ID, names(simLocation)))

  fig1 = plotDecorate(ensdb, treeList, treeListClusters_filter, peaks_GR, query)
  fig1
}

# newFigure <- gatherData("10", "chr12", 1, 3)
# sh sch_mount.sh mount "/sc/arion/projects/epigenAD/Bukola"
# pdf("/sc/arion/projects/epigenAD/Bukola/finalPlayrun.pdf")
# newFigure
# dev.off()

# hw:
# gatherData <- function(clusterSize, chrome, MACscore, LEFscore, rangel, rangeh)

#### Making scatterplots: ------------------------------------------------------
clstScore = list(cleanclst100, cleanclst125, cleanclst150, cleanclst200, cleanclst250, cleanclst300)
names(clstScore) = c("100", "125", "150", "200", "250", "300")

clstScore[["100"]]$quantGroup <- NA
 for(i in 1:length(clstScore[["100"]]$N)){
  if(clstScore[["100"]]$N[i] >= quantile(clstScore[["100"]]$N, probs = seq(0,1,0.20))[[1]] && clstScore[["100"]]$N[i] <= quantile(clstScore[["100"]]$N, probs = seq(0,1,0.20))[[2]]){
    clstScore[["100"]][i,]$quantGroup <- 0.2
  } else if(clstScore[["100"]]$N[i] >= quantile(clstScore[["100"]]$N, probs = seq(0,1,0.20))[[2]] && clstScore[["100"]]$N[i] <= quantile(clstScore[["100"]]$N, probs = seq(0,1,0.20))[[3]]){
      clstScore[["100"]][i,]$quantGroup <- 0.4
  } else if(clstScore[["100"]]$N[i] >= quantile(clstScore[["100"]]$N, probs = seq(0,1,0.20))[[3]] && clstScore[["100"]]$N[i] <= quantile(clstScore[["100"]]$N, probs = seq(0,1,0.20))[[4]]){
      clstScore[["100"]][i,]$quantGroup <- 0.6
  } else if(clstScore[["100"]]$N[i] >= quantile(clstScore[["100"]]$N, probs = seq(0,1,0.20))[[4]] && clstScore[["100"]]$N[i] <= quantile(clstScore[["100"]]$N, probs = seq(0,1,0.20))[[5]]){
      clstScore[["100"]][i,]$quantGroup <- 0.8
  } else if(clstScore[["100"]]$N[i] >= quantile(clstScore[["100"]]$N, probs = seq(0,1,0.20))[[5]] && clstScore[["100"]]$N[i] <= quantile(clstScore[["100"]]$N, probs = seq(0,1,0.20))[[6]]){
      clstScore[["100"]][i,]$quantGroup <- 1
  }
 }

clstScore[["125"]]$quantGroup <- NA
  for(i in 1:length(clstScore[["125"]]$N)){
   if(clstScore[["125"]]$N[i] >= quantile(clstScore[["125"]]$N, probs = seq(0,1,0.20))[[1]] && clstScore[["125"]]$N[i] <= quantile(clstScore[["125"]]$N, probs = seq(0,1,0.20))[[2]]){
     clstScore[["125"]][i,]$quantGroup <- 0.2
   } else if(clstScore[["125"]]$N[i] >= quantile(clstScore[["125"]]$N, probs = seq(0,1,0.20))[[2]] && clstScore[["125"]]$N[i] <= quantile(clstScore[["125"]]$N, probs = seq(0,1,0.20))[[3]]){
       clstScore[["125"]][i,]$quantGroup <- 0.4
   } else if(clstScore[["125"]]$N[i] >= quantile(clstScore[["125"]]$N, probs = seq(0,1,0.20))[[3]] && clstScore[["125"]]$N[i] <= quantile(clstScore[["125"]]$N, probs = seq(0,1,0.20))[[4]]){
       clstScore[["125"]][i,]$quantGroup <- 0.6
   } else if(clstScore[["125"]]$N[i] >= quantile(clstScore[["125"]]$N, probs = seq(0,1,0.20))[[4]] && clstScore[["125"]]$N[i] <= quantile(clstScore[["125"]]$N, probs = seq(0,1,0.20))[[5]]){
       clstScore[["125"]][i,]$quantGroup <- 0.8
   } else if(clstScore[["125"]]$N[i] >= quantile(clstScore[["125"]]$N, probs = seq(0,1,0.20))[[5]] && clstScore[["125"]]$N[i] <= quantile(clstScore[["125"]]$N, probs = seq(0,1,0.20))[[6]]){
       clstScore[["125"]][i,]$quantGroup <- 1
   }
  }

  clstScore[["150"]]$quantGroup <- NA
    for(i in 1:length(clstScore[["150"]]$N)){
     if(clstScore[["150"]]$N[i] >= quantile(clstScore[["150"]]$N, probs = seq(0,1,0.20))[[1]] && clstScore[["150"]]$N[i] <= quantile(clstScore[["150"]]$N, probs = seq(0,1,0.20))[[2]]){
       clstScore[["150"]][i,]$quantGroup <- 0.2
     } else if(clstScore[["150"]]$N[i] >= quantile(clstScore[["150"]]$N, probs = seq(0,1,0.20))[[2]] && clstScore[["150"]]$N[i] <= quantile(clstScore[["150"]]$N, probs = seq(0,1,0.20))[[3]]){
         clstScore[["150"]][i,]$quantGroup <- 0.4
     } else if(clstScore[["150"]]$N[i] >= quantile(clstScore[["150"]]$N, probs = seq(0,1,0.20))[[3]] && clstScore[["150"]]$N[i] <= quantile(clstScore[["150"]]$N, probs = seq(0,1,0.20))[[4]]){
         clstScore[["150"]][i,]$quantGroup <- 0.6
     } else if(clstScore[["150"]]$N[i] >= quantile(clstScore[["150"]]$N, probs = seq(0,1,0.20))[[4]] && clstScore[["150"]]$N[i] <= quantile(clstScore[["150"]]$N, probs = seq(0,1,0.20))[[5]]){
         clstScore[["150"]][i,]$quantGroup <- 0.8
     } else if(clstScore[["150"]]$N[i] >= quantile(clstScore[["150"]]$N, probs = seq(0,1,0.20))[[5]] && clstScore[["150"]]$N[i] <= quantile(clstScore[["150"]]$N, probs = seq(0,1,0.20))[[6]]){
         clstScore[["150"]][i,]$quantGroup <- 1
     }
  }

  clstScore[["200"]]$quantGroup <- NA
    for(i in 1:length(clstScore[["200"]]$N)){
     if(clstScore[["200"]]$N[i] >= quantile(clstScore[["200"]]$N, probs = seq(0,1,0.20))[[1]] && clstScore[["200"]]$N[i] <= quantile(clstScore[["200"]]$N, probs = seq(0,1,0.20))[[2]]){
       clstScore[["200"]][i,]$quantGroup <- 0.2
     } else if(clstScore[["200"]]$N[i] >= quantile(clstScore[["200"]]$N, probs = seq(0,1,0.20))[[2]] && clstScore[["200"]]$N[i] <= quantile(clstScore[["200"]]$N, probs = seq(0,1,0.20))[[3]]){
         clstScore[["200"]][i,]$quantGroup <- 0.4
     } else if(clstScore[["200"]]$N[i] >= quantile(clstScore[["200"]]$N, probs = seq(0,1,0.20))[[3]] && clstScore[["200"]]$N[i] <= quantile(clstScore[["200"]]$N, probs = seq(0,1,0.20))[[4]]){
         clstScore[["200"]][i,]$quantGroup <- 0.6
     } else if(clstScore[["200"]]$N[i] >= quantile(clstScore[["200"]]$N, probs = seq(0,1,0.20))[[4]] && clstScore[["200"]]$N[i] <= quantile(clstScore[["200"]]$N, probs = seq(0,1,0.20))[[5]]){
         clstScore[["200"]][i,]$quantGroup <- 0.8
     } else if(clstScore[["200"]]$N[i] >= quantile(clstScore[["200"]]$N, probs = seq(0,1,0.20))[[5]] && clstScore[["200"]]$N[i] <= quantile(clstScore[["200"]]$N, probs = seq(0,1,0.20))[[6]]){
         clstScore[["200"]][i,]$quantGroup <- 1
     }
  }

  clstScore[["250"]]$quantGroup <- NA
    for(i in 1:length(clstScore[["250"]]$N)){
     if(clstScore[["250"]]$N[i] >= quantile(clstScore[["250"]]$N, probs = seq(0,1,0.20))[[1]] && clstScore[["250"]]$N[i] <= quantile(clstScore[["250"]]$N, probs = seq(0,1,0.20))[[2]]){
       clstScore[["250"]][i,]$quantGroup <- 0.2
     } else if(clstScore[["250"]]$N[i] >= quantile(clstScore[["250"]]$N, probs = seq(0,1,0.20))[[2]] && clstScore[["250"]]$N[i] <= quantile(clstScore[["250"]]$N, probs = seq(0,1,0.20))[[3]]){
         clstScore[["250"]][i,]$quantGroup <- 0.4
     } else if(clstScore[["250"]]$N[i] >= quantile(clstScore[["250"]]$N, probs = seq(0,1,0.20))[[3]] && clstScore[["250"]]$N[i] <= quantile(clstScore[["250"]]$N, probs = seq(0,1,0.20))[[4]]){
         clstScore[["250"]][i,]$quantGroup <- 0.6
     } else if(clstScore[["250"]]$N[i] >= quantile(clstScore[["250"]]$N, probs = seq(0,1,0.20))[[4]] && clstScore[["250"]]$N[i] <= quantile(clstScore[["250"]]$N, probs = seq(0,1,0.20))[[5]]){
         clstScore[["250"]][i,]$quantGroup <- 0.8
     } else if(clstScore[["250"]]$N[i] >= quantile(clstScore[["250"]]$N, probs = seq(0,1,0.20))[[5]] && clstScore[["250"]]$N[i] <= quantile(clstScore[["250"]]$N, probs = seq(0,1,0.20))[[6]]){
         clstScore[["250"]][i,]$quantGroup <- 1
     }
  }

  clstScore[["300"]]$quantGroup <- NA
    for(i in 1:length(clstScore[["300"]]$N)){
     if(clstScore[["300"]]$N[i] >= quantile(clstScore[["300"]]$N, probs = seq(0,1,0.20))[[1]] && clstScore[["300"]]$N[i] <= quantile(clstScore[["300"]]$N, probs = seq(0,1,0.20))[[2]]){
       clstScore[["300"]][i,]$quantGroup <- 0.2
     } else if(clstScore[["300"]]$N[i] >= quantile(clstScore[["300"]]$N, probs = seq(0,1,0.20))[[2]] && clstScore[["300"]]$N[i] <= quantile(clstScore[["300"]]$N, probs = seq(0,1,0.20))[[3]]){
         clstScore[["300"]][i,]$quantGroup <- 0.4
     } else if(clstScore[["300"]]$N[i] >= quantile(clstScore[["300"]]$N, probs = seq(0,1,0.20))[[3]] && clstScore[["300"]]$N[i] <= quantile(clstScore[["300"]]$N, probs = seq(0,1,0.20))[[4]]){
         clstScore[["300"]][i,]$quantGroup <- 0.6
     } else if(clstScore[["300"]]$N[i] >= quantile(clstScore[["300"]]$N, probs = seq(0,1,0.20))[[4]] && clstScore[["300"]]$N[i] <= quantile(clstScore[["300"]]$N, probs = seq(0,1,0.20))[[5]]){
         clstScore[["300"]][i,]$quantGroup <- 0.8
     } else if(clstScore[["300"]]$N[i] >= quantile(clstScore[["300"]]$N, probs = seq(0,1,0.20))[[5]] && clstScore[["300"]]$N[i] <= quantile(clstScore[["300"]]$N, probs = seq(0,1,0.20))[[6]]){
         clstScore[["300"]][i,]$quantGroup <- 1
     }
  }

# clstScore[["100"]]$quantGroup <- NA
# quantiles <- quantile(clstScore[["100"]]$N, probs = seq(0,1,0.20))
# clstScore[["100"]]$quantGroup  <- factor(findInterval(clstScore[["100"]]$N, quantiles))
# clstScore[["100"]]$quantGroup[clstScore[["100"]]$quantGroup == 1] <- 2
# c("100", "125", "150", "200", "250", "300")

### SCATTER PLOTS of  mac and lef values per mean cluster size ------------------
pdf("/sc/arion/projects/epigenAD/Bukola/DNAMethPlots/newSizes/viewScatterPlots.pdf")

  ggplot(clstScore[["100"]], aes(clstScore[["100"]]$mean_abs_corr, clstScore[["100"]]$LEF, color = clstScore[["100"]]$quantGroup)) + geom_point() + labs(title = "For Cluster Size 100",
    x = "MAC", y = "LEF") + labs(color = "Quantile Group")

  ggplot(clstScore[["125"]], aes(clstScore[["125"]]$mean_abs_corr, clstScore[["125"]]$LEF, color = clstScore[["125"]]$quantGroup)) + geom_point() + labs(title = "For Cluster Size 125",
    x = "MAC", y = "LEF") + labs(color = "Quantile Group")

  ggplot(clstScore[["150"]], aes(clstScore[["150"]]$mean_abs_corr, clstScore[["150"]]$LEF, color = clstScore[["150"]]$quantGroup)) + geom_point() + labs(title = "For Cluster Size 150",
    x = "MAC", y = "LEF") + labs(color = "Quantile Group")

  ggplot(clstScore[["200"]], aes(clstScore[["200"]]$mean_abs_corr, clstScore[["200"]]$LEF, color = clstScore[["200"]]$quantGroup)) + geom_point() + labs(title = "For Cluster Size 200",
    x = "MAC", y = "LEF") + labs(color = "Quantile Group")

  ggplot(clstScore[["250"]], aes(clstScore[["250"]]$mean_abs_corr, clstScore[["250"]]$LEF, color = clstScore[["250"]]$quantGroup)) + geom_point() + labs(title = "For Cluster Size 250",
    x = "MAC", y = "LEF") + labs(color = "Quantile Group")

  ggplot(clstScore[["300"]], aes(clstScore[["300"]]$mean_abs_corr, clstScore[["300"]]$LEF, color = clstScore[["300"]]$quantGroup)) + geom_point() + labs(title = "For Cluster Size 300",
    x = "MAC", y = "LEF") + labs(color = "Quantile Group")
dev.off()


### SCATTER PLOTS for MAC and Number of Features -------------------------------
pdf("/sc/arion/projects/epigenAD/Bukola/DNAMethPlots/newSizes/byNviewScatterPlots.pdf")

ggplot(clstScore[["100"]][["chr1"]], aes(clstScore[["100"]][clstScore[["100"]]$chrom == "chr1",]$N,
  clstScore[["100"]][clstScore[["100"]]$chrom == "chr1",]$mean_abs_corr)) + geom_point() + labs(title = "For Cluster Size 100 on Chr1",
  x = "Number of Features", y = "MAC")

ggplot(clstScore[["125"]][["chr1"]], aes(clstScore[["125"]][clstScore[["125"]]$chrom == "chr1",]$N,
    clstScore[["125"]][clstScore[["125"]]$chrom == "chr1",]$mean_abs_corr)) + geom_point() + labs(title = "For Cluster Size 125 on Chr1",
    x = "Number of Features", y = "MAC")

  ggplot(clstScore[["150"]][["chr1"]], aes(clstScore[["150"]][clstScore[["150"]]$chrom == "chr1",]$N,
    clstScore[["150"]][clstScore[["150"]]$chrom == "chr1",]$mean_abs_corr)) + geom_point() + labs(title = "For Cluster Size 150 on Chr1",
    x = "Number of Features", y = "MAC")

  ggplot(clstScore[["200"]][["chr1"]], aes(clstScore[["200"]][clstScore[["200"]]$chrom == "chr1",]$N,
    clstScore[["200"]][clstScore[["200"]]$chrom == "chr1",]$mean_abs_corr)) + geom_point() + labs(title = "For Cluster Size 200 on Chr1",
    x = "Number of Features", y = "MAC")

  ggplot(clstScore[["250"]][["chr1"]], aes(clstScore[["250"]][clstScore[["250"]]$chrom == "chr1",]$N,
    clstScore[["250"]][clstScore[["250"]]$chrom == "chr1",]$mean_abs_corr)) + geom_point() + labs(title = "For Cluster Size 250 on Chr1",
    x = "Number of Features", y = "MAC")

  ggplot(clstScore[["300"]][["chr1"]], aes(clstScore[["300"]][clstScore[["300"]]$chrom == "chr1",]$N,
    clstScore[["300"]][clstScore[["300"]]$chrom == "chr1",]$mean_abs_corr)) + geom_point() + labs(title = "For Cluster Size 300 on Chr1",
    x = "Number of Features", y = "MAC")

dev.off()

# ****** FIRST 5: chr1 clst10 ---------------------------------------------------------
# gatherData <- function(clusterSize, chrome, MACscore, LEFscore, rangel, rangeh)
# high = 3rd quartile = [5]
# low = 1st quartile = [2]

# High MAC Low LEF ***
tenhighmaclowleffirst5 <- gatherData("10", "chr1", as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "mean_abs_corr"])[5]),
  as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[2] + summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[4])/2, 1, 5)
pdf("/sc/arion/projects/epigenAD/Bukola/DNAMethPlots/clst10highMAClowLEF.pdf")
  tenhighmaclowleffirst5
dev.off()

# Low MAC High LEF ** having trouble plotting
tenlowmachighleffirst5 <- gatherData("10", "chr1", as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "mean_abs_corr"])[2]),
  as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[5]), 1, 5)
pdf("/sc/arion/projects/epigenAD/Bukola/DNAMethPlots/clst10lowMAChighLEFonly3.pdf")
  tenlowmachighleffirst5
dev.off()


# Low MAC High LEF
tenlowmaclowleffirst5 <- gatherData("10", "chr1", as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "mean_abs_corr"])[2]),
  as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[2]), 1, 5)
pdf("/sc/arion/projects/epigenAD/Bukola/DNAMethPlots/clst10lowMAClowLEF.pdf")
  tenlowmaclowleffirst5
dev.off()


# high mac high LEF
tenhighmachighleffirst5 <- gatherData("10", "chr1", as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "mean_abs_corr"])[5]),
  as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[5]), 1, 5)
pdf("/sc/arion/projects/epigenAD/Bukola/DNAMethPlots/clst10highMAChighLEF.pdf")
  tenhighmachighleffirst5
dev.off()

# FIRST 2: chr1 clst10 ---------------------------------------------------------

# High MAC Low LEF
tenhighmaclowleffirst2 <- gatherData("10", "chr1", as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "mean_abs_corr"])[5]),
  as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[3] + summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[4])/2, 1, 2)
pdf("/sc/arion/projects/epigenAD/Bukola/DNAMethPlots/clst10highMAClowLEFFIRST2.pdf")
  tenhighmaclowleffirst2
dev.off()


# Low MAC High LEF
tenlowmachighleffirst2 <- gatherData("10", "chr1", as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "mean_abs_corr"])[2]),
  as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[5]), 1, 2)
pdf("/sc/arion/projects/epigenAD/Bukola/DNAMethPlots/clst10lowMAChighLEFFIRST2.pdf")
  tenlowmachighleffirst2
dev.off()

# low mac low LEF
tenlowmaclowleffirst2 <- gatherData("10", "chr1", as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "mean_abs_corr"])[2]),
  as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[2]), 1, 2)
pdf("/sc/arion/projects/epigenAD/Bukola/DNAMethPlots/clst10lowMAClowLEFFIRST2.pdf")
  tenlowmaclowleffirst2
dev.off()


# high mac high LEF
tenhighmachighleffirst2 <- gatherData("10", "chr1", as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "mean_abs_corr"])[5]),
  as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[5]), 1, 2)
pdf("/sc/arion/projects/epigenAD/Bukola/DNAMethPlots/clst10highMAChighLEFFIRST2.pdf")
  tenhighmachighleffirst2
dev.off()


## new hw 7/15/2021:
## FIRST MAKE SURE THAT THE FEATURES = 1 THING FIXES
## 1) Re-run everything following clstScore_filter. Just resave clstScore with the filter version.
## 2) Get CTCf data from Kiran: "Acetylated Chromatin Domains Link Chromosomal Organization to
## Cell- and Circuit-level Dysfunction in Schizophrenia and Bipolar Disorder"
##          Use Figure 3 to compare. Function the code starts at line 38.


## starting the CTCF code here:

## Issues:
## clstScore still has N = 1, no fix (significant problem for my previous graphs) [Remove them myself, Kiran will run and check]
## Going through featurelocation for each mean clst size and giving them methylation markers

CTCF=read.table("/sc/arion/projects/epigenAD/Bukola//h1_neuron_ENCFF372JOV.txt")
colnames(CTCF)=c("chr","start","end")
CTCF$GeneID=paste0("peak_",1:dim(CTCF)[1])
CTCF=as(CTCF,"GRanges")
CTCF=sort(CTCF)
​
newCRD100 = featureLocation
newCRD100$cluster <- NA
newCRD125 = featureLocation
newCRD125$cluster <- NA
newCRD150 = featureLocation
newCRD150$cluster <- NA
newCRD200 = featureLocation
newCRD200$cluster <- NA
newCRD250 = featureLocation
newCRD250$cluster <- NA
newCRD300 = featureLocation
newCRD300$cluster <- NA

group <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
  "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
  "chr19", "chr20", "chr21", "chr22")


## for CRD cluster size 100:
for (m in group){
  df100 = as.data.frame(treeListClusters_filter[["100"]][[m]])
  df100$ID=rownames(df100)
  names(df100) <- c("cluster", "ID")
  newCRD100 = as.data.frame(featureLocation)
  newCRD100 <- newCRD100[newCRD100$seqnames == m,]
  newCRD100$ID <- rownames(newCRD100)
  assign(paste("forc100and", m, sep = ""), merge(newCRD100, df100, by = "ID"))
}

finCRD100 = rbind(forc100andchr1, forc100andchr2, forc100andchr3, forc100andchr4, forc100andchr5,
  forc100andchr6, forc100andchr7, forc100andchr8, forc100andchr9, forc100andchr10, forc100andchr11,
  forc100andchr12, forc100andchr13, forc100andchr14, forc100andchr15, forc100andchr16, forc100andchr17,
  forc100andchr18, forc100andchr19, forc100andchr20, forc100andchr21, forc100andchr22)

## for CRD cluster size 125:
for (m in group){
  df125 = as.data.frame(treeListClusters_filter[["125"]][[m]])
  df125$ID=rownames(df125)
  names(df125) <- c("cluster", "ID")
  newCRD125 = as.data.frame(featureLocation)
  newCRD125 <- newCRD125[newCRD125$seqnames == m,]
  newCRD125$ID <- rownames(newCRD125)
  assign(paste("forc125and", m, sep = ""), merge(newCRD125, df125, by = "ID"))
}

finCRD125 = rbind(forc125andchr1, forc125andchr2, forc125andchr3, forc125andchr4, forc125andchr5,
  forc125andchr6, forc125andchr7, forc125andchr8, forc125andchr9, forc125andchr10, forc125andchr11,
  forc125andchr12, forc125andchr13, forc125andchr14, forc125andchr15, forc125andchr16, forc125andchr17,
  forc125andchr18, forc125andchr19, forc125andchr20, forc125andchr21, forc125andchr22)

## for CRD cluster size 150:
for (m in group){
  df150 = as.data.frame(treeListClusters_filter[["150"]][[m]])
  df150$ID=rownames(df150)
  names(df150) <- c("cluster", "ID")
  newCRD150 = as.data.frame(featureLocation)
  newCRD150 <- newCRD150[newCRD150$seqnames == m,]
  newCRD150$ID <- rownames(newCRD150)
  assign(paste("forc150and", m, sep = ""), merge(newCRD150, df150, by = "ID"))
}

finCRD150 = rbind(forc150andchr1, forc150andchr2, forc150andchr3, forc150andchr4, forc150andchr5,
  forc150andchr6, forc150andchr7, forc150andchr8, forc150andchr9, forc150andchr10, forc150andchr11,
  forc150andchr12, forc150andchr13, forc150andchr14, forc150andchr15, forc150andchr16, forc150andchr17,
  forc150andchr18, forc150andchr19, forc150andchr20, forc150andchr21, forc150andchr22)

## for CRD cluster size 200:
for (m in group){
  df200 = as.data.frame(treeListClusters_filter[["200"]][[m]])
  df200$ID=rownames(df200)
  names(df200) <- c("cluster", "ID")
  newCRD200 = as.data.frame(featureLocation)
  newCRD200 <- newCRD200[newCRD200$seqnames == m,]
  newCRD200$ID <- rownames(newCRD200)
  assign(paste("forc200and", m, sep = ""), merge(newCRD200, df200, by = "ID"))
}

finCRD200 = rbind(forc200andchr1, forc200andchr2, forc200andchr3, forc200andchr4, forc200andchr5,
  forc200andchr6, forc200andchr7, forc200andchr8, forc200andchr9, forc200andchr10, forc200andchr11,
  forc200andchr12, forc200andchr13, forc200andchr14, forc200andchr15, forc200andchr16, forc200andchr17,
  forc200andchr18, forc200andchr19, forc200andchr20, forc200andchr21, forc200andchr22)

## for CRD cluster size 250:
for (m in group){
  df250 = as.data.frame(treeListClusters_filter[["250"]][[m]])
  df250$ID=rownames(df250)
  names(df250) <- c("cluster", "ID")
  newCRD250 = as.data.frame(featureLocation)
  newCRD250 <- newCRD250[newCRD250$seqnames == m,]
  newCRD250$ID <- rownames(newCRD250)
  assign(paste("forc250and", m, sep = ""), merge(newCRD250, df250, by = "ID"))
}

finCRD250 = rbind(forc250andchr1, forc250andchr2, forc250andchr3, forc250andchr4, forc250andchr5,
  forc250andchr6, forc250andchr7, forc250andchr8, forc250andchr9, forc250andchr10, forc250andchr11,
  forc250andchr12, forc250andchr13, forc250andchr14, forc250andchr15, forc250andchr16, forc250andchr17,
  forc250andchr18, forc250andchr19, forc250andchr20, forc250andchr21, forc250andchr22)

## for CRD cluster size 300:
for (m in group){
  df300 = as.data.frame(treeListClusters_filter[["300"]][[m]])
  df300$ID=rownames(df300)
  names(df300) <- c("cluster", "ID")
  newCRD300 = as.data.frame(featureLocation)
  newCRD300 <- newCRD300[newCRD300$seqnames == m,]
  newCRD300$ID <- rownames(newCRD300)
  assign(paste("forc300and", m, sep = ""), merge(newCRD300, df300, by = "ID"))
}

finCRD300 = rbind(forc300andchr1, forc300andchr2, forc300andchr3, forc300andchr4, forc300andchr5,
  forc300andchr6, forc300andchr7, forc300andchr8, forc300andchr9, forc300andchr10, forc300andchr11,
  forc300andchr12, forc300andchr13, forc300andchr14, forc300andchr15, forc300andchr16, forc300andchr17,
  forc300andchr18, forc300andchr19, forc300andchr20, forc300andchr21, forc300andchr22)


finCRD100 <- as(finCRD100, "GRanges")
finCRD125 <- as(finCRD125, "GRanges")
finCRD150 <- as(finCRD150, "GRanges")
finCRD200 <- as(finCRD200, "GRanges")
finCRD250 <- as(finCRD250, "GRanges")
finCRD300 <- as(finCRD300, "GRanges")

allCRDs = list(finCRD100, finCRD125, finCRD150, finCRD200, finCRD250, finCRD300)
names(allCRDs) = c("100","125","150","200", "250", "300")

## continuing with function for CTCF
CTCF_density=function(query_grange,CTCF,CTCF_dist_bins,dist_window,direction){
dist_window=dist_window
CTCF_dist_bins=CTCF_dist_bins
counts=matrix(0,nrow=dist_window)
dist_crds=matrix(0,nrow=dist_window)
dist_crds[1]=0
query_bins_list=list()
if (direction == "left"){
df1=as(data.frame("chr"=seqnames(query_grange),"end"=start(query_grange),"start"=start(query_grange)-CTCF_dist_bins),"GRanges")
counts[1]=uniqueN(as.data.frame(findOverlaps(CTCF,df1))[,1])/length(CTCF)
} else {
df1=as(data.frame("chr"=seqnames(query_grange),"start"=end(query_grange),"end"=end(query_grange)+CTCF_dist_bins),"GRanges")
counts[1]=uniqueN(as.data.frame(findOverlaps(CTCF,df1))[,1])/length(CTCF)
}
query_bins_list[[1]]=df1
for (i in (2:dist_window)) {
   if (direction == "left"){
    df1=as(data.frame("chr"=seqnames(df1),"end"=start(df1),"start"=start(df1)-CTCF_dist_bins),"GRanges")
    end(df1)=end(df1)-1
    counts[i]=uniqueN(as.data.frame(findOverlaps(CTCF,df1))[,1])/length(CTCF)
    dist_crds[i]=-1*i*CTCF_dist_bins
    } else {
    df1=as(data.frame("chr"=seqnames(df1),"start"=end(df1),"end"=end(df1)+CTCF_dist_bins),"GRanges")
    start(df1)=start(df1)+1
    counts[i]=uniqueN(as.data.frame(findOverlaps(CTCF,df1))[,1])/length(CTCF)
    dist_crds[i]=1*i*CTCF_dist_bins
    }
    query_bins_list[[i]]=df1
}
return(list(query_bins_list,dist_crds,counts))
}

## Running CTCF function:

CRD_CTCF_density=lapply(1:length(allCRDs),function(x){
  CRD_left=CTCF_density(allCRDs[[x]],CTCF,1000,200,"left")
  CRD_right=CTCF_density(allCRDs[[x]],CTCF,1000,200,"right")
  density_crd_count=c(rev(CRD_left[[3]]),CRD_right[[3]])
  dist_crd_count=c(rev(CRD_left[[2]]),CRD_right[[2]])
  data.frame("Density_CRD"=density_crd_count,"Distance"=dist_crd_count)
})

​names(CRD_CTCF_density)=names(allCRDs)
CRD_CTCF_density_df = reshape2::melt(CRD_CTCF_density,id=c("Density_CRD","Distance"))

## Printing CRD Densities:

pdf("/sc/arion/projects/epigenAD/Bukola/DNAMethPlots/newSizes/plottingCRDDensitywithlines.pdf", width = 20,
  height = 5)
  ggplot(CRD_CTCF_density_df, aes(Distance, Density_CRD, color = L1)) + geom_point() + geom_line()
dev.off()

# pdf("/sc/arion/projects/epigenAD/Bukola/plottingCRDDensity10.pdf")
#  ggplot(CRD_CTCF_density_df[CRD_CTCF_density_df$L1 == 10,], aes(Distance, Density_CRD, color = L1)) + geom_point()
# dev.off()

table(CRD_CTCF_density_df$L1)

# Making histograms:
clusterGR=unlist(range(split(allCRDs[[1]], ~cluster)))
pdf("/sc/arion/projects/epigenAD/Bukola/DNAMethPlots/newSizes/histfornewSizes.pdf")
  hist(width(clusterGR),100)
dev.off()

# Making histograms:
clusterGR=unlist(range(split(allCRDs[[6]], ~cluster)))
pdf("/sc/arion/projects/epigenAD/Bukola/DNAMethPlots/newSizes/histfornewSizes300.pdf")
  hist(width(clusterGR),100)
dev.off()
