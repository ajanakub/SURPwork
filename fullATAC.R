## Bukola Ajanaku
## July 24, 2021
## fullATAC.R: using ATAC sequence for both control and Schizo
## on screen 3 as fullATAC
## module load udunits proj gdal geos
## module load R/4.0.2
## R
## save(list = ls(), file = "/sc/arion/work/ajanab01/newFinalPlayData.RDATA")

.libPaths(c("/hpc/users/ajanab01/.Rlib", .libPaths()))

suppressPackageStartupMessages(library(decorate))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(BiocParallel))

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

# loading data
load("/sc/arion/projects/epigenAD/Bukola/ATAC_SCZ/Processed.RDATA")


# filtering metaData for only dx of scz and control
keep = metaData$Dx %in% c('SCZ', 'Control')
metaData = metaData[keep,]
chipCounts = chipCounts[,keep]
metaData$Dx = factor(metaData$Dx, c("Control", "SCZ"))

## Processing data:
isexpr = rowSums(cpm(chipCounts)>1) >= 0.2*ncol(chipCounts)
peakLocs2 = peakLocs[which(isexpr)]

# Standard usage of limma/voom
countObj = DGEList( chipCounts[isexpr,] )
countObj = calcNormFactors( countObj )
design = model.matrix( ~ as.character(`ATACSeq_report:Sequencing_Order_ID`) +
`ATACSeq_report:Mean_GC_content`+
`ATACSeq_report:Mapped_Reads` +
`Age of Death` +
`PMI (in hours)` + Sex , metaData)
vobj = voom( countObj, design, plot= FALSE)

# loading some necessary variables for residualizing Data
dcmp = svd(vobj$E, nv=5, nu=0)
frac = dcmp$d^2 / sum(dcmp$d^2) * 100
xlab = paste0('PC1: ', round(frac[1], 1), '%')
ylab = paste0('PC2: ', round(frac[2], 1), '%')

# Residualizing the data.
dsgn = model.matrix( ~ dcmp$v[,1:2] + as.character(`ATACSeq_report:Sequencing_Order_ID`) +
`ATACSeq_report:Mean_GC_content`+
`ATACSeq_report:Mapped_Reads` +
`Age of Death` +
`PMI (in hours)` + Sex , metaData)
fitPC = lmFit(vobj, dsgn)
quantResid = residuals( fitPC, vobj )

vobj2 = voom( countObj, dsgn, plot=FALSE)

fitPC2 = lmFit(vobj2, dsgn)
quantResid2 = residuals( fitPC2, vobj2 )

dsgn = model.matrix( ~ dcmp$v[,1:2] + as.character(`ATACSeq_report:Sequencing_Order_ID`) +
`ATACSeq_report:Mean_GC_content`+
`ATACSeq_report:Mapped_Reads` +
`Age of Death` +
`PMI (in hours)` + Sex + Dx, metaData)

fitDE = lmFit(vobj2, dsgn)

fitDE = eBayes(fitDE)

# Now, for clustering the correlated ATAC sequences
treeList = runOrderedClusteringGenome( quantResid2, peakLocs2)
treeListOriginal = runOrderedClusteringGenome( vobj$E, peakLocs2)
treeListClusters = createClusters( treeList, method='meanClusterSize', meanClusterSize=c(10, 25, 50, 100) )
n_clusters = countClusters( treeListClusters )
clstScore = scoreClusters(treeList, treeListClusters, BPPARAM=SerialParam() )

# Now, dropping all the N = 1s
clstScore = as.data.frame(clstScore)
clearing out all the clusters with only 1 feature
cleanclst10 = clstScore[["10"]][clstScore[["10"]]$N != 1, ]
cleanclst25 = clstScore[["25"]][clstScore[["25"]]$N != 1, ]
cleanclst50 = clstScore[["50"]][clstScore[["50"]]$N != 1, ]
cleanclst100 = clstScore[["100"]][clstScore[["100"]]$N != 1, ]
clstScore = list(cleanclst10, cleanclst25, cleanclst50, cleanclst100)
names(clstScore) = c("10", "25", "50", "100")

# Retaining clusters based on strength
clustInclude = retainClusters( clstScore, "LEF", 0.1 )

# get retained clusters
treeListClusters_filter = filterClusters( treeListClusters, clustInclude )

# Collapse similar clusters
treeListClusters_collapse = collapseClusters( treeListClusters_filter, peakLocs2 )
n_clusters = countClusters( treeListClusters_collapse )

# starting the CTCF code here:
CTCF=read.table("/sc/arion/projects/epigenAD/Bukola//h1_neuron_ENCFF372JOV.txt")
colnames(CTCF)=c("chr","start","end")
CTCF$GeneID=paste0("peak_",1:dim(CTCF)[1])
CTCF=as(CTCF,"GRanges")
CTCF=sort(CTCF)
​
newCRD10 = peakLocs
newCRD10$cluster <- NA
newCRD25 = peakLocs
newCRD25$cluster <- NA
newCRD50 = peakLocs
newCRD50$cluster <- NA
newCRD100 = peakLocs
newCRD100$cluster <- NA

group <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
  "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
  "chr19", "chr20", "chr21", "chr22")

## for CRD cluster size 10:
for (m in group){
    df10 = as.data.frame(treeListClusters_filter[["10"]][[m]])
    df10$ID=rownames(df10)
    names(df10) <- c("cluster", "ID")
    newCRD10 = as.data.frame(peakLocs)
    newCRD10 <- newCRD10[newCRD10$seqnames == m,]
    newCRD10$ID <- rownames(newCRD10)
    assign(paste("forc10and", m, sep = ""), merge(newCRD10, df10, by = "ID"))
}

finCRD10 = rbind(forc10andchr1, forc10andchr2, forc10andchr3, forc10andchr4, forc10andchr5,
  forc10andchr6, forc10andchr7, forc10andchr8, forc10andchr9, forc10andchr10, forc10andchr11,
  forc10andchr12, forc10andchr13, forc10andchr14, forc10andchr15, forc10andchr16, forc10andchr17,
  forc10andchr18, forc10andchr19, forc10andchr20, forc10andchr21, forc10andchr22)


## for CRD cluster size 25:
for (m in group){
  df25 = as.data.frame(treeListClusters_filter[["25"]][[m]])
  df25$ID=rownames(df25)
  names(df25) <- c("cluster", "ID")
  newCRD25 = as.data.frame(peakLocs)
  newCRD25 <- newCRD25[newCRD25$seqnames == m,]
  newCRD25$ID <- rownames(newCRD25)
  assign(paste("forc25and", m, sep = ""), merge(newCRD25, df25, by = "ID"))  }

finCRD25 = rbind(forc25andchr1, forc25andchr2, forc25andchr3, forc25andchr4, forc25andchr5,
  forc25andchr6, forc25andchr7, forc25andchr8, forc25andchr9, forc25andchr10, forc25andchr11,
  forc25andchr12, forc25andchr13, forc25andchr14, forc25andchr15, forc25andchr16, forc25andchr17,
  forc25andchr18, forc25andchr19, forc25andchr20, forc25andchr21, forc25andchr22)

## for CRD cluster size 50:
for (m in group){
  df50 = as.data.frame(treeListClusters_filter[["50"]][[m]])
  df50$ID=rownames(df50)
  names(df50) <- c("cluster", "ID")
  newCRD50 = as.data.frame(peakLocs)
  newCRD50 <- newCRD50[newCRD50$seqnames == m,]
  newCRD50$ID <- rownames(newCRD50)
  assign(paste("forc50and", m, sep = ""), merge(newCRD50, df50, by = "ID"))
}

finCRD50 = rbind(forc50andchr1, forc50andchr2, forc50andchr3, forc50andchr4, forc50andchr5,
  forc50andchr6, forc50andchr7, forc50andchr8, forc50andchr9, forc50andchr10, forc50andchr11,
  forc50andchr12, forc50andchr13, forc50andchr14, forc50andchr15, forc50andchr16, forc50andchr17,
  forc50andchr18, forc50andchr19, forc50andchr20, forc50andchr21, forc50andchr22)

## for CRD cluster size 100:
for (m in group){
  df100 = as.data.frame(treeListClusters_filter[["100"]][[m]])
  df100$ID=rownames(df100)
  names(df100) <- c("cluster", "ID")
  newCRD100 = as.data.frame(peakLocs)
  newCRD100 <- newCRD100[newCRD100$seqnames == m,]
  newCRD100$ID <- rownames(newCRD100)
  assign(paste("forc100and", m, sep = ""), merge(newCRD100, df100, by = "ID"))
}

finCRD100 = rbind(forc100andchr1, forc100andchr2, forc100andchr3, forc100andchr4, forc100andchr5,
  forc100andchr6, forc100andchr7, forc100andchr8, forc100andchr9, forc100andchr10, forc100andchr11,
  forc100andchr12, forc100andchr13, forc100andchr14, forc100andchr15, forc100andchr16, forc100andchr17,
  forc100andchr18, forc100andchr19, forc100andchr20, forc100andchr21, forc100andchr22)

finCRD10 <- as(finCRD10, "GRanges")
finCRD25 <- as(finCRD25, "GRanges")
finCRD50 <- as(finCRD50, "GRanges")
finCRD100 <- as(finCRD100, "GRanges")

allCRDs = list(finCRD10, finCRD25, finCRD50, finCRD100)
names(allCRDs) = c("10","25","50","100")

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

names(CRD_CTCF_density)=names(allCRDs)
CRD_CTCF_density_df = reshape2::melt(CRD_CTCF_density,id=c("Density_CRD","Distance"))

pdf("/sc/arion/projects/epigenAD/Bukola/ATACPlots/Full/CRDDensity.pdf", width = 20,
  height = 5)
  ggplot(CRD_CTCF_density_df, aes(Distance, Density_CRD, color = L1)) + geom_point() + geom_line() + xlim(-2e5, 2e5)
dev.off()

# Making TreeList Plots: -------------------------------------------------------
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

#### Making scatterplots: ------------------------------------------------------
clstScore[["10"]]$quantGroup <- NA
 for(i in 1:length(clstScore[["10"]]$N)){
  if(clstScore[["10"]]$N[i] >= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[1]] && clstScore[["10"]]$N[i] <= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[2]]){
    clstScore[["10"]][i,]$quantGroup <- 0.2
  } else if(clstScore[["10"]]$N[i] >= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[2]] && clstScore[["10"]]$N[i] <= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[3]]){
      clstScore[["10"]][i,]$quantGroup <- 0.4
  } else if(clstScore[["10"]]$N[i] >= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[3]] && clstScore[["10"]]$N[i] <= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[4]]){
      clstScore[["10"]][i,]$quantGroup <- 0.6
  } else if(clstScore[["10"]]$N[i] >= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[4]] && clstScore[["10"]]$N[i] <= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[5]]){
      clstScore[["10"]][i,]$quantGroup <- 0.8
  } else if(clstScore[["10"]]$N[i] >= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[5]] && clstScore[["10"]]$N[i] <= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[6]]){
      clstScore[["10"]][i,]$quantGroup <- 1
  }
 }

clstScore[["25"]]$quantGroup <- NA
  for(i in 1:length(clstScore[["25"]]$N)){
   if(clstScore[["25"]]$N[i] >= quantile(clstScore[["25"]]$N, probs = seq(0,1,0.20))[[1]] && clstScore[["25"]]$N[i] <= quantile(clstScore[["25"]]$N, probs = seq(0,1,0.20))[[2]]){
     clstScore[["25"]][i,]$quantGroup <- 0.2
   } else if(clstScore[["25"]]$N[i] >= quantile(clstScore[["25"]]$N, probs = seq(0,1,0.20))[[2]] && clstScore[["25"]]$N[i] <= quantile(clstScore[["25"]]$N, probs = seq(0,1,0.20))[[3]]){
       clstScore[["25"]][i,]$quantGroup <- 0.4
   } else if(clstScore[["25"]]$N[i] >= quantile(clstScore[["25"]]$N, probs = seq(0,1,0.20))[[3]] && clstScore[["25"]]$N[i] <= quantile(clstScore[["25"]]$N, probs = seq(0,1,0.20))[[4]]){
       clstScore[["25"]][i,]$quantGroup <- 0.6
   } else if(clstScore[["25"]]$N[i] >= quantile(clstScore[["25"]]$N, probs = seq(0,1,0.20))[[4]] && clstScore[["25"]]$N[i] <= quantile(clstScore[["25"]]$N, probs = seq(0,1,0.20))[[5]]){
       clstScore[["25"]][i,]$quantGroup <- 0.8
   } else if(clstScore[["25"]]$N[i] >= quantile(clstScore[["25"]]$N, probs = seq(0,1,0.20))[[5]] && clstScore[["25"]]$N[i] <= quantile(clstScore[["25"]]$N, probs = seq(0,1,0.20))[[6]]){
       clstScore[["25"]][i,]$quantGroup <- 1
   }
  }

  clstScore[["50"]]$quantGroup <- NA
    for(i in 1:length(clstScore[["50"]]$N)){
     if(clstScore[["50"]]$N[i] >= quantile(clstScore[["50"]]$N, probs = seq(0,1,0.20))[[1]] && clstScore[["50"]]$N[i] <= quantile(clstScore[["50"]]$N, probs = seq(0,1,0.20))[[2]]){
       clstScore[["50"]][i,]$quantGroup <- 0.2
     } else if(clstScore[["50"]]$N[i] >= quantile(clstScore[["50"]]$N, probs = seq(0,1,0.20))[[2]] && clstScore[["50"]]$N[i] <= quantile(clstScore[["50"]]$N, probs = seq(0,1,0.20))[[3]]){
         clstScore[["50"]][i,]$quantGroup <- 0.4
     } else if(clstScore[["50"]]$N[i] >= quantile(clstScore[["50"]]$N, probs = seq(0,1,0.20))[[3]] && clstScore[["50"]]$N[i] <= quantile(clstScore[["50"]]$N, probs = seq(0,1,0.20))[[4]]){
         clstScore[["50"]][i,]$quantGroup <- 0.6
     } else if(clstScore[["50"]]$N[i] >= quantile(clstScore[["50"]]$N, probs = seq(0,1,0.20))[[4]] && clstScore[["50"]]$N[i] <= quantile(clstScore[["50"]]$N, probs = seq(0,1,0.20))[[5]]){
         clstScore[["50"]][i,]$quantGroup <- 0.8
     } else if(clstScore[["50"]]$N[i] >= quantile(clstScore[["50"]]$N, probs = seq(0,1,0.20))[[5]] && clstScore[["50"]]$N[i] <= quantile(clstScore[["50"]]$N, probs = seq(0,1,0.20))[[6]]){
         clstScore[["50"]][i,]$quantGroup <- 1
     }
  }

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

  ###* SCATTER PLOTS of  mac and lef values per mean cluster size ------------------
  pdf("/sc/arion/projects/epigenAD/Bukola/ATACPlots/Full/ATACviewScatterPlots.pdf")
    ggplot(clstScore[["10"]], aes(clstScore[["10"]]$mean_abs_corr, clstScore[["10"]]$LEF, color = clstScore[["10"]]$quantGroup)) + geom_point() + labs(title = "For Cluster Size 10",
      x = "MAC", y = "LEF") + labs(color = "Quantile Group")

    ggplot(clstScore[["25"]], aes(clstScore[["25"]]$mean_abs_corr, clstScore[["25"]]$LEF, color = clstScore[["25"]]$quantGroup)) + geom_point() + labs(title = "For Cluster Size 25",
      x = "MAC", y = "LEF") + labs(color = "Quantile Group")

    ggplot(clstScore[["50"]], aes(clstScore[["50"]]$mean_abs_corr, clstScore[["50"]]$LEF, color = clstScore[["50"]]$quantGroup)) + geom_point() + labs(title = "For Cluster Size 50",
      x = "MAC", y = "LEF") + labs(color = "Quantile Group")

    ggplot(clstScore[["100"]], aes(clstScore[["100"]]$mean_abs_corr, clstScore[["100"]]$LEF, color = clstScore[["100"]]$quantGroup)) + geom_point() + labs(title = "For Cluster Size 100",
      x = "MAC", y = "LEF") + labs(color = "Quantile Group")

  dev.off()

  ### SCATTER PLOTS for MAC and Number of Features -------------------------------
  pdf("/sc/arion/projects/epigenAD/Bukola/ATACPlots/Full/ATACbyNviewScatterPlots.pdf")

  ggplot(clstScore[["10"]][["chr1"]], aes(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1",]$N,
    clstScore[["10"]][clstScore[["10"]]$chrom == "chr1",]$mean_abs_corr)) + geom_point() + labs(title = "For Cluster Size 10 on Chr1",
    x = "Number of Features", y = "MAC")

  ggplot(clstScore[["25"]][["chr1"]], aes(clstScore[["25"]][clstScore[["25"]]$chrom == "chr1",]$N,
      clstScore[["25"]][clstScore[["25"]]$chrom == "chr1",]$mean_abs_corr)) + geom_point() + labs(title = "For Cluster Size 25 on Chr1",
      x = "Number of Features", y = "MAC")

    ggplot(clstScore[["50"]][["chr1"]], aes(clstScore[["50"]][clstScore[["50"]]$chrom == "chr1",]$N,
      clstScore[["50"]][clstScore[["50"]]$chrom == "chr1",]$mean_abs_corr)) + geom_point() + labs(title = "For Cluster Size 50 on Chr1",
      x = "Number of Features", y = "MAC")

    ggplot(clstScore[["100"]][["chr1"]], aes(clstScore[["100"]][clstScore[["100"]]$chrom == "chr1",]$N,
      clstScore[["100"]][clstScore[["100"]]$chrom == "chr1",]$mean_abs_corr)) + geom_point() + labs(title = "For Cluster Size 100 on Chr1",
      x = "Number of Features", y = "MAC")

  dev.off()

  # FIRST 5: chr1 clst10 ---------------------------------------------------------
  # gatherData <- function(clusterSize, chrome, MACscore, LEFscore, rangel, rangeh)
  # high = 3rd quartile = [5]
  # low = 1st quartile = [2]

  peaks = as.data.frame(peakLocs2)
  peaks$PeaksID = rownames(peaks)
  peaks$names = peaks$PeakIDs
  peaks_GR = as(peaks, "GRanges")

  library(EnsDb.Hsapiens.v75)

  ensdb= EnsDb.Hsapiens.v75

  # High MAC Low LEF
  tenhighmaclowleffirst5 <- gatherData("10", "chr1", as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "mean_abs_corr"])[5]),
    as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[2]), 1, 5)
  pdf("/sc/arion/projects/epigenAD/Bukola/ATACPlots/Full/ATACclst10highMAClowLEF.pdf")
    tenhighmaclowleffirst5
  dev.off()

  # Low MAC High LEF ** having trouble plotting
  tenlowmachighleffirst5 <- gatherData("10", "chr1", as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "mean_abs_corr"])[1]),
    as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[2]), 1, 3)
  pdf("/sc/arion/projects/epigenAD/Bukola/ATACPlots/Full/ATACclst10lowMAChighLEFonly3.pdf")
    tenlowmachighleffirst5
  dev.off()

  # Low MAC High LEF
  tenlowmaclowleffirst5 <- gatherData("10", "chr1", as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "mean_abs_corr"])[2]),
    as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[2]), 1, 5)
  pdf("/sc/arion/projects/epigenAD/Bukola/ATACPlots/Full/ATACclst10lowMAClowLEF.pdf")
    tenlowmaclowleffirst5
  dev.off()


  # high mac high LEF
  tenhighmachighleffirst5 <- gatherData("10", "chr1", as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "mean_abs_corr"])[5]),
    as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[4]), 1, 5)
  pdf("/sc/arion/projects/epigenAD/Bukola/ATACPlots/Full/ATACclst10highMAChighLEF.pdf")
    tenhighmachighleffirst5
  dev.off()

  # FIRST 2: chr1 clst10 ---------------------------------------------------------

  # High MAC Low LEF
  tenhighmaclowleffirst2 <- gatherData("10", "chr1", as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "mean_abs_corr"])[5]),
    as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[2] + summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[4])/2, 1, 2)
  pdf("/sc/arion/projects/epigenAD/Bukola/ATACPlots/Full/ATACclst10highMAClowLEFFIRST2.pdf")
    tenhighmaclowleffirst2
  dev.off()


  # Low MAC High LEF
  tenlowmachighleffirst2 <- gatherData("10", "chr1", as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "mean_abs_corr"])[1]),
    as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[4]), 1, 2)
  pdf("/sc/arion/projects/epigenAD/Bukola/ATACPlots/Full/ATACclst10lowMAChighLEFFIRST2.pdf")
    tenlowmachighleffirst2
  dev.off()

  # low mac low LEF
  tenlowmaclowleffirst2 <- gatherData("10", "chr1", as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "mean_abs_corr"])[2]),
    as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[2]), 1, 2)
  pdf("/sc/arion/projects/epigenAD/Bukola/ATACPlots/Full/ATACclst10lowMAClowLEFFIRST2.pdf")
    tenlowmaclowleffirst2
  dev.off()


  # high mac high LEF
  tenhighmachighleffirst2 <- gatherData("10", "chr1", as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "mean_abs_corr"])[5]),
    as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[5]), 1, 2)
  pdf("/sc/arion/projects/epigenAD/Bukola/ATACPlots/Full/ATACclst10highMAChighLEFFIRST2.pdf")
    tenhighmachighleffirst2
  dev.off()

  ### NOW, testing differential signals: comparing cases and controls:
  # get total number of clusters
  n_clusters = countClusters( treeListClusters_collapse )

  # Evaluate Differential Correlation between two subsets of data
  param = SnowParam(6, "SOCK", progressbar=TRUE)

  ecdBox = evalDiffCorr( quantResid2, metaData$Dx, peakLocs2, treeListClusters_collapse, BPPARAM=param, method="Box.permute" )

  df = summary( ecdBox )
  # df = head(df[df$id == 10,])
  # print results
  # head(df)
  # extract peak ID's from most significant cluster
  peakIDs = getFeaturesInCluster( treeListClusters_collapse, df$chrom[1], df$cluster[1], df$id[1])
  # get location of peaks in this cluster
  query = range(peakLocs2[names(peakLocs2) %in% peakIDs])
  # expand window to include adjacent clusters
  window = 1e5
  start(query) = start(query) - window
  end(query) = end(query) + window

  pdf("/sc/arion/projects/epigenAD/Bukola/significantClustersATAC.pdf")
    plotDecorate( ensdb, treeList, treeListClusters_collapse, peakLocs2, query, showGenes=TRUE)
  dev.off()

  # Now to compare the level of correlation in this cluster for cases vs control.
  main = paste0(df$chrom[1], ': cluster ', df$cluster[1])
  pdf("/sc/arion/projects/epigenAD/Bukola/plotComparingCorrATAC.pdf")
    plotCompareCorr( quantResid2, peakIDs, metaData$Dx) + ggtitle(main)
  dev.off()


  # Examining differential accessability signal for these peaks
  topTable(fitDE, coef='DxSCZ', number=Inf)[peakIDs,]
  df_results = combineResults( ecdBox, clstScore, treeListClusters, peakLocs2)
  head(df_results)

  #444444444 Trying to pull my own data, using clst size 100:
  df = summary( ecdBox )
  df = head(df[df$id == 100,])
  peakIDs = getFeaturesInCluster( treeListClusters_collapse, df$chrom[1], df$cluster[1], df$id[1])
  query = range(peakLocs2[names(peakLocs2) %in% peakIDs])
  window = 1e5
  start(query) = start(query) - window
  end(query) = end(query) + window

  pdf("/sc/arion/projects/epigenAD/Bukola/significantClustersATAC100.pdf")
    plotDecorate( ensdb, treeList, treeListClusters_collapse, peakLocs2, query, showGenes=TRUE)
  dev.off()

  # Now to compare the level of correlation in this cluster for cases vs control.
  main = paste0(df$chrom[1], ': cluster ', df$cluster[1])
  pdf("/sc/arion/projects/epigenAD/Bukola/plotComparingCorrATAC100.pdf")
    plotCompareCorr( quantResid2, peakIDs, metaData$Dx) + ggtitle(main)
  dev.off()

  #444444444 Trying to pull my own data, using clst size 50:
  df = summary( ecdBox )
  df = head(df[df$id == 50,])
  peakIDs = getFeaturesInCluster( treeListClusters_collapse, df$chrom[1], df$cluster[1], df$id[1])
  query = range(peakLocs2[names(peakLocs2) %in% peakIDs])
  window = 1e5
  start(query) = start(query) - window
  end(query) = end(query) + window

  pdf("/sc/arion/projects/epigenAD/Bukola/significantClustersATAC50.pdf")
    plotDecorate( ensdb, treeList, treeListClusters_collapse, peakLocs2, query, showGenes=TRUE)
  dev.off()

  # Now to compare the level of correlation in this cluster for cases vs control.
  main = paste0(df$chrom[1], ': cluster ', df$cluster[1])
  pdf("/sc/arion/projects/epigenAD/Bukola/plotComparingCorrATAC50.pdf")
    plotCompareCorr( quantResid2, peakIDs, metaData$Dx) + ggtitle(main)
  dev.off()

  #444444444 Trying to pull my own data, using clst size 25:
  df = summary( ecdBox )
  df = head(df[df$id == 25,])
  peakIDs = getFeaturesInCluster( treeListClusters_collapse, df$chrom[1], df$cluster[1], df$id[1])
  query = range(peakLocs2[names(peakLocs2) %in% peakIDs])
  window = 1e5
  start(query) = start(query) - window
  end(query) = end(query) + window

  pdf("/sc/arion/projects/epigenAD/Bukola/significantClustersATAC25.pdf")
    plotDecorate( ensdb, treeList, treeListClusters_collapse, peakLocs2, query, showGenes=TRUE)
  dev.off()

  # Now to compare the level of correlation in this cluster for cases vs control.
  main = paste0(df$chrom[1], ': cluster ', df$cluster[1])
  pdf("/sc/arion/projects/epigenAD/Bukola/plotComparingCorrATAC25.pdf")
    plotCompareCorr( quantResid2, peakIDs, metaData$Dx) + ggtitle(main)
  dev.off()

  #444444444 Trying to pull my own data, using clst size 10:
  df = summary( ecdBox )
  df = head(df[df$id == 10,])
  peakIDs = getFeaturesInCluster( treeListClusters_collapse, df$chrom[1], df$cluster[1], df$id[1])
  query = range(peakLocs2[names(peakLocs2) %in% peakIDs])
  window = 1e5
  start(query) = start(query) - window
  end(query) = end(query) + window

  pdf("/sc/arion/projects/epigenAD/Bukola/significantClustersATAC10.pdf")
    plotDecorate( ensdb, treeList, treeListClusters_collapse, peakLocs2, query, showGenes=TRUE)
  dev.off()

  # Now to compare the level of correlation in this cluster for cases vs control.
  main = paste0(df$chrom[1], ': cluster ', df$cluster[1])
  pdf("/sc/arion/projects/epigenAD/Bukola/plotComparingCorrATAC10.pdf")
    plotCompareCorr( quantResid2, peakIDs, metaData$Dx) + ggtitle(main)
  dev.off()



  ###
  # peakIDs = getFeaturesInCluster( treeListClusters_collapse, df_fdr$chrom[4], df_fdr$cluster[4], df_fdr$id[4])

  # temp=as.numeric(as.character(df_results$pValue))
  # pdf("/sc/arion/projects/epigenAD/Bukola/anothertest.pdf")
  #   plotDecorate( ensdb, treeList, treeListClusters_collapse, peakLocs2, query, showGenes=TRUE)
  # dev.off()

  # pdf("/sc/arion/projects/epigenAD/Bukola/histogramm.pdf")
  #   hist(temp, 100)
  # dev.off()
