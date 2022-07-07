## Bukola Ajanaku
## August 3, 2021
## finalizing.R : finalized code for my poster plots
## on screen 1 as forFINALf
## module load udunits proj gdal geos
## module load R/4.0.2
## R
## save(list = ls(), file = "/sc/arion/work/ajanab01/newFinalPlayData.RDATA")

.libPaths(c("/hpc/users/ajanab01/.Rlib", .libPaths()))
# New paths: .libPaths(c("/sc/arion/projects/psychgen/scratch/temp_from_kiran/RLib",.libPaths()))

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
  # CONTROL
    treeListCON = runOrderedClusteringGenome( quantResid2[,metaData$Dx=='Control'], peakLocs2)
    treeListOriginalCON = runOrderedClusteringGenome( vobj$E[,metaData$Dx=='Control'], peakLocs2)
    treeListClustersCON = createClusters( treeListCON, method='meanClusterSize', meanClusterSize=c(10, 12, 15, 17, 20, 25, 30, 35, 40, 45, 50, 75, 100) )
    n_clustersCON = countClusters(treeListClustersCON)
    clstScoreCON = scoreClusters(treeListCON, treeListClustersCON, BPPARAM=SerialParam())
  # DISEASE
    treeListSCZ = runOrderedClusteringGenome( quantResid2[,metaData$Dx=='SCZ'], peakLocs2)
    treeListOriginalSCZ = runOrderedClusteringGenome( vobj$E[,metaData$Dx=='SCZ'], peakLocs2)
    treeListClustersSCZ = createClusters( treeListSCZ, method='meanClusterSize', meanClusterSize=c(10, 12, 15, 17, 20, 25, 30, 35, 40, 45, 50, 75, 100) )
    n_clustersSCZ = countClusters(treeListClustersSCZ)
    clstScoreSCZ = scoreClusters(treeListSCZ, treeListClustersSCZ, BPPARAM=SerialParam())
  # ALL
    treeList = runOrderedClusteringGenome( quantResid2, peakLocs2)
    treeListOriginal = runOrderedClusteringGenome( vobj$E, peakLocs2)
    treeListClusters = createClusters( treeList, method='meanClusterSize', meanClusterSize=c(10, 12, 15, 17, 20, 25, 30, 35, 40, 45, 50, 75, 100))
    n_clusters = countClusters( treeListClusters )
    clstScore = scoreClusters(treeList, treeListClusters, BPPARAM=SerialParam() )


# Now, dropping all the N = 1s
  # CONTROL
    originalClstScoreCON = clstScoreCON
    # clearing out all the clusters with only 1 feature
    cleanclst10CON = clstScoreCON[["10"]][clstScoreCON[["10"]]$N != 1, ]
    cleanclst12CON = clstScoreCON[["12"]][clstScoreCON[["12"]]$N != 1, ]
    cleanclst15CON = clstScoreCON[["15"]][clstScoreCON[["15"]]$N != 1, ]
    cleanclst17CON = clstScoreCON[["17"]][clstScoreCON[["17"]]$N != 1, ]
    cleanclst20CON = clstScoreCON[["20"]][clstScoreCON[["20"]]$N != 1, ]
    cleanclst25CON = clstScoreCON[["25"]][clstScoreCON[["25"]]$N != 1, ]
    cleanclst30CON = clstScoreCON[["30"]][clstScoreCON[["30"]]$N != 1, ]
    cleanclst35CON = clstScoreCON[["35"]][clstScoreCON[["35"]]$N != 1, ]
    cleanclst40CON = clstScoreCON[["40"]][clstScoreCON[["40"]]$N != 1, ]
    cleanclst45CON = clstScoreCON[["45"]][clstScoreCON[["45"]]$N != 1, ]
    cleanclst50CON = clstScoreCON[["50"]][clstScoreCON[["50"]]$N != 1, ]
    cleanclst75CON = clstScoreCON[["75"]][clstScoreCON[["75"]]$N != 1, ]
    cleanclst100CON = clstScoreCON[["100"]][clstScoreCON[["100"]]$N != 1, ]
    clstScoreCON = list(cleanclst10CON, cleanclst12CON, cleanclst15CON, cleanclst17CON,
      cleanclst20CON, cleanclst25CON, cleanclst30CON, cleanclst35CON, cleanclst40CON,
      cleanclst45CON, cleanclst50CON, cleanclst75CON, cleanclst100CON)
    names(clstScoreCON) = c("10", "12", "15", "17", "20", "25", "30", "35", "40", "45", "50", "75", "100")

  # DISEASE
    originalClstScoreSCZ = clstScoreSCZ
    # clearing out all the clusters with only 1 feature
    cleanclst10SCZ = clstScoreSCZ[["10"]][clstScoreSCZ[["10"]]$N != 1, ]
    cleanclst12SCZ = clstScoreSCZ[["12"]][clstScoreSCZ[["12"]]$N != 1, ]
    cleanclst15SCZ = clstScoreSCZ[["15"]][clstScoreSCZ[["15"]]$N != 1, ]
    cleanclst17SCZ = clstScoreSCZ[["17"]][clstScoreSCZ[["17"]]$N != 1, ]
    cleanclst20SCZ = clstScoreSCZ[["20"]][clstScoreSCZ[["20"]]$N != 1, ]
    cleanclst25SCZ = clstScoreSCZ[["25"]][clstScoreSCZ[["25"]]$N != 1, ]
    cleanclst30SCZ = clstScoreSCZ[["30"]][clstScoreSCZ[["30"]]$N != 1, ]
    cleanclst35SCZ = clstScoreSCZ[["35"]][clstScoreSCZ[["35"]]$N != 1, ]
    cleanclst40SCZ = clstScoreSCZ[["40"]][clstScoreSCZ[["40"]]$N != 1, ]
    cleanclst45SCZ = clstScoreSCZ[["45"]][clstScoreSCZ[["45"]]$N != 1, ]
    cleanclst50SCZ = clstScoreSCZ[["50"]][clstScoreSCZ[["50"]]$N != 1, ]
    cleanclst75SCZ = clstScoreSCZ[["75"]][clstScoreSCZ[["75"]]$N != 1, ]
    cleanclst100SCZ = clstScoreSCZ[["100"]][clstScoreSCZ[["100"]]$N != 1, ]
    clstScoreSCZ = list(cleanclst10SCZ, cleanclst12SCZ, cleanclst15SCZ, cleanclst17SCZ,
      cleanclst20SCZ, cleanclst25SCZ, cleanclst30SCZ, cleanclst35SCZ, cleanclst40SCZ,
      cleanclst45SCZ, cleanclst50SCZ, cleanclst75SCZ, cleanclst100SCZ)
    names(clstScoreSCZ) = c("10", "12", "15", "17", "20", "25", "30", "35", "40", "45", "50", "75", "100")

  # ALL
    originalClstScore = clstScore
    # clearing out all the clusters with only 1 feature
    cleanclst10 = clstScore[["10"]][clstScore[["10"]]$N != 1, ]
    cleanclst12 = clstScore[["12"]][clstScore[["12"]]$N != 1, ]
    cleanclst15 = clstScore[["15"]][clstScore[["15"]]$N != 1, ]
    cleanclst17 = clstScore[["17"]][clstScore[["17"]]$N != 1, ]
    cleanclst20 = clstScore[["20"]][clstScore[["20"]]$N != 1, ]
    cleanclst25 = clstScore[["25"]][clstScore[["25"]]$N != 1, ]
    cleanclst30 = clstScore[["30"]][clstScore[["30"]]$N != 1, ]
    cleanclst35 = clstScore[["35"]][clstScore[["35"]]$N != 1, ]
    cleanclst40 = clstScore[["40"]][clstScore[["40"]]$N != 1, ]
    cleanclst45 = clstScore[["45"]][clstScore[["45"]]$N != 1, ]
    cleanclst50 = clstScore[["50"]][clstScore[["50"]]$N != 1, ]
    cleanclst75 = clstScore[["75"]][clstScore[["75"]]$N != 1, ]
    cleanclst100 = clstScore[["100"]][clstScore[["100"]]$N != 1, ]
    clstScore = list(cleanclst10, cleanclst12, cleanclst15, cleanclst17, cleanclst20,
      cleanclst25, cleanclst30, cleanclst35, cleanclst40, cleanclst45,
      cleanclst50, cleanclst75, cleanclst100)
    names(clstScore) = c("10", "12", "15", "17", "20", "25", "30", "35", "40", "45", "50", "75", "100")

# Retaining clusters based on strength and getting these retained clusters. Then collapsing similar clusters.
  # CONTROL
    clustIncludeCON = retainClusters( clstScoreCON, "LEF", 0.1 )
    treeListClusters_filterCON = filterClusters( treeListClustersCON, clustIncludeCON )
    treeListClusters_collapseCON = collapseClusters( treeListClusters_filterCON, peakLocs2 )
    n_clustersCON = countClusters( treeListClusters_collapseCON )

  # DISEASE
    clustIncludeSCZ = retainClusters( clstScoreSCZ, "LEF", 0.1 )
    treeListClusters_filterSCZ = filterClusters( treeListClustersSCZ, clustIncludeSCZ )
    treeListClusters_collapseSCZ = collapseClusters( treeListClusters_filterSCZ, peakLocs2 )
    n_clustersSCZ = countClusters( treeListClusters_collapseSCZ )

  # ALL
    clustInclude = retainClusters( clstScore, "LEF", 0.1 )
    treeListClusters_filter = filterClusters( treeListClusters, clustInclude )
    treeListClusters_collapse = collapseClusters( treeListClusters_filter, peakLocs2 )
    n_clusters = countClusters( treeListClusters_collapse )

### NOW FOR CRD DATA:

# starting the CTCF code here:
  CTCF = read.table("/sc/arion/projects/epigenAD/Bukola/h1_neuron_ENCFF372JOV.txt")
  colnames(CTCF) = c("chr","start","end")
  CTCF$GeneID = paste0("peak_",1:dim(CTCF)[1])
  CTCF = as(CTCF,"GRanges")
  CTCF = sort(CTCF)
​
  newCRD10 = peakLocs
  newCRD10$cluster <- NA
  newCRD12 = peakLocs
  newCRD12$cluster <- NA
  newCRD15 = peakLocs
  newCRD15$cluster <- NA
  newCRD17 = peakLocs
  newCRD17$cluster <- NA
  newCRD20 = peakLocs
  newCRD20$cluster <- NA
  newCRD25 = peakLocs
  newCRD25$cluster <- NA
  newCRD30 = peakLocs
  newCRD30$cluster <- NA
  newCRD35 = peakLocs
  newCRD35$cluster <- NA
  newCRD40 = peakLocs
  newCRD40$cluster <- NA
  newCRD45 = peakLocs
  newCRD45$cluster <- NA
  newCRD50 = peakLocs
  newCRD50$cluster <- NA
  newCRD75 = peakLocs
  newCRD75$cluster <- NA
  newCRD100 = peakLocs
  newCRD100$cluster <- NA

group <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
  "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
  "chr19", "chr20", "chr21", "chr22")

# ADDING THE METHYLATION NAMES PER CLST SIZE (FOR BOTH DISEASE AND CONTROL: ALL)
  # 10:
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

  # 12:
  for (m in group){
      df12 = as.data.frame(treeListClusters_filter[["12"]][[m]])
      df12$ID=rownames(df12)
      names(df12) <- c("cluster", "ID")
      newCRD12 = as.data.frame(peakLocs)
      newCRD12 <- newCRD12[newCRD12$seqnames == m,]
      newCRD12$ID <- rownames(newCRD12)
      assign(paste("forc12and", m, sep = ""), merge(newCRD12, df12, by = "ID"))
  }

  finCRD12 = rbind(forc12andchr1, forc12andchr2, forc12andchr3, forc12andchr4, forc12andchr5,
    forc12andchr6, forc12andchr7, forc12andchr8, forc12andchr9, forc12andchr10, forc12andchr11,
    forc12andchr12, forc12andchr13, forc12andchr14, forc12andchr15, forc12andchr16, forc12andchr17,
    forc12andchr18, forc12andchr19, forc12andchr20, forc12andchr21, forc12andchr22)

  # 15:
  for (m in group){
      df15 = as.data.frame(treeListClusters_filter[["15"]][[m]])
      df15$ID=rownames(df15)
      names(df15) <- c("cluster", "ID")
      newCRD15 = as.data.frame(peakLocs)
      newCRD15 <- newCRD15[newCRD15$seqnames == m,]
      newCRD15$ID <- rownames(newCRD15)
      assign(paste("forc15and", m, sep = ""), merge(newCRD15, df15, by = "ID"))
  }

  finCRD15 = rbind(forc15andchr1, forc15andchr2, forc15andchr3, forc15andchr4, forc15andchr5,
    forc15andchr6, forc15andchr7, forc15andchr8, forc15andchr9, forc15andchr10, forc15andchr11,
    forc15andchr12, forc15andchr13, forc15andchr14, forc15andchr15, forc15andchr16, forc15andchr17,
    forc15andchr18, forc15andchr19, forc15andchr20, forc15andchr21, forc15andchr22)

  # 17:
  for (m in group){
      df17 = as.data.frame(treeListClusters_filter[["17"]][[m]])
      df17$ID=rownames(df17)
      names(df17) <- c("cluster", "ID")
      newCRD17 = as.data.frame(peakLocs)
      newCRD17 <- newCRD17[newCRD17$seqnames == m,]
      newCRD17$ID <- rownames(newCRD17)
      assign(paste("forc17and", m, sep = ""), merge(newCRD17, df17, by = "ID"))
  }

  finCRD17 = rbind(forc17andchr1, forc17andchr2, forc17andchr3, forc17andchr4, forc17andchr5,
    forc17andchr6, forc17andchr7, forc17andchr8, forc17andchr9, forc17andchr10, forc17andchr11,
    forc17andchr12, forc17andchr13, forc17andchr14, forc17andchr15, forc17andchr16, forc17andchr17,
    forc17andchr18, forc17andchr19, forc17andchr20, forc17andchr21, forc17andchr22)

  # 20:
  for (m in group){
    df20 = as.data.frame(treeListClusters_filter[["20"]][[m]])
    df20$ID=rownames(df20)
    names(df20) <- c("cluster", "ID")
    newCRD20 = as.data.frame(peakLocs)
    newCRD20 <- newCRD20[newCRD20$seqnames == m,]
    newCRD20$ID <- rownames(newCRD20)
    assign(paste("forc20and", m, sep = ""), merge(newCRD20, df20, by = "ID"))  }

  finCRD20 = rbind(forc20andchr1, forc20andchr2, forc20andchr3, forc20andchr4, forc20andchr5,
    forc20andchr6, forc20andchr7, forc20andchr8, forc20andchr9, forc20andchr10, forc20andchr11,
    forc20andchr12, forc20andchr13, forc20andchr14, forc20andchr15, forc20andchr16, forc20andchr17,
    forc20andchr18, forc20andchr19, forc20andchr20, forc20andchr21, forc20andchr22)

  # 25:
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

  # 30:
  for (m in group){
    df30 = as.data.frame(treeListClusters_filter[["30"]][[m]])
    df30$ID=rownames(df30)
    names(df30) <- c("cluster", "ID")
    newCRD30 = as.data.frame(peakLocs)
    newCRD30 <- newCRD30[newCRD30$seqnames == m,]
    newCRD30$ID <- rownames(newCRD30)
    assign(paste("forc30and", m, sep = ""), merge(newCRD30, df30, by = "ID"))  }

  finCRD30 = rbind(forc30andchr1, forc30andchr2, forc30andchr3, forc30andchr4, forc30andchr5,
    forc30andchr6, forc30andchr7, forc30andchr8, forc30andchr9, forc30andchr10, forc30andchr11,
    forc30andchr12, forc30andchr13, forc30andchr14, forc30andchr15, forc30andchr16, forc30andchr17,
    forc30andchr18, forc30andchr19, forc30andchr20, forc30andchr21, forc30andchr22)

  # 35:
  for (m in group){
    df35 = as.data.frame(treeListClusters_filter[["35"]][[m]])
    df35$ID=rownames(df35)
    names(df35) <- c("cluster", "ID")
    newCRD35 = as.data.frame(peakLocs)
    newCRD35 <- newCRD35[newCRD35$seqnames == m,]
    newCRD35$ID <- rownames(newCRD35)
    assign(paste("forc35and", m, sep = ""), merge(newCRD35, df35, by = "ID"))  }

  finCRD35 = rbind(forc35andchr1, forc35andchr2, forc35andchr3, forc35andchr4, forc35andchr5,
    forc35andchr6, forc35andchr7, forc35andchr8, forc35andchr9, forc35andchr10, forc35andchr11,
    forc35andchr12, forc35andchr13, forc35andchr14, forc35andchr15, forc35andchr16, forc35andchr17,
    forc35andchr18, forc35andchr19, forc35andchr20, forc35andchr21, forc35andchr22)

  # 40:
  for (m in group){
    df40 = as.data.frame(treeListClusters_filter[["40"]][[m]])
    df40$ID=rownames(df40)
    names(df40) <- c("cluster", "ID")
    newCRD40 = as.data.frame(peakLocs)
    newCRD40 <- newCRD40[newCRD40$seqnames == m,]
    newCRD40$ID <- rownames(newCRD40)
    assign(paste("forc40and", m, sep = ""), merge(newCRD40, df40, by = "ID"))  }

  finCRD40 = rbind(forc40andchr1, forc40andchr2, forc40andchr3, forc40andchr4, forc40andchr5,
    forc40andchr6, forc40andchr7, forc40andchr8, forc40andchr9, forc40andchr10, forc40andchr11,
    forc40andchr12, forc40andchr13, forc40andchr14, forc40andchr15, forc40andchr16, forc40andchr17,
    forc40andchr18, forc40andchr19, forc40andchr20, forc40andchr21, forc40andchr22)

  # 45:
  for (m in group){
    df45 = as.data.frame(treeListClusters_filter[["45"]][[m]])
    df45$ID=rownames(df45)
    names(df45) <- c("cluster", "ID")
    newCRD45 = as.data.frame(peakLocs)
    newCRD45 <- newCRD45[newCRD45$seqnames == m,]
    newCRD45$ID <- rownames(newCRD45)
    assign(paste("forc45and", m, sep = ""), merge(newCRD45, df45, by = "ID"))  }

  finCRD45 = rbind(forc45andchr1, forc45andchr2, forc45andchr3, forc45andchr4, forc45andchr5,
    forc45andchr6, forc45andchr7, forc45andchr8, forc45andchr9, forc45andchr10, forc45andchr11,
    forc45andchr12, forc45andchr13, forc45andchr14, forc45andchr15, forc45andchr16, forc45andchr17,
    forc45andchr18, forc45andchr19, forc45andchr20, forc45andchr21, forc45andchr22)

  # 50:
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

  # 75:
  for (m in group){
    df75 = as.data.frame(treeListClusters_filter[["75"]][[m]])
    df75$ID=rownames(df75)
    names(df75) <- c("cluster", "ID")
    newCRD75 = as.data.frame(peakLocs)
    newCRD75 <- newCRD75[newCRD75$seqnames == m,]
    newCRD75$ID <- rownames(newCRD75)
    assign(paste("forc75and", m, sep = ""), merge(newCRD75, df75, by = "ID"))  }

  finCRD75 = rbind(forc75andchr1, forc75andchr2, forc75andchr3, forc75andchr4, forc75andchr5,
    forc75andchr6, forc75andchr7, forc75andchr8, forc75andchr9, forc75andchr10, forc75andchr11,
    forc75andchr12, forc75andchr13, forc75andchr14, forc75andchr15, forc75andchr16, forc75andchr17,
    forc75andchr18, forc75andchr19, forc75andchr20, forc75andchr21, forc75andchr22)

  # 100:
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
  finCRD12 <- as(finCRD12, "GRanges")
  finCRD15 <- as(finCRD15, "GRanges")
  finCRD17 <- as(finCRD17, "GRanges")
  finCRD20 <- as(finCRD20, "GRanges")
  finCRD25 <- as(finCRD25, "GRanges")
  finCRD30 <- as(finCRD30, "GRanges")
  finCRD35 <- as(finCRD35, "GRanges")
  finCRD40 <- as(finCRD40, "GRanges")
  finCRD45 <- as(finCRD45, "GRanges")
  finCRD50 <- as(finCRD50, "GRanges")
  finCRD75 <- as(finCRD75, "GRanges")
  finCRD100 <- as(finCRD100, "GRanges")

  # allCRDs = list(finCRD10, finCRD20, finCRD25, finCRD30, finCRD35, finCRD40, finCRD45,
  #   finCRD50, finCRD75, finCRD100)
  # names(allCRDs) = c("10", "20", "25", "30", "35", "40", "45", "50", "75", "100")

  allCRDs = list(finCRD10, finCRD12, finCRD15, finCRD17, finCRD20, finCRD25, finCRD30)
  names(allCRDs) = c("10", "12", "15", "17", "20", "25", "30")

      # MAIN CODE FOR CTCF Density:
      CTCF_density = function(query_grange,CTCF,CTCF_dist_bins,dist_window,direction){
      dist_window = dist_window
      CTCF_dist_bins = CTCF_dist_bins
      counts = matrix(0,nrow=dist_window)
      dist_crds = matrix(0,nrow=dist_window)
      dist_crds[1] = 0
      query_bins_list = list()

      if (direction == "left"){
        df1 = as(data.frame("chr"=seqnames(query_grange),"end"=start(query_grange),"start"=start(query_grange)-CTCF_dist_bins),"GRanges")
        counts[1] = uniqueN(as.data.frame(findOverlaps(CTCF,df1))[,1])/length(CTCF)
        } else {
          df1 = as(data.frame("chr"=seqnames(query_grange),"start"=end(query_grange),"end"=end(query_grange)+CTCF_dist_bins),"GRanges")
          counts[1] = uniqueN(as.data.frame(findOverlaps(CTCF,df1))[,1])/length(CTCF)
          }

      query_bins_list[[1]] = df1

      for (i in (2:dist_window)) {
         if (direction == "left"){
            df1 = as(data.frame("chr"=seqnames(df1),"end"=start(df1),"start"=start(df1)-CTCF_dist_bins),"GRanges")
            end(df1) = end(df1)-1
            counts[i] = uniqueN(as.data.frame(findOverlaps(CTCF,df1))[,1])/length(CTCF)
            dist_crds[i] = -1*i*CTCF_dist_bins
            } else {
              df1 = as(data.frame("chr"=seqnames(df1),"start"=end(df1),"end"=end(df1)+CTCF_dist_bins),"GRanges")
              start(df1) = start(df1)+1
              counts[i] = uniqueN(as.data.frame(findOverlaps(CTCF,df1))[,1])/length(CTCF)
              dist_crds[i] = 1*i*CTCF_dist_bins
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

  names(CRD_CTCF_density) = names(allCRDs)
  CRD_CTCF_density_df = reshape2::melt(CRD_CTCF_density,id=c("Density_CRD","Distance"))
  CRD_CTCF_density_df$L1 = factor(CRD_CTCF_density_df$L1, levels = c("10", "12", "15", "17", "20", "25", "30"))

  pdf("/sc/arion/projects/epigenAD/Bukola/FinalPlots/finalATACCRDDensity.pdf", height = 6, width = 8)
    ggplot(CRD_CTCF_density_df, aes(Distance, Density_CRD, color = L1)) + geom_point() + geom_line() + theme_classic() +
      scale_colour_discrete("Mean Cluster Size") +
      ggtitle("CTCF Protein Enrichments for Cis-Regulatory Domains (CRDs)") + xlab("Distance from CRD Border") + ylab("CTCF Protein Density") +
      theme(plot.title = element_text(size=16, hjust = 0.25), axis.title.x = element_text(size= 14), axis.text=element_text(size=12),
      axis.title.y = element_text(size= 14))
  dev.off()

# Making Scatter Plot for Cluster Size : 10 --------------------------------------------------------
 # MAC vs LEF
      clstScore[["10"]]$quantGroup <- NA
       for(i in 1:length(clstScore[["10"]]$N)){
        if(clstScore[["10"]]$N[i] >= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[1]] && clstScore[["10"]]$N[i] <= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[2]]){
          clstScore[["10"]][i,]$quantGroup <- "Only 2"
        } else if(clstScore[["10"]]$N[i] >= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[2]] && clstScore[["10"]]$N[i] <= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[3]]){
            clstScore[["10"]][i,]$quantGroup <- "4-3"
        } else if(clstScore[["10"]]$N[i] >= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[3]] && clstScore[["10"]]$N[i] <= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[4]]){
            clstScore[["10"]][i,]$quantGroup <- "7-5"
        } else if(clstScore[["10"]]$N[i] >= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[4]] && clstScore[["10"]]$N[i] <= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[5]]){
            clstScore[["10"]][i,]$quantGroup <- "15-8"
        } else if(clstScore[["10"]]$N[i] >= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[5]] && clstScore[["10"]]$N[i] <= quantile(clstScore[["10"]]$N, probs = seq(0,1,0.20))[[6]]){
            clstScore[["10"]][i,]$quantGroup <- "243-16"
        }
       }

      clstScore[["10"]]$quantGroup = factor(clstScore[["10"]]$quantGroup, levels = c("243-16", "15-8", "7-5", "4-3", "Only 2"))

      pdf("/sc/arion/projects/epigenAD/Bukola/FinalPlots/finalATACviewScatterPlots.pdf", width = 8, height = 8)
       ggplot(clstScore[["10"]], aes(clstScore[["10"]]$mean_abs_corr, clstScore[["10"]]$LEF, color = clstScore[["10"]]$quantGroup)) + geom_point() + labs(title = "For Cluster Size 10",
           x = "MAC", y = "LEF") + theme_classic() + labs(color = "Number of Features") + ylim(0, 0.8) +
           ggtitle("Correlation Trends for Mean Cluster Size 10") + xlab("Lead Eigenvalue") + ylab("Mean Absolute Correlation") +
           theme(
           plot.title = element_text(size=20, hjust = 0.5),
           axis.title.x = element_text(size= 18),
           axis.text=element_text(size=14),
           axis.title.y = element_text(size= 18),
           legend.position = "bottom",
           legend.text = element_text(size = 14),
           legend.title = element_text(size = 16))
      dev.off()

  # MAC vs Number of Features

    pdf("/sc/arion/projects/epigenAD/Bukola/FinalPlots/finalATACbyNviewScatterPlots.pdf")
      ggplot(clstScore[["10"]][["chr1"]], aes(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1",]$N,
        clstScore[["10"]][clstScore[["10"]]$chrom == "chr1",]$mean_abs_corr)) + geom_point() + labs(title = "For Cluster Size 10 on Chr1",
        x = "Number of Features", y = "MAC")
    dev.off()

# Ending Scatter Plots: --------------------------------------------------------
  # NOW for TREELIST PLOTS:
    # My function:
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
  # Showing first 2 Genes:
      peaks = as.data.frame(peakLocs2)
      peaks$PeaksID = rownames(peaks)
      peaks$names = peaks$PeakIDs
      peaks_GR = as(peaks, "GRanges")

      library(EnsDb.Hsapiens.v75)
      ensdb= EnsDb.Hsapiens.v75

# ** Need to find a nice treeList plot. I changed it to chr22, check for chr1.
      tenhighmachighleffirst2 <- gatherData("10", "chr22", as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr22", "mean_abs_corr"])[5]),
        as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr22", "LEF"])[5]), 1, 5)
      pdf("/sc/arion/projects/epigenAD/Bukola/FinalPlots/ATACclst10highMAChighLEFFIRST2.pdf")
        tenhighmachighleffirst2
      dev.off()

    # SCZ
    tenhighmachighleffirst2SCZ <- gatherData("10", "chr1", as.numeric(summary(clstScoreSCZ[["10"]][clstScoreSCZ[["10"]]$chrom == "chr1", "mean_abs_corr"])[5]),
      as.numeric(summary(clstScore[["10"]][clstScore[["10"]]$chrom == "chr1", "LEF"])[5]), 1, 2)
    pdf("/sc/arion/projects/epigenAD/Bukola/FinalPlots/SCZATACclst10highMAChighLEFFIRST2.pdf")
      tenhighmachighleffirst2SCZ
    dev.off()

# Test differential signals:
  n_clusters = countClusters( treeListClusters_collapse )
  param = SnowParam(6, "SOCK", progressbar=TRUE)
  ecdBox = evalDiffCorr( quantResid2, metaData$Dx, peakLocs2, treeListClusters_collapse, BPPARAM=param, method="Box.permute" )

  df = summary( ecdBox )
  peakIDs = getFeaturesInCluster( treeListClusters_collapse, df$chrom[1], df$cluster[1], df$id[1])
  query = range(peakLocs2[names(peakLocs2) %in% peakIDs])
  window = 0
  # 1e5
  start(query) = start(query) - window
  end(query) = end(query) + window

  pdf("/sc/arion/projects/epigenAD/Bukola/FinalPlots/significantClustersATAC2.pdf")
    plotDecorate( ensdb, treeList, treeListClusters_collapse, peakLocs2, query, showGenes=TRUE)
  dev.off()

# Now to compare the level of correlation in this cluster for cases vs control.
  main = paste0(df$chrom[1], ': cluster ', df$cluster[1])
  pdf("/sc/arion/projects/epigenAD/Bukola/FinalPlots/plotComparingCorrATAC.pdf")
    plotCompareCorr( quantResid2, peakIDs, metaData$Dx) + ggtitle(main)
  dev.off()

# Examine differential accessability signal for these peaks
topTable(fitDE, coef='DxSCZ', number=Inf)[peakIDs,]
df_results = combineResults( ecdBox, clstScore, treeListClusters, peakLocs2)
head(df_results)
