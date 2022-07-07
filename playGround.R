## Bukola Ajanaku
## June 29, 2021
## playGround.R - playing around with the epigenAD tutorial and challenging
## Gabriel Hoffman's process. Special focus on understand treeList!
## Notes -----------------------------------------------------------------------
## Correlation: This is important for the 3D structure of the genome because the
## structural organization of the genome's chromatin is dependent on the correlation
## of the epigenetic histone peaks. Coordinated structures of histone peaks
## (since this is genomic data: specifically when adjacent to one anther) indicates
## regions of high correlation in these CRDs are areas of generally similar expression
## due to the influence of the methylation in these regions.
## -
## Adjacent Correlation: is a more strict approach to merging features when
## clustering in which only features that are adjacent to one another are
## allowed to be clustered. Maintaining the ordinal integrity of the genetic data.
## -
## Mean absolute correlation: takes the absolute value of the correlations measured
## from the matrix produced through the adjacent correlation measure.
## -
## Q: What part of the brain is this? A: PFC
## -
## Start Up Instructions:
## Log into Minerva first! Then, run:
## module load udunits proj gdal geos
## module load R/4.0.2
## R
## For saving objects: save(list = ls(), file = "/sc/arion/work/ajanab01/newRunner.RDATA") (but first I need to save the list, check code on slack)
## pathname for bigger folder: /sc/arion/work/ajanab01/
## example for particular variables of interest: save(clustscore, file="SCZ_sclustscore.RDATA")
## Latest: runningCode.RDATA
## Set the location for R to find its packages. Especially since I'm using
## Minerva.
.libPaths(c("/hpc/users/ajanab01/.Rlib", .libPaths()))

## Load packages and supress their start-up messages. Just for no extra output.
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

## Sets the style of output for the tables that are made later in the code.
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

# Loads gene location information.
response = readRDS("/sc/arion/projects/epigenAD/Bukola/DNA_methylation/response.RDS" )
metadata = readRDS("/sc/arion/projects/epigenAD/Bukola/DNA_methylation/metadata.RDS")
featureLocation = readRDS("/sc/arion/projects/epigenAD/Bukola/DNA_methylation/featureLocation.RDS")

# features must have width of at least 2
end(featureLocation) = end(featureLocation) + 1

# Conditions must be factors
metadata$Dx = factor(metadata$Dx, c("Control", "Schizo"))

# Computing residuals ----------------------------------------------------------

# Creates a model matrix for the columns called from the metadata. It makes
# "dummy" variables (using relative coefficients) that just describes the
$ relative differences among the values. Note that it does this for all each
# sample based on the methylation markers that are recorded.
design = model.matrix(~ Dx + Age + Race + negControl_PC1 + negControl_PC2 +
  negControl_PC3 + negControl_PC4, metadata)

# Takes the sequenced data and fits the variables of the response matrix for each
# sample to match that of the design matrix. So each sample is fit into relative
# "dummy" variables relative to one another. These coefficients are attached to
# each methylation marker.
fit = lmFit(response, design)

# Now that each sample is properly weighed, the methylation markers are now able
# to be processed through the culmination of each sample. Here we are getting
# relative clusters (their graphable statistics) and their statistics based on
# the weighed data for each sample of each methylation maker.
fit = eBayes(fit)

# Residuals helps get the coefficients for the best line of fit for each sample's
# methylation marker. All of which are relative to one another. Fits it with the
# same response matrix using statistical data obtained from the fit matrix. Cleaned.
residValues = residuals( fit, response)

# Now that we have found the residual values for each sample's methylation markers,
# our goal is to find methylation markers that are highly correlated with one
# another (using statistical processing of fit through this limma package). The
# following code gives the top, most highly correlated methylation markers.
topTable(fit, coef='DxSchizo')

# Now it's time to look at the genomes. ----------------------------------------
# Starting with the Covariates corrected read counts matrix.

# Preserves sequential order of the genome and only clusters peaks of methylation
# markers of the CRDs for ONLY features that are adjacent to one another. Important
# for creating the 3D structure. This is for residualized data (cleaned response).
treeList = runOrderedClusteringGenome(residValues[,metadata$Dx=='Control'], featureLocation, method.corr="spearman")

# Repeat using the response (sequenced) data. Trying to pick out significant areas
# of clustering based on the coordinated structure of methylation peaks. Will later compare
# to that of the expected residual values (based on the weighted fit from earlier)
# aka the residValues. Does so for each chromosome.
treeListOriginal = runOrderedClusteringGenome(response[,metadata$Dx=='Control'], featureLocation, method.corr="spearman")

# Using clustering inforation of highly correlated markers, this line evaluates
# the decay of correlation versus distance between features. It is found that the
# closer genes are, the more correlated they are. Around 1 megabase, this
# correlation is restricted solely 1 megabase around a gene. Do not need to
# dispute.

dfDistOriginal = evaluateCorrDecay( treeListOriginal, featureLocation, "chr22" )

pdf("CorrDecayforChr22")
plotCorrDecay( dfDistOriginal, method="R", xlim=c(500, 1e6), outlierQuantile=1e-5 )
dev.off()

dfDist = evaluateCorrDecay( treeList, featureLocation, "chr22" )

pdf("CorrDecayforChr22part2")
plotCorrDecay( dfDist, method="R", xlim=c(500, 1e6), outlierQuantile=1e-5 )
dev.off()

# Now to make the genome 3D through grouping up clusters! ----------------------

# Makes clusters of the gene. I got from the above graphs that these clusters
# are from sizes 10-100. Residualized data.
treeListClusters = createClusters(treeList, method='meanClusterSize',
  meanClusterSize=c(10, 25, 50, 100))

# Calculates the number of clusters.
n_clusters = countClusters( treeListClusters )

# For bar graph:---------------------------------------

clusterCounter = c(10, 25, 50, 100)
numberClstersinClst = as.data.frame(clusterCounter)
numberClstersinClst$counting = NA

numberClstersinClst$clusterCounter = factor(numberClstersinClst$clusterCounter, levels = clusterCounter)

for (i in clusterCounter){
 numberClstersinClst[clusterCounter == i, ]$counting = countClusters(treeListClusters[names(treeListClusters) == i,])
}

pdf("plotforNumberCluster")
  plotforNumberCluster = ggplot(data= numberClstersinClst, aes(x= clusterCounter, y= counting )) +
    geom_bar(stat="identity")
  plotforNumberCluster
dev.off()

# End bar graph:----------------------------------------

#* Making for loop bar graph for number of clusters of varying cluster sizes:

# Will score clusters based on their correlation in the structure.
clstScore = scoreClusters(treeList, treeListClusters)

# Graphically presents cluster strength!
df_LEF = do.call('rbind', clstScore )
df_LEF$id = factor(df_LEF$id, sort(unique(df_LEF$id)))

pdf("clustersLEFplot.pdf")
  clustIdeaLEF = ggplot(df_LEF, aes(LEF, color=id)) + geom_density() + theme_bw(17) + theme(aspect.ratio=1)
  clustIdeaLEF
dev.off()
# xdg-open clustersLEFplot.pdf
# scp -r XX@minerva.hpc.mssm.edu:/hpc/users/minerva/plotname.pdf ~/Desktop/
#

# Filtering the clustScores based on the strength of correlation.
clustInclude = retainClusters( clstScore, "LEF", 0.05 )

# Gets "retained" clusters that meet the cut off value based on LEF. This means
# that we are looking for structures with specific markers where the correlation
# is the greatest. Since this is a matrix, we must also be sure to collapse the
# similar clusters (the repeats).
treeListClusters_filter = filterClusters( treeListClusters, clustInclude )

treeListClusters_collapse = collapseClusters( treeListClusters_filter, featureLocation, jaccardCutoff=0.9 ) #**

# Now, we are looking at a different signal.
# Here we are taking advantage of the collapse function in order to identify
# individual clusters of the highest correlation.
n_clusters = countClusters( treeListClusters_collapse )

# Differential Correlation function: identifies correlated epigenetic features
# and finds clusters of features that are differentially correlated between two
# or more subsets of the data (in this case schizo vs control)
ecdBox = evalDiffCorr( residValues, metadata$Dx, featureLocation,
  treeListClusters_collapse, npermute, method = "Box.permute", method.corr="spearman")

# Summarize results. Looks at pairs of features (schizo vs control) for each
# cluster size. This summary returns the top features, across cluster sizes,
# that have the greatest correlative difference between the schizo and control
# groups.
df = summary( ecdBox )
head(df)

# * # save(list = ls(), file = "savingData.RDATA")

# Merges info for each feature pair all into one dataframe.
df_results = combineResults( ecdBox, clstScore, treeListClusters, featureLocation)
head(df_results)

# Creates a histogram of LEF values (higher LEF ~ higher similarity between features)

pdf("plotLEFhist.pdf")
  ggplot(df_results, aes(LEF, fill=id)) + geom_histogram(alpha=0.7) +
    theme_bw(17) + xlim(0,1) +
    theme(aspect.ratio=1, legend.position="bottom",
      plot.title = element_text(hjust = 0.5)) +
    scale_fill_discrete(name = "Requested mean cluster size") +
    xlab("Lead eigenvalue fraction (LEF)") +
    ggtitle("Summarize LEF")
dev.off()

# Creates a histogram of mean absolute correlation values. Mean absolute correlation
# computes the correlations between all pairs of features for each class and then
# averages these values for all pairs and all classes. Simply helps correct any significant
# variation in the graphs.

pdf("plotMAChist.pdf")
ggplot(df_results, aes(mean_abs_corr, fill=id)) + geom_histogram(alpha=0.7) + theme_bw(17) + xlim(0,1) +
  theme(aspect.ratio=1, legend.position="bottom",
    plot.title = element_text(hjust = 0.5)) +
  scale_fill_discrete(name = "Requested mean cluster size") +
   xlab("Mean absolute correlation") +
   ggtitle("Summarize absolute correlation")
dev.off()

# My ideas: 1) We compare the areas of major difference between the LEF and
# MAC graphs or 2) We compare the lower regions of LEF/MAC where there is a
# peak among the smaller cluster sizes.

# Splitting up df_results by the mean cluster size (= id) and making it into a
# genomic range object.

newGRange10 = sort(unlist(makeGRangesListFromDataFrame(df_results[df_results$id == 10,])))
newGRange25 = sort(unlist(makeGRangesListFromDataFrame(df_results[df_results$id == 25,])))
newGRange50 = sort(unlist(makeGRangesListFromDataFrame(df_results[df_results$id == 50,])))
newGRange100 = sort(unlist(makeGRangesListFromDataFrame(df_results[df_results$id == 100,])))

# make plots of the correlation matrix using function plotDecorate
# plot for high lef/low mac & vice versa
# http://deepfigv.mssm.edu/img/software/decorate/decorate_example.html
# DONE: download oxfuse (double version for mac)

###### -------------------------------- Correlation Matrices with plotDecorate -

# Loads the gene locations which is in a library.
library(EnsDb.Hsapiens.v75)

ensdb = EnsDb.Hsapiens.v75

# 1) get grange object with a specified MAC score.
# 2) fig1 = plotDecorate( ensdb, treeList, treeListClusters, simLocation, query)
# what about specified LEF values?

# newFigure <- gatherData("50", "chr12", 0.17)

gatherData <- function(clusterSize, chrome, MACscore) {

  df <- clstScore[[clusterSize]][clstScore[[clusterSize]]]$chrom == chrome,]

  MACscoreLL <- MACscore - 0.05
  MACscoreHL <- MACscore + 0.05

  newdf <- df[df$mean_abs_corr >= MACscoreLL & df$mean_abs_corr <= MACscoreHL,]

  otherdf <- as.data.frame(treeListClusters[[clusterSize]][[chrome]])
  colnames(otherdf) <- "numberCluster"
  otherdf$ID <- rownames(otherdf)

  probe <- unique(newdf$cluster)
  otherdf <- otherdf[otherdf$numberCluster %in% probe,]

  simLocation <- featureLocation[otherdf$ID]
  query <- range(simLocation)
  print(identical(otherdf$ID, names(simLocation)))

  fig1 = plotDecorate(ensdb, treeList, treeListClusters, simLocation, query)

# simLocation$names = names(simLocation)
  fig1
}

fromOther <- unique(otherdf$numberCluster)
fromNew <- unique(newdf$cluster)


names(treeListClusters_collapse[["50"]][["chr12"]]) # features in TreeListClusters


looking <- as.data.frame(simLocation)
looking$features <- rownames(looking)

glob <- as.data.frame(names(treeListClusters[["50"]][["chr12"]]))
bleh <- as.data.frame(looking$features)

sum(glob$tester %in% bleh$tester) # == 0

treeList = runOrderedClusteringGenome(residValues, featureLocation, method.corr="spearman" )

--
look <- rownames(residValues)
lookingfeat <- as.data.frame(featureLocation)
lookingfeat$features <- rownames(lookingfeat)

glob <- as.data.frame(look)
bleh <- as.data.frame(lookingfeat$features)
names(bleh) <- "tester"

sum(glob$look %in% bleh$tester) # == 456513
---

gatherER <- function(clusterSize, chrome, MACscore) {
  df <- clstScore[[clusterSize]][clstScore[[clusterSize]]]$chrom == chrome,]
  otherdf <- as.data.frame(treeListClusters[[clusterSize]][[chrome]])

  trying <- as.data.frame(featureLocation)
  trying$ID <- rownames(trying)

  names(otherdf)[1] <- "cluster"
  otherdf$ID <- rownames(otherdf)

  merged <- merge(otherdf, trying, by = "ID")


  names(trying)[12] <- "id"
  forThings <- df[df$mean_abs_corr >= MACscoreLL & df$mean_abs_corr <= MACscoreHL,]

  merged <- df[df$id %in% unique(otherdf$cluster),]

  mergedd <- merged[merged$cluster %in% unique(forThings$cluster),]

  featureLocation$ID <- names(featureLocation)
  simLocation <- featureLocation[mergedd$ID]
  query <- range(simLocation)

  figure1 = plotDecorate(ensdb, treeList, treeListClusters, simLocation, query)

}
