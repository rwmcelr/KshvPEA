## Library calls/ set script variables -------------------------------------
setwd("/Users/rmcelroy/Desktop/PEA/")
library("clusterProfiler")
library("DOSE")
library("pathview")

wtMut <- read.csv("wtMut.csv")
wmGeneList <- wtMut[,2]
names(wmGeneList) <- as.character(wtMut[,1])
geneList = sort(wmGeneList, decreasing = TRUE)
kegg_geneList = sort(wmGeneList, decreasing = TRUE)

data(wmgeneList, package="DOSE")
wmGene <- names(wmgeneList)[abs(wmgeneList) > 1]
kk <- enrichKEGG(gene = wmGene,
                 universe = names(kegg_geneList),
                 organism = 'hsa',
                 pvalueCutoff = 0.05)
browseKEGG(kk, 'hsa04110')

hsa04110 <- pathview(gene.data  = wmgeneList,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(wmGene=max(abs(geneList)), cpd=1))

kk2 <- gseKEGG(geneList     = kegg_geneList,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
# -------------------------------------------------------------------------

# wtVec <- read.csv("wtVec.csv")
# wvGeneList <- wtVec[,2]
# names(wvGeneList) <- as.character(wtVec[,1])
# wvGeneList <- sort(wvGeneList, decreasing = TRUE)
# 
# mutVec <- read.csv("mutVec.csv")
# mvGeneList <- mutVec[,2]
# names(mvGeneList) <- as.character(mutVec[,1])
# mvGeneList <- sort(mvGeneList, decreasing = TRUE)

