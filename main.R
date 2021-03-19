## Library calls/ set script variables -------------------------------------
setwd("/Users/rmcelroy/Desktop/PEA/")
source("/Users/rmcelroy/Desktop/PEA/R/format.R")
source("/Users/rmcelroy/Desktop/PEA/R/PEA_ORA.R")
source("/Users/rmcelroy/Desktop/PEA/R/PEA_GSEA.R")

## Main Body --------------------------------------------------------------
wtMut <- read.csv("WtMut_PathwayData.csv")
wtMut_geneList <- formatGeneList(wtMut)
wtMut_genes <- formatGenes(wtMut)
PEA_ORA("Wt_Vs_Mutant", wtMut_geneList, wtMut_genes)

wtVec <- read.csv("WtVec_PathwayData.csv")
wtVec_geneList <- formatGeneList(wtVec)
wtVec_genes <- formatGenes(wtVec)
PEA_ORA("Wt_Vs_Vector", wtVec_geneList, wtVec_genes)

MutVec <- read.csv("MutVec_PathwayData.csv")
MutVec_geneList <- formatGeneList(MutVec)
MutVec_genes <- formatGenes(MutVec)
PEA_ORA("Mutant_Vs_Vector", MutVec_geneList, MutVec_genes)