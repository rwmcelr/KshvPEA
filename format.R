## Library calls/ set script variables -------------------------------------
setwd("/Users/rmcelroy/Desktop/PEA/")
library("clusterProfiler")
library("pathview")

WtMut <- read.csv("WtMut_PathwayData.csv")
wtMut <- WtMut[,c("entrez","log2FoldChange","pvalue")]
wtMut <- wtMut[complete.cases(wtMut), ]
wtMut <- wtMut[order(-wtMut$log2FoldChange),]

wtMut_geneList <- wtMut$log2FoldChange
names(wtMut_geneList) <- wtMut$entrez

wtMut_sig <- subset(wtMut, pvalue < 0.05)
wtMut_genes <- wtMut_sig$log2FoldChange
names(wtMut_genes) <- wtMut_sig$entrez
genes <- names(wtMut_genes)

go_enrich <- enrichGO(gene = genes,
                      universe = names(wtMut_geneList),
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'ENTREZID',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05)

kk <- enrichKEGG(gene = genes,
                 organism = 'hsa',
                 pvalueCutoff = 0.2)

wtMut_geneList <- wtMut_geneList[!duplicated(names(wtMut_geneList))]