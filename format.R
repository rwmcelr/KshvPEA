## Library calls/ set script variables -------------------------------------
setwd("/Users/rmcelroy/Desktop/PEA/")
library("clusterProfiler")
library("pathview")
library("ggplot2")
library("enrichplot")
library("ReactomePA")

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

## GO (Gene Ontology) Analysis --------------------------------------------
# Over-Representation Analysis
# BP = Biological Process, MF = Molecular Function
go_enrichBP <- enrichGO(gene = genes,
                      universe = names(wtMut_geneList),
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'ENTREZID',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.2)
write.csv(go_enrichBP, file="Wt Vs Mut GO BP enrichment.csv")
dotplot(go_enrichBP) + ggtitle("Wt Vs Mut GO BP Enrichment (ORA) Dotplot")
go_enrichBP <- pairwise_termsim(go_enrichBP)
p <- emapplot(go_enrichBP) + ggtitle("Wt Vs Mut GO BP Enrichment Map")
p
goplot(go_enrichBP)

go_enrichMF <- enrichGO(gene = genes,
                      universe = names(wtMut_geneList),
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'ENTREZID',
                      readable = T,
                      ont = "MF",
                      pvalueCutoff = 0.2)
write.csv(go_enrichMF, file="Wt Vs Mut GO MF enrichment.csv")
dotplot(go_enrichMF) + ggtitle("Wt Vs Mut GO MF Enrichment (ORA) Dotplot")
go_enrichMF <- pairwise_termsim(go_enrichMF)
p1 <- emapplot(go_enrichMF) + ggtitle("Wt Vs Mut GO MF Enrichment Map")
p1
goplot(go_enrichMF)

## KEGG Pathways Analysis --------------------------------------------
# Over-Representation Analysis
kk <- enrichKEGG(gene = genes,
                 organism = 'hsa',
                 pvalueCutoff = 0.2)
write.csv(kk, file="WtVsMut_KEGG_ORA.csv")
dotplot(kk) + ggtitle("Wt Vs Mut KEGG Dotplot")
kk2 <- pairwise_termsim(kk)
p2 <- emapplot(kk2) + ggtitle("Wt Vs Mut KEGG Enrichment Map")
p2

## WikiPathways Analysis --------------------------------------------
# Over-Representation Analysis
wiPath1 <- enrichWP(genes, organism = "Homo sapiens")
write.csv(wiPath1, file="WtVsMut_WikiPathways_ORA.csv")
dotplot(wiPath1) + ggtitle("Wt Vs Mut WikiPathways ORA Dotplot")
wiPath1 <- pairwise_termsim(wiPath1)
p3 <- emapplot(wiPath1) + ggtitle("Wt Vs Mut WikiPathways ORA Map")
p3

# Gene Set Enrichment Analysis
wiPath2 <- gseWP(wtMut_geneList, organism = "Homo sapiens")
write.csv(wiPath2, file="WtVsMut_WikiPathways_GSEA.csv")
dotplot(wiPath2) + ggtitle("Wt Vs Mut WikiPathways GSEA Dotplot")
wiPath2 <- pairwise_termsim(wiPath2)
p4 <- emapplot(wiPath2) + ggtitle("Wt Vs Mut WikiPathways GSEA Map")
p4

## Reactome Pathways Analysis --------------------------------------------
# Over-Representation Analysis
reactome1 <- enrichPathway(gene=genes, pvalueCutoff = 0.05, readable=TRUE)
write.csv(reactome1, file="Wt Vs Mut ReactomePathways ORA.csv")
dotplot(reactome1) + ggtitle("Wt Vs Mut Reactome ORA Dotplot")
reactome1 <- pairwise_termsim(reactome1)
p5 <- emapplot(reactome1) + ggtitle("Wt Vs Mut ReactomePathways ORA Map")
p5

# Gene Set Enrichment Analysis
reactome2 <- gsePathway(wtMut_geneList,
                        pvalueCutoff = 0.2,
                        pAdjustMethod = "BH",
                        verbose = FALSE)
write.csv(reactome2, file="Wt Vs Mut ReactomePathways GSEA.csv")
dotplot(reactome2) + ggtitle("Wt Vs Mut Reactome GSEA Dotplot")
reactome2 <- pairwise_termsim(reactome2)
p6 <- emapplot(reactome2) + ggtitle("Wt Vs Mut ReactomePathways GSEA Map")
p6
