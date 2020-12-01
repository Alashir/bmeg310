library(pheatmap)
library(ggplot2)
library(DESeq2)

# Import metadata
colData = read.csv("GSE37704_metadata.csv", row.names=1)
colData

# import as dataframe
countData = read.csv("GSE37704_featurecounts.csv",row.names=1)
# convert dataframe to matrix
countData = as.matrix(countData) 

countData = countData[,-1]
countData = countData[rowSums(countData)>1, ]
head(countData)

sampleDists = dist(t(countData),upper = TRUE)
sampleDists

annot_col = data.frame(colData$condition)
row.names(annot_col) <- rownames(colData)

sampleDistMatrix = as.matrix( sampleDists )
rownames(sampleDistMatrix) = colnames(countData)
colnames(sampleDistMatrix) = colnames(countData)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE,
         annotation_col=annot_col)

pca_res <- prcomp(t(countData), scale. = TRUE)
score <- pca_res$x

score = as.data.frame(score)
score$color <- as.factor(colData$condition)


plt <- ggplot(score, aes(x=PC1, y=PC2,  color=color)) + geom_point(size = 4)
plt

dds = DESeqDataSetFromMatrix(countData=countData,
                             colData=colData,
                             design=~condition)

dds = DESeq(dds)
dds

res <- results(dds)
res

res = results(dds, contrast=c("condition", "hoxa1_kd", "control_sirna"))
mcols(res, use.names = TRUE)

summary(res)

res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)

resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

res <- res[order(res$pvalue),]
summary(res)

sum(res$padj < 0.1, na.rm=TRUE)

sum(res$pvalue < 0.05, na.rm=TRUE)

sum(!is.na(res$pvalue))

sum(res$padj < 0.06, na.rm=TRUE)

resSig <- subset(res, padj < 0.06)
head(resSig[ order( resSig$log2FoldChange ), ])

head(resSig[ order( resSig$log2FoldChange, decreasing=TRUE), ])

plotMA(res, ylim=c(-2,2))

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
res$symbol = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")

res$entrez = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")

res$name =   mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")

head(res, 10)

library(pathview)

library(gage)

library(gageData)

data(kegg.sets.hs)
data(sigmet.idx.hs)

# Focus on signaling and metabolic pathways only
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]

# Examine the first 3 pathways
head(kegg.sets.hs, 3)

foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)

# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs)

attributes(keggres)

# Look at the first few down (less) pathways
head(keggres$less)

pathview(gene.data=foldchanges, pathway.id="hsa04110")

## Focus on top 5 upregulated pathways here for demo purposes only
keggrespathways <- rownames(keggres$greater)[1:5]

# Extract the 8 character long IDs part of each string
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

pathview(gene.data=foldchanges, pathway.id=keggresids, species="hsa")

#1
sum(res$pvalue < 0.03, na.rm=TRUE)
sum(!is.na(res$pvalue))

#2
sum(res$padj < 0.0458, na.rm=TRUE)

#3
resSig <- subset(res, padj < 0.0458)
head(resSig[ order( resSig$log2FoldChange ), ])
head(resSig[ order( resSig$log2FoldChange, decreasing=TRUE), ])

#4
library(pheatmap)
library(ggplot2)
library(DESeq2)

# Import metadata
colData = read.csv("GSE37704_metadata.csv", row.names=1)
colData

# import as dataframe
countData = read.csv("GSE37704_featurecounts.csv",row.names=1)
# convert dataframe to matrix
countData = as.matrix(countData) 



ENSG00000171587 <- countData[19684,]
ENSG00000155011 <- countData[4840,]
ENSG00000101306 <- countData[17247,]
ENSG00000179855 <- countData[17705,]
ENSG00000188581 <- countData[16171,]
ENSG00000128052 <- countData[4601,]
ENSG00000141668 <- countData[17020,]
ENSG00000004799 <- countData[7449,]
ENSG00000109321 <- countData[4694,]
ENSG00000162892 <- countData[1765,]

sampleData <- rbind(ENSG00000171587, ENSG00000155011, ENSG00000101306, ENSG00000179855, ENSG00000188581, ENSG00000128052, ENSG00000141668, ENSG00000004799, ENSG00000109321, ENSG00000162892)

sampleDists = dist(t(sampleData),upper = TRUE)

annot_col = data.frame(colData$condition)
row.names(annot_col) <- rownames(colData)

sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = colnames(sampleData)
colnames(sampleDistMatrix) = colnames(sampleData)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE,
         annotation_col=annot_col)

#5
keggrespathways <- rownames(keggres$greater)[1:5]
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

keggrespathways <- rownames(keggres$less)[1:5]
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

