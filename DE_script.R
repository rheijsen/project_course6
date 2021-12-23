#Running the DEseq analysis
library("DESeq2")

#It is absolutely critical that the columns of the count matrix and the rows 
#of the column data (information about samples) are in the same order. DESeq2 
#will not make guesses as to which column of the count matrix belongs to which 
#row of the column data, these must be provided to DESeq2 already in consistent 
#order.

samples <- c("SRR11068823", "SRR11068824","SRR11068825", "SRR11068826", "SRR11068827", "SRR11068828")
condition <- c(rep("untreated", times = 3), rep("treated", times = 3))


#make a experimental table
experiment.data <- data.frame(samples = samples, condition = condition)
experiment.data

#read in count table
counts <- as.matrix(read.csv("orginele_counts.txt", sep = "\t", col.names = samples))
head(counts)
#make dseq matrix
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = experiment.data,
                              design = ~ condition)
#Note on factor levels
#By default, R will choose a reference level for factors based on alphabetical 
#order. Then, if you never tell the DESeq2 functions which level you want to 
#compare against (e.g. which level represents the control group), the 
#comparisons will be based on the alphabetical order of the levels.

dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
dds$condition <- relevel(dds$condition, ref = "untreated")


dds_s <- DESeq(dds)
res <- results(dds_s)
res
#order on pvalues
#resOrdered <- res[order(res$pvalue),]
#head(resOrdered)
#summary(res)

################# plotting some results ###############
#plot MA plot
plotMA(res, ylim=c(-2,2))

#shrinking for visualization and ranking
resLFC <- lfcShrink(dds_s, coef="condition_treated_vs_untreated", type="apeglm")
plotMA(resLFC, ylim=c(-2,2))

#plot count
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

# plot PCA which needs some transformation of the data
vsd <- vst(dds_s, blind=FALSE)
rld <- rlog(dds_s, blind=FALSE)

plotPCA(vsd, intgroup=c("condition"))

#dispersion plot
plotDispEsts(dds_s)

#clustering
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 

library(pheatmap)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#volcano plot (see: https://github.com/kevinblighe/EnhancedVolcano)
library(EnhancedVolcano)
EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue')

######################### Filter results on FC en p-vals ########################
res
res_filtered <- res[!is.na(res$padj) & res$padj <= 0.1,] #
res_filtered

res_up <- as.data.frame(res_filtered[res_filtered$log2FoldChange >= 1,])
res_up

res_down <- as.data.frame(res_filtered[res_filtered$log2FoldChange <= -1,])
res_down

######################### Annotation #################
annotation <- read.csv("Annotation_Celegans.txt", sep = "\t", row.names = "Your.Input")
rownames(annotation)

#upregulated
res_up_anno <- merge(res_up, annotation, by='row.names')
res_up_anno_sorted <- res_up_anno[order(res_up_anno$log2FoldChange, decreasing = TRUE),]
head(res_up_anno_sorted)

#downregulated
res_down_anno <- merge(res_down, annotation, by = 'row.names')
res_down_anno_sorted <- res_down_anno[order(res_down_anno$log2FoldChange, decreasing = FALSE),]
head(res_down_anno_sorted)

#write to new file
concat_res = rbind(res_up_anno_sorted, res_down_anno_sorted)

write.table(concat_res, "DEGs_Celegans_course6B.txt", append = FALSE, sep = "\t", dec = ",",
            row.names = FALSE, col.names = TRUE)
