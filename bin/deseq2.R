library("DESeq2")
library("RColorBrewer")
library("pheatmap")

# Create empty pdf 
pdf("deseq2.pdf")

# Path to directory of htseq outputs 
args <- commandArgs(trailingOnly = TRUE)
directory <- args[1]

# Retrieve htseq count files
sampleFiles <- grep("treated",list.files(directory),value=TRUE)

# Retrieve condition (treated/untreated) for each file
sampleCondition <- sub(".*_(treated|untreated).*", "\\1", sampleFiles)

# Create a table of sample metadata 
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
# Level                          
sampleTable$condition <- factor(sampleTable$condition)

# Read htseq datasets 
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)

# Specify the reference level:
ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "untreated")

# DE analysis
ddsHTSeq <- DESeq(ddsHTSeq, parallel=TRUE)
res <- results(ddsHTSeq)

# Shrinkage
resLFC <- lfcShrink(ddsHTSeq, coef="condition_treated_vs_untreated", type="apeglm")

# Export
resOrdered <- res[order(res$pvalue),]
# sum(res$padj < 0.1, na.rm=TRUE) ## adjusted p-values were less than 0.1

vsd <- varianceStabilizingTransformation(ddsHTSeq, blind=FALSE)
# PCA
plotPCA(vsd, intgroup=c("condition"))

# Heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Dispersion estimate
plotDispEsts(ddsHTSeq)

# MA plot
plotMA(resLFC, ylim=c(-2,2))

# P value
h1 <- hist(res$pvalue, breaks=0:50/50, plot=FALSE, color='powderblue')
plot(h1, main="P-value Histogram", xlab="P-value", ylab="Frequency")

dev.off()
