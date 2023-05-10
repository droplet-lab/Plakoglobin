library("pheatmap")
library(tximport)
library("DESeq2")
library(readr)
library(genefilter)
library(gplots)
library(RColorBrewer)
library(org.Mm.eg.db)
library(GOstats)


# Import data from featureCounts
countdata <- read.table("count_matrix.txt", header=TRUE, row.names = 1)

# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))

# Convert to matrix
countdata <- as.matrix(countdata)

# Assign condition 
condition <- factor(c("Bulk2D", "Bulk2D", "Bulk2D", "Bulk3D", "Bulk3D", "Bulk3D"))
coldata <- data.frame(row.names=colnames(countdata), condition)


# Import in  DESeq
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~ condition)

# Run the DESeq pipeline
dds <- DESeq(dds)
# Print size factors for each sample
sizeFactors(dds)
res <- results(dds)

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])


rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))

  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
}

# order by pvalue
resOrdered <- res[order(res$padj),]
write.csv(as.data.frame(resOrdered),file="Ordered_genes.csv")
write.csv(as.data.frame(counts(dds,normalized=TRUE)),file="Normalized_counts_all_genes.csv")
write.csv(as.data.frame(fpkm(dds)),file="FPKM_all_genes.csv")
# get differentially expressed gene matrix
sig <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj<0.01 & abs(resOrdered$log2FoldChange)>=0.1,]


# select genes
selected <- rownames(sig);selected

# save selected genes
write.csv(as.data.frame(counts(dds,normalized=TRUE)[rownames(dds) %in% selected,]), file="All_selected_genes.csv")

)
