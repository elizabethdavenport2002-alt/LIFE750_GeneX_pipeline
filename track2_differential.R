# load library 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)

install.packages("pheatmap")
library(pheatmap)

# load dtaa 
dir()
counts <- read.table('data/gene_counts.tsv',header=TRUE, row.names=1, sep="\t")
meta <- read.table('data/sample_metadata.tsv',header=TRUE, row.names=1, sep="\t")

head(counts)
head(meta)

# makesure colum order in counts matches row order in meta 
counts <- counts[,rownames(meta)]
head(counts)

# Build DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = meta,
                              design    = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "mutant", "normal"))
res_df <- as.data.frame(res)

write.csv(res_df, "track2_results/DESeq2_results.csv")


# Significant DE genes (padj < 0.05, |log2FC| > 1)
sig <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(sig, "track2_results/DESeq2_significant.csv")

cat("Significant DE genes:", nrow(sig), "\n")
# 54 

vsd <- vst(dds, blind=TRUE)
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png("track2_results/PCA_plot.png", width=800, height=600)
plot(pcaData$PC1, pcaData$PC2,
     col = ifelse(pcaData$condition == "mutant", "tomato", "steelblue"),
     pch = 19, cex = 1.5,
     xlab = paste0("PC1: ", percentVar[1], "% variance"),
     ylab = paste0("PC2: ", percentVar[2], "% variance"),
     main = "PCA of RNA-seq samples")
legend("topright", legend=c("mutant","normal"),
       col=c("tomato","steelblue"), pch=19)
dev.off()


# MA plot
png("track2_results/MA_plot.png", width=800, height=600)
plotMA(res, ylim=c(-5,5), main="MA plot: mutant vs normal")
dev.off()


# Volcano plot
png("track2_results/volcano_plot.png", width=800, height=600)
with(res_df, {
  plot(log2FoldChange, -log10(pvalue),
       pch=20, cex=0.5, col="grey70",
       main="Volcano plot: mutant vs normal",
       xlab="log2 Fold Change", ylab="-log10(p-value)")
  points(log2FoldChange[padj < 0.05 & log2FoldChange >  1],
         -log10(pvalue[padj < 0.05 & log2FoldChange >  1]),
         col="tomato", pch=20)
  points(log2FoldChange[padj < 0.05 & log2FoldChange < -1],
         -log10(pvalue[padj < 0.05 & log2FoldChange < -1]),
         col="steelblue", pch=20)
  abline(v=c(-1,1), h=-log10(0.05), lty=2, col="black")
  legend("topright", legend=c("Up in mutant","Down in mutant","NS"),
         col=c("tomato","steelblue","grey70"), pch=20)
})
dev.off()


# Heatmap of top 30 DE genes
anno <- data.frame(condition = meta$condition)
rownames(anno) <- rownames(meta)
anno_colours <- list(condition = c(mutant="tomato", normal="steelblue"))

png("track2_results/heatmap_top30.png", width=800, height=600)
pheatmap(mat,
         annotation_col  = anno,
         annotation_colors = anno_colours,
         cluster_rows    = TRUE,
         cluster_cols    = TRUE,
         show_rownames   = TRUE,
         main            = "Top 30 DE genes (VST normalised)")

dev.off()


# up and down regulated csv with bound genes 
sig <- read.csv("track2_results/DESeq2_significant.csv", row.names=1)

up   <- subset(sig, log2FoldChange >  1 & padj < 0.05)
down <- subset(sig, log2FoldChange < -1 & padj < 0.05)

cat("Upregulated in mutant:  ", nrow(up), "\n")
cat("Downregulated in mutant:", nrow(down), "\n")

# check GeneX-bound regulation 
bound_genes <- c("GENE_0019","GENE_0034","GENE_0040","GENE_0079","GENE_0093",
                 "GENE_0116","GENE_0124","GENE_0184","GENE_0216","GENE_0223",
                 "GENE_0259","GENE_0270","GENE_0292","GENE_0361","GENE_0394",
                 "GENE_0408","GENE_0442","GENE_0477","GENE_0478","GENE_0489",
                 "GENE_0525","GENE_0549","GENE_0575","GENE_0587","GENE_0599",
                 "GENE_0695","GENE_0713","GENE_0735","GENE_0774","GENE_0806",
                 "GENE_0886","GENE_0930","GENE_0939","GENE_0978","GENE_0981")

bound_de <- sig[rownames(sig) %in% bound_genes, c("log2FoldChange","padj")]
bound_de$direction <- ifelse(bound_de$log2FoldChange > 0, "up", "down")
print(bound_de)

cat("\nBound targets up:  ", sum(bound_de$direction=="up"), "\n")
cat("Bound targets down:", sum(bound_de$direction=="down"), "\n")
#Bound targets up:   16 
#Bound targets down: 19 

# GENE_0079, GENE_0292 -very close to 1 (−1.001,1.001)

# save 
write.csv(bound_de, "track3_results/GeneX_bound_DE_genes.csv")
write.csv(up,   "track2_results/upregulated_genes.csv")
write.csv(down, "track2_results/downregulated_genes.csv")