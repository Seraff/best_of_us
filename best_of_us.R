###############################################################################
## Description: Determine sugnificant genes for differntial expression analysis
## Author: Masoud Nazarizadeh, Serafim Nenarokov
###############################################################################

Sys.setenv(LANG = "en")

library(DESeq2)
library(edgeR)

library(VennDiagram)
library(pheatmap)
library(ggplot2)
library(argparse)
library(glue)

output_path_for <- function (filename) {
  return(file.path(args$output, filename))
}

parse_arguments <- function () {
  parser <- ArgumentParser()
  parser$add_argument("-c", "--count_table", type="character", required=TRUE,
                      help="Path to the table with counts data")
  parser$add_argument("-m", "--metadata", type="character", required=TRUE,
                      help="Path to the table experiment metadata")
  parser$add_argument("-o", "--output", type="character", required=TRUE,
                      help="Path to the output folder (will be created if doesn't exist)")
  parser$add_argument("--padj_cutoff", type="double", default="0.05",
                      help="Cut-off for adjusted p-value; all the genes with padj <= padj_cutoff will be ignored")
  args <- parser$parse_args()
  return(args)
}

args <- c()

if (Sys.getenv("RSTUDIO") == "1") {
  library(rstudioapi)
  my_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  
  args$count_table = file.path(my_dir, "data/results/counts.txt")
  args$metadata = file.path(my_dir, "data/results/metadata.csv")
  args$output = file.path(my_dir, "data/results/plots")
  args$padj_cutoff = 0.05
} else {
  args <- parse_arguments()
}

# Reading count matrix
count_data <- read.table(args$count_table, header=TRUE, row.names=1)
count_data <- count_data[ ,6:ncol(count_data)]
colnames(count_data) <- gsub("\\.[sb]am$", "", colnames(count_data))
colnames(count_data) <- gsub("(\\w+\\.)+(\\w+)", "\\2", colnames(count_data))
colnames(count_data) <- gsub("_mapping", "", colnames(count_data))
count_data <- as.matrix(count_data)

# Reading metadata
metadata <- read.csv(args$metadata, row.names = 1)

# Creating the output folder if it doesn't exist
system(glue("mkdir -p {args$output}"))

## DESeq2 analysis ##
metadata_deseq <- metadata[,1:ncol(metadata)]
metadata_deseq = t(metadata_deseq)
deseq_data <- DESeqDataSetFromMatrix(countData = count_data,
                                     colData = metadata_deseq,
                                     design = ~ state)
deseq_data <- DESeq(deseq_data)
result_deseq2 <- results(deseq_data)

## edgeR analysis ##
metadata_edgeR <- c(t(metadata))
dge <- DGEList(counts = count_data, group = metadata_edgeR)
dge <- calcNormFactors(dge)
design <- model.matrix(~metadata_edgeR)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = 2)
result_edgeR <- topTags(qlf, n = Inf)
result_edgeR_df <- as.data.frame(result_edgeR)

## Find intersecting genes ##
sig_genes_deseq2 <- rownames(result_deseq2[which(result_deseq2$padj < args$padj_cutoff),])
sig_genes_edgeR <- rownames(result_edgeR_df[which(result_edgeR_df$FDR < args$padj_cutoff),])
intersect_genes <- intersect(sig_genes_deseq2, sig_genes_edgeR)

## Plotting plots ##

# Venn diagram
venn_data <- list(
  DESeq2 = sig_genes_deseq2,
  edgeR = sig_genes_edgeR
)
venn.plot <- venn.diagram(venn_data, filename = output_path_for("vienn.png"),
                          disable.logging = TRUE,
                          imagetype = "png")

# Heatmap
counts_intersect <- count_data[intersect_genes, ]
log_counts_intersect <- log2(counts_intersect + 1)

pheatmap(log_counts_intersect, 
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = metadata$state,
         filename = output_path_for("heatmap.png"))

# PCA plot
pca_data <- t(assay(rlog(deseq_data[intersect_genes, ])))
pca_result <- prcomp(pca_data)
pca_df <- data.frame(PC1 = pca_result$x[, 1],
                     PC2 = pca_result$x[, 2],
                     Condition = t(metadata))

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = state)) +
            geom_point(size = 4) +
            theme_bw() +
            xlab("PC1") +
            ylab("PC2") +
            ggtitle("PCA Plot of Overlapped DEGs")

ggsave(output_path_for("pca.png"), plot = pca_plot)

# Volcano plot
volcano_data <- as.data.frame(result_deseq2[intersect_genes, ])

volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
                geom_point(alpha = 0.5) +
                theme_bw() +
                xlab("log2 Fold Change (DESeq2)") +
                ylab("-log10 Adjusted p-value (DESeq2)") +
                ggtitle("Volcano plot (Overlapped DEGs)")

ggsave(output_path_for("volcano.png"), plot = volcano_plot)

# MA plot
dge_intersect <- dge[intersect_genes,]
M_intersect <- result_edgeR_df[intersect_genes, ]$logFC
A_intersect <- log2(rowMeans(dge_intersect$counts + 0.5))

ma_plot <- ggplot(data.frame(A = A_intersect, M = M_intersect), aes(x = A, y = M)) +
           geom_point(alpha = 0.5) +
           theme_bw() +
           xlab("A") +
           ylab("M") +
           ggtitle("MA plot (Overlapped DEGs)")

ggsave(output_path_for("ma.png"), plot = ma_plot)

# List of genes
write.table(intersect_genes, output_path_for("significant_genes.txt"),
          quote=FALSE,
          row.names=FALSE,
          col.names=FALSE, na='')
