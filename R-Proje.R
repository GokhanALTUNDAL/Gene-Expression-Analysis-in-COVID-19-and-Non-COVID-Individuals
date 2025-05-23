
options(scipen = 999)

# Gerekli kutuphaneler
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "EnhancedVolcano", "pheatmap", "clusterProfiler", "org.Hs.eg.db"))

library(limma)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

# Veri yukleme
raw_counts <- read.delim("GSE157103_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, sep = "\t", check.names = FALSE)
rownames(raw_counts) <- raw_counts$GeneID
raw_counts <- raw_counts[, -1]

# Etiket verisi
meta_lines <- readLines("GSE157103_family.soft")
sample_lines <- grep("Series_sample_id", meta_lines, value = TRUE)
sample_ids <- gsub("!Series_sample_id = ", "", sample_lines)
disease_lines <- grep("disease state", meta_lines, value = TRUE)
disease_labels <- gsub("!Sample_characteristics_ch1 = disease state: ", "", disease_lines)
sample_info <- data.frame(SampleID = sample_ids, Group = disease_labels)

# Ortak ornekler
common_samples <- intersect(colnames(raw_counts), sample_info$SampleID)
raw_counts <- raw_counts[, common_samples]
sample_info <- sample_info[sample_info$SampleID %in% common_samples, ]

# log2 donusumm ve kontrol
log_counts <- log2(raw_counts + 1)
if (max(log_counts, na.rm = TRUE) > 50) stop("log2 donusumu eksik olabilir!")

# limma analizi
group <- factor(sample_info$Group)
design <- model.matrix(~ group)
fit <- lmFit(log_counts, design)
fit <- eBayes(fit)
results <- topTable(fit, coef = 2, number = Inf, sort.by = "P")

# Volcano plot
png("volcano.png", width = 800, height = 600)
EnhancedVolcano(results,
                lab = rownames(results),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'COVID-19 vs non-COVID-19')
dev.off()

# Heatmap
top_gene_ids <- rownames(results)[1:50]
heat_data <- log_counts[top_gene_ids, ]
annotation_col <- data.frame(Group = sample_info$Group)
rownames(annotation_col) <- sample_info$SampleID
png("heatmap.png", width = 800, height = 600)
pheatmap(heat_data,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         show_colnames = FALSE,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Top 50 Differentially Expressed Genes")
dev.off()

# GO analizi
deg_ids <- rownames(results[results$adj.P.Val < 0.05, ])
ego <- enrichGO(gene = deg_ids,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable = TRUE)

# GO dotplot
pdf("go_dotplot_adjusted.pdf", width = 10, height = 8)
dotplot(ego, showCategory = 15, title = "GO: Biological Processes", font.size = 9)
dev.off()

# CSV Ciktilarr
write.csv(as.data.frame(ego), "go_results.csv", row.names = FALSE)

# Top 10 gen
top10 <- head(results, 10)
top10_symbols <- mapIds(org.Hs.eg.db,
                        keys = rownames(top10),
                        column = "SYMBOL",
                        keytype = "ENTREZID",
                        multiVals = "first")
top10_table <- data.frame(
  EntrezID = rownames(top10),
  Symbol = top10_symbols,
  logFC = top10$logFC,
  adj.P.Val = top10$adj.P.Val
)
top_named <- na.omit(top10_table)
write.csv(top_named, "top10_genes_named.csv", row.names = FALSE)

# Tam sonuclari kaydet
write.csv(results, "top_genes.csv")
