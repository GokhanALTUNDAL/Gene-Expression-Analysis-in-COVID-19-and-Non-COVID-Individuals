# 💡 Differential Gene Expression Analysis in COVID-19 Patients

This project aims to identify genes that are differentially expressed between COVID-19 patients and non-COVID individuals using RNA-Seq data (GSE157103). Statistical analysis and biological interpretation were performed using the **limma** package in R.

---

## 📁 Contents

- `R-Project4_FIXED_RELATIVE.R` – Full R analysis script
- `GSE157103_raw_counts_GRCh38.p13_NCBI.tsv` – Raw gene count data
- `GSE157103_family.soft` – Metadata including sample group labels
- `top_genes.csv` – All genes with fold change and adjusted p-values
- `top10_genes_named.csv` – Top 10 significant genes with symbols
- `go_results.csv` – Gene Ontology enrichment results
- `volcano.png` – Volcano plot of differential expression
- `heatmap.png` – Heatmap of top 50 genes
- `go_dotplot_adjusted.pdf` – Dotplot of enriched GO terms

---

## 🧪 Methods

1. **Data Preprocessing**: Raw count matrix filtered and log2-transformed
2. **Differential Expression Analysis**: Using `limma` with empirical Bayes adjustment
3. **Visualization**:
   - Volcano Plot
   - Heatmap (Top 50 genes)
4. **GO Enrichment Analysis**: via `clusterProfiler::enrichGO` using biological process (BP)

---

## 🔍 Key Findings

- Genes like **CDC6, MCM10, TOP2A** showed significant downregulation in COVID-19 samples.
- GO analysis revealed enrichment in processes such as:
  - **neuron apoptotic process**
  - **positive regulation of phosphorus metabolism**
  - **regulation of body fluid levels**

---

## 📦 Dependencies

- R (≥ 4.3)
- limma
- EnhancedVolcano
- pheatmap
- clusterProfiler
- org.Hs.eg.db

Install them via:

```r
BiocManager::install(c("limma", "EnhancedVolcano", "pheatmap", "clusterProfiler", "org.Hs.eg.db"))
```

---

## 📌 Dataset Source

GSE157103 – [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157103)

---

## 📃 License

MIT License

---

## ✍️ Author

Gökhan Altundal – 2025
