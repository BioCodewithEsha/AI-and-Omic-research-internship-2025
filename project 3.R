
# Load filtered data (from your saved CSV)
filtered_data <- read.csv("GSE124226_RMA_Normalized_Expression.csv", row.names = 1)

# (Optional) If you filtered it already earlier and saved filtered version:
# filtered_data <- read.csv("GSE124226_Filtered_Expression.csv", row.names = 1)

# Check
dim(filtered_data)
head(filtered_data[, 1:4])




#### Step 1: Map probe IDs to gene symbols and report statistics ####

# Install/load required Bioconductor package for GPL570 (HG-U133 Plus 2.0)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("hgu133plus2.db", quietly = TRUE)) BiocManager::install("hgu133plus2.db")
library(AnnotationDbi)
library(hgu133plus2.db)

# Use filtered_data from Assignment 1
if (!exists("filtered_data")) stop("filtered_data not found. Load filtered_data before running this.")

# Get probe IDs
probe_ids <- rownames(filtered_data)
length(probe_ids)  # total probes present after filtering

# Query annotation: get SYMBOL for each probe (select returns possibly multiple rows if multi-mapping)
ann <- AnnotationDbi::select(hgu133plus2.db,
                             keys = probe_ids,
                             columns = c("PROBEID","SYMBOL"),
                             keytype = "PROBEID")

# Show a small sample
cat("\n--- Example mappings (first 12 rows) ---\n")
print(head(ann, 12))

# Summary statistics
total_probes <- length(probe_ids)
probes_with_symbol_rows <- ann[!is.na(ann$SYMBOL) & ann$SYMBOL != "", ]
unique_probes_mapped <- length(unique(probes_with_symbol_rows$PROBEID))
unique_genes_mapped <- length(unique(probes_with_symbol_rows$SYMBOL))
probes_without_mapping <- setdiff(probe_ids, unique(probes_with_symbol_rows$PROBEID))
num_unmapped_probes <- length(probes_without_mapping)

cat("\n--- Mapping summary ---\n")
cat("Total probes (filtered_data):", total_probes, "\n")
cat("Probes with >=1 SYMBOL mapping (unique probes):", unique_probes_mapped, "\n")
cat("Number of unique gene symbols mapped:", unique_genes_mapped, "\n")
cat("Probes with NO mapping:", num_unmapped_probes, "\n")

# How many probes map to the same gene? List top genes by probe count
probe_counts_per_gene <- table(probes_with_symbol_rows$SYMBOL)
probe_counts_per_gene_sorted <- sort(probe_counts_per_gene, decreasing = TRUE)

cat("\n--- Top genes by number of mapping probes (top 20) ---\n")
print(head(probe_counts_per_gene_sorted, 20))

# Show a few example probes that have multiple SYMBOL rows (multi-mapping)
multi_map_probes <- names(which(table(probes_with_symbol_rows$PROBEID) > 1))
cat("\nNumber of probes that returned multiple annotation rows (same probe → multiple SYMBOL rows):", length(multi_map_probes), "\n")
if (length(multi_map_probes) > 0) {
  cat("\nExample probe(s) with multiple annotation rows:\n")
  print(unique(head(multi_map_probes, 10)))
  cat("\nShow full annotation rows for the first example probe:\n")
  print(ann[ann$PROBEID == multi_map_probes[1], ])
}

# Save the annotation table for inspection
write.csv(ann, file = "probe_to_symbol_mapping_table.csv", row.names = FALSE)
cat("\nSaved full mapping table to 'probe_to_symbol_mapping_table.csv' in working directory.\n")

# Prepare a simple one-to-one mapping per probe (choose first SYMBOL if multiple)
probe_to_symbol_first <- tapply(probes_with_symbol_rows$SYMBOL, probes_with_symbol_rows$PROBEID, function(x) x[1])
# probe_to_symbol_first is named vector: names = probeIDs, values = first SYMBOL

# Quick check: how many probes would remain if we drop unmapped probes and collapse by probe->first SYMBOL?
num_probes_after_map <- length(probe_to_symbol_first)
cat("\nProbes retained if using first SYMBOL per probe:", num_probes_after_map, "\n")

# Print out instructions for next step:
cat("\n--- NEXT ---\nRun this and paste the printed mapping summary above (the 'Mapping summary' and 'Top genes' sections).\n")
cat("After you paste those numbers, I'll provide the code to handle duplicates (collapse probes to gene-level)\n")




#### Step 2 (Revised): Collapse probes to gene-level matrix ####
# ---------------------------------------------------------------

library(matrixStats)

# Remove probes with no gene symbol
mapped_probes <- names(probe_to_symbol_first)
filtered_data_mapped <- filtered_data[mapped_probes, ]

# Map probe IDs to gene symbols
gene_symbols <- unname(probe_to_symbol_first[rownames(filtered_data_mapped)])
filtered_data_mapped$Gene <- gene_symbols

# Collapse duplicates: take mean expression per gene
cat("Collapsing duplicate probes by gene symbol (mean expression)...\n")

gene_level_expr <- aggregate(. ~ Gene, data = filtered_data_mapped, FUN = mean)

# Set gene names as rownames and remove redundant 'Gene' column
rownames(gene_level_expr) <- gene_level_expr$Gene
gene_level_expr <- gene_level_expr[, -1]

# Check summary
cat("Gene-level matrix dimensions:", dim(gene_level_expr), "\n")
head(gene_level_expr[, 1:4])

# Save for later use
write.csv(gene_level_expr, "GSE124226_GeneLevel_Expression.csv", row.names = TRUE)
cat("\n✅ Saved gene-level expression matrix as 'GSE124226_GeneLevel_Expression.csv'\n")






#### Step 3: Differential Gene Expression using Limma ####
# --------------------------------------------------------

# Install/load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) BiocManager::install("EnhancedVolcano")

library(limma)
library(pheatmap)
library(EnhancedVolcano)

# --------------------------------------------------------
# Load gene-level expression matrix
# --------------------------------------------------------
if (!exists("gene_level_expr")) {
  gene_level_expr <- read.csv("GSE124226_GeneLevel_Expression.csv", row.names = 1)
}

# --------------------------------------------------------
# Define sample groups (4 Control vs 4 PCOS)
# --------------------------------------------------------
group <- factor(c("Control", "PCOS", "PCOS", "PCOS",
                  "Control", "Control", "PCOS", "Control"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
print(design)

# --------------------------------------------------------
# Fit linear model and compute DEGs
# --------------------------------------------------------
fit <- lmFit(gene_level_expr, design)
contrast_matrix <- makeContrasts(PCOSvsControl = PCOS - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Get all DE results
deg_results <- topTable(fit2, number = Inf, adjust.method = "BH")
head(deg_results)

# --------------------------------------------------------
# Save DEG tables
# --------------------------------------------------------
upregulated <- deg_results[deg_results$logFC > 1 & deg_results$adj.P.Val < 0.05, ]
downregulated <- deg_results[deg_results$logFC < -1 & deg_results$adj.P.Val < 0.05, ]

write.csv(deg_results, "DEG_results_complete.csv")
write.csv(upregulated, "DEG_results_upregulated.csv")
write.csv(downregulated, "DEG_results_downregulated.csv")

cat("✅ DEG tables saved as CSV files.\n")
cat("Upregulated genes:", nrow(upregulated), "\n")
cat("Downregulated genes:", nrow(downregulated), "\n")

# --------------------------------------------------------
# Volcano Plot
# --------------------------------------------------------
png("VolcanoPlot_PCOS_vs_Control.png", width = 1200, height = 900)
EnhancedVolcano(deg_results,
                lab = rownames(deg_results),
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Volcano Plot: PCOS vs Control',
                subtitle = 'Limma DE analysis',
                caption = paste("Total DEGs:", nrow(upregulated) + nrow(downregulated)))
dev.off()
cat("✅ Volcano plot saved as 'VolcanoPlot_PCOS_vs_Control.png'\n")

# --------------------------------------------------------
# Heatmap (Top 25 DEGs)
# --------------------------------------------------------
top25_genes <- rownames(deg_results)[1:25]

# ✅ Fix: ensure annotation matches columns of expression matrix
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(gene_level_expr)

png("Heatmap_Top25_DEGs.png", width = 1000, height = 800)
pheatmap(gene_level_expr[top25_genes, ],
         scale = "row",
         annotation_col = annotation_col,
         main = "Top 25 Differentially Expressed Genes (PCOS vs Control)",
         cluster_cols = TRUE,
         fontsize_row = 8)
dev.off()
cat("✅ Heatmap saved as 'Heatmap_Top25_DEGs.png'\n")

# --------------------------------------------------------
# Summary for Report
# --------------------------------------------------------
cat("\n🧾 Summary:\n")
cat("• Comparison performed: PCOS vs Control\n")
cat("• Multiple probes per gene were collapsed by mean expression.\n")
cat("• Upregulated genes:", nrow(upregulated), "\n")
cat("• Downregulated genes:", nrow(downregulated), "\n")
cat("• DEG tables and plots saved in working directory.\n")




upregulated <- deg_results[deg_results$logFC > 0.5 & deg_results$adj.P.Val < 0.1, ]
downregulated <- deg_results[deg_results$logFC < -0.5 & deg_results$adj.P.Val < 0.1, ]

