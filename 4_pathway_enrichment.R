# Load necessary libraries
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

# Load the functional impact variant data
functional_impact_variants <- read.csv("results/functional_impact_variants.csv")

# Extract the gene symbols for pathway analysis
genes <- unique(functional_impact_variants$Symbol)

# Map Ensembl Gene IDs to Entrez IDs
gene_mapping <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform Reactome pathway enrichment analysis
pathway_enrichment <- enrichPathway(gene = gene_mapping$ENTREZID, pvalueCutoff = 0.05, readable = TRUE)

# Save the pathway enrichment results
write.csv(as.data.frame(pathway_enrichment), file = "results/pathway_enrichment_results.csv", row.names = FALSE)
