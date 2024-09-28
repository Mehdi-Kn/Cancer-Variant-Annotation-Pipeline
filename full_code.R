# full_code.R
# Comprehensive script for the Cancer Variant Annotation Pipeline

# Load necessary libraries
library(VariantAnnotation)
library(biomaRt)
library(dplyr)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)

# ---------------------------
# Step 1: Read Variant Data
# ---------------------------
# Define input VCF file path (make sure to update this if needed)
vcf_file <- "data/variants.vcf"

# Read the VCF file into R
vcf <- readVcf(vcf_file, genome = "hg38")

# Print a summary of the VCF data
cat("Step 1: VCF file successfully read.\n")
print(vcf)


# ---------------------------
# Step 2: Annotate Variants
# ---------------------------
# Connect to Ensembl BioMart
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Prepare variant IDs
variants <- as.data.frame(rowRanges(vcf))
variants$variant_id <- paste0(variants$seqnames, "_", variants$start, "_", variants$REF, "/", variants$ALT)

# Retrieve annotations from Ensembl VEP
vep_results <- getBM(
  attributes = c('Uploaded_variation', 'Consequence', 'Gene', 'Symbol', 'Protein_position',
                 'Amino_acids', 'Codons', 'SIFT_score', 'SIFT_prediction',
                 'PolyPhen_score', 'PolyPhen_prediction'),
  filters = 'Uploaded_variation',
  values = variants$variant_id,
  mart = ensembl
)

# Merge annotations with original variant data
annotated_variants <- left_join(variants, vep_results, by = c("variant_id" = "Uploaded_variation"))

# Save the annotated variants to a CSV file
write.csv(annotated_variants, file = "results/annotated_variants.csv", row.names = FALSE)

cat("Step 2: Variant annotation completed.\n")


# ---------------------------
# Step 3: Predict Functional Impact
# ---------------------------
# Convert SIFT and PolyPhen scores to numeric
annotated_variants <- annotated_variants %>%
  mutate(
    SIFT_score = as.numeric(SIFT_score),
    PolyPhen_score = as.numeric(PolyPhen_score)
  )

# Save the results
write.csv(annotated_variants, file = "results/functional_impact_variants.csv", row.names = FALSE)

cat("Step 3: Functional impact prediction completed.\n")


# ---------------------------
# Step 4: Pathway Enrichment Analysis
# ---------------------------
# Extract the gene symbols for pathway analysis
genes <- unique(annotated_variants$Symbol)

# Map Ensembl Gene IDs to Entrez IDs
gene_mapping <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform Reactome pathway enrichment analysis
pathway_enrichment <- enrichPathway(gene = gene_mapping$ENTREZID, pvalueCutoff = 0.05, readable = TRUE)

# Save the pathway enrichment results
write.csv(as.data.frame(pathway_enrichment), file = "results/pathway_enrichment_results.csv", row.names = FALSE)

cat("Step 4: Pathway enrichment analysis completed.\n")


# ---------------------------
# Step 5: Visualization
# ---------------------------
# Plot the distribution of variant consequences
consequence_plot <- ggplot(annotated_variants, aes(x = Consequence)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  xlab("Variant Consequence") +
  ylab("Count") +
  ggtitle("Distribution of Variant Consequences") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the consequence plot
ggsave("results/variant_consequences.png", plot = consequence_plot)

# Visualize pathway enrichment (ReactomePA format)
dotplot_pathway <- dotplot(pathway_enrichment, showCategory = 10) +
  ggtitle("Top 10 Enriched Pathways")

# Save the pathway enrichment plot
ggsave("results/pathway_enrichment.png", plot = dotplot_pathway)

cat("Step 5: Visualization completed. Check the results directory for plots and CSV files.\n")
