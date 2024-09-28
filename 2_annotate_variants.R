# Load necessary libraries
library(biomaRt)
library(dplyr)

# Load the previously read variant data
source("scripts/1_read_variant_data.R")

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
