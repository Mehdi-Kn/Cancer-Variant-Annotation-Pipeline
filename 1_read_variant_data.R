# Load necessary libraries
library(VariantAnnotation)

# Define input VCF file path (make sure to update this if needed)
vcf_file <- "data/variants.vcf"

# Read the VCF file into R
vcf <- readVcf(vcf_file, genome = "hg38")

# Print a summary of the VCF data
print(vcf)
