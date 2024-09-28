# Load necessary libraries
library(dplyr)

# Load the annotated variant data
annotated_variants <- read.csv("results/annotated_variants.csv")

# Convert SIFT and PolyPhen scores to numeric
annotated_variants <- annotated_variants %>%
  mutate(
    SIFT_score = as.numeric(SIFT_score),
    PolyPhen_score = as.numeric(PolyPhen_score)
  )

# Save the results
write.csv(annotated_variants, file = "results/functional_impact_variants.csv", row.names = FALSE)
