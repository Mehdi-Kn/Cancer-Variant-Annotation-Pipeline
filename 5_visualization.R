# Load necessary libraries
library(ggplot2)
library(ReactomePA)

# Load the annotated variant data
annotated_variants <- read.csv("results/functional_impact_variants.csv")

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

# Load the pathway enrichment results
pathway_enrichment_results <- read.csv("results/pathway_enrichment_results.csv")

# Visualize pathway enrichment (ReactomePA format)
dotplot_pathway <- dotplot(as(pathway_enrichment_results, "enrichResult"), showCategory = 10) +
  ggtitle("Top 10 Enriched Pathways")

# Save the pathway enrichment plot
ggsave("results/pathway_enrichment.png", plot = dotplot_pathway)
