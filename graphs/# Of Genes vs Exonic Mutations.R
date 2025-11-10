# Load file with header
rmsd <- read.delim("/Users/ms201/OneDrive/Documents/projects/UChicago - MAC-RIBS/BIOS10007-Final/data/SRR701471.annovar.hg38_multianno.exonic.txt", sep = "\t", header = TRUE)

print("Column names:")
print(colnames(rmsd))
print("First few rows:")
print(head(rmsd))

# Define chromosome labels of interest (1–22 and X)
chromosomes <- paste0("chr", c(1:22, "X"))

# Initialize result data frame
gene_mutation_data <- data.frame(
  Chromosome = character(),
  NumGenes = integer(),
  NumMutations = integer(),
  stringsAsFactors = FALSE
)
# Loop through each chromosome
for (chr in chromosomes) {
  chr_data <- rmsd[rmsd[[1]] == chr, ]
  
  if (nrow(chr_data) > 0) {
    # Count total mutations for this chromosome
    mutation_count <- nrow(chr_data)
    
    # Count unique genes for this chromosome
    # Remove any empty or NA gene names
    gene_names <- chr_data[[7]]
    gene_names <- gene_names[!is.na(gene_names) & gene_names != "" & gene_names != "."]
    unique_genes <- length(unique(gene_names))
    
    gene_mutation_data <- rbind(
      gene_mutation_data,
      data.frame(
        Chromosome = chr,
        NumGenes = unique_genes,
        NumMutations = mutation_count
      )
    )
  }
}

# Show result
print(gene_mutation_data)

# Plot: Number of Genes vs Number of Mutations
plot(
  gene_mutation_data$NumMutations, gene_mutation_data$NumGenes,
  pch = 19,
  col = "darkblue",
  main = "Number of Genes vs. Number of Exonic Mutations",
  ylab = "Number of Genes",
  xlab = "Number of Exonic Mutations",
  ylim = c(0, max(gene_mutation_data$NumGenes) * 1.1),
  xlim = c(0, max(gene_mutation_data$NumMutations) * 1.1)
)

# Add chromosome labels to points
text(gene_mutation_data$NumMutations, gene_mutation_data$NumGenes, 
     labels = gsub("chr", "", gene_mutation_data$Chromosome), 
     pos = 3, cex = 0.8)

# Add regression line
model <- lm(NumGenes ~ NumMutations, data = gene_mutation_data)
summary_stats <- summary(model)
slope <- coef(model)[2]
conf_int <- confint(model, level = 0.95)

cat("\n=== REGRESSION RESULTS ===\n")
cat("Slope (β₁):", round(slope, 4), "genes per mutation\n")
cat("95% CI:", round(conf_int[2,1], 4), "to", round(conf_int[2,2], 4), "\n")
cat("R²:", round(summary_stats$r.squared, 4), "\n")
cat("p-value:", format.pval(summary_stats$coefficients[2,4]), "\n")
abline(model, col = "red", lwd = 2)

# Add R-squared value
r_squared <- summary(model)$r.squared
text(
  x = max(gene_mutation_data$NumMutations) * 0.6,
  y = max(gene_mutation_data$NumGenes) * 0.95,
  labels = paste("R² =", round(r_squared, 4)),
  col = "red",
  cex = 1.2
)
