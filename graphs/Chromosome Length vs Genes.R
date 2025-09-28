# Load file with header
rmsd <- read.delim("SRR701471.annovar.hg38_multianno.txt", sep = "\t", header = TRUE)

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
    if ("Gene.refGene" %in% colnames(chr_data)) {
      unique_genes <- length(unique(chr_data$Gene.refGene[chr_data$Gene.refGene != "" & !is.na(chr_data$Gene.refGene)]))
    }
    
    mutation_count <- nrow(chr_data)
    
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

# Plot with regression
plot(
  gene_mutation_data$NumGenes, gene_mutation_data$NumMutations,
  pch = 19,
  col = "darkblue",
  main = "Number of Genes vs. Mutations by Chromosome",
  xlab = "Number of Genes",
  ylab = "Number of Mutations",
  xlim = c(0, max(gene_mutation_data$NumGenes) * 1.1),
  ylim = c(0, max(gene_mutation_data$NumMutations) * 1.1)
)

# Add chromosome labels to points
text(
  gene_mutation_data$NumGenes, 
  gene_mutation_data$NumMutations,
  labels = gsub("chr", "", gene_mutation_data$Chromosome),
  pos = 3,
  cex = 0.8,
  col = "darkred"
)

# Fit regression model
model <- lm(NumMutations ~ NumGenes, data = gene_mutation_data)
abline(model, col = "red", lwd = 2)

# Calculate and display R-squared
r_squared <- summary(model)$r.squared
text(
  x = max(gene_mutation_data$NumGenes) * 0.7,
  y = max(gene_mutation_data$NumMutations) * 0.95,
  labels = paste("R² =", round(r_squared, 4)),
  col = "red",
  cex = 1.2
)

print(summary(model))