# Load file with header
rmsd <- read.delim("SRR701471.annovar.hg38_multianno.intronic.txt", sep = "\t", header = TRUE)

# Check column structure first
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

# Create bar plot matching the reference style
par(mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins for better spacing

# Create the bar plot for mutations (matching the reference image)
barplot(
  gene_mutation_data$NumMutations,
  names.arg = gene_mutation_data$Chromosome,
  col = "gray70",
  border = "black",
  main = "Intronic Variants per Chromosome",
  xlab = "Chromosome",
  ylab = "Number of Variants",
  ylim = c(0, max(gene_mutation_data$NumMutations) * 1.1),
  cex.main = 1.2,
  cex.lab = 1.1,
  cex.axis = 1.0,
  las = 2  # Rotate x-axis labels vertically
)

# Add a subtle grid for easier reading (optional)
grid(nx = NA, ny = NULL, col = "white", lty = "solid", lwd = 0.5)

# Redraw the bars on top of the grid
barplot(
  gene_mutation_data$NumMutations,
  names.arg = gsub("chr", "", gene_mutation_data$Chromosome),
  col = "gray",
  border = "black",
  add = TRUE,
  las = 2  # Rotate x-axis labels vertically
)

# Print summary statistics
cat("Total intronic variants:", sum(gene_mutation_data$NumMutations), "\n")
cat("Average variants per chromosome:", mean(gene_mutation_data$NumMutations), "\n")
cat("Chromosome with most variants:", 
    gene_mutation_data$Chromosome[which.max(gene_mutation_data$NumMutations)], 
    "with", max(gene_mutation_data$NumMutations), "variants\n")

# Print the data table
print("Summary by chromosome:")
print(gene_mutation_data[, c("Chromosome", "NumMutations", "NumGenes")])