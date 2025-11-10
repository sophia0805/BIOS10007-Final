rm(list = ls())
# Load exonic mutations file
rmsd_exonic <- read.delim("../data/SRR701471.annovar.hg38_multianno.exonic.txt", sep = "\t", header = TRUE)

# Load total mutations file (assuming you have the complete annotation file)
rmsd_total <- read.delim("../data/SRR701471.annovar.hg38_multianno.txt", sep = "\t", header = TRUE)

# Define chromosome labels of interest (1–22 and X)
chromosomes <- paste0("chr", c(1:22, "X"))

# Initialize result data frame
mutation_percentage_data <- data.frame(
  Chromosome = character(),
  ChromosomeNum = numeric(),
  TotalMutations = integer(),
  ExonicMutations = integer(),
  PercentExonic = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each chromosome
for (i in seq_along(chromosomes)) {
  chr <- chromosomes[i]
  
  # Count total mutations for this chromosome
  total_chr_data <- rmsd_total[rmsd_total[[1]] == chr, ]
  total_mutations <- nrow(total_chr_data)
  
  # Count exonic mutations for this chromosome
  exonic_chr_data <- rmsd_exonic[rmsd_exonic[[1]] == chr, ]
  exonic_mutations <- nrow(exonic_chr_data)
  
  # Calculate percentage only if there are total mutations
  if (total_mutations > 0) {
    percent_exonic <- (exonic_mutations / total_mutations) * 100
    
    # Create numeric chromosome for plotting (1-22, then 23 for X)
    chr_num <- ifelse(chr == "chrX", 23, as.numeric(gsub("chr", "", chr)))
    
    mutation_percentage_data <- rbind(
      mutation_percentage_data,
      data.frame(
        Chromosome = chr,
        ChromosomeNum = chr_num,
        TotalMutations = total_mutations,
        ExonicMutations = exonic_mutations,
        PercentExonic = percent_exonic
      )
    )
  }
}

# Show result
print("Exonic mutation percentages by chromosome:")
print(mutation_percentage_data)

# Set up plot with proper margins and styling
par(mar = c(5, 4, 4, 2) + 0.1)

# Create scatter plot matching the reference style
plot(
  mutation_percentage_data$ChromosomeNum, 
  mutation_percentage_data$PercentExonic,
  pch = 16,
  col = "purple",
  cex = 1.0,
  main = "% of Mutations that are Exonic",
  xlab = "Chromosome",
  ylab = "Percent Exonic",
  xlim = c(0.5, 23.5),
  ylim = c(6, 20),
  xaxt = "n",
  frame.plot = TRUE,
  cex.main = 1.2,
  cex.lab = 1.1
)

# Add custom x-axis with chromosome labels
x_positions <- 1:23
x_labels <- c(paste0("chr", 1:22), "chrX")
axis(1, at = x_positions, labels = x_labels, las = 2, cex.axis = 0.8)

# Add regression line
model <- lm(PercentExonic ~ ChromosomeNum, data = mutation_percentage_data)
abline(model, col = "red", lwd = 2)
summary_stats <- summary(model)
slope <- coef(model)[2]
conf_int <- confint(model, level = 0.95)

cat("Slope (β₁):", round(slope, 4), "# of genes per # of exonic mutations\n")
cat("95% CI:", round(conf_int[2,1], 4), "to", round(conf_int[2,2], 4), "\n")
cat("R²:", round(summary_stats$r.squared, 4), "\n")
cat("p-value:", format.pval(summary_stats$coefficients[2,4]), "\n")
abline(model, col = "red", lwd = 2)
# Add R-squared value in top right
r_squared <- summary(model)$r.squared
text(
  x = 21,
  y = 20,
  labels = bquote(R^2 == .(round(r_squared, 3))),
  col = "red",
  cex = 1.0
)