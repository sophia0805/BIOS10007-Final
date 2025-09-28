# Load file with header
rmsd <- read.delim("SRR701471.annovar.hg38_multianno.txt", sep = "\t", header = TRUE)

# Define chromosome labels of interest (1–22 and X)
chromosomes <- paste0("chr", c(1:22, "X"))

# Initialize result data frame
bp_lengths <- data.frame(
  Chromosome = character(),
  LastRowPosition = numeric(),
  NumMutations = integer(),
  stringsAsFactors = FALSE
)

# Loop through each chromosome
for (chr in chromosomes) {
  chr_data <- rmsd[rmsd[[1]] == chr, ]
  if (nrow(chr_data) > 0) {
    last_row_pos <- chr_data[nrow(chr_data), 2]  
    row_count <- nrow(chr_data)
    bp_lengths <- rbind(
      bp_lengths,
      data.frame(
        Chromosome = chr,
        LastRowPosition = as.numeric(last_row_pos),
        NumMutations = row_count
      )
    )
  }
}

# Show result
print(bp_lengths)

# Plot with regression
plot(
  bp_lengths$NumMutations, bp_lengths$LastRowPosition,
  pch = 19,
  col = "darkblue",
  main = "Chromosome Length vs. Mutations",
  ylab = "Chromosome Length (Last Base Pair Position)",
  xlab = "Number of Mutations",
  ylim = c(0, max(bp_lengths$LastRowPosition) * 1.1)
)


model <- lm(LastRowPosition ~ NumMutations, data = bp_lengths)
abline(model, col = "red", lwd = 2)

# Calculate coefficient of determination
r_squared <- summary(model)$r.squared
text(
  x = max(bp_lengths$NumMutations) * 0.6,
  y = max(bp_lengths$LastRowPosition) * 0.95,
  labels = paste("R² =", round(r_squared, 4)),
  col = "red",
  cex = 1.2
)

