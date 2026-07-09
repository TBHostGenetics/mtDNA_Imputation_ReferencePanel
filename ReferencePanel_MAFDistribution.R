library(data.table)

# --- Read in files ---
legend <- fread("mtDNA_ReferencePanel.legend.gz", header = TRUE)      # id, position, a0, a1, ...
hap    <- fread("mtDNA_ReferencePanel.hap.gz", header = FALSE)         # 0/1 matrix, one row per variant

stopifnot(nrow(hap) == nrow(legend))  # sanity check

freq1 <- scan("freq1.txt")

legend$MAF <- pmin(freq1, 1 - freq1)

# --- Quick look ---
summary(legend$MAF)

library(data.table)
library(ggplot2)

legend$MAF <- pmin(freq1, 1 - freq1)

legend[, MAF_bin := cut(
  MAF,
  breaks = c(-Inf, 0, 1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.5),
  labels = c("Monomorphic", "Singletons", "0.00001 to 0.0001", "0.0001 to 0.001",
             "0.001 to 0.01", "0.01 to 0.05", ">=0.05"),
  right = TRUE
)]

bin_counts <- legend[, .N, by = MAF_bin][order(MAF_bin)]
bin_counts


plot <- ggplot(bin_counts, aes(x = MAF_bin, y = N)) +
  geom_col(fill = "#356920") +
  ylim(0, 6000) +
  theme_minimal() +
  labs(
       x = "MAF bin", y = "Number of variants") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
plot

ggsave("VariantDistribution_byMAF_ReferencePanel.png", plot = plot, height = 6, width = 10)
