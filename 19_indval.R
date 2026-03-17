#!/usr/bin/env Rscript
# =============================================================================
# Indicator Species Analysis (IndVal)
# ISS Chile Pepper Rhizosphere Microbiome Study
#
# Identifies ASVs / genera specifically diagnostic for Space Flight or
# Terrestrial Soil conditions using the IndVal method (Dufrene & Legendre 1997).
# Uses indicspecies::multipatt() with 999 permutations.
#
# Prerequisites:
#   - version-2_integrated/exported_table_clean/feature-table.tsv
#   - version-2_integrated/exported_taxonomy/taxonomy.tsv
#   - version-2_integrated/integrated_metadata.tsv
#
# Run:
#   conda activate qiime2-amplicon
#   Rscript 19_indval.R
#
# Output (version-2_integrated/):
#   - IndVal_results.csv        : all ASVs with IndVal stat, p-value, group
#   - IndVal_significant.csv    : significant indicators (p < 0.05)
#   - Main_Fig_IndVal.png/.pdf  : top indicator taxa per group
# =============================================================================

suppressPackageStartupMessages({
  library(indicspecies)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

set.seed(42)
WORKDIR  <- "version-2_integrated"
N_PERM   <- 999
TOP_N    <- 15   # top indicators to plot per group

cat("============================================================\n")
cat("Indicator Species Analysis (IndVal)\n")
cat("============================================================\n\n")

# ── 1. Load data ──────────────────────────────────────────────
cat("[1] Loading data...\n")

otu_raw <- read.table(
  file.path(WORKDIR, "exported_table_clean/feature-table.tsv"),
  sep = "\t", header = TRUE, row.names = 1, skip = 1, check.names = FALSE,
  comment.char = ""
)

tax <- read.table(
  file.path(WORKDIR, "exported_taxonomy/taxonomy.tsv"),
  sep = "\t", header = TRUE, row.names = 1
)

meta <- read.table(
  file.path(WORKDIR, "integrated_metadata.tsv"),
  sep = "\t", header = TRUE, row.names = 1
)
meta <- meta[meta$Study_Group %in% c("Space_Flight", "Terrestrial_Soil"), , drop = FALSE]

# Samples as rows
otu <- t(otu_raw)
otu <- otu[rownames(otu) %in% rownames(meta), ]
meta <- meta[rownames(otu), , drop = FALSE]

# Relative abundance (% per sample)
otu_rel <- sweep(otu, 1, rowSums(otu), "/") * 100

cat(sprintf("  %d samples × %d ASVs\n", nrow(otu_rel), ncol(otu_rel)))
cat(sprintf("  Space_Flight: %d  |  Terrestrial_Soil: %d\n",
            sum(meta$Study_Group == "Space_Flight"),
            sum(meta$Study_Group == "Terrestrial_Soil")))

# ── 2. Genus-level collapse ────────────────────────────────────
cat("\n[2] Collapsing to genus level...\n")

get_genus <- function(taxon_str) {
  if (is.na(taxon_str)) return("Unassigned")
  parts <- strsplit(taxon_str, ";")[[1]]
  for (p in rev(parts)) {
    p <- trimws(p)
    if (startsWith(p, "g__")) {
      g <- sub("^g__", "", p)
      if (nchar(g) > 0) return(g)
    }
  }
  return("Unassigned")
}

tax$Genus <- sapply(tax$Taxon, get_genus)

# Map ASVs to genus, then sum per genus per sample
asv_ids <- colnames(otu_rel)
shared  <- intersect(asv_ids, rownames(tax))
otu_sub <- otu_rel[, shared]

# Build genus matrix
genus_df <- as.data.frame(t(otu_sub))
genus_df$Genus <- tax[shared, "Genus"]
genus_mat <- aggregate(. ~ Genus, data = genus_df, FUN = sum)
rownames(genus_mat) <- genus_mat$Genus
genus_mat$Genus <- NULL
genus_mat <- t(genus_mat)   # samples × genera

# Filter genera present in ≥5% of samples
prev <- colMeans(genus_mat > 0)
genus_mat <- genus_mat[, prev >= 0.05]
cat(sprintf("  %d genera (prevalence ≥5%% in samples)\n", ncol(genus_mat)))

# ── 3. Run IndVal ──────────────────────────────────────────────
cat(sprintf("\n[3] Running multipatt() — %d permutations...\n", N_PERM))

groups <- ifelse(meta$Study_Group == "Space_Flight", 1, 2)

indval_res <- multipatt(
  as.data.frame(genus_mat),
  groups,
  func      = "IndVal.g",
  control   = how(nperm = N_PERM),
  duleg     = TRUE
)

cat("  Done.\n")

# ── 4. Extract results ─────────────────────────────────────────
cat("\n[4] Extracting results...\n")

summary_df <- as.data.frame(indval_res$sign)
summary_df$Genus <- rownames(summary_df)

# Determine which group is indicated
summary_df$Group <- dplyr::case_when(
  summary_df$s.1 == 1 & summary_df$s.2 == 0 ~ "Space_Flight",
  summary_df$s.2 == 1 & summary_df$s.1 == 0 ~ "Terrestrial_Soil",
  TRUE                                         ~ "Both"
)

# Add mean relative abundance per group
space_mean <- colMeans(genus_mat[meta$Study_Group == "Space_Flight", , drop = FALSE])
terr_mean  <- colMeans(genus_mat[meta$Study_Group == "Terrestrial_Soil", , drop = FALSE])
summary_df$Space_Mean_RA  <- space_mean[summary_df$Genus]
summary_df$Terr_Mean_RA   <- terr_mean[summary_df$Genus]

# Save all results
write.csv(summary_df, file.path(WORKDIR, "IndVal_results.csv"), row.names = FALSE)
cat(sprintf("  Saved: IndVal_results.csv  (%d genera)\n", nrow(summary_df)))

# Significant indicators
sig <- summary_df[!is.na(summary_df$p.value) & summary_df$p.value < 0.05, ]
sig <- sig[sig$Group != "Both", ]
write.csv(sig, file.path(WORKDIR, "IndVal_significant.csv"), row.names = FALSE)
cat(sprintf("  Significant (p<0.05): %d genera\n", nrow(sig)))
cat(sprintf("    Space indicators:       %d\n", sum(sig$Group == "Space_Flight")))
cat(sprintf("    Terrestrial indicators: %d\n", sum(sig$Group == "Terrestrial_Soil")))

# ── 5. Plot ────────────────────────────────────────────────────
cat("\n[5] Plotting...\n")

plot_indval <- function(sig_df, top_n) {
  if (nrow(sig_df) == 0) {
    cat("  No significant indicators to plot.\n")
    return()
  }

  # Top indicators per group by IndVal stat
  top <- sig_df %>%
    group_by(Group) %>%
    slice_max(stat, n = top_n) %>%
    ungroup() %>%
    arrange(Group, stat)

  top$Genus <- factor(top$Genus, levels = top$Genus)
  top$Group_label <- ifelse(top$Group == "Space_Flight", "Space Flight", "Terrestrial Soil")

  p <- ggplot(top, aes(x = stat, y = Genus, fill = Group_label)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("p=%.3f", p.value)), hjust = -0.1, size = 2.8) +
    scale_fill_manual(values = c("Space Flight" = "#2980b9", "Terrestrial Soil" = "#27ae60")) +
    facet_wrap(~ Group_label, scales = "free_y") +
    xlim(0, 1.15) +
    labs(
      title    = "Indicator Species Analysis (IndVal)",
      subtitle = sprintf("Top %d indicators per group, p < 0.05 (%d permutations, seed=42)", top_n, N_PERM),
      x        = "IndVal Statistic",
      y        = NULL,
      fill     = "Group"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "none",
      axis.text.y      = element_text(face = "italic"),
      strip.background = element_rect(fill = "#f0f0f0"),
      strip.text       = element_text(face = "bold")
    )

  ggsave(file.path(WORKDIR, "Main_Fig_IndVal.png"), p, width = 12, height = max(6, nrow(top) * 0.28 + 2), dpi = 300)
  ggsave(file.path(WORKDIR, "Main_Fig_IndVal.pdf"), p, width = 12, height = max(6, nrow(top) * 0.28 + 2))
  cat("  Saved: Main_Fig_IndVal.png\n")
  cat("  Saved: Main_Fig_IndVal.pdf\n")
}

plot_indval(sig, TOP_N)

cat("\n============================================================\n")
cat("IndVal analysis complete.\n")
cat("============================================================\n")
