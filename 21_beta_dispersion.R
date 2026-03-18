#!/usr/bin/env Rscript
# =============================================================================
# Beta-Dispersion Analysis (PERMDISP2)
# ISS Chile Pepper Rhizosphere Microbiome Study
#
# Tests whether within-group beta-diversity dispersion differs between
# Space Flight and Terrestrial Soil.
#
# Ecological interpretation:
#   Higher dispersion in Space → community composition is more variable
#   (less predictable) within the group, consistent with stochastic assembly
#   shown by βNTI analysis (16_bnti_rcbray_analysis.R).
#
# Method:
#   vegan::betadisper() computes distance from each sample to its group
#   centroid in ordination space (Bray-Curtis or UniFrac).
#   vegan::permutest()  tests whether group dispersions differ
#   (H0: within-group distances to centroid are equal across groups).
#
# Prerequisites:
#   - version-2_integrated/exported_table_clean/feature-table.tsv
#   - version-2_integrated/integrated_metadata.tsv
#   - version-2_integrated/exported_tree/tree.nwk  (for UniFrac)
#
# Run in qiime2-amplicon environment:
#   conda activate qiime2-amplicon
#   Rscript 21_beta_dispersion.R
#
# Output (version-2_integrated/):
#   - BetaDispersion_results.csv     : per-sample distances to centroid
#   - BetaDispersion_summary.csv     : group means ± SD, p-value
#   - Main_Fig_BetaDispersion.png/.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(ape)
})

set.seed(42)
WORKDIR   <- "version-2_integrated"
N_PERM    <- 999

cat("============================================================\n")
cat("Beta-Dispersion (PERMDISP2) Analysis\n")
cat("============================================================\n\n")

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
cat("[1] Loading data...\n")

otu_raw <- read.table(
  file.path(WORKDIR, "exported_table_clean/feature-table.tsv"),
  sep = "\t", header = TRUE, row.names = 1,
  skip = 1, check.names = FALSE, comment.char = ""
)
otu <- t(otu_raw)

meta <- read.table(
  file.path(WORKDIR, "integrated_metadata.tsv"),
  sep = "\t", header = TRUE, row.names = 1
)
meta <- meta[meta$Study_Group %in% c("Space_Flight", "Terrestrial_Soil"), , drop = FALSE]
otu  <- otu[rownames(otu) %in% rownames(meta), ]
meta <- meta[rownames(otu), , drop = FALSE]
cat(sprintf("  %d samples: Space=%d, Terrestrial=%d\n",
  nrow(otu),
  sum(meta$Study_Group == "Space_Flight"),
  sum(meta$Study_Group == "Terrestrial_Soil")))

# Relative abundance
otu_rel <- otu / rowSums(otu)
group   <- factor(meta$Study_Group, levels = c("Space_Flight", "Terrestrial_Soil"))

# ---------------------------------------------------------------------------
# 2. Bray-Curtis beta-dispersion
# ---------------------------------------------------------------------------
cat("\n[2] Bray-Curtis beta-dispersion...\n")

bc_dist <- vegdist(otu_rel, method = "bray")
bd_bc   <- betadisper(bc_dist, group, type = "centroid")
perm_bc <- permutest(bd_bc, permutations = N_PERM, pairwise = TRUE)

cat("  Group distances to centroid (mean ± SD):\n")
for (grp in levels(group)) {
  d <- bd_bc$distances[group == grp]
  cat(sprintf("    %-20s %.4f ± %.4f\n", grp, mean(d), sd(d)))
}
cat(sprintf("  Permutation test p = %.4f\n", perm_bc$tab["Groups", "Pr(>F)"]))

# F-statistic and p
F_bc <- perm_bc$tab["Groups", "F"]
p_bc <- perm_bc$tab["Groups", "Pr(>F)"]

# ---------------------------------------------------------------------------
# 3. UniFrac beta-dispersion (if tree is available)
# ---------------------------------------------------------------------------
uf_available <- FALSE
tree_path <- file.path(WORKDIR, "exported_tree/tree.nwk")

if (file.exists(tree_path)) {
  cat("\n[3] UniFrac beta-dispersion...\n")
  tryCatch({
    suppressPackageStartupMessages(library(picante))
    tree <- read.tree(tree_path)
    otu_matched <- match.phylo.comm(tree, otu_rel)
    otu_uf  <- otu_matched$comm
    tree_uf <- otu_matched$phy
    group_uf <- factor(meta[rownames(otu_uf), "Study_Group"],
                       levels = c("Space_Flight", "Terrestrial_Soil"))

    # Weighted UniFrac via picante
    uf_dist  <- unifrac(otu_uf, tree_uf)
    bd_uf    <- betadisper(as.dist(uf_dist), group_uf, type = "centroid")
    perm_uf  <- permutest(bd_uf, permutations = N_PERM, pairwise = TRUE)
    F_uf <- perm_uf$tab["Groups", "F"]
    p_uf <- perm_uf$tab["Groups", "Pr(>F)"]

    cat("  Group UniFrac distances to centroid (mean ± SD):\n")
    for (grp in levels(group_uf)) {
      d <- bd_uf$distances[group_uf == grp]
      cat(sprintf("    %-20s %.4f ± %.4f\n", grp, mean(d), sd(d)))
    }
    cat(sprintf("  Permutation test p = %.4f\n", p_uf))
    uf_available <- TRUE
  }, error = function(e) {
    cat(sprintf("  UniFrac skipped: %s\n", conditionMessage(e)))
  })
} else {
  cat("\n[3] UniFrac skipped (tree not found).\n")
}

# ---------------------------------------------------------------------------
# 4. Save results
# ---------------------------------------------------------------------------
cat("\n[4] Saving results...\n")

# Per-sample distances
res_df <- data.frame(
  SampleID      = names(bd_bc$distances),
  Study_Group   = as.character(group[names(bd_bc$distances)]),
  BC_Distance   = round(bd_bc$distances, 6),
  stringsAsFactors = FALSE
)
if (uf_available) {
  uf_dist_vec <- bd_uf$distances[res_df$SampleID]
  res_df$UF_Distance <- round(ifelse(is.na(uf_dist_vec), NA, uf_dist_vec), 6)
}
write.csv(res_df, file.path(WORKDIR, "BetaDispersion_results.csv"), row.names = FALSE)
cat("  Saved: BetaDispersion_results.csv\n")

# Summary statistics
sum_rows <- lapply(levels(group), function(grp) {
  d <- bd_bc$distances[group == grp]
  row <- data.frame(Study_Group = grp,
    BC_Mean = round(mean(d), 4), BC_SD = round(sd(d), 4),
    BC_F = round(F_bc, 4), BC_p = round(p_bc, 4),
    stringsAsFactors = FALSE)
  if (uf_available) {
    d2 <- bd_uf$distances[group_uf == grp]
    row$UF_Mean <- round(mean(d2), 4)
    row$UF_SD   <- round(sd(d2), 4)
    row$UF_F    <- round(F_uf, 4)
    row$UF_p    <- round(p_uf, 4)
  }
  row
})
sum_df <- do.call(rbind, sum_rows)
write.csv(sum_df, file.path(WORKDIR, "BetaDispersion_summary.csv"), row.names = FALSE)
cat("  Saved: BetaDispersion_summary.csv\n")

cat("\n--- Summary ---\n")
print(sum_df, row.names = FALSE)

# ---------------------------------------------------------------------------
# 5. Plot
# ---------------------------------------------------------------------------
cat("\n[5] Plotting...\n")

plot_df <- res_df %>%
  mutate(Group_Label = ifelse(Study_Group == "Space_Flight",
                              "Space Flight", "Terrestrial Soil"),
         Group_Label = factor(Group_Label, levels = c("Space Flight", "Terrestrial Soil")))

# How many panels?
n_panels <- if (uf_available) 2 else 1
fig_w    <- if (uf_available) 10 else 5.5

# Significance label helper
sig_label <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "n.s."
}

# Bray-Curtis plot
p1 <- ggplot(plot_df, aes(x = Group_Label, y = BC_Distance, fill = Group_Label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.4, color = "gray30") +
  scale_fill_manual(values = c("Space Flight"    = "#3498db",
                                "Terrestrial Soil" = "#2ecc71")) +
  annotate("text", x = 1.5, y = max(plot_df$BC_Distance) * 1.05,
           label = sprintf("F = %.2f, p = %.3f %s",
                           F_bc, p_bc, sig_label(p_bc)),
           size = 3.8, hjust = 0.5) +
  labs(title = "Bray-Curtis Beta-Dispersion",
       subtitle = sprintf("PERMDISP2, %d permutations", N_PERM),
       x = NULL, y = "Distance to Group Centroid") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

if (uf_available) {
  plot_df_uf <- data.frame(
    SampleID   = names(bd_uf$distances),
    Study_Group = as.character(group_uf),
    UF_Distance = bd_uf$distances,
    stringsAsFactors = FALSE
  ) %>%
    mutate(Group_Label = ifelse(Study_Group == "Space_Flight",
                                "Space Flight", "Terrestrial Soil"),
           Group_Label = factor(Group_Label, levels = c("Space Flight", "Terrestrial Soil")))

  p2 <- ggplot(plot_df_uf, aes(x = Group_Label, y = UF_Distance, fill = Group_Label)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.4, color = "gray30") +
    scale_fill_manual(values = c("Space Flight"    = "#3498db",
                                  "Terrestrial Soil" = "#2ecc71")) +
    annotate("text", x = 1.5, y = max(plot_df_uf$UF_Distance) * 1.05,
             label = sprintf("F = %.2f, p = %.3f %s",
                             F_uf, p_uf, sig_label(p_uf)),
             size = 3.8, hjust = 0.5) +
    labs(title = "UniFrac Beta-Dispersion",
         subtitle = sprintf("PERMDISP2, %d permutations", N_PERM),
         x = NULL, y = "Distance to Group Centroid") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))
}

# Combine and save
if (requireNamespace("patchwork", quietly = TRUE) && uf_available) {
  library(patchwork)
  fig <- p1 + p2 + plot_annotation(
    title = "Beta-Dispersion: Community Compositional Variability",
    subtitle = "Higher dispersion = more variable (less predictable) community within each group",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )
} else {
  fig <- p1
}

ggsave(file.path(WORKDIR, "Main_Fig_BetaDispersion.png"), fig,
       width = fig_w, height = 5, dpi = 300)
ggsave(file.path(WORKDIR, "Main_Fig_BetaDispersion.pdf"), fig,
       width = fig_w, height = 5)
cat("  Saved: Main_Fig_BetaDispersion.png\n")
cat("  Saved: Main_Fig_BetaDispersion.pdf\n")

cat("\n============================================================\n")
cat("Beta-dispersion analysis complete.\n")
cat("============================================================\n")
