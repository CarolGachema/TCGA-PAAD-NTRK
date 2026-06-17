
# TCGA-PAAD : NTRK1 (TrkA) / NTRK2 (TrkB) expression, co-expression,
#             survival association, and pathway enrichment
# Data source: UCSC Xena Browser, TCGA-PAAD cohort
#   Expression : HiSeqV2 (IlluminaHiSeq RNAseqV2, gene-level,
#                log2(normalised_count + 1))  — genes × samples
#   Survival   : TCGA-PAAD survival curated file
#                (columns: sample, OS, OS.time, ...)


#Packages 
cran_pkgs <- c("data.table", "dplyr", "tidyr", "tibble",
               "ggplot2", "ggpubr", "ggrepel",
               "survival", "survminer")

new_cran <- cran_pkgs[!cran_pkgs %in% installed.packages()[, "Package"]]
if (length(new_cran)) install.packages(new_cran)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot")
new_bioc  <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[, "Package"]]
if (length(new_bioc))
  BiocManager::install(new_bioc, update = FALSE, ask = FALSE)

library(data.table); library(dplyr);  library(tidyr);    library(tibble)
library(ggplot2);    library(ggpubr); library(ggrepel)
library(survival);   library(survminer)
library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)


# Loading data
survival_raw <- read.delim(
  "data/PAAD_survival.txt",
  check.names = FALSE
)

expr <- read.delim(
  "data/gene_exp.tsv",
  check.names = FALSE
)

expr <- read.delim(
  "data/TCGA_matrix.gz",
  check.names = FALSE
)


message("Expression matrix: ", nrow(expr), " genes × ", ncol(expr) - 1, " samples")
message("Survival rows: ", nrow(survival_raw))


# Output directory + reproducibility seed

if (!dir.exists("figures")) dir.create("figures")
set.seed(42)   # makes GSEA permutations reproducible


# Identify gene column & extract NTRK1 / NTRK2


# Direct assignment since your file is already sample × gene
ntrk_wide <- expr

# Dynamically identify and standardize the sample column name
gene_col <- colnames(ntrk_wide)[1]
colnames(ntrk_wide)[colnames(ntrk_wide) == gene_col] <- "sample"

message("Columns detected: ", paste(colnames(ntrk_wide), collapse = " | "))

# Trim barcodes to standard 15-character primary tumor format
ntrk_wide$sample <- substr(ntrk_wide$sample, 1, 15)


# Safe coercion to numeric values for downstream statistical models
if ("NTRK1" %in% colnames(ntrk_wide)) ntrk_wide$NTRK1 <- as.numeric(ntrk_wide$NTRK1)
if ("NTRK2" %in% colnames(ntrk_wide)) ntrk_wide$NTRK2 <- as.numeric(ntrk_wide$NTRK2)

# Verify data integrity
message("Actual patients in expression data: ", nrow(ntrk_wide))
print(head(ntrk_wide[, c("sample", "NTRK1", "NTRK2")]))


# Tidy survival data + merge

# Force all column names to lowercase to eliminate Xena casing variations
names(survival_raw) <- tolower(names(survival_raw))

# Handle alternate Xena underscore prefixes if they exist
if (!"os" %in% names(survival_raw) && "_os" %in% names(survival_raw)) {
  names(survival_raw)[names(survival_raw) == "_os"] <- "os"
}
if (!"os.time" %in% names(survival_raw) && "_os_time" %in% names(survival_raw)) {
  names(survival_raw)[names(survival_raw) == "_os_time"] <- "os.time"
}

# Clean and isolate my endpoints (using explicit namespace to avoid package conflicts)
surv_clean <- survival_raw %>%
  mutate(sample = substr(sample, 1, 15)) %>%
  dplyr::select(sample, OS = os, OS.time = os.time) %>%   # Forces R to use dplyr's select
  filter(
    !is.na(OS),
    !is.na(OS.time),
    OS.time > 0                                           # removes zero-time entries
  )

# Inner join: matches your already-wide expression data with your clean survival endpoints
df <- inner_join(ntrk_wide, surv_clean, by = "sample") %>%
  filter(!is.na(NTRK1), !is.na(NTRK2))

message("Samples after merging expression + survival: ", nrow(df))


# Sample size safety check
if (nrow(df) < 30) {
  warning(
    "Very few samples matched (n = ", nrow(df), "). ",
    "Check that sample IDs in both files use the same TCGA barcode format. ",
    "Try substr(sample, 1, 12) in both files as an alternative trim length."
  )
}


# Median dichotomisation (High / Low)

med_ntrk1 <- median(df$NTRK1, na.rm = TRUE)
med_ntrk2 <- median(df$NTRK2, na.rm = TRUE)

df <- df %>%
  mutate(
    # Low is the reference level so that HR > 1 means High = worse prognosis
    NTRK1_group = factor(
      ifelse(NTRK1 >= med_ntrk1, "High", "Low"),
      levels = c("Low", "High")
    ),
    NTRK2_group = factor(
      ifelse(NTRK2 >= med_ntrk2, "High", "Low"),
      levels = c("Low", "High")
    )
  )

message(
  "NTRK1 — median: ", round(med_ntrk1, 3),
  " | High: ", sum(df$NTRK1_group == "High"),
  " | Low: ", sum(df$NTRK1_group == "Low")
)
message(
  "NTRK2 — median: ", round(med_ntrk2, 3),
  " | High: ", sum(df$NTRK2_group == "High"),
  " | Low: ", sum(df$NTRK2_group == "Low")
)


# Publication theme + save helper

# Colour palette — consistent throughout all figures
col_high <- "#C0392B"   # deep red    = high expression / poor outcome
col_low  <- "#2C6BAC"   # steel blue  = low expression
col_pt   <- "#4D4D4D"   # dark grey   = scatter points
col_line <- "#E8844A"   # orange      = regression line

pub_theme <- theme_classic(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle    = element_text(size = 10.5, hjust = 0.5, colour = "grey45"),
    axis.title       = element_text(face = "bold", size = 12),
    axis.text        = element_text(size = 11, colour = "grey20"),
    legend.title     = element_text(face = "bold", size = 11),
    legend.text      = element_text(size = 10),
    strip.text       = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.4)
  )

# Save a ggplot object as both PDF and high-res PNG
save_fig <- function(p, name, w = 7, h = 5.5) {
  ggsave(file.path("figures", paste0(name, ".pdf")),
         plot = p, width = w, height = h, device = "pdf")
  ggsave(file.path("figures", paste0(name, ".png")),
         plot = p, width = w, height = h, dpi = 300, device = "png")
  message("  Saved → figures/", name, ".pdf & .png")
}

# Save a survminer ggsurvplot (needs explicit device handling)
save_km <- function(km_obj, name, w = 8, h = 6.5) {
  pdf(file.path("figures", paste0(name, ".pdf")), width = w, height = h)
  print(km_obj)
  dev.off()
  png(file.path("figures", paste0(name, ".png")),
      width = w, height = h, units = "in", res = 300)
  print(km_obj)
  dev.off()
  message("  Saved → figures/", name, ".pdf & .png")
}


# FIGURE 1 — Expression distributions
# Violin + boxplot + jitter for NTRK1 and NTRK2

message("\n[ Figure 1 ] Expression distributions …")

# Reshape to long format for facetting
df_violin <- df %>%
  dplyr::select(sample, NTRK1, NTRK2) %>%   # Added dplyr:: right here!
  pivot_longer(
    cols      = c(NTRK1, NTRK2),
    names_to  = "Gene",
    values_to = "Expression"
  ) %>%
  mutate(Gene = factor(Gene, levels = c("NTRK1", "NTRK2")))
# Per-gene medians for annotation
medians_df <- df_violin %>%
  group_by(Gene) %>%
  summarise(med = median(Expression, na.rm = TRUE), .groups = "drop")


fig1 <- ggplot(df_violin, aes(x = Gene, y = Expression, fill = Gene)) +
  geom_violin(trim = FALSE, alpha = 0.55, colour = NA, width = 0.9) +
  geom_boxplot(
    width        = 0.14,
    outlier.shape = NA,
    colour       = "black",
    linewidth    = 0.75,
    fill         = "white",
    alpha        = 0.85
  ) +
  geom_jitter(
    width   = 0.09,
    size    = 0.9,
    alpha   = 0.30,
    colour  = col_pt
  ) +
  geom_text(
    data     = medians_df,
    aes(x = Gene, y = med, label = paste0("Median: ", round(med, 2))),
    nudge_x  = 0.30,
    size     = 3.4,
    fontface = "italic",
    colour   = "grey30",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(NTRK1 = "#3B82C4", NTRK2 = col_high),
    guide  = "none"
  ) +
  labs(
    title    = expression(italic("NTRK1") ~ "and" ~ italic("NTRK2") ~
                            "Expression in TCGA-PAAD"),
    subtitle = paste0(
      "Pancreatic ductal adenocarcinoma  |  n = ", nrow(df), " tumour samples"
    ),
    x = NULL,
    y = expression(log[2] ~ "(normalised count + 1)")
  ) +
  pub_theme

print(fig1)
save_fig(fig1, "Figure1_expression_distribution", w = 6, h = 5.5)


#FIGURE 2: NTRK1 vs NTRK2 co-expression scatter

message("\n[ Figure 2 ] Co-expression scatter …")

# Pearson correlation test
cor_test <- cor.test(df$NTRK1, df$NTRK2, method = "pearson")
r_val    <- round(cor_test$estimate, 3)
p_raw    <- cor_test$p.value
p_label  <- if (p_raw < 0.001) {
  paste0("Pearson r = ", r_val, ",  p < 0.001")
} else {
  paste0("Pearson r = ", r_val, ",  p = ", signif(p_raw, 3))
}

# Label the 6 samples with the highest combined NTRK1 + NTRK2 expression
top_samples <- df %>%
  mutate(combined = NTRK1 + NTRK2) %>%
  slice_max(combined, n = 6) %>%
  mutate(short_id = substr(sample, 6, 12))   # e.g. "HZ-7289" — readable

fig2 <- ggplot(df, aes(x = NTRK1, y = NTRK2)) +
  geom_point(
    size   = 2.2,
    alpha  = 0.60,
    colour = col_pt
  ) +
  geom_smooth(
    method   = "lm",
    colour   = col_line,
    fill     = "#F5CBA7",
    se       = TRUE,
    linewidth = 1.1,
    alpha    = 0.30
  ) +
  geom_label_repel(
    data         = top_samples,
    aes(label = short_id),
    size         = 2.8,
    colour       = col_high,
    fill         = "white",
    label.padding = unit(0.15, "lines"),
    box.padding  = unit(0.35, "lines"),
    max.overlaps = Inf,
    seed         = 42
  ) +
  annotate(
    "text",
    x        = min(df$NTRK1, na.rm = TRUE),
    y        = max(df$NTRK2, na.rm = TRUE),
    label    = p_label,
    size     = 4,
    fontface = "italic",
    hjust    = 0,
    vjust    = 1,
    colour   = "grey25"
  ) +
  labs(
    title    = expression(italic("NTRK1") ~ "–" ~ italic("NTRK2") ~
                            "Co-expression  |  TCGA-PAAD"),
    subtitle = "Each point = one tumour sample. Labels: top 6 co-expressing samples.",
    x        = expression(italic("NTRK1") ~ "expression  " *
                            log[2] * "(norm + 1)"),
    y        = expression(italic("NTRK2") ~ "expression  " *
                            log[2] * "(norm + 1)")
  ) +
  pub_theme

print(fig2)
save_fig(fig2, "Figure2_coexpression_scatter", w = 6.5, h = 5.5)


# FIGURE 3 — Kaplan-Meier survival curves
# NTRK2 high vs low  (primary — most directly supports paper thesis)
# 3b  NTRK1 high vs low  (secondary)
# Both include:
#  Log-rank p-value
#  HR (95 % CI) from Cox proportional-hazards
#  Median survival lines
#  Risk table


message("\n[ Figure 3 ] Kaplan-Meier curves …")

km_palette <- c("Low" = col_low, "High" = col_high)

# Helper: extract HR and CI from a fitted coxph model for a binary covariate
get_hr <- function(cox_fit, term) {
  hr <- round(exp(coef(cox_fit)[term]), 2)
  ci <- round(exp(confint(cox_fit))[term, ], 2)
  list(
    hr    = hr,
    lower = ci[1],
    upper = ci[2],
    label = paste0("HR = ", hr, "  (95% CI ", ci[1], "–", ci[2], ")")
  )
}

# NTRK2

fit_ntrk2 <- survfit(
  Surv(OS.time / 30.44, OS) ~ NTRK2_group,   # days-months
  data = df
)

cox_ntrk2 <- coxph(Surv(OS.time, OS) ~ NTRK2_group, data = df)
hr2       <- get_hr(cox_ntrk2, "NTRK2_groupHigh")

lr_p2 <- survdiff(Surv(OS.time, OS) ~ NTRK2_group, data = df)
lr_p2 <- round(1 - pchisq(lr_p2$chisq, df = 1), 4)

fig3a <- ggsurvplot(
  fit_ntrk2,
  data              = df,
  palette           = c(col_low, col_high),      # unnamed, ordered Low then High
  pval              = TRUE,
  pval.size         = 4.2,
  pval.coord        = c(2, 0.08),
  conf.int          = FALSE,
  risk.table        = TRUE,
  risk.table.height = 0.26,
  risk.table.fontsize = 3.6,
  risk.table.col    = "strata",
  legend.labs       = c(
    paste0("Low NTRK2  (n=", sum(df$NTRK2_group=="Low"),")"),
    paste0("High NTRK2 (n=", sum(df$NTRK2_group=="High"),")")
  ),
  legend.title      = "",
  xlab              = "Time (months)",
  ylab              = "Overall survival probability",
  title             = expression("Overall Survival by" ~ italic("NTRK2") ~
                                   "Expression — TCGA-PAAD"),
  subtitle          = hr2$label,
  surv.median.line  = "hv",
  ggtheme           = theme_classic(base_size = 13),
  fontsize          = 4
)

print(fig3a)
save_km(fig3a, "Figure3a_KM_NTRK2", w = 8, h = 6.5)

# NTRK1

fit_ntrk1 <- survfit(
  Surv(OS.time / 30.44, OS) ~ NTRK1_group,
  data = df
)

cox_ntrk1 <- coxph(Surv(OS.time, OS) ~ NTRK1_group, data = df)
hr1       <- get_hr(cox_ntrk1, "NTRK1_groupHigh")

lr_p1 <- survdiff(Surv(OS.time, OS) ~ NTRK1_group, data = df)
lr_p1 <- round(1 - pchisq(lr_p1$chisq, df = 1), 4)

fig3b <- ggsurvplot(
  fit_ntrk1,
  data              = df,
  palette           = c(col_low, col_high),
  pval              = TRUE,
  pval.size         = 4.2,
  pval.coord        = c(2, 0.08),
  conf.int          = FALSE,
  risk.table        = TRUE,
  risk.table.height = 0.26,
  risk.table.fontsize = 3.6,
  risk.table.col    = "strata",
  legend.labs       = c(
    paste0("Low NTRK1  (n=", sum(df$NTRK1_group=="Low"),")"),
    paste0("High NTRK1 (n=", sum(df$NTRK1_group=="High"),")")
  ),
  legend.title      = "",
  xlab              = "Time (months)",
  ylab              = "Overall survival probability",
  title             = expression("Overall Survival by" ~ italic("NTRK1") ~
                                   "Expression — TCGA-PAAD"),
  subtitle          = hr1$label,
  surv.median.line  = "hv",
  ggtheme           = theme_classic(base_size = 13),
  fontsize          = 4
)

print(fig3b)
save_km(fig3b, "Figure3b_KM_NTRK1", w = 8, h = 6.5)


# FIGURE 4 — GSEA-KEGG pathway enrichment
# Compute Spearman ρ between every gene and NTRK2 across all samples
# Build a ranked list (high ρ → positively co-expressed with NTRK2)
# Run gseKEGG on that ranked list
# Dot plot of top pathways, with our 4 key pathways highlighted
# Individual enrichment plots for those 4 pathways


message("\n[ Figure 4 ] GSEA pathway enrichment …")

# Build a sample-aligned NTRK2 vector


#Figure 4: 
# Extract the gene symbols from the very first column
gene_symbols <- expr[[1]]

# Locate where NTRK2 lives in the rows
ntrk2_idx <- which(gene_symbols == "NTRK2")
if (length(ntrk2_idx) == 0) {
  stop("NTRK2 not found in the expression matrix rows!")
}

# Grab the patient columns (skipping the first column) and trim to 15 chars
expr_patients <- colnames(expr)[-1]
expr_patients_trimmed <- substr(expr_patients, 1, 15)

# Find which columns match your clean survival dataset (df)
matched_cols_idx <- which(expr_patients_trimmed %in% df$sample)

if (length(matched_cols_idx) < 10) {
  stop("Fewer than 10 samples matched. Check your sample ID structures.")
}

# Isolate the data for matched patients (+1 to account for skipping column 1)
sub_matrix   <- expr[, matched_cols_idx + 1, drop = FALSE]
ntrk2_vector <- as.numeric(expr[ntrk2_idx, matched_cols_idx + 1])

message("Computing Spearman correlations row-by-row for ", nrow(sub_matrix), " genes...")

# Run a fast row-by-row correlation against your NTRK2 vector
gene_rho <- apply(sub_matrix, 1, function(gene_vec) {
  suppressWarnings(
    cor(as.numeric(gene_vec), ntrk2_vector, method = "spearman", use = "complete.obs")
  )
})

# Assign gene names, filter out NAs, and sort descending for GSEA ranking
names(gene_rho) <- gene_symbols
gene_rho        <- sort(gene_rho[!is.na(gene_rho)], decreasing = TRUE)

# 
message("  Genes in ranked list: ", length(gene_rho))



# [ Figure 4 ] Complete GSEA Pipeline (Tall Matrix Optimized)

# Compute gene-level Spearman correlations with NTRK2

# Extract the gene symbols from the very first column of your tall matrix
gene_symbols <- expr[[1]]

# Locate where NTRK2 lives in the rows
ntrk2_idx <- which(gene_symbols == "NTRK2")
if (length(ntrk2_idx) == 0) {
  stop("NTRK2 not found in the expression matrix rows!")
}

# Grab the patient columns (skipping column 1) and trim barcodes to 15 chars
expr_patients <- colnames(expr)[-1]
expr_patients_trimmed <- substr(expr_patients, 1, 15)

# Find which columns match your clean survival dataset (df)
matched_cols_idx <- which(expr_patients_trimmed %in% df$sample)

if (length(matched_cols_idx) < 10) {
  stop("Fewer than 10 samples matched. Check your sample ID structures.")
}

# Isolate data for matched patients (+1 accounts for skipping column 1)
sub_matrix   <- expr[, matched_cols_idx + 1, drop = FALSE]
ntrk2_vector <- as.numeric(expr[ntrk2_idx, matched_cols_idx + 1])

message("Computing Spearman correlations row-by-row for ", nrow(sub_matrix), " genes...")

#Run a fast row-by-row correlation against your NTRK2 vector
gene_rho <- apply(sub_matrix, 1, function(gene_vec) {
  suppressWarnings(
    cor(as.numeric(gene_vec), ntrk2_vector, method = "spearman", use = "complete.obs")
  )
})

# Assign the thousands of gene names, filter out NAs, and sort descending
names(gene_rho) <- gene_symbols
gene_rho        <- sort(gene_rho[!is.na(gene_rho)], decreasing = TRUE)


message("  Genes in ranked list: ", length(gene_rho))


# Convert SYMBOL → Entrez IDs

id_map <- bitr(
  names(gene_rho),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

# One Entrez ID per gene symbol; re-rank
gene_list_df <- data.frame(
  SYMBOL = names(gene_rho),
  rho    = gene_rho,
  stringsAsFactors = FALSE
) %>%
  inner_join(id_map, by = "SYMBOL") %>%
  arrange(desc(rho)) %>%
  distinct(ENTREZID, .keep_all = TRUE)   # drop duplicate Entrez mappings

gene_list_entrez <- setNames(gene_list_df$rho, gene_list_df$ENTREZID)
message("  Genes with Entrez IDs: ", length(gene_list_entrez))


#. Run GSEA

set.seed(42)
message("  Running gseKEGG …  (~1-3 min)")

gsea_res <- gseKEGG(
  geneList      = gene_list_entrez,
  organism      = "hsa",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.25,     # permissive for a full picture; tighten to 0.05
  pAdjustMethod = "BH",     # Benjamini-Hochberg FDR
  verbose       = FALSE,
  eps           = 0         # exact p-values (slower but more accurate)
)

gsea_df <- as.data.frame(gsea_res)
message("  Pathways at padj < 0.25: ", nrow(gsea_df))

# Our 4 mechanistically critical pathways
target_pathways <- c(
  "hsa04722",   # Neurotrophin signalling
  "hsa04360",   # Axon guidance
  "hsa04010",   # MAPK signalling
  "hsa04151"    # PI3K-Akt signalling
)

gsea_df <- gsea_df %>%
  mutate(
    highlight  = ID %in% target_pathways,
    Direction  = ifelse(NES > 0, "Enriched", "Depleted")
  ) %>%
  arrange(desc(abs(NES)))


## Dot plot

# Take top 20 pathways by |NES|
plot_df <- gsea_df %>%
  slice_max(order_by = abs(NES), n = pmin(20, nrow(gsea_df))) %>%
  arrange(NES) %>%
  mutate(Description = factor(Description, levels = Description))

fig4 <- ggplot(
  plot_df,
  aes(x = NES, y = Description, colour = p.adjust, size = setSize)
) +
  geom_point(alpha = 0.85) +
  
  # Orange ring around our 4 key pathways
  geom_point(
    data       = subset(plot_df, ID %in% target_pathways),
    shape      = 21,
    stroke     = 2,
    fill       = NA,
    colour     = "#E67E22",
    show.legend = FALSE
  ) +
  # Label for the ring legend
  annotate(
    "text",
    x      = min(plot_df$NES, na.rm = TRUE),
    y      = 1,
    label  = "○  Key pathway (neurotrophin / axon / MAPK / PI3K)",
    size   = 3.2,
    colour = "#E67E22",
    hjust  = 0
  ) +
  scale_colour_gradient(
    low   = col_high,
    high  = "#AEB6BF",
    name  = "Adjusted\np-value"
  ) +
  scale_size_continuous(
    range  = c(3, 9),
    name   = "Gene set\nsize"
  ) +
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    colour     = "grey55",
    linewidth  = 0.6
  ) +
  labs(
    title    = expression(
      "GSEA-KEGG: Pathways co-expressed with" ~ italic("NTRK2") ~ "in TCGA-PAAD"
    ),
    subtitle = "Genes ranked by Spearman ρ with NTRK2 expression",
    x        = "Normalised Enrichment Score (NES)",
    y        = NULL
  ) +
  pub_theme +
  theme(axis.text.y = element_text(size = 8.5))

print(fig4)
save_fig(fig4, "Figure4_GSEA_KEGG_dotplot", w = 9.5, h = 7)
## ── 9f. Individual enrichment plots for the 4 key pathways ──────────────────

detected_targets <- intersect(target_pathways, gsea_df$ID)

if (length(detected_targets) == 0) {
  message(
    "  NOTE: none of the 4 key pathways reached the current padj threshold.\n",
    "  Try relaxing pvalueCutoff to 0.50 in Section 9d, then re-run.\n",
    "  This is common with small cohort sizes (PAAD n ≈ 150-180)."
  )
} else {
  for (pid in detected_targets) {
    pname <- gsea_df$Description[gsea_df$ID == pid]
    p_enr <- gseaplot2(
      gsea_res,
      geneSetID = pid,
      title     = pname,
      color     = col_high,
      base_size = 12
    )
    fname <- paste0("Figure4_enrichplot_", pid)
    pdf(file.path("figures", paste0(fname, ".pdf")), width = 7, height = 5)
    print(p_enr)
    dev.off()
    png(file.path("figures", paste0(fname, ".png")),
        width = 7, height = 5, units = "in", res = 300)
    print(p_enr)
    dev.off()
    message("  Saved → figures/", fname, ".pdf & .png")
  }
}


# 3. Save the full data tables inside the results folder
if (exists("gsea_df") && nrow(gsea_df) > 0) {
  write.csv(gsea_df, "results/GSEA_full_pathway_results.csv", row.names = FALSE)
}
if (exists("df") && nrow(df) > 0) {
  write.csv(df, "results/matched_clinical_expression_data.csv", row.names = FALSE)
}
if (exists("gene_rho") && length(gene_rho) > 0) {
  rho_df <- data.frame(Gene = names(gene_rho), Spearman_Rho = as.numeric(gene_rho))
  write.csv(rho_df, "results/NTRK2_genome_wide_correlations.csv", row.names = FALSE)
}
