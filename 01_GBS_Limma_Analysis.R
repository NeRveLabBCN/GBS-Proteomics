# ==============================================================================
# Script Name: 01_GBS_Limma_Analysis.R
# Description: Differential abundance analysis of Guillain-Barré Syndrome (GBS) 
#              proteomics data using limma. Includes:
#              1. Paired (Acute vs Recovery)
#              2. Covariate-adjusted (Acute vs Healthy Controls)
#              3. Covariate-adjusted (Recovery vs Healthy Controls)
# Input:       - Data/GBS_Proteomics_Raw_Github.csv
#              - Data/SomaScan_Protein_Dictionary_Github.csv 
# ==============================================================================

# Set seed for reproducible volcano plot label placement
set.seed(1234)

# Load required packages
library(dplyr)
library(limma)
library(writexl)
library(ggplot2)
library(ggrepel)
library(viridis)

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================

# Read the clean dataset and dictionary. 
gbs_proj <- read.csv("Data/GBS_Proteomics_Raw_Github.csv", check.names = FALSE)
dictionary <- read.csv("Data/SomaScan_Protein_Dictionary_Github.csv")

# Identify metadata and apply log2 transformation to protein columns
meta_cols <- c("SampleId", "baseline_hc", "baseline_remission", "matching_ID", "SampleType", "Group", "BaseID", "Age", "Sex")
prot_cols <- setdiff(colnames(gbs_proj), meta_cols)
gbs_proj[, prot_cols] <- log2(gbs_proj[, prot_cols])

# Global Thresholds
fdr_cut <- 0.05
lfc_cut <- 1  

# Note: Take into account that this will identidy differentially abundant aptamers, and there may be more than one aptamer per protein.

# ==============================================================================
# 2. ANALYSIS 1: Acute vs Recovery Guillain-Barré syndrome (paired)
# ==============================================================================

cat("\n--- Running Acute vs Recovery Analysis ---\n")

# 1) Subset to the two groups
dat_AR <- gbs_proj[gbs_proj$Group %in% c("GBS_Acute", "GBS_Recovery"), ]
dat_AR$Group <- droplevels(factor(dat_AR$Group, levels = c("GBS_Acute", "GBS_Recovery")))

# 2) Expression matrix: features (rows) × samples (cols)
expr_AR <- t(as.matrix(dat_AR[, prot_cols]))

# 3) Model design
design_AR <- model.matrix(~ 0 + Group, data = dat_AR)
colnames(design_AR) <- sub("^Group", "", colnames(design_AR))

# Identify paired samples via matching_ID
block_ids <- dat_AR$matching_ID

# 4) Estimate within-patient correlation and fit the linear model
corfit <- duplicateCorrelation(expr_AR, design_AR, block = block_ids)
fit_AR <- lmFit(expr_AR, design_AR, block = block_ids, correlation = corfit$consensus)
C_AR <- makeContrasts(GBS_Acute - GBS_Recovery, levels = design_AR)
fit2_AR <- contrasts.fit(fit_AR, C_AR)
fit2_AR <- eBayes(fit2_AR, trend = TRUE, robust = TRUE)
res_AR <- topTable(fit2_AR, number = Inf, adjust.method = "BH")

# 5) Annotate, reorder columns, save CSV
res_AR_annot <- res_AR %>%
  tibble::rownames_to_column("AptName") %>%
  left_join(dictionary, by = "AptName") %>%
  dplyr::select(Target, AptName, TargetFullName, EntrezGeneSymbol, UniProt, logFC, AveExpr, t, P.Value, adj.P.Val, B)

write.csv(res_AR_annot, "limma_acute_vs_recovery.csv", row.names = FALSE)

# 6) Count differentially abundant aptamers
sig_AR <- res_AR_annot %>% filter(adj.P.Val < fdr_cut & abs(logFC) >= lfc_cut)
cat("Acute vs Recovery -> Upregulated:", sum(sig_AR$logFC >= lfc_cut), 
    "| Downregulated:", sum(sig_AR$logFC <= -lfc_cut), 
    "| Total:", nrow(sig_AR), "\n")

# 7) Volcano Plot
df_AR <- res_AR_annot
df_AR$neg_log10_fdr <- -log10(pmax(df_AR$adj.P.Val, .Machine$double.xmin))
top_fc_AR <- df_AR[order(-abs(df_AR$logFC)), ][1:30, ]
x_limit_AR <- max(abs(df_AR$logFC), na.rm = TRUE) * 1.05

vp_AR <- ggplot(df_AR, aes(x = logFC, y = neg_log10_fdr)) +
  geom_point(aes(color = adj.P.Val), alpha = 0.8) +
  scale_color_viridis(option = "plasma", trans = "log10", direction = -1, name = "FDR") +
  scale_x_continuous(limits = c(-x_limit_AR, x_limit_AR)) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", color = "black") +
  labs(title = "Acute GBS vs Recovery GBS", x = "Log2 Fold-Change", y = "-log10 FDR") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5),
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.80),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    axis.text = element_text(color = "black")
  ) +
  geom_text_repel(
    data = top_fc_AR,
    aes(label = ifelse(is.na(EntrezGeneSymbol) | EntrezGeneSymbol == "", AptName, EntrezGeneSymbol)),
    size = 3.5, max.overlaps = Inf
  )
print(vp_AR)


# ==============================================================================
# 3. ANALYSIS 2: Acute Guillain-Barré syndrome vs Healthy controls (Covariate-adjusted)
# ==============================================================================

cat("\n--- Running Acute vs Healthy Controls Analysis ---\n")

# 1) Subset to the two groups
dat_AC <- gbs_proj[gbs_proj$Group %in% c("GBS_Acute","Healthy_Control"), ]
dat_AC$Group <- droplevels(factor(dat_AC$Group, levels = c("GBS_Acute","Healthy_Control")))

# Ensure covariates are in the right form (Male baseline, Age numeric)
dat_AC$Sex <- droplevels(factor(dat_AC$Sex, levels = c("male","female")))
stopifnot(is.numeric(dat_AC$Age))

# 2) Expression matrix: features (rows) × samples (cols)
expr_AC <- t(as.matrix(dat_AC[, prot_cols]))

# 3) Design matrix (no intercept, covariates included)
design_AC <- model.matrix(~ 0 + Group + Age + Sex, data = dat_AC)
colnames(design_AC) <- sub("^Group", "", colnames(design_AC))

# 4) Fit limma model and contrast: Acute - Healthy
fit_AC <- lmFit(expr_AC, design_AC)
C_AC <- makeContrasts(GBS_Acute - Healthy_Control, levels = design_AC)
fit2_AC <- contrasts.fit(fit_AC, C_AC)
fit2_AC <- eBayes(fit2_AC, trend = TRUE, robust = TRUE)
res_AC <- topTable(fit2_AC, number = Inf, adjust.method = "BH")

# 5) Annotate, reorder columns, save CSV
res_AC_annot <- res_AC %>%
  tibble::rownames_to_column("AptName") %>%
  left_join(dictionary, by = "AptName") %>%
  dplyr::select(Target, AptName, TargetFullName, EntrezGeneSymbol, UniProt, logFC, AveExpr, t, P.Value, adj.P.Val, B)

write.csv(res_AC_annot, "limma_acute_vs_healthy.csv", row.names = FALSE)

# 6) Count differentially abundant aptamers
sig_AC <- res_AC_annot %>% filter(adj.P.Val < fdr_cut & abs(logFC) >= lfc_cut)
cat("Acute vs HC -> Upregulated:", sum(sig_AC$logFC >= lfc_cut), 
    "| Downregulated:", sum(sig_AC$logFC <= -lfc_cut), 
    "| Total:", nrow(sig_AC), "\n")

# 7) Volcano Plot
df_AC <- res_AC_annot
df_AC$neg_log10_fdr <- -log10(pmax(df_AC$adj.P.Val, .Machine$double.xmin))
top_fc_AC <- df_AC[order(-abs(df_AC$logFC)), ][1:30, ]
x_limit_AC <- max(abs(df_AC$logFC), na.rm = TRUE) * 1.05

vp_AC <- ggplot(df_AC, aes(x = logFC, y = neg_log10_fdr)) +
  geom_point(aes(color = adj.P.Val), alpha = 0.85) +
  scale_color_viridis(option = "plasma", trans = "log10", direction = -1, name = "FDR") +
  scale_x_continuous(limits = c(-x_limit_AC, x_limit_AC)) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", color = "black") +
  labs(title = "Acute GBS vs Healthy Controls", x = "Log2 Fold-Change", y = "-log10 FDR") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5),
    legend.position = "inside",
    legend.position.inside = c(0.10, 0.80),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    axis.text = element_text(color = "black")
  ) +
  geom_text_repel(
    data = top_fc_AC,
    aes(label = ifelse(is.na(EntrezGeneSymbol) | EntrezGeneSymbol == "", AptName, EntrezGeneSymbol)),
    size = 3.5, max.overlaps = Inf
  )
print(vp_AC)


# ==============================================================================
# 4. ANALYSIS 3: Recovery Guillain-Barré syndrome VS Healthy controls (covariate-adjusted)
# ==============================================================================

cat("\n--- Running Recovery vs Healthy Controls Analysis ---\n")

# 1) Subset to the two groups
dat_RC <- gbs_proj[gbs_proj$Group %in% c("GBS_Recovery","Healthy_Control"), ]
dat_RC$Group <- droplevels(factor(dat_RC$Group, levels = c("GBS_Recovery","Healthy_Control")))

# Ensure covariates are in the right form (male baseline, Age numeric)
dat_RC$Sex  <- droplevels(factor(dat_RC$Sex, levels = c("male","female")))
stopifnot(is.numeric(dat_RC$Age))

# 2) Expression matrix: features (rows) × samples (cols)
expr_RC <- t(as.matrix(dat_RC[, prot_cols])) 

# 3) Design matrix (no intercept, covariates included)
design_RC <- model.matrix(~ 0 + Group + Age + Sex, data = dat_RC)
colnames(design_RC) <- sub("^Group", "", colnames(design_RC))

# 4) Fit limma model and contrast: Recovery - Healthy
fit_RC <- lmFit(expr_RC, design_RC)
C_RC <- makeContrasts(GBS_Recovery - Healthy_Control, levels = design_RC)
fit2_RC <- contrasts.fit(fit_RC, C_RC)
fit2_RC <- eBayes(fit2_RC, trend = TRUE, robust = TRUE)
res_RC <- topTable(fit2_RC, number = Inf, adjust.method = "BH")

# 5) Annotate, reorder columns, save CSV
res_RC_annot <- res_RC %>%
  tibble::rownames_to_column("AptName") %>%
  left_join(dictionary, by = "AptName") %>%
  dplyr::select(Target, AptName, TargetFullName, EntrezGeneSymbol, UniProt, logFC, AveExpr, t, P.Value, adj.P.Val, B)

write.csv(res_RC_annot, file = "limma_recovery_vs_healthy.csv", row.names = FALSE)

# 6) Count differentially abundant aptamers
sig_RC <- res_RC_annot %>% filter(adj.P.Val < fdr_cut & abs(logFC) >= lfc_cut)
cat("Recovery vs HC -> Upregulated:", sum(sig_RC$logFC >= lfc_cut), 
    "| Downregulated:", sum(sig_RC$logFC <= -lfc_cut), 
    "| Total:", nrow(sig_RC), "\n")

# 7) Volcano Plot
df_RC <- res_RC_annot
df_RC$neg_log10_fdr <- -log10(pmax(df_RC$adj.P.Val, .Machine$double.xmin))
top_fc_RC <- df_RC[order(-abs(df_RC$logFC)), ][1:30, ]
x_limit_RC <- max(abs(df_RC$logFC), na.rm = TRUE) * 1.05

vp_RC <- ggplot(df_RC, aes(x = logFC, y = neg_log10_fdr)) +
  geom_point(aes(color = adj.P.Val), alpha = 0.85) +
  scale_color_viridis(option = "plasma", trans = "log10", direction = -1, name = "FDR") +
  scale_x_continuous(limits = c(-x_limit_RC, x_limit_RC)) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", color = "black") +
  labs(title = "Recovery GBS vs Healthy Controls", x = "Log2 Fold-Change", y = "-log10 FDR") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5),
    legend.position = "inside",
    legend.position.inside = c(0.10, 0.80), 
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    axis.text = element_text(color = "black")
  ) +
  geom_text_repel(
    data = top_fc_RC,
    aes(label = ifelse(is.na(EntrezGeneSymbol) | EntrezGeneSymbol == "", AptName, EntrezGeneSymbol)),
    size = 3.5, max.overlaps = Inf
  )
print(vp_RC)
