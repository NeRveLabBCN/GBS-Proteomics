# ==============================================================================
# Script Name: 02_GSEA_Analysis.R
# Description: Pathway enrichment analysis using fgsea. 
#              Automates rank building, GSEA, and dot plot generation for 
#              three comparisons (AR, AC, RC).
# Input:       - limma_acute_vs_recovery.csv
#              - limma_acute_vs_healthy.csv
#              - limma_recovery_vs_healthy.csv
# Output:      - Combined_GSEA_Results.xlsx
#              - Figure_GSEA_Pyramid_PrintReady.tiff
# ==============================================================================

# Set seed for reproducible permutation and plot placement
set.seed(1234)

# Load required packages
library(dplyr)
library(readr)
library(fgsea)
library(msigdbr)
library(ggplot2)
library(stringr)
library(writexl)
library(patchwork)

# ==============================================================================
# 1. LOAD MSIGDB PATHWAYS (HALLMARK, REACTOME, GO:BP)
# ==============================================================================
cat("Loading MSigDB collections (Hallmark, Reactome, and GO:BP)...\n")

hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)

reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol)

gobp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)

toPathways <- function(df) split(df$gene_symbol, df$gs_name)
all_pathways <- c(toPathways(hallmark), toPathways(reactome), toPathways(gobp))


# ==============================================================================
# 2. AUTOMATION FUNCTION FOR GSEA & PLOTTING
# ==============================================================================
run_gsea_and_plot <- function(limma_file, plot_title) {
  cat("\nProcessing:", plot_title, "...\n")
  
  # 1. Load Data & Build Ranks (using moderated t-statistic)
  res <- read_csv(limma_file, show_col_types = FALSE)
  ranks_tbl <- res %>%
    transmute(gene = EntrezGeneSymbol, rank = t) %>%
    filter(!is.na(gene), gene != "", !is.na(rank)) %>%
    group_by(gene) %>%
    slice_max(order_by = abs(rank), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  ranks <- sort(setNames(ranks_tbl$rank, ranks_tbl$gene), decreasing = TRUE)
  
  # 2. Run fgsea
  fgsea_res <- fgseaMultilevel(
    pathways = all_pathways, stats = ranks, minSize = 10, maxSize = 500
  ) %>% as_tibble()
  
  # 3. Format Labels and Filter for Plotting
  df_plot <- fgsea_res %>%
    filter(!is.na(NES), !is.na(padj), !is.na(size)) %>%
    mutate(
      label = pathway %>%
        str_replace("^(HALLMARK_|REACTOME_|GOBP_)", "") %>%
        str_replace_all("_", " ") %>% str_squish() %>% str_to_title(),
      mlog10padj = -log10(pmax(padj, .Machine$double.xmin))
    )
  
  # Get Top 5 Up and Top 5 Down (FDR < 0.05)
  top_up <- df_plot %>% filter(NES > 0, padj < 0.05) %>% arrange(desc(NES)) %>% slice_head(n = 5)
  top_down <- df_plot %>% filter(NES < 0, padj < 0.05) %>% arrange(NES) %>% slice_head(n = 5)
  
  plot_df <- bind_rows(top_up, top_down) %>%
    mutate(label_wrapped = stringr::str_wrap(toupper(label), width = 28)) %>%
    arrange(NES) %>%
    mutate(label_fac = factor(label_wrapped, levels = unique(label_wrapped)))
  
  # 4. Generate Dot Plot
  p <- ggplot(plot_df, aes(x = NES, y = label_fac)) +
    geom_point(aes(size = size, color = mlog10padj), alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_viridis_c(name = "-log10 FDR") +
    scale_size_continuous(name = "Gene set size") +
    labs(title = plot_title, x = "Normalized Enrichment Score", y = NULL) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      axis.text.y = element_text(size = 9, color = "black"), 
      axis.text.x = element_text(size = 10, color = "black"),
      legend.position = "right", aspect.ratio = 3.0
    )
  
  return(list(results = fgsea_res, plot = p))
}

# ==============================================================================
# 3. EXECUTE PIPELINE
# ==============================================================================

# Run the function for all three comparisons
gsea_AR <- run_gsea_and_plot("limma_acute_vs_recovery.csv", "Acute GBS vs Recovery GBS")
gsea_AC <- run_gsea_and_plot("limma_acute_vs_healthy.csv", "Acute GBS vs Healthy Controls")
gsea_RC <- run_gsea_and_plot("limma_recovery_vs_healthy.csv", "Recovery GBS vs Healthy Controls")


# ==============================================================================
# 4. EXPORT RESULTS
# ==============================================================================

# 1. Save Excel Tables
write_xlsx(
  list(
    "Acute_vs_Recovery" = gsea_AR$results %>% arrange(padj), 
    "Acute_vs_HC"       = gsea_AC$results %>% arrange(padj), 
    "Recovery_vs_HC"    = gsea_RC$results %>% arrange(padj)
  ), 
  path = "Combined_GSEA_Results.xlsx"
)
cat("\nExported full GSEA results to 'Combined_GSEA_Results.xlsx'\n")

# Load required packages
library(dplyr)
library(readr)
library(fgsea)
library(msigdbr)
library(ggplot2)
library(stringr)
library(writexl)
library(patchwork)

# ==============================================================================
# 1. LOAD MSIGDB PATHWAYS (REACTOME & GO:BP)
# ==============================================================================
cat("Loading MSigDB collections (Reactome and GO:BP)...\n")

reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol)

gobp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)

toPathways <- function(df) split(df$gene_symbol, df$gs_name)
all_pathways <- c(toPathways(reactome), toPathways(gobp))


# ==============================================================================
# 2. AUTOMATION FUNCTION FOR GSEA & PLOTTING
# ==============================================================================
run_gsea_and_plot <- function(limma_file, plot_title) {
  cat("\nProcessing:", plot_title, "...\n")
  
  # 1. Load Data & Build Ranks (using moderated t-statistic)
  res <- read_csv(limma_file, show_col_types = FALSE)
  ranks_tbl <- res %>%
    transmute(gene = EntrezGeneSymbol, rank = t) %>%
    filter(!is.na(gene), gene != "", !is.na(rank)) %>%
    group_by(gene) %>%
    slice_max(order_by = abs(rank), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  ranks <- sort(setNames(ranks_tbl$rank, ranks_tbl$gene), decreasing = TRUE)
  
  # 2. Run fgsea
  fgsea_res <- fgseaMultilevel(
    pathways = all_pathways, stats = ranks, minSize = 10, maxSize = 500
  ) %>% as_tibble()
  
  # 3. Format Labels and Filter for Plotting
  df_plot <- fgsea_res %>%
    filter(!is.na(NES), !is.na(padj), !is.na(size)) %>%
    mutate(
      label = pathway %>%
        str_replace("^(REACTOME_|GOBP_)", "") %>%
        str_replace_all("_", " ") %>% str_squish() %>% str_to_title(),
      mlog10padj = -log10(pmax(padj, .Machine$double.xmin))
    )
  
  # Get Top 5 Up and Top 5 Down (FDR < 0.05)
  top_up <- df_plot %>% filter(NES > 0, padj < 0.05) %>% arrange(desc(NES)) %>% slice_head(n = 5)
  top_down <- df_plot %>% filter(NES < 0, padj < 0.05) %>% arrange(NES) %>% slice_head(n = 5)
  
  plot_df <- bind_rows(top_up, top_down) %>%
    mutate(label_wrapped = stringr::str_wrap(toupper(label), width = 28)) %>%
    arrange(NES) %>%
    mutate(label_fac = factor(label_wrapped, levels = unique(label_wrapped)))
  
  # 4. Generate Dot Plot
  p <- ggplot(plot_df, aes(x = NES, y = label_fac)) +
    geom_point(aes(size = size, color = mlog10padj), alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_viridis_c(name = "-log10 FDR") +
    scale_size_continuous(name = "Gene set size") +
    labs(title = plot_title, x = "Normalized Enrichment Score", y = NULL) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      axis.text.y = element_text(size = 9, color = "black"), 
      axis.text.x = element_text(size = 10, color = "black"),
      legend.position = "right", aspect.ratio = 3.0
    )
  
  return(list(results = fgsea_res, plot = p))
}

# ==============================================================================
# 3. EXECUTE PIPELINE
# ==============================================================================

# Run the function for all three comparisons
gsea_AR <- run_gsea_and_plot("limma_acute_vs_recovery.csv", "Acute GBS vs Recovery GBS")
gsea_AC <- run_gsea_and_plot("limma_acute_vs_healthy.csv", "Acute GBS vs Healthy Controls")
gsea_RC <- run_gsea_and_plot("limma_recovery_vs_healthy.csv", "Recovery GBS vs Healthy Controls")


# ==============================================================================
# 4. EXPORT RESULTS
# ==============================================================================

# 1. Save Excel Tables
write_xlsx(
  list(
    "Acute_vs_Recovery" = gsea_AR$results %>% arrange(padj), 
    "Acute_vs_HC"       = gsea_AC$results %>% arrange(padj), 
    "Recovery_vs_HC"    = gsea_RC$results %>% arrange(padj)
  ), 
  path = "Combined_GSEA_Results.xlsx"
)
cat("\nExported full GSEA results to 'Combined_GSEA_Results.xlsx'\n")