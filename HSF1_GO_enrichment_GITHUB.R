# =============================================================================
# HSF1 Target Gene GO Biological Process Enrichment Analysis
# =============================================================================
# Title:   GO Biological Process enrichment of HSF1 target genes
# Journal: Mechanisms of Ageing and Development (Elsevier)
# Authors: [Your name(s)]
# Date:    May 2026
#
# Description:
#   Performs GO Biological Process enrichment on HSF1 target genes from
#   HSF1Base, annotated with HSF1 output class (Output 1/2/3) and aging
#   hallmark category (Lopez-Otin et al., 2023).
#   Produces two publication-quality figures:
#     - Figure 5: Homo sapiens (upregulated + downregulated)
#     - Figure 6: C. elegans (upregulated + downregulated)
#
# Data source:
#   HSF1Base: Kovacs et al. (2019) Int. J. Mol. Sci. 20, 5815
#   File: all_species.xlsx
#
# Output files:
#   Figure5_Homo_sapiens.pdf / .tiff
#   Figure6_C_elegans.pdf / .tiff
#   Suppl_GO_human_upregulated.csv
#   Suppl_GO_human_downregulated.csv
#   Suppl_GO_worm_upregulated.csv
#   Suppl_GO_worm_downregulated.csv
#
# Software:
#   R >= 4.6.0, Bioconductor >= 3.23
#   clusterProfiler v4 (Yu et al., 2012, OMICS 16:284-287)
#   org.Hs.eg.db, org.Ce.eg.db (Bioconductor)
#   ggplot2, ggnewscale, cowplot, patchwork
#
# Reproducibility:
#   Run sessionInfo() at end of script for full package version record
# =============================================================================


# -----------------------------------------------------------------------------
# STEP 1: SET WORKING DIRECTORY
#   Change this path to the folder containing all_species.xlsx
# -----------------------------------------------------------------------------

setwd("C:/Documents/HSF1 review article/GO terms")


# -----------------------------------------------------------------------------
# STEP 2: INSTALL PACKAGES
#   Safe to re-run — skips already installed packages
# -----------------------------------------------------------------------------

if (!require("readxl"))     install.packages("readxl")
if (!require("dplyr"))      install.packages("dplyr")
if (!require("ggplot2"))    install.packages("ggplot2")
if (!require("stringr"))    install.packages("stringr")
if (!require("scales"))     install.packages("scales")
if (!require("patchwork"))  install.packages("patchwork")
if (!require("cowplot"))    install.packages("cowplot")
if (!require("ggnewscale")) install.packages("ggnewscale")

if (!require("BiocManager"))     install.packages("BiocManager")
if (!require("clusterProfiler")) BiocManager::install("clusterProfiler")
if (!require("org.Hs.eg.db"))    BiocManager::install("org.Hs.eg.db")
if (!require("org.Ce.eg.db"))    BiocManager::install("org.Ce.eg.db")


# -----------------------------------------------------------------------------
# STEP 3: LOAD PACKAGES
# -----------------------------------------------------------------------------

library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)
library(patchwork)
library(cowplot)
library(ggnewscale)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Ce.eg.db)


# -----------------------------------------------------------------------------
# STEP 4: COLOR SCHEME
#   Matches HSF1 Output Allocation figure colors exactly:
#   Output 1 = orange  (#D85A30) — Acute inducible proteotoxic defense
#   Output 2 = teal    (#1D9E75) — Basal tissue-maintenance proteostasis
#   Output 3 = purple  (#7F77DD) — Biosynthetic / growth-supporting
#   Unclassified = gray (#AAAAAA)
# -----------------------------------------------------------------------------

output_colors <- c(
  "Output 1"     = "#D85A30",
  "Output 2"     = "#1D9E75",
  "Output 3"     = "#7F77DD",
  "Unclassified" = "#AAAAAA"
)

hallmark_colors <- c(
  "Loss of proteostasis"         = "#C0392B",
  "Deregulated nutrient sensing" = "#8E44AD",
  "Epigenetic alterations"       = "#2980B9",
  "Disabled macroautophagy"      = "#27AE60",
  "Cellular senescence"          = "#E67E22",
  "Genomic instability"          = "#16A085",
  "Mitochondrial dysfunction"    = "#D35400",
  "Altered intercellular comm."  = "#7F8C8D",
  "Other"                        = "#BDC3C7"
)


# -----------------------------------------------------------------------------
# STEP 5: ANNOTATION FUNCTIONS
#   Assigns each GO term to an HSF1 output class and aging hallmark
#
#   Output class assignment based on:
#     Li et al. (2016) — degenerate HSE / EFL-1/DPL-1 / Output 2
#     Morphis et al. (2022) — adult somatic sub-chaperome / Output 2
#     Mendillo et al. (2012) — cancer biosynthetic program / Output 3
#     Labbadia & Morimoto (2015) — inducible HSR / Output 1
#
#   Note: "protein folding" assigned to Output 2 (constitutive chaperones)
#   not Output 1, as constitutive chaperone-mediated folding is a basal
#   tissue-maintenance function (Morphis et al., 2022)
#
#   Hallmark assignment based on:
#     Lopez-Otin et al. (2023) — hallmarks of aging
# -----------------------------------------------------------------------------

assign_output_class <- function(description) {
  desc <- tolower(description)
  dplyr::case_when(

    # Output 1 — Acute inducible proteotoxic defense
    # Canonical HSR: rapid induction upon proteotoxic stress
    str_detect(desc, paste(c(
      "heat", "unfold",
      "stress response", "proteotox",
      "response to heat", "response to temperature",
      "thermotolerance", "response to unfolded"
    ), collapse = "|")) ~ "Output 1",

    # Output 2 — Basal tissue-maintenance proteostasis
    # Constitutive chaperones, autophagy, cytoskeletal, proteostasis
    # Note: protein folding = Output 2 (constitutive, not inducible)
    str_detect(desc, paste(c(
      "protein folding", "protein fold",
      "chaperone", "chaperone-mediated",
      "autophagy", "macroautophagy",
      "cytoskeletal", "cytoskeleton",
      "proteasom", "ubiquitin",
      "mitochondr", "oxidative stress",
      "lysosom", "proteostasis",
      "protein catabol", "protein quality",
      "protein homeostasis",
      "actin", "tubulin",
      "muscle contraction", "muscle fiber",
      "collagen", "cuticle",
      "tissue maintenance",
      "autophagy of mitochondr",
      "protein polyubiquit",
      "vacuole organization"
    ), collapse = "|")) ~ "Output 2",

    # Output 3 — Biosynthetic, reproductive and growth-supporting
    # Translation, metabolism, cell cycle, immune, chromatin, signaling
    str_detect(desc, paste(c(
      "translat", "ribosom", "biosynthet",
      "metabol", "proliferat", "germline",
      "reproduct", "growth", "cell cycle",
      "nuclear division", "dna repair",
      "chromatin", "epigenet", "rna splicing",
      "mrna", "nucleosome", "dna metabol",
      "dna template", "cell adhesion",
      "cell differentiat", "immune",
      "cytokine", "inflammatory",
      "leukocyte", "mononuclear",
      "apoptot", "cell death", "senescen",
      "fat cell", "lipid", "fatty acid",
      "amino acid metabol", "organic acid",
      "small molecule", "carboxylic acid",
      "sulfur compound", "phosphate",
      "neuron", "nervous",
      "signal transduct", "signaling cassette",
      "erk", "mapk", "kinase cascade",
      "dauer", "catabol", "microrna",
      "mirna", "skin epidermis",
      "rhythmic", "circadian",
      "biosynthetic process"
    ), collapse = "|")) ~ "Output 3",

    TRUE ~ "Unclassified"
  )
}

assign_hallmark <- function(description) {
  desc <- tolower(description)
  dplyr::case_when(
    str_detect(desc, paste(c(
      "protein folding", "protein fold",
      "chaperone", "unfold",
      "proteostasis", "proteotox",
      "heat shock", "protein catabol",
      "protein quality", "unfolded protein"
    ), collapse = "|")) ~ "Loss of proteostasis",

    str_detect(desc, paste(c(
      "autophagy", "macroautophagy",
      "lysosom", "mitophagy",
      "vacuole", "protein polyubiquit"
    ), collapse = "|")) ~ "Disabled macroautophagy",

    str_detect(desc, paste(c(
      "metabol", "nutrient", "lipid",
      "fatty acid", "biosynthet",
      "ribosom", "translat",
      "amino acid", "organic acid",
      "small molecule", "carboxylic",
      "sulfur", "phosphate", "catabol"
    ), collapse = "|")) ~ "Deregulated nutrient sensing",

    str_detect(desc, paste(c(
      "chromatin", "epigenet", "histone",
      "methylat", "acetylat",
      "heterochromatin", "nucleosome",
      "dna metabol", "dna template",
      "dna-templated"
    ), collapse = "|")) ~ "Epigenetic alterations",

    str_detect(desc, paste(c(
      "apoptot", "senescen", "cell cycle",
      "nuclear division", "cell death",
      "cell proliferat", "cell differentiat",
      "cell adhesion", "leukocyte",
      "mononuclear", "immune", "cytokine",
      "dauer", "rhythmic", "circadian"
    ), collapse = "|")) ~ "Cellular senescence",

    str_detect(desc, paste(c(
      "dna repair", "genome", "genomic",
      "nucleotide", "dna damage"
    ), collapse = "|")) ~ "Genomic instability",

    str_detect(desc, paste(c(
      "mitochondr", "oxidative",
      "reactive oxygen", "electron transport"
    ), collapse = "|")) ~ "Mitochondrial dysfunction",

    str_detect(desc, paste(c(
      "intercellular", "cytokine",
      "signaling cassette",
      "signal transduct", "erk", "mapk"
    ), collapse = "|")) ~ "Altered intercellular comm.",

    TRUE ~ "Other"
  )
}


# -----------------------------------------------------------------------------
# STEP 6: LOAD DATA
# -----------------------------------------------------------------------------

cat("Loading HSF1Base data...\n")

df <- read_excel("all_species.xlsx")

# Rename key columns
colnames(df)[2]  <- "target_id"   # Ensembl/WormBase gene ID
colnames(df)[11] <- "taxid_B"     # NCBI taxonomy ID of target species
colnames(df)[46] <- "causal"      # MI:2236 = upregulated; MI:2241 = downregulated

cat("Total interactions loaded:", nrow(df), "\n")


# -----------------------------------------------------------------------------
# STEP 7: EXTRACT GENE LISTS
#   Human:     Ensembl IDs (ENSG...)
#   C. elegans: WormBase IDs (WBGene...)
#   Direction: MI:2236 = upregulated; MI:2241 = downregulated
#   Note: entries with missing causal annotation excluded
# -----------------------------------------------------------------------------

clean_id <- function(x) str_replace(x, "^ensembl:", "") %>% str_trim()

human_up <- df %>%
  filter(taxid_B == "taxid:9606", causal == "MI:2236") %>%
  pull(target_id) %>% clean_id() %>% unique() %>% na.omit()

human_down <- df %>%
  filter(taxid_B == "taxid:9606", causal == "MI:2241") %>%
  pull(target_id) %>% clean_id() %>% unique() %>% na.omit()

worm_up <- df %>%
  filter(taxid_B == "taxid:6239", causal == "MI:2236") %>%
  pull(target_id) %>% clean_id() %>% unique() %>% na.omit()

worm_down <- df %>%
  filter(taxid_B == "taxid:6239", causal == "MI:2241") %>%
  pull(target_id) %>% clean_id() %>% unique() %>% na.omit()

cat("\nGene list sizes:\n")
cat("  Human upregulated   :", length(human_up),   "unique genes\n")
cat("  Human downregulated :", length(human_down), "unique genes\n")
cat("  Worm  upregulated   :", length(worm_up),    "unique genes\n")
cat("  Worm  downregulated :", length(worm_down),  "unique genes\n")


# -----------------------------------------------------------------------------
# STEP 8: GO ENRICHMENT FUNCTION
#   Parameters:
#     gene_ids  — character vector of Ensembl or WormBase IDs
#     orgdb     — annotation database (org.Hs.eg.db or org.Ce.eg.db)
#     keytype   — "ENSEMBL" for human; "WORMBASE" for C. elegans
#     label     — label for progress messages
#     top_n     — number of top terms to retain (default 10)
#     p_cutoff  — BH-adjusted p-value cutoff (default 0.05)
#     q_cutoff  — q-value cutoff (default 0.20)
# -----------------------------------------------------------------------------

run_go <- function(gene_ids, orgdb, keytype, label,
                   top_n = 10, p_cutoff = 0.05, q_cutoff = 0.20) {

  cat("\nRunning GO enrichment:", label, "...\n")

  if (length(gene_ids) < 5) {
    message("  Too few genes (n=", length(gene_ids), ") — skipping: ", label)
    return(NULL)
  }

  res <- enrichGO(
    gene          = gene_ids,
    OrgDb         = orgdb,
    keyType       = keytype,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = p_cutoff,
    qvalueCutoff  = q_cutoff,
    readable      = TRUE
  )

  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    message("  No significant terms found: ", label)
    return(NULL)
  }

  # Remove redundant GO terms (semantic similarity cutoff = 0.7)
  res_simplified <- simplify(res, cutoff = 0.7, by = "p.adjust")

  df_out <- as.data.frame(res_simplified) %>%
    arrange(p.adjust) %>%
    head(top_n) %>%
    mutate(
      group        = label,
      log_padj     = -log10(p.adjust),
      GeneRatio_n  = as.numeric(str_extract(GeneRatio, "^\\d+")),
      GeneRatio_d  = as.numeric(str_extract(GeneRatio, "\\d+$")),
      ratio        = GeneRatio_n / GeneRatio_d,
      Description  = str_wrap(Description, width = 22),
      output_class = assign_output_class(Description),
      hallmark     = assign_hallmark(Description)
    )

  cat("  Significant terms:", nrow(as.data.frame(res_simplified)),
      "| Showing top:", nrow(df_out), "\n")
  cat("  Output class distribution:\n")
  print(table(df_out$output_class))

  return(df_out)
}


# -----------------------------------------------------------------------------
# STEP 9: RUN ENRICHMENT FOR ALL FOUR GENE SETS
# -----------------------------------------------------------------------------

go_human_up   <- run_go(human_up,   org.Hs.eg.db, "ENSEMBL",
                         "Homo sapiens — Upregulated")
go_human_down <- run_go(human_down, org.Hs.eg.db, "ENSEMBL",
                         "Homo sapiens — Downregulated")
go_worm_up    <- run_go(worm_up,    org.Ce.eg.db, "WORMBASE",
                         "C. elegans — Upregulated")
go_worm_down  <- run_go(worm_down,  org.Ce.eg.db, "WORMBASE",
                         "C. elegans — Downregulated")


# -----------------------------------------------------------------------------
# STEP 10: PLOT FUNCTION
#   Design:
#     - No figure title (title goes in manuscript caption only)
#     - Dot fill:        -log10(BH-adjusted p-value) gradient
#     - Dot size:        gene count
#     - Outer circle:    HSF1 output class color (thicker ring)
#     - Square marker:   aging hallmark color (right of each row)
#     - base_size = 11:  readable but compact for 19x24cm page
#     - Two independent color scales via ggnewscale
# -----------------------------------------------------------------------------

plot_go_annotated <- function(go_df, title_expr, direction = "up") {

  # Empty panel for missing data
  if (is.null(go_df) || nrow(go_df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = "No significant terms",
                 size = 4, hjust = 0.5, color = "grey50") +
        theme_void()
    )
  }

  # Color gradient direction
  fill_low  <- if (direction == "up") "#FEE08B" else "#C6DBEF"
  fill_high <- if (direction == "up") "#A50026" else "#08306B"

  # Order terms by gene ratio
  go_df <- go_df %>%
    mutate(Description = factor(Description,
                                levels = Description[order(ratio)]))

  # Position for hallmark square markers (right of plot area)
  x_max   <- max(go_df$ratio, na.rm = TRUE)
  x_min   <- min(go_df$ratio, na.rm = TRUE)
  x_range <- x_max - x_min
  strip_x <- x_max + x_range * 0.10

  p <- ggplot(go_df, aes(x = ratio, y = Description)) +

    # Layer 1: outer circle = HSF1 output class (thick ring)
    geom_point(aes(size  = Count + 3,
                   color = output_class),
               alpha = 1) +
    scale_color_manual(
      values = output_colors,
      name   = "HSF1 output class (outer circle)",
      labels = c(
        "Output 1"     = "Output 1 - Acute inducible defense",
        "Output 2"     = "Output 2 - Basal tissue maintenance",
        "Output 3"     = "Output 3 - Biosynthetic / growth",
        "Unclassified" = "Unclassified"
      )
    ) +

    # Layer 2: inner dot = significance (fill gradient)
    geom_point(aes(size = Count,
                   fill = log_padj),
               shape  = 21,
               color  = "white",
               stroke = 0.3,
               alpha  = 0.92) +
    scale_fill_gradient(
      low  = fill_low,
      high = fill_high,
      name = expression(-log[10](p[adj]))
    ) +

    # New color scale for hallmark squares
    new_scale_color() +

    # Layer 3: hallmark square marker
    geom_point(aes(x     = strip_x,
                   color = hallmark),
               size  = 2.5,
               shape = 15) +
    scale_color_manual(
      values = hallmark_colors,
      name   = "Aging hallmark (square)"
    ) +

    # Size scale
    scale_size_continuous(
      name  = "Gene count",
      range = c(2, 8),
      guide = guide_legend(
        override.aes = list(
          color = "grey60", fill = "grey80", shape = 21
        )
      )
    ) +

    # X axis
    scale_x_continuous(
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0.05, 0.18))
    ) +

    coord_cartesian(clip = "off") +
    labs(title = title_expr, x = "Gene ratio", y = NULL) +

    theme_bw(base_size = 9) +
    theme(
      plot.title         = element_text(
                             size   = 9,
                             face   = "bold",
                             hjust  = 0,
                             margin = margin(b = 5)),
      axis.text.y        = element_text(size = 7.5, color = "grey15"),
      axis.text.x        = element_text(size = 7.5),
      axis.title.x       = element_text(size = 8.5, margin = margin(t = 4)),
      panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3),
      panel.grid.major.x = element_line(color = "grey93", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "grey60", linewidth = 0.4),
      legend.position    = "none",   # legends handled by shared bottom legend
      plot.margin        = margin(5, 45, 5, 5)
    )

  return(p)
}


# -----------------------------------------------------------------------------
# STEP 11: BUILD SHARED BOTTOM LEGEND
#   Single shared legend for both figures — output class + hallmark + dot anatomy
# -----------------------------------------------------------------------------

# Output class legend
output_leg_df <- data.frame(
  x     = 1:4,
  y     = rep(1, 4),
  class = names(output_colors),
  stringsAsFactors = FALSE
)

output_leg_plot <- ggplot(output_leg_df,
                           aes(x = x, y = y,
                               fill = class, color = class)) +
  geom_point(size = 5, shape = 21, stroke = 2) +
  scale_fill_manual(
    values = output_colors,
    name   = "HSF1 output class (outer circle)",
    labels = c(
      "Output 1"     = "Output 1 - Acute inducible defense",
      "Output 2"     = "Output 2 - Basal tissue maintenance",
      "Output 3"     = "Output 3 - Biosynthetic / growth-supporting",
      "Unclassified" = "Unclassified"
    )
  ) +
  scale_color_manual(values = output_colors, guide = "none") +
  theme_void() +
  theme(
    legend.position  = "bottom",
    legend.title     = element_text(size = 8, face = "bold"),
    legend.text      = element_text(size = 7.5),
    legend.key.size  = unit(0.5, "cm"),
    legend.direction = "horizontal"
  )

output_leg <- cowplot::get_legend(output_leg_plot)

# Hallmark legend
hallmark_leg_df <- data.frame(
  x        = seq_along(hallmark_colors),
  y        = rep(1, length(hallmark_colors)),
  hallmark = names(hallmark_colors),
  stringsAsFactors = FALSE
)

hallmark_leg_plot <- ggplot(hallmark_leg_df,
                             aes(x = x, y = y, color = hallmark)) +
  geom_point(size = 3.5, shape = 15) +
  scale_color_manual(
    values = hallmark_colors,
    name   = "Aging hallmark (square marker; Lopez-Otin et al., 2023)"
  ) +
  theme_void() +
  theme(
    legend.position  = "bottom",
    legend.title     = element_text(size = 8, face = "bold"),
    legend.text      = element_text(size = 7.5),
    legend.key.size  = unit(0.4, "cm"),
    legend.direction = "horizontal"
  )

hallmark_leg <- cowplot::get_legend(hallmark_leg_plot)


# -----------------------------------------------------------------------------
# STEP 12: BUILD FIGURE 5 — Homo sapiens
# -----------------------------------------------------------------------------

p1 <- plot_go_annotated(
  go_human_up,
  expression(italic("Homo sapiens") ~ "— Upregulated"),
  "up"
)

p2 <- plot_go_annotated(
  go_human_down,
  expression(italic("Homo sapiens") ~ "— Downregulated"),
  "down"
)

# Combine two panels side by side
human_grid <- plot_grid(
  p1, p2,
  ncol           = 2,
  labels         = c("A", "B"),
  label_size     = 12,
  label_fontface = "bold"
)

# Add shared legends below
figure5 <- plot_grid(
  human_grid,
  output_leg,
  hallmark_leg,
  ncol        = 1,
  rel_heights = c(1, 0.10, 0.10)
)


# -----------------------------------------------------------------------------
# STEP 13: BUILD FIGURE 6 — C. elegans
# -----------------------------------------------------------------------------

p3 <- plot_go_annotated(
  go_worm_up,
  expression(italic("C. elegans") ~ "— Upregulated"),
  "up"
)

p4 <- plot_go_annotated(
  go_worm_down,
  expression(italic("C. elegans") ~ "— Downregulated"),
  "down"
)

# Combine two panels side by side
worm_grid <- plot_grid(
  p3, p4,
  ncol           = 2,
  labels         = c("A", "B"),
  label_size     = 12,
  label_fontface = "bold"
)

# Add shared legends below
figure6 <- plot_grid(
  worm_grid,
  output_leg,
  hallmark_leg,
  ncol        = 1,
  rel_heights = c(1, 0.10, 0.10)
)


# -----------------------------------------------------------------------------
# STEP 14: EXPORT FIGURES
#   Format: TIFF at 600 DPI (Elsevier combination figure requirement)
#           PDF for vector quality submission
#   Size:   19 x 24 cm (Elsevier full page double column)
#   No title in figure — title provided in manuscript caption only
# -----------------------------------------------------------------------------

cat("\nExporting Figure 5 (Homo sapiens)...\n")

ggsave("Figure5_Homo_sapiens.tiff",
       plot   = figure5,
       width  = 26,
       height = 18,
       units  = "cm",
       dpi    = 600,
       device = "tiff",
       bg     = "white")

ggsave("Figure5_Homo_sapiens.pdf",
       plot   = figure5,
       width  = 26,
       height = 18,
       units  = "cm",
       device = "pdf",
       bg     = "white")

cat("  Saved: Figure5_Homo_sapiens.tiff / .pdf\n")

cat("\nExporting Figure 6 (C. elegans)...\n")

ggsave("Figure6_C_elegans.tiff",
       plot   = figure6,
       width  = 26,
       height = 18,
       units  = "cm",
       dpi    = 600,
       device = "tiff",
       bg     = "white")

ggsave("Figure6_C_elegans.pdf",
       plot   = figure6,
       width  = 26,
       height = 18,
       units  = "cm",
       device = "pdf",
       bg     = "white")

cat("  Saved: Figure6_C_elegans.tiff / .pdf\n")


# -----------------------------------------------------------------------------
# STEP 15: EXPORT SUPPLEMENTARY ENRICHMENT TABLES
#   Full GO enrichment results for Supplementary Data
#   Includes output class and hallmark annotations for each term
# -----------------------------------------------------------------------------

cat("\nSaving supplementary enrichment tables...\n")

save_csv <- function(go_df, filename) {
  if (!is.null(go_df) && nrow(go_df) > 0) {
    go_df %>%
      dplyr::select(
        ID, Description, GeneRatio, BgRatio,
        pvalue, p.adjust, qvalue, Count,
        output_class, hallmark, geneID
      ) %>%
      write.csv(filename, row.names = FALSE)
    cat("  Saved:", filename, "\n")
  }
}

save_csv(go_human_up,   "Suppl_Table_S1_GO_human_upregulated.csv")
save_csv(go_human_down, "Suppl_Table_S2_GO_human_downregulated.csv")
save_csv(go_worm_up,    "Suppl_Table_S3_GO_worm_upregulated.csv")
save_csv(go_worm_down,  "Suppl_Table_S4_GO_worm_downregulated.csv")


# -----------------------------------------------------------------------------
# STEP 16: SESSION INFO
#   Print full package version information for Methods section
# -----------------------------------------------------------------------------

cat("\n=== Session Info (include in Methods / GitHub README) ===\n")
sessionInfo()

cat("\n=== DONE ===\n")
cat("All files saved to:", getwd(), "\n")
cat("\nFiles produced:\n")
cat("  Figure5_Homo_sapiens.tiff  — submit to journal\n")
cat("  Figure5_Homo_sapiens.pdf   — vector version\n")
cat("  Figure6_C_elegans.tiff     — submit to journal\n")
cat("  Figure6_C_elegans.pdf      — vector version\n")
cat("  Suppl_Table_S1_GO_human_upregulated.csv\n")
cat("  Suppl_Table_S2_GO_human_downregulated.csv\n")
cat("  Suppl_Table_S3_GO_worm_upregulated.csv\n")
cat("  Suppl_Table_S4_GO_worm_downregulated.csv\n")
