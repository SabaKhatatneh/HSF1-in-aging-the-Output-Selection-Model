# HSF1 Target Gene GO Enrichment Analysis

Reproducible R code for GO Biological Process enrichment of HSF1 target genes, as presented in:

> HSF1 in aging: the Output Selection Model  
> Saba Khatatneh, Csaba Sőti, Milán Somogyvári, 2026
> *Mechanisms of Ageing and Development*

---

## Description

This repository contains the complete R script to reproduce Figures 5 and 6 from the manuscript. The analysis performs GO Biological Process enrichment on HSF1 target genes from HSF1Base, annotated with:

- **HSF1 output class** (Output 1/2/3) matching the Output Allocation Model
- **Aging hallmark category** (Lopez-Otin et al., 2023)

Two publication-quality figures are produced:
- **Figure 5** — *Homo sapiens* upregulated and downregulated HSF1 targets
- **Figure 6** — *Caenorhabditis elegans* upregulated and downregulated HSF1 targets

---

## Data

Download `all_species.xlsx` from HSF1Base:

> Kovacs et al. (2019). HSF1Base: A Comprehensive Database of HSF1 (Heat Shock Factor 1) Target Genes. *International Journal of Molecular Sciences*, 20(22), 5815.  
> https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6888049/

Place the file in the same directory as the R script.

---

## Requirements

- R >= 4.6.0
- Bioconductor >= 3.23

All required packages are installed automatically when the script is run for the first time.

| Package | Version | Reference |
|---|---|---|
| clusterProfiler | 4.x | Yu et al., 2012, OMICS 16:284 |
| org.Hs.eg.db | Bioconductor 3.23 | Carlson, 2024 |
| org.Ce.eg.db | Bioconductor 3.23 | Carlson, 2024 |
| ggplot2 | 3.x | Wickham, 2016 |
| ggnewscale | 0.x | Campitelli, 2024 |
| cowplot | 1.x | Wilke, 2024 |
| patchwork | 1.x | Pedersen, 2024 |
| readxl | 1.x | Wickham & Bryan, 2023 |

---

## Usage

1. Clone or download this repository
2. Download `all_species.xlsx` from HSF1Base (see above)
3. Open `HSF1_GO_enrichment_GITHUB.R` in RStudio
4. Set your working directory in **STEP 1** of the script
5. Click **Source** to run the complete analysis

The script will:
- Install all required packages automatically (first run only)
- Load and process HSF1Base data
- Run GO enrichment for all four gene sets
- Export figures and supplementary tables

---

## Output files

| File | Description |
|---|---|
| `Figure5_Homo_sapiens.tiff` | Figure 5 — journal submission (600 DPI) |
| `Figure5_Homo_sapiens.pdf` | Figure 5 — vector format |
| `Figure6_C_elegans.tiff` | Figure 6 — journal submission (600 DPI) |
| `Figure6_C_elegans.pdf` | Figure 6 — vector format |
| `Suppl_Table_S1_GO_human_upregulated.csv` | Full enrichment results |
| `Suppl_Table_S2_GO_human_downregulated.csv` | Full enrichment results |
| `Suppl_Table_S3_GO_worm_upregulated.csv` | Full enrichment results |
| `Suppl_Table_S4_GO_worm_downregulated.csv` | Full enrichment results |

---

## Methods summary

HSF1 target genes were retrieved from HSF1Base (Kovacs et al., 2019) in PSI-MITAB format and stratified by NCBI taxonomy identifier (Human: taxid:9606; *C. elegans*: taxid:6239) and causal statement (upregulated: MI:2236; downregulated: MI:2241). GO Biological Process enrichment was performed using clusterProfiler v4 (Yu et al., 2012) with Benjamini-Hochberg multiple testing correction (adjusted p < 0.05). Redundant GO terms were removed using the simplify function (semantic similarity cutoff = 0.7). Each GO term was annotated with an HSF1 output class and aging hallmark category as described in the manuscript.

---

## Citation

If you use this code, please cite our article and the underlying tools:

> Yu G, Wang LG, Han Y, He QY (2012). clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters. *OMICS*, 16(5), 284-287.

> Kovacs D et al. (2019). HSF1Base: A Comprehensive Database of HSF1 Target Genes. *Int. J. Mol. Sci.* 20(22), 5815.

> Lopez-Otin C et al. (2023). Hallmarks of aging: An expanding universe. *Cell*, 186(2), 243-278.

---

## License

MIT License — free to use and modify with attribution.
