# Kidney and PBMC scRNA-seq Analysis Demo (HS vs LS)

This repository contains a demonstration workflow for analyzing single-cell RNA sequencing (scRNA-seq) data from rat kidney and PBMC samples. The project focuses on comparing High Salt (HS) and Low Salt (LS) treatment conditions.

> **Note:** This is a **demo version** using a downsampled kidney dataset and a simulated PBMC dataset. The code is optimized for quick reproduction and logic verification, not for biological discovery.

## Data Description

The analysis uses two Seurat objects stored in the `data` directory.

* **Organism:** Rat (*Rattus norvegicus*)
* **Sample Groups:**
    * **HS:** High Salt diet group (samples containing "HS" in the sample ID).
    * **LS:** Low Salt diet group (samples containing "LS" in the sample ID).
* **Input Format:** Pre-merged Seurat objects with an `orig.ident` metadata column.
* **Input Files:**
    * `data/kidney_demo_input.rds`: downsampled kidney demonstration dataset.
    * `data/pbmc_demo_input.rds`: downsampled PBMC dataset.

## Analysis Pipeline

The R script `src/analysis_demo.R` runs the complete workflow using the following modules:

1.  **Quality Control and Filtering:**
    - Removes low-quality cells based on mitochondrial content (`percent.mt`), feature counts (`nFeature_RNA`), and total RNA counts (`nCount_RNA`).
    - Uses relaxed thresholds for the small demonstration datasets.
    - Uses **DoubletFinder** on a per-sample basis when the package is available; records cells as `Not_tested` and continues the workflow when unavailable.
    - **Source file:** `src/preprocess.R`

2.  **Clustering:**
    - Performs normalization, variable-feature selection, scaling, PCA, graph-based clustering, and UMAP.
    - Uses **Harmony** for sample correction when available and PCA as the documented fallback.
    - Identifies cluster marker genes using `FindAllMarkers`; annotates kidney tubular, endothelial, stromal, and immune subtypes using a transparent marker reference; adds an `Uncertain_` prefix when marker evidence is insufficient.
    - **Source file:** `src/preprocess.R` (clustering), `src/annotation.R` (cell-type annotation)

3.  **Differential Expression Analysis:**
    - Aggregates counts by biological sample and cell type.
    - Performs sample-aware pseudobulk differential expression with **edgeR** to avoid treating cells as independent biological replicates.
    - Scores osmotic stress, sodium transport, RAAS, inflammation, oxidative stress, fibrosis, endothelial dysfunction, and tissue-injury programs.
    - Ranks pathway changes using effect size, replicate consistency, gene coverage, and differential-expression evidence.
    - **Source file:** `src/mechanism_analysis.R`

4.  **Cell-Cell Communication Analysis:**
    - Compares HS and LS ligand-receptor networks at the biological-sample level.
    - Prioritizes core driver pathways using effect size, replicate consistency, supporting interactions, and statistical evidence.
    - **Source file:** `src/communication_analysis.R`

5.  **Integration of PBMC and Kidney scRNA-seq Datasets:**
    - Generates a joint shared-gene embedding of kidney and PBMC cells.
    - Compares high-salt pathway responses across tissues and shared immune cell classes.
    - **Source file:** `src/pbmc_integration.R`

The workflow is organized into the following source files:

* `src/analysis_demo.R`: complete workflow entry point.
* `src/utils.R`: shared data handling and reproducibility utilities.
* `src/preprocess.R`: QC, doublet detection, clustering, and UMAP.
* `src/annotation.R`: kidney and PBMC cell-type annotation.
* `src/mechanism_analysis.R`: pseudobulk differential expression and high-salt pathway analysis.
* `src/pbmc_integration.R`: PBMC demo generation and PBMC-kidney integration.
* `src/communication_analysis.R`: differential ligand-receptor analysis and core pathway prioritization.

## Dependencies

To run the complete workflow, you need **R (>= 4.0)** and the following packages:

* [Seurat](https://satijalab.org/seurat/) (v5 recommended)
* [Harmony](https://github.com/immunogenomics/harmony)
* [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)
* `SeuratObject`
* `Matrix`
* `edgeR`
* `ggplot2`
* `patchwork`

The pathway definitions, cell-type marker references, and ligand-receptor reference are included in the source code, so no online database download is required.

## How to Run

1.  Clone this repository.
2.  Ensure `kidney_demo_input.rds` and `pbmc_demo_input.rds` are located in the `data` directory.
3.  Run the complete workflow from the repository root in R or RStudio:

```R
source("src/analysis_demo.R")
```
