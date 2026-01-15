# scRNA-seq Data processing for CRL rats (HS vs LS)

This repository contains a demonstration workflow for analyzing single-cell RNA sequencing (scRNA-seq) data from rat kidney samples. The project focuses on comparing High Salt (HS) and Low Salt (LS) treatment conditions.

> **Note:** This is a **demo version** using a downsampled dataset. The code is optimized for quick reproduction and logic verification, not for biological discovery.

## Data Description

The analysis uses a subset of kidney scRNA-seq data stored in `kidney_demo_input.rds`. 

* **Species:** Rat (*Rattus norvegicus*)
* **Sample Groups:**
    * **HS:** High Salt diet group (Samples containing "HS" in ID)
    * **LS:** Low Salt diet group (Samples containing "LS" in ID)
* **Input Format:** Pre-merged Seurat object.

## Analysis Pipeline

The R script `analysis_demo.R` performs the following steps:

1.  **QC & Filtering:**
    * Removes low-quality cells based on mitochondrial content (`percent.mt`) and feature counts (`nFeature_RNA`).
    * *Adjustment:* Thresholds are slightly relaxed for this small demo dataset.

2.  **Doublet Removal:**
    * Uses **DoubletFinder** to predict and remove doublets on a per-sample basis.
    * Includes safety checks for samples with low cell counts.

3.  **Integration (Batch Correction):**
    * Uses **Harmony** to integrate samples and remove batch effects driven by `orig.ident`.

4.  **Clustering & UMAP:**
    * Dimensionality reduction using PCA and UMAP.
    * Graph-based clustering (Louvain algorithm).

5.  **Annotation:**
    * Identifies marker genes using `FindAllMarkers`.
    * Assigns generic cell type labels (Type_0, Type_1, etc.) for demonstration purposes.

6.  **Cell-Cell Interaction Preparation:**
    * Extracts High Salt (HS) condition cells for focused analysis.
    * Converts rat genes to human homologs for CellPhoneDB compatibility.
    * Prepares expression matrix and metadata in CellPhoneDB format.

## Dependencies

To run this code, you need **R (>= 4.0)** and the following packages:

* [Seurat](https://satijalab.org/seurat/) (v5 recommended)
* [Harmony](https://github.com/immunogenomics/harmony)
* [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)
* [homologene](https://github.com/oganm/homologene) - For rat-to-human gene conversion
* [AnnoProbe](https://github.com/jmzeng1314/AnnoProbe) - For gene annotation
* `dplyr`, `ggplot2`, `patchwork`

## How to Run

1.  Clone this repository.
2.  Ensure `kidney_demo_input.rds` is in the root directory.
3.  Run the script in R or RStudio:

```R
source("analysis_demo.R")
