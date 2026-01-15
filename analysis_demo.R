# Analysis of Kidney scRNA-seq Data (HS vs LS) - DEMO VERSION
# Project: Kidney Salt Sensitivity Study
# Description: Standard scRNA-seq pipeline including QC, Doublet Removal, Harmony Integration, and Clustering.
# Dependencies: Seurat v5, DoubletFinder, Harmony, dplyr, ggplot2

library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(patchwork)
library(dplyr)
library(harmony)

# Set working directory to the current script location
# Note: When running in RStudio, you can use Session -> Set Working Directory -> To Source File Location
# setwd("./") 

# ------------------------------------------------------------------------------
# 1. Data Loading 
# ------------------------------------------------------------------------------

# Load pre-processed demo data object instead of raw 10X files
input_file <- "kidney_demo_input.rds"

if(file.exists(input_file)){
  pbmc.merged <- readRDS(input_file)
  message("Success: Demo data loaded.")
} else {
  stop("Error: 'kidney_demo_input.rds' not found. Please ensure the data file is in the working directory.")
}

# Calculate mitochondrial percentage 
# Note: Pattern is '^Mt-' for Mouse data (use '^MT-' for Human)
pbmc.merged[["percent.mt"]] <- PercentageFeatureSet(pbmc.merged, pattern = "^Mt-")

# Subset based on QC metrics
# Demo Adjustment: Lower thresholds are relaxed here to prevent excessive filtering on this small demo dataset.
# Standard cutoff suggestions: nFeature > 200, percent.mt < 10-20% depending on tissue quality.
pbmc <- subset(pbmc.merged, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & 
                 nCount_RNA < 25000 & percent.mt < 30)

# ------------------------------------------------------------------------------
# 2. Doublet Detection (Per-Sample Basis)
# ------------------------------------------------------------------------------

samples <- unique(pbmc$orig.ident)
pbmc_list <- list()

for(sample_name in samples) {
  message("Processing DoubletFinder for: ", sample_name)
  
  # Subset current sample
  seu_sample <- subset(pbmc, subset = orig.ident == sample_name)
  
  # Safety Check: If cell count is too low (<50) after QC, skip DoubletFinder to avoid errors
  if(ncol(seu_sample) < 50) {
    warning(paste("Sample", sample_name, "has too few cells. Skipping DoubletFinder and marking as Singlet."))
    seu_sample$DoubletStatus <- "Singlet" 
    pbmc_list[[sample_name]] <- seu_sample
    next
  }
  
  # Standard Pre-processing for DoubletFinder
  seu_sample <- NormalizeData(seu_sample)
  seu_sample <- FindVariableFeatures(seu_sample, nfeatures = 2000)
  seu_sample <- ScaleData(seu_sample)
  seu_sample <- RunPCA(seu_sample, verbose = FALSE)
  
  # Demo Adjustment: Using dims 1:10 for stability on small data
  seu_sample <- RunUMAP(seu_sample, dims = 1:10, verbose = FALSE) 
  
  # Parameter Sweep
  sweep.res.list <- paramSweep(seu_sample, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # Select optimal pK
  pk_opt <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  # Define expected doublet rate 
  # Demo Adjustment: Dynamically calculated based on cell count (approx. 5%) to prevent errors
  nExp_poi <- round(0.05 * ncol(seu_sample)) 
  if(nExp_poi == 0) nExp_poi <- 1 
  
  seu_sample <- doubletFinder(seu_sample, PCs = 1:10, pN = 0.25, pK = pk_opt, 
                              nExp = nExp_poi, sct = FALSE)
  
  # Standardize metadata column names for merging
  df_col <- grep("DF.classifications", colnames(seu_sample@meta.data), value = TRUE)
  seu_sample$DoubletStatus <- seu_sample@meta.data[[df_col]]
  pbmc_list[[sample_name]] <- seu_sample
}

# Re-merge after doublet classification
pbmc <- merge(pbmc_list[[1]], y = pbmc_list[2:length(pbmc_list)])

# ------------------------------------------------------------------------------
# 3. Post-Doublet Processing & Batch Correction (Harmony)
# ------------------------------------------------------------------------------

# Filter doublets
cat("Total cells before filtering:", ncol(pbmc), "\n")
pbmc_clean <- subset(pbmc, subset = DoubletStatus == "Singlet")
cat("Total cells after filtering:", ncol(pbmc_clean), "\n")

# Normalization and Feature Selection
pbmc_clean <- NormalizeData(pbmc_clean)
pbmc_clean <- FindVariableFeatures(pbmc_clean, nfeatures = 2000)
pbmc_clean <- ScaleData(pbmc_clean)
pbmc_clean <- RunPCA(pbmc_clean, verbose = FALSE)

# Integration using Harmony 
# Integrating based on 'orig.ident' to remove batch effects between samples
pbmc_clean <- RunHarmony(pbmc_clean, group.by.vars = "orig.ident", assay.use = "RNA", plot_convergence = FALSE)

# ------------------------------------------------------------------------------
# 4. Clustering and Dimensionality Reduction
# ------------------------------------------------------------------------------

# Use Harmony embeddings for downstream steps
# Demo Adjustment: dims=1:10 and resolution=0.3 are optimized for this small subset
pbmc_clean <- FindNeighbors(pbmc_clean, reduction = "harmony", dims = 1:10)
pbmc_clean <- FindClusters(pbmc_clean, resolution = 0.3) 
pbmc_clean <- RunUMAP(pbmc_clean, reduction = "harmony", dims = 1:10)

# Save intermediate object
saveRDS(pbmc_clean, file = "kidney_processed_demo.rds")

# ------------------------------------------------------------------------------
# 5. Marker Identification & Annotation
# ------------------------------------------------------------------------------

# Ensure layers are joined (Seurat v5 requirement)
pbmc_clean <- JoinLayers(pbmc_clean)

# Find all cluster markers
# Demo Adjustment: 'min.pct' lowered to 0.25 to ensure marker detection in sparse demo data
cluster_markers <- FindAllMarkers(pbmc_clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
print(head(cluster_markers))

# Cluster Renaming
# Note: Since this is a random subset, cluster IDs may not match the full dataset.
# We assign generic names (Type_0, Type_1...) for demonstration. 
current_ids <- levels(pbmc_clean)
new_cluster_ids <- setNames(paste0("Type_", current_ids), current_ids)

pbmc_clean <- RenameIdents(pbmc_clean, new_cluster_ids)
pbmc_clean$cell_type <- Idents(pbmc_clean)

# Add Treatment Condition metadata (HS = High Salt, LS = Low Salt)
pbmc_clean$Treat <- ifelse(grepl("HS", pbmc_clean$orig.ident), "HS", "LS")

# Final Save
saveRDS(pbmc_clean, file = "kidney_final_demo.rds")

# ------------------------------------------------------------------------------
# 6. Visualization
# ------------------------------------------------------------------------------

# UMAP by Cell Type
p1 <- DimPlot(pbmc_clean, reduction = "umap", label = TRUE, repel = TRUE) + 
  NoLegend() + ggtitle("Cell Types (Demo)")

# UMAP by Treatment Group
p2 <- DimPlot(pbmc_clean, reduction = "umap", group.by = "Treat") + 
  ggtitle("Treatment Condition")

final_plot <- p1 | p2
print(final_plot)

# Save plot
ggsave("Final_UMAP_Results_Demo.png", plot = final_plot, width = 12, height = 5)