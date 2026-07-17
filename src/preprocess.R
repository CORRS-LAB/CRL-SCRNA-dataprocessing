run_optional_doublet_finder <- function(object, enabled = TRUE, seed = 2026L) {
  if (!enabled || !requireNamespace("DoubletFinder", quietly = TRUE)) {
    object$DoubletStatus <- "Not_tested"
    if (enabled) warning("DoubletFinder is not installed; doublet removal was skipped and recorded as Not_tested.")
    return(object)
  }
  set.seed(seed)
  sample_objects <- Seurat::SplitObject(object, split.by = "orig.ident")
  processed <- lapply(names(sample_objects), function(sample_name) {
    x <- sample_objects[[sample_name]]
    if (ncol(x) < 50L) {
      x$DoubletStatus <- "Not_tested_low_cell_count"
      return(x)
    }
    tryCatch({
      x <- Seurat::NormalizeData(x, verbose = FALSE)
      x <- Seurat::FindVariableFeatures(x, nfeatures = min(2000L, nrow(x)), verbose = FALSE)
      x <- Seurat::ScaleData(x, features = Seurat::VariableFeatures(x), verbose = FALSE)
      npcs <- max(5L, min(20L, ncol(x) - 1L, length(Seurat::VariableFeatures(x)) - 1L))
      x <- Seurat::RunPCA(x, npcs = npcs, verbose = FALSE)
      dims <- seq_len(min(10L, npcs))
      sweep <- DoubletFinder::paramSweep(x, PCs = dims, sct = FALSE)
      stats <- DoubletFinder::summarizeSweep(sweep, GT = FALSE)
      bcmvn <- DoubletFinder::find.pK(stats)
      pk <- suppressWarnings(as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])))
      if (!is.finite(pk)) stop("No finite pK was found")
      expected <- max(1L, round(0.05 * ncol(x)))
      x <- DoubletFinder::doubletFinder(x, PCs = dims, pN = 0.25, pK = pk,
                                        nExp = expected, sct = FALSE)
      column <- tail(grep("^DF.classifications", colnames(x[[]]), value = TRUE), 1)
      x$DoubletStatus <- x[[column]][, 1]
      x
    }, error = function(e) {
      warning("DoubletFinder failed for ", sample_name, ": ", conditionMessage(e), ". Keeping cells as Not_tested.")
      x$DoubletStatus <- "Not_tested_error"
      x
    })
  })
  merged <- if (length(processed) == 1L) processed[[1]] else merge(processed[[1]], y = processed[-1])
  join_rna_layers(merged)
}

preprocess_sc_object <- function(object, tissue = c("Kidney", "PBMC"), output_dir,
                                 use_doublet_finder = TRUE, seed = 2026L,
                                 min_features = 200L, max_features = 6000L,
                                 max_counts = 25000L, max_percent_mt = 30,
                                 resolution = 0.4, find_markers = TRUE) {
  assert_packages(c("Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork"))
  tissue <- match.arg(tissue)
  set.seed(seed)
  out <- safe_dir_create(output_dir)
  object <- join_rna_layers(object)
  object <- add_standard_metadata(object, tissue)
  mt_pattern <- if (any(grepl("^MT-", rownames(object)))) "^MT-" else "^Mt-"
  object[["percent.mt"]] <- Seurat::PercentageFeatureSet(object, pattern = mt_pattern)
  before <- ncol(object)
  keep <- object$nFeature_RNA > min_features & object$nFeature_RNA < max_features &
    object$nCount_RNA < max_counts & object$percent.mt < max_percent_mt
  object <- subset(object, cells = colnames(object)[keep])
  if (ncol(object) < 50L) stop("QC retained only ", ncol(object), " cells; review thresholds.")
  write_csv_safe(data.frame(tissue = tissue, cells_before_QC = before, cells_after_QC = ncol(object),
                            min_features = min_features, max_features = max_features,
                            max_counts = max_counts, max_percent_mt = max_percent_mt),
                 file.path(out, "qc_summary.csv"))
  object <- run_optional_doublet_finder(object, use_doublet_finder, seed)
  doublet <- grepl("Doublet", as.character(object$DoubletStatus), ignore.case = TRUE)
  if (any(doublet)) object <- subset(object, cells = colnames(object)[!doublet])
  object <- join_rna_layers(object)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- Seurat::FindVariableFeatures(object, nfeatures = min(2000L, nrow(object)), verbose = FALSE)
  object <- Seurat::ScaleData(object, features = Seurat::VariableFeatures(object), verbose = FALSE)
  npcs <- max(5L, min(30L, ncol(object) - 1L, length(Seurat::VariableFeatures(object)) - 1L))
  object <- Seurat::RunPCA(object, npcs = npcs, verbose = FALSE)
  reduction <- "pca"
  if (requireNamespace("harmony", quietly = TRUE) && length(unique(object$orig.ident)) > 1L) {
    object <- harmony::RunHarmony(object, group.by.vars = "orig.ident", plot_convergence = FALSE, verbose = FALSE)
    reduction <- "harmony"
  } else {
    message("Harmony unavailable: continuing with PCA. Sample-aware downstream tests remain valid.")
  }
  dims <- seq_len(min(15L, ncol(Seurat::Embeddings(object, reduction))))
  object <- Seurat::FindNeighbors(object, reduction = reduction, dims = dims, verbose = FALSE)
  object <- Seurat::FindClusters(object, resolution = resolution, verbose = FALSE)
  object <- Seurat::RunUMAP(object, reduction = reduction, dims = dims, seed.use = seed, verbose = FALSE)
  object <- annotate_clusters(object, tissue, out)
  object <- join_rna_layers(object)
  if (find_markers) {
    markers <- tryCatch(Seurat::FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25,
      logfc.threshold = 0.25, max.cells.per.ident = 250, verbose = FALSE),
      error = function(e) { warning("Marker testing failed: ", conditionMessage(e)); data.frame() })
    write_csv_safe(markers, file.path(out, paste0(tolower(tissue), "_cluster_markers.csv")))
  }
  p1 <- Seurat::DimPlot(object, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
    ggplot2::ggtitle(paste(tissue, "fine annotation"))
  p2 <- Seurat::DimPlot(object, reduction = "umap", group.by = "Treat") +
    ggplot2::ggtitle("Salt condition")
  save_plot(p1 + p2, file.path(out, paste0(tolower(tissue), "_umap.png")), 14, 6)
  saveRDS(object, file.path(out, paste0(tolower(tissue), "_processed.rds")))
  object
}
