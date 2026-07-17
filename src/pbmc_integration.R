generate_pbmc_demo <- function(kidney_input, output_file, cells_per_sample = 120L, seed = 2026L) {
  assert_packages(c("Seurat", "SeuratObject", "Matrix"))
  set.seed(seed)
  kidney_input <- join_rna_layers(kidney_input)
  kidney_counts <- get_assay_matrix(kidney_input, "counts")
  genes <- rownames(kidney_counts)
  base <- Matrix::rowMeans(kidney_counts) + 1e-5
  base <- base / sum(base)
  references <- pbmc_marker_reference()
  cell_types <- c("CD4_T", "CD8_T", "NK", "B", "Classical_monocyte", "Nonclassical_monocyte", "Dendritic")
  proportions <- c(0.25, 0.18, 0.12, 0.18, 0.13, 0.09, 0.05)
  samples <- c("PBLS1", "PBLS2", "PBHS1", "PBHS2")
  hs_genes <- intersect(c("Nfat5", "Sgk1", "Nfkbia", "Ccl2", "Il1b", "Tgfb1", "Spp1", "Hmox1"), genes)
  i_all <- integer(); j_all <- integer(); x_all <- numeric()
  metadata <- list(); barcodes <- character(); column <- 1L
  for (sample in samples) {
    types <- sample(cell_types, cells_per_sample, replace = TRUE, prob = proportions)
    for (idx in seq_len(cells_per_sample)) {
      type <- types[idx]
      prob <- base * 0.55
      markers <- intersect(references[[type]], genes)
      if (length(markers)) prob[match(markers, genes)] <- prob[match(markers, genes)] + 0.40 / length(markers)
      if (grepl("HS", sample) && length(hs_genes)) prob[match(hs_genes, genes)] <- prob[match(hs_genes, genes)] + 0.05 / length(hs_genes)
      prob <- prob / sum(prob)
      library_size <- max(350L, round(stats::rlnorm(1, log(1200), 0.25)))
      draws <- sample.int(length(genes), library_size, replace = TRUE, prob = prob)
      tab <- tabulate(draws, nbins = length(genes)); nz <- which(tab > 0)
      i_all <- c(i_all, nz); j_all <- c(j_all, rep.int(column, length(nz))); x_all <- c(x_all, tab[nz])
      barcode <- paste(sample, sprintf("C%04d", idx), sep = "_")
      barcodes <- c(barcodes, barcode)
      metadata[[column]] <- data.frame(orig.ident = sample, Treat = ifelse(grepl("HS", sample), "HS", "LS"),
        tissue = "PBMC", simulated_cell_type = type, row.names = barcode)
      column <- column + 1L
    }
  }
  counts <- Matrix::sparseMatrix(i = i_all, j = j_all, x = x_all,
    dims = c(length(genes), length(barcodes)), dimnames = list(genes, barcodes))
  meta <- do.call(rbind, metadata)
  object <- Seurat::CreateSeuratObject(counts, project = "PBMC_demo", meta.data = meta, min.cells = 0, min.features = 0)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, output_file)
  object
}

run_joint_kidney_pbmc_analysis <- function(kidney, pbmc, output_dir, seed = 2026L) {
  out <- safe_dir_create(output_dir)
  kidney <- add_standard_metadata(join_rna_layers(kidney), "Kidney")
  pbmc <- add_standard_metadata(join_rna_layers(pbmc), "PBMC")
  joint <- merge(kidney, y = pbmc, add.cell.ids = c("KID", "PBMC"), project = "Kidney_PBMC_HS")
  joint <- join_rna_layers(joint)
  joint <- Seurat::NormalizeData(joint, verbose = FALSE)
  joint <- Seurat::FindVariableFeatures(joint, nfeatures = min(2500L, nrow(joint)), verbose = FALSE)
  joint <- Seurat::ScaleData(joint, features = Seurat::VariableFeatures(joint), verbose = FALSE)
  npcs <- min(30L, ncol(joint) - 1L, length(Seurat::VariableFeatures(joint)) - 1L)
  joint <- Seurat::RunPCA(joint, npcs = npcs, verbose = FALSE)
  dims <- seq_len(min(20L, npcs))
  joint <- Seurat::FindNeighbors(joint, dims = dims, verbose = FALSE)
  joint <- Seurat::FindClusters(joint, resolution = 0.4, verbose = FALSE)
  joint <- Seurat::RunUMAP(joint, dims = dims, seed.use = seed, verbose = FALSE)
  p1 <- Seurat::DimPlot(joint, group.by = "tissue") + ggplot2::ggtitle("Joint shared-gene embedding")
  p2 <- Seurat::DimPlot(joint, group.by = "cell_class", label = TRUE, repel = TRUE) + ggplot2::ggtitle("Shared cell classes")
  p3 <- Seurat::DimPlot(joint, group.by = "Treat") + ggplot2::ggtitle("Salt condition")
  save_plot((p1 + p2) / p3, file.path(out, "kidney_pbmc_joint_umap.png"), 13, 11)

  gene_sets <- high_salt_gene_sets()
  tissue_results <- lapply(c("Kidney", "PBMC"), function(tissue_name) {
    x <- subset(joint, subset = tissue == tissue_name)
    activity <- calculate_pathway_activity(x, gene_sets, cell_type_col = "cell_class")
    comp <- compare_pathway_activity(activity)
    # Always include a tissue-wide response. Shared immune classes add finer
    # concordance when both tissues contain enough cells in every replicate.
    x$integration_group <- "All_cells"
    global_activity <- calculate_pathway_activity(x, gene_sets, cell_type_col = "integration_group")
    comp <- rbind(comp, compare_pathway_activity(global_activity))
    comp$tissue <- tissue_name
    comp
  })
  response <- do.call(rbind, tissue_results)
  kidney_response <- response[response$tissue == "Kidney", ]
  pbmc_response <- response[response$tissue == "PBMC", ]
  concordance <- merge(kidney_response, pbmc_response, by = c("cell_type", "pathway"),
    suffixes = c("_Kidney", "_PBMC"))
  if (nrow(concordance)) {
    concordance$direction_concordant <- sign(concordance$delta_HS_minus_LS_Kidney) == sign(concordance$delta_HS_minus_LS_PBMC)
    concordance$joint_effect <- sqrt(abs(concordance$delta_HS_minus_LS_Kidney * concordance$delta_HS_minus_LS_PBMC))
    concordance$joint_priority <- concordance$joint_effect *
      (concordance$consistency_Kidney + concordance$consistency_PBMC) / 2 *
      ifelse(concordance$direction_concordant, 1, 0.25)
    concordance <- concordance[order(-concordance$joint_priority), ]
  }
  write_csv_safe(response, file.path(out, "tissue_specific_pathway_responses.csv"))
  write_csv_safe(concordance, file.path(out, "kidney_pbmc_response_concordance.csv"))
  saveRDS(joint, file.path(out, "kidney_pbmc_joint.rds"))
  list(object = joint, tissue_response = response, concordance = concordance)
}
