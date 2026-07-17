kidney_marker_reference <- function() {
  list(
    "Podocyte" = c("Nphs1", "Nphs2", "Podxl", "Synpo", "Wt1", "Magi2"),
    "Proximal_tubule_S1" = c("Slc5a2", "Slc34a1", "Lrp2", "Cubn", "Akap12"),
    "Proximal_tubule_S2" = c("Slc22a6", "Slc22a8", "Aqp1", "Acsm2", "Hnf4a"),
    "Proximal_tubule_S3" = c("Slc7a13", "Slc13a3", "Gsta1", "Gsta2", "Havcr1"),
    "Thin_descending_limb" = c("Aqp1", "Slc14a2", "Clcnka", "Cryab"),
    "Thin_ascending_limb" = c("Clcnka", "Clcnkb", "Umod", "Kcnj1"),
    "TAL" = c("Slc12a1", "Umod", "Kcnj1", "Clcnkb", "Casr"),
    "DCT1" = c("Slc12a3", "Pvalb", "Trpm6", "Calb1"),
    "DCT2_CNT" = c("Calb1", "Trpv5", "Slc8a1", "Fxyd4", "Scnn1g"),
    "Collecting_duct_principal" = c("Aqp2", "Aqp3", "Fxyd4", "Scnn1a", "Scnn1g"),
    "Intercalated_A" = c("Atp6v1b1", "Atp6v0d2", "Slc4a1", "Foxi1"),
    "Intercalated_B" = c("Atp6v1b1", "Slc26a4", "Slc4a9", "Foxi1"),
    "Glomerular_endothelial" = c("Ehd3", "Kdr", "Pecam1", "Emcn", "Plvap"),
    "Peritubular_endothelial" = c("Plvap", "Emcn", "Kdr", "Klf2", "Pecam1"),
    "Lymphatic_endothelial" = c("Prox1", "Pdpn", "Flt4", "Ccl21"),
    "Fibroblast" = c("Col1a1", "Col1a2", "Dcn", "Pdgfra", "Lum"),
    "Myofibroblast" = c("Acta2", "Tagln", "Col3a1", "Postn", "Ctgf"),
    "Pericyte" = c("Rgs5", "Pdgfrb", "Cspg4", "Des", "Kcnj8"),
    "Vascular_smooth_muscle" = c("Acta2", "Myh11", "Tagln", "Cnn1", "Smtn"),
    "Macrophage" = c("Adgre1", "C1qa", "C1qb", "Cd68", "Lyz2"),
    "Monocyte" = c("Lyz2", "S100a8", "S100a9", "Ccr2", "Ly6c"),
    "Dendritic" = c("Itgax", "Flt3", "Clec10a", "H2-Ab1", "Cd74"),
    "T_cell" = c("Cd3d", "Cd3e", "Trbc1", "Lck", "Il7r"),
    "NK_cell" = c("Nkg7", "Klrk1", "Prf1", "Gzmb", "Ccl5"),
    "B_cell" = c("Cd79a", "Ms4a1", "Cd37", "Cd74", "Cd19"),
    "Neutrophil" = c("S100a8", "S100a9", "Csf3r", "Mpo", "Retnlg"),
    "Cycling" = c("Mki67", "Top2a", "Cenpf", "Ube2c", "Birc5")
  )
}

pbmc_marker_reference <- function() {
  list(
    "CD4_T" = c("Cd3d", "Cd3e", "Il7r", "Ltb", "Lck"),
    "CD8_T" = c("Cd3d", "Cd3e", "Cd8a", "Ccl5", "Lck"),
    "NK" = c("Nkg7", "Klrk1", "Prf1", "Gzmb", "Ccl5"),
    "B" = c("Cd79a", "Ms4a1", "Cd37", "Cd74", "Cd19"),
    "Plasma" = c("Jchain", "Mzb1", "Sdc1", "Xbp1", "Igha"),
    "Classical_monocyte" = c("Lyz2", "S100a8", "S100a9", "Ccr2", "Ly6c"),
    "Nonclassical_monocyte" = c("Lyz2", "Fcgr3", "LST1", "Ms4a7", "C1qa"),
    "Dendritic" = c("Itgax", "Flt3", "Clec10a", "H2-Ab1", "Cd74")
  )
}

broad_cell_class <- function(labels) {
  x <- as.character(labels)
  out <- rep("Other", length(x))
  out[grepl("T_cell|CD4_T|CD8_T", x)] <- "T_cell"
  out[grepl("NK", x)] <- "NK_cell"
  out[grepl("B_cell|^B$|Plasma", x)] <- "B_cell"
  out[grepl("Mono|Macrophage", x)] <- "Mononuclear_phagocyte"
  out[grepl("Dendritic", x)] <- "Dendritic"
  out[grepl("Neutrophil", x)] <- "Neutrophil"
  out[grepl("endothelial", x, ignore.case = TRUE)] <- "Endothelial"
  out[grepl("Fibroblast|Myofibroblast|Pericyte|smooth_muscle", x)] <- "Stromal"
  out[grepl("Podocyte", x)] <- "Podocyte"
  out[grepl("tubule|limb|TAL|DCT|Collecting|Intercalated", x, ignore.case = TRUE)] <- "Tubular_epithelial"
  out
}

score_cluster_markers <- function(object, marker_reference, cluster_col = "seurat_clusters",
                                  min_detected_markers = 2L, min_margin = 0.05,
                                  min_marker_detection = 0.05) {
  if (!cluster_col %in% colnames(object[[]])) stop("Missing cluster metadata: ", cluster_col)
  object <- join_rna_layers(object)
  groups <- object[[cluster_col]][, 1]
  avg <- average_sparse_expression(object, groups, layer = "data")
  detected <- average_sparse_expression(object, groups, layer = "data", detection = TRUE)
  gene_center <- rowMeans(avg)
  gene_sd <- apply(avg, 1, stats::sd)
  gene_sd[!is.finite(gene_sd) | gene_sd < 1e-8] <- 1
  z <- sweep(sweep(avg, 1, gene_center, "-"), 1, gene_sd, "/")

  rows <- list()
  k <- 1L
  for (label in names(marker_reference)) {
    genes <- intersect(marker_reference[[label]], rownames(z))
    if (!length(genes)) next
    marker_score <- colMeans(z[genes, , drop = FALSE])
    pct <- colMeans(detected[genes, , drop = FALSE])
    for (cluster in colnames(z)) {
      rows[[k]] <- data.frame(cluster = cluster, candidate = label,
                              score = marker_score[[cluster]], marker_detection = pct[[cluster]],
                              n_markers = length(genes), stringsAsFactors = FALSE)
      k <- k + 1L
    }
  }
  scores <- do.call(rbind, rows)
  assignments <- lapply(split(scores, scores$cluster), function(d) {
    d <- d[order(-d$score, -d$marker_detection), , drop = FALSE]
    best <- d[1, ]
    runner_up <- if (nrow(d) > 1) d$score[2] else -Inf
    margin <- best$score - runner_up
    label <- best$candidate
    if (best$n_markers < min_detected_markers || best$marker_detection < min_marker_detection || margin < min_margin) {
      label <- paste0("Uncertain_", best$candidate)
    }
    data.frame(cluster = best$cluster, cell_type = label, confidence_margin = margin,
               marker_score = best$score, marker_detection = best$marker_detection,
               stringsAsFactors = FALSE)
  })
  assignments <- do.call(rbind, assignments)
  list(assignments = assignments, scores = scores)
}

annotate_clusters <- function(object, tissue = c("Kidney", "PBMC"), output_dir = NULL) {
  tissue <- match.arg(tissue)
  reference <- if (tissue == "Kidney") kidney_marker_reference() else pbmc_marker_reference()
  result <- score_cluster_markers(object, reference)
  mapping <- setNames(result$assignments$cell_type, result$assignments$cluster)
  object$cell_type <- unname(mapping[as.character(object$seurat_clusters)])
  object$cell_class <- broad_cell_class(object$cell_type)
  Seurat::Idents(object) <- "cell_type"
  if (!is.null(output_dir)) {
    write_csv_safe(result$assignments, file.path(output_dir, paste0(tolower(tissue), "_cluster_annotations.csv")))
    write_csv_safe(result$scores, file.path(output_dir, paste0(tolower(tissue), "_annotation_scores_all.csv")))
  }
  attr(object, "annotation_table") <- result$assignments
  object
}
