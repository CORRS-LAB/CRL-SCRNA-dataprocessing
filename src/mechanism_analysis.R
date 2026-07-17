high_salt_gene_sets <- function() {
  list(
    "Osmotic_stress_NFAT5" = c("Nfat5", "Akr1b1", "Slc5a3", "Smit1", "TauT", "Hspa1a", "Hspa1b", "Hsp90aa1"),
    "Sodium_transport_ENaC" = c("Scnn1a", "Scnn1b", "Scnn1g", "Sgk1", "Nedd4l", "Fxyd4", "Atp1a1", "Atp1b1"),
    "RAAS_aldosterone" = c("Ren", "Agt", "Ace", "Ace2", "Agtr1a", "Agtr2", "Cyp11b2", "Nr3c2", "Sgk1"),
    "Oxidative_stress" = c("Nox1", "Nox2", "Nox4", "Cyba", "Ncf1", "Nrf2", "Nfe2l2", "Hmox1", "Sod1", "Sod2", "Gpx1"),
    "NFkB_inflammation" = c("Nfkb1", "Rela", "Ikbkb", "Nfkbia", "Tnf", "Il1b", "Il6", "Ccl2", "Cxcl2"),
    "Interferon_response" = c("Stat1", "Stat2", "Irf7", "Isg15", "Ifit1", "Ifit2", "Ifit3", "Mx1", "Oas1a"),
    "TGFb_fibrosis" = c("Tgfb1", "Tgfb2", "Tgfbr1", "Tgfbr2", "Smad2", "Smad3", "Ctgf", "Col1a1", "Col3a1", "Fn1", "Postn"),
    "ECM_remodeling" = c("Col1a1", "Col1a2", "Col3a1", "Col4a1", "Fn1", "Mmp2", "Mmp9", "Timp1", "Lox", "Dcn"),
    "Endothelial_dysfunction" = c("Nos3", "Klf2", "Klf4", "Edn1", "Ednra", "Vcam1", "Icam1", "Sele", "Kdr"),
    "Hypoxia_HIF1" = c("Hif1a", "Epas1", "Vegfa", "Slc2a1", "Ldha", "Pdk1", "Bnip3", "EglN1"),
    "Macrophage_activation" = c("Adgre1", "Cd68", "Lyz2", "C1qa", "C1qb", "Trem2", "Spp1", "Nos2", "Arg1"),
    "Apoptosis_injury" = c("Bax", "Bcl2", "Casp3", "Casp8", "Fas", "Fasl", "Havcr1", "Lcn2", "Krt8", "Krt18")
  )
}

run_edger_one_celltype <- function(counts, sample_info) {
  condition <- factor(sample_info$Treat, levels = c("LS", "HS"))
  design <- stats::model.matrix(~ condition)
  y <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(y, design = design, min.count = 3)
  if (sum(keep) < 10L) return(NULL)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  robust <- requireNamespace("statmod", quietly = TRUE)
  y <- edgeR::estimateDisp(y, design, robust = robust)
  fit <- edgeR::glmQLFit(y, design, robust = robust)
  test <- edgeR::glmQLFTest(fit, coef = "conditionHS")
  tab <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
  tab$gene <- rownames(tab)
  rownames(tab) <- NULL
  tab[, c("gene", setdiff(colnames(tab), "gene")), drop = FALSE]
}

run_pseudobulk_de <- function(object, output_dir, min_cells_per_sample = 10L,
                              cell_type_col = "cell_type") {
  assert_packages(c("Seurat", "SeuratObject", "Matrix", "edgeR"))
  meta <- object[[]]
  required <- c("orig.ident", "Treat", cell_type_col)
  if (!all(required %in% colnames(meta))) stop("Missing metadata for pseudobulk DE: ", paste(setdiff(required, colnames(meta)), collapse = ", "))
  results <- list()
  diagnostics <- list()
  k <- 1L
  for (cell_type in sort(unique(as.character(meta[[cell_type_col]])))) {
    cells <- rownames(meta)[as.character(meta[[cell_type_col]]) == cell_type]
    cell_meta <- meta[cells, , drop = FALSE]
    per_sample <- table(cell_meta$orig.ident)
    valid_samples <- names(per_sample)[per_sample >= min_cells_per_sample]
    info <- unique(cell_meta[cell_meta$orig.ident %in% valid_samples, c("orig.ident", "Treat"), drop = FALSE])
    info <- info[!duplicated(info$orig.ident), , drop = FALSE]
    n_hs <- sum(as.character(info$Treat) == "HS")
    n_ls <- sum(as.character(info$Treat) == "LS")
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(cell_type = cell_type,
      n_cells = length(cells), valid_samples = length(valid_samples), n_HS = n_hs, n_LS = n_ls)
    if (n_hs < 2L || n_ls < 2L) next
    keep_cells <- rownames(cell_meta)[cell_meta$orig.ident %in% valid_samples]
    sub <- subset(object, cells = keep_cells)
    counts <- aggregate_sparse_counts(sub, sub$orig.ident)
    info <- info[match(colnames(counts), info$orig.ident), , drop = FALSE]
    de <- tryCatch(run_edger_one_celltype(counts, info), error = function(e) {
      warning("edgeR failed for ", cell_type, ": ", conditionMessage(e)); NULL
    })
    if (is.null(de)) next
    de$cell_type <- cell_type
    de$n_HS <- n_hs
    de$n_LS <- n_ls
    de$method <- "edgeR_QL_pseudobulk"
    results[[k]] <- de
    k <- k + 1L
  }
  diagnostics <- do.call(rbind, diagnostics)
  de_all <- if (length(results)) do.call(rbind, results) else data.frame()
  write_csv_safe(diagnostics, file.path(output_dir, "pseudobulk_DE_diagnostics.csv"))
  write_csv_safe(de_all, file.path(output_dir, "pseudobulk_DE_all.csv"))
  if (nrow(de_all)) {
    sig <- de_all[is.finite(de_all$FDR) & de_all$FDR < 0.05 & abs(de_all$logFC) >= 0.5, ]
    sig <- sig[order(sig$FDR, -abs(sig$logFC)), ]
    write_csv_safe(sig, file.path(output_dir, "pseudobulk_DE_significant.csv"))
  }
  de_all
}

calculate_pathway_activity <- function(object, gene_sets, cell_type_col = "cell_type") {
  meta <- object[[]]
  key <- paste(meta$orig.ident, meta[[cell_type_col]], sep = "||")
  avg <- average_sparse_expression(object, key, layer = "data")
  gene_mean <- rowMeans(avg)
  gene_sd <- apply(avg, 1, stats::sd)
  gene_sd[!is.finite(gene_sd) | gene_sd < 1e-8] <- 1
  z <- sweep(sweep(avg, 1, gene_mean, "-"), 1, gene_sd, "/")
  key_parts <- strsplit(colnames(avg), "\\|\\|")
  key_meta <- data.frame(sample = vapply(key_parts, `[`, character(1), 1),
                         cell_type = vapply(key_parts, `[`, character(1), 2),
                         stringsAsFactors = FALSE)
  sample_condition <- unique(meta[, c("orig.ident", "Treat"), drop = FALSE])
  key_meta$Treat <- as.character(sample_condition$Treat[match(key_meta$sample, sample_condition$orig.ident)])
  rows <- list(); k <- 1L
  for (pathway in names(gene_sets)) {
    genes <- intersect(gene_sets[[pathway]], rownames(z))
    if (length(genes) < 2L) next
    scores <- colMeans(z[genes, , drop = FALSE])
    coverage <- length(genes) / length(unique(gene_sets[[pathway]]))
    rows[[k]] <- transform(key_meta, pathway = pathway, score = as.numeric(scores),
                           genes_detected = length(genes), gene_coverage = coverage)
    k <- k + 1L
  }
  do.call(rbind, rows)
}

compare_pathway_activity <- function(activity) {
  groups <- split(activity, interaction(activity$cell_type, activity$pathway, drop = TRUE))
  rows <- lapply(groups, function(d) {
    hs <- d$score[d$Treat == "HS"]; ls <- d$score[d$Treat == "LS"]
    if (length(hs) < 2L || length(ls) < 2L) return(NULL)
    delta <- mean(hs) - mean(ls)
    p <- tryCatch(stats::t.test(hs, ls)$p.value, error = function(e) 1)
    data.frame(cell_type = d$cell_type[1], pathway = d$pathway[1],
      mean_HS = mean(hs), mean_LS = mean(ls), delta_HS_minus_LS = delta,
      consistency = pairwise_consistency(hs, ls, sign(delta)), activity_p = p,
      genes_detected = d$genes_detected[1], gene_coverage = d$gene_coverage[1])
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame())
  out$activity_FDR <- stats::p.adjust(out$activity_p, method = "BH")
  rownames(out) <- NULL
  out
}

add_de_pathway_evidence <- function(pathway_table, de, gene_sets) {
  if (!nrow(pathway_table)) return(pathway_table)
  rows <- lapply(seq_len(nrow(pathway_table)), function(i) {
    row <- pathway_table[i, ]
    d <- de[de$cell_type == row$cell_type, , drop = FALSE]
    genes <- intersect(gene_sets[[row$pathway]], d$gene)
    if (length(genes) < 2L || nrow(d) < 20L) {
      return(data.frame(de_gene_count = length(genes), de_sig_count = 0,
                        mean_gene_logFC = NA_real_, enrichment_p = 1))
    }
    statistic <- sign(d$logFC) * -log10(pmax(d$PValue, 1e-300))
    inside <- d$gene %in% genes
    p <- tryCatch(stats::wilcox.test(statistic[inside], statistic[!inside], exact = FALSE)$p.value,
                  error = function(e) 1)
    data.frame(de_gene_count = length(genes),
      de_sig_count = sum(d$FDR[inside] < 0.1 & abs(d$logFC[inside]) >= 0.25, na.rm = TRUE),
      mean_gene_logFC = mean(d$logFC[inside], na.rm = TRUE), enrichment_p = p)
  })
  evidence <- do.call(rbind, rows)
  out <- cbind(pathway_table, evidence)
  fisher_stat <- -2 * (log(pmax(out$activity_p, 1e-300)) + log(pmax(out$enrichment_p, 1e-300)))
  out$combined_p <- stats::pchisq(fisher_stat, df = 4, lower.tail = FALSE)
  out$combined_FDR <- stats::p.adjust(out$combined_p, method = "BH")
  out$priority_score <- ave(abs(out$delta_HS_minus_LS), out$cell_type, FUN = scaled_priority) *
    (0.5 + out$consistency) * sqrt(out$gene_coverage) * (1 + pmin(-log10(pmax(out$combined_FDR, 1e-12)), 6))
  out$direction <- ifelse(out$delta_HS_minus_LS > 0, "Activated_in_HS", "Suppressed_in_HS")
  out[order(-out$priority_score), ]
}

run_mechanism_analysis <- function(object, output_dir, top_n_per_celltype = 5L) {
  out <- safe_dir_create(output_dir)
  de <- run_pseudobulk_de(object, out)
  gene_sets <- high_salt_gene_sets()
  activity <- calculate_pathway_activity(object, gene_sets)
  comparison <- compare_pathway_activity(activity)
  prioritized <- add_de_pathway_evidence(comparison, de, gene_sets)
  write_csv_safe(activity, file.path(out, "pathway_activity_by_sample_celltype.csv"))
  write_csv_safe(prioritized, file.path(out, "high_salt_mechanism_pathways_ranked.csv"))
  if (nrow(prioritized)) {
    top <- do.call(rbind, lapply(split(prioritized, prioritized$cell_type), function(d) head(d, top_n_per_celltype)))
    write_csv_safe(top, file.path(out, "high_salt_mechanism_top_pathways.csv"))
    p <- ggplot2::ggplot(top, ggplot2::aes(x = stats::reorder(pathway, priority_score), y = priority_score,
                                           fill = direction)) +
      ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::facet_wrap(~cell_type, scales = "free_y") +
      ggplot2::labs(x = NULL, y = "Mechanism priority score", title = "Prioritized high-salt mechanisms") +
      ggplot2::theme_bw(base_size = 10)
    save_plot(p, file.path(out, "high_salt_mechanism_priority.png"), 13, 9)
  }
  list(de = de, pathway_activity = activity, pathways = prioritized)
}
