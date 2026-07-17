ligand_receptor_reference <- function() {
  data.frame(
    pathway = c(
      rep("TGFB", 3), rep("CCL", 5), rep("CXCL", 5), rep("TNF", 3), rep("IL1", 2),
      rep("VEGF", 3), rep("ANGPT", 2), rep("NOTCH", 4), rep("WNT", 4), rep("FGF", 3),
      rep("PDGF", 3), rep("MIF", 2), rep("SPP1", 3), rep("CSF", 3), rep("OSM", 2)
    ),
    ligand = c(
      "Tgfb1", "Tgfb2", "Tgfb3", "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl7",
      "Cxcl1", "Cxcl2", "Cxcl9", "Cxcl10", "Cxcl12", "Tnf", "Tnfsf10", "Cd40lg",
      "Il1a", "Il1b", "Vegfa", "Vegfb", "Vegfc", "Angpt1", "Angpt2",
      "Dll1", "Dll4", "Jag1", "Jag2", "Wnt2", "Wnt4", "Wnt5a", "Wnt7b",
      "Fgf1", "Fgf2", "Fgf7", "Pdgfa", "Pdgfb", "Pdgfc", "Mif", "Mif",
      "Spp1", "Spp1", "Spp1", "Csf1", "Csf2", "Csf3", "Osm", "Osm"
    ),
    receptor = c(
      "Tgfbr1", "Tgfbr2", "Tgfbr2", "Ccr2", "Ccr1", "Ccr5", "Ccr5", "Ccr2",
      "Cxcr2", "Cxcr2", "Cxcr3", "Cxcr3", "Cxcr4", "Tnfrsf1a", "Tnfrsf10b", "Cd40",
      "Il1r1", "Il1r1", "Kdr", "Flt1", "Flt4", "Tek", "Tek",
      "Notch1", "Notch1", "Notch2", "Notch2", "Fzd4", "Fzd6", "Fzd2", "Fzd4",
      "Fgfr1", "Fgfr1", "Fgfr2", "Pdgfra", "Pdgfrb", "Pdgfra", "Cd74", "Cxcr4",
      "Cd44", "Itgav", "Itgb1", "Csf1r", "Csf2ra", "Csf3r", "Osmr", "Lifr"
    ), stringsAsFactors = FALSE
  )
}

calculate_sample_lr_edges <- function(object, lr_db = ligand_receptor_reference(),
                                      cell_type_col = "cell_type", min_pct = 0.05) {
  meta <- object[[]]
  key <- paste(meta$orig.ident, meta[[cell_type_col]], sep = "||")
  avg <- average_sparse_expression(object, key, "data")
  pct <- average_sparse_expression(object, key, "data", detection = TRUE)
  available <- lr_db$ligand %in% rownames(avg) & lr_db$receptor %in% rownames(avg)
  lr_db <- lr_db[available, , drop = FALSE]
  if (!nrow(lr_db)) stop("No ligand-receptor genes were found in the RNA assay.")
  key_parts <- strsplit(colnames(avg), "\\|\\|")
  group_info <- data.frame(key = colnames(avg),
    sample = vapply(key_parts, `[`, character(1), 1),
    cell_type = vapply(key_parts, `[`, character(1), 2), stringsAsFactors = FALSE)
  condition_map <- unique(meta[, c("orig.ident", "Treat"), drop = FALSE])
  rows <- list(); k <- 1L
  for (sample in unique(group_info$sample)) {
    gi <- group_info[group_info$sample == sample, , drop = FALSE]
    for (i in seq_len(nrow(lr_db))) {
      lig <- lr_db$ligand[i]; rec <- lr_db$receptor[i]
      lig_expr <- avg[lig, gi$key]; lig_pct <- pct[lig, gi$key]
      rec_expr <- avg[rec, gi$key]; rec_pct <- pct[rec, gi$key]
      score <- sqrt(outer(lig_expr, rec_expr, "*")) * sqrt(outer(lig_pct, rec_pct, "*"))
      score[outer(lig_pct, rec_pct, function(a, b) a < min_pct | b < min_pct)] <- 0
      grid <- expand.grid(sender = gi$cell_type, receiver = gi$cell_type, stringsAsFactors = FALSE)
      grid$score <- as.vector(score)
      grid$sample <- sample
      grid$Treat <- as.character(condition_map$Treat[match(sample, condition_map$orig.ident)])
      grid$pathway <- lr_db$pathway[i]
      grid$ligand <- lig
      grid$receptor <- rec
      rows[[k]] <- grid[grid$score > 0, ]; k <- k + 1L
    }
  }
  do.call(rbind, rows)
}

summarize_lr_edge_differences <- function(edges) {
  samples <- unique(edges[, c("sample", "Treat"), drop = FALSE])
  keys <- unique(edges[, c("pathway", "ligand", "receptor", "sender", "receiver"), drop = FALSE])
  template <- merge(keys, samples, by = NULL)
  full <- merge(template, edges, by = c("pathway", "ligand", "receptor", "sender", "receiver", "sample", "Treat"), all.x = TRUE)
  full$score[is.na(full$score)] <- 0
  groups <- split(full, interaction(full$pathway, full$ligand, full$receptor, full$sender, full$receiver, drop = TRUE))
  rows <- lapply(groups, function(d) {
    hs <- d$score[d$Treat == "HS"]; ls <- d$score[d$Treat == "LS"]
    if (length(hs) < 2L || length(ls) < 2L) return(NULL)
    delta <- mean(hs) - mean(ls)
    data.frame(pathway = d$pathway[1], ligand = d$ligand[1], receptor = d$receptor[1],
      sender = d$sender[1], receiver = d$receiver[1], mean_HS = mean(hs), mean_LS = mean(ls),
      delta = delta, consistency = pairwise_consistency(hs, ls, sign(delta)), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

prioritize_communication_pathways <- function(edges, edge_differences, top_n = 8L) {
  sample_totals <- stats::aggregate(score ~ sample + Treat + pathway, edges, sum)
  all_combinations <- merge(unique(edges[, c("pathway"), drop = FALSE]),
                            unique(edges[, c("sample", "Treat"), drop = FALSE]), by = NULL)
  sample_totals <- merge(all_combinations, sample_totals, by = c("pathway", "sample", "Treat"), all.x = TRUE)
  sample_totals$score[is.na(sample_totals$score)] <- 0
  rows <- lapply(split(sample_totals, sample_totals$pathway), function(d) {
    hs <- d$score[d$Treat == "HS"]; ls <- d$score[d$Treat == "LS"]
    delta <- mean(hs) - mean(ls)
    pooled_sd <- sqrt((stats::var(hs) + stats::var(ls)) / 2)
    if (!is.finite(pooled_sd) || pooled_sd < 1e-8) pooled_sd <- mean(c(hs, ls)) + 1e-6
    p <- tryCatch(stats::t.test(hs, ls)$p.value, error = function(e) 1)
    ed <- edge_differences[edge_differences$pathway == d$pathway[1], , drop = FALSE]
    support <- sum(ed$consistency >= 0.75 & sign(ed$delta) == sign(delta) & abs(ed$delta) > 0, na.rm = TRUE)
    data.frame(pathway = d$pathway[1], mean_HS = mean(hs), mean_LS = mean(ls),
      delta_HS_minus_LS = delta, standardized_effect = delta / pooled_sd,
      consistency = pairwise_consistency(hs, ls, sign(delta)), p_value = p,
      supporting_edges = support,
      affected_senders = length(unique(ed$sender[ed$consistency >= 0.75 & sign(ed$delta) == sign(delta)])),
      affected_receivers = length(unique(ed$receiver[ed$consistency >= 0.75 & sign(ed$delta) == sign(delta)])),
      lr_pairs_detected = length(unique(paste(ed$ligand, ed$receptor))), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$FDR <- stats::p.adjust(out$p_value, "BH")
  out$direction <- ifelse(out$delta_HS_minus_LS > 0, "Increased_in_HS", "Decreased_in_HS")
  # With very small n, an almost-zero within-group SD can create an unstable
  # standardized effect. Preserve the raw value but cap its contribution to rank.
  out$ranking_effect <- pmin(abs(out$standardized_effect), 5)
  out$driver_score <- out$ranking_effect * (0.5 + out$consistency) *
    log1p(out$supporting_edges) * (1 + pmin(-log10(pmax(out$FDR, 1e-12)), 6))
  out <- out[order(-out$driver_score), ]
  out$core_rank <- seq_len(nrow(out))
  core <- out[out$supporting_edges >= 1 & out$consistency >= 0.5, , drop = FALSE]
  core <- head(core, top_n)
  list(screening = out, core = core, sample_totals = sample_totals)
}

run_communication_analysis <- function(object, output_dir, top_n = 8L) {
  out <- safe_dir_create(output_dir)
  edges <- calculate_sample_lr_edges(object)
  edge_diff <- summarize_lr_edge_differences(edges)
  ranked <- prioritize_communication_pathways(edges, edge_diff, top_n)
  write_csv_safe(ranked$screening, file.path(out, "communication_pathway_screening_all.csv"))
  write_csv_safe(ranked$core, file.path(out, "core_driver_pathways_ranked.csv"))
  if (nrow(ranked$core)) {
    core_edges <- edge_diff[edge_diff$pathway %in% ranked$core$pathway & edge_diff$consistency >= 0.75, ]
    core_edges <- core_edges[order(-abs(core_edges$delta)), ]
    write_csv_safe(core_edges, file.path(out, "core_driver_supporting_edges.csv"))
    plot_data <- ranked$core
    plot_data$pathway <- factor(plot_data$pathway, levels = rev(plot_data$pathway))
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(pathway, driver_score, fill = direction)) +
      ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::theme_bw(base_size = 11) +
      ggplot2::labs(x = NULL, y = "Core driver score",
                    title = "Prioritized HS-versus-LS communication drivers")
    save_plot(p, file.path(out, "core_driver_pathways.png"), 8, 5)
  }
  saveRDS(list(core = ranked$core, pathway_screening = ranked$screening,
               supporting_edges = edge_diff, sample_scores = ranked$sample_totals),
          file.path(out, "communication_analysis.rds"))
  ranked
}
