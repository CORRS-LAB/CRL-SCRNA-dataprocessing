options(stringsAsFactors = FALSE)

assert_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required R package(s): ", paste(missing, collapse = ", "),
         ". Install them before running this workflow.")
  }
}

safe_dir_create <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create output directory: ", path)
  }
  # Keep relative paths relative. On Windows, R running in the C locale can fail
  # when a normalized absolute path contains non-ASCII directory names.
  path
}

write_csv_safe <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file = file, row.names = FALSE, na = "")
}

join_rna_layers <- function(object) {
  Seurat::DefaultAssay(object) <- "RNA"
  layers <- SeuratObject::Layers(object[["RNA"]])
  if (length(layers) > 1L || !all(c("counts", "data") %in% layers)) {
    object <- SeuratObject::JoinLayers(object, assay = "RNA")
  }
  object
}

get_assay_matrix <- function(object, layer = c("counts", "data")) {
  layer <- match.arg(layer)
  object <- join_rna_layers(object)
  SeuratObject::LayerData(object, assay = "RNA", layer = layer)
}

infer_condition <- function(sample_id) {
  id <- toupper(as.character(sample_id))
  condition <- ifelse(grepl("(^|[^A-Z])HS|KHS|HIGH", id), "HS",
                      ifelse(grepl("(^|[^A-Z])LS|KLS|LOW|NORMAL|NS", id), "LS", NA_character_))
  factor(condition, levels = c("LS", "HS"))
}

add_standard_metadata <- function(object, tissue = "Kidney") {
  if (!"orig.ident" %in% colnames(object[[]])) stop("Metadata column 'orig.ident' is required.")
  if (!"Treat" %in% colnames(object[[]])) object$Treat <- infer_condition(object$orig.ident)
  object$Treat <- factor(as.character(object$Treat), levels = c("LS", "HS"))
  if (anyNA(object$Treat)) {
    bad <- unique(as.character(object$orig.ident[is.na(object$Treat)]))
    stop("Could not infer HS/LS condition for sample(s): ", paste(bad, collapse = ", "))
  }
  object$tissue <- tissue
  object
}

make_group_matrix <- function(groups) {
  groups <- factor(as.character(groups), levels = unique(as.character(groups)))
  design <- Matrix::sparse.model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  design
}

aggregate_sparse_counts <- function(object, groups) {
  counts <- get_assay_matrix(object, "counts")
  if (length(groups) != ncol(counts)) stop("groups must have one value per cell")
  design <- make_group_matrix(groups)
  out <- counts %*% design
  colnames(out) <- colnames(design)
  out
}

average_sparse_expression <- function(object, groups, layer = "data", detection = FALSE) {
  mat <- get_assay_matrix(object, layer)
  design <- make_group_matrix(groups)
  denom <- Matrix::colSums(design)
  source <- if (detection) (mat > 0) * 1 else mat
  out <- source %*% design
  out <- sweep(as.matrix(out), 2, denom, "/")
  colnames(out) <- colnames(design)
  out
}

scaled_priority <- function(x) {
  x <- as.numeric(x)
  center <- stats::median(x, na.rm = TRUE)
  spread <- stats::mad(x, center = center, na.rm = TRUE)
  if (!is.finite(spread) || spread < 1e-8) spread <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(spread) || spread < 1e-8) spread <- 1
  abs(x - center) / spread
}

pairwise_consistency <- function(hs, ls, expected_sign) {
  differences <- as.vector(outer(hs, ls, "-"))
  if (!length(differences) || expected_sign == 0) return(0)
  mean(sign(differences) == expected_sign, na.rm = TRUE)
}

save_plot <- function(plot, filename, width = 10, height = 6) {
  # Avoid ragg here: it cannot write through some Windows C-locale/non-ASCII paths.
  grDevices::png(filename = filename, width = width * 300, height = height * 300,
                 units = "px", res = 300)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  invisible(filename)
}
