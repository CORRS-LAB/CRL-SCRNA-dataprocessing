# End-to-end kidney/PBMC high-salt scRNA-seq workflow.
# Run from the repository root with:
#   Rscript src/analysis_demo.R
# The same file can also be sourced from RStudio.

locate_this_script <- function() {
  from_source <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(from_source)) return(normalizePath(from_source, winslash = "/", mustWork = FALSE))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
  NA_character_
}

script_file <- locate_this_script()
project_root <- if (file.exists(file.path("src", "utils.R"))) "." else if (!is.na(script_file)) dirname(dirname(script_file)) else getwd()
if (!file.exists(file.path(project_root, "src", "utils.R"))) project_root <- getwd()
if (!file.exists(file.path(project_root, "src", "utils.R"))) stop("Run this workflow from the repository root.")

source(file.path(project_root, "src", "utils.R"))
source(file.path(project_root, "src", "annotation.R"))
source(file.path(project_root, "src", "preprocess.R"))
source(file.path(project_root, "src", "mechanism_analysis.R"))
source(file.path(project_root, "src", "communication_analysis.R"))
source(file.path(project_root, "src", "pbmc_integration.R"))

parse_cli <- function(args) {
  defaults <- list(input = file.path(project_root, "data", "kidney_demo_input.rds"),
                   pbmc = file.path(project_root, "data", "pbmc_demo_input.rds"),
                   outdir = file.path(project_root, "results"),
                   seed = 2026L, generate_pbmc = TRUE, doublets = TRUE, markers = TRUE)
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) next
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", pair[1]); value <- paste(pair[-1], collapse = "=")
    if (!key %in% names(defaults)) next
    if (is.logical(defaults[[key]])) value <- tolower(value) %in% c("true", "1", "yes", "y")
    if (is.integer(defaults[[key]])) value <- as.integer(value)
    defaults[[key]] <- value
  }
  defaults
}

run_full_workflow <- function(config = parse_cli(commandArgs(trailingOnly = TRUE))) {
  set.seed(config$seed)
  assert_packages(c("Seurat", "SeuratObject", "Matrix", "edgeR", "ggplot2", "patchwork"))
  output_root <- safe_dir_create(config$outdir)
  if (!file.exists(config$input)) stop("Kidney input not found: ", config$input)
  message("[1/6] Loading and preprocessing kidney data")
  kidney_input <- readRDS(config$input)
  kidney <- preprocess_sc_object(kidney_input, tissue = "Kidney",
    output_dir = file.path(output_root, "01_kidney_preprocessing"),
    use_doublet_finder = config$doublets, seed = config$seed, find_markers = config$markers)

  message("[2/6] Running sample-aware differential expression and pathway analysis")
  mechanism <- run_mechanism_analysis(kidney, file.path(output_root, "02_high_salt_mechanisms"))

  message("[3/6] Preparing PBMC data")
  if (file.exists(config$pbmc)) {
    pbmc_input <- readRDS(config$pbmc)
  } else if (config$generate_pbmc) {
    message("PBMC input was not supplied; generating a clearly labelled simulation for code validation only.")
    pbmc_input <- generate_pbmc_demo(kidney_input, config$pbmc, seed = config$seed)
  } else {
    stop("PBMC input not found and --generate-pbmc=false: ", config$pbmc)
  }
  pbmc <- preprocess_sc_object(pbmc_input, tissue = "PBMC",
    output_dir = file.path(output_root, "03_pbmc_preprocessing"), use_doublet_finder = FALSE,
    seed = config$seed, resolution = 0.5, find_markers = config$markers)
  if ("simulated_cell_type" %in% colnames(pbmc[[]])) {
    pbmc$marker_score_cell_type <- pbmc$cell_type
    pbmc$cell_type <- pbmc$simulated_cell_type
    pbmc$cell_class <- broad_cell_class(pbmc$cell_type)
    pbmc$annotation_source <- "known_simulation_label_for_code_validation"
    Seurat::Idents(pbmc) <- "cell_type"
  }

  message("[4/6] Integrating kidney and PBMC response analyses")
  integration <- run_joint_kidney_pbmc_analysis(kidney, pbmc,
    file.path(output_root, "04_kidney_pbmc_integration"), seed = config$seed)

  message("[5/6] Prioritizing differential communication drivers")
  communication <- run_communication_analysis(kidney,
    file.path(output_root, "05_cell_communication"), top_n = 8L)

  message("[6/6] Saving reproducibility metadata")
  capture.output(sessionInfo(), file = file.path(output_root, "sessionInfo.txt"))
  config_table <- data.frame(parameter = names(config), value = vapply(config, as.character, character(1)))
  write_csv_safe(config_table, file.path(output_root, "run_config.csv"))
  saveRDS(list(kidney = kidney, pbmc = pbmc, mechanism = mechanism,
               integration = integration, communication = communication),
          file.path(output_root, "analysis_bundle.rds"))
  message("Workflow completed. Results: ", output_root)
  invisible(list(kidney = kidney, pbmc = pbmc, mechanism = mechanism,
                 integration = integration, communication = communication))
}

workflow_result <- run_full_workflow()
