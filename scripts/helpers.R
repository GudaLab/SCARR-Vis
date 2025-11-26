`%||%` <- function(a, b) if (is.null(a)) b else a

# -----------------------------
# Helpers
# -----------------------------
read_10x_any <- function(path_or_zip) {
  stopifnot(file.exists(path_or_zip))
  ext <- tolower(file_ext(path_or_zip))
  if (file.info(path_or_zip)$isdir) return(Seurat::Read10X(path_or_zip))
  if (ext == "zip") {
    td <- tempfile(pattern = "unz_10x_"); dir.create(td)
    utils::unzip(path_or_zip, exdir = td)
    cand <- list.dirs(td, recursive = TRUE, full.names = TRUE)
    has_mtx <- vapply(
      cand,
      function(d) any(grepl("matrix[.]mtx([.]gz)?$", list.files(d, full.names = TRUE))),
      logical(1)
    )
    dirs <- cand[has_mtx]
    if (length(dirs) == 0) stop("Could not find 10x files inside the zip.")
    return(Seurat::Read10X(dirs[[1]]))
  }
  if (ext == "h5") return(Seurat::Read10X_h5(path_or_zip))
  stop("Unsupported file type: ", ext)
}

makeSoupX <- function(raw_counts, filt_counts) {
  sc <- SoupX::SoupChannel(tod = raw_counts, toc = filt_counts)
  if (length(intersect(colnames(raw_counts), colnames(filt_counts))) == 0)
    warning("No overlapping barcodes between raw and filtered matrices. Proceeding anyway.")
  sc
}

map_ensembl_to_symbol <- function(ids, species_genome) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) return(setNames(ids, NULL))
  dataset <- switch(
    species_genome,
    "Human (GRCh38)" = "hsapiens_gene_ensembl",
    "Mouse (GRCm39)" = "mmusculus_gene_ensembl",
    NULL
  )
  if (is.null(dataset)) return(setNames(ids, NULL))
  ids_stripped <- sub("[.].*$", "", ids)
  mart <- NULL
  try({ mart <- biomaRt::useEnsembl("ensembl", dataset = dataset) }, silent = TRUE)
  if (is.null(mart)) try({
    mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = dataset,
                             host = "https://www.ensembl.org")
  }, silent = TRUE)
  if (is.null(mart)) return(setNames(ids, NULL))
  tbl <- tryCatch({
    biomaRt::getBM(c("ensembl_gene_id","external_gene_name"), "ensembl_gene_id",
                   unique(ids_stripped), mart)
  }, error = function(e) NULL)
  if (is.null(tbl) || nrow(tbl) == 0) return(setNames(ids, NULL))
  map <- setNames(tbl$external_gene_name, tbl$ensembl_gene_id)
  unname(ifelse(ids_stripped %in% names(map) & nzchar(map[ids_stripped]),
                map[ids_stripped], ids))
}

collapse_duplicated_rows <- function(m) {
  rn <- rownames(m)
  f <- factor(rn, levels = unique(rn))
  Mt <- as(m, "dgTMatrix")
  i_new <- as.integer(f)[Mt@i + 1]
  res <- Matrix::sparseMatrix(i = i_new, j = Mt@j + 1, x = Mt@x,
                              dims = c(length(levels(f)), ncol(m)))
  rownames(res) <- levels(f); colnames(res) <- colnames(m); res
}

species_defaults <- function(species_genome) {
  if (identical(species_genome, "Mouse (GRCm39)")) {
    list(
      mito_regex      = "^mt-",
      default_nonexp  = c("Erdr1","Ccl5","Ccl4","Lyz2","Prg4")  # example defaults
    )
  } else if (identical(species_genome, "Human (GRCh38)")) {
    list(
      mito_regex      = "^MT-",
      default_nonexp  = c("HBB","HBA1","HBA2","IGKC","PRM1")    # example defaults
    )
  } else {
    list(
      mito_regex      = "^MT-",
      default_nonexp  = character(0)
    )
  }
}

theme_clean <- ggplot2::theme(
  panel.background = ggplot2::element_blank(),
  panel.border     = ggplot2::element_rect(fill = NA),
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank(),
  strip.background = ggplot2::element_blank(),
  axis.text.x      = ggplot2::element_text(color = "black"),
  axis.text.y      = ggplot2::element_text(color = "black"),
  axis.ticks       = ggplot2::element_line(color = "black"),
  plot.margin      = grid::unit(c(1, 1,  1, 1), "lines")
)

write_10x_h5 <- function(mat, file) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) stop("Package 'rhdf5' is required to write .h5")
  mat <- as(mat, "dgCMatrix")
  rhdf5::h5createFile(file)
  rhdf5::h5createGroup(file, "matrix")
  rhdf5::h5write(as.numeric(mat@x), file, "matrix/data")
  rhdf5::h5write(as.integer(mat@i), file, "matrix/indices")
  rhdf5::h5write(as.integer(mat@p), file, "matrix/indptr")
  rhdf5::h5write(as.integer(c(nrow(mat), ncol(mat))), file, "matrix/shape")
  rhdf5::h5write(as.character(colnames(mat)), file, "matrix/barcodes")
  rhdf5::h5createGroup(file, "matrix/features")
  gene_ids <- rownames(mat); gene_names <- rownames(mat)
  rhdf5::h5write(as.character(gene_names), file, "matrix/features/name")
  rhdf5::h5write(as.character(gene_ids),   file, "matrix/features/id")
  rhdf5::h5write(rep("Gene Expression", length(gene_ids)), file, "matrix/features/feature_type")
  rhdf5::H5close()
}

round_sparse <- function(x) {
  x <- as(x, "dgCMatrix")
  if (length(x@x)) x@x <- round(x@x)
  x
}

# ---- Safety helpers ----
safe_rowsums <- function(m) {
  if (is.null(m)) return(NULL)
  dm <- dim(m)
  if (is.null(dm) || any(dm == 0)) return(NULL)
  Matrix::rowSums(m)
}
safe_colsums <- function(m, logical_gt0 = FALSE) {
  if (is.null(m)) return(NULL)
  dm <- dim(m)
  if (is.null(dm) || any(dm == 0)) return(NULL)
  if (logical_gt0) return(Matrix::colSums(m > 0))
  Matrix::colSums(m)
}

# ---- SoupX profile helpers ----
compute_soup_counts <- function(sc, raw_mat = NULL, filt_mat = NULL) {
  gsum <- safe_rowsums(sc$tod)
  if (is.null(gsum)) gsum <- safe_rowsums(raw_mat)
  if (is.null(gsum)) gsum <- safe_rowsums(filt_mat)
  genes <- names(gsum)
  if (is.null(genes) || any(!nzchar(genes))) {
    genes <- rownames(sc$tod)
    if (is.null(genes) || length(genes) == 0) genes <- rownames(raw_mat)
    if (is.null(genes) || length(genes) == 0) genes <- rownames(filt_mat)
  }
  if (is.null(genes) || length(genes) == 0) genes <- "GENE1"
  if (is.null(gsum) || !is.finite(sum(gsum)) || sum(gsum) == 0) {
    gsum <- setNames(rep(1, length(genes)), genes)
  } else if (is.null(names(gsum))) {
    names(gsum) <- genes
  }
  gsum
}
soup_profile_df <- function(sc, raw_mat = NULL, filt_mat = NULL) {
  counts <- compute_soup_counts(sc, raw_mat, filt_mat)
  est <- as.numeric(counts) / sum(counts)
  data.frame(counts = as.numeric(counts), est = est,
             row.names = names(counts), check.names = FALSE, stringsAsFactors = FALSE)
}

# Misc
clean_gene_list <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  unique(trimws(unlist(strsplit(x, "[,;\\s]+"))))
}


# -----------------------------
# FastCAR helpers
# -----------------------------

fastcar_available <- function() {
  suppressWarnings(
    requireNamespace("FastCAR", quietly = TRUE) &&
      requireNamespace("Matrix", quietly = TRUE)
  )
}

# internal: sparse-friendly background computation mimicking FastCAR logic
scarr_fastcar_determine_background <- function(full_mat,
                                               empty_droplet_cutoff,
                                               contamination_chance_cutoff) {
  stopifnot(methods::is(full_mat, "dgCMatrix"))
  
  # which barcodes are "empty" but not completely unused
  cs <- Matrix::colSums(full_mat)
  empties <- cs < empty_droplet_cutoff & cs > 0
  if (!any(empties)) {
    bg <- numeric(nrow(full_mat))
    names(bg) <- rownames(full_mat)
    return(bg)
  }
  
  sub_mat <- full_mat[, empties, drop = FALSE]
  
  # row-wise maxima in empty droplets
  bg_max <- NULL
  if (requireNamespace("qlcMatrix", quietly = TRUE)) {
    bg_max <- as.numeric(qlcMatrix::rowMax(sub_mat))
  } else {
    # fallback: dense; OK for small datasets
    warning("Package 'qlcMatrix' not available; using dense matrix for FastCAR background computation. This may use a lot of memory for large data.")
    bg_max <- apply(as.matrix(sub_mat), 1L, max)
  }
  names(bg_max) <- rownames(full_mat)
  
  # how often each gene shows up in empties
  occurrences <- Matrix::rowSums(sub_mat != 0)
  n_empty     <- sum(empties)
  if (!is.finite(n_empty) || n_empty == 0) {
    bg_max[] <- 0
    return(bg_max)
  }
  
  prob_contam <- occurrences / n_empty
  bg_max[prob_contam < contamination_chance_cutoff] <- 0
  bg_max
}

scarr_fastcar_wrapper <- function(full_mat, cell_mat,
                                  empty_droplet_cutoff = 100L,
                                  contamination_chance_cutoff = 0.05,
                                  do_profile = TRUE,
                                  profile_start = 10L, profile_stop = 500L, profile_by = 10L,
                                  use_recommended_cutoff = TRUE,
                                  seed = 1L) {
  stopifnot(methods::is(full_mat, "dgCMatrix"),
            methods::is(cell_mat, "dgCMatrix"))
  if (!fastcar_available()) {
    stop("FastCAR is not available. Please install the 'FastCAR' package in this R environment.")
  }
  set.seed(seed)
  
  amb_profile <- NULL
  if (isTRUE(do_profile)) {
    amb_profile <- FastCAR::describe.ambient.RNA.sequence(
      fullCellMatrix           = full_mat,
      start                    = profile_start,
      stop                     = profile_stop,
      by                       = profile_by,
      contaminationChanceCutoff = contamination_chance_cutoff
    )
    if (isTRUE(use_recommended_cutoff)) {
      rec <- tryCatch(
        FastCAR::recommend.empty.cutoff(amb_profile),
        error = function(e) NA_integer_
      )
      if (is.finite(rec) && rec > 0) {
        empty_droplet_cutoff <- rec
      }
    }
  }
  
  # use our sparse-friendly background computation (avoids rowQ / Biobase issues)
  background_to_remove <- scarr_fastcar_determine_background(
    full_mat,
    empty_droplet_cutoff = empty_droplet_cutoff,
    contamination_chance_cutoff = contamination_chance_cutoff
  )
  
  # then apply FastCAR's efficient sparse subtraction
  corrected <- FastCAR::remove.background(cell_mat, background_to_remove)
  
  delta <- cell_mat - corrected
  reads_removed_per_cell <- Matrix::colSums(delta)
  reads_removed_per_gene <- Matrix::rowSums(delta)
  
  list(
    method           = "FastCAR",
    params           = list(
      emptyDropletCutoff        = empty_droplet_cutoff,
      contaminationChanceCutoff = contamination_chance_cutoff,
      profiled                  = isTRUE(do_profile),
      profile_grid              = if (isTRUE(do_profile)) list(start = profile_start,
                                                               stop  = profile_stop,
                                                               by    = profile_by) else NULL,
      usedRecommendedCutoff     = isTRUE(use_recommended_cutoff)
    ),
    corrected_counts = corrected,
    aux              = list(
      amb_profile            = amb_profile,
      reads_removed_per_cell = reads_removed_per_cell,
      reads_removed_per_gene = reads_removed_per_gene
    )
  )
}

# build a ggplot / patchwork version of the ambient profile
fastcar_ambient_profile_plot <- function(amb_profile) {
  if (is.null(amb_profile)) return(NULL)
  
  df <- as.data.frame(amb_profile)
  
  # old FastCAR: rownames are cutoffs, 3 unnamed columns
  # new FastCAR: has columns cutoffValue, nEmptyDroplets, genesInBackground, genesContaminating
  if ("cutoffValue" %in% colnames(df)) {
    cutoff            <- df$cutoffValue
    nEmptyDroplets    <- df$nEmptyDroplets
    genesInBackground <- df$genesInBackground
    genesContaminating<- df$genesContaminating
  } else {
    cutoff            <- as.numeric(rownames(df))
    nEmptyDroplets    <- df[[1]]
    genesInBackground <- df[[2]]
    genesContaminating<- df[[3]]
  }
  
  d <- data.frame(
    cutoff            = as.numeric(cutoff),
    nEmptyDroplets    = as.numeric(nEmptyDroplets),
    genesInBackground = as.numeric(genesInBackground),
    genesContaminating= as.numeric(genesContaminating)
  )
  
  p_top <- ggplot2::ggplot(d, ggplot2::aes(x = cutoff, y = nEmptyDroplets)) +
    ggplot2::geom_point() +
    ggplot2::ggtitle("Total number of empty droplets at cutoffs") +
    ggplot2::xlab("empty droplet UMI cutoff") +
    ggplot2::ylab("Number of empty droplets") +
    theme_clean
  
  p_mid <- ggplot2::ggplot(d, ggplot2::aes(x = cutoff, y = genesInBackground)) +
    ggplot2::geom_point() +
    ggplot2::ggtitle("Number of genes in ambient RNA") +
    ggplot2::xlab("empty droplet UMI cutoff") +
    ggplot2::ylab("Genes in empty droplets") +
    theme_clean
  
  p_bot <- ggplot2::ggplot(d, ggplot2::aes(x = cutoff, y = genesContaminating)) +
    ggplot2::geom_point() +
    ggplot2::ggtitle("number of genes to correct") +
    ggplot2::xlab("empty droplet UMI cutoff") +
    ggplot2::ylab("Genes identified as contamination") +
    theme_clean
  
  # use patchwork (already used elsewhere in the app) -> no grid.arrange needed
  patchwork::wrap_plots(p_top, p_mid, p_bot, ncol = 1)
}

