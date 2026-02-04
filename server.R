
# -----------------------------
# Server
# -----------------------------
server <- function(input, output, session) {
  
  ###### session info
  
  sess_txt <- reactive({
    paste(capture.output(utils::sessionInfo()), collapse = "\n")
  })
  
  output$sess <- renderPrint({
    cat(sess_txt())
  })
  
  output$download_sess <- downloadHandler(
    filename = function() {
      paste0(
        "SCARR-Vis_session-info_",
        format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
        ".txt"
      )
    },
    content = function(file) {
      writeLines(sess_txt(), con = file, useBytes = TRUE)
    }
  )
  
  # Timeout
  observeEvent(input$timeOut, { 
    print(paste0("Session (", session$token, ") timed out at: ", Sys.time()))
    showModal(modalDialog(
      title = "Timeout",
      paste("Session timeout due to", input$timeOut, "inactivity -", Sys.time()),
      footer = NULL
    ))
    session$close()
  })
  
  # view counter
  session$onFlushed(function() {
    current <- increment_count()
    output$view_count <- renderText(format(current, big.mark = ","))
  }, once = TRUE)
  
  
  # hide analysis tabs until Run is clicked
  analysis_tabs <- c(
    "QC (Pre)",
    "Estimation",
    "QC (Post)",
    "Cluster counts (Pre vs post)",
    "Cells Table",
    "Top Genes",
    "UMAP/TSNE (Pre vs Post)",
    "Heatmap (Pre vs Post)",
    "Feature Plot (Pre vs Post)",
    "Session Info"
  )
  
  # Hide analysis tabs on initial load
  observe({
    lapply(
      analysis_tabs,
      function(tb) shinyjs::hide(selector = sprintf('a[data-value="%s"]', tb))
    )
  })
  
  # Show all analysis tabs after Run is clicked
  observeEvent(input$run, {
    lapply(
      analysis_tabs,
      function(tb) shinyjs::show(selector = sprintf('a[data-value="%s"]', tb))
    )
    # jump to QC (Pre) after run:
    updateTabsetPanel(session, "tabs", selected = "QC (Pre)")
  }, ignoreInit = TRUE)
  
  
  # ---- upload validators + reset ----
  is_valid_10x <- function(fname, kind = c("raw","filtered")) {
    kind <- match.arg(kind)
    pat <- if (kind == "raw") "^raw_feature_bc_matrix\\.(h5|zip)$"
    else                      "^filtered_feature_bc_matrix\\.(h5|zip)$"
    grepl(pat, basename(fname), ignore.case = FALSE)
  }
  
  same_ext <- function(a, b) tools::file_ext(a) == tools::file_ext(b)
  
  hard_reset <- function(session, msg = "Resetting the app…") {
    showNotification(msg, type = "error", duration = 3)
    later::later(function() {
      if (is.function(session$reload)) {
        session$reload()
      } else {
        session$sendCustomMessage("force-reload", TRUE)
      }
    }, 2.5)
  }
  
  rv <- reactiveValues(
    raw=NULL, filt=NULL, sc=NULL, auto=NULL, rho=NULL,
    adj_counts=NULL, seu_pre=NULL, seu_post=NULL,
    cells_df=NULL, cluster_counts_df=NULL,
    upload_type=NULL, log=character(),
    decontx_sce=NULL,
    GCGs=NULL, sccdc_quant=NULL,
    scCDC_pdf_abs = NULL,
    scCDC_pdf_url = NULL
  )
  
  # Example files
  ex_dir   <- file.path("www", "example", "GSM7681687")
  ex_raw   <- file.path(ex_dir, "raw_feature_bc_matrix.h5")
  ex_filt  <- file.path(ex_dir, "filtered_feature_bc_matrix.h5")
  
  clear_data <- function() {
    rv$raw <- rv$filt <- rv$sc <- rv$auto <- rv$rho <- rv$adj_counts <- NULL
    rv$seu_pre <- rv$seu_post <- rv$cells_df <- rv$cluster_counts_df <- NULL
    rv$decontx_sce <- NULL
  }
  
  # Serve scCDC PDF reports
  reports_dir <- file.path(tempdir(), "scCDC_reports")
  dir.create(reports_dir, showWarnings = FALSE, recursive = TRUE)
  shiny::addResourcePath("reports", reports_dir)
  
  clear_scCDC_pdfs <- function(reports_dir) {
    old_reports <- list.files(reports_dir, pattern = "[.]pdf$", full.names = TRUE)
    if (length(old_reports)) try(unlink(old_reports, force = TRUE), silent = TRUE)
    mis_pat  <- "^scCDC_reportsdefault_SE-plot.*[.]pdf$"
    mis_pdfs <- list.files(tempdir(), pattern = mis_pat, full.names = TRUE)
    if (length(mis_pdfs)) try(unlink(mis_pdfs, force = TRUE), silent = TRUE)
  }
  
  log_msg <- function(...) {
    rv$log <- c(rv$log, paste0(format(Sys.time(), "%H:%M:%S"), " - ", paste(..., collapse = " ")))
  }
  
  # Reset to Status whenever method changes
  observeEvent(input$method, ignoreInit = TRUE, {
    updateTabsetPanel(session, "tabs", selected = "Status")
  })
  
  observeEvent(input$species_genome, ignoreInit = TRUE, {
    defs <- species_defaults(input$species_genome)
    updateTextInput(session, "nonexp", value = paste(defs$default_nonexp, collapse = ","))
    log_msg("Species set to:", input$species_genome)
  })
  
  observeEvent(input$filt, ignoreInit = TRUE, {
    if (!is.null(input$filt$name)) {
      rv$upload_type <- if (grepl("[.]h5$", input$filt$name, ignore.case = TRUE)) "h5" else "mtx"
    }
  })
  
  output$dl_adjusted_ui <- renderUI({
    label <- if (identical(rv$upload_type, "h5")) "Download Cleaned (.h5)" else "Download Cleaned (MTX zip)"
    downloadBttn("dl_adjusted_data", label)
  })
  
  # Data source toggle (example vs upload)
  observeEvent(input$data_mode, ignoreInit = TRUE, {
    if (identical(input$data_mode, "example")) {
      clear_data()
      if (!file.exists(ex_raw) || !file.exists(ex_filt)) {
        showNotification(paste0(
          "Example files not found. Please ensure they exist at: ",
          normalizePath(ex_dir, mustWork = FALSE)
        ), type = "error", duration = 6)
        return()
      }
      log_msg("Using example data from:", normalizePath(ex_dir, mustWork = FALSE))
      ok <- TRUE
      try({
        rv$raw  <- read_10x_any(ex_raw)
        rv$filt <- read_10x_any(ex_filt)
        rv$upload_type <- "h5"
      }, silent = TRUE)
      if (is.null(rv$raw) || is.null(rv$filt)) ok <- FALSE
      if (!ok) {
        showNotification("Failed to read example matrices.", type = "error")
        return()
      }
      log_msg("Example RAW loaded:",  nrow(rv$raw),  "genes x", ncol(rv$raw),  "cells")
      log_msg("Example FILTERED loaded:", nrow(rv$filt), "genes x", ncol(rv$filt), "cells")
    } else {
      clear_data()
      rv$upload_type <- NULL
      log_msg("Switched to Upload mode; please select your files.")
    }
  })
  
  # Upload handling
  observeEvent({
    input$raw
    input$filt
  }, {
    req(input$raw, input$filt)
    raw_name  <- input$raw$name
    filt_name <- input$filt$name
    
    if (!is_valid_10x(raw_name,  "raw") ||
        !is_valid_10x(filt_name, "filtered")) {
      return(hard_reset(session,
                        "Upload error: files must be named raw_feature_bc_matrix.(h5|zip) and filtered_feature_bc_matrix.(h5|zip)."))
    }
    
    if (!same_ext(raw_name, filt_name)) {
      return(hard_reset(session,
                        "Upload error: file formats must match (both .h5 or both .zip)."))
    }
    
    rv$upload_type <- if (tools::file_ext(filt_name) == "h5") "h5" else "mtx"
    
    log_msg("Loading RAW matrix from", raw_name)
    log_msg("Loading FILTERED matrix from", filt_name)
    ok <- TRUE
    try({
      rv$raw  <- read_10x_any(input$raw$datapath)
      rv$filt <- read_10x_any(input$filt$datapath)
    }, silent = TRUE)
    if (is.null(rv$raw) || is.null(rv$filt)) ok <- FALSE
    
    if (!ok) {
      return(hard_reset(session, "Upload error: failed to read one or both matrices."))
    }
    
    log_msg("RAW loaded:",  nrow(rv$raw),  "genes x", ncol(rv$raw),  "cells")
    log_msg("FILTERED loaded:", nrow(rv$filt), "genes x", ncol(rv$filt), "cells")
  }, ignoreInit = TRUE)
  
  # -----------------------------
  # Main run
  # -----------------------------
  observeEvent(input$run, {
    if (identical(input$data_mode, "upload") &&
        (is.null(input$raw$datapath) || is.null(input$filt$datapath))) {
      showModal(modalDialog(
        title = "Upload required",
        "Please upload BOTH Raw and Filtered files from 10x matrices before running.",
        easyClose = TRUE, footer = modalButton("OK")
      ))
      return()
    }
    
    req(rv$raw, rv$filt)
    common_genes <- intersect(rownames(rv$raw), rownames(rv$filt))
    if (length(common_genes) == 0) {
      showNotification("No overlapping genes between raw and filtered matrices.", type="error")
      return()
    }
    raw  <- rv$raw [common_genes, , drop=FALSE]
    filt <- rv$filt[common_genes, , drop=FALSE]
    
    if (isTRUE(input$use_gene_symbols)) {
      rn <- rownames(raw)
      looks_ens <- grepl("^(ENSG|ENSMUSG)", rn)
      if (any(looks_ens)) {
        log_msg("Converting Ensembl IDs to symbols via biomaRt for", input$species_genome)
        conv <- tryCatch(map_ensembl_to_symbol(rn, input$species_genome), error=function(e) rn)
        rownames(raw)  <- conv
        rownames(filt) <- conv
      } else if (any(grepl("_", rn))) {
        conv <- sub("^.*_", "", rn)
        rownames(raw)  <- conv
        rownames(filt) <- conv
      }
    }
    
    if (isTRUE(input$collapse_dups)) {
      if (any(duplicated(rownames(raw))) || any(duplicated(rownames(filt)))) {
        log_msg("Collapsing duplicated gene names by summing counts...")
        raw  <- collapse_duplicated_rows(raw)
        filt <- collapse_duplicated_rows(filt)
      }
    } else if (any(duplicated(rownames(raw))) || any(duplicated(rownames(filt)))) {
      log_msg("Forcing unique gene names (suffixing .1, .2...). Consider enabling 'Collapse duplicate gene names'.")
      rownames(raw)  <- make.unique(rownames(raw))
      rownames(filt) <- make.unique(rownames(filt))
    }
    
    log_msg("Creating Seurat object (pre-clean) and computing clusters/UMAP/TSNE...")
    seu <- Seurat::CreateSeuratObject(counts = filt, min.cells = input$min_cells %||% 3)
    defs <- species_defaults(input$species_genome)
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = defs$mito_regex)
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu)
    seu <- ScaleData(seu, verbose=FALSE)
    seu <- RunPCA(seu, npcs=30, verbose=FALSE)
    seu <- FindNeighbors(seu, dims=1:20)
    seu <- FindClusters(seu, resolution=0.5)
    ## identities = seurat_clusters (not orig.ident)
    if (!is.null(seu$seurat_clusters)) {
      Seurat::Idents(seu) <- seu$seurat_clusters
    }
    suppressWarnings({
      if (ncol(seu) >= 50)
        seu <- tryCatch(RunUMAP(seu, dims=1:20), error=function(e) seu)
    })
    suppressWarnings({
      if (ncol(seu) >= 50)
        seu <- tryCatch(RunTSNE(seu, dims = 1:20, check_duplicates = FALSE), error = function(e) seu)
    })
    rv$seu_pre <- seu
    
    log_msg("Decontamination method:", input$method)
    
    # ---- SoupX ----
    if (identical(input$method, "SoupX")) {
      log_msg("Building SoupChannel and setting clusters...")
      sc <- makeSoupX(raw, filt)
      try({ sc <- setClusters(sc, setNames(seu$seurat_clusters, colnames(seu))) }, silent = TRUE)
      
      user_nonexp <- unique(trimws(unlist(strsplit(input$nonexp, ","))))
      user_nonexp <- user_nonexp[user_nonexp != ""]
      nonexp <- unique(c(user_nonexp, species_defaults(input$species_genome)$default_nonexp))
      genes_present <- intersect(nonexp, rownames(filt))
      if (length(genes_present) > 0) {
        log_msg("Marking non-expressing genes:", paste(head(genes_present,5), collapse=", "),
                ifelse(length(genes_present)>5,"…",""))
        try({ sc <- setNonExpressingGenes(sc, genes_present) }, silent = TRUE)
      } else {
        log_msg("No provided non-expressing genes found in data; proceeding without them.")
      }
      
      sr <- c(as.numeric(input$soupRange[1]), as.numeric(input$soupRange[2]))
      if (!all(is.finite(sr))) sr <- c(1, 100)
      sr[1] <- max(1, floor(sr[1]))
      sr[2] <- max(sr[1] + 1, ceiling(sr[2]))
      
      tod_ok <- !is.null(sc$tod) && !is.null(dim(sc$tod)) && all(dim(sc$tod) > 0)
      toc_ok <- !is.null(sc$toc) && !is.null(dim(sc$toc)) && all(dim(sc$toc) > 0)
      if (!tod_ok || !toc_ok) {
        log_msg("WARNING: tod/toc not available or empty; using global soup profile fallback.")
        prof_df <- soup_profile_df(sc, raw_mat = raw, filt_mat = filt)
        sc      <- SoupX::setSoupProfile(sc, prof_df)
      } else {
        drop_umis <- safe_colsums(sc$tod)
        idx <- which(drop_umis >= sr[1] & drop_umis <= sr[2])
        if (length(idx) < 50) {
          log_msg(sprintf("Only %d droplets in soupRange (%d,%d) — using global soup profile.", length(idx), sr[1], sr[2]))
          prof_df <- soup_profile_df(sc, raw_mat = raw, filt_mat = filt)
          sc      <- SoupX::setSoupProfile(sc, prof_df)
        } else {
          log_msg(sprintf("Estimating soup profile: soupRange=(%d,%d), keepDroplets=%s",
                          sr[1], sr[2], as.character(input$keepDroplets)))
          sc <- tryCatch({
            SoupX::estimateSoup(sc, soupRange = sr, keepDroplets = isTRUE(input$keepDroplets))
          }, error = function(e) {
            log_msg("estimateSoup failed:", e$message, "— using global soup profile fallback.")
            prof_df <- soup_profile_df(sc, raw_mat = raw, filt_mat = filt)
            SoupX::setSoupProfile(sc, prof_df)
          })
        }
      }
      
      md <- tryCatch(sc$metaData, error = function(...) NULL)
      if (is.null(md) || !nrow(md)) md <- data.frame(row.names = colnames(sc$toc %||% filt))
      md$nUMIs  <- safe_colsums(sc$toc)
      md$nDrops <- safe_colsums(sc$tod, logical_gt0 = TRUE)
      if (is.null(md$nUMIs))  md$nUMIs  <- rep(0, length(rownames(md)))
      if (is.null(md$nDrops)) md$nDrops <- rep(0, length(rownames(md)))
      md$nUMIs [!is.finite(md$nUMIs )] <- 0
      md$nDrops[!is.finite(md$nDrops)] <- 0
      sc$metaData <- md
      
      if (isTRUE(input$do_auto)) {
        log_msg("Estimating contamination via autoEstCont...")
        auto <- tryCatch({ SoupX::autoEstCont(sc, doPlot = input$do_plot) }, error=function(e) e)
        if (inherits(auto, "error")) {
          showNotification(paste("autoEstCont failed:", auto$message), type="error")
          log_msg("autoEstCont failed:", auto$message, "— falling back to ρ = 0.05 uniform.")
          md <- sc$metaData
          md$rho <- rep(0.05, nrow(md))
          sc$metaData <- md
          rv$auto <- NULL
          rv$sc   <- sc
        } else {
          rv$auto <- auto
          rv$sc   <- auto
        }
        rv$rho <- tryCatch(rv$sc$metaData$rho, error=function(...) NULL)
      } else {
        log_msg(sprintf("Manual SoupX mode: using uniform rho = %.3f", input$manual_rho %||% 0.05))
        md <- sc$metaData
        md$rho <- rep(as.numeric(input$manual_rho %||% 0.05), nrow(md))
        sc$metaData <- md
        rv$sc <- sc
        rv$auto <- NULL
        rv$rho <- md$rho
      }
      
      log_msg("Adjusting counts with SoupX::adjustCounts (defaults; integers rounded)...")
      log_msg("Completed SoupX")
      adj <- tryCatch({ SoupX::adjustCounts(rv$sc, roundToInt = FALSE) }, error=function(e) e)
      if (inherits(adj, "error")) {
        showNotification(paste("adjustCounts failed:", adj$message), type="error")
        log_msg("adjustCounts failed:", adj$message)
        return()
      }
      rv$adj_counts <- round_sparse(adj)
      
      # ---- DecontX ----
    } else if (identical(input$method, "DecontX")) {
      if (!requireNamespace("celda", quietly = TRUE) ||
          !requireNamespace("SingleCellExperiment", quietly = TRUE) ||
          !requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        showNotification("Packages 'celda', 'SingleCellExperiment', and 'SummarizedExperiment' are required for DecontX.", type = "error")
        return()
      }
      log_msg("Preparing SingleCellExperiment for decontX...")
      sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = as(filt, "dgCMatrix"))
      )
      z <- NULL
      if (isTRUE(input$decontx_use_clusters)) {
        cl <- as.character(seu$seurat_clusters)
        cl <- cl[match(colnames(sce), names(cl))]
        if (all(!is.na(cl))) {
          log_msg("Using Seurat clusters as priors for decontX.")
          z <- as.integer(as.factor(cl))
        } else {
          log_msg("Cluster prior could not be matched to all cells; running decontX without priors.")
        }
      }
      delta_vec <- c(10,10)
      if (!is.null(input$decontx_delta) && nzchar(input$decontx_delta)) {
        toks <- as.numeric(unlist(strsplit(input$decontx_delta, "[, ]+")))
        if (length(toks) >= 2 && all(is.finite(toks[1:2]))) delta_vec <- toks[1:2]
      }
      log_msg(sprintf("Running decontX (maxIter=%d, delta=(%.3g,%.3g), estimateDelta=%s, convergence=%.4g, iterLogLik=%d, varGenes=%d)...",
                      as.integer(input$decontx_maxiter %||% 500),
                      delta_vec[1], delta_vec[2],
                      as.character(isTRUE(input$decontx_estimateDelta)),
                      as.numeric(input$decontx_convergence %||% 0.001),
                      as.integer(input$decontx_iterLogLik %||% 10),
                      as.integer(input$decontx_varGenes %||% 5000)))
      log_msg("Completed DecontX")
      dx <- tryCatch({
        celda::decontX(
          x = sce, z = z,
          maxIter       = as.integer(input$decontx_maxiter %||% 500),
          delta         = delta_vec,
          estimateDelta = isTRUE(input$decontx_estimateDelta),
          convergence   = as.numeric(input$decontx_convergence %||% 0.001),
          iterLogLik    = as.integer(input$decontx_iterLogLik %||% 10),
          varGenes      = as.integer(input$decontx_varGenes %||% 5000),
          verbose       = TRUE
        )
      }, error = function(e) e)
      if (inherits(dx, "error")) {
        showNotification(paste("decontX failed:", dx$message), type = "error")
        log_msg("decontX failed:", dx$message)
        return()
      }
      rv$decontx_sce <- dx
      adj <- SummarizedExperiment::assay(dx, "decontXcounts")
      if (!inherits(adj, "dgCMatrix")) adj <- as(adj, "dgCMatrix")
      rv$adj_counts <- round_sparse(adj)
      contam <- tryCatch(SummarizedExperiment::colData(dx)$decontX_contamination, error = function(...) NULL)
      nUMIs  <- Matrix::colSums(SummarizedExperiment::assay(sce, "counts"))
      if (!is.null(contam)) {
        rv$rho <- as.numeric(contam)
        rv$sc  <- list(metaData = data.frame(rho = rv$rho, nUMIs = as.numeric(nUMIs),
                                             row.names = colnames(sce)))
        rv$auto <- NULL
      } else {
        rv$rho <- NULL
        rv$sc  <- list(metaData = data.frame(nUMIs = as.numeric(nUMIs),
                                             row.names = colnames(sce)))
        rv$auto <- NULL
      }
      
      # ---- FastCAR ----
    } else if (identical(input$method, "FastCAR")) {
      if (!fastcar_available()) {
        showNotification("Package 'FastCAR' (and 'Matrix') is required for FastCAR decontamination.", type = "error")
        log_msg("FastCAR aborted: package not available.")
        return()
      }
      
      empty_cutoff  <- as.integer(input$fastcar_empty_cutoff %||% 100L)
      contam_cutoff <- as.numeric(input$fastcar_contam_cutoff %||% 0.05)
      do_profile    <- isTRUE(input$fastcar_do_profile)
      prof_start    <- as.integer(input$fastcar_profile_start %||% 10L)
      prof_stop     <- as.integer(input$fastcar_profile_stop %||% 500L)
      prof_by       <- as.integer(input$fastcar_profile_by %||% 10L)
      use_reco      <- isTRUE(input$fastcar_use_recommended)
      
      log_msg(sprintf(
        paste(
          "FastCAR parameters:",
          "empty_droplet_cutoff=%d, contamination_chance_cutoff=%.3f,",
          "do_profile=%s, profile_start=%d, profile_stop=%d,",
          "profile_by=%d, use_recommended_cutoff=%s"
        ),
        empty_cutoff, contam_cutoff, do_profile,
        prof_start, prof_stop, prof_by, use_reco
      ))
      
      log_msg("Completed FastCAR")
      
      fc <- tryCatch(
        scarr_fastcar_wrapper(
          raw,
          filt,
          empty_droplet_cutoff        = empty_cutoff,
          contamination_chance_cutoff = contam_cutoff,
          do_profile                  = do_profile,
          profile_start               = prof_start,
          profile_stop                = prof_stop,
          profile_by                  = prof_by,
          use_recommended_cutoff      = use_reco
        ),
        error = function(e) e
      )
      
      if (inherits(fc, "error")) {
        showNotification(paste("FastCAR failed:", fc$message), type = "error")
        log_msg("FastCAR failed:", fc$message)
        return()
      }
      
      adj <- fc$corrected_counts
      if (!inherits(adj, "dgCMatrix")) adj <- as(adj, "dgCMatrix")
      rv$adj_counts <- round_sparse(adj)
      
      nUMIs <- tryCatch(Matrix::colSums(filt), error = function(...) NULL)
      if (is.null(nUMIs)) nUMIs <- rep(NA_real_, ncol(filt))
      md <- data.frame(
        nUMIs = as.numeric(nUMIs),
        row.names = colnames(filt),
        stringsAsFactors = FALSE
      )
      
      rv$rho  <- NULL
      rv$auto <- NULL
      rv$sc   <- list(
        metaData = md,
        params   = c(
          empty_droplet_cutoff        = empty_cutoff,
          contamination_chance_cutoff = contam_cutoff,
          do_profile                  = do_profile,
          profile_start               = prof_start,
          profile_stop                = prof_stop,
          profile_by                  = prof_by,
          use_recommended_cutoff      = use_reco
        ),
        fastcar_aux = fc$aux
      )
      
      # ---- scCDC ----
    } else if (identical(input$method, "scCDC")) {
      if (!requireNamespace("scCDC", quietly = TRUE)) {
        showNotification("Package 'scCDC' is required for this method.", type = "error")
        log_msg("scCDC aborted: package not available.")
        return()
      }
      seu <- rv$seu_pre
      if (is.null(seu) || is.null(seu$seurat_clusters)) {
        showNotification("scCDC needs clustered Seurat object; retry after QC (Pre) clustering.", type = "error")
        log_msg("scCDC aborted: seu_pre missing or clusters missing.")
        return()
      }
      restriction_factor <- as.numeric(input$sccdc_restriction %||% 0.5)
      min_cell           <- as.integer(input$sccdc_min_cell %||% 100)
      percent_cutoff     <- as.numeric(input$sccdc_percent_cutoff %||% 0.2)
      
      out_dir <- normalizePath(reports_dir, winslash = "/", mustWork = TRUE)
      if (!grepl("/$", out_dir)) out_dir <- paste0(out_dir, "/")
      
      clear_scCDC_pdfs(reports_dir)
      
      log_msg(sprintf(
        "Running scCDC (restriction_factor=%.2f, min.cell=%d, percent.cutoff=%.2f) with out_path.plot=%s",
        restriction_factor, min_cell, percent_cutoff, out_dir
      ))
      log_msg("Completed scCDC")
      
      GCGs <- tryCatch(
        scCDC::ContaminationDetection(
          seu,
          restriction_factor = restriction_factor,
          min.cell           = min_cell,
          percent.cutoff     = percent_cutoff,
          out_path.plot      = out_dir
        ),
        error = function(e) e
      )
      if (inherits(GCGs, "error")) {
        showNotification(paste("ContaminationDetection failed:", GCGs$message), type = "error")
        log_msg("scCDC ContaminationDetection failed:", GCGs$message)
        return()
      }
      rv$GCGs <- GCGs
      
      mis_pat <- "^scCDC_reportsdefault_SE-plot.*[.]pdf$"
      mis_pdf <- list.files(tempdir(), pattern = mis_pat, full.names = TRUE)
      if (length(mis_pdf)) {
        dst <- file.path(reports_dir, sub("^scCDC_reports", "", basename(mis_pdf[1])))
        file.copy(mis_pdf[1], dst, overwrite = TRUE)
        rv$scCDC_pdf_abs <- dst
        rv$scCDC_pdf_url <- paste0("reports/", basename(dst))
        log_msg("Moved misplaced scCDC PDF to:", dst)
      } else {
        pdfs <- list.files(reports_dir, pattern = "[.]pdf$", full.names = TRUE)
        if (length(pdfs)) {
          newest <- pdfs[which.max(file.info(pdfs)$mtime)]
          rv$scCDC_pdf_abs <- newest
          rv$scCDC_pdf_url <- paste0("reports/", basename(newest))
          log_msg("scCDC PDF report:", newest)
        } else {
          rv$scCDC_pdf_abs <- NULL
          rv$scCDC_pdf_url <- NULL
          showNotification("scCDC ran, but no PDF was found. Check permissions or filenames.", type = "warning")
          log_msg("No scCDC PDF found under:", reports_dir)
        }
      }
      
      gcg_genes <- tryCatch(rownames(GCGs), error=function(...) NULL)
      if (is.null(gcg_genes)) gcg_genes <- tryCatch(as.character(GCGs$gene), error=function(...) NULL)
      if (is.null(gcg_genes)) {
        showNotification("scCDC: Could not extract GCG gene list.", type = "error")
        return()
      }
      gcg_genes <- unique(gcg_genes[gcg_genes %in% rownames(seu)])
      rv$sccdc_quant <- tryCatch(scCDC::ContaminationQuantification(seu, gcg_genes), error = function(e) NULL)
      
      seu_corr <- tryCatch(scCDC::ContaminationCorrection(seu, gcg_genes), error = function(e) e)
      if (inherits(seu_corr, "error")) {
        showNotification(paste("ContaminationCorrection failed:", seu_corr$message), type = "error")
        log_msg("scCDC ContaminationCorrection failed:", seu_corr$message)
        return()
      }
      corrected <- NULL
      if ("Corrected" %in% names(seu_corr@assays)) {
        DefaultAssay(seu_corr) <- "Corrected"
        corrected <- tryCatch(
          Seurat::GetAssayData(seu_corr, layer = "counts", assay = "Corrected"),
          error=function(e) NULL
        )
      }
      if (is.null(corrected)) {
        showNotification("scCDC: Could not retrieve corrected counts from 'Corrected' assay.", type = "error")
        log_msg("scCDC: 'Corrected' assay not found or no counts slot.")
        return()
      }
      if (!inherits(corrected, "dgCMatrix")) corrected <- as(corrected, "dgCMatrix")
      rv$adj_counts <- round_sparse(corrected)
      
      rv$rho <- NULL
      rv$sc  <- list(metaData = data.frame(
        nUMIs = Matrix::colSums(rv$adj_counts),
        row.names = colnames(rv$adj_counts)
      ))
      rv$auto <- NULL
    }
    
    do_adjust_and_update()
    updateTabsetPanel(session, "tabs", selected = "QC (Pre)")
  })
  
  # -----------------------------
  # Post-adjustment
  # -----------------------------
  do_adjust_and_update <- function() {
    req(rv$seu_pre, rv$adj_counts)
    
    pre_umis  <- tryCatch(sum(rv$seu_pre@assays[[DefaultAssay(rv$seu_pre)]]@counts),  error=function(...) NA_real_)
    post_umis <- tryCatch(sum(rv$adj_counts), error=function(...) NA_real_)
    if (is.finite(pre_umis) && is.finite(post_umis)) {
      log_msg(sprintf("UMIs pre: %.0f; post: %.0f; removed: %.0f (%.2f%%)",
                      pre_umis, post_umis, pre_umis - post_umis,
                      100 * (pre_umis - post_umis) / pre_umis))
    }
    
    defs <- species_defaults(input$species_genome)
    seu2 <- Seurat::CreateSeuratObject(counts = rv$adj_counts, min.cells = input$min_cells %||% 3)
    if (ncol(seu2) == 0) {
      showNotification("No cells remained after decontamination. Try less aggressive settings.", type = "error")
      return()
    }
    seu2[["percent.mt"]] <- PercentageFeatureSet(seu2, pattern = defs$mito_regex)
    seu2 <- NormalizeData(seu2)
    seu2 <- FindVariableFeatures(seu2)
    seu2 <- ScaleData(seu2, verbose=FALSE)
    seu2 <- RunPCA(seu2, npcs=30, verbose=FALSE)
    seu2 <- FindNeighbors(seu2, dims=1:20)
    seu2 <- FindClusters(seu2, resolution=0.5)
    ## identities = seurat_clusters here as well
    if (!is.null(seu2$seurat_clusters)) {
      Seurat::Idents(seu2) <- seu2$seurat_clusters
    }
    suppressWarnings({
      if (ncol(seu2) >= 50)
        seu2 <- tryCatch(RunUMAP(seu2, dims=1:20), error=function(e) seu2)
    })
    suppressWarnings({
      if (ncol(seu2) >= 50)
        seu2 <- tryCatch(RunTSNE(seu2, dims = 1:20, check_duplicates = FALSE), error = function(e) seu2)
    })
    rv$seu_post <- seu2
    
    pre_counts  <- tapply(rv$seu_pre$seurat_clusters,  rv$seu_pre$seurat_clusters,  length)
    post_counts <- tapply(rv$seu_post$seurat_clusters, rv$seu_post$seurat_clusters, length)
    # Ensure numeric cluster ordering: 0, 1, 2, ... 10, 11 (not 0, 1, 10, 11, ...)
    all_clusters <- sort_clusters(c(names(pre_counts), names(post_counts)))
    cluster_counts_df <- data.frame(
      cluster = all_clusters,
      n_pre   = as.integer(pre_counts [all_clusters]),
      n_post  = as.integer(post_counts[all_clusters]),
      stringsAsFactors = FALSE
    )
    cluster_counts_df$n_pre [is.na(cluster_counts_df$n_pre )] <- 0L
    cluster_counts_df$n_post[is.na(cluster_counts_df$n_post)] <- 0L
    cluster_counts_df$total_change <- cluster_counts_df$n_post - cluster_counts_df$n_pre
    rv$cluster_counts_df <- cluster_counts_df
    
    pre_df <- data.frame(
      cell            = colnames(rv$seu_pre),
      cluster_pre     = as.character(rv$seu_pre$seurat_clusters),
      nUMIs_pre       = rv$seu_pre@meta.data$nCount_RNA,
      nFeature_pre    = rv$seu_pre@meta.data$nFeature_RNA,
      percent.mt_pre  = rv$seu_pre@meta.data$percent.mt,
      stringsAsFactors = FALSE
    )
    post_df <- data.frame(
      cell            = colnames(rv$seu_post),
      cluster_post    = as.character(rv$seu_post$seurat_clusters),
      nUMIs_post      = rv$seu_post@meta.data$nCount_RNA,
      nFeature_post   = rv$seu_post@meta.data$nFeature_RNA,
      percent.mt_post = rv$seu_post@meta.data$percent.mt,
      stringsAsFactors = FALSE
    )
    cells_df <- merge(pre_df, post_df, by = "cell", all = TRUE)
    md <- tryCatch(rv$sc$metaData, error = function(e) NULL)
    if (!is.null(md) && "rho" %in% colnames(md)) {
      rho_named <- setNames(md$rho, rownames(md))
      cells_df$rho <- rho_named[match(cells_df$cell, names(rho_named))]
    } else {
      cells_df$rho <- NA_real_
    }
    cells_df$delta_nUMIs    <- with(cells_df, nUMIs_post    - nUMIs_pre)
    cells_df$delta_features <- with(cells_df, nFeature_post - nFeature_pre)
    cells_df <- cells_df[
      order(-cells_df$nUMIs_post %||% 0, -cells_df$nUMIs_pre %||% 0),
    ]
    rv$cells_df <- cells_df
  }
  
  ##########################
  ## Plots
  ##########################
  
  plot_pre_nFeature <- reactive({
    req(rv$seu_pre)
    ggplot(data.frame(nFeature_RNA = rv$seu_pre@meta.data$nFeature_RNA),
           aes(nFeature_RNA)) + geom_histogram(bins = 50) +
      ggtitle("QC (Pre): Detected genes per cell") +
      xlab("Detected genes per cell (pre)") + theme_clean
  })
  plot_pre_nCount <- reactive({
    req(rv$seu_pre)
    ggplot(data.frame(nCount_RNA = rv$seu_pre@meta.data$nCount_RNA),
           aes(nCount_RNA)) + geom_histogram(bins = 50) +
      ggtitle("QC (Pre): UMIs per cell") +
      xlab("UMIs per cell (pre)") + theme_clean
  })
  plot_pre_pctMT <- reactive({
    req(rv$seu_pre)
    ggplot(data.frame(percent.mt = rv$seu_pre@meta.data$percent.mt),
           aes(percent.mt)) + geom_histogram(bins = 50) +
      ggtitle("QC (Pre): Mitochondrial %") +
      xlab("Mitochondrial % (pre)") + theme_clean
  })
  
  plot_umap_pre <- reactive({
    req(rv$seu_pre)
    if (!"umap" %in% names(rv$seu_pre@reductions)) return(NULL)
    DimPlot(rv$seu_pre, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
      ggtitle("UMAP (Pre)")
  })
  plot_umap_post <- reactive({
    req(rv$seu_post)
    if (!"umap" %in% names(rv$seu_post@reductions)) return(NULL)
    DimPlot(rv$seu_post, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
      ggtitle("UMAP (Post)")
  })
  
  plot_tsne_pre <- reactive({
    req(rv$seu_pre)
    if (!"tsne" %in% names(rv$seu_pre@reductions)) return(NULL)
    DimPlot(rv$seu_pre, reduction = "tsne", group.by = "seurat_clusters", label = TRUE) +
      ggtitle("TSNE (Pre)")
  })
  plot_tsne_post <- reactive({
    req(rv$seu_post)
    if (!"tsne" %in% names(rv$seu_post@reductions)) return(NULL)
    DimPlot(rv$seu_post, reduction = "tsne", group.by = "seurat_clusters", label = TRUE) +
      ggtitle("TSNE (Post)")
  })
  
  plot_est_density <- reactive({
    r <- rv$rho
    if (is.null(r)) return(NULL)
    r <- r[is.finite(r)]
    if (!length(r)) return(NULL)
    
    d <- density(pmin(pmax(r, 0), 1), from = 0, to = 1, n = 512, na.rm = TRUE)
    df <- data.frame(x = d$x, y = d$y)
    m  <- mean(r, na.rm = TRUE)
    q  <- stats::quantile(r, c(0.1, 0.9), na.rm = TRUE)
    
    ggplot(df, aes(x, y)) + geom_line() + geom_vline(xintercept = m, color = "red") +
      ggtitle("Estimation: rho density") +
      labs(x = "Contamination Fraction", y = "Probability Density",
           subtitle = sprintf("ρ mean %.3f (10th=%.3f, 90th=%.3f)  |  red = mean", m, q[[1]], q[[2]])) +
      theme_clean
  })
  
  plot_rho_vs_counts <- reactive({
    req(rv$sc)
    md <- rv$sc$metaData
    if (!all(c("rho","nUMIs") %in% colnames(md))) return(NULL)
    ggplot(md, aes(nUMIs, rho)) + geom_point(alpha = 0.4, size = 0.6) +
      ggtitle("Estimation: rho vs nUMIs") +
      scale_x_log10() + xlab("nUMIs (log10)") + ylab("rho") + theme_clean
  })
  
  plot_auto_plot <- reactive({
    req(rv$auto, input$do_plot)
    md <- rv$auto$metaData
    md$q <- cut(
      md$nUMIs,
      breaks = quantile(md$nUMIs, probs = seq(0,1,0.1), na.rm = TRUE),
      include.lowest = TRUE
    )
    df <- aggregate(rho ~ q, data = md, FUN = function(x) mean(x, na.rm = TRUE))
    ggplot(df, aes(x = q, y = rho, group = 1)) + geom_line() + geom_point() +
      ggtitle("Estimation: autoEstCont diagnostic") +
      xlab("UMI decile") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_clean
  })
  
  output$pre_nFeature <- renderPlot({ req(plot_pre_nFeature()); plot_pre_nFeature() })
  output$pre_nCount   <- renderPlot({ req(plot_pre_nCount());   plot_pre_nCount()   })
  output$pre_pctMT    <- renderPlot({ req(plot_pre_pctMT());    plot_pre_pctMT()    })
  
  output$umap_pre     <- renderPlot({ p <- plot_umap_pre();  if (!is.null(p)) p })
  output$umap_post    <- renderPlot({ p <- plot_umap_post(); if (!is.null(p)) p })
  
  output$tsne_pre     <- renderPlot({ p <- plot_tsne_pre();  if (!is.null(p)) p })
  output$tsne_post    <- renderPlot({ p <- plot_tsne_post(); if (!is.null(p)) p })
  
  output$est_density  <- renderPlot({
    if (identical(input$method, "scCDC")) return(NULL)
    p <- plot_est_density(); if (!is.null(p)) p
  })
  output$rho_vs_counts<- renderPlot({
    if (identical(input$method, "scCDC")) return(NULL)
    p <- plot_rho_vs_counts(); if (!is.null(p)) p
  })
  output$auto_plot    <- renderPlot({
    if (!identical(input$method, "SoupX")) return(NULL)
    p <- plot_auto_plot(); if (!is.null(p)) p
  })
  
  plot_post_nFeature <- reactive({
    if (is.null(rv$seu_post)) return(NULL)
    ggplot(data.frame(nFeature_RNA = rv$seu_post@meta.data$nFeature_RNA),
           aes(nFeature_RNA)) + geom_histogram(bins = 50) +
      ggtitle("QC (Post): Detected genes per cell") +
      xlab("Detected genes per cell (post)") + theme_clean
  })
  plot_post_nCount <- reactive({
    if (is.null(rv$seu_post)) return(NULL)
    ggplot(data.frame(nCount_RNA = rv$seu_post@meta.data$nCount_RNA),
           aes(nCount_RNA)) + geom_histogram(bins = 50) +
      ggtitle("QC (Post): UMIs per cell") +
      xlab("UMIs per cell (post)") + theme_clean
  })
  plot_post_pctMT <- reactive({
    if (is.null(rv$seu_post)) return(NULL)
    ggplot(data.frame(percent.mt = rv$seu_post@meta.data$percent.mt),
           aes(percent.mt)) + geom_histogram(bins = 50) +
      ggtitle("QC (Post): Mitochondrial %") +
      xlab("Mitochondrial % (post)") + theme_clean
  })
  
  output$post_nFeature <- renderPlot({ p <- plot_post_nFeature(); if (!is.null(p)) p })
  output$post_nCount   <- renderPlot({ p <- plot_post_nCount();   if (!is.null(p)) p })
  output$post_pctMT    <- renderPlot({ p <- plot_post_pctMT();    if (!is.null(p)) p })
  
  plot_cluster_counts_bar <- reactive({
    if (is.null(rv$cluster_counts_df) || !nrow(rv$cluster_counts_df)) return(NULL)
    df <- rv$cluster_counts_df[, c("cluster", "n_pre", "n_post")]
    df$cluster <- as.character(df$cluster)
    df$cluster <- factor(df$cluster, levels = sort_clusters(df$cluster))
    long <- reshape2::melt(df, id.vars = "cluster", variable.name = "set", value.name = "n")
    ggplot(long, aes(x = cluster, y = n, fill = set)) +
      geom_col(position = position_dodge(width = 0.8)) +
      ggtitle("Cluster counts: Pre vs Post") +
      labs(x = "Cluster", y = "Cell count", fill = NULL) + theme_clean
  })
  output$cluster_counts_bar <- renderPlot({
    p <- plot_cluster_counts_bar(); if (!is.null(p)) p
  })
  
  # DecontX contamination plot + panel
  plot_decontx_contam <- reactive({
    if (!identical(input$method, "DecontX")) return(NULL)
    if (is.null(rv$decontx_sce)) return(NULL)
    tryCatch(celda::plotDecontXContamination(rv$decontx_sce), error = function(e) NULL)
  })
  output$decontx_contam_plot <- renderPlot({
    p <- plot_decontx_contam(); if (!is.null(p)) p
  })
  output$decontx_contam_panel <- renderUI({
    if (!identical(input$method, "DecontX")) return(NULL)
    tagList(
      h3("DecontX UMAP contamination plot"),
      plotOutput("decontx_contam_plot") %>% withSpinner(type = 4),
      div(class="dl-row", actionBttn("open_dl_decontx_contam", "Download plot",
                                     style = "unite", color = "primary", icon = icon("download")))
    )
  })
  observeEvent(input$open_dl_decontx_contam, {
    showModal(modalDialog(
      title = "Download: DecontX contamination",
      numericInput("decontx_contam_h", "Figure height", 6),
      numericInput("decontx_contam_w", "Figure width",  8),
      numericInput("decontx_contam_dpi", "Figure resolution", 300),
      selectInput("decontx_contam_fmt", "Image format",
                  choices = c(".jpg", ".tiff", ".pdf", ".svg", ".bmp", ".eps", ".ps")),
      downloadBttn("download_decontx_contam", "Download"),
      easyClose = TRUE
    ))
  })
  output$download_decontx_contam <- downloadHandler(
    filename = function() {
      ext <- input$decontx_contam_fmt %||% ".jpg"
      paste0("decontx_contamination", ext)
    },
    content = function(file) {
      p <- plot_decontx_contam()
      if (is.null(p)) stop("No DecontX contamination plot available")
      h <- as.numeric(input$decontx_contam_h %||% 6)
      w <- as.numeric(input$decontx_contam_w %||% 8)
      dpi <- as.numeric(input$decontx_contam_dpi %||% 300)
      ggplot2::ggsave(filename = file, plot = p,
                      width = w, height = h, dpi = dpi, limitsize = FALSE)
    }
  )
  
  # scCDC diagnostics
  plot_scCDC_delta <- reactive({
    if (!identical(input$method, "scCDC")) return(NULL)
    if (is.null(rv$adj_counts) || is.null(rv$seu_pre)) return(NULL)
    pre  <- rv$seu_pre@meta.data$nCount_RNA
    post <- Matrix::colSums(rv$adj_counts)
    if (length(post) != length(pre)) return(NULL)
    df <- data.frame(nUMIs_pre = pre, nUMIs_post = as.numeric(post))
    df$delta <- df$nUMIs_post - df$nUMIs_pre
    ggplot(df, aes(nUMIs_pre, nUMIs_post)) +
      geom_point(alpha = 0.4, size = 0.6) +
      scale_x_log10() + scale_y_log10() +
      ggtitle("scCDC: UMIs per cell (Post vs Pre)") +
      xlab("Pre (log10)") + ylab("Post (log10)") + theme_clean
  })
  output$scCDC_diag_plot <- renderPlot({
    if (!identical(input$method, "scCDC")) return(NULL)
    p <- plot_scCDC_delta(); if (!is.null(p)) p
  })
  output$scCDC_diag_panel <- renderUI({
    if (!identical(input$method, "scCDC")) return(NULL)
    tagList(
      h3("UMIs per cell before and after Ambient RNA removal"),
      plotOutput("scCDC_diag_plot") %>% withSpinner(type = 4),
      div(class="dl-row", actionBttn("open_dl_scCDC_diag", "Download plot",
                                     style = "unite", color = "primary", icon = icon("download")))
    )
  })
  observeEvent(input$open_dl_scCDC_diag, {
    showModal(modalDialog(
      title = "Download: scCDC diagnostic (UMIs Post vs Pre)",
      numericInput("scCDC_diag_h", "Figure height", 6),
      numericInput("scCDC_diag_w", "Figure width",  8),
      numericInput("scCDC_diag_dpi", "Figure resolution", 300),
      selectInput("scCDC_diag_fmt", "Image format",
                  choices = c(".jpg", ".tiff", ".pdf", ".svg", ".bmp", ".eps", ".ps")),
      downloadBttn("download_scCDC_diag", "Download"),
      easyClose = TRUE
    ))
  })
  output$download_scCDC_diag <- downloadHandler(
    filename = function() {
      ext <- input$scCDC_diag_fmt %||% ".jpg"
      paste0("scCDC_diag", ext)
    },
    content = function(file) {
      p <- plot_scCDC_delta()
      if (is.null(p)) stop("No scCDC diagnostic plot available")
      h <- as.numeric(input$scCDC_diag_h %||% 6)
      w <- as.numeric(input$scCDC_diag_w %||% 8)
      dpi <- as.numeric(input$scCDC_diag_dpi %||% 300)
      ggplot2::ggsave(filename = file, plot = p,
                      width = w, height = h, dpi = dpi, limitsize = FALSE)
    }
  )
  
  # scCDC extra panel: Top GCGs
  output$sccdc_panel <- renderUI({
    if (!identical(input$method, "scCDC")) return(NULL)
    tagList(
      h3("Mean count before and after Ambient RNA removal"),
      plotOutput("sccdc_top_gcgs_prepost", height = 350) %>% withSpinner(type = 4),
      div(class="dl-row", actionBttn("open_dl_sccdc_top_gcgs_prepost", "Download plot",
                                     style = "unite", color = "primary", icon = icon("download")))
    )
  })
  
  sccdc_gcgs_genes <- reactive({
    if (!identical(input$method, "scCDC")) return(character(0))
    G <- rv$GCGs; if (is.null(G)) return(character(0))
    genes <- tryCatch(rownames(G), error=function(...) NULL)
    if (is.null(genes)) genes <- tryCatch(as.character(G$gene), error=function(...) NULL)
    genes <- unique(genes)
    genes[genes %in% rownames(rv$seu_pre)]
  })
  
  plot_top_gcgs_prepost <- reactive({
    req(rv$seu_pre, rv$adj_counts)
    genes <- sccdc_gcgs_genes()
    if (!length(genes)) return(NULL)
    pre_counts  <- Seurat::GetAssayData(rv$seu_pre,  layer = "counts", assay = DefaultAssay(rv$seu_pre))
    post_counts <- rv$adj_counts
    genes <- intersect(genes, intersect(rownames(pre_counts), rownames(post_counts)))
    if (!length(genes)) return(NULL)
    pre_m  <- as.numeric(Matrix::rowMeans(pre_counts[genes, , drop=FALSE]))
    post_m <- as.numeric(Matrix::rowMeans(post_counts[genes, , drop=FALSE]))
    df <- data.frame(gene = genes, pre = pre_m, post = post_m)
    df <- df[order(-df$pre), , drop = FALSE]
    df_top <- head(df, 20)
    long <- reshape2::melt(df_top, id.vars = "gene", variable.name = "set", value.name = "mean_expr")
    ggplot(long, aes(x = reorder(gene, -mean_expr), y = mean_expr, fill = set)) +
      geom_col(position = position_dodge(width = 0.8)) +
      theme_clean + xlab("") + ylab("Mean counts") +
      ggtitle("Top GCGs: Pre vs Post (mean counts)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  output$sccdc_top_gcgs_prepost <- renderPlot({
    p <- plot_top_gcgs_prepost(); if (!is.null(p)) p
  })
  
  observeEvent(input$open_dl_sccdc_top_gcgs_prepost, {
    showModal(modalDialog(
      title = "Download: Top GCGs (Pre vs Post mean)",
      numericInput("sccdc_top_h", "Figure height", 6),
      numericInput("sccdc_top_w", "Figure width",  10),
      numericInput("sccdc_top_dpi", "Figure resolution", 300),
      selectInput("sccdc_top_fmt", "Image format",
                  choices = c(".jpg", ".tiff", ".pdf", ".svg", ".bmp", ".eps", ".ps")),
      downloadBttn("download_sccdc_top", "Download"),
      easyClose = TRUE
    ))
  })
  output$download_sccdc_top <- downloadHandler(
    filename = function() {
      ext <- input$sccdc_top_fmt %||% ".jpg"
      paste0("scCDC_top_GCGs_pre_post", ext)
    },
    content = function(file) {
      p <- plot_top_gcgs_prepost()
      if (is.null(p)) stop("No plot available")
      ggplot2::ggsave(filename = file, plot = p,
                      width  = as.numeric(input$sccdc_top_w   %||% 10),
                      height = as.numeric(input$sccdc_top_h   %||% 6),
                      dpi    = as.numeric(input$sccdc_top_dpi %||% 300),
                      limitsize = FALSE)
    }
  )
  
  # scCDC PDF panel
  output$scCDC_pdf_panel <- renderUI({
    if (!identical(input$method, "scCDC")) return(NULL)
    if (is.null(rv$scCDC_pdf_url)) {
      return(tags$div(class = "text-muted",
                      "Run scCDC to generate a PDF report."))
    }
    tagList(
      tags$h3("scCDC PDF report (from ContaminationDetection)"),
      tags$iframe(src = rv$scCDC_pdf_url,
                  style = "width:100%;height:600px;border:1px solid #ddd;"),
      br(),
      downloadBttn("dl_scCDC_pdf", "Download scCDC report (PDF)")
    )
  })
  output$dl_scCDC_pdf <- downloadHandler(
    filename = function() sprintf("scCDC_report_%s.pdf", format(Sys.Date(), "%Y%m%d")),
    content = function(file) {
      req(rv$scCDC_pdf_abs)
      file.copy(rv$scCDC_pdf_abs, file, overwrite = TRUE)
    }
  )
  
  # FastCAR estimation plots
  fastcar_profile_obj <- reactive({
    if (!identical(input$method, "FastCAR")) return(NULL)
    aux <- tryCatch(rv$sc$fastcar_aux, error = function(e) NULL)
    amb <- tryCatch(aux$amb_profile,    error = function(e) NULL)
    if (is.null(amb)) return(NULL)
    fastcar_ambient_profile_plot(amb)
  })
  
  output$fastcar_profile_plot <- renderPlot({
    if (!identical(input$method, "FastCAR")) return(NULL)
    p <- fastcar_profile_obj()
    if (!is.null(p)) p
  })
  
  fastcar_removed_hist_obj <- reactive({
    if (!identical(input$method, "FastCAR")) return(NULL)
    aux <- tryCatch(rv$sc$fastcar_aux, error = function(e) NULL)
    v <- aux$reads_removed_per_cell
    if (is.null(v)) return(NULL)
    df <- data.frame(removed = as.numeric(v))
    ggplot(df, aes(removed)) +
      geom_histogram(bins = 50) +
      ggtitle("FastCAR: reads removed per cell") +
      xlab("UMIs removed") +
      theme_clean
  })
  
  output$fastcar_removed_hist <- renderPlot({
    if (!identical(input$method, "FastCAR")) return(NULL)
    p <- fastcar_removed_hist_obj()
    if (!is.null(p)) p
  })
  
  # Heatmap + download
  heatmap_plot_obj <- reactive({
    req(rv$seu_pre, rv$seu_post)
    src  <- input$heatmap_source %||% "post"
    topn <- as.integer(input$heatmap_topn %||% 50)
    topn <- max(10L, min(200L, topn))
    obj  <- if (identical(src, "pre")) rv$seu_pre else rv$seu_post
    feats <- head(VariableFeatures(obj), topn)
    if (length(feats) < 10) {
      obj   <- FindVariableFeatures(obj)
      feats <- head(VariableFeatures(obj), topn)
    }
    suppressWarnings(DoHeatmap(obj, features = feats, group.by = "seurat_clusters"))
  })
  output$heatmap_plot <- renderPlot({
    req(heatmap_plot_obj())
    heatmap_plot_obj()
  })
  observeEvent(input$open_heatmap_download, {
    showModal(modalDialog(
      title = "Download Heatmap",
      numericInput("heatmap_plot1_height", "Figure height", 8),
      numericInput("heatmap_plot1_width",  "Figure width",  8),
      numericInput("heatmap_plot1_dpi",    "Figure resolution", 300),
      selectInput("heatmap_plot1_type", "Image format",
                  choices = c(".jpg", ".tiff", ".pdf", ".svg", ".bmp", ".eps", ".ps")),
      downloadBttn("downloadoutput_heatmap1", "Download"),
      easyClose = TRUE
    ))
  })
  output$downloadoutput_heatmap1 <- downloadHandler(
    filename = function() {
      ext <- input$heatmap_plot1_type %||% ".jpg"
      paste0("heatmap_", input$heatmap_source %||% "post", ext)
    },
    content = function(file) {
      p <- heatmap_plot_obj()
      if (is.null(p)) stop("No heatmap available")
      h <- as.numeric(input$heatmap_plot1_height %||% 8)
      w <- as.numeric(input$heatmap_plot1_width  %||% 8)
      dpi <- as.numeric(input$heatmap_plot1_dpi  %||% 300)
      ggplot2::ggsave(filename = file, plot = p,
                      width = w, height = h, dpi = dpi, limitsize = FALSE)
    }
  )
  
  # Feature Plot
  output$featplot_warning <- renderUI({
    genes <- clean_gene_list(input$featplot_genes)
    if (length(genes) == 0)
      return(tags$div(class = "text-muted", "Enter one or more gene symbols above."))
    NULL
  })
  featplot_combined_obj <- reactive({
    req(rv$seu_pre, rv$seu_post)
    genes <- clean_gene_list(input$featplot_genes)
    if (!length(genes)) return(NULL)
    present <- genes[genes %in% rownames(rv$seu_pre) & genes %in% rownames(rv$seu_post)]
    if (!length(present)) return(NULL)
    rows <- lapply(present, function(g) {
      p_pre  <- tryCatch(
        FeaturePlot(rv$seu_pre,  features = g, order = TRUE) +
          ggtitle(paste0(g, " — Pre")),
        error=function(e) NULL
      )
      p_post <- tryCatch(
        FeaturePlot(rv$seu_post, features = g, order = TRUE) +
          ggtitle(paste0(g, " — Post")),
        error=function(e) NULL
      )
      if (is.null(p_pre) || is.null(p_post)) return(NULL)
      patchwork::wrap_plots(p_pre, p_post, ncol = 2)
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) return(NULL)
    Reduce(`+`, lapply(rows, function(r) r + patchwork::plot_layout(ncol = 1)))
  })
  output$featplot_combined <- renderPlot({
    p <- featplot_combined_obj(); if (!is.null(p)) p
  })
  observeEvent(input$open_dl_featplot, {
    showModal(modalDialog(
      title = "Download: Feature Plot (Pre vs Post)",
      numericInput("featplot_h", "Figure height", 6),
      numericInput("featplot_w", "Figure width", 12),
      numericInput("featplot_dpi", "Figure resolution", 300),
      selectInput("featplot_fmt", "Image format",
                  choices = c(".jpg", ".tiff", ".pdf", ".svg", ".bmp", ".eps", ".ps")),
      downloadBttn("download_featplot", "Download"),
      easyClose = TRUE
    ))
  })
  output$download_featplot <- downloadHandler(
    filename = function() {
      ext <- input$featplot_fmt %||% ".jpg"
      paste0("featureplot_pre_post", ext)
    },
    content = function(file) {
      p <- featplot_combined_obj()
      if (is.null(p)) stop("No Feature Plot available")
      h <- as.numeric(input$featplot_h  %||% 6)
      w <- as.numeric(input$featplot_w  %||% 12)
      dpi <- as.numeric(input$featplot_dpi %||% 300)
      ggplot2::ggsave(filename = file, plot = p,
                      width = w, height = h, dpi = dpi, limitsize = FALSE)
    }
  )
  
  # Tables & downloads
  output$top_genes <- renderDT({
    req(rv$seu_pre, rv$seu_post)
    pre_counts  <- Seurat::GetAssayData(rv$seu_pre,  layer = "counts", assay = DefaultAssay(rv$seu_pre))
    post_counts <- Seurat::GetAssayData(rv$seu_post, layer = "counts", assay = DefaultAssay(rv$seu_post))
    if (nrow(pre_counts) == 0 || ncol(pre_counts) == 0 ||
        nrow(post_counts) == 0 || ncol(post_counts) == 0) {
      return(datatable(
        data.frame(message = "Counts are empty; cannot compute Top Genes."),
        options = list(pageLength = 5),
        rownames = FALSE
      ))
    }
    pre_sum  <- Matrix::rowSums(pre_counts)
    post_sum <- Matrix::rowSums(post_counts)
    common <- intersect(rownames(pre_counts), rownames(post_counts))
    df <- data.frame(
      gene       = common,
      total_pre  = as.numeric(pre_sum[common]),
      total_post = as.numeric(post_sum[common])
    )
    df$delta <- df$total_post - df$total_pre
    df <- df[order(-df$total_pre), ]
    datatable(df,
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE, selection = "none")
  })
  
  output$cells_dt <- renderDT({
    req(rv$cells_df)
    datatable(rv$cells_df,
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE, selection = "none")
  })

  # Total cell count (Pre vs Post)
  total_cell_counts_df <- reactive({
    req(rv$seu_pre, rv$seu_post)
    data.frame(
      total_pre    = ncol(rv$seu_pre),
      total_post   = ncol(rv$seu_post),
      total_change = ncol(rv$seu_post) - ncol(rv$seu_pre),
      stringsAsFactors = FALSE
    )
  })

  # output$total_cell_counts_dt <- renderDT({
  #   df <- total_cell_counts_df()
  #   datatable(
  #     df,
  #     options = list(
  #       dom = 't',
  #       paging = FALSE,
  #       searching = FALSE,
  #       info = FALSE,
  #       scrollX = TRUE
  #     ),
  #     rownames = FALSE,
  #     selection = "none"
  #   )
  # })
  
  output$cluster_counts_dt <- renderDT({
    req(rv$cluster_counts_df)
    df <- rv$cluster_counts_df
    df$cluster <- as.character(df$cluster)
    lev <- sort_clusters(df$cluster)
    df <- df[order(factor(df$cluster, levels = lev)), , drop = FALSE]
    datatable(
      df,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE,
      selection = "none"
    )
  })
  
  output$dl_topgenes_csv <- downloadHandler(
    filename = function() sprintf("top_genes_%s.csv", format(Sys.Date(), "%Y%m%d")),
    content = function(file) {
      req(rv$seu_pre, rv$seu_post)
      pre_counts  <- Seurat::GetAssayData(rv$seu_pre,  layer = "counts", assay = DefaultAssay(rv$seu_pre))
      post_counts <- Seurat::GetAssayData(rv$seu_post, layer = "counts", assay = DefaultAssay(rv$seu_post))
      if (nrow(pre_counts) == 0 || nrow(post_counts) == 0)
        stop("Counts are empty.")
      pre_sum  <- Matrix::rowSums(pre_counts)
      post_sum <- Matrix::rowSums(post_counts)
      common <- intersect(rownames(pre_counts), rownames(post_counts))
      df <- data.frame(
        gene       = common,
        total_pre  = as.numeric(pre_sum[common]),
        total_post = as.numeric(post_sum[common])
      )
      df$delta <- df$total_post - df$total_pre
      df <- df[order(-df$total_pre), ]
      utils::write.csv(df, file, row.names = FALSE)
    }
  )
  
  output$dl_cells_csv <- downloadHandler(
    filename = function() sprintf("cells_table_%s.csv", format(Sys.Date(), "%Y%m%d")),
    content = function(file) {
      req(rv$cells_df)
      utils::write.csv(rv$cells_df, file, row.names = FALSE)
    }
  )

  # output$dl_total_cell_counts_csv <- downloadHandler(
  #   filename = function() sprintf("total_cell_counts_%s.csv", format(Sys.Date(), "%Y%m%d")),
  #   content = function(file) {
  #     utils::write.csv(total_cell_counts_df(), file, row.names = FALSE)
  #   }
  # )
  
  output$dl_cluster_counts_csv <- downloadHandler(
    filename = function() sprintf("cluster_counts_%s.csv", format(Sys.Date(), "%Y%m%d")),
    content = function(file) {
      req(rv$cluster_counts_df)
      df <- rv$cluster_counts_df
      df$cluster <- as.character(df$cluster)
      lev <- sort_clusters(df$cluster)
      df <- df[order(factor(df$cluster, levels = lev)), , drop = FALSE]
      utils::write.csv(df, file, row.names = FALSE)
    }
  )
  
  # output$dl_adjusted_rds <- downloadHandler(
  #   filename = function() sprintf("%s_seurat_adjusted_%s.rds",
  #                                 tolower(input$method %||% "soupx"),
  #                                 format(Sys.Date(), "%Y%m%d")),
  #   content = function(file) {
  #     req(rv$seu_post)
  #     saveRDS(rv$seu_post, file)
  #   }
  # )
  
  output$dl_adjusted_rds <- downloadHandler(
    filename = function() sprintf(
      "%s_seurat_adjusted_%s.rds",
      tolower(input$method %||% "soupx"),
      format(Sys.Date(), "%Y%m%d")
    ),
    content = function(file) {
      req(rv$seu_post)
      
      # make a copy so we don't change rv$seu_post in the app
      seu_save <- rv$seu_post
      
      # if orig.ident exists, set Idents to that before saving
      if ("orig.ident" %in% colnames(seu_save@meta.data)) {
        Seurat::Idents(seu_save) <- seu_save$orig.ident
      }
      
      saveRDS(seu_save, file)
    }
  )
  
  
  # *** IMPORTANT: cleaned .h5 / MTX zip download ***
  output$dl_adjusted_data <- downloadHandler(
    filename = function() {
      tag <- tolower(input$method %||% "soupx")
      if (identical(rv$upload_type, "h5")) {
        sprintf("%s_adjusted_%s.h5", tag, format(Sys.Date(), "%Y%m%d"))
      } else {
        sprintf("%s_adjusted_mtx_%s.zip", tag, format(Sys.Date(), "%Y%m%d"))
      }
    },
    contentType = "application/octet-stream",  # always binary
    content = function(file) {
      req(rv$adj_counts)
      
      mat <- rv$adj_counts
      if (!methods::is(mat, "dgCMatrix")) {
        mat <- as(mat, "dgCMatrix")
      }
      
      ## make sure 10x-like dimnames exist
      if (is.null(rownames(mat))) {
        if (!is.null(rv$filt) && !is.null(rownames(rv$filt))) {
          rownames(mat) <- rownames(rv$filt)[seq_len(nrow(mat))]
        } else {
          rownames(mat) <- sprintf("gene_%s", seq_len(nrow(mat)))
        }
      }
      if (is.null(colnames(mat))) {
        if (!is.null(rv$filt) && !is.null(colnames(rv$filt))) {
          colnames(mat) <- colnames(rv$filt)[seq_len(ncol(mat))]
        } else {
          colnames(mat) <- sprintf("cell_%s", seq_len(ncol(mat)))
        }
      }
      
      if (identical(rv$upload_type, "h5")) {
        ## ---- H5 cleaned output ----
        if (!requireNamespace("rhdf5", quietly = TRUE)) {
          showNotification("Package 'rhdf5' is required to write .h5 output.", type = "error")
          stop("rhdf5 not installed")
        }
        write_10x_h5(mat, file)
        
      } else {
        ## ---- MTX zip cleaned output ----
        if (!requireNamespace("zip", quietly = TRUE)) {
          showNotification("Package 'zip' is required for MTX zip download.", type = "error")
          stop("zip package not installed")
        }
        
        td <- tempfile(pattern = "mtx_")
        dir.create(td, recursive = TRUE, showWarnings = FALSE)
        
        # write matrix.mtx(.gz), barcodes.tsv(.gz), features.tsv(.gz)
        write_cleaned_10x(mat, td, gzip = TRUE)
        
        oldwd <- getwd()
        on.exit(setwd(oldwd), add = TRUE)
        setwd(td)
        
        # make a proper zip archive
        zip::zip(
          zipfile = file,
          files   = c("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")
        )
      }
    }
  )
  
  
  
  output$log  <- renderText(paste(rv$log, collapse = "\n"))
  
  # reusable download helper
  open_simple_download <- function(id_btn, id_dl, title_txt, plot_fun, fname_prefix) {
    observeEvent(id_btn(), {
      showModal(modalDialog(
        title = paste("Download:", title_txt),
        numericInput(paste0(id_dl, "_h"), "Figure height", 6),
        numericInput(paste0(id_dl, "_w"), "Figure width",  8),
        numericInput(paste0(id_dl, "_dpi"), "Figure resolution", 300),
        selectInput(paste0(id_dl, "_fmt"), "Image format",
                    choices = c(".jpg", ".tiff", ".pdf", ".svg", ".bmp", ".eps", ".ps")),
        downloadBttn(id_dl, "Download"),
        easyClose = TRUE
      ))
    })
    output[[id_dl]] <- downloadHandler(
      filename = function() {
        ext <- input[[paste0(id_dl, "_fmt")]] %||% ".jpg"
        paste0(fname_prefix, ext)
      },
      content = function(file) {
        p <- plot_fun()
        if (is.null(p)) stop("Plot is not available")
        h <- as.numeric(input[[paste0(id_dl, "_h")]]   %||% 6)
        w <- as.numeric(input[[paste0(id_dl, "_w")]]   %||% 8)
        dpi <- as.numeric(input[[paste0(id_dl, "_dpi")]] %||% 300)
        ggplot2::ggsave(filename = file, plot = p,
                        width = w, height = h, dpi = dpi, limitsize = FALSE)
      }
    )
  }
  
  # Wire up buttons
  open_simple_download(reactive(input$open_dl_pre_nFeature), "dl_pre_nFeature",
                       "QC (Pre): nFeature",      plot_pre_nFeature,     "qc_pre_nFeature")
  open_simple_download(reactive(input$open_dl_pre_nCount),   "dl_pre_nCount",
                       "QC (Pre): nCount",        plot_pre_nCount,       "qc_pre_nCount")
  open_simple_download(reactive(input$open_dl_pre_pctMT),    "dl_pre_pctMT",
                       "QC (Pre): %MT",           plot_pre_pctMT,        "qc_pre_pctMT")
  
  open_simple_download(reactive(input$open_dl_est_density),   "dl_est_density",
                       "Estimation: rho density", plot_est_density,      "est_rho_density")
  open_simple_download(reactive(input$open_dl_rho_vs_counts), "dl_rho_vs_counts",
                       "Estimation: rho vs nUMIs",plot_rho_vs_counts,    "est_rho_vs_nUMIs")
  open_simple_download(reactive(input$open_dl_auto_plot),     "dl_auto_plot",
                       "Estimation: autoEstCont", plot_auto_plot,        "est_auto")
  
  open_simple_download(
    reactive(input$open_dl_fastcar_profile), "dl_fastcar_profile",
    "FastCAR: ambient profile",
    fastcar_profile_obj,
    "fastcar_ambient_profile"
  )
  
  open_simple_download(
    reactive(input$open_dl_fastcar_removed), "dl_fastcar_removed",
    "FastCAR: UMIs removed per cell",
    fastcar_removed_hist_obj,
    "fastcar_reads_removed_per_cell"
  )
  
  open_simple_download(reactive(input$open_dl_post_nFeature), "dl_post_nFeature",
                       "QC (Post): nFeature",     plot_post_nFeature,    "qc_post_nFeature")
  open_simple_download(reactive(input$open_dl_post_nCount),   "dl_post_nCount",
                       "QC (Post): nCount",       plot_post_nCount,      "qc_post_nCount")
  open_simple_download(reactive(input$open_dl_post_pctMT),    "dl_post_pctMT",
                       "QC (Post): %MT",          plot_post_pctMT,       "qc_post_pctMT")
  
  open_simple_download(reactive(input$open_dl_cluster_counts_bar), "dl_cluster_counts_bar",
                       "Cluster counts",          plot_cluster_counts_bar,"cluster_counts_bar")
  
  open_simple_download(reactive(input$open_dl_umap_pre),  "dl_umap_pre",
                       "UMAP (Pre)",              plot_umap_pre,         "umap_pre")
  open_simple_download(reactive(input$open_dl_umap_post), "dl_umap_post",
                       "UMAP (Post)",             plot_umap_post,        "umap_post")
  
  open_simple_download(reactive(input$open_dl_tsne_pre),  "dl_tsne_pre",
                       "TSNE (Pre)",              plot_tsne_pre,         "tsne_pre")
  open_simple_download(reactive(input$open_dl_tsne_post), "dl_tsne_post",
                       "TSNE (Post)",             plot_tsne_post,        "tsne_post")
}
