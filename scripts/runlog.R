# =============================================================================
# SCARR-Vis - Run log
# -----------------------------------------------------------------------------
# A live, per-session process log shown in the "Run log" tab and saved as
# SCARR-Vis_run_log.txt (standalone download and inside the Bulk Download zip).
#
# What is recorded (tag shown in the log):
#   [UI]    menus (navbar tabs), Analysis sub-tabs, every input change, button
#           clicks including the Run (submit) button, and download requests
#   [CMD]   R commands run by the app (calls wrapped with run_cmd() in server.R)
#   [OUT]   everything printed to the R console (print(), cat(), Rcpp output)
#   [MSG]   messages (message(), package start-up messages, verbose output)
#   [WARN]  warnings
#   [ERROR] errors that reach the console
#   [APP]   SCARR-Vis status / progress messages
#   [SYS]   session and start-up information
#
# How it works:
#   While an observer, render function or download handler of a session runs,
#   stdout is tee'd with sink(split = TRUE), so the RStudio console still shows
#   everything, and message / warning / error conditions are observed with
#   withCallingHandlers() without muffling them. New lines are pushed to the
#   browser immediately with session$sendCustomMessage(), so the tab updates
#   while a long analysis is still running.
#   Text that C code or txtProgressBar(file = stderr()) writes straight to
#   stderr bypasses R's condition system and is not mirrored.
#
# Files: only this script (sourced first in global.R). The tab's JavaScript and
# CSS are embedded below (.rl_js / .rl_css), so nothing is needed in www/.
# =============================================================================

.scarr_runlog_global <- new.env(parent = emptyenv())
.scarr_runlog_global$startup <- list(t = character(), k = character(), x = character())
.scarr_runlog_global$startup_at <- NULL

.rl_max_lines <- 200000L

# Inputs that are internal widget state and are not logged as user actions.
# "<id>_bttn" is the button shinyWidgets::downloadBttn() wraps around a download
# link; downloads are logged by the download handler itself.
.rl_skip_input_re <- paste0(
  "^runlog_|_bttn$|",
  "_(rows_current|rows_all|rows_selected|columns_selected|cells_selected|",
  "cell_clicked|cell_info|cell_edit|cell_edit_info|state|search_columns|",
  "row_last_clicked|row_edit)$"
)

# Human-readable names for SCARR-Vis inputs (unknown ids are shown as-is)
.rl_input_labels <- c(
  main_nav               = "Menu",
  tabs                   = "Analysis sub-tab",
  species_genome         = "Species / Genome build",
  data_mode              = "Data source",
  raw                    = "Raw matrix (raw_feature_bc_matrix)",
  filt                   = "Filtered matrix (filtered_feature_bc_matrix)",
  collapse_dups          = "Collapse duplicate gene names",
  use_gene_symbols       = "Convert Ensembl IDs to gene symbols",
  method                 = "Decontamination method",
  min_cells              = "Min cells to keep (pre-filter)",
  do_auto                = "Run autoEstCont",
  manual_rho             = "Manual contamination fraction (rho)",
  do_plot                = "Show estimation plots",
  soupRange              = "soupRange (UMIs in empty droplets)",
  keepDroplets           = "keepDroplets",
  nonexp                 = "Optional markers (non-expressed genes)",
  decontx_min_cells      = "DecontX: min cells to keep (pre-filter)",
  decontx_use_clusters   = "DecontX: use Seurat clusters as priors",
  decontx_maxiter        = "DecontX: maxIter",
  decontx_delta          = "DecontX: delta",
  decontx_convergence    = "DecontX: convergence",
  decontx_iterLogLik     = "DecontX: iterLogLik",
  decontx_varGenes       = "DecontX: varGenes",
  decontx_estimateDelta  = "DecontX: estimateDelta",
  sccdc_restriction      = "scCDC: restriction_factor",
  sccdc_min_cell         = "scCDC: min.cell",
  sccdc_percent_cutoff   = "scCDC: percent.cutoff",
  fastcar_empty_cutoff   = "FastCAR: empty droplet UMI cutoff",
  fastcar_contam_cutoff  = "FastCAR: contamination chance cutoff",
  fastcar_do_profile     = "FastCAR: profile ambient RNA",
  fastcar_profile_start  = "FastCAR: profile start (UMI)",
  fastcar_profile_stop   = "FastCAR: profile stop (UMI)",
  fastcar_profile_by     = "FastCAR: profile step",
  fastcar_use_recommended = "FastCAR: use recommended cutoff",
  run                    = "Run (submit)",
  open_bulk_download     = "Bulk Download Tables + Images",
  bulk_download_fmt      = "Bulk download: image format",
  heatmap_source         = "Heatmap: show heatmap for",
  heatmap_topn           = "Heatmap: top variable genes",
  open_heatmap_download  = "Download heatmap",
  featplot_genes         = "Feature plot: gene symbols",
  open_dl_featplot       = "Download plot (Feature plot)",
  timeOut                = "Inactivity timeout"
)

.rl_plot_names <- c(
  pre_nFeature = "QC (Pre) detected genes", pre_nCount = "QC (Pre) UMIs",
  pre_pctMT = "QC (Pre) mitochondrial %", post_nFeature = "QC (Post) detected genes",
  post_nCount = "QC (Post) UMIs", post_pctMT = "QC (Post) mitochondrial %",
  est_density = "rho density", rho_vs_counts = "rho vs nUMIs",
  auto_plot = "autoEstCont diagnostic", fastcar_profile = "FastCAR ambient profile",
  fastcar_removed = "FastCAR UMIs removed", cluster_counts_bar = "cluster counts",
  umap_pre = "UMAP (Pre)", umap_post = "UMAP (Post)", tsne_pre = "TSNE (Pre)",
  tsne_post = "TSNE (Post)", decontx_contam = "DecontX contamination",
  scCDC_diag = "scCDC diagnostic", sccdc_top_gcgs_prepost = "scCDC top GCGs",
  featplot = "Feature plot"
)

.rl_input_label <- function(id) {
  if (!is.na(.rl_input_labels[id])) return(unname(.rl_input_labels[id]))
  if (grepl("^open_dl_", id)) {
    key <- sub("^open_dl_", "", id)
    nm  <- if (!is.na(.rl_plot_names[key])) .rl_plot_names[[key]] else key
    return(sprintf("Download plot (%s)", nm))
  }
  fields <- c("_h$" = "figure height", "_height$" = "figure height",
              "_w$" = "figure width",  "_width$"  = "figure width",
              "_dpi$" = "figure resolution",
              "_fmt$" = "image format", "_type$" = "image format")
  for (pat in names(fields)) {
    if (grepl(pat, id)) return(sprintf("Download dialog %s: %s", sub(pat, "", id), fields[[pat]]))
  }
  if (grepl("_search$", id)) return(sprintf("Table search (%s)", sub("_search$", "", id)))
  id
}

.rl_fmt_bytes <- function(b) {
  b <- suppressWarnings(as.numeric(b))
  if (!length(b) || is.na(b[1])) return("unknown size")
  b <- b[1]
  if (b >= 1024^3) return(sprintf("%.2f GB", b / 1024^3))
  if (b >= 1024^2) return(sprintf("%.1f MB", b / 1024^2))
  if (b >= 1024)   return(sprintf("%.1f KB", b / 1024))
  sprintf("%d bytes", as.integer(b))
}

# Value shown for a UI change ("on"/"off", "0 - 100", file names, ...)
.rl_fmt_value <- function(v) {
  if (is.null(v)) return("(none)")
  if (is.data.frame(v)) {
    if ("name" %in% names(v)) {
      return(paste(sprintf("%s (%s)", v$name, vapply(v$size, .rl_fmt_bytes, "")), collapse = ", "))
    }
    return(sprintf("<table %d x %d>", nrow(v), ncol(v)))
  }
  if (is.list(v)) return(sprintf("<list of %d>", length(v)))
  if (is.logical(v) && length(v) == 1L && !is.na(v)) return(if (v) "on" else "off")
  v <- as.character(v)
  if (!length(v)) return("(empty)")
  num <- suppressWarnings(as.numeric(v))
  s <- if (length(v) == 2L && !anyNA(num)) paste(v, collapse = " - ") else paste(v, collapse = ", ")
  if (!nzchar(s)) s <- "\"\""
  if (nchar(s) > 200L) s <- paste0(substr(s, 1L, 197L), "...")
  s
}

# Value shown in the Run settings summary (R syntax, easy to copy)
.rl_fmt_param <- function(v) {
  if (is.null(v)) return("NULL")
  if (is.data.frame(v) && "name" %in% names(v)) return(paste0("\"", v$name, "\"", collapse = ", "))
  if (!is.atomic(v)) return(sprintf("<%s>", class(v)[1]))
  attributes(v) <- NULL
  paste(deparse(v, width.cutoff = 200L, control = c("keepNA", "niceNames")), collapse = "")
}

# Split text into display lines; for "\r" progress updates keep the final state
.rl_clean_lines <- function(text) {
  text <- as.character(text)
  if (!length(text)) return(character())
  text <- paste(text, collapse = "\n")
  text <- sub("\n$", "", text)
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  if (!length(lines)) lines <- ""
  lines <- vapply(lines, function(l) {
    if (!grepl("\r", l, fixed = TRUE)) return(l)
    segs <- strsplit(l, "\r", fixed = TRUE)[[1]]
    segs <- segs[nzchar(trimws(segs))]
    if (length(segs)) segs[length(segs)] else ""
  }, character(1), USE.NAMES = FALSE)
  sub("[[:space:]]+$", "", lines)
}

.rl_condition_text <- function(prefix, cond, context = NULL) {
  msg  <- tryCatch(conditionMessage(cond), error = function(e) "")
  call <- tryCatch(conditionCall(cond), error = function(e) NULL)
  ctxt <- ""
  internal <- c("withCallingHandlers", "eval", "doTryCatch", "force")
  if (!is.null(call) && is.call(call) && !is.function(call[[1]]) &&
      !(is.name(call[[1]]) && as.character(call[[1]]) %in% internal)) {
    ctxt <- tryCatch(paste(deparse(call, nlines = 1L, width.cutoff = 80L), collapse = " "),
                     error = function(e) "")
  }
  out <- if (nzchar(ctxt)) sprintf("%s in %s : %s", prefix, ctxt, msg) else sprintf("%s: %s", prefix, msg)
  if (!is.null(context) && nzchar(context)) out <- sprintf("%s   [%s]", out, context)
  out
}

.rl_ignorable_error <- function(e) {
  inherits(e, c("shiny.silent.error", "shiny.output.cancel", "shiny.output.progress", "validation"))
}

# Lines printed to stdout since the last call (complete lines only)
.rl_take_stdout <- function(st) {
  if (is.null(st$tc)) return(character())
  v <- tryCatch(textConnectionValue(st$tc), error = function(e) character())
  n <- length(v)
  if (n <= st$tc_read) return(character())
  out <- v[seq.int(st$tc_read + 1L, n)]
  st$tc_read <- n
  out
}

# Core capture: tee stdout and observe conditions while `expr` runs.
# Nothing is muffled, so the R console shows exactly what it showed before.
.rl_capture <- function(st, emit, label, expr) {
  if (st$depth > 0L) return(expr)   # an outer capture is already active
  st$depth <- 1L
  n_sinks <- sink.number()
  tc <- textConnection(NULL, open = "w")
  sink(tc, split = TRUE)
  st$tc <- tc
  st$tc_read <- 0L
  on.exit({
    if (isTRUE(tryCatch(isIncomplete(tc), error = function(e) FALSE))) {
      try(cat("\n", file = tc), silent = TRUE)   # complete a trailing partial line
    }
    while (sink.number() > n_sinks) sink()
    emit(NULL)                                    # flush remaining stdout lines
    try(close(tc), silent = TRUE)
    st$tc <- NULL
    st$depth <- 0L
  }, add = TRUE)
  withCallingHandlers(
    expr,
    message = function(m) emit("MSG", tryCatch(conditionMessage(m), error = function(e) "")),
    warning = function(w) emit("WARN", .rl_condition_text("Warning", w, label)),
    error   = function(e) {
      if (!.rl_ignorable_error(e)) emit("ERROR", .rl_condition_text("Error", e, label))
    }
  )
}

# -----------------------------------------------------------------------------
# Start-up capture (used in global.R): records what global.R prints so it can
# be shown at the top of every session's Run log.
# -----------------------------------------------------------------------------
scarr_runlog_startup <- function(expr) {
  g <- .scarr_runlog_global
  g$startup <- list(t = character(), k = character(), x = character())
  g$startup_at <- Sys.time()
  st <- new.env(parent = emptyenv())
  st$depth <- 0L; st$busy <- FALSE; st$tc <- NULL; st$tc_read <- 0L
  emit <- function(kind, text = character()) {
    if (st$busy) return(invisible(NULL))
    st$busy <- TRUE
    on.exit(st$busy <- FALSE)
    out <- .rl_clean_lines(.rl_take_stdout(st))
    new_k <- rep("OUT", length(out)); new_x <- out
    if (!is.null(kind)) {
      x <- .rl_clean_lines(text)
      new_k <- c(new_k, rep(kind, length(x))); new_x <- c(new_x, x)
    }
    if (length(new_x)) {
      g$startup$t <- c(g$startup$t, rep(format(Sys.time(), "%H:%M:%S"), length(new_x)))
      g$startup$k <- c(g$startup$k, new_k)
      g$startup$x <- c(g$startup$x, new_x)
    }
    invisible(NULL)
  }
  invisible(.rl_capture(st, emit, "global.R", expr))
}

# -----------------------------------------------------------------------------
# UI: content of the "Run log" navbar tab
# -----------------------------------------------------------------------------
# Styles and script for the Run log tab. They are embedded here (not served from
# www/) so the tab works on any server and the browser never uses a stale copy.
.rl_css <- c(
  "/* SCARR-Vis Run log styles */",
  ".rl-head h3 { margin-top: 18px; }",
  ".rl-head p { max-width: 900px; }",
  "",
  ".rl-toolbar {",
  "  display: flex;",
  "  flex-wrap: wrap;",
  "  align-items: center;",
  "  gap: 10px 18px;",
  "  margin: 8px 0 10px;",
  "}",
  ".rl-toolbar .btn-group .badge { margin-left: 2px; }",
  ".rl-check {",
  "  font-weight: normal;",
  "  margin: 0;",
  "  cursor: pointer;",
  "  user-select: none;",
  "}",
  ".rl-check input { margin: 0 5px 0 0; vertical-align: -2px; }",
  ".rl-status { color: #6b7280; font-size: 13px; white-space: nowrap; }",
  ".rl-dot {",
  "  display: inline-block;",
  "  width: 8px; height: 8px;",
  "  border-radius: 50%;",
  "  margin-right: 6px;",
  "  background: #9ca3af;",
  "  vertical-align: 1px;",
  "}",
  ".rl-status[data-state=\"idle\"] .rl-dot { background: #22a06b; }",
  ".rl-status[data-state=\"busy\"] .rl-dot { background: #e2a400; animation: rl-pulse 1s ease-in-out infinite; }",
  ".rl-status[data-state=\"off\"]  .rl-dot { background: #c9372c; }",
  "@keyframes rl-pulse { 50% { opacity: .35; } }",
  ".rl-download { margin-left: auto; }",
  "",
  ".rl-console {",
  "  background: #1e1f22;",
  "  color: #d4d4d4;",
  "  font-family: \"Cascadia Mono\", Consolas, \"SFMono-Regular\", Menlo, \"Liberation Mono\", monospace;",
  "  font-size: 12.5px;",
  "  line-height: 1.5;",
  "  border: 1px solid #2b2d31;",
  "  border-radius: 6px;",
  "  height: calc(100vh - 300px);",
  "  min-height: 380px;",
  "  overflow: auto;",
  "  padding: 8px 0;",
  "}",
  ".rl-empty { color: #8b949e; padding: 4px 14px; font-style: italic; }",
  ".rl-line {",
  "  display: flex;",
  "  gap: 12px;",
  "  padding: 0 14px;",
  "  white-space: pre-wrap;",
  "  word-break: break-word;",
  "}",
  ".rl-console.rl-nowrap .rl-line { white-space: pre; word-break: normal; }",
  ".rl-line:hover { background: rgba(255, 255, 255, 0.04); }",
  ".rl-t { color: #6b7280; flex: none; }",
  ".rl-tag {",
  "  flex: none;",
  "  width: 5ch;",
  "  font-size: 11px;",
  "  font-weight: 600;",
  "  letter-spacing: .02em;",
  "  padding-top: 1px;",
  "  color: #8b949e;",
  "}",
  ".rl-x { flex: 1 1 auto; min-width: 0; }",
  "",
  ".rl-k-UI    .rl-tag { color: #4ec9b0; }",
  ".rl-k-UI    .rl-x   { color: #bdeee3; }",
  ".rl-k-CMD   .rl-tag { color: #569cd6; }",
  ".rl-k-CMD   .rl-x   { color: #9cdcfe; }",
  ".rl-k-OUT   .rl-x   { color: #d4d4d4; }",
  ".rl-k-MSG   .rl-tag { color: #b392f0; }",
  ".rl-k-MSG   .rl-x   { color: #cfcfcf; }",
  ".rl-k-WARN  .rl-tag { color: #e5c07b; }",
  ".rl-k-WARN  .rl-x   { color: #f2d58f; }",
  ".rl-k-ERROR { background: rgba(244, 135, 113, 0.10); }",
  ".rl-k-ERROR .rl-tag { color: #f48771; }",
  ".rl-k-ERROR .rl-x   { color: #ffb4a8; }",
  ".rl-k-APP   .rl-tag { color: #dcdcaa; }",
  ".rl-k-APP   .rl-x   { color: #ecebc6; }",
  ".rl-k-SYS   .rl-x   { color: #9ca3af; font-style: italic; }",
  "",
  "/* A new Run starts a visually separated block */",
  ".rl-mark { border-top: 1px solid #3a3d44; margin-top: 6px; padding-top: 6px; }",
  "",
  "/* Filters */",
  ".rl-console[data-filter=\"console\"] .rl-line:not(.rl-k-CMD):not(.rl-k-OUT):not(.rl-k-MSG):not(.rl-k-WARN):not(.rl-k-ERROR) { display: none; }",
  ".rl-console[data-filter=\"ui\"]      .rl-line:not(.rl-k-UI) { display: none; }",
  ".rl-console[data-filter=\"app\"]     .rl-line:not(.rl-k-APP):not(.rl-k-SYS) { display: none; }",
  ".rl-console[data-filter=\"issues\"]  .rl-line:not(.rl-k-WARN):not(.rl-k-ERROR) { display: none; }"
)

.rl_js <- c(
  "/* SCARR-Vis Run log - live console view (see scripts/runlog.R) */",
  "(function () {",
  "  'use strict';",
  "",
  "  var MAX_DOM_LINES = 20000;   // older lines stay in the downloadable .txt",
  "  var RUNLOG_TAB = 'Run log';",
  "  var pending = [];",
  "  var issues = 0;",
  "",
  "  function arr(a) { return a == null ? [] : (Array.isArray(a) ? a : [a]); }",
  "  function root() { return document.getElementById('scarr-runlog'); }",
  "  function panel() { var el = root(); return el ? el.closest('.tab-pane') || document : document; }",
  "  function q(sel) { var p = panel(); return p ? p.querySelector(sel) : null; }",
  "",
  "  function autoScrollOn() {",
  "    var cb = q('.rl-autoscroll');",
  "    return !cb || cb.checked;",
  "  }",
  "  function nearBottom(el) {",
  "    return el.scrollHeight - el.scrollTop - el.clientHeight < 48;",
  "  }",
  "  function makeLine(t, k, x, mark) {",
  "    var d = document.createElement('div');",
  "    d.className = 'rl-line rl-k-' + k + (mark ? ' rl-mark' : '');",
  "    var st = document.createElement('span');",
  "    st.className = 'rl-t';",
  "    st.textContent = t;",
  "    var sk = document.createElement('span');",
  "    sk.className = 'rl-tag';",
  "    sk.textContent = k;",
  "    var sx = document.createElement('span');",
  "    sx.className = 'rl-x';",
  "    sx.textContent = x;",
  "    d.appendChild(st);",
  "    d.appendChild(sk);",
  "    d.appendChild(sx);",
  "    return d;",
  "  }",
  "",
  "  function apply(el, msg) {",
  "    var follow = autoScrollOn() && nearBottom(el);",
  "    var t = arr(msg.t), k = arr(msg.k), x = arr(msg.x), m = arr(msg.m);",
  "    var empty = el.querySelector('.rl-empty');",
  "    if (empty && x.length) empty.remove();",
  "",
  "    if (msg.op === 'replace_last') {",
  "      var last = el.lastElementChild;",
  "      if (last && last.classList.contains('rl-line') && x.length) {",
  "        last.querySelector('.rl-t').textContent = t[0];",
  "        last.querySelector('.rl-x').textContent = x[0];",
  "      }",
  "    } else {",
  "      var frag = document.createDocumentFragment();",
  "      var newIssues = 0;",
  "      for (var i = 0; i < x.length; i++) {",
  "        frag.appendChild(makeLine(t[i], k[i], x[i], m[i] === true));",
  "        if (k[i] === 'WARN' || k[i] === 'ERROR') newIssues++;",
  "      }",
  "      el.appendChild(frag);",
  "      var extra = el.childElementCount - MAX_DOM_LINES;",
  "      while (extra-- > 0 && el.firstElementChild) el.removeChild(el.firstElementChild);",
  "      if (newIssues) {",
  "        issues += newIssues;",
  "        var ic = q('.rl-issue-count');",
  "        if (ic) ic.textContent = String(issues);",
  "      }",
  "    }",
  "    var c = q('.rl-count');",
  "    if (c && msg.n != null) c.textContent = Number(msg.n).toLocaleString();",
  "    if (follow) el.scrollTop = el.scrollHeight;",
  "  }",
  "",
  "  function handle(msg) {",
  "    var el = root();",
  "    if (!el) { pending.push(msg); return; }",
  "    while (pending.length) apply(el, pending.shift());",
  "    apply(el, msg);",
  "  }",
  "",
  "  function setState(text, cls) {",
  "    var s = q('.rl-state');",
  "    var st = q('.rl-status');",
  "    if (s) s.textContent = text;",
  "    if (st) st.setAttribute('data-state', cls);",
  "  }",
  "",
  "  function register() {",
  "    if (!window.Shiny || !Shiny.addCustomMessageHandler) return false;",
  "    Shiny.addCustomMessageHandler('scarr-runlog', handle);",
  "    return true;",
  "  }",
  "  if (!register()) {",
  "    document.addEventListener('DOMContentLoaded', register);",
  "  }",
  "",
  "  $(document).on('shiny:connected', function () { setState('Live', 'idle'); });",
  "  $(document).on('shiny:busy', function () { setState('Running\\u2026', 'busy'); });",
  "  $(document).on('shiny:idle', function () { setState('Live', 'idle'); });",
  "  $(document).on('shiny:disconnected', function () { setState('Disconnected', 'off'); });",
  "",
  "  // Filters",
  "  $(document).on('click', '.rl-filters [data-rl-filter]', function () {",
  "    var el = root();",
  "    if (!el) return;",
  "    $(this).addClass('active').siblings().removeClass('active');",
  "    el.setAttribute('data-filter', this.getAttribute('data-rl-filter'));",
  "    if (autoScrollOn()) el.scrollTop = el.scrollHeight;",
  "  });",
  "",
  "  // Line wrapping",
  "  $(document).on('change', '.rl-wrap', function () {",
  "    var el = root();",
  "    if (el) el.classList.toggle('rl-nowrap', !this.checked);",
  "  });",
  "",
  "  // Re-enable auto-scroll: jump to the newest line",
  "  $(document).on('change', '.rl-autoscroll', function () {",
  "    var el = root();",
  "    if (el && this.checked) el.scrollTop = el.scrollHeight;",
  "  });",
  "",
  "  // Opening the tab: show the newest lines",
  "  $(document).on('shown.bs.tab', 'a[data-toggle=\"tab\"], a[data-bs-toggle=\"tab\"]', function (e) {",
  "    if (e.target.getAttribute('data-value') !== RUNLOG_TAB) return;",
  "    var el = root();",
  "    if (el && autoScrollOn()) el.scrollTop = el.scrollHeight;",
  "  });",
  "})();"
)

scarr_runlog_ui <- function() {
  filter_btn <- function(value, label, active = FALSE) {
    tags$button(
      type = "button", class = paste("btn btn-default btn-sm", if (active) "active"),
      `data-rl-filter` = value, label
    )
  }
  tagList(
    tags$head(
      tags$style(HTML(paste(.rl_css, collapse = "\n"))),
      tags$script(HTML(paste(.rl_js, collapse = "\n")))
    ),
    fluidRow(
      column(
        width = 10, offset = 1,
        div(
          class = "rl-head",
          h3("Run log"),
          p(
            class = "text-muted",
            "Live process log for this session: every menu, sub-tab, setting and button you use ",
            "(including Run), the R commands SCARR-Vis runs, and everything they print to the R console ",
            "- output, messages, warnings and errors. It is also saved as ",
            tags$code("SCARR-Vis_run_log.txt"), " in the Bulk Download zip."
          )
        ),
        div(
          class = "rl-toolbar",
          div(
            class = "btn-group rl-filters", role = "group", `aria-label` = "Filter log lines",
            filter_btn("all", "All", active = TRUE),
            filter_btn("console", "R console"),
            filter_btn("ui", "UI actions"),
            filter_btn("app", "Status"),
            tags$button(
              type = "button", class = "btn btn-default btn-sm", `data-rl-filter` = "issues",
              "Warnings & errors ", tags$span(class = "badge rl-issue-count", "0")
            )
          ),
          tags$label(class = "rl-check",
                     tags$input(type = "checkbox", class = "rl-autoscroll", checked = NA,
                                `data-shiny-no-bind-input` = NA),
                     "Auto-scroll"),
          tags$label(class = "rl-check",
                     tags$input(type = "checkbox", class = "rl-wrap", checked = NA,
                                `data-shiny-no-bind-input` = NA),
                     "Wrap lines"),
          tags$span(class = "rl-status", tags$span(class = "rl-dot"),
                    tags$span(class = "rl-state", "Connecting"), HTML(" &middot; "),
                    tags$span(class = "rl-count", "0"), " lines"),
          div(class = "rl-download",
              shinyWidgets::downloadBttn("download_runlog", "Download run log (.txt)",
                                         style = "unite", color = "primary", size = "sm"))
        ),
        div(
          id = "scarr-runlog", class = "rl-console", role = "log", `data-filter` = "all",
          div(class = "rl-empty", HTML("Waiting for activity&hellip;"))
        ),
        br()
      )
    )
  )
}

# -----------------------------------------------------------------------------
# Server: call once at the top of server(). Returns helpers used in server.R.
# -----------------------------------------------------------------------------
scarr_runlog_start <- function(input, output, session) {
  st <- new.env(parent = emptyenv())
  st$depth <- 0L; st$busy <- FALSE; st$tc <- NULL; st$tc_read <- 0L
  st$t <- character(1024L); st$k <- character(1024L); st$x <- character(1024L)
  st$m <- logical(1024L)
  st$n <- 0L; st$dropped <- 0L
  st$started <- Sys.time()
  st$last_key <- NULL; st$last_n <- -1L; st$last_time <- 0; st$last_old <- NULL

  store <- function(t, k, x, m) {
    add <- length(x)
    if (!add) return(invisible(NULL))
    need <- st$n + add
    if (need > length(st$x)) {
      len <- max(need, 2L * length(st$x))
      length(st$t) <- len; length(st$k) <- len; length(st$x) <- len; length(st$m) <- len
    }
    idx <- seq.int(st$n + 1L, need)
    st$t[idx] <- t; st$k[idx] <- k; st$x[idx] <- x; st$m[idx] <- m
    st$n <- need
    if (st$n > .rl_max_lines) {
      drop <- st$n - as.integer(.rl_max_lines * 0.9)
      keep <- seq.int(drop + 1L, st$n)
      st$t <- st$t[keep]; st$k <- st$k[keep]; st$x <- st$x[keep]; st$m <- st$m[keep]
      st$n <- length(keep); st$dropped <- st$dropped + drop
    }
    invisible(NULL)
  }

  send <- function(op, t, k, x, m) {
    try(session$sendCustomMessage("scarr-runlog", list(
      op = op, t = I(t), k = I(k), x = I(x), m = I(m), n = st$n + st$dropped
    )), silent = TRUE)
  }

  # Add lines of one kind (plus any stdout printed since the last call) and push
  # them to the browser. `change` = list(key, label, old, new) for an input
  # change: repeated changes of the same input within 3 s (typing, slider
  # drags) are collapsed into one line that keeps the original old value.
  emit <- function(kind, text = character(), mark = FALSE, change = NULL) {
    if (st$busy) return(invisible(NULL))
    st$busy <- TRUE
    on.exit(st$busy <- FALSE)
    now <- format(Sys.time(), "%H:%M:%S")
    out <- .rl_clean_lines(.rl_take_stdout(st))
    if (length(out)) {
      store(rep(now, length(out)), rep("OUT", length(out)), out, rep(FALSE, length(out)))
      send("append", rep(now, length(out)), rep("OUT", length(out)), out, rep(FALSE, length(out)))
    }
    if (is.null(kind)) return(invisible(NULL))
    if (!is.null(change)) {
      t_now <- as.numeric(Sys.time())
      same <- identical(change$key, st$last_key) && st$last_n == st$n && (t_now - st$last_time) < 3
      st$last_time <- t_now
      if (same) {
        x <- sprintf("%s: %s -> %s", change$label, st$last_old, change$new)
        st$t[st$n] <- now; st$x[st$n] <- x
        send("replace_last", now, kind, x, FALSE)
        return(invisible(NULL))
      }
      st$last_key <- change$key
      st$last_old <- change$old
      text <- sprintf("%s: %s -> %s", change$label, change$old, change$new)
    } else {
      st$last_key <- NULL
    }
    x <- .rl_clean_lines(text)
    if (!length(x)) return(invisible(NULL))
    ms <- c(isTRUE(mark), rep(FALSE, length(x) - 1L))
    store(rep(now, length(x)), rep(kind, length(x)), x, ms)
    send("append", rep(now, length(x)), rep(kind, length(x)), x, ms)
    st$last_n <- st$n
    invisible(NULL)
  }

  capture <- function(label, expr) .rl_capture(st, emit, label, expr)

  # ---- full log as text (header + one line per entry) ----
  text_lines <- function() {
    n <- st$n
    body <- if (n) sprintf("%s %-7s %s", st$t[seq_len(n)], paste0("[", st$k[seq_len(n)], "]"), st$x[seq_len(n)]) else character()
    c(
      "SCARR-Vis run log",
      sprintf("Generated:       %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      sprintf("Session started: %s", format(st$started, "%Y-%m-%d %H:%M:%S %Z")),
      sprintf("R: %s | shiny %s | %s", R.version.string, as.character(utils::packageVersion("shiny")),
              R.version$platform),
      sprintf("Lines: %d%s", n, if (st$dropped) sprintf(" (the %d earliest lines were dropped)", st$dropped) else ""),
      paste("Legend: [UI] menus, sub-tabs, inputs, buttons (incl. Run/submit) and downloads |",
            "[CMD] R commands run by SCARR-Vis | [OUT] console output | [MSG] messages |",
            "[WARN] warnings | [ERROR] errors | [APP] SCARR-Vis status | [SYS] session / start-up"),
      strrep("-", 100),
      body
    )
  }
  write_txt <- function(path) {
    writeLines(enc2utf8(text_lines()), path, useBytes = TRUE)
    invisible(path)
  }

  # ---- run_cmd(): echo an R command, then run it. In the echo, input$... and
  # simple local values (numbers, text, TRUE/FALSE) are replaced by their values
  # so the logged command shows the parameters actually used. ----
  simple_value <- function(v) {
    if (is.null(v) || !is.atomic(v) || length(v) > 6L || !is.null(attributes(v))) return(NULL)
    if (!(is.numeric(v) || is.character(v) || is.logical(v))) return(NULL)
    if (is.character(v) && any(nchar(v) > 200L)) return(NULL)
    v
  }
  resolve_values <- function(e, env) {
    if (is.name(e)) {
      nm <- as.character(e)
      if (!nzchar(nm) || !exists(nm, envir = env, inherits = TRUE)) return(e)
      v <- tryCatch(simple_value(get(nm, envir = env, inherits = TRUE)), error = function(err) NULL)
      return(if (is.null(v)) e else v)
    }
    if (!is.call(e)) return(e)
    head <- e[[1]]
    if (length(e) == 3L && (identical(head, as.name("$")) || identical(head, as.name("[["))) &&
        identical(e[[2]], as.name("input"))) {
      nm <- if (is.name(e[[3]])) as.character(e[[3]]) else if (is.character(e[[3]])) e[[3]] else NULL
      if (!is.null(nm)) {
        v <- tryCatch(shiny::isolate(input[[nm]]), error = function(err) NULL)
        if (!is.null(v) && is.atomic(v) && length(v) <= 6L) { attributes(v) <- NULL; return(v) }
      }
      return(e)
    }
    skip <- integer()
    if (is.name(head) && as.character(head) %in% c("<-", "=", "<<-")) skip <- 2L   # assignment target
    if (is.name(head) && as.character(head) %in% c("$", "@")) skip <- 3L           # element name
    for (i in setdiff(seq_along(e)[-1L], skip)) {
      if (identical(e[[i]], quote(expr = )) || is.null(e[[i]])) next   # empty arg, e.g. x[, 1]
      e[[i]] <- resolve_values(e[[i]], env)
    }
    e
  }
  cmd <- function(expr) {
    e <- substitute(expr)
    env <- parent.frame()
    shown <- tryCatch(resolve_values(e, env), error = function(err) e)
    txt <- tryCatch(deparse(shown, width.cutoff = 100L, control = c("keepNA", "niceNames", "showAttributes")),
                    error = function(err) "<command>")
    emit("CMD", c(paste0("> ", txt[1]), if (length(txt) > 1L) paste0("+ ", txt[-1L])))
    eval(e, env)
  }

  # ---- wrappers for observers / renders / downloads (assigned in server.R) ----
  wrap_expr <- function(expr, label) as.call(list(capture, label, expr))
  deparse1_short <- function(e) {
    s <- tryCatch(paste(deparse(e, width.cutoff = 60L), collapse = " "), error = function(err) "?")
    if (nchar(s) > 80L) paste0(substr(s, 1L, 77L), "...") else s
  }

  observeEvent_w <- function(eventExpr, handlerExpr, event.env = parent.frame(), event.quoted = FALSE,
                             handler.env = parent.frame(), handler.quoted = FALSE, ...) {
    if (!event.quoted) eventExpr <- substitute(eventExpr)
    if (!handler.quoted) handlerExpr <- substitute(handlerExpr)
    lbl <- sprintf("observeEvent(%s)", deparse1_short(eventExpr))
    shiny::observeEvent(eventExpr, wrap_expr(handlerExpr, lbl),
                        event.env = event.env, event.quoted = TRUE,
                        handler.env = handler.env, handler.quoted = TRUE, ...)
  }

  observe_w <- function(x, env = parent.frame(), quoted = FALSE, ..., label = NULL) {
    if (!quoted) x <- substitute(x)
    if (is.null(label)) label <- sprintf("observe(%s)", deparse1_short(x))
    shiny::observe(wrap_expr(x, label), env = env, quoted = TRUE, ..., label = label)
  }

  wrap_render_fun <- function(rf, what) {
    if (!is.function(rf)) return(rf)
    out <- function(...) {
      nm <- tryCatch(list(...)[["name"]], error = function(e) NULL)
      capture(if (is.character(nm)) paste0("output$", nm) else what, rf(...))
    }
    attrs <- attributes(rf)
    attrs$srcref <- NULL
    attributes(out) <- attrs
    out
  }
  make_render <- function(pkg, fun_name) {
    force(pkg); force(fun_name)
    function(expr, ..., env = parent.frame(), quoted = FALSE) {
      if (!quoted) expr <- substitute(expr)
      fun <- getExportedValue(pkg, fun_name)
      wrap_render_fun(fun(expr, ..., env = env, quoted = TRUE), fun_name)
    }
  }

  downloadHandler_w <- function(filename, content, ...) {
    force(filename); force(content)
    wrapped_content <- function(file) {
      fname <- tryCatch(if (is.function(filename)) filename() else as.character(filename),
                        error = function(e) "download")
      capture(paste("download", fname), {
        emit("UI", sprintf("Download requested: %s", fname))
        res <- content(file)
        emit("APP", sprintf("Download ready: %s (%s)", fname,
                            .rl_fmt_bytes(tryCatch(file.size(file), error = function(e) NA))))
        res
      })
    }
    shiny::downloadHandler(filename = filename, content = wrapped_content, ...)
  }

  # ---- UI action log (menus, sub-tabs, inputs, buttons) ----
  run_settings <- function() {
    i <- function(id) input[[id]]
    f <- function(ids) paste(vapply(ids, function(id) sprintf("%s=%s", id, .rl_fmt_param(i(id))), ""),
                             collapse = ", ")
    method <- i("method") %||% "?"
    src <- if (identical(i("data_mode"), "example")) {
      "example data (GSM7681687)"
    } else {
      sprintf("uploaded files (raw: %s; filtered: %s)",
              if (is.data.frame(i("raw"))) i("raw")$name[1] else "none",
              if (is.data.frame(i("filt"))) i("filt")$name[1] else "none")
    }
    params <- switch(
      method,
      SoupX   = c("min_cells", "do_auto", if (!isTRUE(i("do_auto"))) "manual_rho", "do_plot",
                  "soupRange", "keepDroplets", "nonexp"),
      DecontX = c("decontx_min_cells", "decontx_use_clusters", "decontx_maxiter", "decontx_delta",
                  "decontx_convergence", "decontx_iterLogLik", "decontx_varGenes", "decontx_estimateDelta"),
      scCDC   = c("min_cells", "sccdc_restriction", "sccdc_min_cell", "sccdc_percent_cutoff"),
      FastCAR = c("min_cells", "fastcar_empty_cutoff", "fastcar_contam_cutoff", "fastcar_do_profile",
                  if (isTRUE(i("fastcar_do_profile")))
                    c("fastcar_profile_start", "fastcar_profile_stop", "fastcar_profile_by",
                      "fastcar_use_recommended")),
      character()
    )
    c(
      sprintf("  Data source: %s | Species: %s | Method: %s", src, i("species_genome") %||% "?", method),
      sprintf("  Options: %s", f(c("collapse_dups", "use_gene_symbols"))),
      if (length(params)) sprintf("  %s parameters: %s", method, f(params))
    )
  }

  log_input_change <- function(id, old, new) {
    label <- .rl_input_label(id)
    if (identical(id, "main_nav")) {
      emit("UI", sprintf("Menu -> %s", .rl_fmt_value(new)))
    } else if (identical(id, "tabs")) {
      emit("UI", sprintf("Analysis sub-tab -> %s", .rl_fmt_value(new)))
    } else if (inherits(new, "shinyActionButtonValue")) {
      n_new <- suppressWarnings(as.integer(new)); n_old <- suppressWarnings(as.integer(old))
      if (length(n_new) && !is.na(n_new) && (!length(n_old) || is.na(n_old) || n_new > n_old)) {
        if (identical(id, "run")) {
          emit("UI", sprintf("Run (submit) clicked - run #%d", n_new), mark = TRUE)
          emit("UI", run_settings())
        } else {
          emit("UI", sprintf("Clicked: %s", label))
        }
      }
    } else if (identical(id, "timeOut")) {
      emit("UI", sprintf("Session timed out after %s of inactivity", .rl_fmt_value(new)))
    } else if (is.data.frame(new) && "name" %in% names(new)) {
      emit("UI", sprintf("File selected - %s: %s", label, .rl_fmt_value(new)))
    } else {
      emit("UI", change = list(key = id, label = label,
                               old = .rl_fmt_value(old), new = .rl_fmt_value(new)))
    }
  }

  prev_inputs <- NULL
  shiny::observe(label = "runlog: UI actions", priority = 1e6, {
    cur <- shiny::reactiveValuesToList(input)
    cur <- cur[!grepl(.rl_skip_input_re, names(cur))]
    shiny::isolate({
      if (is.null(prev_inputs)) {
        if (!is.null(cur$main_nav)) emit("UI", sprintf("Opened menu: %s", .rl_fmt_value(cur$main_nav)))
      } else {
        ids <- names(cur)[names(cur) %in% names(prev_inputs)]   # new inputs (e.g. a dialog opened) are a baseline
        ids <- c(intersect(c("main_nav", "tabs"), ids), sort(setdiff(ids, c("main_nav", "tabs"))))
        for (id in ids) {
          if (!identical(prev_inputs[[id]], cur[[id]])) log_input_change(id, prev_inputs[[id]], cur[[id]])
        }
      }
      prev_inputs <<- cur
    })
  })

  # ---- standalone download of the run log ----
  output$download_runlog <- shiny::downloadHandler(
    filename = function() sprintf("SCARR-Vis_run_log_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) {
      emit("UI", "Download requested: run log (.txt)")
      write_txt(file)
    },
    contentType = "text/plain"
  )

  # ---- session header + start-up output ----
  emit("SYS", sprintf("SCARR-Vis run log - session started %s",
                      format(st$started, "%Y-%m-%d %H:%M:%S %Z")))
  emit("SYS", sprintf("%s | shiny %s | %s", R.version.string,
                      as.character(utils::packageVersion("shiny")), R.version$platform))
  su <- .scarr_runlog_global$startup
  if (length(su$x)) {
    emit("SYS", sprintf("---- R console output from app start-up (global.R, %s) ----",
                        format(.scarr_runlog_global$startup_at, "%Y-%m-%d %H:%M:%S")))
    store(su$t, su$k, su$x, rep(FALSE, length(su$x)))
    send("append", su$t, su$k, su$x, rep(FALSE, length(su$x)))
    emit("SYS", "---- end of start-up output ----")
  }

  list(
    capture         = capture,
    cmd             = cmd,
    app             = function(...) emit("APP", paste(..., collapse = " ")),
    ui              = function(...) emit("UI", paste(..., collapse = " ")),
    sys             = function(...) emit("SYS", paste(..., collapse = " ")),
    text            = text_lines,
    write_txt       = write_txt,
    observeEvent    = observeEvent_w,
    observe         = observe_w,
    renderPlot      = make_render("shiny", "renderPlot"),
    renderUI        = make_render("shiny", "renderUI"),
    renderText      = make_render("shiny", "renderText"),
    renderPrint     = make_render("shiny", "renderPrint"),
    renderDT        = make_render("DT", "renderDT"),
    downloadHandler = downloadHandler_w
  )
}
