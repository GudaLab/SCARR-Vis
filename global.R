# global.R — dependencies + common setup

# # ---------- Install helpers ----------
# install_if_missing <- function(pkgs, repo = c("CRAN","Bioc"), github = FALSE) {
#   repo <- match.arg(repo)
#   to_install <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
#   if (!length(to_install)) return(invisible(TRUE))
#   if (github) {
#     if (!requireNamespace("remotes", quietly = TRUE) && !requireNamespace("devtools", quietly = TRUE)) {
#       install.packages("remotes")
#     }
#     install_fun <- if (requireNamespace("remotes", quietly = TRUE)) remotes::install_github else devtools::install_github
#     install_fun(to_install, upgrade = "never", dependencies = TRUE)
#   } else if (repo == "Bioc") {
#     if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#     BiocManager::install(to_install, update = FALSE, ask = FALSE)
#   } else {
#     install.packages(to_install, dependencies = TRUE)
#   }
# }
# 
# # ---------- Ensure all packages are installed ----------
# # CRAN
# install_if_missing(c(
#   "shiny","shinythemes","shinyWidgets","shinycssloaders","DT",
#   "ggplot2","patchwork","SoupX","reshape2","zip","Seurat","Matrix"
# ), repo = "CRAN")
# 
# # Bioconductor
# install_if_missing(c(
#   "celda","SingleCellExperiment","SummarizedExperiment","biomaRt","rhdf5"
# ), repo = "Bioc")
# 
# # Some users prefer decontX as a standalone Bioc pkg; install if available
# if (!requireNamespace("decontX", quietly = TRUE)) {
#   if (requireNamespace("BiocManager", quietly = TRUE)) {
#     try(BiocManager::install("decontX", update = FALSE, ask = FALSE), silent = TRUE)
#   }
# }
# 
# # GitHub (scCDC)
# if (!requireNamespace("scCDC", quietly = TRUE)) {
#   if (!requireNamespace("remotes", quietly = TRUE) && !requireNamespace("devtools", quietly = TRUE)) {
#     install.packages("remotes")
#   }
#   install_fun <- if (requireNamespace("remotes", quietly = TRUE)) remotes::install_github else devtools::install_github
#   install_fun("ZJU-UoE-CCW-LAB/scCDC", upgrade = "never")
# }
# 
# # ---------- Load libraries ----------
# suppressPackageStartupMessages({
#   library(shiny)
#   library(shinythemes)
#   library(shinyWidgets)
#   library(shinycssloaders)
#   library(DT)
#   library(ggplot2)
#   library(patchwork)
#   library(SoupX)
#   library(Seurat)
#   library(Matrix)
#   library(reshape2)
#   library(zip)
#   library(biomaRt)
#   library(rhdf5)
#   # Use celda::decontX even if decontX package is present
#   suppressWarnings(requireNamespace("celda", quietly = TRUE))
#   # Load optional packages if present (don’t hard fail)
#   suppressWarnings(suppressMessages(requireNamespace("decontX", quietly = TRUE)))
#   suppressWarnings(suppressMessages(requireNamespace("SingleCellExperiment", quietly = TRUE)))
#   suppressWarnings(suppressMessages(requireNamespace("SummarizedExperiment", quietly = TRUE)))
#   suppressWarnings(suppressMessages(requireNamespace("scCDC", quietly = TRUE)))
#   # Base-recommended; loaded implicitly but referenced directly in code
#   # (no install needed)
#   # tools, grid
#   library(tools)
#   library(grid)
# })
if (!require("shiny")) install.packages("shiny", dependencies = TRUE)
if (!require("shinythemes")) install.packages("shinythemes")
if (!require("shinyWidgets")) install.packages("shinyWidgets")
if (!require("shinycssloaders")) install.packages("shinycssloaders")
if (!require("shinydashboard")) install.packages("shinydashboard")
if (!require("shinyjs")) install.packages("shinyjs")
if (!require("DT")) install.packages("DT")
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require("devtools")) install.packages("devtools")
if (!require("Matrix")) install.packages("Matrix")
if (!require("reshape2")) install.packages("reshape2")
if (!require("zip")) install.packages("zip")
if (!require("tools")) install.packages("tools")
if (!require("grid")) install.packages("grid")
if (!require("filelock")) install.packages("filelock")
if (!require("Seurat")) install.packages("Seurat")
if (!require("BiocManager")) install.packages("BiocManager", update = FALSE)
if (!require("biomaRt")) BiocManager::install("biomaRt", update = FALSE)
if (!require("rhdf5")) BiocManager::install("rhdf5", update = FALSE)
if (!require("celda")) BiocManager::install("celda", update = FALSE)
if (!require("decontX")) BiocManager::install("decontX", update = FALSE)
if (!require("SingleCellExperiment")) BiocManager::install("SingleCellExperiment", update = FALSE)
if (!require("SummarizedExperiment")) BiocManager::install("SummarizedExperiment", update = FALSE)
if (!require("patchwork")) install_github("thomasp85/patchwork", update = FALSE)
if (!require("SoupX"))install_github("constantAmateur/SoupX", update = FALSE)
if (!require("scCDC"))install_github("ZJU-UoE-CCW-LAB/scCDC", update = FALSE)
if (!require("FastCAR"))install_github("https://git.web.rug.nl/P278949/FastCAR", update = FALSE)

# ---------- App options ----------
options(shiny.maxRequestSize = 500 * 1024^2)

# ---------- Source helpers ----------
# Keep your custom helpers (read_10x_any, species_defaults, theme_clean, etc.) here
if (file.exists("scripts/helpers.R")) {
  source("scripts/helpers.R", chdir = TRUE)
} else if (file.exists("scripts/helpers.r")) {
  source("scripts/helpers.r", chdir = TRUE)
}

# (Optional) Global infix for null-coalesce used throughout the app
`%||%` <- function(a, b) if (is.null(a)) b else a

#views
count_file <- "view_counter.rds"
lock_file  <- "view_counter.lock"
if (!file.exists(count_file)) saveRDS(0L, count_file)

read_count <- function() {
  readRDS(count_file)
}
increment_count <- function() {
  lock <- lock(lock_file, timeout = 5000)   # wait up to 5s for the lock
  on.exit(unlock(lock), add = TRUE)
  n <- readRDS(count_file)
  n <- n + 1L
  saveRDS(n, count_file)
  n
}
