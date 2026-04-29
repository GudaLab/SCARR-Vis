if (!require("shiny")) install.packages("shiny", dependencies = TRUE)
if (!require("remotes")) install.packages("remotes", dependencies = TRUE)
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


## --- BiocFileCache: set a consistent cache directory -----------------------

# Use R's recommended cross-platform cache location for BiocFileCache
cache_path <- tools::R_user_dir("BiocFileCache", which = "cache")

# Create directory if it doesn't exist
if (!dir.exists(cache_path)) {
  dir.create(cache_path, recursive = TRUE, showWarnings = FALSE)
}

# Tell BiocFileCache to use this path
# (must be set before BiocFileCache is created/used)
Sys.setenv(BFC_CACHE = cache_path)
# Alternatively, if BiocFileCache is already loaded:
# BiocFileCache::setBFCOption("CACHE", cache_path)

## --- Download GSM7681687 example data into www/example/GSM7681687 ----------

# URL of the zip file
zip_url <- "https://www.gudalab-rtools.net/SCARR-Vis/example/GSM7681687.zip"

# Define target directory
target_dir <- file.path("www", "example", "GSM7681687")

# Check if example data already exists
if (!dir.exists(target_dir)) {
  # Create www/example/GSM7681687 if it doesn't exist
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Path for the downloaded zip file
  zip_file <- tempfile(fileext = ".zip")
  
  # Download the zip file
  download.file(zip_url, zip_file, mode = "wb")
  
  # Extract the zip file into the target directory
  unzip(zip_file, exdir = target_dir)
  
  # Remove the zip file after extraction
  unlink(zip_file)
  
  cat("Files extracted to:", normalizePath(target_dir), "\n")
} else {
  cat("Example data folder already exists at:",
      normalizePath(target_dir), "- skipping download.\n")
}

