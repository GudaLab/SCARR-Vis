# SCARR-Vis

SCARR-Vis is an R/Shiny application for interactive assessment and correction of ambient RNA contamination in single-cell and single nucleus RNA-seq data. The interface follows a typical workflow: upload 10x matrices, choose a decontamination method, inspect pre- and post-correction QC, and explore clustering and gene expression.

## Features

Estimation Method and Parameters
1. SoupX: models background 'soup' RNA and adjusts counts.
2. DecontX: from celda infers cell-specific contamination fractions.
3.scCDC: identifies contamination-causing genes (GCGs) and optionally corrects them.
4. FastCAR: profiles ambient RNA using empty droplets and estimates per-cell contamination.

## Use SCARR-Vis Online
ScRDAVis is deployed online and accessible at:  
**[https://www.gudalab-rtools.net/SCARR-Vis](https://www.gudalab-rtools.net/SCARR-Vis)**

## Launch SCARR-Vis Locally

### Prerequisites
Ensure the following software is installed:
- **R** (>= 4.5.2): [Download R](https://www.r-project.org/)
- **RStudio** (>= 2025.09.2): [Download RStudio](https://posit.co/download/rstudio-desktop/)
- **Bioconductor** (>= 3.22)
- **Shiny** (>= 1.11.1)

**Note:** SCARR-Vis has been tested with these versions. Using older versions of R may cause errors during package installation. Updating to the latest R version is recommended.

### Installation

Run the following commands in an R session to install required packages:

```R
if (!require("BiocManager")) install.packages("BiocManager", update = FALSE)
if (!require("devtools")) install.packages("devtools", update = FALSE)

# Function to check and install CRAN packages
install_cran_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

# Function to check and install Bioconductor packages
install_bioc_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      BiocManager::install(pkg, update = FALSE, dependencies = TRUE)
    }
  }
}

# Function to check and install GitHub packages
install_github_packages <- function(packages) {
  for (pkg in names(packages)) {
    if (!require(pkg, character.only = TRUE)) {
      devtools::install_github(packages[[pkg]], upgrade = "never", dependencies = TRUE)
    }
  }
}

# List of CRAN packages
cran_packages <- c("shiny", "shinythemes", "shinyWidgets", "shinycssloaders", "shinydashboard", "shinyjs", "DT",  "ggplot2", "devtools", "Matrix", "reshape2", "zip", "tools", "grid", "filelock", "Seurat", BiocManager)


# List of Bioconductor packages
bioc_packages <- c("biomaRt", "rhdf5", "celda", "decontX", "SingleCellExperiment", "SummarizedExperiment")

# List of GitHub packages
github_packages <- list(
  "patchwork" = "thomasp85/patchwork",
  "SoupX" = "constantAmateur/SoupX"",
  "scCDC" = "ZJU-UoE-CCW-LAB/scCDC",
  "FastCAR" = "https://git.web.rug.nl/P278949/FastCAR"
)

# Install all packages
install_cran_packages(cran_packages)
install_bioc_packages(bioc_packages)
install_github_packages(github_packages)
```
## Start the App

To launch SCARR-Vis, follow one of these methods:

#### Option 1: Run Directly from GitHub
1. Open an R session in **RStudio**.
2. Execute the following lines of code:

```R
library(shiny)
shiny::runGitHub('SCARR-Vis', 'GudaLab')
```

#### Option 2: Download the source code from GitHub and run:
```
library(shiny)
runApp('/path/to/the/SCARR-Vis-master', launch.browser = TRUE)
```
Replace /path/to/the/SCARR-Vis-master with the actual path to the downloaded folder
## Usage

A detailed user manual is available under the "Manual" tab at: [https://www.gudalab-rtools.net/SCARR-Vis](https://www.gudalab-rtools.net/SCARR-Vis)

## Example Datasets

To ensure seamless analysis and reproducibility, **SCARR-Vis** includes one reference dataset [GSM7681687](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7681687). These datasets, sourced from **NCBI**, have been pre-tested with the tool and allow users to explore its functionalities and understand the analysis workflow effectively.


These datasets are ideal for:
- Demonstrating SCARR-Vis functionalities.
- Familiarizing users with the tool's analysis workflow.
- Testing the application.


## Tested Platforms

This application was tested on: Linux (Red Hat and Ubuntu) and Windows (10 and 11)

