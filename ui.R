# -----------------------------
# UI (navbarPage)
# -----------------------------
timeoutSeconds <- 36000

inactivity <- sprintf("function idleTimer() {
var t = setTimeout(logout, %s);
window.onmousemove = resetTimer; // catches mouse movements
window.onmousedown = resetTimer; // catches mouse movements
window.onclick = resetTimer;     // catches mouse clicks
window.onscroll = resetTimer;    // catches scrolling
window.onkeypress = resetTimer;  //catches keyboard actions

function logout() {
Shiny.setInputValue('timeOut', '%ss')
}

function resetTimer() {
clearTimeout(t);
t = setTimeout(logout, %s);  // time is in milliseconds (1000 is 1 second)
}
}
idleTimer();", timeoutSeconds*1000, timeoutSeconds, timeoutSeconds*1000)

ui <- navbarPage(
  theme = shinytheme("cerulean"),
  "",
  useShinyjs(),
  # head block so server's session$sendCustomMessage('force-reload') works
  tags$head(
    tags$script(HTML("
      Shiny.addCustomMessageHandler('force-reload', function(_) {
        window.location.reload(true);
      });
    "))
  ),
  
  # --- Introduction tab ---
  tabPanel(
    "SCARR-Vis",
    tags$script(inactivity),
    mainPanel(
      column(
        width = 10,
        br(),
        h3("SCARR-Vis - Single Cell Ambient RNA Removal and Visualization"),
        p("SCARR-Vis removes ambient RNA from 10x Genomics single-cell data using",
          strong("SoupX"), ", ", strong("DecontX"), ", ", strong("scCDC"), ", or ", strong("FastCAR"),
          "and lets you compare pre vs post cleanup across QC, clustering, UMAP, heatmaps, and per-cell tables."
        ),
        h4("Data upload"),
        p("Provide both the raw/droplet matrix and the filtered/cell matrix. You can upload as:"),
        tags$li("A zip folder containing matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz; or"),
        tags$li("HDF5 (.h5) from Cell Ranger / other tools."),
        h4("Estimate contamination"),
        p("Go to Step 2. Estimate Contamination and choose a method."),
        tags$li("SoupX: estimates sample-level ambient contamination (global ρ) from empty droplets; supports parameter tuning and cluster-aware count adjustment."),
        tags$li("DecontX: Infers per-cell Bayesian contamination; uses user-provided clusters or re-clusters; outputs decontaminated counts."),
        tags$li("scCDC: gene-specific contamination detection and correction (with doublet-aware adjustments)."),
        tags$li("FastCAR: fast ambient RNA correction using an empty-droplet UMI cutoff and gene-level contamination probability threshold, with optional ambient profiling to suggest a cutoff."),
        h4("Compare results"),
        tags$li("Check QC (Pre) vs QC (Post), Estimation, Cluster Counts (Pre vs Post), UMAP (Pre vs Post), Top Genes, and the Cells table."),
        h4("Outputs and Visualization"),
        tags$li("Download publication-quality plots in seven formats: JPG, TIFF, PDF, SVG, BMP, EPS, and PS. Summary tables are also generated in .csv format for easy visualization and download, cleaned Seurat/SCE objects, summaries, and a zip bundle."),
        hr(),
        h4("Tips"),
        tags$ul(
          tags$li("If SoupX errors on empty droplets, widen soupRange or use the automatic global soup profile fallback."),
          tags$li("If your dataset contains erythroid/testis cells, adjust 'Optional Markers'.")
        ),
        hr(),
        HTML("<h3> use SCARR-Vis online</h3>
          <p>SCARR-Vis is deployed at: <a href='https://www.gudalab-rtools.net/SCARR-Vis'>https://www.gudalab-rtools.net/SCARR-Vis</a></p>
          <hr>
          <h3> Launch SCARR-Vis using R and GitHub </h3>
          <p> SCARR-Vis were deposited under the GitHub repository: <a href='https://github.com/GudaLab/SCARR-Vis'>https://github.com/GudaLab/SCARR-Vis</a><br>
          Before running the app, users must have the following versions installed: R (>= 4.5.2), RStudio (>= 2025.09.2), Bioconductor (>= 3.22) and Shiny (>= 1.11.1) (Tested with this version).<br>
          Note: SCARR-Vis has been tested with these versions. If users are running an older version of R, they may encounter errors during package installation. Therefore, it is recommended to update R to the latest version first.<br>
          Once R is open in the command line or in RStudio, users should run the following command in R to install the shiny package.<br><br></p>
          <pre>install.packages('shiny')<br>library(shiny)</pre>
          <hr>
          <h3>Start the app</h3>
          Start the R session using RStudio and run these lines:<br><br>
          <pre>shiny::runGitHub('SCARR-Vis','GudaLab')</pre>
          or
        Alternatively, download the source code from GitHub and run the following command in the R session using RStudio:
          <pre>library(shiny)<br>runApp('/path/to/the/SCARR-Vis-master', launch.browser=TRUE)</pre>
          <hr>
          <h3>Usage</h3>
          <p>Please refer our Manual tab.</p>
          <hr>
          <h3> Developed and maintained by</h3>
          <p>SCARR-Vis was developed by Sankarasubramanian Jagadesan and Babu Guda. We share a passion for developing a user-friendly tool for biologists, particularly those who do not have access to bioinformaticians or programming expertise.</p>
          <hr>
        "),
        tags$head(
          tags$style(HTML("
            #view_count { color: #dc2626; font-weight: 700; }  /* red, bold number */
            .views-center { text-align: center; }
          "))
        ),
        tags$div(
          class = "views-center",
          tags$p("Total number of views: ", textOutput("view_count", inline = TRUE))
        ),
        hr()
      )
    )
  ),
  
  # --- Analysis tab (all existing tabs live here) ---
  tabPanel(
    "Analysis",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("1) Upload 10x Matrices"),
        helpText("Provide BOTH raw and filtered counts. Each may be a .zip of a 10x folder or a single .h5 file."),
        selectInput(
          "species_genome", "Species / Genome build",
          choices = c("Human (GRCh38)", "Mouse (GRCm39)"),
          selected = "Human (GRCh38)"
        ),
        radioButtons(
          "data_mode", "Data source",
          choices = c("Upload files" = "upload", "Use example data" = "example"),
          inline = TRUE, selected = "upload"
        ),
        conditionalPanel(
          condition = "input.data_mode == 'upload'",
          fileInput("raw",  "Raw (raw_feature_bc_matrix)",           accept = c(".zip", ".h5")),
          fileInput("filt", "Filtered (filtered_feature_bc_matrix)", accept = c(".zip", ".h5"))
        ),
        conditionalPanel(
          condition = "input.data_mode == 'example'",
          tags$div(
            class = "well",
            tags$b("Using built-in example data  (Organism: Human)"),
            tags$p("raw_feature_bc_matrix.h5 and filtered_feature_bc_matrix.h5 from GSM7681687")
          )
        ),
        checkboxInput("collapse_dups", "Collapse duplicate gene names by summing", TRUE),
        checkboxInput("use_gene_symbols", "Convert Ensembl IDs to gene symbols (if possible)", TRUE),
        div(class = "small-note", "If features are already symbols, they are kept."),
        hr(),
        h4("2) Estimate Contamination"),
        radioButtons(
          "method", "Decontamination method",
          choices = c("SoupX", "DecontX", "scCDC", "FastCAR"),
          inline = TRUE, selected = "SoupX"
        ),
        
        # SoupX controls
        conditionalPanel(
          condition = "input.method == 'SoupX'",
          numericInput("min_cells", "Min cells to keep (pre-filter)", value = 3, min = 1, step = 1),
          checkboxInput("do_auto", "Run autoEstCont (recommended)", TRUE),
          conditionalPanel(
            condition = "!input.do_auto",
            sliderInput(
              "manual_rho", "Manual contamination fraction (ρ)",
              min = 0, max = 1, value = 0.05, step = 0.005
            ),
            div(class = "small-note", "Applied uniformly to all cells when autoEstCont is off.")
          ),
          checkboxInput("do_plot", "Show estimation plots", TRUE),
          hr(),
          h5("Soup profile options"),
          sliderInput(
            "soupRange", "soupRange (UMIs in empty droplets)",
            min = 0, max = 2000, value = c(0, 100), step = 10
          ),
          checkboxInput("keepDroplets", "keepDroplets (uses memory)", FALSE),
          h4("Optional Markers (non-expressed in most cells)"),
          textInput("nonexp", "Comma-separated genes (eg: HBB,HBA2,IGKC)", value = "")
        ),
        
        # DecontX controls
        conditionalPanel(
          condition = "input.method == 'DecontX'",
          numericInput("min_cells", "Min cells to keep (pre-filter)", value = 3, min = 1, step = 1),
          checkboxInput("decontx_use_clusters", "Use Seurat clusters as priors (recommended)", TRUE),
          numericInput("decontx_maxiter", "maxIter", value = 500, min = 50, step = 50),
          textInput("decontx_delta", "delta (two numbers, e.g., 10,10)", value = "10,10"),
          numericInput("decontx_convergence", "convergence (tolerance)", value = 0.001, min = 1e-6, step = 1e-4),
          numericInput("decontx_iterLogLik", "iterLogLik (logLik interval)", value = 10, min = 1, step = 1),
          numericInput("decontx_varGenes", "varGenes (features to use)", value = 5000, min = 100, step = 100),
          checkboxInput("decontx_estimateDelta", "estimateDelta", TRUE)
        ),
        
        # scCDC controls (dropdowns)
        conditionalPanel(
          condition = "input.method == 'scCDC'",
          h5("scCDC parameters"),
          selectInput(
            "sccdc_restriction", "restriction_factor",
            choices = c(0.2, 0.3, 0.4, 0.5, 0.6), selected = 0.5
          ),
          selectInput(
            "sccdc_min_cell", "min.cell",
            choices = c(50, 100, 200, 300, 500), selected = 100
          ),
          selectInput(
            "sccdc_percent_cutoff", "percent.cutoff",
            choices = c(0.1, 0.2, 0.3, 0.4, 0.5), selected = 0.2
          )
        ),
        
        # FastCAR controls
        conditionalPanel(
          condition = "input.method == 'FastCAR'",
          h5("FastCAR parameters"),
          numericInput("fastcar_empty_cutoff", "Empty droplet UMI cutoff", value = 100, min = 10, max = 5000, step = 10),
          sliderInput("fastcar_contam_cutoff", "Contamination chance cutoff", min = 0, max = 0.5, value = 0.05, step = 0.005),
          checkboxInput("fastcar_do_profile", "Profile ambient RNA to suggest cutoff", value = TRUE),
          conditionalPanel(
            condition = "input.fastcar_do_profile",
            fluidRow(
              numericInput("fastcar_profile_start", "Profile start (UMI)", value = 10,  min = 1, max = 2000, step = 1),
              numericInput("fastcar_profile_stop",  "Profile stop (UMI)",  value = 500, min = 50, max = 10000, step = 10),
              numericInput("fastcar_profile_by",    "Profile step",        value = 10,  min = 1, max = 100,   step = 1)
            ),
            checkboxInput("fastcar_use_recommended", "Use recommended cutoff from profile", value = TRUE)
          )
        ),
        hr(),
        actionBttn("run", "Run", style = "unite", color = "primary", icon = icon("download")),
        hr()
      ),
      
      mainPanel(
        width = 9,
        tabsetPanel(
          id = "tabs",
          
          tabPanel(
            "Status",
            br(),
            conditionalPanel(
              condition = "input.run == 0",
              box(
                h3("Instructions for Uploading Sample Files"),
                HTML("<ol>
               <li><strong>H5 Files (Cell Ranger Output)</strong></li>
               <ul>
               <li>Cell Ranger file: raw_feature_bc_matrix.h5, filtered_feature_bc_matrix.h5.</li>
               </ul>
               
               <li><strong>Cell Ranger Matrix Files</strong></li>
               <ul>
               <li>Cell Ranger files: filtered_feature_bc_matrix (matrix.mtx.gz, feature.tsv.gz, barcode.tsv.gz) and filtered_feature_bc_matrix (matrix.mtx.gz, feature.tsv.gz, barcode.tsv.gz).</li>
               <li>Compress the folder to .zip format.</li>
               </ul>
               <br>
               <img src='images/folder_image.jpg' width='600' height='400' alt=''/>
               </ol>
                "),
                h4("Clarification example file format"),
                h5("Users can download this example dataset to better understand the required structure. Following this reference will help ensure that your files are correctly prepared and fully compatible with our tool"),
                tags$b("H5 File (Cell Ranger Output)"),
                br(),
                a(
                  href = "https://www.gudalab-rtools.net/SCARR-Vis/example/GSM7681687/raw_feature_bc_matrix.h5",
                  "Raw_feature_bc_matrix H5 File", style = "color:red;", download = NA, target = "_blank"
                ),
                br(),
                a(
                  href = "https://www.gudalab-rtools.net/SCARR-Vis/example/GSM7681687/filtered_feature_bc_matrix.h5",
                  "Filtered_feature_bc_matrix H5 File", style = "color:red;", download = NA, target = "_blank"
                ),
                br(),
                tags$b("or"),
                br(),
                tags$b("Cell Ranger Matrix Files"),
                br(),
                a(
                  href = "https://www.gudalab-rtools.net/SCARR-Vis/example/GSM7681687/raw_feature_bc_matrix.zip",
                  "Raw_feature_bc_matrix Matrix Files", style = "color:red;", download = NA, target = "_blank"
                ),
                br(),
                a(
                  href = "https://www.gudalab-rtools.net/SCARR-Vis/example/GSM7681687/filtered_feature_bc_matrix.zip",
                  "Filtered_feature_bc_matrix Matrix Files", style = "color:red;", download = NA, target = "_blank"
                ),
                br()
              )
            ),
            verbatimTextOutput("log") %>% withSpinner(type = 4),
            br(),
            conditionalPanel(
              condition = "input.run > 0",
              downloadBttn("dl_adjusted_rds", "Download Seurat (RDS)"),
              br(),
              br(),
              uiOutput("dl_adjusted_ui")
            )
          ),
          
          tabPanel(
            "QC (Pre)",
            fluidRow(
              column(
                6,
                h3("Detected genes per cell"),
                plotOutput("pre_nFeature") %>% withSpinner(type = 4),
                div(class = "dl-row", actionBttn("open_dl_pre_nFeature", "Download plot", style = "unite", color = "primary", icon = icon("download")))
              ),
              column(
                6,
                h3("UMIs per cell"),
                plotOutput("pre_nCount") %>% withSpinner(type = 4),
                div(class = "dl-row", actionBttn("open_dl_pre_nCount", "Download plot", style = "unite", color = "primary", icon = icon("download")))
              )
            ),
            h3("Mitocondrial percentage"),
            plotOutput("pre_pctMT") %>% withSpinner(type = 4),
            div(class = "dl-row", actionBttn("open_dl_pre_pctMT", "Download plot", style = "unite", color = "primary", icon = icon("download")))
          ),
      
          tabPanel(
            "Estimation",
            
            ## ---- SoupX & DecontX: rho density + rho vs nUMIs ----
            conditionalPanel(
              condition = "input.method == 'SoupX' || input.method == 'DecontX'",
              h3("rho density"),
              plotOutput("est_density") %>% withSpinner(type = 4),
              div(
                class = "dl-row",
                actionBttn(
                  "open_dl_est_density", "Download plot",
                  style = "unite", color = "primary", icon = icon("download")
                )
              ),
              br()
            ),
            conditionalPanel(
              condition = "input.method == 'SoupX' || input.method == 'DecontX'",
              fluidRow(
                column(
                  6,
                  h3("rho vs nUMIs"),
                  plotOutput("rho_vs_counts") %>% withSpinner(type = 4),
                  div(
                    class = "dl-row",
                    actionBttn(
                      "open_dl_rho_vs_counts", "Download plot",
                      style = "unite", color = "primary", icon = icon("download")
                    )
                  )
                ),
                column(
                  6,
                  ## This UI only shows anything when method == "DecontX"
                  uiOutput("decontx_contam_panel")
                )
              )
            ),
            
            ## ---- scCDC: diagnostics + top GCGs + PDF ----
            conditionalPanel(
              condition = "input.method == 'scCDC'",
              fluidRow(
                column(
                  12,
                  uiOutput("scCDC_diag_panel"),
                  uiOutput("sccdc_panel"),
                  uiOutput("scCDC_pdf_panel")
                )
              )
            ),
            
            ## ---- SoupX-only: autoEstCont diagnostic ----
            conditionalPanel(
              condition = "input.method == 'SoupX'",
              h3("Auto estimation contamination diagnostic"),
              plotOutput("auto_plot") %>% withSpinner(type = 4),
              div(
                class = "dl-row",
                actionBttn(
                  "open_dl_auto_plot", "Download plot",
                  style = "unite", color = "primary", icon = icon("download")
                )
              )
            ),
            
            ## ---- FastCAR-only: ambient profile + removed-reads histogram ----
            conditionalPanel(
              condition = "input.method == 'FastCAR'",
              fluidRow(
                column(
                  6,
                  h3("Ambient profile plot"),
                  plotOutput("fastcar_profile_plot", height="600px") %>% withSpinner(type = 4),
                  div(
                    class = "dl-row",
                    actionBttn(
                      "open_dl_fastcar_profile", "Download plot",
                      style = "unite", color = "primary", icon = icon("download")
                    )
                  )
                ),
                column(
                  6,
                  h3("UMIs removed per cell"),
                  plotOutput("fastcar_removed_hist") %>% withSpinner(type = 4),
                  div(
                    class = "dl-row",
                    actionBttn(
                      "open_dl_fastcar_removed", "Download plot",
                      style = "unite", color = "primary", icon = icon("download")
                    )
                  )
                )
              )
            )
          ),
          
          tabPanel(
            "QC (Post)",
            fluidRow(
              column(
                6,
                h3("Detected genes per cell"),
                plotOutput("post_nFeature") %>% withSpinner(type = 4),
                div(class = "dl-row", actionBttn("open_dl_post_nFeature", "Download plot", style = "unite", color = "primary", icon = icon("download")))
              ),
              column(
                6,
                h3("UMIs cell"),
                plotOutput("post_nCount") %>% withSpinner(type = 4),
                div(class = "dl-row", actionBttn("open_dl_post_nCount", "Download plot", style = "unite", color = "primary", icon = icon("download")))
              )
            ),
            h3("Mitocondrial percentage"),
            plotOutput("post_pctMT") %>% withSpinner(type = 4),
            div(class = "dl-row", actionBttn("open_dl_post_pctMT", "Download plot", style = "unite", color = "primary", icon = icon("download")))
          ),
          
          tabPanel(
            "Cluster counts (Pre vs post)",
            fluidRow(column(12, downloadBttn("dl_cluster_counts_csv", "Download Cluster Counts (CSV)"))),
            br(),
            h3("Cell counts table before and after ambinent RNA removal"),
            DTOutput("cluster_counts_dt") %>% withSpinner(type = 4),
            br(),
            h3("Cell counts plot before and after ambinent RNA removal"),
            plotOutput("cluster_counts_bar", height = 400) %>% withSpinner(type = 4),
            div(class = "dl-row", actionBttn("open_dl_cluster_counts_bar", "Download plot", style = "unite", color = "primary", icon = icon("download")))
          ),
          
          tabPanel(
            "Cells Table",
            fluidRow(column(12, downloadBttn("dl_cells_csv", "Download Cells Table (CSV)"))),
            br(),
            h3("Cell stats table"),
            DTOutput("cells_dt") %>% withSpinner(type = 4)
          ),
          
          tabPanel(
            "Top Genes",
            fluidRow(column(12, downloadBttn("dl_topgenes_csv", "Download Top Genes (CSV)"))),
            br(),
            h3("Top genes"),
            DTOutput("top_genes") %>% withSpinner(type = 4)
          ),
          
          tabPanel(
            "UMAP/TSNE (Pre vs Post)",
            fluidRow(
              column(
                6,
                h3("UMAP plot before ambient RNA removal"),
                plotOutput("umap_pre", height = 450) %>% withSpinner(type = 4),
                div(class = "dl-row", actionBttn("open_dl_umap_pre", "Download plot", style = "unite", color = "primary", icon = icon("download")))
              ),
              column(
                6,
                h3("UMAP plot after ambient RNA removal"),
                plotOutput("umap_post", height = 450) %>% withSpinner(type = 4),
                div(class = "dl-row", actionBttn("open_dl_umap_post", "Download plot", style = "unite", color = "primary", icon = icon("download")))
              ),
              column(
                6,
                h3("TSNE plot before ambient RNA removal"),
                plotOutput("tsne_pre", height = 450) %>% withSpinner(type = 4),
                div(class = "dl-row", actionBttn("open_dl_tsne_pre", "Download plot", style = "unite", color = "primary", icon = icon("download")))
              ),
              column(
                6,
                h3("TSNE plot after ambient RNA removal"),
                plotOutput("tsne_post", height = 450) %>% withSpinner(type = 4),
                div(class = "dl-row", actionBttn("open_dl_tsne_post", "Download plot", style = "unite", color = "primary", icon = icon("download")))
              )
            )
          ),
          
          tabPanel(
            "Heatmap (Pre vs Post)",
            fluidRow(
              column(4, radioButtons("heatmap_source", "Show heatmap for:", c("Pre" = "pre", "Post" = "post"), inline = TRUE)),
              column(4, numericInput("heatmap_topn", "Top variable genes", 50, min = 10, step = 10)),
              column(4, actionBttn("open_heatmap_download", "Download heatmap", style = "unite", color = "primary", icon = icon("download")))
            ),
            h3("Heatmap before or after ambient RNA removal"),
            plotOutput("heatmap_plot", height = 600) %>% withSpinner(type = 4)
          ),
          
          tabPanel(
            "Feature Plot (Pre vs Post)",
            fluidRow(
              column(8, textInput("featplot_genes", "Gene symbols (comma-separated)", value = "")),
              column(4, actionBttn("open_dl_featplot", "Download plot", style = "unite", color = "primary", icon = icon("download")))
            ),
            br(),
            uiOutput("featplot_warning"),
            h3("Feature plot before ambient RNA removal"),
            plotOutput("featplot_combined", height = 600) %>% withSpinner(type = 4)
          ),
          
          tabPanel(
            "Session Info",
            tags$h4("R Session info"),
            downloadBttn("download_sess", "Download session-info.txt"),
            br(),
            withSpinner(verbatimTextOutput("sess"), type=4)
          )
        )
      )
    )
  ),
  
  tabPanel(
    "Manual",
    fluidRow(
      column(
        8, offset = 2,
        br(),
        h2("SCARR-Vis: Manual & Parameter Reference"),
        p(
          "SCARR-Vis is an R/Shiny application for interactive assessment and correction of ambient RNA contamination in single-cell and single nucleus RNA-seq data. The interface follows a typical workflow: upload 10x matrices, choose a decontamination method, inspect pre- and post-correction QC, and explore clustering and gene expression."
           ),
        div(
          style = "text-align:center; margin-bottom:20px;",
          tags$img(src = "images/folder_image.jpg", class = "img-responsive", style = "max-width:100%;"),
        ),
        h3("Overview of the Example Dataset (GSM7681687)"),
        p(
          "Throughout this manual we use the single-cell RNA-seq sample ",
          strong("GSM7681687"), " from NCBI GEO. FASTQ files were processed with ",
          strong("Cell Ranger"), " to generate both ", strong("raw"),
          " and ", strong("filtered feature-barcode matrices"), ". ",
          "Because raw matrices are typically not submitted to GEO, SCARR-Vis bundles ",
          "this dataset as ", strong("example data"), " so users can play with the ",
          "full pipeline, including ambient RNA estimation that relies on the raw matrix."
        ),
        p(
          "To quickly test SCARR-Vis, select ",
          strong("“Use example data”"),
          " in the upload panel and run the pipeline with the default parameters."
        ),
        
        hr(),
        h3("Upload 10x Matrices"),
        p(
          "In step 1, SCARR-Vis expects both a raw and a filtered 10x feature-barcode matrix. ",
          "You can either use the bundled GSM7681687 example or upload your own matrices."
        ),
        tags$ul(
          tags$li(strong("Species / genome build:"), " choose the appropriate reference."),
          tags$li(strong("Raw matrix:"), " a 10x HDF5 file ."),
          tags$li(strong("Filtered matrix:"), " a zipped directory contain barcode, matrix feature file show in above image."),
          tags$li("Optional: collapse duplicate gene names; convert Ensembl IDs to gene symbols.")
        ),
        div(
          style = "text-align:center; margin-bottom:20px;",
          tags$img(src = "images/1.jpg", class = "img-responsive", style = "max-width:100%;"),
          tags$small(
            em("Figure 1. Upload panel and estimation method selector. ",
               "Sub-panels show: (a) upload matrices, ",
               "(b) SoupX, (c) DecontX, (d) scCDC, (e) FastCAR parameter panels.")
          )
        ),
        p(
          strong("Note:"),
          " when starting from FASTQ files, you must first run Cell Ranger (or a similar ",
          "pipeline) to generate the raw and filtered matrices required by SCARR-Vis."
        ),
        
        hr(),
        h3("Choose Estimation Method and Parameters"),
        p(
          "In step 2 (Estimate Contamination), SCARR-Vis provides adapters for four methods:"
        ),
        tags$ul(
          tags$li(strong("SoupX:"), " models background 'soup' RNA and adjusts counts."),
          tags$li(strong("DecontX:"), "from celda infers cell-specific contamination fractions."),
          tags$li(strong("scCDC:"), " identifies contamination-causing genes (GCGs) and optionally corrects them."),
          tags$li(strong("FastCAR:"), " profiles ambient RNA using empty droplets and estimates per-cell contamination.")
        ),
        p(
          "Each method has its own parameter panel (see Figure 1b–1e). ",
          "After setting parameters, click ", strong("Run"), " to start estimation. ",
          "Progress and status messages are shown in the log."
        ),
        
        hr(),
        h3("1. QC (Pre)"),
        p(
          "The ", strong("QC (Pre)"), " tab summarizes per-cell metrics before any correction:"
        ),
        tags$ul(
          tags$li("Detected genes per cell."),
          tags$li("UMIs per cell."),
          tags$li("Mitochondrial percentage (species-aware MT gene pattern).")
        ),
        div(
          style = "text-align:center; margin-bottom:20px;",
          tags$img(src = "images/2.jpg", class = "img-responsive", style = "max-width:100%;"),
          tags$small(
            em("Figure 2. QC (Pre) for GSM7681687: (a) detected genes per cell, ",
               "(b) UMIs per cell, (c) mitochondrial %. ")
          )
        ),
        
        hr(),
        h3("2. Estimation Diagnostics"),
        p(
          "The ", strong("Estimation"), " tab displays method-specific diagnostic plots."
        ),
        
        h4("2.1 SoupX"),
        p(
          "SoupX output includes the distribution of contamination fractions (ρ), ",
          "their relationship with UMI counts, and other model diagnostics."
        ),
        div(
          style = "text-align:center; margin-bottom:20px;",
          tags$img(src = "images/3.1.jpg", class = "img-responsive", style = "max-width:100%;"),
          tags$small(
            em("Figure 3.1. SoupX diagnostics: (a) ρ density, (b) ρ vs nUMIs (c) auto estimation contamination diagnostic.")
          )
        ),
        
        h4("2.2 DecontX"),
        p(
          "DecontX produces a contamination density plot, ρ vs nUMIs, and a UMAP colored ",
          "by estimated contamination."
        ),
        div(
          style = "text-align:center; margin-bottom:20px;",
          tags$img(src = "images/3.2.jpg", class = "img-responsive", style = "max-width:100%;"),
          tags$small(
            em("Figure 3.2. DecontX diagnostics: contamination density, ρ vs nUMIs, ",
               "and contamination on UMAP.")
          )
        ),
        
        h4("2.3 scCDC"),
        p(
          "scCDC focuses on genes driving contamination. SCARR-Vis shows UMIs per cell ",
          "post vs pre, mean counts of top GCGs, and entropy vs mean expression to highlight ",
          "putative contamination genes."
        ),
        div(
          style = "text-align:center; margin-bottom:20px;",
          tags$img(src = "images/3.3.jpg", class = "img-responsive", style = "max-width:100%;"),
          tags$small(
            em("Figure 3.3. scCDC diagnostics: (a) UMIs per cell (post vs pre), ",
               "(b) top GCGs (pre vs post), (c) entropy vs mean expression.")
          )
        ),
        
        h4("2.4 FastCAR"),
        p(
          "FastCAR scans a grid of empty-droplet UMI cutoffs to profile ambient RNA and ",
          "identifies an appropriate threshold for empty droplets and contamination. ",
          "SCARR-Vis displays the number of empty droplets and genes in ambient RNA at each ",
          "cutoff, as well as the distribution of UMIs removed per cell."
        ),
        div(
          style = "text-align:center; margin-bottom:20px;",
          tags$img(src = "images/3.4.jpg", class = "img-responsive", style = "max-width:100%;"),
          tags$small(
            em("Figure 3.4. FastCAR diagnostics: (a) empty-droplet profile, ",
               "(b) reads removed per cell.")
          )
        ),
        
        hr(),
        h3("3. QC (Post)"),
        p(
          "After decontamination, the ", strong("QC (Post)"), " tab repeats the same metrics ",
          "as QC (Pre) but on the corrected counts, allowing a direct comparison."
        ),
        div(
          style = "text-align:center; margin-bottom:20px;",
          tags$img(src = "images/4.jpg", class = "img-responsive", style = "max-width:100%;"),
          tags$small(
            em("Figure 4. QC (Post) histograms for GSM7681687: detected genes, UMIs, and mitochondrial %. ")
          )
        ),
        
        hr(),
        h3("4. Cluster counts (Pre vs Post)"),
        p(
          "The ", strong("Cluster counts (Pre vs Post)"), " tab compares the number of cells ",
          "in each cluster before and after correction, making it easy to see whether ",
          "ambient RNA disproportionately affected particular clusters."
        ),
        div(
          style = "text-align:center; margin-bottom:20px;",
          tags$img(src = "images/table_5.jpg", class = "img-responsive", style = "max-width:100%;"),
          tags$small(
            em("Figure 5. Cluster counts per cluster (pre vs post).")
          )
        ),
        
        h4("5. Cell table"),
        p(
          "The ", strong("Cells Table"), " tab lists per-cell metrics such as cluster ID, ",
          "UMIs, QC statistics, and contamination estimates. Users can sort and filter ",
          "rows for detailed inspection."
        ),
        div(
          style = "text-align:center; margin-bottom:20px;",
          tags$img(src = "images/able_6.jpg", class = "img-responsive", style = "max-width:100%;"),
          tags$small(
            em("Figure 6. Cell-level summary table.")
          )
        ),
        
        h4("6 Top genes"),
        p(
          "The ", strong("Top Genes"), " tab highlights the genes most affected by ambient RNA ",
          "and compares their pre- and post-correction expression."
        ),
        div(
          style = "text-align:center; margin-bottom:20px;",
          tags$img(src = "images/table_7.jpg", class = "img-responsive", style = "max-width:100%;"),
          tags$small(
            em("Figure 7. Top genes affected by ambient RNA (pre vs post).")
          )
        ),
        
        hr(),
        h3("7. UMAP / tSNE (Pre vs Post)"),
        p(
          "Global structure is visualized in the ", strong("UMAP/TSNE (Pre vs Post)"),
          " tab. SCARR-Vis recomputes embeddings on both the original and corrected counts using Seurat."
        ),
        tags$ul(
          tags$li("UMAP (Pre) and UMAP (Post)."),
          tags$li("tSNE (Pre) and tSNE (Post).")
        ),
        div(
          style = "text-align:center; margin-bottom:20px;",
          tags$img(src = "images/8.jpg", class = "img-responsive", style = "max-width:100%;"),
          tags$small(
            em("Figure 8. UMAP and tSNE embeddings for GSM7681687 before and after correction.")
          )
        ),
        
        hr(),
        h3("8. Heatmap (Pre vs Post)"),
        p(
          "The ", strong("Heatmap"), " tab visualizes expression of top variable genes across clusters. ",
          "Users can switch between pre- and post-correction matrices to see how ambient removal ",
          "changes gene-level patterns."
        ),
        div(
          style = "text-align:center; margin-bottom:20px;",
          tags$img(src = "images/9.jpg", class = "img-responsive", style = "max-width:100%;"),
          tags$small(
            em("Figure 9. Heatmaps of top variable genes for pre- and post-correction data.")
          )
        ),
        
        hr(),
        h3("9. Feature plots (Pre vs Post)"),
        p(
          "The ", strong("Feature Plot"), " tab displays per-gene expression over UMAP/tSNE. ",
          "Users provide one or more comma-separated gene symbols (e.g. ", code("FTH1"), "), and ",
          "SCARR-Vis plots paired pre- and post-correction feature maps."
        ),
        div(
          style = "text-align:center; margin-bottom:20px;",
          tags$img(src = "images/10.jpg", class = "img-responsive", style = "max-width:100%;"),
          tags$small(
            em("Figure 10. Example feature plots for FTH1 (Pre and Post).")
          )
        ),
        
        hr(),
        h3("10. Reproducibility and Session Info"),
        p(
          "SCARR-Vis provides a reproducibility summary and full R session information. ",
          "The reproducibility table records the selected method, key parameter values, ",
          "and dataset-level statistics. The Session Info tab shows R version, platform, ",
          "and package versions. These should be included in reports or manuscripts so ",
          "that analyses can be fully reproduced."
        ),
        div(
          style = "text-align:center; margin-bottom:20px;",
          tags$img(src = "images/table_1.jpg", class = "img-responsive", style = "max-width:100%;"),
          tags$small(
            em("Figure 11. Reproducibility summary table with selected parameters.")
          )
        ),
        
        hr(),
        p(
          strong("Pipeline summary: "),
          "The Seurat-based processing (normalization, variable feature selection, PCA, ",
          "neighbors, clustering, and UMAP) is run twice: first on the uploaded filtered ",
          "counts (Pre) and again on the decontaminated counts (Post). The downstream ",
          "visualization tabs always reflect this paired design."
        ),
        
        # -------------------------------------------------------------------
        # EXISTING PARAMETER REFERENCE TABLES (unchanged from your code)
        # -------------------------------------------------------------------
        hr(),
        h3("Parameter reference"),
        
        h4("General"),
        tags$table(
          class = "table table-condensed",
          tags$thead(tags$tr(tags$th("Name"), tags$th("Default"), tags$th("Min"), tags$th("Max"), tags$th("Notes"))),
          tags$tbody(
            tags$tr(
              tags$td(code("min_cells")), tags$td("3"), tags$td("1"), tags$td("—"),
              tags$td("Minimum cells per gene to retain when creating Seurat objects.")
            )
          )
        ),
        
        h4("SoupX"),
        tags$table(
          class = "table table-condensed",
          tags$thead(tags$tr(tags$th("Name"), tags$th("Default"), tags$th("Min"), tags$th("Max"), tags$th("Notes"))),
          tags$tbody(
            tags$tr(tags$td(code("do_auto")), tags$td("TRUE"), tags$td("FALSE"), tags$td("TRUE"),
                    tags$td("If TRUE, uses autoEstCont to estimate ρ per cell.")),
            tags$tr(tags$td(code("manual_rho")), tags$td("0.05"), tags$td("0"), tags$td("1"),
                    tags$td("Used only if do_auto = FALSE; uniform ρ.")),
            tags$tr(tags$td(code("soupRange")), tags$td("c(0, 100)"), tags$td("0"), tags$td("2000"),
                    tags$td("UMI range of empty droplets to build soup profile.")),
            tags$tr(tags$td(code("keepDroplets")), tags$td("FALSE"), tags$td("FALSE"), tags$td("TRUE"),
                    tags$td("Keeps droplet table in memory; uses more RAM."))
          )
        ),
        
        h4("DecontX"),
        tags$table(
          class = "table table-condensed",
          tags$thead(tags$tr(tags$th("Name"), tags$th("Default"), tags$th("Min"), tags$th("Max"), tags$th("Notes"))),
          tags$tbody(
            tags$tr(tags$td(code("decontx_use_clusters")), tags$td("TRUE"), tags$td("FALSE"), tags$td("TRUE"),
                    tags$td("If TRUE, uses Seurat clusters as priors.")),
            tags$tr(tags$td(code("decontx_maxiter (maxIter)")), tags$td("500"), tags$td("50"), tags$td("10000"),
                    tags$td("Maximum EM iterations.")),
            tags$tr(tags$td(code("decontx_delta")), tags$td("10,10"), tags$td(">0,>0"), tags$td("—"),
                    tags$td("Dirichlet prior hyperparameters as two numbers.")),
            tags$tr(tags$td(code("decontx_estimateDelta")), tags$td("TRUE"), tags$td("FALSE"), tags$td("TRUE"),
                    tags$td("Estimate delta during fitting.")),
            tags$tr(tags$td(code("decontx_convergence")), tags$td("0.001"), tags$td("1e-6"), tags$td("0.1"),
                    tags$td("EM tolerance for convergence.")),
            tags$tr(tags$td(code("decontx_iterLogLik")), tags$td("10"), tags$td("1"), tags$td("1000"),
                    tags$td("Iterations between log-likelihood checks.")),
            tags$tr(tags$td(code("decontx_varGenes")), tags$td("5000"), tags$td("100"), tags$td("30000"),
                    tags$td("Number of variable genes used by decontX."))
          )
        ),
        
        h4("scCDC"),
        tags$table(
          class = "table table-condensed",
          tags$thead(tags$tr(tags$th("Name"), tags$th("Default / Option"), tags$th("Notes"))),
          tags$tbody(
            tags$tr(tags$td(code("restriction_factor")), tags$td("0.5 (dropdown)"),
                    tags$td("Controls aggressiveness of GCG detection.")),
            tags$tr(tags$td(code("min.cell")), tags$td("100 (dropdown)"),
                    tags$td("Minimum cells per gene for estimation.")),
            tags$tr(tags$td(code("percent.cutoff")), tags$td("0.2 (dropdown)"),
                    tags$td("Threshold for ambient fraction filtering."))
          )
        ),
        
        h4("FastCAR"),
        tags$table(
          class = "table table-condensed",
          tags$thead(tags$tr(tags$th("Name"), tags$th("Default"), tags$th("Min"), tags$th("Max"), tags$th("Notes"))),
          tags$tbody(
            tags$tr(tags$td(code("fastcar_empty_cutoff")), tags$td("100"), tags$td("10"), tags$td("5000"),
                    tags$td("Maximum UMIs to call a droplet 'empty'. Higher values can over-correct lowly expressed genes.")),
            tags$tr(tags$td(code("fastcar_contam_cutoff")), tags$td("0.05"), tags$td("0"), tags$td("0.5"),
                    tags$td("Contamination chance cutoff used for background detection; lower is more conservative.")),
            tags$tr(tags$td(code("fastcar_do_profile")), tags$td("TRUE"), tags$td("FALSE"), tags$td("TRUE"),
                    tags$td("If TRUE, runs describe.ambient.RNA.sequence to profile ambient RNA over a grid of empty-droplet cutoffs.")),
            tags$tr(tags$td(code("fastcar_profile_start")), tags$td("10"), tags$td("1"), tags$td("2000"),
                    tags$td("Lower bound of UMI cutoff grid for ambient profiling.")),
            tags$tr(tags$td(code("fastcar_profile_stop")), tags$td("500"), tags$td("50"), tags$td("10000"),
                    tags$td("Upper bound of UMI cutoff grid for ambient profiling.")),
            tags$tr(tags$td(code("fastcar_profile_by")), tags$td("10"), tags$td("1"), tags$td("100"),
                    tags$td("Step size of UMI cutoff grid for ambient profiling.")),
            tags$tr(tags$td(code("fastcar_use_recommended")), tags$td("TRUE"), tags$td("FALSE"), tags$td("TRUE"),
                    tags$td("If TRUE, uses FastCAR's recommended empty-droplet cutoff based on the ambient profile."))
          )
        ),
        
        h4("Seurat processing (defaults used in app)"),
        tags$table(
          class = "table table-condensed",
          tags$thead(tags$tr(
            tags$th("Step"), tags$th("Key parameters (value)"), tags$th("Notes")
          )),
          tags$tbody(
            tags$tr(
              tags$td("Mito %"),
              tags$td(HTML(paste0(
                code("PercentageFeatureSet"), " pattern = ",
                code("^MT- (human) / ^mt- (mouse)")
              ))),
              tags$td("Species-aware mitochondrial regex.")
            ),
            tags$tr(
              tags$td("NormalizeData"),
              tags$td(HTML(paste(
                code('normalization.method=\"LogNormalize\"'), ",",
                code("scale.factor=10000")
              ))),
              tags$td("Standard log-normalization.")
            ),
            tags$tr(
              tags$td("FindVariableFeatures"),
              tags$td(HTML(paste(
                code('selection.method=\"vst\"'), ",",
                code("nfeatures=2000")
              ))),
              tags$td("Top 2,000 HVGs (Seurat default).")
            ),
            tags$tr(
              tags$td("ScaleData"),
              tags$td(HTML(paste(
                code("center=TRUE"), ",", code("scale=TRUE"), ",", code("verbose=FALSE")
              ))),
              tags$td("Centers and scales features before PCA.")
            ),
            tags$tr(
              tags$td("RunPCA"),
              tags$td(HTML(paste(
                code("features=VariableFeatures(object)"), ",",
                code("npcs=30"), ",", code("verbose=FALSE")
              ))),
              tags$td("PCA on HVGs; 30 PCs kept.")
            ),
            tags$tr(
              tags$td("FindNeighbors"),
              tags$td(HTML(paste(
                code("reduction='pca'"), ",",
                code("dims=1:20"), ",",
                code("k.param=20")
              ))),
              tags$td("SNN graph on first 20 PCs; k=20.")
            ),
            tags$tr(
              tags$td("FindClusters"),
              tags$td(HTML(paste(
                code("resolution=0.5"), ",",
                code("algorithm=1")
              ))),
              tags$td("Louvain (algorithm 1) at res=0.5.")
            ),
            tags$tr(
              tags$td("RunUMAP"),
              tags$td(HTML(paste(
                code("reduction='pca'"), ",",
                code("dims=1:20"), ",",
                code("n.neighbors=30"), ",",
                code("min.dist=0.3"), ",",
                code('umap.method=\"uwot\"'), ",",
                code('metric=\"cosine\"')
              ))),
              tags$td("UMAP via uwot; first 20 PCs.")
            )
          )
        ),
        
        hr(),
        p("See the Estimation, QC, and visualization tabs for diagnostics and plots after you run the pipeline.")
      )
    )
  )
  
)
