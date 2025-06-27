suppressMessages(library(shiny))
suppressMessages(library(DT))
suppressMessages(library(shinyFiles))
suppressMessages(library(shinyjs))

fluidPage(
    useShinyjs(),
    tags$head(
        tags$title("CaMutQC: Cancer Somatic Mutation Quality Control and Filtration"),
        tags$style(HTML("
      body {
        background-color: #FFF5E6;
      }
      .title-container {
        text-align: center;
        margin-bottom: 30px;
      }
      .title-text {
        font-weight: bold;
        display: inline-block;
        vertical-align: middle;
        margin-left: 20px;
      }
      /* Run Filtration button styles */
      #run_filter {
        background-color: #FFB6C1 !important;  /* Light red */
        border-color: #FFB6C1 !important;
      }
      #run_filter:hover, #run_filter:active {
        background-color: #FF6347 !important;  /* Tomato red */
        border-color: #FF6347 !important;
      }
      /* Other button styles */
      .btn-default {
        background-color: #ADD8E6 !important;  /* Light blue */
        border-color: #ADD8E6 !important;
      }
      .btn-default:hover, .btn-default:active {
        background-color: #1E90FF !important;  /* Dodger blue */
        border-color: #1E90FF !important;
      }
    "))
    ),
    
    # Centered logo and title
    div(class = "title-container",
        tags$img(src = "CaMutQC_logo.png", height = "120px", 
                 style = "display: inline-block; vertical-align: middle;"),
        div(class = "title-text", h1("CaMutQC: Cancer Mutation Quality Control"))
    ),
    
    sidebarLayout(
        sidebarPanel(
            fileInput("vcf_file", "Upload VCF file(s)", multiple = FALSE, accept = ".vcf", placeholder = "example.vcf"),
            downloadButton("download_example", "Download VCF"),
            checkboxInput("PON_filter", "Filter based on PON file", value = FALSE),
            conditionalPanel(
              condition = "input.PON_filter == true",
              fileInput("PON_file", "Upload one PON file, txt or vcf", multiple = FALSE, accept = c(".vcf", ".txt"),
                        placeholder = "PON_example.txt"),
              selectInput("PON_file_type", "PON file type", choices = c("txt", "vcf")),
            ),
            selectInput("panel", "Filtration Panel", choices = c("Customized","MSKCC", "WES")),
            numericInput("VAF", "VAF cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
            numericInput("VAFratio", "VAF ratio cutoff", value = 0, min = 0, max = Inf, step = 0.5),
            numericInput("tumorDP", "Tumor Depth (DP) cutoff", value = 20, min = 0, max = Inf, step = 1),
            numericInput("normalDP", "Normal Depth (DP) cutoff", value = 10, min = 0, max = Inf, step = 1),
            numericInput("tumorAD", "Tumor Allele Depth (AD) cutoff", value = 5, min = 0, max = Inf, step = 1),
            numericInput("normalAD", "Normal Allele Depth (AD) cutoff", value = 500, min = 0, max = Inf, step = 1),
            numericInput("dbsnpCutoff", "Normal Allele Depth (AD) cutoff for dbSNP variants", value = 19, min = 0, max = Inf, step = 1),
            numericInput("nonCutoff", "Normal Allele Depth (AD) cutoff for non-dbSNP variants", value = 8, min = 0, max = Inf, step = 1),
            selectInput("sb", "Strand Bias method", choices = c('SOR', 'Fisher')),
            numericInput("sbscore", "Strand Bias cutoff", value = 3, min = 0, max = Inf, step = 1),
            numericInput("maxIndelLen", "Max length of indel will be filtered when next to SNV", value = 50, min = 0, max = Inf, step = 1),
            numericInput("minInterval", "Minimum interval between the indel will be filtered and its cloest SNV", value = 10, min = 0, max = Inf, step = 1),
            textInput('tagFILTER', 'tag used in the FILTER column', value = "PASS"),
            numericInput("dbVAF", "Database VAF cutoff", value = 0.01, min = 0, max = 1, step = 0.01),
            checkboxInput("ExAC", "ExAC", value = TRUE),
            checkboxInput("Genomesprojects1000", "Genomesprojects1000", value = TRUE),
            checkboxInput("gnomAD", "gnomAD", value = TRUE),
            checkboxInput("dbSNP", "dbSNP", value = FALSE),
            checkboxInput("ESP6500", "ESP6500", value = TRUE),
            checkboxInput("keepCOSMIC", "keepCOSMIC", value = TRUE),
            selectInput('keepType', 'keepType (kept during filtration)', choices = c('nonsynonymous', 'exonic', 'all')),
            checkboxInput("report", "Generate Filter Report", value = FALSE),
            # add a folder to save the report, only display when report is set to true
            conditionalPanel(
                condition = "input.report == true",
                shinyDirButton("dir_selector", "Select a folder to store report", "Browse", 
                               buttonClass = "btn-default",
                               icon = icon("folder-open")),
                verbatimTextOutput("dir_path")
            ),
            checkboxInput("TMB", "Calculate TMB (without filter)", value = FALSE),
            conditionalPanel(
                condition = "input.TMB == true",
                selectInput("assay", "TMB assay", choices = c('MSK-v3', 'MSK-v2', 'MSK-v1', 'FoundationOne', 
                                                              'Pan-Cancer Panel',
                                                              'Customized'))
            ),
            conditionalPanel(
              condition = "input.assay === 'Pan-Cancer Panel' || input.assay === 'Customized'",
              selectInput("mutType", "mutType (kept during TMB calculation)", choices = c('nonsynonymous', 'exonic')),
              
            ),
            conditionalPanel(
              condition = "input.assay === 'Customized'",
              textInput('genelist', 'Target genes during TMB calculation (separated by comma)', value = NULL),
              
            ),
            checkboxInput("bedFilter", "bedFilter", value = FALSE),
            conditionalPanel(
              condition = "input.bedFilter == true || input.assay === 'Customized' ",
              fileInput("bed_file", "Select a bed file", multiple = FALSE, 
                        accept = c(".bed", ".rds"), placeholder = "Select a bed file")
            ),
            conditionalPanel(
              condition = "input.bedFilter == true || input.assay === 'Customized'",
              checkboxInput("bedHeader", "bed has a header", value = FALSE)
            ),
            actionButton("run_filter", "Filter NOW"),
            downloadButton("download_maf", "Download Labeled MAF")
        ),
        mainPanel(
          verbatimTextOutput("error_message"),
            DTOutput("labeled_table"),
            verbatimTextOutput("tmb_value")
        )
    )
)

