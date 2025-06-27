library(shiny)
library(CaMutQC)
library(DT)

function(input, output, session) {
    # Set root directories
    volumes <- c(Working = getwd(), Home = "~", Root = "/", getVolumes()())
    
    # Initialize folder selector (only relevant when report is checked)
    shinyDirChoose(
        input, 
        "dir_selector",
        roots = volumes,
        defaultPath = getwd(),
        defaultRoot = "Working"
    )
    
    # Reactive value for report directory
    report_dir <- reactive({
        if (input$report) {
            if (!is.integer(input$dir_selector)) {
                path <- parseDirPath(volumes, input$dir_selector)
                message("[Report] User selected folder: ", path)  # Print to console
                return(path)
            } else {
                default_path <- getwd()
                message("[Report] Using default working directory: ", default_path)  # Print to console
                return(default_path)
            }
        } else {
            return(NULL)
        }
    })
    
    # Display selected path in UI
    output$dir_path <- renderText({
        if (input$report) {
            if (!is.null(report_dir())) {
                paste("Report will be saved to:", report_dir())
            } else {
                "No folder selected - reports will be saved to working directory"
            }
        }
    })
    
    # Reactive expression for PON file path
    pon_file_path <- reactive({
        if (input$PON_filter) {
            if (!is.null(input$PON_file)) {
                return(input$PON_file$datapath)
            } else {
                default_path <- file.path("dom", "PON_example.txt")
                if (!file.exists(default_path)) {
                    showNotification("Default PON file not found in dom/PON_example.txt", type = "error")
                    return(NULL)
                }
                return(default_path)
            }
        } else {
            return(NULL)
        }
    })
    
    
    bed_file_path <- reactive({
        if (!is.null(input$bed_file)) {
            message("Using uploaded BED file: ", input$bed_file$name)
            return(input$bed_file$datapath)
        } else if (input$bedFilter) {
            # when no bed file is uploaded
            stop("BED filter enabled but no file uploaded")
        } else {
            # when bedFilter set to FALSE
            return(NULL)
        }
    })
    
    # Reactive expression for VCF file path
    vcf_file_path <- reactive({
        if (!is.null(input$vcf_file)) {
            # Return the uploaded file path
            input$vcf_file$datapath[1]
        } else {
            # Return default file path
            default_path <- file.path("dom", "example.vcf")
            # Check if default file exists
            if (!file.exists(default_path)) {
                showNotification("Default VCF file not found in dom/example.vcf", type = "error")
                return(NULL)
            }
            return(default_path)
        }
    })
    
    maf_data <- reactiveVal(NULL)
    tmb_val <- reactiveVal(NULL)
    results <- reactiveVal(NULL)
    
    # Download example VCF
    output$download_example <- downloadHandler(
        filename = function() {
            if (!is.null(input$vcf_file) && nrow(input$vcf_file) > 0) {
                input$vcf_file$name[1]  # Name of the uploaded file
            } else {
                "example.vcf"
            }
        },
        content = function(file) {
            if (!is.null(input$vcf_file) && nrow(input$vcf_file) > 0) {
                # Download the first uploaded VCF file
                file.copy(input$vcf_file$datapath[1], file)
            } else {
                # Download the example VCF from www/
                file.copy("dom/example.vcf", file)
            }
        }
    )
    
    observeEvent(input$run_filter, {
        req(vcf_file_path())
        if (input$PON_filter) {
            req(pon_file_path())
        }
        # If multiple VCFs, use multiVCF=TRUE
        maf <- vcfToMAF(vcf_file_path())
        # process gene list
        if (input$genelist == ""){
            genes <- NULL
        }else{
            genes <- unlist(strsplit(input$genelist, split = ","))
        }
        # handle ponformat when PON filter is set to FALSE
        if (!(input$PON_filter)){
            PON_format_final <- '.vcf'
        }else{
            PON_format_final <- input$PON_file_type
        }
        # Apply filtration
        maf_filtered <- mutFilterCom(
            maf,
            panel = input$panel,
            report = input$report,
            reportDir = selected_dir(),
            normalAD = input$normalAD,
            normalDP = input$normalDP,
            VAF = input$VAF,
            VAFratio = input$VAFratio,
            tumorDP = input$tumorDP,
            tumorAD = input$tumorAD,
            dbsnpCutoff = input$dbsnpCutoff,
            nonCutoff = input$nonCutoff,
            PONfile = pon_file_path(),
            PONformat = PON_format_final,
            progressbar = FALSE,
            SBmethod = input$sb,
            SBscore = input$sbscore,
            maxIndelLen = input$maxIndelLen,
            minInterval = input$minInterval,
            tagFILTER = input$tagFILTER,
            dbVAF = input$dbVAF,
            ExAC = input$ExAC,
            Genomesprojects1000 = input$Genomesprojects1000,
            gnomAD = input$gnomAD,
            dbSNP = input$dbSNP,
            keepCOSMIC = input$keepCOSMIC,
            keepType = input$keepType,
            bedFilter = input$bedFilter,
            bedFile = bed_file_path(),
            bedHeader = input$bedHeader,
            mutFilter = FALSE,
            ESP6500 = input$ESP6500,
            mutType = input$mutType,
            genelist = genes,
            # do not run TMB here
            TMB = FALSE
        )
        maf_data(maf_filtered)
        # Calculate TMB if requested
        if (input$TMB) {
            tmb_val(calTMB(maf, assay = input$assay, 
                           bedFile = bed_file_path(), bedHeader = input$bedHeader,
                           genelist = genes, mutType = input$mutType, 
                           bedFilter = input$bedFilter))
        } else {
            tmb_val(NULL)
        }
        filtered_results <- maf_data()
        rownames(filtered_results) <- 1:nrow(filtered_results)
        filtered_results <- filtered_results[,c(125,1,5,6,7,9,10,11,13)]
        results(filtered_results)
    })
    
    output$labeled_table <- renderDT({
        req(maf_data())
        req(results())
        datatable(results(), options = list(pageLength = 20))
    })
    
    output$download_maf <- downloadHandler(
        filename = function() {
            paste0("labeled_", Sys.Date(), ".maf")
        },
        content = function(file) {
            write.table(maf_data(), file, sep = "\t", row.names = FALSE, quote = FALSE)
        }
    )
    
    output$tmb_value <- renderPrint({
        req(tmb_val())
        paste("TMB value (before filtering):", tmb_val())
    })
}