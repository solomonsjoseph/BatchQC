### UPLOAD TAB ###

## Obtain the counts matrix and count table location
observeEvent(input$counts, {
    req(input$counts)
    reactivevalue$counts_location <- input$counts$datapath
    reactivevalue$counts <- read.table(reactivevalue$counts_location,
        header = TRUE,
        row.names = 1,
        sep = get.delim(reactivevalue$counts_location,
            n = 10,
            delims = c('\t', ',')),
        check.names = FALSE)
    output$counts_header <- renderDT((datatable(reactivevalue$counts)))
    output$counts_dimensions <- renderText(paste(dim(reactivevalue$counts),
        c('observations and', 'samples')))
})

## Obtain the metadata matrix and metadata table location
observeEvent(input$md, {
    req(input$md)
    reactivevalue$metadata_location <- input$md$datapath

    reactivevalue$metadata <- read.table(reactivevalue$metadata_location,
        header = TRUE, row.names = 1,
        sep = get.delim(reactivevalue$metadata_location,
            n = 10,
            delims = c('\t', ',')))

    output$metadata_header <- renderDT(datatable(reactivevalue$metadata))
})

## Obtain count matrix and metadata from the summarized experiment input
observeEvent(input$se, {
    req(input$se)
    reactivevalue$se_location <- input$se$datapath
    se <- readRDS(input$se$datapath)
    reactivevalue$se <- se
    output$counts_header <- renderDT(datatable(assays(reactivevalue$se)[[1]]))
    output$counts_dimensions <- renderText(paste(dim(reactivevalue$se),
        c('observations and', 'samples')))

    reactivevalue$metadata <- as.data.table(colData(reactivevalue$se))
    output$metadata_header <- renderDT(datatable(reactivevalue$metadata))
})

## Obtain count matrix and count location for example data
observeEvent(input$exampleData, {
    req(input$exampleData)
    if (input$exampleData == "proteinData") {

        data(protein_data)
        reactivevalue$counts <- protein_data
        output$counts_header <- renderDT(datatable(reactivevalue$counts))
        output$counts_dimensions <- renderText(paste(dim(reactivevalue$counts),
            c('observations and', 'samples')))

        data(protein_sample_info)
        reactivevalue$metadata <- protein_sample_info
        output$metadata_header <- renderDT(datatable(reactivevalue$metadata))
    }else if (input$exampleData == "signatureData") {
        data(signature_data)
        reactivevalue$counts <- signature_data
        output$counts_header <- renderDT(datatable(reactivevalue$counts))
        output$counts_dimensions <- renderText(paste(dim(reactivevalue$counts),
            c('observations and', 'samples')))

        data(batch_indicator)
        reactivevalue$metadata <- batch_indicator
        output$metadata_header <- renderDT(datatable(reactivevalue$metadata))
    }else if (input$exampleData == "bladderData") {
        bladder_data <- bladder_data_upload()
        reactivevalue$counts <- assays(bladder_data)$counts
        output$counts_header <- renderDT(datatable(reactivevalue$counts))
        output$counts_dimensions <- renderText(paste(dim(reactivevalue$counts),
            c('observations and', 'samples')))

        reactivevalue$metadata <- as.data.frame(colData(bladder_data))
        output$metadata_header <- renderDT(datatable(reactivevalue$metadata))

    }else if (input$exampleData == "TBData") {
        TB_data <- tb_data_upload()
        reactivevalue$counts <- assays(TB_data)$counts
        output$counts_header <- renderDT(datatable(reactivevalue$counts))
        output$counts_dimensions <- renderText(paste(dim(reactivevalue$counts),
            c('observations and', 'samples')))

        reactivevalue$metadata <- as.data.frame(colData(TB_data))
        output$metadata_header <- renderDT(datatable(reactivevalue$metadata))
     }
})

## Create summarized experiment object and set up plot options
observeEvent(input$submit, {
    withBusyIndicatorServer("submit", {
        if (input$uploadChoice == "countFile" &
                !is.null(reactivevalue$counts_location) &
                !is.null(reactivevalue$metadata_location)) {
            se <- summarized_experiment(reactivevalue$counts,
                reactivevalue$metadata)
            se <- se[which(rownames(se) != 'NA')]
        } else if (input$uploadChoice == "seObject" &
                !is.null(input$se$datapath)) {
            se <- readRDS(reactivevalue$se_location)
            se <- se[rowSums(se@assays@data[[1]]) > 0, ]
        } else if (input$uploadChoice == "example") {
            se <- summarized_experiment(reactivevalue$counts,
                reactivevalue$metadata)
            if (input$exampleData == "proteinData") {
                assayNames(se) <- "mass_to_charge_ratio"
            }else if (input$exampleData == "signatureData") {
                assayNames(se) <- "log_intensity"
            }else if (input$exampleData == "bladderData") {
                assayNames(se) <- "log_CPM"
            }
        }

        reactivevalue$se <- se
        reactivevalue$metadata <- data.frame(colData(reactivevalue$se))
        output$se_download <- renderPrint(reactivevalue$se)
        output$downloadData <- downloadHandler(
            filename = function() {
                paste("se", ".RDS", sep = "")
            },
            content = function(file) {
                saveRDS(reactivevalue$se, file)
            }
        )

        # Display metadata table
        output$metadata <- renderDataTable(data.table(data.frame(
            colData(reactivevalue$se)),
            keep.rownames = TRUE))

        # Add options to input selections
        setupSelections()
    })
})

## Complete NB Check
observeEvent(input$nb_check, {
    req(input$nb_test, input$counts_matrix,
        input$condition_of_interest)
    withBusyIndicatorServer("nb_check", {
        check_res <- goodness_of_fit_DESeq2(reactivevalue$se,
            input$counts_matrix,
            input$condition_of_interest,
            input$nb_variables,
            input$num_genes)
        output$recommendation <- renderText(check_res$recommendation)
        output$nb_histogram <- renderPlot(check_res$res_histogram)
        output$reference <- renderText(check_res$reference)
    })
})

## Compute Lambda Statistic
observeEvent(input$lambda_stat, {
    req(input$lambda_matrix,
        input$batch_ind,
        input$group_ind)
    withBusyIndicatorServer("lambda_stat", {
        lambda_res <- compute_lambda(dat = assays(reactivevalue$se)[[input$lambda_matrix]],
            batchind = colData(reactivevalue$se)[, input$batch_ind],
            groupind = colData(reactivevalue$se)[, input$group_ind])
        output$lambda_table <- renderDataTable(data.table(lambda_res))
    })
})
