## Update batch effect corrected assay name
observe( {
    req(input$correction_assay, input$correction_batch, input$correction_method)
    if (!is.null(input$correction_covariates)) {
        if (!is.null(input$group_for_batch)) {
            updateTextInput(session = session,
                inputId = 'corrected_assay_name',
                'Name for the corrected assay',
                value = paste(input$correction_assay,
                    input$correction_batch,
                    input$group_for_batch,
                    input$correction_method,
                    paste(input$correction_covariates,
                        collapse = '_'),
                    sep = '_'))
        }else {
            updateTextInput(session = session, inputId = 'corrected_assay_name',
                'Name for the corrected assay',
                value = paste(input$correction_assay,
                    input$correction_batch,
                    input$correction_method,
                    paste(input$correction_covariates,
                        collapse = '_'),
                    sep = '_'))
        }
    }else {
        if (!is.null(input$group_for_batch)) {
            updateTextInput(session = session, inputId = 'corrected_assay_name',
                'Name for the corrected assay',
                value = paste(input$correction_assay,
                    input$correction_batch,
                    input$group_for_batch,
                    input$correction_method,
                    sep = '_'))
        }else if (input$correction_method == 'sva') {
            updateTextInput(session = session, inputId = 'corrected_assay_name',
                            'Name for the corrected assay',
                            value = paste(input$correction_assay,
                                          input$correction_batch,
                                          input$group_for_batch,
                                          input$correction_method,
                                          input$psva,
                                          sep = '_'))
        }else if (input$correction_method == 'svaseq') {
            updateTextInput(session = session, inputId = 'corrected_assay_name',
                            'Name for the corrected assay',
                            value = paste(input$correction_assay,
                                          input$correction_batch,
                                          input$group_for_batch,
                                          input$correction_method,
                                          input$num_sv,
                                          sep = '_'))
        }else {
            updateTextInput(session = session, inputId = 'corrected_assay_name',
                'Name for the corrected assay',
                value = paste(input$correction_assay,
                    input$correction_batch,
                    input$correction_method,
                    sep = '_'))
        }
    }
})

## Complete NB Check
observeEvent(input$nb_check, {
    req(input$nb_test, input$correction_assay,
        input$condition_of_interest)
    withBusyIndicatorServer("nb_check", {
        check_res <- goodness_of_fit_DESeq2(reactivevalue$se,
            input$correction_assay,
            input$condition_of_interest,
            input$nb_variables,
            input$num_genes)
        output$recommendation <- renderText(check_res$recommendation)
        output$nb_histogram <- renderPlot(check_res$res_histogram)
        output$reference <- renderText(check_res$reference)
    })
})

## Run AIC

observeEvent(input$compute_aic, {
    req(reactivevalue$se, input$correction_assay,
        input$aic_batch, input$aic_covar, input$aic_nb_maxit,
        input$aic_zero_filt_percent)
    withBusyIndicatorServer("compute_aic", {
        aic_res <- compute_aic(reactivevalue$se,
            input$correction_assay,
            input$aic_batch,
            input$aic_covar,
            input$aic_nb_maxit,
            input$aic_zero_filt_percent)
        output$total_aic <- renderTable({
            total_aic_data <- aic_res[["total_AIC"]]
            df <- data.frame(
                Model = c("NB_AIC", "Lognormal_AIC", "Voom_AIC"),
                AIC_Value = round(as.numeric(total_aic_data), 3),
                stringsAsFactors = FALSE
            )
            return(df)
        }, rownames = FALSE)

        output$min_aic <- renderTable({
            min_aic_data <- aic_res[["min_AIC"]]
            return(min_aic_data)
        }, rownames = FALSE)
    })
})

## Run batch effect correction
observeEvent(input$correct, {
    req(input$correction_assay, input$correction_batch, input$correction_method)
    tryCatch({{
        msg <- sprintf('Start the batch correction process')
        withProgress(message = msg, {
            setProgress(0.5, 'Correcting...')
            reactivevalue$se <- batch_correct(reactivevalue$se,
                input$correction_method,
                input$correction_assay,
                input$correction_batch,
                group = NULL,
                input$correction_covariates,
                input$corrected_assay_name)
            setProgress(1, 'Complete!')
        })
    }},
        error = function(error) {
            showNotification('Confounding', type = "error")
            print(error)
        })
    setupSelections()
    showNotification('Batch Correction Completed', type = "message")
})

## Update normalized assay name
observe({
    req(input$normalization_method, input$correction_assay)
    updateTextInput(session = session, inputId = 'normalized_assay_name',
        'Name for the normalized Assay',
        value = paste(input$correction_assay,
            input$normalization_method, sep = '_'))
})


## Normalize a selected assay
observeEvent(input$normalize, {
    req(input$normalization_method, input$correction_assay,
        input$normalized_assay_name)
    withBusyIndicatorServer("normalize", {
        reactivevalue$se <- normalize_SE(reactivevalue$se,
            input$normalization_method,
            input$log,
            input$correction_assay,
            input$normalized_assay_name)
        setupSelections()
    })
})

