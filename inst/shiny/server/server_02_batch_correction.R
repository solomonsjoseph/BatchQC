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

