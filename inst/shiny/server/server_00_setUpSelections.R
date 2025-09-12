### UPDATING SELECTIONS ###

setupSelections <- function() {
    # Experimental design
    updateSelectizeInput(session = session, inputId = "design_batch",
        choices = names(colData(reactivevalue$se)),
        selected = NULL)
    updateSelectizeInput(session = session, inputId = "design_covariate",
        choices = names(colData(reactivevalue$se)),
        selected = NULL)

    # Negative Binomial Check
    updateSelectizeInput(session = session, inputId = "counts_matrix",
        choices = assayNames((reactivevalue$se)),
        selected = NULL,
        options = list(placeholder = "Please select an option below",
            onInitialize = I('function() { this.setValue(""); }')))

    updateSelectizeInput(session = session, inputId = "condition_of_interest",
        choices = names(colData(reactivevalue$se)),
        selected = NULL,
        options = list(placeholder = "Please select an option below",
            onInitialize = I('function() { this.setValue(""); }')))

    updateSelectizeInput(session = session, inputId = "nb_variables",
        choices = names(colData(reactivevalue$se)),
        selected = NULL,
        options = list(placeholder = "Please select an option below",
            onInitialize = I('function() { this.setValue(""); }')))

    # AIC Computation
    updateSelectizeInput(session = session, inputId = "aic_batch",
                        choices = names(colData(reactivevalue$se)),
                        selected = NULL,
                        options = list(placeholder = "Please select an option below",
                                        onInitialize = I('function() { this.setValue(""); }')))

    updateSelectizeInput(session = session, inputId = "aic_covar",
                        choices = names(colData(reactivevalue$se)),
                        selected = NULL,
                        options = list(placeholder = "Please select an option below",
                                        onInitialize = I('function() { this.setValue(""); }')))

    # Lambda Statistic
    updateSelectizeInput(session = session, inputId = "lambda_matrix",
        choices = assayNames((reactivevalue$se)),
        selected = NULL,
        options = list(placeholder = "Please select an option below",
            onInitialize = I('function() { this.setValue(""); }')))
    updateSelectizeInput(session = session, inputId = "batch_ind",
        choices = names(colData(reactivevalue$se)),
        selected = NULL)
    updateSelectizeInput(session = session, inputId = "group_ind",
        choices = names(colData(reactivevalue$se)),
        selected = NULL)

    # Batch Correction
    updateSelectizeInput(session = session,
        inputId = "correction_batch",
        choices = (names(colData(reactivevalue$se))),
        selected = NULL,
        options = list(placeholder = "Please select an option below",
            onInitialize = I('function() { this.setValue(""); }')))
    updateSelectizeInput(session = session,
        inputId = "correction_assay",
        choices = (assayNames((reactivevalue$se))),
        selected = NULL,
        options = list(placeholder =
                'Please select an option below',
            onInitialize =
                I('function() { this.setValue(""); }')))
    updateSelectizeInput(session = session,
        inputId = "correction_covariates",
        choices = (names(colData(reactivevalue$se))),
        selected = NULL,
        options = list(placeholder =
                'Please select an option below',
            onInitialize =
                I('function() { this.setValue(""); }')))

    # Heatmap
    updateSelectizeInput(session = session, inputId = 'heatmap_assay_name',
        choices = assayNames((reactivevalue$se)),
        selected = NULL)
    updateSelectInput(session = session, inputId = 'variates_to_display',
        choices = colnames(colData(reactivevalue$se)),
        selected = NULL)
    updateNumericInput(session = session, inputId = 'top_n_heatmap',
        value = dim(reactivevalue$se)[1], min = 2, max = dim(reactivevalue$se)[1])

    # Dendrogram
    updateSelectizeInput(session = session, inputId = 'dend_assay_name',
        choices = assayNames((reactivevalue$se)),
        selected = NULL)
    updateSelectInput(session = session, inputId = 'dend_batch_to_display',
        choices = colnames(colData(reactivevalue$se)),
        selected = NULL)
    updateSelectInput(session = session, inputId = 'dend_category_to_display',
        choices = colnames(colData(reactivevalue$se)),
        selected = NULL)

    # PCA
    updateSelectizeInput(session = session, inputId = 'pca_assays',
        choices = assayNames((reactivevalue$se)),
        selected = NULL)
    updateNumericInput(session = session, inputId = 'top_n_PCA',
        value = dim(reactivevalue$se)[1], min = 2, max = dim(reactivevalue$se)[1])
    updateSelectizeInput(session = session, inputId = 'variates_shape',
        choices = colnames(colData(reactivevalue$se)),
        selected = NULL)
    updateSelectizeInput(session = session, inputId = 'variates_color',
        choices = colnames(colData(reactivevalue$se)),
        selected = NULL)
    updateSelectizeInput(session = session, inputId = 'variates_batch',
        choices = colnames(colData(reactivevalue$se)),
        selected = NULL)
    updateSelectizeInput(session, inputId = 'batch',
        choices = colnames(colData(reactivevalue$se)),
        selected = NULL)
    updateSelectizeInput(session, inputId = 'group',
        choices = colnames(colData(reactivevalue$se)),
        selected = NULL)

    # umap Analysis
    updateSelectizeInput(session = session, inputId = 'umap_assay',
        choices = assayNames((reactivevalue$se)),
        selected = NULL)
    updateSelectizeInput(session = session, inputId = 'batch_variable',
        choices = colnames(colData(reactivevalue$se)),
        selected = NULL)
    updateNumericInput(session = session, inputId = 'num_neighbors',
        value = 15, min = 2,
        max = dim(reactivevalue$se)[1] - 1)
    updateNumericInput(session = session, inputId = 'distance',
        value = 0.1, min = 0.000000000001)
    updateNumericInput(session = session, inputId = 'spread',
        value = 1, min = 0.000000000001)

    # Variation Analysis
    updateSelectizeInput(session = session, inputId = "variation_assay",
        choices = names(assays(reactivevalue$se)),
        selected = NULL)
    updateSelectizeInput(session = session, inputId = "variation_batch",
        choices = names(colData(reactivevalue$se)),
        selected = NULL)

    # Differential expression analysis
    updateSelectizeInput(session = session, inputId = "DE_assay",
        choices = names(assays(reactivevalue$se)),
        selected = NULL)
    updateSelectizeInput(session = session, inputId = "DE_batch",
        choices = names(colData(reactivevalue$se)),
        selected = NULL)

    # kBET Analysis
    updateSelectizeInput(session = session, inputId = "kbet_assay",
                         choices = names(assays(reactivevalue$se)),
                         selected = NULL)
    updateSelectizeInput(session = session, inputId = "kbet_batch",
                         choices = names(colData(reactivevalue$se)),
                         selected = NULL)
}
