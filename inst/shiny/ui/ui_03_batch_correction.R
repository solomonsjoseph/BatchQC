tabPanel('Batch Correction',
    titlePanel("Batch Correction"),
    sidebarLayout(
        sidebarPanel(
            selectizeInput('correction_method', 'Choose correction method',
                multiple = FALSE,
                choices = c('ComBat-Seq', 'ComBat', 'limma', 'sva'),
                selected = NULL,
                options = list(placeholder =
                        'Please select an option below',
                    onInitialize = I(
                        'function() { this.setValue(""); }'
                    ))),
            selectizeInput('correction_assay',
                'Choose the assay on which to do correction',
                multiple = FALSE,
                choices = c(''),
                selected = NULL,
                options = list(placeholder =
                        'Please select an option below',
                    onInitialize = I(
                        'function() { this.setValue(""); }')
                )),
            selectizeInput('correction_batch',
                'Select the variable that represents batch (or the experimental variable for sva)',
                multiple = FALSE,
                choices = c(''),
                selected = NULL,
                options = list(placeholder =
                        'Please select an option below',
                    onInitialize = I(
                        'function() { this.setValue(""); }')
                )),
            selectizeInput('correction_covariates',
                'Choose the covariates you would like to preserve (or for sva, include as adjustment variables)',
                multiple = TRUE,
                choices = c(''),
                selected = NULL,
                options = list(placeholder =
                        'Please select an option below',
                    onInitialize = I(
                        'function() { this.setValue(""); }')
                )),
            textInput(inputId = 'corrected_assay_name',
                'Name for the corrected assay'),
            actionButton(inputId = 'correct', label = 'Correct')
        ),
        mainPanel(
            h4(strong("Usage")),
            conditionalPanel(condition = "input.correction_method == 'ComBat-Seq'",
                h5("ComBat-Seq uses a negative binomial regression to model batch effects. It requires untransformed, raw count data to adjust for batch effect. Please use this option with a counts assay")
                ),
            conditionalPanel(condition = "input.correction_method == 'ComBat'",
                h5("ComBat corrects for Batch effect using a parametric empirical Bayes framework and data should be cleaned and normalized. Therefore, please select a normalized assay to run this on.")
                ),
            conditionalPanel(condition = "input.correction_method == 'limma'",
                h5("limma batch correction fits a linear model to the data, including batch and regular treatments, then removes the component due to the batch effects. Please run limma on log expression data.")
            )
        )
    )
)
