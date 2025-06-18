tabPanel('Batch Correction/Normalization',
    titlePanel("Batch Correction/Normalization"),
    sidebarLayout(
        sidebarPanel(
            selectizeInput('correction_assay',
                'Choose the assay to modify',
                multiple = FALSE,
                choices = c(''),
                selected = NULL,
                options = list(placeholder =
                        'Please select an option below',
                    onInitialize = I(
                        'function() { this.setValue(""); }')
                ))
        ),
        mainPanel(
            tabsetPanel(
                tabPanel('Batch Correction',
                    h4(strong("Usage")),
                    conditionalPanel(condition = "input.correction_method == 'ComBat-Seq'",
                        h5("ComBat-Seq uses a negative binomial regression to model batch effects. It requires untransformed, raw count data to adjust for batch effect. Please use this option with a counts assay")
                    ),
                    conditionalPanel(condition = "input.correction_method == 'ComBat'",
                        h5("ComBat corrects for Batch effect using a parametric empirical Bayes framework and data should be cleaned and normalized. Therefore, please select a normalized assay to run this on.")
                    ),
                    conditionalPanel(condition = "input.correction_method == 'limma'",
                        h5("limma batch correction fits a linear model to the data, including batch and regular treatments, then removes the component due to the batch effects. Please run limma on log expression data.")
                    ),
                    conditionalPanel(condition = "input.correction_method == 'sva'",
                        h5("sva correction identifies suurogate variables to correct for unknown batch effects. This is the two-step implementation available as psva in the sva package.")
                    ),
                    selectizeInput('correction_method', 'Choose correction method',
                        multiple = FALSE,
                        choices = c('ComBat-Seq', 'ComBat', 'limma', 'sva'),
                        selected = NULL,
                        options = list(placeholder =
                                'Please select an option below',
                            onInitialize = I(
                                'function() { this.setValue(""); }'
                            ))),
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
                tabPanel('Normalization',
                    h4(strong("Usage")),
                    conditionalPanel(condition = "input.normalization_method == 'CPM'",
                        h5("CPM calculates the counts mapped to a feature relative to the total counts mapped to a sample times one million.")
                    ),
                    conditionalPanel(condition = "input.normalization_method == 'DESeq'",
                        h5("DESeq calculates the counts mapped to a feature divided by sample-specific size factors. Size factors are determined by the median ratio of gene counts relative to the geometric mean per feature.")
                        ),
                    conditionalPanel(condition = "input.normalization_method == 'edgeR'",
                        h5("edgeR calculates scale factors using a trimmed mean of M-values between each pair of samples and multiplies the scale factors with the original library size to get the normalized library size.")
                        ),
                    selectizeInput('normalization_method',
                        'Choose normalization method',
                        multiple = FALSE,
                        choices = c('CPM', 'DESeq', 'edgeR', 'none'),
                        selected = NULL,
                        options = list(placeholder =
                                'Please select an option below',
                            onInitialize = I(
                                'function() { this.setValue(""); }'
                            ))),
                    textInput(inputId = 'normalized_assay_name', 'Name for the normalized assay',
                        value = ''),
                    checkboxInput('log', 'log(x+1) transform the results'),
                    withBusyIndicatorUI(actionButton(inputId = 'normalize',
                        label = 'Normalize')),
                    br()
                )
            )
        )
    )
)
