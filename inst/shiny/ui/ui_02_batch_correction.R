tabPanel('Batch Correction/Normalization',
    titlePanel("Batch Correction/Normalization"),
    sidebarLayout(
        sidebarPanel(
            selectizeInput('correction_assay',
                'Choose the assay of interest',
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
                tabPanel('Negative Binomial Check',
                    h4(strong("Usage")),
                    h5("This features allows you to check your data (must be counts/whole number data) to see if it conforms to the required negative binomial assumption needed for various downstream analysis. If The negative binomial assumption is not met, you should normalize your data or perform other preprocessing step and/or use other more appropraite analysis tools."),
                    selectizeInput('nb_test',
                        'Choose the test to perform',
                        multiple = FALSE,
                        choices = c('nb_DESeq2'),
                        selected = NULL,
                        options = list(placeholder =
                                'Please select an option below',
                            onInitialize = I(
                                'function() { this.setValue(""); }'
                            ))),
                    selectizeInput('condition_of_interest',
                        'Select the variable you are interested in analyzing',
                        multiple = FALSE,
                        choices = c(''),
                        selected = NULL, options = list(placeholder =
                                'Please select an option below',
                            onInitialize = I(
                                'function() { this.setValue(""); }'
                            ))),
                    selectizeInput('nb_variables',
                        'Select other variables you would like to include in your analysis',
                        multiple = TRUE,
                        choices = c(''),
                        selected = NULL,
                        options = list(placeholder =
                                'Please select an option below',
                            onInitialize = I(
                                'function() { this.setValue(""); }'
                            ))),
                    checkboxInput('nb_advanced_options', 'Advanced Options',
                        value = FALSE),
                    conditionalPanel(condition = "input.nb_advanced_options == 1",
                        numericInput('num_genes',
                            'Number of genes to analyze (downsampling)',
                            value = 500,
                            min = 2)),
                    withBusyIndicatorUI(actionButton(inputId = 'nb_check',
                        label = 'Check Distribution')),
                    textOutput('recommendation'),
                    plotOutput('nb_histogram'),
                    textOutput('reference'),
                    br()
                ),
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
                    conditionalPanel(condition = "input.correction_method == 'svaseq'",
                                     h5("svaseq correction is a variant of sva correction for sequencing data")
                    ),
                    selectizeInput('correction_method', 'Choose correction method',
                        multiple = FALSE,
                        choices = c('ComBat-Seq', 'ComBat', 'limma', 'sva', 'svaseq'),
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
                    conditionalPanel(condition = "input.correction_method == 'svaseq'",
                                     checkboxInput('num_sv',
                                                   'Uncheck this if the number of samples is small (the number of latent factors that need to be estimated (n.sv) is set to 1);
                                                   otherwise, svaseq function will estimate n.sv for you.',
                                                   value = FALSE)
                    ),
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
