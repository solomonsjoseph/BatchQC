
tabPanel(
    "lambda/Variation Analysis",

    # Application title
    titlePanel("lambda/Variation Analysis"),

    sidebarLayout(sidebarPanel(
        h3("Variation Analysis"),
        selectizeInput('variation_assay', 'Select Assay Name', choices = "",
            options = list(placeholder = 'Please select an option below',
            onInitialize = I('function() { this.setValue(""); }'))),
        selectizeInput('variation_batch', 'Select Batch Variable', choices = "",
            options = list(placeholder = 'Please select an option below',
            onInitialize = I('function() { this.setValue(""); }'))),
        selectizeInput('variation_condition', 'Select Condition Variable(s)',
            choices = "",
            multiple = TRUE,
            options = list(placeholder = 'Please select an option below',
            onInitialize = I('function() { this.setValue(""); }'))),
        withBusyIndicatorUI(actionButton('variation', label = 'Here we go!'))
    ),
    mainPanel(
        tabsetPanel(
            tabPanel('Lambda Statistic',
                h4(strong("Usage")),
                h5("The lambda statistic should be applied to uncorrected data
                    sets to aid in determining if a batch correction should be
                    applied to the uncorrected data. It is not appropriate to
                    use on corrected data sets."),
                textOutput('lambda_rec')
            ),
            tabPanel("Explained Variation - Individual Variable",
                h5("The boxplot and p-value table display the individual, or raw
                    variation, explained by each variable."),
                     plotOutput('EV_show_plot'),
                     dataTableOutput('EV_show_table')
            ),
            tabPanel("Explained Variation - Residual",
                h5("The boxplot and p-value table display the residual
                    varaition."),
                plotOutput('EV_residual_show_plot'),
                dataTableOutput('EV_residual_show_table')
            ),
            tabPanel("Individual Variation Variable/Batch Ratio",
                h5("The boxplot and table display the individual, or raw
                    variation, divided by the batch. A ratio less that 0
                    indicates that batch has a stronger affect than the
                    variable of interest."),
                conditionalPanel(condition = "input.variation_condition == ''",
                    h5(strong("Condition is required for ratio"))
                ),
                conditionalPanel(condition = "input.variation_condition != ''",
                    plotOutput('EV_ratio_plot'),
                    dataTableOutput('EV_ratio_table')
                    )
            ),
            tabPanel("Residual Variation Variable/Batch Ratio",
                h5("The boxplot and table display the residual divided by the
                    batch. A ratio less that 1 indicates that batch has a
                    stronger affect than the variable of interest."),
                conditionalPanel(condition = "input.variation_condition == ''",
                    h5(strong("Condition is required for residual ratio"))
                ),
                conditionalPanel(condition = "input.variation_condition != ''",
                    plotOutput('EV_residual_ratio_plot'),
                    dataTableOutput('EV_residual_ratio_table')
                )
            )
        )
        )
    )
)
