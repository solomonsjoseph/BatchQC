tabPanel("umap",

    # Application title
    titlePanel("umap Analysis"),

    # Sidebar
    sidebarLayout(
        sidebarPanel(
            # List of assays to plot from se object
            selectizeInput('umap_assay',
                'Assay to plot',
                choices = c(),
                multiple = FALSE),

            # Options for plots
            selectizeInput('batch_variable',
                'Select batch varaiable (to display color)',
                choices = c(),
                multiple = FALSE,
                selected = NULL,
                options = list(placeholder = 'Please select an option below',
                    onInitialize = I('function() { this.setValue(""); }'))),
            selectizeInput('covar_variable',
                'Select biological varaiable (to display shape)',
                choices = c(),
                multiple = FALSE,
                selected = NULL,
                options = list(placeholder = 'Please select an option below',
                    onInitialize = I('function() { this.setValue(""); }'))),
            numericInput('num_neighbors',
                'Choose how many nearest neighbors to use',
                value = 15,
                min = 0,
                max = 500),
            numericInput('distance',
                'min distance between points',
                value = 0.1,
                min = 0,
                max = 1),
            numericInput('spread',
                'spread',
                value = 1,
                min = 0,
                max = 100),
            checkboxInput('exploratory',
                'create an exploratory matrix',
                value = FALSE),
            actionButton('umap_plot', label = 'Here we go!')
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotOutput('umap')
        )
    )
)
