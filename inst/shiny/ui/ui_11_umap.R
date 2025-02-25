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
            actionButton('umap_plot', label = 'Here we go!')
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotOutput('umap')
        )
    )
)
