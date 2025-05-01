tabPanel("kBET",
         
         # Application title
         titlePanel("k-Nearest Neighbor Batch Effect Test"),
         
         # Sidebar with a slider input for number of bins
         sidebarLayout(
           sidebarPanel(
             # List of assays to plot from se
             selectizeInput('kbet_assay',
                            'Assay',
                            choices = c(),
                            multiple = FALSE,
                            options = list(placeholder ='Please select an option below',
                                           onInitialize = I('function() { this.setValue(""); }'))),
             selectizeInput('kbet_batch',
                            'Select batch',
                            choices = c(),
                            multiple = FALSE,
                            selected = NULL,
                            options = list(placeholder = 'Please select an option below',
                                           onInitialize = I('function() { this.setValue(""); }'))),
             checkboxInput('kbet_advanced_options', 'Show Advanced Options',
                           value = FALSE),
             conditionalPanel(condition = "input.kbet_advanced_options == 1",
                              numericInput('knet_k0',
                                           'Number of nearest neighbors to test on (for NULL, set as 0)',
                                           value = 0,
                                           min = 0)),
             conditionalPanel(condition = "input.kbet_advanced_options == 1",
                              selectizeInput('knet_doPCA',
                                             'Perform a pca prior to knn search?',
                                             choices = c("Yes", "No"),
                                             multiple = FALSE,
                                             selected = "Yes")),
             conditionalPanel(condition = "input.kbet_advanced_options == 1",
                              numericInput('knet_dimPCA',
                                           'Number of dimensions to consider if performing a pca prior to knn search',
                                           value = 50,
                                           min = 1)),
             conditionalPanel(condition = "input.kbet_advanced_options == 1",
                              selectizeInput('knet_heuristic',
                                             'Compute an optimal neighbourhood size k?',
                                             choices = c("Yes", "No"),
                                             multiple = FALSE,
                                             selected = "Yes")),
             conditionalPanel(condition = "input.kbet_advanced_options == 1",
                              numericInput('knet_alpha',
                                           'Significance level',
                                           value = 0.05,
                                           min = 0,
                                           max = 1)),
             actionButton('run_kBET', label = 'Run kBET')
           ),
           
           # Show a plot of the generated distribution
           mainPanel(
             h4(strong("Usage")),
             h5("The k-Nearest Neighbor Batch Effect Test (kBET) is a statistical method used to quantify and test for batch effects in high-dimensional data. evaluate whether samples from different batches are well-mixed in local neighborhoods of the data. A higher rejection rate indicates a more prevalent batch effect in the data."),
             plotOutput('kBET_plot')
           )
         )
)
