accepted <- c("text/csv",
    "text/comma-separated-values",
    "text/plain",
    ".csv")


tabPanel("Upload Data",
    useShinyjs(),
    tags$style(appCSS),
    tags$div(
        class = "jumbotron",
        tags$div(
            class = "container",
            fluidRow(column(7, h1("BatchQC"))),
            tags$p("Batch Effects Quality Control Software"),
            uiOutput("tab")

        )
    ),
    # Application title
    titlePanel("Upload Data"),

    # Place for uploading data
    sidebarLayout(
        sidebarPanel(
            h4("Select the type of input you would like to provide:"),
            radioButtons("uploadChoice", "",
                c("Feature File and Metadata File" = "countFile",
                    "Summarized Experiment Object" = "seObject",
                    "Example Data" = "example"
                )),
            #Only show panel if uploading count and metadata files
            conditionalPanel(condition = "input.uploadChoice == 'countFile'",
                h4("Upload feature and metadata table"),
                tags$div(tags$p(
                    'Metadata file must be a table with headers and sample names.'
                )),
                fileInput(
                    "counts",
                    "Feature table",
                    multiple = FALSE,
                    accept = accepted
                ),
                fileInput("md", "Metadata",
                    multiple = FALSE,
                    accept = accepted)),
            conditionalPanel(condition = "input.uploadChoice == 'seObject'",
                fileInput(
                    "se",
                    "Summarized Experiment",
                    multiple = FALSE,
                    accept = c(accepted, ".RDS")
                )),
            conditionalPanel(condition = "input.uploadChoice == 'example'",
                selectInput("exampleData",
                    "Example Data",
                    choices = c("", "proteinData", "signatureData",
                        "bladderData", "TBData", "No Selection"),
                    selected = ""), ),
            withBusyIndicatorUI(actionButton(inputId = 'submit',
                label = 'Upload'))
        ),

        # Show a table of the inputted data
        mainPanel(
            tabsetPanel(
                tabPanel('Input Summary',
                    h4(strong("Usage")),
                    h5("The table displayed here is a preview of the data selected to be uploaded to BatchQC. Please take a moment to verify that this is the intended data to upload, then select the Upload button to complete uploading the data to BatchQC."),
                    conditionalPanel(condition = "input.uploadChoice == 'example'",
                        h4(strong("Example Data Description")),
                        conditionalPanel(condition = "input.exampleData == 'proteinData'",
                            h5("This data set is mass-spectromery data for protein expression in an unpublished data set.")
                        ),
                        conditionalPanel(condition = "input.exampleData == 'signatureData'",
                            h5("This data set is RNA-seq data from signature data captured when activating different growth pathway genes in human mammary epithelial cells (GEO accession: GSE73628)."),
                            h6("Reference: Rahman M, et al. Activity of distinct growth factor receptor network components in breast tumors uncovers two biologically relevant subtypes. Genome Med. 2017 Apr 26;9(1):40. doi: 10.1186/s13073-017-0429-x.")
                        ),
                        conditionalPanel(condition = "input.exampleData == 'bladderData'",
                            h5("This data set is from bladder cancer data and contains Affymetrix gene expression microarray data. This data set is from the
                                    bladderbatch package which must be installed to use this data example set"),
                            h6("Reference: Leek J., et al. Tackling the widespread and critical impact of batch effects in high-throughput data. Nat Rev Genet 11, 733–739 (2010). https://doi.org/10.1038/nrg2825")
                        ),
                        conditionalPanel(condition = "input.exampleData == 'TBData'",
                            h5("This data set is from a TB cohort study from Pondicherry India and contains RNA-seq counts data."),
                            h6("Reference: Johnson WE, et al. Comparing tuberculosis gene signatures in malnourished individuals using the TBSignatureProfiler. BMC Infect Dis. 2021 Jan 22;21(1):106. doi: 10.1186/s12879-020-05598-z. PMID: 33482742; PMCID: PMC7821401.")
                        )),
                    h4(strong("Data Preview")),
                    DTOutput('counts_header'),
                    textOutput('counts_dimensions'),
                    DTOutput('se_counts'),
                    textOutput('se_dimensions'),
                    br()
                ),
                tabPanel('Full Metadata',
                    h4(strong("Usage")),
                    h5("The table displayed here shows all metadata associated with the uploaded data. Please refer to this table to view sample condition information."),
                    h4(strong("Metadata Table")),
                    DTOutput('metadata_header'),
                    DTOutput('se_meta'),
                    br()
                )
            )
        )
    )
)
