globalVariables(c("protein_sample_info", "protein_data", "batch_indicator",
    "signature_data", "bladder_meta", "bladder_data", "merged_IDs"))

#' Batch and Condition indicator for protein expression data
#'
#' This data consists of two batches and two conditions
#' corresponding to case and control for the protein expression data
#'
#' @name protein_sample_info
#' @docType data
#' @format A data frame with 24 rows and 2 variables:
#' \describe{
#'     \item{batch}{Batch Indicator}
#'     \item{category}{Condition (Case vs Control) Indicator}
#' }
#' @keywords datasets
#' @usage data(protein_sample_info)
"protein_sample_info"

#' Protein data with 39 protein expression levels
#'
#' This data consists of two batches and two conditions
#' corresponding to case and control. The columns are case/control
#' samples, and the rows represent 39 different proteins.
#'
#' @name protein_data
#' @docType data
#' @format A data frame with 39 rows and 24 variables
#' @keywords datasets
#' @usage data(protein_data)
"protein_data"


#' Batch and Condition indicator for signature data
#'
#' This dataset is from signature data captured when activating different growth
#' pathway genes in human mammary epithelial cells (GEO accession: GSE73628).
#' This data consists of three batches and ten different conditions
#' corresponding to control and nine different pathways.
#'
#' @name batch_indicator
#' @docType data
#' @format A data frame with 89 rows and 2 variables:
#' \describe{
#'     \item{batch}{batch}
#'     \item{condition}{condition}
#' }
#' @keywords datasets
#' @usage data(batch_indicator)
"batch_indicator"

#' Signature data with 1600 gene expression levels
#'
#' This data consists of three batches and ten conditions.
#' The columns are samples, and the rows represent
#' 1600 different genes.
#'
#' @name signature_data
#' @docType data
#' @format A data frame with 1600 rows and 89 variables
#' @keywords datasets
#' @usage data(signature_data)
"signature_data"

#' Bladder data upload
#' This function uploads the Bladder data set from the bladderbatch package.
#' This dataset is from bladder cancer data with 22,283 different microarray
#' gene expression data. It has 57 bladder samples with 3 metadata variables
#' (batch, outcome and cancer). It contains 5 batches, 3 cancer types (cancer,
#' biopsy, control), and 5 outcomes (Biopsy, mTCC, sTCC-CIS, sTCC+CIS, and
#' Normal). Batch 1 contains only cancer, 2 has cancer and controls, 3 has only
#' controls, 4 contains only biopsy, and 5 contains cancer and biopsy
#'
#' @usage bladder_data_upload()
#' @return a SE object with counts data and metadata
#'
#' @examples
#' library(bladderbatch)
#' se_object <- bladder_data_upload()
#'
#' @export

bladder_data_upload <- function() {
    if (!requireNamespace("bladderbatch")) {
        stop("You need to install the 'bladderbatch' package to use this
            data set.")
    }
    data(bladderdata, package = "bladderbatch", envir = environment())
    bladderEset <- bladderEset
    pheno <- pData(bladderEset) %>% select(-sample)
    edata <- exprs(bladderEset)
    se_object <- BatchQC::summarized_experiment(edata, pheno)
    colData(se_object)[['batch']] <- as.factor(colData(se_object)[['batch']])
    return(se_object)
}

#' BMI and matched sample names for TB data
#'
#' This is support data for the TB data set that contains the BMI data and ID
#' numbers from both the curatedTBData database and the original study the data
#' was used in
#'
#' @name merged_IDs
#' @docType data
#' @format A data frame with 91 rows and 3 columns
#' \describe{
#'     \item{subjectID_curatedTBData}{Subject ID found in curatedTBData}
#'     \item{subjectID_TB_Paper}{Subject ID in the original paper}
#'     \item{BMI}{subject's BMI from the original study}
#' }
#' @keywords datasets
#' @usage data(merged_IDs)
"merged_IDs"


#' TB data upload
#' This function uploads the TB data set from the curatedTBData package.
#'
#' @usage tb_data_upload()
#' @return a SE object with raw counts data and metadata
#' @import SummarizedExperiment
#' @import dplyr
#'
#' @examples
#' library(curatedTBData)
#' se_object <- tb_data_upload()
#'
#' @export

tb_data_upload <- function() {
    if (!requireNamespace("curatedTBData")) {
        stop("You need to install the 'curatedTBData' package to use this
            data set.")
    } else if (!requireNamespace("MultiAssayExperiment")) {
        stop("You need to install the 'MultiAssayExperiment' package to use this
            data set.")
    }
    curatedData <- curatedTBData::curatedTBData(c("GSE152218", "GSE101705"),
        dry.run = FALSE, curated.only = FALSE)

    batch1_se <- MultiAssayExperiment::experiments(
        curatedData$GSE152218)$object_raw
    batch1_data <- SummarizedExperiment::assays(batch1_se)$assay_raw
    batch1_metadata <- MultiAssayExperiment::colData(curatedData$GSE152218)
    batch1_metadata$Experiment <- rep("GSE152218", length(batch1_metadata[, 1]))
    batch1_metadata <- batch1_metadata %>% as.data.frame() %>%
        select("TBStatus", "HIVStatus", "BMI", "Experiment")
    batch2_data <- MultiAssayExperiment::experiments(
        curatedData$GSE101705)$assay_reprocess_hg38 %>%
        as.data.frame() %>%
        select(-"GSM2712712")
    batch2_metadata <- MultiAssayExperiment::colData(curatedData$GSE101705)
    batch2_metadata <-
        batch2_metadata[-which(rownames(batch2_metadata) == "GSM2712712"), ]
    batch2_metadata$Experiment <- rep("GSE101705", length(batch2_metadata[, 1]))

    batch2_metadata <- BMI_data(batch2_metadata)
    batch2_metadata <- batch2_metadata %>% as.data.frame() %>%
        select("TBStatus", "HIVStatus", "BMI", "Experiment")

    all_data <- merge(batch1_data, batch2_data, by = 0)
    rownames(all_data) <- all_data$Row.names
    all_data <- all_data %>% select(-1) %>% as.matrix() %>% round()
    all_metadata <- rbind(batch1_metadata, batch2_metadata)

    se <- summarized_experiment(all_data, all_metadata)
    colData(se)[['Experiment']] <- as.factor(colData(se)[['Experiment']])
    colData(se)[['TBStatus']] <- as.factor(colData(se)[['TBStatus']])
    colData(se)[['HIVStatus']] <- as.factor(colData(se)[['HIVStatus']])
    colData(se)[['BMIcat']] <- as.factor(vapply(colData(se)[['BMI']],
        function(x) {
            if (x < 16) {
                return("mal")
            }else {
                return("well")
                }
            }, FUN.VALUE = character(1)))
    return(se)
}

#' This function returns BMI data that comes form the data in "Comparing
#' tuberculosis gene signatures in malnourished individuals using the
#' TBSignatureProfiler" paper. Subject IDs were matched as shown on
#' "github.com/jessmcc22/BatchQCv2_Manuscript/blob/devel/R/subjectID_match.R"
#'
#' @param meta dataframe; metadata that needs to be matched to BMI
#' @usage BMI_data(meta)
#' @return dataframe provided as input with BMI info added
#'
BMI_data <- function(meta) {
    data("merged_IDs", envir = environment())
    subjectID_curatedTBData <- merged_IDs$subjectID_curatedTBData
    BMI <- merged_IDs$BMI
    IDs <- data.frame(subjectID_curatedTBData, BMI)
    meta$BMI <- rep(NA, length(meta$Experiment))
    for (i in seq_len(nrow(meta))){
        position <- which(IDs$subjectID_curatedTBData == rownames(meta)[i])
        meta$BMI[i] <- IDs$BMI[position]
    }
    return(meta)
}
