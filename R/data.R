globalVariables(c("protein_sample_info", "protein_data", "batch_indicator",
    "signature_data", "bladder_meta", "bladder_data"))

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
    all_data <- all_data %>% select(-1) %>% as.matrix()
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
    subjectID_curatedTBData <- c("GSM4609427", "GSM4609429", "GSM4609424",
        "GSM4609396", "GSM4609417", "GSM4609387", "GSM4609415", "GSM4609408",
        "GSM4609426", "GSM4609412", "GSM4609400", "GSM4609423", "GSM4609416",
        "GSM4609425", "GSM4609420", "GSM4609386", "GSM4609422", "GSM4609406",
        "GSM4609419", "GSM4609414", "GSM4609389", "GSM4609432", "GSM4609428",
        "GSM4609397", "GSM4609421", "GSM4609410", "GSM4609413", "GSM4609411",
        "GSM4609433", "GSM4609431", "GSM4609430", "GSM4609418", "GSM4609407",
        "GSM4609402", "GSM4609404", "GSM4609398", "GSM4609394", "GSM4609392",
        "GSM4609399", "GSM4609391", "GSM4609388", "GSM4609409", "GSM4609403",
        "GSM4609395", "GSM4609393", "GSM4609390", "GSM4609405", "GSM4609401",
        "GSM2712704", "GSM2712705", "GSM2712706", "GSM2712677", "GSM2712707",
        "GSM2712696", "GSM2712698", "GSM2712699", "GSM2712708", "GSM2712700",
        "GSM2712709", "GSM2712710", "GSM2712711", "GSM2712678", "GSM2712713",
        "GSM2712714", "GSM2712679", "GSM2712715", "GSM2712680", "GSM2712681",
        "GSM2712702", "GSM2712716", "GSM2712717", "GSM2712682", "GSM2712718",
        "GSM2712683", "GSM2712684", "GSM2712692", "GSM2712685", "GSM2712686",
        "GSM2712719", "GSM2712693", "GSM2712703", "GSM2712687", "GSM2712694",
        "GSM2712697", "GSM2712688", "GSM2712689", "GSM2712690", "GSM2712695",
        "GSM2712691", "GSM2712676", "GSM2712701")
    BMI <- c(13.00, 13.61, 13.72, 13.87, 13.88, 13.92, 13.95, 13.95, 14.04,
        14.11, 14.21, 14.38, 14.43, 14.77, 14.80, 14.81, 15.01, 15.02, 15.02,
        15.02, 15.09, 15.17, 15.59, 15.61, 15.66, 15.68, 15.72, 15.75, 15.81,
        15.81, 15.82, 15.86, 19.09, 19.52, 21.01, 21.26, 22.31, 22.44, 22.81,
        22.87, 24.15, 26.31, 27.20, 27.44, 27.81, 28.12, 30.83, 33.16, 28.80,
        24.80, 24.30, 21.50, 26.30, 15.10, 14.30, 15.60, 23.50, 15.10, 24.10,
        38.70, 23.80, 19.40, 26.00, 21.80, 21.70, 22.30, 19.20, 24.90, 14.20,
        22.00, 25.30, 20.10, 25.00, 20.40, 19.20, 14.30, 18.90, 20.30, 21.70,
        17.90, 15.60, 18.60, 16.50, 11.50, 23.20, 26.00, 19.30, 21.40, 19.00,
        20.90, 13.80)
    IDs <- data.frame(subjectID_curatedTBData, BMI)
    meta$BMI <- rep(NA, length(meta$Experiment))
    for (i in seq_len(rownames(meta))){
        position <- which(IDs$subjectID_curatedTBData == rownames(meta)[i])
        meta$BMI[i] <- IDs$BMI[position]
    }
    return(meta)
}
