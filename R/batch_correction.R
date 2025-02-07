#' Batch Correct
#' This function allows you to Add batch corrected count matrix to the SE object
#' @param se SummarizedExperiment object
#' @param method Normalization Method ("ComBat-Seq", "ComBat", "limma")
#' @param assay_to_normalize Which assay use to do normalization
#' @param batch The batch
#' @param group The group variable
#' @param covar list of covariates
#' @param output_assay_name name of results assay
#' @usage batch_correct(se, method, assay_to_normalize, batch, group = NULL,
#' covar, output_assay_name)
#' @return a summarized experiment object with normalized assay appended
#' @import SummarizedExperiment
#' @import sva
#' @examples
#' library(scran)
#' se <- mockSCE()
#' se <- BatchQC::batch_correct(se, method = "ComBat-Seq",
#'                                     assay_to_normalize = "counts",
#'                                     batch = "Mutation_Status",
#'                                     covar = "Treatment",
#'                                     output_assay_name =
#'                                         "ComBat_Seq_Corrected")
#' se <- BatchQC::batch_correct(se, method = "ComBat",
#'                                     assay_to_normalize = "counts",
#'                                     batch = "Mutation_Status",
#'                                     covar = "Treatment",
#'                                     output_assay_name =
#'                                         "ComBat_Corrected")
#' se
#'
#' @export
batch_correct <- function(se, method, assay_to_normalize, batch, group = NULL,
    covar, output_assay_name) {
    se <- se
    batch <- data.frame(colData(se))[, batch]
    if (method == 'ComBat-Seq') {
        se <- ComBat_seq_correction(se, assay_to_normalize, batch, group, covar,
            output_assay_name)
    } else if (method == 'ComBat') {
        se <- ComBat_correction(se, assay_to_normalize, batch, covar,
            output_assay_name)
    } else if (method == 'limma') {
        se <- limma_correction(se, assay_to_normalize, batch, covar,
            output_assay_name)
    }
    return(se)
}

#' ComBat-Seq Correction
#' This function applies ComBat-seq correction to your summarized experiment
#' object
#' @param se SummarizedExperiment object
#' @param assay_to_normalize Assay that should be corrected
#' @param batch The variable that represents batch
#' @param group The group variable
#' @param covar list of covariates
#' @param output_assay_name name of results assay
#' @usage ComBat_seq_correction(se, assay_to_normalize, batch, group, covar,
#' output_assay_name)
#' @return SE object with an added ComBat-seq corrected array
#' @import SummarizedExperiment
#' @import sva

ComBat_seq_correction <- function(se, assay_to_normalize, batch,
    group, covar, output_assay_name) {
    if (is.null(covar)) {
        assays(se)[[output_assay_name]] <- ComBat_seq(as.matrix(
            assays(se)[[assay_to_normalize]]), batch = batch)
    } else {
        if (length(covar) == 1) {
            cov <- data.frame(colData(se))[, covar]
            cov <- as.factor(cov)
            cov <- as.numeric(cov)
            cov <- as.matrix(cov)
            rownames(cov) <- rownames(data.frame(colData(se)))

            if (!is.null(group)) {
                assays(se)[[output_assay_name]] <- ComBat_seq(
                    as.matrix(assays(se)[[assay_to_normalize]]),
                    batch = batch, covar_mod = cov, group = group,
                    full_mod = TRUE)
            } else {
                assays(se)[[output_assay_name]] <- ComBat_seq(as.matrix(
                    assays(se)[[assay_to_normalize]]),
                    batch = batch, covar_mod = cov, group = group)
            }
        } else {
            cov <- data.frame(colData(se))[, covar]
            for (i in seq_len(ncol(cov))) {
                cov[, i] <- as.factor(cov[, i])
                cov[, i] <- as.numeric(cov[, i])
            }

            if (!is.null(group)) {
                assays(se)[[output_assay_name]] <- ComBat_seq(as.matrix(
                    assays(se)[[assay_to_normalize]]),
                    batch = batch, covar_mod = cov, group = group,
                    full_mod = TRUE)
            } else {
                assays(se)[[output_assay_name]] <- ComBat_seq(as.matrix(
                    assays(se)[[assay_to_normalize]]),
                    batch = batch, covar_mod = cov, group = group)
            }
        }
    }
    return(se)
}

#' ComBat Correction
#' This function applies ComBat correction to your summarized experiment object
#' @param se SummarizedExperiment object
#' @param assay_to_normalize Assay that should be corrected
#' @param batch The variable that represents batch
#' @param covar list of covariates
#' @param output_assay_name name of results assay
#' @return SE object with an added ComBat corrected array
#' @import SummarizedExperiment
#' @import sva

ComBat_correction <- function(se, assay_to_normalize, batch,
    covar, output_assay_name) {
    if (is.null(covar)) {
        assays(se)[[output_assay_name]] <-
            ComBat(dat = assays(se)[[assay_to_normalize]], batch = batch)
    } else {
        if (length(covar) == 1) {
            cov <- data.frame(colData(se))[, covar]
            cov <- as.factor(cov)
            cov <- as.numeric(cov)
            cov <- data.frame(cov)
            colnames(cov) <- covar
            rownames(cov) <- rownames(data.frame(colData(se)))

            model <- stats::model.matrix(stats::as.formula(
                paste0('~', colnames(cov))), data = cov)
            results <- ComBat(dat = assays(se)[[assay_to_normalize]],
                batch = batch,
                mod = model)
            results[is.na(results)] <- 0
            assays(se)[[output_assay_name]] <- results
        } else {
            cov <- data.frame(colData(se))[, covar]

            for (i in seq_len(ncol(cov))) {
                cov[, i] <- as.factor(cov[, i])
                cov[, i] <- as.numeric(cov[, i])
            }

            cov <- data.frame(cov)
            rownames(cov) <- rownames(data.frame(colData(se)))
            colnames(cov) <- covar

            linearmodel <- stats::as.formula(paste0('~',
                paste(colnames(cov),
                    sep = '+')))
            model <- stats::model.matrix(linearmodel, data = cov)

            results <- ComBat(dat = assays(se)[[assay_to_normalize]],
                batch = batch,
                mod = model)
            results[is.na(results)] <- 0
            assays(se)[[output_assay_name]] <- results

        }
    }
    return(se)
}

#' Limma Correction
#' This function applies limma batch correction to your provided assay
#' @param se SummarizedExperiment object
#' @param assay_to_normalize Log assay that should be corrected
#' @param batch Factor containing batch information
#' @param covar list of covariates
#' @param output_assay_name name of results assay
#' @return SE object with an added limma corrected array
#' @import SummarizedExperiment
#' @importFrom limma removeBatchEffect
#'
limma_correction <- function(se, assay_to_normalize, batch, covar,
    output_assay_name) {
    if (is.null(covar)) {
        limma_corrected <- limma::removeBatchEffect(
            assays(se)[[assay_to_normalize]],
            batch = batch
            #design = model.matrix(~)
        )
    }else {
        cov <- data.frame(colData(se))[, covar]

        for (i in seq_len(ncol(cov))) {
            cov[, i] <- as.factor(cov[, i])
            cov[, i] <- as.numeric(cov[, i])
        }

        cov <- data.frame(cov)
        rownames(cov) <- rownames(data.frame(colData(se)))
        colnames(cov) <- covar

        limma_corrected <- limma::removeBatchEffect(
            assays(se)[[assay_to_normalize]],
            batch = batch,
            covariates = cov
        )
    }
    assays(se)[[output_assay_name]] <- limma_corrected
    return(se)
    }
