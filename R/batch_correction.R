#' Batch Correct
#' This function allows you to Add batch corrected count matrix to the SE object
#' @param se SummarizedExperiment object
#' @param method Normalization Method ("ComBat-Seq", "ComBat", "limma", "sva",
#' svaseq)
#' @param assay_to_normalize Which assay use to do normalization
#' @param batch The batch
#' @param group The group variable
#' @param covar list of covariates
#' @param output_assay_name name of results assay
#' @param ... Arguments to be passed to specific methods, such as `num_sv` for
#' `svaseq_correction`
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
    var_of_interest <- batch
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
    } else if (method == 'sva') {
        se <- sva_correction(se, assay_to_normalize, var_of_interest, covar,
            output_assay_name)
    } else if (method == "svaseq") {
        se <- svaseq_correction(se, assay_to_normalize, var_of_interest, covar,
            output_assay_name, num_sv = FALSE)
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

#' sva Correction
#' This function applies sva correction to a summarized experiment object
#' (implementation adapted from sva::psva)
#' @param se SummarizedExperiment object
#' @param assay_to_normalize string; name of assay that should be corrected
#' @param var_of_interest string; name of  experimental variable of interest
#' @param covar list; sting list  of covariates to include in sva analysis
#' @param output_assay_name string; name of results assay
#' @return SE object with an added sva corrected array
#' @import SummarizedExperiment
#' @import sva

sva_correction <- function(se, assay_to_normalize, var_of_interest,
    covar, output_assay_name) {
    if (is.null(covar)) {
        mod0 <- model.matrix(~1, data = colData(se))
        n.sv <- sva::num.sv(assays(se)[[assay_to_normalize]], mod,
            method = "leek")
        sva_assay <- sva::psva(dat = assays(se)[[assay_to_normalize]],
            batch = data.frame(colData(se))[, var_of_interest],
            mod0 = mod0,
            n.sv = n.sv)
    }else {
        str_formula <- ""
        for (i in seq_len(covar)) {
            if (i < length(covar)) {
                str_formula <- paste0(str_formula,
                    "as.factor(", covar[i], ") + ")
            }else {
                str_formula <- paste0(str_formula, "as.factor(", covar[i], ")")
            }
        }
        form_mod0 <- formula(paste0("~", str_formula))
        mod0 <- model.matrix(form_mod0, data = colData(se))
        mod_formula <- formula(paste0("~as.factor(", var_of_interest, ") + ",
            str_formula))
        mod <- model.matrix(mod_formula, data = colData(se))

        psva.SV <- sva::sva(assays(se)[[assay_to_normalize]], mod = mod,
            mod0 = mod0) #, n.sv = n.sv)
        colnames(psva.SV$sv) <- paste('sv', seq_len(ncol(psva.SV$sv)))
        psva.fit <- lmFit(assays(se)[[assay_to_normalize]],
            cbind(mod, psva.SV$sv))
        sva_assay <- sweep(psva.fit$coefficients[,
            paste('sv', seq_len(ncol(psva.SV$sv)))] %*% t(psva.SV$sv), 1,
            psva.fit$coefficients[, "(Intercept)"], FUN = "+")
    }
    colnames(sva_assay) <- colnames(assays(se)[[assay_to_normalize]])
    assays(se)[[output_assay_name]] <- sva_assay
    return(se)
}

#' svaseq Correction
#' This function applies sva correction to a summarized experiment object
#' with count based RNA-seq data
#' @param se SummarizedExperiment object
#' @param assay_to_normalize string; name of assay that should be corrected
#' @param var_of_interest string; name of  experimental variable of interest
#' @param covar list; sting list  of covariates to include in sva analysis
#' @param output_assay_name string; name of results assay
#' @param num_sv boolean; Default is FALSE: the number of estimated latent
#' factor is set to 1 for a small number of samples. If set to TRUE, svaseq
#' function will estimate the number of latent factors for you.
#' @return SE object with an added sva corrected array
#' @import SummarizedExperiment
#' @import sva

svaseq_correction <- function(se, assay_to_normalize, var_of_interest,
    covar, output_assay_name, num_sv = FALSE) {
    dat <- assays(se)[[assay_to_normalize]]
    if (is.null(covar)) {
        mod0 <- model.matrix(~1, data = colData(se))
        mod1 <- model.matrix(~as.factor(get(var_of_interest)),
                            data = colData(se)
                            )
    } else {
        mod1 <- model.matrix(
            ~as.factor(get(var_of_interest)) + as.factor(get(covar)),
            data = colData(se))
        mod0 <- model.matrix(
            ~as.factor(get(covar)),
            data = colData(se)
            )
    }
    if (num_sv){
        batch_unsup_sva <- svaseq(dat, mod1, mod0)$sv
    } else {
        batch_unsup_sva <- svaseq(dat, mod1, mod0, n.sv = 1)$sv
    }
    colnames(batch_unsup_sva) <- paste('sv', seq_len(ncol(batch_unsup_sva)))
    mod1Sv <- cbind(mod1, batch_unsup_sva)
    psva.fit <- lmFit(dat, mod1Sv)
    
    sv_coef <- psva.fit$coefficients[, colnames(batch_unsup_sva)]
    sv_effects <- sv_coef %*% t(batch_unsup_sva)
    sva_assay <- dat - sv_effects
    assays(se)[[output_assay_name]] <- svaseq_assay
    return(se)
}