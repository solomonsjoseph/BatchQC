globalVariables(c(glm.nb, AIC, glm, gaussian, glm.control))
#' Compute the AIC for lognormal (ComBat) model, negative binomial (ComBat-seq)
#' model and the Voom model
#'
#' This function calculates the AIC based on lognormal distribution,
#' negative-binomial distribution as well as the Voom transformation.
#' It then compares the AICs of the three models across different genes.
#'
#' @param se SummarizedExperiment object
#' @param assay_of_interest The assay name from se that you are interested in
#' analyzing. This assay need to be a counts assay containing only non-negative
#' integers.
#' @param batchind Factor or numeric vector of length = ncol(dat);
#' batch indicator for each sample.
#' @param groupind Factor or numeric vector of length = ncol(dat);
#' biological group label/indicator for each sample.
#' @param maxit Integer giving the maximal number of IWLS iterations. Default is
#' 25.
#' @param zero_filt_percent Numeric value between 0 and 100, the percentage of
#' zeros allowed for each gene to be included in the AIC calculation. Genes with
#' more than this percentage of zeros will be filtered out. Default is 100.
#'
#' @description
#'   \describe{
#'     \item{nb_result}{A vector contains the AIC based on negative binomial
#'     model for individual genes.}
#'     \item{lognormal_result}{A vector contains the AIC based on lognormal
#'     model for individual genes.}
#'     \item{voom_result}{A vector contains the AIC based on voom transformation
#'     for individual genes.}
#'     \item{total_AIC}{The sum of AICs across all genes for the three models
#'     in comparison.}
#'     \item{min_AIC}{The number of minimum AIC across the three models in
#'     comparison for individual genes.}
#'   }
#'
#' @return A list with the following two elements:
#' \describe{
#'     \item{total_AIC}{The sum of AICs across all genes for the three models
#'     in comparison.}
#'     \item{min_AIC}{The number of minimum AIC across the three models in
#'     comparison for individual genes.}
#'   }
#' @examples
#' library(scran)
#' se <- mockSCE()
#' compare_aic <- compute_aic(se, assay_of_interest = "counts",
#'                             batchind = "Cell_Cycle", groupind = "Treatment")
#' print(compare_aic["total_AIC"])
#' print(compare_aic["min_AIC"])
#'
#' @importFrom limma voom
#' @importFrom stats AIC gaussian glm glm.control
#' @importFrom MASS glm.nb
#' @import SummarizedExperiment
#' @export
compute_aic <- function(se, assay_of_interest, batchind,
                        groupind, maxit = 25, zero_filt_percent = 100) {
    dat <- assays(se)[[assay_of_interest]]
    batchind <- as.factor(colData(se)[[batchind]])
    groupind <- as.factor(colData(se)[[groupind]])
    if (!all(dat == floor(dat)) || any(dat < 0)) {
        stop("Counts must be non-negative integers only.")
    }
    dat <- dat[rowSums(dat) != 0, ]
    dat <- dat[rowSums(dat == 0) <= zero_filt_percent / 100 * ncol(dat), ]
    nb_result <- apply(dat, 1, function(x) {
        tryCatch({
        nb_model <- glm.nb(x ~ 1, control = glm.control(maxit = maxit))
        nb_AIC <- AIC(nb_model)
        return(nb_AIC)
    }, error = function(e) {
        return(NA)
    })})
    lognormal_result <- apply(dat, 1, function(x) {
        tryCatch({
        lognormal_model <- glm(log(x + 1e-100) ~ 1, family = gaussian)
        lognormal_AIC <- AIC(lognormal_model)
        return(lognormal_AIC)
    }, error = function(e) {
        return(NA)
    })})
    voom_dat <- voom(dat, design = model.matrix(~ groupind + batchind))
    voom_dat <- voom_dat$E
    voom_result <- apply(voom_dat, 1, function(x) {
        tryCatch({
        voom_lm_model <- lm(x ~ 1)
        voom_lm_AIC <- AIC(voom_lm_model)
        return(voom_lm_AIC)
    }, error = function(e) {
        return(NA)
    })})
    aic_matrix <- cbind(nb_result, lognormal_result, voom_result)
    colnames(aic_matrix) <- c("NB_AIC", "Lognormal_AIC", "Voom_AIC")
    total_AIC <- colSums(aic_matrix, na.rm = TRUE)
    min_model <- apply(aic_matrix, 1, function(x) {
        if (all(is.na(x))) return(NA)
        which.min(x)
    })
    min_AIC <- table(factor(min_model, levels = seq_along(seq(1, 3))))
    names(min_AIC) <- c("NB", "Lognormal", "Voom")
    return(list(total_AIC = total_AIC, min_AIC = min_AIC))
}