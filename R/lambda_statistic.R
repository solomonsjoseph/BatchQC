#' Compute the lambda index for determining a need for batch correction
#'
#' This function calculates the proportions of variation explained by batch,
#' group, and residual for each gene using two-way ANOVA and computes the lambda
#' index based on these three proportions.
#'
#' @param dat Numeric matrix of dimension (genes x samples) where each row
#'   represents one gene's expression across samples.
#' @param batchind Factor or numeric vector of length = ncol(dat); batch
#'   indicator for each sample
#' @param groupind Factor or numeric vector of length = ncol(dat); biological
#'   group label/indicator for each sample.
#' @importFrom stats lm anova
#' @return Named numeric vector with elements:
#'   \describe{
#'     \item{BatchV}{Proportion of total variance explained by batch effects.}
#'     \item{GroupV}{Proportion of total variance explained by group effects.}
#'     \item{ResidV}{Proportion of total variance that is residual noise.}
#'     \item{lambda_raw}{Raw lambda index = total SS_batch / total SS_group.}
#'     \item{lambda_adj}{Adjusted lambda = lambda_raw * ResidV/(1-ResidV).}
#'   }
#' @examples
#'
#' library(scran)
#' se <- mockSCE()
#' res <- BatchQC::compute_lambda(assays(se)[["counts"]],
#'   colData(se)$Mutation_Status,
#'   colData(se)$Treatment)
#' print(res)
#'
#' @export

compute_lambda <- function(dat, batchind, groupind) {
    dat <- t(dat)
    n <- ncol(dat)
    expvar <- t(apply(dat, 2, function(x) {
        fit <- lm(x ~ factor(batchind) + factor(groupind))
        av  <- anova(fit)
        av$`Sum Sq`
    }))

    ss <- colSums(expvar)
    ss_list <- ss / sum(ss)
    lambda_raw <- ss[1] / ss[2]
    lambda_adj <- lambda_raw * ss_list[3] / (1 - ss_list[3])

    res <- c(BatchV = ss_list[1], GroupV = ss_list[2], ResidV = ss_list[3],
        lambda_raw = lambda_raw, lambda_adj = lambda_adj)

    res <- t(data.frame(res))
    return(res)
}
