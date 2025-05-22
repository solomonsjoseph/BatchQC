#' kBET rejection rate
#' This function runs the k-nearest neighbor batch effect test (kBET) to
#' evaluate whether the data has detectable batch effect.
#'
#' @param se SummarizedExperiment object
#' @param batch character string of column name that represents batch
#' @param k0 integer representing number of nearest neighbors to test on
#' (neighborhood size)
#' @param knn n x k matrix of nearest neighbors for each cell (optional)
#' @param testSize integer representing number of data points to test
#' @param do.pca Boolean, if TRUE, perform a pca prior to knn search
#' (defaults to TRUE)
#' @param dim.pca Boolean, if do.pca=TRUE, choose the number of dimensions to
#' consider (defaults to 50)
#' @param heuristic Boolean, if true, compute an optimal neighborhood size k
#' (defaults to TRUE)
#' @param n_repeat numeric representing 'n_repeat' subsets to evaluate in order
#' to create a statistics on batch estimates
#' @param alpha numeric for significance level
#' @param addTest Boolean, if TRUE, perform an LRT-approximation to the
#' multinomial test AND a multinomial exact test (if appropriate)
#' @param verbose Boolean, if TRUE, display stages of current computation
#' (defaults to FALSE)
#' @param adapt Boolean, if TRUE, frequencies will be adapted (defaults to TRUE)
#'
#' @returns list object from kBET() function
#'    \enumerate{
#'    \item \code{summary} - a rejection rate for the data,
#'         an expected rejection rate for random
#'         labeling and the significance for the observed result
#'    \item \code{results} - detailed list for each tested cells;
#'         p-values for expected and observed label distribution
#'    \item \code{average.pval} - significance level over the averaged
#'         batch label distribution in all neighbourhoods
#'    \item \code{stats} - extended test summary for every sample
#'    \item \code{params} - list of input parameters and adapted parameters,
#'    respectively
#'    \item \code{outsider} - only shown if \code{adapt=TRUE}. List of samples
#'         without mutual nearest neighbour: \itemize{
#'     \item \code{index} - index of each outsider sample)
#'     \item \code{categories} - tabularised labels of outsiders
#'     \item \code{p.val} - Significance level of outsider batch label distribution
#'         vs expected frequencies.
#'     If the significance level is lower than \code{alpha},
#'     expected frequencies will be adapted}
#'    }
#' @import SummarizedExperiment
#' @import kBET
#' @import ggplot2
#'
#' @examples
#' library(scran)
#' se <- mockSCE()
#' kBET_result <- BatchQC::run_kBET(
#'   se = se,
#'   assay_to_normalize = "counts",
#'   batch = "Treatment"
#' )
#'
#' BatchQC::plot_kBET(kBET_result)
#'
#' @export
run_kBET <- function(se, assay_to_normalize, batch, k0 = NULL, knn = NULL,
                     testSize = NULL, do.pca = TRUE, dim.pca = 50,
                     heuristic = TRUE, n_repeat = 100, alpha = 0.05,
                     addTest = FALSE, verbose = FALSE, adapt = TRUE) {
  # run kBET
  batch.estimate <- kBET::kBET(
    df = as.matrix(assays(se)[[assay_to_normalize]]),
    batch = data.frame(colData(se))[, batch],
    k0 = k0, knn = knn,
    testSize = testSize, do.pca = do.pca,
    dim.pca = dim.pca, heuristic = heuristic,
    n_repeat = n_repeat, alpha = alpha,
    addTest = addTest, verbose = verbose,
    plot = FALSE, adapt = adapt
  )

  return(batch.estimate)
}


#' kBET Rejection Plotter
#' This function generates a boxplot of observed and expected rejection rates
#' for the provided kBET output list object
#'
#' @param kBET_res list object output from kBET function
#'
#' @returns ggplot object containing kBET rejection boxplot
#' @export
plot_kBET <- function(kBET_res) {
  # create ggplot object for plotting kBET's rejection rate
  plot.data <- data.frame(
    class = rep(c("observed", "expected"),
      each = length(kBET_res$stats$kBET.observed)
    ),
    data = c(
      kBET_res$stats$kBET.observed,
      kBET_res$stats$kBET.expected
    )
  )
  g <- ggplot(plot.data, aes(class, data)) +
    geom_boxplot() +
    labs(x = "Test", y = "Rejection rate", title = "kBET test results") +
    theme_bw() +
    scale_y_continuous(limits = c(0, 1))

  return(g)
}
