#' kBET rejection rate
#' This function runs the k-nearest neighbor batch effect test (kBET) to 
#' evaluate whether the data has detectable batch effect.
#'
#' @param se SummarizedExperiment object
#' @param batch The variable that represents batch
#' @param k0 number of nearest neighbors to test on (neighborhood size)
#' @param knn n x k matrix of nearest neighbors for each cell (optional)
#' @param testSize number of data points to test
#' @param do.pca if TRUE, perform a pca prior to knn search
#' (defaults to TRUE)
#' @param dim.pca if do.pca=TRUE, choose the number of dimensions to consider
#' (defaults to 50)
#' @param heuristic heuristic compute an optimal neighborhood size k
#' (defaults to TRUE)
#' @param n_repeat to create a statistics on batch estimates, evaluate 
#' 'n_repeat' subsets
#' @param alpha significance level
#' @param addTest perform an LRT-approximation to the multinomial test AND a 
#' multinomial exact test (if appropriate)
#' @param verbose if TRUE, display stages of current computation
#' (defaults to FALSE)
#' @param adapt In some cases, a number of cells do not contribute
#' to any neighbourhood and this may cause an imbalance in observed and 
#' expected batch label frequencies.
#' 
#' @import SummarizedExperiment
#' @import kBET
#' @import ggplot2
#' 
#' @export

run_kBET <- function(se, assay_to_normalize, batch, k0 = NULL, knn = NULL, testSize = NULL,
                     do.pca = TRUE, dim.pca = 50, heuristic = TRUE,
                     n_repeat = 100, alpha = 0.05, addTest = FALSE,
                     verbose = FALSE, adapt = TRUE){
  
  # run kBET
  batch.estimate <- kBET::kBET(df = as.matrix(assays(se)[[assay_to_normalize]]),
                               batch = data.frame(colData(se))[, batch],
                               k0 = k0, knn = knn,
                               testSize = testSize, do.pca = do.pca,
                               dim.pca = dim.pca, heuristic = heuristic,
                               n_repeat = n_repeat, alpha = alpha,
                               addTest = addTest, verbose = verbose,
                               plot = FALSE, adapt = adapt)
  
  return(batch.estimate)
}
plot_kBET <- function(kBET_res){
  # create ggplot object for plotting kBET's rejection rate
  plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                    each=length(kBET_res$stats$kBET.observed)), 
                          data =  c(kBET_res$stats$kBET.observed,
                                    kBET_res$stats$kBET.expected))
  g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
    labs(x='Test', y='Rejection rate',title='kBET test results') +
    theme_bw() +  
    scale_y_continuous(limits=c(0,1))
  
  return(g)
}