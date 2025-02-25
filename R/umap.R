#' Create a umap plot; wrapper function for umap package pplus custom plotting
#' @param se_object se_object; containing data of interest
#' @param assay_of_interest string; the assay in the se_object to plot
#' @param batch string; representing batch
#' @param neighbors integer; number of nearest neighbors, default 15 per umap;
#'   lower values prioritize local structure, higher values will represent
#'   bigger picture but lose finer details
#' @param min_distance numeric; how close points appear in final layout; higher
#'   values puts less emphasis on global structure
#' @import umap
#' @return umap plot
#' @examples
#' library(scran)
#' se <- mockSCE()
#' se$Treatment <- as.factor(se$Treatment)
#' se$Mutation_Status <- as.factor(se$Mutation_Status)
#' umap_plot <- BatchQC::umap(se_object = se, assay_of_interest = "counts",
#' batch = "Treatment")
#' umap_plot
#'
#' @export

umap <- function(se_object, assay_of_interest, batch, neighbors = 15,
    min_distance = 0.1) {
    plot_data <- t(assays(se_object)[[assay_of_interest]])
    umap_data <- umap::umap(plot_data, n_neighbors = neighbors,
        min_dist = min_distance)

    df <- data.frame(x = umap_data$layout[, 1],
                    y = umap_data$layout[, 2],
                    batch = colData(se_object)[[batch]])

    plot <- ggplot(df, aes(x, y, color = batch)) +
        geom_point()

    return(plot)
}
