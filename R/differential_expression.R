globalVariables(c("chosen", "P.Value", "adj.P.Val", "effects", "pval"))

#' Differential Expression Analysis
#'
#' This function runs DE analysis on a count matrix (DESeq, edgeR) or a normalized log
#' or log-CPM matrix (limma) contained in the se object
#' @param se SummarizedExperiment object
#' @param method DE analysis method option ('DESeq2', 'limma', or 'edgeR')
#' @param batch metadata column in the se object representing batch
#' @param conditions metadata columns in the se object representing additional
#'   analysis covariates
#' @param assay_to_analyze Assay in the se object (either counts for DESeq2 or edgeR or
#'   normalized data for limma) for DE analysis
#' @return A named list containing the log2FoldChange, pvalue and adjusted
#'   pvalue (padj) for each analysis returned by DESeq2 or limma or edgeR
#' @import SummarizedExperiment
#' @import DESeq2
#' @import scran
#' @import edgeR
#' @importFrom stats model.matrix as.formula t.test aov coef
#' @importFrom limma lmFit eBayes topTable makeContrasts contrasts.fit
#' @examples
#' library(scran)
#' se <- mockSCE()
#' differential_expression <- BatchQC::DE_analyze(se = se,
#'                                                 method = "DESeq2",
#'                                                 batch = "Treatment",
#'                                                 conditions = c(
#'                                                 "Mutation_Status"),
#'                                                 assay_to_analyze = "counts")
#' pval_summary(differential_expression)
#' pval_plotter(differential_expression)
#'
#' @export
DE_analyze <- function(se, method, batch, conditions, assay_to_analyze) {
    data <- assays(se)[[assay_to_analyze]]
    rownames(data) <- names(se)
    analysis_design <- as.data.frame(colData(se)[c(conditions, batch)])
    res <- list()

    if (method == 'DESeq2') {
        # Check if the assay contains counts (e.g. non negative integer data),
        for (item in data){
            if (round(item) != item) {
                stop("Data contains non-integers")
            }else if (item < 0) {
                stop("Data: data contains negative integers")
            }
        }
        colnames(data) <- rownames(analysis_design)
        data[is.na(data)] <- 0
            dds <- DESeqDataSetFromMatrix(countData = data,
                                        colData = analysis_design,
                                        design = stats::as.formula(paste(" ~ ",
                                            paste(colnames(analysis_design),
                                            collapse = "+"))))
        dds <- DESeq(dds)
        for (covar in DESeq2::resultsNames(dds)){
            imp_data <- data.frame("log2FoldChange" =
                    DESeq2::results(dds, name = covar)$log2FoldChange,
                "pvalue" =  DESeq2::results(dds, name = covar)$pvalue,
                "padj" = DESeq2::results(dds, name = covar)$padj,
                row.names = rownames(DESeq2::results(dds, name = covar)))
            res[[covar]] <- imp_data
        }
    }else if (method == 'limma') {
        design <- stats::model.matrix(stats::as.formula(paste(" ~",
            paste(colnames(analysis_design), collapse = "+"))),
            data = analysis_design)

        fit <- limma::lmFit(data, design)
        eBayes_res <- limma::eBayes(fit)

        for (i in seq_len(length(colnames(eBayes_res$coefficients)))){
            results <- limma::topTable(eBayes_res, coef = i, number = Inf) %>%
                select(c(1, P.Value, adj.P.Val))
            colnames(results) <- c("log2FoldChange", "pvalue", "padj" )
            res[[colnames(eBayes_res$coefficients)[[i]]]] <- results
        }
    }else if (method == 'edgeR') {
      # define design matrix
      design <- stats::model.matrix(
        stats::as.formula(
          paste(
            "~", paste(colnames(analysis_design)),
            collapse = "+"
          )
        ),
        data = analysis_design
      )
      
      # Test for significant DE in each gene using QL F-test and get results
      fit <- edgeR::glmQLFit(data, design)
      
      for (i in 1:length(colnames(design))){
        qlf <- edgeR::glmQLFTest(fit, coef = i)
        results <- edgeR::topTags(qlf, n = Inf, adjust.method = "BH")$table %>%
          select(logFC, PValue, FDR)
        colnames(results) <- c("log2FoldChange", "pvalue", "padj" )
        res[[qlf$comparison]] <- results
      }
    }else {
        stop("Please select a method: 'DESeq2', 'limma', or 'edgeR'")
    }
    return(res)
}


#' Returns summary table for p-values of explained variation
#'
#' @param res_list Differential Expression analysis result (a named list of
#'   dataframes corresponding to each analysis completed with a "pvalue" column)
#' @return summary table for p-values of explained variation for each analysis
#'
#' @examples
#' library(scran)
#' se <- mockSCE()
#' differential_expression <- BatchQC::DE_analyze(se = se,
#'                                                 method = "DESeq2",
#'                                                 batch = "Treatment",
#'                                                 conditions = c(
#'                                                 "Mutation_Status"),
#'                                                 assay_to_analyze = "counts")
#' pval_summary(differential_expression)
#'
#' @export
pval_summary <- function(res_list) {

    pval_sum_table <- vector()
    for (res_table in res_list){
        pval_sum_table <- as.data.frame(cbind(pval_sum_table, res_table$pvalue))

    }

    colnames(pval_sum_table) <- names(res_list)
    rownames(pval_sum_table) <- rownames(res_list[[1]])

    return(pval_table = pval_sum_table)
}


#' P-value Plotter
#' This function allows you to plot p-values of explained variation
#' @param DE_results Differential Expression analysis result (a named list of
#' dataframes corresponding to each analysis completed with a "pvalue" column)
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @importFrom data.table data.table
#' @return boxplots of pvalues for each condition
#' @examples
#' library(scran)
#' se <- mockSCE()
#' differential_expression <- BatchQC::DE_analyze(se = se,
#'                                                 method = "DESeq2",
#'                                                 batch = "Treatment",
#'                                                 conditions = c(
#'                                                 "Mutation_Status"),
#'                                                 assay_to_analyze = "counts")
#' pval_summary(differential_expression)
#' pval_plotter(differential_expression)
#'
#' @export
pval_plotter <- function(DE_results) {
    pval_table <- data.frame(row.names = row.names(DE_results[[1]]))
    for (covar in DE_results){
        pval_table <- cbind(pval_table, covar$pvalue)
    }

    colnames(pval_table) <- names(DE_results)

    pval_table <- pivot_longer(pval_table, 2:length(colnames(pval_table)),
        names_to = "effects",
        values_to = "pval")

    covar_boxplot <- ggplot(pval_table,
        aes(x = effects, y = pval, fill = effects)) +
        geom_violin(width = 1.4) +
        geom_boxplot(width = 0.1) +
        scale_x_discrete(name = "") +
        scale_y_continuous(name = "P-Values", limits = c(0, 1)) +
        coord_flip() +
        labs(title =
                "Distribution of Batch and Covariate Effects (P-Values)
                Across Genes") +
        theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
    return(covar_boxplot = covar_boxplot)
}
