globalVariables(c("chosen", "P.Value", "adj.P.Val", "effects", "pval"))

#' Differential Expression Analysis
#'
#' This function runs DE analysis on a count matrix (DESeq), a normalized log (ANOVA), a normalized log
#' or log-CPM matrix (limma), or an edgeR TMM-normalized matrix (edgeR)
#' contained in the se object.
#' @param se SummarizedExperiment object
#' @param method DE analysis method option ('DESeq2', 'limma', 'edgeR', or 'ANOVA')
#' @param batch metadata column in the se object representing batch
#' @param conditions metadata columns in the se object representing additional
#'   analysis covariates
#' @param assay_to_analyze Assay in the se object (either counts for DESeq2 or
#'   normalized data for limma or edgeR) for DE analysis
#' @param padj_method correction method for adjusted p-value from p.adjust.methods
#' @return A named list containing the log2FoldChange, pvalue and adjusted
#'   pvalue (padj) for each analysis returned by DESeq2, limma, edgeR, or ANOVA
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
DE_analyze <- function(se, method, batch, conditions, assay_to_analyze, padj_method) {
    data <- assays(se)[[assay_to_analyze]]
    rownames(data) <- names(se)
    analysis_design <- as.data.frame(colData(se)[c(conditions, batch)])
    res <- list()
    design <- stats::model.matrix(stats::as.formula(paste(" ~",
            paste(colnames(analysis_design), collapse = "+"))),
            data = analysis_design)
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
        fit <- limma::lmFit(data, design)
        eBayes_res <- limma::eBayes(fit)

        for (i in seq_len(length(colnames(eBayes_res$coefficients)))){
            results <- limma::topTable(eBayes_res, coef = i, number = Inf) %>%
                select(c(1, P.Value, adj.P.Val))
            colnames(results) <- c("log2FoldChange", "pvalue", "padj" )
            res[[colnames(eBayes_res$coefficients)[[i]]]] <- results
        }
    }else if (method == 'edgeR') {
        res <- edgeR_DE(data, design)
    }else if (method == 'ANOVA') {
        feature_list <- datatable_DE(se, data, assay_to_analyze, analysis_design)
        res <- anova_DE(feature_list, padj_method, assay_to_analyze, analysis_design)
    }else {
        stop("Please select a method: 'DESeq2', 'limma', or 'edgeR'")
    }
    return(res)
}


edgeR_DE <- function(data, design) {
    fit <- edgeR::glmQLFit(data, design)
    res <- list()

    for (i in seq_len(length(colnames(design)))){
        quasi_likelihood <- edgeR::glmQLFTest(fit, coef = i)
        results <- edgeR::topTags(quasi_likelihood,
            n = Inf, adjust.method = "BH")$table %>%
            select(logFC, PValue, FDR)
        colnames(results) <- c("log2FoldChange", "pvalue", "padj" )
        res[[quasi_likelihood$comparison]] <- results
    }
    return(res)
}

datatable_DE <- function(se, assay_to_analyze, conditions, batch) {
    data <- assays(se)[[assay_to_analyze]]
    features <- rownames(data)
    analysis_design <- as.data.frame(colData(se)[c(conditions, batch)])

    assay_dt <- as.data.table(data, keep.rownames = "features")
    design_dt <- as.data.table(analysis_design, keep.rownames = "samples")

    assay_long <- data.table::melt(assay_dt, id.vars = "features", variable.name = "samples", value.name = assay_to_analyze)

    merged_dt <- assay_long[design_dt, on = "samples"]
    feature_list <- split(merged_dt, by = "features", keep.by = FALSE)
    return(feature_list)
}

anova_DE <- function(se, feature_list, padj_method, assay_to_analyze, batch, conditions) {
    analysis_design <- as.data.frame(colData(se)[c(conditions, batch)])
    model <- stats::as.formula(paste(assay_to_analyze, "~", paste(colnames(analysis_design), collapse = "+")))
    res <- list()
    all_res <- list()
    for (feature in names(feature_list)) {
        feature_dt <- data.table::as.data.table(feature_list[[feature]])
        
        anov_model <- aov(model, data = feature_dt)
        model_summary <- anova(anov_model)
        
        result_vars <- setdiff(rownames(model_summary), "Residuals")
        
        for (var_name in result_vars) {
            var_levels <- as.character(unique(feature_dt[[var_name]]))
            if (length(var_levels) > 1) {
                pval <- model_summary[var_name, "Pr(>F)"]
                fval <- model_summary[var_name, "F value"]
                
                means_dt <- feature_dt[, .(mean_val = mean(get(assay_to_analyze))), by = var_name]
                
                
                ref_level <- var_levels[1]
                ref_mean <- means_dt[get(var_name) == ref_level, mean_val]
                non_ref_results <- means_dt[get(var_name) != ref_level, .(feature = feature,
                                                                       log2FoldChange = log2(mean_val / ref_mean),
                                                                       fvalue = fval,
                                                                       pvalue = pval,
                                                                       var = var_name,
                                                                       reflevel = ref_level,
                                                                       currentlevel = get(var_name))]
                all_res[[length(all_res) + 1]] <- non_ref_results
            }else{
                stop("Each factor needs to have more than two levels!")
            }
        }
    }
    if (length(all_res) > 0) {
        combined_dt <- rbindlist(all_res)
        
        combined_dt[, padj := p.adjust(pvalue, method = padj_method), by = var]
        
        combined_dt[, comparison := paste0(currentlevel, ":", reflevel)]
        combined_dt[, comparison := paste0(var, ":", comparison)]
        
        for (i in unique(combined_dt$comparison)) {
            var_data <- combined_dt[comparison == i]
            var_df <- as.data.frame(var_data[, .(log2FoldChange, fvalue, pvalue, padj)])
            rownames(var_df) <- var_data$feature
            res[[i]] <- var_df
        }
    }
    
    return(res)
}

kw_DE <- function(feature_list, padj_method, assay_to_analyze, batch) {
    res <- list()
    for (n in names(feature_list)){
        feature_df <- feature_list[[n]]
        kw_model <- kruskal.test(assay_to_analyze ~ batch, data = feature_df)

        var_name <- batch
        var_levels <- as.character(unique(feature_df[[var_name]]))
        pval <- model_summary[var_name, "Pr(>F)"]
        if (length(var_levels) > 1) {
            ref_level <- var_levels[1]
            ref_mean <- mean(
                feature_df[assay_to_analyze][feature_df[var_name] == ref_level],
                na.rm = TRUE
            )
            for (j in 2:length(var_levels)) {
                current_level <- var_levels[j]
                current_mean <- mean(
                    feature_df[assay_to_analyze][feature_df[var_name] == current_level],
                    na.rm = TRUE
                )
                
                log2FC <- log2(current_mean / ref_mean)
            }
        }
    }
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

    if ("(Intercept)" %in% colnames(pval_table)) {
        pval_table <- pval_table |>
            select(-"(Intercept)")
    }

    pval_table <- pivot_longer(pval_table, 1:length(colnames(pval_table)),
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

