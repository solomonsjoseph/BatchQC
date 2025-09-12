#' @importFrom stats pchisq pnorm
# a wrapper for kBET to fix a neighbourhood size
scan_nb <- function(x, df, batch, knn) {
    res <- kBET(
        df = df, batch = batch, k0 = x, knn = knn, testSize = NULL,
        heuristic = FALSE, n_repeat = 10, alpha = 0.05,
        addTest = FALSE, plot = FALSE, verbose = FALSE, adapt = FALSE
    )
    result <- res$summary
    result$kBET.observed[1]
}

# the residual score function of kBET
residual_score_batch <- function(knn.set, class.freq, batch) {
    # knn.set: indices of nearest neighbours
    # empirical frequencies in nn-environment (sample 1)
    # ignore NA entries (which may arise from subsampling a knn-graph)
    if (all(is.na(knn.set))) { # if all values of a neighbourhood are NA
        return(NA)
    } else {
        valid.knn <- knn.set[!is.na(knn.set)]
        freq.env <- table(batch[valid.knn]) / length(valid.knn)
        full.classes <- rep(0, length(class.freq$class))
        full.classes[class.freq$class %in% names(freq.env)] <- freq.env
        exp.freqs <- class.freq$freq
        # compute chi-square test statistics
        sum((full.classes - exp.freqs)^2 / exp.freqs)
    }
}

# which batch has the largest deviance (and is underrepresented)
max_deviance_batch <- function(knn.set, class.freq, batch) {
    # knn.set: indices of nearest neighbours
    # empirical frequencies in nn-environment (sample 1)
    if (all(is.na(knn.set))) { # if all values of a neighbourhood are NA
        return(NA)
    } else {
        valid.knn <- knn.set[!is.na(knn.set)]
        freq.env <- table(batch[valid.knn]) / length(valid.knn)
        full.classes <- rep(0, length(class.freq$class))
        full.classes[class.freq$class %in% names(freq.env)] <- freq.env
        exp.freqs <- class.freq$freq
        # compute chi-square test statistics
        allScores <- (full.classes - exp.freqs) / exp.freqs
        batch[which(allScores == min(allScores))]
    }
}


# the core function of kBET
chi_batch_test <- function(knn.set, class.freq, batch, df) {
    # knn.set: indices of nearest neighbours
    # empirical frequencies in nn-environment (sample 1)
    if (all(is.na(knn.set))) { # if all values of a neighbourhood are NA
        return(NA)
    } else {
        freq.env <- table(batch[knn.set[!is.na(knn.set)]])
        full.classes <- rep(0, length(class.freq$class))
        full.classes[class.freq$class %in% names(freq.env)] <- freq.env
        exp.freqs <- class.freq$freq * length(knn.set)
        # compute chi-square test statistics
        chi.sq.value <- sum((full.classes - exp.freqs)^2 / exp.freqs)
        result <- 1 - pchisq(chi.sq.value, df) # p-value for the result
        if (is.na(result)) { # I actually would like to now when 'NA' arises.
            return(NA)
        } else {
            result
        }
    }
}

lrt_approximation <- function(knn.set, class.freq, batch, df) {
    # knn.set: indices of nearest neighbours
    # empirical frequencies in nn-environment (sample 1)
    if (all(is.na(knn.set))) { # if all values of a neighbourhood are NA
        return(NA)
    } else {
        # observed realisations of each category
        obs.env <- table(batch[knn.set[!is.na(knn.set)]])
        # observed 'probabilities'
        freq.env <- obs.env / sum(obs.env)
        full.classes <- rep(0, length(class.freq$class))
        obs.classes <- class.freq$class %in% names(freq.env)
        # for stability issues (to avoid the secret division by 0): introduce
        # another alternative model where the observed probability
        # is either the empirical frequency or 1/(sample size) at minimum
        if (length(full.classes) > sum(obs.classes)) {
            dummy.count <- length(full.classes) - sum(obs.classes)
            full.classes[obs.classes] <- obs.env / (sum(obs.env) + dummy.count)
            pmin <- 1 / (sum(obs.env) + dummy.count)
            full.classes[!obs.classes] <- pmin
        } else {
            full.classes[obs.classes] <- freq.env
        }
        exp.freqs <- class.freq$freq # expected 'probabilities'
        # compute likelihood ratio of null and alternative hypothesis,
        # test statistics converges to chi-square distribution
        full.obs <- rep(0, length(class.freq$class))
        full.obs[obs.classes] <- obs.env

        lrt.value <- -2 * sum(full.obs * log(exp.freqs / full.classes))

        result <- 1 - pchisq(lrt.value, df) # p-value for the result
        if (is.na(result)) { # I actually would like to now when 'NA' arises.
            return(NA)
        } else {
            result
        }
    }
}

# truncated normal distribution distribution function
ptnorm <- function(x, mu, sd, a = 0, b = 1, alpha = 0.05, verbose = FALSE) {
    # this is the cumulative density of the truncated normal distribution
    # x ~ N(mu, sd^2), but we condition on a <= x <= b
    if (!is.na(x)) {
        if (a > b) {
            warning("Lower and upper bound are interchanged.")
            tmp <- a
            a <- b
            b <- tmp
        }

        if (sd <= 0 || is.na(sd)) {
            if (verbose) {
                warning("Standard deviation must be positive.")
            }
            if (alpha <= 0) {
                stop("False positive rate alpha must be positive.")
            }
            sd <- alpha
        }
        if (x < a || x > b) {
            warning("x out of bounds.")
            cdf <- as.numeric(x > a)
        } else {
            alp <- pnorm((a - mu) / sd)
            bet <- pnorm((b - mu) / sd)
            zet <- pnorm((x - mu) / sd)
            cdf <- (zet - alp) / (bet - alp)
        }
        cdf
    } else {
        return(NA)
    }
}

# from EMT package (function is needed to run ExactMultinomialTest.R)
findVectors <- function(groups, size) {
    if (groups == 1) {
        mat <- size
    } else {
        mat <- matrix(rep(0, groups - 1), nrow = 1)
        for (i in seq_len(size)) {
            mat <- rbind(mat, findVectors(groups - 1, i))
        }
        mat <- cbind(mat, size - rowSums(mat))
    }
    invisible(mat)
}

# ExactMultinomialTest
# @description from EMT package (solemnly adapted to stay quiet)
#
# @param observed
# @param prob
# @param size
# @param groups
# @param numEvents
# @param verbose
#' @importFrom stats dmultinom
ExactMultinomialTest <- function(
    observed, prob, size,
    groups, numEvents, verbose) {
    pObs <- stats::dmultinom(observed, size = size, prob)
    eventMat <- findVectors(groups, size)
    if (nrow(eventMat) != numEvents) {
        stop("Wrong number of events calculated. \n This is probably a bug.")
    }
    eventProb <- apply(eventMat, 1, function(x) {
        dmultinom(x,
            size = size, prob = prob
        )
    })
    p.value <- sum(eventProb[eventProb <= pObs])
    if (round(sum(eventProb), digits = 2) != 1) {
        stop("Wrong values for probabilities. \n This is probably a bug.")
    }


    if (verbose) {
        head <- paste("\n Exact Multinomial Test, distance measure: p\n\n")
        tab <- as.data.frame(cbind(
            numEvents, round(pObs, digits = 4),
            round(p.value, digits = 4)
        ))
        colnames(tab) <- c("   Events", "   pObs", "   p.value")
        warning(head)
        print(tab, row.names = FALSE)
    }

    invisible(list(
        id = "Exact Multinomial Test", size = size,
        groups = groups, stat = "lowP", allProb = sort(eventProb,
            decreasing = TRUE
        ), ntrial = NULL, p.value = round(p.value,
            digits = 4
        )
    ))
}

# from EMT package (function needed to run ExactMultinomialTest)

multinomial.test <- function(
    observed, prob, useChisq = FALSE,
    MonteCarlo = FALSE, ntrial = 1e+05, atOnce = 1e+06, verbose = FALSE) {
    if (!is.vector(observed, mode = "numeric")) {
        stop(
            " Observations have to be stored in a vector, ",
            "e.g.  'observed <- c(5,2,1)'"
        )
    }
    if (!is.vector(prob, mode = "numeric")) {
        stop(
            " Probabilities have to be stored in a vector, ",
            "e.g.  'prob <- c(0.25, 0.5, 0.25)'"
        )
    }
    if (round(sum(prob), digits = 1) != 1) {
        stop("Wrong input: sum of probabilities must not deviate from 1.")
    }
    if (length(observed) != length(prob)) {
        stop(" Observations and probabilities must have same dimensions.")
    }
    size <- sum(observed)
    groups <- length(observed)
    numEvents <- choose(size + groups - 1, groups - 1)
    res <- ExactMultinomialTest(
        observed, prob, size,
        groups, numEvents, verbose
    )

    invisible(res)
}

# wrapper for the multinomial exact test function
multiNom <- function(x, y, z) {
    z.f <- factor(z)
    tmp <- multinomial.test(as.numeric(table(z.f[x])), y)
    tmp$p.value
}

# significance test for pcRegression (two levels)
correlate.fun_two <- function(rot.data, batch, batch.levels) {
    # rot.data: some vector (numeric entries)
    # batch: some vector (categoric entries)
    a <- stats::lm(rot.data ~ batch)
    result <- numeric(2)
    result[1] <- summary(a)$r.squared # coefficient of determination
    result[2] <- summary(a)$coefficients[2, 4] # p-value (significance level)
    t.test.result <- t.test(rot.data[batch == batch.levels[1]],
        rot.data[batch == batch.levels[2]],
        paired = FALSE
    )
    result[3] <- t.test.result$p.value
    result
}

# significance test for pcRegression (more than two levels)
correlate.fun_gen <- function(rot.data, batch) {
    # rot.data: some vector (numeric covariate)
    # batch: some vector (categoric covariate)
    a <- stats::lm(rot.data ~ batch)
    result <- numeric(2)
    # coefficient of determination
    result[1] <- summary(a)$r.squared
    F.test.result <- aov(rot.data ~ batch)
    F.test.summary <- summary(F.test.result)

    # p-value (significance level)
    result[2] <- summary(a)$coefficients[2, 4]
    # p-value of the one-way anova test
    result[3] <- F.test.summary[[1]]$"Pr(>F)"[1]

    result
}

initialize.kbet <- function(
    df, batch, k0, knn, testSize, do.pca, dim.pca, heuristic, n_repeat,
    alpha = 0.05, addTest, verbose, adapt) {
    dof <- length(unique(batch)) - 1 # degrees of freedom
    if (is.factor(batch)) batch <- droplevels(batch)
    frequencies <- table(batch) / length(batch)
    batch.shuff <- replicate(3, batch[sample.int(length(batch))])
    class.frequency <- data.frame(
        class = names(frequencies), freq = as.numeric(frequencies))
    inputs <- validate.inputs(df, batch, verbose)
    if (!inputs$valid) return(NA)
    stopifnot(is(n_repeat, "numeric"), n_repeat > 0)
    k0_info <- determine.k0(k0, heuristic, class.frequency, knn,
        inputs$dim.dataset, verbose)
    if (is.null(k0_info$k0)) return(NA)
    k0 <- k0_info$k0
    do_heuristic <- k0_info$do_heuristic
    knn <- knn %||% find.knn(inputs$dataset, do.pca, k0, verbose, dim.pca,
        inputs$dim.dataset) # find KNNs
    knn <- knn.graph.backward.compatibility(knn, verbose)
    testSize <- set.number.tests(testSize, inputs$dim.dataset, verbose)
    if (adapt) {  # decide to adapt general frequencies
        adapt.freq.res <- adapt.freq(inputs$dim.dataset, knn, k0, inputs$batch,
            dof, class.frequency, alpha, verbose)
    }
    is.imbalanced <- if (adapt) adapt.freq.res$is.imbalanced else NULL
    new.class.frequency <- if (adapt) adapt.freq.res$new.class.freq else NULL
    outsider <- if (adapt) adapt.freq.res$outsider else NULL
    p.out <- if (adapt) adapt.freq.res$p.out else NULL
    if (do_heuristic) {
        k0 <- calculate_nb_size(scan_nb, k0, inputs$dataset,
            inputs$batch, knn, verbose)
    }
    rejection <- initialize.result.list(inputs$dim.dataset)
    # get average residual score
    env <- as.vector(cbind(
        knn[, seq_len(k0 - 1)], seq_len(inputs$dim.dataset[1])))
    cf <- if (adapt && is.imbalanced) new.class.frequency else class.frequency
    score <- k0 * residual_score_batch(env, cf, batch)
    rejection$average.pval <- 1 - pchisq(score, dof)
    return(list(
        dim.dataset = inputs$dim.dataset, testSize = testSize, k0 = k0,
        is.imbalanced = is.imbalanced, batch = inputs$batch, knn = knn,
        new.class.frequency = new.class.frequency, outsider = outsider,
        class.frequency = class.frequency, batch.shuff = batch.shuff,
        kBET.expected = numeric(n_repeat), kBET.observed = numeric(n_repeat),
        kBET.signif = numeric(n_repeat), rejection = rejection, dof = dof,
        p.out = p.out))
}

validate.inputs <- function(df, batch, verbose) {
    dim.dataset <- dim(df)
    # check the feasibility of data input
    if (dim.dataset[1] != length(batch) && dim.dataset[2] != length(batch)) {
        msg <- paste(
            "Input matrix and batch information do not match.",
            "Execution halted.",
            sep = " "
        )
        stop(msg)
    }

    if (dim.dataset[2] == length(batch) && dim.dataset[1] != length(batch)) {
        if (verbose) {
            warning(
                "Input matrix has samples as columns. ",
                "kBET needs samples as rows. Transposing...\n"
            )
        }
        df <- t(df)
        dim.dataset <- dim(df)
    }

    # check if the dataset is too small per se
    if (dim.dataset[1] <= 10) {
        if (verbose) {
            warning(
                "Your dataset has less than 10 samples.",
                "Abort and return NA.\n"
            )
        }
        return(list(valid = FALSE))
    }

    list(dataset = df, batch = batch, dim.dataset = dim.dataset, valid = TRUE)
}

determine.k0 <- function(k0, heuristic, class.freq, knn, dim.dataset, verbose) {
    do_heuristic <- FALSE
    if (is.null(k0) || k0 >= dim.dataset[1]) {
        do_heuristic <- heuristic
        k0 <- floor(mean(class.freq$freq) * dim.dataset[1] *
            ifelse(heuristic, 0.75, 0.25))

        if (verbose) {
            warning(paste0("Initial neighborhood size set to ", k0, ".\n"))
        }
    }

    if (k0 < 10 & (heuristic | is.null(knn))) {
        if (verbose) {
            warning(
                "Your dataset has too few samples to run a heuristic.\n",
                "Return NA.\n",
                "Please assign k0 and set heuristic = FALSE."
            )
        }
        return(NA)
    }

    list(k0 = k0, do_heuristic = do_heuristic)
}

find.knn <- function(dataset, do.pca, k0, verbose, dim.pca, dim.dataset) {
    if (!do.pca) {
        if (verbose) {
            warning("finding knns...")
            tic <- proc.time()
        }
        # use the nearest neighbour index directly for further use in the
        # package
        knn <- get.knn(dataset, k = k0, algorithm = "cover_tree")$nn.index
    } else {
        dim.comp <- min(dim.pca, dim.dataset[2])
        if (verbose) {
            warning("reducing dimensions with svd first...\n")
        }
        data.pca <- svd(x = dataset, nu = dim.comp, nv = 0)
        if (verbose) {
            warning("finding knns...")
            tic <- proc.time()
        }
        knn <- get.knn(data.pca$u, k = k0, algorithm = "cover_tree")
    }
    if (verbose) {
        warning("done. Time:\n")
        warning(proc.time() - tic)
    }

    return(knn)
}

knn.graph.backward.compatibility <- function(knn, verbose) {
    if (is(knn, "list")) {
        knn <- knn$nn.index
        if (verbose) {
            msg <- "KNN input is a list, extracting nearest neighbour index.\n"
            warning(msg)
        }
    }
    return(knn)
}

set.number.tests <- function(testSize, dim.dataset, verbose) {
    if (is.null(testSize) ||
        (floor(testSize) < 1 || dim.dataset[1] < testSize)) {
        test.frac <- 0.1
        testSize <- ceiling(dim.dataset[1] * test.frac)
        if (testSize < 25 && dim.dataset[1] > 25) {
            testSize <- 25
        }
        if (verbose) {
            warning(
                "Number of kBET tests is set to ",
                testSize,
                ".\n",
                sep = ""
            )
        }
    }

    return(testSize)
}

initialize.result.list <- function(dim.dataset) {
    rejection <- list()
    rejection$summary <- data.frame(kBET.expected = numeric(4),
                                    kBET.observed = numeric(4),
                                    kBET.signif = numeric(4))

    rejection$results <- data.frame(tested = numeric(dim.dataset[1]),
                                    kBET.pvalue.test = rep(0, dim.dataset[1]),
                                    kBET.pvalue.null = rep(0, dim.dataset[1]))

    return(rejection)
}

adapt.freq <- function(
    dim.dataset, knn, k0, batch, dof, class.frequency, alpha, verbose) {
    outsider <- which(!(seq_len(dim.dataset[1]) %in% knn[, seq_len(k0 - 1)]))
    is.imbalanced <- FALSE # initialisation
    new.class.frequency <- NULL
    p.out <- 1
    if (length(outsider) > 0) {
        p.out <- chi_batch_test(outsider, class.frequency, batch, dof)
        if (!is.na(p.out)) {
            is.imbalanced <- p.out < alpha
            if (is.imbalanced) {
                new.frequencies <- table(batch[-outsider]) /
                    length(batch[-outsider])
                new.class.frequency <- data.frame(
                    class = names(new.frequencies),
                    freq = as.numeric(new.frequencies)
                )
                if (verbose) {
                    percent <- length(outsider) / length(batch)
                    outs_percent <- round(percent * 100, 3)
                    warning(paste(
                        sprintf(
                            paste0(
                                "There are %s cells (%s%%) that do ",
                                "not appear in any neighbourhood."
                            ),
                            length(outsider), outs_percent
                        ),
                        paste0(
                            "The expected frequencies for each category ",
                            "have been adapted."
                        ),
                        "Cell indexes are saved to result list.",
                        "", sep = "\n"
                    ))
                }
            } else {
                if (verbose) warning(paste0("No outsiders found."))
            }
        } else {
            if (verbose) warning(paste0("No outsiders found."))
        }
    }
    return(list(
        is.imbalanced = is.imbalanced, outsider = outsider, p.out = p.out,
        new.class.freq = new.class.frequency
    ))
}

calculate_nb_size <- function(scan_nb, k0, dataset, batch, knn, verbose) {
    if (verbose) {
        warning("Determining optimal neighbourhood size ...")
    }
    opt.k <- bisect(scan_nb,
        bounds = c(10, k0), known = NULL, dataset, batch,
        knn
    )
    # result
    if (length(opt.k) > 1) {
        k0 <- opt.k[2]
        if (verbose) {
            warning(paste0(
                "done.\nNew size of neighbourhood is set to ",
                k0, ".\n"
            ))
        }
    } else {
        if (verbose) {
            warning(
                "done.",
                "Heuristic did not change the neighbourhood.",
                sprintf(
                    "If results appear inconclusive, change k0 = %s.",
                    k0
                ),
                "",
                sep = "\n"
            )
        }
    }
    return(k0)
}

run.kbet.addTest <- function(initialize.kbet.res, adapt, alpha, n_repeat) {

    # initialize result list
    kbet.addTest.res <- initialize.kbet.addTest(
        rejection = initialize.kbet.res$rejection,
        k0 = initialize.kbet.res$k0,
        dof = initialize.kbet.res$dof,
        batch = initialize.kbet.res$batch,
        dim.dataset = initialize.kbet.res$dim.dataset,
        n_repeat = n_repeat
        )

    # run kBet test (with addTest)
    kbet.addTest.res <- kbet.addTest(
        rejection = kbet.addTest.res$rejection,
        n_repeat = n_repeat,
        dim.dataset = initialize.kbet.res$dim.dataset,
        testSize = initialize.kbet.res$testSize,
        k0 = initialize.kbet.res$k0,
        knn = initialize.kbet.res$knn,
        is.imbalanced = initialize.kbet.res$is.imbalanced,
        class.frequency = initialize.kbet.res$class.frequency,
        new.class.frequency = initialize.kbet.res$new.class.frequency,
        dof = initialize.kbet.res$dof,
        batch = initialize.kbet.res$batch,
        batch.shuff = initialize.kbet.res$batch.shuff,
        kBET.expected = initialize.kbet.res$kBET.expected,
        kBET.observed = initialize.kbet.res$kBET.observed,
        kBET.signif = initialize.kbet.res$kBET.signif,
        lrt.expected = kbet.addTest.res$lrt.expected,
        lrt.observed = kbet.addTest.res$lrt.observed,
        lrt.signif = kbet.addTest.res$lrt.signif,
        exact.expected = kbet.addTest.res$exact.expected,
        exact.observed = kbet.addTest.res$exact.observed,
        exact.signif = kbet.addTest.res$exact.signif,
        adapt = adapt,
        alpha = alpha
    )

    return(list(
        rejection = kbet.addTest.res$rejection,
        kBET.expected = kbet.addTest.res$kBET.expected,
        kBET.observed = kbet.addTest.res$kBET.observed,
        kBET.signif = kbet.addTest.res$kBET.signif,
        lrt.expected = kbet.addTest.res$lrt.expected,
        lrt.observed = kbet.addTest.res$lrt.observed,
        lrt.signif = kbet.addTest.res$lrt.signif
    ))
}

initialize.kbet.addTest <- function(
    rejection, k0, dof, batch, dim.dataset, n_repeat) {

    # initialize result list
    rejection$summary$lrt.expected <- numeric(4)
    rejection$summary$lrt.observed <- numeric(4)

    rejection$results$lrt.pvalue.test <- rep(0, dim.dataset[1])
    rejection$results$lrt.pvalue.null <- rep(0, dim.dataset[1])

    lrt.expected <- numeric(n_repeat)
    lrt.observed <- numeric(n_repeat)
    lrt.signif <- numeric(n_repeat)

    # decide to perform exact test or not
    if (choose(k0 + dof, dof) < 5e5 && k0 <= min(table(batch))) {
        exact.expected <- numeric(n_repeat)
        exact.observed <- numeric(n_repeat)
        exact.signif <- numeric(n_repeat)

        rejection$summary$exact.expected <- numeric(4)
        rejection$summary$exact.observed <- numeric(4)
        rejection$results$exact.pvalue.test <- rep(0, dim.dataset[1])
        rejection$results$exact.pvalue.null <- rep(0, dim.dataset[1])
    }

    return(list(
        rejection = rejection, lrt.expected = lrt.expected,
        lrt.observed = lrt.observed, lrt.signif = lrt.signif,
        exact.expected = exact.expected, exact.observed = exact.observed,
        exact.signif = exact.signif
    ))
}

kbet.addTest <- function(
    rejection, n_repeat, dim.dataset, testSize, k0, knn, is.imbalanced,
    class.frequency, new.class.frequency, batch, dof, alpha, batch.shuff,
    kBET.expected, kBET.observed, kBET.signif, lrt.expected, lrt.observed,
    lrt.signif, exact.expected, exact.observed, exact.signif, adapt) {
    for (i in seq_len(n_repeat)) {
        idx.runs <- sample.int(dim.dataset[1], size = testSize)
        env <- cbind(knn[idx.runs, seq_len(k0 - 1)], idx.runs)
        pka.res <- perform.kbet.addTest(
            env, adapt, is.imbalanced, class.frequency, new.class.frequency,
            batch, dof, alpha, batch.shuff)
        kBET.expected[i] <- sum(pka.res$p.val.test.null < alpha, na.rm = TRUE) /
            sum(!is.na(pka.res$p.val.test.null))
        kBET.observed[i] <- sum(pka.res$is.rejected, na.rm = TRUE) /
            sum(!is.na(pka.res$p.val.test))
        kBET.signif[i] <- 1 - ptnorm(
            kBET.observed[i], mu = kBET.expected[i],
            sd = sqrt(kBET.expected[i] * (1 - kBET.expected[i]) / testSize),
            alpha = alpha
        )
        rejection$results$tested[idx.runs] <- 1
        rejection$results$kBET.pvalue.test[idx.runs] <- pka.res$p.val.test
        rejection$results$kBET.pvalue.null[idx.runs] <- pka.res$p.val.test.null
        compute.LRT.res <- compute_LRT(
            rejection, lrt.expected, lrt.observed, lrt.signif, adapt,
            is.imbalanced, class.frequency, new.class.frequency, batch, dof,
            batch.shuff, env, alpha, i, testSize, idx.runs
        )
        rejection <- compute.LRT.res$rejection
        lrt.expected <- compute.LRT.res$lrt.expected
        lrt.observed <- compute.LRT.res$lrt.observed
        lrt.signif <- compute.LRT.res$lrt.signif
        if (exists(x = "exact.observed")) {
            exact.observed.res <- run.kbet.exact.observed(
                adapt, is.imbalanced, env, new.class.frequency, batch,
                batch.shuff, alpha, testSize, exact.expected, exact.observed,
                exact.signif, rejection, idx.runs, i, class.frequency
            )
            rejection <- exact.observed.res$rejection
            exact.expected <- exact.observed.res$exact.expected
            exact.observed <- exact.observed.res$exact.observed
            exact.signif <- exact.observed.res$exact.signif
        }
    }
    return(list(
        rejection = rejection, kBET.expected = kBET.expected,
        kBET.observed = kBET.observed, kBET.signif = kBET.signif,
        lrt.expected = lrt.expected, lrt.observed = lrt.observed,
        lrt.signif = lrt.signif))
}

perform.kbet.addTest <- function(
    env, adapt, is.imbalanced, class.frequency, new.class.frequency, batch, dof,
    alpha, batch.shuff) {

    # perform test
    if (adapt && is.imbalanced) {
        p.val.test <- apply(env, 1,
            FUN = chi_batch_test,
            new.class.frequency, batch, dof
            )
    } else {
        p.val.test <- apply(env, 1,
            FUN = chi_batch_test,
            class.frequency, batch, dof
        )
    }

    is.rejected <- p.val.test < alpha

    p.val.test.null <- apply(apply(
        batch.shuff, 2,
        function(x, freq, dof, envir) {
            apply(envir, 1, FUN = chi_batch_test, freq, x, dof)
        },
        class.frequency, dof, env
    ), 1, mean, na.rm = TRUE)

    return(list(p.val.test = p.val.test, p.val.test.null = p.val.test.null))
}

compute_LRT <- function(
    rejection, lrt.expected, lrt.observed, lrt.signif, adapt, is.imbalanced,
    class.frequency, new.class.frequency, batch, dof, batch.shuff, env, alpha,
    i, testSize, idx.runs) {

    cf <- if (adapt && is.imbalanced) {
        new.class.frequency
    } else {
        class.frequency
    }
    p.val.test.lrt <- apply(env, 1,
        FUN = lrt_approximation, cf, batch, dof
    )
    p.val.test.lrt.null <- apply(apply(
        batch.shuff, 2,
        function(x, freq, dof, envir) {
            apply(envir, 1, FUN = lrt_approximation, freq, x, dof)
        },
        class.frequency, dof, env
    ), 1, mean, na.rm = TRUE)

    lrt.expected[i] <-
        sum(p.val.test.lrt.null < alpha, na.rm = TRUE) /
            sum(!is.na(p.val.test.lrt.null))
    lrt.observed[i] <-
        sum(p.val.test.lrt < alpha, na.rm = TRUE) /
            sum(!is.na(p.val.test.lrt))

    lrt.signif[i] <-
        1 - ptnorm(lrt.observed[i],
            mu = lrt.expected[i], sd = sqrt(
                lrt.expected[i] * (1 - lrt.expected[i]) / testSize
            ),
            alpha = alpha
        )

    rejection$results$lrt.pvalue.test[idx.runs] <- p.val.test.lrt
    rejection$results$lrt.pvalue.null[idx.runs] <- p.val.test.lrt.null

    return(list(
        rejection = rejection, lrt.expected = lrt.expected,
        lrt.observed = lrt.observed, lrt.signif = lrt.signif))
}

run.kbet.exact.observed <- function(
    adapt, is.imbalanced, env, new.class.frequency, batch, batch.shuff, alpha,
    testSize, exact.expected, exact.observed, exact.signif, rejection,
    idx.runs, i, class.frequency) {
    if (adapt && is.imbalanced) {
        p.val.test.exact <- apply(
            env, 1, multiNom, new.class.frequency$freq, batch
        )
    } else {
        p.val.test.exact <- apply(
            env, 1, multiNom, class.frequency$freq, batch
        )
    }
    p.val.test.exact.null <- apply(apply(
        batch.shuff, 2, function(x, freq, envir) {
            apply(envir, 1, FUN = multiNom, freq, x)
        }, class.frequency$freq, env), 1, mean, na.rm = TRUE)
    exact.expected[i] <- sum(
        p.val.test.exact.null < alpha,
        na.rm = TRUE
    ) / testSize
    exact.observed[i] <- sum(
        p.val.test.exact < alpha,
        na.rm = TRUE
    ) / testSize
    # compute the significance level for the number of rejected data points
    exact.signif[i] <-
        1 - ptnorm(exact.observed[i],
            mu = exact.expected[i],
            sd = sqrt(
                exact.expected[i] * (1 - exact.expected[i]) / testSize
            ),
            alpha = alpha
        )
    # p-value distribution
    rejection$results$exact.pvalue.test[idx.runs] <- p.val.test.exact
    rejection$results$exact.pvalue.null[idx.runs] <- p.val.test.exact.null

    return(list(
        rejection = rejection, exact.expected = exact.expected,
        exact.observed = exact.observed, exact.signif = exact.signif))
}

run.kbet.only <- function(initialize.kbet.res, adapt, alpha, n_repeat) {
    rejection <- initialize.kbet.res$rejection
    dim.dataset <- initialize.kbet.res$dim.dataset
    testSize <- initialize.kbet.res$testSize
    k0 <- initialize.kbet.res$k0
    knn <- initialize.kbet.res$knn
    batch <- initialize.kbet.res$batch
    dof <- initialize.kbet.res$dof
    batch.shuff <- initialize.kbet.res$batch.shuff
    kBET.expected <- initialize.kbet.res$kBET.expected
    kBET.observed <- initialize.kbet.res$kBET.observed
    kBET.signif <- initialize.kbet.res$kBET.signif
    for (i in seq_len(n_repeat)) {
        # choose a random sample from dataset
        idx.runs <- sample.int(dim.dataset[1], size = testSize)
        env <- cbind(knn[idx.runs, seq_len(k0 - 1)], idx.runs)
        # perform test
        kbet.res <- perform.kbet.test(
            adapt, initialize.kbet.res$is.imbalanced,
            initialize.kbet.res$new.class.frequency,
            initialize.kbet.res$class.frequency,
            batch.shuff, env, batch, dof, alpha
        )
        # summarise test results
        kBET.expected[i] <- mean(apply(
            kbet.res$p.val.test.null, 2,
            function(x) sum(x < alpha, na.rm = TRUE) / sum(!is.na(x))
        ))
        kBET.observed[i] <- sum(kbet.res$is.rejected, na.rm = TRUE) /
            sum(!is.na(kbet.res$p.val.test))
        # compute significance
        kBET.signif[i] <- 1 - ptnorm(
            kBET.observed[i],
            mu = kBET.expected[i],
            sd = sqrt(kBET.expected[i] * (1 - kBET.expected[i]) / testSize),
            alpha = alpha
        )
        # assign results to result table
        rejection$results$tested[idx.runs] <- 1
        rejection$results$kBET.pvalue.test[idx.runs] <- kbet.res$p.val.test
        rejection$results$kBET.pvalue.null[idx.runs] <- rowMeans(
            kbet.res$p.val.test.null,
            na.rm = TRUE
        )
    }
    return(list(
        rejection = rejection, kBET.expected = kBET.expected,
        kBET.observed = kBET.observed, kBET.signif = kBET.signif
    ))
}

perform.kbet.test <- function(
    adapt, is.imbalanced, new.class.frequency, class.frequency, batch.shuff,
    env, batch, dof, alpha) {
    cf <- if (adapt && is.imbalanced) {
        new.class.frequency
    } else {
        class.frequency
    }
    p.val.test <- apply(env, 1, chi_batch_test, cf, batch, dof)
    is.rejected <- p.val.test < alpha
    p.val.test.null <- apply(
        batch.shuff, 2,
        function(x) {
            apply(
                env, 1, chi_batch_test, class.frequency,
                x, dof
            )
        }
    )
    return(list(
        p.val.test = p.val.test, p.val.test.null = p.val.test.null,
        is.rejected = is.rejected
    ))
}

mean_ci <- function(x, probs) {
    c(Mean = mean(x, na.rm = TRUE), quantile(x, probs, na.rm = TRUE))
}

summarize_kbet_results <- function(
    rejection, kBET.expected, kBET.observed, kBET.signif, lrt.expected = NULL,
    lrt.observed = NULL, lrt.signif = NULL, exact.observed = NULL,
    exact.expected = NULL, exact.signif = NULL, n_repeat = 1, addTest = FALSE) {
    if (n_repeat > 1) {
        # summarize chi2-results
        CI95 <- c(0.025, 0.5, 0.975)
        rejection$summary$kBET.expected <- mean_ci(kBET.expected, CI95)
        rownames(rejection$summary) <- c("mean", "2.5%", "50%", "97.5%")
        rejection$summary$kBET.observed <- mean_ci(kBET.observed, CI95)
        rejection$summary$kBET.signif <- mean_ci(kBET.signif, CI95)
        if (!addTest) {
            rejection$stats$kBET.expected <- kBET.expected
            rejection$stats$kBET.observed <- kBET.observed
            rejection$stats$kBET.signif <- kBET.signif
        } else {
            # summarize lrt-results
            rejection$summary$lrt.expected <- mean_ci(lrt.expected, CI95)
            rejection$summary$lrt.observed <- mean_ci(lrt.observed, CI95)
            rejection$summary$lrt.signif <- mean_ci(lrt.signif, CI95)
            # summarize exact test results
            if (!is.null(exact.observed)) {
                rejection$summary$exact.expected <- mean_ci(
                    exact.expected, CI95)
                rejection$summary$exact.observed <- mean_ci(
                    exact.observed, CI95)
                rejection$summary$exact.signif <- mean_ci(exact.signif, CI95)
            }
        }
        if (n_repeat < 10) {
            warning("Warning: The quantile computation for", n_repeat,
                "subset results is not meaningful.\n")
        }
    } else {
        rejection$summary$kBET.expected <- kBET.expected[1]
        rejection$summary$kBET.observed <- kBET.observed[1]
        rejection$summary$kBET.signif <- kBET.signif[1]
        if (addTest) {
            rejection$summary$lrt.expected <- lrt.expected
            rejection$summary$lrt.observed <- lrt.observed
            rejection$summary$lrt.signif <- lrt.signif
            if (!is.null(exact.observed)) {
                rejection$summary$exact.expected <- exact.expected
                rejection$summary$exact.observed <- exact.observed
                rejection$summary$exact.signif <- exact.signif
            }
        }
    }
    return(rejection)
}

plot_kbet_helper <- function(
    kBET.observed, kBET.expected, lrt.observed = NULL, lrt.expected = NULL,
    exact.observed = NULL, exact.expected = NULL, n_repeat) {
    if (!is.null(exact.observed)) {
        plot.data <- data.frame(
            class = rep(c(
                "kBET", "kBET (random)", "lrt", "lrt (random)",
                "exact", "exact (random)"
            ), each = n_repeat),
            data = c(
                kBET.observed, kBET.expected, lrt.observed,
                lrt.expected, exact.observed, exact.expected
            )
        )
    } else if (!is.null(lrt.observed)) {
        plot.data <- data.frame(
            class = rep(c("kBET", "kBET (random)", "lrt", "lrt (random)"),
                each = n_repeat
            ),
            data = c(
                kBET.observed, kBET.expected, lrt.observed,
                lrt.expected
            )
        )
    } else {
        plot.data <- data.frame(
            class = rep(c("observed(kBET)", "expected(random)"),
                        each = n_repeat),
            data = c(kBET.observed, kBET.expected)
        )
    }

    g <- ggplot(plot.data, aes(class, data)) +
        geom_boxplot() +
        theme_bw() +
        labs(x = "Test", y = "Rejection rate") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_y_continuous(limits = c(0, 1))

    return(g)
}