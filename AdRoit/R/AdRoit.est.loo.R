#' Estimate the proportions for a single sample which is left out in training data
#'
#' This function estimate the proportions of a single bulk sample while being left out when estimating gene cross-sample variations and cross-platform biases.
#'
#' @param n index of sample (column) to be estimated.
#' @param bulk.sample a matrix or data.frame with the rows being genes and columns being samples (or spatial spots for spatial transcriptomics data).
#' @param single.ref the reference object built by the function `ref.build()`.
#' @param silent whether to print out messages. Default is FALSE.
#' @importFrom parallel detectCores
#' @importFrom doParallel registerParallel stopImplicitCluster
#' @importFrom stats optim
#'
#'
#' @return the cell type proportions for the nth sample.
#' @export

AdRoit.est.loo <- function(n, bulk.sample, single.ref, silent = FALSE){

    i <- NULL

    # Calculate the number of cores
    no_cores <- detectCores() - 1

    # Initiate cluster
    registerDoParallel(no_cores)

    genes = rownames(single.ref[[1]])

    x = single.ref[[1]]
    z = single.ref[[2]]
    w0 = single.ref[[3]]

    nb = ncol(bulk.sample)
    if (nb < 3) {
        w = single.ref[[4]]
    } else {
        tmp = foreach(i = genes, .combine = rbind) %dopar%
            negbin.est(as.integer(bulk.sample[i, -n]))
        rownames(tmp) = genes
        w = 1/(1 + tmp[, 2]/tmp[, 1])
    }
    w[which(is.infinite(w) | is.na(w))] = 0

    M = nnls::nnls(x, rowMeans(bulk.sample[genes, -n]))
    ptheta = M$x/sum(M$x)
    msf = abs(rowMeans(bulk.sample[genes, -n])/(x %*% ptheta))
    msf[which(is.infinite(msf) | is.na(msf))] = median(msf, na.rm = T)

    ns = ncol(x)
    theta0 = diff(sort(c(runif(ns - 1), 0, 1)))

    if (silent == FALSE) {
        message(colnames(bulk.sample)[n])
    }

    y = bulk.sample[genes, n]
    fn = function(theta) {
        (w0 * w) %*% (y - log(msf + 1) * (x %*% theta))^2 + sum(theta^2)
    }
    gn = function(theta) {
        -2 * t(x) %*% ((w0 * w) * (y - log(msf + 1) * (x %*% theta))) + 2 * theta
    }

    op <- options(show.error.messages = FALSE)
    on.exit(op)
    esti = try(optim(theta0, fn, gr = gn, lower = rep(0, ns),
                     upper = rep(Inf, ns), method = "L-BFGS-B"),
             silent = T)
    if (class(esti) == "try-error") {
        res = rep(NA, length(theta0))
    } else {
        res = esti$par
    }

    res = as.matrix(data.frame(thetas = res/sum(res)))

    colnames(res) = colnames(bulk.sample)[n]
    rownames(res) = colnames(x)
    options(show.error.messages = T)

    ## stop nodes ##
    stopImplicitCluster()

    return(res)
}
