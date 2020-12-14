#' Estimate the proportions of cell types given by the reference
#'
#' This is the main function for estimating cell type proportions.
#'
#' @param bulk.sample a matrix or data.frame with the rows being genes and columns being samples (or spatial spots for spatial transcriptomics data).
#' @param single.ref the reference object built by the function `ref.build()`.
#' @param silent whether to print out messages. Default is FALSE.
#' @importFrom parallel detectCores
#' @importFrom doParallel registerParallel stopImplicitCluster
#' @importFrom stats optim
#'
#' @return a matrix with rows being cell types and columns being bulk samples, and entries are estimated proportions.
#' @export

AdRoit.est <- function(bulk.sample, single.ref, silent = FALSE){

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
        negbin.est(as.integer(bulk.sample[i, ]))
        rownames(tmp) = genes
      w = 1/(1 + tmp[, 2]/tmp[, 1])
    }
    w[which(is.infinite(w) | is.na(w))] = 0

    M = nnls::nnls(x, rowMeans(bulk.sample[genes, ]))
    ptheta = M$x/sum(M$x)
    msf = abs(rowMeans(bulk.sample[genes, ])/(x %*% ptheta))
    msf[which(is.infinite(msf) | is.na(msf))] = median(msf, na.rm = T)
    r = log(msf + 1)

    ns = ncol(x)
    theta0 = diff(sort(c(runif(ns - 1), 0, 1)))
    res <- NULL

    for (i in seq_len(nb)){

        if (silent == FALSE) {
            message(colnames(bulk.sample)[i])
        }

        y = bulk.sample[genes, i]
        fn = function(theta) {
            (w0 * w) %*% (y - r * (x %*% theta))^2 + sum(theta^2)
        }

        gn = function(theta) {
            -2 * t(sweep(x, 1, r, `*`)) %*% (w0 * w * (y - r*(x %*% theta))) + 2 * theta
        }

        op <- options(show.error.messages = FALSE)
        on.exit(op)
        esti = try(optim(theta0, fn, gr = gn,
                         lower = rep(0, ns),
                         upper = rep(Inf, ns), method = "L-BFGS-B"),
               silent = T)
        if (class(esti) == "try-error") {
            res = cbind(res, rep(NA, length(theta0)))
        } else {
            res = cbind(res, esti$par)
        }
    }

    res = sweep(res, 2, colSums(res), `/`)

    colnames(res) = colnames(bulk.sample)
    rownames(res) = colnames(x)
    options(show.error.messages = T)

    ## stop nodes ##
    stopImplicitCluster()

    return(res)
}
