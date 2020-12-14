#' Estimate gene mean and dispersion by fitting negative binomial distribution
#'
#' This function takes a vector of gene counts, fits the negative binomial distribution, and estimate the size and mean parameters.
#'
#' @param X a vector of gene count to fit negative binomial distribution.
#' @return estimated size and mean.
#' @export

negbin.est <- function(X){

    defaultW <- getOption("warn")
    options(warn = -1)
    options(show.error.messages = FALSE)
    if (sum(X) == 0){
        par = c(size = 0, mu = 0)
    } else {
        fit = tryCatch(fitdistrplus::fitdist(X, "nbinom"), error = function(e){})
        if (is.null(fit)){
            par = c(size = 0, mu = 0)
        } else {
            par = fit$estimate
        }
    }
    options(warn = defaultW)
    options(show.error.messages = TRUE)
    return(par)
}
