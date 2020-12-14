#' plot ROC curves to compare two methods
#'
#' This function plots ROC curves to compare two methods given the estimated proportions from each one, and the true proportions.

#' @param method1.ests a matrix or data frame of estimated proporitons with row being cell types and columns being samples (or spatial spots).
#' @param method2.ests a matrix or data frame of estimated proporitons with row being cell types and columns being samples (or spatial spots).
#' @param prop.real a matrix or data frame of the ground truth proportions.
#' @param labels a vector of two elements labeling the two methods. Default is `c("Method1", "Metthod2")`.
#'
#' @return a ggplot object. Use `print` to render the plot.
#' @export


ROC_Plot <- function(method1.ests,
                     method2.ests,
                     prop.real,
                     labels = c("Method1", "Method2")){

    TPR <- FPR <- Method <- NULL

    tmp1 = cbind(as.numeric(method1.ests), as.numeric(prop.real))
    tmp2 = cbind(as.numeric(method2.ests), as.numeric(prop.real))

    get.Frate <- function(dat, cf = 0.0001){

        fpn = sum(dat[,1] > cf & dat[,2] == 0)
        fnn = sum(dat[,1] < cf & dat[,2] > 0)

        tpn = sum(dat[,1] > cf & dat[,2] > 0)
        tnn = sum(dat[,1] < cf & dat[,2] == 0)

        fpr = fpn/(fpn + tnn)
        tpr = tpn/(tpn + fnn)

        return(c(fpr, tpr))
    }

    met1.c <- met2.c <- NULL
    for (n in seq(0, 1, by = 0.0001)){
        met1.c <- rbind(met1.c, get.Frate(tmp1, n))
        met2.c <- rbind(met2.c, get.Frate(tmp2, n))
    }

    met1.c = as.data.frame(met1.c) %>% mutate(method = labels[1])
    met2.c = as.data.frame(met2.c) %>% mutate(method = labels[2])
    colnames(met1.c) = colnames(met2.c) = c("FPR", "TPR", "Method")

    tmp = rbind(met1.c, met2.c) %>% arrange(TPR)

    pout <- ggplot2::ggplot(tmp, aes(FPR, TPR, color = Method)) +
                geom_line(size = 1, alpha = 0.7)+
                scale_color_manual(values = c("blue", "orange")) +
                labs(x = "False Positive Rate (1-Specificity)",
                    y = "True Positive Rate (Sensitivity)") +
                theme_bw() +
                theme(legend.title = element_blank(),
                      legend.position = "top")
    return(pout)
}
