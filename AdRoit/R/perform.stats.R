#' Computate the statistics to evaluate deconvolution performance
#'
#' This function computes four evaluation statistics of the estimated proportions given the the ground truth.
#' The 4 statistics include mean absolute difference (mAD), rooted mean squared deviation (RMSD), Pearson and Spearman correlations.
#' It also summarizes the number of false postives and false negatives based on a given cutoff of detection.
#'
#' @param prop.ests a matrix or data frame of estimated proporitons with row being cell types and columns being samples (or spatial spots).
#' @param prop.real a matrix or data frame of the ground truth proportions.
#' @param false.cutoff a numeric value (between 0 and 1) used to determine whether a cell type is detected. Default value is 1e-4.
#'
#' @return a list of elements including mAD, RMSD, Pearson and Spearman correlations, number of false positives and false negatives.
#' The function summarizes the overall performance in two ways:
#' 1) the four statisitcs are computated per sample then the means are summarized (stored in `Stat_means`),
#' and 2) the samples were pooled first then the four statistics are computated over all values (stored in `Stat_all`).
#' @export
#'

perform.stats <- function(prop.ests, prop.real, false.cutoff = 1e-4){


  if (!(all(rownames(prop.ests) %in% rownames(prop.real)) &
        all(colnames(prop.ests) %in% colnames(prop.real)))) {
    stop("Some cell types (rows) or Samples (columns) do not exsit in the reference!")
  } else {
    prop.real = prop.real[rownames(prop.ests), colnames(prop.ests)]
  }

  nsub = ncol(prop.ests)

  if (nsub == 1){
    message("The first column is used for evaluation.")
  }

  mAD <- RMSD <- Pearson <- Spearman <- Eucd <- n.FN <- n.FP <- NULL
  for (i in seq_len(nsub)){
    mAD <- c(mAD, mean(abs(prop.ests[,i]-prop.real[,i])))
    RMSD <- c(RMSD, sqrt(mean((prop.ests[,i]-prop.real[,i])^2)))
    Pearson <- c(Pearson, cor(prop.ests[,i],prop.real[,i]))
    Spearman <- c(Spearman, cor(prop.ests[,i],prop.real[,i], method = "spearman"))
    n.FN <- c(n.FN, sum(prop.ests[,i] < false.cutoff & prop.real[,i] > 0))
    n.FP <- c(n.FP, sum(prop.real[,i] == 0 & prop.ests[,i] > false.cutoff))
  }

  stat.mean = c("mAD" = mean(mAD),
                "RMSD" = mean(RMSD),
                "Pearson" = mean(Pearson),
                "Spearman" = mean(Spearman),
                "n.FN" = mean(n.FN),
                "n.FP" = mean(n.FP))

  stat.all = c("mAD" = round(mean(abs(c(as.matrix(prop.real)) - c(as.matrix(prop.ests)))), digits = 5),
               "RMSD" = round(sqrt(mean((as.matrix(prop.real) - as.matrix(prop.ests))^2)), digits = 5),
               "Pearson" = round(cor(c(as.matrix(prop.real)), c(as.matrix(prop.ests))), digits = 4),
               "Spearman" = round(cor(c(as.matrix(prop.real)), c(as.matrix(prop.ests)), method = "spearman"), digits = 4),
               "n.FN" = sum(n.FN),
               "n.FP" = sum(n.FP))


  return(list("mAD" = mAD,
              "RMSD" = RMSD,
              "Pearson_correlation" = Pearson,
              "Spearman_correlation" = Spearman,
              "n.FN" = n.FN,
              "n.FP" = n.FP,
              "Stats_means" = stat.mean,
              "Stats_all" = stat.all))
}
