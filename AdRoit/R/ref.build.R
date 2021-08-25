#' Build the single cell reference
#'
#' This function takes a vector of gene counts, fits the negative binomial distribution, and estimate the size and mean parameters.
#'
#' @param counts a matrix/data.frame/sparse matrix (dgCMatrix) of single cell counts with the rows being genes and columns being cells.
#' @param annotations a vector annotates the cell types of cells.
#' @param genes the set of genes selected and to be used for deconvolution.
#' @param samples a vector annotates which sample ID the cells come from. Default is `NULL`.
#' @param normalize specify how to normalize the single cell raw counts. It should be chosen from "None", "Total" and "Median". Default is "None".
#' @param scale.factor the total number of reads to be normalized to if `normalize = "Total"`. Default is 1e+5.
#' @param multi.sample.bulk specify whether the bulk to be deconvoluted has multiple samples. Default is TRUE.
#' @param multi.sample.single specify whether the single cell reference has multiple samples. Default is TRUE.
#' @param silent whether to print out messages. Default is FALSE.
#' @importFrom parallel detectCores
#' @importFrom doParallel registerParallel stopImplicitCluster
#'
#' @return a list of elements including `mus`, `sizes`, `cell.specificity.w` and `cross.sample.w`.
#' \itemize{
#'   \item{mus    }{estimated means for selected genes.}
#'    \item{sizes    }{estimated size parameters for selected genes.}
#'    \item{cell.specificity.w    }{cell-type specificity weights.}
#'    \item{cross.sample.w    }{cross-sample varaibility weights.}
#'    }
#' @export

ref.build <- function(counts,
                       annotations,
                       genes,
                       samples = NULL,
                       normalize = "None",
                       scale.factor = 1e+05,
                       multi.sample.bulk = TRUE,
                       multi.sample.single = TRUE,
                       silent = FALSE){

  no_cores <- detectCores() - 1
  registerDoParallel(no_cores)
  negbin.par <- list()
  cell.types = unique(sort(annotations))
  for (c in cell.types) {
    cells.id = colnames(counts)[which(annotations == c)]
    if (normalize == "Total") {
      temp = round(scale.factor * (sweep(counts[, cells.id, drop = F],
                                         2, colSums(counts[, cells.id, drop = F]), `/`)))
    } else if (normalize == "Median") {
      medcnt = median(colSums(counts[, cells.id]))
      temp = round(medcnt * (sweep(counts[, cells.id, drop = F],
                                   2, colSums(counts[, cells.id, drop = F]), `/`)))
    } else if (normalize == "None") {
      temp = counts[, cells.id, drop = F]
    } else {
      message("normalize has to be one of Total, Median and None")
    }

    if (silent == FALSE) {
      message(paste("Estimate means and dispersions for cell type",
                    c, sep = ": "))
    }
    dim(temp)
    negbin.par[[c]] <- foreach(i = genes, .combine = rbind) %dopar%
      negbin.est(temp[i, ])
    rownames(negbin.par[[c]]) = genes
    rm(temp)
  }
  mus <- sizes <- NULL
  for (i in 1:length(negbin.par)) {
    mus = cbind(mus, negbin.par[[i]][, 2])
    sizes = cbind(sizes, negbin.par[[i]][, 1])
  }
  vars = mus + mus^2/sizes
  colnames(mus) = colnames(sizes) = cell.types
  midx = apply(mus, 1, which.max)
  w0 <- NULL
  for (i in 1:length(midx)) {
    w0 = c(w0, mus[i, midx[i]]/vars[i, midx[i]])
  }
  w0[which(is.na(w0) | is.infinite(w0))] = 0
  if (multi.sample.bulk == TRUE) {
    ref.est = list(mus = mus, lambda = sizes, cell.specificity.w = w0)
  } else {
    if (multi.sample.single == TRUE) {
      if (is.null(samples)) {
        stop("Please provide sample ID for each cell.")
      }
      sid = samples
    } else {
      nc = ncol(counts)
      subn = floor(nc/10)
      idx = NULL
      vc = seq_len(nc)
      sid = rep("N", 10 * subn)
      for (i in 1:10) {
        if (i == 10) {
          sid[vc] = paste("sample", i, sep = "_")
        } else {
          idx = sample(vc, subn, replace = F)
          sid[idx] = paste("sample", i, sep = "_")
          vc = vc[-idx]
        }
      }
    }
    single.pool = synthesize.bulk(counts, annotations, SampleID = sid)
    smv = foreach(i = genes, .combine = rbind) %dopar% negbin.est(single.pool[[1]][i,
    ])
    rownames(smv) = genes
    subv = smv[, 2] * (1 + smv[, 2]/smv[, 1])
    ref.est = list(mus = mus, lambda = sizes, cell.specificity.w = w0,
                   cross.sample.w = subv)
  }
  stopImplicitCluster()
  return(ref.est)
}
