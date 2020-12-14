#' Normalize single cell count matrix
#'
#' The function normalize the total counts of cells to the medians or means of counts in each cell cluster,
#' or to the same total number of reads given by `scale.factor`.
#'
#' @param counts the raw single cell count matrix, sparse matrix (dgCMatrix).
#' @param annotations a vector of length equal to the number of cells, each element annotates the cell type of a cell in the count matrix.
#' @param scale.factor total number of reads per cell (column) to be scaled to. Default is 1e+5.
#' @param normalize it has to one of "median", "mean" and "total".
#'
#' @return a normalized count matrix.
#' @export

normalize.cell <- function(counts, annotations, scale.factor = 1e5, normalize = "Total"){

  if (class(counts) == 'dgCMatrix'){
      counts = as.matrix(counts)
  }

  if (normalize %in% c("Median", "median", "Mean", "mean")) {
    message(paste0("scale.factor will be ignored. do ", normalize," normalization."))
    tcnt = colSums(counts)
    cell.cluster = unique(sort(annotations))
    normed.counts <- NULL
    for (i in cell.cluster){
      idx = which(annotations == i)
      if (normalize %in% c("Median", "median")){
        m.factor = floor(median(tcnt[idx]))
      } else if (normalize %in% c("Mean", "mean")){
        m.factor = floor(mean(tcnt[idx]))
      }
      sub.cnt = round(m.factor * sweep(counts[, idx], 2, colSums(counts[, idx]), `/`))
      normed.counts = cbind(normed.counts, as.matrix(sub.cnt))
    }
  } else if (normalize %in% c("Total", "total")){
    normed.counts = scale.factor*sweep(counts, 2, colSums(counts), `/`)
  } else {
    stop("Now only support Median, Mean, and Total normalization.")
  }
  return(normed.counts)
}
