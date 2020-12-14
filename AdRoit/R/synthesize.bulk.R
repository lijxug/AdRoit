#' Synthesize bulk by pooling single cell, and with user specified percentages
#'
#' The function creates simualted bulk samples by pooling the reads from each single cell type per sample.
#' It also can be used to create synthetic bulk samples with user-defined proportions.
#'
#' @param counts single cell UMI count matrix or sparse matrix (dgCMatrix), with rownames being genes and colnames being cell names.
#' @param annotations cell type annotation, a vector of length equal to the number of columns in counts.
#' @param SampleID sample identifications, a vector of length equal to the number of columns in counts.
#' @param prop either a string `original` or a vector of porportions with the length equal to the number of cell types and add up to 1.
#' Default is `original`. `original` pools all the cells in the provided single cell matrix. Alternatively, the cells will be sampled
#' to match the user-defined proportions, then pool.
#'
#' @return simulated bulk samples,  the ground truth cell counts and proportions per cell type per sample.
#' @export

synthesize.bulk <- function(counts, annotations, SampleID, prop = "original"){

  clusters = unique(sort(annotations))
  if (class(counts) == "dgCMatrix"){
      counts = as.data.frame(counts)
  }

  if (prop == "original"){
    s.SCcounts = counts
    s.annotations = annotations
    s.SampleID = SampleID

  } else if (length(prop) == length(clusters) & sum(prop) == 1) {
    N = ncol(counts)
    ns = round(prop*N)

    sampled.cells <- NULL
    for (i in 1:length(ns)){
      cluster.cells = colnames(counts)[which(annotations == clusters[i])]
      sampled.cells = c(sampled.cells, sample(cluster.cells, ns, replace = T))
    }

    s.SCcounts = counts[, sampled.cells]
    s.annotations = annotations[which(colnames(counts) %in% sampled.cells)]
    s.SampleID = SampleID[which(colnames(counts) %in% sampled.cells)]

  } else {
    message("prop needs to be either \"original\" or a vector of length
        	equal to number of cluster and sums to 1" )
  }

  tsingle = as.data.frame(transpose(s.SCcounts))
  tsingle$SampleID = s.SampleID
  tsingle$cluster = s.annotations

  simbulk = tsingle %>% group_by(SampleID) %>%
    summarise_if(is.numeric, sum) %>% transpose()
  colnames(simbulk) = simbulk[1,]
  simbulk = apply(simbulk[-1,],2, as.numeric)
  rownames(simbulk) = rownames(counts)


  # compute the real proportion
  count.table = table(s.SampleID, s.annotations) %>%
    as.data.frame.matrix() %>% t()
  prop.table = sweep(count.table, 2, colSums(count.table), `/`)
  rownames(simbulk) = rownames(counts)

  return(list(simbulk, prop.table, count.table))
}
