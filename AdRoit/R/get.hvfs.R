#' Get highly variable genes
#'
#' This function extract the highly variable genes from Seurat object. The function provides the option to balance the cluster sizes
#' before compute variations.
#'
#' @param ref.obj reference single cell matrix, data.frame or Seurat object.
#' @param annotations a vector annotates the cell types of cells.
#' @param n.top number of top highly variable genes to be selected.
#' @param balance whether to balance the sizes of the cell clusters. If TRUE, the function balances the number of cells in each cell
#' cluster to the median size of all by sampling cells from each cluster. Default is FALSE.
#'
#' @return a vector of top highly variable genes.
#' @export

get.hvfs <- function(ref.obj, annotations, n.top = 2000, balance = FALSE){

  if (balance == TRUE){
    cell.cluster = unique(sort(annotations))
    ncl = floor(median(table(annotations)))

    all.sid <- NULL
    for (i in cell.cluster){
      idx = which(annotations == i)
      if (ncl >= length(idx)){
        sidx = sample(idx, ncl, replace = T)
      } else {
        sidx = sample(idx, ncl, replace = F)
      }
      all.sid = c(all.sid, sidx)
    }

    if (class(ref.obj) == "Seurat"){
      umi.mat = as.data.frame(GetAssayData(ref.obj, slot = 'counts'))
      b.umi.mat = umi.mat[, all.sid]
    } else {
      b.umi.mat = ref.obj[, all.sid]
    }

    b.seurat.obj <- Seurat::CreateSeuratObject(b.umi.mat)
    b.seurat.obj <- Seurat::NormalizeData(object = b.seurat.obj, normalization.method = "LogNormalize",scale.factor = 1e+5)
    b.seurat.obj <- Seurat::FindVariableFeatures(b.seurat.obj, nfeature = n.top)
    hvfs <- VariableFeatures(b.seurat.obj)

  } else {
    if (class(ref.obj) == "data.frame" | class(ref.obj) == "matrix"){

      seurat.obj <- Seurat::CreateSeuratObject(ref.obj)
      seurat.obj <- Seurat::NormalizeData(object = seurat.obj, normalization.method = "LogNormalize",scale.factor = 1e+5)
      seurat.obj <- Seurat::FindVariableFeatures(seurat.obj, nfeature = n.top)
      hvfs <- VariableFeatures(seurat.obj)

    } else if (class(ref.obj) == "Seurat"){

      hvfs <- VariableFeatures(ref.obj)
      if (length(hvfs) == 0 | length(hvfs) < n.top){
        ref.obj <- Seurat::NormalizeData(object = ref.obj, normalization.method = "LogNormalize",scale.factor = 1e+5)
        ref.obj <- Seurat::FindVariableFeatures(ref.obj, nfeature = n.top)
        hvfs <- VariableFeatures(ref.obj)
      } else {
        hvfs = head(hvfs, n.top)
      }

      message("High variable genes are directly extracted Seurat object pre-computed!")

    } else {
      stop("ref.obj needs be counts matrix or data.frame or seurat object.")
    }
  }

  return(hvfs)
}
