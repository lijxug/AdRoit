#' Function to find marker genes per cell cluster using Seurat
#'
#' This function is a wraper to find markers for cell clusters using Seruat
#'
#' @param ref.obj a Seurat object after normalization and scaling.
#' @param annotations a vector of length equal to number of cells, each element annotates the cell type of a cell in the `ref.obj`.
#' @return a matrix containing a ranked list of putative markers, and associated statistics (e.g., p-values, fold change, etc.)
#' @export


get.Markers.Seurat <- function(ref.obj, annotations){

    if (class(ref.obj) == "Seurat") {
        Idents(ref.obj) = annotations
        ref.markers <- Seurat::FindAllMarkers(ref.obj)
    } else {
        stop("ref.obj needs be processed Seurat object")
    }

    return(ref.markers)
}
