#' Select the top marker genes from each cell cluster
#'
#' This function selects the top marker genes given marker table computated using function `get.Markers.Seurat` or user provided.
#' The function filters genes given p value cutoff, ranks them according to the absolute log fold change and select top genes.
#'
#' @param marker.table a data.frame or matirx contains gene names, effect size (e.g., fold change), cell cluster and p.value columns.
#' @param gene.col the column name of gene names in marker.table.
#' @param cluster.col the column name of cell clusters in marker.table.
#' @param fold.change.col the column name of fold change in marker.table.
#' @param p.value.col the column name of p values marker.table.
#' @param mito.genes a vector of mitochondria genes. These genes will be excluded from selected genes. Default is `NULL`.
#' @param pval.cutoff a numeric value to specify cutoff of p values.
#' @param n.top specify how many genes to be selected from each cell cluster. Default is 200.
#'
#' @return a vector of unique gene names to be used for deconvolution.
#' @export

get.TopMarkers <- function(marker.table,
                           gene.col,
                           cluster.col,
                           fold.change.col,
                           p.value.col,
                           mito.genes = NULL,
                           pval.cutoff = 0.01,
                           n.top = 200){
    FC <- p.value <- cluster <- alFC <- NULL

    ref.markers = data.frame("gene" = marker.table[, gene.col],
                             "cluster" = marker.table[, cluster.col],
                             "FC" = marker.table[, fold.change.col],
                             "p.value" = marker.table[, p.value.col])

    colnames(ref.markers) = c("gene", "cluster", "FC", "p.value")

    TopMarkers <- ref.markers %>% mutate(alFC = abs(log(FC))) %>%
        filter(p.value < pval.cutoff) %>% group_by(cluster) %>%
        top_n(n.top, alFC)

    u.feat = TopMarkers$gene %>% unique() %>% as.character()

    if (is.null(mito.genes)){
        out = u.feat
    } else {
        out = u.feat[which(!(u.feat %in% mito.genes))]
    }

    message("Total unique markers selected: ", length(out))
    if (length(out) < 1000) {
        warning("Less than 1K features are selected. Consider increase n.top to get better estimation.")
    }

    return(out)
}
