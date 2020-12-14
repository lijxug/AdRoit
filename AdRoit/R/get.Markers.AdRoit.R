#' Function to get top marker genes for each cell cluster using AdRoit simplified way
#'
#' This function first normalizes the single cell count matrix (gene by cell) by equalizing the total counts then take `log1p`,
#' and perform Wilcoxon test to identify top differential expressed genes in each cell type.
#'
#' @param SC.counts the raw single cell count matrix, sparse matrix (dgCMatrix). The rows are genes and the columns are cells.
#' @param annotation a vector of length equal to number of cells, each element annotate the cell type of a cell in the count matrix.
#' @return a matrix containing a ranked list of putative markers, and associated statistics (p-values, fold change, etc.)
#' @export

get.Markers.AdRoit <- function(SC.counts, annotation){

    j <- FC <- P.value <- afc <- Gene <- Cluster <- NULL

    if (class(SC.counts) == "dgCMatrix"){
        SC.counts = as.matrix(SC.counts)
    }

    do.wilcoxon.test <- function(SC.mat, cell.id1, cell.id2, gene){

        a = as.numeric(SC.mat[gene, cell.id1])
        b = as.numeric(SC.mat[gene, cell.id2])

        FC = (exp(mean(a))+1e-9)/(exp(mean(b))+1e-9)
        w.test = wilcox.test(a, b)
        pvalue = w.test$p.value
        return(c(FC, pvalue))
    }

    genes = rownames(SC.counts)
    n.SCcounts = log1p(100000*sweep(SC.counts, 2, colSums(SC.counts), `/`))

    markers <- NULL
    for (i in unique(annotation)){

        message(paste("Do testing for ", i, " ......"))

        no_cores <- detectCores() - 1
        registerDoParallel(no_cores)

        idx1 = which(annotation == i)
        idx2 = which(annotation != i)

        cell.idx1 = colnames(SC.counts)[idx1]
        cell.idx2 = colnames(SC.counts)[idx2]

        tmp = foreach(j = genes, .combine = rbind) %dopar%
            do.wilcoxon.test(n.SCcounts, cell.idx1, cell.idx2, j)%>%
            as.data.frame()

        colnames(tmp) = c("FC", "P.value")
        tmp$Cluster = i
        tmp$Gene = genes

        tmp = tmp %>%
            arrange(desc(FC)) %>%
            filter(!is.na(P.value), FC > 0) %>%
            mutate(afc = ifelse(FC < 1, 1/FC, FC)) %>%
            top_n(2000, afc) %>%
            dplyr::select(Gene, Cluster, FC, P.value)

        markers = rbind(markers, tmp)
        rm(tmp)

        # stop node #
        stopImplicitCluster()
    }
    return(markers)
}
