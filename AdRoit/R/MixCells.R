#' Simualte spatial spots using single cell reference data

#' This function samples cells from single cell reference data and mix them to simulate spatial spots.
#'
#' @param counts single cell count matrix (gene by cell) where cells are sampled from.
#' @param annotation cell type annotation for the cells in the count matrix.
#' @param celltypes a vector of cell type names to be included for sampling. The names must exsit in the annotation.
#' @param percent a vector of the same length as `celltypes`, specifying the percentages of corresponding cell types to be simulated.
#' @param batchsize specify how many cells per spot. Default is 10.
#' @param nspots specify how many spots to be simulated. Default is 100.
#' @param normalize whether to normalize the count matrix before mixing. Default is `FALSE`. If `TRUE`, the cells will be normalized to have total counts of 1e+5.
#' @param replace whether to sample cell with replacement. Defualt is `FALSE`.
#'
#' @return a matrix that simulated spatial spots.
#' @export

MixCells<-function(counts,
                   annotation,
                   celltypes,
                   percent,
                   batchsize=10,
                   normalize=FALSE,
                   nspots=100,
                   replace=F){

    set.seed(42)
    mixed=c()
    if(normalize){
        bases=colSums(counts)
        counts=sweep(counts, 2, bases, FUN="/")*1e+5
    }
    for(r in 1:nspots){
        tmp=c()
        for(i in 1:length(celltypes)){
            cell=celltypes[i]
            p=percent[i]/sum(percent)

            selall=names(annotation)[which(annotation==cell)]
            sampled=sample(x=selall,size = batchsize,replace=replace)
            tmp=cbind(tmp,rowSums(counts[,sampled])*p)
        }
        mixed=cbind(mixed,rowSums(tmp))
    }
    return(mixed)
}
