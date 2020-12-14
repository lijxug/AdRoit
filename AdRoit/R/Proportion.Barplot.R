#' Make bar plot to display percentages of cell types per sample
#'
#' This function makes bar plot to visualize the percentages of cell types per sample or spatial spots.

#' @param prop.ests a matrix or data frame of estimated proporitons with row being cell types and columns being samples (or spatial spots).
#'
#' @return a ggplot object. Use `print` to render the plot.
#' @export


Proportion.Barplot <- function(prop.ests){

    sample.id <- percent <- cell.type <- NULL

    tmp = reshape2::melt(prop.ests)
    colnames(tmp) = c("cell.type", "sample.id", "percent")
    pout <- ggplot2::ggplot(tmp, aes(sample.id, percent, fill = cell.type)) +
        geom_bar(stat = "identity") + theme_bw() + ylab("Cell proprotion") + xlab("") +
        theme(axis.text.x=element_text(angle = 45, size = 8, hjust = 1),
            text = element_text(size=12),
            legend.title = element_blank())
    return(pout)
}
