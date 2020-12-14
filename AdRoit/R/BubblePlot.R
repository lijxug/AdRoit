#' Make bubble plot to compare estimated proportions with truth
#'
#' This function makes overlaid bubble plots that help visualize the difference between estimated proportions and the true proportions.

#' @param prop.ests a matrix or data frame of estimated proporitons with row being cell types and columns being samples (spatial spots).
#' @param prop.real a matrix or data frame of the ground truth proportions.
#' @param labels a vector of two elements, labels the first two arguments. Default is `c("Estimated", "True")`.
#'
#' @return a ggplot object. Use `print` to render the plot.
#' @export


BubblePlot <- function(prop.ests,
                        prop.real,
                        labels = c("Estimatd", "Truth")){

    SampleID <- cell.type <- Percent <- method <- NULL

    ord = rowMeans(prop.real) %>% sort() %>% names() %>% rev()
    prop.ests = reshape2::melt(prop.ests) %>%
      mutate(method = labels[1])

    prop.real = reshape2::melt(prop.real) %>%
      mutate(method = labels[2])

    tmp = rbind(prop.ests,  prop.real)
    colnames(tmp)[1:3] = c("cell.type", "SampleID", "Percent")

    tmp$cell.type = factor(tmp$cell.type, levels = ord)

    pout = ggplot2::ggplot(tmp, aes(x=SampleID, y=cell.type, size=Percent, color=method)) +
        geom_point(alpha=0.4) +  scale_colour_manual(values = c("red", "blue")) +
        scale_size(range = c(.1, 24), name="Proportion") +  xlab("") + ylab("Cell type") + theme_classic() +
        theme(axis.text.x = element_text(angle = 35, vjust = 0.5, size = 15),
            axis.text.y = element_text(size = 17),
            axis.title.y = element_text(face = "bold", size = 17),
            legend.position = "top")

    return(pout)
}
