#' @title Plot proportion
#'
#' @param object ENIGMA object
#'
#' @importFrom reshape2 melt
#' @importFrom magrittr %>%
#' @import ggplot2
#' @export
#'
plot_proportion <- function(object) {
    object@result_cell_proportion %>%
        melt(varnames = c("sample", "cell_type"), value.name = "Freq.") %>%
        ggplot(aes(sample, Freq., fill=cell_type)) +
        geom_bar(stat="identity",position="fill") +
        labs(x = "Sample", y = "Proportion", fill = "Cell type") +
        coord_flip() +
        theme_bw() +
        theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12, color = 'black')) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

