#' @importFrom plyr alply
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SingleCellExperiment colData
#' @importFrom tidyr separate
#' @importFrom S4Vectors DataFrame
res2sce <- function(res_array) {
    dat_matrix = res_array %>%
        alply(3, .dims = T) %>%
        Map(function(x, i) {
            colnames(x) = paste( colnames(x), i, sep = ":" )
            return(x)
        }, ., names(.)) %>%
        do.call(cbind, .)

    sce = SingleCellExperiment(assays=list(counts = dat_matrix))
    colData(sce) = data.frame( label=colnames(dat_matrix) ) %>%
        separate(label, into = c("sample", "cell_type"), sep = ":", remove = F) %>%
        DataFrame()
    return(sce)
}
