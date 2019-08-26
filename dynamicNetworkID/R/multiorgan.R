

#' Multi-organ gene expression data
#'
#' Data are from Anderson et al., 2017, Plos Computational Biology.
#'
#' @docType data
#'
#' @usage data(multiorgan)
#'
#' @keywords datasets
#'
#' @references Moore et al. (2013) Genetics 195:1077-1086
#' (\href{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005627})
#'
#'
#' @examples
#' data(multiorgan)
#' times <- unique(multiorgan[,1])
#' organs <- unique(multiorgan[,2])
#' genes = colnames(multiorgan)[2:ncol(multiorgan)]
"multiorgan"
