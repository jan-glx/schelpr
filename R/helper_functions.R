
#' @export
seq2vec <- function(seq) seq %>% str_replace_all("A", "4000") %>% str_replace_all("C", "0400") %>% str_replace_all("G", "0040") %>% str_replace_all("T", "0004") %>% str_replace_all("N", "1111") %>% tstrsplit("") %>% lapply(as.integer) %>% do.call(cbind, .) %>% `/`(4)


#' @author Aaron Lun
#' @export
setGeneric(".my_colVars", function(x, ...) standardGeneric(".my_colVars"))

#' @author Aaron Lun
#' @importFrom Matrix t rowSums
#' @export
setMethod(".my_colVars", "ANY", function(x, center=NULL) {
  if (is.null(center)) center <- Matrix::colMeans(x)
  y <- t(x) - center
  Matrix::rowSums(y^2)/(ncol(y)-1)
})

#' @author Aaron Lun
#' @author Jan Gleixner
#' @importFrom Matrix t colSums
#' @export
setMethod(".my_colVars", "dgCMatrix", function(x, center=NULL) {
  if (is.null(center)) center <- Matrix::colMeans(x)
  nzero <- diff(x@p)
  expanded <- rep(center, nzero)
  x@x <- (x@x - expanded)^2
  (Matrix::colSums(x) + (nrow(x)-nzero) * center^2)/(nrow(x)-1)
})


#' @export
setGeneric(".my_rowVars", function(x, ...) standardGeneric(".my_rowVars"))

#' @importFrom Matrix t rowSums
#' @export
setMethod(".my_rowVars", "ANY", function(x, center=NULL) {
  if (is.null(center)) center <- Matrix::rowMeans(x) else stopifnot(length(center) == nrow(x))
  y <- x - center
  Matrix::rowSums(y^2)/(ncol(y)-1)
})

#' @importFrom Matrix t colSums
#' @export
setMethod(".my_rowVars", "dgCMatrix", function(x, center=NULL) {
  if (is.null(center)) center <- Matrix::rowMeans(x) else stopifnot(length(center) == nrow(x))
  nzero <- base::tabulate(x@i+1L, nrow(x))
  x@x <- (x@x - center[x@i+1L])^2
  (Matrix::rowSums(x) + (ncol(x)-nzero) * center^2)/(ncol(x)-1)
})

#' @export
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# .my_aggregate.matrix <- function (x, groupings = NULL, form = NULL, fun = "sum", ...)
# {
#   if (!is(x, "Matrix"))
#     x <- Matrix(as.matrix(x), sparse = TRUE)
#   if (fun == "count")
#     x <- x != 0
#   groupings2 <- groupings
#   if (!is(groupings2, "data.frame"))
#     groupings2 <- as(groupings2, "data.frame")
#   groupings2 <- data.frame(lapply(groupings2, as.factor))
#   groupings2 <- data.frame(interaction(groupings2, sep = "_"))
#   colnames(groupings2) <- "A"
#   if (is.null(form))
#     form <- as.formula("~0+.")
#   form <- as.formula(form)
#   mapping <- dMcast(groupings2, form)
#   colnames(mapping) <- substring(colnames(mapping), 2)
#
#
#   if (fun == "var") {
#     center <- .my_aggregate.matrix(x, groupings = groupings, form = form, fun = "mean", ...)
#     t(mapping)
#   }
#
#   result <- t(mapping) %*% x
#   if (fun == "mean")
#     result@x <- result@x/(aggregate.Matrix(x, groupings2,
#                                            fun = "count"))@x
#   attr(result, "crosswalk") <- grr::extract(groupings, match(rownames(result),
#                                                              groupings2$A))
#   return(result)
# }
#
# if (is.null(center)) center <- Matrix::colMeans(x)
# nzero <- diff(x@p)
# expanded <- rep(center, nzero)
# x@x <- (x@x - expanded)^2
# (Matrix::colSums(x) + (nrow(x)-nzero) * center^2)/(nrow(x)-1)
