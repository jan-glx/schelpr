
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



#' shuffel, i. e. randomly reorder observations
#'
#' @export
shuffle <- function(object) UseMethod("shuffle", object)

#' @method shuffle data.frame
#' @export
shuffle.data.frame <- function(df) df[sample(nrow(df)),]

#' @method shuffle vector
#' @export
shuffle.vector <- function(x) x[sample(length(x))]

#' Create data.table from sparse matrix
#' @export
#' @examples
#' sparse2long(Matrix::spMatrix(1:3,2:4,3:5, ncol=4, nrow=3))
sparse2long <- function(mat, value_name = "value") {
  mat <- as(mat, "dgTMatrix")
  dt <- data.table(row = mat@i + 1, col = mat@j + 1, value = mat@x)
  if (!is.null(mat@Dimnames[[1]])) dt[, row := mat@Dimnames[[1]][row]]
  if (!is.null(names(mat@Dimnames)[1])) setnames(dt, "row", names(mat@Dimnames)[1])
  if (!is.null(mat@Dimnames[[2]])) dt[, col := mat@Dimnames[[2]][col]]
  if (!is.null(names(mat@Dimnames)[2])) setnames(dt, "col", names(mat@Dimnames)[2])
  setnames(dt, "value", value_name)
  dt[]
}


#' Create sparse matrix from long data.table
#' @export
#' @examples
#' long2sparse(1:3, 2:4, runif(3))
long2sparse <- function(rows, cols, values, dimname_rows = base::deparse1(substitute(rows)), dimname_cols = base::deparse1(substitute(rows))) {
  force(dimname_rows)
  force(dimname_cols)
  rows <- factor(rows)
  cols <- factor(cols)
  dimnames <- list(levels(rows), levels(cols))
  names(dimnames) <- c(dimname_rows, dimname_cols)
  Matrix::sparseMatrix(i = as.integer(rows), j = as.integer(cols), x = values, dimnames = dimnames)
}
