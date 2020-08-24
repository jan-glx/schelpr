
#' @export
print_encoded <- function(obj) {
  b <-  capture.output(a <- knitr::knit_print(obj))
  if (identical(b, character(0))) {
    if (!identical(a, obj)) {
      cat(paste0(sprintf("%s\n", a), collapse = ""))
    }
  } else {
    cat(paste0(sprintf("%s\n", b), collapse = ""))
  }
  invisible(obj)
}
