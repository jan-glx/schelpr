
#' @export
print_encoded <- function(obj) {
  b <-  capture.output(a <- knitr::knit_print(obj))
  if (identical(b, character(0))) {
    if (!identical(a, obj)) {
      cat(paste0(a, collapse="\n"))
    }
  } else {
    cat(paste0(b, collapse="\n"))
  }
  invisible(obj)
}
