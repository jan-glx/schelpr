#'
#' @export
rc <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
