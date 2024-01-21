
H <- function(p) -sum(p * log2(p))

H_max <- function(p, k) k*H(p/k)


find_best_partitioning_recursive <- function(p, k) {
  n <- length(p)
  best_H <- 0
  best_partitioning <- integer(0)

  if(k>0) {
    for(i in seq_len(if(k==1) pmin(n, k) else n)-1) {
      this_best_left <- find_best_partitioning_recursive(p[seq_len(i)], k-1)
      this_H_right <- H(sum(p[seq_len(n-i)+i]))
      this_H <- this_best_left$H + this_H_right
      this_partitioning <- c(this_best_left$partitioning, n)
      if(best_H <= this_H) {
        best_partitioning <- this_partitioning
        best_H <- this_H
      }
    }
  }


  list(partitioning=best_partitioning, H=best_H)
}


find_best_partitioning_non_recursive <- function(p, k) {
  n <- length(p)
  k <- pmin(n, k)
  best_H <- 0
  best_partitioning <- integer(0)

  best_partitionings <- vector(mode = "list", length = n) # element i contains a list of the best partitions of of the first i elements of p for different values of k
  empty_partitioning <- list(partitioning=integer(0), H=0)
  if(k<1) return(empty_partitioning)
  if(n<1) return(empty_partitioning)

  for (j in seq_len(k)) {
    best_partitionings_jm1 <- if(j-1==0) replicate(n, empty_partitioning, simplify=FALSE) else best_partitionings[[j-1]]
    #H_jm1 <- sapply(best_partitionings_jm1, `[[`, "H")
    #end_point_jm1 <- sapply(best_partitionings_jm1, `[[`, "H")
    best_i_jm1 <- 0
    for(i in if(j==k) n else seq_len(n)) {
      best_i_jm1 <- 0
      best_partitioning_i_jm1 <- integer(0)
      best_H <- 0
      for(i_jm1 in if(j==1) 0 else seq_len(i)-1) {
        current_partionin_jm1 <- if(i_jm1 == 0) empty_partitioning else best_partitionings_jm1[[i_jm1]]
        this_H <- current_partionin_jm1$H + H(sum(p[seq_len(i-i_jm1)+i_jm1]))
        if(this_H >= best_H) {
          best_H <- this_H
          best_i_jm1 <- i_jm1
          best_partitioning_i_jm1 <- current_partionin_jm1$partitioning
        }
      }
      best_partitionings[[j]][[i]] <- list(partitioning = c(best_partitioning_i_jm1, i) , H = best_H)
    }
  }
  best_partitionings[[k]][[n]]
}

find_best_partitioning <- function(weights, k) find_best_partitioning_non_recursive(weights/sum(weights), k)$partitioning



#' Bin a vector into a predefined number of as-equally-sized-as-possible intervals
#'
#' This function bins a vector of comparable values into a specified number of as-equally-sized intervals.
#' It uses the empirical cumulative distribution function (ECDF) and an optimal (in terms of entropy) partioning approach
#' to determine the bin edges, ensuring a more informative binning process compared to simple range division.
#' It can be considered a better version of base:cut
#'
#' @param x Numeric vector to be binned.
#' @param nbins Integer indicating the number of bins to create. Defaults to 5.
#' Must be a positive integer with a length of 1.
#'
#' @return An ordered factor with levels corresponding to the bins.
#'
#' @details The function first checks for input validity, ensuring that `x` is non-empty
#' and `nbins` is a valid integer greater than or equal to 1. It then computes the ECDF
#' of the input vector and applies an optimal binning algorithm (using find_best_partioning)
#' to find the best bin edges based on the desired number of bins. If `nbins` is set to 1,
#' the function returns an empty matrix.
#'
#' @examples
#' set.seed(123)
#' data <- rpois(100, rexp(100))
#' binned_data <- bin_it(data, nbins = 5)
#' binned_data
#' table(binned_data)
#'
#' @export
#' @importFrom stats ecdf
bin_it <- function(x, nbins = 5L) {
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  #if(length(x)==0) return(factor(x[0], ordered=TRUE))
  if(!is.wholenumber(nbins) || (length(nbins)>1)) stop("Illegal argument error: Number of bins (nbins) must be an length-1 integer ≥ 1, but is `", nbins, "`!")

  ECDF <- environment(ecdf(x))
  vx <- ECDF$x
  vy <- ECDF$y
  weights <- diff(c(0, vy))
  res <- find_best_partitioning(weights/sum(weights), nbins)
  nbins <- length(res)
  breaks <- vx[c(1, res)]
  labels <- paste0(rep(c("[", "("), c(1, pmax(0, nbins-1))), breaks[-(nbins+1)],",", breaks[-1],"]", recycle0 = TRUE)

  factor(findInterval(x,  breaks, left.open = TRUE, rightmost.closed = TRUE), levels= seq_len(nbins), labels = labels, ordered = TRUE)
}
