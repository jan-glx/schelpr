
partition_path <- function(weights, k) {
  if(k==0) return(integer(0))
  n <- length(weights)
  if (k >= n) {
    return(seq_len(n))
  }

  cumsum_weights <- cumsum(weights)
  total_weight <- cumsum_weights[n]
  target_weight <- total_weight/k
  ends <- c(seq_len(k-1), n)
  block_weights <- diff(c(0, cumsum_weights[ends]))

  i <- k-1
  first_correct <- k
  moved <- FALSE
  while(i>0) {
    this_weigth <- weights[ends[i]+1]
    if(block_weights[i+1]-this_weigth >= block_weights[i] && (ends[i+1] > ends[i]+1)) {
      ends[i] <- ends[i]+1
      block_weights[i+1] <- block_weights[i+1] - this_weigth
      block_weights[i] <- block_weights[i] + this_weigth
      moved <- TRUE
    } else {
      if(moved && i+1 < k) {
        i <- i+1
      } else {
        i <- i-1
      }
      moved <- FALSE
    }
  }
  ends
}

#' Bin a vector into a predfined number of as-equally-sized-as-possible intervals
#'
#' This function bins a vector of compareble values into a specified number of as-equally-sized intervals.
#' It uses the empirical cumulative distribution function (ECDF) and an optimal matching approach
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
#' of the input vector and applies an optimal binning algorithm (using the Hungarian method)
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
  if(length(x)==0) return(factor(x[0], ordered=TRUE))
  if(length(x)>0 && (nbins<1 || !is.wholenumber(nbins)) || (length(nbins)>1)) stop("Illegal argument error: Number of bins (nbins) must be an length-1 integer ≥ 1, but is `", nbins, "`!")

  #ideal_quantiles <- seq(0, 1, length.out=1+nbins)[-1]
  ECDF <- environment(ecdf(x))
  vx <- ECDF$x
  vy <- ECDF$y
  weights <- diff(c(0, vy))
  res <- partition_path(weights, nbins)
  breaks <-vx[c(1, res)]
  labels <- paste0(rep(c("[", "("), c(1, nbins-1)), breaks[-(nbins+1)],",", breaks[-1],"]")
  factor(findInterval(x,  breaks[-1], left.open = TRUE, rightmost.closed = TRUE), levels= seq_len(nbins), labels = labels, ordered = TRUE)
}

#'
#'
#'
#' #' @title StatBin2dBinned
#' #' @description A 2d binning statistic that calculates bins by rounding to integer values
#' #' @importFrom ggplot2 ggproto StatBin2d
#' #' @export
#' StatBin2dBinned <- ggproto("StatBin2dBinned", StatBin2d,
#'                            setup_data = function(data, params) {
#'                              cat("params")
#'                              print( params)
#'                              nbins <- params$nbins
#'                              if(length(nbins)==1) nbins <- rep(nbins, 2)
#'                              cat("data in:")
#'                              print(data)
#'                              data$x <- bin_it(data$x, nbins[1])
#'                              data$y <- bin_it(data$y, nbins[2])
#'                              cat("data out:")
#'                              print(data)
#'                              data
#'                             },
#'                            extra_params = c("na.rm", "nbins")
#'   #                           compute_group = function(data, scales, nbins = 30, drop = TRUE) {
#'   #                             if(length(nbins)==1) nbins <- rep(nbins, 2)
#'   # print(data)
#'   #
#'   #                             binned_x <- bin_it(data$x_cont, nbins[1])
#'   #                             binned_y <- bin_it(data$y_cont, nbins[2])
#'   #
#'   #                             # Count the number of points in each bin
#'   #                             bin_data <- as.data.frame(table(binned_x, binned_y, data$group, data$PANEL, exclude = NULL))
#'   #                             names(bin_data) <- c("x", "y", "group", "PANEL", "count")
#'   #
#'   #                          #   bin_data$width <- 1
#'   #                         #    bin_data$height <- 1
#'   #
#'   #                             bin_data <- rbind(data.frame(
#'   #                               x= bin_data[1, "x"],
#'   #                               y= bin_data[1, "y"],
#'   #                               group= bin_data[1, "group"],
#'   #                               PANEL= bin_data[1, "PANEL"],
#'   #                               count = 0
#'   #                            #   width = Inf,
#'   #                            #   height= Inf
#'   #                             ), bin_data)
#'   #                             bin_data$fill <- bin_data$count
#'   #                             bin_data$x_cont <- bin_data$x
#'   #                             bin_data$y_cont <- bin_data$y
#'   #                             print(bin_data)
#'   #                             bin_data
#'   #                           },
#'                            # required_aes = c("x", "y"),
#'                            # default_aes = aes(x=after_stat(x), y=after_stat(y), fill = after_stat(fill))
#' )
#'
#' #' @title stat_bin_2d_binned
#' #' @description A function to use the StatBin2dBinned statistic
#' #' @param mapping The aesthetic mapping, usually constructed with aes or aes_.
#' #' @param data The data to be displayed in the layer.
#' #' @param geom The geometric object to use display the data
#' #' @param position Position adjustment, either as a string, or the result of a call to a position adjustment function.
#' #' @param na.rm If FALSE, the default, missing values are removed with a warning. If TRUE, missing values are silently removed.
#' #' @param show.legend logical. Should this layer be included in the legends?
#' #' @param inherit.aes If FALSE, overrides the default aesthetics, rather than combining with them.
#' #' @param bins number of bins
#' #' @param ... other arguments passed on to layer.
#' #' @examples
#' #' library(ggplot2)
#' #' ggplot(mtcars, aes(mpg, disp)) +
#' #'   stat_bin_2d_integer(aes(fill = after_stat(count)), nbins = 30) +
#' #'   scale_y_log10()
#' #' @return ggplot layer
#' #' @importFrom ggplot2 layer
#' #' @export
#' stat_bin_2d_binned <- function(mapping = NULL, data = NULL, geom = "tile",
#'                                 position = "identity", na.rm = FALSE, show.legend = NA,
#'                                 inherit.aes = TRUE, nbins = 30, ...) {
#'   ggplot2::layer(
#'     stat = StatBin2dBinned, data = data, mapping = mapping, geom = geom,
#'     position = position, show.legend = show.legend, inherit.aes = inherit.aes,
#'     params = list(nbins = nbins, na.rm = na.rm, ...)
#'   )
#' }
#'
#'
#' #' @title StatCountBinned2d
#' #' @description A 2d binning statistic that calculates bins by rounding to integer values
#' #' @importFrom ggplot2 ggproto StatBin2d
#' #' @export
#' StatCountBinned2d <- ggproto("StatCountBinned2d", StatSum,
#'                            setup_data = function(data, params) {
#'                              cat("params")
#'                              print( params)
#'                              nbins <- params$nbins
#'                              if(length(nbins)==1) nbins <- rep(nbins, 2)
#'                              cat("data in:")
#'                              print(data)
#'                              data$x <- bin_it(data$x, nbins[1])
#'                              data$y <- bin_it(data$y, nbins[2])
#'                              cat("data out:")
#'                              print(data)
#'                              browser()
#'                              data
#'                            },
#'                            extra_params = c("na.rm", "nbins"),
#'                            compute_layer = function(self, data, scales, ...)
#'                            #                           compute_group = function(data, scales, nbins = 30, drop = TRUE) {
#'                            #                             if(length(nbins)==1) nbins <- rep(nbins, 2)
#'                            # print(data)
#'                            #
#'                            #                             binned_x <- bin_it(data$x_cont, nbins[1])
#'                            #                             binned_y <- bin_it(data$y_cont, nbins[2])
#'                            #
#'                            #                             # Count the number of points in each bin
#'                            #                             bin_data <- as.data.frame(table(binned_x, binned_y, data$group, data$PANEL, exclude = NULL))
#'                            #                             names(bin_data) <- c("x", "y", "group", "PANEL", "count")
#'                            #
#'                            #                          #   bin_data$width <- 1
#'                            #                         #    bin_data$height <- 1
#'                            #
#'                            #                             bin_data <- rbind(data.frame(
#'                            #                               x= bin_data[1, "x"],
#'                            #                               y= bin_data[1, "y"],
#'                            #                               group= bin_data[1, "group"],
#'                            #                               PANEL= bin_data[1, "PANEL"],
#'                            #                               count = 0
#'                            #                            #   width = Inf,
#'                            #                            #   height= Inf
#'                            #                             ), bin_data)
#'                            #                             bin_data$fill <- bin_data$count
#'                            #                             bin_data$x_cont <- bin_data$x
#'                            #                             bin_data$y_cont <- bin_data$y
#'                            #                             print(bin_data)
#'                            #                             bin_data
#'                            #                           },
#'                            # required_aes = c("x", "y"),
#'                            # default_aes = aes(x=after_stat(x), y=after_stat(y), fill = after_stat(fill))
#' )
#'
#' #' @title stat_bin_2d_binned
#' #' @description A function to use the StatBin2dBinned statistic
#' #' @param mapping The aesthetic mapping, usually constructed with aes or aes_.
#' #' @param data The data to be displayed in the layer.
#' #' @param geom The geometric object to use display the data
#' #' @param position Position adjustment, either as a string, or the result of a call to a position adjustment function.
#' #' @param na.rm If FALSE, the default, missing values are removed with a warning. If TRUE, missing values are silently removed.
#' #' @param show.legend logical. Should this layer be included in the legends?
#' #' @param inherit.aes If FALSE, overrides the default aesthetics, rather than combining with them.
#' #' @param bins number of bins
#' #' @param ... other arguments passed on to layer.
#' #' @examples
#' #' library(ggplot2)
#' #' ggplot(mtcars, aes(mpg, disp)) +
#' #'   stat_bin_2d_integer(aes(fill = after_stat(count)), nbins = 30) +
#' #'   scale_y_log10()
#' #' @return ggplot layer
#' #' @importFrom ggplot2 layer
#' #' @export
#' stat_count_binned_2d <- function(mapping = NULL, data = NULL, geom = "tile",
#'                                position = "identity", na.rm = FALSE, show.legend = NA,
#'                                inherit.aes = TRUE, nbins = 30, ...) {
#'   ggplot2::layer(
#'     stat = StatCountBinned2d, data = data, mapping = mapping, geom = geom,
#'     position = position, show.legend = show.legend, inherit.aes = inherit.aes,
#'     params = list(nbins = nbins, na.rm = na.rm, ...)
#'   )
#' }
