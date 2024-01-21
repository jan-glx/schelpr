
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
