

#' @title StatBin2dInteger
#' @description A 2d binning statistic that calculates bins by rounding to integer values
#' @importFrom ggplot2 ggproto StatBin2d
#' @export
StatBin2dInteger <- ggproto("StatBin2dInteger", StatBin2d,
                            compute_group = function(data, scales, nbins = 30, drop = TRUE) {
                              if(length(nbins)==1) nbins <- rep(nbins, 2)

                              # Create helper functions
                              with_untrans <- function(trans, f, x) trans$trans(f(trans$inverse(x)))
                              round_to_halves <- function(x) round(x-0.5)+0.5
                              floor1p <- function(x) floor(x+1)
                              ceiling_minus_0.5 <- function(x) ceiling(x)-0.5

                              # Compute bins separately for x and y
                              compute_bins <- function(data, scale, nbins) {
                                binwidth <- diff(range(data)) / nbins
                                wu <- function(f, x) with_untrans(scale$trans, f, x)
                                bins <- unique(seq(wu(floor, min(data)), wu(floor1p, max(data)), by = binwidth))
                                bins_lower <- wu(ceiling_minus_0.5, bins)
                                bins <- bins[!duplicated(bins_lower)]
                                bins_lower <- bins_lower[!duplicated(bins_lower)]
                                binned_data <- factor(findInterval(data, bins_lower), levels = seq_len(length(bins_lower)-1))
                                list(bins, bins_lower, binned_data)
                              }

                              list_x <- compute_bins(data$x, scales$x, nbins[1])
                              list_y <- compute_bins(data$y, scales$y, nbins[2])

                              bins_x <- list_x[[1]]
                              bins_x_lower <- list_x[[2]]

                              bins_y <- list_y[[1]]
                              bins_y_lower <- list_y[[2]]

                              # Count the number of points in each bin
                              bin_data <- as.data.frame(table(list_x[[3]], list_y[[3]]))
                              names(bin_data) <- c("x", "y", "count")

                              # Convert x and y back to numeric
                              bin_data$x <- as.numeric(as.character(bin_data$x))
                              bin_data$y <- as.numeric(as.character(bin_data$y))


                              get_out_domain <- function(trans) trans$transform(trans$trans$domain)

                              # Compute the coordinates of the rectangles
                              bin_data <- within(bin_data, {
                                xmin <- bins_x_lower[x]
                                xmax <-bins_x_lower[x + 1]
                                ymin <- bins_y_lower[y]
                                ymax <- bins_y_lower[y + 1]
                                width <- xmax - xmin
                                height <- ymax - ymin
                                total_width <- max(xmax) - min(xmin)
                                total_height <- max(ymax) - min(ymin)
                                density <- count / (width/total_width) / (height/total_height)
                              })

                              # Compute the centers of the bins
                              bin_data$x <- bins_x[bin_data$x]
                              bin_data$y <- bins_y[bin_data$y]

                              bin_data <- rbind(data.frame(
                                x= NA, y= NA, count = 0,
                                xmin = -Inf, #get_out_domain(scales$x)[1],
                                xmax = Inf, #get_out_domain(scales$x)[2],
                                ymin = -Inf, #get_out_domain(scales$y)[1],
                                ymax = Inf, #get_out_domain(scales$y)[2],
                                width = Inf,
                                height= Inf,
                                total_width = Inf,
                                total_height = Inf,
                                density = 0
                              ), bin_data)

                              bin_data
                            },
                            required_aes = c("x", "y"),
                            default_aes = aes(weight = 1, fill = after_stat(density))
)

#' @title stat_bin_2d_integer
#' @description A function to use the StatBin2dInteger statistic
#' @param mapping The aesthetic mapping, usually constructed with aes or aes_.
#' @param data The data to be displayed in the layer.
#' @param geom The geometric object to use display the data
#' @param position Position adjustment, either as a string, or the result of a call to a position adjustment function.
#' @param na.rm If FALSE, the default, missing values are removed with a warning. If TRUE, missing values are silently removed.
#' @param show.legend logical. Should this layer be included in the legends?
#' @param inherit.aes If FALSE, overrides the default aesthetics, rather than combining with them.
#' @param bins number of bins
#' @param ... other arguments passed on to layer.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, disp)) +
#'   stat_bin_2d_integer(aes(fill = after_stat(count)), nbins = 30) +
#'   scale_y_log10()
#' @return ggplot layer
#' @importFrom ggplot2 layer
#' @export
stat_bin_2d_integer <- function(mapping = NULL, data = NULL, geom = "rect",
                                position = "identity", na.rm = FALSE, show.legend = NA,
                                inherit.aes = TRUE, nbins = 30, ...) {
  ggplot2::layer(
    stat = StatBin2dInteger, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(nbins = nbins, na.rm = na.rm, ...)
  )
}
