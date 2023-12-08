
#' @export
#' @import ggplot2
StatSplitSina <- ggproto(
  "StatSplitSina",
  ggforce::StatSina,
  required_aes = c("x", "y", "split"),
  compute_panel =  function(self, data, ...) {
    new_data <- self$super()$compute_panel(data, ...)
    side_factor <- (-1) ^ (new_data$split)
    new_data$x_diff <- abs(new_data$x_diff) * side_factor
    new_data$xmin <- new_data$x + pmin(0, new_data$width * side_factor)
    new_data$xmax <- new_data$x + pmax(0, new_data$width * side_factor)
    new_data
  }
)

#' @export
#' @import ggplot2
geom_split_sina <- function (mapping = NULL, data = NULL, stat = "splitSina", position = "dodge",
                             ..., na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomPoint,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, ...))
}

#' @export
#' @import ggplot2
stat_split_sina <- function (mapping = NULL, data = NULL, geom = "Point", position = "dodge",
                             scale = "area", method = "density", bw = "nrd0", kernel = "gaussian",
                             maxwidth = NULL, adjust = 1, bin_limit = 1, binwidth = NULL,
                             bins = NULL, seed = NA, ..., na.rm = FALSE, show.legend = NA,
                             inherit.aes = TRUE)
{
  method <- match.arg(method, c("density", "counts"))
  layer(data = data, mapping = mapping, stat = StatSplitSina, geom = geom,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(scale = scale, method = method, bw = bw,
                      kernel = kernel, maxwidth = maxwidth, adjust = adjust,
                      bin_limit = bin_limit, binwidth = binwidth, bins = bins,
                      seed = seed, na.rm = na.rm, ...))
}
