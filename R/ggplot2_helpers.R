
#' @export
scale_color_d3 <- ggsci::scale_color_d3
#' @export
scale_fill_d3 <- ggsci::scale_fill_d3


#' @export
geom_densityjitter <- purrr::partial(ggbeeswarm::geom_quasirandom, method = "pseudorandom")


#' @export
quantile_ci <- function(y, q=0.5, alpha=0.95, na.rm=FALSE) {
  if (na.rm) y <- y[!is.na(y)]
  n <- length(y)
  data.frame(
    y = quantile(y, probs = q, na.rm=na.rm, type = 8),
    ymin = quantile(y, probs = qbinom((1-alpha)/2, size=n, prob = q) / n, na.rm=na.rm, type = 1, g=0),
    ymax = quantile(y, probs = qbinom(1-(1-alpha)/2, size=n, prob = q) / n, na.rm=na.rm, type = 1, g=1)
  )
}


#' @export
fix_plot_limits <- function(p) p + coord_cartesian(xlim=ggplot_build(p)$layout$panel_params[[1]]$x.range, ylim=ggplot_build(p)$layout$panel_params[[1]]$y.range)



#' @export
#' @import ggplot2
#' @importFrom ggforce StatSina
StatSplitSina <- ggproto("StatSplitSina", StatSina, required_aes = c("x", "y", "split"), compute_panel =  function(self, data, ...) {
  new_data <- self$super()$compute_panel(data, ...)
  side_factor <- (-1) ^ (new_data$split)
  new_data$x_diff <- abs(new_data$x_diff) * side_factor
  new_data$xmin <- new_data$x + pmin(0, new_data$width * side_factor)
  new_data$xmax <- new_data$x + pmax(0, new_data$width * side_factor)
  new_data
})

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
