
#' @export
#' @import ggplot2
gg_background <- function(background_value = 0, ..., mappping = aes(fill=background_value), data=data.frame(), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, show.legend = FALSE, inherit.aes = FALSE)
  geom_rect(mapping = mappping, data = data, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, show.legend = show.legend, inherit.aes = inherit.aes, ...)
