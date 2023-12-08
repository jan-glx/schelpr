
#' @export
scale_color_d3 <- ggsci::scale_color_d3
#' @export
scale_fill_d3 <- ggsci::scale_fill_d3

#' @export
geom_densityjitter <- purrr::partial(ggbeeswarm::geom_quasirandom, method = "pseudorandom")

#' @export
fix_plot_limits <- function(p) p + coord_cartesian(xlim=ggplot_build(p)$layout$panel_params[[1]]$x.range, ylim=ggplot_build(p)$layout$panel_params[[1]]$y.range)

#' @export
geom_cell <- function(..., size = 0.2, raster.dpi = 600) ggrastr::rasterize(ggplot2::geom_point(..., size = size), dpi = raster.dpi)

#' @export
print_by <- function(x, by, plot) {
  x <- as.data.table(x)
  split(x, by=by) %>% lapply(function(dt) print(plot(dt)))
  invisible(x)
}

#' @export
fix_width_lab <- function(lab, width=40) paste0(paste0(rep(" ", width), collapse = ""), "\n", lab)
