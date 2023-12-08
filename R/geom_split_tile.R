

#' Provides a diagonally split version of GeomTile
#'
#'
#' based on ggplot2::GeomTile by the ggplot2 authors https://github.com/tidyverse/ggplot2/graphs/contributors
#' @import ggplot2 rlang
#' @export
draw_key_split_tile <- function(data, params, size) {

  data$width <- data$width %||% params$width %||% 1
  data$height <- data$height %||% params$height %||% 1
  data$width[is.na(data$width)] <- 1
  data$height[is.na(data$height)] <- 1
  if (isTRUE(data$split)) {
    x <- c(0, 1, 0)
    y <- c(0, 1, 1)
  } else {
    x <- c(0, 1, 1)
    y <- c(0, 1, 0)
  }
  x <- 0.5 + (x-0.5) * data$width
  y <- 0.5 + (y-0.5) * data$height

  grid::polygonGrob(
    x = x,
    y = y,
    default.units = "npc",
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha)
    )
  )
}


#' @export
geom_split_tile <- function(mapping = NULL, data = NULL,
                            stat = "identity", position = "identity",
                            ...,
                            linejoin = "mitre",
                            na.rm = FALSE,
                            show.legend = NA,
                            inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSplitTile,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = rlang:::list2(
      linejoin = linejoin,
      na.rm = na.rm,
      ...
    )
  )
}

#' @export
GeomSplitTile <- ggproto(
  "GeomSplitTile", GeomPolygon,
  extra_params = c("na.rm"),

  setup_data = function(data, params) {
   data$width <- data$width %||% params$width %||% resolution(data$x, FALSE)
   data$height <- data$height %||% params$height %||% resolution(data$y, FALSE)
   data$split <- data$split %||% params$split %||% FALSE %>% as.factor

   K = 3
   n <- nrow(data)
   new_data <- data.frame(
     x = rep(data$x, each=K) + rep(3-as.integer(data$split)*2, each=K) * rep(c(-1,  1, 1), n) * rep(data$width / 2, each=K),
     y = rep(data$y, each=K) + rep(3-as.integer(data$split)*2, each=K) * rep(c(-1, -1, 1), n) * rep(data$height / 2, each=K),
     group = rep(seq_len(n), each=K)
   )
   new_data <- cbind(new_data, data[rep(seq_len(n), each = K), setdiff(colnames(data), c("x", "y", "group")), drop = FALSE])
   new_data
  },

  default_aes = aes(fill = "grey20", colour = NA, linewidth = 0.1, linetype = 1,
                   alpha = NA, width = NA, height = NA),

  required_aes = c("x", "y", "split"),

  draw_key = draw_key_split_tile
)

#' @export
scale_split <- function(..., scale_name="scale_direction", palette = function(n) if(n>2) error(paste0(scale_name, " can handle at most 2 levels")) else c(FALSE, TRUE) ) discrete_scale(aesthetics = "split", scale_name=scale_name, palette = palette, ... )


#df <- reshape2::melt(UCBAdmissions)
#ggplot(data = df, aes(x=Dept , y= Admit, fill = value, split=Gender, width=0.9)) +geom_split_tile(size=2) + scale_split()
#
#df <- reshape2::melt(HairEyeColor)
#ggplot(data = df, aes(x=Hair , y= Eye, fill = value, split=Sex)) + geom_split_tile() + expand_limits(fill=0) +scale_fill_viridis_c() + scale_split()
