#' Convert "R" colors to colorspaces' sRGB colors
#'
#' This function takes R color strings and converts them into colorspaces::sRGB values.
#' It utilizes grDevices::col2rgb for conversion.
#' @param col A character vector of color names or hexadecimal color codes.
#' @return An object of class ["sRGB"][colorspace::sRGB] containing the red, green, and blue components of the color.
#' @export
#' @examples
#' col2sRGB("red")
#' col2sRGB("#FF0000")
col2sRGB <- function(col) {
  col <- grDevices::col2rgb(col)/255
  colorspace::sRGB(col["red",], col["green",], col["blue",])
}

#' `r lifecycle::badge('deprecated')`
#' @export
col2RGB <- function(col) {
  lifecycle::deprecate_soft("2024-02-20", "col2RGB()", "col2sRGB()") ;
  col2sRGB(col)
}

#' Mix two "R" colors
#'
#' `mixcolorr` mixes two R color strings in a given proportion. It is a wrapper
#' around [colorspace::mixcolor()] that accepts and returns [R color strings][grDevices::col2rgb] like
#' `"red"` or`"#FF0000"`. It is vectorized.
#' @param alpha A numeric vector of length 1 or n with values between 0 (only col1) and 1 (only col2) indicating the mixing coefficient.
#' @param col1 The first color to mix, as a name or hexadecimal code.
#' @param col2 The second color to mix, as a name or hexadecimal code.
#' @param ... other arguments passed to [colorspace::mixcolor()]
#' @return A hexadecimal color code representing the mixed color.
#' @export
#' @seealso [colorspace::mixcolor()] [grDevices::colors()] [grDevices::col2rgb()]
#' @examples
#' mixcolorr(0.8, "red", "gray")
mixcolorr <- function(alpha, col1, col2, ...) {
  col1 <- col2sRGB(col1)
  col2 <- col2sRGB(col2)
  ret <- colorspace::mixcolor(alpha, col1, col2, ...)
  colorspace::hex(ret, fixup = TRUE)
}














if(FALSE) {

observed_celltypes <- dt[, sort(unique(cell_type))]



tmp <- str_split(observed_celltypes, " & ")
names(observed_celltypes) <- observed_celltypes
colors_of_cell_types <- sapply(observed_celltypes, function(cell_types) {
  cell_types <- str_split(cell_types, " & ")[[1]]
  colors <- colors_of_cell_type[cell_types]
  colors <- lapply(colors, colorspace::hex2RGB)
  color <- Reduce(x=colors, f= \(a,b) list(colorspace::mixcolor(alpha=1/a[[2]], a[[1]], b, where="polarLUV"), a[[2]]+1), init=list(colorspace::hex2RGB("#aaaaaa"),1))[[1]]
  colorspace::hex(color, fixup = TRUE)
} )

colorspace::hex(colorspace::mixcolor(0.5, colorspace::hex2RGB("#ce5a02"), colorspace::hex2RGB("#e6ab02")), fixup = TRUE)


labels_of_cell_type <- c(
  "pseudo bulk",
  "dead cells",
  "Stem cells",
  "TA cells",
  "EC prog.",
  "Enterocytes",
  "EECs",
  "Paneth cells",
  "Goblet &\nPaneth cells",
  "Goblet cells",
  "Tuft cells")
legend_dot_size <- 2

single_cell_types <- (function(x) (!grepl(pattern = "&",x))|(grepl(pattern = "Goblet",x) & grepl(pattern = "Paneth",x)) )(names(colors_of_cell_types))

scale_color_cell_type <- function (
    name = "cell type",
    breaks = names(colors_of_cell_types)[single_cell_types],
    labels =  str_wrap(names(colors_of_cell_types), width=20)[single_cell_types], #labels_of_cell_type, #names(colors_of_cell_types), #
    values = colors_of_cell_types, limits = names(colors_of_cell_types), na.value = "#808080",
    guide = guide_legend(override.aes = list(size = legend_dot_size)),
    ...)
  scale_color_manual(name = name, breaks = breaks, limits = limits,
                     labels = labels, values = values, na.value = na.value, #guide = guide,
                     ...)

ggplot(dt[], aes(x=UMAP_1, y=UMAP_2, color = cell_type)) +
  geom_point() +
  scale_color_cell_type() #labels=\(x) str_wrap(x, width=20)) +
NULL






GeomPointHCL <- ggproto("GeomPointHCL", GeomPoint)
GeomPointHCL$default_aes$hue <- 1
GeomPointHCL$default_aes$chroma <- 1
GeomPointHCL$default_aes$luminance <- 0.6

geom_point_hcl <- function(mapping = NULL, data = NULL, stat = "identity", position = "identity",
                           ..., na.rm = FALSE, show.legend = NA, inherit.aes = TRUE)
{
  layer(data = data, mapping = mapping, stat = stat, geom = GeomPointHCL,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list2(na.rm = na.rm, ...))
}

ggplot(dt[], aes(x=UMAP_1, y=UMAP_2, luminance = 0.6,
                 chroma=(1-(1-cell_type_score^2.5)*0.999)*0.998,
                 color = stage(cell_type, after_scale = colorspace::hex(colorspace::polarLUV(luminance*100, chroma*132, as(colorspace::hex2RGB(colorspace::darken(color,0.0001)), "polarLUV")@coords[, "H"]), fixup=TRUE))
)) +
  geom_point_hcl() +
  sir::scale_color_cell_type()+
  # scale_color_discrete(labels=\(x) str_wrap(x, width=20)) +
  NULL

}
