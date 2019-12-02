
#' @importFrom magrittr `%<>%` `%>%` 
#' @import rlang
NULL


#source: https://github.com/lonsbio/fastqe
fastq_emoji_map_binned= c(
  '!'= 'no_entry_sign',
  '"'= 'no_entry_sign',
  
  #2–9 6
  '#'= 'skull',
  '$'= 'skull',
  '%'= 'skull',
  '&'= 'skull',
  '\''= 'skull',
  '('= 'skull',
  ')'= 'skull',
  '*'= 'skull',
  
  #10–19 15
  '+'= 'poop' ,
  ','= 'poop' ,
  '-'= 'poop' ,
  '.'= 'poop' ,
  '/'= 'poop' ,
  '0'= 'poop' ,
  '1'= 'poop' ,
  '2'= 'poop' ,
  '3'= 'poop' ,
  '4'= 'poop' ,
  
  #20–24 22
  '5'= 'warning',
  '6'= 'warning',
  '7'= 'warning',
  '8'= 'warning',
  '9'= 'warning',
  
  
  #25–29 27
  ':'= 'smile',
  ';'= 'smile',
  '<'= 'smile',
  '='= 'smile',
  '>'= 'smile',
  
  
  #30–34 33
  '?'= 'laughing',
  '@'= 'laughing',
  'A'= 'laughing',
  'B'= 'laughing',
  'C'= 'laughing',
  
  #35–39 37
  'D'= 'sunglasses',
  'E'= 'sunglasses',
  'F'= 'sunglasses',
  'G'= 'sunglasses',
  'H'= 'sunglasses',
  
  #≥ 40 40
  'I'= 'heart_eyes',
  'J'= 'heart_eyes'
)
fastq_emoji_map_binned <- sapply(fastq_emoji_map_binned, emo::ji)

#' @export
fastq2emoji <- function(x) do.call(paste0,lapply(tstrsplit(x, split=""), function(x) fastq_emoji_map_binned[x]))


#' @export
seq2vec <- function(seq) seq %>% str_replace_all("A", "4000") %>% str_replace_all("C", "0400") %>% str_replace_all("G", "0040") %>% str_replace_all("T", "0004") %>% str_replace_all("N", "1111") %>% tstrsplit("") %>% lapply(as.integer) %>% do.call(cbind, .) %>% `/`(4)


#' @author Aaron Lun
#' @export
setGeneric(".my_colVars", function(x, ...) standardGeneric(".my_colVars"))

#' @author Aaron Lun
#' @importFrom Matrix t rowSums
#' @export
setMethod(".my_colVars", "ANY", function(x, center=NULL) {
  if (is.null(center)) center <- Matrix::colMeans(x)
  y <- t(x) - center
  Matrix::rowSums(y^2)/(ncol(y)-1)
}) 

#' @author Aaron Lun
#' @author Jan Gleixner
#' @importFrom Matrix t colSums
#' @export
setMethod(".my_colVars", "dgCMatrix", function(x, center=NULL) {
  if (is.null(center)) center <- Matrix::colMeans(x)
  nzero <- diff(x@p)
  expanded <- rep(center, nzero)
  x@x <- (x@x - expanded)^2
  (Matrix::colSums(x) + (nrow(x)-nzero) * center^2)/(nrow(x)-1)
})


#' @export
setGeneric(".my_rowVars", function(x, ...) standardGeneric(".my_rowVars"))

#' @importFrom Matrix t rowSums
#' @export
setMethod(".my_rowVars", "ANY", function(x, center=NULL) {
  if (is.null(center)) center <- Matrix::rowMeans(x) else stopifnot(length(center) == nrow(x))
  y <- x - center
  Matrix::rowSums(y^2)/(ncol(y)-1)
}) 

#' @importFrom Matrix t colSums
#' @export
setMethod(".my_rowVars", "dgCMatrix", function(x, center=NULL) {
  if (is.null(center)) center <- Matrix::rowMeans(x) else stopifnot(length(center) == nrow(x))
  nzero <- base::tabulate(x@i+1L, nrow(x))
  x@x <- (x@x - center[x@i+1L])^2
  (Matrix::rowSums(x) + (ncol(x)-nzero) * center^2)/(ncol(x)-1)
})

#' @export
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

FindLowNoiseFeatures <- function (object, ...) UseMethod(generic = "FindVariableFeatures", object = object)

FindLowNoiseFeatures.Seurat <- function (object, assay = NULL, selection.method = "vst", loess.span = 0.3, 
          clip.max = "auto", mean.function = FastExpMean, dispersion.function = FastLogVMR, 
          num.bin = 20, binning.method = "equal_width", nfeatures = 2000,  verbose = TRUE, 
          ...) 
{
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  assay.data <- FindLowNoiseFeatures(object = assay.data, nfeatures = nfeatures, verbose = verbose, ...)
  object[[assay]] <- assay.data
  object <- LogSeuratCommand(object = object)
  return(object)
}
#' 
#' eurat:::FindVariableFeatures.Assay
#' function (object, selection.method = "vst", loess.span = 0.3, 
#'           clip.max = "auto", mean.function = FastExpMean, dispersion.function = FastLogVMR, 
#'           num.bin = 20, binning.method = "equal_width", nfeatures = 2000, 
#'           mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf), verbose = TRUE, 
#'           ...) 
#' {
#'   if (length(x = mean.cutoff) != 2 || length(x = dispersion.cutoff) != 
#'       2) {
#'     stop("Both 'mean.cutoff' and 'dispersion.cutoff' must be two numbers")
#'   }
#'   if (selection.method == "vst") {
#'     data <- GetAssayData(object = object, slot = "counts")
#'     if (IsMatrixEmpty(x = data)) {
#'       warning("selection.method set to 'vst' but count slot is empty; will use data slot instead")
#'       data <- GetAssayData(object = object, slot = "data")
#'     }
#'   }
#'   else {
#'     data <- GetAssayData(object = object, slot = "data")
#'   }
#'   hvf.info <- FindVariableFeatures(object = data, selection.method = selection.method, 
#'                                    loess.span = loess.span, clip.max = clip.max, mean.function = mean.function, 
#'                                    dispersion.function = dispersion.function, num.bin = num.bin, 
#'                                    binning.method = binning.method, verbose = verbose, ...)
#'   object[[names(x = hvf.info)]] <- hvf.info
#'   hvf.info <- hvf.info[which(x = hvf.info[, 1, drop = TRUE] != 
#'                                0), ]
#'   if (selection.method == "vst") {
#'     hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized, 
#'                                decreasing = TRUE), , drop = FALSE]
#'   }
#'   else {
#'     hvf.info <- hvf.info[order(hvf.info$mvp.dispersion, decreasing = TRUE), 
#'                          , drop = FALSE]
#'   }
#'   selection.method <- switch(EXPR = selection.method, mvp = "mean.var.plot", 
#'                              disp = "dispersion", selection.method)
#'   top.features <- switch(EXPR = selection.method, mean.var.plot = {
#'     means.use <- (hvf.info[, 1] > mean.cutoff[1]) & (hvf.info[, 
#'                                                               1] < mean.cutoff[2])
#'     dispersions.use <- (hvf.info[, 3] > dispersion.cutoff[1]) & 
#'       (hvf.info[, 3] < dispersion.cutoff[2])
#'     rownames(x = hvf.info)[which(x = means.use & dispersions.use)]
#'   }, dispersion = head(x = rownames(x = hvf.info), n = nfeatures), 
#'   vst = head(x = rownames(x = hvf.info), n = nfeatures), 
#'   stop("Unkown selection method: ", selection.method))
#'   VariableFeatures(object = object) <- top.features
#'   vf.name <- ifelse(test = selection.method == "vst", yes = "vst", 
#'                     no = "mvp")
#'   vf.name <- paste0(vf.name, ".variable")
#'   object[[vf.name]] <- rownames(x = object[[]]) %in% top.features
#'   return(object)
#' }
#' 
#' #' @importFrom rlang %||%
#' #' @import Seurat
#' #' @importFrom magrittr %>% %<>%
#' #' @export
#' FindLowNoiseFeatures <- function(object, nfeatures = 5000, assay = NULL, subset=NULL){
#'   subset <-  subset %||% TRUE
#'   #assay <- assay %||% DefaultAssay(object)
#'   
#'   mf <- GetAssay(object, assay = assay)@meta.features
#'   X <- GetAssayData(object, assay=assay)[, subset, drop=FALSE]
#'   
#'   mf$mean_normalized <- Matrix::rowMeans(X)
#'   mf$variance_normalized <- .my_rowVars(X, center = mf$mean_normalized)
#'   mf$is_variable <- FALSE
#'   mf[GetAssay(object, assay = assay)@var.features, "is_variable_Seurat"] <- TRUE
#'   
#'   x <- seq(0.00, max(mf$mean_normalized, na.rm=TRUE), length.out = 100)
#'   y <- (Vectorize(function(.) var(log1p(rpois(n = 1000, lambda = rlnorm(n = 1000, meanlog = log(exp(.)-1), sdlog=sqrt(1/40))))), SIMPLIFY = TRUE))(x)
#'   expected_variance <- function(z) predict(loess(y~x), z)
#'   mf$expected_variance <- expected_variance(mf$mean_normalized)
#'   variable_genes <- rownames(mf)[order(-mf$variance_normalized / mf$expected_variance)[seq_len(nfeatures)]]
#'   variable_genes <- variable_genes[!is.na(variable_genes)]
#'   mf$is_variable_log1p <- FALSE
#'   mf[variable_genes, "is_variable_normalized"] <- TRUE
#'   
#'   object[[assay]]@var.features <- variable_genes
#'   object %<>% Seurat::AddMetaData(object[[assay]], mf)
#'   object
#' }

#' @import Seurat
#' @export
WeightData <- function(object, assay = NULL){
  assay <- assay %||% Seurat::DefaultAssay(object)
  
  X <- Seurat::GetAssayData(assay, slot="scale.data")
  
  X <- X * mf$variance_normalized / mf$expected_variance / mean( mf$variance_normalized / mf$expected_variance)

  Seurat::SetAssayData(object, slot="scale.data", new.data=X, assay = assay)

}

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




