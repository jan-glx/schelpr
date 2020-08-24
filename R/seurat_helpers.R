
#' extract and join meta information and reductions (PCA, TSNE etc) into a single, randomly ordered data.table
#' @export
meta_and_reductions <- function(si) {
  dt <- if("cell" %in% colnames(si@meta.data)) data.table(si@meta.data) else data.table(si@meta.data, keep.rownames = "cell")
  dt <- cbind(dt, do.call(cbind, setNames(lapply(si@reductions, function(reduction) reduction@cell.embeddings[, 1:2] %>% data.table), NULL)))
  shuffle(dt)
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
#
# eurat:::FindVariableFeatures.Assay
# function (object, selection.method = "vst", loess.span = 0.3,
#           clip.max = "auto", mean.function = FastExpMean, dispersion.function = FastLogVMR,
#           num.bin = 20, binning.method = "equal_width", nfeatures = 2000,
#           mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf), verbose = TRUE,
#           ...)
# {
#   if (length(x = mean.cutoff) != 2 || length(x = dispersion.cutoff) !=
#       2) {
#     stop("Both 'mean.cutoff' and 'dispersion.cutoff' must be two numbers")
#   }
#   if (selection.method == "vst") {
#     data <- GetAssayData(object = object, slot = "counts")
#     if (IsMatrixEmpty(x = data)) {
#       warning("selection.method set to 'vst' but count slot is empty; will use data slot instead")
#       data <- GetAssayData(object = object, slot = "data")
#     }
#   }
#   else {
#     data <- GetAssayData(object = object, slot = "data")
#   }
#   hvf.info <- FindVariableFeatures(object = data, selection.method = selection.method,
#                                    loess.span = loess.span, clip.max = clip.max, mean.function = mean.function,
#                                    dispersion.function = dispersion.function, num.bin = num.bin,
#                                    binning.method = binning.method, verbose = verbose, ...)
#   object[[names(x = hvf.info)]] <- hvf.info
#   hvf.info <- hvf.info[which(x = hvf.info[, 1, drop = TRUE] !=
#                                0), ]
#   if (selection.method == "vst") {
#     hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized,
#                                decreasing = TRUE), , drop = FALSE]
#   }
#   else {
#     hvf.info <- hvf.info[order(hvf.info$mvp.dispersion, decreasing = TRUE),
#                          , drop = FALSE]
#   }
#   selection.method <- switch(EXPR = selection.method, mvp = "mean.var.plot",
#                              disp = "dispersion", selection.method)
#   top.features <- switch(EXPR = selection.method, mean.var.plot = {
#     means.use <- (hvf.info[, 1] > mean.cutoff[1]) & (hvf.info[,
#                                                               1] < mean.cutoff[2])
#     dispersions.use <- (hvf.info[, 3] > dispersion.cutoff[1]) &
#       (hvf.info[, 3] < dispersion.cutoff[2])
#     rownames(x = hvf.info)[which(x = means.use & dispersions.use)]
#   }, dispersion = head(x = rownames(x = hvf.info), n = nfeatures),
#   vst = head(x = rownames(x = hvf.info), n = nfeatures),
#   stop("Unkown selection method: ", selection.method))
#   VariableFeatures(object = object) <- top.features
#   vf.name <- ifelse(test = selection.method == "vst", yes = "vst",
#                     no = "mvp")
#   vf.name <- paste0(vf.name, ".variable")
#   object[[vf.name]] <- rownames(x = object[[]]) %in% top.features
#   return(object)
# }
#
# # @importFrom rlang %||%
# # @import Seurat
# # @importFrom magrittr %>% %<>%
# # @export
# FindLowNoiseFeatures <- function(object, nfeatures = 5000, assay = NULL, subset=NULL){
#   subset <-  subset %||% TRUE
#   #assay <- assay %||% DefaultAssay(object)
#
#   mf <- GetAssay(object, assay = assay)@meta.features
#   X <- GetAssayData(object, assay=assay)[, subset, drop=FALSE]
#
#   mf$mean_normalized <- Matrix::rowMeans(X)
#   mf$variance_normalized <- .my_rowVars(X, center = mf$mean_normalized)
#   mf$is_variable <- FALSE
#   mf[GetAssay(object, assay = assay)@var.features, "is_variable_Seurat"] <- TRUE
#
#   x <- seq(0.00, max(mf$mean_normalized, na.rm=TRUE), length.out = 100)
#   y <- (Vectorize(function(.) var(log1p(rpois(n = 1000, lambda = rlnorm(n = 1000, meanlog = log(exp(.)-1), sdlog=sqrt(1/40))))), SIMPLIFY = TRUE))(x)
#   expected_variance <- function(z) predict(loess(y~x), z)
#   mf$expected_variance <- expected_variance(mf$mean_normalized)
#   variable_genes <- rownames(mf)[order(-mf$variance_normalized / mf$expected_variance)[seq_len(nfeatures)]]
#   variable_genes <- variable_genes[!is.na(variable_genes)]
#   mf$is_variable_log1p <- FALSE
#   mf[variable_genes, "is_variable_normalized"] <- TRUE
#
#   object[[assay]]@var.features <- variable_genes
#   object %<>% Seurat::AddMetaData(object[[assay]], mf)
#   object
# }

#' @import Seurat
#' @export
WeightData <- function(object, assay = NULL){
  assay <- assay %||% Seurat::DefaultAssay(object)

  X <- Seurat::GetAssayData(assay, slot="scale.data")

  X <- X * mf$variance_normalized / mf$expected_variance / mean( mf$variance_normalized / mf$expected_variance)

  Seurat::SetAssayData(object, slot="scale.data", new.data=X, assay = assay)

}


#' @export
mean_scores <- function(data, features) {
  lapply(features, function(features) {
    colMeans(data[features %>% .[. %in% rownames(data)], ])
  }) %>% as.data.frame(row.names = colnames(data)) %>%
    setNames(names(features))
}

#' @export
pca_scores <- function(data, features) {
  gene_vars <- schelpr::.my_rowVars(data)
  usefull_genes <-  names(gene_vars)[gene_vars>0]
  features <- lapply(features, function(x) x[x %in% usefull_genes])

  lapply(features, function(features) {

    dat <- t(data[features %>% .[. %in% rownames(data)], ])
    tryCatch({
      res <- tryCatch(irlba::irlba(dat, nv=1),
                    warning = function(w) svd(dat, nu=1, nv=1),
                    error=function(e) svd(dat, nu=1, nv=1))
      with(res, u*sign(mean(v)))
    }, error = function(e) rep(NA_real_, ncol(data)))
  }) %>% as.data.frame(row.names = colnames(data)) %>%
    setNames(names(features))
}
