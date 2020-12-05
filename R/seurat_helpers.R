
#' extract and join meta information and reductions (PCA, TSNE etc) into a single, randomly ordered data.table
#' @export
meta_and_reductions <- function(si) {
  dt <- if("cell" %in% colnames(si@meta.data)) data.table(si@meta.data) else data.table(si@meta.data, keep.rownames = "cell")
  dt <- cbind(dt, do.call(cbind, setNames(lapply(si@reductions, function(reduction) reduction@cell.embeddings[, 1:2] %>% data.table), NULL)))
  shuffle(dt)
}

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
