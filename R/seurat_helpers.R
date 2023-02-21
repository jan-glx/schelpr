
#' extract and join meta information and reductions (PCA, TSNE etc) into a single, randomly ordered data.table
#' @export
meta_and_reductions <- function(si, max_n_dim_embeddings = 2) {
  dt <- if("cell" %in% colnames(si@meta.data)) data.table(si@meta.data) else data.table(si@meta.data, keep.rownames = "cell")
  dt <- cbind(dt, do.call(cbind, setNames(lapply(si@reductions, function(reduction) reduction@cell.embeddings[, seq_len(min(ncol(reduction@cell.embeddings), max_n_dim_embeddings))] %>% data.table), NULL)))
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
  features <- lapply(features, intersect, usefull_genes)
  features <- lapply(features, intersect, rownames(data))

  lapply(features, function(features) {
    dat <- t(data[features, ])
    nn <- pmax(1, pmin(5, ncol(dat)-2)) # number of PCs to look at
    mean_scores <- Matrix::rowMeans(dat)
    tryCatch({
      res <- tryCatch(irlba::irlba(dat, scale = sqrt(gene_vars[features]), center = Matrix::colMeans(dat), nv=nn),
                      warning = function(w) {svd(scale(dat), nu=nn, nv=nn)},
                      error = function(e) {svd(scale(dat), nu=nn, nv=nn)})
      scores <- with(res, t(t(u)*sign(colMeans(v))))
      cors <- cor(scores, mean_scores)
      if(which.max(cors)!=1) message("The first principal component has not the highest correlation with the average of features (", paste0(features, collapse = ", "), ") average (correlations:", sprintf(" %.3f", cors), ")")
      scores[, 1]
    }, error = function(e) {warning("pca_score_computation failed with \"", e, "\", returning NAs.") ;rep(NA_real_, ncol(data))})
  }) %>% as.data.frame(row.names = colnames(data)) %>%
    setNames(names(features))
}

#' @export
aggregate_log_normalized_expression <- function(x, group, group_name = deparse1(substitute(group))) {
  group <- as.factor(group)
  total_counts_of_cell <- Matrix::colSums(x)
  x <- Matrix::t(x)
  mean_counts_of_gene <- Matrix::colMeans(x)

  #subset to detected genes
  x <- x[, mean_counts_of_gene>0]

  total_counts_of_group_x_gene <- as.matrix(Matrix.utils::aggregate.Matrix(x, group))
  n_cells_of_group <- tabulate(group)
  total_counts_of_group <- Matrix::rowSums(total_counts_of_group_x_gene)

  # log1p normalize
  scaling_factor <-  mean(total_counts_of_cell)
  size_factor <- scaling_factor / total_counts_of_cell
  x <- log1p(x * size_factor)


  x_mean <- as.matrix(Matrix.utils::aggregate.Matrix(x, group) / n_cells_of_group)
  x_var <- as.matrix(aggregate_var.Matrix(x, group))
  n_genes <- ncol(x_mean)

  agg_dt <- data.table(
    group = rep(rownames(x_mean), n_genes),
    gene = rep(colnames(x_mean), each = nrow(x_mean)),
    expression_mean = as.vector(x_mean),
    expression_variance = as.vector(x_var),
    count_sum = as.vector(total_counts_of_group_x_gene),
    total_count_sum = rep(total_counts_of_group, n_genes),
    n_cells_of_group = rep(n_cells_of_group, n_genes),
    scaling_factor = scaling_factor
  )

  agg_dt[, c("expression_log1p_lower", "expression_log1p_upper", "expression_log1p_mean") :=
           data.table(binom::binom.exact(count_sum, total_count_sum, conf.level = diff(pnorm(c(-0.5, 0.5)))))[
             , data.table(log1p(cbind(lower, upper, mean) * scaling_factor))], by=.(group, total_count_sum)]
  agg_dt[, expression_counting_variance := (expression_log1p_upper - expression_log1p_lower)^2 * n_cells_of_group , by=.(group, total_count_sum)]
  setnames(agg_dt, "group", group_name)
  setnames(agg_dt, "n_cells_of_group", paste0("n_cells_of_", group_name))
  agg_dt[]
}
