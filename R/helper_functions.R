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
#' @importFrom magrittr %>%
seq2vec <- function(seq) seq %>% str_replace_all("A", "4000") %>% str_replace_all("C", "0400") %>% str_replace_all("G", "0040") %>% str_replace_all("T", "0004") %>% str_replace_all("N", "1111") %>% tstrsplit("") %>% lapply(as.integer) %>% do.call(cbind, .) %>% `/`(4)


#' @export
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' @importFrom rlang %||%
#' @export
FindLowNoiseFeatures <- function(object, nfeatures = 5000, assay = NULL, subset=NULL){
  subset <-  subset %||% TRUE
  assay <- assay %||% Seurat::DefaultAssay(object)
  
  mf <- Seurat::GetAssay(object, assay = assay)@meta.features
  X <- Seurat::GetAssayData(object, assay=assay)[, subset, drop=FALSE]
  
  mf$mean_normalized <- Matrix::rowMeans(X)
  mf$variance_normalized <- Matrix::rowMeans(X^2) - Matrix::rowMeans(X)^2
  mf$is_variable <- FALSE
  mf[object[[assay]]@var.features, "is_variable_Seurat"] <- TRUE
  
  x <- seq(0.00, max(mf$mean_normalized, na.rm=TRUE), length.out = 100)
  y <- (Vectorize(function(.) var(log1p(rpois(n = 1000, lambda = rlnorm(n = 1000, meanlog = log(exp(.)-1), sdlog=sqrt(1/40))))), SIMPLIFY = TRUE))(x)
  expected_variance <- function(z) predict(loess(y~x), z)
  mf$expected_variance <- expected_variance(mf$mean_normalized)
  variable_genes <- rownames(mf)[order(-mf$variance_normalized / mf$expected_variance)[seq_len(nfeatures)]]
  variable_genes <- variable_genes[!is.na(variable_genes)]
  mf$is_variable_log1p <- FALSE
  mf[variable_genes, "is_variable_normalized"] <- TRUE
  
  object[[assay]]@var.features <- variable_genes
  object %<>% Seurat::AddMetaData(object[[assay]], mf)
  object
}

#' @export
WeightData <- function(object, assay = NULL){
  assay <- assay %||% Seurat::DefaultAssay(object)
  
  X <- Seurat::GetAssayData(assay, slot="scale.data")
  
  X <- X * mf$variance_normalized / mf$expected_variance / mean( mf$variance_normalized / mf$expected_variance)

  Seurat::SetAssayData(object, slot="scale.data", new.data=X, assay = assay)

}




#' @export
render_separately <- function(...)  callr::r(function(...) rmarkdown::render(..., envir = globalenv()), args = list(...), show = TRUE)

#' @export
knit_all_formats <- function(rmd_file, params_=list(), figure_formats=c("png", "pdf", "svg")) {
  rmd_file_plain <- fs::path_ext_remove(fs::path_file(rmd_file))
  
  params_string <- gsub(" +", " ", paste(deparse(params_[seq_along(params_)]), collapse=""))
  params_hash <- substr(digest::digest(params_string), 1, 8)
  results_dir <- fs::path("results", fs::path_sanitize(rmd_file_plain), fs::path_sanitize(params_hash))
  fs::dir_create(results_dir)
  cat(params_string, file = fs::path(results_dir, "params.txt"))
  cat(params_string, file = fs::path(results_dir,  paste0(fs::path_sanitize(params_string), ext="txt")))
  
  lapply(
    figure_formats,
    function(figure_format) render_separately(
      input = rmd_file,
      params = c(params_, results_dir=results_dir),
      output_file = paste0(rmd_file_plain, "_", figure_format),
      output_format = rmarkdown::html_document(dev=figure_format, keep_md=TRUE),
      output_dir = fs::path(results_dir), 
      encoding = 'UTF-8'
    )
  )
  results_dir
}











