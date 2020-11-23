
#' @export
fread_mtx <- function(mtx_file) {
  dt <- fread(mtx_file, header = TRUE)
  meta <- as.integer(colnames(dt))
  meta <- as.list(setNames(meta, c("n_row", "n_col", "row_count")))
  setnames(dt, c("row", "col", "val"))
  with(meta, {
    dt[, {
      assertthat::are_equal(.N, row_count)
      assertthat::assert_that(min(row) >= 1L)
      assertthat::assert_that(min(col) >= 1L)
      assertthat::assert_that(max(row) <= n_row)
      assertthat::assert_that(max(col) <= n_col)
    }]
  })
  dt[]
}

#' @export
fread_cellSNP <- function(input_dir) {
  base_vcf <- vcfR::read.vcfR(fs::path(input_dir, "cellSNP.base.vcf.gz"))
  SNPs_dt <- data.table(base_vcf@fix)


  dt <- fread_mtx(fs::path(input_dir, "cellSNP.tag.DP.mtx"))
  setnames(dt, "val", "DP")
  ad <- fread_mtx(fs::path(input_dir, "cellSNP.tag.AD.mtx"))
  dt[, `:=`(AD = 0, OTH = 0)]
  dt[ad, AD := val, on = .(row, col)]
  oth <- fread_mtx(fs::path(input_dir, "cellSNP.tag.OTH.mtx"))
  dt[oth, OTH := val, on = .(row, col)]

  cells_dt <- fread(fs::path(input_dir, "cellSNP.samples.tsv"), header = FALSE, col.names = "cell")
  cells_dt[, I := .I]
  SNPs_dt[, I := .I]
  dt <- cells_dt[dt, on = .(I = col)]
  dt[, I := NULL]
  dt <- SNPs_dt[dt, on = .(I = row)]
  dt[, I := NULL]
  dt[]
}
