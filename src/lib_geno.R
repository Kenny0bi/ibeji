# Shared helpers for reading ibeji genotypes and models.
suppressPackageStartupMessages({
  library(data.table)
  library(pgenlibr)
})

project_root <- function() {
  f <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
  normalizePath(file.path(dirname(f), ".."))
}

geno_file <- function(root, ext) file.path(root, "data", "processed", "geno", paste0("geuv445.", ext))

# Variant lookup without loading all 17M pvar rows: IDs are chr:pos:ref:alt, so each
# requested ID is resolved against a per-chromosome index (global pgen row, POS, ID)
# that is loaded once per chromosome and cached.
load_pvar <- function(root) {
  env <- new.env()
  env$root <- root
  env$cache <- list()
  env
}

variant_rows <- function(pvar, var_ids) {
  chr <- as.integer(sub(":.*$", "", var_ids))
  out <- vector("list", length(unique(chr)))
  for (k in seq_along(unique(chr))) {
    c_now <- unique(chr)[k]
    key <- as.character(c_now)
    if (is.null(pvar$cache[[key]])) {
      d <- fread(file.path(pvar$root, "data", "processed", "geno", "chr_index", sprintf("chr%d.tsv", c_now)),
                 col.names = c("idx", "POS", "ID"))
      setkey(d, ID)
      pvar$cache[[key]] <- d
    }
    out[[k]] <- pvar$cache[[key]][J(var_ids[chr == c_now]), nomatch = NULL]
  }
  rbindlist(out)
}

# pgenlibr needs a sorted 1-based sample subset; rows come back in file order.
sample_subset <- function(root, ids) {
  ps <- fread(geno_file(root, "psam"))
  setnames(ps, "#IID", "IID")
  idx <- match(ids, ps$IID)
  stopifnot(!anyNA(idx))
  o <- order(idx)
  list(idx = idx[o], ids = ids[o])
}

set_ids <- function(root, set) fread(file.path(root, "data", "processed", "samples", paste0(set, ".txt")))[["#IID"]]

# ALT-allele dosage matrix (samples x variants), columns named and ordered as var_ids.
read_dosage <- function(pg, pvar, var_ids) {
  rows <- variant_rows(pvar, var_ids)
  rows <- rows[order(idx)]
  x <- ReadList(pg, rows$idx, meanimpute = FALSE)
  colnames(x) <- rows$ID
  x[, var_ids[var_ids %in% rows$ID], drop = FALSE]
}

load_models <- function(root, set, r2_min = 0.01, p_max = 0.05) {
  d <- file.path(root, "results", "models", set)
  s <- rbindlist(lapply(list.files(d, "summary.tsv$", full.names = TRUE), fread))
  w <- rbindlist(lapply(list.files(d, "weights.tsv$", full.names = TRUE), fread))
  s[, significant := !is.na(cv_r2) & cv_r2 > r2_min & cv_pval < p_max & n_model > 0]
  list(summary = s, weights = w[gene %in% s[significant == TRUE, gene]])
}

load_residuals <- function(root, set) {
  r <- fread(file.path(root, "data", "processed", "expr", paste0(set, ".residuals.tsv.gz")))
  m <- as.matrix(r[, -1])
  rownames(m) <- r$gene
  m
}

safe_cor <- function(x) {
  r <- suppressWarnings(cor(x))
  r[is.na(r)] <- 0
  diag(r) <- 1
  r
}
