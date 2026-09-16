# Train elastic-net GReX models for one ibeji training set.
#
# Usage: Rscript src/03_train_grex.R <set> <chromosomes: e.g. 22 or 1-22> <cores>
#
# For each gene: cis SNPs within +/-500 kb of the gene coordinate with MAF >= 0.05
# in the training set; elastic net (alpha = 0.5); lambda by 5-fold CV;
# performance by nested 5-fold CV (R^2 = squared Pearson correlation of held-out
# predictions). Final weights come from a fit on all samples. Weights are per
# ALT allele dosage on the original genotype scale.
#
# Why 500 kb, MAF 0.05 and 5 inner folds: a pilot with +/-1 Mb, MAF 0.01 and 10 inner
# folds ran more than 18 minutes on chr22 (633 genes) for YRI87 on 4 physical cores,
# which projects to over 10 hours per training set. Yoruba 1 Mb windows hold ~10,000
# SNVs. MAF 0.05 is also the more defensible floor at n = 87, where a 1% variant is
# carried only 2 to 8 times. The same settings apply to every training set.
suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(pgenlibr)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
set <- args[1]
chr_arg <- if (length(args) >= 2) args[2] else "1-22"
cores <- if (length(args) >= 3) as.integer(args[3]) else 4L
chroms <- if (grepl("-", chr_arg)) do.call(seq, as.list(as.integer(strsplit(chr_arg, "-")[[1]]))) else as.integer(strsplit(chr_arg, ",")[[1]])

root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
gdir <- file.path(root, "data", "processed", "geno")
edir <- file.path(root, "data", "processed", "expr")
sdir <- file.path(root, "data", "processed", "samples")
odir <- file.path(root, "results", "models", set)
dir.create(odir, showWarnings = FALSE, recursive = TRUE)

WINDOW <- 5e5
MAF_MIN <- 0.05
ALPHA <- 0.5
OUTER <- 5L
INNER <- 5L
SEED <- 2026L

pvar_path <- file.path(gdir, "geuv445.pvar")
pgen_path <- file.path(gdir, "geuv445.pgen")
# Per-chromosome variant index (global pgen row, POS, ID), written once by awk from
# geuv445.pvar. Loading one chromosome at a time keeps memory small: the full 17M-row
# table held in 7 forked workers pushed a 16 GB machine into heavy swapping.
load_chr_index <- function(chr) {
  d <- fread(file.path(gdir, "chr_index", sprintf("chr%d.tsv", chr)), col.names = c("idx", "POS", "ID"))
  stopifnot(!is.unsorted(d$POS))
  d
}
# Positions are sorted, so a gene's cis window is a binary search.
cis_window <- function(d, coord) {
  lo <- findInterval(coord - WINDOW - 1, d$POS) + 1L
  hi <- findInterval(coord + WINDOW, d$POS)
  if (hi < lo) d[0] else d[lo:hi]
}

psam <- fread(file.path(gdir, "geuv445.psam"))
setnames(psam, "#IID", "IID")
ids <- fread(file.path(sdir, paste0(set, ".txt")))[["#IID"]]
sample_idx <- match(ids, psam$IID)
stopifnot(!anyNA(sample_idx))
# pgenlibr returns rows in file order for a sorted sample subset
ord <- order(sample_idx)
sample_idx_sorted <- sample_idx[ord]
ids_sorted <- ids[ord]

resid <- fread(file.path(edir, paste0(set, ".residuals.tsv.gz")))
annot <- fread(file.path(edir, "gene_annotation.tsv"))
annot <- annot[gene %in% resid$gene]

fit_gene <- function(y, x) {
  n <- length(y)
  set.seed(SEED)
  outer <- sample(rep(seq_len(OUTER), length.out = n))
  pred <- rep(NA_real_, n)
  for (k in seq_len(OUTER)) {
    tr <- outer != k
    inner <- sample(rep(seq_len(INNER), length.out = sum(tr)))
    cv <- tryCatch(cv.glmnet(x[tr, , drop = FALSE], y[tr], alpha = ALPHA, foldid = inner),
                   error = function(e) NULL)
    pred[!tr] <- if (is.null(cv)) mean(y[tr]) else as.numeric(predict(cv, x[!tr, , drop = FALSE], s = "lambda.min"))
  }
  r <- suppressWarnings(cor(pred, y))
  if (is.na(r)) { r2 <- 0; pval <- 1 } else {
    r2 <- r^2
    pval <- suppressWarnings(cor.test(pred, y)$p.value)
  }
  inner_all <- sample(rep(seq_len(INNER), length.out = n))
  cv_all <- tryCatch(cv.glmnet(x, y, alpha = ALPHA, foldid = inner_all), error = function(e) NULL)
  if (is.null(cv_all)) return(list(r2 = r2, pval = pval, lambda = NA_real_, w = numeric(0)))
  b <- as.matrix(coef(cv_all, s = "lambda.min"))[-1, 1]
  list(r2 = r2, pval = pval, lambda = cv_all$lambda.min, w = b[b != 0])
}

run_chunk <- function(genes_chunk) {
  pg <- NewPgen(pgen_path, sample_subset = sample_idx_sorted)
  on.exit(ClosePgen(pg))
  buf_summary <- vector("list", nrow(genes_chunk))
  buf_weights <- vector("list", nrow(genes_chunk))
  for (i in seq_len(nrow(genes_chunk))) {
    g <- genes_chunk[i]
    win <- cis_window(chr_index, g$coord)
    y <- as.numeric(unlist(resid[gene == g$gene, ..ids_sorted]))
    if (nrow(win) < 2) {
      buf_summary[[i]] <- data.table(gene = g$gene, n_window = nrow(win), n_cis = 0L, n_model = 0L, cv_r2 = NA_real_, cv_pval = NA_real_, lambda = NA_real_)
      next
    }
    x <- ReadList(pg, win$idx, meanimpute = FALSE)
    af <- colMeans(x) / 2
    maf <- pmin(af, 1 - af)
    keep <- maf >= MAF_MIN
    x <- x[, keep, drop = FALSE]
    colnames(x) <- win$ID[keep]
    if (ncol(x) < 2) {
      buf_summary[[i]] <- data.table(gene = g$gene, n_window = nrow(win), n_cis = ncol(x), n_model = 0L, cv_r2 = NA_real_, cv_pval = NA_real_, lambda = NA_real_)
      next
    }
    f <- fit_gene(y, x)
    buf_summary[[i]] <- data.table(gene = g$gene, n_window = nrow(win), n_cis = ncol(x), n_model = length(f$w),
                                   cv_r2 = f$r2, cv_pval = f$pval, lambda = f$lambda)
    if (length(f$w)) buf_weights[[i]] <- data.table(gene = g$gene, varID = names(f$w), weight = unname(f$w))
  }
  cat(sprintf("[%s] %s chunk of %d genes done (first %s)\n", set, format(Sys.time(), "%H:%M:%S"),
              nrow(genes_chunk), genes_chunk$gene[1]), file = file.path(root, "logs", sprintf("progress_%s.log", set)), append = TRUE)
  list(summary = rbindlist(buf_summary), weights = rbindlist(buf_weights))
}

setkey(resid, gene)
for (chr_now in chroms) {
  chr <- chr_now
  if (file.exists(file.path(odir, sprintf("chr%d.weights.tsv", chr_now)))) {
    message(sprintf("[%s] chr%d already done, skipping", set, chr_now))
    next
  }
  t0 <- Sys.time()
  genes <- annot[annot$chr == chr_now]
  if (!nrow(genes)) next
  chr_index <- load_chr_index(chr_now)
  chunks <- split(genes, ceiling(seq_len(nrow(genes)) / 25))
  res <- mclapply(chunks, run_chunk, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) stop(sprintf("chr%d: %d chunks failed: %s", chr, sum(bad), as.character(res[bad][[1]])))
  fwrite(rbindlist(lapply(res, `[[`, "summary")), file.path(odir, sprintf("chr%d.summary.tsv", chr)), sep = "\t")
  fwrite(rbindlist(lapply(res, `[[`, "weights")), file.path(odir, sprintf("chr%d.weights.tsv", chr)), sep = "\t")
  message(sprintf("[%s] chr%d: %d genes in %.1f min", set, chr, nrow(genes), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
