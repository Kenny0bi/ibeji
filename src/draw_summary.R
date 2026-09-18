## Every number in the paper that depends on how many European draws have been trained.
## Also re-tests the margin claim: the smallest noise floor must still exceed the
## largest frequency plus LD sum, and it reports PASS or FAIL.
##
## Usage: Rscript src/draw_summary.R

suppressPackageStartupMessages(library(data.table))

find_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f[1L])) else getwd()
  normalizePath(file.path(d, ".."), mustWork = FALSE)
}
ROOT   <- find_root()
MODELS <- file.path(ROOT, "results", "models")
DEC    <- file.path(ROOT, "results", "decomposition")
YRI    <- "YRI87"

complete_draws <- function() {
  dirs <- list.dirs(MODELS, recursive = FALSE)
  dirs <- dirs[grepl("EUR87_r", basename(dirs))]
  keep <- vapply(dirs, function(d) length(list.files(d, pattern = "^chr.*\\.summary\\.tsv$")) == 22L, logical(1))
  sort(basename(dirs[keep]))
}

yield_row <- function(set_name) {
  files <- sort(list.files(file.path(MODELS, set_name), pattern = "^chr.*\\.summary\\.tsv$", full.names = TRUE))
  d  <- rbindlist(lapply(files, fread, sep = "\t"), fill = TRUE)
  cv <- d[cv_r2 > 0.01 & cv_pval < 0.05]
  usable <- cv[n_model > 0]
  data.table(set = set_name, attempted = nrow(d), cleared_cv = nrow(cv), usable = nrow(usable),
             empty_final = nrow(cv) - nrow(usable), median_r2 = median(usable$cv_r2))
}

load_dec <- function(path) {
  d <- fread(path, sep = "\t")
  d[degenerate == FALSE | degenerate == "False"]
}

shares <- function(path) {
  d <- load_dec(path)
  t <- abs(d[, .(phi_w, phi_D, phi_R)])
  s <- rowSums(t)
  o <- as.data.table(t / s)
  setnames(o, c("s_w", "s_D", "s_R"))
  o[, gene := d$gene]
  o[, dom := max.col(t, ties.method = "first") - 1L]
  o[]
}

draws <- complete_draws()
cat(sprintf("%d complete European draws: %s\n\n", length(draws), paste(draws, collapse = ", ")))

cat(strrep("=", 78), "\n", sep = "")
cat("YIELD  (paper: the draws span X to Y usable, and A to B clear cross-validation)\n")
cat(strrep("=", 78), "\n", sep = "")
y <- rbindlist(lapply(c(draws, YRI), yield_row))
print(y)
eur <- y[grepl("^EUR87", set)]
yri <- y[set == YRI]
cat(sprintf("\n  usable span across draws: %s to %s (width %d)\n",
            format(min(eur$usable), big.mark = ","), format(max(eur$usable), big.mark = ","),
            max(eur$usable) - min(eur$usable)))
cat(sprintf("  Yoruba usable: %s  -> inside the span: %s\n",
            format(yri$usable, big.mark = ","),
            min(eur$usable) <= yri$usable && yri$usable <= max(eur$usable)))
cat(sprintf("  cleared cross-validation, all sets: %s to %s\n",
            format(min(y$cleared_cv), big.mark = ","), format(max(y$cleared_cv), big.mark = ",")))
frac <- y$empty_final / y$cleared_cv
cat(sprintf("  empty final fit as a share of cleared: %.0f%% to %.0f%%\n", 100 * min(frac), 100 * max(frac)))

cat("\n", strrep("=", 78), "\n", sep = "")
cat("NOISE FLOOR  (paper: the floor lies between X and Y)\n")
cat(strrep("=", 78), "\n", sep = "")
floors <- list()
if (length(draws) >= 2L) {
  cb <- combn(draws, 2L)
  for (i in seq_len(ncol(cb))) {
    a <- cb[1L, i]; b <- cb[2L, i]
    p <- file.path(DEC, sprintf("%s_vs_%s.tsv", a, b))
    if (!file.exists(p)) { cat(sprintf("  %s vs %s: MISSING, not yet computed\n", a, b)); next }
    d <- load_dec(p)
    floors[[length(floors) + 1L]] <- list(name = sprintf("%s_vs_%s", a, b), n = nrow(d),
                                          v = median(abs(d$phi_w)))
  }
}
for (f in floors) cat(sprintf("  %-28s n=%4d  median |phi_w| = %.3f\n", f$name, f$n, f$v))
fv <- vapply(floors, function(f) f$v, numeric(1L))
cat(sprintf("\n  %d of %d possible pairs present\n", length(floors), choose(length(draws), 2L)))
cat(sprintf("  floor range %.3f to %.3f, median %.3f\n", min(fv), max(fv), median(fv)))

cat("\n", strrep("=", 78), "\n", sep = "")
cat("CROSS-ANCESTRY COMPONENTS  (paper: weights term lies between X and Y)\n")
cat(strrep("=", 78), "\n", sep = "")
cw <- numeric(0); sums <- numeric(0)
for (s in draws) {
  p <- file.path(DEC, sprintf("%s_vs_%s.tsv", s, YRI))
  if (!file.exists(p)) { cat(sprintf("  %s: MISSING\n", s)); next }
  d  <- load_dec(p)
  w  <- median(abs(d$phi_w)); dd <- median(abs(d$phi_D)); r <- median(abs(d$phi_R))
  cw <- c(cw, w); sums <- c(sums, dd + r)
  cat(sprintf("  %-12s n=%4d  |phi_w|=%.3f  |phi_D|=%.3f  |phi_R|=%.3f  D+R=%.3f\n",
              s, nrow(d), w, dd, r, dd + r))
}
cat(sprintf("\n  weights term across draws: %.3f to %.3f\n", min(cw), max(cw)))
cat(sprintf("  frequency plus LD sum:     %.3f to %.3f\n", min(sums), max(sums)))

cat("\n", strrep("=", 78), "\n", sep = "")
cat("THE MARGIN CLAIM  (paper: the floor exceeds frequency and LD combined,\n")
cat("                   whichever control pair and whichever draw is used)\n")
cat(strrep("=", 78), "\n", sep = "")
lo_floor <- min(fv); hi_sum <- max(sums)
cat(sprintf("  smallest floor      %.3f\n", lo_floor))
cat(sprintf("  largest D+R sum     %.3f\n", hi_sum))
cat(sprintf("  margin              %+.3f\n", lo_floor - hi_sum))
cat(sprintf("  VERDICT: %s\n", if (lo_floor > hi_sum) "PASS, the claim holds for every pairing"
                               else "FAIL, the sentence must be rewritten"))

cat("\n", strrep("=", 78), "\n", sep = "")
cat("DRAW DEPENDENCE  (paper: pooling all pairs gives N comparisons over G genes,\n")
cat("                  the largest component changes for X percent)\n")
cat(strrep("=", 78), "\n", sep = "")
per <- list()
for (s in draws) {
  p <- file.path(DEC, sprintf("%s_vs_%s.tsv", s, YRI))
  if (file.exists(p)) per[[s]] <- shares(p)
}
nms <- names(per)
merged <- list(); pair_names <- list()
if (length(nms) >= 2L) {
  cb <- combn(nms, 2L)
  for (i in seq_len(ncol(cb))) {
    mm <- merge(per[[cb[1L, i]]], per[[cb[2L, i]]], by = "gene", suffixes = c("_1", "_2"))
    merged[[length(merged) + 1L]] <- mm
    pair_names[[length(pair_names) + 1L]] <- paste(cb[1L, i], "vs", cb[2L, i])
  }
}
m <- rbindlist(merged)
M <- matrix(0L, 3L, 3L)
for (i in 0:2) for (j in 0:2) M[i + 1L, j + 1L] <- sum(m$dom_1 == i & m$dom_2 == j)
agree <- mean(m$dom_1 == m$dom_2)
cat(sprintf("  %d pairs, %s pooled comparisons, %d distinct genes\n",
            length(merged), format(nrow(m), big.mark = ","), uniqueN(m$gene)))
cat(sprintf("  agreement %.1f%%, so the largest component changes for %.1f%%\n", 100 * agree, 100 * (1 - agree)))
cat(sprintf("  matrix diagonal: weights %d, frequency %d, LD %d; switched %d\n",
            M[1, 1], M[2, 2], M[3, 3], sum(M) - sum(diag(M))))
cat(sprintf("  weights-share Spearman %.2f, median absolute change %.2f\n",
            cor(m$s_w_1, m$s_w_2, method = "spearman"), median(abs(m$s_w_1 - m$s_w_2))))
cat("  per-pair agreement:\n")
for (i in seq_along(merged))
  cat(sprintf("    %s: %d genes, %.1f%%\n", pair_names[[i]], nrow(merged[[i]]),
              100 * mean(merged[[i]]$dom_1 == merged[[i]]$dom_2)))

cat("\n  median share of the total, per draw (Figure 13 panel c):\n")
for (s in nms) {
  d <- per[[s]]
  cat(sprintf("    %-12s weights %.3f  frequency %.3f  LD %.3f\n",
              s, median(d$s_w), median(d$s_D), median(d$s_R)))
}
sw <- vapply(per, function(d) median(d$s_w), numeric(1L))
sd_ <- vapply(per, function(d) median(d$s_D), numeric(1L))
sr <- vapply(per, function(d) median(d$s_R), numeric(1L))
cat(sprintf("    ranges: weights %.3f to %.3f, frequency %.3f to %.3f, LD %.3f to %.3f\n",
            min(sw), max(sw), min(sd_), max(sd_), min(sr), max(sr)))
