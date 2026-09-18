## Export the real numbers the Manim animation shows, and prove they match the analysis.
## Refuses to write if any component disagrees with the decomposition by more than 1e-6.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

find_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f[1L])) else getwd()
  normalizePath(file.path(d, ".."), mustWork = FALSE)
}
ROOT <- find_root()
LOC  <- file.path(ROOT, "results", "figdata", "locus")
GENE <- "ENSG00000172404.4"
NAME <- "DNAJB7"

detail <- fread(file.path(LOC, sprintf("%s_snps_detail.tsv", NAME)), sep = "\t")
snps   <- as.character(detail$varID)

dosages <- function(pop) {
  raw <- fread(file.path(LOC, sprintf("%s_%s.raw", NAME, pop)), sep = "\t")
  first <- which(names(raw) == "PHENOTYPE") + 1L
  x <- raw[, first:ncol(raw), with = FALSE]
  setnames(x, sub("_[^_]*$", "", names(x)))
  as.matrix(x[, snps, with = FALSE])
}

ingredients <- function(x) {
  p <- colMeans(x) / 2
  d <- 2 * p * (1 - p)
  r <- suppressWarnings(cor(x))
  r[!is.finite(r)] <- 0
  diag(r) <- 1
  list(d = d, R = r)
}

eur <- ingredients(dosages("EUR358"))
yri <- ingredients(dosages("YRI87"))
wE  <- as.numeric(detail$w_EUR)
wY  <- as.numeric(detail$w_YRI)

logV <- function(w, d, R) {
  s <- sqrt(d) * w
  log(as.numeric(t(s) %*% R %*% s))
}

players <- c("w", "D", "R")
bits_grid <- expand.grid(w = 0:1, D = 0:1, R = 0:1)
vertex <- list()
for (i in seq_len(nrow(bits_grid))) {
  b <- as.integer(bits_grid[i, ])
  w <- if (b[1L]) wY else wE
  d <- if (b[2L]) yri$d else eur$d
  R <- if (b[3L]) yri$R else eur$R
  vertex[[paste(b, collapse = "")]] <- logV(w, d, R)
}

## Shapley value for each player: average marginal change over the 6 orders.
orders <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
phi   <- numeric(3)
paths <- list()
for (o in orders) {
  state <- c(0L, 0L, 0L)
  steps <- list()
  for (k in o) {
    before <- vertex[[paste(state, collapse = "")]]
    state[k] <- 1L
    after <- vertex[[paste(state, collapse = "")]]
    steps[[length(steps) + 1L]] <- list(player = players[k], change = after - before)
    phi[k] <- phi[k] + (after - before) / length(orders)
  }
  paths[[length(paths) + 1L]] <- list(order = players[o], steps = steps)
}

dec <- fread(file.path(ROOT, "results", "decomposition", "EUR87_r1_vs_YRI87.tsv"), sep = "\t")
row <- dec[gene == GENE][1L]
check <- c(phi_w = abs(phi[1L] - row$phi_w), phi_D = abs(phi[2L] - row$phi_D),
           phi_R = abs(phi[3L] - row$phi_R),
           delta = abs((vertex[["111"]] - vertex[["000"]]) - row$delta))
worst <- max(check)
if (worst > 1e-6)
  stop(sprintf("animation numbers disagree with the decomposition output: %s",
               paste(sprintf("%s=%.2e", names(check), check), collapse = ", ")), call. = FALSE)

out <- list(gene = NAME, eur_set = "EUR87_r1", yri_set = "YRI87", n_snps = length(snps),
            n_nonzero_eur = sum(wE != 0), n_nonzero_yri = sum(wY != 0),
            vertex_logV = vertex, delta = vertex[["111"]] - vertex[["000"]],
            phi = list(w = phi[1L], D = phi[2L], R = phi[3L]), paths = paths,
            weights_eur = wE, weights_yri = wY,
            max_abs_difference_vs_decomposition = worst)
write(toJSON(out, auto_unbox = TRUE, digits = NA, pretty = TRUE),
      file.path(ROOT, "animation", "dnajb7_numbers.json"))

cat(sprintf("8 vertices exported; delta %+.3f; phi w %+.3f D %+.3f R %+.3f; max difference vs decomposition %.2e\n",
            out$delta, phi[1L], phi[2L], phi[3L], worst))
for (p in paths)
  cat(paste(vapply(p$steps, function(s) sprintf("%s %+.3f", s$player, s$change), character(1L)),
            collapse = " -> "), "\n")
