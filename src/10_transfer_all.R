# Cross-population prediction accuracy for every usable model, not only genes usable in both.
#
# Usage: Rscript src/10_transfer_all.R <source_set> <target_pop: EUR358 or YRI87> <cores>
#   e.g. Rscript src/10_transfer_all.R EUR87_r1 YRI87 2
#
# For each usable model in the source set: predict expression in the target individuals from
# their genotypes and correlate with their adjusted expression (signed r^2). Genes not expressed
# in the target set are recorded as not tested. This is the comparable benchmark to published
# cross-population GEUVADIS results, which average over all models rather than a selected subset.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "lib_geno.R"))
suppressPackageStartupMessages(library(parallel))

args <- commandArgs(trailingOnly = TRUE)
src_set <- args[1]
target <- args[2]
cores <- if (length(args) >= 3) as.integer(args[3]) else 2L
root <- project_root()
odir <- file.path(root, "results", "transfer")
dir.create(odir, showWarnings = FALSE, recursive = TRUE)

pvar <- load_pvar(root)
ref_T <- sample_subset(root, set_ids(root, target))
resid_T <- load_residuals(root, target)[, ref_T$ids, drop = FALSE]
m <- load_models(root, src_set)
genes <- m$summary[significant == TRUE, gene]
w_all <- split(m$weights, by = "gene")

one <- function(g, pg) {
  w <- w_all[[g]]
  if (!(g %in% rownames(resid_T))) {
    return(data.table(gene = g, tested = FALSE, n_snps = nrow(w), n_snps_found = NA_integer_, r = NA_real_, r2_signed = NA_real_, p = NA_real_))
  }
  x <- read_dosage(pg, pvar, w$varID)
  w <- w[match(colnames(x), varID)]
  pred <- as.numeric(x %*% w$weight)
  if (sd(pred) == 0) {
    return(data.table(gene = g, tested = TRUE, n_snps = nrow(w_all[[g]]), n_snps_found = ncol(x), r = 0, r2_signed = 0, p = 1))
  }
  ct <- cor.test(pred, resid_T[g, ])
  r <- unname(ct$estimate)
  data.table(gene = g, tested = TRUE, n_snps = nrow(w_all[[g]]), n_snps_found = ncol(x), r = r, r2_signed = r^2 * sign(r), p = ct$p.value)
}

run_chunk <- function(gs) {
  pg <- NewPgen(geno_file(root, "pgen"), sample_subset = ref_T$idx)
  on.exit(ClosePgen(pg))
  rbindlist(lapply(gs, one, pg = pg))
}

res <- mclapply(split(genes, ceiling(seq_along(genes) / 50)), run_chunk, mc.cores = cores, mc.preschedule = FALSE)
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) stop(as.character(res[bad][[1]]))
out <- rbindlist(res)
src_r2 <- m$summary[, .(gene, source_cv_r2 = cv_r2)]
out <- merge(out, src_r2, by = "gene", all.x = TRUE)
fwrite(out, file.path(odir, sprintf("%s_into_%s.tsv", src_set, target)), sep = "\t")
t <- out[tested == TRUE]
message(sprintf("%s -> %s: %d usable models, %d tested; median signed r2 %.3f (mean %.3f); nominal p<0.05 in target: %.1f%%; median source cv R2 %.3f",
                src_set, target, nrow(out), nrow(t), median(t$r2_signed), mean(t$r2_signed), 100 * mean(t$p < 0.05), median(t$source_cv_r2)))
