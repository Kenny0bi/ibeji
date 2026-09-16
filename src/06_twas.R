# Summary-statistic TWAS (S-PrediXcan form) for one ibeji model set, all five traits.
#
# Usage: Rscript src/06_twas.R <set> <cores>
#
#   z_g = sum_j w_j (sigma_j / sigma_g) z_j,   sigma_g^2 = w' Sigma w
# with Sigma and sigma_j from the EUR358 genotypes, because all five GWAS are European.
# SNPs absent from a GWAS are dropped from that trait's sum and from sigma_g.
# Coverage is recorded three ways: share of model SNPs used, share of |weight| used,
# and share of the model's European-reference variance retained.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "lib_geno.R"))
suppressPackageStartupMessages(library(parallel))

args <- commandArgs(trailingOnly = TRUE)
set <- args[1]
cores <- if (length(args) >= 2) as.integer(args[2]) else 4L
chrom_filter <- if (length(args) >= 3) as.integer(strsplit(args[3], ",")[[1]]) else NULL
root <- project_root()
odir <- file.path(root, "results", "twas", set)
dir.create(odir, showWarnings = FALSE, recursive = TRUE)

TRAITS <- c("ASD", "SCZ", "BIP", "MDD", "PTSD")

pvar <- load_pvar(root)
ref_E <- sample_subset(root, set_ids(root, "EUR358"))
m <- load_models(root, set)
model_snps <- unique(m$weights$varID)

gwas <- lapply(TRAITS, function(t) {
  d <- fread(file.path(root, "data", "processed", "gwas", paste0(t, ".tsv.gz")))
  d <- d[varID %in% model_snps]
  setNames(d$z, d$varID)
})
names(gwas) <- TRAITS

twas_gene <- function(g, pg) {
  w <- m$weights[gene == g]
  x <- read_dosage(pg, pvar, w$varID)
  w <- w[match(colnames(x), varID)]
  S <- cov(x)
  sd_j <- sqrt(diag(S))
  full_var <- as.numeric(crossprod(w$weight, S %*% w$weight))
  rbindlist(lapply(TRAITS, function(t) {
    z <- gwas[[t]][w$varID]
    use <- !is.na(z) & sd_j > 0
    if (!any(use)) {
      return(data.table(gene = g, trait = t, n_model = nrow(w), n_used = 0L, z = NA_real_, p = NA_real_,
                        frac_snps_used = 0, frac_absw_used = 0, frac_var_used = 0))
    }
    wu <- w$weight[use]
    Su <- S[use, use, drop = FALSE]
    sg <- sqrt(as.numeric(crossprod(wu, Su %*% wu)))
    zg <- if (sg > 0) sum(wu * sd_j[use] * z[use]) / sg else NA_real_
    data.table(gene = g, trait = t, n_model = nrow(w), n_used = sum(use),
               z = zg, p = if (is.na(zg)) NA_real_ else 2 * pnorm(-abs(zg)),
               frac_snps_used = mean(use), frac_absw_used = sum(abs(wu)) / sum(abs(w$weight)),
               frac_var_used = if (full_var > 0) sg^2 / full_var else NA_real_)
  }))
}

run_chunk <- function(gs) {
  pg <- NewPgen(geno_file(root, "pgen"), sample_subset = ref_E$idx)
  on.exit(ClosePgen(pg))
  rbindlist(lapply(gs, twas_gene, pg = pg))
}

genes <- m$summary[significant == TRUE, gene]
out_name <- "twas_all_traits.tsv"
if (!is.null(chrom_filter)) {
  annot <- fread(file.path(root, "data", "processed", "expr", "gene_annotation.tsv"))
  genes <- intersect(genes, annot[chr %in% chrom_filter, gene])
  out_name <- sprintf("twas_chr%s.tsv", paste(chrom_filter, collapse = "_"))
}
chunks <- split(genes, ceiling(seq_along(genes) / 40))
res <- mclapply(chunks, run_chunk, mc.cores = cores, mc.preschedule = FALSE)
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) stop(as.character(res[bad][[1]]))
out <- rbindlist(res)
out[, p_bonf := pmin(1, p * sum(!is.na(p))), by = trait]
out[, fdr := p.adjust(p, "BH"), by = trait]
fwrite(out, file.path(odir, out_name), sep = "\t")
print(out[, .(genes = sum(!is.na(z)), bonferroni_hits = sum(p_bonf < 0.05, na.rm = TRUE),
              fdr05_hits = sum(fdr < 0.05, na.rm = TRUE), median_frac_var_used = median(frac_var_used, na.rm = TRUE)), by = trait])
