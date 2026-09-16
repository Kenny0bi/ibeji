# Export one model set, restricted to chosen chromosomes, in PredictDB form so the
# official S-PrediXcan (MetaXcan) can be run on exactly the same models.
#
# Usage: Rscript src/09_export_predictdb.R <set> <chromosomes, e.g. 22>
# Writes results/validation/<set>_chr<..>/weights.tsv, extra.tsv, covariance.txt.gz
# Covariances are computed from the EUR358 genotypes, the same reference 06_twas.R uses.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "lib_geno.R"))

args <- commandArgs(trailingOnly = TRUE)
set <- args[1]
chroms <- as.integer(strsplit(args[2], ",")[[1]])
root <- project_root()
odir <- file.path(root, "results", "validation", sprintf("%s_chr%s", set, paste(chroms, collapse = "_")))
dir.create(odir, showWarnings = FALSE, recursive = TRUE)

annot <- fread(file.path(root, "data", "processed", "expr", "gene_annotation.tsv"))
m <- load_models(root, set)
genes <- intersect(m$summary[significant == TRUE, gene], annot[chr %in% chroms, gene])
w <- m$weights[gene %in% genes]
parts <- tstrsplit(w$varID, ":", fixed = TRUE)
w[, `:=`(ref_allele = parts[[3]], eff_allele = parts[[4]])]
fwrite(w[, .(rsid = varID, gene, weight, ref_allele, eff_allele)], file.path(odir, "weights.tsv"), sep = "\t")

s <- m$summary[gene %in% genes]
s[, qval := p.adjust(cv_pval, "BH")]
fwrite(s[, .(gene, genename = gene, `n.snps.in.model` = n_model, `pred.perf.R2` = cv_r2,
             `pred.perf.pval` = cv_pval, `pred.perf.qval` = qval)], file.path(odir, "extra.tsv"), sep = "\t")

pvar <- load_pvar(root)
ref_E <- sample_subset(root, set_ids(root, "EUR358"))
pg <- NewPgen(geno_file(root, "pgen"), sample_subset = ref_E$idx)
cov_rows <- vector("list", length(genes))
for (i in seq_along(genes)) {
  wg <- w[gene == genes[i]]
  x <- read_dosage(pg, pvar, wg$varID)
  S <- cov(x)
  ij <- which(upper.tri(S, diag = TRUE), arr.ind = TRUE)
  cov_rows[[i]] <- data.table(GENE = genes[i], RSID1 = colnames(S)[ij[, 1]], RSID2 = colnames(S)[ij[, 2]], VALUE = S[ij])
}
ClosePgen(pg)
fwrite(rbindlist(cov_rows), file.path(odir, "covariance.txt.gz"), sep = " ")
message(sprintf("exported %d genes, %d weights to %s", length(genes), nrow(w), odir))
