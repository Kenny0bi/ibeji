# Prepare expression phenotypes for each ibeji training set.
#
# Within each set, independently (so EUR87 draws get exactly the treatment YRI87 gets):
#   - autosomal genes expressed (RPKM > 0.1) in at least 20% of the set
#   - rank-based inverse normal transform per gene
#   - expression PCs from the transformed matrix
#   - residualize on sex, 3 genotype PCs, and K expression PCs
# K = 15 for n >= 300, K = 5 otherwise (the same K for EUR87 and YRI87).
suppressPackageStartupMessages(library(data.table))

root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
raw <- file.path(root, "data", "raw")
sdir <- file.path(root, "data", "processed", "samples")
gdir <- file.path(root, "data", "processed", "geno")
out <- file.path(root, "data", "processed", "expr")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

expr <- fread(file.path(raw, "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz"))
expr <- expr[Chr %in% as.character(1:22)]
psam <- fread(file.path(raw, "phase3_corrected.psam"))
setnames(psam, "#IID", "IID")

read_pcs <- function(pop) {
  pc <- fread(file.path(gdir, paste0("pca_", pop, ".eigenvec")))
  setnames(pc, "#IID", "IID")
  pc[, .(IID, PC1, PC2, PC3)]
}
pcs <- rbind(read_pcs("EUR358"), read_pcs("YRI87"))

int_transform <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))

sets <- sub("\\.txt$", "", list.files(sdir, pattern = "^(EUR358|YRI87|EUR87_r[0-9])\\.txt$"))
log_rows <- list()

for (set in sets) {
  ids <- fread(file.path(sdir, paste0(set, ".txt")))[["#IID"]]
  n <- length(ids)
  y <- as.matrix(expr[, ..ids])
  rownames(y) <- expr$TargetID

  expressed <- rowMeans(y > 0.1) >= 0.2
  y <- y[expressed, , drop = FALSE]
  y_int <- t(apply(y, 1, int_transform))

  k <- if (n >= 300) 15L else 5L
  epcs <- prcomp(t(y_int), center = TRUE, scale. = FALSE)$x[, seq_len(k), drop = FALSE]

  cov <- data.table(IID = ids)
  cov <- merge(cov, psam[, .(IID, SEX)], by = "IID", sort = FALSE)
  cov <- merge(cov, pcs, by = "IID", sort = FALSE)
  stopifnot(identical(cov$IID, ids))
  cmat <- cbind(1, sex = cov$SEX, as.matrix(cov[, .(PC1, PC2, PC3)]), epcs)

  qr_c <- qr(cmat)
  resid <- t(apply(y_int, 1, function(g) qr.resid(qr_c, g)))
  colnames(resid) <- ids

  dt <- data.table(gene = rownames(resid), resid)
  fwrite(dt, file.path(out, paste0(set, ".residuals.tsv.gz")), sep = "\t")
  log_rows[[set]] <- data.table(set = set, n = n, genes_expressed = nrow(resid), expression_pcs = k)
  message(sprintf("%s: n=%d, genes=%d, K=%d", set, n, nrow(resid), k))
}

annot <- expr[, .(gene = TargetID, chr = as.integer(Chr), coord = as.integer(Coord))]
fwrite(annot, file.path(out, "gene_annotation.tsv"), sep = "\t")
fwrite(rbindlist(log_rows), file.path(out, "expression_prep_log.tsv"), sep = "\t")
