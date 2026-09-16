# Decompose cross-ancestry disagreement between a European and a Yoruba GReX model.
#
# Usage: Rscript src/05_decompose.R <eur_set> <yri_set> <cores>
#   e.g. Rscript src/05_decompose.R EUR87_r1 YRI87 7
#
# For every gene with a significant model in both sets:
#   V(w, D, R) = w' D^{1/2} R D^{1/2} w, with D = diag(2p(1-p)) and R the LD correlation.
#   Reference populations are always the full EUR358 and YRI87 genotypes, so the
#   frequency and LD ingredients do not depend on which European draw trained the model.
#   Delta = log V(wY, DY, RY) - log V(wE, DE, RE) is split exactly into
#   phi_w + phi_D + phi_R by the 3-player Shapley value.
# Also records weight agreement on shared SNPs, cross-population prediction accuracy
# against observed expression, and frequency and LD divergence at the models' SNPs.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "lib_geno.R"))
suppressPackageStartupMessages(library(parallel))

args <- commandArgs(trailingOnly = TRUE)
eur_set <- args[1]
yri_set <- args[2]
cores <- if (length(args) >= 3) as.integer(args[3]) else 4L  # 4, not 7: seven workers each holding genotype buffers thrashed swap on a 16 GB machine
root <- project_root()
odir <- file.path(root, "results", "decomposition")
dir.create(odir, showWarnings = FALSE, recursive = TRUE)

pvar <- load_pvar(root)
# Noise-floor mode: when both sets are European draws, both "populations" are EUR358,
# so phi_D = phi_R = 0 and delta = phi_w measures weight disagreement from training noise
# and tag choice alone. (Transfer accuracy in this mode is partly in-sample, because the
# draws overlap EUR358; it is not interpreted.)
same_ancestry <- startsWith(yri_set, "EUR")
second_pop <- if (same_ancestry) "EUR358" else "YRI87"
ref_E <- sample_subset(root, set_ids(root, "EUR358"))
ref_Y <- sample_subset(root, set_ids(root, second_pop))

mE <- load_models(root, eur_set)
mY <- load_models(root, yri_set)
genes <- intersect(mE$summary[significant == TRUE, gene], mY$summary[significant == TRUE, gene])

resid_E <- load_residuals(root, "EUR358")[, ref_E$ids]
resid_Y <- load_residuals(root, second_pop)[, ref_Y$ids]

wE_all <- split(mE$weights, by = "gene")
wY_all <- split(mY$weights, by = "gene")

V <- function(w, d, R) {
  s <- sqrt(d) * w
  as.numeric(crossprod(s, R %*% s))
}

# v(S): log V with the ingredients in S taken from YRI and the rest from EUR.
shapley3 <- function(vals) {
  # vals is a named list over subsets of c("w","D","R") keyed "S", "Sw", "SD", "SR", "SwD", "SwR", "SDR", "SwDR"
  # (an "S" prefix, because R cannot look up a list element named "")
  players <- c("w", "D", "R")
  key <- function(s) paste0("S", paste0(players[players %in% s], collapse = ""))
  out <- setNames(numeric(3), players)
  for (p in players) {
    others <- setdiff(players, p)
    subsets <- list(character(0), others[1], others[2], others)
    wts <- c(2, 1, 1, 2) / 6
    out[p] <- sum(vapply(seq_along(subsets), function(i) {
      wts[i] * (vals[[key(c(subsets[[i]], p))]] - vals[[key(subsets[[i]])]])
    }, numeric(1)))
  }
  out
}

decompose_gene <- function(g, pgE, pgY) {
  wE <- wE_all[[g]]
  wY <- wY_all[[g]]
  snps <- union(wE$varID, wY$varID)
  xE <- read_dosage(pgE, pvar, snps)
  xY <- read_dosage(pgY, pvar, snps)
  snps <- colnames(xE)
  w_E <- setNames(numeric(length(snps)), snps); w_E[wE$varID] <- wE$weight
  w_Y <- setNames(numeric(length(snps)), snps); w_Y[wY$varID] <- wY$weight

  pE <- colMeans(xE) / 2; pY <- colMeans(xY) / 2
  dE <- 2 * pE * (1 - pE); dY <- 2 * pY * (1 - pY)
  RE <- safe_cor(xE); RY <- safe_cor(xY)

  pick <- function(s) list(w = if ("w" %in% s) w_Y else w_E,
                           d = if ("D" %in% s) dY else dE,
                           R = if ("R" %in% s) RY else RE)
  combos <- list("", "w", "D", "R", c("w", "D"), c("w", "R"), c("D", "R"), c("w", "D", "R"))
  keys <- c("S", "Sw", "SD", "SR", "SwD", "SwR", "SDR", "SwDR")
  rawV <- vapply(combos, function(s) { k <- pick(s); V(k$w, k$d, k$R) }, numeric(1))
  names(rawV) <- keys
  degenerate <- any(rawV <= 0)
  phi <- if (degenerate) c(w = NA_real_, D = NA_real_, R = NA_real_) else shapley3(as.list(log(rawV)))

  shared <- wE$varID[wE$varID %in% wY$varID]
  w_r <- if (length(shared) >= 3) suppressWarnings(cor(w_E[shared], w_Y[shared])) else NA_real_
  same_sign <- if (length(shared)) mean(sign(w_E[shared]) == sign(w_Y[shared])) else NA_real_

  acc <- function(x, w, y) {
    pr <- as.numeric(x %*% w)
    if (sd(pr) == 0) return(c(r2 = 0, p = 1))
    ct <- cor.test(pr, y)
    r <- unname(ct$estimate)
    # signed r^2: a model that predicts in the wrong direction counts as negative accuracy
    c(r2 = r^2 * sign(r), p = unname(ct$p.value))
  }
  # Each training set applies its own expression filter, so a gene usable in EUR87 and YRI87
  # can be absent from the EUR358 (or YRI87) residual matrix. Transfer accuracy is then NA;
  # the variance decomposition itself does not need expression and is still computed.
  na_acc <- c(r2 = NA_real_, p = NA_real_)
  eur_in_yri <- if (g %in% rownames(resid_Y)) acc(xY, w_E, resid_Y[g, ]) else na_acc
  yri_in_eur <- if (g %in% rownames(resid_E)) acc(xE, w_Y, resid_E[g, ]) else na_acc

  # Tag-invariant agreement: do the two models predict the same genetic value in the
  # same people? Elastic net picks one tag among correlated SNPs somewhat arbitrarily,
  # so models with no shared SNPs can still capture one signal. High prediction
  # correlation with low shared-SNP overlap means the weights differ only in tag choice.
  pred_cor <- function(x) {
    a <- as.numeric(x %*% w_E); b <- as.numeric(x %*% w_Y)
    if (sd(a) == 0 || sd(b) == 0) NA_real_ else cor(a, b)
  }
  agree_in_E <- pred_cor(xE)
  agree_in_Y <- pred_cor(xY)

  a <- (abs(w_E) + abs(w_Y)); a <- a / sum(a)
  pbar <- (pE + pY) / 2
  ok <- pbar > 0 & pbar < 1
  # Nei's Fst for two equally weighted populations: Var(p) / (pbar (1 - pbar)) with
  # Var(p) = (pE - pY)^2 / 4. The factor of 4 was missing until 2026-09-15, which made
  # every stored value exactly four times the true Fst (per-SNP term capped at 4, not 1).
  # 'a' is normalised on line above, so this is a weighted mean, not a sum of ratios.
  fst_w <- sum(a[ok] * (pE[ok] - pY[ok])^2 / (4 * pbar[ok] * (1 - pbar[ok])))
  aa <- tcrossprod(a); diag(aa) <- 0
  ld_div <- if (length(snps) >= 2 && sum(aa) > 0) sum(aa * abs(RE - RY)) / sum(aa) else NA_real_

  data.table(gene = g, n_snps_E = nrow(wE), n_snps_Y = nrow(wY), n_union = length(snps), n_shared = length(shared),
             V_E = rawV[["S"]], V_Y = rawV[["SwDR"]], delta = log(rawV[["SwDR"]]) - log(rawV[["S"]]),
             phi_w = phi[["w"]], phi_D = phi[["D"]], phi_R = phi[["R"]], degenerate = degenerate,
             weight_r_shared = w_r, same_sign_shared = same_sign,
             pred_agreement_in_EUR = agree_in_E, pred_agreement_in_YRI = agree_in_Y,
             V_E_empirical = var(as.numeric(xE %*% w_E)), V_Y_empirical = var(as.numeric(xY %*% w_Y)),
             r2_eur_model_in_yri = eur_in_yri[["r2"]], p_eur_model_in_yri = eur_in_yri[["p"]],
             r2_yri_model_in_eur = yri_in_eur[["r2"]], p_yri_model_in_eur = yri_in_eur[["p"]],
             fst_weighted = fst_w, ld_divergence = ld_div,
             frac_weight_private_E = sum(abs(w_E)[dY == 0]) / sum(abs(w_E)),
             frac_weight_private_Y = sum(abs(w_Y)[dE == 0]) / sum(abs(w_Y)))
}

run_chunk <- function(gs) {
  pgE <- NewPgen(geno_file(root, "pgen"), sample_subset = ref_E$idx)
  pgY <- NewPgen(geno_file(root, "pgen"), sample_subset = ref_Y$idx)
  on.exit({ ClosePgen(pgE); ClosePgen(pgY) })
  rbindlist(lapply(gs, decompose_gene, pgE = pgE, pgY = pgY))
}

chunks <- split(genes, ceiling(seq_along(genes) / 40))
res <- mclapply(chunks, run_chunk, mc.cores = cores, mc.preschedule = FALSE)
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) stop(as.character(res[bad][[1]]))
out <- rbindlist(res)
out[, check_sum := delta - (phi_w + phi_D + phi_R)]
fwrite(out, file.path(odir, sprintf("%s_vs_%s.tsv", eur_set, yri_set)), sep = "\t")
message(sprintf("%s vs %s: %d genes decomposed (%d degenerate); max |Shapley sum error| = %.2e",
                eur_set, yri_set, nrow(out), sum(out$degenerate), max(abs(out$check_sum), na.rm = TRUE)))
