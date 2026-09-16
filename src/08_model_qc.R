# Model QC across ibeji training sets.
#
# Reports, per set: genes attempted, genes with a fitted model, significant models
# (cv R^2 > 0.01 and p < 0.05), median and mean R^2 among significant models,
# and R^2 for positive-control genes with strong, well-known LCL eQTLs.
# Writes results/qc/model_yield.tsv and results/qc/positive_controls.tsv.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "lib_geno.R"))

root <- project_root()
qdir <- file.path(root, "results", "qc")
dir.create(qdir, showWarnings = FALSE, recursive = TRUE)

# Ensembl IDs (GRCh37, version-stripped) for strong LCL eQTL genes
controls <- c(ERAP2 = "ENSG00000164308", GSTM3 = "ENSG00000134202", CHI3L2 = "ENSG00000064886",
              HLA_DQA1 = "ENSG00000196735", PEX6 = "ENSG00000124587", ZNF880 = "ENSG00000221923")

# Only sets with all 22 chromosomes. A set that is still training already has a directory but
# few or no summary files, and load_models would then rbindlist an empty list into a
# data.table with no columns and die with "object 'cv_r2' not found". This step first ran on
# 2026-09-16 at 00:38, while EUR87_r3 was training, and failed twice for exactly that reason.
# Partial sets are skipped rather than reported, because QC on 3 of 22 chromosomes would be
# a misleading yield number rather than a missing one.
all_sets <- list.dirs(file.path(root, "results", "models"), full.names = FALSE, recursive = FALSE)
n_summ <- vapply(all_sets, function(s)
  length(list.files(file.path(root, "results", "models", s), "summary.tsv$")), integer(1))
sets <- all_sets[n_summ == 22L]
if (any(n_summ != 22L)) {
  message("skipping sets that are not yet complete: ",
          paste(sprintf("%s (%d/22)", all_sets[n_summ != 22L], n_summ[n_summ != 22L]), collapse = ", "))
}
stopifnot(length(sets) > 0)
yield <- list()
ctrl <- list()
for (s in sets) {
  m <- load_models(root, s)
  sm <- m$summary
  cv_sig <- sm[!is.na(cv_r2) & cv_r2 > 0.01 & cv_pval < 0.05]
  sig <- sm[significant == TRUE]   # cv-significant AND non-empty final model
  yield[[s]] <- data.table(set = s, chromosomes_done = length(list.files(file.path(root, "results", "models", s), "summary.tsv$")),
                           genes_attempted = nrow(sm), genes_fit = sum(!is.na(sm$cv_r2)),
                           cv_significant = nrow(cv_sig), cv_significant_empty_final = sum(cv_sig$n_model == 0),
                           usable_models = nrow(sig), median_r2_usable = median(sig$cv_r2), mean_r2_usable = mean(sig$cv_r2),
                           median_snps_in_model = median(sig$n_model))
  sm[, gene_base := sub("\\..*$", "", gene)]
  ctrl[[s]] <- data.table(set = s, control = names(controls),
                          cv_r2 = sm$cv_r2[match(controls, sm$gene_base)],
                          n_model = sm$n_model[match(controls, sm$gene_base)])
}
yield <- rbindlist(yield)
ctrl <- dcast(rbindlist(ctrl), control ~ set, value.var = "cv_r2")
fwrite(yield, file.path(qdir, "model_yield.tsv"), sep = "\t")
fwrite(ctrl, file.path(qdir, "positive_controls.tsv"), sep = "\t")
print(yield)
print(ctrl)
