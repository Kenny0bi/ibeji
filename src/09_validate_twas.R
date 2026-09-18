## Validate the ibeji TWAS implementation against the official S-PrediXcan.
## Builds the PredictDB file with the sqlite3 binary, since RSQLite is not installed.
##
## Usage: Rscript src/09_validate_twas.R <set> <chromosomes> <trait>

suppressPackageStartupMessages(library(data.table))

find_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f[1L])) else getwd()
  normalizePath(file.path(d, ".."), mustWork = FALSE)
}
ROOT <- find_root()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("usage: Rscript src/09_validate_twas.R <set> <chromosomes> <trait>", call. = FALSE)
set_name <- args[1L]; chroms <- args[2L]; trait <- args[3L]
tag  <- gsub(",", "_", chroms, fixed = TRUE)
vdir <- file.path(ROOT, "results", "validation", sprintf("%s_chr%s", set_name, tag))

db <- file.path(vdir, "model.db")
if (file.exists(db)) unlink(db)

weights <- fread(file.path(vdir, "weights.tsv"), sep = "\t")
extra   <- fread(file.path(vdir, "extra.tsv"),   sep = "\t")
w_csv <- file.path(vdir, ".weights_import.csv")
e_csv <- file.path(vdir, ".extra_import.csv")
fwrite(weights, w_csv)
fwrite(extra,   e_csv)

sql <- c(".mode csv",
         sprintf(".import --skip 1 '%s' weights", w_csv),
         sprintf(".import --skip 1 '%s' extra",   e_csv))
## sqlite3's .import infers no types, so declare the tables first and import into them.
create_w <- sprintf("CREATE TABLE weights (%s);",
                    paste(sprintf("%s %s", names(weights),
                                  ifelse(vapply(weights, is.numeric, logical(1)), "REAL", "TEXT")),
                          collapse = ", "))
create_e <- sprintf("CREATE TABLE extra (%s);",
                    paste(sprintf("%s %s", names(extra),
                                  ifelse(vapply(extra, is.numeric, logical(1)), "REAL", "TEXT")),
                          collapse = ", "))
script <- c(create_w, create_e, sql,
            "CREATE INDEX weights_rsid ON weights (rsid);",
            "CREATE INDEX weights_gene ON weights (gene);")
system2("sqlite3", db, input = script)
unlink(c(w_csv, e_csv))

model_snps <- unique(as.character(weights$rsid))
g <- fread(file.path(ROOT, "data", "processed", "gwas", paste0(trait, ".tsv.gz")), sep = "\t")
g <- g[varID %in% model_snps]
parts <- tstrsplit(g$varID, ":", fixed = TRUE)
g[, A2 := parts[[3L]]][, A1 := parts[[4L]]]
gwas_path <- file.path(vdir, sprintf("gwas_%s.tsv.gz", trait))
fwrite(g[, .(SNP = varID, A1, A2, Z = z)], gwas_path, sep = "\t")

out <- file.path(vdir, sprintf("spredixcan_%s.csv", trait))
status <- system2("python3",
  c(file.path(ROOT, "tools", "MetaXcan", "software", "SPrediXcan.py"),
    "--model_db_path", db, "--covariance", file.path(vdir, "covariance.txt.gz"),
    "--gwas_file", gwas_path, "--snp_column", "SNP", "--effect_allele_column", "A1",
    "--non_effect_allele_column", "A2", "--zscore_column", "Z", "--keep_non_rsid",
    "--output_file", out, "--throw"))
if (status != 0L) quit(status = status)

theirs <- fread(out)
ours   <- fread(file.path(ROOT, "results", "twas", set_name, sprintf("twas_chr%s.tsv", tag)), sep = "\t")
ours   <- ours[trait == ..trait]
m <- merge(theirs, ours, by = "gene", suffixes = c("_spredixcan", "_ours"))
m <- m[is.finite(zscore) & is.finite(z)]

r    <- cor(m$zscore, m$z)
diff <- abs(m$zscore - m$z)
report <- data.table(set = set_name, chromosomes = chroms, trait = trait,
                     genes_compared = nrow(m), pearson_r = r,
                     max_abs_diff = max(diff), median_abs_diff = median(diff),
                     n_snps_used_match = mean(m$n_snps_used == m$n_used))
fwrite(report, file.path(vdir, sprintf("validation_%s.tsv", trait)), sep = "\t")
print(report)

worst <- m[order(-diff)][1:min(5L, nrow(m))]
print(worst[, .(gene, zscore, z, n_snps_used, n_used, diff = abs(zscore - z))])
