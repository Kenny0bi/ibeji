## Build the training sample lists for ibeji.
## Matches GEUVADIS expression samples to 1000 Genomes phase 3 genotypes, then writes
## PLINK 2 --keep files for the European and Yoruba sets plus five European draws.
##
## This will not reproduce the existing downsamples. R's sample() and pandas'
## .sample(random_state=rep) are different generators, so it draws five valid but
## different European sets and overwrites the lists the trained models depend on.
## ALL, EUR358 and YRI87 are deterministic and do match.

suppressPackageStartupMessages(library(data.table))

find_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f[1L])) else getwd()
  normalizePath(file.path(d, ".."), mustWork = FALSE)
}
ROOT <- find_root()
RAW  <- file.path(ROOT, "data", "raw")
OUT  <- file.path(ROOT, "data", "processed", "samples")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

EUR_POPS          <- c("CEU", "FIN", "GBR", "TSI")
N_DOWNSAMPLE_REPS <- 5L

## Read only the header of the expression matrix: the sample identifiers begin at
## column 5, after the four annotation columns.
expr_header  <- fread(file.path(RAW, "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz"),
                      sep = "\t", nrows = 0L)
expr_samples <- names(expr_header)[-seq_len(4L)]

psam <- fread(file.path(RAW, "phase3_corrected.psam"), sep = "\t")
setnames(psam, sub("^#", "", names(psam)))
setkey(psam, IID)

matched <- expr_samples[expr_samples %in% psam$IID]
missing <- expr_samples[!expr_samples %in% psam$IID]

pops <- psam[match(matched, IID), Population]
names(pops) <- matched
eur <- sort(names(pops)[pops %in% EUR_POPS])
yri <- sort(names(pops)[pops == "YRI"])

write_keep <- function(name, ids) {
  d <- data.table(IID = ids)
  setnames(d, "IID", "#IID")
  fwrite(d, file.path(OUT, paste0(name, ".txt")), sep = "\t")
}

write_keep("ALL",    sort(c(eur, yri)))
write_keep("EUR358", eur)
write_keep("YRI87",  yri)

rows <- list()
for (rep in seq_len(N_DOWNSAMPLE_REPS)) {
  set.seed(rep)
  draw <- sort(sample(eur, length(yri)))
  write_keep(sprintf("EUR87_r%d", rep), draw)
  comp <- as.list(table(pops[draw]))
  rows[[rep]] <- c(list(set = sprintf("EUR87_r%d", rep), n = length(draw)), comp)
}

summary_rows <- c(
  list(c(list(set = "EUR358", n = length(eur)), as.list(table(pops[eur])))),
  list(list(set = "YRI87", n = length(yri), YRI = length(yri))),
  rows)
summary <- rbindlist(summary_rows, fill = TRUE)
for (j in names(summary)) if (is.numeric(summary[[j]])) set(summary, which(is.na(summary[[j]])), j, 0)
fwrite(summary, file.path(OUT, "set_composition.tsv"), sep = "\t")

cat(sprintf("expression samples: %d; matched to 1KG phase 3: %d\n", length(expr_samples), length(matched)))
cat(sprintf("missing from phase 3 (%d): %s\n", length(missing), paste(missing, collapse = ", ")))
print(summary)
