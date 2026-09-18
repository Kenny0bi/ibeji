## Harmonize the five GWAS to the ibeji genotype variant IDs.
## Every z ends up aligned to the ALT allele; palindromic SNPs are dropped, because
## strand cannot be resolved from alleles alone.

suppressPackageStartupMessages(library(data.table))

find_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f[1L])) else getwd()
  normalizePath(file.path(d, ".."), mustWork = FALSE)
}
ROOT <- find_root()
RAW  <- file.path(ROOT, "data", "raw")
OUT  <- file.path(ROOT, "data", "processed", "gwas")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

COMP <- c(A = "T", T = "A", C = "G", G = "C")

TRAITS <- list(
  ASD  = list(file = "iPSYCH-PGC_ASD_Nov2017.gz", chrom = "CHR", pos = "BP",
              ea = "A1", oa = "A2", effect = c("OR", "SE")),
  SCZ  = list(file = "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz", chrom = "CHROM",
              pos = "POS", ea = "A1", oa = "A2", effect = c("BETA", "SE")),
  BIP  = list(file = "pgc-bip2021-all.vcf.tsv.gz", chrom = "#CHROM", pos = "POS",
              ea = "A1", oa = "A2", effect = c("BETA", "SE")),
  MDD  = list(file = "pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz", chrom = "#CHROM", pos = "POS",
              ea = "EA", oa = "NEA", effect = c("BETA", "SE")),
  PTSD = list(file = "eur_ptsd_pcs_v4_aug3_2021.vcf.gz", chrom = "#CHROM", pos = "POS",
              ea = "A1", oa = "A2", effect = c("Z"))
)

load_reference <- function() {
  ## Read geuv445.pvar compactly: skip the '##' metadata block, keep chrom/pos/ref/alt.
  ## Variant IDs are chr:pos:ref:alt by construction, so they are rebuilt only for
  ## matched rows.
  path <- file.path(ROOT, "data", "processed", "geno", "geuv445.pvar")
  pvar <- fread(path, sep = "\t", skip = "#CHROM",
                select = c("#CHROM", "POS", "REF", "ALT"))
  setnames(pvar, c("chrom", "pos", "ref", "alt"))
  pvar[, chrom := as.integer(chrom)][, pos := as.integer(pos)]
  pvar[, ref := as.character(ref)][, alt := as.character(alt)]
  pvar[]
}

count_meta_lines <- function(path) {
  ## PGC VCF-style tsv files open with '##' metadata lines before the header.
  con <- gzfile(path, "rt"); on.exit(close(con))
  n <- 0L
  repeat {
    line <- readLines(con, n = 1L)
    if (!length(line) || !startsWith(line, "##")) break
    n <- n + 1L
  }
  n
}

harmonize <- function(trait, spec, ref) {
  path <- file.path(RAW, spec$file)
  cols <- c(spec$chrom, spec$pos, spec$ea, spec$oa, spec$effect)
  g <- fread(path, sep = "\t", skip = count_meta_lines(path), select = cols)
  setnames(g, c(spec$chrom, spec$pos, spec$ea, spec$oa), c("chrom", "pos", "ea", "oa"))
  n_raw <- nrow(g)

  g[, chrom := suppressWarnings(as.integer(sub("chr", "", as.character(chrom), fixed = TRUE)))]
  g <- g[!is.na(chrom) & chrom >= 1L & chrom <= 22L]
  g[, pos := as.integer(pos)]
  g[, ea := toupper(ea)][, oa := toupper(oa)]
  g <- g[ea %in% names(COMP) & oa %in% names(COMP)]

  if (identical(spec$effect, c("OR", "SE")))        g[, z := log(OR) / SE]
  else if (identical(spec$effect, c("BETA", "SE"))) g[, z := BETA / SE]
  else                                              g[, z := Z]
  g <- g[is.finite(z)]

  palindromic <- COMP[g$ea] == g$oa
  n_pal <- sum(palindromic, na.rm = TRUE)
  g <- g[!palindromic]

  m <- merge(g, ref, by = c("chrom", "pos"), allow.cartesian = TRUE)
  m[, varID := paste(chrom, pos, ref, alt, sep = ":")]

  same   <- m$ea == m$alt & m$oa == m$ref
  swap   <- m$ea == m$ref & m$oa == m$alt
  same_c <- COMP[m$ea] == m$alt & COMP[m$oa] == m$ref
  swap_c <- COMP[m$ea] == m$ref & COMP[m$oa] == m$alt
  same_c[is.na(same_c)] <- FALSE; swap_c[is.na(swap_c)] <- FALSE

  m[, z_alt := fifelse(same | same_c, z, fifelse(swap | swap_c, -z, NA_real_))]
  flipped_flag <- swap | swap_c
  strand_flag  <- same_c | swap_c
  keep <- is.finite(m$z_alt)
  m <- m[keep]; flipped_flag <- flipped_flag[keep]; strand_flag <- strand_flag[keep]

  ## drop any varID appearing more than once, keeping none of them
  dup <- m$varID %in% m$varID[duplicated(m$varID)]
  m <- m[!dup]; flipped_flag <- flipped_flag[!dup]; strand_flag <- strand_flag[!dup]

  fwrite(m[, .(varID, z = z_alt)], file.path(OUT, paste0(trait, ".tsv.gz")), sep = "\t")
  row <- data.table(trait = trait, rows_raw = n_raw, palindromic_dropped = n_pal,
                    matched = nrow(m), flipped = sum(flipped_flag),
                    strand_complement = sum(strand_flag),
                    reference_variants = nrow(ref), coverage = nrow(m) / nrow(ref))
  fwrite(row, file.path(OUT, paste0(trait, ".report.tsv")), sep = "\t")
  print(row)
  row
}

if (sys.nframe() == 0L) {
  reference <- load_reference()
  for (t in names(TRAITS)) {
    if (file.exists(file.path(OUT, paste0(t, ".report.tsv")))) {
      cat(t, ": already harmonized, skipping\n", sep = ""); next
    }
    harmonize(t, TRAITS[[t]], reference)
  }
  reports <- lapply(names(TRAITS), function(t) {
    f <- file.path(OUT, paste0(t, ".report.tsv"))
    if (file.exists(f)) fread(f, sep = "\t") else NULL
  })
  report <- rbindlist(Filter(Negate(is.null), reports))
  fwrite(report, file.path(OUT, "harmonization_report.tsv"), sep = "\t")
  print(report)
}
