## Verify every reference in references.bib against Crossref, and every citation
## against the .bib. Exits non-zero on any failure, so it can gate the compile.
##
## Usage: Rscript paper/verify_refs.R [references.bib] [ibeji.tex]

suppressPackageStartupMessages({
  library(jsonlite)
  library(curl)
})

HERE  <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1L])) else getwd()
})
args  <- commandArgs(trailingOnly = TRUE)
BIB   <- if (length(args) >= 1L) args[1L] else file.path(HERE, "references.bib")
TEX   <- if (length(args) >= 2L) args[2L] else file.path(HERE, "ibeji.tex")
CACHE <- file.path(HERE, ".crossref_cache.json")
UA    <- "ibeji-reference-check (mailto:obidelek19@gmail.com)"

parse_bib <- function(text) {
  entries <- list()
  ## (?s) so "." crosses newlines: a BibTeX entry spans many lines, and without
  ## the flag this regex silently matches nothing and the file reads as empty.
  ms <- gregexpr("(?s)@(\\w+)\\{([^,]+),(.*?)\\n\\}", text, perl = TRUE)[[1L]]
  if (ms[1L] == -1L) return(entries)
  lens <- attr(ms, "match.length")
  for (i in seq_along(ms)) {
    blk  <- substr(text, ms[i], ms[i] + lens[i] - 1L)
    kind <- tolower(sub("(?s)^@(\\w+)\\{.*$", "\\1", blk, perl = TRUE))
    if (kind == "ieeetranbstctl") next   # bibliography style settings, not a reference
    key  <- trimws(sub("(?s)^@\\w+\\{([^,]+),.*$", "\\1", blk, perl = TRUE))
    fm   <- gregexpr("(\\w+)\\s*=\\s*\\{((?:[^{}]|\\{[^{}]*\\})*)\\}", blk, perl = TRUE)[[1L]]
    flen <- attr(fm, "match.length")
    fields <- list()
    if (fm[1L] != -1L) for (j in seq_along(fm)) {
      piece <- substr(blk, fm[j], fm[j] + flen[j] - 1L)
      nm    <- tolower(sub("^(\\w+)\\s*=.*", "\\1", piece, perl = TRUE))
      val   <- sub("^\\w+\\s*=\\s*\\{(.*)\\}$", "\\1", piece, perl = TRUE)
      fields[[nm]] <- trimws(gsub("\\s+", " ", val))
    }
    entries[[key]] <- fields
  }
  entries
}

norm_str <- function(s) {
  s <- tolower(gsub("[{}\\\\$]", "", if (is.null(s)) "" else s))
  s <- sub("^the\\s+", "", s)
  trimws(gsub("[^a-z0-9]+", " ", s))
}

similar <- function(a, b) {
  a <- norm_str(a); b <- norm_str(b)
  if (!nzchar(a) && !nzchar(b)) return(1)
  1 - adist(a, b)[1, 1] / max(nchar(a), nchar(b), 1L)
}

fetch_json <- function(url) {
  h <- new_handle(); handle_setheaders(h, "User-Agent" = UA)
  r <- tryCatch(curl_fetch_memory(url, handle = h), error = function(e) NULL)
  if (is.null(r) || r$status_code != 200L) return(NULL)
  tryCatch(fromJSON(rawToChar(r$content), simplifyVector = FALSE), error = function(e) NULL)
}

crossref <- function(doi, cache) {
  if (!is.null(cache[[doi]])) return(cache[[doi]])
  j <- fetch_json(paste0("https://api.crossref.org/works/", URLencode(doi, reserved = TRUE)))
  if (is.null(j)) return(list(error = "lookup failed"))
  m <- j$message
  dp <- NULL
  for (k in c("published-print", "published-online", "issued"))
    if (is.null(dp) && !is.null(m[[k]]$`date-parts`)) dp <- m[[k]]$`date-parts`[[1L]][[1L]]
  years <- unique(unlist(lapply(c("published-print", "published-online", "issued"),
                                function(k) if (!is.null(m[[k]]$`date-parts`)) m[[k]]$`date-parts`[[1L]][[1L]])))
  rec <- list(title = if (length(m$title)) m$title[[1L]] else "",
              journal = if (length(m$`container-title`)) m$`container-title`[[1L]] else "",
              first_author = if (length(m$author) && !is.null(m$author[[1L]]$family)) m$author[[1L]]$family else "",
              year = dp, years_all = sort(unlist(years)))
  cache[[doi]] <<- rec
  Sys.sleep(0.2)
  rec
}

main <- function() {
  entries <- parse_bib(paste(readLines(BIB, warn = FALSE), collapse = "\n"))
  tex     <- paste(readLines(TEX, warn = FALSE), collapse = "\n")
  ## Structural guard. On 2026-09-15 an edit that replaced the span from the
  ## availability section to \end{document} silently deleted these two lines; the
  ## paper still compiled, but with 65 undefined citations and no reference list.
  ## This turns that into a hard stop.
  structural  <- c("\\bibliographystyle{", "\\bibliography{")
  missing_cmd <- structural[!vapply(structural, function(c) grepl(c, tex, fixed = TRUE), logical(1))]
  if (length(missing_cmd)) {
    cat("FAIL  manuscript is missing: ", paste(missing_cmd, collapse = ", "), "\n", sep = "")
    quit(status = 1L)
  }
  groups <- regmatches(tex, gregexpr("\\\\cite\\{[^}]*\\}", tex, perl = TRUE))[[1L]]
  cited  <- unique(trimws(unlist(strsplit(sub("^\\\\cite\\{", "", sub("\\}$", "", groups)), ","))))
  cited  <- cited[nzchar(cited)]

  cache <<- if (file.exists(CACHE)) fromJSON(CACHE, simplifyVector = FALSE) else list()
  failures <- 0L
  for (key in sort(names(entries))) {
    f   <- entries[[key]]
    doi <- f$doi
    if (is.null(doi) || !nzchar(doi)) {
      cat(sprintf("FAIL  %s: no DOI\n", key)); failures <- failures + 1L; next
    }
    cr <- crossref(doi, cache)
    if (!is.null(cr$error)) {
      cat(sprintf("FAIL  %s: DOI %s did not resolve (%s)\n", key, doi, cr$error))
      failures <- failures + 1L; next
    }
    problems <- character(0)
    if (similar(f$title, cr$title) < 0.9)
      problems <- c(problems, sprintf("title differs: bib '%s' vs Crossref '%s'", f$title, cr$title))
    bib_first <- norm_str(sub(",.*", "", strsplit(if (is.null(f$author)) "" else f$author, "\\s+and\\s+")[[1L]][1L]))
    if (nzchar(cr$first_author) && nzchar(bib_first) && norm_str(cr$first_author) != bib_first)
      problems <- c(problems, sprintf("first author differs: bib '%s' vs Crossref '%s'", bib_first, cr$first_author))
    if (!identical(as.character(f$year), as.character(cr$year)) &&
        !(suppressWarnings(as.integer(f$year)) %in% cr$years_all))
      problems <- c(problems, sprintf("year differs: bib %s vs Crossref %s",
                                      f$year, paste(cr$years_all, collapse = ",")))
    if (!is.null(f$journal) && nzchar(f$journal) && nzchar(cr$journal) && similar(f$journal, cr$journal) < 0.8)
      problems <- c(problems, sprintf("journal differs: bib '%s' vs Crossref '%s'", f$journal, cr$journal))
    if (length(problems)) {
      failures <- failures + 1L
      cat(sprintf("FAIL  %s: %s\n", key, paste(problems, collapse = "; ")))
    } else {
      cat(sprintf("OK    %s: %s %s, %s\n", key, cr$first_author, cr$year, cr$journal))
    }
  }
  write(toJSON(cache, auto_unbox = TRUE, pretty = TRUE), CACHE)
  missing <- sort(setdiff(cited, names(entries)))
  uncited <- sort(setdiff(names(entries), cited))
  if (length(missing)) {
    failures <- failures + length(missing)
    cat("FAIL  cited but not in .bib: ", paste(missing, collapse = ", "), "\n", sep = "")
  }
  cat(sprintf("\n%d entries, %d cited keys; uncited entries (%d): %s\n",
              length(entries), length(cited), length(uncited),
              if (length(uncited)) paste(uncited, collapse = ", ") else "none"))
  cat(if (failures == 0L) "ALL REFERENCES VERIFIED\n" else sprintf("%d FAILURE(S)\n", failures))
  quit(status = if (failures) 1L else 0L)
}

if (sys.nframe() == 0L) main()
