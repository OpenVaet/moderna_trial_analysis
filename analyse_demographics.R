# Moderna demographics, recruitment and BLA Figure 1, v1.15.1 -- R 4.3+.
# Save as analyse_demographics.R in the project root.
# Open in RStudio, adjust SETTINGS, select all, and Run (or use Source).
# No options() or separate function call is required. The final line starts it.
# Diagram numbers are computed from ADSL subjects, not copied from the paper.
# The paper's reference counts are used only to report MATCH/MISMATCH.
# This reconstructs sponsor population flags; it does not independently derive
# eligibility from raw visits or establish that every exclusion was appropriate.
# Base/recommended R packages suffice for CSV input, PDF, SVG and PNG output.
# haven is needed only for XPT/SAS7BDAT input; pdftools only for P1 PDF input.
# P301 JMP support is included: no JMP software, Python or reference CSV required.
# All study inputs are discovered locally in prod_data; no network downloads occur.
# The map uses clinical_trial_sites.csv and screening_map_template.html when present.
# Missing map assets do not prevent Figure 1 or recruitment charts.
# First/last recruitment dates by phase are included in HTML and Recruitment_periods.csv.
# P1 PDF input needs pdftools (or an installed pdftotext command); P2 XPT needs haven.

# ---- SETTINGS: edit these values, then run the whole script ----
PHMPT_ROOT <- ""                   # Empty = detect project; or "D:/Projects/moderna_trial_analysis"
PHMPT_ADSL_FILE <- ""             # Empty = find local ADSL; or supply its path
PHMPT_OUTPUT_DIR <- ""            # Empty = a new dated folder under analysis/
PHMPT_EXCEL_FR <- TRUE             # Also write semicolon-delimited subject CSV
PHMPT_STRICT <- TRUE               # Stop after saving results if counts disagree
PHMPT_SHOW_DIAGRAM <- interactive() # Show the diagram in RStudio's Plots pane
PHMPT_PNG_DPI <- 180L
PHMPT_SITES_FILE <- "clinical_trial_sites.csv"       # In the project root
PHMPT_MAP_TEMPLATE <- "screening_map_template.html"  # In the project root
PHMPT_P1_FILE <- "" # Empty = discover the local DMID 20-0003 16.1.6 PDF
PHMPT_P2_FILE <- "" # Empty = discover local P201-d57 ADSL (original or converted)
# P1/P2/P3 mean clinical trial phases, not P301 Parts A/B/C.
# Set P1/P2 paths explicitly to choose other validated source snapshots.

# ---- FUNCTIONS ----

bla1_project_root <- function(root = "") {
  if (!is.null(root) && nzchar(root)) {
    root <- normalizePath(root,winslash="/",mustWork=TRUE)
    if(!dir.exists(root))stop("PHMPT_ROOT is not a directory: ",root)
    return(root)
  }
  paths <- getwd()
  for(frame in sys.frames()) if(!is.null(frame$ofile)) paths <- c(paths,dirname(frame$ofile))
  arg <- grep("^--file=",commandArgs(),value=TRUE)
  if(length(arg))paths <- c(paths,dirname(sub("^--file=","",arg[1])))
  if(interactive() && requireNamespace("rstudioapi",quietly=TRUE) && rstudioapi::isAvailable()) {
    document <- tryCatch(rstudioapi::getActiveDocumentContext()$path,error=function(e)"")
    if(nzchar(document))paths <- c(paths,dirname(document))
  }
  paths <- unique(normalizePath(paths[dir.exists(paths)],winslash="/",mustWork=TRUE))
  for(base in paths) {
    p <- base
    for(i in seq_len(8L)) {
      if(dir.exists(file.path(p,"prod_data")))return(p)
      parent <- dirname(p);if(parent==p)break;p<-parent
    }
  }
  stop("Cannot locate prod_data. Set PHMPT_ROOT to your project root and run extract_and_categorize.R first.",call.=FALSE)
}


# Read the fixed-width character and Float64/date columns used in the PHMPT
# P301 ADSL JMP. This is a deliberately limited reader, not a general JMP importer.
# Binary-layout references: JMPReader.jl and jmpio-python, revision
# 5e9d49f07f052e266c183f9e19a8ce087367324a (metadata.py / column.py).
# https://github.com/martyblaber/jmpio-python
# Unsupported encodings, truncated data and invalid text stop the read. No
# placeholder values are substituted. The original JMP is never modified.
bla1_read_jmp <- function(path, fields = bla1_fields()) {
  size <- file.info(path)$size
  if (is.na(size) || size < 64 || size > .Machine$integer.max)
    stop("Unsupported JMP file size: ", path)
  b <- readBin(path, "raw", n = size)
  take <- function(buf, at, n) {
    if (length(at) != 1L || length(n) != 1L || !is.finite(at) || !is.finite(n) ||
        at < 1 || n < 0 || at + n - 1 > length(buf)) stop("JMP byte range outside file/payload.")
    if (n == 0) raw() else buf[seq.int(at, length.out = n)]
  }
  uint <- function(buf, at, width) {
    z <- sum(as.integer(take(buf, at, width)) * 256^(seq_len(width) - 1L))
    if (!is.finite(z) || z > 2^53 - 1) stop("JMP integer exceeds exact R representation.")
    z
  }
  txt <- function(z) {
    nonzero <- which(z != as.raw(0))
    if (!length(nonzero)) return("")
    s <- rawToChar(z[seq_len(tail(nonzero, 1L))])
    if (!validUTF8(s)) stop("Invalid UTF-8 in JMP string.")
    Encoding(s) <- "UTF-8"; s
  }
  if (!identical(b[1:8], as.raw(c(255,255,0,0,7,0,0,0))))
    stop("Unsupported JMP signature. Supply an exported CSV/XPT/SAS7BDAT instead.")
  header <- grepRaw(as.raw(c(7,0,8,0,0,0)), b, fixed = TRUE)[1]
  if (is.na(header) || header < 33) stop("JMP dimensions not found.")
  nr <- uint(b, header - 32, 8); nc <- uint(b, header - 24, 4)
  if (nr < 1 || nr > .Machine$integer.max || nc < 1 || nc > 10000) stop("Invalid JMP dimensions.")
  offsets <- NULL
  # Column directory precedes column payloads; validate every pointer and name.
  starts <- grepRaw(as.raw(c(255,255)), b, all = TRUE, fixed = TRUE)
  starts <- starts[starts > 2]
  for (at in starts) {
    p <- at
    while (p <= length(b) && b[p] == as.raw(255)) p <- p + 1L
    p <- p + 10L
    if (p + 15 > length(b)) next
    if (uint(b, p, 4) + uint(b, p + 4, 4) != nc) next
    p <- p + 16 + 6 * nc + 28
    found <- tryCatch({
      while (identical(take(b, p, 2), as.raw(c(253,255))) ||
             identical(take(b, p, 2), as.raw(c(254,255))) ||
             identical(take(b, p, 2), as.raw(c(255,255)))) p <- p + 10 + uint(b, p + 2, 8)
      if (uint(b, p, 4) != nc) stop("Directory column count mismatch.")
      off <- vapply(seq_len(nc), function(i) uint(b, p + 4 + (i - 1) * 8, 8) + 1, numeric(1))
      if (any(diff(off) <= 0) || any(off <= p) || any(off + 8 > length(b))) stop("Invalid JMP column pointers.")
      nm <- vapply(off, function(o) txt(take(b, o + 2, uint(b, o, 2))), character(1))
      if (anyDuplicated(nm) || any(!nzchar(nm))) stop("Invalid JMP column names.")
      list(offsets = off, names = nm)
    }, error = function(e) NULL)
    if (!is.null(found)) { offsets <- found$offsets; column_names <- found$names; break }
  }
  if (is.null(offsets)) stop("JMP column directory not found or unsupported.")
  missing <- setdiff(fields, column_names)
  if (length(missing)) stop("JMP missing required columns: ", paste(missing, collapse = ", "))
  selected <- c(fields, intersect("COUNTRY", setdiff(column_names, fields)))
  values <- setNames(vector("list", length(selected)), selected)
  for (name in selected) {
    i <- match(name, column_names); p <- offsets[i] + 2 + uint(b, offsets[i], 2)
    dt <- as.integer(take(b, p, 6)); p <- p + 6
    end <- if (i < nc) offsets[i + 1] - 1 else length(b)
    if (!dt[1] %in% c(9L,10L)) stop("Unsupported JMP encoding in ", name, "; export this table to CSV.")
    h <- take(b, p, min(2048, end - p + 1))
    marker <- grepRaw(as.raw(c(239,190,254,202)), h, fixed = TRUE)[1]
    if (is.na(marker)) stop("JMP compression header missing in ", name)
    q <- p + marker - 1 + 4
    packed <- uint(b, q, 8); unpacked <- uint(b, q + 8, 8)
    if (unpacked > .Machine$integer.max || q + 16 + packed - 1 > end)
      stop("Invalid JMP compressed block length in ", name)
    payload <- withCallingHandlers(memDecompress(take(b, q + 16, packed), type = "gzip"),
      warning = function(w) stop("JMP decompression warning in ", name, ": ", conditionMessage(w)))
    if (length(payload) != unpacked) stop("JMP decompressed size mismatch in ", name)
    if (dt[1] == 9L && dt[2] %in% c(1L,2L) && all(dt[3:4] == 0L) && dt[5] > 0L) {
      width <- dt[5]
      if (length(payload) != nr * width) stop("Unsupported fixed-width JMP payload in ", name)
      values[[name]] <- vapply(seq_len(nr), function(j) txt(take(payload, (j - 1) * width + 1, width)), character(1))
    } else if (dt[1] == 10L && dt[6] == 8L && dt[4] == dt[5] && dt[4] %in% c(0L,3L,99L,114L)) {
      if (length(payload) != nr * 8) stop("JMP numeric row count mismatch in ", name)
      z <- readBin(payload, "double", n = nr, size = 8, endian = "little")
      if (any(is.infinite(z))) stop("Non-finite JMP number in ", name)
      missing <- is.na(z)
      if (dt[4] == 114L) {
        if (any(abs(z[!is.na(z)] / 86400 - round(z[!is.na(z)] / 86400)) > 1e-8))
          stop("JMP date includes a time component in ", name)
        z <- format(as.Date(z / 86400, origin = "1904-01-01"), "%Y-%m-%d")
      } else z <- as.character(z)
      z[missing | is.na(z)] <- ""; values[[name]] <- z
    } else stop("Unsupported JMP column type in ", name, ": ", paste(dt, collapse = ", "))
  }
  as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
}

bla1_local_candidates <- function(root, phase) {
  dirs <- file.path(root, "prod_data", c("csv_data", "xpt_data", "sas7bdat_data", "jmp_data", "pdf_data", "drops_data"))
  files <- unique(c(list.files(root, full.names = TRUE),
    unlist(lapply(dirs[dir.exists(dirs)], list.files, recursive = TRUE, full.names = TRUE))))
  files <- files[file.exists(files) & !dir.exists(files)]
  name <- basename(files)
  if (phase == "P1") {
    keep <- grepl("dmid[-_]200003[-_]16[-_]1[-_]6(__[0-9a-f]{32})?\\.pdf$|^P1_recruitment\\.csv$", name, ignore.case = TRUE)
  } else {
    keep <- grepl("(^|[_ -])adsl([_. -]|$)", name, ignore.case = TRUE) &
      grepl("\\.(csv|xpt|sas7bdat|jmp)$", name, ignore.case = TRUE)
    if (phase == "P2") keep <- keep & grepl("p201[-_]d57", name, ignore.case = TRUE)
    else keep <- keep & !grepl("p201|p101|dmid", name, ignore.case = TRUE)
  }
  files <- normalizePath(files[keep], winslash = "/", mustWork = TRUE)
  # Prefer original source formats to derivatives, and categorized copies to drops.
  priority <- if (phase == "P1") c("pdf","csv") else c("jmp","xpt","sas7bdat","csv")
  files[order(match(tolower(tools::file_ext(files)), priority), grepl("/drops_data/", files, fixed = TRUE), files)]
}

bla1_content_md5 <- function(x, phase) {
  fields <- if (phase == "P3") bla1_fields() else
    c("STUDYID","USUBJID","SITEID","COUNTRY","RECRUIT_DATE","DATE_BASIS","POPULATION","PHASE","ARM")
  x <- x[order(x$USUBJID), fields, drop = FALSE]
  x[] <- lapply(x, function(z) { z <- enc2utf8(as.character(z)); z[is.na(z)] <- ""; z })
  rownames(x) <- NULL
  tmp <- tempfile("adsl_content_"); on.exit(unlink(tmp))
  saveRDS(x, tmp, version = 2, compress = FALSE)
  unname(tools::md5sum(tmp))
}

bla1_resolve_input <- function(root, explicit = "", phase = "P3", output_dir) {
  chosen <- !is.null(explicit) && length(explicit) == 1L && nzchar(explicit)
  candidates <- if (chosen) bla1_path(root, explicit) else bla1_local_candidates(root, phase)
  audit <- data.frame(phase=character(),path=character(),format=character(),source_md5=character(),
    status=character(),rows=integer(),content_md5=character(),detail=character())
  audit_file <- file.path(output_dir, paste0("Input_discovery_", phase, ".csv"))
  reader <- switch(phase, P3 = function(p) bla1_classify(bla1_read(p)), P1 = bla1_p1, P2 = bla1_p2)
  hashes <- character(); contents <- character(); successful <- list()
  for (p in candidates) {
    message("Checking ", phase, " source: ", p)
    hash <- if (file.exists(p) && !dir.exists(p)) unname(tools::md5sum(p)) else ""
    if (is.na(hash)) hash <- ""
    same <- if (nzchar(hash)) match(hash, hashes) else NA_integer_
    if (!is.na(same)) {
      j <- same
      audit <- rbind(audit, data.frame(phase=phase,path=p,format=tolower(tools::file_ext(p)),source_md5=hash,
        status="IDENTICAL_COPY",rows=audit$rows[j],content_md5=audit$content_md5[j],detail=paste("Same bytes as",audit$path[j])))
      hashes <- c(hashes,hash); next
    }
    parsed <- tryCatch({
      if (!file.exists(p) || dir.exists(p)) stop("Input file not found.")
      x <- reader(p)
      if (!nrow(x)) stop("No eligible participant records.")
      list(data=x,fingerprint=bla1_content_md5(x,phase))
    },error=function(e)e)
    ok <- !inherits(parsed,"error")
    audit <- rbind(audit,data.frame(phase=phase,path=p,format=tolower(tools::file_ext(p)),source_md5=hash,
      status=if(ok)"VALID" else "REJECTED",rows=if(ok)nrow(parsed$data) else NA_integer_,
      content_md5=if(ok)parsed$fingerprint else "",detail=if(ok)"" else conditionMessage(parsed)))
    hashes <- c(hashes,hash)
    if(ok){successful[[length(successful)+1L]] <- list(path=p,data=parsed$data);contents <- c(contents,parsed$fingerprint)}
  }
  bla1_csv(audit,audit_file)
  if (!length(successful)) stop("No usable ",phase," source found under ",root,". ",
    if(phase=="P3") paste0("Expected the P301 BLA ADSL with CUTOFFDT=2021-03-26, including the PHMPT JMP file in prod_data/jmp_data or prod_data/drops_data. ",
      "The known release is md-production-050125.zip. ") else "",
    "See ",audit_file,if(nrow(audit))paste0("\n",paste(paste(basename(audit$path),audit$detail,sep=": "),collapse="\n")) else " Check PHMPT_ROOT and run extract_and_categorize.R first.",call.=FALSE)
  country_conflict <- FALSE
  if (phase == "P3") {
    country_sources <- successful[vapply(successful,function(s)"COUNTRY" %in% names(s$data),logical(1))]
    if (length(country_sources) > 1L) {
      countries <- lapply(country_sources,function(s)bla1_country(s$data$COUNTRY[order(s$data$USUBJID)]))
      country_conflict <- any(!vapply(countries,function(z)identical(z,countries[[1L]]),logical(1)))
    }
  }
  if (length(unique(contents)) > 1L || country_conflict) stop("Several non-equivalent ",phase," datasets passed validation. ",
    "No snapshots were combined. Set ",if(phase=="P3")"PHMPT_ADSL_FILE" else paste0("PHMPT_",phase,"_FILE"),
    " to the intended source; see ",audit_file,call.=FALSE)
  selected <- successful[[1L]]
  audit$status[audit$status=="VALID"] <- "EQUIVALENT_FORMAT"
  audit$status[audit$path==selected$path] <- "SELECTED"
  bla1_csv(audit,audit_file)
  message("Selected ",phase,": ",selected$path)
  selected
}


bla1_csv <- function(x, path, sep = ",") {
  con <- file(path, open = "wb")
  on.exit(close(con))
  writeBin(as.raw(c(239, 187, 191)), con) # UTF-8 BOM for Excel
  write.table(x, con, sep = sep, dec = ".", row.names = FALSE,
              col.names = TRUE, quote = TRUE, qmethod = "double", na = "")
  invisible(path)
}

bla1_fields <- function() strsplit(paste(
  "STUDYID USUBJID SUBJID SITEID AGE SEX RANDFL RANDNUM FASFL SAFFL",
  "MITTFL PPROTFL DOSE1FL DOSE2FL RTPCRBL ELECSBL SCOV2BL TRT01P TRT01A RANDDT",
  "TR01SDT DOSE2DT AP01SDT AP01EDT AP01SEDT AP02SDT UBLDDT PDVDT EUAVACDT",
  "EFFCODT CUTOFFDT FASEREAS SAFEREAS MITTERS PPEREAS EOSSTT EOSDT",
  "EOT01STT EOT01DT DCT01RS DCSREAS"), " ")[[1]]

bla1_read <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "jmp") {
    x <- bla1_read_jmp(path)
  } else if (ext == "csv") {
    x <- read.csv(path, colClasses = "character", na.strings = "NA",
                  check.names = FALSE, fileEncoding = "UTF-8-BOM")
  } else if (ext %in% c("xpt", "sas7bdat")) {
    if (!requireNamespace("haven", quietly = TRUE)) stop("Install haven for XPT/SAS input.")
    x <- as.data.frame(if (ext == "xpt") haven::read_xpt(path) else haven::read_sas(path))
    x[] <- lapply(x, function(z) {
      if (inherits(z, "Date")) format(z, "%Y-%m-%d")
      else if (inherits(z, "POSIXt")) format(z, "%Y-%m-%d", tz = "UTC")
      else as.character(z)
    })
  } else stop("Supported ADSL inputs: JMP, CSV, XPT and SAS7BDAT.")
  if (anyDuplicated(names(x))) stop("Duplicate input column names.")
  missing <- setdiff(bla1_fields(), names(x))
  if (length(missing)) stop("ADSL is missing required columns: ", paste(missing, collapse = ", "))
  if (!"source_row" %in% names(x)) x$source_row <- as.character(seq_len(nrow(x)))
  x <- x[c("source_row", bla1_fields(), intersect("COUNTRY", names(x)))]
  x[] <- lapply(x, function(z) { z[is.na(z)] <- ""; enc2utf8(z) })
  if (any(!validUTF8(unlist(x, use.names = FALSE)))) stop("Input has invalid UTF-8 text; export it as UTF-8.")
  nz <- x$AGE != ""
  age <- suppressWarnings(as.numeric(x$AGE[nz]))
  if (any(!is.finite(age))) stop("Invalid numeric AGE.")
  x$AGE[nz] <- as.character(age) # Normalize 45 vs 45.0 across file formats.
  x
}

bla1_classify <- function(x) {
  if (!nrow(x) || any(x$USUBJID == "") || anyDuplicated(x$USUBJID)) stop("Missing/duplicate ADSL subject keys.")
  if (any(x$STUDYID != "mRNA-1273-P301")) stop("Unexpected study.")
  if (any(x$CUTOFFDT != "2021-03-26")) stop("This reconstruction requires the 26 March 2021 BLA snapshot (CUTOFFDT).")
  if (any(!x$TRT01P %in% c("Placebo", "mRNA-1273"))) stop("Missing/unrecognized randomized arm.")
  flags <- c("RANDFL", "FASFL", "SAFFL", "MITTFL", "PPROTFL", "DOSE1FL", "DOSE2FL")
  for (n in flags) if (any(!x[[n]] %in% c("Y", "N"))) stop("Missing/invalid ", n)
  if (any(x$SAFFL == "Y" & !x$TRT01A %in% c("Placebo", "mRNA-1273"))) stop("Missing safety treatment.")
  if (any(x$FASFL == "Y" & x$RANDFL != "Y") ||
      any(x$SAFFL == "Y" & x$RANDFL != "Y") ||
      any(x$MITTFL == "Y" & x$FASFL != "Y") ||
      any(x$PPROTFL == "Y" & x$MITTFL != "Y")) stop("Inconsistent population nesting.")
  if (any(x$RANDFL == "Y" & x$FASFL == "N" & (x$DOSE1FL == "Y" | x$DOSE2FL == "Y")))
    stop("A randomized FAS exclusion received a dose: cannot label this branch 'no dose'.")
  x$IN_FIGURE1 <- x$RANDFL
  x$SAFETY_ARM <- ifelse(x$SAFFL == "Y", x$TRT01A, "")
  x$SAFETY_ARM_DIFFERS <- ifelse(x$SAFFL == "Y", ifelse(x$TRT01P != x$TRT01A, "Y", "N"), "")
  ex_mitt <- x$FASFL == "Y" & x$MITTFL == "N"
  x$MITT_EXCLUSION <- ""
  x$MITT_EXCLUSION[ex_mitt & x$SCOV2BL == "Positive"] <- "BASELINE_POSITIVE"
  x$MITT_EXCLUSION[ex_mitt & x$SCOV2BL == ""] <- "BASELINE_MISSING"
  x$MITT_EXCLUSION[ex_mitt & x$MITT_EXCLUSION == ""] <- "UNRESOLVED"
  ex_pp <- x$MITTFL == "Y" & x$PPROTFL == "N"
  phrases <- c(INCORRECT_VACCINE = "Received incorrect vaccination in blinded phase",
               DOSE2_OUTSIDE_WINDOW = "Received Vaccination #2 out of window",
               NO_SECOND_DOSE = "Did not receive Vaccination #2 per schedule",
               OTHER_MAJOR_DEVIATION = "Had major protocol deviations")
  indicators <- sapply(phrases, function(p) ex_pp & grepl(p, x$PPEREAS, fixed = TRUE))
  for (k in names(phrases)) x[[paste0("PP_HAS_", k)]] <- ifelse(indicators[, k], "Y", "N")
  x$PP_REASON_COUNT <- rowSums(indicators)
  x$PP_EXCLUSION <- ifelse(ex_pp, "UNRESOLVED", "")
  # Display rule inferred from the figure, not a submitted sponsor priority variable.
  # Preserve all overlapping reasons. The 11 observed overlaps are other-deviation
  # plus no-dose2 (10) or out-of-window (1); specific reasons take display priority.
  for (k in rev(names(phrases))) x$PP_EXCLUSION[indicators[, k]] <- k
  x$PP_CATEGORY_BASIS <- ifelse(ex_pp, ifelse(x$PP_REASON_COUNT > 1,
    "INFERRED_DISPLAY_PRIORITY; ALL_REASONS_RETAINED", "SINGLE_ADSL_REASON"), "")
  x$FIGURE1_FOOTNOTE_GROUP <- ""
  not_rand <- x$RANDFL == "N"
  duplicate_ids <- c("US3902125", "US3522027") # ADRG section 3.3, page 12
  x$FIGURE1_FOOTNOTE_GROUP[not_rand & x$SUBJID %in% duplicate_ids] <- "DUPLICATE_RANDOMIZATION"
  x$FIGURE1_FOOTNOTE_GROUP[not_rand & !x$SUBJID %in% duplicate_ids] <- "OTHER_ALL_SET_EXCLUSION"
  x$FIGURE1_BRANCH <- ifelse(not_rand, "EXCLUDED_BEFORE_FIGURE1",
    ifelse(x$FASFL == "N", "NO_DOSE",
    ifelse(x$MITTFL == "N", x$MITT_EXCLUSION,
    ifelse(x$PPROTFL == "Y", "PER_PROTOCOL", x$PP_EXCLUSION))))
  # Retain every source date; EFFCODT is subject-specific and is NOT CUTOFFDT.
  date_fields <- c("RANDDT", "TR01SDT", "DOSE2DT", "AP01SDT", "AP01EDT", "AP01SEDT",
                   "AP02SDT", "UBLDDT", "PDVDT", "EUAVACDT", "EFFCODT", "CUTOFFDT", "EOSDT", "EOT01DT")
  for (n in date_fields) {
    v <- x[[n]]; nz <- v != ""
    d <- suppressWarnings(as.Date(v[nz], "%Y-%m-%d"))
    if (any(is.na(d)) || any(format(d, "%Y-%m-%d") != v[nz])) stop("Invalid ISO date in ", n)
  }
  x$DOSE_INTERVAL_DAYS <- as.integer(as.Date(x$DOSE2DT) - as.Date(x$TR01SDT))
  x
}

bla1_nodes <- function(x) {
  masks <- list(); meta <- list()
  add <- function(id, label, arm, basis, expected, mask) {
    masks[[id]] <<- mask
    meta[[id]] <<- data.frame(node_id = id, figure1_label = label, arm = arm,
      arm_basis = basis, expected_n = expected, observed_n = sum(mask), stringsAsFactors = FALSE)
  }
  add("ORIGINAL", "Original randomization records (footnote)", "All", "ADSL records", 30423, rep(TRUE, nrow(x)))
  add("FIGURE1", "Randomization population shown in Figure 1", "All", "RANDFL", 30415, x$RANDFL == "Y")
  add("EXCLUDED_BEFORE", "Excluded before Figure 1 and all analysis sets", "All", "RANDFL", 8, x$RANDFL == "N")
  add("DUPLICATE_RANDOMIZATION", "Duplicate randomization exclusions (footnote)", "All", "ADRG identifiers", 2,
      x$FIGURE1_FOOTNOTE_GROUP == "DUPLICATE_RANDOMIZATION")
  add("OTHER_ALL_SET_EXCLUSION", "Other six all-set exclusions (major deviations in footnote)", "All", "RANDFL + ADRG identifiers", 6,
      x$FIGURE1_FOOTNOTE_GROUP == "OTHER_ALL_SET_EXCLUSION")
  exp_p <- c(15206, 40, 15166, 15162, 421, 337, 84, 14745, 581, 7, 119, 425, 30, 14164)
  exp_v <- c(15209, 29, 15180, 15184, 434, 347, 87, 14746, 459, 6, 102, 334, 17, 14287)
  keys <- c("RANDOMIZED", "NO_DOSE", "FAS", "SAFETY", "MITT_EXCLUDED", "BASELINE_POSITIVE",
            "BASELINE_MISSING", "MITT", "PP_EXCLUDED", "INCORRECT_VACCINE", "DOSE2_OUTSIDE_WINDOW",
            "NO_SECOND_DOSE", "OTHER_MAJOR_DEVIATION", "PP")
  labels <- c("Assigned treatment", "Did not receive a dose", "Full analysis population", "Safety population",
              "Excluded from modified intention-to-treat", "SARS-CoV-2 positive at baseline", "Missing baseline data",
              "Modified intention-to-treat population", "Excluded from per-protocol",
              "Received incorrect vaccine", "Received second dose outside planned time frame",
              "Discontinued without receiving second dose", "Other major protocol deviation", "Per-protocol population")
  for (a in c("Placebo", "mRNA-1273")) {
    arm <- x$TRT01P == a; prefix <- if (a == "Placebo") "P" else "V"
    expected <- if (a == "Placebo") exp_p else exp_v
    mm <- list(x$RANDFL == "Y" & arm,
      x$RANDFL == "Y" & x$FASFL == "N" & arm, x$FASFL == "Y" & arm,
      x$SAFFL == "Y" & x$TRT01A == a, x$FASFL == "Y" & x$MITTFL == "N" & arm,
      x$MITT_EXCLUSION == "BASELINE_POSITIVE" & arm, x$MITT_EXCLUSION == "BASELINE_MISSING" & arm,
      x$MITTFL == "Y" & arm, x$MITTFL == "Y" & x$PPROTFL == "N" & arm,
      x$PP_EXCLUSION == "INCORRECT_VACCINE" & arm, x$PP_EXCLUSION == "DOSE2_OUTSIDE_WINDOW" & arm,
      x$PP_EXCLUSION == "NO_SECOND_DOSE" & arm, x$PP_EXCLUSION == "OTHER_MAJOR_DEVIATION" & arm,
      x$PPROTFL == "Y" & arm)
    for (i in seq_along(keys)) add(paste0(prefix, "_", keys[i]), labels[i], a,
      if (keys[i] == "SAFETY") "TRT01A (actual blinded treatment)" else "TRT01P (randomized treatment)", expected[i], mm[[i]])
  }
  add("SWITCH_P_TO_V", "Placebo assigned, vaccine safety arm", "Placebo -> mRNA-1273", "TRT01P -> TRT01A", 7,
      x$SAFFL == "Y" & x$TRT01P == "Placebo" & x$TRT01A == "mRNA-1273")
  add("SWITCH_V_TO_P", "Vaccine assigned, placebo safety arm", "mRNA-1273 -> Placebo", "TRT01P -> TRT01A", 3,
      x$SAFFL == "Y" & x$TRT01P == "mRNA-1273" & x$TRT01A == "Placebo")
  counts <- do.call(rbind, meta); rownames(counts) <- NULL
  counts$difference <- counts$observed_n - counts$expected_n
  counts$status <- ifelse(counts$difference == 0, "MATCH", "MISMATCH")
  list(counts = counts, masks = masks)
}

bla1_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub('"', "&quot;", x, fixed = TRUE)
}

bla1_diagram <- function(counts) {
  # Diagram construction deliberately has no dependency on expected_n. Every
  # population/reason number below is looked up in the computed observed_n.
  if (anyDuplicated(counts$node_id)) stop("Duplicate count nodes.")
  n <- function(id) {
    j <- match(id, counts$node_id)
    if (anyNA(j)) stop("Missing diagram node: ", paste(id[is.na(j)], collapse = ", "))
    unname(counts$observed_n[j])
  }
  fmt <- function(z) format(z, big.mark = ",", scientific = FALSE, trim = TRUE)
  num <- function(id) fmt(n(id))
  boxes <- list()
  add <- function(id, ids, label, x, y, w, h, kind = "main", size = 18) {
    mismatch <- any(counts$status[match(ids, counts$node_id)] != "MATCH")
    fill <- switch(kind, main = "#edf4fb", safety = "#eaf5f1", exclusion = "#f5f6f8", root = "#173a5e")
    border <- if (mismatch) "#b45309" else if (kind == "root") "#173a5e" else "#53677b"
    if (mismatch) fill <- "#fff3df"
    text_color <- if (kind == "root" && !mismatch) "#ffffff" else "#172c42"
    lines <- strsplit(label, "\n", fixed = TRUE)[[1]]
    # Keep unexpected/unresolved lines visible without overflowing their box.
    size <- min(size, (h - 12) / (length(lines) * 1.18))
    boxes[[id]] <<- data.frame(node_id = id, source_node_ids = paste(ids, collapse = "|"),
      label = label, x = x, y = y, width = w, height = h, fill = fill,
      border = border, text_color = text_color, font_size = size,
      line_height = size * 1.18, stringsAsFactors = FALSE)
  }
  add("FIGURE1", "FIGURE1", paste0(num("FIGURE1"), " participants underwent\nrandomization"),
      530, 100, 340, 55, "root", 19)
  edges <- list()
  edge <- function(x, y) edges[[length(edges) + 1L]] <<- list(x = x, y = y)
  for (prefix in c("P", "V")) {
    left <- prefix == "P"; arm_name <- if (left) "placebo" else "mRNA-1273"
    id <- function(k) paste0(prefix, "_", k)
    nn <- function(k) num(id(k))
    mx <- if (left) 335 else 715; cx <- mx + 175
    sx <- if (left) 35 else 1065
    safety_x <- if (left) 35 else 1095
    add(id("RANDOMIZED"), id("RANDOMIZED"), paste0(nn("RANDOMIZED"), " assigned to receive\n", arm_name),
        mx, 195, 350, 60)
    add(id("NO_DOSE"), id("NO_DOSE"), paste0(nn("NO_DOSE"), " did not receive ", arm_name),
        sx, 285, 300, 45, "exclusion", 17)
    add(id("FAS"), id("FAS"), paste0(nn("FAS"), " received at least one dose\nand were included in the full\nanalysis population"),
        mx, 365, 350, 78)
    add(id("SAFETY"), id("SAFETY"), paste0(nn("SAFETY"), " received at least one\ndose and were included\nin the safety population"),
        safety_x, 365, 270, 78, "safety", 16.5)
    unresolved_mitt <- n(id("MITT_EXCLUDED")) - n(id("BASELINE_POSITIVE")) - n(id("BASELINE_MISSING"))
    mitt_text <- paste0(nn("MITT_EXCLUDED"), " excluded from the modified\nintention-to-treat population\n",
      nn("BASELINE_POSITIVE"), " SARS-CoV-2 positive\nat baseline\n", nn("BASELINE_MISSING"), " missing baseline data",
      if (unresolved_mitt > 0) paste0("\n", fmt(unresolved_mitt), " unresolved exclusion reason") else "")
    add(id("MITT_EXCLUDED"), id(c("MITT_EXCLUDED", "BASELINE_POSITIVE", "BASELINE_MISSING")),
        mitt_text, sx, 475, 300, 145, "exclusion", 17)
    add(id("MITT"), id("MITT"), paste0(nn("MITT"), " included in the modified\nintention-to-treat population"),
        mx, 660, 350, 60)
    reasons <- c("INCORRECT_VACCINE", "DOSE2_OUTSIDE_WINDOW", "NO_SECOND_DOSE", "OTHER_MAJOR_DEVIATION")
    unresolved_pp <- n(id("PP_EXCLUDED")) - sum(n(id(reasons)))
    pp_text <- paste0(nn("PP_EXCLUDED"), " excluded from the\nper-protocol population\n",
      nn("INCORRECT_VACCINE"), " received incorrect vaccine\n",
      nn("DOSE2_OUTSIDE_WINDOW"), " second dose outside\nplanned time frame\n",
      nn("NO_SECOND_DOSE"), " did not receive second dose\nper schedule\n",
      nn("OTHER_MAJOR_DEVIATION"), " other major protocol deviations",
      if (unresolved_pp > 0) paste0("\n", fmt(unresolved_pp), " unresolved exclusion reason") else "")
    add(id("PP_EXCLUDED"), id(c("PP_EXCLUDED", reasons)), pp_text, sx, 755, 300, 200, "exclusion", 17)
    add(id("PP"), id("PP"), paste0(nn("PP"), " included in the\nper-protocol population"),
        mx, 995, 350, 60)
    edge(c(700, 700, cx, cx), c(155, 175, 175, 195))
    edge(c(cx, cx), c(255, 365))
    edge(c(cx, if (left) sx + 300 else sx), c(307.5, 307.5))
    edge(c(if (left) mx else mx + 350, if (left) safety_x + 270 else safety_x), c(404, 404))
    edge(c(cx, cx), c(443, 660))
    edge(c(cx, if (left) sx + 300 else sx), c(547.5, 547.5))
    edge(c(cx, cx), c(720, 995))
    edge(c(cx, if (left) sx + 300 else sx), c(855, 855))
  }
  good <- sum(counts$status == "MATCH")
  footnotes <- c(
    paste0(num("ORIGINAL"), " source records; ", num("EXCLUDED_BEFORE"),
      " excluded before Figure 1/all analysis sets: ", num("DUPLICATE_RANDOMIZATION"),
      " duplicate-randomization records and ", num("OTHER_ALL_SET_EXCLUSION"), " other exclusions."),
    paste0("Safety arm: ", num("SWITCH_P_TO_V"), " placebo-assigned participants counted with mRNA-1273; ",
      num("SWITCH_V_TO_P"), " vaccine-assigned participants counted with placebo."),
    "Counts use submitted ADSL population flags. PP reason priority is inferred where reasons overlap; all raw reasons are retained.",
    "The no-second-dose branch uses the ADSL reason 'Did not receive Vaccination #2 per schedule'; it does not establish a withdrawal date.",
    "Reference: El Sahly et al. (2021), DOI 10.1056/NEJMoa2113017, Figure 1. All displayed counts are calculated from the input subjects.")
  list(width = 1400, height = 1265, boxes = do.call(rbind, boxes), edges = edges,
    title = "Figure 1. Randomization and analysis populations",
    subtitle = "Moderna mRNA-1273-P301 | Production ADSL | Data cutoff: 26 March 2021",
    status = paste0(good, "/", nrow(counts), " published counts match"),
    status_color = if (good == nrow(counts)) "#176a50" else "#b45309",
    footnotes = footnotes, counts = counts)
}

bla1_draw_diagram <- function(d) {
  grid::grid.newpage()
  # One native unit represents 1/100 inch in the exported PDF/PNG.
  # Scale the entire figure together in a smaller RStudio plotting pane.
  scale <- min(grDevices::dev.size("in") / c(d$width / 100, d$height / 100))
  grid::pushViewport(grid::viewport(width = grid::unit(d$width / 100 * scale, "inches"),
    height = grid::unit(d$height / 100 * scale, "inches"),
    xscale = c(0, d$width), yscale = c(d$height, 0)))
  on.exit(grid::popViewport())
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  txt <- function(label, x, y, size, color = "#172c42", just = "left", bold = FALSE) {
    grid::grid.text(label, x = grid::unit(x, "native"), y = grid::unit(y, "native"),
      just = c(just, "centre"), gp = grid::gpar(fontsize = size * 72 / 100 * scale,
        fontfamily = "sans", fontface = if (bold) "bold" else "plain", col = color))
  }
  txt(d$title, 35, 35, 25, bold = TRUE)
  txt(d$subtitle, 35, 65, 16, "#53677b")
  for (e in d$edges) grid::grid.lines(x = grid::unit(e$x, "native"), y = grid::unit(e$y, "native"),
    arrow = grid::arrow(length = grid::unit(0.085 * scale, "inches"), type = "closed"),
    gp = grid::gpar(col = "#53677b", fill = "#53677b", lwd = scale))
  for (i in seq_len(nrow(d$boxes))) {
    b <- d$boxes[i, ]; lines <- strsplit(b$label, "\n", fixed = TRUE)[[1]]
    # Sizes use positive viewport fractions because the native y scale runs
    # downward. A native height with that reversed scale would flip the box.
    grid::grid.rect(x = grid::unit(b$x + b$width / 2, "native"),
      y = grid::unit(b$y + b$height / 2, "native"),
      width = grid::unit(b$width / d$width, "npc"), height = grid::unit(b$height / d$height, "npc"),
      gp = grid::gpar(fill = b$fill, col = b$border, lwd = scale))
    ys <- b$y + b$height / 2 + (seq_along(lines) - (length(lines) + 1) / 2) * b$line_height
    for (j in seq_along(lines)) txt(lines[j], b$x + b$width / 2, ys[j], b$font_size,
                                    b$text_color, "centre")
  }
  grid::grid.lines(x = grid::unit(c(35, 1365), "native"), y = grid::unit(c(1090, 1090), "native"),
                   gp = grid::gpar(col = "#d7e0e8"))
  for (i in seq_along(d$footnotes)) txt(d$footnotes[i], 35, 1116 + (i - 1) * 23, 13.5, "#53677b")
  txt(d$status, 35, 1240, 15, d$status_color, "left", TRUE)
  invisible(d)
}

bla1_svg <- function(d) {
  esc <- bla1_escape
  text <- function(label, x, y, size, color = "#172c42", anchor = "start", bold = FALSE) {
    sprintf('<text x="%s" y="%s" font-size="%s" fill="%s" text-anchor="%s" dominant-baseline="central"%s>%s</text>',
      x, y, size, color, anchor, if (bold) ' font-weight="700"' else "", esc(label))
  }
  s <- c(sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %s %s" role="img" aria-labelledby="title desc">', d$width, d$height),
    paste0('<title id="title">', esc(d$title), '</title>'),
    '<desc id="desc">Population flow reconstructed from subject-level ADSL flags. Published values are used only for comparison.</desc>',
    '<rect width="100%" height="100%" fill="white"/>', '<g font-family="sans-serif">',
    text(d$title, 35, 35, 25, bold = TRUE), text(d$subtitle, 35, 65, 16, "#53677b"))
  for (e in d$edges) {
    s <- c(s, sprintf('<polyline points="%s" fill="none" stroke="#53677b" stroke-width="1.4"/>',
      paste(paste(e$x, e$y, sep = ","), collapse = " ")))
    # Explicit arrow polygons also work in SVG viewers without marker support.
    tip <- c(tail(e$x, 1), tail(e$y, 1))
    delta <- c(tail(diff(e$x), 1), tail(diff(e$y), 1))
    direction <- delta / sqrt(sum(delta^2))
    normal <- c(-direction[2], direction[1])
    points <- rbind(tip, tip - 9 * direction + 3.5 * normal, tip - 9 * direction - 3.5 * normal)
    s <- c(s, sprintf('<polygon points="%s" fill="#53677b"/>',
      paste(paste(points[, 1], points[, 2], sep = ","), collapse = " ")))
  }
  for (i in seq_len(nrow(d$boxes))) {
    b <- d$boxes[i, ]; lines <- strsplit(b$label, "\n", fixed = TRUE)[[1]]
    s <- c(s, sprintf('<g id="%s"><rect x="%s" y="%s" width="%s" height="%s" fill="%s" stroke="%s" stroke-width="1.4"/>',
      b$node_id, b$x, b$y, b$width, b$height, b$fill, b$border))
    ys <- b$y + b$height / 2 + (seq_along(lines) - (length(lines) + 1) / 2) * b$line_height
    for (j in seq_along(lines)) s <- c(s, text(lines[j], b$x + b$width / 2, ys[j], b$font_size, b$text_color, "middle"))
    s <- c(s, '</g>')
  }
  s <- c(s, '<line x1="35" x2="1365" y1="1090" y2="1090" stroke="#d7e0e8"/>')
  for (i in seq_along(d$footnotes)) s <- c(s, text(d$footnotes[i], 35, 1116 + (i - 1) * 23, 13.5, "#53677b"))
  s <- c(s, text(d$status, 35, 1240, 15, d$status_color, bold = TRUE))
  paste(c(s, '</g></svg>'), collapse = "\n")
}

bla1_export_diagram <- function(d, output_dir, png_dpi = 180L) {
  if (length(png_dpi) != 1L || !is.finite(png_dpi) || png_dpi < 72 || png_dpi > 600)
    stop("PHMPT_PNG_DPI must be between 72 and 600.")
  path <- function(ext) file.path(output_dir, paste0("BLA_Figure1_diagram.", ext))
  svg <- bla1_svg(d)
  writeLines(enc2utf8(svg), path("svg"), useBytes = TRUE)
  draw_file <- function(open_device) {
    open_device(); device <- grDevices::dev.cur()
    on.exit(grDevices::dev.off(device))
    bla1_draw_diagram(d)
  }
  draw_file(function() grDevices::pdf(path("pdf"), width = d$width / 100, height = d$height / 100,
                                      family = "Helvetica", useDingbats = FALSE, title = d$title))
  # PDF and SVG remain available if a platform cannot open a bitmap device.
  tryCatch(draw_file(function() grDevices::png(path("png"), width = d$width / 100,
      height = d$height / 100, units = "in", res = png_dpi,
      type = if (capabilities("cairo")) "cairo" else if (.Platform$OS.type == "windows") "windows" else "Xlib")),
    error = function(e) { if (file.exists(path("png"))) unlink(path("png"));
      warning("PNG unavailable; PDF/SVG were saved: ", conditionMessage(e), call. = FALSE) })
  bla1_csv(d$boxes, file.path(output_dir, "BLA_Figure1_diagram_nodes.csv"))
  cols <- c("node_id", "arm", "observed_n", "expected_n", "difference", "status")
  rows <- apply(d$counts[cols], 1, function(r) paste0("<tr>", paste0("<td>", bla1_escape(r), "</td>", collapse = ""), "</tr>"))
  html <- paste0('<!doctype html><html lang="en"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1">',
    '<title>Moderna demographics and recruitment</title><style>body{font:16px/1.5 Arial,sans-serif;color:#172c42;background:#f1f5f9;margin:0;padding:24px}',
    'main{max-width:1400px;margin:auto}nav,section{padding:20px;background:white;margin:18px 0;border-radius:8px}',
    'a{color:#175e9c}svg{display:block;width:100%;height:auto;background:white}table{border-collapse:collapse;width:100%}',
    'th,td{text-align:left;padding:8px;border-bottom:1px solid #d7e0e8}th{background:#edf4fb}.table{overflow-x:auto}</style></head><body><main>',
    '<nav><b>Subject-derived reconstruction</b> | <a href="BLA_Figure1_diagram.pdf">PDF</a> | ',
    '<a href="BLA_Figure1_diagram.svg">SVG</a> | <a href="BLA_Figure1_subjects.csv">Participants CSV</a> | ',
    '<a href="BLA_Figure1_counts.csv">Counts CSV</a> | <a href="Input_discovery_P3.csv">Input discovery</a><br>Counts come from ADSL records. Reference values are comparison targets only.</nav>',
    svg, '<!-- RECRUITMENT --><section><h2>Calculated counts versus the published figure</h2><p>',
    'A mismatch is retained and highlighted; it is never replaced by the published value.</p><div class="table"><table><thead><tr>',
    paste0('<th>', c("Node", "Arm", "Calculated", "Published", "Difference", "Status"), '</th>', collapse = ""),
    '</tr></thead><tbody>', paste(rows, collapse = "\n"), '</tbody></table></div></section></main></body></html>')
  writeLines(enc2utf8(html), file.path(output_dir, "demographics.html"), useBytes = TRUE)
  invisible(path("pdf"))
}

# ---- RECRUITMENT: P1 DMID 20-0003, P2 P201, P3 P301 ----
# Separate study cohorts, not P301 Parts A/B/C. P1 first-dose dates are a
# recruitment proxy; P2/P3 use randomization dates. Screen failures are excluded.

bla1_path <- function(root, path) {
  if (!grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)) path <- file.path(root, path)
  path
}

bla1_iso_date <- function(x, field) {
  x <- trimws(as.character(x)); x[is.na(x)] <- ""
  good <- x == "" | grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", x)
  d <- suppressWarnings(as.Date(x, "%Y-%m-%d"))
  if (any(!good | (x != "" & (is.na(d) | format(d, "%Y-%m-%d") != x)), na.rm = TRUE))
    stop(field, " must contain complete ISO dates (YYYY-MM-DD); no date origin is guessed.")
  d
}

bla1_table_read <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    x <- read.csv(path, colClasses = "character", na.strings = "NA",
                  check.names = FALSE, fileEncoding = "UTF-8-BOM")
  } else if (ext %in% c("xpt", "sas7bdat")) {
    if (!requireNamespace("haven", quietly = TRUE)) stop("Install haven to read ", basename(path))
    x <- as.data.frame(if (ext == "xpt") haven::read_xpt(path) else haven::read_sas(path))
    x[] <- lapply(x, function(z) if (inherits(z, "Date")) format(z, "%Y-%m-%d")
      else if (inherits(z, "POSIXt")) format(z, "%Y-%m-%d", tz = "UTC") else as.character(z))
  } else stop("Unsupported recruitment table: ", basename(path))
  if (anyDuplicated(names(x))) stop("Duplicate columns in ", basename(path))
  x[] <- lapply(x, function(z) { z[is.na(z)] <- ""; enc2utf8(z) })
  x
}

bla1_need <- function(x, fields) {
  missing <- setdiff(fields, names(x))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
}

bla1_country <- function(x) {
  x <- trimws(x); x[is.na(x) | x == ""] <- "Unknown / unmapped"
  x[toupper(x) %in% c("US", "USA", "UNITED STATES", "UNITED STATES OF AMERICA")] <- "United States"
  x
}

bla1_p1 <- function(path) {
  if (tolower(tools::file_ext(path)) == "pdf") {
    if (requireNamespace("pdftools", quietly = TRUE)) {
      pages <- pdftools::pdf_text(path)
    } else {
      command <- unname(Sys.which("pdftotext"))
      if (!nzchar(command)) stop('Install pdftools to read the local P1 PDF: install.packages("pdftools")')
      tmp <- tempfile(fileext=".txt"); on.exit(unlink(tmp),add=TRUE)
      diagnostic <- system2(command, args=vapply(c("-layout","-enc","UTF-8",normalizePath(path,winslash="/"),tmp),shQuote,""),stdout=TRUE,stderr=TRUE)
      status <- attr(diagnostic,"status")
      if ((!is.null(status) && status != 0L) || !file.exists(tmp)) stop("pdftotext failed: ",paste(diagnostic,collapse=" | "))
      pages <- strsplit(paste(readLines(tmp,warn=FALSE,encoding="UTF-8"),collapse="\n"),"\f",fixed=TRUE)[[1]]
    }
    if (!grepl("20-0003", paste(pages, collapse = " "), fixed = TRUE) ||
        !grepl("1st Treatment", paste(pages, collapse = " "), fixed = TRUE))
      stop("P1 PDF is not the expected DMID 20-0003 first-treatment listing.")
    rows <- list()
    months <- setNames(sprintf("%02d", 1:12), c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"))
    for (p in seq_along(pages)) for (line in strsplit(pages[p], "\n")[[1]]) {
      if (!grepl("^\\s*COV\\.[0-9]+\\b", line, perl = TRUE)) next
      id <- regmatches(line, regexpr("COV\\.[0-9]+", line))
      dates <- regmatches(line, gregexpr("[0-9]{2}[A-Z]{3}[0-9]{4}", line))[[1]]
      if (!length(dates)) stop("No first-dose date on P1 subject row: ", id)
      dt <- dates[1]; mm <- unname(months[substr(dt, 3, 5)])
      if (is.na(mm)) stop("Unrecognized English month: ", dt)
      rows[[length(rows) + 1L]] <- data.frame(STUDYID = "DMID 20-0003", USUBJID = id,
        SITEID = "", COUNTRY = "United States", RECRUIT_DATE = paste(substr(dt, 6, 9), mm, substr(dt, 1, 2), sep = "-"),
        DATE_BASIS = "First dose (recruitment proxy)", POPULATION = "First-dose recipients in original CSR", SOURCE_PAGE = as.character(p))
    }
    if (!length(rows)) stop("No P1 participants parsed; scanned PDFs require OCR first.")
    x <- do.call(rbind, rows)
  } else {
    x <- bla1_table_read(path)
    bla1_need(x, c("STUDYID", "USUBJID", "SITEID", "COUNTRY", "RECRUIT_DATE", "DATE_BASIS", "POPULATION", "SOURCE_PAGE"))
    if (any(x$STUDYID != "DMID 20-0003")) stop("Unexpected P1 STUDYID.")
  }
  if (!nrow(x) || any(!nzchar(x$USUBJID)) || anyDuplicated(x$USUBJID)) stop("P1 missing/duplicate subject keys.")
  x$RECRUIT_DATE <- bla1_iso_date(x$RECRUIT_DATE, "P1 RECRUIT_DATE")
  x$COUNTRY <- bla1_country(x$COUNTRY)
  x$PHASE <- "P1"; x$ARM <- ""
  x
}

bla1_p2 <- function(path) {
  x <- bla1_table_read(path)
  bla1_need(x, c("STUDYID", "USUBJID", "SITEID", "COUNTRY", "RANDFL", "RANDDT", "TRT01P"))
  if (any(x$STUDYID != "mRNA-1273-P201") || any(!x$RANDFL %in% c("Y", "N"))) stop("Invalid P201 study or randomization flags.")
  if (!nrow(x) || any(!nzchar(x$USUBJID)) || anyDuplicated(x$USUBJID)) stop("P2 missing/duplicate ADSL subject keys.")
  x <- x[x$RANDFL == "Y", , drop = FALSE]
  data.frame(STUDYID = x$STUDYID, USUBJID = x$USUBJID, SITEID = x$SITEID,
    COUNTRY = bla1_country(x$COUNTRY), RECRUIT_DATE = bla1_iso_date(x$RANDDT, "P2 RANDDT"),
    DATE_BASIS = "Randomization date", POPULATION = "Randomized P201 initial cohort (RANDFL=Y)",
    SOURCE_PAGE = "", PHASE = "P2", ARM = x$TRT01P)
}

bla1_plot_files <- function(draw, stem, output_dir, width, height, dpi) {
  device <- function(open) { open(); id <- grDevices::dev.cur(); on.exit(grDevices::dev.off(id)); draw() }
  device(function() grDevices::pdf(file.path(output_dir, paste0(stem, ".pdf")), width = width, height = height, useDingbats = FALSE))
  if (capabilities("cairo")) device(function() grDevices::svg(file.path(output_dir, paste0(stem, ".svg")), width = width, height = height))
  tryCatch(device(function() grDevices::png(file.path(output_dir, paste0(stem, ".png")), width = width,
    height = height, units = "in", res = dpi, type = if (capabilities("cairo")) "cairo" else if (.Platform$OS.type == "windows") "windows" else "Xlib")),
    error = function(e) warning("PNG unavailable for ", stem, ": ", conditionMessage(e), call. = FALSE))
}

bla1_js_string <- function(x) {
  # JSON-compatible string quoting, including HTML script-closing protection.
  vapply(enc2utf8(x), function(s) paste0('"', paste(vapply(utf8ToInt(s), function(n) {
    if (n == 34L) '\\"' else if (n == 92L) '\\\\'
    else if (n < 32L || n %in% c(60L, 62L, 38L, 8232L, 8233L)) sprintf("\\u%04x", n)
    else intToUtf8(n)
  }, character(1)), collapse = ""), '"'), character(1))
}

bla1_map <- function(site_counts, template_file, output_dir, total) {
  template <- paste(readLines(template_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  token <- "[---DATAPOINTS---]"
  pos <- gregexpr(token, template, fixed = TRUE)[[1]]
  if (pos[1] < 0 || length(pos) != 1L || !grepl("currentMouseOver", template, fixed = TRUE))
    stop("Map template must contain one [---DATAPOINTS---] token and currentMouseOver panel.")
  point <- function(i) {
    r <- site_counts[i, ]
    number <- function(v) if (is.na(v) || !is.finite(v)) "null" else format(v, scientific = FALSE, trim = TRUE, digits = 12)
    paste0('{"id":', bla1_js_string(r$SITEID), ',"name":', bla1_js_string(r$site_name),
      ',"country":', bla1_js_string(r$country), ',"lat":', number(r$latitude), ',"lon":', number(r$longitude),
      ',"n":', r$subjects_recruited, ',"vaccine":', r$mRNA_1273, ',"placebo":', r$Placebo, '}')
  }
  points <- paste(vapply(seq_len(nrow(site_counts)), point, character(1)), collapse = ",\n")
  js <- paste0('var recruitmentSites = [', points, '];\nvar recruitmentTotal = ', total, ';\n',
    paste(bla1_map_js(), collapse = "\n"))
  # Literal insertion avoids gsub replacement backreference interpretation.
  template <- paste0(substr(template, 1, pos - 1), js, substr(template, pos + nchar(token), nchar(template)))
  template <- sub("<title>.*?</title>", "<title>Moderna P3 recruitment by site</title>", template, perl = TRUE)
  template <- gsub("attributionControl: false", "attributionControl: true", template, fixed = TRUE)
  # The inherited example's unused draw plugin and jQuery callback are left
  # compatible; our own UI uses native DOM methods and safe text nodes.
  style <- paste0('<style>body{font:15px/1.45 Arial,sans-serif;color:#20354a;background:#f3f6fa}',
    'body>div:first-child{width:70%!important}body>div:nth-child(2){width:30%!important;max-height:100vh;overflow:auto}',
    '#currentMouseOver{padding:20px}h1{font-size:22px;margin:0 0 8px}h2{font-size:17px}',
    '.site-item{display:block;width:100%;text-align:left;background:white;border:1px solid #dce3ea;border-radius:5px;padding:9px;margin:6px 0;cursor:pointer}',
    '.site-item:focus,.site-item:hover{border-color:#2378aa;background:#edf5fb}input{box-sizing:border-box;width:100%;padding:9px;border:1px solid #99aabc}',
    '.muted{color:#52677a;font-size:13px}.site-detail{border-bottom:1px solid #dce3ea;padding:10px 0}',
    '@media(max-width:800px){body>div:first-child,body>div:nth-child(2){width:100%!important;max-height:none}#map{height:62vh!important}}</style>')
  template <- sub("</head>", paste0(style, "</head>"), template, fixed = TRUE)
  writeLines(enc2utf8(template), file.path(output_dir, "Recruitment_by_site_map.html"), useBytes = TRUE)
}

bla1_map_js <- function() c(
  '(function(){',
  'const host=document.getElementById("currentMouseOver"); host.replaceChildren();',
  'function el(tag,text,parent,cls){const x=document.createElement(tag); if(text!==null)x.textContent=text; if(cls)x.className=cls; if(parent)parent.appendChild(x); return x;}',
  'const fmt=n=>n.toLocaleString("en-US");',
  'el("h1","P3 recruitment by site",host);',
  'el("p",fmt(recruitmentTotal)+" randomized participants across "+recruitmentSites.filter(s=>s.n>0).length+" sites",host);',
  'el("p","P301 Figure 1 population (RANDFL=Y). Includes both randomized arms.",host,"muted");',
  'el("p","Bubble area represents participants at that location. Sites sharing coordinates appear together. Coordinates refer to sponsor-listed addresses and can represent campuses or offices.",host,"muted");',
  'el("p","Click a bubble or select a site below. Map tiles require an internet connection.",host,"muted");',
  'const details=el("div",null,host); el("h2","Select a site",details);',
  'const label=el("label","Find a site",host); label.htmlFor="site-search";',
  'const search=el("input",null,host); search.id="site-search"; search.placeholder="Site ID, name or country";',
  'const list=el("div",null,host); const groups=new Map(); const lookup=new Map();',
  'recruitmentSites.filter(s=>s.n>0&&Number.isFinite(s.lat)&&Number.isFinite(s.lon)).forEach(s=>{const k=s.lat+","+s.lon;if(!groups.has(k))groups.set(k,{lat:s.lat,lon:s.lon,n:0,sites:[]});const g=groups.get(k);g.n+=s.n;g.sites.push(s);});',
  'function content(g){const box=document.createElement("div");g.sites.forEach(s=>{const d=el("div",null,box,"site-detail");el("strong",s.id+" | "+s.name,d);el("div",s.country,d);el("div",fmt(s.n)+" randomized participants",d);el("div","mRNA-1273: "+fmt(s.vaccine)+"; placebo: "+fmt(s.placebo),d);});return box;}',
  'const maxN=Math.max(1,...Array.from(groups.values()).map(g=>g.n));',
  'const markers=[]; groups.forEach(g=>{const marker=L.circleMarker([g.lat,g.lon],{radius:28*Math.sqrt(g.n/maxN),color:"#174c78",weight:1.3,fillColor:"#3283b3",fillOpacity:0.66}).addTo(map);',
  'marker.bindTooltip(g.sites.map(s=>s.id).join(" / ")+": "+fmt(g.n));marker.bindPopup(content(g),{maxWidth:380});',
  'marker.on("mouseover",()=>{details.replaceChildren(content(g));});marker.on("click",()=>{details.replaceChildren(content(g));});g.marker=marker;g.sites.forEach(s=>lookup.set(s.id,g));markers.push(marker);});',
  'if(markers.length)map.fitBounds(L.featureGroup(markers).getBounds().pad(0.12));',
  'const mapped=Array.from(groups.values()).reduce((a,g)=>a+g.n,0); if(mapped!==recruitmentTotal)el("p",fmt(recruitmentTotal-mapped)+" participants have no usable map coordinates; retained in the CSV and plots.",host);',
  'function renderList(){const q=search.value.toLowerCase();list.replaceChildren();recruitmentSites.filter(s=>s.n>0&&(s.id+" "+s.name+" "+s.country).toLowerCase().includes(q)).sort((a,b)=>b.n-a.n||a.id.localeCompare(b.id)).forEach(s=>{const b=el("button",s.id+" | "+fmt(s.n)+" | "+s.name,list,"site-item");b.type="button";b.onclick=()=>{const g=lookup.get(s.id);if(g){map.setView([g.lat,g.lon],Math.max(map.getZoom(),8));g.marker.openPopup();details.replaceChildren(content(g));}else{details.replaceChildren();el("h2",s.id+" | "+s.name,details);el("p","Coordinates unavailable. "+fmt(s.n)+" participants.",details);}};});}',
  'search.addEventListener("input",renderList);renderList();',
  'const legend=L.control({position:"bottomleft"});legend.onAdd=function(){const box=L.DomUtil.create("div");box.style.cssText="background:white;padding:8px 12px;border:1px solid #cbd5e1;font:13px Arial";box.textContent="Bubble area proportional to recruitment; largest location: "+fmt(maxN);return box;};legend.addTo(map);',
  'map.attributionControl.addAttribution("Site coordinates include &copy; OpenStreetMap contributors (ODbL)");',
  '})();'
)

bla1_recruitment <- function(subjects, root, output_dir, adsl_file, sites_file, template_file, p1_file, p2_file, dpi) {
  sites_file <- bla1_path(root, sites_file); template_file <- bla1_path(root, template_file)
  has_sites <- file.exists(sites_file) && !dir.exists(sites_file)
  has_template <- file.exists(template_file) && !dir.exists(template_file)
  map_available <- has_sites && has_template
  sites <- if (has_sites) bla1_table_read(sites_file) else
    data.frame(SITEID=character(),site_name=character(),country=character(),latitude=character(),longitude=character())
  if (!map_available) message("Map assets incomplete: site counts and other charts will still be generated.")
  bla1_need(sites, c("SITEID", "site_name", "country", "latitude", "longitude"))
  if (any(!nzchar(sites$SITEID)) || anyDuplicated(sites$SITEID)) stop("Site asset has missing/duplicate SITEID values.")
  for (nm in c("latitude", "longitude")) sites[[nm]] <- suppressWarnings(as.numeric(sites[[nm]]))
  bad_geo <- !is.finite(sites$latitude) | !is.finite(sites$longitude) | abs(sites$latitude)>90 | abs(sites$longitude)>180
  sites$latitude[bad_geo] <- NA_real_; sites$longitude[bad_geo] <- NA_real_
  sites$country <- bla1_country(sites$country)
  cohorts <- list(); sources <- list()
  fields <- c("STUDYID","USUBJID","SITEID","COUNTRY","RECRUIT_DATE","DATE_BASIS","POPULATION","SOURCE_PAGE","PHASE","ARM")
  for (phase in c("P1", "P2")) {
    input <- ""
    result <- tryCatch({
      selected <- bla1_resolve_input(root, if (phase == "P1") p1_file else p2_file, phase, output_dir)
      input <- selected$path; z <- selected$data
      if (!nrow(z)) stop("Selected source contains no eligible participants.")
      z[fields]
    }, error = function(e) e)
    if (inherits(result, "error")) {
      sources[[phase]] <- data.frame(phase=phase, status="UNAVAILABLE", participants=NA_integer_, dated=NA_integer_, missing_dates=NA_integer_,
        date_basis="", population="", source_file=input, source_md5=if(nzchar(input) && file.exists(input)) unname(tools::md5sum(input)) else "", detail=conditionMessage(result))
      warning(phase, ": ", conditionMessage(result), call. = FALSE)
    } else {
      cohorts[[phase]] <- result
      sources[[phase]] <- data.frame(phase=phase, status="AVAILABLE", participants=nrow(result), dated=sum(!is.na(result$RECRUIT_DATE)), missing_dates=sum(is.na(result$RECRUIT_DATE)),
        date_basis=paste(unique(result$DATE_BASIS),collapse="; "), population=paste(unique(result$POPULATION),collapse="; "),
        source_file=normalizePath(input,winslash="/"), source_md5=unname(tools::md5sum(input)), detail="")
    }
  }
  p3_country <- if ("COUNTRY" %in% names(subjects)) subjects$COUNTRY else rep("",nrow(subjects))
  missing_country <- is.na(p3_country) | p3_country == ""
  p3_country[missing_country] <- sites$country[match(subjects$SITEID[missing_country],sites$SITEID)]
  cohorts$P3 <- data.frame(STUDYID=subjects$STUDYID, USUBJID=subjects$USUBJID, SITEID=subjects$SITEID,
    COUNTRY=bla1_country(p3_country), RECRUIT_DATE=bla1_iso_date(subjects$RANDDT,"P3 RANDDT"),
    DATE_BASIS="Randomization date", POPULATION="Figure 1 randomization population (RANDFL=Y)", SOURCE_PAGE="", PHASE="P3", ARM=subjects$TRT01P)
  sources$P3 <- data.frame(phase="P3",status="AVAILABLE",participants=nrow(subjects),dated=sum(!is.na(cohorts$P3$RECRUIT_DATE)),missing_dates=sum(is.na(cohorts$P3$RECRUIT_DATE)),
    date_basis="Randomization date",population="Figure 1 randomization population (RANDFL=Y)",source_file=adsl_file,source_md5=unname(tools::md5sum(adsl_file)),detail="P3 country uses ADSL COUNTRY when present, otherwise the site asset; missing country remains unknown.")
  source_table <- do.call(rbind,sources); rownames(source_table) <- NULL
  source_table$first_recruitment_date <- as.Date(rep(NA_character_, nrow(source_table)))
  source_table$last_recruitment_date <- as.Date(rep(NA_character_, nrow(source_table)))
  for (i in seq_len(nrow(source_table))) {
    cohort <- cohorts[[source_table$phase[i]]]
    if (is.null(cohort)) next
    dates <- cohort$RECRUIT_DATE[!is.na(cohort$RECRUIT_DATE)]
    if (length(dates)) {
      source_table$first_recruitment_date[i] <- min(dates)
      source_table$last_recruitment_date[i] <- max(dates)
    }
  }
  periods <- source_table[c("phase", "status", "participants", "dated", "missing_dates",
    "first_recruitment_date", "last_recruitment_date", "date_basis", "population", "source_file", "source_md5")]
  people <- do.call(rbind,cohorts); rownames(people) <- NULL
  if(anyDuplicated(paste(people$STUDYID,people$USUBJID,sep="|"))) stop("Duplicate recruitment study/subject keys.")
  people$WEEK_START <- people$RECRUIT_DATE - ((as.integer(format(people$RECRUIT_DATE,"%w"))+6L)%%7L)
  dated <- people[!is.na(people$WEEK_START),,drop=FALSE]
  if(!nrow(dated)) stop("No valid recruitment dates available.")
  weeks <- seq(min(dated$WEEK_START),max(dated$WEEK_START),by="7 days")
  weekly <- expand.grid(PHASE=c("P1","P2","P3"),WEEK_START=weeks,stringsAsFactors=FALSE)
  weekly$WEEK_START <- as.Date(weekly$WEEK_START,origin="1970-01-01")
  t <- table(dated$PHASE,as.character(dated$WEEK_START))
  weekly$subjects_recruited <- vapply(seq_len(nrow(weekly)),function(i){
    p<-weekly$PHASE[i]; w<-as.character(weekly$WEEK_START[i]); if(!p %in% names(cohorts)) return(NA_integer_)
    if(p %in% rownames(t) && w %in% colnames(t)) as.integer(t[p,w]) else 0L
  },integer(1))
  weekly$ISO_WEEK <- format(weekly$WEEK_START,"%G-W%V")
  weekly$DATE_BASIS <- source_table$date_basis[match(weekly$PHASE,source_table$phase)]
  country <- as.data.frame(table(COUNTRY=people$COUNTRY,PHASE=factor(people$PHASE,levels=c("P1","P2","P3"))),stringsAsFactors=FALSE)
  names(country)[3] <- "subjects_recruited"
  country$subjects_recruited[!country$PHASE %in% names(cohorts)] <- NA_integer_
  totals <- aggregate(rep(1L,nrow(people)),list(COUNTRY=people$COUNTRY),sum);names(totals)[2]<-"subjects_recruited"
  totals <- totals[order(-totals$subjects_recruited),,drop=FALSE]
  totals$scope <- paste(names(cohorts),collapse=" + ")
  unknown <- setdiff(unique(subjects$SITEID),sites$SITEID)
  if(length(unknown)) sites <- rbind(sites[c("SITEID","site_name","country","latitude","longitude")],
    data.frame(SITEID=unknown,site_name="Name unavailable",country=bla1_country(p3_country[match(unknown,subjects$SITEID)]),latitude=NA_real_,longitude=NA_real_))
  site_counts <- sites[c("SITEID","site_name","country","latitude","longitude")]
  site_counts$subjects_recruited <- as.integer(table(factor(subjects$SITEID,levels=site_counts$SITEID)))
  site_counts$mRNA_1273 <- as.integer(table(factor(subjects$SITEID[subjects$TRT01P=="mRNA-1273"],levels=site_counts$SITEID)))
  site_counts$Placebo <- as.integer(table(factor(subjects$SITEID[subjects$TRT01P=="Placebo"],levels=site_counts$SITEID)))
  site_counts$mapped <- is.finite(site_counts$latitude)&is.finite(site_counts$longitude)
  stopifnot(sum(site_counts$subjects_recruited)==nrow(subjects),sum(totals$subjects_recruited)==nrow(people))
  for(p in names(cohorts)) stopifnot(sum(weekly$subjects_recruited[weekly$PHASE==p])==sum(!is.na(cohorts[[p]]$RECRUIT_DATE)))
  for(pair in list(list(people,"Recruitment_subjects.csv"),list(source_table,"Recruitment_phase_sources.csv"),list(periods,"Recruitment_periods.csv"),
    list(weekly,"Recruitment_by_week.csv"),list(country,"Recruitment_by_country_phase.csv"),list(totals,"Recruitment_by_country.csv"),list(site_counts,"Recruitment_by_site.csv")))
    bla1_csv(pair[[1]],file.path(output_dir,pair[[2]]))
  asset_paths <- c(sites_file,template_file); present <- c(has_sites,has_template)
  asset_hashes <- rep("",2L); asset_hashes[present] <- unname(tools::md5sum(asset_paths[present]))
  bla1_csv(data.frame(asset=c("clinical_trial_sites.csv","screening_map_template.html"),path=asset_paths,status=ifelse(present,"AVAILABLE","MISSING"),md5=asset_hashes),file.path(output_dir,"Recruitment_assets.csv"))
  if (map_available) bla1_map(site_counts,template_file,output_dir,nrow(subjects))
  colors <- c(P1="#ad6532",P2="#23836b",P3="#236ba4")
  weekly_draw <- function(){
    old<-graphics::par(mfrow=c(3,1),mar=c(2.7,5.1,3.2,1.3),oma=c(2.1,0,2.5,0),family="sans");on.exit(graphics::par(old))
    for(p in names(colors)){
      w<-weekly[weekly$PHASE==p,,drop=FALSE];s<-source_table[source_table$phase==p,,drop=FALSE]
      ymax<-if(all(is.na(w$subjects_recruited))) 1 else max(1,w$subjects_recruited,na.rm=TRUE)
      xx<-as.numeric(w$WEEK_START)
      graphics::plot(xx,rep(0,length(xx)),type="n",ylim=c(0,ymax*1.17),xlim=c(as.numeric(as.Date(format(min(weeks),"%Y-%m-01")))-4,as.numeric(max(weeks))+4),xaxt="n",xlab="",ylab="Participants",bty="l",yaxs="i")
      graphics::abline(h=graphics::axTicks(2),col="#e5ebf0",lwd=0.8)
      if(s$status=="AVAILABLE") graphics::rect(xx-2.5,0,xx+2.5,w$subjects_recruited,col=colors[p],border=NA)
      else graphics::text(mean(range(xx)),0.5,"Data unavailable (not zero recruitment)",col="#755b44")
      ticks<-seq(as.Date(format(min(weeks),"%Y-%m-01")),max(weeks)+7,by="month")
      graphics::axis(1,at=as.numeric(ticks),labels=format(ticks,"%Y-%m"),cex.axis=0.88)
      graphics::title(main=paste0(p," | ",if(s$status=="AVAILABLE") paste(format(s$participants,big.mark=","),"participants") else "source unavailable"),adj=0,line=1.65,cex.main=1.05)
      graphics::mtext(if(s$status=="AVAILABLE") paste0(s$date_basis,"; ",s$dated," dated; ",s$missing_dates," missing dates") else "See Recruitment_phase_sources.csv",side=3,line=0.35,adj=0,cex=0.75,col="#52677a")
    }
    graphics::mtext("Weekly recruitment by clinical trial phase",outer=TRUE,side=3,line=0.7,adj=0,cex=1.2,font=2)
    graphics::mtext("Monday-start calendar weeks. Separate vertical scales. P1 uses first dose as a recruitment proxy.",outer=TRUE,side=1,line=0.7,cex=0.8)
  }
  bla1_plot_files(weekly_draw,"Recruitment_by_week",output_dir,12,9,dpi)
  country_draw <- function(){
    old<-graphics::par(mar=c(6,10,4.5,3.5),family="sans");on.exit(graphics::par(old))
    countries<-rev(totals$COUNTRY);mat<-matrix(0,nrow=3,ncol=length(countries),dimnames=list(names(colors),countries))
    for(i in seq_len(nrow(country))) if(!is.na(country$subjects_recruited[i])) mat[as.character(country$PHASE[i]),as.character(country$COUNTRY[i])]<-country$subjects_recruited[i]
    tt<-colSums(mat); mids<-graphics::barplot(mat,horiz=TRUE,col=colors,border=NA,las=1,ylim=c(-0.2,max(2,length(countries)*1.2+0.7)),xlim=c(0,max(tt)*1.23),xlab="Participants counted in the available study cohorts",cex.names=0.95)
    graphics::text(tt,mids,labels=format(tt,big.mark=","),pos=4,cex=1,font=2)
    graphics::title("Recruitment by country",adj=0,line=2.7)
    graphics::mtext("P1 first-dose recipients; P2/P3 randomized participants",side=3,line=1.2,adj=0,cex=0.85,col="#52677a")
    graphics::legend("topright",horiz=TRUE,legend=vapply(names(colors),function(p)paste0(p,": ",if(p %in% names(cohorts))format(nrow(cohorts[[p]]),big.mark=",") else "unavailable"),character(1)),fill=colors,bty="n",cex=0.85)
    graphics::mtext("Totals count study participation records; people are not linked across different studies.",side=1,line=4.6,cex=0.75,col="#52677a")
  }
  bla1_plot_files(country_draw,"Recruitment_by_country",output_dir,12,max(5,2.8+0.55*nrow(totals)),dpi)
  writeLines(c("RECRUITMENT METHODS", "",
    "P1 = DMID 20-0003; P2 = mRNA-1273-P201; P3 = mRNA-1273-P301. These are clinical phases, not P301 Parts A/B/C.",
    "P3 map population: RANDFL=Y in the same BLA ADSL used by Figure 1. Screen failures and the 8 pre-figure excluded records are not included.",
    "P2 population: RANDFL=Y in the selected P201 initial-cohort ADSL. Later booster/variant cohorts are not merged.",
    "P1 population: participants listed with a first treatment in the original CSR 16.1.6. First-dose dates are a proxy, not screening or randomization dates.",
    "The P1 synopsis states that study centers are in the United States; country is assigned from that study-level source. P1 site identifiers are unavailable in this listing.",
    "P1 source: 125752_S3_M5_5352_dmid-200003-16-1-6.pdf. Country/scope: dmid-200003-synopsis.pdf, pages 1 and 5.",
    "Known P2 PHMPT source: 125752_S3_M5_mrna-1273-p201-d57_A_D_adsl.xpt, md-production-050125.zip. Snapshot CUTOFFDT=2020-11-05.",
    "Missing phase sources are UNAVAILABLE, not zero. Missing dates remain in country/site totals, are omitted from weekly bars, and are counted in phase_sources.csv.",
    "Recruitment_periods.csv gives the earliest and latest non-missing recruitment dates for each available cohort (YYYY-MM-DD), also included in Recruitment_phase_sources.csv and the HTML summary.",
    "These are observed date ranges in the selected source cohorts, not necessarily the full enrollment periods of every later study extension. P1 uses first dose as a recruitment proxy; P2/P3 use randomization.",
    "Missing dates are not imputed. A phase with no dated records has blank endpoints; incomplete date coverage may understate its actual recruitment period.",
    "Weekly bins start Monday and are calculated from full calendar dates. Plots share a horizontal calendar axis but use separate vertical scales.",
    "Each study/USUBJID is counted once. Country totals count study participation records; cross-study people cannot be deduplicated with these identifiers.",
    "The map is P3 only because clinical_trial_sites.csv describes P301. It uses the supplied root HTML template, adds proportional-area markers and actual SITEID labels, and retains map attribution.",
    "Coincident coordinates are grouped into one bubble; all distinct sites and their individual counts remain accessible in the popup and searchable sidebar.",
    "Coordinate estimates may describe campus or administrative addresses; consult the clinical_trial_sites source bundle for details.",
    "HTML map libraries and tiles load online. CSV/PDF/SVG/PNG outputs do not require an internet connection to read.",
    "Source and asset hashes are saved in Recruitment_phase_sources.csv and Recruitment_assets.csv."),file.path(output_dir,"Recruitment_README.txt"),useBytes=TRUE)
  list(sources=source_table,periods=periods,weekly=weekly,country=country,totals=totals,sites=site_counts,map_available=map_available)
}

bla1_recruitment_html <- function(r,output_dir) {
  esc<-bla1_escape
  summary <- r$sources[c("phase","status","participants","dated","missing_dates","first_recruitment_date","last_recruitment_date","date_basis")]
  for (nm in c("first_recruitment_date", "last_recruitment_date")) {
    summary[[nm]] <- format(summary[[nm]], "%Y-%m-%d")
    missing <- is.na(summary[[nm]])
    summary[[nm]][missing] <- ifelse(summary$status[missing] == "UNAVAILABLE", "Unavailable", "No dated records")
  }
  rows<-apply(summary,1,function(z)paste0("<tr>",paste0("<td>",esc(z),"</td>",collapse=""),"</tr>"))
  picture<-function(stem){ext<-if(file.exists(file.path(output_dir,paste0(stem,".svg"))))"svg" else if(file.exists(file.path(output_dir,paste0(stem,".png"))))"png" else "";
    if(nzchar(ext))paste0('<img style="width:100%;height:auto" src="',stem,'.',ext,'" alt="',gsub("_"," ",stem),'"><p><a href="',stem,'.pdf">PDF</a> | <a href="',stem,'.csv">CSV</a></p>') else paste0('<p><a href="',stem,'.pdf">Open plot PDF</a></p>')}
  paste0('<section id="recruitment"><h2>Recruitment across P1, P2 and P3</h2>',
    '<p>P1 uses first-dose dates as a recruitment proxy. P2 and P3 use randomization dates. The phases are separate study cohorts.</p>',
    '<div class="table"><table><thead><tr>',paste0('<th>',c('Phase','Source status','Participants','Dated','Missing dates','First recruitment date','Last recruitment date','Date basis'),'</th>',collapse=''),'</tr></thead><tbody>',paste(rows,collapse=''),'</tbody></table></div>',
    '<p>Dates are the earliest and latest observed in each selected cohort (YYYY-MM-DD). Missing dates are excluded from the range and counted above. These periods do not include later study extensions absent from the selected sources. <a href="Recruitment_periods.csv">Recruitment periods CSV</a></p>',
    '<h3>P3 recruitment by site</h3><p><a href="Recruitment_by_site.csv">Site counts CSV</a></p>',
    if(r$map_available) paste0('<p><a href="Recruitment_by_site_map.html">Open interactive map</a></p>',
      '<iframe title="P3 recruitment by site" src="Recruitment_by_site_map.html" style="width:100%;height:680px;border:1px solid #d7e0e8" loading="lazy"></iframe>') else
      '<p>Map unavailable: place clinical_trial_sites.csv and screening_map_template.html in the project root. Site counts are available above.</p>',
    '<h3>Weekly recruitment</h3>',picture('Recruitment_by_week'),'<h3>Total recruitment by country</h3>',picture('Recruitment_by_country'),
    '<p>Country totals count study participation records. Missing phases are reported as unavailable. Coordinates identify sponsor-listed addresses and can represent campuses or offices.</p>',
    '<p><a href="Recruitment_README.txt">Definitions and source notes</a> | <a href="Recruitment_phase_sources.csv">Phase sources</a> | <a href="Recruitment_subjects.csv">Subject-level recruitment CSV</a></p></section>')
}


analyse_demographics <- function(root = "", adsl_file = NULL,
  output_dir = NULL, excel_fr = TRUE, strict = TRUE, show_diagram = interactive(), png_dpi = 180L,
  sites_file = "clinical_trial_sites.csv", template_file = "screening_map_template.html", p1_file = "", p2_file = "") {
  root <- bla1_project_root(root)
  if (is.null(output_dir) || !nzchar(output_dir)) {
    analysis <- file.path(root, "analysis")
    if (!dir.exists(analysis) && !dir.create(analysis, recursive = TRUE)) stop("Cannot create: ", analysis)
    output_dir <- tempfile(paste0("demographics_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = analysis)
  } else if (!grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", output_dir)) output_dir <- file.path(root, output_dir)
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)))
    stop("Choose an empty output directory; existing results are not overwritten: ", output_dir)
  if (!dir.exists(output_dir) && !dir.create(output_dir,recursive=TRUE)) stop("Cannot create: ",output_dir)
  selected <- bla1_resolve_input(root,adsl_file,"P3",output_dir)
  adsl_file <- selected$path; x <- selected$data
  x$SOURCE_FILE <- basename(adsl_file)
  x$SOURCE_MD5 <- unname(tools::md5sum(adsl_file))
  x$SOURCE_RECORD <- paste0(x$SOURCE_FILE, ":", x$source_row)
  nodes <- bla1_nodes(x)
  x$FIGURE1_NODES <- ""
  for (id in names(nodes$masks)) {
    j <- which(nodes$masks[[id]])
    x$FIGURE1_NODES[j] <- ifelse(x$FIGURE1_NODES[j] == "", id, paste(x$FIGURE1_NODES[j], id, sep = "|"))
  }
  leading <- c("SUBJID", "USUBJID", "SITEID", "TRT01P", "SAFETY_ARM", "IN_FIGURE1",
    "FASFL", "SAFFL", "MITTFL", "PPROTFL", "FIGURE1_BRANCH", "MITT_EXCLUSION", "PP_EXCLUSION",
    "PP_REASON_COUNT", "PP_CATEGORY_BASIS", "SAFETY_ARM_DIFFERS", "RTPCRBL", "ELECSBL", "SCOV2BL",
    "RANDDT", "TR01SDT", "DOSE2DT", "DOSE_INTERVAL_DAYS", "FASEREAS", "SAFEREAS", "MITTERS", "PPEREAS")
  x <- x[order(x$SUBJID), c(leading, setdiff(names(x), leading))]
  subjects <- x[x$IN_FIGURE1 == "Y", , drop = FALSE]
  excluded <- x[x$IN_FIGURE1 == "N", , drop = FALSE]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  bla1_csv(subjects, file.path(output_dir, "BLA_Figure1_subjects.csv"))
  if (excel_fr) bla1_csv(subjects, file.path(output_dir, "BLA_Figure1_subjects_Excel_FR.csv"), sep = ";")
  bla1_csv(excluded, file.path(output_dir, "BLA_Figure1_excluded_before_randomization_population.csv"))
  bla1_csv(nodes$counts, file.path(output_dir, "BLA_Figure1_counts.csv"))
  bla1_csv(x[x$PP_REASON_COUNT > 1, , drop = FALSE], file.path(output_dir, "BLA_Figure1_overlapping_reasons.csv"))
  bla1_csv(x[x$SAFETY_ARM_DIFFERS == "Y", , drop = FALSE], file.path(output_dir, "BLA_Figure1_safety_arm_switches.csv"))
  diagram <- bla1_diagram(nodes$counts)
  bla1_export_diagram(diagram, output_dir, png_dpi)
  if (show_diagram) bla1_draw_diagram(diagram)
  recruitment <- bla1_recruitment(subjects, root, output_dir, adsl_file,
    sites_file, template_file, p1_file, p2_file, png_dpi)
  report_file <- file.path(output_dir, "demographics.html")
  report <- paste(readLines(report_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  marker <- "<!-- RECRUITMENT -->"; at <- regexpr(marker, report, fixed = TRUE)[1]
  if (at < 1L) stop("Report insertion marker missing.")
  report <- paste0(substr(report, 1, at - 1), bla1_recruitment_html(recruitment, output_dir),
    substr(report, at + nchar(marker), nchar(report)))
  writeLines(enc2utf8(report), report_file, useBytes = TRUE)
  # A dictionary covers every column, including all preserved source columns.
  dictionary <- data.frame(column = names(x), origin = ifelse(names(x) %in% c(bla1_fields(),"COUNTRY"), "ADSL", "DERIVED/PROVENANCE"),
    definition = ifelse(names(x) %in% c(bla1_fields(),"COUNTRY"), "Original ADSL value; see the corresponding PHMPT define.xml and ADRG.", ""))
  definitions <- c(IN_FIGURE1 = "Y iff RANDFL=Y. Main CSV contains only these records; pre-figure exclusions are exported separately.",
    SAFETY_ARM = "TRT01A only when SAFFL=Y; efficacy grouping remains TRT01P.",
    FIGURE1_BRANCH = "One terminal branch of the efficacy population flow; safety membership is separate.",
    MITT_EXCLUSION = "FAS=Y and MITT=N, classified using SCOV2BL as positive or missing.",
    PP_EXCLUSION = "Display category among MITT=Y/PPROT=N: incorrect vaccine, dose2 outside window, no second dose, then other major deviation. Priority inferred from published counts.",
    PP_REASON_COUNT = "Number of the four recognized PP reasons among MITT=Y and PPROT=N; raw PPEREAS is retained.",
    PP_CATEGORY_BASIS = "Single source reason or inferred display priority for overlapping reasons; not a sponsor-supplied priority field.",
    SAFETY_ARM_DIFFERS = "Y iff SAFFL=Y and TRT01P differs from TRT01A; blank outside safety population.",
    FIGURE1_FOOTNOTE_GROUP = "Two invalid randomization IDs identified by ADRG; other RANDFL=N records are the remaining all-set exclusions.",
    DOSE_INTERVAL_DAYS = "DOSE2DT minus TR01SDT, with missing dates left blank.",
    FIGURE1_NODES = "Pipe-separated node IDs, matching BLA_Figure1_counts.csv. Safety nodes use actual treatment.",
    source_row = "One-based data row in the source ADSL/export, excluding header.",
    SOURCE_FILE = "Input basename.", SOURCE_MD5 = "MD5 of the exact input file.", SOURCE_RECORD = "Input basename and one-based source data row.")
  dictionary$definition[match(names(definitions), dictionary$column)] <- unname(definitions)
  ii <- grepl("^PP_HAS_", dictionary$column)
  dictionary$definition[ii] <- "Y if this reason appears in PPEREAS among MITT=Y/PPROT=N; multiple reasons can be Y."
  bla1_csv(dictionary, file.path(output_dir, "BLA_Figure1_dictionary.csv"))
  writeLines(c("Moderna demographics, recruitment and BLA Figure 1, v1.15.1", "",
    paste("Input:", normalizePath(adsl_file, winslash = "/")),
    paste("Main CSV records:", nrow(subjects)), paste("Pre-figure exclusion records:", nrow(excluded)),
    paste("Published count comparisons:", nrow(nodes$counts)),
    paste("Mismatches:", sum(nodes$counts$status != "MATCH")),
    "", "CSV: comma-delimited UTF-8 with BOM. Excel_FR: semicolon-delimited version of the same records.",
    "BLA_Figure1_diagram.pdf/.svg/.png: diagram redrawn from calculated ADSL counts. No source figure image is used.",
    "demographics.html: figure, recruitment map, weekly/country plots, and observed-versus-published count table.",
    "The published-count match statement appears at the end of the summary below the figure.",
    "Recruitment_by_site_map.html: interactive P3 map generated from the two assets in the project root.",
    "Recruitment_by_week and Recruitment_by_country: CSV data and PDF/PNG/SVG plots.",
    "Recruitment_periods.csv: first and last observed recruitment dates by phase, with cohort, date basis and missing-date counts; also shown in demographics.html.",
    "Recruitment_README.txt and Recruitment_phase_sources.csv: cohort definitions, date proxies, source status and hashes.",
    "BLA_Figure1_diagram_nodes.csv: diagram text and the count-node identifiers used for each box.",
    "Every displayed population/reason count uses observed_n. Published expected_n values are comparison targets only.",
    "Use TRT01P for efficacy/randomization and SAFETY_ARM/TRT01A for safety.",
    "Use population flags, not positive-test or case flags, to reconstruct enrollment/analysis populations.",
    "PP_EXCLUSION is a display reconstruction. All original overlapping reasons remain in PPEREAS and PP_HAS_*.",
    "No second dose uses the ADSL reason 'Did not receive Vaccination #2 per schedule'; it does not establish a withdrawal date.",
    "CUTOFFDT=2021-03-26 identifies this snapshot. EFFCODT can be earlier for individual subjects.",
    "Count agreement validates population bookkeeping; it does not establish the appropriateness of each exclusion.",
    "Reference figure: El Sahly et al., NEJM 2021, DOI 10.1056/NEJMoa2113017, Figure 1 and footnote.",
    "All study data are read from local PHMPT files; no reference CSV bundle or online download is required.",
    "Known BLA source: 125752_S3_M5_mrna-1273-p301_A_D_adsl.jmp in md-production-050125.zip.",
    "JMP import reads only the analysis columns and COUNTRY; unsupported unrelated JMP columns are not imported.",
    "Input_discovery_P1/P2/P3.csv record candidate files, rejection reasons, duplicate copies and selected sources.",
    "Equivalent file formats are compared on the fields used in this analysis; conflicting valid snapshots require an explicit choice.",
    "Original JMP/XPT/SAS formats are preferred when readable. An equivalent converted CSV is also supported.",
    "clinical_trial_sites.csv and screening_map_template.html are optional presentation assets; only the map depends on them.",
    "The input path and SOURCE_MD5 identify the exact source analyzed in this run."),
    file.path(output_dir, "README.txt"), useBytes = TRUE)
  capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  files <- list.files(output_dir, full.names = TRUE)
  bla1_csv(data.frame(file = basename(files), md5 = unname(tools::md5sum(files))), file.path(output_dir, "output_checksums.csv"))
  bad <- any(nodes$counts$status != "MATCH") || any(x$FIGURE1_BRANCH == "UNRESOLVED")
  message("Wrote ", nrow(subjects), " Figure 1 participants; ", nrow(excluded), " pre-figure exclusions. ",
    sum(nodes$counts$status == "MATCH"), "/", nrow(nodes$counts), " published counts match.\nRecruitment: ",
    paste(paste0(recruitment$sources$phase, "=", ifelse(is.na(recruitment$sources$participants), "unavailable", recruitment$sources$participants)), collapse = "; "),
    "\nObserved recruitment periods: ", paste(vapply(seq_len(nrow(recruitment$sources)), function(i) {
      s <- recruitment$sources[i, , drop = FALSE]
      paste0(s$phase, "=", if (s$status != "AVAILABLE") "unavailable" else if (is.na(s$first_recruitment_date)) "no dated records" else
        paste0(format(s$first_recruitment_date, "%Y-%m-%d"), " to ", format(s$last_recruitment_date, "%Y-%m-%d"), " (", s$date_basis, ")"))
    }, character(1)), collapse = "; "),
    "\n", normalizePath(output_dir, winslash = "/"))
  if (bad && strict) stop("Reconstruction has discrepancies. CSVs and the diagram were saved for review; no count was forced to match.")
  invisible(list(subjects = subjects, excluded = excluded, counts = nodes$counts, recruitment = recruitment, output_dir = output_dir))
}

# ---- RUN: starts automatically when the whole script is executed ----
analyse_demographics(root = PHMPT_ROOT, adsl_file = PHMPT_ADSL_FILE,
  output_dir = PHMPT_OUTPUT_DIR, excel_fr = PHMPT_EXCEL_FR, strict = PHMPT_STRICT,
  show_diagram = PHMPT_SHOW_DIAGRAM, png_dpi = PHMPT_PNG_DPI,
  sites_file = PHMPT_SITES_FILE, template_file = PHMPT_MAP_TEMPLATE,
  p1_file = PHMPT_P1_FILE, p2_file = PHMPT_P2_FILE)
