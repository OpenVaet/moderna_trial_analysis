# Moderna demographics, recruitment and BLA Figure 1, v1.16.0 -- R 4.3+.
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
# Leaflet and overview boundaries are embedded: the generated map works offline.
# Missing map assets do not prevent Figure 1 or recruitment charts.
# Audit exports: V01 weekly recruitment, V02 populations, V03 static site map.
# V03 needs clinical_trial_sites.csv coordinates; it does not need the HTML template.
# All three PNGs are refreshed in audit_assets only after they render successfully.
# First/last recruitment dates by phase are included in HTML and Recruitment_periods.csv.
# P1 PDF input needs pdftools (or an installed pdftotext command); P2 XPT needs haven.

# ---- SETTINGS: edit these values, then run the whole script ----
PHMPT_ROOT <- ""                   # Empty = detect project; or "D:/Projects/moderna_trial_analysis"
PHMPT_ADSL_FILE <- ""             # Empty = find local ADSL; or supply its path
PHMPT_OUTPUT_DIR <- ""            # Empty = a new dated folder under analysis/
PHMPT_EXCEL_FR <- TRUE             # Also write semicolon-delimited subject CSV
PHMPT_STRICT <- TRUE               # Stop after saving results if counts disagree
PHMPT_SHOW_DIAGRAM <- interactive() # Show the diagram in RStudio's Plots pane
PHMPT_PNG_DPI <- 180L              # Existing report figures
PHMPT_AUDIT_DPI <- 300L            # Print PNGs: 8.5 inches wide = 2550 px
PHMPT_AUDIT_ASSET_DIR <- "audit_assets" # Stable V01/V02/V03 files, relative to project root
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
  # Preserve the template layout, but replace the example's executable code
  # and remote stylesheets. No OSM tiles, plugins, jQuery or CDN requests remain.
  template <- gsub("(?is)<script\\b[^>]*>.*?</script\\s*>", "", template, perl = TRUE)
  template <- gsub("(?is)<link\\b[^>]*\\brel\\s*=\\s*['\"]stylesheet['\"][^>]*>", "", template, perl = TRUE)
  if (grepl(token, template, fixed = TRUE)) stop("Map placeholder must be inside a script element.")
  template <- sub("<title>.*?</title>", "<title>Moderna P3 recruitment by site</title>", template, perl = TRUE)
  insert <- function(html, marker, value) {
    at <- regexpr(marker, html, fixed = TRUE)[1]
    if (at < 0L) stop("Map template is missing ", marker)
    paste0(substr(html, 1L, at - 1L), value, substr(html, at, nchar(html)))
  }
  js <- paste0('var recruitmentSites = [', points, '];\nvar recruitmentTotal = ', total, ';\n',
    'var map = L.map("map", {attributionControl:true, minZoom:2, maxZoom:10, zoomAnimation:false}).setView([38,-97],4);\n',
    'var recruitmentCountries = ', bla1_map_asset("countries"), ';\n',
    'var recruitmentStates = ', bla1_map_asset("states"), ';\n',
    paste(bla1_map_background_js(), collapse = "\n"), '\n',
    paste(bla1_map_js(), collapse = "\n"))
  style <- paste0('<style>body{font:15px/1.45 Arial,sans-serif;color:#20354a;background:#f3f6fa}',
    'body>div:first-child{width:70%!important}body>div:nth-child(2){width:30%!important;max-height:100vh;overflow:auto}',
    '#map{background:#e7f1f7}.state-label{color:#697a88;font:10px Arial;text-align:center;white-space:nowrap;pointer-events:none}',
    '#currentMouseOver{padding:20px}h1{font-size:22px;margin:0 0 8px}h2{font-size:17px}',
    '.site-item{display:block;width:100%;text-align:left;background:white;border:1px solid #dce3ea;border-radius:5px;padding:9px;margin:6px 0;cursor:pointer}',
    '.site-item:focus,.site-item:hover{border-color:#2378aa;background:#edf5fb}input{box-sizing:border-box;width:100%;padding:9px;border:1px solid #99aabc}',
    '.muted{color:#52677a;font-size:13px}.site-detail{border-bottom:1px solid #dce3ea;padding:10px 0}',
    '@media(max-width:800px){body>div:first-child,body>div:nth-child(2){width:100%!important;max-height:none}#map{height:62vh!important}}</style>')
  template <- insert(template, "</head>", paste0('<style>', bla1_map_asset("leaflet_css"), '</style>', style))
  template <- insert(template, "</body>", paste0('<script>', bla1_map_asset("leaflet_js"), '</script>\n<script>', js, '</script>'))
  writeLines(enc2utf8(template), file.path(output_dir, "Recruitment_by_site_map.html"), useBytes = TRUE)
}

bla1_map_background_js <- function() c(
  'map.createPane("overview"); map.getPane("overview").style.zIndex="200"; map.getPane("overview").style.pointerEvents="none";',
  'L.geoJSON(recruitmentCountries,{pane:"overview",interactive:false,style:{color:"#9baab7",weight:0.7,fillColor:"#f5f6f2",fillOpacity:1}}).addTo(map);',
  'L.geoJSON(recruitmentStates,{pane:"overview",interactive:false,style:{color:"#96a8b6",weight:0.8,fillColor:"#f9faf6",fillOpacity:1}}).addTo(map);',
  'const stateLabels=L.layerGroup(); recruitmentStates.features.forEach(f=>{const p=f.properties; L.marker([p.lat,p.lon],{pane:"overview",interactive:false,keyboard:false,icon:L.divIcon({className:"state-label",html:p.postal,iconSize:[24,12]})}).addTo(stateLabels);});',
  'function updateLabels(){if(map.getZoom()>=4){if(!map.hasLayer(stateLabels))stateLabels.addTo(map);}else if(map.hasLayer(stateLabels))map.removeLayer(stateLabels);}',
  'map.on("zoomend",updateLabels);updateLabels();',
  'map.attributionControl.addAttribution("Overview boundaries: <a href=https://www.naturalearthdata.com/>Natural Earth</a> (public domain)");'
)

bla1_map_js <- function() c(
  '(function(){',
  'const host=document.getElementById("currentMouseOver"); host.replaceChildren();',
  'function el(tag,text,parent,cls){const x=document.createElement(tag); if(text!==null)x.textContent=text; if(cls)x.className=cls; if(parent)parent.appendChild(x); return x;}',
  'const fmt=n=>n.toLocaleString("en-US");',
  'el("h1","P3 recruitment by site",host);',
  'el("p",fmt(recruitmentTotal)+" randomized participants across "+recruitmentSites.filter(s=>s.n>0).length+" sites",host);',
  'el("p","P301 Figure 1 population (RANDFL=Y). Includes both randomized arms.",host,"muted");',
  'el("p","Bubble area represents participants at that location. Sites sharing coordinates appear together. Coordinates refer to sponsor-listed addresses and can represent campuses or offices.",host,"muted");',
  'el("p","Click a bubble or select a site below. The embedded country/state overview works offline; it does not include street-level detail.",host,"muted");',
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
  'function fitSites(){if(markers.length)map.fitBounds(L.featureGroup(markers).getBounds().pad(0.12),{maxZoom:7});}fitSites();',
  'const reset=L.control({position:"topleft"});reset.onAdd=function(){const b=L.DomUtil.create("button","site-item");b.type="button";b.textContent="All sites";b.title="Show all mapped recruitment sites";b.style.cssText="margin:0;padding:6px 10px;width:auto;box-shadow:0 1px 5px #0004";L.DomEvent.disableClickPropagation(b);L.DomEvent.on(b,"click",fitSites);return b;};reset.addTo(map);',
  'if(window.ResizeObserver){new ResizeObserver(()=>map.invalidateSize({pan:false})).observe(document.getElementById("map"));}', 
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
    "The map is P3 only because clinical_trial_sites.csv describes P301. It preserves the supplied root HTML template layout, replaces its legacy scripts with embedded assets, and adds proportional-area markers, actual SITEID labels and source attribution.",
    "Coincident coordinates are grouped into one bubble; all distinct sites and their individual counts remain accessible in the popup and searchable sidebar.",
    "Coordinate estimates may describe campus or administrative addresses; consult the clinical_trial_sites source bundle for details.",
    "The HTML map is self-contained: Leaflet, country boundaries and US state boundaries are embedded. No map tiles, CDN requests, API keys or local web server are required. The overview has no street-level detail.",
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


# ---- AUDIT ASSETS: V01/V02/V03, print and Substack PNGs ----
# Uses the same calculated tables as the report. No screenshots or downloads.
bla1_audit_number <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

bla1_audit_png <- function(path, draw, width, height, dpi) {
  if (length(dpi) != 1L || !is.finite(dpi) || dpi < 150 || dpi > 600)
    stop("PHMPT_AUDIT_DPI must be between 150 and 600.")
  type <- if (capabilities("cairo")) "cairo" else if (.Platform$OS.type == "windows") "windows" else "Xlib"
  success <- FALSE
  grDevices::png(path, width = width, height = height, units = "in", res = dpi,
    pointsize = 13, bg = "white", type = type)
  device <- grDevices::dev.cur()
  on.exit({
    if (device %in% grDevices::dev.list()) grDevices::dev.off(device)
    if (!success && file.exists(path)) unlink(path)
  })
  draw()
  grDevices::dev.off(device)
  con <- file(path, "rb"); on.exit(close(con), add = TRUE)
  sig <- readBin(con, "raw", n = 8L)
  if (!identical(sig, as.raw(c(137,80,78,71,13,10,26,10)))) stop("Invalid PNG: ", path)
  seek(con, 16L); dims <- readBin(con, "integer", n = 2L, size = 4L, endian = "big")
  if (any(abs(dims - round(c(width, height) * dpi)) > 1L)) stop("Unexpected PNG dimensions: ", path)
  success <- TRUE
  dims
}

bla1_audit_weekly <- function(r) {
  colors <- c(P1 = "#A65B22", P2 = "#16816F", P3 = "#286AA0")
  english_date <- function(d) paste(format(d,"%d"),month.abb[as.integer(format(d,"%m"))],format(d,"%Y"))
  weeks <- sort(unique(r$weekly$WEEK_START))
  left <- as.Date(format(min(weeks), "%Y-%m-01"))
  ticks <- seq(left, max(weeks) + 7, by = "month")
  old <- graphics::par(mfrow = c(3,1), mar = c(2.7,4.8,3.2,0.8),
    oma = c(3.2,0,3.7,0), family = "sans", cex = 1)
  on.exit(graphics::par(old))
  for (p in names(colors)) {
    s <- r$sources[r$sources$phase == p, , drop = FALSE]
    w <- r$weekly[r$weekly$PHASE == p, , drop = FALSE]
    available <- nrow(s) == 1L && s$status == "AVAILABLE"
    ymax <- if (all(is.na(w$subjects_recruited))) 1 else max(1, w$subjects_recruited, na.rm = TRUE)
    graphics::plot(as.numeric(weeks), rep(0,length(weeks)), type = "n",
      ylim = c(0,ymax * 1.12), xlim = c(as.numeric(left)-4,as.numeric(max(weeks))+7),
      xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
    yy <- pretty(c(0,ymax), n = 4); yy <- yy[yy >= 0 & yy <= ymax * 1.12]
    graphics::abline(h = yy, col = "#DCE3E8", lwd = 0.7)
    if (available) graphics::rect(as.numeric(w$WEEK_START)-2.5,0,
      as.numeric(w$WEEK_START)+2.5,w$subjects_recruited,col=colors[p],border=NA)
    else graphics::text(mean(graphics::par("usr")[1:2]),0.5,
      "Source unavailable\n(not zero recruitment)", col = "#765A36", cex = 1)
    graphics::axis(2,at=yy,labels=bla1_audit_number(yy),las=1,tck=-0.015,col="#74838C",cex.axis=0.9)
    graphics::axis(1,at=as.numeric(ticks),labels=format(ticks,"%Y-%m"),
      tck=-0.025,col="#74838C",cex.axis=0.83)
    graphics::mtext("Participants", side=2, line=3.6, cex=0.9)
    graphics::mtext(paste0(p,"  |  ", if(available) paste0(bla1_audit_number(s$participants)," participants") else "unavailable"),
      side=3,line=1.85,adj=0,font=2,cex=1.13,col=colors[p])
    period <- if (!available || is.na(s$first_recruitment_date)) "No observed recruitment dates" else
      paste(english_date(s$first_recruitment_date),"to",english_date(s$last_recruitment_date))
    coverage <- if(available) paste0("  |  ",s$dated," dated; ",s$missing_dates," missing") else ""
    graphics::mtext(paste0(period,coverage),side=3,line=0.45,adj=0,cex=0.78,col="#4E606B")
  }
  graphics::mtext("V01  Weekly recruitment by clinical trial phase",outer=TRUE,side=3,line=2,at=0.035,adj=0,cex=1.3,font=2)
  graphics::mtext("P1: first-dose proxy. P2 and P3: randomization dates.",outer=TRUE,side=3,line=0.5,at=0.035,adj=0,cex=0.9,col="#4E606B")
  graphics::mtext("Monday-start weeks; a common calendar axis and separate vertical scales.",outer=TRUE,side=1,line=0.8,at=0.035,adj=0,cex=0.85,col="#4E606B")
  graphics::mtext("Local PHMPT records; missing dates are excluded from weekly bars.",outer=TRUE,side=1,line=2,at=0.035,adj=0,cex=0.85,col="#4E606B")
}

bla1_audit_flow <- function(counts) {
  # A two-column layout replaces the original four-column diagram only in V02.
  # Every number still comes from observed_n; status never substitutes reference counts.
  n <- function(id) {
    i <- match(id,counts$node_id)
    if(anyNA(i)) stop("Missing audit diagram node: ",paste(id[is.na(i)],collapse=", "))
    counts$observed_n[i]
  }
  num <- function(id) bla1_audit_number(n(id))
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(xscale=c(0,100),yscale=c(100,0)))
  on.exit(grid::popViewport())
  txt <- function(label,x,y,size=12,col="#20313E",bold=FALSE,just="left")
    grid::grid.text(label,grid::unit(x,"native"),grid::unit(y,"native"),just=c(just,"centre"),
      gp=grid::gpar(fontfamily="sans",fontsize=size,col=col,fontface=if(bold)"bold" else "plain",lineheight=1.16))
  box <- function(ids,label,x,y,w,h,kind="main",size=12,bold=FALSE) {
    mismatch <- any(counts$status[match(ids,counts$node_id)] != "MATCH")
    fill <- if(mismatch) "#FFF1D9" else switch(kind,root="#234D67",exclusion="#F3F5F6",main="#EAF2F8")
    grid::grid.rect(grid::unit(x+w/2,"native"),grid::unit(y+h/2,"native"),
      grid::unit(w/100,"npc"),grid::unit(h/100,"npc"),
      gp=grid::gpar(fill=fill,col=if(mismatch)"#B86616" else "#6D8492",lwd=0.8))
    txt(label,x+w/2,y+h/2,size,if(kind=="root" && !mismatch)"white" else "#20313E",bold,"centre")
  }
  arrow <- function(x,y1,y2) grid::grid.lines(x=grid::unit(c(x,x),"native"),y=grid::unit(c(y1,y2),"native"),
    arrow=grid::arrow(length=grid::unit(0.07,"inches"),type="closed"),gp=grid::gpar(col="#617886",fill="#617886",lwd=0.9))
  txt("V02  Randomization and analysis populations",4,2.3,17,bold=TRUE)
  txt("Moderna P301 | BLA cutoff: 26 March 2021",4,5,11,col="#4E606B")
  box("FIGURE1",paste0(num("FIGURE1")," randomized participants"),21,7,58,4.5,"root",14,TRUE)
  for(pre in c("P","V")) {
    x <- if(pre=="P")4 else 53; center <- x+21.5
    id <- function(k) paste0(pre,"_",k)
    nn <- function(k) num(id(k))
    arm <- if(pre=="P") "Placebo" else "mRNA-1273"
    grid::grid.lines(x=grid::unit(c(50,50,center,center),"native"),y=grid::unit(c(11.5,12.6,12.6,14),"native"),
      gp=grid::gpar(col="#617886",lwd=0.9))
    box(id(c("RANDOMIZED","NO_DOSE")),paste0(arm,": ",nn("RANDOMIZED")," assigned\nNo study dose: ",nn("NO_DOSE")),x,14,43,6,"main",13,TRUE)
    arrow(center,20,22)
    box(id(c("FAS","SAFETY")),paste0("Full analysis: ",nn("FAS"),"\nSafety (actual treatment): ",nn("SAFETY")),x,22,43,6,"main",12.5)
    arrow(center,28,30)
    unresolved <- n(id("MITT_EXCLUDED"))-sum(n(id(c("BASELINE_POSITIVE","BASELINE_MISSING"))))
    label <- paste0("Excluded before modified ITT: ",nn("MITT_EXCLUDED"),
      "\n",nn("BASELINE_POSITIVE")," positive at baseline\n",nn("BASELINE_MISSING")," missing baseline data",
      if(unresolved != 0)paste0("\nUnresolved: ",unresolved) else "")
    box(id(c("MITT_EXCLUDED","BASELINE_POSITIVE","BASELINE_MISSING")),label,x,30,43,10,"exclusion",12)
    arrow(center,40,42)
    box(id("MITT"),paste0("Modified intention-to-treat\n",nn("MITT")," participants"),x,42,43,6,"main",13,TRUE)
    arrow(center,48,50)
    reasons <- c("INCORRECT_VACCINE","DOSE2_OUTSIDE_WINDOW","NO_SECOND_DOSE","OTHER_MAJOR_DEVIATION")
    unresolved <- n(id("PP_EXCLUDED"))-sum(n(id(reasons)))
    label <- paste0("Excluded before per-protocol: ",nn("PP_EXCLUDED"),
      "\n\n",nn("INCORRECT_VACCINE")," received incorrect vaccine\n",
      nn("DOSE2_OUTSIDE_WINDOW")," second dose outside window\n",
      nn("NO_SECOND_DOSE")," no second dose per schedule\n",
      nn("OTHER_MAJOR_DEVIATION")," other major deviations",
      if(unresolved != 0)paste0("\nUnresolved: ",unresolved) else "")
    box(id(c("PP_EXCLUDED",reasons)),label,x,50,43,16,"exclusion",11.8)
    arrow(center,66,68)
    box(id("PP"),paste0("Per-protocol population\n",nn("PP")," participants"),x,68,43,6,"main",13.5,TRUE)
  }
  notes <- c(
    paste0(num("ORIGINAL")," ADSL records; ",num("EXCLUDED_BEFORE")," excluded before the figure: ",num("DUPLICATE_RANDOMIZATION")," duplicates and ",num("OTHER_ALL_SET_EXCLUSION")," other exclusions."),
    paste0("Safety uses actual treatment: ",num("SWITCH_P_TO_V")," placebo-assigned participants enter the vaccine safety arm;"),
    paste0(num("SWITCH_V_TO_P")," vaccine-assigned participants enter the placebo safety arm."),
    "Full analysis and efficacy populations use randomized treatment; safety is a separate population.",
    "Counts use ADSL flags. Overlapping per-protocol reasons use inferred display priority.",
    "All raw reasons remain in the CSV. No-second-dose classification does not establish withdrawal.",
    "Reference: El Sahly et al. (2021), Figure 1. All counts are calculated from the input records.")
  for(i in seq_along(notes)) txt(notes[i],4,78+(i-1)*2.15,10.3,col="#4E606B")
  good <- sum(counts$status=="MATCH")
  txt(paste0(good,"/",nrow(counts)," published counts match"),4,96,12,
    col=if(good==nrow(counts))"#1C7056" else "#B86616",bold=TRUE)
}

bla1_audit_albers <- function(lon,lat,lon0=-96,lat0=23,parallels=c(29.5,45.5)) {
  rad <- pi/180; p <- parallels*rad
  n <- (sin(p[1])+sin(p[2]))/2
  cc <- cos(p[1])^2+2*n*sin(p[1])
  rho <- 6371*sqrt(pmax(0,cc-2*n*sin(lat*rad)))/n
  rho0 <- 6371*sqrt(cc-2*n*sin(lat0*rad))/n
  theta <- n*((lon-lon0+180)%%360-180)*rad
  cbind(x=rho*sin(theta),y=rho0-rho*cos(theta))
}

bla1_audit_map <- function(r) {
  s <- r$sites[r$sites$subjects_recruited>0,,drop=FALSE]
  total <- sum(s$subjects_recruited)
  valid <- is.finite(s$latitude)&is.finite(s$longitude)
  # Explicit geographic coverage: never silently drop a mapped site off the chart.
  region <- rep(NA_character_,nrow(s))
  region[valid & s$longitude>=-125 & s$longitude<=-66 & s$latitude>=24 & s$latitude<=50] <- "CONUS"
  region[valid & (s$longitude<=-129 | s$longitude>=170) & s$latitude>=50 & s$latitude<=72] <- "AK"
  region[valid & s$longitude>=-161 & s$longitude<=-154 & s$latitude>=18 & s$latitude<=23] <- "HI"
  mapped <- valid & !is.na(region)
  if(!any(mapped)) stop("V03 requires valid US coordinates in clinical_trial_sites.csv; no sites could be plotted.")
  key <- paste(region[mapped],sprintf("%.4f",s$longitude[mapped]),sprintf("%.4f",s$latitude[mapped]),sep="|")
  rows <- split(which(mapped),key)
  points <- do.call(rbind,lapply(rows,function(ii)data.frame(region=region[ii[1]],lon=s$longitude[ii[1]],
    lat=s$latitude[ii[1]],n=sum(s$subjects_recruited[ii]),sites=length(ii))))
  geometry <- bla1_audit_boundaries()
  max_n <- max(points$n,600)
  old <- graphics::par(no.readonly=TRUE);on.exit(graphics::par(old))
  panel <- function(fig,new=TRUE) {
    graphics::par(fig=fig,mar=rep(0,4),new=new,family="sans")
    graphics::plot.new();graphics::plot.window(c(0,1),c(0,1),xaxs="i",yaxs="i")
  }
  panel(c(0,1,0.88,1),FALSE)
  graphics::text(.035,.77,"V03  Phase 3 recruitment by site",adj=0,cex=1.35,font=2,col="#20313E")
  graphics::text(.035,.43,paste0(bla1_audit_number(total)," randomized participants | ",nrow(s)," recruiting sites"),adj=0,cex=1.05,col="#405663")
  graphics::text(.035,.12,"Circle area represents participants; coincident locations are combined.",adj=0,cex=.9,col="#405663")
  draw_region <- function(which,fig) {
    g <- geometry[if(which=="CONUS")!geometry$state %in% c("AK","HI") else geometry$state==which,,drop=FALSE]
    project <- switch(which,CONUS=function(x,y)bla1_audit_albers(x,y),
      AK=function(x,y)bla1_audit_albers(x,y,-154,50,c(55,65)),
      HI=function(x,y)cbind(x=x*cos(20*pi/180),y=y))
    xy <- project(g$lon,g$lat)
    graphics::par(fig=fig,mar=c(0,0,0,0),new=TRUE)
    graphics::plot.new();graphics::plot.window(range(xy[,1]),range(xy[,2]),asp=1,xaxs="r",yaxs="r")
    for(ii in split(seq_len(nrow(g)),g$polygon)) graphics::polygon(xy[ii,1],xy[ii,2],col="#EEF1F3",border="#99A8B2",lwd=.65)
    pts <- points[points$region==which,,drop=FALSE]
    if(nrow(pts)) {
      pts<-pts[order(-pts$n),,drop=FALSE]; pp<-project(pts$lon,pts$lat)
      graphics::symbols(pp[,1],pp[,2],circles=sqrt(pts$n),inches=.22*sqrt(max(pts$n)/max_n),
        add=TRUE,bg=grDevices::adjustcolor("#286AA0",alpha.f=.55),fg="#174E76",lwd=.8)
    }
    if(which!="CONUS") graphics::mtext(if(which=="AK")"Alaska" else "Hawaii",side=1,line=-.5,adj=.05,cex=.85,col="#405663")
    if(which=="CONUS") {
      labs <- unique(g[c("state","label_lon","label_lat")])
      labs <- labs[!labs$state %in% c("DC","DE","RI","CT","MA","NJ","MD"),,drop=FALSE]
      ll <- project(labs$label_lon,labs$label_lat)
      graphics::text(ll[,1],ll[,2],labs$state,cex=.65,col="#4D5F6A")
    }
  }
  draw_region("CONUS",c(.025,.975,.29,.88))
  draw_region("AK",c(.035,.235,.13,.31))
  draw_region("HI",c(.255,.405,.13,.26))
  panel(c(.48,.97,.12,.31))
  graphics::text(.02,.88,"Participants per location",adj=0,font=2,cex=.95,col="#20313E")
  nn<-c(100,300,600); xx<-c(.15,.47,.8)
  graphics::symbols(xx,rep(.53,3),circles=sqrt(nn),inches=.22*sqrt(max(nn)/max_n),add=TRUE,
    bg=grDevices::adjustcolor("#286AA0",alpha.f=.55),fg="#174E76",lwd=.8)
  graphics::text(xx,.12,bla1_audit_number(nn),cex=.9,col="#405663")
  panel(c(0,1,0,.12))
  coverage <- paste0(sum(mapped),"/",nrow(s)," sites mapped; ",bla1_audit_number(sum(s$subjects_recruited[mapped])),"/",bla1_audit_number(total)," participants represented.")
  if(any(!mapped)) coverage<-paste0(coverage,"  Unmapped/outside extent: ",sum(!mapped)," sites.")
  graphics::text(.035,.78,coverage,adj=0,cex=.83,col=if(all(mapped))"#405663" else "#9B4B0A")
  graphics::text(.035,.48,"P301 randomized population (RANDFL=Y). Site coordinates: clinical_trial_sites.csv.",adj=0,cex=.83,col="#405663")
  graphics::text(.035,.18,"Natural Earth 1:110m boundaries. Alaska/Hawaii insets are not at the main-map scale.",adj=0,cex=.83,col="#405663")
}

bla1_export_audit_assets <- function(r,counts,root,output_dir,asset_dir="audit_assets",dpi=300L) {
  if(is.null(asset_dir)||!nzchar(asset_dir))asset_dir<-"audit_assets"
  asset_dir<-bla1_path(root,asset_dir)
  # Render a complete run-specific set first. Failed runs do not publish partial assets.
  stage<-file.path(output_dir,"audit_assets")
  if(!dir.exists(stage)&&!dir.create(stage,recursive=TRUE))stop("Cannot create ",stage)
  files<-c("V01_weekly_recruitment_by_phase.png","V02_randomization_analysis_populations.png","V03_recruitment_by_site.png")
  dims<-rbind(c(8.5,10),c(8.5,12.8),c(8.5,7.5))
  draw<-list(function()bla1_audit_weekly(r),function()bla1_audit_flow(counts),function()bla1_audit_map(r))
  pixels<-matrix(NA_integer_,3,2)
  for(i in seq_along(files)) pixels[i,]<-bla1_audit_png(file.path(stage,files[i]),draw[[i]],dims[i,1],dims[i,2],dpi)
  manifest<-data.frame(visual=c("V01","V02","V03"),file=files,width_px=pixels[,1],height_px=pixels[,2],dpi=dpi,
    generated_at=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),run_directory=normalizePath(output_dir,winslash="/"),
    source_adsl_md5=r$sources$source_md5[match("P3",r$sources$phase)],
    md5=unname(tools::md5sum(file.path(stage,files))),stringsAsFactors=FALSE)
  bla1_csv(manifest,file.path(stage,"V01_V03_demographics_manifest.csv"))
  captions<-c("AUDIT ASSETS: DEMOGRAPHICS (V01-V03)","",
    "PNG files are generated from the same calculated counts as demographics.html, at 8.5 inches wide and the configured audit DPI (default 300 = 2550 pixels).",
    "Opaque white backgrounds, large sans-serif text and dedicated narrow-column layouts. Upload the original PNGs to Substack; do not screenshot the HTML.","",
    "V01: Weekly recruitment by phase. P1 uses first dose as a recruitment proxy; P2/P3 use randomization dates. Monday-start weeks and independent vertical scales. Unavailable phases are labelled, not zero-filled.",
    "V02: Independently drawn randomization/analysis diagram from observed ADSL counts. Safety follows actual treatment; other populations follow randomized treatment. The count-match statement is last.",
    "V03: Static P301 recruitment map. Marker area represents participants; coincident coordinates are combined. Alaska and Hawaii have separate insets. Unmapped/out-of-extent counts are disclosed on the image. Coordinates are not jittered.","",
    "V03 uses the embedded Natural Earth geometry and clinical_trial_sites.csv. No screenshot tool, web browser, API key, network access or screening_map_template.html is needed for this PNG.",
    "The existing interactive map still uses screening_map_template.html. Existing report/data/PDF/SVG outputs are preserved.",
    "These stable V01-V03 filenames are refreshed on successful runs. Other audit_assets files are not changed. A run-specific copy and manifest remain under analysis/.")
  writeLines(captions,file.path(stage,"V01_V03_demographics_README.txt"),useBytes=TRUE)
  if(!dir.exists(asset_dir)&&!dir.create(asset_dir,recursive=TRUE))stop("Cannot create ",asset_dir)
  for(name in c(files,"V01_V03_demographics_manifest.csv","V01_V03_demographics_README.txt")) {
    src<-file.path(stage,name);dst<-file.path(asset_dir,name)
    if(normalizePath(src,winslash="/")!=normalizePath(dst,winslash="/",mustWork=FALSE)) {
      if(!file.copy(src,dst,overwrite=TRUE))stop("Cannot refresh audit asset: ",dst)
      if(unname(tools::md5sum(src))!=unname(tools::md5sum(dst)))stop("Audit asset copy verification failed: ",dst)
    }
  }
  message("Audit PNGs: ",normalizePath(asset_dir,winslash="/")," (V01, V02, V03; ",dpi," dpi)")
  invisible(list(directory=asset_dir,manifest=manifest))
}


# Compact base-R boundary table derived from the existing embedded Natural Earth
# 5.1.2 1:110m state geometries (public domain); no extra package is required.
bla1_audit_boundaries <- function() {
  hex <- paste0(c(
    "1f8b0800000000000203bd9dcbce65b9919dc7c7ef523f78bf0c0b52db2aabab64b464c9ed8921030d4f04b7016be2b73757c4fe330f79f65a27",
    "bb25e44048e5c9afc8cd5b3022180cfe9f7ffdcbfffb5ffffabf1ffff7af7ffeebbf3cfeb2fedf5ffefcd7f5bffff92f7ff91ffe37fb7f7ffeeb",
    "7f888f9f7f79fc30e6c7acbd3f4aff48a3cdc70f337fe4561ea57d843ad345cdf011731c46cdd43995d75f9d1a1c2a56500e5954d76a7168460e",
    "8ddafd9b62a2d5c58fd022a0183b67522bd60721c5cca9d286fd3a53e5d0fa4bc6afa33751548fad19957ae3d428113ff6beeaa5d06ce8d045a5",
    "495b983e428cd328d10fe923f56254ab2d0aeaaa31cc9a0435bd8961f28e481fb5e4d5f3f563aec1e454af068d312787d63f0e503dd4c1a93666",
    "025523efd455a14dc0ba9adaa3a25a7d94f231a7eaae3e5b07d5f950af7e083502aa238b26e651adc28a71e263dd5202557212936b5a491963c9",
    "27e0fa76506955c8a996ba35308c2aa8345056fe98b98b85b16acaa056516abde25f17d5a2aa308d04a88698d40434267289b51a16a243995657",
    "3f725d93d3aa0bb4baf6516aead7470d41e5597c04aba06a0b0d9334cfcca18e2216547a10458d5c6c85b55455592318d59b2aabd6b59160b5f6",
    "5404953101d7afaa1f6a1e268fe29c95532d44fb35d7282a6cb6761a46517617e411e4641335f6524d02ce5214d5a79535b9d4c20e10cab5b126",
    "45f90655526b822a2581aaa5a9b2a655d8cba09fb576c258ada8594be2d4d20a1e65ac3ff8925e501915d4d22086a0da5aa58b5a63d80495201e",
    "c6472b292a2a588dbd8ecaa9148b7dfd0c7c9e7648913505e7920161f2c51feb6a1bca9a5c8e805ae3b2caca7ceb291f23f6ea105fd6657540b0",
    "8e5025e534ad81ad072a4b970ccd90228b0a53ecd3eb5fad8145a83e2d4302623a4ca1da8dec5f9e7a6043b8b452fba41043e6c892b7e53d55a6",
    "335474a3a4946c75cd16aaf8a4988d1a718afa7a329db4d165b398511c2a8d6e72ef14f3f4f8f90f8f1f625c7f2d235d13b4e097b07e713575b6",
    "fa85cbf8cfaa735170eb576c14c6658a8525203e6bc57265d8da3a67ba4aab9aebdda4525c8257721080d0176b7b539eed7aa10481550872dbab",
    "96f695546b97446d5edee229373f22648995d7791f2f9ba447c7ca6c926bdd34b95aa2e646750ec68ee0e22cce35d1dcc52dc161da6197736559",
    "440641420a68ed20c6d58cd5c5b9debdd23287f8b8a504e7e6aa7213630bae25535e0bd461c1b519e2556f57dce8307ba03287a8175ab2f25a9c",
    "5371f12acf6c0dc1e53eac1d7d606fe1dca7d2385a925c0dc9e67c4c2529ae8d6096522a53b6b7cfb529c0a26a5561d3b65b9867490cdb5ab997",
    "c2dab3eae602a5c39a516d9fe45cce4ba57505587235ba25bab4d724b9e4dc486ad8d6740a3d7f355a05375d7b8da5cbf6dae82ed320ea6e49a3",
    "989adbe0b051cd18d5ccf356a2ec9651b2f933d624141c66537417ca54ddb266d318d976b51a24b72c32f3dbe452a3e27a707fd22a6ebee1acde",
    "36ca50dc32d6872bbe4a5abddb76f3e3975fdfa88f31c047967cefefcfd8a68b526c576c05f6ac25736c53b939b6e9ef1c7bb60528b5db15043b",
    "6c144e6df60ec736e34960cf9618c536ab8e539b85c8b1cddca4d86ebb726c338429b65bd5027bb2d009b50c82124d5ccf96699df83534f70aac",
    "d5d138973ed6dabab828ca7b510439b72b969cdb1555da8e96e6fcaa1ed3f52c977d79fce6a7f56b5d5f9fd7d7c709e51cadaddda6558ad815e6",
    "3317460497b32d08c2adb1eab982ebd8e315874d6c7123da32a4f58e1aacded96d55f3ef83b44e013a8d2cae45c312f45ac1f5519d2ba69251ae",
    "84645ce8a61a512ea562dd32a76c4568295baf54d37808b7f6c4b09ab9b81ad4e72d2e47abb6d49cf5e72d89b3b894a76c6e853cc76c19b92aae",
    "8db537c5b59bc426bb59cfbefa955ba628bab9378161f2adee5d582b122b71fdfbc2aaf907395747f35a4756d8b2510c1bad682e2e73617133b5",
    "a1b8b52936fc1c62aa8a4bb5dbdc9bb773fe1b3baf7dfdbcbac66cfd1cfbdde73d715077f179ed6eee7d724b6f4ad5ca0bb30fc52d3dccca4bd1",
    "4cd17fe7f7f58b5bfb7fc2a8c5a52ddc8dee572ead5d083faf5d26490e078ee0fadd9afcca85abbd3d768ead5f6149af9f1bfcf0825bff6e9f97",
    "534a8a836b2dd9016291186c3870b772f45b7b6f5cdc92506b9b7a241cf5d4acb8a5400d70318fa4b83a97a4585cca4d72b9613227e8ab5372d0",
    "dad7cf53575b5ab6e68e91fffdad9d8f9f7e7daf93c348b4d3a4be8cb52fdcab8ecfb8d36610dc668330eeb4691877da488c3b6dae7beed58663",
    "dc691332eeb431197798ac0c3b2d60ca1d16356dc661a153eeb0f869330e0f02e34e8f04e34e0fc73df7ea3161dce18061d8e9cf61dce91f62dc",
    "e96f62dce9bf62dce90f63dce95fbbe75efd758c3bfd7f8c3bfd89943bfc93f7dcabbf9371a7ff94719b339641a7679771a7a79871a7e7997265",
    "5a8c4019bab1a16292269c5234bd86f07d111abda817113d10a5ab3c09450b025886bd280b713fed8a73507db28a0bdda326a2eae38e3eb6314b",
    "7d46cd21a408d1006a6c3b0e316c6c5b6baa19e39a9f6b33988a6b339bd88b439757eb3079514a298a2bd5cf646a2f5596b716b989bdd46479cd",
    "a5ad19636ad4a27b3124139a41b1f721b9392c222bc39a97b3aed9c6dc53d2239663b98e8db3aa572b36313cfef4e36d81cbecc81efe5646dbc1",
    "fd0b19f8d26401ee7d28c0af03c2a0d5803e6cb39ad3f63e06c2237285868d2cb93a936da633a724c13e2e2da306d10ef83aeaf0795a70bc2ec0",
    "16a6b9aa7acd538223776b74eb6a5060d495e6e13ab957092644590274958983159e19d349d290601fe3524aa2e0127673578666d560319f36c0",
    "24b956e3d5e63735f774a957aeae09b0269f6379e8121760333b4dc56578fd7c09b89ecdc1d25de3bdcc410ef6decc3f9c4692dc84ee625c569f",
    "0877e8f4c0538bfb136086dbc204845af70b6ca3db29421c26fe39585bb31382dc67978df1538950e498d4d1e365fb4c09668f246e7a124e0f6f",
    "8923e859880860bd88434030c917d7af10c15a56c7c78fffdd830a42593a56b61387e95ad2ccf591d7182d39b0834b18e6b89a1c05082d2bad19",
    "0d10617f1ccc08e87ee4a53d65551e629656398bab163e28c09e569b33d6a809570e56180b0063d32596eee07c07e2dc21c37a78d396066b7481",
    "4b7476dd1658d519f6cdd425d635a000d73e50755b6088aedf43ebba2d881ecf501cdd30a45cae6bb2e672051aa99a61d1677492abd3bc2d0888",
    "5ebf0f9748a21b533270b637ddd8ca92c2192133efaa1e6b4165b35f8b2e113b4486425adf7c635d3b77c69632e63bd0be71d4f2667ee7255e33",
    "d496f86636f605ac351dfdec55cd097ce33780c9264f732120270fc2c20106dd8d01d1eda8b9b63760ca5660bc0c7c0e9a3881202b6aa897388e",
    "cb82cab0a4e4508725b75b7250ad98f732343d7ef5e3fd7484e618f15f2e7d66078f094ec173c950f05c84143c9735050f39c16b3e040f6fcb21",
    "ca28780a47de8d87b8e5e02ebf795b8e0d81b7e5d862785b8e4d8b80b6ea23403b42e2a00568f87e996b9660827683567b840307e14ac0b8a4a8",
    "b8015fac75e3325aa6042bbc6b9812a9ab46e330d266777883953c86adb732c4cc3125a81b584a4f126c305db15cc26812eca598505e824c7ca3",
    "9d1d1517642d8ab908e5b32cbb6b899d5ae3d06034f1d4476e929b3821ce4bff8c5915b8c66d29978b9b70f5286ea990ebf7892b4a629c4dcf5f",
    "8aefeaf59825163b9c3cf12396942498d00273cc490931ed08c7b8a4472f587baf130e3a5d6d7bcc88bb2e0a6c1f198b38e34f2d8cdfec03f9f1",
    "abdf79c08a7b2ca1f79b740838a6cffecd736c60815a82107a8f1c1260bcf6aa24b88c10948bab0a5c63822b49063605565ca19b0e7a10170161",
    "a23c6da60c7bd9750558b3b973edee14e5063a7b5e5c951fe8fee6b72042b0db55a22d51de89018776069af2a2866f7c23f8397372fe46f07e8a",
    "95c72f7fbcf176e3af100125cdf9157a759d736ed7d62877a87f943bf449c61dfa29c34e7d9771a702cdb8532367dca9e2d3f20e9b41709b1142",
    "b9c3aa61dc6926d1f20e0177cfbd08cc7bec5500df72a738bf2fec656fb8c73ea73bb85ede73b7cba23e7ef9f9d67912102457a05fe042e70eee",
    "428d80a79024d8abd025e0ab1827e0ebc620c0684badb80bea1d987121787e03183f9d5a1c9cc811f01e8447fa1bb0762d8f6f01e3524c004617",
    "cf14acc1abeebe8019383ef0ab61f30d967c7ea5f8667ec989d81ebffba7db839461d7ab47cef919da4fc938779cba71f038c61325eee782e213",
    "f783460e1e27970c3c4f4239779cac2a703baa1535ef67bfbc2dc761b228713b9c16dc93907b477d11990c3c4530e15ef57d0e569c2c154b9a20",
    "b988f581068721be1026cef4f95553d6e0cde1c83df87ada72cfdd1cdfdc8337074277e0dd09d33d77736445c0d743b07bf0e658ed1e7c3da7bb",
    "e76e0efeeec0bb93c47beee668f21ebc39ecbc076f8e4fefc0db03d97bf0e688f71ebc3934be5f7f37c7d0f7e0cdc1f63df87a527ec7dd1ebd13",
    "11c6f79ffef8af7fb8b73610cf54a07fadddac7f055f237f048870371c62f5c9c117a352702fbee1db9aefbccdb7e09dfffa16bcf388535076e3",
    "78fce99f6f03b27077b95abfa2751bb8df2d16e07e5999802f979f19f772999a812fb7b31978def6165c4d36494308518387b14f3ff1741fd046",
    "9f0e0902beba3868771f3e133e7efb7a51e0b60205b807f309504fc6f9f8f19f3cb5424b66028584d93d1394007302f7daf23387eb4db0bd8ba2",
    "7aaf66297773a0716ee2ea43b60d684a6e0cb3a8ab05c7130e39d210bf0c0bddb248702c172bce34448e2dfbdded733b75105c5bf21ce75025ca",
    "6a4bf1c3a0a564cbf28af58a9d7a4cc5d551eccca8d41215d7a018c2f13fbbfc3e4b9886c3aaa6db31879fbb855a45bdd006e057811ed98ae4a0",
    "3e2e6e2de4acb868071d38a852fd173dc90b8e6cc4d7d91acb7eb093452faf7f9fc50b53730a9d016b1a5c8a92c385b9c555b3e33997b21fdf55",
    "3bb8e05c4938084179a1480ef5e17cea4d71b8b80a4b375ac61ace65f7ad9591b3e25a1ae6dfa8b631f05e2ed15d8e6b7bc86ad472e8c9b9dce5",
    "9c82fde99c5c1b4b57b0f24a8b724d86e4670f19370ca4048ae6734a598b205c8c80ab332589e5e8dc5bb12cc5770a8f9f7e24c9cf062ec4c065",
    "9ebe7277a9d4eeb89bbc6c77d86b8ab73bea2e5d1ce162b68f5bf254961746b6c6a621b1d8b3559b96bda839575c92fcb86540fafefa252526f9",
    "384b78097bb244ce8565730d2bafc16c115c2bc18cfb1a6a515c1d50c0604ccea4b8d28a71f94bea4bc2456b4672f143b104f1b87e8e53511157",
    "c941592882280df7999642357291e5e559cdd7d2c29b4634d3e3967a51651f4f24ea8a4b96ea191510f118d1c54d4f15e752ed4d55bb34356b6d",
    "ec6afd045cc47705379537eb0217f291bd3576f97d989ac5f407595e897e38b6d4812839cb3b86f88221fb6535c3beaf54bdd47a1dc6b5a89650",
    "f276e2fb662a5a98256bef979c75f7d89870290153420fca6fbb38898de6b5f6a6e672c55c376e962c45f7a8d566418c2d2b6e5aefad510945cc",
    "026c296e88d554b2e4dcbaaab535852ddbce6ca13efa909c4d02244955bddc10cd324d64943834d74c00a554647905b9ec20f7c6ac92abd5b82f",
    "391c19964d3ed61164714be536c9d76694dd57b39b6f36caaa5e4c4e6c57bab86229e910f05de4e779ac77f6e4b0926bef77f0b76a487cfcf6f7",
    "37811873e00eccc4f165dda1e7b00e4a6d47f8bcac11df514825f2aeac39ddc5a1a1e6b2d520b345efa08a83dff9165a967e863b6a3456dbb250",
    "abc7240dcffc4820e47541e052aa91d786384538bf32eb739464d17813e74f9d42c3a2cd907ac3d28ede430df7bb17145575cd94e581c40cbc0b",
    "5a4a4f87c2a475f5cb39ef60102203e7f87a6a4c2641c424f87210cce653fc525f898c422289cf0ab39ac074c5a4c7cfbfbbb7162c3c6f20b6b4",
    "7f6277c6c73d769a32f7d48b61748fbdd859a4b4d36cbbc36eacc03becc6a8bcc36e6cd43becc6e42518b2cc7e9d88bcb46dbef2d2b6b94fb17d",
    "1d31ec58930c3bd6f73df6222b6869bbdca1d82ec368a59b3864d4215919b6f62b53e552acaace167cf3c82527817524a28162a82bddf5473a54",
    "bb36cab05db7650be1d094d9b23a146f861d7a3cc30eb3809686139d82a52f66118cae61c14523058525a412c848295c051590b2cc26784d421a",
    "8d12ed2428c73004d6cbb429992c890ac51a32bc2c2cf4a24a5be6c2b048aed095405d0a63c5cfa355555a36dff42a2d07f56d11cf52ac9f7b9b",
    "456d09d183a4eb904d585b810551963c55a5b94ddbc2d2a7ef92549a46bc22b85569c13c97fd63763aa4387140123c506d6481b56a51e8ddb2a4",
    "51aae24c1a589c822a48deb5a8d2e59755d87519e9a59a2a2de3901907a0830b9985a5110c8b76a04ab16514b91a5393c472b108ce51f9b60c2c",
    "fa71ef9a444d7d5b37c5a245be462dbdf830ac8626bb2df9a79530557f54e4bd839ad2d4806a052a3f7ef9875b9b635e6924a30792dd811688ac",
    "c02db49982e7f9a902b7135955f5749ff01b10bbc7ec6fc1a59a56cb70f29ecbd94ce3d14d3b23982d55338dfb344705e5c688e63118790ec5b5",
    "52b2732d2a2ed9f56b5c4cd65c48d83571dbc90436e19a473bc215601b3ac7c270876fb58d8e72c8a8687ee1d415b6fb3338b6bb472877b85b04",
    "b7b96f2877b88338b7bb9728b77bab38b63bbf08f7e24ca3dce19ce3dceeece3dcee3ca4dcae4b726cd74c297728ba943bf466ca1d6a389fcbbb",
    "b784702fae172eaa0e47ce3bb9fb6f91e4f76079fceeb77721e296bed9af0dac4efbca9d063dc14ee700c15e3c0db7dcabdbe21e7b7181dc632f",
    "ee947becc5d4a5d86e38f3d236339c61c74136c5f673718a6da7ec8c9a38be422cbf65efa3dd16a7df149d963a92624b727a697e04c7b0a5b25b",
    "a000623e3935bb87138c3e7961481438ac9d7df629b01c9b7d5af7d7ec18565ab24091116a1658c77b76565a54df368b0562cc281a0a85b7f8a7",
    "7531544bcbee5ee7c8b23433d910e79cbba02612d6a230858da5b27b3cd1d2a6abc06ab51919d58c1c50c5ed8e6e6ca5a9058f3341c4fee4aa30",
    "cbac675783b3a02a0e2c1170e42f3032acc3f8b6f8a53e94cc0a263e708d5f49ca4fae7ec619bce13eefdc08911afcde95e62c92a8bee7cebb3e",
    "8a533b437dfc9e64874ffe3461e9313e63db63709cdb9f96e3dcfe509d286f7bf64e7cdff31b7a1cdb1fe413dc76b424aa7d3ea952d8f3c197aa",
    "f6f9204d70cfe772a2daed988f73fbb1a1fabce76348c53d5901a2da679b424c81cd4411dcb3c5c3b0d380e2dc6e90716e37f038b71b8c9cdb0d",
    "50c61df62cc776f39870aff636054f0b9e82af51da1cdcc3be39f7f25604f9c4d7c72748a35f5fb3b8efc5d7d731d88a53d2b43dfef1c7bb30cf",
    "e8c139790d42eae50b77468d12ec2508957047482ba502dcad48bc522df68c97d6ab71c3c3a278793918d72c105e703d59d695d2cd0d4a38b366",
    "81e5d0ba6e06f6c9f491ec653bd10c73f5e365e52acbcb78533a473c2a21cb2bd1b9261bbb849e51250c596bc34e9091efb7cb41abc333c18562",
    "3614ed3c7b65015cacbc1588bdc7f530e372951c3c24198717a25e3c9568ceebb0ec58b3f579791317691005d686e2dada52f173b1fc37824bd0",
    "82d63aec557e5f83cb25e3e19326b93abd1db185f2861bcebda1fa7baa20bf799a784c49f6dd5a178b2a63c8d252ec565a9a31cad2b0ac1717eb",
    "9423d1edb59089e3efa8665eac6b1da23c7b67837303a91bc1598221b532d60c49b84bd2956c2cf839e1f42465256a13ee4b274ce83c14b784a3",
    "733137c52dc3b1daa0a5a8cbc356b5b8d99414cd1fcdb28f99eba569aed802ca21c8ada5e3ce24b8d6e4f779fa0bdc8cac55714b1336ae47b507",
    "65e8842ea8fa94eda8b8bd08c1179b6c47830b0a82d4d28788f6dab9d89a36bdca7eb6606088f934653bcc0b85f40325897adf6ef5fdf187ff76",
    "97a660e2743498a88e393e737bde03c16d891418b7e76560d499e681737bde8877dc67228a77dc176b97702fd633e15eac71c2bd58f76fb82fde",
    "8237dc17efc32df7eacc60d8e11a61d8ee6761d4e1b4b9c75e5c400cdbfd498c3a9c53f7d88bab8b61bbdf8c5187138e61874bef1e3bfd838c3a",
    "9c8d0c3b5c970c3b1ca10c3bdcaaf7d88b939661bbc7975187fb98618733fa1e7b716d336cf793536abfdbc6b0e3aa1cc38e9b7714dbc5fe3df6",
    "b28b30ecd8941876ec710c3bb64c861d3b30c38e0d9d96b6eb07f4db767583b674d75e6869bb32c4b14db7a2a3b0ab6a744c77cd8f96b62b926c",
    "22799545cd369ca442671eb8e6a56496192f0bcba1672501d777f58f91a52c4585eb4fcbb32eb03537124e5c359601a528217c158eaa4547f48f",
    "253e0165be194f2470c2c737acd129f6c5901c138b64facc4af067cbb26af2fee2ebd28e551116b1c62706a5142114b20e70cd1288716e6d9f36",
    "c5c486079da8c1119bcc922b4ac7b24998cc372075b608e9bd3e6f36a51361139866c125c94119724baff6a4b00e9181552cd69da57ac0cbc6d9",
    "9eb78c8a6bd955d9e0162155b4fd926cfcc870a00a4e2bf8e3f1ab3f3c7ee84b9e62d7c62516a47bf8a1e36ee7b4a8905646fbc225a4d1f12c0f",
    "101f0c8b7873cf8b8b70bc700e5b07424f2cd59628ae7a9a28555842ba630b78f191e5dce87e6532a5206a5d02ba748bec9855d50b6e64e7ca37",
    "965786e4ec2e24721520f2447078dc120139265839577af5f65a800fe76acd16ed912c0392284fcf95f9f8f9c7bb49709967762f37cf2fdc39a7",
    "08f6522de5d27592517255b51668f9e0123c488c8bc8a4e107374362d19eec4ed8ae45b56b503f4f5006d623e706c28ec099e926b8e45cceb2da",
    "328be585193e07281786df4eeeb6f5732ec25d6fb16fbcda36d14acfef86f018c1cde653d4b2e1f35a5bf6bc72a5573d18f6183b5ef854adc5a0",
    "f5e913b737c5e5e09265098e22b9ee17e883450e724eaf8c1c1ebffce67e12987b0547bf8b7fe25e26d52d7733490977cef95bec660911ae426f",
    "b643c8307479c39a315ba89243a8360eb6edf931c141a25812855a1497b161e0674b5ac3b9d4829d6a5a4a7fc1c5eec920861db3702e202978b1",
    "7867d92fab199e3f32b426b960f92873286a7007f43c6076b58573755ca7886d48acf8c7f9bd1bce35cb00583cb6429517b2650a0cf62680e26a",
    "746ec8395f2e2e598655ced901246239c29bb556fc687be92c597223daa4cac8a12ab83ebb4fd238bbe4aafdaa5610b4cbeac7e9789554706f04",
    "4b7cfcd34faf8a17fa7cb8984e253d51bb1647b117a947b81729cab8532a13ee45ca33eedc35583b54a7a4c71fff70271ca1e6c42b67ebd2923e",
    "b9535de158f64089ead3897219694ac199378673052eeec5354b7f23ca8319079967e1d7a2bce459dd42fc54be4979c9655eb25b2aaa1dfe84ac",
    "4d06d90e7ba6b7f5264bc35528e4fe0d4573a5f74bf864cebd081fc16dc28c728770e4e5edc296739beca6d8b11510eed85838b56f539cdbb73d",
    "ca1ddb28e58e6d9972c7364fb9436de0dcae86706e576b2877a849a23c2959f2e3c77f7cfc3006ee9d586a3ff3108c0675ce5e225a73237ec1f0",
    "f4eafce422e72c68d20e5996593a1597afb7967a816dc0b98814437818094f18500ef5350b95b1eb241c1b7358c04faaba1901e9bf117963978e",
    "c5e75d279efea0b82aaf1a577b97d80cc5436f70495555dbcde99e2c5da428cf32b85a62415d2fd4998c444caa979161b45a50c81c652aaec5ee",
    "e5a52ecb6bd9dc5dc39e1ae258f1c01b7f854960b55b008c853909ae424459e04d11bd8787ac3c8028db255fc195cf636f595cead7294e545498",
    "d58e5947680ad3cbb63cfee36d53a35f41811fb7b6f1c41d5dc7b96d2428b68f2bc7f66942b963da11ee651a13ae4037b3eff3bd80713807f15e",
    "6916a1c538dc016c1e36519bc2a21d53aff964d1c29c0b764dc802c88ae6fc40ca92b5312ce2f10f1bb4613b1ae7bac533e016ffc88a2b571c5c",
    "47dc85e072f408b2800bfa945ba3d5ddd71ecd9bc4b98a2080f5734130a6e4f0a3e54c1254a8addae9c890d4b043a2a578492cc21f8fd39dd824",
    "07357961c96e7273ac598f549c5bca118bb86f6fe5c5284702765bc2b1811eb11e9df3932c31dfb39f9de9199070ec61474b908f021b88c9c7b1",
    "de50dd929010c2ce082b0e6205e7417f51adec8e64fa382fcdad2949e1970ec3652a73c9136da24c1cf727c525c4a4821b4d96d77bf6635fc4a9",
    "530e97fd8371250429dfdfec03f5f19f7ebc93a0b09cba4752d7d29fb95db1a0dca1a8f0f276c587739b1e25aaddd43251dca6e611ee456da4dc",
    "ae85f2cfdb955aca1d4a32e50ea59b37c3466d71bd492e5fcfceaab1cd4b3ab5abb4a9b0888f824283180c8ee5e2b9b97b68f55b6aad45cd9425",
    "0482e7075e5b78555c8797db62c772545c9d1eaab3fa78282e590e9c7c654ce15cc445257b4eb588a987f09bcfdce15172c3c3b41282ca0556f1",
    "4813827f526c8a2bd967e8cc8aca78cd06f72746aab2b478bd5a1bc6545c4a7ebf63097ad9d88888685cb4884963f66b4c517e9db9244c0aa865",
    "16e007f597ea7517e7dcfcfe8425e9179db2eb51943bf4323e03763d8fcebc5d6d14d8a685f289bc6bb57c616c4a325ddf87ce2d64d9a6c37311",
    "af77b4f6f8f9f72f17907075200edb19d6986fd496a19f737bc67fc16d2f0830ee7c918073fb0b079cdb5f4ce0dcfe0203e7f6171d04b7bd10c1",
    "b9fdc509c16d2f58100e374b6ca30717793fbf98d514db8d748a6d163fa576f701c7766f046daa6f2c486ed914f67cbf4752e39ba8a79b479cdb",
    "6f32096ebb1925b8eda615e7f69b5b9cdb6f8289f2b69b6582db6eaaf159bcdf7ce3d262bf49c7b9fd669e286fbbe9c7b9e77b839cda6f210a19",
    "badd6ae4dc7e4b52c8c6e74b97a2b8ed0ea7e0b63ba1a219db1d53d58ce73baba23cb501f5c7ef7f75a383e02cf7ba8650fbf882bda83494db35",
    "24813deb5b143bd437ca1dea20e50ef59297f7acac52ea507d29b76bd2143b1473c2bd28fa943b0c07ca1d8608e50ec3867287a144b817c38b72",
    "872147b9dd2e7c5bad5b99bc1586d959abecbd02870f5ec61a21cac9128b3ddc157bec72a1152fceb55ab11e877d1d1e7ba41842b4ba3d3e1643",
    "2e8aebf600f7ea947758763505996519d7714fd5a7c0088ab30b54feb39af2abbcd05cab2dd00b04977086ba7e8e1636c0db912d8900d2714cc9",
    "69e1381e3ffde34d3e5ea84c76676dae25dd9eb1edd50f8aed8f88506c7f938460e713271cdb5e4ce1dfb63dc0c2bfedf939175ae7fe380cc5b6",
    "a76628b53f5c43b1fd1d1c8aedcfeaf0d29edee8e1dffffcde8f68e5f3eb4114db1f23e25dbbbd6d748f995e1ce25b6c2dd3eebd5143a6958e0e",
    "e3c5bfcdb2b97c03d6a6c00ae4b861559516b2c7b01684a50bec8a1cce814fc96179a3acb4280b334dadd88ead5a90a0bf5aa4581802cb4899fd",
    "1eabc9dfd87c8b455b7b5dac03c37c55bdc5fa25b24294d8b03cdd196f9373ac154b341ea3ecb66a09c44309f2c39081c19274cfaab068a9bc7b",
    "aef2f3e1985e58f3c30cfaf93d5a9d350c3587962d67a98dd30c4d2caa108765235fda635718368bbc4a8d539686f850cbab5d051571496a51cd",
    "8f1f6961d051327a2faa3a53f1d24a1e0a5bb6a5e5f2f6535b4a217f17d26ae72645516b963f3c54318706148f61a515be690c4b180f2a85a0a8",
    "58bc3ba25d15e1d896f19b627bfe708aedd9c829b6e736a7d89e299d97f69c769def1acf29dc29b5e783a73bd5965c9e6fb45ba67aae046c79ef",
    "39b665d117bac2734e7ea52b3c65f8e7956eef05f00d7e7b7d406cf0cf6f19706c7b198162fb3b0b14db5f6de0dadff60604c5b60725b8d2bcbd",
    "4e41b1fdad0b5e9ad4d4e7e3a75fee846ec3eb2ff6ad6329ddcfdc2ec309f7b22550eed86128776c589cdbb63f8eed9b29e79e7766de886d97e7",
    "85ed3a83e0361544709b4623b84d4112dca66f51ee50df28776883846b9e65e31bb891aebceddfca59220dc6e1e1faabde1ac578984d9fdf9767",
    "4719a624b7debbe6a6dd4f5d3bd3949c457ec0f73e757978920af3cfdd459c1b2ecc429d45717dfa9e326a4db2bc926d53e923ca7e0ef15a945d",
    "2db68a932393196deaf2d66a4c569ec5650b2e557fd62ccb6a6bf36daab4a4b0666f8840a2296ae0fa08a83ae41cf50c22d015c254dc5ae2d677",
    "a1b4a6b89c7c138d4d8dad653c6a5eef94f5b6e84acaa8a3ebb5e6aad1ec4dcae58024b5b97b9e198ec51c4d3fed96a6857396e917d56aac21a5",
    "2d6a4d4d8af935e1bcb8a0779799aecfb35c44629794bb69098fdffef3cdbb34988b1009506ee77ce6b6576e38b6bf9923b8ed091ec1edfa3de7",
    "76738172bbf141b017538672876544b9c3d0e2dc66b651ec30022977d894943b4c545ede66f0726c9f7b847b99cb9cdbd706e58eb546b97de952",
    "ec900494db050bc15ee414e50eb947b9438e52ee90cb943be43ce58e7d83702ffb10e5b65d8d52fb1649b17dc3e5d8be7f53eed0072877e81784",
    "7bd1572877e83fbcbc5d9fe2dcae9f516ea9a3661f048f9da05c0efe92d81c6d2a2ee431dd4c17cdcd58b06e5a798036e55a77e7c09a79922bd9",
    "3faff9491de5d6d4f4d922870d79af3fcb53822ff97081b36c7c9cebdd7d214b4e26c55d9f576295c555eb949cc29bc2427b72e750ae588a48cb",
    "8d28cb2b280772ca720a722ed9bbaab0b36b575cbcb85a55274738adecfbaa92ca69ed7dd5f748b93470f6934c8ca63c14d62c00a423c8227e4f",
    "ae63aa43c328417223567f6e289721a772cfa6c88d949ae486bf70bccad5e5590c62c375eba2b8d8877176c54d703907cbc76a4b5d71a3b8425a",
    "5bd522239ba259e794922a577f4aaab6ae047dcad1b016a6deaee2f45abbdcc3c3f4b79ad7d63ba44ed0bb95d7b2d694526c5e9e5247d7922de9",
    "52e625d55d47d75d0c2dd31fc46ab548e536b479a5daedd2d4786392c4c72f3741251d9e348be258e22b7da1f66b10144a763504e943a150532c",
    "e09995f5736a516216e4899b1ffcc390c9d47e2d96ba90523e7f31a15a17d8b2041c43a2578e25bb005e3fec1d3e4a45b89a3332375b6a847b2c",
    "2df53d589d0317401556ede731661798b955edb1d2aa4aab335bacf0d4df96529efedc7b5595067b2f0def314241211882f2ddfa4d30cf38d62c",
    "1a176974c5b7c1095a3c457413230f1fdef465306a534db0e39786ecb3f56f2cad4f3cdfe75849f4db7a47e48abf589f916b9f60c861ec8f02ae",
    "1d9df65bafc8abec580d4d60a30d17d135a6bf1d5bfb8baf98cafb6d3521e341168baae292a14334f5e60feff13a61187bd01f2c324eed015a74",
    "a8f6782f32f067f898c09e82d128b587b6d1c9b605ca71ea39ea8e527b0c1f5d067b482017465b84219580726b498fdffde6ceaa820c8bfe40e8",
    "5a45cfdcee45e7dcee9517dce6e5a7dc6796cc36cc8a241c6eeb7c0b6539c4c1d99b9c9c8bc58e84a6abfb1cc325138be799920bb85865c16316",
    "794d388458b7eb8844755e721bd33aaf46c5d939f53770096690a5511ce2fbd01fc583e0c6988aebcd8fc0962622cbb340e8afd172944b982408",
    "850bb25ac40f832a5db416972a92d59ae6ac8aab96e56859571e2fcd3984523837bf917b53deb09f5b89ba1dd58faeeb0c4372f06717440a66cd",
    "b5684ffea66c5629e53afcc56b9bed6683536ca038480c386615874987e3fc5a921add58837125b528274b6b76329fbbaa17774d1c2b536235a5",
    "e4c5c526e77cb15e89736a0cc1f2f04ec5295b3ba23944a6e569162b32b8af708c243b65460f2498334b0984784efbba911416ed65079c31c4a4",
    "054be84fce29ca65849dc24b5482fcbcdd9944b1c337c5e5edeeebe2f27bf79d51eef0c5f16d68f7ed51eef015f2cd6ff73dd2cd74776552ecf0",
    "8c52eef0b4f2bd5eeb18f9f1875feeee50e2c9fa68aa1a5e0679e2a2bdf480e40a4562c52d5dbffd7c8f21de0c9b1994e8cf14568cab6e2af610",
    "a7e23a365bfc6ccf88736e59c5f6f32c926a69ba4999a3e22ab26fc3bc6b5d62c9fd2125e86a0fff05e15efc219cdbfd2b847bf1d770eed9fb43",
    "a9dd9544b017cf14e50e4f17e14ec719c50e3f1ce14eb71ec15ebc84943bbc8e9cdbbd987f7379dff67db8a6d69cbb423508b6fb0c18773a0d38",
    "b77b3444799b838473bbbf45b46373df706ef706716e772e716ef755716e777d29eed993c6a7c1ee98e3dce6e6e3d8ee347c37fb3e7d909cdb3c",
    "9a7c4d6eee51816dce56ceedbe5b210a9ef3e150eac8ae4305d04b8a4c22f75e526e12217fa602247bcb910281eedf6fd481f2f8e38f3782aa23",
    "b424470ff108f989db8f6f38b71f07096e3b5e12dc765cc5b8f3f88b73fb711ae7f6e3b9efc76da7911cdb0f3719779c950acc22fd3b5ebb8a8a",
    "1b7958a44b9a45563ba27f5eea55d6db5ab6cf8bf3f2f913aee4ece599b2ccb95cae76c859f5f7e72cc9cefabc72f9756fb1b00cc99aafe22437",
    "90bfdeee97cc28b92be82d0f355b1697fda4bee0e130c1f57495a6a0b5c54f9f7a497e5cc5db51c6b5efca150fc6b207d505b6768de4a1715362",
    "e67258d8a7f79971f6f62392a186a6abb591e869c4ef87c5da3da22308e183c39f562c122b64d177e0f0a4002249ba903e8b1bd14de6e489dd29",
    "d7f14e2022586a979cbdaa642e82d2bf2757aff696d9657b97f8b4fe2b62afc285f9e417a1ec9d538ec5ee7ea25aabfc3a3fab181f4d8e59b0c8",
    "57b82fc45ce97853cddd3ffeb8e0f7e3c64cc323e314d52f67d2d2aa92e25af2f8b9a59d75c5d5e2fed36039550587dc3fe0def45e499727d32e",
    "b9712ee3453c77a0ca76e4e2d77f4a103bc6df9db3505673f0c63a2596afd6aa990cf994ac97734a45716dba7f37c796145711418b9b7363c87a",
    "8b2577c18bf4b97c4faee2f513f44b9ab2979144c25d7b8acae1f27a1735973b6ee15d771dcbdf010bc8fc80755b843a00ce8ef8b16e755bbf99",
    "b36c5b0337c564af84e99753dbe851712979735b2cb2bdd9c676ad901965bd297a1465cebade50fcfbf208babdc17783046720e5966133ae5bb6",
    "7976c565cb338998f93e15972cebda2a177b16e5f0887334d7ce983d2bee25eae23b712fb123acffce5014321e2f912d4c073a0365988973c6dd",
    "308bee0ce321dce6e1e376f3ee30e4dcee80149cf613d4cfe1c8f6e6372ee9d6a238bf8091629b9c6a5044ede774dfda2fa5f5e4d7d973ff362a",
    "2aac05bfe51dfe1eed7cc3b5c79f7eba528dfa0338355cf923ab3dee844791c27ce66ae9fe9653d158b7f76a901e557296ef257f3a3c29d7f0b6",
    "833d3b5365bd036f74e35d203fb7a6dcb4278b72bb2ec8df525813c91e862a3935c5e5eb09a455e9549c3fe554f0cf925b8db0d6f62eabedb8fc",
    "55e0b36b5173abd3f02c5deb926b038fa1c1cd3a8ae2067254151ca655d97d08c029e67cd6bd02eb6b71ae2b736e6da46b6c7192d3b2e430aac6",
    "d524398c5631d123bb6519adc6f518e55c9e57bd3d4dc90dc4e134c4ebc829daa3377724350b22dee7cac6f5269750c1a3ee8b9b57424ab6221b",
    "16788783571687138f026d4eb43558ce346029a6a1389ba15035675458194ea9e563993ea6553a659d115e06c3126f839d1d77eb91349aac3644",
    "ffbc2b8a84733d16efbaa8bf0f5e151b09dddc88b8390c6c6eb2f3123c5298286524c94d9fee2de85eaed5e667eb5956eb39daecb8574e9502c5",
    "b6344f162dcb6bf67d55cef760af14039b7ad472ce5edcccf2f322b288827bd32d269d16674fd2f34995e0872f781326096ee00d772b2e209249",
    "706becbb73a1292e7a2b42e8536226f32aee0b72216547febeb5f42ed6dac0155613f12d8ba5619101ddaaad7526c9c178c403aab375c5adb969",
    "3b5a1e5d721da990b143263159eca929ab367ac65981256b468ca12a6e9902fe64635393057e73df70d5be8c880938a7d08aa10723c2fdbfb8f2",
    "99b3ea9683c992be89eb388dfe06ae2015334623e5a4b8d4ae76c8d1401c4cb45915b32eaf4230a29b8b9040f6fa59b1272a47694571030e547b",
    "1933eafe8ba3bb9a19e5ec0b39d963964196d73d540b4ab59f6ff3f1f557835bed5363f57abe35cad917116887d785536b7a5a41926624f7d2cd",
    "40801dde941c29ea6e39923cb2d571e68c6442e34c41f9962b59c9c897049944d4bfe4db643ad099be93703dfb33bead8e371a5ab1c0ebfef9ea",
    "02e162f23790a7d4d19025d64637a4a4299fa44be3d35f3772f7974f65a5e99a54a534cd69fbb53ffef4c79bcbbbab33a3bbfd5b49e90bf67217",
    "9872c7dd62ca1d779539b7df7da6dc71979a73cf37b34569cf91b9143b027d79a7ec81c3bc93f74064caed71cd047b0993e6dc1e764db93d8a9b",
    "62475038e7b618738a1d21eb94db23e0297604d453ee08d027dc4bc03fe7f60b0494dbef2350ecb8dec0b9fdba04e7f6eb17943bae73680ecdf5",
    "5b2ef75c474ace5ebf891bedeb410fc796e263c5e59086e262f3e3b4a296901d96669b7ca5f6a9b8913c75646da929aefa7955b573378e156fc5",
    "32d5aac212f4a46c170f452b90d6e7aa76aa5e8130f6f3af52e654dc7e7c28ca7b3e8d14d8f3d9266fed7e54fa1db9ed68580ddaf351b39803db",
    "d1b5e0b6a370ceed47eb9cdb8feac55c7e3af8572be3398ce07b7167d884e09e8230841cd8423a38b7458808f1f31c6fc2b13d7c85737b38ccf7",
    "e3f6f01fceede1449cdbc39304b7853bf1cd650f9fe2dc168df5ddb02df04c705b209baaf6292e8e635b949dd897b7a0bdefc7ed418a42fdf91a",
    "f12874a92d7c52705b38a6e0b6f04eceede1a242c77c8e3ee5d81eccfafdb83d7897737b30b0b20b9e838b05b7052b0b7b640b7e169c345dc7e3",
    "d7ff6027bea1783abb82fd143f94e8912ff606d857ae47df1f4bb15031ced9ad0aa8ac4869a7caf3a0a3112cc88a722d99f6334a97d85253bcda",
    "5965b5cb423775a5d9c387824bdd0d1ccfb043b91c2ff5a25ae002e582c5a5406d90d5be198df9f8f5afeec2933a4e1d5d419818dd7f3bb7874f",
    "516e8fc612c55db6f5b06813c2c1e5198d1b7dfc0dcda8e1f1f3afefc217f01fba3320238becbf9d7b8eada0d01ea7f10d18a71a8c169b4fd9dc",
    "039cbb6c87aa9865739ba2d23d950a2fcbf6ee011f5b56dc11d345b8971831ceed3167943b62d82877c4c4f1f2f6183bfe7d7bcc1ee7f61840d1",
    "de7316937e7e5915a4bc975576cfbd2cdabf093be23b29b7458b726a8f3da5dc11cbfaddb823769772472c30e58ed862caedce06ceedce0bc69d",
    "ce10ce6dbe158e6d9e1a8eed8e1fceed8e24c69d8e29ceed8e2ece6d7e33816d6e38d12b182c7071be1b8ccbfd17f5e44b393b57e5a6d19abf7f",
    "a2bfef4537e3dcaeeb51eed05638573d67450925fe3db6e5375c7cfcf29fef94b8f251e04641328ad5b2676e570a29772899944bc1bde3c39276",
    "512c267761cd32457178acd4bdd961d89d14ca75ac22bc6f635e6a5e2f92c5c3e92dabadfedeb83daa62bd2c5a6b6f4cd6d1642747e4c840b6b0",
    "d275b523d9cfa34ddd7b19f182f62ae554bd327bb50455c1823a442f077f7b255bf279ceb56ec5e5aae640f6e3517c9e674de45cf5e6f61ae4e7",
    "253c238e51ebb9c8f25ab2e1281607c8cb8b8827c45c0eb6c679f7854b2978838de0dbbd6f1c6a69f82b5265c849f566e9a6c72fff6ccd2d4879",
    "57ecd0cfedbeb5cfd9c9ee12c1e389abc8b05e2c82782a6e4d4e4b079782ddb1a35c9b9e5dce5f68139c3dad1ae07e8e92b350263cc55a86e646",
    "beb8a6b86ecfc3ce1c45ad489f16adb54bf513bd87a477c527b36c2b3226e1f04b4213c7a8a052970371ac0bceedeb8c70e7b2a5d8c0c355f62c",
    "6d2b9c4306fd6069fe960ced92b388b78813c124b932e6d7979538173d95e2f0ebc214f310f3f831d56c6f70655b0c49b0f7451837b06a2d7625",
    "70682202347e33b4f4fe28b96ec13c6b96bfc33eb9a91a30f1844641d6c5dc2417fdb1de25c592aa37582854ea8ab10ca08845696a454c379711",
    "8b927a52df3660708173339d721dd95910d09594c8415c9a5363286ad9f0033fafed424ea5b5fd5b9fcc5464ad760521f8cd7e0a216dab45b835",
    "b918721c165937dce8e28bba350bb05cd3bc6801d6bb735aa8e776d51bb410b6082d881e496544305be8d14c9283ff1d97552cef97a815170aac",
    "93a7dc4a32628e31b2f610ba6a456c7e8727ebdeb3102a040856c999fb085cd33b40f2b8b525a3e4cefe4603c88ffff2e34d6c064c76e4453345",
    "33c48d7b4eddf94ddc549c3d7e04ae5d9e813bee45e2716c17a06fb9c0a15dae13e8659320dcb9e7506cdfc128766c889cdb3758ceed1b36e15e",
    "1400caedfa042f6eb72a28775829bcf776ab873777b7a27879bb5546b9ddc8e3d5ee36a318b5671394f6ca61d1f25ede2d64dedadde216ad7d36",
    "e02976f80328b79f6151ec3812e3dce9c6210bf7c52d74cbddb89998b438dd56847b71831111f5e256a392f6395aee0dc7257c79fc8c4336fcbb",
    "2b52c152a83618ad7e9f2396f184758f77c7cce30cf2e85bf8329619c526e62562e26dcbe358f19dcb354682414614bbdf910394545a5ab35b11",
    "9ee3904239d89d8d89cb689c0ad86c70290bb7d709b67e4d38a02970ae775a5a4388d575f36d148525bcb988fb7b66ea120caff105bb1f176c7f",
    "15985fc6e921a9d262be7e850aca30a42128dd35403eee0da13ac96f0ea7aafa2df8a5663502112a2e9430e86a7cc8e5ccae8f9f6faf28aeddc7",
    "bc0f96ac383f714b31b5812acd22a228d71d6b16e34ba9a53f5aad235a0420e56a4ac62d8bbe2a0eb76fe003b698574e21e5af7155d69af08881",
    "71164f48b988c408e04a17e58da5cd43bb1e3882107ddc716da4786b2d809b7229f8f5ce618e33c6358f3f5fdcf59e37e55a4f565ef354c894b3",
    "6b3205ef7ab721b935998ccb53726b555b7b8b054c702e9425b7b0faed810cc6e1ca616b7ecb5773cba0b77a935f6ea05ca9c3ea8db21db81def",
    "524736038a4eb00b99aab082cbd8265fa76c04b26bb7ebb26deb8aeb302521afedc49c737e971e57de467fc3a56fe15af03bf2c52f6e522e67bf",
    "57ba0c64d92fb124eb971a53925cf2cbeac59f3ee51c422d71dbb687a9b830fd0e7f3253967176b5ceee48c7a80469bef45824fe9558f3d18d29",
    "7589459f2def6a2d88bfc0e63e8b2caf4ebf7f3ade9467763b02592415fb8c7e7957b502f7c73ed347f4a9b8de5c452966a170aed660adc8efb8",
    "eb3a6bb2c348ce95516deb0c4ade26dc34334d60c62a9b5b82fbb42c7a5870d9f22e2039481a924b9ecd23fb2d50ca59568035f367d2d8e5ab5a",
    "db64535c9c907b788620c96e597b7d7bba2d2aca73e5a7f422bfaf5cd7092d944960e6fb4a52dc268859732d95acf602e4af766f70e8534f7998",
    "95e09a1267e7534362dd6e2f1771e9f3fc0e129565c7ab4a7c6fd95f69125cbf1e2e6a72c3ed88eefae2c0e11bfdc0c590b7e53d712d480563c6",
    "ab5ec9411170cfcc3b2ef63afce1a72c15b9887bcc78e04889bd7e196f09b94ba4be775c67e6dc7e3d9a7fde7edd9a73fbf56dceedd7c105f77c",
    "bb5c74cb76595d74cb76f95dcd82e7cbf462566d97f329775cf6a7dc913c80abcb7b3202ae2eefc90d78797bb204de8e3df902efbf3d99839856",
    "f59bb823d9043787b6dc151cdb536150ee48ad21b8e74c1d143b127ff0e2f64422943b1293706e4f74c2bf6fcb9bc28bdbd3b0501bf648ebc2b9",
    "2d4b0c3789f7a433943b92d870937d4f8a73cfbde4d8a1d891b28772470a20c13d6714e2d89ea0887247c22351de964089737b4226caedf99d38",
    "b6a78be2bdbca79fe2dc9ece8a0fee9e1e8b737bba2d3aa9def8c6dae3c7df3e7e8835e31a537f54386e267e88307e1f0dcf8484f9955bf606f2",
    "bc55bb8e572418e1b65ee0da8aa200915363cd946ab27708706d92482058d7d865188002ac71ed6a0b9ca66c081029012b1c82386c139fd85bb2",
    "460f64a553a0eec5ee5cc3533d7d3c5af0c3720e56bc755300a6999b042b7c2415d72f65c51181d51599094a91a05d9d4179252a107bd75aa0f8",
    "c298e3dfd2e6e1605f2a6c5e4bbd214603c9fa041871e3a199c5af4a5cc20673b021cc4d14d826ded7ace04abe9db39f20ee504cab79d9675582",
    "1defc52d305ac0bda8ba269498ae976f445b124c4c94386f47f02bd8204a1698fbedbafada8d6b00ad77325e069203737563096f4a7c3382d341",
    "8b235c6dc1d15afabb60c8eeb8665a78435df330b422419c7537bbe33b38e777f6c1895a718ba7ac0da7dadd812840bb345b1f15c74aa2ebf0da",
    "579f1d05167b228e821531a0bee44757edc0d650c1a52e6b5e26386ca48a4b7c5135257b3441456a1d1869145cd65e5ebf575c3b49aa2d78ea6c",
    "b5a122dfb4ec451cf101c44960158d4e784311581a551667765085b7200df9810117bd2bdecb1eba254b454e00cb94a3923e301be0e2953d5d33",
    "a47f7f576d41a4828d5d4c39ca591390670ba0e9fc621eb6e425a698929cd80354952b69608bc022c9aa5370049eaaed4e3d7525609002caf7cf",
    "308a5aec08d01db63ea7bd7d49c1eaafca35dcc1bf95e89f60bf76a7f5ada3881e2cabe7f0b26583c0ce920b7005e113fbd060ed61b82e50545b",
    "acfbacd16bc58b4facb8d13f4cbbc888c1e0a0e5ed2d00971d3d2438f07616aa2e55720547d20db90264814b9c2366b6ad159873d36d819dd510",
    "3bafa6043e1197e0aceaae1a8d49366c0083796728583c4b16bab12a81083d72d9b42e628b2e3165df2d62519b594574054ac48d3db5b42069c2",
    "3075bca4a4c10a0f67f5b3408ac1ea86a083fa27c76520534136913dd5ac5dbdd78283d58e0d19889d1b0f3f563cc4a8aac65448566096bdbd94",
    "b3947d5709497de2527c42c3ec2ed013558905f9d31ac0126489d89c9778ad7690a2ab76c51d027c16d998617d632ff689ae19f8b0c58d59655f",
    "570c6f85d7644409861e0c2c5181153a7631a5204435cc9838cd2d4d7fcf8082d76d5598a45d6215abb32271a092b3b5bb5714aa5a54bb2ed600",
    "8e2cd124a98d607bb4c532dee86a6b0d04ec3ed5f250e9c56279dbf18d3916094e3c898d4fe8498f74aebef786aed4bf66c1a4d668ec468a9b0d",
    "5b2fc0a8971f727741416d6a9c3167e11586a0905b41c3e1e9b5a3ea659a91fcc8f58d296777c1a3d7c98c86a61413fc1d021beb26bc03919b10",
    "fb5a577a3e4c7bc4312e5017d771ada9d91b965ade60ab5adcea43dde265d098cdac772018837589ed66c185e58df05c8b007f562973907fd44d",
    "e15aa29e1035bba3a048ad045b7df7c674cb4623c089f4b30db717929edc36c51a047c7a5735cc52ec08436f5619d11c00bb524b56f7946a9f58",
    "eb9b8da057af39675933fa1973b65c7e7b35797008894697a6e76c49fe8d6dc8ad1ceea88a9fc35092b641dbc09480df4c16d79189169e82a144",
    "3c5a80ec69cddcd95a42b4eabe87b7a08d2faad6c36727d703608c722e22b463da37f65cf5f87584f5222a34caa58a275f22b8863b215236cdb5",
    "532230b4c8c90d773b162092cc0d3d21eca26a437ada37e39231c3b0f5d62adbbc2685f98eb2deadd2da208381cb80d7dbdac4a95ab78393ac65",
    "04742b781f8b32a1ab399012c061995ea42ebd0610609266014e5556f7746cc1528fc08963373016092ee180606e7815a734349621ed2d910653",
    "8a06650d5953b52dd78b75707e63e502aa595aeb761d0cee55a526e130d2da271d914bfe25cca439a5c7a1c2878cb26ad23e8c609b1266a6dae6",
    "deb849ff3f6fae45b8736e0100"
  ), collapse = "")
  starts <- seq.int(1L,nchar(hex),by=2L)
  text <- rawToChar(memDecompress(as.raw(strtoi(substring(hex,starts,starts+1L),16L)),type="gzip"))
  utils::read.table(text=text,header=TRUE,sep="\t",stringsAsFactors=FALSE)
}

analyse_demographics <- function(root = "", adsl_file = NULL,
  output_dir = NULL, excel_fr = TRUE, strict = TRUE, show_diagram = interactive(), png_dpi = 180L,
  sites_file = "clinical_trial_sites.csv", template_file = "screening_map_template.html", p1_file = "", p2_file = "",
  audit_asset_dir = "audit_assets", audit_dpi = 300L) {
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
  writeLines(c("Moderna demographics, recruitment and BLA Figure 1, v1.16.0", "",
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
    "clinical_trial_sites.csv supplies coordinates for V03. screening_map_template.html is needed only by the interactive HTML map.",
    "audit_assets: V01_weekly_recruitment_by_phase.png, V02_randomization_analysis_populations.png and V03_recruitment_by_site.png; print/Substack layouts at PHMPT_AUDIT_DPI (default 300).",
    "A complete run-specific PNG set and its manifest are also retained under this results directory.",
    "The input path and SOURCE_MD5 identify the exact source analyzed in this run."),
    file.path(output_dir, "README.txt"), useBytes = TRUE)
  capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  audit_assets <- bla1_export_audit_assets(recruitment, nodes$counts, root, output_dir,
    asset_dir = audit_asset_dir, dpi = audit_dpi)
  audit_links <- paste0('<section><h2>Audit PNGs for print / Substack</h2><p>',
    paste0('<a href="audit_assets/', audit_assets$manifest$file, '">',
      audit_assets$manifest$visual, ' PNG</a>', collapse = ' | '),
    '</p><p>White backgrounds; 8.5-inch export width; ', audit_dpi,
    ' dpi. Stable copies are refreshed in ', bla1_escape(audit_assets$directory),
    '.</p></section>')
  report <- sub('</main>', paste0(audit_links, '</main>'), report, fixed = TRUE)
  writeLines(enc2utf8(report), report_file, useBytes = TRUE)
  files <- list.files(output_dir, full.names = TRUE, recursive = TRUE)
  bla1_csv(data.frame(file = substring(files, nchar(output_dir) + 2L), md5 = unname(tools::md5sum(files))), file.path(output_dir, "output_checksums.csv"))
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
  invisible(list(subjects = subjects, excluded = excluded, counts = nodes$counts, recruitment = recruitment, audit_assets = audit_assets, output_dir = output_dir))
}

# Embedded offline map assets. Leaflet 1.9.4 (BSD-2-Clause; full license in
# the generated HTML) and Natural Earth 5.1.2 1:110m (public domain).
# Sources: https://unpkg.com/leaflet@1.9.4/dist/
# https://github.com/nvkelso/natural-earth-vector/tree/v5.1.2/geojson
# Natural Earth countries and US states; coordinates rounded to 4 decimals.
# Gzip hexadecimal payloads use only base R to unpack; no runtime downloads.
bla1_map_asset <- function(name) {
  hex <- switch(name,
    leaflet_js = paste0(c(
      "1f8b0800000000000203ccbd6b5bdcc6b230fa3dbf62d0ce21123473013b713428f362c0366b8161018e93b07879c44cc3281ea489a4318c81ff7eaa",
      "aaef9206c8de6b9fe7ecbd62467dbf54575755d7a5b3d2dae7f1d58497ad5efbe7f6abd62419f2b4e0e1776f4f765aeb6bdb937856f0d6be48fdeebb",
      "ed6c3acf93eb71d9f287416bbddbebaead77d73758ebd76c928de637f3bcb5751d5f65e997246d2edcebb1d6f6249b8d0ee211ff6e6b32695191a295",
      "f382e75ff9a8fddd77c77c9414659e5cceca244b5b713a6ae12892b45564b37cc829e53249e37cdebacaf29b82b56e9372dcca72fa9bcdcaef6eb251",
      "72950c636c80b5e29cb7a63cbf49ca928f5ad33cfb9a8ce047398e4bf88743239349769ba4d7ad61968e12ac5450a51b5e86df7dd76bb7dc2115adec",
      "4a8d65988da0dcac286106650c63c406e3cbec2b66a90548b3129690415e527cd76ac132437968c3ee2e1d55c6023d0e277172c3735893f5fa18a02f",
      "6b11d4186076a3198ceb8961e00870247f77182d39bb51369cddf0b4a4d5c5c6a05207163f83ccbc7513973c4fe24961169a76876a5a1380499d7ed8",
      "3b699d1cbe3bfdbc75bcdb82df47c787bfeeedeceeb4defede3afdb0dbda3e3cfafd78effd87d3d687c3fd9ddde393d6d6c71d48fd787abcf7f6d3e9",
      "2124785b2750d3fb0e33b63efededafdede878f7e4a47578dcda3b38dadf83c6a0f5e3ad8fa77bbb27acb5f7717bffd3cedec7f7ac050db43e1e9eb6",
      "f6f70ef64ea1d8e921c34ebfab576b1dbe6b1dec1e6f7f80cfadb77bfb7ba7bfd340deed9d7ec4bede41675bada3ade3d3bded4ffb5bc7ada34fc747",
      "8727bb2d98d6773b7b27dbfb5b7b07bb3b6de81d7a6cedfebafbf1b475f2616b7fbf71963876678e6f7761905b6ff777bfa39e60963b7bc7bbdba738",
      "1df36b1b560ec6b7cf5a2747bbdb7bf863f7b75d98ccd6f1ef4cb679b2fbaf4f500832bfdbd93ad87a0f73f39f5912d893ed4fc7bb0738665887934f",
      "6f4f4ef74e3f9deeb6de1f1eeee0427f77b27bfcebdef6ee49bfb57f7842abf5e96497b576b64eb7a8636802960ab2e1f7db4f277bb4687b1f4f778f",
      "8f3f1d9dee1d7e0c60e69f61598ebfdbde82aa3bb4ba871f69aab04287c7bf63a3b806b4f8acf5f9c32ea41fe382d24a6de1129cc08a6d9f5ac5be83",
      "fe60014fad39b63eeebedfdf7bbffb717b17477388ad7cde3bd90d60aff64eb0c09ee8f6f316f4f989a68c7b04a3fa8e7e5a10cb68275b7bef5a5b3b",
      "bfeee1b06561d8fb933d0927b464db1f5a62b901e8573adf75565aff672a11de77ad0a0e0674d5fac7099ccdcb5c9e6d387770a4e26199c069be89a7",
      "45bb352ecb6911763a1351f3cfa23dcc6eb029072fb77e9dc4a3e426b1b0327391b1858b5b30b0a5ab593ac4c3e9978c07f75e76f9271f965e1495f3",
      "290744c1efa6595e16cbcbde2c1df1ab24e5236f496502ce9d4df880fbb254107aaa39d382a8b5bc2cfeb6e39bd140fcf4cf3c59cf3b87be43eefb65",
      "d4d4cdf524bb8c27a780bf06e667583e3c147c7215b4e58a44f78fc1a38f588e9939c18cf022411c04b3eaab8cd604b360a1fdaf71dee22c6169d463",
      "5914e7d784e70a6834bd2ec7fd7433eba7abab0116a5fb283165ced2f3a03ce3e75102fff4e12e98e5802d1fb1c5e3e890d6b13dcc3920c787077b44",
      "b2e4c73660cb32c33946254bf96debe3a319e0473fb87fd45fb1d81d6c1a87ba95e7f1dc546f177889eb11b4e18e180dc49f763c9d4ee6503b6d0fe3",
      "c9c4d78367bd2008fd24aaa5af0766f9f4584bd90eac945c994102f097c275eb575b8076cd6fd8121cf44ed435531b9b45f02ee4e65d24230fef9a87",
      "07bf6c5b89d1eaea4ec09c24b32a7fe2aab044ac4bca325644d6d0d36809b67479d9cfe5e0139605b0c74bbde091e54ec981b5f3a15feaf2661eace0",
      "e529dc884069f8056c056cc2523778546b9e9b417d700615f1b31e0037fce99e43dfe99ad9a6288ad2e5e564508600f76b3cf87fb2d50cffe1a6ad99",
      "de81a59e494d0430a80c688763233c3a88cb717b9addfabd2efb9a25a35697f27e0c61bc94976770bafc7285071d1e98063f5b5059b6e1a8dc0cc41f",
      "3f08cb3610189378c8fdcefffd77b1fa00ff7ddfb9669e67d57f67d5c7b6dac57492947e07ca76ac6243316e75ea123c4ef29818501ec7c5e16d7a94",
      "6740bd9573015a25f3b229112b5e400022bf22fd6b706c528310100174a413ce927358ffc41c51956186f6a90249d1d9791fc79912011464ede9ac18",
      "fb3c45d2ef135c39d9cd344b0130fc6490b6cbecd31446bb1d171c162c0d56bdc85b6d285b22c208e4287cc0886bbd25d89f361c537e7778e57b032f",
      "1878cb5e883f56b3f69f5992faf02d8ed0af51e7dff7ad15ffecdfb717adb5f3d5a0b5f2efc7ceb539567fc11c126b1fd5befdca5c14cf05be4aae7c",
      "03234139ceb3db1622a1dd3c87897b1fb3d6d77832e386a4c36b090692c49713def2564b3593168f1ab03ede0a09e07489004612672505fd6d4288de",
      "99b87b5a54e21c9a8a6ad05166270097e9b5848bc04298ef2bd0051827d92c151e4f0087c38c4b84059caf1c7a22e7b0d6a351fe33f246711987c94d",
      "7ccd3bd7c955ff1236f5c757ecb83b797fb833196ffd6bebedd6d64e676bfb768bfe8fbeb7b6b7768ac8ba5d7eb70e0450d4a3ecf6ccbbe5975f9212",
      "d6edfce141a5dd64df2a09057ed3587eb351e63fb041baaaa255dca41db85458220efc4d7ce77759efc7359faffda6e10beaf3d5848986db16f64220",
      "a10ebe8f6466ceff9af1a2dc4a61dad8dbbb3cbe812beb77df3b6ecac133f80ff62f557b18a7433e69a8bcdd90817575d6c2d66dd850bd00fecff514",
      "ec7dbf534717b677297978f81eced43fd4fe7e2f004534c2c42d1af44b916aa3c01cfb2a9797ffe55428c54a95405b5c5c10145e5c84e96c3261fcae",
      "e4e9289c3071c187c70c6fdb3066d740d44de2a2dc1b99db73e79115657c330dc70ccf59594e78f827bbcde3e9c7d94df8815d01ebc4dfa5e18c215b",
      "1797989a3044c0e16746b8f473968f8af01dde428702778543ecea288635134722fcc44a7e03471e46f31793e72c1c31895dc2f70c72cbf91e82f6a7",
      "7c12fe93c97d878ebf676217e1e7bf98050eb423e11d337b2c52726bfd7989940a2fdb6249227bf34c21283324da0c2e73f817701e709ff124f90668",
      "b092a0481624e4ac0b984ae1e6ec41c90f59f6a5f0c5ee2041c4db1717c50c90f0c545440535ca005c7eeca7016173bc71fc0ce916a00767c332cb01",
      "13585458469d042fbb937078490074b12e07046b9440821e8a4cf0f1f2a161c1df00116fd92e90951e02590d60c8f427523a493a9ccc46bc10e71def",
      "389582359bc8e3fde565fc5ffb20b94bd2e0be88467e110c8af0ac109718b6930336c93735599b033a2cce72c485b25e7bf72baef3f232ae4e36e1ed",
      "db182e296f8703d702541ee07f390c941eb87542215a98caf5018c3799b42e792be737d957aad8ba9ac141c01438c6052f586b4a3f206bccf3a46c5d",
      "e5d90db44aed518da2e4f1a8ed315f5f49012ed2f04bf038910042e7f02c3b5774a8b574c1a33c79131ff6346023e8b78491ab753609aa0acb145100",
      "fba57f03756c088cd42630d8c4946286f400fab27d912800052a02be1d98b5c94e445804aca6c63694e523a048dd5acbcb95040383f260549b00ea54",
      "ef7d19758102b507a6c0a0dce4fd1240c1ce3b2bcfadd61f1f1987ffe9b5758eb7b890dcf3a696a2af77a0721ecbc0ac17ac7573ed888b69ddf0fc9a",
      "4b8ce7dbcb8cb9342ebb44d4c0654d16f440e3508dc4a3915ada860936725b75760aeee3063aa804d2dcda73c202b0c2168ee3868da80cd582a48539",
      "0f0f00640b7305dd9ac8a9f6c57ceeb334b449427979d658ff32501024c960016950296548cbc2c8fb1c6eae962a9601a0151130f1c810040ac8b2cd",
      "a29f0190e9da251c5aead59ef623cbaeae9a865565c983bf33d4ababca58efc5e0fa66c8c83f55fb8059748145d4d463b199f70b3c2703d36e09e815",
      "d823e79b06fd486b22318cc8e78427dde95ed4768101f23610b46428691725dfe6597a4d825c9e7238dd502a447a5c960e42620945c7a254219713af",
      "a2e8fe2a0d391b9677611225b2e040700261f2c8522cd4068cca91c165f6f823fbe3e10190a09d00501d55be0d6c9a34c949018a66174d3bae5979ba",
      "27adba30aeb4da4180c0d0b4810439a2f05582d4d13670be02480490a655e84cf126b94aa359bf61e7a02bdad5ffccee0055e0f3c60dc21d4aa214d8",
      "33561b3c66d1006bab9ec26c0827f980185362bee1374709cf23831678d3b1a2464cefd0b714584c7cd858ceee69e4252b6169794970cec4fbc7a948",
      "01c20d50a195f2f040801d54374e349c55f70e8b6502213af38c6a29abbd87879e3eb1e26866d5a3499dc45186a7701cc5b052fd98c058d2977444d9",
      "98c56d00fd808d25fec62f3170387d8fb5aed7d61e1f13d9023200d3f81a6821225380a87c740eb45ccca6435d10856e80a6ac90593fa812042f2d05",
      "cd287c0534c7473f047db1861cf16b1314d2e9585a42a6308bc469760faf9a84393bb80190ac8eccf2324166152c33040d296aea629554a35af93aa4",
      "5b056684f6dbdd7f990ac466db01765818d3ae1672b18be665b40eb40d556913cac1512e29167f69491f7585ede890c935eabb77570d2de084256688",
      "225865fa009881af443199d9a3193e42dc7fe87265807fff47f72b35e08028d039bb664fc20672a9be71513da98afbd516a360172e01495709f2ff6f",
      "f7870f055504ec74209baf1c47a7712dd85f00a10d4d73d8634095406412950a587012cf791ec2376137a6bb1bbd036645a73f222189a2e047969480",
      "d7db6a8df7e51510e16dca78db5a0e2b8be41a5b9389935e609dab2b468d1da61ce51ef536871c0ae098292b12bf1971abb5d6e4a1629a4dc71bd8f0",
      "f05305a9b430775132b0c5c54128b6b43d7733f0c9863870e89d9281a94e878d6f1dddcd52d4bc9a6419ca5a42fa1af2a422c5bb512351f041bc208a",
      "1fe0cc4c8196c6f10c90298431a34cbd3ceb011986bc201c2dc8af1f34e0c8efe86d017ecce987aadfbe0356640eb5c5a7a0def450ae844011f14e60",
      "a48a7c7006c9e7b024299cc62c4aaa2f45b45272911394f89a162f1cc13d5e98f6ecae60f43812d1afa955fcc7c7716d8fc35de2420ea2a80ce2ab85",
      "cb92e263fc1106f3f0207e0128d5c4c77be9d7780218763f2ef7811e129b12b67ca084563d0604115ff5022fe80b32242ea355096393f43a5ae5f2f9",
      "42ca4f28239e40217b7f6e9f0295af125460cf254098cb17c066b08177817acac2717f35e004ffac034cad3f514440dc736007f332f006b571a2cc83",
      "198ad412270b886402f8a1c4f905a1f56a233a101505644e2d51d5fd106af1b0fe52276199ceb13cb50121fe85d8975af283f60514f26ff06601e4fa",
      "5485bbd5a8d48de3efb944c9c5ecb2c447ebe7bb5225757fcf56bd5bb33a5db33a1d25f83cf176fe7ca7aa24621c76f16cbdbb4ea4f05e475f6a37b3",
      "499900ebfe92fe4c59d1e30beadeade83e57749f0590aa0b066aeff68a599f9592767d96beb06ac754ed88aa84e11bc0cb9da1b807706e4f96bf8bec",
      "eb84520275a35473e64a2c4337c5b32310f7098ee0c9f27204f2f669188195a34780b7d3b303a02b0cfb7faab4ec5e5c770dbd9b0cdd39dda5cff64e",
      "a5a8fb27cbdf45f8fce876ac52748fa8a985d8f3340beb3230203ae9a8c23954b805cfa0f8982b2a932653fc95973e5fe1abe50a823dff6b164f8a06",
      "08344d2a71c51d30293032f935871dc85254b02b1acf0bd516efdaf1254a07ef82cdc87c8ad92e2f5b05e6b50208eaea4db1be78de5196c03279ab89",
      "6a4edc5d89aa8c7718107e573662962f43b525640931c6e22277a89a87072f9ddd5cf2dcba3de0a2816441b7043c4a68b2821b80fab41f17b4783749",
      "0a3c5f89af806c893f3c2c69b604c728ee572823785c2c35f0559a024bf80974ab4400942e81044aeb32f19d9fe832981e98f273bb9db9956eb533b7",
      "db995be9a886a26a10414ca0ad0b006123932aec36d0dedbf85a903701c78d6fcd72d51a74679d999cf9aa350cc841708566df666599ddecf3abb2e1",
      "3cdd580ddbb3a09aa7d9f418d52e9fa887a3b05748d55bd09d2a690d6c510faa7d2a7a927ce34f9469eb0b5875102c386d0a76f599f517006b8552bf",
      "092f02045087c045b1970db1c04000680379d3befb2532ebbabc0ca0b61999054351d3dc2a31c71273ab04e00ad2982b800673c72fce49df7aaa90fd",
      "cbaa403e8b2195f210654867fc02907887af2e308c44213b4c9c63e21c1335d2cb006b01e3ff95e79378fa1fe95cf76dbad63dd73a4e8a5f91d2ae6f",
      "f6927c5da2d3bfa4ba847d9ec60de869c1f034ba8431adc1708215a04aecd439a4ce21550dea02483a2c8b2c2f5c0f01bb41bcb10acde2710b9eba10",
      "969650c829d6cc1c11511cb0b739273e0a4af5dedb05ac13028590292f5e869b99922415d9ac1c7fe645c994c432cdf272bc1b1765bf8ab9bf06c087",
      "45a5785840a1939b5d041a150f6c6eec9678a76b3c1d2457ed13462fedbef17c988e1b10bbc285b0b5d9c04f8985b29030321a9488026164a9ecbcf4",
      "5a24e263a5a9478819eb65a25e66d51379508f1235c636038e04a3223ae6a2f1caeac922a28f44b423298f2700d25d13b745173ca1d9356abc06a2d0",
      "d59ae8d080a9e071ad5113c062495c323352025c4816a0db70e738f4f457bfba30d44665dc34487d11d965b19f6a591876675df47da24a2ec2eca62d",
      "2aff51b5b2b0bceec7945fd0bee444b1962aa96e69f8c64a7e608d7241af6e2b54d26a052ba9565e364d5c1ed3e94bcac762aa2f5b15dd3ccdf725e5",
      "b1f9c67bf48517e757a0fab4d8e036bc0ef4fdf1421cd52cd421fc6243901f108eb18104565e5fc730845f22422078d3c29fcd28135f7448282fbda6",
      "bcf49af270a1165ec1d7955bf099535d1b6a591b2add94344a2ec6958a510af482e569985c0c3315c3a4d1bee8dafe5f1cb0335e77b8ce681b065b66",
      "6fdf66778bb89533f730b2a74f19ab1ce573a96ecabca62bda125bca4bfa5a5fd2d6e9b2af627b39b8baafada36297b5d7899e869f256a4ca78ab631",
      "4dc3b52f9ecbcae87e42f2c7d38c38b9ca844aad5581a23bb44f923f716a72e34980e273e71da50d94735a08553dac75a1bf85d0944db1b3d34cc83e",
      "c3aaf6add32e9383a834394bdd36ab0a226ac0b3d40cf991c9df0ba54b8d334531d1dfa8e8f648d368aab6fefac715ad8e2e38ab6f5976d354948a4d",
      "3200a80ed40a3af4bdff719d10ef91e88c8fdea288a8813571562649af50f5450a4f7d5edbe04b6a452dbad80479e1cbe7c5ca46986de082509740ff",
      "5431209ccb80a0580e275cf2bf3a64e8f3c7eb56cb3690f872851648ad0822825532907c1194d5666469770f7a7cede7900e565dd861f4af05bc1a71",
      "07a12657e481346029c41e8be445723fe084ab0242cde956487a51ebb5e12285e29823872054cf10df35ef7a19f5de7457caceab6ef7a7d7ddde4f78",
      "5a05dc0cb342acd6d15e87cac889e81374ed9f9945d44f0e6bfc9ce974eb2962959f2363fc84acfdabb5567a85e47b0520a28099692da271a9447a3d",
      "f8e097a6ba4cc4a7c5b0b4af02ab53d9b4a8a80620135545a4c7e5ebc2a3359645c749dc5b82c8153482d5266e0d99be58d436dab6180a5b0db24bb6",
      "2d5d699992bce882b449f24492e489d81c45929732b994c982c32b85ce4a56b27bb96ae1d91a6c3e83ffced971f8e3c64f3dbea1e1b572eec4739a05",
      "33749f433f2b09dcd95cfe2a448902ae48c92f8813b8920071ceb89bab8fa1c82da362a558d5e099062bfa7716acf01554da581769318c70dd376254",
      "75c4e9a307bc897b131caf2408627919c11cdff4367ec29ff7c7615eb283addf2ef6b7d06c7367377cf3badd7dddebadbff9e9a79fdf345e1362f3ed",
      "55905b6fb7635b1c688632611217ad25812a800b91ac986bcb7a6b385ea1a559910a9ec72b1afdfbbdd524e8c02c9340703ccdf79218298cb02347eb",
      "1e0cdf5a4931487e874fbdf38ee82f08d6d42cd7031c45fb6e85abbc47262e08f90e7bb6e6c37ae6e58aac0173cc4b4014b0bcf0d77eb646f563a5ed",
      "432fd492b78b23f164481f9791785da48f6184af8cf263041f1be79aad8ed5231054e1ba7ca20ba7d62be8d8ead95a086b408f71695f3efaaa6abc7f",
      "244763511f5a22cbc533c3d3b5db77a8dbc61f1e7ac18a9a0d3e4c49bef612957fe7115779437ca9927923c8c34d7fa60b014b28eeefa89ed654e31d",
      "d923a30de76baa61993e5494e144a20c20deefd14c2bf4768f4ede871b6f5effe431432ee039722ffa10961b2ab75f77d42db392970038acfd9aad4d",
      "4af813c0799cc9e6276ef33f77bb3ff736bc474bf182c828353365642f8d4677271cbf3e9ef81edafe869dceeded6dfb76a39de5d79df56eb7db29be",
      "5e7b485ae9e64665c5088a912a25cb23cf6371d46563a3d31a6fc27fa83847961111a23d1f7894b3f85c6b16259b29194ee5ab919f0cbc7d2ff40ebc",
      "60d52fa20c2d1960573d200de8fd671530fd651b0634f0be41b13bf8cf5332aa1c78da836eabeb91b2c84519e989aa1f72aaeda29c4f389b9691b745",
      "a6cfbf09530ce486a545ce4d194d81425ad24d54d56d601ade4db11fc3928c3fe509564de3afc9755c6639d4f33d55f300360673d577c0aeca68ee2b",
      "3bad805dd3679c8e7220a4e0fbabfddd5a474325eb7b034acc61703130c07b69e9773ef3cb7f26e5bf3bfe5977ede7f3d5e0e1fb0ea0223ef4f570da",
      "b382e75bd7d837e2855e975ab886e941bbefb3ec7ac23d609ce6e5e6eb8d9f68ec5bb351927d9403172b12b0bb325a5a926652688111b35b4849a995",
      "e138cfd0a68a5dd2e0aff9f04b866d2e5de122dee13fd3921d41f95bd16d115fc5b06a01dba70ad3719c96d98d87ecab77788a6781bc4660ff17253b",
      "2d232430cd8cd0fc080f8bb169fc9ca450fb406c9b57565bf8067b2d566afbe4048e549edc99b941859b5e8fb6100ebd9c62a534cce06b097829f20e",
      "b26f473c2f5063124047b6bf8d531315f72f76f64ed0a5c2c5c60ed0d807c0397e83ff2e5049640918eca5fd92fdd16c809ee589727e419b7e935d26",
      "b839eca48cfe80795d956c57fcfa56b22fa64be23a794e00babc2c130f4eec647608c5fde60a4b5f002ab76048196cc26c3806fa25b70e0370beaade",
      "29668bf6deda33fe7878717af869fb03cc770be67a08eded8981de95ec58fcba2cd94720a837a5d15b7bc4bf26437e94dcf1c9314e59592702bb9473",
      "9ecafcdf76e0067533e0124f80a3c29c80ed94b6a98b20da977afd329fcb1b5c5a5989b556b6558835bd695c14b883ec1ef58c6de38968a90bf44e5f",
      "765b3dfabe57027179246a53cec96c8afe053c3643936d6dee5853c47ba6e6e3302e8764b5ae9f201fe142fc130f5e33d2869317a75fe3c20b889c86",
      "fde3b0e01f68ab17a179d84eb80e3c44ea81cc3bf9f5fd312c52c03e63571f9035f43f5bf8b3d2e728f9ea05019c3d98d287d383fdc8db84d63abf78",
      "eca91b040e30348ada8345b93d4e26a3e565e7b39dc637bc98c643b466b6b506e78e825fd480d7da65b69fdd2ad3688d14a459e565749ff010f04fc2",
      "27e5cfe10d5c5da36b1ea64ca0e0104e95c4aee1b5feb9be117ed51f276536fc12ce4b46782f842516082fbc2d19e1ba10805ba0b4f0a864129b85fb",
      "b2426f3dcc102ac2531cc3c6283c2865dff0fb9b6c027e5e608773f8b15d3271f4c33fd4afcf62ac27ee3794dd8594429ee7f04b298442f0f310680a",
      "3cade15bf9e3638cf82adc522d1cd25cf6d4e77b9ac731909cbc4cd238fc084d59805a843b2513c016fe0993fd7a1d7e28d9d79b4948f0629f1e79f4",
      "9e811f86ca05160cfdb0f9352cc631aaa58ffe8cbc9ed7f9e50766c387b10317f777fb928f011432203a66f9c4ff2f38e1f16c52fed7af07fb81c7b8",
      "d1c833b6e27094ff34874c6b53e3214bd209e0871398d667589178183e77e11cc4439804549add3d5b761f4ba1c0f35d195db6f5760d3c8da377b2db",
      "146819b97923fc629f16153e00bc620a2396f1d8af8b0a7f9a9aa2b3a9c7fe5a5450184b9bc2c218d763ef81d533d742f84e8213761b7e921ff8daf9",
      "abfc2dea857f0191fd4fb76a45428937d421f9eab938fdfd6857dd20bc2d07708a7c043a2e682e7788e252cea52cd48c89733326f55b8e89f347f63b",
      "3a7061bfd12d61ecde4b57c9d3b3ee41f49df0f0f01bba0c59480dfab02adf932c842d2e038bf5afe7cae022f267cafca5cbe02ce0cf3fd13bcc0088",
      "6bfa41fe5894dd306ba8ff1e4ba1da3eda120233d86458443aeb2d247292ab848f420f2fb699c5027c4f1cc5efe559a9b66b6f741e95a6c0bf9a0ac0",
      "c552ab62355aa294f45eeae5570a9a5272cff125dc8115b277aa02cbc1e1a7935d600c6e32b82de0d672dd82a0a707041532a1fdbd0cf4b7301ffbbd",
      "244b6ade1e0246bfe6a35355180da7d0ae5f1886f3082e39034c29171e321a96de1b5d4e869364f8c543b344e504c4aacb8dc099ac64806e29e36432",
      "48f5cf504e85b486add93f3c94d24c6a3b9ec680d181fce50510774dc9a4c35f186aae00e0062ef723275daf02ae36bf49466f3ff621d2de8f2ff944",
      "a275d4a98b4be1830d3bc8f2472478f9cb5bdb4ba733c5a2559e2b31fb04c062a8f21f0331627450d14eb35b3f08d6b2cd28e10354a75e5d458f377e",
      "b3fec5fda336d2279d37dc639692f316443ee4d488cb431484dc6824a501d96c44d626a67253a275f89914a7f90c367a840e81d2f64591dcccd04d02",
      "1a49b3f491f43084fb276e29993500896c1c4edcbdea2a4c1834b7a3beb880bc8cb382c385cd62cec69c014f3be3d12d3abcd2e204a0c904a9716aa5",
      "1cda1fc0cfd89f3785f93a0fd85036683522f82a663359cc61da4c9bea5b362a3ecfd1503eaab7081b370430769aa2b4c190af7abbe908ee26d327a0",
      "77cb13ca05b71e18a4919b651bad9129be4b0a087a3bdf13362706b14cb92daa95620238eb78b486b31c0d784e300941ddfe366eb9bc78566644e926",
      "e4a2813ed1ba00886e238d1084caaf09bf1d24a16fcb294c8ea0e76fe044f0117582beb5669349100c107d93494020de9f123383236d60b388f8322f",
      "1428710302ef2350dd78c1798092f0a13d9e4e616989e2a2a757d3f8a92dbd9f923d13ca08fa746c04b7a3aa994a375c194af56d62ae1f946e159bd0",
      "b36c63f8537da2ff13aa81981213dcb15b86284f36637a56ed00e6012ef12d87611348340fee2b770486da94442e2c9ee781f5bbadf435d0a0a9bb09",
      "fb7ec72dfbbae565143e1cf3ebddbba9eff9fff7e1dfff2e02b25ff1e1d7c3f701706bc83d3a133bb000d6785a72c7600c00a377f8acb2c002d11e29",
      "5a68a005a2d40496137d7898e30f1f605b8c7c90a0602ef4bc60d5b6e1f92607f5dc7a88fdc7d5a0763ffb3e8af9a8696c37d09ea53c61c7d3024402",
      "c981d5d59cdb7d45ba2f04ea36fa54fa359e0c1c400f1b4a44962bb43b0b91b474b33ebe230d3338f2c5143d9aa6d7f23c0d9a93c332787a204dc330",
      "83d8d6048e9701439c94c29e4c20a440fe6dcbac886b7d6def2a990039601796c80c6fa4c8dbf98ddcf26804df3e488679566457657b6b025cab4772",
      "13447da2a1a29d94fc062d9035d3240de195d9eb23b72d2b7addee0a405802b4687b3745075ea308d290826649fb508d171f0cc5144437ab91875e5d",
      "ae8101f756530078353369bef568196371c7eef23911ef62ef5c48fc694547e3b4cdf67f77c915362504291e04baac8b8e9dc4d530e3e7917fd9469e",
      "7e20ae26bcef7d1836caaca7770c7fcdf157a0ee2e2cb0316a2a020da3d03b19782da130e0ad2638fbd0f181f787d429319e0aa719ba31b96c93e860",
      "20068d4f3c728927fc0a8d36a9338fa9d4329b46b25fbbf523eeb8e8b3fbb017c054d8478f4b27ca8f9537cae36bc135b143abd42996faf25ca90373",
      "55ac916b8976195fee210f0db7858bb7e534809df8864f46131ef96340ed813a19b392b8798281b22131f252a0203c6646fe85cf05d7fd8ddb08069b",
      "bf1f634763ded4ce04a9ae08fa16a8827d79aec16d33c925df77a715a0356ec1cbcfc908ef83527e7ee0a8670c4408ac8986f7cb6c340ffae632b720",
      "84575ed9e90d1e3d0d4c50ba8c823e5fd5bbbf0b81efc3ee3aa5ddf9c3438fcd216b4c5d77aa23e989f74ba74da449638e32e48208752d4456431ef8",
      "39b765b566eded0a0811b663d02f0b4b018c1792384541a0600fe08a1092b14f76d2a1f305c4a9f37d53589f409cbac3245156bfe0a42ef30cbac990",
      "9b289039c47f2590b9f3a1869e6e86ea66d49af41ef018f4a775df70a7c75b1f4fde1d1e1f84332e3ef6d0e37138b4bf2e763fee84238e3a4ee105fd",
      "21c00da75cb9943b9216f4e1a9f0dc16dea0d0e45d9ea1221bfe7c1b03c7710d1c465c6ce39d157ee568e8297e1fc8dae2eb1b7a8e133fe7d499f87d",
      "c7854739c2eae1362329789e4dc35bcad0f7517849df479920f1c33f48354b7d1d518be8e91b6f96537e578a3d737465724efa5c2d9e3e5128e6a456",
      "8425e842dc018414ee73262a99a453cea63989420ec5a10f0f90dd2a4aa00b55ca37aeec604647fa2487db2295d4d6fee0d60bfe89794647f9459394",
      "52136c198928822f9cfc637074b790545c35082f25c008236917549d95504d8eae4ab043c7470356de8d8c5f5ce1d9c2e2a8be38036d747e73424a6d",
      "ca8dc2d9eeb9311c7bc9bc0e1be645260a44a6ae2fea91d56509e4e0f43d7ab8c3a3d2b040bcba32878b57c66c156f70248bd3543e4112ed8076d409",
      "cebae7c84e2c7120971e1ea8fd943836f14f228545bb3cba27018ed0ef17c21cd40ef618fd84fd80a3289367809e6ec79c4fc225a00453fa693dc9e2",
      "bb259613e916907d41b73ca5d83c5a707a9b2f57d145d7aa9f0ebc0b0f7ea70192167d0e3302a607fe3dcb81d5c557f726076789f064933e3c40d3fa",
      "ad8ec00630dcd265db7a64585ebe54a229d2ce22e78252164ec5bc609045ff2869945910caca301f2d5ac14a5026c5994009af2a25c165e0035b4e8b",
      "35808dd58260952056477e58eb453dd4fd68f8bb9cfcbb20682ec1b49c4790e5e57bf91d02a518c87da2ad747aa05d7cba07e6ae72195516957da64d",
      "844d2e480576a917840d6d95b0b532ab2c6329ce4360f156716919ee6b84ff903315b9cb51665124875c9d74287f9f4519ec7033a4d0796239aa6c11",
      "cce0bf7078fbe83be11900e0750040375b0095097b0fb765c4cf0769e393e97b727896d3245f26a946cfe0cd10c5b1d384f939928a8defb38e684ffe",
      "068aa0b9b02e6989e6a06fafa1b0f2c4d0d4ceae90317121920f85a0d7d9494e927cb1da11497eccf66db9243bdc4dd323e94206b207b514f2b59de5",
      "c97592c6c21dcca0f2ddbec02a53925b229f4a6f276f67977035a224d345926fedeee9be91078e6df1a05276cf2e8b17a1382c4826b7cc416ea975c7",
      "e021f89e0d37c04c36673812797d5f50c9fab00e9d459197f88e10af0daa09d2fd3816fe157d5123b7ec36776c8f1cdb66b4ec9552243e27bb6b8c99",
      "003c139204706b19819b95ea1b9f501c5f1f9248b9fae927fd806b8f8a095c3316e7a3de404daf3bdc71daaf4b20acfb49f40727caa04eb033ade136",
      "a4b4dfd6526216834ed2be5be33215ed20992ef3fb1aba219f6291b92e720a09dab98d6a8ce92ae211fe4f1e5db6e9b5147183783e1f2cd00481d37b",
      "130f071b2b8bf2bb9b0b7206eb0beb58af7d1fecddbc6ca33600c00481ed0e9f94f1ef9d753a85f853dd5ff485ba5c833595d5f9935bc57a9562ebdd",
      "155dd22ab65e29f66343b1df50e04cbffe18e00934231bf8ce38b19cf90ce4b0f141c2f147209282cd8df59f7e7cad07468f49fa8149fde8ac51a995",
      "1fbb61d7f2dcef48c651da4b0f1ba7125ed123b671b28632246270519e8cf2bf7ee00271c313bc2232b038c1cb4d9ddf0116e084bc777e61159c166e",
      "7145cf9f0cf36c32b1f3deeabc6dc41576d69ea6ef2526080fa9e9f0984ba31055144f6cf891520f106769a66487d23eeb4d083f00b22e76efe0e603",
      "7c2a7d187e268649217c9885b80674c29747f68e4789769a759fcf9a9c760a6312189faf7dee4db40a6f828ea3af81392914366c5f8c66398d1edf1f",
      "daebaf55a5b840eee5083566a25e47ab5a0369d77ecddaebdade0ab1314c343ae2c6084b080200045ccb7e5dd8a98c0ec52c47f29443decf3c41b2a9",
      "c231b966e7e8608496bfe2c8d69e9c72d3041df0a96f5c87226a453ec4176e854483f58630636f14ddf94ec7ccf6704c0d93af20f8d1a012aee7b356",
      "9929a0ea1edf587197becf3713e9cd153695dc9afbce36f8bc8397bb56cad6a57acd5353f9cd46667a1f48866f6f59db711314f4e9154b79d6617fa8",
      "414db4e991da283ef5c8ef8c1c85639fe2dbab1a344062cf6e0b9fecb02939f326739ede9a092ea2ec661c80255b972dfbac4847cbe1fd302fc249c9",
      "8682b3920239b20693bf6f92f40ffb33beb33fc90f5f11c2550c19d26e456601e21af1dc6d54871308972a09a76398fb389b8cc257ec2a1e71a7e44d",
      "9c7fe1f91fd506f47bedc6ce7e729394e19b8d376f7eecbea1964fd2781af6e8a7c0313d2c3ffc72cc0b74bfb1d4451330e5d6be66073854ee9fe56a",
      "8ee3743441777dda69ae9879a49d2d6247b400fb950cec6e5ba8421824835d6f8ba72de4476c47e1501fe32868dc918a1147b1ef26b0aa8771c16da1",
      "31405b6f8634ed04503e5049be951d185f6ebc8d33d098023fb4bf51585b5c7c5f9409b0070130cbcbb5faaa3f7c8ff56f7d5512eb6001768f81b64a",
      "5cffe6b001d6728acd86651b126f44b2fae565609a2c5d3ce54156c073db81a986a6f4fc841c0f93e1dcddcda1df135f7b92bd9bb391b6edc06bd73c",
      "aeefa623b1f01a0d8fe48efbce400478206e164bd1e477d457b24a1125480f36acadbb36451589d2048cec10e53855bfd6d62678831930b4afc04906",
      "470cd73269d376c02ff2e7873e8e8d77bfb6c4f4e4dd980062e2dfab5b42e7c2391260c1905c491bcb3085d92149fd7c14e5037bc797cc420041a436",
      "ed0f0358d564b19aa2ffb056e228c6e596dd288ec28d56a28f28a6e4622c8e19d7052d1041b3e80b1ac3277ff8d34e33d42a3cc1db43b9cb80b27fb8",
      "b6b435fb1db1f803e7a028836f6ddc57b27b0205fea82f39712499b4bfc4afbd74413f2817f1d5eb56ed80487c28175e8ed8ea64556807e24ff7c6f9",
      "cf75b126ba90195b150f754e902eb53458f0445966673a1985c87ed0d6fe03c9a6cfafb9265586a06470ae31aff02146ea5b8a32cb02fbd2efadf53a",
      "69a09c3b0c9d7aca88dccf04dd50b1ff535b9ba2c36edacbe4112f72f5b85488bd6e80167a39c5836b3d451503eb37f0deda09027a1f03c6174690e0",
      "d4c807cfc3834e797838ebb2ee394ec12a67b9e2692a6b5657f4274f1a10268950733091b2b2a8e63983ee17ac32b01cdbc824144b6601e0bb5ea73b",
      "b89754876bd64a4b953da285786ab62571b73875acc7eb4e0d326d1eeb14717c3e044c0f000b1ae3f99426e9f648499c6ac90192fff6b266ab2b7510",
      "a41786a02dbd25f8c6df90c4e40d7040b5998b18e4fdc94a717d7274cb5af1e72a9aa108b96956b6c8b96bdb13e3fb9ce593c5de42f504fcb3b3b59f",
      "bb6c8dac72cf7e16e6b9e7e41a00109d63455ec567063b9ae3cdeea116c517c3ea8e8f4ba37aa323940968972e07a55604ba1e444e7d3e007e8564a1",
      "f22a91cef2ccc9d71a42c6b213fac42b00c3d1d9dfe4e5e89dbad865621b06754f5c8b22afe0d630f7fd09e430d48a6ecc056ae0519103dcbe0ee428",
      "05098f6cabe2dd964458393999817f208778134fa165ce3c29ab5b83e1ad61294f433242cc812806bc8a6fe12cb36a95b901df53693f8171aa1b58f0",
      "b7bc8dacc23ea0b5384f4ae3f2ef228f6f8f70f37c97b9a16782c0fa200e2550d7af755dbaa7ca3d8d95fb2e502894d9e8de97dfb2ed4a97e886746e",
      "8326da530a592245d9f0e98eba7f0c0ce82cc90b2b6846d458bd6f934b7d19e6e9c9a1e32b915d20d50e4a3490b238b20ec738f2d3e816cb658600cc",
      "06719819af12e8f391a30346d8fdf14afd0e8c11134da2bc6d3c3df845805a06b3a8d77eb5ce86d16c6566c47723e9e7285949d6c62be355bf1cacf5",
      "e0825e19c2ff4f562641c75f5f81c4241c63da24506eba84a9fb4ab98aa6bcfa82db44af15d0c29bd07808719d696b4fa2c6ce5b5a77e36fb496efac",
      "3fda4ece9b2aac562b90942b1af95dcbd4e9cac26ae3157feadf042bfe05a293d5d94a1974b0a160ed029283ce905ab8b654a0d9d7c81ff930b99ba0",
      "3363739433cae33120c184fe0cf1f3eb4afb8d4b1edea823ef2373e792844647dfd7be5a4dcf6bd741670edbeb3bfc3b67bdf6eb2058f9dae79b514f",
      "613482741261a018c661fd481fb072d68424233718a270a89a2bb8dc3a1380add8f244846045876e8cb6db09934b09cb87eb88fe1cb0f8bd3873c8be",
      "a9134f03c07730f11370228a981eab61b0e4716dbc309fbb139d25a7661a2e462227359bdbec96d6784652f12e343e31a8730fb0fb88eb96022b4ac9",
      "0b4a63983875e70ffc66e64c0b201523b6a83585d1611a2fe85aa1ee7a7f240b969c03ae9214ea2c707aa2eb8b62ceea57f26a33e1421fd8b73036ee",
      "cf847fe59342988478552cbfd9d47230707907016e6a935f347c51ac79f822ef3f32fc5f9a5a5e347cbd6dcd0741c8eed2ab2c1f26ca130cc68b93ea",
      "f53536b1412c90d884d8b5cd9724dab15720091422ee2cb2b3de758f55c6dde8a2e5c6f7251d17fceff0228b6fdfb2522060855e277acad22729892e",
      "fc33e1b85410f8ace26f370dceb176622e6eb3721699a9e9cafa3e917f1ac3b4b823c5071d29802d038b82b5b034cbd0833fba97ed0ed0a56b885e66",
      "33f4e48ffe65296d1e6a1fcd72f75cdc8ffa0b8b7753e1e244c480403b6657fe69873e1487c2f5596dcb77003890ce47592a11e9e5c064751f43d74d",
      "9f2186e0ca7b42344a060ab466066b99cac2c391cd9829ca17c36c34a5a7361b9919b737c85fa4c85f9446ce858aa5e940ad2db493925230897c2a44",
      "71da4014337c10c447e3213f10c87ae03f2f6eb212222b16706c7010156386f265ebdd6ed04c143ba3ca4950ecb1fb6c32a26de6f87e4dbf12757fd7",
      "5f8f2ab4b125b93112498361548fb6f40765ef0e07f415e86b620cbcc011473eb2493674de9c6a5ead35bf2c4aaa68920085a55828a0cb5eb15b94d2",
      "a27a11bec478d73ca3d2d29585b63246579b4a8c2e24faef4dc963b24728945c3d595c92d86e4583b5a9eb81354628f119d3f646967db3352451433d",
      "86fac893c1edd05c14f187309aaa9457c457f3e07ce1daa6cb6e7851c4d73cf4ac12242428840f053e6a7b8f5a76097bb25fd90fe3daa7697ccbcbcd",
      "c326a0a74578be84dfb87446526dedba96de3ba98a81b3b0db826569f050a85e0ba560cfe8ca24238ac888fa27f816803ec0c55a02f627638e8137e5",
      "f94d5214c45ef134e1238fe2d2608edcacd60c661f27137c4847630601b15ee3e4d434d0fab44e1b2a518eef1c70353b8eb3f3a44323debceb5486c2",
      "1dae7ab8e941d332a94350bd109e5826db4c5479821b6618131be5ad4939c3e5d329597a4d4988af512d465eceeb2bba483c1c02c7359ceb7bdc5d26",
      "d4dcb5d6cae74db2caac22454b59d120943489013a6064e89014a3fea4cad358c670c74478f0b2ad7f3f069631ac1c765077202c73ce12b4e0ce29cc",
      "b449927206771bafb05b8fe17b2dea3e7c108f8ecdb237aee74e612661e1b9137358bd58aa989e0e96862a68b22754b77dedeb5beab4bb61906b0f8c",
      "5a7fa0c666fc5d46c905aabd910e71d80c6835c1e7413c6de9a2ada4685d72a0755a394705b9d6e5bc1503aa1b63967c0cc0a056f9fcde0962d7d817",
      "6b2eb237326a384fd456862db5da32e351477f461032367f4d8850bda86adcecbb57283b75257c7a5911c1e243479e4de0ee50f44b25593726488563",
      "115e1ee1d56f486f2c2c742bed5e25dc16b5e7460be46729261a4c289e4bed0f8c472bad1dc5792b4dc04080265e04a7bef5455133dd7736f110af5f",
      "eca9187e3a7bab04a34ea2d25d50c186e8b1188bd55e6a8e846716476ecb3d14af792d950428d79865622b50dc0b563d5116b58499126beb5dd453d1",
      "779e9a65045bcc9f76772ff762cc875ff6aef669f1cd5e68f25a5f341409de0f06d502cabba07a3c9b571fbeb4d04694dfd72594a378974baf8c0f29",
      "45193ac5e584954baa462ece0d12509531b79d10310ab3554ba85030148682d4d3aaf290ba356b93886260039a6ce3e1a11b369515fd540517cff523",
      "ef2c37577528337b9d6ed8901e36b564ad77ede15198c093840c19fa44b1e0fdfa3bac0f3c70d77a269473f76dde5b0ec327f978db8e60605cccea70",
      "04c69d6d035b8caa5e1878c461f39962c0544a8e416b2d2e338f16be4e4bdda018c3c7507835e4cc29ba1a74c4075a081fa3ab174d2fe01733eb6004",
      "a6a5650d9353b05fcbbc37edf879a787fcda8afc811cecc084164b3b79b0928756a4339162b99856ef0230080c851c3c1b4588384a7dbe2d365b8b2e",
      "48a5c872fd6a5d644205599a34d62e3199adac1abb41939a53cf4e556844f80c376779b10c588a8ee4eb3c737dd56208ced268ea19c980d5c121a9e6",
      "ff4dcc3835354d4b44743f2151b6416b980b74e5ba4537bebe65a4197d864219d3aa72b13ce192c3bd09c252572f16ce950a8beb426de1e275512534",
      "f23ea9f890b73c7d54666e9e739b759a507546f974ef24dac7fe23b38fd28b7a42fd8ea77af9466772c5f461094065745481339b9cf2db14f6c23e6a",
      "bbee4432907a590b9c24fff73a70a21650604eeaa07e27576265d08bba7d56ace3e1571456dce8016a4ae6527f2e0201cddb049eb402893677feb8c8",
      "fdfaa2b365793ab71cd62ff499fe7c33f2805e8bc616f8205fd48cf6a04fdb7deb930bffba86d053ab47bb53d109af3cf25736d9d55c5ad8a6d163ae",
      "37d7acc454819b05ea4e9abebba9a93b35918716183d3ff05a2b552d2d5b83cb1a0a8102f4446651c4a2bea4373200aa60bd4a232f80fc2756488440",
      "5b382601726e774f1e8445e477f38865e45e47cdb8698bcde423f206456629e5739c360a104956d1f6a846f937d9f478020cd568de320ad8d41259b7",
      "15648be269ed67619be2bc78db0cf55884ef358ad30bd80853897c688a4775ad670e348b73be1d0d74ad7bccd0739066f47483c0ed49cb498be71306",
      "9be4e60772857b502b5b24e87c72b26a65679351c275ae70926a654b47d022bf3e17ab24260b35222c1c30f2d6a5c59390d2f7e2cb229bcc4a4e3e1b",
      "97973d32534217bff2fb0a90f6487d146532fc32975fda5d8a6a2f32956dc5742250fca0a2faeeca25dc545fed698db6b1b75430f4f78f7df3792c19",
      "774bff5e72d412a5694edef76486574303dab043c9068c1317a5b36e35532613fae535e4890856f345d9c5381e65b78b7285d5c3a2dc32cb509f6351",
      "f6349bce9c4cc308d66c29802900c86e9bfe0c9093aeee3819613358c68cb8b10cd9daf0452aef4f2cabe4311dc937f097befd8d4f843c6ad28d771f",
      "9aa6b0ffce5b93d22c57ba5d4665278361290d1ab8eb8d024b566b5535993a622c21c4c2f8e1aacd05d483530b07a41402f9425541f994e18a64b55d",
      "9b211d8d105c3cc889f85a913bf77e4df4a2a357580fae659d15325cd9477e6b1370189a6940deefa6498a0f3d95097ae427d5cf1e1e4ca1a0b154ed",
      "2115d29ce9efa6cd2aacb5d6c4eb68b3b2e0c58277cebca66255d1a2b4b47348a9d2b2a44070970fc2ce00ab90feacfea6540c27960b9fd4177068b6",
      "5865ad267c8136ea62f6ba395f4529c06882d87a290b1489c8a8cde69e2b8f05aee2408d1238e165eb06a90102b7569c8e5ab86d2d7215482acbd643",
      "834baf08db0e320d259c2f350b065fc2933e57a13d65fed9b82ad408ce098e039b0e10ae088ce5bcb1aba75fb3694bfb196929d722e20742956d68df",
      "fac2e753b4dd6b490756f817fd46db4f313b870722088045caa8f5b52cd39697b9f6dea4deef1b0dbfb44154b525c7260e6fe9817a959122c13430ba",
      "79acfa5093a507e2c4e166a82e9bec175ff03270e75781cf553bf1ef2bea12a457f8a814d534f5578761238e12f4e269368deaa22a9187a2a4a84bcd",
      "1d5491894b50540e68bf21c61bc65677a2d607bf44cf6c40a35eb4a3c4e5aa3cc348af0000044d2f00ba72a9582164f03d238b8ceb9b881c6f5b1e72",
      "3001a5b2d24f02793cce87d27716cbd1a377d117260b4acfdc9ca22238479330cf380381b601d2cdb70e1e89fee9aef139f1805e130083dfe7a8b876",
      "0944c99747689e6e01a57549aeb8a9570ceef199fb05ea345059642852edc9210b740345547b1e94358aa8b02de595b050fbeb7c78c81f1e32158cd2",
      "19c3200d29b6dc39a29e828cdd7784e1bb7be50867ac30d9fa1b655fbaaf6876b36199f83b7eeb1ef5f3bb420c4d4a300bf6adefe059b861eb9d73e9",
      "a2daec5d499e9835c1ed4e16d61ac840e35b84fcf172590736f980fbeaa2bec0ab550d9a005210408a0f2dc233e5e5c5720d633b6db23c3579b6bf92",
      "73e6b45d2320eb73f1fd58847342c9937034adc1b369b4318ba9185218f6bb46edc409b13382ad396e191eb622eaf60b0c002b7d63a177ccb3e2dc85",
      "aae5e54c42306405fd34ca5054308c51c5f091e05bfae6ba77e62fa8b94365da95b3384a3144d938ba779ccd84a570c7eda97bc75b324ba29d2856d2",
      "e03eb252fc3c8a11e7085102fa228a51b96d94cc0066f4cfcd08031f8d2b728e281f3c654e67b5eb2b3da5270521d081916e3c2b7aaa8e465447c50d",
      "1897d3f7c2774bbbbf20905b9a3a5b0a7b84bb4a542467638a1c305ee0f207100b6a045105750d5ca2eb1fa0b10eccb1585e962ed71417a273f00849",
      "d7b070925c4cda40f892074c2a95e0d699df52a183de724341ef3416930fbe92f2bbccee2c0b5bf9a5ca20bd67bfe9db57a7d662880c63a6d54ee462",
      "969bbc5fa2a7642717dff625a2c22e6ec71c088d78347f89d9aca05d24dbc4ee057a14aa8c41a895e58931635ab15990d7e6726f20b08f784591c271",
      "de2a189117530e864359d2d194c56baa1889fdc2d53869bfb46deb070b9930745bdb2c5e7f5eb05c6f6ff1fb8b7996bdb0946c9b421f4b4edab6dcbc",
      "582c8fd6cf05301c854f604c8db2576938bf90234d160463ae8f487726b81ad365c3cbe4cbfbbc90e16b175aa55646c11618a8beb8b4fd7cfec2b69b",
      "4b2bb0acea712ce2439f937b37438bd3c92139587942dc5d13f437be95346a9ed8e610cd56e52c53beab13207baad6c4ec390bf32c122fd18e663913",
      "660559a02b1b7b26315b943761db9a81c8809d803bd6723c95014b01290a81d7ac923319eb9366585bc3ba712d46b81163f59366859a401b44a0e207",
      "1944685b4cf721c9990a4e046f77dba2bd613cc2137845798357c1b10aafbc7a7212329fa030df1610e06740418e1d330e1465de35c5b705268c343b",
      "7db4aa58e3c8cb55d3e7983e9782cccb9a5b022ba65ab92ab53804022bd738fa12d32c639719158f5205b6b53284aa0737e0fa944593a3629334a8d8",
      "a4cf29bbe8d52097ccf6b83b69b0925aca269cd9e18303c187d76ca26b2cb96d0751af5261babfbdc0ea798114cf75efb15887c4c632883aa02214d5",
      "57f3924f1e4f7c0a6621ed3f02a531d3685c1e284d3c61194248628928838a2f9985cf1598d9a03088c92d47982ec733f202eb6903a570b45ab56823",
      "8ae6c1bab478ac0e4633aecdc8a9431da7a57fc91dc737957bc331b2744515c210599a5ff6e880468d5d2896534c0b0850dba14b4660ed408aaf4440",
      "cc21e75a1549153627053acef3209697eaa5caa51b20913c9beb1dd202a66a86bd75a7f6b2a8a1a066f3cbc6e2e893521bd2c5dad332285b2ec42b3a",
      "774f6f95b6fa727da5e0ce208cd6dc183588572b1bd4dd8c283a3545ffa4a01bc61fad0e77143cb39117a8824de6805b75ff72f22c56e57626e050f1",
      "76bead427ef8de824322790d173934aa371a09b53355ca50ee9296ea5ea3248397547c50d4a6862c95bec995ab3b3a38bf2cf652a53d9f294951afff",
      "84e39bf4090ca72914e1b586def22b8e9c16e2b73430ddd7c5b6ae957be27aba58ea05da29a0f6c94498511c332500332e061bf646bb6a5477023a9d",
      "d2bacfce6e195341d9dc695679bed2e97f8877ae854e363420550ce70c26d52e629854cd4ab34f53945f87a972607651f29be93ba887bd1143afb52c",
      "6b39eeab5cc368eb56fdd55936cd51ba71232945dd84b0e9684a2cb5feba2b6ff85a89663f90167af02bfbf5ed05abcc1af7b3f73f9ef482e5ae3f37",
      "561073ad469351a71e9b746d80ee148dfb894fb6635c2436dfaa18ba6f23dee07151294b845e994d73d4a2f516b8220ceea51b42a527aa6a2ee0ca14",
      "7a513d3c3a7121169096b85d7d6353a48392095d0cbfe1b14a6b7b944c504fa351bdf4dfd73885661c3f4312ea940988d1e6884adbacd8d26092d4c2",
      "9664650cbba34c2785de39d580d16e67798ab2281dc7a875e07357bf074a010090d82e3197df25191878c0a3b941cfb81b7a372c1d4a8d5b5368d768",
      "14314d662f60dd1acb753c124fe10c9ed65e3a15d1744cd5d5f3a0fc3447d4190d92354f0f47acbc31f521fee82a1bce8ac3f4209e36101334be12a9",
      "88520621ffcdfef85dd10d38001b56e056a2767df25afa2b8ffc2d58e6e16436e2febd81b546c9419b80c8afafe293951c10b3b4e6842e5295504b2c",
      "b8335084fa4669a481c7d322005d4c41a96200d2554fe659aa64355da4bee54346f8163b03c6f3dc6aa5c4206fabe92a867c782c5131684adcc515e9",
      "90a86f8168284182af5d4627c9624aeadbbc0a8ef55ad3620436545a270d4db46a591a6e6b2677a6e6823c853b004ede36b8b6cd2693785a707c5566",
      "0d589761a0c93f28561696400da613b89226fc6d5c906f8322cb4be1db537dbd938bd040bc280df3cd141d2ba59bc9a017761f173a9965d6bb528aeb",
      "3834d1878d2d8f5c7f0af25af73d6b7d17a598887e042749bb705d8faefed50d2b3d779bd4c0f562ea9767e9399ab3e0d8521172a65282530941e011",
      "b26d38fb0d7e6c674435c1a76f21f240f3ad162f45fa25ea5d52acbf404496efff6e9f6fda6ba19e1bb87e6a904684fc5cbcf5504730879638ec46db",
      "813a10862bb297975f50b2e4db36b95cc7873581802cf50aad19c1ef80891fed5d7dcc00b748c0f469098f6b78be8299fff6fa94b03e65e3fa94d5f5",
      "29f5fa604f2f5a205a093c2394f3d4938d812a66df7e031722a41b0a287c28d4255fdaa445f92d6c553a8daf8c146301be7cc27d4bac243a1f3beaef",
      "ea26934b8df18d86eaa695698a760070e83d3b66012a36401c546ff91a99b2263a5a1375b9b1ac2d441c2fa9a92b62d3d13ddeb76553f212965ce17c",
      "ad661426fca19f66d355e0188c1b38b70fdb34ccb83694b90b872cb4742ee3fce9319722ee62f8edbfdd6cc3c1d197beba2c1ac8ad6fffddc57f462f",
      "7d412b1eab53b6eab63764adb64451b804a80f9fa3c3832d154f1ca03b4fe2b5715c9032b047c7658fa4086fb92147e53a621ff227f4b3eaada14283",
      "476f191623ae68575bc3420d41528b274007dbd1ba6cf477125ff1c9dc0edae5b4f0e8bad616ebb19fa45f7070b11856995dc3158d81be837ed61e03",
      "111a79ff05e456bb4c4a8c51b92f1731abac052c30da5dc347898aef30cc8cdd4b55090733f43648bf04b2b6b3915697b1c78f685b8414b7eb1dfa15",
      "5caf0b2b696622e546225bdf8c970a9552845fb3d73059ccf2b4a52b6cd534ced1958a5348a75a25a5e279bd459541655d51722adfe7eac8f265f7ae",
      "f615a1afdee5e5b1df781ba307e03268b0d6871228269248be6e9e6c71162fbfd11d70127a39f7f411962cc5a810403acb9b2779ac30bb8606ac6279",
      "c8f063ff795b319b6ef44b317fc6e5dfb28d0340dfb1a877ecc2bf6ac090a9386bf4c8243ffc2af5b7bacaac02b5ece0395a445e448d0e400c4b72af",
      "a371b8901b309d6183df9304ad2159041d1da15579b79f35021946dc4e1c13788c82a8e9823d0cbb9cd07b206ac4c93130f202bb64bed3d5487f0cba",
      "61afef2cb6cb05e83576672a2fa65102577d3c2747b3cbcbbdcd3418785e280398568f6ca50e221653fad18ea5c85c186e20af1de25e09a31511e108",
      "476c8a456a1292cb3cae57008f90d1a81b2893154c0dd5873c6161b52cae0a6d84748c482187faa5cdd993244d387b97ef64c7f128c954c8efeaf3dd",
      "0f9b0982464ba8f2a1125ae6b528e877b4f09ea750a580fd5a788622ef07e0887f807f7d3ef8a145573e1f459efce1fd4046535ee717bc6675485931",
      "323928812d359183af1e00f61f4e0ff62337acbdc05408768daa9b8b5a9fc4971cc54aa9457bc1252dc95a8920f13a2ff526a177c545cdd182795afd",
      "91267a99dd7918a44387707f76f1d02db408e7b42d972c0d424d88d436cebcc5200cc8d628eca1c0b8e9e253af1d12491c288cf8d4b4917a70a80b38",
      "0c548f605dbdad2d588a02d09a2459526bd7503e2211ed33358dbb47fb62ccd02d6345949755ae4e66ed561d07868df8d26922799a4a4d082d989568",
      "b6bbb0d97b496632f3f05bdf0ba53f0e88b8516cd0ed9b90b04ae575add7ef6e4645bf585b0bca08155c590ddb946a6b037dcfc9033848d5fe873aed",
      "e1215389d49fd0c7cc6c7dcc869392a16aad8d6a2c564f645a8da5cf3486ba9b81e583a7adf94ccae93f2954b105a144ee35ed61a300edb9cd719833",
      "a9ea91557622eb677227b297ed8454bc54cea8c81f7dc559cdf272ba594b7c7868aa104be5d1f4175e7736d348693c2dd15e6a666fe4463ba4b36e5f",
      "50da8bf42e14a3a8ed7e5d21589791ddb1c43987f2183a1d057deb51cd89b3ee54e435f9da23cab003f6178fde3ef31044b25819d4050370873f6c22",
      "4e6a09162e198d781a79653ee3de2fab9b1dccfae507551c199fd0c3156f25a9a7e2b63cd3caf27fddadaff7d6fb4e6358cd6a0d35f49b847c82c2a8",
      "dd29d29e4ff31b7cd518220b3e3c7598d7ba4de25efa96f833e7de114980d1cdea30fd81a3c58ed670e6b6e7e3bdd436ff84893dd7b25c3266be74db",
      "64a990b8cdb9724d6d46814a73964cb355731edd584d31884de24023acfaefb4a963b62f3c751a1de867f5fa9484582a7d613bbda7dbb9a8852e0aee",
      "97dc366cb44ec5371d4ca8f5deec72a2555b4e518d47b44206ccc955f94f3e1f6c844231a61ee5e8f9d1fce28e4669e6554773a83dfebe7c383660d6",
      "5e1930d6742265217608f6c4215113231549a454849cef3c290e499a4447825245d9ed1e99ff9cf889c677c7dcfdceb4f827a91070f6f56800bd021d",
      "8b9037ac1cb3508dda12afffcdafa30d1ca99dae4fbd113fdbc59ba6ac3b60de553c29b8d788455e58b30e4746f5c3811c6181ffdcb85fda3d21f8e0",
      "05ddab6354e9de5ecebfbd6ab2737c937b4f6fb7371cb84ee9bcd5a74058ea319662006e21a9a5e2ffd5148fec4323ab299edc4a220de7bfa4997545",
      "1dc12e07c36a7c293437b178081597315031e4042dec75d19132cc7888ef84c9cd94c39c69f82fbf1bc9179573393e7d1bc22c48e7abc0f861248d4d",
      "b9b0956fd38d254ecfe7314ff7461360c5953e62e8d9c23051cab8aa5676367e43f6d3978fb31b7fbf73c12a8b19d5a57a6db1ba46ca4c05cd52097f",
      "006db5f2ba5cd250ae518a55c728a5fddcd15927c510edd1a95c14e9edaccbf879c016e757ade60980b08e6221c4d8e4d69666b4b5a5a91c00b5425a",
      "b822aa1d50b2d10296a5cd4a55caefc90cbb6bd1c642fd2058a763d414ff3843a74ff569387bc6f8668f6f0c90f2bbf142de810ff8f905c842deb1ba",
      "54c368929bc0566cb4d7df74df6cfcfcf34ad97f0d3f37cb819f340c0786d8c17c9716b34725208425389e04ce4f87639094a4716a2f6904f142d229",
      "83cac63580b410f7dde2fe3beaf64d00b29204f4b88530ae2f722ea4f16a800dfb63a20675996f191594d09ad68a5deb91d653477b00e14077443dd4",
      "ed1df4bae16bfcfb3adcc03f1be13afe59075284d8967f3ec3b6086429952980f3b94aee80e1885b82fc1897e5b4083b1d8909ff4466eec66bc9679a",
      "add63fe2aff1c9304fa6656b925ce6713e6f015f0c1c4cc9d18a23f9cad15544e1fdf20379484a11039e7cbd1e004bf3f5ba81a369dddd4cd242f40b",
      "dddededeb66f37da597edd59ef76bb1da8e4b5c47678bd75af255f15bd375e0bddbbbccdee22afdbeab67aebad373591632cef3b98f6dad524bef67e",
      "d99cc27ab7ae92c90408ad57db3fbddded79ad51e41d401be3defad7571fbadfbc8e5becddbb9dd7ddae2af60a8b6d3414dbedbedd36c57ec2623d59",
      "0c67f14b4bca30f7c5e0363bf12f2fd0205472723313d4596abac1f46b4f0b0d0cdb560d75dd565d69d48d1ad4b56755f6888aac6bab19e7cff633d1",
      "b52132c80597badaad44ff890a55df8646fd44b067540fa5d9e6f9c0aaddec6dec69e6ec252d8a5bd04aae34243431aa53f79be7de54981c1e6b8597",
      "21850011577295b012e92f6a51bd4a09755271cc9ff3b328b08156c836abaf543d162d42cdedb50db067e5b956adaea447dd2608870c7c1373c6e02a",
      "063e3f8ec656170d6f6d6d416f4f3dada136a612d6aa37314bb7ceee2258301a2e04a83a26ac12e8bafb417e98a85c3d8fe4f17463e852bcfd275037",
      "bec75a5a3dded20831375522cafdd05a286c7a9052a6d60f70b33cd2e3469531a8a399bfc31fd46ba31f66e00afec9035b2994e87f19a1fd570ebf49",
      "01fd2ffc2568c9f7f8d3028ae89f9c7d922f154554bde62dadef5fb94c7cc4f2e460ab019e04a342d4179422a66051b1f7ba9835b94585ff290a5bfa",
      "e60bee034bf9ee59898ef4266085bb72054518a242046d5750fe9cac49b6a04f8f6eb167ab1d571a95851aec86969c568866b92a2361297b9a45cd6f",
      "f36d1359c37715aeefa5fb120e50f23b8fb45347fa435ed98c8f28e275d46ff61b8f9206557f12839c66d013f995dd78421b145548e55dad05975cbc",
      "4e692c8afe1bc8e646f828c167eed2d5eb3c9c9548264549d3de56375505bcaf340b53d72f6f3b303be7323420d02c59ac6ef36f5c7a064ad26b691d",
      "a86d32d2a418ef401e0514f9f23f198c6d3f42dec7e9b16ca7aad9a3024b35c3a2aaccb496bc5cff4643163c14023e38cca8273cab884f8946072f9a",
      "7b10dac5c8cb8f14093e3c88566fc7097ad7131f425e475e25646f301077952571e6420439f0716615b07df2926f9956645f690876ec4302ef40f821",
      "12fd0df4afb3ee798837d676ad65159f04b651f8741136d8895410fc8da95fbf57ca16d1d182c684e32481a3ff20edb52472bd1509a50476e2ab375e",
      "960c3cedae0da32fe1b05d55a10363e4d0506f3655b550ea4e3f8678902defb49fa4443308a419d541953cb4c34c592067c1ce660d721c80ec864bbe",
      "72d52f3c32e99a5103d4b9db13e825e77ac9ab8eb9ad8dc2a8db005cbc3db70d19db773210aff89a57e275ba580ee0072a74a2dab68980c64d197366",
      "e9ced1ac2b4008a06d3bc5344bc30ef4aeb52fb3d1dc1c55751c3c3bee89c49b0b7ca109af7bed935fdfcb943d1935c8687ea9267440a1ecaaf54435",
      "bfde75250178e23ca7985b18f4f253a13412026db368caba935b137330941940081e20f73ca9b0f5f60977432d4ab344876c35065396a0a3c9fa4cc8",
      "d5aacea9dcb61fed4857809570f0a4c9fa4705cddad370ac31550d3a609fa64fdc3735ecfac8cc97732cbffd37000736f4db4b77a50e75c276ea8b41",
      "330639b59e464db53ab369eb6568899d52f0f1aa56181d1ffd5a4540d17721a4c79c3b05beb9634b89396492709f667bc0070051037781f6a06f6fa5",
      "1da1dc05cee0d1b5a4fc07af59c6b08c152c67311bb3099b45673df68aadb337e76c1875d9084eb254eb186e8efa4354ed381b9e239332e25191f8f8",
      "c51ca59257a40042ed8fa319eaae9006cc90f40e7dd360b0d6a336b3085bcd236a29863fd939cb4507cbe341ac7ee1253c89b2c48f61b063218836c3",
      "9890b05ff05413943fea7ae8c0eec5f5d4cf1c85c951aa75064d14f5efb9e50286598b27dcba3f3c50141235c59a7756e145ad20f7eed3b828c801c5",
      "9e08780b7c5f91a14033ce53531290b870060f094badc3743207429f0b77ae2d0c5e02d87132695df2d6ac106af1b084dd7363a932896ea55f2cd85d",
      "0a3d348cfc59c513920d41b38a67986065567126532b6db5146cf67eea7661d527d1bfb87471235703a00900e102a0e46273d8bf4028c1014e6180e5",
      "d9050c79a4986213eee26c8a7e8ad626f82f83dfe935fe4eafcf116fe2142fa23c8aa3314057120d01a4b0e924c2c6d36804ad02d08dce12742a98b6",
      "e72b59fb6e2d83bf295c95f12af23177abe82168a56063fa9cafa27b20f8cc57a38d95a23f8b7043f3c108effab3b80330d4c9cf6125b9e537e8c69f",
      "1961188c7986a35d15639ee16857e5980d1cd1dad85adf30fcaead9c6b6be4de0bf5305ca60c9689af46e4bd8925f423bd66e9eaeaa3e99e77529674",
      "d27361225d26ec1a2e908b0bb2a0bab8200d5254ae9f1e6593f935dc36ffe06c2a7e4a7f4edf738636f979968cc27ff147833e78222d24af7c0ae5ab",
      "5c2c4aa3cd68810fd5e80c411211410fa747be52e5f4d2cda20fa30ffc26b76840d0ac61ec28743f34875f73fd1abf92acf215602285311d2dcb2fe4",
      "cf400a5f20955c4805b094c55acf4865cee0cb44b14b1eb1bb089a42818952fa12bea77c1912e813106d6fb6f23c9e2f4975ae55cf1b98d490fe0d7c",
      "ec0b6619656729740133b5cc4aa56a012bc4cc24b68d84ee5f1c65abbd7ebc19c1d0fa31acc538ca139f9fc5805849e10c7128327093cdb1f02a0955",
      "c7012c1d4a6092b31c3b337da0c0577ee40c658a7e02a99c0174adf5c42585113ea5104ceaee11ce468dbee5e5dc76aba91c786b4cd82ec8f2ca82e4",
      "24a9886a888e2dfeca4b3f9759e489d6d448135bdb420e0897244a076512e2ad82cfb01344d0d8b090d625d1a42ffcc02ef9e3878972a671064d9d23",
      "f21d2f4f8c7f8d18ebe688f6b1a72282eb63021d66e4070bbddc8e077e19c1618e62b82c38fc9ac0af4733c6ac798c7590c451a3fa1e05c1be89efd4",
      "92bd594e063e02fadd6abc8288050b071d0e4b0f1f41f8cac91fdbf963cc5f17f988a830b8db1ce01db115761e7460d0bd6534198299d905c6a68074",
      "fd8fa3ceacb52f12fb8075b504bf7db7295c7fdd0d9287a81742c22fc2e7175a4140d23a22f2b92c34c742afa0d05c169a8b426ff08d5bf7a5b73f15",
      "fd6511572b4721e00115c3d013983a87952d568ad57c25ef6bef5f31b4d9db84e3e1a34fe8b50cb0f22afa835ee318ba0da638f05151135b4c70c544",
      "8d0c70e24acc0047e62b31791e2d65475495a503d95128162843c5793de43d23f45b1af974913e3c78d925a2796f4947098374f89fe53e52a69886e2",
      "c4121fba37fb0e47ffb9e8f406afec16b47781973b6b4d81cc2de81e6fedb7f793947f2a93493b29de412eb1433cc6202e0c07697a1a278d14c9ff6f",
      "e991b1a14726821e9945fee4097a6452a347264fd0239345f4c8b84e8f0c2954299c82d1e6ac3f5217ed055db42318f2b0811eb9207a84fce1b20ba2",
      "47c60e3d9244b24140e6d824dcd14368cc1e227cafc225d459c73dc2fd4982180a75cffb7c027b4f76e1d556ae00d6b11d8067519d259b7e0a401e65",
      "76d32826b98fa3b302601ef09f9fae2501863af40b227ceee03cc029c8f17b8e8450702ebd874f2a044dec103444caac8a491329b32a274d14c6d73a",
      "79512437d34972350f79c26460b5137e8dfcd58e626092048890ac804d22b9c8612a0b840b5e00143611968153553a4de80d7f7774cdf7907891aac8",
      "6126d2dfa29f8911874b855d147f6d377698274c1cb2708fd1510ce384282214eb4992689c3cb239ccb31e6dce792498e21b1f1064f8d2572e084fe7",
      "54c050de738638fb5185c416be24cfd67a6fba6cede7ee393bc35f3f9363db3b18c271f8e3c64f6f7a1b3fb1e38b83bd8f87f8fdfac79f5eafb7377a",
      "afd65fbdeefdf473b5a9f56e77e3a7d7dd37ed8d57eb90bbd67bfdeae71f5fffd46dffb4f1f34febd043a540efcd4f3fbefaf1f58fedf58ddec61b38",
      "a94db3b034248ef63a3846a9f4718c6a57b8002b5cf9e89023edc0151c19faa0b796ad6488a2b3159198a468e75188220026be6afdd55ada590f3a5a",
      "1d03aa1641c7efad02039e410ef4b896884626d9b5af3d2d16acc7d77adda0eaaa12376985afe0d5b46893346d0ed382c975e450944ad971756aa91a",
      "b73db55c24f1bba98ff77b076617eb255b5f5b17438e71aa40b721173389dabdfe78b3f77a791986fed3a6960a4e82fe183001b0ad66b16234c5b317",
      "85d3a200355be0a200873359d41b90bc6b71df81c418d603207105780734b1b8ad1feb7d349b2bc379c90e783e449bbbf0ae6427d331cf93613cd189",
      "79f9c88e4ae136be28d93db2dba1b77b74f23edcd8f8f9b5a7e0090f2a34a05ddd91b7b6705cfa9765d47edd51dbbf7257b68f03d67ecdd62e4bf883",
      "6f8c93a4b1f9571beb3f3acdcf9b9aef11b8f6e020d00fd1e2be1c70062d3edb0090d36bf00f2c744539c92583715bd649aef6adea7bd42e86305b4a",
      "f0deffb8be3050a2628bf0da118806a929bc90e8bcf56b14b86494027c1643315dc9c5936f16f959d9de8de19a8c6005f137ac1eee4d74643edfbcfe",
      "299ae8cf9fbbdd9f7b1bd14c27e05a4793043f4f10dff3681f08f7267f5c187dd589d4c5ac0757015a0ddedc49fd73a14b90d21810bddca593c87d97",
      "a30b3df54c6bd9249d665b238cca628a2dd258704ca12a3eb916049b6d392aed1443ac747cc992ee45e90695c6a513ee40f68cd69610773edd851da8",
      "4786e651f37a514bd508ee9506cd6c1bd53b9a5556ac4d27f7e0b0768dbabd56508e7e62196ee232db7e1222db604487194cdc6fedb1538508f08dd3",
      "4f9184b67fa8b2c45d97a80b747b125246d24ab78e3939394673a36da1be12b62f8399691526651e4f4ef58532e02ca93ae2aafb09106f5d32188a5c",
      "c04007bb1dfbae2d38aa6da1e0b492845e8134b847427bad7d491ed62001a15b7fd8dee92cdd66d3392be5ecfb35e6e214f80240a4e83173d4123c14",
      "86a644ee216ed1d4281ed522ff3182be208d9b663706e6f587e2160887059613b68abf2b5d8fb9b5ec77302aa3365d6d113a2c90f9328bb808a1e261",
      "2da4edba4d016dd3d1c22925ee6c1e198f87e3266f3f5a88d6aa540954980677c3937367b92c570fc522871370dcd0be7e803efd06658861afcfce03",
      "ed776233313e9f2c2756fc5c2adba17e0fc56172da9771a0ab068b707e6b39c24423d04a30787ac9dbb918b6c260cefb19754a065abef4248e1bd33c",
      "9452005173e3a468e6c049bdc4933d57536b8f777047f73ab8c86bf8c7f22eab42d1f941551fad3284a01a6ed01d9e5a496645bbced4b20e78a89d8d",
      "0339ae571bc6631716bb33488c77f584e964b57602c4a44d4b84bea8714a03d14e983880280d6f220ea5ac429ca5da2babb304353fa28ef95dc04c88",
      "f78a7ab7ede9da19a04cd466bdbf2c2e53c8320d2516f6accc781b265ded79737199c69ee516a11a5ff6b4b69925f8c14b41aa3771f712a040ad8110",
      "4e74056f268e76b209b51a8e76722e488e453859cd6e5fd83d2fc6d0fad00a7f3275d3efd225db1a9c88456515f1a9b80d4f0cc0ea6e628e79b31d7b",
      "a558f3a5519e3f8fd77d2f9ddd5c523c3825276c1c6a1df79323c8ba0dbb3d257d3df8161dbb2fcccd052192a45fb32fbcd1dc238de88dc2f2a14752",
      "7c717f0099450fdfc02a0b3cc46b378d5f710903a89f163441e4134fa793b96fdb4a8a955ae0b8d09a88060362909a35bfab559c8907ff8bd766a397",
      "a626282fcf4de18638ca67e7fd05fb48c245321b20b56fe926735177627bd181982ce9299fc512b2c231d290c3249a251a67349e61bb55435b07c48a",
      "846257882e3e22fd21494941a38e074641bed59d302e227391c65d78c8178dc8f71b8fbfbbfaa881251a7e66d456ef4f0fbc4eff29b781b0f8276874",
      "f4925da28262932e7372610ffce4e2783faaa65dd63375dfc6ae439227aa62514f80463d3eada248904c2faa64873c218e0f5ebdce7daec02a69eba6",
      "07d66f3f0813272e9d7a34e478918d9205aebaa7b3e9563a1c677948af024c46dd76d286795614327ed652ef25bebc8591f7de7081cd81e5a100cbf8",
      "5e3224a77cbae60945e17e595d11b15baa4f35765c89f485940f96f894a395a00881a85ea41def097b37d718056879d9db3b784fb10b812dbe467743",
      "40d61976031d705083047305059561152f1f6605a1356fa942cd58d9408627f677b454a37dac5c74b0b5280f684c9c9b585b724757631113b10ac413",
      "16a4f4d7c279b4646b2dbfe09cde9746d9b008dac665144e1e034edbf35ebce08a404ea3e48caf7a6894eaa13ac28d9f46b55b3b1d9ca52c3d0f5314",
      "93dde8ed15912313199f5d80273a41c329a82f58dbd48a6245efdefdb2c949948804bfe6ad72f22d4d8abbb2d0c38307f46e462aacc2c4110ac37252",
      "905d7cc451868c762646e7c5071de1c193a5566d61919756eb293fa5aa8e81dd9bebe67879e4f56da18fac9b6bf49055e44374e1250c2b258c379da2",
      "4b40c76592c66e6ce5b372d53ba674a8e69d4bc99c9549c9a8d1460f7549344aea2845c254e8c935c6eff634bdf618fed2cd3bf96beb77a288d85d3b",
      "5ba498061078c2b3f5d7ec153ec2e9cd0fcf7aeb94e460b41e5bdb7855c369bd1fd9dafa9b73d99968f0550f6b3fb76e5e512292374fd117008037f1",
      "353f8ad1aac8b73f231599a7e4436a1213d1e6c3af1832cbf20f0f76ed6075645f9ed6b09ccb9382bee7c9b43a5ca332e43ee32110c196f1a108f7c4",
      "81f07ad42137b16ce7ffcef2c9bf7dffec07ef1c2880f66af0efdebf83ef3b6c9dca5309bfbd12589bf76fdc1c28d113f186ece9d6e9b19a1da5f4c9",
      "26c0008d443de6e8ade213cb14fbf52ee162bd261be3355a263451a8e6ec890c85a69c9654dc0715554a0922f50a62ac3240ecb012bade5f339ecf4f",
      "a40b39ff8749927e3943e3dfef352a690f8bc23bff21c05809988341da0494f85d2693a49db253c5182f879e8737f4340164f042db265c7b346f52f6",
      "498b3d0160c9f6052e6ebf12c55acb7b740a5125bfe18a96843a59251f85b4f75a553dd4562bcaa08549dd679d739473cc644ea24ec168e876e26e3a",
      "92f25c19e5d5271d715645da7a3c467a595b84dac8afaefe3f1aba0efacaea5b60479031a94f4fb01694d57549a48ad682944b1d689211fe392b4a37",
      "945c35400af64be107c97944ea828feda6159a3999723ed2d11f9acb1c0164c2116064eec2c53c03d236aa066764b195a842bb42e2859f93ae9331e5",
      "8803151f3227a52727c78a5219f46313f5a940b926d211e6513a167a55ac68df056bf223e8f8a2d13b9db2e66bc15d2c54b3740dfc1035f0974e61d5",
      "3ee6a8eea1fa98eb3ee63aa5d687a9811fba8fb94e09da3733bccb2673a07090c493c1fa32762f430721855e3fba4a55fd422e612d5fdb53a8127fa8",
      "7d5bd4968946479608962113dcf280388b32521e7a3508b62f13e922021f2ea4159b3a84f5039c4d46829d7181cde272dc53461a2d47480118715325",
      "a35ac5044292063875931c11ec51208485e8b87a025055d4af2ec8ff64895e747a05490307d83e76993cd74d81c23188daf2f21f244042f77e32ce78",
      "265468d4cfea3ee86fe12d56459282013719b4480ce9d2530d6be348fe4c170bb68bc265568d354a92dade2451d64c98924ace45c22c1f1ce8f1e80b",
      "9f5f66714e71544ae11ed063f1a40cbd03ead663dff6447c030affda65d9341e26e53cecb13c29f861fa015fd731860a7d8a52ebafbb4c3cbd8bb18b",
      "9777b143f446ed99df5ef3037c8f49683a4cdfa187311c9f3e8756b6f55362def00c7a7f0dbcbb8db4c35e77a181ae66df2de935edffad5f3e1107c5",
      "79faad27d1a3a07e2cc055f8c38e0cc816d470427d3a81d1282c9a6d228bd32139804c506e099e92685a21e3ed4f135c5ebd4aaaa11b12495996e93a",
      "15836e1b86559945b3541e0f5f324da93a604d54260949892f044e126e6ad40239250badf561e81126e7052fedd44729d0c463b790e4106041c2b86a",
      "c90a621225fb0db509a8dcfdaa85a563f71a05849c897aa13384c012dc3e136c5b7b59b30aabd711033192c97a4e8f02110975fda47ccb2eed3ec4f8",
      "4f43aebe24ac7b0cef04717d59797639d59da50753f3405ed948318d05aec4fab5a2d6030e455836cf394d91ccad7d0aa463284b5806d4835f064e0c",
      "91eac25b6c8c795f75edd4bd55bf7eb2069e0e9b1a7ae8eb5e1a779318c2121b9a89a1f622e0cf24f213fd362a26ecd767ef9cc394d081f04049d62f",
      "c2f5934c0998121b264a6c4885e05e8122f02f499a90c7a190e91490836e21d1547c29e25d795def39df961a1722a8a112a17529e9f078be08594217",
      "9574d66dc9ba45c4128ce2a2e6a98fd66343b40673252d2f933f4c7280e919f246e66ac77cd6f24b8c2581813ec8a3f952af6f3640a4eb2d70711d16",
      "eec2b642b65c3cb1aa9e094a4025692de455bdd973fdc41d8a649f84753a8e226a85b95ecb6d38b1ceada21f080e602059a591b26d2ef7a606e5bcb5",
      "ea4415fc9d156fd84ee467ff97f6f38b35e7e736f6b4be3e0b94dc9a1612c1554a92ed0daef3010a1a4eddd57337dbc8a49d687e26641b1de33fecd9",
      "956e132657d2d13a5fdc1b64e7b3f03e61b595b6228e343c2f2ec034523eac7ab4ba5f257aba39c2b0890c4fca820a257fe4b70bb0b20c78ce540074",
      "8da5fb152cadf0b381f8fa65e15e7866ebe980dad0a4f0b755469dd9261d4b0768a689b9776a8419ab1076aa609db6abd06f5a5e53c910fe2812a580",
      "d74c2a0af701b054129b3c4b0c4864d4a49851454b8e57ad5afb4d17a46abd6f43d6f693f0be5d8577e871d173a9a38565c5df71b0149d0481d61a1e",
      "d51b1b41473b166669782ed521043092af800225bcb5a92c4d03b5d543c1e006b860f51184377e1775c651066e5e0e7419f1a94b9164652f2d808470",
      "0fcefd54b058a7d914df84803895096fc95be431beeb84a59148f1805edf2f284eaf799c780191690524d6f54443a7ce8bc6df69ca790b514f39578d",
      "0c7351e6a8d20234ce309b4037de7f6d6cbc797375e5b15b7abc0a372c3e988c75e269e8111af1e8fb1f807574c2282ec6c2629714cef15312edf48d",
      "ae21918dc5bfdbd49d4e56c0df5ea7cfe319f2e7e8fb274315870a2bbf48915dabd436e0e01ce6cd7314a9e3fd7d2cbfa43fb52aeb5bad2408027adc",
      "7119360c5a6fbe5461001553f6a9e098a68ab81bed5a24acbe7d3acc44b515e9e3c41de7d33a1555d79aaa318b2c735ba7969a226a0b502267bd87a4",
      "896cbd6d8de3e2f0363dcab329cfcbb97cde629e80312f70f1a3121babe1bf44b7c31a7775d8767567595ea0f7f17cb358bb6a9df01c57864f600a7d",
      "d6814259085663ce523863c7155bad7dbf61431cfb87b658f1ce7ad80d56ab1b6cd087f682d42559db75f2ff56f7eddd6d1bc99effefa79078f76a00",
      "0ba44925994c48215c3fe2c413cbf6b19538898e8e0f44c2528f298001205bb4a5efbef5e8371a24e5ccdd397b3263118d7ea3bbbaaabaea57e9bbc0",
      "2d306fe7618201ae80a67fa9e2492d3caac5a5f89cf66fd343f096fa523dc49dd426af784e36f10af6a0f50cb89da586fdeafc556ae6a9bdd1e53491",
      "470de553de38fca4e6ec9d7d17ad2a0a5871e9c12944263cd97de73ebd9609386d7b354212a2039a43eab0b572a64f070fe2c7df115d505bfc38db07",
      "d7417ad2ec8370b92f54e8a8e5353792b25fa7350a73d81b153fa5e35d4ad301d9dec101b2998001a5f2aed1c384f791a8668af262d0a0cba5cb2c06",
      "56840a4464932d1a1e33f0e8c85b47eec089ccc89b35fa4061c7b0161c12cd457ce84cfe7e78da1547f241a4e76293f135da355866d76d5ba27c4a2e",
      "8522f92ce9507e1b8f8374467a45b4f75fdcb69ae2c9dee1f73bb3ac40f3a9b37c076ae829f1e9f2d57a92b5810698f277dff0b2e83a0bc4f4c4d916",
      "fea6b01e4f6d1fd53ab2776beb2ea9634358d1a0379592db25b6a95588f404a98a15641351181ce24fce9a959794350878a543265ca4996546574f2e",
      "f4424ed3b436cb7a1a158edbb5be09e779bf0f59d11db981ead4717d52ed438f4e51f4b013fb2ab191287696f11a6a2d330b13208fa9c3b2e96c06b2",
      "395f1ee3aff25e11f7b55b72054fc64959c0537cdfe4356ff189dec6d0dd5dde0105f05b88a8817688455aded719ad11dfab6287ca813865aee25b77",
      "f9de19ae9a990ec739dd9debe91049dd2f4ee3c1b5bb1e15d24b3c8ecaf4c29a910b5dd4d536ebbe9c38df05212fdc6e67db1e39ce49cc2b157b1e28",
      "5dc66807103ca824615b85f9a5fab22c9b8b2719da19813455948f1662d96d6c1b629b3423c45a1bcd92d41b78129b875a43919cea7daa24ea1fc247",
      "ceaed392347aba55d01366e63a9ce3d8738b24c122ad8482af52009b8d899c7b581356972a5ad99910c6294b47b0c32b953f3bbc20b027a21a8be42a",
      "99a538ad8bb43ac9fa234418ab10016a77184f66873979adceb0279827b932c6d58294f59a5018d7ef3ca68033f015248445172ebaf24a14ca1c8b0d",
      "e25e5f644bdb8a6869d3a6b8eded797455373b182c9a48ec4e53626c891d967a77ae6a4486d17d89e24e0375e793333f207d913c16db029c26c3d470",
      "e7f970d5d8ed895dafda04edd5f4b6633ddaa5d984de39c96b1355fe435e399bc1e95860bc4f9d7d5fc753e7d161f7117696f920ab918ee57b720a6b",
      "06f17c24b49c85bb5612ee9a9846f949714a90378899169a1fcc1003db82f9c203c40c3af8c96df7e1c88cabe4505044c500d58923593ee331bad391",
      "582512ff2b8afa57a04e7469ded80f4a82cc3ebe547c7313a48cdbf2f0010e5dc28924dabdd46a4df7c063db4f5a39c9b4aacdc13b39b26bc99c9c52",
      "7fddd96a31a51ad993c0972c28dd0f0484a5fcffde4535af168c4f8b4b04d64545eba2341fbb4b34a2f5a24f9292d680dc68657c6b4119995ac35f1a",
      "21fc18d8e72dc3120acf86c011a43c6161a26919515a6b3de9cfe03e3bd2451c0b2ff0161f7cb15da5b5f80c4293f6adb3583dca0e133c04623f4c2e",
      "ec826aca2b20fb953c26a4eb7569e73ba91ccfebfe887caf61120b01539c9f26f00fa23fe177461359c480837322c57f6e6e60f0f8833f4281b05b49",
      "71323a45175f2ab6b797036b25fa07c05cd91947a70816b34f2b4b4137b53f838ea9eb0cd7d383db7c44d2e56f7a224e53c4924460ce7ca3701ae9ed",
      "27d786b9c977faea29a5a4b9fb1d055904c5d4626c97ac194480eda012d602ed5a91da5814d6a372d30d321902980c21570f137404f9ca1c5643e805",
      "8408bbb8f3cab46008b1089566c85ba3173ec26965c498e0ce3d4c2f147aa2020504de0f58c695181026569a09621faf45baea56b7dd9511bbb9d9f5",
      "8f36cda26dcdb9fceffcff35e7b2eef0e51db2729c27dcecaede2a37fb438ee7e0103d8d7cc2cd2983fccfab6c8144334720533cf29625d99476732d",
      "6e5f4c26b71f3e076261ab336f73e23c9e6e666a22af4aa4471bf89a58f236773e067c32c4aae444a1fae724a528bd997bde92e98b3a5f816bfe4f9e",
      "2842f169ed9343716c401dfea9163c1f1988804b87818e2e64930e4673fd623298501c92355a37e7205440e2644b830094ce44e9dd1ea47c1abc9519",
      "92619285e5ac4ccb59356402fa77a12f63b5b46593c01a8fdc14016e459a331e790ec72cb02eabef41acdf451738fc817be9fa90405605424213d6a8",
      "2023f842feda1704870ac35b687e777173d3daec66aaec3d465329a5f08f229d89cdd80e6b001d4cc0b6c759939190217f777d9e5402bc3483777906",
      "7dcf89792af9606766a4f4515f2254c2a0c3ffe03c2f294226c64329d4d30a7fabdaf0f7ac2cabb9286005d5b1d7430f204001110702b446bb6867bd",
      "804d7573a37eb1bf187cf2333cb8aa389e46b56a39fd03d156937a20e9918c78951acb3d29df29f51db0bd09866ffb219b5d3ce13a10f8d84980566a",
      "6319c23ee9b5f60437758514160632641a466fd07a44e234c6060d8721e81a6f248e8605aa91d79cfe15e76a21174a1cba6dc9bb1025da12aedf4eae",
      "80aed421a3de79cbb5d10a52f8543df5ce5287931e038dd3e497232025fd55f2d88d56701682b64d7bf203998834d346afc731aa54ab69652f44bec8",
      "bf4042bec093158414a5004618892b4ea302b5d20adfdc3c23e7e5dd6c6f6fb78a2dafe349fd5134b38ba8a2a6618b6675dea3bdde1bcb6c2f31b000",
      "deec5c1142633ca13c47e83723332a1257a30ed0552051292461c905d36e555b92bb309633115dc8aa112cf8357b868e4d5b76aa2a971e632404bb00",
      "c2dbf32c0ec7a3e48a81a357028322cbda254ebc5db54ef2ea55e9a6d2d1f840557a6d57faa3fc608fcac5829116ed59a92c72e34f102d8719d281cf",
      "faabdbf991a7c5a6c77aa1201e24deedc3cb315dfccb875beacdde9e9ce89956aa79d32bebe9eaaaa27fc18eceb1a3260f01bac793b96e741e68546e",
      "fe718b537d5a7c40edc60e4cde3f5fbf782e41d51045cde041bfb450b7d57e9f12227c4cce2797027d6dd0519dfd21eaa705c2772a52032f2c70e967",
      "a28d577b82e2e70922359c1cd8210d8e4520b6086e3a3ea635975b834852537c867c8a65f048066616fb1709dc7731a521f6118ba77ad597a6b123e1",
      "6800351077c41abf18cd71a72742c2f1a2812e45df6cd44f780d3f4fc7e12c1690f7276b3ad5b8d68d4a763a9f7e92437bcaa3c1eb061c24d4333e92",
      "af0abda177115a606f6f78586a1e4e5684310514e83e4c89e9d823e105dd934b6c2a4f10f998984d82778f7052da78e77f58dfb74d576f6e024bdfa2",
      "bae3cfddfbec330a708a26dfde1e37e9e7a694eb3674643e62f57322ab643a99d854fc482aa81de73b0a70634e8cd7eea4e08afdc849b79742432b1e",
      "03b3f0c17d3c771f57e6b1a3dbac5fddf5a5a6497038b01c986ef6645c61437e9d217e129ede13e9e730c1e87b1cc927b9fed28ea17b1e1c6501910c",
      "e4a256bb627a30966d0f93a66b54c21d9522fece909a34674c41eeff95d37f73243a6c8a0eb0ba995f51713d067a32f0a252ad3c871d85d6bd45669d",
      "c84e9f055ad40466d76245a41cfa4e318ff6936ebe2399768f1c9bdd07c32fda3303ad6238f5d081890538bad23a8ece4c0d5b253b53055fba5053a8",
      "fa178fc3db1e0a13bd8879fb733989ff55188a03fcacaa32c64917533dede6901e17f28c0e0c0bc8944b58ac77896a645cdc6a61eaa8495f8be487a0",
      "0da9311025774acf5c73e4410b25391eb42f18a998504094e72554a0115a10302159631222b590d5425f38c82ba3739c7129d3896ecb4e4297d0f008",
      "64fc8e299167d3b8c9b1c4cbde61168f357fb15d3c1656edac755fa18cae656ad8f6f438502238826d7c3da88ebf641f4f55acb7905f0b0aa62a3470",
      "8faab87e75274b525288bfcbdb93b49dc528153f0f1687ae7560f5b497b39c15bb1a46feb12af3efd29ae0dd6e6af9ed3b56c3c678b9c3b154064634",
      "8ea55cc8772ce5d45bd7f8d2f5848dd81305fd5ad29613ec96f07cae3f8e0f562b9d71cc70b6bcf7dec68912275fb9c8e0eff61c292f40fd1d9537a0",
      "3647a44ad4e18339c62f25a2d3e428b2a903e6eb93124a7bbab973e964d7c17cc91bd2dbcc9a9a72a1f03bbcba8704c2bf213084f48a52745c47f98c",
      "58ce6966909ca93ae80824f7b80e22ec2acb5b8945ffa220c14ae6a62cbdbf849c96ffbb90d3028bcaa542664d99af96d22e046e9337a34a86f193cb",
      "a7e373b758ebcbe5c23b2c15942d059395db254e020e5fbcae8ddb97bc62b75775dbefeb52149333872625ac726a1ba1db57009c33778c61d79af14a",
      "7b012f3e5072a792148788a30ec5c45353223001513cf903cf7f1c4d9cb8d86ba2037b4d48ecb5755e568fdc79099c55dd3e7e0e3f6147ab0a2d2ed4",
      "6d75bd73097dcb3dd0cd4e981bcefe6af588addbe57e9b041cc93c366cd2e8c55f512775975aa7923c873aae0f03547660ddeda184f25ea43f042e38",
      "d1497509bd419f8245592e19b5235f3ea897409b5e219a04265d5e01a9434612f3d64f0bf442a27bd1b504fad7a78f7f78f105241a6d35cb1edd58fd",
      "ff43a7f339acd5ac9b5a3776e448da5ff20cac1fae8e7956a25e5d5e5533645149fce19899c28e9959a83086e5292e09659f83eb65a8734e8bf109ad",
      "9853b257f93c37abcac0d1639113fdf3d43f1fbc45b0a57f91dc46898c2df74490a3914a1ee8d4b48737eb7486a91598eeeeb6dca7311db2e0c2f45f",
      "631a92255c98fe3b4a44eb38b358fd2cd62b1dc28d83469ab5aa54604ae75921ca9efa4293ca3d8d48dde7880815ec3b25cabd7076dfba4de36c0b55",
      "f9c4c16bdb661b345fb20d9ace6dd0b4b641e36e0335d00761882057382dd93f5182e05ae227e3fa30414ff0c21159c4355229d05d24997650f50f37",
      "3764eb3e8ddace08fada311e47be0f1e7fd7d497cd64278cb50e3fa7a14c28fac1e77f11c22e899a542392ab5daa58796e9a6d3fbcb03018683e180f",
      "880c7db771540ce1c2ab4a9af2fc7cb13ea60f13666a2dc2f8d3ad31d89307034aed046d02992fb30ab880bc30b23ece53c46146b4d3e1766048a138",
      "38ea761aaf568d7a01460c7c4d44514dde0185361be451e41543cf64c2ab3f169708efe0e063c8b42d1402babe160acfe61e8c6419c777722b6587a9",
      "e30b151ea6db9df84e6effa7517006039396c24e57739a45c732ec44d29aae83e110b6e571c73c86c7ffe94ee35faf5471a6e0afa235ad778974ccff",
      "945d4a2d24e889075b675c7291a19324b1ab57921051b75a799b70661fae697b975a3d6521b4a380a1a3ca2e99810fa2166762816a2904149ae745cf",
      "512fc80178e68f664faf9bbe756df50253bcb532c66dafa595715fff45f58ca89164064c0e778354de3d39be48f1e66dbb2f51beb5aab0c9bf1b4a09",
      "435e0be7c0886d0348b6b512f6e13e1326949a745e718f50362432f13e0a345b2b4ef253aef833148307191c163b204c63ce5126e86d83c58511a5e2",
      "2615b660c5566e2a0f130295475de6511e1aac1554a0edfcf80b019fa0952664dae15e48cb4dae0951e265d5cac2df91f1ec6fa2a952480ba9373ac9",
      "6fce5e0b9aa42aeec6751ae0c4e7e51c65b9804d8c936dea3c45f64cf3791d8f9d1c84af5f2b430e6565133740ff6165fd747cf42ccd795ef14b4f1a",
      "5cfb740c637f60622671e3a2520f2890333d2072bd7d70e7f1ad25bccb0ef0aca00142e2ede9a067e4c4a676f99d7cafd34bd7ac9279e338b1622a30",
      "a6481462f5a77fb48e408a63120934e932fe8f9ca4a56e9d9ba155528af82efc978ccc4fae578c61e4bd7f835aa0fb07f17ea37d0b7dba7bc6f5e752",
      "4f14ccb3c0668446ebd7036e131c12185099f1d00d2948cc9e54ab78b7565690766d6e61539486608c711f3651110fccb9c98e6668af603693404b05",
      "8a5e75b7a645b8e929082125e2c396868366bbb9922eb50b8c86c13d93f118edfea9ce3c105d914909dd8651508dc4f5ed6972995dd3a71b7f351c26",
      "97a2e0876f86f8e22746a021e2aa4150871e08aa42e8b173bd6c63f5045e8f4fbe49be39250dd3d3e257383ce9ae12a58b8784b9a7daa288d6849083",
      "3f5e143fd4b36c99ff9caf38cdbeb0ec14b976bf50e64283590d12e9c1400e74e75ac152654e62f5e91747a07c601bd1725fbd48036da9c72d83af3d",
      "4c0849aca819acb39770c4191977d3116935ff2587bcb628daf23aa5ed05fb4e28e14a568578b270ba936b482f7990774a0fde8038c79a31d157ffc2",
      "4105cade7154881feb0f6b0d7be80c4ee733a35366265158512d57f8236c6ddaf9661cf0fad0d8db320ff9360c54cf5323b97b4294d97d5440823ea7",
      "1e57ccbe17468c0ee87795b04553dd6b9f2e3a2a45c3d160c2da2582eb840c5dfa2a752e7dacf0d8762aedcbe3ba2fdff5d008b1c5a0844af412a09f",
      "4fd13637799847ad2271f21a35c99472dd0005b9a28520175123968f4263ecc38bbe1e7d2fb10bf8f97a81ba5aba394d17e143450651452763a5991c",
      "17a6f6257229aab8d6439bb6a04fb34a64fd4576962f201751b81dfeac9815a36da4bdff925b4b586cd8df0e6b386c76a8340b7169afa9aef2def77b",
      "ff7530fa6672781fdf7fffb784d04ce596b249c30bcda24a25938a3dac39afaee5176244e54d14fa1faa5b2a90f7e0e78568f2d7cb6c86602625ae16",
      "bcb5e253918e3f5c663a668048fc689c3a4fa4b25128513f74a6cca62fc8f647ccf8b81dc0677955d633dde0933756228dd538bf494a6fc3f5eb5985",
      "1644f3de2442c3d743747ad6f18fb8e1a378fc290676a46c09c6d4557fc3da33b2e11af5afc2626a6a6131b8933637ab38d8381c70a3854e829e3a93",
      "100cab1d609b00579185b6c56795062ba834d10ce816222b0ae05da6ad94747734b6119b31a44e56d5a8638a96795b39c531ad9845827d351ac63737",
      "c336d97416c47eae1d86dc6f97c052947eda6db63de9e7fd20ab0f74810251bc6ccbeb8ce76c23c968c2c4df5370e4b0d0ec4af62ea6e0636b3248ce",
      "f1e6a6dc5897c54952f69a5104f9529abc91f1fe39db2f07d7dfd7e4bf54d18d74d687274c45c275ddaffa05ba3f65e9305659fa05bf5cede7906f05",
      "a557f02aa37beb1c4aaf309532f433c8bb3a1caad7f8142751757393b5b0ecdd53b5bd548686c791e92a3006471a3981957b1ac7f1061968e7327219",
      "21979f75c130a7eb5e4631dfc07094f6a724525de6d5792e4de6a3cf2d0e03af7c41e230e20f729094a1652ce85af31821297a2812cea22f01fc6b0d",
      "bfc26dee5098f746a65c1e235ca12bac69b4f730de8655531aeeb78d116f4e774af8292be68bbcaa8159cfe7ea9215d1c16716fba6a70b236300a356",
      "d73adc11c8362f3141465e1a9b23910b98440a404489b7dd3d4815cb7b9b5c1581718706ad57348260ffa7ba3db2df39b1e7436b6dcd305cbda1bef6",
      "6671ce91b71dbc7fe7b6aa51a28154d76817542913eaf5ebde63855671b8a77675cec295777377ae848b397b4ad4bc85d7ea920355b1fa39624353aa",
      "62e39542a0164b63d158770b5b0c0ced6d829f9c34c093408bac7f7b9547caff1ccf32bc426de8a62709ac0090c64000f104c26970ade4e6ea52868c",
      "90a17ce278dca185577469ea6d8b281e77d605c3d63b257017aa6755297855410a09a4b6a35370f41539b7946483972d483ec5a0088f8067565367f5",
      "245eaf559240c7be5e69789acce1706320733ae17a09486497900d560c6cebba0182b2a2eb7f853bfcdd17293f7c83e390ed983e6a656fa5a6433ead",
      "d32804c244ebdb8924a07358d30039afff758d88ac53e935b61843386ef436c3083622c7f16f5081b846c57a792071eed05b6cab82901dbea39661db",
      "4827019d42a7aa2370051912c465877b9bb28bb915d45316827e4bb0e975527297c4a4633978f719c6abde92680a7fa0c0b60b7dade10b08f685148a",
      "13629d38d1c48987e3a90908b0f7852d8b82085438525172d1756fb208dc9b8814a67b8986883588c86975ff20c9e2718f2f279ce4613c865f3d9651",
      "f90d248f7b156166d3f370dcc31b0bf9720cf2c421081bd3a84e652eaaa44e395752ed1fdc8b2e40225980d01127d9fd03bcc6b3803c2fd18795712e",
      "5030bbe03fc0f57c8a8af697e7263a5e528b1def70023a5e9d2981f428f4b6b75f63d0c28283186cb80cbbc3d597130f43c815d8e9d4b2d69b452fce",
      "9059cdff68708f8e41fc3501ce0942a0e0840d1c021ea3d3f5855cd1ee957b5b864795ccbdbdd4f674a3d4d6aed27217da209485bae375462e47395d",
      "a2964598455537db85559956e1c89269785076165ba8c37cb2262b444bedd7da3ecc1c66b4cb3a43c707d79f2272c5b5d68c744d47b4bebbfae6c3fe",
      "46ada9b1c5acaeaa4207c564570f5756e58bc168dbdd03da0c4729eaa545fab92528cab67c515126df6e9aefa94082991a36564d1d1c473a9253fb75",
      "a2de5e3569bb2f18b75ef3224e31cf32129569145ee599a881378091c3deb37bd47adfb580786313eca9b1ab6d4d07973ec94f511bd73df3e9aeb246",
      "5db72d032ba9436a569d0dc9cdfa9d23391ba02d33504f5cd6e694d2d06caa89e9834a648ff37a064c503e7fb87a51585b67ec3a05cbf5e741d17457",
      "226f17d6932ecbfa45f63e0ecd9975654162f7d69bd6ad2528b35bf46d537daea42e0b6e94d53b3ad321af6f3934e239bd05dfb20a6c7df156912ff8",
      "e05d7568809fae0c0107b1903993d56bd4329a27b4f082b54ef78532a85bcb0cb2b57d9236f997f75e7cf178b6b8aa7a6dc26daec63ad7b7c7dedc75",
      "289412b8199cab86ce563d7fdf1bd145cc63a9b509911d5b6fe3ad3edb4e144f4b2f92ac938617d6022deb74a405a62cb3fcc9223b57ba1b3b0dc349",
      "4ab3ac58dea12386a50e776c7fb0dc2f3a4a72e70440060f1d083abe6c58e9e4945f43fea74a9d239d2062ad0fea9cd0541599ac3f573c0590176bf4",
      "5287756acb696e41d74172e6e475b97d73c7d715a220f6cf10a3d7ca490d3517e178d3780b343e191d24a383d3e4a2b9a4105767e718a090ac8e2c13",
      "2125cf8074dec7a23d20fa3a5c69e85a736faff7f8e9af1210861db3a6cd780e5b1bf7880cb6a9760c09fdb107d3a831a307d835fb4855742fba2425",
      "a5678948d9519b685db0ef8e5045496fa6fc674c17c83458b236c4386bf464fc33cfb2d9fb73b2d9536249dacfd95b73a7b70f3fd941d3c25fc3c920",
      "77ff1a5d87689e48fb6287160d107f8a0d09db61021fea31c33ba56f1951f579d013a8118b9c3edec1377fb7a29cb158f9e6222f9e0213333e83f578",
      "0639ad7414ddd0800cb63227126ffa215b8c0f86430366c1ce8fbc08203f097c64d846bfa4f3013c3ecfd003c04914453bb128df54d9921c1f59f909",
      "9d62cda7e3b684777f0fafdebd03f27bd0e1aca4410f9b35001950506f3ecda42ff20ff982c111a5c5087c25f348e6e778ed8847484734363a1e715c",
      "cfc4a568ba63a45991dbf0c583c5e2181b8b4c4052dbbc5b65f32af68572b6d4367dc7ec29cfef5d8cd4a32e2b7566f6ae9a52fa8e2bfb8cf8ce46ec",
      "d1f99d9a1085ef9b4145b6f197f80bd01d9d681dff462807513f2b33b2d7ecf47fe1f7a190798138d7f6f1de5a5a9101a2cf2e97b83aa296873ebaa5",
      "eb18c66a11190943ad2a7750cf68e3185f783710f87a970f74ec0079c6f6ed108c52075540e7c3781cb80d131fba831225aba1ef4d8f50415fcc6f83",
      "8a71971a6a218c4b798fe9bf22afce245097229692814cd8ea2fb54b863d49377ba9f009811f33b058d69e9a2c58a803617d0856756ee8f9b24f5590",
      "c5c72664c146dc004b5df897b1037c63f6f5f8010e1509c7b650da77e3653843e6a0ca8b040ef126ea63e816f87f2c63b6081f9e314fd939dcea4e82",
      "297afb58c347c1a548d1be7c9f10a2eb277804e514ad280a518d621f7ba05d153dea11af437cb021ce81733f1b887cd17c6790202c8d6de8466f220d",
      "bdd27dfcd28fa10d9829389731c63cf9fbe08ad8b14e48697caf11e729f1a4399d8048755555c402a31b397aeec3684b63fe370261a1af5ec5f7d133",
      "31791481c8b9403bbaf270348596876348613fc22942a2e89d0d43fff38af603c2f1aa3c29dda1155a6429ca97d555a1b5bf4b7c50041163df28b116",
      "9d2e9f54ec83ed26a4d751e040501222dd069bae8caf1297bb58b32b3481695f6e69669a1d93d65fb6c5dd27ccd6585f77f1b8b5eec288f0775d4ea8",
      "d3c2c756974c22fa1e19b260813698f5c51c592cd2e71417cfa8e2f8c5893885b5a2f7ad8e092008092f0a64756847ce2e37d919b0e27d5d7559fc62",
      "c646ee3dd171a02adf071756d4033e41ad9a98ed533525f31c3e69bee355c63b4eddffc9e4e6d48465cb961ab617964cd4ca092c6a0c3d6a2f1efc00",
      "c64cbac3eebbed215cb4670ad25836055143052d2bad506665286a9a8ec80ea53fd95c43cdbcc6719515357cf54b208ba57da3c94f921d49ae22ea90",
      "6343ac26f8119d763cc145ecac8eb4e0bd697d4ddc9cce47e104ab124c3024a23370846d0c9b0678285cdac5f7a1457f73531cfac6cc74c2877936e3",
      "7bd7a2bada4c5dd15b98ee1c3f219c5492ee4eb62da8e9f4ae1868a7ef1c2b22506bbdd4b17a69e800522e08c7081af209ffdfff26d6711729db23b9",
      "29fd8cfb0770fe05fbe5764a0ec654aab619d9cd10b693b7f1c2c7fcba56f25339c0c1a79440815a6d491429fba30423da845b713b9db88ced867a24",
      "e9f3484f63931e97be481289aa9e109d695410d89658d021303a9fbbc3e3aee1b3fcdda2c40f7aff00f533564a7e9fa2400a60644a6942bddf200784",
      "81674a98746135fe48e2a7ff9cafa212a18a9d555a6a30d7125d05253f50aa258f5c02bd306c867905f4e750988f6b2d623598580f572ddbc0800d7a",
      "f3c1bd06b840f877ffc089d957c39b7c52c39b1cde5870363cf6120302546954e1c8f7471d6347fbf0ca197b0563c7f00272d0953d687aa1075dd983",
      "86260e0b77d87a53a2013f76a1600e468b4ffe9d3d3b3a896276818ad6778bd571695dc59332c4217f9a82078862d26c301170ab9537ffca12007184",
      "0745c9d49c437d4921b60b78ced3116a568320d68c060a66e8d04b9a7acf63bbacadd242dc7f27e1b0997a29631958c79b5dd76fd572fda5cd119696",
      "e41102cbbce368e928c6270c140b9e3c532e35f6550365cb6aa8a523842a5bba8209b2deb5b9cf5464a5547e1e6765d548a5875a9c769a8f4221550b",
      "b60eeec70a03f4257ab4a5cbd3523464756e38fcbe2311a4bbbb22ee60486a19b3a2fdc2b39ce8e6593b381d8ff74b3ada0907e4eb025414898254cc",
      "b5f9e9e01211a0978bd5c315f95bebe880c60c07810cf38f36d7864d69ae6d7236c88ad557f329622ba24c9623f6fc1ff2b78657c4afb1c63c89bb64",
      "e2672536f7818a062d852b0545a45db38c0451ba43a73ebf29abc5dc018854d9e349a9f506e78bf22ca3b3ee555668a82f1540e8b8d46fe8d0310e90",
      "bfc15ec7bfcf7099ee7a01995057bdb777621f7e260cf030d1254f86a788d23fb8be8f1181128e459c8b4557f691c9be8a4faddefcae7a933577ec8d",
      "2e8711178654fd6a4d6f4cf6919d1d7a43ccb4d47c75e0cfec5a8b8ba93c6c6849631dc5201b6ee1d4cfe95306b07e3d045182ca742b9d6aff40fd26",
      "3f2e59e0b4e487b1f5a000379dade32e9d04a341ebb8ccfe5b39bfb8d12ce44e2beaf53de1c458c160616ea030e365e743f9b4078d628486b3f7e9b3",
      "3d46234aa712dc92bb674c44dda3c10c4846fab47c04fdaf82e4b45cb36d723c275cf1ad424c454fe837373768634a535332c40ada44a2af5a6447be",
      "ce937e7e2a05c1e37249de69115b6cca7734e8dd48abd24a8aca768da109dda49597945db77241d22a8edb58310f9a26bf5c82948c7031c865ed6405",
      "50fa775470a7208e7b076ffaf180e9b15c5304d563170e475728c16372c1a287a76ecf4c6033661d2f4088bac03e3a072b57c4b25bba3b8a114a6774",
      "a8951aa2ef7decd863b180da1b0193e248a772ca268bc354cecb646171b85732c3f5e44a65809f3aa84e3297acee55b2882773189b47c7a5829f83ec",
      "926034470d68347345d220633c8f4fe3e9cc8c7738ae54081c1a7935a8818d88c2a68d3ac2f57119d5713f779f6fe304f7860aaead889abc7bd1732e",
      "9fd1cac142b291a93d29d7e244bc4d3d7dfc63f9f8a4caced90024b94a8730893a9ef7950e6f0b6b9c66a63ab93a4dde2a769b5806d483d89ab9b784",
      "7c8992a599d1b538ccd6414cdb271fa8d5ec9219efd49479f5791861a43bc12b010583ebef85dc573737bbd6498521f064be15e65bc97c2b13adefd6",
      "8eeea9ba7726c3fd69881f7b2d18b2741e054ac503060f5ed692cebecf57aad01a23ad4003f4824a1febb714932271b23ffff83adf707485599c418d",
      "47109c1a020fa6426207c9e3034e62a35ec3e8189fd0dddd2411cf47a772bbeb41cb6bbfcb6416c141c823198e6964a2707bb742cc646827505c4af2",
      "236bd06e5671789b366ca47cbddf1bf7f61b348ea0bf9fd45732f31cb21719d4cb8568222813ebc09cfbd4f57da7ef39d09d7d0c2e95e48e9ea8f323",
      "e91b0b547b1dc3c98a3a164789a272d8fb1ed3ae089118fd78b07e2cc8816094df667bf5dcc6b132fdf5bbe480c9628512db3ac4224f5c84f0bc0321",
      "5c19a06cc6954de43d9136986edd1b3cf2616b2dfc7049b342218e5dcee26559478dcbe0db34def01ae6b633326cb09ec524b38ec1577936e7fb9804",
      "76e7c4b2d3c6d20a91c1542849eee1c1dede75d4555142c6142813fc119589b0d56478eca69ff3c5b8d41f3b91a7124175bbe8c4aa7db36870c9b04b",
      "ba5c37a69e5bb587a82b2d2950070ad27531067af299fe8ef384ea13567dcee5427bba279e36b8c098c091902a25eb0ed062ed7d7c5281cb7e086cdf",
      "5fbb4a433b727d8797b4c577d8f470221c7173ee1d07f7b6179c6935c9c2de9cb7b1d10650f5c725ea1f222b6ab13aeb47fe590fcdc8bda2e3d1b6a7",
      "450dd40c4086ea74505ac37992836f861a9e406e9a30215507489b3eb4247ee61f583780559bed14208b36dc0509c3d39f22833d4729683c6727fc4e",
      "5956760a66593954b9c13b0044aa6ccb115dc4d9a57a9e68858cc555a1662137e2190764765e91a46b8b18a304ce0c9c0a67096ca5de37a189f57585",
      "bcc4d620972a18bac2eb7e2cd2e7011bcbb6e9dce81f09f4715e5e22fb3fee6567b31e4795c24eaa905290f1057b04c34eb9acd16e0ed35ec147ae6a",
      "82fe9ec3e1356b5ee5202d67813855550e92588576a20b314367e175d0df56a404e0c80cb8f7c06e646fef0cb5cff46b78982bad2442e568558f7d37",
      "a053e9da6060751f8b9811f6fb1c1483b488da7440d76f5eee8f9086d845f7f7139d2fb5b405b28079d91f2175f15a81ac4393488a04bb93db750a4b",
      "6dd3833869c3700ecc3a60687ffd98da0f8a2dd286008816673327eab688f61a1df771205c53cb110c3f386913b466816f050dc61b2f09892c3a6003",
      "356dda163257b2b8822e83258e1a243b4258563c08fb9c7e51e04ee5839a54960c7a25cf4237a715a087b3fe85e83ce2df159da7f5ad9daceedea4e3",
      "d84d4ad76447381d8cd68391ea4ce0044945e18b93765c1b84f911bb98f07eaec66a2f4f7bffe7e01ae333268aa545ad975a7b58d9351d042ba4f5c9",
      "27930717f593b2c216e3db8091e26e5028d5e2a86535efc9a22c40f6f559a30dd72e7997ac308a587ed2ebaf7aa7f8f34f13862359780148b42cc72b",
      "ca5ba8f2949fba4736af25e20e81c119c22117c9075d9717b4a6adc277c2d4c8ef801727a4b4b33c42e0fbf5d01cb3e0fd0f5f13d10223c18dd93bda",
      "3343c67ea8b056643a9eb8dfe44b4d76343e906bb36693c53e0889f17eeb3d5364ee875e3f9ea06acc71409c41dee1bf5d8f5c43f2a4d22464d36972",
      "61d0d0c4b95c0a630a6f638da0ec0442b6b12d830a941fd5ec5fc13654d1c3f0f7acbc5ca26049e404bfe8cf892f862aab8be396cd5097184a63d48c",
      "aee272f3db388eef280193c659736c6807d45a90c9cf71f2dcc6a1b01af01057d78832ae9562b42b6e6e44c7eaff59698b76dc6675e576ab54bb1305",
      "fe5fed98be8f654c5f9c829f44fad830663274f69bcbfa6506120b7068755e7d10b37cdc7b73f4ba077cd39f5779dd60dcd1e608c10d19089d082479",
      "75e02fbcc1cb200fc58fb9ffaf657e0e87305eed2dd99802f82fdc2c04e0321ac07f84ae2b3942a082ec557185409b183f7c6bee6c22a3cba61cc219",
      "5ff8e3611d255d57e6b1b9b594adc3a22c4ec4699a6bb3b426dd86e39b1e8c474a707719f6422a28cac1f53d90a9955ea21cacee496dc247d5b5b4e8",
      "0e4c0293e21fb5352a14ed8b4529795cd6bff2dc3244e013d8855206d22d0de4ec2bf5ca68f0d561ea179ff6a05238f96af877e2963fc94f53dd2fd8",
      "adf33c79bc1e4c67c381dbada113a621f8fd3612fa6e2aa7d8cbf6338a33095a6d017749982040b5e1390a0e4f5dc9e1cca6e9424c4fd848eb1ae1f8",
      "10d0ef14435bb1e996e0d4d378f0af120efd5e4260b18f3d2c1a393c67dc8a428bfd5fbc6f9008df7c5eaef6d8b375d52fa6bdbd870f5ffc96c237d9",
      "3b3b2bafd35ebc9f3366166fd5202fbbf0db6de20e8e1568060c097679fa9348fe25b0481aac1229c84f9282d00a7a13f4c05a4a24eec168b3a75272",
      "e1b8f4305149ed879b9bcfb79dfe4c6d9362dfc3694d041b0f06f80b82efbc75429a90fc2171fddde02170bcd749a74f94c21e4154e80dc565fa1c28",
      "7155aeac716ef63569b992d86e24254b65ae0f09b74ed96cc712ccd9722cd9e8d45112f4a9f6e9b09e03548f5bb68c475c9b28e601ddb29d45b7b5d2",
      "3237d361631457abdb3245c92de796d8378e95b7e59665cae01b97535c2a7853dbe043dfc6f3e2e329c8b176c7caa52ff87eb8dcd6da0547e35ab978",
      "5b04a82819bcf8c9c1d0a2016591bb2d3abe89faa4d6b4d9fb3eb61f903934e1afcd22dc644bda51939a580b6f08f7d817d6661977b4a314053c8de4",
      "b74ef2d022b1b1b1d6f93ee7ce0a682ceb25274a342b0c4522d8fec22e34da3fb86795d35e95f459d275fb866ddadb5b08fda99f88f44d40ddd7948b",
      "bcc26be7f1702da97ab335cedbc076a0d34406eff75f56c01d2381633a6327b5e9c5b2ac9b25c829bf986590e25d46ebc079b386d1d174994ed5db35",
      "9e303eb0387bbf74e9846659f12143fb8ad77843a62fad8c66eb0893d03649433b3412037d677eb6e01f546c5e7e2cf8d7d572c7019dd71e0598d9ae",
      "4681d7a82c1704ff422dbeb86ad479840a763a4edfce458d0186dee6f4c9cc55c6acb96660e547dc6cd43b98531c16ff1473bcfbb42f27cee92b163e",
      "3ca1102a0e790dbf6f27b9c58c936cd7e6376aeec0f208d206eee743b9e98c3bf2269a31b18b475634a00e04f2c9264b336d41e0ac5759a9b3625327",
      "b8740b6f90b0ba0d592a5247e6f9030f3c0ed9acae628b7b82d866850f7f4f105f1d8ce79c87e339274693ae3f31dff54407c9416c16d480a4ca054e",
      "619f6f43ae13f963e55c86a9304bc193eb8d27ccc3ebf69e0eac007d5716241ea336f316e9bb6f7cec64781e67f5c583aacae87e12045067d500a31c",
      "a35e0ff66a35079af199d2c74d0234f0c3d8d09f67590df226ec3312a66f27ee1bdd759530c0ac696353304cd62e4b98f004435ba5deb332345625f8",
      "32bc637c527bf08a9779639433adfc0a71844749768bd843b23bc4914ec454d08f34f7069de6095e4ed078c4d8efbdd0444055ed3a8b58939c74f6d9",
      "7cd2c018f9cc7b65110206e030dc063e78ac5167230495b1d54ad95493ceebc7c66b2be30d4ea62ae3f90db733c84b98fb27c9cee9fefdeeb0d2d860",
      "aedd194f4a14e1ebe7d973b455e2237d22234f63b834ddce5bdd505a50a8e99de0bb40c76889595312984b1b34a0e3ebb50fa13490767373ed9c5706",
      "45a95d6d08b0c81898d245536486f331e7f004c358fbcd38474d3be9e68678be4066c59559cd11c5b48d4fd1be74dba2d9b5324bc552bca77da886d0",
      "045ac7a3dd8685df60358c1d7cab6fb703efa1176ff98a5b9f0f1865377279b2a4e398267f1ac81e386e03733bc9a7cac2581e8be644a26a5e91424a",
      "1e48fc779590bd009afe2948273ad5b20f6e79bc6cd112d22819c27ff46fb80d7cede303d019db4abd90515f4c2d70d0356595f3a114c4d6480c44af",
      "33fc56f78574970c4dc8597e8efe66cd45e4368e367c7286a4856422f56dee5011b2ccfe886c80aa4d4d0bff409a000501c118e97f0c8b8ce1a01252",
      "b03b1bcc5ec314c4b746e3307478a67786c66bee4cb73e0acda28559bb58792a03adec97554825828a2483e7dc32ab9a3ac9b40d6ea263e634d7782b",
      "913109be70e65300811587d944206525e3eb740854b742176d49748bc37252c0eb9a524f8ad3e4e2a498f6301efd71d91bf718f2ab771ad530f935da",
      "abe47b7b1712429e5ab995464c62b1780d8cfafb3cba48a40faa8ce929aa59e86e8543e5b89387179f40652f97e874cf312c917f820f6054bdc0cf17",
      "c658c00e0b09eb309b8b2b60ca46f0bfdd342ad3686dc6df312790c2f87ec1765d6ac14a767244166ec29dd44156cd94dbeefd12c6304c0eee51ed2f",
      "9f26bbdc72b9b727ccb74f5a53241236943549412d52aebdf4c400b392ca8c2f7c1f2c9617180b0693b5a1d880eb039e40be79542eca0aaf8e66f843",
      "6688f8ddab2bc474e9a12056ce49b309c3a6beeceda1b9b790c70a358a8068b026904b90c87cea31d29dc40bdaf6a97b7373728a3c8ddbed52771997",
      "da1bc9f7738364a589fd5023519dc7ac8fb225628df02f99f6cf92dcbad54f5d9e6997145abb604e1c2d9c0d3fe742c8952d42524ecab4644202ebb4",
      "644a12874384e79ad6d65c23adb7480622638837e0aee073e8807f3a3176a24311fa204ad028cfa3af0cde94a6b936ec4481860610edee02dd2527ba",
      "583a27492d4080b730e6068a610f221dda2f5dc1d232465f3f8f2d2dc14f688fceca46e8a6ab3ee8bc0bbac042f99caa9edc3d04bb354b27397a191a",
      "3586eea1d582cd8b5cea3e1f5f8054d82c28588cd7757adf26f2bb9d5558ae924814d7ae35b161ad097fade10a02d22fe249a1efe7ede16996ca9d7c",
      "32448677475f3cb58535b5d05cc7dcc6ed922a1a8cdbcb9376e2a9897310fc30c3c4b14e69fbe505bfa68240fdea80542dba5fad5b7adbbc15323d7e",
      "71c4fdcfd1c336a7ddcbbb2f0c2f679f845a989d90358b6451c8f101e5d8041779c4c22d7c49906e4982cdc7b923abb3148bca7fcee909fc41a9be88",
      "e55fb3c8b5b05c74888db1352817cf6e8b31456a504403d548e2a91add58b4f40f306435248eb72e1512eec07d65842cd0d251ac19956d19f18bb0ac",
      "89cf06ac669da2dcf404dfb0e60407fbab48ed9555ad5ab063457699d718759102bb44bdc5874b8c6f0979c6f5ec22bfcceafea5401bb8f25dd39fa1",
      "4b3dbc8f938055739712f810ab447791fdbfed10fe51ca8d7cff3718d42c6b66a8b4f8aca2b3dfa562aef3fa7251409dc11e0f548fdb2d03f3f3a949",
      "3f6f0dfad446e981aaac58a6e1eb9335a7d2460567580b18d4c269ad93e9eeaf22ead517d9129d518e1031d9ee36bdd8884e954b8326b400ee8d7646",
      "3dd20641db5439fe40871ef79e9933b83767c42d9923d6550876e8ddda239af892a11ba33e69c2870e79ee3e3529c7043a6c03c46cd3f631b7c00542",
      "b575abe3d6a8c47443cc1592868a386fa250723878e8da5d29250f09e7c16e217f272531ce328d6b90af30700b4122aa567865d04ffcc8a53f8fb964",
      "768114f20fa5f22666175299e9cd35ee66a179e6c2a8af8066e26f6694e791f5269e5a0fd23065a70744cb4aadf2e5026852743fdab997ecdc8bef9f",
      "2794c7ae94407ea1e333e0bb0bc577eb9214c71736ebbb45d6d02ac58648692f33233f1ed32155cacf6aad2435551cf68167740aa783e089c4479a46",
      "fcd19e44740bd4b365093b6aea843575969414d3f1e2f646506fa83deacb46d9558ba51ac34ba46149b408a6ff0e32991db405454b947b94e83bed1d",
      "0d7786207f3f7806f483dcd712fcbba200d6829e0af80da2e757df7c77f0f7e1b0277132bc6d4626704c2d061fd27c1d2b82b8b6f67dd89a13fedccf",
      "7a9bfc29d2b30110bde9af623c6b921f9d0bdeedc9ff9fb86b3e9c6f0c4745939f577dbe4684e928cac2f0a05559363fc2942fa9be506d6dab1d5d66",
      "d3a5e35fb94c749375936e328c1fd564eb2ff9fe7d577ca849766ef9e2809d81878aa9faa8a306f1e30018aa6c51b354eabc49f356a46f52442275b9",
      "6e47016775648fd47de86228d435a29f11aff71f9640384fd4055fe36a0b4f0dedfbc2c39e8ee23ff5513c69da27f9dede11dd74b65e749e970eaf10",
      "12a3befc4cf7d69692d9c3a0de3a5b077f113ed6d5fbeec3fdd8aea3eb3cd759369feaed5176dc986dc30760a378b3adb56b39ebb8d449ee5b75f329",
      "8ed1b6f060c1632e94a12f8f1bf22656f11f3bb2cab5dfc8b3bf331f1ea1330e122f8fdeb55971a9f72ca518ce8ccd2b040b62860c33f4ecdcc80548",
      "2cbc35255403ecbeb0ae05c6a7ec39f937b521cb906f59f89368b22f79073f1f710ef2ad640f3abf22e6b1bfa1cd33843357578b5ce56c6953c7e1be",
      "c81eafbd0ff0f982b9dc2a78056083e36ce44e1457b256515ea4bd4cf1145ba9cb05c6af1cee8c92e14e8ff0ab5a5ccb11d404f4b78f3925df52ec1f",
      "dcc336a00cfdeecb871613946f6664dc799df7b8cc06c64611acf53c0de7b2b5003fbb5a0038cf6e6e98d3414dc08f461340897b7b3f0a1d53ee53e3",
      "04290792f70a486c5eb5a2e3441631821dc17914e196a6c6aaa875ed8bb120944eb6abb47a3627b24ad1a6f328eeebda633746aa0917a71a820341a5",
      "24d261c9ee5e40b9ccf5f518776445612b762dcf49fe6d009ab4227969554b40ac7c2114ea35c7046e4c446faf28aac41816d129b73148c192bc071f",
      "91e6676f8f14423737b42294e3f2ef22bd368c6ea71bcab5cd84995cb63780c581a9c834645d20dd513740b4e888350a99a55d13ad7e373550dd092c",
      "460a89435c6079d55cbc0166980e59b4112eab40c20f9949a0229c708a93f4a3901a25e0a1402c181061821ee0fe4ee74df211ef85cbcbbca956f2ae",
      "203d13983a93be1edcd5f45920b14e8f29d5c481c5b7e9919558ebd44f82db6a9e406faeaa3c7d440959ad9eff10b05b2f73604d5ef0e78f3e9f95ec",
      "75be3b949ffb2d0aec6bbfb6a58d6a070675740cd65a65aea4c8eb81b54b6c50c7d70d790b93223b1d123a09f0b08e07b39459940701ecd39fcaf2bd",
      "63b1f8baadcad7769f9eb9e863485255f119ddaaedfd97d58675cdbb03a0e05b238185042f227b2ecf8849ca768fa62a645d129846d324a9fdd9cae3",
      "317b0ccf5f056b332871ad1a813e5105eac2a1235bdcfd61ad9b3a9c32df166b1798d40bf1aef9395fdddc60402660202fc4ec626f8f1f501b531606",
      "65c1b274690f2a6af703c1e9f09f67e62521b1d03dd29651bad02e58e99293cf9611f1f8559e68e3e471db3659da1dbb6f7e5926eff315ae289d0e83",
      "c7b9b9d511f5baae36adafab0f3ff9ad879adc5eb735cee4f1738642654baa5f7317461eec179930d75c38532c6c42655829a3dac4c666f40eb32add",
      "0d25b48729dffa4c8c72aa25f8c91ffa3cb85642b44ee80049f25e7b604974235688fa22a081c7d92590285d07c67ade66cee224c345778cffbcffcf",
      "2ca05fdca8769dfb4b2f269e06e32824477fd75da7767f00fbc664925792430b956c8d030887b56bad8c2d8bd0c28add65dc687347774d53dc42754e",
      "4a1ea5919ed57aaa9d793df8964c1be0ab3c2ae7e6fad29bccbb102d722979d460746952578806cfaaa887472047a385dd2dbbd84bde1257ee9ef3f3",
      "f2ea6c91939d883eef1373d087ce520ba0018e62e54761cebdc7a6ca0de7a855133af76d55959bba111ad640bc16c6a888bce41ee78b2663a35d3bba",
      "a23e67a6a25f8cc57e31b182d19b1abc796305800cf331d60f0f941ce9ae3af6448f933f367c39af911e7ceac01794b62a18980f5b00c6ccfaf9389f",
      "e5e85a44c1e2bffa7aa85f1c65d7af9739f02318ce28cfea1ccd9ab20a63060d0e928f88a1fca85caefe797549d1f74032e635feaba86718d27035de",
      "bc52381099be8466f31d2422d647f2ded216ff8d245278cbcc60c0b154e7c735485340bb5d53bac790f21a53127ce724538aedb08989e8b09978ae4a",
      "a601b44dc204b3305f563916a3a07bdad5471b25db9367dfe86f5123c2299a0a21a326359e7fa9cef4f1222f08ed200ae688cdd94d2af3c061745e65",
      "673a084c535ecd2efad4b7f65ce405fe892cdf0f8a6d59a351bbf255bfccf1b173cf7fdaa2333d737676655bdb4de96a156de4b97511a5ce3775287e",
      "9c2dc0beb4068ad067afc76e436ac475cef1d42a978e07ae05b5a24d74bb5feaed396569b63ba78e1642ba465ac9e9dbc80e17ab444b2f18ad2b0c3b",
      "6ebbfd91ddf3753558c2b25f03198a045c3ff5b9fc418dd18e186661736d333d711c8fdbe367c01e79cee302908cab4cd174866e7ccdf5061154dbeb",
      "a96b57c47a3184ec8626ed8e9bba73edf55f13b3e4846133af5e9675da5a88d9590de99647927a0155c00b7f37b39f89056d099d57be2736fea20a6d",
      "4b80ba2dd65f7a62da73d7639b30b770cb3c75323af43bc436db7b7bdf0c0f9bbed5ab93e1e924f633d3011eb9bd97695dbed15e6c822082bb1f4c99",
      "d7b665e379324c86a736cca8c8381601ebfc199b4dba9fc706231173b019f0fa5006d6e5e4e01a0642bb00a3a9c3c2ed806beb474d3f8fef79bb8626",
      "c13ec0023469e2955134c7da2c1612576b45f98ef65606c995d7663aad3ad117e450fa87309a746a75440e96102765268aa4b99245560c2c1d2eb252",
      "fe2654e4fa7bc60bdbdc0aa157532bb2c8e65618c83ae9989bb4734a98f04901c2e207ba56aa593a2abc92643739bc52d7424c8acecf768d01988a7e",
      "be2fe2ff6ef6f33e9a3246c53e3ff7f1d902e32a21d580ca17f0342dc7c5a483fa74b6395b94451e754e177e1c4d359d6debba53e83093c85456a9d6",
      "642759ba5b293a8a976045f9543f58144261ff4e6a8b6031e31527bb167577299fa6c2c6635892617f0f680245f83c99715f657aee91b5f8fe28ff0a",
      "3e6b357078732a07228b7562e6f733c45f32b0f54c87f40d189c8e7afc8ae94f32eaae5d8d806ae06b8bfb26b72d3bdc830146a55ba45fc042d3f808",
      "88f85ec2c29f42361c17ee085e7018552ac47f5cdb66c974a9039596c9e7f99594570a4f30c993a2448507735220e448741a925a41a01ad7d6d983df",
      "8ec5acd79bc42c2941f540206bcb5720ae9f955935c7f6d46f904b487c1cffc31582a4685f8f3f2ff277cdf8e4ab6f4f930a7548f0f3bbd384943127",
      "5f0f4f93ab25a4fce394c0639e42dae81fdf26a3e1b7c9df47c9e8db11a7a31b00bcf80e5e7c977cf335bcf8eab41331a8ad7aa79b44eea67573e68f",
      "c08b74e467d78272875edd3f386d73ca41939d5198c8c3744834533da7bd618f4114e0bc9f5dd55a127b824fc9d9e2aad2490fe121d1caf5714bb7ee",
      "4a6d5229e154abba6dd76bc926b75be928ac9791b17af28493ff91d1bc7bf765c3096bd3db76545437e9892decd5b312f1cecdb3fa216da3094cac9e",
      "55e562715ca29d8d79207c317ec2d82de61d3e25a1691b500f60563f8a026645574568963c109ac9f647915d374a6dc37cd2ab1e97c6495f5378d456",
      "5a4301635528777ad089cadcd7fe9cafeaf4f3ad3a76152148d8f1105dbddf35beeb617ec2e927c5e9697ad21fdd6b80704f8cbba21810ed0815a317",
      "5cae55089756a80ca67311982fb7c8d53254e06aa9b263df4e4d38b2cd33827463d39430ed0b35cc6fa8f1c6e9a7248b5d65e01515ea4b0bb1c035a0",
      "d6aff7a472dc2856a4ca566fa0cecbbf3bd41156045380a206e178e93eab19cc9a6a217f5ee64d063f63d79b444e9e13fb1763c91a8415b906e3927e",
      "22781886bd54bf910f7c5995e7555ed7b4d59d42e4afa4159fe80f9aa21ed33eebbf222be4802e82bcb29cd39e221d07f23a35383ab229b9119ac021",
      "5e08612798b00e5d25c3a250ad300e2018402ec6a56423e0058530da716749ad4b9826a5afb5620aef4796faf7abf148894eaa144641e0c04850ebc1",
      "b708c07b73b38bb35c2ed1de50ffb4ec21cb1a817feb59b6ccf1b34a2c87523a126366104b5fc9e8b530981f36302bea00ef015bd366569878beb9c8",
      "f38554e7271ff1e1718ef60f33ba71197f2d135f5ebfcc2bccc6218effbe59a1fbba432b47d5995d40edbfa6ae78b87c4037d2e1c67bf5bb564f1bcd",
      "4a0ddc0afce4a066da2a96d6fc905db0e9ef7e9a278655a733d588339bef84adbb28acdc9807ab94d41221b4124bf48d64d1f70ae05d583c09dcb7a3",
      "ec50d97a8f2a74b3b6cc2bc450d061b7d0ce26792523115b2f37c56eb483e7b56715b7cceb02bd6a872843369e629366f67ef4f5bd8eefe12ecd1825",
      "93afd9d57c519e4707f7a3d13e3de5d7cbc8446b2f801ac4f7e9f1d9f3036445a6269861715fc4f704c814e49d47048be3e5ed47c343ab575331ee03",
      "f5e82b1d885cb4fe8794ee784000adfb21c3387b1b71da687a93ef171800c3bb216a2db004f3b1fcf27e0349f0daea010d695386265bfe542ee6e3b3",
      "01a9ce39102d22f5d6d93b90aef0d765792616b0d8b3e5b186a31b7df3c54481da91f16ccc85ded6c63577ad4e3e06ccf0025be502a6421ba58cd8e3",
      "1cebcf8dae912018642285c1b22f95d5fe972a1d8e8b92232a07ac85df12f5eb77e3f8abdb4bd73be6ce70da177aa388fa385b520c378285b01918ee",
      "1b6a285e38a627fac50efde00ab549073dd195e2cb8a7c4834891270dca3ffe307c6d6b341f088116762f1f7e170736b21083d7ea311f4bc0a7cdc3e",
      "05d947a68bad3eeb8fbcd3781c993329ef374f8aa9df5e839b56ccfbbf3afef7db8e5f4757f563b7594b7312588b8d5e8b8d5e8b140f502da7ce1b2d",
      "a976b3144aeeca8f0f03e4dea6192827d86bc9b765a61e1ee94313d5006757676788910eec0a4f136a00f109bd5cc62c182640e6f2bcf80d6fd7f997",
      "4cf95da7fc9ec8418fcd569483d729bfdfa2238bee21498f39f41ffd3270d04bf418961d63f2fb6203f99594b507743a4076f12bd1712a096fc28cc6",
      "03a2fea436aeb731ba38daea3e144f5d56ad6c4d428f318dafcad7d3e54f77e9c1fbbfd40327b1cbc35ecb41bb7a2fdcdc1cecb649393a46faded186",
      "afe7389ee8b5b98e8fb3761bc764dc2af3e8d4458d656bc2c2bac8796bdf2f59fb4c9afe16ebcd95ac5ae344b1212895987da9969fd15f8b424eeca6",
      "36722961399d8c9d6e3e1608ba60d30ae1bca7f04685c52aba76a789f3197027169a4fdc7042d0fa70905db73bff64496340f142dfb8e82a7d595d7f",
      "d0bdbd8320abe00c236c87b4f5e22ab6cd8c8b8bd43d66f28bf8bef771261a618b4089199a0d119e58f4f5bf559cece64ed8549750e94584951d5265",
      "471cb00a9993f270646faceff93d0788c1f7a3c332766ac0c8dd860fb780aecd62766caef462b680c31c38666beb2053574a69fb964576982b540a61",
      "d47b7b454331f7a2c6de55387f58d3de1e3a6808bc2fe41f5a9277ba610762cd3d7c727fe3d9c8dee6ae48d820d20a624c59f24aab14aa2482ddb23b",
      "f27694156c11099ec2059e1469a60a2779d2812f9e7ca62ee2b14b4209851761e71475ab6d559a5e47cc7a27d8ea0b3a2c136b7b75dac83a9b651a79",
      "2460141cc0267ec9a706efbf941a84055906eec7d81f2ddfe3dc5abbce8c529ee0ea4e768749b7c01c5ba174d9d6936372dfa5056ddb62a6953899df",
      "c4064e466d30e09e899779c8a6a3e9db061e1ebbd688e9234aacb2f3f40ffcf5b3544ca5aff1e9b52b93a63f60e231734ae97b7a508da52ff0de5f82",
      "4ebec39f55f9b1466f18f8fde8d5ebb4c4f7ec8a943e41876876ff4b3f98df4759f51e4a9c530a41b1e454088827f4237d08bf1f8b0f4fe1a44b5f09",
      "7e78c1fe26e9037a2e2f89dea6970d3ffdd28845baa4076d9af81bda225236d86b64ac201d68d8e37e8615fd9897ff7cfde279fa911e2a21b19e9ee3",
      "a39c6844e1195057e698fa148301a9cefc8029f260fea07fcab9a92901ab2bd52f6ef98a4ac12141bdfe803d3bca96e903fa4b137389398ec4b52852",
      "bc8619900fd23b4c649e6449bf162baae0bc914fe7d0c76b211fd0b5365df1d3121a7dc83f6bbd39d227383f2f99e8e1f347ac07e138b3e21ca6ef77",
      "418fd207ef0d3ebdfef5c7f447f9434dc10b7cc698313cd4c7f45896a8894e9fd2838201e556336c85fadde0af5f81ae97aaaef798ff0cc48734c31f",
      "3c8d6f5124e2d5f40bbe9ff16af221aeaca02e1f844cbcd5d9e5bc76c6813997716012b6fdc545f80b4ecf5caec280371b167b451e76908f65817481",
      "4eb6f622eb6c70661a3c976bf0b5900f353478849373ae176447fbcf55fb624d27e73a93bd74d7ccdf0ff6fcb1e554fa51ff94abfb9c12f49aee1ce8",
      "9519285e7777e67b60655bffad2e4d8dec997249bf78fd7796bab64bc9edd1997965675eae1bde4393b3d25ba733f7ef26379c1452ea4c670404995d",
      "2ed30bfcf5e13cfd59f08f2d3ed60bfb63357a1ffe0b6b68e43eececcf53d39fc6dda6846e29e36ba5bdd1e0bbc1d7e841fec1deae6b7af5def44ade",
      "133f4339be40fef81d1c494ddad667a88ce96f42c61ed5291800fb7ffd5f770e99c1b4450200"
    ), collapse = ""),
    leaflet_css = paste0(c(
      "1f8b0800000000000203bd5b6d6fe3b811febcfb2b78595c2f39588ee5d87971706daff7826eb15b1cdaeba1df0acaa22d36b2a88a52ec64e1ffde19",
      "92924889b2bcd8a21b246b91a3e170382fcf90f4f5b7a460ffa978c16222cb979449f2edf5dbb7d394d14dcaca20a7199bb48f254fedc71d2d9e5811",
      "f0b5c8faad32a1b1d8775e0e80b4a43c63c5c41d84fc9ec8e76dbf714db3672aadf657217641240e5613dfd12d0b52fae27055cfe4d3db37b990bce4",
      "225b111a499156257b7cfb26659b724566f0a914b9fe706cdf6da4c4f7c5332b36a9d8af48c2e398652ee9e7ab0479067b163df132a82436b394ad41",
      "9a4c64281a2124d889d7a13ef5cfdfe7708d0bba6dfa8e6fafbf25bf14ec9965a524ef7f0a43b229c40e66b44d52f82d79b6253815497846a2b46268",
      "07ce24572b3d1e68126710d1f5d3b6105516af4859d04ce6b400e68f6fd5507fa71b5a70b0ad2c668544298282c120540d428085798c5859829af7bc",
      "4c48997039215155921f12108e112ec95e14d295456ad68e6c2890b6023d22cc66d56843e425dff1576d7b059565ad8f04a60083d292e4b566923d51",
      "7623c9852c41c47502ac2ec83e611949058d514d19db1b558d89e5dad19ec765b222e1ed6c961f408684a1e2ad865a60a5ce8d28768128f89683e1ce",
      "ba067aa699c55ce6309f15aca8583fd5f3f69839f8de0a9644b2123f921d3d044a5c12b375ca0aaa565d263ccf21508085fc05bc30a55f91cb7fab0f",
      "53516cafc8cdf4802af18fc07730c28ee6b8a851219e40a1c0e7e79f09df901751c1623cb3765cd0ca6cf635da495fd32dcfa6095d14e6a983064e00",
      "e66ef14227205ff15d2e8a92a28daace5afffddee3c9c18c96d55830abc94962bd1667122b9b39490aed23616f80ed176ac4bc46ab52c0534ee358b9",
      "98b6ca3141ebf1312e300661b42c73b9babe8eaaad9caed1d579b54313baceafebc76b2e65c5e475cc8061fa071e7f075e12ce6768096f76fc104429",
      "783ac4c91818e669250315c6588101c823502b8ca8d689ca222a0eefa469a16b9d2340fdc141fd7d51c961b8ef383e0cc6e081617866c478548af999",
      "a6298654027e8f5eb14ff81abc4f30997d033e59e5b816d65b4a0f2e4b13ea87c73953e0b355e51d50370ee5d226c4d13c68720f90a4a2e8a5110f07",
      "3acaa3d846f472194e48780f7fe6f38709994d17570ebfda1e373c058359411c4a2065a0993f73c9239ef2f26538d9079807200a7e72c95b26c73e58",
      "b1c2ffcc8afcf859a57a2009247f553e158902321836d5f97eb0f735e090ed0e2b723f9bd5c11d12e6135579599992e56baf3c85480dc3a9ffd1dd64",
      "22f6ff829ee97acbd1c5eeefef6fc20727d8fa22eb203a397690638d56c827d288ba0051894dd8063d9770ae09fd92f4197a23ae4db7ecd0d961dca6",
      "bbedd09542a425cf6b116dcaa54b998bbcf2d1dd75a70c69509369746bd3863d29f356f75e05b5b4cfbb14679f331b6db850433d217ba055288e41de",
      "e5e83855915ebe8bd9865669f9eeb78f1fc0655af8c0b314fc2fa8518417518324687f0a648994d43488989c92a226709079c152c018cfac6fd4b9e0",
      "197869a001da8a289f4bd92f149be3470243beff893c04e1ac89960a4a60aa5241b2cbc0e4306785732b8146a22c75c4f34ed35e299f7c8d2bd8ec91",
      "9bafce287059b0b36803c2d1278afed423c02a06bbdb6ae6e85534942f1408900c6800d1d102434999f8c439c14011f467e77b03dc0ba06ba0661dce",
      "1abb736776e2c57ac6fd77d5a44fbca995d17f6f7876e645b308cd9bca9e55faa3184d21eac307405d68139d3a193b03ec246e3450f5634ed72a45cc",
      "ba18df58972180443597041d8d1656197816a5fa374e793c2971136b4ecc21f42438ad149d12c76a18a2a6254ff63773f1d64010088787df436e0bd6",
      "09cdb66cd5f2703161f31239731a46a5bda9a1729792acab88af8388bd72565cce26f003cd93f06a600d9db6b3d9f416b8a7a87136c7933ae86c6420",
      "caed130ce9a6b375e1efeacdc1031b861607d018eb422e0ba1e9d4531552141dcf54d119c1e9b37a5fd3002ed661dbd5cab6a0914d54cf14db1fdbe6",
      "7a96fd66f5cf34db91b81012123377eab4ba8d8c89da50ba4c5bb431f1c6fcfaf57ebe4390bfc5bcec4c7c3246e015f4e45b0662a90e1aa5ceac7642",
      "657b9faa236032a0ee7e57a372d37554a6a04726bf230635e2b6562d3386451be27ede2ea6bfea6eb730a13c2c93c15d4f0d2b06108357a411a55be2",
      "8c507644f453dbc17584b5e2e29bcfff0ea1e14a02b78aa6a4dc33fa24fb48b2d1aabd1ff92e8e636020aa52c156b1d948e6a05f7f61694ac877b3d9",
      "ddfdf7f7c3959caec10080e70798518929e3ddcdfde6d1154255a2f3e57252ffcea6cbab36586d198c8b137bc905986f9ebc0c6c7221fa820748d83b",
      "9e42d4bbf8334b9f59c9d794fc9555ec6242be2f384d27a4699f100931168ab382a350ea6da81f21238673a502ab6536bd5b160c372194aa9a4261ba",
      "f4880a9550440bef4101b62bd51c8ccb60c206859325fc2a55a85c044ab8555a30756c41635ec1922f7a00919a4569555a97f8ef369b4dcba1c18930",
      "0e80740e6bb15eafdb2dabf9ad5301994767b2a6adb7595ab24319d0946f2155ad994917aa31666ba137451bf735c24529351badee54fa315a7b16ee",
      "a96cb73a345a536d4b8fe5ec6bfc758c2b2858cea8aabbcdc747ef5eaf2bc22ac16038e9b66ec4ba92439a5ee08f8fd7861712a692f0346e5d02d1be",
      "82df9d75b57a15c61e59f6554a3dbcf5429f606f087c23746ca55fa3a971db4cc625a6abd8ce57a6347e3cada7ba218aa24e4857bb64a46fdfc64c6f",
      "dc3301f3e898e94dafa2f1f13c7765e6a757667ed658672fd47c74a1ea019b9aab8632becd031d8e793619e881c85fc74dac7341bcf01ee2c3371faa",
      "358f29f901ea3791b26f26e4a3c8e85a4c008d6442428dc36a0fc7223fb372c6a022ba1291b30851c026b2eb483c9fdb65a739853aa90443735ed85d",
      "5c75f2931b46eb6558fa72e458b45248a171abdaa2ddc06b1eed7a5c9f007ec1485687d1210c5287f3e3e892f50632922f168ee48bc5699d0c324eb9",
      "2c87e33e3b0022c38de471b99ad0de0f5c9fcb1465f29e0c7af7e2ce18e893751c84dac7dd13f541af421d106f6e6efcf6373c8484c2274d0db0a80f",
      "e1031059773c5a8d07ab1834d2d4bb396306adb7b045d1ddb3d2dee8d58aded31a318a94462cf56adac663375d3c164e67f7a02b85c84e099d53c01f",
      "5a6afb34a18de70e22d278584f4fa984046a9dd4a75b2bf2fea8b31cc1ca83fce36f1fdc93579303555d1268f88fe8be926007224b5ff044179b836d",
      "c5a4c4523061550116c7d7804919231fa6efb1a4a94741f03fe8dc2610d2b22c7854a90368dc86936b0a4e61f4315810f41467b3f9e40f843eec4eda",
      "3fb3e9fd95a544ff66ab3d8ec7f595f0016675c76f66c64a5d67e900f2c5e8801a500c01d4b1777bf0d04fd5c0c5de3815debc40a1ddc1acd7834daa",
      "cf433b070bee59334881554c5a43ef88828ff28113e950794ae302d3dbdbdb87aef378378cf55af4b68d9767ee54f75eafb1e560bc71d7deae218d93",
      "dedddd753cd8ac5cc71042fb047e6ebc5887a37dc24b162820836fefa1aefcc223c633bd4219438d40422dcf8908dfea629589f2d2c6ac572e681dd6",
      "5007ccdbb13b988f2f427f60ddd242daab1ea6f508730632f4878493986472127777009fefc4e9cb78774dd30192f3ab4e0d94f21caf46e8acabadc6",
      "c46e7d7ce05e28ac4f14bc876bde8abbe360f35e19a477427196f04a80469f9b9db6da4bb477d8dccd515807fd86f393bc5b77d7a99bcc1768eaea93",
      "a7600ba7379f93ee7730c9fe41ad578edc91e40ec69f59a909c665bb6e82d20cf020db77296de1969f6636276e4eea2d09c7e9b4cc4e340d0c23cf15",
      "caf173d34662fb38fbce3dcfd68fee3ab78ad02047ed33fa8e6af53ea3e7cc0a429d2869c92e17cb986d7b07398324eec9969fc439b3e9938c99f564",
      "403f769c5659a08fbd9d1a116d365c788b44ff066957ac54481640542b35a0f2bb737de45d74116a7b79c8e7f1f5d69d5b8799475dd487805aaf95fb",
      "fd4a13b1a313f21b2b629a75373f6b1ddc2df1e7c4fedde07dd6cfd78717499df15a03ad6aa197f7f8e3b30a5d02e9f2ab75ae76e3bcbd3593c6bc7b",
      "a0eb0b9568b9b833d13be1f532709dd22c4c838c5bc702a6f525ab8bbc105b1eaf7efce77b84f9bfd62e30fdc8f1804b6ccae9470a89f270f9310cbf",
      "9b4def6677e1ecf6ee7e423e8673f7791e7e17b80d36c1d5c5637bb5ebff36e869b5db7b3f9351aa5e9e3e671127e72d59ed7e5671f8f0f0d0a6ec98",
      "3febe2cfc9dad0aa8abea1d2a9cf147078cbf4577d85ca5c0efe13e0f9fa0c01afa89509235084ef30b1a9dbd8099550ca94ed4bddbb588301c7de89",
      "781c3d3b7004f6ee86dd3865d97c3e1fc6d8272ef40fde9743fff0b5fbda06f2657fe7ef662ca81b259e7f427ece0d2a7d470e037ec460519d2f8d98",
      "4e03e006fb156418ecd52944770fafff9092f472df36cbed46f8e1d86f9c0cc2d745bb47c20bfd3588ceee833bcfee7e92b248bfcabaa4c110ed9806",
      "5bf5b797c21ca06650996783d4c7c0be74d601e141d83b3f703d6c4c7ae7469c0b23e71e1c79eb393f181db0be1de761751c30b01ef910f5e719ab9a",
      "9907320f2e84c5be733dd0bdabd6590775e632aa97ae74eded4157532eefc2be66dd30d7df2b02afabef79fe71c7624e21f1725d2fb5df3bd26db837",
      "aabe7a54b09d78c6b7ba3b8092884db3c137559b84befb2d4dc8556cb560018dff5d49909d1d20a081806f4e761e7106ff058d5b8fb803370000"
    ), collapse = ""),
    countries = paste0(c(
      "1f8b0800000000000203acbddbb22649729df72a6373bd77599c0fbc832093441325938994e982c68b16d0045a047a60838649208def2eff96bbe7ae",
      "aafd8fc41a716ea6b757fe999191111e7e58befc3ffdfeb77ffe879f7fff2f7effdffdfcd36ffff4c79ffff20f7ff7773fffd56fbffce1d7dfbffdfe",
      "dfbbec1f7fff2ffeed7ffaee3afbd77ff8e31ffee1e73ffef60bfffe9f7effeb4f7faf7ffee5fffce5f7fff9edf77ff3f31ffefee7dffef8cffc53fc",
      "f27ffaa7bffbed97ffe50f7ff7cf7fa39bffd51ffef0c7bffee5d79f7ed3fded7ff594b7f7babe94b5ebbf7b7bfe9c7336fedcf74b5fa34a764a1d92",
      "9d2fbb4d64fb4ba92d44f3aea3cb56bff9d372d7926cf411a251fb94a8effbf503f5fc7f67827fabbbd536971e30cb38f184bebbcb7abf3707528f0f",
      "64b5337324b36d93d94d665992ed2fb787ac9d13a27eeaf4cbd6e821b39be86ebb8d1da2b5cbf5a79e5ae309df8ece47fdcedbeddbbbbf4e2b0c51c2",
      "5bf7f0492df548f8ed24bf7ff70d5edcc9fef79ffff3db7fe972f8373ffdfa1f7ffaf5979ffec492f8d3aba1f72fb7749ba8f2e5329b7d7c29bbf1aa",
      "5fcabc8cbddb84dc635774fbb8db057bad8b606d0dbfdf2fadb4f6f63e4cb2f423fb4ea3d8cc4ebbfdb921b9d7d68b2d8d3db7ff6a0cfe3e23feb287",
      "ec2f71c37afdcf5dfa8d474c7bc2f952ca3e71c95913c93873ba64f6693fb27fa92dee79c7e5658a8d5d8f19e54bafcb45bd6af42f44dcaa5549ce5d",
      "31fcd1361755964cf379388f689d98abb16b976886687dd99b9123ba63ba6856edb06acb4eaba64f7b741b128dae6d625f62facdfd69c39eab115556",
      "38df6d94c3cbaedbb74b76bf7afd116fd1be6c7b2492d6b52f4c526fbd4cdaed455364b79fab23d9abb9a0ce7d10d86ff5a3f285d56cbb71bcfaabf1",
      "950aafda6c2c8d0faa01db5fa369010cff6beaaf11ffd607cb65fa1f7ba222c697315bdcb55606658bf196ea9259e664cdd9bce6a8266bae7fe9f3f6",
      "1c59d73aed2304a6176c469a29b21d3f1aeb4a625ae97e929cbcf1e948dad9cd25a72ebe173b217f651389a476e90d7eb5eed5bea93e6df5cb59cc87",
      "49da1ef1cd3ef6da0fedeffffdcbeffef54f7ffbd31f7f7c83db975acb1652b36d3c17c3704945328f363d92637bd124fd6a6f4ab2fb5b9b5f8e2b41",
      "16a26d8581c814d54851dfedadd937d87ae9f76a73bd8744ed0cd77dddbee6396fcd96e3ae2daebacde6afd52fa6a5fdf6e8015b2288ec7fae144d3f",
      "da3a69365ff79e14b532b86ab4f65c55ba4b8ac65e877d41dd7c96b242b2baed96c63709356dfba9d5c6a87a2d2d067a6ead8c7ddd1ae3340dfed6ec",
      "86bba4a4e84e8bffbf310b76263544b58c9c2b3b1217a27362ae6c314e530d26bab79c10d5c38cae6754b615faecdf4a06bf628b9df856db941a921a",
      "ef7b6c79da8efdfe0b7f7cf31f5a6a7ff9d3af3ffdf54f7fa66d61f3630b60bc0d9f19fbe3b21a861906a54d970d3b1e8b64a69263d14cd31c63bd4d",
      "f6d50a99e9d23e253bbd3cb27b4d0d4ea67c9e949d690b60da876c37aeb323a0d96133d18073a4ac175b60131552fd3a3484990413ddba3f44dd36bb",
      "c94e69e74fca1ac37f9bf38b5b21ba88a309d1cdb1990ed89271f2375f1a6865d3bd6f9373b5f93c759418ef6aa74c2dbeac51f8db56a3c9ce8a2d68",
      "323b42dea629f41d536227879da707d93e3b1ecb196b4acd6463f578c446bb4e6e774b8c18256ea68e5d36578ecebe54bf6fabf81eb177b2d5b74ccb",
      "fbb92ed1bdf64a26b3af323e6476ef659be2119cf5b66c4cb1e13f4b38c09b99088b83bfe618e758e36dd9186fbe897ddb85e48406e8ac157b61fb5d",
      "aff39997615a17d9725dab79b6e96dc8ece8e8f93d46af03d91d63e687dcd81adb76a39fa4be7e8a9d1c0c775fff6db3d9b275c3756602df94d92c36",
      "64e3f4673d0e7b21647def96eb7b73caf2ae27d7b28db98da3e7ce3352d6ced573cde08db1d87966df5fefb1727ccdd473d7dc99ed13ef61efb69b9e",
      "61dfdec75279cf76753f7b1197d9bab3030959a9d5dfcd946f6b3ef16754bf5fb5a56d47a33e8699822e43436ec9ba9f84bace96b3c96c6f96b6f2ba",
      "71edf830995941a112ed84b41d63229b011715866747243f8dab8a8dc0f61822b3be43739a1578cd38e2a727be2ed7b1326d20353673c1aeb7d3c864",
      "6bc6aa2fa66d27929d7f56fd19fad7863e4c9d98c054827faa32ecfcf1574cdbbd748e8a2a51ce6cb1d99ec3df702cd7d4f636a55e7f433f45afede7",
      "d1f536fbf8302fdecad18a1e6edabe5fbbc7f2619ae5e557d940ebd5285a970125519bba7dbb7ea65cdbfb1c201a989f05d74ebaa35bd9d17942d27a",
      "d37e2931507d31ad28bb76c7cd65b8b2b863a399a8731aee8a2d7ee287cd8e2144b7b912b9ac08047bf901793100d8127604d6e337b783b8dee53bcc",
      "579c5d35cad2c2f44b4cc3b269d820773fa2e95fcb149c4e35537bb6e3b9c814801f7d7c25368c7df350068725b1a73ec3f6793ab698a62faa7ae387",
      "0bf555b5f49aaff8336dfa97becc8e830951930a39c715e419a6198a66d33eb2dfdd36a3ed1e4da7198e12e12414e9992aab5a92567c71ae30b130a8",
      "86be9e29ab99579931cba8e6ddf9c38e9fc3d8abdb0dc74e8b81b6b38f5663ec03c361206a73f718fbc65e5b9c3d312c399cfa61315d125353a411ed",
      "d3d694f4c69ac5e8992ba6d4d41b92e27af9e0e23274939801179f6bf3c6663b853985c4be8ca9e115822acb119129af1acbc7bc0d24c5fd6f2dd7c1",
      "67c7f53acf9a6e57c750cd7bb3cc870e4c5b9b276ed5aaceb71ddfca16d95e5cb4312d77881a0723e6941f3a97e76c89cc60ddf176b63327a7f699f5",
      "5967fdea20b765ba9e7932db86f3febab9c8fcdaf186c834f189af605a5656410fedc8e75b65c976186e28d9126adbf4b489ea886f653ba49b6e946d",
      "d2525439c1ed4577aea0c1623433c766517752f0813b55b6ea0851354d6592d95de16c5b84682afb5d1efc5ca4b7b1f53d6afeee984185682d7fa0fd",
      "b099a5a471f6a595b76d2aefd4db9ceeca79db8ab37bc8eae9f2f7102d1615b31cd3403ca6d879651fac846d4424a79a56c0e429e3b9aa2f9937e675",
      "a7c8ced5cd72d8f9d6b685a44e30e66f5e65a69496d69caeabb6590e97e386d98eb9e1dd8e4435fc9c6d1aad5fddcb74b73ee26e4431a6ae2abe704d",
      "bd11fc785bb8e7ae07d078dbdcbb8553e89b9733fbb22b2a67f808d16ac56daa56f32a9b54bdb6795f7e15b1af317d267c5c3a695882f60daacf17bb",
      "bd5489f6f25dcf3ec68d62beb61f09b649a78c053348ba6f68dba536c16e200e9f55e6c9f620f7aa11945a78314bfbe0ae1857e510df887a9c4b0b0f",
      "da268745bfdd169141797d858747641ffee2944d59b03344e6e869c9ade68b8927751640b74dec76f0e4c83e15d18a73608635a58da0280622338886",
      "76c2700782ad5accfec357889998f84c4b9bc3dc7717d9d4ed218fa2785090172abd571c0fdb11275fbbd9a79d1c8fbe9898c2c56eb4ab7cf7dbd4f7",
      "d5bb4b4e8b0f644bc96fe52623a7b75d85f733c370e1d8af76d0e224ad938be49ac7f1366ca9878ae018b63df1364c2935574a769e36732b10996d5d",
      "e35ef6c186ad9e78bcfd3bde98adc89e8344f523c22a88f5516d250f42767e0cd94576ae0c4465bb1bca6970ed2bf32c37884d62f6a69e6eba2ee7a4",
      "6e9b81310945c5c44df352919c3090173e977dd081dfeedb9725348edecd74dd8a256466e0e587666ae4da2b6c4c13d972e9f140db7c767b028f6e7f",
      "2cce0f3bc54d6457cde7c398bd33ccf31b333f8c19d603d10a2b05d1aabad75833ef653332f5c49bbf33f525c18c69d8ecf681c8fe7f86689b9da369",
      "08775122737678e9b5728bdb0806223bb96f88e2a21a13ca57329376c8a738f1d10b13b5ae7f50dc8ca2215e37054cb034c250a0366bc716f898727f",
      "4334fc27cd27681329b1cd69123b4e5a68c16e4f74d14e5d79d676c9a3503b910f9b32b36e67287a9b0524e5ba85c169a0bfe39d360b5813ef8179e9",
      "732c6e896e1e041b2529d17c0ea8cde731d762c59df16010b43c9c6cbb0fec8d997f9b8e427256f7a3972dbd25eaebc330db76fe9b2802f83a79cddc",
      "18a8f8b029cdc4bb3ceb2b91b96584184c74622399a812eca8a8ba9d120d290c371b1b7f8d8f135ee3cd7fb405b4b44a77cdc79a2f72b424d3d4b4c3",
      "4c5ba0a7516c975fdf28f60d6e8806ea973571470e6411c21cf815cf55c74eb56f447c35b30bd8873b8c5bec81a61f4e733b436447f246348a7f4213",
      "2d7c4044ddcf5f13997fa4ab660f9b1b8bf24a6398239be6f4d93eae5b7c039bc5b4a4b4b62ce6b0ab70134d32674f2375608f99c84f21acfe7da4c5",
      "ccaa48235546d5908f7042b4fcef34fd3a7fb6745396749e1fd4662a2edd6fa4a3b37ad1dd56f1ad6f4668278a89687fb85fba87ebf1cbc6e09b9aff",
      "7dd2f7aa735597f83a9088dba43328134b21b388679622ad7ed3412dc3c76d76e74d0fb93c1790ecf8f8a33fffb5f85d7bfe94e1fffcd9ca57313b0f",
      "e04592cbd689697c1952a63c73a1b759640f991b50631b9db0f6c3724372f11acc0aff705cec845bb2a3a62b5413d9f9eba23473bf7b608c830c1981",
      "aedd7432c476d6ea415999259122b34c1099864c519f5d57ed795355145c2fae32d734159704e661ec1098ea197a5e73958d318adb8328c22e2429af",
      "a925446683855ab2c3ac23b203b585e89bb1e7d492709a7288ccbc5c69bfb72b51ff104d59a83dc23f26b0a3edf8349e149957393134579cf4879c97",
      "6981a5f0748acca5a87a5ee8ebef87102363b5d766031ea19d6c25dbbabd08eccbe5ba5de645bdc9868ec802babde8aadb7646324ccf9bf33f0941e5",
      "bdce34cd833d3edc50b6e75db4f0e6a45dcf86fb18c2c7b86c71da7798e906369c0f8cfd8c2c708c9b8bc735339c4cb4121fc4fcc0ee8bf372b237fd",
      "f0e40036860ae9df1efbf2ab67e5e3972d21160006f7cc98c8c465df846e4f4ecce94522b3b572579bf2387a9ea7a8090499c7e5571d5fcf8a0d99c6",
      "e25ed33f998996d9c826317ff5b9153601225fe0df8f2a5797d2e18309ecf1c94cff9947ac4f36abaf2e2231ade92af3dfcee337dbcee5939d95a273",
      "ccf6dad871a12e393fbabdd2c2e76af1e2179bd844239cc88bc9cd2e628c353f87b93695b9a8e114b93b2f51b96e135e8e413632f9533fef51d065ba",
      "c46d09420c36ebfe403f774c74f5fdf1c547bcf46966a832f630a20ec96362b4f389561c1687bff44cd3409005fdd014671e87fb5cade67543839163",
      "40754cc24619fea9e00a10f90bb219d99f4c69f729c5bd96a662af3c5b76907b65476514a72ac0816884758ef66dc7b762a47b16d3367551c6c0bef9",
      "f4895750a6da574a8d9078c1b1364f13c444586c558929363b290a3f82c87b9901a09d56fcc589f2da8b6b5fedde5a5eb7c9e468ed67e2f0dba73e83",
      "e14b999edff2f6478ec68cbbfaa60845cbc11c42bd6c8a0860f294895522d9fe90111f151c2462ec24cf75ff99b7fae68931129cc752e42946fe5f5e",
      "27ee64210d1a827dedee762c9fda1ecf74b89f68935dc30d1dbb4ee5b6baeb07e222cd7db01921277ce03edd0670a535f1a56c2f71d1f453c49cdbb1",
      "75ea9bc9ea56b1bc62d92fb3fa1a3249396ea5ac48fe98077ca75f55a72b5713ade11ed7ece17b73128729361fd1248086281419c1205beef239cbce",
      "31dc858908fec6953011aae16ecb9e3135d30353fcf0a4cf7e59a2187ae97aa3ebabcc3a338ed2676f9c0cdc3d42a24409407e30f8ede607e9af25e7",
      "781537078963d4ab2198969ae9d88345c2bc9b35c312b6d17c029babcea9e35b130fc2271f388f7cfd9df3601fdf34a1826a9180e0b39a02999cb7ed",
      "f9f6768010663b2767e19b75f56150990e5d845c6dbdafd42187c0ec005a93b1bd3514851d718c1da01a5d2156d391e9cf1ca046c461e7ce03bd16a2",
      "3d1c6879113108138c5a536599d7dc14accd0027868fe9198201a5de8f50b79e37430d281eee43283334e204e67074afbed33847dff2bb511f4d3adb",
      "f178b1ef8b43ca89b4d1607fdcbc95e96026c6a63bc3a7764674448f461cb857c723d63ba3df4a1b6262d4b434bf99e3b42189f694851d66e36961d3",
      "1d292cacb5081df2dddc5ab3037a87ef6b3a7779a2c377c256345cf991de4bc609db68ca7dcc38f8ed3f951e8cf001c1e9a5dce06c192b30bf5da991",
      "e6d3be88b62933523f823a07f77cc92bce38841d789ea55a116bc43d6fca1d9cf00df9e49a63bb65a4a9cd1a3c5893000c22886856242686e21f25c3",
      "2ab629b4ac6ec4e09722d39709cd80a16d0a0014c4b8fc632965a62c705b6e2390da18f1f922adb408bddff83227e333767e541e677a6a64c846e912",
      "3bc0c2dee05f95cfe89c69378359b5cb105e7be4c4b4ed9ec998f599f3aee5b8230b4072111b8144443ccf5ea2b9b3d24ea94f3c68ca7ebe81f4e043",
      "ecedc6792f19072ecd731ab7babbc62a502e868146d0c6661d97804d7232bcb23da552ef135de99e5259fbe60a3b95b033931022d3b19a28db32e97f",
      "98fed63632f57873452bf94ac4e766f8bdd4e3a2f06eec13cfab5d34a667f758e3b7eb139f5a739877f802ead5cf547bbb8521b0942f49d16d9e4b33",
      "9573f3fdcce3f4f5f9f182a32b653b67fc905578ab67f7c7e37535070ad48c786d695de50ad7f3d675f955a77d447be6d65eabe18de0880996402c3c",
      "f442250a26d11e3dc319f68eca87db7991a99da27465cbf4d73a55088111dff810bec5e82b019b538e6ce1c1d886afee5360781e4f859a3e4cebb401",
      "42f13552e32af34e965ccff8a8dcabf0430cacfde846db55884ee9a941e76a12d5481b90f3dbeea0da49f14458c03029a05732f0643edf7287380d4f",
      "3ba9aa94e08763b8d9f0262aab6696e51bed994ed100a321a7c2a6e43e4e850499c995f363630febc62ee06861e584b14d26f086d7bc1fe76ef9005a",
      "786d1cb1f79b315d70585dfef78dd7bbb1d8ed71233252fcf04e1fc16ae9727c33eab44a097fcc21bfa3c654bc960927296fe4ac80925460235d7673",
      "06cb88b3cc5daabba8e7815b54d030d8d7633eb2dea75c99b377feb492a326821a17d9761a45de9419742d1f6a364f97c3337bcbc1c9d7b82f879faf",
      "0a12491999321eb492c05089d0b1691ee441db07520b78e856e6b09e8f9fcdbd3d0b196e123fe5df91ad1e2e44277849d2a626e6a94cb2b3fdc96223",
      "6b6058272ba53ed791e2c2147df05339ea7c93c2f7b5b30fec7a493808f6ec95acae7476cc5ec29894031028289315f33590f5163887eb18707ca2b2",
      "67023ebabeada2892e62cdb1aec151c5ccdb48aa4252268b0c1611bb568bfc1fd31023af1bc2441c8edafec8be7e8b6749da9b5437fb1325040249b1",
      "ce1d5042107453999f93383bfbfc53fe43c0580015eaef980c3085dd1d9a9968a6422803efe53ea2311c98377b62220f06b29ca3bd3f80796603f975",
      "91a6e579a05ac0f495404a0a071539b36794231ca60f78602bb2c8cd8ae90981baad4dcf90e5eb7d3529cf4455e94b77964f0636b100e540b7917027",
      "c58696f6fe09a81ee0fce32141334f6a429b662d521b369696d79905a800630bbb952dd8a4d3e39a464a14b57fbefbbb9fbe1e50d76d6e19e6820576",
      "eab899f17cc0e6e09a135e28b2a9c323023a4c97b9665276a5b5e746d7639be70371b9a5eaec77fbd1105f4fd5b393cc78a87778b024906565c54050",
      "aced814e99572e1d13e73ba25da476c6becf2f3bd51f9b9cd147c04143499fdc5ea6e1f113bbd81957b0c3f978f462f40779b6ab221a6486d2bbdf77",
      "283e537359a046ab077bca4894324653d13bd9277b3e63e79ba0891f1cdb28ae25cdf0f88829dce2973d8b079fcb9f6a4ab5e53b51c9e1efde9fb885",
      "3070dcaf0576b11004a8d2e88c3c616b93e42710899d9a801cbf47cd5a8959f9fecb3c5f0cbc7dd3695e02c38a22294d0b658513881252a27b0bc43e",
      "52667e93208891d445a76d3b2c64eb3c5837b09b4597d949dc531d9a192fb4a16d9bf9fc9694821097234585b216ccb90066201bdda1507be6726277",
      "0a6a58f30dae1b81275fbe3637e56a7d26dda15a4fa4a70f99922b4c73d68cad3e3d28e31e5a476eff655a5340470dba3d9ac09ee71ed83acf29de1d",
      "36baf21af338641296301ab98f9d8572ee5a78c25aca0a7f027d8bb3b1526ae04e608f5cb07deb2f5d5ee0737b5b12c7fdc4f2ace25bfc63cd91f12d",
      "53cfbd0d5723f759d9f658f74b1f1940f0e9c6e99d8fdd304af8aae7b99d599432a67af9808b9a7f2f536d3f30d00970dd6dca1a8681bd81d0973263",
      "ebf3563639b2e9567e69de7ec8d2bc89ccaed41ab935389e85431ea474b710eb03172d65f9b2ce1201409f0763cf96f52c4f360ce7405662aa30e23c",
      "4b4ab19c544e5869dc0ea363a5befa7a2f3d5bcc36433fdfce816d80a9fbe54fc932f913cccd4940a879d69e8b8a8a0982f73a4f2239c2b2f6e32582",
      "fd05d8072f1d20edbbf418805a193bd7375b9943e8f3fa770defcf1eda5bd3a9d1233c2c13a7fadedffd3c1b5d331936059bd2f3538fa5f4cd2b7fa5",
      "6d343d8fa2f1d95a8f42d17c8ec718eabae97ea6d70fa6f1fc1d374b57e27cf1f87e8bd9d954002813b772c276191e6c0ff54e1e85807cbcfc899494",
      "b9c0cfcbdb02516aa9ccaf64ab5dcf4065925500a8ea0980c07223c3ac22c7d1639b83fc3d9e60e8810005f90bc20a153e225701965717cd731ee06f",
      "bcdb4736c9cefef126efebc9a975c209c4eccaee4f8eebb8683f22a13dc12e94c89601a835334708d2707ac0a1d94be2c9b59ea9903111cc486193f5",
      "a8087a264150093796174ea7f29d97d7bf71452552b98591ec0f5c722a401f0553caac178ff8b748345c4f9086b19c09b6ab9cc2493f1409257cd8ec",
      "b5679aa90268630c1199b944007d54bbe52e29b2ba2fab35374ab1af443874eef3f8855f4d777e0502b6fb4de9ea48fcd96fde8e12a2e18d5ffb2848",
      "6a061c6d216f9b29e136ef0315e6ef1521f083f72a4186570f7ff50cc00efdfa64aadc1e880517618b6e870a821a852fda0112ec404a92593a473f8a",
      "dc345bb7eaa216b9702181a62eaaeb6688ceef5356c2054f3323450f3f11fb5b1ae97dc2857a8be6434540908857892a1ac2535d8bd18ccb93a1c1e9",
      "22fb4009935ad3d60d41ec5233c637cc4572912f38c52eb5d203fda4da3c3e58243c40c11d3e7386a740e869c126d094c02f024c8c98ba4b2907b31a",
      "693676389521acb2009b2bec84f10aa4e8e655b6c0953e9acd9d0461045c143e1bb1b673943d3b01c66081eb10a556f864067fcfabab5ab94fc25efe",
      "1fd8bd9ee8825bddeb2c1fb1699db3ca78465466120fc720dd27814855319f45242f5667ad9e156d3777ed724334b5c02ddb2ddd968a41f820861360",
      "dc4399ef5602bc865d6a3bc173e6cf4d3d6fb8cb1349db9e677f203f8eaa3861c7d94f943d5269f723aa3374c0c92443df479ae2de95d1af42d11b1f",
      "f1e325a71c1f4ae69e7ceded6cb3f2a47cc86184f24c3834e828ce083443e254c0c5907c5de3013d494be608e7951a357321359e9d060b911993097f",
      "6f52be2b546fc749f4944dcd2288afd44c262448622c22a7e666b75876677999470d1f566b644954221cbf05235549cc8c021b0247a69f544e304e26",
      "32ccd891019c9a195fc9f4a997b7b8ad85a8742f3a3879d1d7a3fa38ae385abc46230120c3aba4a8964ab44794b204048ee296eabf7a700f83af84b1",
      "3f0335a59225bf735c72e3373b51104d21bd42ca75a6e914a3794a66660856eafbafc61b2f419cbc0a6c4575c8037fed5e8f180953f20294000f8e93",
      "95e0d3832d423af67ca8c2e250da4485933eb87eaf12134e16ec785065940f40f65763f81899294665564b4f94b659deca990690b58168ec82d955b7",
      "4305e5717cde18bd3e48eeaeac6d09b373a92a5788d41b2ecd12b04f29da7e9f0c8baa5607464a6490be1dd48f5592fe6fbffef2dbcf7ffdbb7ffd1b",
      "f5a1bffbc3bfffdd5ffcfdcf7ffce5affe6bd596963f0165fb1ee8d65f02e3be81ccbdc4d47d05ba7b01cafb0eb5f702d6f719fbf7091cf8093ef802",
      "64f81d0cf13b90e2b710c66ff18daf00902f60929fa0949fd1962f20992f809b2fe19d2f40a02fa0a22f00a52f60a79fc0a92f20ac2f81ae2fe0b09f",
      "40b32fa0b52fe0b72f20badf8178bf46f87e07ff7d850ffe04227e01357e0148fe8c5afe0c6d7e0180fe8491fe0e42fd0264fd028afd02b0fd09d4fd",
      "19f5fd3d2efc0578fc33c2fc050cfd0558fd05a0fd13e8fd3b58fc0bd8fc37b8fa17d0fbcff0fc1710fe1740ff97e5009f8a0648fb7b29c3091b9133",
      "ba374d5a0fac8d4e51e6f9663180ad729ff828682aa4942692724bbef9910e6719b835be9507d15ae9f9f6647a5887e5990eadcba7ea804cf4a0bc26",
      "2b13b07b59cb2dde89e098242baa9f341604b19f30de216f40e24c284a8c2f04e1b1707eea1e2d53e9a665f40b533c2b5604c52a0553af3ea08baadd",
      "51ebcdbb5226af8b9e55a4bf57fe7d597985e5fa64c3cdc51b94949507f131598c85906e2ec6bb5083455c30719519c7e3ad6329c4c3cc0ab66fd4b5",
      "f972595f926162b559697335fe6e8177904b72908c7b9eda8aa6bba42957cc6246b01f23ad30409869a2c8ca31b61351a9336f7309f3f58d4997e6a4",
      "1d975da208abe04e365d1415f0b80ec40ea008a80f86c05c1fde613e457e4d8f378d915e6a27e3c6e34fcfe28fabbbb4ac0c24b3a507d57d1fd7d68e",
      "d6cefa7f26e7ea5d97800c8f552c09241acf907957f0b72507dd41b274b6d0c9eac46e3e09a25953290200808ba5cd0f088529f14e72b83e9525cb3c",
      "a44e1d7b4f3041838980746494eb1336811e81c466b929a25897cc694428411c003e867a6094acf21d5405c219f35c637a00418940353ef0796b1900",
      "c7f9c0fc241f1721008910940f84a93d5a94232938222a31c7362515fc5bd3a997a26ee640a3506f255ade26da04250113c41cfc8afaf14608027c0a",
      "be5cf7dc59ca2a4a939d58830de00156849ba7aa394da23431433ac3218cfc8499cb21daf5ee59216067067fdf88c1c93f63b2d23914740dfaa1f404",
      "a1e5d82e0a1095a0115d3fab91a1c4d0e1efac119fdd3f4adb23dddd01b14a6792c67d80c3a6debb38a31209020cba073193bc6b3e66232fffe170eb",
      "8546b80eaa64dfbac60ed0070602e683d79e37af827de24dcc254fb5722569ce55753f26208c0988c6798cc4b6fcf577e22d0660356898224045bd34",
      "9c0ff05f44699e10efc745b52750de3c4f0dfdc1fc172d3b1091e97741628345b91fd74c64420187b73db8c52d743e42e04a0e3f7076add2fad56f73",
      "35ddad453ed3643eb1b8c295bc842b1ba64586beb5fe3328ebeb6e3e315fb87dcc72ae4f847bce836465c2974c1751aaa660e3ce78ef42ebb7afb840",
      "3c5a8ee4c11af8bd03ea509481e5ebee27d136f5b9c3c3543aae4c6df9da32c3679b0eaba0736e3e39d5c596ed183419942619896eae5fc693b5285b",
      "179961fc243c6ccc2eb91f92aa9ff5284b514605a58e28400d558036ddac3bb59527d24a68b548d055ac3689664b70caf05c17008d16b0848acded2a",
      "72d5f624b26a635b01529a4f46be01af061432d693b4ba665c758ccac4ba881f474a788f645f0202c02133e0d97a784f8aa93591d93d92d14d4998a8",
      "c72e212d6f6e87642329681a3641d3dd5644ef90edc1012127ef66fe7cc338d35535fac084e6e048c2734ad00519169dac3b89a020459a7eb667909b",
      "94bd9d667ede072217999d253213a2b2884cbeb96a327fc86ca6cc41da9565f1d0c8a082b0ceee7de866263420836ccc2372270b0dfbc14053088452",
      "511b21060de41eb946b3f58febf6f572e435cec7e86c4d03040f3711d93a98d18400c7f35bb33bdc09dd1f2c37555795679226a03fece87c83462c11",
      "5bfbf97bfa7d5fd6ad5501b239ba0a43f58bc09c37330521c689bd6ca2026e161cd58cd52da2277b7a25931dd00b6404d85949e117fbedcc70020092",
      "30a739bd36416bd00dc7aaf27453ee95c0562c6f10f1e0feb86e0684628a2e51bfed1134b01567c3ef7aee0e975d32ead119737e7410e900fd78b5c8",
      "13e97e449a00c1ece711d049216a811c46b6cf74d938cf6b7c337b5fcdeae523b32922eba711724e7de82d8aefc91c9968468819d9c4bde397916504",
      "585faa4427aaef25638fc3fad602c132a949804f0d76bd407aa856606b80f7ab97fb7a70cf98c9471c71cbd5fd71cb4d5c8e5b26ba266b0d9025f801",
      "c68ae9fc722d110cdfdfef798e2ddff8fdaecf2defa81ab9ede547463adcd9ef42134f129b5b4f79e6109e8ba39fb6e7835234b6745d4f7d0a21c4ec",
      "baddec91f047769bae9b3b17f5b7a37b064d309b3933c3adc65931156717855e301954e712d31d4f7ebd2b141dbf2c89b083866c38d91fcc58299bf0",
      "6ef084be1ed9374fcdd1006c5f412f766265831e27c30165c68dcdbf20ebb9e2b988d703af2d6eb8ab6ab8941de27a5c15056995a836f62777abfdf9",
      "ed374f7de6a607df09798b981a52c2104600c7885b8a7502aa363205615e5061437418ce89a48303b448e842f524e3512b93ac24d41709ce821163eb",
      "762369a9a8fd594d8f3dc990f7dde872d0ff25e46b2f19db5ed0ba7dcffef62d3ddc4b02b9173c73afe8e85ed1d6bd60b77b4582f78a2cef15a9de2b",
      "f2bd17247daf98fcfe8b28005fca30d4a66e571345d6217adaa2057a305536e026669530007ad7e6d572797e27f3452f5f2375c9bc1518e790254eb8",
      "43aee3d7b5841d434e8b60e6673bac6f3e541479f2758fe90d7d8188f6db5162bbca775ac9d39a60d33cfac8374a1eab5839a88dc0860e2b89b058f5",
      "df9e2052aa04a6a9713559f0d35405a3192cc7d479647397e39b3936fdd0b0f4d81d6c6a95723c7c24298f542e6833f60f7e53ee33b2d657b21df630",
      "942fe001200d6a7959010378c5ced343aff00871bd51ea958a989cdd3cfedb3c3299a9a3372bab3edbb654d7497dc6c2648f9eed443eeb63cbb7ee5b",
      "c44cc8e7e4c27de2fbefb462591de548858c08b2e91025d180ba38cf9150595de8b07c0b726c50d128e4fd685ccae590cda831ab6c6db02cd8272335",
      "2bb1365dd673c0a4b5baef9ad25229ab827088bde784e943f9da99628c1a25af635337b1fc9825f0fcd6d5328b394c1fa5a8f4d4b026955e737ead44",
      "6c411d63afc8ab961284817c008a1264719ce7a4139f2d345f355e952959aea33310a613b6b9529dcfd13c45267a540f94b232cc4899826adf67d221",
      "8d9ce29d7b3ec4d487381f2a6d62dc1dd755ede343ec1af7eb7982a083aa9e3b03d3c29bf5e91bb444b68149e9532b6725f29bfaa475aecbeaf3c108",
      "7279a1e47da6132a5696412e763285b9179f6fd8293a72b590a48c2abd0223299b71e601dc55740798e82b19e59b6ca82000d2b18c570813505eb489",
      "c632c6fc1d00e52b8e83642f4036a9258753283709255cb3ab8caad53b3ed6a1d7329a193d9ff50ac39a57203eef3f555fd5e5f6e70cefedcf309fe8",
      "594f777a6556cb4dc76f99278a3153d1ebb75b258d3ba23afa3aa4a8449d178a422c16badd9c1f4b7d4fff69eff9d31ed4a7baee63ebaca2ebc658cf",
      "748ee6f75b2777b1a893a9db2a27d62720008c68994e7991e9d1ae133f0d38513851fab4603e7cbeff9a6e197c2de31df5dbe7bd44807654365873ee",
      "9a43e37856308de9b76268e478c8ef388000883b3310155a3ab7a83aac8c9c62a0f64e90d8e779debfc15c46c9d8576f0125b3d0c381f764d5448158",
      "7ff60e61c5e24563f33cfb097e7a61711352ca0ae80eb41d238ed829bb4778fd73ebfa50bc5718e0d352b5e3e394e3e8cb7d3e1ccc2d5996b47104d8",
      "62d4fdd6cd23c00e5a7f409e304d20e1544360cd3e20c2533cec8cb17f9c970207f73c8e5754a8251c77405ae0a5a6fe37a5160235df3ce92774a262",
      "a46d8f455076d4d00553caffa781f97447006e0a2285e513f76312aa6f435be6fe6b259faa2a9457ee6002b7be3347cfef8756f702d07653735092c6",
      "8e035f975b89c09f53a8dc1ddb8659b8ae114c718c942d4a0aa8d3dca111b6037624abfb91ed18de281fd77df36a3f068cf81f7ffa8f3ffd87bffdc7",
      "df7efaf587e9dc619d74148a824ea0a0aeb95b58505e8c20dacfd31c31a0d122699e92f46233155593b586818db9153ad1afd01e20e13c23e0e3610d",
      "d1506c255dbbce98ad6ac94a9ce83a2852a90fcf45466918f986437a99fa507f100788d28ed3ef81d5a17cf95d71d78ab140aa724aa9502247b19152",
      "d68ac892bea276503977de0f106f572afc38c731a93bd547b150b47b84509eddb3a23524a6aa140e9b5e9d442a925227eee31833810d97b2fbbb689b",
      "81b3275644dadd93ec288e733dc7af1b93a155e42ddc571524344f3e7a9c81ac31803249148b17fff1f4109e8e66345f55dab37a5e176d8c3e27a2e7",
      "4ce562a2ebc3637c2325fe46d7518d9c1f626120b82885a76d3225594e5921fc50d7cccc225b99c35be0a791cb9d637ff0b793666268f8daf03fcc5e",
      "d2c298ce1533a36c89957115419c8e3ad2a4784c154f08d36688ef56d7c0d93f35b3ce05321560743c43537a9332c4e1600b674e9c426c4fcf526be1",
      "4e8a0bfc8b15b7d42940bc4d6b70382c5fbfaabaf1f62aaaa988b8c800abdf465a5e5350648ee029f4abd5de7c0e041f10cbe02a32c4e186d84781d3",
      "e96114d13ef84c79fd2a6182e559ffe1f1568204c34127cda70a97be38c8c5a98a74491338e6cc7847cc627172f815d59b7ab09d3d7a074d25c530c0",
      "19b4b60551abcdc15e8a7ba2308454a388b8ba6061d17217a7fdc5539bbb3b09a5720fe887ea5c1cdb83c34076960392ba3b1a1a451754aa0ab330c4",
      "fce33828bf2f2999280c9c3a37e5ffede3803899ff90518a8d148f6fc627d8602e61dcb8321044f87fb63703f0f7568c5bfc9ea5697b12d914a3acdc",
      "3a5d232cad2e718f07535fb49d04025d825e6a22f29cae093090099eb3303c7782a947ae866139fbef12f1873842b6e238fc68e1d261a04b21625980",
      "b6d63295d25c3ac569531037dd5e6d1b17ef203e5db1315528ae2892cca6a520834a76fb746d33c5c525afccd3b55882d4d8433ae364a9c2901c3956",
      "dd3bd96c65ade58f5567fa90ea2bfa5575bb089885ed433dcbed42aad38b3fbdf89755e5c754f9ee70ba074af800e24c99b7c7d5b79a6ce0237a8008",
      "2401714995f86ed7f9a714b1c00e5706478971ef05e1d90601a1294245d14ac7eb88f24b3c0802cff289ef7a4320de1d2e392ea1e870e92b9a9531e2",
      "6cacfe24f3d1e2b49cc297dc0c2d1dc04e57943ea7692abe39957f0c23f91fff8f9fffc32f7f9629f052ad7eaf7a3fabe7afb5f7779afd85eaff7c3c",
      "7c3e423e1f339f4ea25787d5e703edbb33efd3a9f8f9e0fc7cb87e77fc7e3a9f5f1ce13a846474544f9763906a5e67fc3becda9e60d3403015b663f4",
      "828a8482b139050d3a2b2ed97b68ac66771f9714c1ff9873ad6a0c509d4498e88a3d6ec5d374cd70849e606b00d74ad2a8e3d20824553252e49d1c74",
      "4d7144852a1e971047d3115e80d0781612877809dddd945e34fbbd86aad88e5cf29a64b579f044a5e22394670837e420b0a5a3e220e87e88e8ab1ee1",
      "7bba43ab900cb31341fc6c29fbe5e072c1965a982995629a0e0b88721a380755682834b32fc32efc3220179de284334ed3972e8ecd40ff426a844de1",
      "40223153af008f7517004c1200cdb18c8a422e1903cd43acec8fa605c3a98580d4ed72ec9b87a5c552ddf424672b9e4a520ec7b1b949f1dddefc219d",
      "f0bffcf40ffff4d3effee79fffafdffdf7fff4cbaf3fffb978e9cadcc03ff5dec44345273b5eab4f75cd6aae02f10231b0df857e52cf3eac1110f1ef",
      "4384e6213a0055df05fb9d21ba875e62402aa5fee5734234f64e8c217eb85dddbf93c1f3bbab76ea22d97e242b74acb6643ba884142586df9f7e6f43",
      "d61fa28d51437bb7e2bdfa862aa5d404adecb8bb425e0b91996112c91a3a6a8cb61c9b2d671934c03b06873630571da06974656a8e82533018af51cd",
      "9b9c8648776bcc572d0a94e4ed4eebde41eef61ec385d6cfef776485f19eb7fa65d5dbafe9d529e57b17574d4ed0a47884266dceb584c84c3b35bb2b",
      "cec4ce67db4410de77163ce3bf9f6b4edabb38a57a889a9d98fcb079393f227a9fe9f6a7d758154bef846a592744050ec8f7078ea4e50461a1faebed",
      "9b22b8c3dfc93d3f92afd75cf47714059d3a06ea5618b55a843212f8b33eede0b8563569ef22431b215afa8443ec1a211a7c1a1cf0eb5fa63aaee21d",
      "f072ab21028eef3de166ae06756583eeabce90d024afc5d329efd43f9ff8f3e896251edaf41aed1956be55be281b9e9701b25472c5889ed377c50991",
      "d952ea95176106ada1ea82d85fa8c3e29b29560160cca31db71daf8268d001804dd89f1ba95ef49d23bfe73e3197764954662e59e81e11d9dd6b0caa",
      "520880c86d478d133e7e7ee80c8fbe4fb63a3b467c81772e84fdde67025a10adc11ac6addbf92dc438f1ae42f113f35769ee80c8bb7820eaf1a57b8c",
      "a1a9918c548f1fcffad2a7f95b7b159c3e1be5188cd459cc3e7d8afc422c2b84dea0d4e11ffe8acbd7be101cf40d24d618df7069c1ae2f337f74d4ca",
      "f2b4f873f1d7cdbfec057571f74dc12329385437cc7242e43d2041807b5b546c7622afccff8805f6d5587fec0cf997bffef51f7efdf91f7ff9af7a78",
      "bcdceddfab049568b9c6f11636dee9cb5b4a3a9687c4f1512f56b1e4e877b419a5308d7e9ed58f85aea62bd267d35961112de2ba9c140e31aadd3940",
      "f5b5bc86827b5d3d50259c9ac9ce697575a6993d324354690f43dfc9e5870058332d06b2b127451d5e2f7a513a3d903a8b1de92bef9fa5ecf3f48ea8",
      "dbebbf5dd48fda533a999ad2c0849c6847e960767120419e48f7497711116d1738664963baaeb19adb5d4a1e43fb20911f41749becde6d72d59876f4",
      "9a7a8642123143a44e4c74a0ecbe97d4a692fdcc5937f3f6945a23ea1e90651017de8c77b03b31ef233e34e671bd2112ab32dad35d2ee5f94154be6b",
      "1afd85284bed6ab0d97d778119803bf39d8c7d9ff1553b9d7a24f293571d5fcda47e17156d0fd1050ef3aeeaf7faff7200c181b3d0223232350e30ce",
      "b20a7c010af2cc194bef3561ebd4a2adfbf1efd605a443de84557e007f8f30199e7fd71d6efc756563f8e2fa7e0831b22e7628e9f9e3dc1c9a5aa8af",
      "df85f49c216ae447d0fd716a33db107ca305475e44631caee97b3d9221e5d4c21cfafe79310c783e305855e1c475509ec09ed255e322c9f12a393546",
      "dd71cdb195ff56d50937aeb131ed37f80bfcc0e09a6ac61d2b67a464903b8799c9df0e42269c82770e173f41e063521be291cf429f76df1a8eab13d5",
      "89ec4d1162e55585d80187f274c950cc40dbb5a4642cef683c6e5c446671ba289a595748ee25a17ba924b289aa2e7287591444c5d7f2bd23452a3fc4",
      "6618d745b47ff1a1873aaaca5aab5dec740f52dc40a432de091b8e14a9e4ed1db7d037620120da985fb35e465e4480e7c1f0068932dfa9383b84489e",
      "a1dc2168e5e725162de74f112a3ec6ad8246c24efe2bf12b996ae05bc6bb4933f437fae7b8cd49bf5778246b425bc4b053cc6daf49cf2a1a1852c0e8",
      "a1f80253f4ca2842b7c1c4b4650f1753410bc911a399c810629114451fd5b92d97d6c79acd9dae949f37fd6da908bd2b33d17f37c068e94b8e0faf67",
      "7a1b716040546ef17543890b8befde92e956dff1c44d8b2f94a688831a90028fe5aaf8bcea0b39a5aace0855f1edb072b0b0fd0d3f395c7fd145f6f8",
      "c1b1bce24d44607effea384051b7151f85dd38456d5cb77ee773abaf6f9e8f546b7426bf7acf5ef1be1152ae4f4bf5e6cd9fdf70a9dd5e01252f042b",
      "9c22f98e6031a8b57f7e54e16c2d493ca34771ba718d1be454091c3f45f2273a8ce2ea450519973ba79268e5ee9d6f181e35277b8139348d13bbf9bb",
      "f7c9b76cea40f046878899cabc6c5fcd4ef22d8d0f00972c540f8198acd51eea3c2780765bdba9d317072186b13bbdb435c6d12b99a9041baeca2cd5",
      "d6f8b3c984e99ae6b5714814f97fd789bc4374b568802ab9f1a19e9ef378c7e83c71bafc1e38315c77d16a591e3ca88bd5624ca735ef34edbb1bf0b4",
      "4eb7aa166b71f309bce89d1cb81f208dc4dfd0691de4c3889ad492e814673c50662010849e43a804d564a9969c165613a2d5478a445af00ed42ece44",
      "c81c6fd8607e3a349114c9523b71cc3751eea843f8cc358b9de3f65c3db1fed431b0f9ce89e58632734369c501d58a33c04a17d74774bdd9f8730d9a",
      "4a3e4173131ff4b51627808ad0e1628190e7607ecd0dd182c049adcc6f5e35c8d0eab019292a9c7f885c855435703ba11de21a1027ef204b7db5d6eb",
      "d4409c0fb3a708926bd504f55c55f0e7101c5fcfa410ae7dc8c03473d89ceca99e1fe1ebbd923b48248d6edc38c64f87bb864915d189115cea81309a",
      "56cd310cfaf3126e59ee2600c7de6118ad67b1ef3ac376f23af04f4fcc81e86493dd6473bb63c49403cb5958cf7abbd4b2c8c918b90bf6d4ef689c17",
      "efdea63b196681ee7c833674abe384182e2a8a900cefbe16c3adf26a62c97f37aa34a38e438fdec5029366d3d91ebc69a7e4434d3b4bb44bcca47854",
      "348efb71e0a9c81b13d5ab5dc43e398a8260a5fb16e25884df8ad17a7591e8de0884307e3f1db8bb4fcdcc05f7cd3063f4e538efdc3b28a316f6866d",
      "eca968dd0e43ac3854f55d4cb569ceec7965669a73d5c23058708ee116fb96a5e9fb6627e0ce852538c4b58868a7d1211649770d1de2c65573f957db",
      "5e57a306c3544b7395b3a488ab112e34fdd0d520a480ade9abf5d0088ae42dcd4d75be4c8c23d3ee5b8388900d6f4d0b22826bc39504e553240eb9ca",
      "e193d45d99b6947f1ae5df88ec1b473020440029a766d0f4d988ab0a29d177411bfdaafda5b97a0927acd07baf7944c81308df7f9dfc664469a5fc8f",
      "bba790b075ed23a2922344870e93b85cc3173995659a1ba8a8dcbea5b00c12d777759c3e71d5814385704c448e1051bba828c0737b116c22da7e9416",
      "c2c947218fe26488225353306e643c1e91290157d9de484a3cbb542ae2435edfb545cec474c72f7e5844a5c029e8924bdc57277575f7080e998d921b",
      "da26f778ae054f45ce03bc34784f82e7351754f35adfa875d37cd036881a3a5e44b7980a56bc29917a5c40b9f34c57fdcee85e890f11cf1878ac0200",
      "741f4557db0d1bac1c882bc0cb7a030d95ef01aa569c4ef9a68b66aa1069fa89424122758b4f8c9d6a443a6d31737e9e14f1171e9d03b1e4e86a6a46",
      "023e979ff545f4d8b224760dc7013faa2862501c84f8ddeafab118d45ffcf16f7efef5b75f7efdb3095f94de635911fefde89ac836091669d25afee7",
      "4916051405c9e5f9b4bbd4dfd9c580bfec7b7c7d31e08ba46020bfc8e5e7ae646a98547f4a96ac0cb04fba24f06adf0f341bcba8490a0e23193e473a",
      "52c0e146390770b42a0508aec0b2d02fd9bfc42c6d36a33022795d1fd765b6dab3838969af2ad92d4fef13736cb99f8096d947a5e958ee585f273b9d",
      "e8dd088e0531b54a5398528a36a23c46a552b257d6c3564c56ccdd0d71a1dd8ffe2a57b2a75145e5c01c92dd76928da8632411058c4ec64bad8d191d",
      "ca78d4e732be1e21c5fd74a9ad0ad981dd5bc980b4073f25f2136f062c6feba7649357b6395597097e7ba2c79b1a9dc64fb3fb2bdb1b25480664e62f",
      "7be1f30ab876f301b45dd055d12344c86245232811ead941b4ebcc133a2d689d604fd553157ff8e877a1ebc6d7fc781354c63b4c27d1344a09d690dc",
      "a4cc6bd49f239b81d01493975e5f7c0bcf2f951c129553cffe19f44545648679f687bd478f3889fc105b897226c208ed941d7876dec92dac87af6403",
      "c879276db0461296d42e3700b28396bbe5c8f6536bde9a5d3a8807cf80068aa87087a4dd9b5d718767369ab8ab92ea436606b028af0886bd167806b7",
      "8f5e2588aac6908dac9b38787c58a33e5c3394d120e3a3858c388b5eb3d7934c8b4b019e21c3ee79001454c85a32b909e53334bd3b3ae46ab8ca6d41",
      "d330f2151464184a5decec2e0d133bcbe1ee9e3200e99245244d2da73d2d0b62e36159d9573e889d32653f0439fea5bb1a0d3d43d1ab89102f5f435d",
      "5e58ae6d3e8d87cde73ab121eaf8f8edf6bdd49d4f8031037d660b9faf8879e0ffd6eeafeb838d664b3467feb256d91dd4617c103d36e9123a549da7",
      "d5ad22797da98e2e65c7b55acd6548512c734c1835ea0525d34fa13a0a460f627058fe5da8d0a4ac995d1a1b06bf64ba510c0eca87be923ea7e8b508",
      "eaac9d64958bc29c77283aa28e582cd3ac8aa6f866b6006a57fee752df865cec5759472af957f69aa1e5d4bbe241b925860c2feaa7ebc9dddf95b955",
      "ddf6ccaecf971e3fefb844673d5bb828502e00c07d5ae52876599f8fa89203c5d8559d711f1d163f8d8a101a04c912a54034a7985817262c3fedc957",
      "b7a44cf0fc47b67a3e8032de292a0f8a056019245d788760c152df7385cb28a5784e30aa2a4c26bef73c4d17ed3f902d8f1cbfab26e7f8acd7d6b3f3",
      "37384dc9ba9bdc6a107e342700129f36658b522f5d971dd554e1e15f67b59ab2451906cfddc1832bac9c2473f8169be2c3eafe116302d4c71c25d970",
      "117c9e1c0be723190e807c9f33422ae249a939e276714c1b2ded5a1ee2dfda183f66a3fde5dffef2773fffd7b5cf3e5b479f6da8d796565302974ab8",
      "67af406529d98a96443c604832cfc3d3d5fc77e33efbeeea24661283ec1ef85591c6c724ec4fbf2e08bafd0384420131c1eea11c237aaaa9c369f3c3",
      "e8045fd796415d755d193b4f86d187ff767b78defbb6cb5d22fe38527fccabac7d7fcac7b6d756e919b7d73c3cfbd09855a8d5ff5f4c4d8a1ee0147a",
      "27959213539508735f3f4fd9ab7c23919c53f3d4160eec5d84f01f3db3a6b61a47e1d3cbabc1b3c19e1fc1948c15d08eab9a00d2fc2955f34a257d56",
      "5daf54dc2b45f84a5d7e56abafd4ef0b35fd27d4f90badfffdd9f0ea087975d4bc3a925e9c5c2f0fb85707e1ab03f3e5b9fae9f47d7d48bf3acc3f1f",
      "faaf8d83cf46c46b63e39551f2ca787965e4bc32865ed84c9f4dab5726d86b53edb349f7caf47b6522be32255fd89b9f6cd2cf96eb2b03f74f19c2b0",
      "40bb2c579d27fb5115d1338c25018b83d456946e4b0d0af783feeacf3ad187901bbdf20d26a057b5058d2a55a8f74a3ce3f4475d7e9651f2bb5cd692",
      "614f9d0278fdd59e7e70772f9765ab16ef1d73727ec7d38df1b8f3b17bdeac0affa34284684e46db2a6e07be304e64d1fbb33d2984180f63dfa2f30e",
      "4b27084c00a94f42b98314717f98c551a743fd376346366c60f20293b21b88796c89e486d78cf0710659c399adf5ba3c63f0782dcf9412aee7be239d",
      "0d5b2e559bf366e7704a5d859fb8e2c90bd9a4f8847d1d8914c98e72d3f0223f3487e6b76d97b5e7c853811b6a6706d9aa4818797f0f21a48a19f2db",
      "3add4a4e6ef6e57973ea184aaee20e14fa5d746c3377e23837cc94e8f0ac0e25320509a996946d85f5e6d356030034d3ae5e488f957ea5c46b622bc5",
      "4649612847d69e8f2ab6936ae9688b12722ff41c1295e438fef69cfc3163e9bffdf9efbffcee7ffdf91fbefcee2ffff0ebdffce18711fb9c2f0fceb1",
      "899f0cf046f3bf8603c486ffb51c6ec6bf61d9116456d0b1ab624471ef577f75c0d6ca52c85c681c81cd53195293226e1b12cca6e8935a7c5dc7702d",
      "69f5a69de68059ff2e5c73978359f7b921b2f9267be8aa3d7ea7305685e156078ede18ed0572c0157513d5bc24c17481444016875c34bfd5a5b5b444",
      "39a8d9baff6edd13cfab9335450f29af9a97d9cea9caddbdd5255e802f96aa7ef588566c0206da94945367b27c1d1d14ea4c76fca2ee0d4071466825",
      "8aa8399f319c5782ada9c9b5741dec5b329e6837ad0f0810a8b69038c016e4d05c8e253a6e733419973e84e2949ae463a55dab3a345517550ae57495",
      "1f73785725d0ccae4af0e836ab8b35a6131389f40de9a513923b6e2447fc2bebd09fc20f7bbd26129d775b8e0712e146b64b14816aaa62f0448febdd",
      "a62deb4910e724575ed2532c51048eb159b79e55bd0f10393121b77554379794ba9472b9c725476962244e31528952b9e0064e0817aa3b52dab39680",
      "f39574f2eeaf15c6b2a29f3407eb832c12de1dffce734230997a3ac6eeead82a8a843d6513896a30cc43e9dfeb4c2135c1f580ff23374f4180a0bfd5",
      "4b6dab9a24068ecb93c6380d4bf99484a78306147a31e86159f2b000be0fb1e9b9e4cadc65453aea4028004f9b1f07c779f769d44d3c1bc49e67d1e3",
      "4f653907fb3704937cf2fbc8960ee0743cd3df9d769fecdc9449098967f3a97283a2cb3af66bee6e8ea1768e2ada59c99baf49dba78fa055aca8b47f",
      "2696aad2115e6da82ce651826204f272036c92a4b8290c02880f05c6ded3800409062094e8cc5085ce21dd6dc6d9f125532e749c22a0724127824a02",
      "d97363b806508552f31402086299afbc853ac2c0e1e89088ebfc51846565012a954d2192e92ab9596233833a9c05115b74893f5e0c39b9d7e1b52618",
      "d25cb00b559c498a21863b9a25926dabae44c0d8bea90835b4d1a865be3923aa0b54d84d614073a556893353f99f8246809a0f5d4211aea95b34afe9",
      "69cb7510e487ae02b66a3210382806a52b6a113acf8480f2afac72e5f4a0a4269b5038e3269313adcd9df39437ede7c6b1b5296862051ffffbb2ecd4",
      "aa61fb4977001ef3896efb1382bd058233318352e7521b36ca5229167e02948ebaef1bc7a7ed3be03fd371fadd0baac8227abc094e57f038f8ef3b0e",
      "ea43e3f977f5f58d4b26d4b3007b4e5e6396b4f2b1dd9b083aedab903e2b142dd69dd71db4facc8914a44cfc30115ce1704ee434b5ed1508cd698abf",
      "b2307ec8b2f9d77ff8fb9ffeee4f82c5ffb445438d9c227cf8954a488948ac0bc939a587c40d01e1200e1df34ed914c7764f3a1d8c677a3e438ba4d4",
      "8f382a2f0ba439552fc586cb53af530723aea09880f0f0a3a878430e01942dea9b15dc3b9ca951cd7c89d7c0de58f723686f548049b7b9c03e3f036e",
      "7113fa4b620078c0438223c9f07ed82f2464a26d0154c11dbb4b36ef071cc12d9e29848bae59de3fcad3236f404cbd0f3c55ddb586e4a9f3668963a1",
      "0c1d87f29ca0a454966ba6c4f42ac549ae69a7908e0d4973b24defd13ade60237fcac501a35273e553c324d9ce50897916688b9a4d65f031dd839acf",
      "24e251d1129cb82277468067aeb4b5a317a8a55f0210168f5caada5ef8415fdecaa43906f83ae5f810c8cf01718ed9ff7eb9fd18ebc6cfbffef38faf",
      "70fcae227f918ca3b48210200e31dec26651bb701d50ec7a58cc6d724648a4eb47b8567ebce9c4806cfceaf5abf30710199262d1f901fd2ce80f114c",
      "4a7b2907f806546fc6e5138b9cbad31e8f34d75987d61454b08b6c902f345cff11b22e1c16c71165248221be2386e67f1e9aeef08db5ba52203223fe",
      "5ea9a5c790de879f7dceaf0542804a09bb854e6d06bce6f8a6d3ff5ed0f979930e171c780851dbf20ea84f9d024098035c5da0e6d50ab71ed7271b2d",
      "2d2bdbd744155bbca84e5fea974f1ae8938e7aa5c608220973e2dc35ea1ee0189ae1c014ea7697b0dfd331d44364ce5e5d1973268a4c41c6976341fb",
      "f5fe27b2a9b4a8bf5f613fa6b8ffe9afff8cd27131f06eb4650b5b43a0edf376d4b033fd1354c2756821b85362ae28b719b6c7c4cbb9f962120ca99e",
      "12a73891b12d85babd72050b66b946ad5ee2a4c0f1721deb3dc6a1f5a53f97ebc6b86640e25f45eb15a611560c76f158e1a85c2ac90081b9aba02c94",
      "ffc80ce57878bbe07dbb3cd7f0a5a80706fd7f565c035d73050bdae3361d58210d5ac3736b5e9e5987e848fccd51d16f150bde7d5c925ba62290d412",
      "fedec103aa4ffd3592f9d68a6c2ed1eceab0e0bf9fd3deff937aa3f8afe51c3bfe077d48ce7652dc131b926ecca2d3edb23dba48199d28b884f21a80",
      "0fea49c4154c91308b6206afdcc6c0fbdf000d3a7f4e57ef5904c3a922b92fd90a208a0e759032b0c157523b0a2c747501bc48b0765c6266b3a6f778",
      "18acab39b2a637f73b59494a153a208a1e7a88b3be3ae56268ae054e5b9b3394dbc1436219d5106c9af1729857ad993ebc1d1847a317cca2930945f0",
      "4914c040602f25f6c0fc9b640bd6bdbe345adbde1136fc117738703d8bbdefba80461a1408794a832440612510ebcaa71672a7d807a3d6afaf2122a2",
      "9d43f518b555ac78793ef0a7d3bdd9ad8190c0b1f7a6584b8f1ff5ea9b6b794b776277681f850a5cc9a31bcd2f118653b628c98c69d639866c792ce9",
      "cd4edf1ea8576104fba63c4100352ff49d3e42303591e1492aec54b8e9721ec8a6c413577427146edbabc0bf161c50411f3f61e96ab29793e4d15362",
      "367d0c1d89382666e2bfa9a5c176d745d40e2aaf5e19d5a1809b5edc55732dccfff597f16832411d7a0b3e1320a67429c760b68227fd9ae7e0f096dc",
      "ded896aba6dbf5b53afdc104ec4f7ffde3fafba536f9a471bed74a9f35d74bedf6bd02fcac233febd1ef75ed2b7dfc9dca7ea5d53f6bfecfa7c3e713",
      "a4e54edf35e358945bb1c8d78ee321d6a36701d0b66085ef533746a5151bbd78fa1dccb64d0359524759b3185801cb238128cb63cb8800bc0764b603",
      "2788366484035a5235018d20d300ea4a74cb1fb19c920e04e3882bcc83e00a9b9708a4341c47b4ae03344196d2a533299088cf3459ee66b245088726",
      "c0d42a644467d0f1978094976872a0b389a6b7632126447a89bddc8242815a09f9452d680bc80a3039f746593dd18e2bade70134c2095d9bea447c4a",
      "d495fa4c618ea2f36fa8b41ab71d402dbc822c9e3c25e9045cc733dc29c9ec8fc4fc3afedf39e33828c819b086fd6bf36fd2c28aa88d98492aba69db",
      "de1c8d0ccac08fb41ba53d5c3db5a7a20f2481ac62c7154970cf63f2e6c739e7778c908c5515537a89601dbd4ae99922ea8f98632a59deb476e3e187",
      "9c1646963b0c9ccfc3d9e80388f0fd16ff21a5f23ffcf4cb6fbffcb0565196b636b53d08834d391f820e8acede4c1729c145d1ca79400e974ca3ac8a",
      "f5248b0ee387b460ede787473f2ce3c9939202aa27f80cd5f18b941d81dd319f442c4549ba6867b6eb1257e54e3d738e804095fc79126c7dddab31ad",
      "913f1b308f55253fb27bd842fb13f31bf33c09e2113f7c207fe65be963d400c06e79136ab7304eff2a638c640538930c1c25dfaca93a33133c2143e7",
      "8b07b0e9fb69ffb19cd41ffefe975f7ff9ab9f7e5562eacffbe6df7d94975feec5f77db10a5ead95cf2b6a9eedfd36a2ab12896a202312954cf151b0",
      "acc96ced8164c02e8ce849f05167aef61b911a040ed4bcfb467648598aade981a5ce04161d541a6a3dda230033a9dbbf7a90a83be1a85eb12530c693",
      "645a2ead25b470f952cc96411a688b8bce93b1ec5d0d40da78529d93e07a556794fef428f4078eb6f22a3a92ea87e7235d4f609f88f97d20178b76f3",
      "44bfe707e2e2abeffa436bea7ffda77ffcffc11ca18e91f483a78b8657969cc25f330a16f94b501329595a521627aacd54c6f777708cd34b12c335a6",
      "08158f33e2c140b4a7f8cd96633148fd7b97c9e58161e80701d1794c3124e37a3bc4ebec5dc4bb2010573f448f217c9628e5e464521e9381eb0820cd",
      "1023c208c2c4ed8453cd4d450d6755e78ed269a5aea7de7ad113a400138293aa3825108491bb74272b6d114d9bea44da82831796feeecd46dd02037c",
      "33fd022fe105b6403f2541104a4666e95e0350e1b4887a805583edcd4b0c4122783749cfff10c0394e83e95ee00c5acce5472620252afdc479a9d021",
      "0e2787a8482f47dca37535316e1ea2247cb6bd4becf20e6fe224db4e4639dcddb91e5a19723847dc86ca341894bcd01b89e2515b5cda2ed958617c99",
      "36e3ce1b56c0e145db11a2b9536d739b973613a2d94392e88c2666b52282caee45d948e0d1c35da911c9a2c5fb102345bc6611c3dda3fb6145a3a3fc",
      "54979198ad66161fe1e018de14959aa8279d69bbaba7a3182cb7a3c90805f2d9609a2c3b427d761a3af764c6fef004a1dc710625da5e718081b8f76a",
      "263c5562bc22dc9172c3e9acb4aec0a03ee93fd25c7aca76df2ea9eaac001542fa8f944b80582a0a1ebf124ce7b81ceee157356585ded1e7b257f704",
      "28ea71762b5cd942bcb967161f079878e9c412ac239ddb298917ab7571327704ee3d7411370f04d7cb86bbfc08515656a725224d34b724c711a6bd78",
      "f74970a2e14653186baef6548545647d204f52dfa41b608e4ec60aa0ea7e323a4436d469c4b118a2d7db742a8091337c5e31fce3d97ac88cf4a57a14",
      "8483b69979b50519deca822ba0179aea2a1bd9b6d18ada7f6ccf14887ea0aae14a8f04daf610282d58bc79f20b89324e6a91519cc7116280ed0cffda",
      "fddfff2dc8b137cbd8f105aa78b7fbdcf9e52912011aee4038ef7029beeee96900926e0342f791793982099880ea1c71763e678bcddec318c4173009",
      "c195c644a9ebd510a7fccde91e843039d8574ecb143d3f8c1923f28a1de409c7ba7776e3156c2e972a9ceb23f846e2d01211a7970813a9bc19b8aed7",
      "eed0a7b446cffb1eeaa2998bc4f0ae2b99a1d4c07d533ba21d716b7a5cc183e9684df5df99a2d80faa4a5871e84d405b8f4c235cc8b19656500c8646",
      "6b4896d7db101b0368b866b6ea46c21e5df2272338764757e782e13c31c4e1d84db0337a73d3bee52dbfa97bc70dfd8f9a511f831129023b10961a20",
      "ac918142ff7275f4887f4ffc401a1d789740c82e2177fdea47d875cdc7e7e8709d85f0b213f35971e041d7b7d4f623728466c86bb68a6f2f1d92adfb",
      "e4b8790097d2d12ae9cedd37c41fa0ef324b9c9b83e65ccb7bd2f8010d13fdf12a41d0a5ad765f6891e7028a2a90bb27a9c87c49503c4d8fc142ad1c",
      "9fc4b7814e9f21c9f1b886caad7c20b3b5a0986e94252f4561b753416f123258af730749f6a05fb7ec592d4e51ef4e0d3ffab58a7c6c796705cfffc1",
      "17d4f4a3e120260a01fd2ed1f55be517cd9b2c680ea281966efbc115abfe0f244065b7a037baf649f3bf66b40d0846cfa20d311d2bbf80a91c6dabe6",
      "f110261c9094c6d9824add6c458d6bd5b992cc78e83ed53954979856d4c5a03b12109dd07d0be3363a4f6983a3860deb013c56049e0a36a75b737c1b",
      "e8af3139e56e301c2aed717dbad39d7afbe637b9a83d98946997d0bd69e156cbc2c9bf9bfd169cc9f4481475b2568238e68f6eeaa60f00e5ed4f8976",
      "ddee6a7a470aef25b6d5714cfaa6f88e01284a006639f5994b00a4b158b6a32b41526b33f4b8ebf09d501dcebf4542a60ba69fb0f2a9b5c3e24c0351",
      "0bf129d3e70e8a5a4e6b65ef55e29281c3b2c4d6e502b3707c5538bfe85688dc75a1fbbb8cbeb81a0e520e5c755be36a6e11ed05ba6b5d266ac80ea2",
      "37fca5b3469321e812e1c2791587ad830b9ebae4c80606ff7bbabe889d83cb251024835a0dca6aba5a5f3dc87cc3eb975c289345f23b9c8e9a4e28fa",
      "8cde3191c605360f02757b5f033be8f841364258b5e8fa795a9255734eaa0d81a2a147edbf76942d1fde047fc7a9939262ba494249b94b0a6d9015e8",
      "b8413a4d6deb1b6feda0ba830f6b0ec896b655cdba1a2eece8d57c0571d2dfdb730657ec5ef62a108ceadd2eb129fed64176b1e7cdfed8337bab5057",
      "df18fbcaa63494a86fa2a9887ad4e7430cc81b91b08ae2724ac5f8805b1db1b2b89cbc014d1e82e00f6a0412055b9e5356a553cbbbd53636f913ae76",
      "0dffef603cf813864690540915fabdcd467746822a85a52105a10a61c6aa3769ad251f85e084b41b1f41684660ddec57f563dfc903b6e93db805f7ba",
      "412b5109f1f17576b04398fe79d31e8a0b96fe2c49a455a1afd92a91480ab20294884fee397f44b62617a21e8c6764a4ec9800ee7d9d5f02ea2cb597",
      "91417f924ec63490aecabb2b9dc8b66ac919547c8ddefd219a711f7f5883f0403759412725ae6409e6074f54d5834a0f36acebdc166c8e6023c0125b",
      "52a427aec05466c73939494755b1e88800ec243804efb5c56ad88381f010e260eb9e959c8e933c94a21fc19ca7e0b97ed8e303cadd1cdacf2be82155",
      "b2ac863ce3ec1455f84ad193dea01522cb41e440eafe26b9b0944912a1765225dcb73c34b56ac5b59f36016a1e86a38882f7320afa06a9d5a258f77d",
      "ad4cf1b375b5f371360a75092319428ccb6f45430022051c6c41a941a32b401f4b1d0b2419a8685db4bc451ead95eef0c3d7b778158ff7d4bd93d698",
      "ffdc186a623e71b0251545cbcd652760475468c0a72aaa96570d7daa12e6459596ebfaddf1c6ac1585e367c00906b36dcfd1882220a4ca347f087f81",
      "408dbb29ed705c34f2b7a426654276df04a61c46756b757959135785fd1a85ac1e693ade3ecc97cbc66dc36821f2e3c8df4d3c77c98b992b20a5b872",
      "475799575ae3799d163a5c15e98c4d1290b3b03ee497cc04449a4bcd5ffdf638d06f2b5040ea41748f04dd2df1aa5654de0cd10c24cf1431f0a2ab66",
      "ec4d7c823bbd7da6c39c718ed48ecdf96a9698bae503eee0b95aea31a9d66ec7d9f2ab3aa96d7981cd636f62fcdfea44d782164ddde4d486f00455ad",
      "da28d9cea3f8d409fcd47e8d308540422d7e56d15e8842bf2e058d6734dbdc79550d07dcbd6b35be43c97993ac1ea2d2aa3a4bb46043a6b541a74b04",
      "49477f1bc64c100c7ee0d192b87713819fa2900b8987467690f55e4a318934788f31753ea6a658fd316a5e44309a9bad4816450b6071688c7cd4a925",
      "e6a9c578b44ae8bc1b84597c9f49bb42e2d39107c4d05073bd1e2f8b7947cb34dce9482672ae7435cd9b0ed7d2a7a3ba85a66b6eb533978b8fc07a0b",
      "ded02548922f1d2f7653ebddd27ca1faf6e4655ad10347707dd2a3325a03f29ec9940c9a8b4e9be361ea863f63aa19db0dd2e5465f5f35412c29eadd",
      "bb76ce321ebe74bfd30a94b61a580dadfaea4e2004f3838961bc9ece8795bc5eef631af941004a476d33cbf18f453da5826a225d5ec138ab36d3ea8a",
      "12fcb250cef7880db94aea722eb4c6f64c125a7aecf0bbe6ae841867b58239d85d43c1f10fe69278d679a8d2bbba132b9476e2380133ac185d8f1f62",
      "167605fb4a92b3ab9de855f7aeb1920579cea1e8a839fd29229cac48780962640c2d20e932fd92df7852c088637c929ff7800f22e0eba80bef3bbbd4",
      "efc03e65d2ffde126dafee790887e945233af99888602b54c039087a39bfbd95c29ff3a787bcfb9d7fea6fc0af43346f613554024748be312370e9f3",
      "55aab7bc165c3849980f2e018dccdc8d1009f3f5de66357a0d10d159deddbe065ba273f22a70fc702eabb2e228de1c7cb08acd32770762b3bcbd9dd6",
      "cb3bd5470e9b38fe52787bc79729914c20aae7c51a5d3d07742b5bc7695bf5eb7d9246d22cd21f6dab17cd0edef3b6a3e52c49ea30b9549b23add767",
      "90fe09b0abf5b8b3ec044361a93dd1fe8a2359bd8208cdd6a438063030c5767392d77144db9f5e931dafe044a1af83f44c3de8b782c93bc7e0cd78df",
      "54debe921753ec8f6a9a34ea439e49145a449d33acd94673b0a9394aeac506c1e51433d87dd88961db029d54833e6f51ed00694db0462aa53fd5afca",
      "7caff402ee6ad56f5552d461001071bfe339c02ddcad1fcee07805ed006c8434c17eccfe25418d260db80574b6124d7f4f5a37804d481c6a8a3723f6",
      "b421445592cd1d2adb18415283c142bfd4ab28564311d3eef28b6612b1adbd347257c705db7a2aa1306bfc4c0ddc3da3e0edd8f0e90e949bf3e10bbd",
      "02b92985e031cdabeacea95c800775aed22cfbe3b0bd3a9035e8e151ad1bb11861e2b5a27038abe6e37aa9d5255ee90334bb73b8c4ac737f51873b5d",
      "e20343999371f5a217558bf62ec9fb7b819d7649cccf3bee1f1f756dc312d67d0e09536fee146cf59f1b2c7edfebe9733fa8cf3da33ef795fadc7bea",
      "537baa4f0dac3ef7b8fad407eb73afaccfedb43e77dcfadc95eb73e7aecfddbd3e7700fbdc25ec5327b1cfddc65e7424fbdcb5ec7367b3cfddcf5e34",
      "48fba67ddad79dd53e755e7bd19dedfb066e9f7bbc7dee03f7a257dce77e729f5ace7dee4af7b973dde7ee769f3ae07dee92f7a293dee76e7b2f1af2",
      "7d6ad9f7b9abdff79dff3e77077cd541f0bb2683dfa5f03db17fbbb845df4ecd52972b97f522696e33a25aa0c53e7cd69d149054f86c8c49e9563414",
      "0d4cb6c2d6d77fb53188f118a7a6f5b2e808896141ea46970e0d36f5387e9eb93555c2b7e459d55fedb2d0c1ff159533fbafbe197392713607cfe232",
      "063010fe4f9c5fc6507b8a0a09235cd492ca75509db481a3250f6613a05001f2d4d2b423c569f5f4f6a7c7654b01aa9eaae263b5473f1f98a4ab027b",
      "2508e72191d3fd94af58d1544831ac931190558900a964643ff6f654206ba58da5209c0264a6d0b39fc3c42f206a386e76e9f86648315292eee4b289",
      "ae05cba9cc96aa206639d1fa89e0849d083c35ca3c415f505bb0a94e39d937c9636ae88e1e928e77b09f3cd2a707e68ce112774568bb439778fb439b",
      "76e2616bd6a72952ab1eb7f5e812e922a82850a0fbe9b9649698442dea5fd57e6e2b90d78a573b0f75a55248aea789f9dd2002f6a2d967e94dd7ccec",
      "5f1852d90a0e78f5c6405312cd051b0f1a15059e4754a5154a458f9a11658fd0ab4ba24c668a2df72259ee68a8e26ce947a63ab2940d62dcf3f41f50",
      "9ba2a1cdd21df58656a26287ede38e143a08d00bc79763acd108c507d83c2825add1f40ed128fdbb17f7d910a9901f66dd4b1144fbdd74def530cc89",
      "7d151d15237a1bb9bda7208367526182d6d90a43540958b382e7b8b9de3dab459f76852656a08f6d29e8e91138fb7e3c3ecad99d5a9db5e3a9524208",
      "70d86ea16ba74b169c205ae7b211a6ac3e0fd77bca950880e2eacd991d090390bf54a4fd46c7bde6dbea781d998814b7c7c6ddcf54bab029c63f47a6",
      "da2a2cc5dcb87b6b458a8cb774c176d355a7619fbead4bf4921c18f0580fce1242233c4a19149f8eb677a37a4aa6ef1a92c54adc6a2510cdf226c147",
      "b17f68e6b9cd989efd6991e85cdbbb9bef9592e1b1d656562c4d259188e1fa41a70e509efd29d180b7799f11192e9a4f1dffa3798a6ec4c13b9aefcc",
      "fc2edf7dbbd4548d1890e026dbfbeeaab91a30def960ac2422efc1a93ddd7f1f4afac8888888b4e05214761f05d4e3a20ac4c43bdaadd019f63a3a4e",
      "6defe40369e7a8debffdb9fb243887f2ebd18189e16c45079e61094a86877bd7731e5c019856a69c7451a5d5f6ce03d5af1a32044e44f508c190e1c3",
      "fc756b4122bde09d0107a7abc990c776a3eddaa8ded01dfcced38a7031022634dab5110d302da5784b8b80083b6e782c703de3fce64b3cede641a11f",
      "a5275bd0c91078b6c5a0c4e73cd9f71db8820305220928d904f443d27e655b7a6a65b67e7b83099730f91160448cb94ffb7ad9d74b8e413c03162632",
      "e8eacb761f19297172a6e5a4e892f85baa74de2ea332b42ae0dd3c26ade1ade9f08c15dc70bcda658312180f846c65432b322e1859dc4f3c175bc30b",
      "525ec5c6053ea14196f6e07be030052e8914b264422b9237aecda7605fb7c488373baa4ab2d1fd19c323e4efdf86f1dfbf8df16b54d34115f982df7c",
      "bbe79b7a3a42b912e774fa163dfafe2db4d447b75d2fcc40c222331badbb6ccee7f9c0ead0284153884cf4125b4d4ece33035f0dc0c7d5bbfa42b19d",
      "76967d2dfa940015a9371078034af42126b7961569b1076754ecdadc9e37115165c5ad1036435a7827d6afe81a3b78e29aa958d5a809b8d15f01bf01",
      "88a34896839fa99453835f47f090cb6cea51de3d6e4f70aa2d5d311d24dfe5130a415b7a569a0d6a461432db810b5cdd4735bd54bbabc457d71c6f9b",
      "fadd1cfd1800f9bff9e96f7ffafb9ffef1cfe5a564299f37eac51c430db2e44ddd70f4e756e1a5fd39f2cfc69f81fea7d2873f8319eeab5bc5a2dcc4",
      "e8dec482e1bfd78f6fdc0ba6eab73615790a5157213c0bc1b993114930f34f6252e2039df99be709f9d0839d72b9734b8e324044086abe9619e6f011",
      "e56b4dd642239a1d607c1e24418d5703770d5dd1592919f8fcd4b3ed28bef8eeb93ff61dffbb9ffeee3ffcdd4fbffef5effee53ffe19a50970878b1e",
      "2e681f8bfef0ea0ef204537f27e3f9f43f6b52c5c6bfe63f17ff3b89ca9d64afc55f677ef39c2f4eb9d993c03cfff187defd7ffec31fffaf9ffef9cf",
      "05d14f1100e2d8ad88f6a1be6d69c98676a37639a002c9ccaab446e3645c6267aaa1c7cfe87289ef4c53b8b42d7734bdbee320326d84159c39139b7e",
      "8b81c705154a5d3733a3ca894351688a9db55a02edb1ba4fc94a32a2e8acbb3ea356abaa2aee88ca441152b130c9675f110fc5469b6f0eb219cfaf4e",
      "b80fc1c6b47178b62cfd18cf3713168afa7b34e767c0e74b5068219ca7ea9295985b8084ce3bd7bcf6b46ebfcff1425f78736830c97ddc85a7b05459",
      "76f59f595ed3b5f01fb1c89d2d14a74360c103b54fd072b531349e1e8c3ca4427c3cd40746f922f5111af34c0747884e21c25a1424de2e8cd88c96c1",
      "dbbb902e07a6478d22780549facd1a44c791f5e8cd272e21e149ab73cc54514c08b19b700da8171c3cea31b0dabd971dc9f711cdc8287b72c98aee2c",
      "e0f416d9f8da22a4ae0a43d290279898447c2549bbd147cdebf656793a8e42377a94d51ed19aacfa12272d198b5e3dc994751deedeab62922bbac764",
      "d0b54769cb12142d041595da8c02560a063d8bb9bcf5c27846ebe59aec08739f971adf28764c15c578535f000fd7d816de9ab871f265a6e3794fc01d",
      "c0c0f2c29492e4371a7dc677ec418b8699fe01afa0d2701fb94b6d26e15c814c0b49c9c548df53c7a1b8bfbbbef81541c80b689c3e350aea4a6d76cd",
      "fa08544a1222952949f55f7db7c3c249df227475cd9485d8b70d05eba6373cf6f2646de158e485240817848b82a9bd15ee3ad10e171c517349a42908",
      "d64b31748f61d07c6a4a49ce282bc783d90a4604281d560a7c54c21c0ebc86e98ff402410d3fa2bf7f81782d35d35490eef87ca88619dccece4e514a",
      "6939a26dc7761e40e009ecf50c57d00e739fcc2953164a5f1c4e5fa7d5058a75aa2436f62841ffe6d13f7614fdf77ffcf9e75f39877ffc0cc6895cc4",
      "5d5ba0a79ce9d70c9a232a539da02aa161bac09ebb8b41a1107578077352872ea62c3d880efeb057d351e20fae858bbc3c0ac2e9a3d82b817bf76008",
      "bae0109da68e5b1281afe78194a23bc32edd75f577fcd1e9d67644a734e35987c2e183ffea54e290a3a93aba86227bf726625dbff3f27727c2e7e954",
      "977805a3badd51fc2dcc680e7b51010d10b48553e4663da2eb71a97701078b443b3855abb81bf5c3767da41486db47f6c057083cb8e76cf16a1de851",
      "b9eafc439099528f2c51bbe17f7db9b9b9c4935aa687c16b0fc90250e247fc8a5f8c4d205c64ddf95ad003eca7cc0216566af311ad168e2de0b5513d",
      "b4ec04d5142d5f85b682c15013b41c485ae7ccdf75780237ad6dddda6c0a1c348982b056bdf8866365bd6f23a2bd1c1c3b963bf5b01450f24d24ca7b",
      "4cf395e62268dc33a4f2ae101f460b8bd9dd4991195ca1faaa479b59162dc246d1765a5d0822ace5d10bee3e862425567c535988a30abd76131e02ba",
      "01aaacd379ab9be016d2a83b1f46a741dccb20b255ff00b4aef8c36e881acdab88877988851dd0616d41c3b795cb72b71dca3967a53bcabd7aeb46ee",
      "54e859a59a5867ffc5b8195e55104951b87c5590479cd0572ea554a0c924f1e04517624447952d24bf8a1a1eaf7308b71f8f54d514b881bee2bbb24d",
      "2eeab177e1929d8aa02cb79da452540d40ff68b7cb87aa63393dfdedf8fb3832d08d1775fb8198632928ef3f62945728a9e56405880615541cdc41bc",
      "cd55e4e955d43042cbe94c5d82356a94a224df8252653b1efc7d220674b5985e9d4b828116e212c54049e174590aa6394f88aecc980f2e69e5d57b63",
      "f02b8a7f4503cee93ab2ea424e0bacd94a90ba7e54feb3aa06a83825ebbba26c5e1d7382aa1ea0135957e8f94fdefd9ea20fdf235432d56f54a2c04e",
      "899dfb9428b9888618d4c3d5e6f514ae6dc000c0e62618e9c8db8b417c29afd6620c763afbea8b9af9a926a2babdfba0c0b6140b6ee769abd1c9cbb2",
      "459a2b2f079a748f03b98e24477e81849333777717cc0a45de0a1a8f6cdba14eb78a467bf8ce2d9e5156f6ff38ddcb40d674d5c280aae38067b42e20",
      "970dca47e13b579e02393ac0ff3ecdb53811b6d83a7b3aa9d701fb250ab7655e4ea9c559f7e370ae80e09736a2dd4abf9ebcab7b66ef323a56793e6d",
      "640b04f84aa58783755b28e02195def753913f8a5b28a53c1d16b4b9b7865ce3f6827a83c2cf72723ae9f97931934d810eed182949a3dfbd3a144ba6",
      "3f7cdb5de70c6195a7ebd4567dc45541f2d3feaaebbcdad1e3009fa4e910acdecf5e8d25c06a28a935b3cf972af4393d6bf616a35d3992de6af61123",
      "76a653d73d2a5c7ca0363ae77bcf86668d96a2d807d7ebf8056b18217a3696292f3df038244354f5f0b81c787c5da50eb150ca1c099f124d312971c1",
      "1a383e7fc4dc9a6eb5a28501103f50f9fec0fcddd7f6d78f455ffef8e577fffacbeffee2d7df7efae35ffdf6cb5ffdee5f9909f88f3f6c0312b685f2",
      "97fc7e8b42a8e984fc0a05b0e5a0711ff24643b0f4b707f1b599f5b7ce3c96dbf0eb1df0236de3d73b6a156f58bd50c83b74177c3d841f9a857ff3cb",
      "dfffe18fefffeae77ffcede71f7e75fa3f5f51de8ba7f081c44507d135a3933279ee2d2a6d5bf7d97b7b8d68501a504b027ce2053feabb94fdadd5f2",
      "99aba2b5f7f29ed9221f7960796d3a417b80fe18c2513fddee66aa44de29f646dfea6fc6fd8364b0fff4dbdffeee2ffefd1f7ff9ab1fe7cbacea572a",
      "eeff19953040a88624e504e703e0626f10d0e74c4fee6cff9d195f33989bc93e7bb7a3f57880a2c0bd0f64560cf7ded0855ecfe93932f72a3ca9191f",
      "f1be4e439d96dd11ab209511dd1a12d33e4b0d09ce3a415cbed6f2b632633ff966a89fd469c6fba5631417b57b07489d245914ff23694e69d25ad000",
      "f33d3fe8ad484cc8e62bc9e3d5ab5e87bc6bb28d6df54b6ccad68437df6a349b5a4ed52dcb6e7a77ac1db464c4308637c7a2716c10a22dcd03f6c308",
      "d15523e30659b57bebcb8b4825f2be95bc196511b28ea3950041c1ebcd1aa60729f9c475473f20af966a0acbb8a8dca43f3edd7b8c2573037c79cd45",
      "cd3b076083aa459408ef6b54ef57eaf97555f0a4e28e64c3881905dc0a2bc9b25e0f0fc0d40a641e46842d049cd1a7f632077c46ec2fcda91383001f",
      "86e9dabb7c2517c01dd7db7715678091b5bcfc5b7777801045cfa20047099cdcbcb1583b39525699f729dbcebc2867d3ef6edf3cb22f877256bf57f2",
      "41ce537ab40b2b211a2b7641901be25876df76dbe3a894e3eb9a2b26c478e786ae6e42a5c6fb89cd9f82fa7c1575e4a150bf3d51200280ef826097f8",
      "f05737ea72c45cd22aad5ca095d841e8679b37ba97b6e8c330c49edd7bd491b018c9004bb262116f7554a6d47d78c084b5cb6a01a67d57746658a4f3",
      "7455b08cab442e2e4af651fa08f8455e31a94d395c741c6fac90b95e78a48f841638b44a50566f67087657efe8b46fc660575dde75f5785a8c50dded",
      "7ed578a2ede27c9028aa5ba0fb8989b951078947b4fd5ed52d7851078d1943dd27e3fdc57bbf36c79088cf0c35408eaed69e3fd4dc5004ecd1585050",
      "7c65f64f800b54bfee1daf764beafca2ce19ed66572489da0d51b64af85adf1329932a56df3bc1c5927bc1b50cd1f0fdb0498c1efabd46a07df6e69d",
      "b8dac8357245d6fb3c510eabbe6cc94812b14bf2de88a2cd03eb8facbcfae2b9ffc173c6f11f9687a5439807ee1e3b9514812bd2679b7c7a9f1f3a4e",
      "cdf4f8c36f7ffb67748b79318b2fc6f6e20d5ebce7cbd9f834672f66f6d5f47ffa462fbee4cbeffdff6316ffa79fffef5ffeeac727f19d2a00b3b3ce",
      "1b5ab37bc615d03f8d34deb449dcc702e2bfcd8940b623c481ecd4a99fb6fc29aeb544a547a44df9063ba075820cc72194abc85917e7bf87bbe8206f",
      "ffedb2402f48f6a6036ac45520e3f60ad94c19350dc87cac657a43100e28672892ac8056d6b78bcb5439bd1155c7194b46ad869a12c41850ac6fea25",
      "907f4368c8d7f2bf4d6d1cfe3de014f41f87f10cc90e50062dda6dfb234b163464026479bb879a32bfb7bb66573d90d47f33feec1f89f9eb0cc8b8ae",
      "fa5353c5b1ab615c96a112ed23fea449fa8c1c3fff4ae39a3692b412d1a61592c27b1e0942c4a0a9206a79d70358b739502c44d008d240e73e37ef47",
      "9251f3d615379473a2c7ab0885bf4469e784a088da0d6abae67ee815258bd8dd9c29e09de20728b5aad2b33b44c2e68b12d1bfec25aa7eb613c5f9f7",
      "bff0112effa13304232a80fa280777c8dafbad628cd003b77f93eb305291e84504f4c20edc7df09e3b44640a413cd53b7af5c12d40ceb139451da223",
      "3a2cd1580616f8fda8778acf9697b2bc1f6a964c3f218a20a489caacba6a7a8dcf3b9517e6b593dc8a16e64860f4518e201ec8a7e89af9e66d847477",
      "2091a209f4b8af8986d86855019d574d6a67701f8637c73cf81bdb39fa023171e022299ae719bddae158d0dff1c720ea544f1ad8884cf976a8f88e97",
      "81322bb6f034e9a5c4ab5c755311fd684092108d2ed6ede394b1cc7011612fd940df69b7ea1e2edaeb11353db139c25f9f59d4db885cd57015e53ba0",
      "06a28fa60a64e0028793eef9a6b45c41349cfd465f9ee58ec8f96c245a435cdf9029c7edc54fedf40db1fa60438685354ac0b420a9e88422758de7aa",
      "1294e03deece0f69a0414d842f0f30fab0ab89dd78e616508babaa5281161b655d7fc516b7e2d95025aecc31b30b0b595975320acd403bb02dd19abe",
      "f8eef6c0ae600e1ea9bf4a92eb5681e3810a23586567bc333541dc49e4bf2d55de11c539567ead1f2a74e95bafd4fb55df452ba2a6b61480d589d877",
      "e4820a804fadb8e62c39aec0c986b1cf47e4728a149fd360e6dd40c18b06f37ad18f44e6d788307ef4f15ca622c3268af09d32fa3368e705ff27b206",
      "2a83523fa7bb7119ee253bbbb5e77e4b2cba9439d70f1999651cb915a728cc1f14bc80a8de79a62d81b17afaecef54c6a9e1a8f2d9e73958075860f7",
      "0b7a9ebfb04a239bce1926d9c09de3946839ed6ac5b304f81aed911df618319c1d5fb11627d05337b97833717ad43765165d1d610ad0d40e9939d52d",
      "657073bd291a1176054b1cd499fa2b8c4746a29183787d5c77207de92a7faf69805422822aa08d6f818cf3403d8e9ecb0484c4e95b7172abe3d7d465",
      "f32bab87f24071e5050e92fa473a92298bdae7c775a68ef10f6b7f6eb7002d635ec459a542495362d8063340a455d8c08bacb8af21d9209f8dc1d872",
      "8ae99aa13931c3e77cbc9726ae46da128a71727fff0f696fb26d47961d89fd4acc34c1cb75fa6648d5a25454d55a92aa1bd4d2e405038c8418019008",
      "80a9ccafd731b37d1c78be4f30792f47116fc3ddaf37a7d98d6d336a4ad5b46d1d3cfe89c0e5cb56aa407435d56f6f7dad5a0cb76d60836a3c367d6a",
      "ab6fc396d078c0a0df00b8d404c2a480664ab96c9d1ec6f87610001470302cbd4f5e986e4e4733e86980df1409c41bd1c608bb4235d48d4e40b68c49",
      "82aa7abd86f0001613f9151384c66f044e62e45cf23534d79ac3eb956fa6894630c82d87cb716d806341366c7c1b86a87610bc34be990aa757e9fb4c",
      "c270887c0cf1db5700bf0e1cb6747dd18204546a1be1aad100f12fbcf1bc071c38a23042da66a8e3402a48cf2572ccc5edb607107ce243a76f030e2c",
      "e4b4cdfd1bc02f86429b29f292da07810c7dcf6fa6c601b796eff86de8db24ecfbee509b1972999586a46d0c39db71ec60016d1a6be7a78042abdbd6",
      "41bf80383b1b44f91e903c1406fdf7cf5f7ffefabbb0c27f250e723af12a5a4d6acb5285a85b696b520b114f313474d0c381f59036b56ea3365229de",
      "4af2ceb24b8c7d28e31376750d0d4f5379a19efa55af84c21cb316611f96a938cac67bab74820ccbd229d934e26b11e9006c538c2abc3be8ba311363",
      "708c0a1cb6255946bd0a723133b5545058af661bb1cf9d9fd958d115d5c6ad86bb6db94c4b5d0985055bcc49b6a0d616daa880ca5750aedfe8c85923",
      "583327d27d8cc7d0d19f5ffff2e1972786c0f92d1f3ec6e99bb94f7b1a01a791721a518781771e9f501cb1bc829667761627253f8744575ed4d9dc2d",
      "215a76f5b7614b7c21ae3aef01da123310f007cbbc8a7ccac09b42344ad240572aeb7b15a14b62c5823fbf5f4062558a70e930bed94cd37ed7b34195",
      "c40cf5305a78b498254a49712dde35e512f80c5438d8afae52fc9919f6b4df70e461a86ad67d398a660b4c3aedcd15104f33eb3eaeb73990cd60fa7c",
      "d8330096d677763ea46bc0f227a849d7f7c06e4cff876fb80361d1a8c27c7d43ea1c1aea089fb44b91798aa21726299c239815751bc1da19e94fa420",
      "eb35091375ea11c1d8668d4948896b14940cf8449622286e514d436f97143bfc0d2dc866ab941346e9ddcad00819911283cd7a79606b3caefcc10294",
      "86c077269a2cb61701aec489bbc64ca3186fa4a974536907ca8b6ac16caa4f568b8e5c5ce0d9881c1db611b26c2118853cfeb94937d95c2bb0ac32bd",
      "4c79d86465f3124cb4b88a0348a5746a4ba263a5ee8a7b01b90e6cd576bf469db78292a1f9b880eda0f3e3056dd41a9428b9f7ac13b34157509887e2",
      "21b5a0ad2988b4a0714b176b2d4055977552d8c4df43f56d6006248d6c7d0b40d8f037c07fbb85c199dee0a95b901cdc0f9177575404a3f678a0de34",
      "f29c732b35c451f4f82bcedb520d409d500c395caaeaa9b15e5ada25b41e1058a0846acd3f9d28669575a5124de10944b02f54abdd47ada0a3c06421",
      "142dc88933d4be84204c60785ee00950f82c53b0b1013409a26a48034b27946a18140a45ac552e5e7f267721967249ab435ae70582ab264e11282a0f",
      "5332c01f3bd732c57d93f949cde2a117b449b47d5449815aa5f19bbe41a52a7db0660a7c30b045bdc0f7b1869520c87cb83ae0c8bebabe4db8026142",
      "ce685971d9c686a4f57ede850b6f889b04e11266d3dc820bcc5a403e62d89064932c4910aa592a6a79cc766da448c211d75fa141d5d5285830d007e8",
      "11d14b22ec156952a9eb8a0f67d3a3228d43bdf36b628187094cadca053736c1bf8b246a30037d4c6c88a6595190a0c5315670229607e466c87b282f",
      "d8a875b8bcdfad51f542bcd73bea23d99f05fade1843360751db1c904decd31e1299ea46f560f351717788a7c19d628d2d99a0152a0c876d01ce15d7",
      "e9f260b81815a82bee95081b0858a4aeb8175c1aa814213dd376c74c23ad03c6fb3e86356b0017c35e8113d95858c5b2be9b09fa0eac18fb24000ad9",
      "7c6407f44601cab237f0899410941d2d6f88de1dcc5ab4dcb7ebaa8097a077bbeefe1e368b20fab3fd97a9c177f10aa0c127062e39a8c7f5bdd9accb",
      "128a152faf10599cb560985f834022532d38f4bd9131a903a4e176f3327a6541eca24f403a3a690e1b22b492321b3f352cfe26255de531fb6690265e",
      "d30b68e6642705891b272b55c8cb82a546c5dd701ec07289a25fdc170e70bbd21685a70f97a8da9a63bc1072d864519aaf1bf2d43b5f5ffe0638ac55",
      "6719f48cedd93c26d88c02b810de76fae66825864ff87589ca12cb888c54c61e1e3706121df7a0e737c8203f0965b2838d0112e6601034cb30807b25",
      "204389d275dbf04ac404c06c6b58c002884e2476dace6a98ac28885d402e9687b19b18721378a3404dda641857098170854cea947f91ba84149fedf7",
      "99119496adc6365a30917ae15a1b379eb4c8bb43ebf6d890b011696ac24a106e0f86f2172245a3c15c970f0e07ed0fea22025a76763abb502a8886a9",
      "ad4098bea04b5a5339a313187e0b52165a7ba0bf08a99a97ba59de81d3afc8ffc1248d3098d6e3a3658ef01b83fa1204d348d7659615fd7317548044",
      "22be41bcd314ae04a696b9edaecbe40b0d4c47075fe6bacf36b055511031ee7e820284209d24d132b2c7a06f5bbda0c46b7536df6c6a2c827b642479",
      "7f4d12617863b1d1fb03ee4f5f3b53b939d3667a05b09106010e6b9208f80be9e9831cdf2a62208091c1284007799a9012f455e876510d4f39e3c2b2",
      "9e1ceedc0df6074ce090b33e2dc352288e213f7fb60b80dc4a572891b5b6a12757980e42213664b074da3097a78d50c1a6c0ba6e67a2745105471a76",
      "fdc17c26a2a3660b0479ec4a378c8ae116014d8b2af947eb1c05508f883414fbba0df929a127168a6d5be5a43394c11c654fe80a913f81ccfbb53440",
      "0705016752f7146c9dd0a84ce7b81d13078f05d29f7ef9f02f4fa89e1ffdd8a3bf7bf68bbdff7cf2b30feeb8f7d97fc7b33f45002e5238041487b0e3",
      "149d9ca29853b4738a8a0ec1d331c43a0462a778ed10d6b9e8ef1425fe5e34e9a2ce53747a8e620fc1ee2d1cbec7cba7a89a34913a6dec240d945964",
      "da013e9a100ad748b041d9a7a24393696b39ef4f95884d23d7e8659a781ba8dc0859478c3405b0c9e538f787e7ca4c75cfed2e92ef850b9afa31398a",
      "bad6a96edd3fa07dd0729336f73246735242226eb67e8e7ad060f25c0395c3bde7a641e9cfb481ed29120ac84c53deb38d605a8833f6b4676029801a",
      "e28d6fc1360acf0e7d971177941121a0a6353dec38a342e2989f79c7b913651b0c901a2e65b70c3661ed24e19b245cd020cc69073f39f3d3234d99b7",
      "98dc448e1bb66281296c696ace5c4111c9cb22a7e09873db5aabb2d576c5736f97a08756beffebfde7af4f2d7bf718f210699ee2511fb41e22db43fc",
      "7b88920f91f421daf621f9216e3f44f7871c80cb14f874c229e9704a4e1c7218e754c73d23f23b1b8e1f15e7d1731a657e349e46ed7970bb39709e2b",
      "f71935d912c5791766dad99f8c86d9172a77ee0e958cca3ebb1aad12cca63a9972d36244314f3ae294f0d00b86a38ac46a64215abf40e18c645b89bc",
      "4ef067252eb210f4d6a39265cb5e93a1007a97222fbb34ade0d479efdabeb3792b2098c9c2fa4b3f10168a4fbd20bfaa08add36b6c30455bc4d0c992",
      "3b876eb5ce56321ec260c47d3866a2d7e8052573c5e803c023a4c5e171db51c84a607c94cbf740bf68607d029c3cd34c840ae2a1e556a3e519a4f538",
      "4f0a6434711515ecc84c81d96f908b185508c45b2a0faad6fc33515ce6a5d658d87c222b58e679d57ac7d8e48ce451be22954edc0d7cc23f58ae00cd",
      "a58c37401165a42dc8c730eaaa5643e6b7ce4c434d71acc294596a08466409cbfa820cc48265d021d151224d757f9a826c17062690347bdeb73a15e6",
      "49eb86ab83b0e67fb0be4c0ac010270ff844da2b5da7cfffbd066a89c4ef63a71e5bf034b2740037770b7d86c28d14f0dbbdb6b699f8aeba6113a80a",
      "a4a0ce6036f735f931f1f14fbf7cfaf5c7a73c609f437369b67b16ce67e9ee993c97ec3be5037dced0e7156fb9479f9f3c24315da2f39c0c3da44c0f",
      "9be2ed2b9d3fe5e183bb6171183c87217618887eb41e86f461e01fa6c761129d661a413350cab475bc49f30e5491cd0c155d85e1dbd341ba10e7ac7d",
      "679318b5c6cf3bd245a9045d0858c2dc744e44d08315396ccb9c4cdfb5699f179261c090e10b5c16744850c4695b48fdc766f066bf35a10884ec4bdb",
      "16e4d3de8126ebda12c0e1410dbc7d4881186421c5e7e677028417f437e3e2625a6fab6e41281e82ac4e25dbc1be0c3bacf10cd7ba886450bf7ae141",
      "a403ec1a90ac17a55321ad2abaf8f759a9a0bf1395f27c5d19449ac8f55d0454b80e62751b086214a4c0fdb8b64d34fc507c610f03caa2620ae4fdcd",
      "d960344d4802868266612061549a05b11b58fe2291ca66596f7dc012c35e9123b84923727dd75a5ba059487aa0ad9bdd880a444e6f7b5d09242ea08d",
      "b14a2825b890f74b5b2287de041a2f63faae953562d7a506d3750ce0d4316ee5362a5e23b70488ccde82884b6cbaf1cbd4226a0218af571107bb2058",
      "8c6ad84e6cc079f3daadbb30e4f39a3bc0ee0ebc4709c6f3ba0dc418ac449ba5602b05d4edb2e4068e644cdded4017488bf7adfc480bb9942e98305c",
      "553288aee5245cdaed40ca01e87f558f1a9a62a0a9d1c7378f34e39816bef302eb3bf69cef286d82b602e0cbab5298e0a0421226ef88af4d5a906833",
      "4ba6b020aa2a5739079019e82ab61d1376f4525634adc42b4ac4945baf781fd2d0d0c21afbbe4c86423290d257b899038b1dd900b0b84c00b61bf9d0",
      "6e968234311c8f7005a9e062c1926c2592b75bdf6371d7ebfa9f27b6dbc38c3dcdeac3cc77ab835f41eeab0cc26d348d438f315f2be5728dc1a2b659",
      "f128c900faa4b92d25578cd86cd9c4be1960507dd93e634c2b98603e785bc81130b62aeae5252bcbbe9d56848503ed79db43cd28a200023bb6af5ba0",
      "1449b8f63e26c027476b41dfd749ec7247cfd976ad0b34a5fbee39846560d1c3cbd06b1e9469c50c1aa29b8205ca2f78a9e5b2542ab7c1dfd81e7a87",
      "48285a6c74085967f0fd6237983aa44d387dd34cdb32a0e33b087ad9960a43b0306820445efb3cdab7ea3620eb36ae1219af824673b048e46d412f30",
      "c40ffab6f4c8ef39b7879fc08f9fb002f5b28f9960a40764beed1faf08c2e60597c3538d355b265be7cd024d6a7ee2ba5fc520c10d0456af9703d958",
      "06d93b9e812e9512a83be8291d9661e41f089590d444262ceeb1431a58209cd38eb02a724ce8c38b3b0c0bf8c21866e675419e6e605b337150b243a2",
      "0743cede0ef1a09081fea6bedd928eeadfb8a0d3b709f9a0d7fddb97d71ffecb33ddd5e78fe03f94ff98fe83fb41e1078e1f5cb7f147b264bc061383",
      "84a52359c9eef8b18f8183845e7e6d1ea08002748e5281db32b112038d2da773801827f3430d1b5b2c7f347e5e6b9900b7068f9952bf8565f6ca0d79",
      "0dbabe4d2b3ec4f62b878d721dfce0d5de4555e2833c46d52c0d151ea4de6dcd01c6aad3b3d97c1b3445fe56319818ce83246924cbfa36b1bd5ac9f8",
      "eb5a71b9bdf4896c01332943f83b86e0c67d075eca6895f9023ac004382fed8743033d8f9a5a8907bc29f43c84ab00364029bc16924872b2fd51001f",
      "e389d66431b26a0484bad77d6200273a34f8f2610d788c10730df8cfaf3f7f7d66d41fefd63dd3e1c98fefe7f016efaffaf03d0e5f0dbd2483a612f2",
      "36559035442c7b3623d0b8dde86e96dd7b149b1cd062158ac1b887a619aed6a60a657ad45dac04826ea7b559f0a8b82f95c16504cb28bbbb2a04800e",
      "e0762af186c629f45003616fd99aa106739a9af2303021d24261c600f87c1af43f090360261214d014afd780091e0b250dec650d76cbc0216bfbfdad",
      "6895d79a9636c45b86ae1e4eec7b4ec55669c9d6180153196a07525b3e4c01984a986adc4bc14416187d0ac6683d50e6666351b489cf04607b475ed2",
      "6bc81436c54838c94c547bc495ca35d6b4aa14f6c7eef99306af35863dcd9a3f93bd4d065f051d5ee9bc85b51bee2ba56ced4ebb678e9a6a34814660",
      "9f08f70a6f21ca311d24c0e768d82562b1e6f11356ab10c304aa7fd4e552ff76571a5946b04a13f6498cbf38e3f5388d815336680017f2c97067cc7c",
      "1dc5ae3548a9ee830648cd11baf5bcaf3e88c64098944f6bcd43abc67ffcf4f1a7af9f5f7f7b6ad138bc7ef7890e9fd17decc390380e9cc3f07243d0",
      "0dd3c3603e0ef9c3c4f0b3e730c5dc443c4cd7e3a43e4cfdc302715846dc627358928e0b571940429157aa9a69ade4fcc5bcc301b862991361ddd6d5",
      "b8c9ee406c093608079291596fa2ecc64d48d0b3d7a4d6ab353572b64ceb1b62f766e6cc8b25ecbe4c2df29b597ba00a9bf9fed0c374f56ed6aa393c",
      "cc02764c7efbed8421de474e03dd8f96a307bf3e5e32d468dabe03c86ba33dd1b0f76350f19ced8963ec5745960eb4a1b4ebbd346cd530f5aba99682",
      "2e91becbde6732c4a27154b9b6a7c90eccba252b602a905f82c91256c3c47570e2bcb6baa0deca10ea1e910dddb8e89aac79fb56a5b0b5728ccbdb8a",
      "e0038aa4fad85bf2841e28dbf295f61b40a6e9692ce5388879e4c587a15f06c884f1b5ead655e7aaa556ceebdb6431f147f62bced3c2fcd082f4b7bf",
      "fcf05f5f7ff997d79f3e7d7e624d3a8e9ddbf83a8cc1e348f5c3d98df9c3cc38ce9fc32c9ba1d014cddb458b73e6d563fbee3e9326ba4a8d18f52059",
      "8df9626d1c2472e48daef5a4edb6e659b94274cb54b0d3b9f0c43176d37447c06792ec66aa593b652a654f33628ca4cfd84ed3f331a6dcafeb5bfdfa",
      "facb133eeaa96df9d8dc7c6a81be374afb6eea43cbf5b131fbd0beed9abc0fade0e78671df56ee9acf7d8bfaa1915dedee0326d461ae057246ad0c6d",
      "0f315237a9b96eafe38d05a08adcca1e7615554c2e0db5fcee9ae916d6c3f27b5ea4dd4a7e98b2875177189bf711ec87f9792e34644c2263e2612f3e",
      "654da2295a2d981ad250b8ba75d0de47e06358b8f7bf7cf8cbfba756b4f3d7be931b380a840351c2914ee10ddd42a2e732f849f6e618504d207ce0db",
      "7e99d90e1f0c548513016267bfbd7d43085b68e02eaf635f9e15733468e68bbd226535ef87ed3e91b38e3dfebb4086215935c77697ef32757eb14b62",
      "eb38bc0f93e03855debce287beebff78fff1fd5fbebe7f66453b82f83dd0ffde0ce01b065c4f816f3bf0ad09be7de1d4e2e0db207cabc49b5e0adf6c",
      "e11b327cd3c6a9b1e3d6fbe1db437c0b896f333995d1eff586534dc2d72d7c6dc3d73f7c8dc495517ca5c557637cc5e65ed5f1951f5f1d7215245f65",
      "3a54a27cb5ca57b44e552f5f193b54cf5c85cd57e16e85ba7b25cfd7fa0e05c153d5d095160f05c8639932d350af5a1ebb1cc0bf91f7112d0565c04c",
      "e8aa93699aa9b352da7e1da133273f366c0ac5f4ca1cad716a502a8b09f768080e9202f39c7d2fa42641a74fbcd00b2134fb71fba291c57b82bc763f",
      "5902dc0d8986728d8d49d732b1d977a333a0928047abb35d6037594abef0d639f365ef70038ce5d0a20340dadc18005366639eb0966ff05474e492d2",
      "c3589191cbad3c6a2d0c7b2d19eb6332bbdcaf0527b32f0835a9792162b119901d64f77c856e69e1103756171ab534d932df8082d5fb2bc9263a9a3e",
      "06c6d7b4f84622c5c84c6ff81488c6b1ad603ac56b95eeac3e550b16088a6669a919cca3926919c748489050ea0a691e907a6dc0357846502fdb5d5e",
      "097d49e46f9f17fc1a2b81c963e2ee22f2dfed62fc820558186c8365f35b2f7735bfab174cecb6cf3ce8b0fff9f5e3137bdba1edcab766f9f62dd7e1",
      "e59bc04e8d626f5bc97cb3d9a11fcdb7acf9b636dffae6bae36e0d7487163bd783776bd33b7801fe0b9ebef27d24b8c1721a4ff711771b936ed8562a",
      "4ae3223bbf5f19a8723bb85af9a0ab825f3650af547eb03d0d83dc429414a0476cd0d68f0fde11880891bd600f071026a2fde0e20c9f509784637275",
      "fc0d6c4f78c1edfa2d6ca94855c8bdc02002397fbed802c13d3e22514639ef81c6099b28d074e8137c8cd0f9ebe70ffcdfc7a7896f10f44d84a746c3",
      "7b33e2b161f16d4fa36f7bf473d4bf96d3abbbbfded327b87fa6d3a7f49ffb3e2428708d99d0ad57010202056aa8c8806df69490e9d504f3222a67e8",
      "3be6dbea669700931e7bea761b66a6a02424a4db255db01c6d24ac376f058278940ae6c53dd1d872d90d4c71fb7a0f92c7bf7efcfbf7cf2a301e1a39",
      "7db3e7a921f4de34ea1b4b4fcda7f70655dfc47a6a74f56fcbbfd1fb3bf75fe5fee5dcc705b731b4a7b1c6577bcc81f20e7ab2e4f9425e109d99c8f7",
      "5f7a16dfbf40095035aa6e417ab74897af814831c1a04e04e05f260f0852d98290df5acdcb3081665030af6819fd804ac6018cbed6bd821843fad1f0",
      "5c1b0cb5443334c91027217b9a7846a14ca064be98b0297e4af79f3d00117faf9d50678cad052a6c23165ee89222ab33741b01a4506867cc5947f44a",
      "0d5231aca380004425ba4c8a3d186a2d05d90f163f00755ade1ad496a7bd9b8aac5aa1bc1c55d2514fc5df598d11900781426ba6f218086ad1ddf40e",
      "f230a2e207b1d7c0df59ee25bd371ed0055785e4390fa89b2b4aaae6b8a4462b5b92f91b6b77114d1472ac09376ef4e14017c3fd8552ace55ab14a56",
      "bef2b0f58da09c9bf955ba366150de447e278bc2d05a83600e992f659b5047815233847d8def4c6c676cf42cfb06a9283daf9002c8ea49cb5a3e4929",
      "8da7a4bef37ad564b506aa9f7acfe2d7c4320f0721ee503d89df1392d79a15f806ebf6619034a4c9c5a2e5baebdba00a47551c3a26288f4d0eec1535",
      "16c9f1a1a1bbb00f8632dd901ca30c7551fbdb6d7268c60c96de0bc53330ac2098c0cfa33f2aff2daaad8fb99875457c21cea58e3260a10e7a939c1f",
      "b4c30bc1d9d2f7fbfee20f5611fefeeb73158423f2fad00871ef9638b454f8be8b4373c6a185e3d0e8e1ba410e2d23e7c69243fb89ef51b937b20025",
      "d8c9808f55b26e38d860a7e3e562e0c4a6d6f966187f00cb1a5ba8d736d476534c6f6ca105b6691fc4170aecdf051eac9029c23b0e6ddf7b263336d4",
      "4df2b73e99e5d9e377bfc110210e132ec22d80d4402f8f69967e0fcfee31ef275cfc0d3aefd1f50e80ef31fa6f47d36358d6afcbf40938b64f4ff58f",
      "24560191cdb7660f7851318bd4b8eda45a1519b2ad9e50ca64e58e8ca83b74af732adb3ce2b872684547190570233896bfb736af7ab59f1ab96fbaa0",
      "bd6feeeaa1f7f17fbcfefafa14a60fed01d0138c4c7f5f78fcd964b2ac24349c880040d370d89f3ca2eb8b29f4eb9b4fcb610fa35c042e383426dad7",
      "fe75c1809b1d14370a32a38f97afba5d3060598a11915397ba17e398de177f73ef8f0121bffef81c0432414a90bc98e3ea884386192693d182bc273a",
      "2c600a1b0cda20c0826ee366f93c705e479af2bc80ccd8f760aa31ed679cc84b529ca46f532c8d14d66b01da6f1eada5e4beb61427bfcfda2260b211",
      "8695790eb2e16e3240f65580151d81b62ab2e89ac8994cd7cdfadad1242156de340cb88a86a346feed604c28380f4ca8e8304f796f050d8a665408ce",
      "fbb632d46171d486c502178eb2d7249ec74c01e8637003c74b53bca08b29511f70ef0f111ce3bcf9ab3fa740ea048f98c61e4b3d8a44bdce0b0c9ee0",
      "b4a3916aeeb53981c6181de2e9fa88d4f9640bbae19922c94b608ac6ec0e18327983f3b88eb10f6d991ce05f7be5496d5c4001b016f23a6dd40b848a",
      "210397a06cbc4e00a00a0def368c067742dea5519f12745a8c22be5de820a4ec5154d860b404854e0cab76819156dc422ed73ec6c62aa4a6e79de302",
      "5c02bd00b95c1329bd8fff87a6dcfffcf0eb8faf3ffee9f154c0497ce720d1e3857cbcdacf491268a08dbcb3c5df6861212a91f883d1e49249fb6a3d",
      "ff6b648e6d32da95222a96c4ec6e153380245541f7a82d1ef8e26c9abfa306f52277e132a1fb53a66801d4124d0da6629dc82a6f81dd962ddc583ba5",
      "4643ae6f9146acc9674a3701a8179ab2dd27a73d99258470261d347b951b78a14d676450e90b2dcd639ab4c67a2fea720ea2f2050f6d6cd6805d1839",
      "8333b31b11c672e04ca4276775420f31d4921c3a8b40632811425ada6657ef267934aa983272947a122081c34c9a71d0242a5dedde5d6b034d53ef65",
      "f66947b5caae7880e4f7892d1a41c375f54e2e2348d42a11cef3a6bee99a32d9148fc23093ca6b2286e1e888b1d5933cd4835c305f7efbd333a9e4f3",
      "6876325807b1ac83a4d64178eb20cf7514f13a487d3941b0836cd8415cec204176102a3bca99dd35cf0ec26807f9b4a3c8da418aeda6d77610759349",
      "c3628dc7b47f103c0f3045bdf9c0d55c042082bb7023379a27616213b9706ca63793e5ae5de339a9d4c217da35d1878abc20836563314d5b5f7cf9db",
      "3aaa891be6b4b29c961fb7461d56b2c37a775c15dddae917d8c32a7c1fdd8f41fe5f7ffdf04c7bf9ef7d5927eae7a4ff0e02810719c1a3d8a09324f4",
      "ba855edb9084c79275b3781530ae88ae71ec0252db02bcd134e3967b5a64eac5c4e6b2200591d00551790d2e7180178dbab97ef57b488019dd4e1479",
      "3cb044129b43c576ca92a9c8896fdd4db552c92cd3fd9664a460462f9440d10fc26fb2153b489781edb636507388129704b976a569dd43b7cbc72232",
      "a529ee1672bf27f1304d91b801700ca4234e2c724fa94830355bb212e624e5c96261db86598dfca38b0600a230bd27db60e9972550c5e8da552cf8a0",
      "3fd316d8990833993831c1a97dcb2667d70d54e1664fd3f530a90f53dfaf0f7e0df1ebcc6131f24bd661617393e4b19ad0fb8fef7f7e7d82455c7b3e",
      "316cd5b29ee8896f6a57e8799b9010a2a6d6b78386b07633ce2d74df83407acd2a47f485004d0001cab644e2d0c154680cfa95fc7430656b6958a675",
      "bca4580c5b4511804ec1966ace340620f2f738518e0599c9ecc49085f68f545c147630287e20116120d86e6dcb639b129a09a0d9678a0d11e83ca15f",
      "b7fa09c29b5c84a3dd4209f2fba097b005bed0409a692a765b1826832d1a9b5226124cca3ea2d2ea6542f501476d1991c4159fa66cf78099d1d9015f",
      "2d6712490ac9a6906a2d979168259a5a1adb3462e3e5ab6114f035287993ae742080a55ddd4c7998aa0a12fe9118e368885640bbb3fa10620efb9b41",
      "909cb21f611f548955cd7caefd610df95c8d8c1f832aaac965136e6019e97af369abbd1018adf3f66941d0cb6e9c9c1caf422b374b48e203ab81a2e5",
      "eb9e0496454ed4541db0660fdee668d7510d128ebc54deafeafb79f298d8ddeb2f1f9e91ba3b0ca3d358f323d28ddbc3e83ecc0170b909f4dae716c3",
      "3365a34a325533e52290e325e382c831108b1bf3969d43c7264dd9d41926062901ffa9993019f000444a2f8745b534b0ca70ae56ab71acf321fd2011",
      "2855d7d6d143b33e29bf0297551997560cda03fcf3a05e9cc113d0a983a40f837b2bd454169f51ffe14db36a977899a8980e0a37f83b4ae02f89c894",
      "59154e33bcf4c6cc4b137f1e66c8a05859909e63213f090cc049c9d0d56631b43b18b919a6a289d90235c1f754cd5546f12cb2bc5855c0199586acb1",
      "89dc53254e782a4788324c215279aa5686a476538fd0347a01aaaeb32764c55bfab00452f318cb2c25a07bd5c96485acb5c70d81a2db55c8ea360b8d",
      "c5156d629a26c6fa864e81c4e98594a4be4580140869249355bbd250fbc33ab458410c0e2e57262bce22eb46045d37e9f90accef30f0998d84d20da0",
      "66f475786d455db0c6cfca659050adb469809071130ecc986d31da0416ab06126fc82f11cc6583be29ca61bbf01e74c32e63f7d7a5b50b83a96c7460",
      "1f787f695a1bb421a9014bb37a05604afaa969d0734b63a30067f874e8fea9c53558a105256041e5364901c88789a78bd930f3582005d633bac7612a",
      "83f1c2490cc337a6ab2d16340fb621585b36140fbbf629c307a1f6a6856945c542af831090973172a2c8802f0917b8c5c950cf16e58b0958c1c50ddc",
      "dc4c43d9f857f96b295e2b5c9aba81d0f7d255e054d272ad66df2f9f0f2ed85f3f7ff8f2faf1192a2d6aeae6623a8ae3da40ab891f0ac789463940c5",
      "61b2d64e2ce151f9ca6e34c9f05d80e34c6c2add8f45890f84f57b2143a96b522e73666b3ac2be273dcda11a01ebc8d450da1526bf22fa55d3afaca7",
      "d5d7afd0f755dcaff47e37386c19a78de5b0fd1c3729b7951d9d3ee71a1e1cc8839b7970460f2eebc1b13db8bf271fb942949cae58de2e58e2c2c65e",
      "f3cb7553f34434f43229721ba53f01bf3053a62705f6c63cae4ba9f4108ca88dbd0f8d597f4c63f3efdf8ce2079b4f3e7ef8f878e64d728f6dcb7403",
      "b95981e1b3dd17680783fde92f34506171929435de22e9ad8afd09f9a1890cac500f010f4d0137fabe20f8163b82b556856b699ed300236bab642f74",
      "148d1580c9e2b5aac222622a12aedd848162a181674c8990c27b0f5c9592703048abcac34f59d951f41df09a5ba41ed24022cf121286dd1022461000",
      "0e0851e18a73de1e44d086112466855d94486e13bd25b17c16d25b088f8ed2165a0e955b446b34788586b6c8b75fe24156889fdf3f0e8340c7568d2c",
      "c949dd880cd08d3598acd18df9cea50820756e4c5c4cba15da8259d8fb8ad499c999a385b55a6dac5b3e27032744657589b35712c370baa950c5aefc",
      "642dba0cc4e8ce04f56b8b723e52d749cdd99a6cb42846b16404968959cdf5096629b4800a68ee940921c7e8aab1b40a3b46b06f0988c09fa845a1a2",
      "d4ddd918222fcb726019590e4527c46770bb560fb8313462dcb29f0cec777c7062839a1a45f5e0c444d10d139c02f385c46d08e17808103b65ea01a2",
      "6037c039bd534aaa13ba352c32cbd2fe05ee5aade7c1306923c9d98bb10a1d564a953f3899df4297b13a58d7122040919a85c1b76888a248e16390bb",
      "468df909b27fecab952f20a36550d1aeaaa26eaab9c97898ae6e4247704a32aad692109208ea2c0b8bd523eb8357fd5de8dea1d377086496d82058f6",
      "e474ceba73e79dc37f0b085cc0e0420a1774b8b0041584ced0a6a90c8a9c685289572d247809d806a1412932752a81461533d5d9cc1ee7ccc26c51ed",
      "e436b31f5f4b9e70c1ee1b895fdeeecbdfdbe5f1be7c1ed657bf02dfd768b78abb75de0dcfdb007643dc4d02374ddc447253ed3e17dd6475d3d94d78",
      "bf269c968dfbcae2161fbf3eb925ecb4cadd5742bf5ad237532820b021b5478d5932d08b8605cd16b82bab34e08ba9a54762aff871d42819b6043b24",
      "361156491f929faac2072982e7322f933b39f2b681650874a2563304a02d81bce1fba7a33ad9df306d96f5b2667f03faa9c8b04e85c1309e8c0b3b1a",
      "e482faa5f93dd82211d4a88e2f9630fe4bb54f4a929ac6363a011dc14454d944a7a90e571c92c0c940c5094d25457216008e000f4cf157a194d164c2",
      "2a6733e82650cb1054d1a85ce126a0ce63682d0fea3769826aa3c36a2d807468f75429e9dd25eff0ef703efec3facfe74f9f3e3ee37ff86d19c4188a",
      "3aad9ec34e4526dcaa0dbd154035518145f35a989f9bd4b0b543c224a3994aede4972197db142f0fc66f8fca3da4582e37829c5cd2ebc1802e644c9c",
      "9ae52c406190d93481300f7b1ff3486628b9912752ab00b932a63a28f7c35026baa146b57f034176337e56f2c84ff44614ed6b4c02ad918caf19ece9",
      "4b1f6846325565d2cc7358c4ed8c0514ea50fe8adb5f0b55bd1dd55ed7586fe31df0bbb96effadb38540d9b3483de229965ff3e726a94c36350013ca",
      "ecb855e68682d5e4a54ee6dd110f3e7489b8672d3b72b3987e986e47f80c0492cdc95e01415a473433a4ced686d4cce7616d17c29b8ccf404657f1a8",
      "a1d3f1836045c193f66c336e9635653259380c6a5cd0659db51d03623530658a68166f33f23e61ef13fa3ee1dd8ae0d70cb7aa1cd69db70b935bb90e",
      "6b9b5bfd6eeba35b414f8bec6d1df64bf56939f74bbedf169c9ffd66ce3fb4d2fcb74f3f3f8e927511a68b415d94eac3d85b947b0f82df86c8f7003a",
      "2ce70ebd4c6958c49b1b0e4ffa8b22db2011296a8f90d3545413a6ba111962b9c803e9cdbf9b862830d2e4d854530cfeec7cd4689a66a0479f433db3",
      "f411b0e0a94d3844f38cdfbc9ec7fa46fff80cd6c7ddc2e92e6f0fe21ef4fe1eeeafe9ed4bbcbfe3efbfc0fdfb20478f3699ba7d33a6e8b997866009",
      "322c7b6885ebc6b48dc14c5285691ace28e507348066a3dc8793c6c6cd54f72913558d8ede8d2e436d090f918ce68cb0a7ca699bead5bf328d5a329b",
      "a588297c9ad21470b399f33618ed0950f3d178e39a95254ac6f407cbc2c87e503cb6fbff3f5f43f887f2e5fd0f3ffd2f7ff72f9f3e7c7e86fec465e3",
      "7dc6de67f55de2ff541b70e5035f61b857217ca5c25733b0eb8a7c209a5a3a760d72af1a1b0176c490c8b3aa7b21eb941afdd3aed04492f14eeb1c39",
      "7c7037246e83e630ac6e03cf0d4db8f7ec5e9e8a525917228b7eec1a31a4e083d3119b497f91d395d730f5abf5fab87376538802cdf5e0d69a8da05e",
      "e565044f06b456eb0c78cb0c518fbe54de693492054c53b48af6b03f397e0f2b69ba2a2940de0110beab2f6be06007969f42ed985a39d3a3d5700a68",
      "54987eb74b14edafcd40ef288191d36218431a583be0e745e37e005d0ef7dbcdd845114d305c6c8c37167ed0022c6fa15b71a94f6c1fb95ff70d2df8",
      "09e44ddd3f222a836814906f67c083ddfb1f3ebe7fa66472c02f1c500e072cc41131e1701507f4c5b1b6ee2b48becce46a51878ad5a9ac75af7db9fa",
      "d8a984e6cb6cf7529c2fd79d4a7abeece74b83be7ce84b8cbe0ce94b95be9ce946981b83b7517a18c7b7917e980a6eaedc67d37dba75d0d96375306e",
      "280833c160fa2a688b1b602fc8b6a8c3a9cfb8424dcd0c19bd1248248f5dd342e717d9fdae9216203768811abba29532a93b8a89d960b40451ed97b1",
      "c70f74f1e848d76f75b082e79bd6b1020ba45ea871bacf5a0e3bd7f96f45b088760bdc8511a36048467dd5f00d0752037cca51c78621954417672887",
      "494532f4f3a22d30ecc2642711fb162b54ee87416ad8482c74e9c151da8813109b72dbd9739c3b13d366a0a5be0a6bb38b74a66ca04fe36cc23b68bb",
      "fcb61e5591b771ecb2fc1638a483f10611fc53e95d54635923045379b96a839a846b917e4aeb57217692760fa6b11f8f004fceaab81f7024ad2ad110",
      "0a388a092b80d97b3a41b39e584f5ffed70fbffdf6faf52968df01cae5005f0758d8013c765ca20f4f7b7a27fecdb9f77bf80ac76fe5bee8e1bb133c",
      "9f442d745546c3d0fabbbcad8de7ca2649336ca562815754c72b30dcf5d316b3adee619e00720f7dd1fffce1c7a752da87e5cbb913cee1f02ec91b9f",
      "c5393507b7c73946ce7572ce1595bac9d32e62c3e5d31566038d97013a0f19d7a8b619613f4127248666d835fe92999de8e35ab722647fd0e359f756",
      "4cec15de40db9033525c01433c7f7fb1bd2dc86ecd3eacea6eddbfed0cf7ade3b0b77cf7f51ec3f57e78fff9f3eb0ffff9fda78fef9f72ac6eabb25b",
      "b74f4bbb5ffeef5b84df460e7b8ddf904e9b96dfd8fce6e73748ff5defdfde8f8fd318f2e36c3973c8990c63dce6eb18f898c948d291d62a4a848dbd",
      "c48112178dbc79facdf031e8c4d7cfffb8bee00fffdbeb6f4fb4619fc06e1e10e7417377609d07df9d007a77109f07fa7930a0070c9e40856f61871e",
      "98e8c18b0ee0e8caaa87c2ab2bcdbeadddbae2eebdfaebcac37740c8013172c3941c72823e27754b921c12293ed9e21332a7a4cd2dcc778900972a70",
      "c9049f6db8e7236e83f2b19ccefb8f5f3ebffef2c3dffcc3e70f7ffffaf187fff2fe9ffef0783d18598a828d6d0507f36ae844a6cdd866d0a906a697",
      "4aa65bebde4232a652b7de3abcbae52682359a44651ea66854d2dae2422575d0febb445e62bda9612d2dd0c55b4e4c98bb0d0fb87a66efa20c3d70f7",
      "6c863f4854be268feabcfada26ab2ce9ea7d4b1427cce2544047f90c4809a28ca6be92028e0c9463d8951129cfc5ea82018a0765d7d1bba35a0fe2ac",
      "8c72450dc51aac48509ab77e01297541fed4e51e1240d67146d4103d166ede96760ec51f571eba1590ee15a64309ca17a95c19eb56e7ba15c25ca9ac",
      "493418bdf2aa1f8186432c76d20422968f4774491eb3110607543128b393983bc6901218582b501ebf1a8bf04922f91cf33e60308d6ca42244f41529",
      "40d9470ee8a945b19ecf015494a95e09e11585718844f4d9504a108fc382d0f9648464c84f5eceba757cceae4dbc2bfa4aa4e7d246ffad2baa76edea",
      "d15a4e3324c327396265805632dabc8a9d310255a65a6bbba712e813660692cd2549d40d614a39b94ae42b0c66e893df6df669ddb6096e24a092cdfa",
      "b60b12f3281072bf4723a5399ad55a26bf5f051e94f8faf844bdc7cfa3810c0c6a916b7d354380b303f65a816c304606854e45d745d1ec75d768f3bc",
      "8e582f1abc30437021b4b5659373e5961a291048cbfaae7615d0fbc1b282311ba303ea02d83944fb85611c3ab5707b09068c5beb15d975ba84cb290e",
      "01f56a6a6a5badb95768dc5239c60076119729bb6a8139070001b686be4b6fc904afd5f44adc100987aa013f1ab428f0b7be3f4509062d43732a0109",
      "c253ba10cedc958374eab3951f33b4859145969b48b7958af4215b5d57caecc0fb8d0dc700550fc889e4edc34714c350573a1187b04b159c32d67959",
      "26fb3ad1b06bb5d2daf5728d6a142e62609b34bcc90bd4c77ed4abf8070210760caf7f19061f581b33899f72de8842a4bc246a9cac843a4a062bd0ae",
      "f552901515f128925ee45da0791f37f29b6d55d472cef604a48c0211eab0f5c715cd5d59fd5e77778579b703b83de2be89b85de6edfc792c75f2fae3",
      "3338907ba1df43011c58c0c109dcbb3bbcdddbfb775fc87dc3c357be0f84d360b90f283fe8fcc0f483d70f7037074eb3c44da4fb6403a0aef132c665",
      "8d509b1c07d056cb423e4c727945e58e283a484def2e568a493430973a0177d0481f499b25ae83498e77f0b84763891bcc1751e7c2b01de0c5213958",
      "f7608fc768e0fef90f3f3c5904b9434a4e43f1ed8dde9fc43f6906236b04b4623ad0ca633433af4fc9ca676c37ef90782a46e8010e4c541f26b34559",
      "b2ebb028f39f190c4e18a600d32014c17c43d2453e5d86e65aa4a5a96b7f590a1beb594cad76109b52985c120f5c060692a42a60b1103b0ae0279844",
      "1498e62e98b1af80010099127d892c96a1774a780aab89ea54a1c59ca90cc806d5d8cbc6651d585c0e5c2f274698236f8c63973970d078a29a039b8d",
      "a7bc39b2571c382eee9df587f6fb4393beefe46735598fbcbc9d2a7f75c5e0436fe1221b99ba54c251fbea908c600a269aaf3941becd3c56dba400a1",
      "4f99e2dc7dfb33f001b13da5bc2fc5170a2d96be8f42321363666ab54fe4673053b5a82ed355c1792dd9416ba424bb940856d09c04f2076459f28e28",
      "7b4d1a7e12605b1fa2d6ac83d6a26cdf742deac5c6e8303736031bc6c7510fd41a0c916b2bcb9fc9d8952ab08830b5b9d93f66e6f043366a6af43524",
      "1b685137044c51a7c53eeca05cf8bd28a960bfd74ad25d4d6381424d713f8df2b818a2931f870a19d34eec703b277963ec20f20a4d3ca85956d8c429",
      "bff63163980ac85cab0269dc470557861c00f7b1ef1796473bb15ffff4782ff67115eac88dbf10dc1665591ed0e0c30a9c9091347ba7e2affeac8d6f",
      "8c2369ff6b12f867ff2b30f5fa7386ae45670d6c2e85e4efa165483e23a33f8526ca97c8d27b349e28b90e9972b822f3884aec67627235db2027a8df",
      "83de808e12412a4ca800bfb0f566ec47885c044031bf1f12020e760f76a920521b5401a6ad979d94290cc8b99c1c96d0c33a7b588d0f6bf67965f7eb",
      "ff6197b86f25b70ffdd8e0faf417ec93fffcf5093ab6dbd840c9c5665cc9a2df6af0b5a2ad61a2fb6aa0c1d3c25385a9c990918f7a15b57195ce0818",
      "7bb3154b26943fb64954aa1918524d5e1b3725c0c36c34650df17fb38969995dd83153eeb65c34ba90059989a04fdb32d300304d1b4a2a4614f87955",
      "bb6a134b14afbe4765603ea1283961235ce46810044f4669666f0b3c4c435b9a1ca44c5c5ed4d63418c6e4c6948c48cea482973993c4f132a4bd93ab",
      "346ac0b5b3bc6f9baf62ff02f119372b5802a622599bec4ab137f1c024c15532095b8bf1c0049bd459a461f08265d07649ce26337067248be0755227",
      "711123733b6d4d864a9309d967e4c18638b7a21065580a4a90a95c379ea211d10cb1966566e98cdfa876a375abe04b261f99dae63205768a88914441",
      "85091c8d432c569bd272e1e19a08979a81af89225482fca83cc3a997806c5b34d3c85d34698885645ae3bf184d5a311355b970f3b94df33b5b0fe26a",
      "1380fe447ee8f9e78e2475072abb13e19d63c53b50e71d08f60e347c8eaccf33fa9d58ff1c35a0e70f3c920c3a27d62fd387b5fcb0e29fb685c3de71",
      "d861dc3e74d8ad8e7bda61e7f3dbe37d0b3decb36e273eedd35ca91fda1edeffd73fbd7ef9f0f129efe3303fc0a3475632839d65ace549dc6549e54c",
      "307a8e66532ddb17862a354d5905b80ce277a333330efb1c551fe1dc56d12033aba0156077b01d27a49bb6f79b7fe89dfdcdc79f3ffdaed6dd5f93ad",
      "38e70821b98554a4b5ec00263890124cd6550cb7a275589a9614f8c60c73eaf66b5d1ed1671a0fc9c837f7223e7bea736301e3c66577d8217c00e685",
      "b6b375010d78b05ccdc8b96310522e7833d921290c8b1817916d2313389a9a53db195f24553bd0c4968043e804ed859d5a2e9dd17810ca04dc02648b",
      "0317b7922324e4a34532cfa84c30a2c521c378cc9643d060e9795b0ab0012f907be11443fd020c11b008c14736027c1930ce65235213175f434ddf6a",
      "5d15426b788239af3207a60143602364eb28bde37e523072b909103e7187d58e19bc32eaf3c3ca27839ed4dc3238a4fa1df4e687a08d64c5adf2618c",
      "6f9de53973f1e7dcb1edd89e4e90a82aeb214d470da56acea1e629203d84ad87e0f61002dfe2e443307d08b90f81b98fde3d59df81d2cfd3fe1da801",
      "0f0482079ac10319e181b2f0406c88327cd089c66f4da93dfe22ea497b624645ce508dddb585c8f945acc0b469491114ee0a56cac74ce5b5c081b67b",
      "cd3272de2f6a3bb26c7093f80038af7767694741dff0c5766293f7880ee2b89b5f47570493fbbe16895c2103d6ac87760ddba2687a5b482b83a12d31",
      "07c4fab93158ae12b7c531f832c81f35cbef76a0f630d3d2ee914ba9d39266b165f1cd4af530aae5ebc79f3e3c93402c8d7b1348cb6d4bab898c996b",
      "0ef7b673012c4905a1a46100c937fa0ba442910d4184b2955636783ab149a9c26872a715b2e0b8a81a8a4c086b50fa83e958aab0644b94685d4e6c94",
      "b3648ac855818eb26326630e88ce0c23d97ef3440fbdc6bffbedf3ebfbc7391e73155d3b7c5549b4670a8a7458447283d004b4de707a25a999450bf4",
      "8e9ec68e1e12d292f4ba4bdec7840a8b51ea64c2210a2c961184ab06107d062bed3427aca0fc072f559b1ecf5ad11ad2a3573a425c6c18ebd18e4940",
      "09c2cd31bff7df66a1281532cd65d81d168849657219db9376c48af4549239af139ac09900999d319990a6458f3a813b786300fdc13277585656c4a5",
      "63769cb6a24318522be6b78e1479d21ac476e181e9899f32f206e02f97eb80af15edbdbff97c8f014edffff8faf1894ad8f14efdd3b8073ebc138078",
      "911cb3b74d7ee9f5b6e193077be089681a96522cb151d669b0544db6dcd83e0e8b815432c5d4f94586d29bf75b7e3031f8d3ebcfafbffdfdebe3e435",
      "103102a88bab7e23fa6a9946d8bbadf835509d6fda2daaa0b1356cca590cc061a6d29bf6a724826d6a2a75c5722b8c8f7622da83b9730ba8851f6c64",
      "9d45e6b99869c59a32f5a8740a8eca3ab1b0d24ea1a8a2985372f2105c9ae411eea6b4014baeb2980a1e14c94a14ff3038cb601ac28c2bdda10c121d",
      "04b11b1b953575c8b87c238a4cb2844d7b3b256c5328db528d109dcb6611d989c28fc0b44801a478539d4b93082a5c4ccae25ac20994f53b10a1647e",
      "43e96c6a8af5b1b9dba34cb9183573efba79e4c88b31651b9723b8c1f0fac8e3c6afb34c03206ae60d5acedb441e782ab956dd57e61e052f574dcc85",
      "740aca17adb76b4751151e5982c92802cf53948cca39da2dcc9695b232b5395c3c6ea7adeedf2bc5320229ee4bcd5254408ae2772f88fb9302f2b5d4",
      "9b69583cbedca3aa6f91e5a922fa91610c33241b569ddd771cc8756af0418899690a6cb0d18e1ae44ecec0ab4c1d44a560267c873e04726c557ed65a",
      "539299da9e5e43f262133d52e654953d48df4ec2c7b48f5e7f79ffdb8acddf3fb34edeb737bf05fa6df2bc95deb7dbfb967cd8b59721f332631ff2e6",
      "761ec44c3c55e53d51481ff8738f2cbb8e8bd713f61e597defd4bf077660c7207ce0193eb0111f398b3db1b1633f3e7024ff3b1883ffdbd78f1f7efb",
      "f00c8aa00098816c9fbc6310d156fa366bae8aec0505058e2a65cf61c0e0044bd290a8648b7087b0261611c054e44691e5530bf4c06be37e2d681d1b",
      "0cd7882741c5a481cd9bd8c77336820bdc1622b82694c69ae41991a07020897f2523138a60464275236dae0a68cf65ec14c9c2385d4c87f3d862d444",
      "15b10d5c0a533b0cec0880c5d8f4a1c2d6e9498c7cb16514792459ddda24769a746cba74c77116627858842b91b63bfd8fb58819c511593c32e34133",
      "00688bf92a7d5bdccdf26ce86a4adc0cb70f0a80ccc2885ecd045526968fdeeca196ebd0f051a74af1108fc69f5571c7eda33f9660fbe5e7671b92ee",
      "7ca2c07742600ac594ddd8b92c0d9656779f26751f8170b0860050924c3afab54feba60f7065e0fae7cd019d8a1e300837091820b2cb0ceca2a963b6",
      "190a2c6337a437e02a981def06660ff6a6d3e64668a85de16b4ceb02c62be49a7b757b0089ca635a35527b14e4e5b09be66687e6265c5ce9b2bcb026",
      "4eef556039e0991a3c7836b2970dd85f5e167dfb623d0deced8517bbdb116ac0fa8fae34a53da0f64a430b461e191b146edb2e8ea2605339efc4b20a",
      "c6b8c0a9d255d46840130d4ca51843df7c604506de19b2769da70c0dfdfbd43dcc6d37fbdf2e0f6efd702b8c5b83dc2ae5d6b1c34a3728403d017a35",
      "43c3489bfaa36133057442dd1f08f71b476b93a7880ea84983b4d351cf01276adbc271930517186cb4a0936db9b7f03f5b157c6d42ab1cde67b65563",
      "6d0f24d4cd822d23bfd341995798693362340451896f671e89f61c179f63eb737c7e8e88dc51953b32f33bdbb9a343f78ce9070ee1dbbaf098e4e25a",
      "635e9f8a559d5f744f788013a77040658316005a88c519c974ceb70cafb2721446d5bed8bbc8cb0479bad88242e2e25c838aa5e8a2c47280d93aac04",
      "df509a61c54d17068c3d37ad5fd982d75ed68087989ab40b51de5f9f5229906d096bd780652d2b3b25a3119e45cf7cc89bf8dcca29ffe2bdc6ef5fe0",
      "43dfecbf7ffcf0e5fd4f3ffccde7d71f7ff8db5f3f7ce65779f40bd628f5294a67712a5578c69d1325c9c9ab54b3a1254446d9501c47688c6ab67214",
      "959ab234a8db8d92df9c803171ebacec2687a16bb65590a54d4c62035a50db1ee7808f45f7c255f31d85b738c82b8a674d44d9896510e8c48f9664e1",
      "648105bd3c89c2e44d9689ae0edc8ce2c265a9104a660a3fb76da994b9acca6a54f42b6c2dcc688630c8a8d9d55f02497a80fcc8bac9e0b5a2850e8b",
      "0aa0c75c2020408eea2bdf5e31cb9b77fed057ffbf5fbf3c9129a96497e03b59614f572aa3efd73fe8ed55d44db29659ee501575cd15cc526d49af1f",
      "1ace931f644801154f171a97e22426811ad5629348b867178ec0d424a6042d27f3e6761e7a03ffe9eb9f5e3f7c793c59c40669cefc35a92d71128752",
      "9d556c71089e67e6cccf867019d0fca5dfa4260b86dc83fb594dca7d7414dcc73ba1e72c74af08a66099968fa1f27266f7909df4e66e1ecb467f7efd",
      "e7c7d76ab7c8ba75b8a097a931fb67f8096a78d3c36849591d8c8949cf272b7356588c8ee61cf1699173878f887ae7340b50e38c3b24c985744ac18f",
      "77f4f9599e64b9f395c14688964c21c2032e92e2bfc24596cedb54b51562efc0eb676ecc4596b5048d774c1c2b2f04c8efa01bb49c9c68d7898d61c6",
      "6c3bb912d1e845bc0d9d567c33e8bd66a635761e05386e2cdf615bdac0820ed6d53277faa57459aa2556d6bce7b72fd24747c2ae71bc5805f23034fd",
      "8071634a5d84dca50420c597d006642832a642941f7efbf11f1a6dffe7afbfeb19fc351cc071c5748baa5f77fddaecd76fbfc6dff701bf5734e5b8b8",
      "f0d2bfab48e665dece50d2a1a2f20efe7372fd4e59960f48dd87ed044e04ca529715cf5005ea4d0e5a93627985573c68a9f614b0345ab25a90715691",
      "28ae75fdd529a61f2809aa63b6220927f67569d2575256524d38d90128bfd12002c9ca46352a08a855b5c23102e12202bb26431f833e69508f100ec1",
      "1e8f0ca8e0f11559819265b1f789f6598a564b78b3ca43872548cc0d1bf6409b35289f9b7ddfd6a22cca34620834e98cb79cf6a060cb78df3dee155e",
      "0295195648de654820720072426f0a29e4406d881e5816ae68a801255543c018b519771299742c48ddb6673eb56dded462d877d9a87d1c0e3bbd3023",
      "276fc43b2c689c25c6876b616dc2bec212b6979301814bd5448e6f977d6c72fe8fd78f5fd701cfe274a0e6bedb22868981362c7365773c4814a081dc",
      "cdc06fd582258ac32badddecefb7d732a44d6b08ca95b76e82c1c2d604c2a36c63df3656d049d5abfa38a5e5bb206f531382bfa3dcb931d572c3b0eb",
      "ebc3f8df7c98cef8c74f3f3da3b719c0d223d2ee28c6339832386823db75939956e8426285942ccd05b05ba4a9486864ad11e079154183494d06b6aa",
      "903ac1e465602a53240d6b42e817b1cd260aaea4a04c5e40e46f895a153d602a18886439130f78002351204b42b5dec175f954c40b6539301cb5623c",
      "f219986023ee3e67752b9bcc4c0ca042175349b326f1c0ea24091562b538fdf6be1ecbd5fef1f5c32faf1f7f7ae2331d5eede1031c3ed3f163ba4f1e",
      "c97e46be0a437805d4a4e2108585a01e010de183dfa4a8991ca63045d65e8211ba07f251887a8ef917b49a0fbd7ecbe922f990d1088e12950cc9e85e",
      "2c5129d17a1842df974d5d9df0eaaf8165f9b1990df8733f069c2170b254bb9988b2427cc7287a8f6e66f2d1d71ff759d0da7cc7b134cc122b1bc0db",
      "149f00aedc499a9442b9ee0789452c0c6aefc4830f1223afa7ac7a8886e1835d9c65848964f62417a651e80f561f4146c35f9ea073a6e040162df81c",
      "ca2f8342a6d819910437f647424a149c38629786d79fc198274f72fde0ba9c5af9451e885b40e10e5d976a4bc5c7e07c4c24f8b7afa15a476afb6bec",
      "af0c1a10fd10f84a2991a494d7409b2b9595a27071b8978c4c53ddc22c932c37d4425a01a63d6f415d0f35912c4b47be36620b3481a1896e7c50b80c",
      "4a16c8d2793b630f4dbe884cff00a43cf65b0060478275ecc1d1178bacd51ef584f6d01189d720ab05c1528077bd475943233b00fc46d210889f7d47",
      "ba8c740dbc6cf718c61e79cb13a61f51c5f0cba1179148831093205781e01e4acc641b45f0d1307dd886d76d194a68728fe24334d3044800e039ab67",
      "a068d5e89314530709645ce1fb2e4540b7b58e80f602b5ab1ef2b69492a5d89ae7691d7f0c46f2fae9b727d6b7c30e70d8270ebbc961613c3ed1fdb1",
      "0fefe6f006fd6b3e7c8bc3173b7d56f7ed0f23e4388eee83ed30228fe3d60feeccbc7520efe75ee8869cf76b23585f015d27c19e0408d82c77be97bd",
      "1682a19461c1d83b0a025ea6955ad9cfdbd1c8001a1a9147e9d565de646ff673003dd777127ddf6ec540a23292fc711f144766caba292ec2a5d86913",
      "4925b2bffa8a87f982dbe58e4c446af89ea6ab1dd048298f5b4d191c783151d86d1aaf4bc0d253289bb3c6623a39288f0186fefcfaf1d7277260fe13",
      "fa35cc2f737e25f4abe56945bdafba7e653eadde7e85f7bbc06da3f07bc961bbf13b92dfb5ee3bdb7def032744171be8f56e7a130b57514d7192dd4f",
      "a4b8fba41a2447b302c5628f4496605409f543205ae91c0c43f5c209806c97305fe05235d99ecb686f8a2c61226d5ff9428b68a06665dfd13ba2ae38",
      "b52688a1062d6bd25559222ac51419d487012f9122d62422bc0946a7c20f9c449034c9a28631b0a6b5191ae05c1815aa89a11f6c64d574383d274bcf",
      "1c6d4dfc95b050230d8c0daa5aa1b499b2a4288b1990f142182a4de1894aa3920649802cfc540b169bcadb03d4020533a84165bb9d9c2a932141690d",
      "3c38f3b4608b911b896c6c625e7a7d4f7bf0c61b44b1940be8042d574d4a09eb0b931f88216f552e6136260b51774aa2a2c4c7ca34ace5d49c104092",
      "50f19324e444da25b13e3d5408c11c4a9d79d6247da16599c8f990e9d8a6502a8d97e9623d9d2879e11148f86453a803a9adbe1b1bc60dfd7990e814",
      "ad1bc61b0acbecb4cff6530d493bbc1cd54971165f32faabf6481ed89590b81283232651e6ba3d8bad2609442c685788c326de443684ad11696f1a2b",
      "baa0a99ad822bd6318463def21877de6ed4af658c6e0c3fb2feb7f9f583e0f31de21123cc68b2eaa3cc49ec708f510c79e3693fb96e3f7a5c3e6e577",
      "b8c33678d82c8f5bea7ddf3d6cce872dfcbcd1876e99a76197aaa1549bf063ff5e6f1c6a23cebd5f67cb52e69ae2f5dea34668b7f80dd09e5a54482e",
      "d71e0ef9b3445599fd46d76b6316725d40270e6824248ddbcba1ecb14aa2d6b87202556f6dfbbc6ea2b5a472f8d8398c825ec0c85479dc9fb05125af",
      "a9e08adf5bbb1717fc209d0398064920912db40f36495927c1f56696440e4b2ad3ee5bcf4dbcec61fb77dd8837a12952f69021e11f2675b8d21c20af",
      "0769cadcaf9d6a43f3f25ec1ccd1282a156de4dd26c963b27b9f3e7ff9e30fffe9d3e7f74f37dea10563bc2be822610bcabff94f4bd6a14db1c02432",
      "d7d3dfb72bf05f8d90695228eb5d41e397027d7466001f515027b44e40b852eb2d162c1ad78911ad01cb649c273085f536612aa62c04528bf572618a",
      "c394901023749a82e9f46033c051993893b94d50bccbe0f5cfdb9441e39ce7dee4615a5b788629888b12bf88c4f1bb3c809b49664aa1d2945594d289",
      "3428350d43e8b8d2d82496716d486d7937a8650d934e42a7261ad6602a7625f8828528a15e4cfc09dd902a9d491e2d516c92c0a13dc748915979a535",
      "e62f135ab07071935f2534db2ea4952519303353c9385f073558eaa6a00262364ebe0353414b1485a669a8db1c26283cf015dbc5b12045bef46a2c92",
      "f8bd0ca4d5dc5592c8024ee579c6fa1c4977b33615f214ccfdaa2284ee80ad90a86102d00a749c911ceafb9396e562174044f2b78f3538fa4a331631",
      "2c6489d79ab602a2fe0b26e642fee53d90c10e589000d4b07a33111e6492fef4f5af4ceb7f65d73d0d89e3c071c3eb3008dd403d0ce6c39067404118",
      "68919f0a5369c08d010057d2651aaac4e6b4677ac02e8f52ec48fb1723927a0081ce6bfa64c00c507596ea041e8755130043a719d6134f18a6097091",
      "c407d558b60bd73d4cd0d741006badd73811106fcceb5a6fdee8834c23e88b7ea23230305de31a52935cb4c0960de1e6d0e5d1150905e196d1799285",
      "37057768a165c8dd414080d65434aca8c105ae3dd80b605187d4241a71c212446730e10a14fc94d5cd27eaee73c0d20563807b8d04352e63a11a5941",
      "d78f230d6c870cd4f2504513c8971e6ee4214d6da481ae6085a546db80816f5b0b3ef024756ed3fa40bcc31ac34ec13574e0e0a677fea8a5c8e74afd",
      "72d9e072f11eaf14d6404057c5a9609bfb804a1d1e6da3010786164dc9066200e35d2afa1e0211077a2a6871c8ea0887c38aba7381408a482ad156b4",
      "6255353d886e0f741789f79d04a18c54ac5cfb4d05a24f848a94ccae3c717d2399d01c3bf809d613543391929bf7a061ae13d7b8415dda140451f982",
      "ce270299b64f442732100a26802d1ab00ad3f2a9048a467b39488ef1dfb64d687ae0a58aaa12f0ddbaae1e84968089e16e817f96ae13d15883a69450",
      "f72f96b4c6028eea76145cc43e61cad6a30e2f1d647a85e425c95e5741b055c82da1f78c6ee435df0b23f869af9eda9a052416aa34207d8125ae80e9",
      "39ed4f96414b88461ca33ba62b886f866058fb025914e77a6c4a64ec13d77ad5612ac1b427038558b1e80feb88867bb9ee02a62a840c865b07b21cbb",
      "c788db978c50de29200e2b3b3f8ac8027bd334783c40ca20838698a2b9a5ec623653dfd1207c12fede55281960f92d801ae4b1cb2095775015a320be",
      "0606a610d63b15a567e8a0c1220aa6b524748c5b7460e9519040014ebb101e614b4b1b83961554642d3f05f23ee86aba9218110deef808523a9ec0da",
      "63e855e425cc52b1ea17c243bb1d83d2399b9a861d9346181ab0f4c8b040921e1be3956ecd0096065b3c7d721e735b561fc3517d7ca6c67b4c55f874",
      "86cb78dc72223e6de2532b3efde25334a734ce3dd5e3d341a794914f2bdd534f3e3d85ea1fa13acd928fe8a29a531661b581966a02cb28f29e4878f0",
      "906a5939b8670691ec8ce326793ea2217e2c3506912a019e258f0d58ea641e27da6e0455a2a4b7a392f3a0dcafe0894d92cb20ef5cc33e91e623ca92",
      "41ea87b7acbce100408b00f422adc9013f693289d4d52b3cc0f117f9e3c63c3886e8d1c0ddd4da1ead43b8d9aaf4ffb234e44be0f5e66c969a26e14f",
      "4df40318f540b421b6173dd1406e96617c1013350ec1a6874cc2e43ec81b1cca34c9f71da071a8ac3d6c695b88d20a56de451a3350129a040aad97cd",
      "63aa46419c9b166800a817559f49fa7b055fc4d4342d8a438cbc84e6c815a725329eaf95df6a44ee1cefc89da4eb50cb8c499b6a7f67aa9d033c32cd",
      "1220c14aa1213540407d2433dc37851458308ec9f44ba84aa7620245828280859dea298cf56bb7437291de63d6e6d2f1426ba6f2333f5d47404a9af7",
      "a4f68d7544ea9ddcea7c751d511489e253b5bf6b9e20529fea48e800ca2079b09e26ebef68925b46a1d2ab0853b0656903e800c32853d145f1de3198",
      "0a51035d65b58eee6d75e54d69c175f49c660a6eb4281d71a4d08a5eee54b309e076c25c4d09d67674d82305424d9db42da52aabc8c5b583b4277054",
      "650ea18eaa55cd4a74a9a1840aba3c694892be916092cb8611d7b7c9f67266b92a732b6d0a214a1d060ed736c0d0a695649a61051684116721203bb8",
      "917292451f055b61d9e0e328cb72f0886a5e0b41922512cc8fb1ce37d8b05d16e29c856fefe8b3c00a05261dfe788f2a49c012f3f572d05082ae164e",
      "950e7219754bad05a9d887006517b09f92d18425a8dbb4d9c782aec5603bcdda43a21db3261f5b2644ccd9b1801702734ba147d1e1a74c863f39d299",
      "eb14814f86699d369c6a67f74f15a6a2a3b33c76f50371a1edd0f753ab53e196df217f90d83195a5ba80a910d5cb5fc39e3f0988675242727dc4ec00",
      "9d5e26d9cad059a50e92293425113bd8450547ade29feb0835f47ad62e69ab02d8a56189435b3edbf08972ed820e63150006076ca5826c0cf230b332",
      "5084641e40f5e15b29e81b20aa2cdc7d2dbb37b83972afcd7c35039f2172afedd2e01ae0e0d2a81086662047d8abb07c2c15613d078f138692706663",
      "2894449429650c1c53329b994ad1120569c4c08dde7a6d075e0db65aec0fcdaebcc233de71d272014b56079fe0be6b13e8d0a2c71d9bbb467a1a6e53",
      "4371f304fb866e7979524dbb6640750b1611b862c78632106eb0eebdb621fd81ec80ba9e269abff41046e806b70f94f124a40d5650426b24d959ba79",
      "8f057a6308fad53d38198ef363167135a26887f80d5d88f231502d1843df321eea408f91ddbc7efcf997d79fdefff6c7c79dba936be38a73be7ee76b",
      "7c61a8616ced89769940df0b3b73b4f73db296c3289111f83f506e482471d7a7ddcbe34e0c803d926c0e7a910877796fc3b6c6095a307a59216fafa1",
      "2b61bfe24e331408f48291aa69cc9f1c8b9bef71f24fbc0fe3fd1cef0b797fc9fb5437b7cb7b6627efede6e1792ff0e0299ebcc9bbc779f24a9de3ea",
      "7ddbbbffeb5ce4db507b6c88fff1eb97271af58e73fcbe0ef8b5c2af276ecdf1eb925fbafcea06f7b272da87b2077488b224e95a201946e703abaebd",
      "ec54d57963d413f7a77aacd8f2fe9f5e1f67743aaee68715dfef0ab79dc3ef2e6e03f27b94dbc6fc4e77db0afd5ee9f753bfe7fa7d39aa4b233339d8",
      "f543995dd3284573e80faacf452ef74d0b03eaa9eaa612abf0c09a516cfd9fb6110796b2c7c67d0e900969ab819cb48532507b4b6c93c9873df541ca",
      "917ffcf0db3393c87b4fdec1f24e9877d44ecedcdde1f34ee1d171bc3997de01f54eaa7764bdb3eb1ce293cf7c73ab8faef7dd3d772ebc77f37d2480",
      "899599cf0085832cb1546e2f10608205f581c4cda448c8aa01d7138a361cbce406eda76a77c3269386d0b8f390204aa4868462e13ce596dfb238b631",
      "05ed7e3345ac3807a35db683b18c2257fcbeb0a4aea5ac32d26b510cb9895b679625b31b168d860c285bb09f869e9f9e1b1a90367d2257df0606ab90",
      "df59a9454f19e57365715837bc473c3644cdf13758a61b2f523abf0a5e15c24efc90487d1a75293895e1b2cb32cd79cf1ac60d9d6e59dc3a4a30b4ce",
      "5e743af872631a1063c6f320d9da467d561e6322a788fc4862801427ab820d608b54c538b003b4363a439028eab28ee440166f8184851b606941d406",
      "621bc5681b91242053fa5c8856491d52b67e292c518d8092f9eb68ce85d011326dd942bfc2b672ea10dafc58379851420a826c765086375aaab6b80e",
      "f647d59d340fd156c5a29354de3b40b7554c2aea4eeb4417b0b7315bec0f5f9b6d942a55226380bc0a69d0c510f176017a8c6ee41f7efee3ebc7e716",
      "be863eafc186ce2ce50df4dc74dd7e5671115d432c49aa91adb15cd1f4c0d9e6f500c13f4e91c450232e5127a9618bcbcce065459edd0974e265aa0a",
      "085866a0f187ea5dae66c96c01c54f0e4b4a0cb465836d7f282b14d5b6c162e81e0a29eb7666cd764c2d91df684aa61663a10c565543db0990b5a489",
      "2227efac499c592c1d4ae224e39541fda45ba62035d5e2ad64d7b194e9c2455d131d7577fd7d1a1c7ef8f821f67608fa217a1ac67ea8fbe9e066cc7d",
      "52b989779c9c7e02df26f97121b82f167e41f18b8e5f98fce2e517b8d322785f28fd627a5f6edd7aec976cbfacfba55f0940dd0c87628bdc3a444ec4",
      "e45c4395254e23d88c769d96d9b16c81650353911235537a97388bc7a0dd96400d1c93063b9f9b5c045c0771210af20268350402618f94a8a7486018",
      "049f9a048fb1bf8e428b40090d7a7359bc3a82a0b44c594d0ec1bd5b55bce24cae3bdbb63b1467488bc3127843af8148a8b22aa70d8953ad2d4dca4a",
      "ad89d011c7885cf8be643dd657f6faff7e78d247f4cb2084b8b134315a927383bfa79201746e00c1992815db42da432432a92afac7102753e9245459",
      "13835e33f537dade57d3acb0c0f930278ea55d22ba6da12c90722ccc2ed8be0a4e475892006cf0d010abc1528a2da6c84003a0634311809ec90b479b",
      "4d18dca06b018937a3e98e068ac9c76aa1ec942a87c7a4ccab16c1d605c15af193e5ad27f9c606b1dde63d57ec78033d1af5bb8c38212bc38ed12b9f",
      "d32fa47ea1f58bb15fb0dd9aee97fdd3d6e0b70fbfc5dcb721bf559db633bfe5ddb645b7731e3657bf01fb4dfafbe1fb182bc79f3ffffce7bf3c1757",
      "51299b45eea4c414361994750a6599bfff02a8c5abc702a1cd98aa7a8760815503bd4e419b69b47d93c8a542e64173b3b0fcb1a4df2d490e6a3341bd",
      "869cac8654360ca22d45ba1b731427d55d1e6a060328d2b847218a7fef34768caab647010e3b959e3b0100d2ebee40427581da943b5ab73201f42d94",
      "5ddf452474f263c68ad0ac57694a131c972c0c04712a2cb9357b33eb2738f9c634c36c42dd9de79d9f9b6efefa39eed701b7541c5613844285771b95",
      "2d59b7131040010edab68fb55640bc86b5765e8e7c27ae340a7583a58cafaedafdb36f13afbb755b7fde0ca90729203fffe3afef9ff495418d14123f",
      "22d0a4e21c98456042e1562aa50404261490187407503828480916e346a2e1824c8103a3090399958babc8af35821bb3340a2abc1a611b47ba781e00",
      "3b21cedce82d5a2278a4abe25b951910a4379b1b11883989572c0deaf4c69b49c1c2d75a347693c61c226740f50b251a2c0c26f52ef6a8ca39005f28",
      "727864f554c05b428e1a3852f1683404e7dc1f8bc8e61a65a6855aeda51f76f6e3eeef3d04ef45784fc37b23ce61f13e8d737bbc67e4bd27808d3251",
      "896546f30173973f5f8d3509cddb82290a8c03e20faab2770287f479c52138b64a3d383ab0f060f34822030174563c8353c2d19544ab7c396bf8d8f0",
      "1ba0b5e60b4c46e735938e89931b55857642addadec4ec81aed1c0dd7a0a884c060e227e93600915ee7ed78a20085f854aa2f0dac68a02cafeca6162",
      "09555cb71a183cefd9911bd7bc59e3c52dc6115aab48b6a0651cbb2db8695fb7110d6cf4ee9572161ca286a4bbcfd54749981e5f1d3cf5cf891de8ce",
      "20e459864e4c4477b622cf68e4598f3c3392674ff20c4b9e85c93335bdb16c92eda4a0d6a42ccb4e570e4abfc9827e5b6e4d6918b138964a80f9f711",
      "1370a7ccce0f7b86523b01cbc338c51bb70958ba6a7985cd823ca66b1dd43b153a5d7c82449969f4da0e5d400b5009575fa1d2f5b582ae33eda4b69c",
      "293e823d38a024996b9a25bdc9993dcd5913650f36bfaeec89c81b81ef85d027c2f99e8c620f3b3f3c4e31af625a24518cce3d4fee73d4cf633fd7fd",
      "7ae0d70cbfaef8b5e7be3cf915ecb0caf995d0c79aa778f410b3bab8d6c7bef7f8d8c7d03ece3ec5e2f754ad4fe7fa94af4f0bbbccb14f2edfd3cf2e",
      "415da91a4bb2a03c827d19f16781fcc6be1d68572820368d82681a4bd1543116a444682b26e36a3606a69e72d731c59c14757735023e8de271e6a8b2",
      "56b2f513623b221e16e7163c32fc1630d8d57806474d7c59c338ae82524a9a04c620df95569962be2d83824364f23469cbb70be763fd177f7e8667f9",
      "a8dde2f441bc84c85dbcc3ab659c1435eeaa1b5e9903dd8d4afd9a28692c53c9e1249db826ad2ffca0caa3f821f638358cca6eeca7e4e4c627ecc699",
      "ba6e8fab4a93c2e6a00fcd65460496d90afee4324e734b96a666ecc6a63cdad86cd277e805d7308a5cdcd28d075a404f1d78a217740c843792c203a1",
      "ac279dfdf748b9fccd6784038f0fa1c3b6e3b726cc53456359d06858149deddd36b16f0b86683b1e2833b02796f31629072b984e0628c01982557a5c",
      "00bd630b45aa48dd0905828be6bc8be79df06e7946e9a24ca4426509bb250a2a164d1d7a6b3d300723330e43701a4d3d644056fe5b82092e071f669a",
      "6f0136408499c968180b3aa633a3cdd3d6fdd8ccffd3fb9fde3feea591f77e7d0370ef59eb5454731318fb760bd4da0ed6d337ac667d5e1287f11d36",
      "b969bd21c068d232b489471349451c91b6e45945b75b43425f0a8660b68bb294ad2436a38e31820969d9745834e2019245a75b43d1c878d260598790",
      "4ad294093b6ab68d28fa29dd41521bd2a2861aa4d1328fa94d8d0a93d488b014899c260273f33b642105c201000af116924e020d0264046eab36b6f2",
      "06f571d7006a703c4c9dafa25bb491a7716bcaaff1004b10fb153152cb3d68da9f7461346cc212922911aefdb9f1e544bd7674b0e039337261f65403",
      "0e17766595eea95f87b39091536f07525f6bcf69d4053071bdf572fb3b06c4d1ded758de3506c6ac6ac56ef41360e9e2462603ce3a0be1b9306d8034",
      "f7b556c3f5b276f422ba4150ff297715e1f92ecf1281bf2461211657d6c659416b11fb96b65c5b32180443dfaa7aebd2f8ad52d4377b1bba8f818bde",
      "fff2faf9ebe3744389596d3e4e54a3066081484757b22e9844db5aea1b2cad6fddb7810c65a5baa749b40d74f9c16bad5b2115e40d957a25754bccad",
      "5d0eae6edceac419122415d4d3dd847a3bda0ee120cbbd8605dfb666f3db1126a1f76e19da164b053bfb3b46ac712b00231b87a036331791e98ef198",
      "a8b9c1f2113adb40ff6186b5968377b26ccd60cc9d0ab055ec47037428bbf5d4716541a7695407dd8a434c3bb60207076ae5b8df540449778d1bd783",
      "370e76011cd3e82b02e181fa31299acdd0507dc0498969b2441e265aaa0a5a848315fed45ac9f9e11ad8d4793762eb64636fa7654aca80c4ce992759",
      "271b30734177bc7eb2d86c07b0933fb5052e57e8c5074f7b4188a0dfa75bd965584b637ac7a441b415628018be52a4a85fd7e5672a36ae9018c62350",
      "ead3eea6544c2e6aa1753d034b2aa01515fbb56088b4485f156d38e092ab68b5d1ba07b7b9d0d2846781038d040287301de9fb54788cd1fd1f3fbf3e",
      "2585e406151453270762300de8a45e36849512b9a0175bf9fe53167fbe390615a41cd51495f2547b67978e4a468193ad9c55893f78b1e85a6273e716",
      "8d46f61906d38984447a67b367efdbbd056c80ed9f22dda1ffd9d8d69aec3acce6b1af526590ac4d8a4dacaaa6a14e0497844141353d7aa4c2e7261b",
      "60812cb0a132db0c85381b0a0958dcb329d60fb0fb43144b7b3b5ca98cfc12f0cf636b16a00453d812db64499893d4556b767f2b3466c76554af2eab",
      "4e3a260812058fb9cea94ecab145ee1179c152b6b40053d76d33b0e18d82e084ae6237f9ee013591c2419b4c4a71ed93e5122fa67b0cefac510f4c9f",
      "b7617ae1c291db22f066a862d4ef14c613daff0b8569f7ca09416a7470e6ad340dfe2758baf67aea800e3e4215ac3b5b2e80fedd96f46c41ad9a26da",
      "03ce07349a237f23881bd627ec937023e92d6179cabac122a910e018d90d5bd1b2bf4fd223d81933eb6d26c50f8922663094bc9541971bc327aa82df",
      "91413bf024130ba1f825ffb6c536a4c27f2f692b90aeb193780df5ff20304593049b5f8b3d7105955869bbe5061880a8e1904b330ba8263888d41685",
      "7b49f8d6ccf0c5fd5670cc00ad905912c73daaca1c8c8969195a4ae31e88201d9957e4b06c3b609a87dd9429055be122d20bb24c2d8293bda5949d8e",
      "b22410f0722a680926a69b9303640eb6bc825e8367455359467f080d752fed2ca18160bf99e27107295061cf6492a59530f458c12c6168661653614e",
      "d6243a300f8b2927335b3c71ac0948775634265d213b66451feab4dfce68496a3637aa35f89e2959fbbb8623340db968a9c1f9b083f95dceef847eb7",
      "741baadf73fdbe7cdabb6fdbbbf300bc9770f424bcb771f748deee2f8fa16d3f3dc5d07bd8c44f1bbd7706bcc3e07c0ae7751c3eebedc3df87861f3e",
      "7e881d87e16da8e2dbf561141016490d710264a53ad98d577548567b395293fca5a4dc5824be5aa4044a4c81b31db47b12101db20c0c4c8a8c2acc24",
      "34859be51cc3a2a4861e713030140b77eadadcf8982a3653ca5e740ba6f246d8bfd80fb42993cfb065e391101b373c0c5e66531100c0c8edbf9b983d",
      "b638fa19685db60829e8abf422be04f4a26b548e6921d35ad2392c9b3e0be409a1a350afb430c8196aa6e39fc4198c406bad7a0c05aafdd20888b3a0",
      "891eec0932742f97a50a108d37c31100f5337d04d45103cf6aaab5c2c27006b7a5989c34378d16716ce0eb4ef846520eb55182a4012c6a0f81ee7b91",
      "679a423bf8b38fa5aebefef6e599ec2755e0b17ba128902d7b30c1a58c75dfa8080090935364894134fc82fd0a922836ca10e1ca9b19c61b86661ea3",
      "8930a57634ea74eec8cde8d30a3864aa766dc5b598eb8d96b5615ae2247759bab1cac1ebcddc7662344e0af4a29265c0146e111f177938dd58890210",
      "58bc9f21de01c01752d1c61db6e299b97512016c603fe2366e0268cbd94edcc5553067642bdf50597afc886d913b19038445e1fe679e34724c047be0",
      "32d95200b1445aac279f29002349e83b0fb5ee961fa1d86bc848cff0fe5ade2f06da0bdce8c7a54d0f86b772f585482b936ec6904b8eb9038c2452ef",
      "923603af09e162e81d10792dba12b21639b59de1f3925a034d8e3b1155a61c8fb5be5a8a6492f003febbcdb8c14dbdd85f6fc6de4303fe3f7efdf8f3",
      "eb1ae40f6f38275fc3fb23de6709737eff3db0a4270dcb2982f884940c7dc51cf72ab0a6362d51733e92818a06756761b94e7285a7f19f9c2ca58a10",
      "a54bfd958a1a819659e34e4f0d4da42c871e53d4a6440d3b8595d1c25dda96863c4e5a3fb1fde4f70bc47d1101fde2d01b6cdd2ce84fe45995719e16",
      "527bcbd1f2709d57d9d29d93e4c9728963d94bad466632ba45b8f0729293a02478ed9036844b9c382fb0b7f6cc6f5e957fc6021d22476633573faa9b",
      "90671d7cd207d9a37ef9e9d3bf3cbe1a1f7df7bb7fef63001f27f858c2c71ba798c4c72d3eb6f1f1cf2d463ac651dfc559c738ec1eabb970ee14f1bd",
      "8d087dc4e8a3cae5e9159d53ec4545dd981c2afc9dc45623e91ebced1ab9e683365216d6779922900788b753920dff7c88c11e1a3cffe5d3afafcf54",
      "a18ecf7b0bb47d2c7e8ad7b5b05570aa5f67e135811075671847eee40b9aaa3050370d3e6f06ada9bdb949785f46d9c82c09e007a03f342f9b08c328",
      "12962d8aaa28c2a39a3586455a4c9893c269bbfa19a2bf9481b740909abf306843c5928e76480abe172b30942eba281579c97f53497054bb9dc39018",
      "725bdb658b5091c2397dd89a5f1baf6a5932b488a38c85a299883cc88f51f83ed72b2bb2f49692b23aa51df705bf75f8ddc5ef407e97f2d1b70bd07d",
      "0cefe37c9f0bf0f9029f5338ad5df779e1e78e9f5f6e0afa39ea67f19b91ff98f0c2872f7ffcfadc84f379629f4b76e9669f91f659eb53507c0b167c",
      "3c41fd645aaa6837700cf60e4633cdce5ab125ab4f41cdbc28b3f1fec0f01cf658ac9555ac202600246320bc847a94f8f124a39659d7eab9ef511498",
      "0fcf3dd9280a01915323a5fb7ef4acf258c88744fb8382195ffee5998fd691184f2ca5f562cb92ca78496879605a80c9adcce19e92fdc782801f0cee",
      "71fd1bf16fcdbf59fff6fd1742413f14def2d8eecc5a8e553114b573226b5c86a58b3b100b297219a8605a35172aae41354d09212546b0bccedaaaad",
      "80b2eebec1b2d6771be32d206a66f5b5dacb60648d97ac55f0f6da1ffad47ffbdb974f1f9ffbd66b918bace6f76ab771b7a037022909540e822d4905",
      "9aded4032c9676ed0054f225e865fa51e41fd9bf16ffeafceb2d5451c1afe7b12da145969e53dd35b615413658d63f8cbd27a2f43cd96ca0b31af943",
      "81b0a896f41d4d96264d4676794fbd0d81f4fffa1b7be8c3fdefef3f2f5fe6f120ed94bcf1f91d9f03f27922974aba679b0e09a90ae71464a821eddc",
      "920c51f85afc74550d6b26519fe2f50f25bab251120384a7049ae41097854ccdc8a099da1e6a75622b9d9bce79cdd128d6d16ef88f90457bb8a1003e",
      "7ef731fe290fe07305f77cc229e770cf4bf8dc85cf6fb81c884b93b844ca10d00c7e8158d9070081744a9ac01e9dccce34680784764866b4a87a157c",
      "9ec937b53e28219ec832890b561da4d66b44da0a1d9048af8a90465d9e62f8007051a86480a49bb1f266c9c24fe8692289aeb60a6cf8b3c8a04b0c82",
      "39d256b8ec845a114bc07009c173e5dfc5fa83b03d10c4d004e8216e865883aeae34e2753919aa5afe40288c1159769730e5ebe50c1807209a77987f",
      "5c51bc525b50438227904560044c31d39ad3189ad6af4ee63057f824ca5d6872eaa2225f8005c81ec1f2bbc19caa72a15160e5c8fc383d95624cd3e8",
      "d6129223884ee936bb1f83ae7cfde5e7d767129dc790e01e36b8c8c2071f3e40f1418c8b737c28e4c3a5c1cc9a5ad8b205fa0d2899421176db75a611",
      "980755b9e187e4aa7e08318a881bdb5a83cc2f8f90a365fbd0b08d60bd7e75f5e4ed8794999abad6b2553e1ac62ddb90a67981b326b521ed206b906c",
      "1ae8b76999b40cbd5c363c05b380d105fd732d5ad002bd133ed4d887b4205edd31778d85c5e9f52662b4d8a7020989972ef2185c86556f0ec57808f9",
      "1edba73ebf7ffff7ef9fd462008a7b02b86a6e20185550ca0ff66681af44b1df76e5aeb6a2ab1a5923ffdedb381b8afa3e77f91f3464054b805c081b",
      "bc5decde448d6e913afba1d8d66f89897d67927df0dfd07f663f14fc70f143ea30eae87e94b8fb177933457dcab9d9bb0a551d8d43fd71a07501a6ae",
      "104c6eef807a2e384bc4d57847e8bd0714d96251c4b4345cc112fb1795c0c4a42e597fdb40196a8556521f0327b52ce0ec8e9206a9f197a3640100d0",
      "62ecddde0b47cfecee9da28323e691edcc763ffc4c43e219f8937d18220ca522b6faabecf68aba94e32e2656eb5eb383a1a809d077da5989428c76d9",
      "7902a2a8c58c8f686ef2f4b6ff8e6a51899a6bf04f04b1b53c0478fb32590c2c7b110da55b0dee1909058e1bfdd9f9a7ad73cb452206780725d2fa90",
      "6a12ef4c2da133479be86a004577e7dc21097a18d8a91576998c8017f2350cbfec3cdee1f9fecf4fceeb63b3d35dbadecbdbe30b4ac723abddf78887",
      "bf63e61daede63ef0ff87c8fe177307fdf0970ea16f01d05f7ae0304ba8910feb42f135bd02f89542403cd19787fb555c38cf568edf493839f4c2162",
      "bf114e0a2bbbe84ca281c8d8bcc3739ae2712017d8ff849433f741c08e66d4cdb48de3e135d47c0cceb3ac76ced6779aa981c697773bcc128a35d25b",
      "38dc28ccc8796b7b4aa3901e27ba2a5ad4940843b8f5dc7624a8c66d4b8c8a78c192e51ab97a568c7ba9e62cefa919f4b470b077ad88e87c8aa28117",
      "334aaefac26c180e862c9ba06ee17ca88647cbb6684ea9b460ec44ade156ca2ce4f3e3bd84484f963d9585d7a959cd38490811adfcbfd32f70ea29b8",
      "f71dbced4bf05d0bbef5cfb707fa1642df66f86676da06e7f621ef2079270a5592ca076f29ed62c3909ecc9ab16619414dd862b1c18852c3c3946a14",
      "7e69263eb991752736b770f9ac86cca2dc168e508872d8014fbba4df49bf7fc8078102bffcf85cb2d3afefdf2dfef79dc16d1cdf6f2b7192ba032acf",
      "fb4f09036dd84ba173810115ac0c98e38cc67992f631995fac18f21f21336792fd7beeda4646df95d7ef2c760671df1819b31a8ea3238ac74430150b",
      "a878c466dc13f62881b52bf8af6367f3471edf3bbde86c1427813a46609841fdfab596dfb3ac15459f395ecd16d4cf8a9bfce5ed37784ce9fef3a7d7",
      "2fcfe1437c19f95e6af6e5685fb2f6656dd6fc68a9697ff98cde37d4a27233cb5a452a2b4f02bf45aa95b29ab2d6e8ddf36118d920520810d4330181",
      "ebec7baea0e285c534456cb1678f55b1ee8d092cf57767e1eef9db129424be6976a964440341142e7e65b309a3141e9b4295eb6e122b25207d54de07",
      "815d215b495055154fc551d6ea7e7bd219f17fa15a2b8612fbb5ac4069236e2a7fa9aa8307881b040454fac86621eb4c2129a5e5c4660ab474e1c305",
      "ff603965edc951960462e2efde0d149e5396704bb21453ab3c24e69d146b5d2f74a7c0b27a4df1c9951247ce2e26d5239ba5f5d680e1212def24def2",
      "b9748c812608b55415536a25b41495d8a27dbad62387d2c8c51e7cb902aa1fe7110f0889073bce3e7cf9cbfbcf4f611e0f30220734ba01916e40258f",
      "65a2c261e631c5c44e0566644d51baa5ec98e611daba81ee518d3de88b41bbf83b20f950df32ca94cd7270bd7258403f807fa72e45982edd49c42106",
      "3457e1b071dd141280be265043a919e85d7f430ab950c94b49397d26b8555d67f42ce5a2aad627eeb3aa1816919fde53832e79e8d28bb74ff058e1e9",
      "ebfff7fed71f3f7dfdfcf3e3ac613ec1e85290b724251a2b0d39aa8d1723bb28479cb8b4d233179a7386ecd3988fb686fdfce1ebaf4f3cd92d53ea1e",
      "d5dda67b9022502f81aaea3c85e8875e8e620eb26011a6aa38807565c14b856f10a5130c537e3ffc69a2a9a3640e51f222083aabe0f1d7fec6c2ae82",
      "416a22d740af1991d65d2fbf59b120db787ef31e1e64aafef2472d29bf3df1fe5d5af99e77be25a6ef99eb436adb25bfdf3eebed55b857f5d7dff580",
      "3c4225c7a30d80c0b6bb15036b2d08482120f12d7eeddb533e0850fffce5ebcf4f5081bf00f3a8306dc8e179a1d6f434ddaf2e4b6acad52691f9bdb0",
      "f996de5b9f4c12bd74e67064e196f982ac41940e98ba275fc00a3b9a7e8b51e90bdc2c096066a1e25e00098bb444d5f160510896f3fea9607c625199",
      "a11726069ae2b0b86f8712a070cbe943e31839e641c2782fe42d89e2e0cb715f783285045e585922790fa04899edc295ddfe9d954cfba92ad54a68d8",
      "cb32a0dc9c2f0efa17d43b93684f44a4f142256e46f1c0a2cbb2a99c5a8d76d650e428f7e105750ef1081affe80b79ccd4899f19eebd602955d6cb3e",
      "c3442a204b49961bf20b327b929b8dca00e2a7b92f01cdb307013eac78dd8a1d33679030a9e2a0fbd079ccc5f8a7d70f8ff7b41fdfb5ff1efe9bf9ef",
      "7af8f6f7e17118416e94dd47a21baca7f1ecc7fc7d5ef8b9e3e7979f83a7797a9fcba7f93e217c81a04fba1c2fd39a2ae0c2377bf4f5ba934a457b2e",
      "77163fd0ffcefdeaa5128c414ba757fa0276dd3ea446a8e7c262d77956a9acdbbf444363a3eaaf5f87f26de2593a29b0d3888522eea4811a24647fb3",
      "d096782bc6a845c92e9558808ad3521db2523f43af2e710d56ee8a214090b63b13487c4d08c623d3e6adf37983e4e89886d7a0029fe5d04cd325a436",
      "3076a0f1120861e5fc158fd14b30065fe4eb648968b3a0a548c2f185421b55d93dc6152f59bd33caee4d7baf6d580650b360ed56d0f406e388c4135f",
      "d8aecaeac92cf4165e80901019d39a0d36ee56282a9e1285732f80b98624d2d0745add1e64177bff540cf142fc1d7ba286caf72f70a58990304a678c",
      "3f50e191612ad87c585139dda1d6b4ac5128567d6cda835f20eb85922eb86e98c57c0183f3604316ea7db656ab4ba6b5bdbba07f8d28b21c6d3ac06b",
      "673979af273a046816ddddc45b635938db9b7df3480ffa4e7ffae13fbcfef2fea7a7804211bc13e3dd0b122cca17ac71c33f15f77794515e00086df6",
      "af6b73c0dfc5fea4881c0f48d6e88068b5145e634d7e9d569403812d16239128d05758264474d61551d45901db9a72fbd4629664970f4967b5fd73fa",
      "f7617f7e7b9c0735ac7ff9f4eba78f3ffcdd6f7f78569abe81ff72be7b41902a20224c19eb0d6c6b21d25b5b2fd70cf6678e783910ce3369eb86b51d",
      "6f9af16edb977a7375d3b307e7d73a790d553b71c514057f77edcdb12148589b3cd6721542605a8bc68b01fef1e78a3ebefd098a225c21ee3ff70fec",
      "9f649132c136ba7546e07bac40f385ade64a3b40b0c96e458c24c05ff1efb6ffece84b7ba1e4ac5d65c5f23014cb9280ab0453102709abef7edd6eca",
      "2e0e42b7fdeb6b33c143c9678745ef1464d3639b22a0f52fc8835bb26440c02bbe7b41c8af8f01f859c43420c5f234d30aae30d5f74983ffae31fbfd",
      "adecbb0311d93bfa049a4915bec9a0298b242c92946ed2546c1c80d92cf0ca5179864866335c0ae982ba4d15b1341668c1116e3ff7f882f23fdfbffe",
      "2b8bf35f9d08908d8122e50be90a94ab6cec4a85a5d9c3819838aec0f8459c583a0c79362c01284908cb0fdbda5e276d25e9992112b9e2219daac417",
      "4e4dbdf1274a9dc90e9bc04fbce479655f491e5c68aadf9990d37b2149b9dd30e070cb197a21f9962534c122992b6d418824de099640963d951b0529",
      "332614eb95972963f266e663c7fed1e566c26432f4240a2e3891922e3a91dcddcb04316521da7058e83415ebe7841e267cc217c22fda7ef6cca50057",
      "53059e371b75676bd8147b6d191b161e2005fb05c02af02a49756897ab14aad6534d7b6fc81c6376823551d97c7cd58e743e6cc3f089bdc39519b4cd",
      "66644424bb283a575a3eb055c8c4c1d6a65d0fc493a5f39dd791865d6fa2a18bdfcb20851db7dfed13d63dc0a63e7429da9c3a1462b1e48012a30b19",
      "791fac7b9143c315f692f2ffb7f62e3bd224d796deabf0013208bb5f8682ba05a16f03a935d040839f6481ac3e3c55449d62377800bdbbfc5b6b5b64",
      "4446164ee70f8deacf5dee1eeee6e666fbb2f65af5de1d08693412ae37f99f4ec942dd5c341549b9c4fd51b16252e07f07c9d4651351fb4d750ef716",
      "4e9ccdcb11bb89f8759c173eb56388c4be9d178e18940ecbf749509a4ea51c7aa6cf2609a4c3cc60ada9d8a6bfa7daee1f458113557777a6cfa2939b",
      "6728ef572bf44d618bc1a3ed82d983a0f674fa995b5b29cbb64c01cf711038db66ae0e9ebe4192785342dec53071a06be41aa9be1c232cf1d7db43b9",
      "64d067c572a636aae8afdc06acddc42a121b9524167ddcc87eb9680dc09a89adecd840541cdfb69de3d4759bf5bbd752709c9fd4342a087f9bc40ab9",
      "035a00b9e71c8d83503932cf8044d4fbed3d4f9eefe869fef6d7dff4e7fecd0558e4afbb6a2aacd05e0f4cac6d86386711c0eead977c2d1bc7069dae",
      "6ccdcb0bdd2fb9f92d5f5b610f5b07d2adf961d2964c567ea6aae3865b7eb05d5fa86d7d05151839d2e521ebf16952d3d012d94454b5c256767c5bc1",
      "a69ca9d6e4e5a5fff289eeb656fd0907f692ceaf1d96edf9c6f4192dce7cbfe3edb90aea399d277b1ebc5811a83eeb0b46a1ccbb38b8b4c5e9700d99",
      "ce20035e03b5816d470d50c0282d9ef508ca6610685d3b07a04b4f4d41b4b24cdb5ca2d9047159cb7d72d09b8576d20ec35e66e07312acd0a75ac83a",
      "53c405292c5bd4b35043deba1304dcbcbea8d8c81b674731896ececbac843731278f1536e05dda64865d76d5f87ad8acc897f3308b5bd87c39d0a3cd",
      "3f515774838be06ec97640b4d89a061435c97a8eebade7183c3b7314b572f3a305216d56e94b0b3e4268e99cab3ab6de59badfcbf57eb1c16a52ce71",
      "9763c9b04830c4d723709e32a5ea0f41b77785171429bcb103ac67cb16fb6b0c7bb3b81cb69ee35941c4b28f977b632ab6eb59afdf2ca2b6f610c3fd",
      "81df4a43550be4b95a1c64ca3986a94200d8e2b07dce6472d37c51bced52aa44f703db36a0941f2dda0185b758e7e60aa55b8e9bd6e8d02f90e5b881",
      "438ced9eba26b488d8469c5acd6185e97267cf3374c5686234b85f8dcd598c88314ab83fb5c9e6623a07c91351cb62cca5663c962252678574bbbe8d",
      "a28cd1fd3544e4b75c54e165b5a11b91427c09db6ed3d1691a01e5969be620729a5e8c79ad6400b65182eb521c84ddc735bbe319a86384acd5992d6c",
      "0bee666cc53a00d8f69a0e77536be7f32c9def091c48771c25dc452fb6cd009613f00cdbc6aa3356800a240edbe530afb00da86f6eb493277b2cf450",
      "b53c65cb3bba1d480132ce9063ee756c95e6436cf119f35ec052dda8f5bb2e8aad9170d361234e65e767b9047de051b17e6cd1a9fdc0e38180e3db41",
      "9f923d418baaeb55b65a8263b451e3e28e5932ec15080bab11a030efaf582d452c5aea75f1d75e4cbd285bcfe91c37e4de69e6fb33e64bd20a959585",
      "acc7a6288acfd20b8a9a1a3502aaf97b902587b3658bc083b6b9c44ec82f449856b4b14ed996e9eab3b055ac95d8e2559040ce36ed11af710bda235b",
      "8baab714ba96ef640426032c20a94b6cb547d38cf54c65cb2e1de66ad0d5652b6ac0942d9bc90d5b0e8e1251486a50ca9dc203942254683ace3c0498",
      "809bdf32d43e35ce5cb4bd632a33501208c5c461d793fb6acc283d19fce47127920f5bb2a1e4211bdd1cd57712044dd8c65c3e7787c357451ceca728",
      "060062934c878e3312340bf0c96aa9ebed63bbdce1aee342e802db46935423d08f492cfa377516dd0fab579ca641be42b2732b1b2c8e5e863102593a",
      "74d32f6d46864368d21ab3319e0c50fbf25499c1ea52bb71ff37ef41f33cc5eefe08505c934db3b1c976795e2b6cc86ae9a30a55d25c9795bdb1ad14",
      "242ac453fdd81cf4123c258d2817b61b8c3800243fd8a211398bed58e7d220e8dfb55be9736b7475e155cdec8fdec84a4c85f28abe79abe96509e12c",
      "3fda6ae185a1fc50e29b6fed1cd7c17ae8fb1ed1fb764dc31a5f5539ce1f21d6f6477a7d55e7560620dd9b6867edbd90f56f9a2cf099e6f3bb63d9d4",
      "62a4a8b494987b075527ad9be5f953bd340869593df592210b5941a13e20c4caa27f8f084f2b97547f5ad816500e7f067177c2cb6f4fe5488762eb2b",
      "a668e9e7b87e799c9ea3d59125364596b805e10ba91c1373745abf18dbf6b948c344a0c0e0564fc86bc52e615b9bc8573025cf6f5c6c45e64cd2754c",
      "65fa9d5d73ba1cbf7e75aff165f5e3ebfb0748d4c7abe5bbc51912c9ca386142864547b32c1882091d7c2ee18c977801f2709ad8657b3c04dda89a79",
      "fbd0fb13d9647f906cb87626e82654285204ad2e619b797aa31fd6a3cb34d3811ab8a96f3c1d1b1520392235825a501e3372e4d945abac5648f95cac",
      "9df6c3bb209d3eb725af04000eb69c894a36d5a94556595c7d89b4db9ba2ad7386d3d5bab74daa16bddb872b413c23625ff975ac447e3651fb66fb7f",
      "23c870843b2876444bace9400f6097b1c3ea71161c81c92266be7e6c95ee2df9c47d9ddf4843e149423cf26eb3d34d722258b1687ccae1885b5b085b",
      "eb72d7ed857409c0d8f94fe77ec96534473f914dcd664e557891d319cda9a01528f9ace7ed54da016fd27a895c2e6fac3b29d6dc30c99bddd4db4944",
      "e5488a61d3676bd2d31a36885e6403947d626ab9032a32c7b94a0e3a8955b2dd0181da18cd77e93fcde43c9d9f4a914bd007b423017602633c089b5a",
      "9e27e04d4863e96783c10880f9ea3e2eadb8653c0859ae9b3b4bd9446c5a3f109971e96931b5d939d6591793f284eb50026609182ddbae9dbdc676b1",
      "77dcc735b177d896cf1d92950b5b9f27ff69fe5b6c092804af67c6aa2ad11625f6e017880d7d91d173d4990da7cb6a5e9831056205e5ade4edd79ddc",
      "0fc136baded3a2fe9025f2ab9c30409c79b6d66b9375a07cf950c77645ea9e8ca3ddcfbd02bfeaa8f378115b61ad5a06ce8fd6592288df0e93701922",
      "10e7a2c7a1499a9e625e99c75623d63f671260567f0025d6182ebcd82c90eab0a00136c114b15dcb6709bf4c500d6c2962c22acae4acef7007d71fc4",
      "df3487ddc409de8f67dac1f563eb3b5c6e883c15382bbfd9c321ce2aece8c2e5ddbd7e489e7cb138f7cb8fbffb4fdf7efaa7af173817fbe97a53d7b6",
      "faad7117aeaf18b4a8dcf7454bd27a5309365bd9b0b221488dae590a71d1c9a8bc760f6d441c31269c9b9043ed689f5e86cbb0c8776f518cf847ae35",
      "a1be518656f0f77c5b5f03b4ffe57abcef4d0c2691199437878346d42ec246997a7a37116f814b3780831333ad08e2e6745a220952fea6626d241c28",
      "350e9952f849592bb3cebb066dc7511db89be2dfbba932c9614633e3d0cb9d3a09f7aa4249f69914a465989734160546c90e0ae7962890199b744f2b",
      "aba352e450722a167bce129cb7576523974e576fbf5becc55f2160adbb2116720d3ea0f538624abd10d0ba3e9745ec6bd2bfeeccfbc23136ad69f79f",
      "b39be2a06c7d9d4bb01d9fb2b5382f584d9bb885927dea9d8220981289d91c93e58fc8be9a4f6e930cec66cfb6e40ac7b4ee0e83250f6403fb318639",
      "3b4ad8559e9246cae9974d68bdb671fffa8e37845a631b579fc33211b5a62860b59a0d719ed51fa7a9c437a39485c29996d4de346a7523770c4ac989",
      "ee8e2e3cd218c1b00083b0fa497af0ed257ea518fc14f16d6a52ecb2c9bb57229a763fd7f6ee9888f8770db0930f9a107f2f3718d8f74faa44bbeb25",
      "476e3499565e7aaa9142cb023b18ba149951659d86315091b5cbea8350878459ce3169696972dbda39719a852c47374156a3b8a6678a960d9c7881b7",
      "9a6a7091f06a6298ece23f39b9ade50e82d32da4d4d6dc6e59b1739655c5363b69782c7905751781603f296165ed45e31e47213a5234f1526c309c38",
      "0c519f6d1e531a66e0b876fca71cb166ec3c29b76a06bed5ef49ded97cf559eeb96013dbad2367c87959e82c0af7ae2f906d9e26140e426a7e704de3",
      "a9a35b254be04e24beb9dd336d659a0b3844eb48b4655a3df1b90d63ccea1b5ee2fa0f7939998a65057ab69fc5890213c13b654a50f28ba90ab15b82",
      "8655b9b82ef68716c5c322251c533ef7281f905c5fc51a1c2783b59b81497347164a1d1c82201f980794362ca138c7e71ea631a12267d927573520cc",
      "6ef012daa941b7e40eb47726b66a4c4c681ab01e4933bab5601e6f85deb0a18e05e356b3f4a94c37990d56ccaa018c6aea53cf40551ebc82f558ecf1",
      "b1f459231d574eb2688ea00ccde35c3e6d2fdc976b9b4fea2978e86a096f865eb12e3e841d2e1f2081251abf1cd52bb5a43675b385dbaa504fabcf4e",
      "778b256ca38683636b95db40e9cb218a8eca08ff44d166815a0e3e1e52724b25dd1af1162af52a16a814ac3980438dfcbc629a937eccd6c3cde1934a",
      "266f6fcb4fed7bb569747711469a8d22d29846a7d650da61d9b3e47a74916a82427d2d69da794cdda8d66a896da599e7580608a79351ce02a48afdfa",
      "94a308dfdc87ebb49eb8ff7c54dff5bea454f73baec8c950782ae61728d1cfc4fad12c7f59424a88ceaf5a660099c7398af42a718775b8bc882d0938",
      "5ecf734fba8736640bd7340b73b42c1a7ebfd555dcdbbe227f44b141804dc85cbd8214e53325af7bc54a679c0dab44f02a4ecba945af7564b893321d",
      "75dc795cb5b236698b5e83bcc294a102842b23edb31415b8e00966ce774a2522358bc08d530528655b12350ed22698dfd44cee0c02151114caf88203",
      "50a01c4693fae0889c42d1540cd1c2580e8b3b29c4735ace0cd9a096a02dcdf7a203086166f9291c081131df4430e6848f922bd71640f6bcd753fde8",
      "64c3c5dc914e7d71d0148cc640b040b1772c9c4c4a4fed6c0b57bc745dbed0577336c3eb3116a631eec5b54c71a548f5f9d4a5a4442f6a9675ea5297",
      "1b93dfc4d090ef3b323e0003d2fda9c8ade5564953cee3fc369e9a602187eb70794afb4d4db1e57813836e3d6d055efa386ab4b8921328492d5aeb4d",
      "73ce49e7c4ec4ec5b7e04f2c4d7750f3328d2fc7a5e9144078c21c97573359d3d0c4ea933a279437517778194dac95d74e216e0ffb1c89789391a7fe",
      "ecf444a286d3348033fa5f5351154f96752cf4a063ba3efce3b3410aa48252dc15e270c86a66911a87a968604ad4db1258100c2d94c930689cfcdd26",
      "dcb0d13d24faed4db35495657ac2eeed50b428355bc3cbac7e763702aab162be992e74da3270a344352757fff25507d407455448cd16e8f024521a1d",
      "8682ce4b9274bb1ab349191649874ebb0e1c4321029e4f8fd62654be6293226fa9f8ca9b7703485181d0065178ed707c4d43bedfa6516e4978d53c4c",
      "9ba214311655c63064c84c948073d7a494baa5251a28a54d4a776008d14e7c7df81620caf5a6bfd9cdab2453976984376bc8d0134cf79f6e249096a4",
      "4e832f63e35a5c71109611c148c2b1e709b616d8a5406ef924dddf9252aa06a218b0b858da936e272d77642aed2379d950cd5e54679ad690ebbd44cc",
      "b58a1fa224b73bda67f663f9b7c48a14c3a5b575f162a5c51af2c58b14e6d43ad79c165f45dca25a0ccd16bab24568ab76bd698b41fd2eb939d29ca2",
      "73482a416111e0df0abd592f78121e2e09aa0e7f1e932dd94abf810c996c7d4bcbfae543349f757dc4125ded46f74db636e85700f9693d98b4784e73",
      "2469959ff8cac554262ef403e1dbe62019f6442e4b36dd44d6e4072f69c5e3e2bf9636418d1b7fd135acbe184d688064dd3d4791c0048206cea50a3d",
      "b26c19d3ea9f97d35be2183b27c37f2dc4b4e415345f949e5ab5c998566f52c79fcb5e8fd6ab2958a2a98c5c3606ea274e832c315a8f108490ee52d2",
      "3b045d487d46341fda1bc91f6447aad1b9f1985ff822e3ceb71fffc777e8a5cbd798282b53f46c672fcef41b6b194fc71f20c27a93ea64781b2253ce",
      "5eb3db3d80413996b0a28fbbd7b059ce1ab9ddb33f933c2d0257ddb7751c10613de6f1659e6eeb6bbd1dbf7efbeb3fbe3b45f5d2e1fdaa68f1aa7af1",
      "aa8cf1aa9ef1aab0e1f2a03abff338ddfd50698a4c3d3af74d014e463abc2e6a08b2d4830080b1d4e1f9890c58bf453f300255066dc6565492763e0a",
      "a73002091717683475df284d125522f629357fedde0f570399392ce348195ce191baf12e3f2a380116f83575e31db606a44af914e6910bedd93ab6d7",
      "4f8704015ea53e3b53e943bfb0ccba91c7d16328bbea5b5e51f4853a7b990367ce2374a1b04dafb11cb9503743cd83e9a5ae680ddc48fee01716ad2d",
      "3b1dad85011a81f0242a063025a999ef0abf76501d2c5388edf047b0e82ad73a1cc22b0c3a91434f77cd58b7fb6d77f00849b614e2f476f830d05915",
      "934f737c2a0c94695da2e18142417513580aa0e3b5f224f74906860a9f7579013a4ab824e8d45bd68fda895bc7a872857b4b86637b09eac12adab32d",
      "3bf491aed116eba21945616dd7bc326dce7221959cdd08432b8e7e4bacc0975395d444370db0019e6d22c7405e4239b01d30c7fa0a42d0d77057c410",
      "4308c1fdead1743da6b351d5fdfdd7dc3131455ae3338684170e858f240b2f2c0c2f2c0d1f581c5e781e5e98205ec9223e2e3781f524a030815c6e35",
      "265b57f4b5ce281a32581c5d86544c4d7b7b373dd41bd6febeabcd44c01844a3ec7dc34c4ef1f6c9fdb9cd37244700e0246da839da87c4c8a95b4b51",
      "9a7ebed74894ff7e8a37c8e362a2c312ca3dfa1fa6a71dd4789515b05bc887a36fa1b8bc25f2a7edef728c3399aa4fc9c1574b779178b60e69ab599c",
      "7cfeb98dafed20ffee879ffef9db2ffff49d7bc8271cafcf24b02f2cb1f236c5611e321c90eabea9b3467f4090c89f2b387b89c568b2f1c5af554cfc",
      "d7d7e260c364264bb338aa27f45b41755d0f7130daa0f3ae2806e3932e7f26a998d1bd0780f92f6241bf3cc798b20297c382be8ec293cfc87edd43fc",
      "ec6dbed2dd1e28338238e6638fb404a06b53b8fb2f2b7f1dc1659a5e24d37cbea29dac627c849c9f2ef7454dd89f7efcf5873ffdee3ffef8d39ffff4",
      "f33f7fe71bffb42bf4639be5c746ccdfecd47c6fe68483a54a9e3916fe5b478747838372dca7dd9b0cf18d5c2567522470ab68bd13535c96e67ede6d",
      "fd65c3676e94e69b6e7958341b8b48869de4baa1a5d0f40453dec98d90d46cfcc9f1e20d1019cce620c87c5541aafdb29affeef4a1f33ea3ebbcc39a",
      "69e6fd6854cfdd94fa71816b47339d7b5ed12e6c9a78ed2fb7a6a600feeee62de35153f70171dfd73bd015dc49c3015d9fcc388dde4048f9db20e21b",
      "691b7d11d387b303a2b56d78e14d1a4516df8e3f93d51d9685fd6e595ca38c5c390dd7c4a95d0877778e67b7f63ab50529aa28d5f5855eb791459f3e",
      "23f90b4b15ea90fa77ea49adc42baa3f8011a44a591cc224e1f3947677fbcb0d9f2012f11e61719976afc24c836e3dc1e1c85b1361366b05fba59112",
      "2bd6633419d0f50e587544a11e2f690e491066935bd2096e52927ba7b8b5339b2afc3735d5cf309c46f1e5a242b7e77d1396484425db4e2bef9a621b",
      "0acae958a646ee7e9f9ef4fb3ee71fbe822f86117ffcce1671406ae97239c778a71280e4eb726986dad8343ef867694f4cd9ed5242b241f5309a5a45",
      "64a2e0718dca20908d1395a72c989a731c866b5fae0647d95dbd89d1b562d986bc0bc38e1f3e2c452313c86259427046b831dec2384cc937a530a72c",
      "d716eb6b13f673072cd8b9c495c61573f3cc25e6bcb2fed73b6718ac2803a2adebae04b0c971eba95f9f37274e333a04309ca3b6a9d985da66d3190a",
      "bd570cd6d3287fadc5e95f7ff8e50fdf7efc6fbf1932ff9b24b542a5ab5e3204ef93dca2b99957915b4bf55ec555be017d13141c9181c12777286961",
      "9f6a4245131f818b6daad744b1956a23101f359baa4c01f22f2f1762bc2f80225b26678e8e493ed56d86d4a057a61e391c75650778a81bec5e4dc962",
      "2aa84dbe5af99064d40008c53e9465b9eebd876555d31b3b85859a1b9038bcd360ec557ba58b263b0c5740e15acb5af1a057c49215f9381fa2cea76e",
      "4bd2f7ac58b9fb32ce3bd87576f0662c8b1cf0163497ba30b3521529ffa1f28320c39ad00409ddc3d0bdaf116d0f179096215e94f8e656c47528b140",
      "2598e4a4b8fc2c4ec3e9c25a3ed7d923e85f9dcfd1fde7e1bc508b271a52f0222ccb3171723b1c9a7e80e7b9e51d5d8755f3e918802e1abea9c159e6",
      "e3e026bbab57d5f0bd1716581a74f57ae570abc36d69b09afd69e221a7cbce01b33a7f56630e7fb88faf0a93fcfccb9fbf8328a20a322bcc433d42f4",
      "531888ee50564da74541648dd944e9b489bea5981f5424ba666771c9472d93d02002f70c72dbaefc0bbda6678c69291411e9995d961d5679aabdbeae",
      "df78a3af6ffde3cc789d3daf33ec7516be5205bfd209bf520e13f0d7c8e394a0c66ebd8aeff272c6624405a8a7d9d86c7e0459c58438b91cc2e387d7",
      "f2355ab0bffcf8d71ffff6b71f7ffae15fbe375547e994a239f0f913c0728b577024c8ca38d546a1b430059a80320e5b5a165278dc6b92d7fcc96ad4",
      "3b19cea7cb9fe0a598a28102d0a9cd2e12f88451d15dae920f3f402e3e9fa6a406e092e825700e92bdb836359de882a360d6d744a195a59603d2d830",
      "0664b1fadc7b3759a5b79442cea508fef94cd729946ff46aa8719773539581d05a7020210f8ff2de67aaa23b159f76fa951a8e031992a850c323703d",
      "df8480e55832b05376f176d002ea9ba6ca125d638470d7808e7bab8dbaa3ae07a12177df8f013c4594709aaf0afc243800112bf2d5b2f193eebc1f43",
      "6a9d88ea3eb41076bc5195a8a701ac208603e1ff3a05e62b78e1b10e1ccadaa3602dd7be97aa37ad1ce2beccefad6aeded9effd0bb96a11e3c8e1a45",
      "e71bbcf1f5349f91341158f5b49e5d7b5f079119cc1f2854509827ad746f27bb16652ebcd7191b81c11ec033c35df2fbc0eaf4f6c03aa3051c45f3e7",
      "f7196f5970cbebe79830e3d4baf31448b406801fbe5c22b34d51e3f41bd635985463b5fede5958df9403886e4be19e3a73b88e7dc05b031039b987de",
      "ea1d7c712d5319c2f23aef4dbfeff7f4febd4192a256c7765e9cee341bfb7c66329e1599dc7d7f2fa03adeb23cd23b52849c2869ec4880823bb85c0e",
      "5265a51d4801bdd46f24cc1cfdf8527cbbed5eeb56153a571d55ef2d869bc08b14d9698ba639b16bed29070bd3d8022aabca3ec097e6a60d77c8cdfb",
      "6c5d5aa042c19869469f1679fb7e47480c26084582107e6606536fd43db4773cc4d051e35c1d2cc7d251abcfbb094f416d4c0768b601bbe16797deee",
      "d59caa87be36d2735e0770c738ec7bcd275195279ab14c072f7ba3e64e6e3bdd67c9a215cd62f0e75ad75d754c29e075426ef056fbde77f81b5d496e",
      "103f08934aafb0a74fbe0345ae7d9c4993ef2da654d1292bcc337a8ff3ea7db6254a0a086ee57e2684405e59489633cc4d377fafe5f1c632fc6432d5",
      "f312334833d58bfafd95d5a90f61463b09bf98ae7d5fa9e0fbdd0fef1e214cf9725be76e49b3702ea8a872168919777b7d4b07c1974860f0410600ad",
      "a85d435bcf9efd20ffcaf4afd65089472631651d151c7d5ea6784c61ccc3d2facedeebf6bd8996aa4056a9ef6c6cb03ff8fec7b14c0dcedcf725b0d0",
      "8596d5e379bf291009fa34fa79c0a767fea242f4b7bf7efbc7bf7c3fb54522d0ea6c472d0a38298902e06d9c9a3fe01114b7d80b53642cb3e564aeed",
      "2eaa6fa05232aa90f75a569252027b61097a8744623367b6bee8664ba8db2149066cb91f44f2c28565a10943a32946cbcc388754244480b6ce730c4d",
      "28f295f71d6e736d6225327b18aec067be695d3d40e7d2704dae199c0f6ca768a564873f5741604ce485356eb803c9d48ccef1e075a7c12f8db4ee68",
      "1c423a82c71a430affce1b1b7ddc0e1b07ae3035a3734ca6e6431ab19c17716d7df5edaea6f2f1659dad4fbabd0cd081f34d402358aef30f47c3b510",
      "755b0ebb4327b1092c769f5e7eba60ca516a110900ee3128a2e06ba0e4cb62851fdf035dbd3a0bfbf5617886ab7acbae4192e440b7e7f551bef1f5b7",
      "030acb7c62e9de94afaa5cc3324d122a6c57e5c39422e2018001a2667db99f55f181f21de2aea5899d8d9adc3828719ad8c9c8078e91fb019d85bfb5",
      "e39e2fdf723018d13b288679b080ef55f3c6ec7fd3e671507283efa4a9961716176bc1a51f5e8d8a665d53ededa0b025174da071c6ddb2bd92ac0a0c",
      "36cac76f8a4cefc0733aec2584bc8e7753f495ee12103da82eaf172885a7834d1f2c5200d9e7f16e24f1a232e2f164b2f45d69cfbdc327b3ee270562",
      "9a3d467a055408cf7546a64c7547213f4fc22f7285fff2f79f7ef8f1eb588f4f46eeb3d1fdf8065edfd2676ff2f56573d54f5ed8d73aaffefaf37fff",
      "e1bb98255fd1160db6ff6e69fb23569009dc9ae84be7517288ee8475cafed7aeeee6840043e24f275d67b58381f8a8d5f12a7af09930c28b76c24779",
      "85570986cf641a3e2839bc8a3d7c822a79059e3c0dd7975ed2fff6e34fdf952887a086dcdc20c0b0e8189d95d78209cfd538524e62841ed46a56684e",
      "09393e40df1ca1aa6b7d8de47a9cd4c1c78e26d2331f02cd2796e254816875f829f9d3cb3a51fdfa82df8600ecf32847f177b597ac1f6243ceae3c82",
      "17bcae99845dfe9ff91b2688eb7d4257799720bfd6080cd76718626d0bee36d1578e506e2cb44c0e758887c67425d705b5a613b084318d547b3ed42b",
      "c2d7efc2cd8f505b22fb3378de1c229bc5d883216a9475541293ea01db1abc0ac7af994291c2a55ae1b6540de8ee6422ec49aca32023c651e583a4e4",
      "b2d4d0fe43a7ef0af68634567d335bb50816c671a48cdab5a78c2d105988eae9a190104a311455950f50393914dbaf599f6571fbbc951ea7aedc94db",
      "519281a1d827452635b1addf5a6ea506640aa1db046d64453130e6be9fe955fce374fdf242f6ed9fbe4f3fd58aa53427f7234cae163c683f7b8ccbf5",
      "ffd4f872ed65f1fa17b5e416acaae4b20857c9d43bb908856d34cb74e7cd09bfa484ca3947aa6982d080bb2d95906a1a234bf3f872d5e31860d56aba",
      "89623ec97a0904a7755489861ae244c918d2417d0a2e77051f81c96a74487148641fe0a3dfbacc1548f52380e44e9ef8ffea40e3ad04888b8c0fa312",
      "d48896405dbf0fa8d822b66c40497c39d0cfd5967df47f44b54ded231e244f1bae18e428044913fafc454983cea5bbd69250801846a8266dc874b044",
      "3fbd70bf53cd61c166a249e8becc19fa5362d5d6493964bd967bccaa59de4864a4689a9a29042c9fa6c9d77a7afff5873ffee57b36d917116fd0bf96",
      "e35eed2eddb4868bcb07a197a7cbcdd75c8d2d7692b0eee90e48b1f2808e09fa4ac11d547f360657de9bfbe722c65ac10c4f32f3883bbdbe8e8f2fec",
      "f975bebeeee709f13a653e9b56cfd3ee7952becc5a8880b67bdd2287334c438da51ef0d7b54ae94e53c4432da8e985236c07deb5244f7ee8423f132c",
      "ff286afea27bfe9934faab7cfa4789f55719f657a9f64fd4dc3f0abe3f4fa72f4de17fffcb8fbffef2c3771496a490d4dfc476a59925a2fa5ab02c4b",
      "bf4239db494341816479487983ca546d175701d121c52766324564c285b63735bc5a456729c4c3b2edca5409570c2eb34200949efda17ce436d05f68",
      "57929d10542a73621cadb270cd557b095976d279f7125021e5cb7b716d49f05f25fc92e978a8ea903710f24adba40554955769867e0a049a65e9edd4",
      "d4129805c2d8e10223680da7129bbb7954dca6731e9240c1faa924e1e5487f7c452529a9c31fa2b41d0f5e5bd231dd9f2325e8ed8ce48cc1215b5e94",
      "86ed1683d4180f651fd190b3059612a76a6bbccfa737fca569f51fbefdedbb2102e2775a4dcd06f9508f2691720a2d19c4290aeaa5b3909b9736752d",
      "6ec92a1c3e25e27cb0c65d6a428703aba35c4b10b30ecd169df374658c83399f768d3145e215721320c8f5de33251a13b577acb60f134852ddb8dd09",
      "f4e8fb6bc92d1fd3415e05f54e7b261ea01761f5e90e5f2b9633983ef8fbdef95a493cc2eb11c524da00d1de52e7ec3e1dc08574144745c25d1e3cd5",
      "6eea2ff796d805604ddd2cf76ee20634abbada7388eaf87a214d49ed5c6be17dd676afb2a993b7a909724629cab42592016fb5de69e0d2d2a04f43f4",
      "b2a87cdcf8d22373ce7a71ed7aef6ad26278e14317fb4b70e5b13b4d7587ae48f10b903e0ca3b8139d5db18e656fab6b7f5579f366f882d91da9cba6",
      "6ea8f73e2c62306ff21f4376442296438e7be7637a9a9e077b1cf42e4dc8a1c38e75392d4ba0eec30ed64510aaea7aecc860c0e98a5587f23a9c5f59",
      "ad0959222087c887ea1085e62867aa066e10cce50d1e7a1e5a0e2cbeba8f6949a0a5f0e087f3ad8a1d81fda71e3eb62aea0c729bf7cb0f0bd4cf590e",
      "19db1e1110f7e05d534f9e70eba1792572326ddef090b5f3884f6313435695a82dfea6825f88724dd317940e4f2115607d7857dc3ac2727d38d96db6",
      "f10191a5583a2f0761365352ea28ef84a64c5c622ea6798abef34a7bbea6792fc752c82ceabf7bdc29109b7abb66f40d70efa3aa016c472df7e3e37c",
      "b1b0ffed976f7ffefbb77f7c1d37d7c16e0f93a165f718dfbadaa2a7594ffdad62db755a182499baf8065673891d55a43fc7d60060e9383745dcc0a1",
      "8ad055a1e40ed3c8c93222c3cd273723b1cd85667a464c2d058bec7e3ff58a8bb6e9d1cc427303f449338908686ba0232560eec30c7beb21698ea99b",
      "3a115b134db3c845e319da21bd65c9d9e71726f18048d44ac031292d408e4473ffb2f20ecfbf6033a709b4059876aaefdbc799118771a71b50b62bd6",
      "3a633ca44d41f49ceeb6396739f757e3dc95e05525eeb6d77c239b5162504ca385e9fab23d9ecb7c7837521322b754696e87094a4bf3e3f6c040d3cc",
      "648adb9e43b046649dd594b4356c0402d5ec7ab5c544d982e1ca74ded83a14a5f4420722f9e3bcfbd27cffbf7ff8e71fbede7807891e144409a89e5a",
      "c90b8e61f5de03fb06bce059adfcd9c754e832455ba9dc1178e07db9cb2a66af6381bb835a857b8771a5c10b498fb60774b0cf263f2bbaa849094057",
      "9245f01d70be422c9145f91550c65d5d004f56b805dc489ca328a3040e6fc2b02b4e6643c5ba349fed86ce03cd2b510f776b615345b8bbd23dc2829c",
      "bc0ae4661201f5b4f0e0de7d555cdf26cff45ad267c89bcfee42a15d4c519deb9c3a02ec55a02be0b78fa63adc287254ed2ec9299e7e4c73bbc3bb03",
      "9d5c96dc6eb8c0971ba1527671d72ac7248db1919738d2382f6eedda7144a97e5539cfbbb32d7042753d5adcf65478d51d1b964c1e80f294f7536d6a",
      "49961e8f347756d0e117e7fb56cce14d4c10b2215c44f30ec980eb84f893f2bb0f28774bc431cb7f4f66ea70434d93e6ab0ecf592bb4c956547dcf66",
      "dc115e75a8b23e7cccfb64ff5aceecdbdffff4e3effe975fbefde17b10838404ea39af11a141ddb9d4799d0dbca8c2e0152cdd59533c37b07c7896c3",
      "16f1f00dbaaa9373fb74470d3093127f9b118540b5056f8d41d3554ad8459dd746b212cd24efbcc944edb85d80e7f0a24700eee4f1aec0e3b9733bc7",
      "1730d3d67d6593dd90a8a7c2501425c59ba1b288c1f4c964d792dbe5bbb5e4c89d91e1c6d20fb897c6baa2ca5780840b24a085176d08a4724c1b8b57",
      "2220fe948a87303db614e232d2db330c59edfd5007b538a2d2f5ce66b26c1002afa8e6119601cfaf9a8077609a17f1323b4b97d7a9e4013fdd4e8331",
      "744304b6e0a9ca59f210e37b13264c41227c4a622868383f25e6235b391ce3769fba846b44ced1dc04d1458658836ea2da3298fb7813718aa8e69e26",
      "f867d3ff934fe4f5337afaca3e7e82af1fe987aff8f5337f5d0a3e2e16af8bc9eb72f371497a5db68ad2ec58ae0835964c541cb98c9125ca4624a178",
      "da6c077eda9b1e201e10c4004bb1982d2243b02bbf047fee08c86f2100d426dd2267b04818157d893bbe2934aca45091c392587a293898daa10acea6",
      "177b4d94198919a0f21413cce12581a52a2a9834cf374ecba6a6d53ab9874c6c068ea59fafbe24b7b62783afc4ce9387e77d2c306221021ce9dd025e",
      "6d72a6402a9dca25a26c49fc14d5e57ec5f8a998e64291a8c28462c68af8719cc2299e8b28f8c0ec09252b74b97dc6318f6be2d73a267efaf5db2f7f",
      "fcf5c73f7e2f74e6d696bff01bda4b5682e4cbe8920a93499eb1145b701f4529118d67dd63fd741cd9ae64dbe55d1c1b2cd15349f011962b5468b275",
      "57b46f112ebedd96d2503e93e67968cb97b82047dc097c17b2752f64dc310222b285d2c3cd04cb9723bbc4fd659f12a0874fdd2342996616476cc1eb",
      "7c53f3dd92a984c39f9db5e686b7cbf4fcc015b457d9865bcb6f4a85c3f30b83c78a3eb71d423966ec286123ef7a33f5c78867f8e4453c99a28b90a6",
      "9cedbb33abe40d9144006cb2edf0dab31b423476bb1c5bae7e153bc42b11d3633dd75394753cf42e665a5e458b504e9baf0725e47938374bd0604928",
      "a79ec0a0f7781523a437a11c5ebede3633c98d66238950703de7049e6de7379e1f368600be1096c5db0421ee00ec376c050e196cd79c7643274cd950",
      "944fb5d5cb44e8c2182435f91d13e55c6c3529e697ad41963ef6a132d1456863c4966af7b99105bc0df18ff939208281a07ac06b96c3d4441f8e3254",
      "a9356c49ea6913bca6df10847bdda64057626bd2769ba42e6312c01854fdfccdd8226cdd6f17b8d0ba5f6f31e3a7889acedd29037a934e97e579a700",
      "735bb64050e87feb4d4eb096d6799e3845cbe75e0375c63d49da66d2a7e60f72a2e8e35f886f0f66995e6dbb0218ff420f8591091e2a7ffe66631608",
      "8be7ef887bb086f06f1845a6e0717059fd26fe28182775a0dfdc06384d708a329bf7ecdba69aa4f9c25d862c2c68f9e953b7bb5721654afa09956f4a",
      "d8b6941da6aafa3b7e22effb10fa301289e5e9286e4e88518fbe836eb3669d413408f593e73da343e604564fc66ff478e8df3236ada9534a15fe31b6",
      "6b05ea08095a315790617adce26cdf00a5e5a2058ee8cecb8a35a58a8dab46dba48ad9faa860a65eebfcce809141c668497cbdcdfb43c9ff8a1f1b39",
      "2ef1a971507cf5fdc3b16523c94528e71fef552991eca75aae8ec8a88e9de707e86ebcc778cdeffcf9af9f7b1d81bf63fdee3b5ed76f1973eade364b",
      "2c05d25c95c6d5647bf3328cd2a908c3300ec3366e888bb615fb696acea9613471f8a4b2eaed996b660962b2dcd7f838060c13a49f60998af58f5b4a",
      "0403fcd02ef1ae3042eeaabdec379e27d2af2b5dab76bb77b5f0e74eee947dfee7fb41b7aced1505168c25df8dd7d84c1971246d44930846f5d5a422",
      "788cb9a530e6194628119a4f4fa616bb49a55252895cd3aeb88c97f7b0c3784e27bd1eb714e09a1b1287191e92e723c162224bc00fc5dcd0c24aee0c",
      "86b1155f168bfcf3a9acd3f89237519379dd97d6a0c40096b291f15a9af968315eb755ce5bddd2525eaa3bf9c574704eb8718b1ee1e8d885b460a391",
      "082b9ad33637d40eae51cf3206fdd30dc9826b0fb3d11587dbbb2007973441f44d49fa24e7035281146dc04ad1d4e77b6fa2d98f07aab99e1f17478e",
      "1e7dc7b07147ca562e25e3fd94484dace1b1bc4282798c2d6d1f694e47d9122d5cb7458ac73b2122120352d49b98bdcf252999b2c6ad72e4592e2388",
      "49129910b2473692812b92dd5d77b0fc0d95dd55251f9c0f55b34ebf4661d97846ae4bd5f3b84ae7d741572dff5028126bdc25c4b544b1e11f870257",
      "5e354e6a7cc9202673f87cd7c9fb18a7c416f9904ab4a1aba7b7d8ed1b3dee9294a3d73ce28e153744ad67791d293bbe34e810765cf324c0190e0bc9",
      "4ee1aeea31eed967ac4d399e5caa60d5abd83ca3498494237428fd3e61ab047ca5bebacb39dd0a42b0b9596c49e3bea7dc6c8a7931e5de6f691e6094",
      "6e5e14d2182f9f689c97b18b2494278985799e3ddc2772f067769160ca0e7bb26963349e535aa922988b95b1899ec7c62bae9be7755c6b4b93f18a1a",
      "d6dd28c160c8f84cadedb7b92531dd0f5a4346c15c31d6158b546b462f612cef366fabe0ff628ca428a4df063c9cc3489f007313376bc597410c2f37",
      "a89f96b2dbbd102a1fcd4c5f37cb5f0c6fdf50a0d838a23e8238728e898440456e3ebd9e77818045cd3ebd0679892a644d73a66913b2912062960fa7",
      "133c5074e1c898ef24097b0b5b8f092bf98938d9c03b39036b27df5073e3b1767338179e6f5d0929b93d0f3f437533fc9666753cf92de1e33e9c6d2e",
      "e1f17436f913a93ae301af5835c9bb8c1d5e4f59be21a9730d3f4e69b1644b1233fb92c930bc9b641c9b1e9cedd7f42862ae7db96431c9bf4663e563",
      "bbbe6edfd0ccb118d1dc3a727e1ea19c9c96d191964cbd41453b73f31bdf393cc304186c7bae66e7596e62d86fdbb33a9fcf3c49e9ce3fb45becc7b4",
      "c348cf5b6c93ab1fcf364934f3e97441526dace114d176228560f828677af7a0e3d67b3d2e3f2c6b653ebf0ada615a4cac80cce986942db5a7e73483",
      "eeb2e8f3a907f226e358d9f1d4c8f37eeb5b725108243b4faa5f1f1272c453382e4f52b9fcf9c84dc66c3a96d8c66f2a6e38ae7749fd8426236e2782",
      "407416a49789df6c0cef6d4b2ccdeea9736bb70d7149898bb931eeb6e505ca940d5abf6d78605278cb514985cc75ad19b1821daf05418523a46eccfe",
      "6d817af58f264bcaded67d27b8b746de1660f7e6474851215df251cad30dcb7df283ba7ac38eedfd7caa85d2bf909d15bf59bc5d2642faf9fcf4ca92",
      "745f2e7b75d18637f6d32fa0285eba6dd730d8c6921f973bbf3ad51ce6e3ae6b9c60779d80c2c19b6f29029a11c1347dd04fa34bd42b75225e8c31cf",
      "4a27c4435d6bea89eb778d5b3bc30b8cbb0fc7dc33185a8665f722d88c0c8da8e81cce1acc80a9f4087b875590b1f5d11ce296162986796635f2a3f1",
      "13d42595f14b0894e5939e5891b2c86ef4502a029e2a6542f23db5d15031c056cd14aedf986c323c63afe7fefa42024c17c9e7fe2e876dc9861c6cd8",
      "067965fe3183968a7fd28ba1b108969afbf30ec8a43cd7af7f8a7d817f44fb28ff140d8eb25039525475f2b10ff1e6b6486f75cd6b98775ce5c19661",
      "1ac59667a4dac4fdc82fdcb5cdb12530a3d8aed77d52684969302eec86173277a26ac2564ef26dcaf125f967248c9274f47d73d96b2ef4536ea7e680",
      "2d946ca9ad772ab5d8420344980eda28b095480be1bb95eae3fac9a14e8579b2458add580576b86186f5f8dda4e855bc492b4c5322da43ec60e59e7f",
      "cc7e8c936723f0ddd383d2dd4da241513034d8899cb61bf4412a1f4730730e1bd24f1e0885dcf10b5d0a5dbcb3c08d0cd87d658322e19eba5c891505",
      "f2a2331b8f7b32ee6d44bafbdd7c6eac9dbc5169b9f17199219e5fa3ccf8e6d91bb009aa551a9479eae89a52da00d45fe3859d875c9ac8123fda712b",
      "5b1b341fc67d4665c9cf0d2dc00799214d54928a395606114676d9aea13b388f32b56526886a0f8263a5b0cdd83c4636bc53e9c075b2c8757667b552",
      "ecff83ffeb9d2d901032a5ea5d61ce788b20ce537dda6440a194b3a4ee5e9e7fe27d4b26c895d0e783a7c85c95fe299b7cac5b0f0ee52efbbc9e2af4",
      "cf141ce0249b5738dc90ca9e9571c76133e297fb12451a2105151eacb5845e5af96bacd14338a4492dc9ee0fb94fca97b21dba416f11cfc791931a32",
      "4d337c2a8469a93ddbe8929c0e45523bc9d2bae574491221c7ce5356546aba21c3b1e764d9908288e3e694dfb14e6fb3c2144721eb8e397adf05dfe3",
      "376c4e89924bdebe3d55e756943e22d50e847a8da8f3b4f318a356a7a5423d57cf01eb94d3fcb18356371a616b799cdd9104a9e3e01d848d90cbd7f2",
      "1459b325484492c03ac7bca05db439d86eb39dafbb741f764574f7f9ae12d17b42a06f3588cb565d3d65c513b7bf520ceee8c7568a32452ed1c66a99",
      "d28ae382f84fa4997a0abecb408ae1ae6da72c4a84e344e3024f712fe95e4c9a233be37039702e26056bae6c29aa64ca0e39df70dd9fcb49ac6fcbcf",
      "914c80834d2aef4b4ab866fc03d61c3f71005a805fd1c57acc545023a5fead3c87f13ed6598c42dc76c3d84dd8e0e40adb8ccb258bf6b886e7512e92",
      "445b4ff90cb4097a3b052b2f0dc8634bbb5b19cc20b57bb09520d643eda4b4f594cba877571d5d04bf350abba0203ed8a26447cdd2cb34e8e2e9298f",
      "2c5d982e8fa23f7d2ee041a5f448a1340219c91b166726ae05bec7106c9701c06cedf3185ef51f9215681448fef261110139a87d4fdfabdf4639bab8",
      "73e85b928d0ed5e6e34a7c92e0249d00a6c45d0e2b619380f04346423232fac43bf986431478adf74e298c58d1a5a7531d8eb518972c09c3e7b04dca",
      "b2e33990e48b3c91f689aad5c23d9f633198027a7df2a625185f56e4dc2306a597793c3bcaf924e388a4dc4c2f39eeecbd6ae6c87696007ee227a7c8",
      "5ce5533b7bdfe6e0ef285147ea275b718f2688a03d7afbf7a7ea568b97243ca4ec13879b4a6f80404654ba8ce7877d76cefe74a2baae5c87dc35e09d",
      "ec4ede912d420a7d6b1d2d4a7f25385fab9c484ce6688278f338f427f6857eddbf974d76754bb491db540c86b88971743c4607e9dc94342c0559a6e5",
      "2a8aaace3d503f5026241d766fa5eff363f4008220e289617f7dc45af61e4d4c2906c952fc56a7d45de4ce2c7393136f4c5f391b53bb71586d493907",
      "c572d102c3a50f0b879a6c7d472101a8829e4d2324689913ad3ede8168cb8aabbe7b9668ae92fe323f684269d1178df8c190f4ea4a793c1d45967fd8",
      "61dba758220d1c3fe07477042b61d48ae39901d8cc189878c422beb8281ebb0391ca407dbc05e5b2942655f75f39fdbac373ac97144dc7abe5f63858",
      "d4bee6ac4f26113bcdc7b705f665c76baf0610838febc9972f562d539f778e87aef26f5810673cf40c70200d0ac921e3304ca22a49921d31ba004d0f",
      "0cc84e55dbdd8cc522ac722ff1a2bf07e9682d7bd9d704af06f0d4149737e00c49ae614331a005f490febedce1832a4cddbfdf5da0aacb4d8fba7107",
      "931232999ea121b7a0e6aae4a37280038586f05dee5aa38b6cfb35af2312218affe550a08d0323cbc271107bfbb39516fb34e4e0babf80ec4e09de0f",
      "22bb71e737143001c0790ec8708370e32164064dd39a4d751c8ce2465f0b53720401ac45d573bc664f3f3a0da12b7c08ab555a52c807476c0de49f20",
      "6f0aaaddda446d264744155a6072889ee26cca324a9760d26606f8df41319866c52b7dc0682acb1541ca4244b6f2d3cf4dbb7d03b0a4bc05c0e64ab8",
      "f17c0610e3e135e51ca49828a4e336fc4db9072331086894d26344ada341a8eec4c43dc683ac4062f4bc1d679c14baa5fef87608e4b56f3e5c1edc65",
      "f17d8da2a05c4f16a1e6b43237a40b8ea169e156fb00ffda39024daf6ceac67f32908851666089a8d8a6a155ec3da0640af7eeb93ccc93c257b5236d",
      "135faf9ea24438e9ec243fd4f67e5c68c63ee9e4a4fe2c5f6b97edad243979c644b7a199d6fafd4a8261c5c3cc1cc88a6cdd2256c6b19d872d6ea2e3",
      "9aa505c0c4f049366eb59fbcefba6038d2a88fbbee548dd05b1e94b916258259ff313826035a63a9b32a1174ef6d3f2edcf26d62a58bd501ef492e18",
      "267974f8767224df83765cc516976ae32847e5da7c547277ca94f442ac17fb4853df1787ea0d1c18d7aece1c24b783e09fd7eeb9b62cec4b5c24b00a",
      "f90a3361ac1265baf795601533adeab370d265d1531ad7aa5edf574344b03c1dd5e9e3e98f4bc802511bb33bb2260b38ea7c5c4116b1e6f00a52bc70",
      "2e688bf5d4e3a090900f2b61dae6e05fcc5792d7ef57dfece9296ebeac10348be4d07d19db34132cdf69b072ec4aafce787c9eadde86f674a25aa5e3",
      "17bdffa2842c40116363a688cdaa134f5dc2598172618f27d3d2daf47479bed7551e8f52f56035672db765c9249abc6b7bcc64527ae8539fff38002f",
      "6a14abea46fabd6b1112b2a9dd06e73cbd738c699db82faa5253f6180dc9a21d3de595bdaeeee0c1495335f0c7f72229437d498f367c84351e726ad2",
      "1297b23ccbf4de8719cbf0b8f7149d98b09aaf160a0d22f41ac9b680a46531bbb7f9b859408864381bc7790514b39277529a7ac661c01a76311e6e45",
      "fff4c88f203805c5adbcc8c3db80de00e8fbc3a410cd622d5e8c4b8853422731f3e3a41332aac7a3c5d72082c4128fd6c79d8cb6945c1fbf07315a96",
      "62db28777e41a9d328d79ddf2962d3c7ebcd93947eb2390e7d785c6253bb137cc8c1e8b823cdf230f96a3a1fc1fb5bab78017b3f0d33dd8eaa6a3cda",
      "aa99cef5f9876412a8ea3edb7b5a5eadd1f627a5c1fede2e9d9c1daed61b525ff09ef60f8ec233eddb4670c270198da2ac35c32fb2997919dbc83e6c",
      "1ed5e269c1edc7995cb5397a90d7eea73d78a7f5fc22e9bf4f3b6ca7edbc2893f26c0372757e37e8e9a8f7cfa70c34d09cb6e33d9ea907d260f7f1b8",
      "dc65692dedf9e80272dc407c469b4148904a733b30acd951a6a02975db6108dd389a4a6af1613dc4c7ba36127bacd78265a22d51759780c47a0dc1c5",
      "ab2b6c3b84b7c85eb6f17cbd2ea6ce47075f5013d5ba71a7a32557523cf9293610d044e059816d83d042094a3be7f13177ed15f331430ed26ea5e1ed",
      "bd045d06d94165d2899442364d1d93c37eceb41f00286eb9fc797716406e6d781865f3969b59b3728be3ecb965eda2274c8ea0542e518f603a6850e4",
      "1854fb2825985e04113e51bf039fecffbd1e3d9eec647e79748c74ee08006b3360045b7606e490ac2a3a8b6c7d69eddc4951d9e95edbccaae2f5807f",
      "66afacda02a7f33dd70b3d03707d04eb019e90a99128fc5086ba9ea1ebcb10956b061f5ba9cda097a3933e241ce904fe0ada4a4a57ce55dc135ed8da",
      "8c6419dcb9c7b6f67e4cbe6579e4d3b64347ad025c7532ef5a41ceab9d86f4d05454ce938d14b0ff9a63aa20f31530cfd94d4749817057e7d7770a96",
      "a57b40be7f6fed03669efacd1f5058ccd06da4cb164558d82684a04aa5a77466723db0fa20a654296e452b44331bacb0a14ae7168bfc60196a91263d",
      "bc63b20bd0816738ca1938e8cd9c1e0f625915e4ba53dcc1367c1f8fa5c4508f69bd544b94dcea38b33a29d3ccbbf4178143bfaa7199d7255ad8e428",
      "bea328016516e19085fef4b248f672f5382eae861e46403fdbf48d3ce250bfd681f35f7effbbfff51f7ffbe5effff2f55ec8228e41718534a73616fc",
      "bcfadbcd4c22acc770ed86d3798d016b4215a2b338add1a9338bf02247ee63c36f5cfb29db6301a8c486186d49c4299017a86e1a675992552eecb6a5",
      "5341d0d6bae2d725da8bc5dd1e58880db9652bff908c2188a8fd701f4332b24592923df778aabc7d3b467c3c8fc4d748b2fe7f1bfb8ff7f472db9f3d",
      "d9ebd3bf8ed0eb28be8ef4ebdbf8f8c65a74aa42dfbdce30429aa019e0ae3abe4d93c144a15f38ba317d879e2fdf3fd4fff9e75f7efee31f7ffe3a83",
      "44b118a9062ef09f26bb01671b1839f0842288d91634b9d142ba68c35d8162850a7e8855a6077450810c9662b58c1b59a4a1639255736f9292a06177",
      "96a3b7362052adca64c72103bf09f6992391b588c5d4376c8c025f52c9ee3636e292aee5a49ebd6931eddb34b32380c4b94fc61f1fadacc3488a05c4",
      "127dbd3d004bce3f3f59e66eeeeb75bee7b659cc30a42827f1a3d77a0ba4f1c08b93241f3e983a1d42852d2710b2903c5f711c47ad409f6587948a02",
      "4e8f0399e9852932d4987ab4fccea88354a797c81cef40c5128384c67bc0f268874c5270bfdf159ebb148b7b4a47562a21f22da4d93c16a97ce7c39b",
      "7813ed20a2ed5419ee05a50ed81ed981035fa62099a5f53ee63e2756c87be8293e6d3fec088c02ccdb310a1d8a1fc91af77d4a4909a40a3dcd070b49",
      "45d58aef398ab1e21f486ea0aee31cd5d4530da968a00645b9a5011df92069a984e82d8f761fd10cca5353e1fe26288aa89db41e7c29803f758f5fde",
      "783405c13a672df38348041530dc85ee83b6b970e16df53e7b3b2c4db4afcd33cf6811ae22b5cfb654300d958a818189e4cea7899ca22f50b72b569a",
      "71204db4c888c62947370f8e57abb1b18560e2365f53c8d4c903e26f77c0a318e7c5b79ab45bdf6db2c57de837313ef8a2f308fe3dad335f63aafbf3",
      "3ffef6ebf7f0d42dfcb4124bebf6bf18f8305e0371ff2f8ca772378a38038af703beaa94bcbab41664a95b245434dfbbce22a20b2d644c4b1763a476",
      "5d45c5544c9fda241b0686df94a4725a0716a41882c9173d38917af520cb55f956e45bbef5019115f3223ab3491eea9436dde52c6e340bd83b64a9d9",
      "12cb7cb0fb146e8617d5986ff029eca2cb9468ba2ed291f3bdb4d861279e35678dd31a4d455567395ff1a965c10c58848cc851eaa1f22806093d378d",
      "d91443cdcb10673532ae65059b1bd315ea14f21c6ec0aad5e5197da467ebeea48b0ac1cdd9f0731dfa44afebc783b6b89d659298aa8c4df70da67636",
      "e1a6b3a6f3fce2715bde083cbfb91f645758998d0be496298dd2c2e6cf9a8d1ede69e946980683fc46d17ab72d2b004f1a253a963b53658b134deb5f",
      "f0c493d8a0b803e38f19e31ee6f597bea2fff4e31ffef11d64ba0fdf88083d5c925cfc19ff943a0792740e72a061dede34f69140bfb6189105b4b483",
      "937b69ffc38308664ba4904d26910f1727794836072b1b53fcdd7005d43bb537df57324dc53ebabd82c7b0358cd038b8d6c92ace81ea7d872eb6a5dd",
      "c3224b5476b908c09feeeaf415126abf0e4f85d4eed06edd73b14cf5a0550a3fa0857e30082e71822c9f318889b413349fb186680586ff48261568ee",
      "e7dcee86871fc5f4f014aee5e4d4908900ea3afcf129d4978e95d852ea3caac3a3863b358ec6f06eb16704ab303a0ba58a103004e2c9062df111464c",
      "56ac80c057b04bf09ca765cb0af6391cc1a5ad671e61ed0211b5bcbdf3ca67f6fa71348ed860f3f0da55825a3695aca79e410107a067688d29214526",
      "ee258d442c898c04772399c175c8809396d6e926f222ed15dd8d7766b2b72c1fe2701bebb07ffb99da29bc930f940b5bf5fd92b765f5d0e2684e61c9",
      "8f8af7ed7afc43e45d663745cd880d63eba1f211c4fa643378dd303edb543e6e3c4ffbd2c71deb8b6bc1bffff52f3ffefcb7efe00822f38833bb6229",
      "a0e1893c1731550e42a8414dc14f14042843d20ea51f82a82c797b5244417f2ac9274dab7de85a89830838a23e3f075cd7bf3762091a13a0672c3e9a",
      "92e85f76c9b31c4a5a92176c69f8b2c15832508651e7451c81009a768d2030c920f9eb7d011fbead7644c10120e0f9bc1bba12c148acccf9f0379e91",
      "265e955824da33e96c0912dde8ee7c67bf44fe8acecb119b63c20527231ca1a6c0eb54ffebd997f1e866283d8bb4b3a154569de1aed50a6dc0417739",
      "d12a7f8f750eb8961014c302bdae9bb088d11510c45d5c8193e48976102fc18b02ad2d88cd78f2cb5310bf5872a69e676722108e58569e34fd4c62b5",
      "8d02d627ccc1afdcb3affcb4af1cb6af3cb72f5cb8af7cb99f70eabef2eebe72f3bef2f77ec2f12b69259f3fa8076449fe944334bc3d74bd1fea6120",
      "c52a97a56004aa202cb0b8ed5bf29a4dfa53979b9e0f8fd0b0186239cc42bd2f494df579be2ca45a241067260ff888366277dd7b689326f27a13014a",
      "88c36a37848f6b0429d6c307fea545e5dffdb71ffff0f3df7ffdbacccaa7acc91f99955fd997410a4e894dd69d0e3d1dfb079ec26c877d49da67a27b",
      "8eb348524b56ebce3df7c9c07e1cfc8f2fe8f5257ef6a23f4e86f789f2e191bfc6f2f6f33f7ffbeb8fdf25200256a8ad625531cb037f6e91a1daffc1",
      "5097840777c829cb2239cd59fbdd80d062eb3918bcc4a449e4eb7de0e3cef132f35ee6e6c7c9fb3abd3ffb043ebecd9717aefac29460a8d76ded5943",
      "12732b149295b85e56983b3ab743a4e2584a08dd5e33403f35f691ac45c456c0c5117c724a3066113b047d9cfa45b15850b109b9af5952e7518bfe37",
      "ded19766cbfff5e76ba27c071ba0444ae6db2d99610d8f46658f2cea09c79a93981e8b2996609c26314b42d0ea77445a4584a220d9428cc625f70c06",
      "f658e0e1ba413664e11512201c934e37251692b7e9f074e218c137904fa983dfc6479098e0b261d14e0ac6de1b49162bf15b7910b461dd7f83a36547",
      "382dfd5abc06ef9e2f861cae2e45f308afc5a2ab26cfd8a2879c581ce80381dcfaf3a48ae5c8b41349b62991b46e1e5df38763b81e3c9c8a64393893",
      "1a112ec3ffc4c38c138c6a60f2a91eeceb78f5f69d02c3fb7bfcd2ccf93ffec7f7cd9c975950f545c8727d2a31d2d7eb84a497ea6439ef8f0a2c7a0a",
      "e7053e5802ac0a712b966a8c4d519f01f4b7343fee93df019b058066c7dcc96a8c04477e724045e46860ffb5566341a9863b3444e1d379fa612e7f7c",
      "d2af6992fdfc2f3ffdf8ed77d710ffee7fffe1977ffdfdd7b544a8c28b153c5434c6b5e9497cda52d3701fd2a9c782d74e10b6c1d06ac13bda22d774",
      "6e3076a7282ba33f82bf4ec2359f60aef3f7154586648970b050283ab34126974cb5b4d2cdfb8ddc317703d7408ebb49f4f33f9c8572cd3575392bee",
      "0696cdeaebac23f5727d2ffe75f78e7c6a215093a186be0e482f286f8b2baa0613e829031bc146a5a1e901f1d952bb97c55425fc104eafa4baa31249",
      "fca8634a328ee5e1157cad32f9f32fbffee577fff9db1f7ff8d3cf3f7d9fc69149db95bf3a3258244c441cbfc3b241f035091cda3035e06a0c084bdf",
      "e691cfae5143d5d9adc49ebd4401749a973f207ec983644ffa7bb909ab04ce1656fa3522e2be669c2c3d3e3056ee24cbb247f7a96566dfcd6a47d867",
      "c2edd0ca014ea3e4c5aadedcbf17b1fc1c7a0bd7d21c32581b5e650627593cebc3707dcdddfae197efa1d3454c87e80f508f9b333fb13031071a72d9",
      "292d48b8603fc162e4de7d08e01eb1c01be9e0cb1501556425ea22508fa4f2f2ca47de6c98c1df991e4a3174fc36092a860416c27122ca2df1eaa7a4",
      "0d295ac6ec40929d53bcfb89eeddc200c518420ea1eb882fcbbd8b4c2980fbac315677651e2e7d47d1575c541e586fe2215af1d3adfa1df675a6b344",
      "81ca6990f964c6bfbee6cfa642bfd6eae7b34496d9ca698942be9aa89bcf787942e3b1570b2ecc38a481f3e78e9d57a4054b5f4511e38b2db959e920",
      "b91f013de8eef522bb5ac7abcbcd2b889515f8766ab72a44a490127d9e1e1d4b473029480f61c94778ac420af63e8290bf758fb27d4611bdfbb26739",
      "6b1296a87c40f5b796bc9755f175e17c5d5c3f2ebfff534b36007b6d0623ae0bc8530213e548493d7d2d5fac87fff4eb0f3ffdf0e75fbe5e12b718f7",
      "f45c5c476f2c65cfce58f729f3ae9832eb3c51f284cd268ec332bddac2e91e6337102ea391e43c22895fe6ebbacba579d6af0fdbcba7dbcfeb16f5e9",
      "3b7d79ef1f67c6f3d4799d5d9fcec0a751fad2dbf98f3fffcbcffffdbbdecc277b48d1875ac2fb4310bdaec7effde5757efa30af9fdceb67f9fae97e",
      "fcbc5f9780cf9689d7a5e465b5f964417a59b45e17b6d7c5efe35ef9ba9f7eb6e73e8df297deec7ffde5c79f7efcd3b73fc9adfdaf3fffe1db9fbf03",
      "9302b24db1f73cf41a68ca1156efc349b0fd773fea19f48702b4387c07128c3fc404fa5fe9d0230c09c7d7e06a79f8a9af3905bfffddfff9f73f7dfb",
      "e97ba2a3e79092a086c62c5a166ac4cc938f9a26b3755af49477ae512e903a2429f5e64c13154e66c789a5a881ae3775d736ff59a7b4d48be5f120c5",
      "a44422319f284f03e2a501d2937df880cb8b6f23fa0ba76457c1359d1a3778be7996440cf47b51a2f545a9c14d123ec5bdf012322123340c2302050f",
      "55070ca5a7dadeabced8a66aa20a02e0867480555b0915b65e6099f669b9b1244bcb396e648fa24c4d36e2b788e0ad9249ea251e6f76b257cd946a22",
      "eabf16e67d2ffe88669d2c790860a9aa9cb8936b612c519b876cf3ddc026a723ce29a2cf7e838ea24640090086a9b84eebe5de9e8d21eda1b9314981",
      "4d6b86547dc10f0350253bae5455903452d3a6524426d25abeeae977022f0893a94f519d5662f4141206dde724f0adf74d718cb54bea8591382800b4",
      "c858fa9da385b4865255cd8d93f763d241f6831340fd8306d31e272d005a62b5def7120583bf0eb240d421bc8dd333fab164f152d3f850f5782d8b7c",
      "a89b7c2cac7c5279792ecdbc546e86dcf2b348bee4719e923c2f29a09724d14b1ae979556029fa7ffedfff0f9f633c6e593c0300"
    ), collapse = ""),
    states = paste0(c(
      "1f8b0800000000000203ad7d6b8f654791ed5fb1fcb969e5fbc137cb0c83076c1078e072111f0abbb04b6e77595ded613c88ff7e63c523f73e27a3ee",
      "f8144896dcbd3af7e3e4ce8c58f1ccbf7ff8fe87efee3ffce9873fbfbf7bfffdbbfb8f1fdfbcb9ffe2fdc3e3db0f5f7df857c19e3efce99ffe7e358e",
      "fef5bb778fdfddbf7bff807ffffb876fefbec53f7ffaf0f6edfdd3e3fb3b0c787c7a7ff706e067f4b73777ef3ffc6969af439d89fe460ff8e94f667e",
      "9d5bf9c7ab0fbfba7ffcf6fefdbb1f70277dd06f1edffcf015bfc6178f8fefbe7c787bf79edfe44f7ffac998af67edfd55e9afd368f3cfaffef49319",
      "5ec71c0743337585323d4ca0a148e1213964bbaad522c88c8a8cdae53e31c955f1756811488c5d81d40a3f3ea498152a6dd09ddaeb99aa22f4433390",
      "d19b0deab13586526f0a8d1281f43e8a22b3e12d094a539e975e87182743f60ae975ea85a1565b34482f0cb32683a63c304c7d87f4ba964cef5e5f4f",
      "faa10af5cac818732a427f1a807aa843a136e8c31154a3be295dc7b357e9c13d2ea8d557a5bc9e73bd569fad03eafa9be915428d40eac8f6c03c2a5f",
      "57e967d98f6e29012a39d9644d1e436bc6becca09b034a749d422d757e5c18d5a034302abf9eb9db97a00b32201ab43e6ac547cdaf5b5cd7a59180d4",
      "10d39a3d06a2ae287a4c888264b9aabece754eb92ac855ed75a9a9eb8d86417916f975d5a0da42c374e69915e9f48f404a0f3668e4c21faba5ba468d",
      "c0506f6b54ad6df227eda91894317bb414d62bd43c7809c539ab422d445ec5b946bbaef19769f885c76b6109613936bbb097ca6b6f96b2a03e79d4d4",
      "5585651f8a6ed1b420d94925b566502909502d6d8d9a7c5d2f436e453b2f561e346b490a9100785506fd4f3f3221655440241986412d2640f4fb9a41",
      "094b61bc6e25c50505beb0d751154ab1f0ed67d019ed581e347f93964098faed631d83474d5d20805ac2a8acdba6bc1eb15741f443177a76e0775863",
      "729afcb8d683ac4f5a9719cb83a0306d2fd31ff971c5644ecb587b988469126e64b975ea817f1e494fbe4d8821ebdf69c1962ba84c016455634c4afc",
      "a1660bd56e133343234ebbac27969d4d3e0a01a308529aecb36bd9fde73fffe31faf7eb476797cfbfeeeeda56ef9fca45b48e854d32d3186d7a1dcae",
      "5c62a42d5246d20fcb3f3fc60cbd55058b8ad1fd695b0ac6bf3f065a9776a52caf184842cca4a3ea81f5ce2b3ad2165a586abcfd676da771bcff4309",
      "0a55ec34deda2443933d95b64a9371750a365f475aa832aecbfb921eeb51a032dbc25a67295b4b3cb051050b22c481c559046bfa58c2688db240ee6b",
      "4e480f329064a302a00dca58cdb27080f52e179639f466a4017213a5d0f4b7026b89c53b498235aecd10f5da6ed8e8d07cd00c211e1f2cf1b816e734",
      "2cea38282ec3721ffc8c3e64cb0233e13c5a5a580d89e73da6920c6b23b04e4c65aee7f6592aebc9560d9a2c3ba04d93fe34fab22ada7bb6572e105d",
      "fc88aa820158ce35ab0a58588da2c249c8a78525c146b29f4653147a3e2b76c6a648f858fa7a2eff5a5249f17895340a4bfd96e3588f1895c9442b71",
      "bdca2899490e4da26298a128ec68daabd00c8d9179d3d7b030d2a2ccab72a9d1b01e84a2d1b079c2f85a6258c330e21443e4bfadaaab6d7b9360f9ec",
      "f1ddfbaf3ff8d9dd3757ccf5b39f2de90211d60fe912c02ed3cdd2c5d31aae6e713490a7a71c6de6e8bc4d316eda7357b1be1e76b4f5aed377c5efb1",
      "839d43384cc3e1230e6bb9a636a4104be4dd45b457961b7dafd08413919250619a88d21ba612ec5939ecc8eb4dae933991e63c6984ebcf7dd3a2fcc5",
      "dddfee1e1ececbf1179fd8724c111b79aee5583b66db5f8e9f7effe6fdc3f36b925ebcd28fcdf4a3e3849a93890016460496b34af34a534c94141889",
      "b20383e8216c44958c74ed20aa0b6c7659747c3fecf7142085d7b01619220b625dda4715aca88027ac90a5052c7415d284a554f855e65c4f209321f3",
      "9b5495c7952454a04f47580d763bc272e44b4bcdf9b81d991b84a53cd763697bf0af082357c3da20291369b3c6b65ef962f6fe4ca8cd2a711dbc766f",
      "e598c09901b5b2a0124b015495ac03aba3c995231b445a96a1d1ca8145e28484cdd48661241308232d1e53358ccc1d9ebf799afbf3cbad776e903809",
      "97c77edc927463e65b369b4392daa9f2b830fb308ca43b8f4b5109c6f5fdd673481c25fcc24882ca7e357154dadcc06843a785913dcf58b76f4abc54",
      "9fdde37a9d49fc0550ab6961f467be5d4e4a1f08233ecd6fb8e6bfc3f6e64bd35aab576fb7de9a5616edf65709c669d5b1133273008b59d94c05a526",
      "d64f58ca6d61b961c213c4fd5c58c16227cbf6b8b4b4cc8f1e633de1e2a9b749914fbebcfbfaf12c443e39741a9115a20b07630603c9ff22c6eca968",
      "5f957b2adfa3061e85d8a98647497ceae2511c87097984c923561e01f3889a47e83ce2b713448f48ba7cd3a3a51e7df568ae47873ddabcd36b8f86fb",
      "74dda3f53bfdf7cc04df9cd86c0ecf30f10c18cfd0a155c8fe23dad4c7b88a094d30fe4e041ef78b503f7a2dfc8158a2346e0191fd454488a2b1d499",
      "9afaa7ec3d6858e8e2c68af6be1defcbbf2b99dc630c2e45388decb776d881fc5b9bb295081f91cc256d8069589b99975d1cc7b85a07af8952541846",
      "b8aac4e4acbd2c625de983f3b24b6d8d6bb262a15dd72f8bc2c2d6df43632076d510c0e6607f67cee1f8f53c232487533a7e558e451d29f99fe7f77f",
      "b87bfafae1ed57ef59802d49f8878fceec7e9c7c07a415690a5f2009bd17f77ea03711cf4cd8c5a4d27dfbe0bd3da789001034f5dfaac600db998925",
      "c5ccb6c5484f0d955a55fc84f899a4ebe4d316f133016b613201ee55970f4847eefc58e2698a61a7347145e65e0d4b21f3ba682619a1b941f25892a9",
      "9881e61e4325992c6ea84ff668d3ca9b756185cd2060c9a056a33ef5b8b4279595265381d5241395c71a477fe6f94dbaa552063f9739376d4158e922",
      "f2172120acf7c6f64a52a940107d4285540224180a533cfceaa3069641b4f8f3eba7268c88145b6e71e8d64b904c8d0d355206dd9e21f69d2a3cfa05",
      "75f4a80a6f1a962584d1d6a44df10ac611d6ac2138b13e6608f0e51d26e8b66e6fda5c1fbd7bf89fc74bafdc47ffd77656a69f4febe4e018243073bd",
      "7d678509e7567a95d9ac5c4eae806046269a56e212fe648e0cc672e94babd132ca2497f352d764ac744055fde8c03a313360f4a14e4e0f9a4c60f160",
      "13a50b364f181995af32d4e2f18806e64018ed89be1e017292a1664fcc81640f30da77753d8218c32b9abcd0fac112c8c8c8500887172957da33b92c",
      "f72a5f0a5e93f146791196da68d910366cf9e2f552628c6ce4756d2bb469329c86a76b4798c0884d2c6283cd078cecf9e37e95a875c6c61cf384f1fd",
      "462dc72c67da221912ece062f47a19df36aa2343a600f7bbc4124f4be3ef7e90276213c0c27abd8028142eadedc0888a33b6880d7478a8b2a4ca222c",
      "850c9557193a792e476f019b674cbfc6f572bc69c37c7cf7e6e1af8fefde3e5cec998f9736ca7098ccd39e81387f8936f2bfa7f7ddbdf5e1ad236fbd",
      "39cbd25bbdde2af77683b76b9ccde5ed416fafee7b9a575e04064b79512b840d210f725d7a3b9190e5e7aa070c1858197e4532f636602df0eb910e9d",
      "869101c31345927698fe853b25b34fdf10e27103502ea39a70ce309b30c5455d22d0c8441cf85384d10ceba5f036a2357e68d558648d37358ca1314a",
      "af58bbb5aa99032cf272eee3501ab3e25390d2202b76e93c521304913eae06114ba21f3111e79d877e2375432faedb9990d8412c897397b4345e8a50",
      "bee0cbf6fde1ba9902a5f5b3023f71997383b73f5eb69d486aa68f89dfb9c8cdd54abf6d3b3ebe797c77f7e58591fcf1afd7669460dee1f88557e805",
      "4672486a4540fdea520456204811c133a5062caa64d1004f863f512195b781263404c57475d2ab0dcc28637d85a4cad8b5a6a3496b66734683d6310c",
      "bcef54a8aef8165b47670ce1a4a6e3f2e11c855f833195acfcbbc68ed99ce4bc6398a7dbfcf8f7ff75f7e5a507fff7eb43d2768377764955f0f29751",
      "fc4b2bcf33053dd5f28c0a725495a3d13cc5e729484f91ba0ad751ccae027714bd4b083645e3ec526733ef3b7e970a6bca81f5b21bdb37ae91bf7df0",
      "e9fd7f3f7c71b1e13ffbf4c45891b9716c7838b8cb0b19eb153b75f7debe43bd7deced774f2e647652002b66812e2c235f665e62f1b054099b698e4b",
      "0c66db25d274515c617156c6a26dee46f328d7765b6424545251641e4892494ad1a5f5377ddd5fbfbbffead2cafff56f4ffece9149cb9dadfcd15ee4",
      "efbc767538de10cf69e23a571c278ce7acd99d3a8eefc7771179aea4dde5f48c6beada83e579b9b6dde9a8768701c0e6cd08a271ded8220f08bae091",
      "6194453ba64c524db93c6f7ebb56ba67ccef46bfe31bf05d089eabc17349389e8bddc1e1f8417c7789e756d9dd2f9e9bc677e76c6e1fcf3df48c1b69",
      "f336794ea9d316b969fffee7fbbbafcfbbf73f3f3fe96f92cbfdec49a82fcaeff1f4b5e7e20506d73a3c0afd799ee48874cf62742dcb6768c2f9f56e",
      "f373fef0f8edc3dbaf2e9c9c7f5ce22ff11b1e8a0d36c14b18d0337e732f91e83ae1c8c94b72d397bc3427271baa200cc1a99dc415e3f3dcd3e3a81e",
      "97f538effec9bd95f1dc0af28209e7c9bbd1d3f6cddddba7bba70b57db6f4fc4a5d7b6d4db4c90a4e325d9d5ad25661fc45757ca20e2fee090960fd8",
      "7b6562d8d54ae1c4be90190a4277008dc154b16a181479daa50ba394f80d676e171ea49a6822f43f854daad10ea891c482d7a244bbae14f13c90925c",
      "f9ddfc0a6cf55bfe631d857d11a5164bf16e73b26723cf6ef74292367b31dabafd1ce2270955c4351282218060e98776e408677608d0fab2dce2c886",
      "3e3c18cdf2943909149e044b8aa64592c5b5a029d0f4c75964884e0c9e0c360648337e906e9cf8de55891c41298b13a5aa354f5049b0f9312a584269",
      "e9937d19ed1854e119419052534c411ac42a2823afccd434984d56dd74f4a2258a79425bcf125173e849a07ca453c7aa90cd3b09411e555ab40f1690",
      "e70a1b21f795c88fdc1d82525e2b24d5c2464d4a86e428d079899e57ed6d31f2c7bf5d188d9f7c749299a1a47c2e58182f70c57929e04ea2f89e4dbe",
      "259c3b49e948b7c87c23da10ab6461647e5c1a6b09f6ccd72522400b12d199ec46447e5884cd73ed0352bec085d6f6e979f0a8568501604795c034ae",
      "865ad6be836406119ab6f54b2b0c652b7608d08940922d547cebc9b75a1b98f60e5f15d521883131b0841fb9d8a83c2b73d0168e7b3716db247beada",
      "f848a18eb419d6b404444f22deb2ad491028552d5808d09bfcbcd85725c9ac4d94472ac79c23490c1524ab5084566904d4e6fa7e258a7783a48a8995",
      "c229da70fd0d7b09ba3bdf8b968d7d9c5e07434d2d2c946930f10ca069e52c20000d9322c49dab20ab44a123519a21434693cb7a2b96cb4ec3199a25",
      "db8a1db5f28f8eea609dc8921b3281410b29b00145ffd6542c8f2f88a6ad75e50d92a2660d4aacd3ea1122ff665479644bf6cb05bc3fb14c5850e3f5",
      "9152b15134a56c90d0deacabe2a23274a433969a79a9d51156c261acbcb2da8c6da512565edbf8a576e1c8bc07c31a5438171cc1b995ce2871b9cc75",
      "250b6a17fbf44a06dc24a17eb9e9fd5ffeeee4a15c91eb3910509e2f734f5ebb1d9ff14eba5ecc112f31a4c75e8e9b1c6fbc18864f30c5796c69f510",
      "8b70765f42449932bccea359d9cdace2871e6a380142de265cd349f575e5f803a09a97b614dfff84b96d5a6fb0c71d997fc57463434e1241f1b8b0b1",
      "c619c8811b4b37a60bcf0d3db12e27cd2a5549f889275fce04070ef1d26d1380a50b4b2020b7ce2ecd9e73f8b6e28987a7a7c7efdf5d24947e7a7673",
      "d732fa893c9232f95770474757bb0addd1fa1e37b866100ecf70d888c3591c66836fdcebc50775beb1bf12f6f5e2acaa6dedf92bf47a1d3bab7ddf12",
      "eebe2181c3723ec5b5255a90ad9a4b4ea60d48e8b336382ebc56108e16d9558da38f1cade5ea364f03c25a2ef8e43a3350d2443cf09b5398ab1233f1",
      "64196d679bb3c92457e31da380e94e5a88cb9ce865f2f425cde30439407a2641a197b1680df61f12f397a1400aa2021aaddaa8cc76088dca211f9c85",
      "05676f73f19a28d1b43ad6ed69a973ccad58c51f2600347dd07b65231fb400a346e66c54609ede5fcfbe6c2f24b20369c36ab15ae5b85fd7346ad470",
      "d12a6368d96705d9c88494beee54636328a766a332322fe1b818d1ccbf34024351fd1b03422d8b50ab6941b4ca000d8d3d028ae22ca1896976afcec2",
      "a0c55857f5d960884ce4f55a496e55c232412b12c521325affa74d80cfeefff2eeeee99bcbd8d1bf2d33205ec4044882a356e15f1704f482605bb0cc",
      "736c780e9004b5312f316cc9d9cf183185ca19911750ce4c884657c1c90b8a09519fc54a1ac6884cca469e5662d24ac90235ab4d499ce184cca30585",
      "045980f87c5f6520552e6c8bf2f530c494a8c1a81bb2fad9b848fd797ab7734087293a7cd2619d0e37dd09ecc6721d2eec32668f57efec7ba7e88e98",
      "7684b923f237c5e0d0ae9d9c791ccee77aff4cc4f6d7dfbcb9fbfaf1db8b7df7eb5f2e56c2bf676d3b9432bd88dc7a813287f03a846be7650e7df348",
      "9e43059f21130ee5d88889e398f1dc379b8f67c25046dc2f2faf439c92ad32b55a016b268ac36a98158e5af2c91ea4316c15cd2e3ea5d1e7da27d84d",
      "998343b60d736c7cab6e6d09c0d5123bcdc8ee5a05e4684bc0a3d68a9c855d5d338e55dc46df8707f56cc56db9cb7523af51acb4110c3a24166a5d30",
      "c820626a439c87a48faa09bfcab3176df606941427d9c4a6848d3e342c7a7805733588f3d4393dc7c87c857300de44ebf380ec9858c5e7d8d75e0abc",
      "1c906a75707ec1ea72ef1cd8397c1b795f5d62ec41ac7e98f71a3badfc9bf6e6ef1ebf77ab227f77549060cd11a73d5545a617e845bf3bc1dec2c069",
      "74e0b543d89b26ec9d159cf60b4e9386dd8477ec7cc71bb0bb0c1cbf82e37d707c14bba6dbb5a1a332aff5aaa37c5d15ed287247dd6fa460670e3bbd",
      "f05888cf569e09eb38d19fad54d32be9dc4b3f9d0ad1ab3578d346f9d5e3f70f4f0f579d097e75e473d2ef06c73535c611f01738919d80c01e357042",
      "0b5bf801aeca22d993351dfee1cad008ab130d99190cb5900eeb2a71ea64e9c39cb113c41c4677ebebee104344bcb40b023b5561baa1934db5513945",
      "ce0e1975b95e4b14a8ad1e3ac4e401142ddee4fe2e85134d92662460d4908c6fa25d662171b124a058cdcce8e83cc450363b0304ef558681b8ac98c1",
      "3612997e239be9d12762b4f0d8b6b16cabd101154df9049420674900f6ba2cb050785454cdc0ad36e4f6b185724043a003e817404119599aa8d1b677",
      "a23927a08cb18ca8d8794c9ad1022234690550ac73aca8dde40b73e82b28534be3515a1a0a7b19f56380342d95679db874421cd09a0615a8bd04a352",
      "e30948ab68fcc03e94d992122ef20e3d660b1f918aaffc6b525ca3b0670922b952acdb09672433fd6a0b2afc1db2d92bc40d7217c842805933db90f6",
      "50ab4224c619a297b0605b1992e8847a0d1b8524062cacd8ecf664a472a25329cb91c4338825daabbdeae0f8252de434edf6e07c9c0d5552f7827937",
      "8997cfeffffbd2f7fbf9ff59a2257275ffe1fd25d5d85f4090bd5c2d3fa76bcffdda33c4bc3c323fdfecf9bcb4ff8dec78a4c8234fff3f9275903187",
      "b3f9cc6ea77f1b477498a44737374eea30578fde3a1c7863ca3b9d7639b7c3cc1dfebeb17cc716f00c06c7aad86c0fc7427143d57b3cdb897a3bb171",
      "47616e5bd3ddc0ce3677848123323cc1b28b1f474879a26c17788e587484a72762af05312a2d5887f415bc68d8156920c0bb782bf41541399cba6abd",
      "4ac40af35a8fb888feafd5660cd11c2478190e2803487101b80b1c1fd3c22444730118c9445e286ecc1d03a7f9289240fa0126cf0c8091d7989ae47d",
      "ea627f11ae21fa0d31989081ef99d80e616d25d322a4db798a2cdb8f644c03d94eac378bc9225e9e89f9c4925911ab986e379bc9182cf2c92a332d0c",
      "c245d468ed2ba9a987ca1f39eab7e26ca50c91ca4d22560653cb223e832a57086324f8709eecd48ac16b417e6321c4dbb768fef8c5f7ef2f6a213e3f",
      "b9411baa7854dd746405ccdb53f03aed92d9d9a6401e155e1cb7ca49d2a7647df688ae0732280a152308bb120e364d78c5a02a49a53a042d0fc52f67",
      "d3de61944826404ab269e8055ae9ec139b552f0434b240651f253e7b401cf7477690f8d700a143057c7cba47f0fb7a9527aa279020b2e3d8a9983485",
      "f37a166e0be6dd3d3ddd7df1f5f74ff7efdf5f50834fcfa92bb48a56ba5f17c57bfbc772bec3feb19c2f8adc1331124bb66f5ca0e00025a19f74f354",
      "a658a5c390c8cd8212e49a5c47936dd6e510b54bd088956f355483034a02e56cd79559389b6fd86720cd3f2475a4abd82328c21862a7305fd7269e23",
      "99ca49cca306e92a5f548be6e832e9b5876492ba5e9d9b21a1378f3e0fbfa6f3285a4ab6647390154aabb218d425f126a8abfa7ad66fae70f8c5ddb7",
      "df3d7dfd20835640e317a75478925ef9b438c0eb6e5f1cee77d8be96f34df70fefac0e54e3a155227c0f61ac5183ef6e8d080121cc04a78cd6d7030a",
      "b949fe93e8818eae13f079d30fd7dc40c88a16d845812241856297f4aaa11667c7929e9c791f4ab597a0bb4b467d5041048aca79f7391c1f10eb0335",
      "93c9c457d52ea79aa705a4c88d2c628a2fc199e385fd58362a644e260f5a3bc8508d020d9be4a250d2e28aae763e1c5fe1f83a45dc3224436d2fcd11",
      "796672ebb62788ecc9fcc569eb1d099f70c90c5bffa58a2f670a3bb95e0c37add9df7efdf8e5fd079f3cbdb97bfbe579c9fef693b3f241a9c8b1646b",
      "18e945f2ec5265b85a65573dce167536b2b3dd1da1707a839b66e9f7f7efbe7d7c7ba19d7ffff9c9194b36c58576cef17663d0db88bb2ccfdc278bb7",
      "932d1ae44666495f54d20cb500638e9bea6683721e52dfa1a12d8c4a92261dd444c3a8249b2e6974976f2f1d5cf025d6edb9e14ceba68632d2845125",
      "139612cfa577dd3cf9f9cde36c316f23eedb75dfd3dbc6df84832f407631e3082347643982cd117f9e90dc45a92b704fabe1b67cf537777fb9bb8ced",
      "7df4abe5ba6067f8f25d0c34331b2f68063e10c1e5cc7025c1038d54a641e20ce300015be0c4dba64259cbcabb76251ce88e34d973d9abd8e9031734",
      "f65c6ab476702a0cbb45535d770f0d6ebd04725eec56ea3cb1a6493caa3254c58906648622fe4dedbc8beb3a1ba3d6390ca36ae151519701a050d980",
      "9c535f14e50c95dd867388f827a8c52ea352b7512db36531b4a69b9022de4dab2d07523b3b29e18e55a866bc28bc9bea8244d9bb384fb326b4002ae6",
      "30b241a9ab991e1508b3b29b650815bffa7c372daf9fbf797cf77059eefbf3b5bc60efc1d8b2e51539f87afbf2727eba3b41fb2cee33ed7c8eeda339",
      "9fb6405ff0bd6cd30fd8f3f20acddccdc87368eca64eeaf226f9c18e28fa441a181a505f69a8e7ba2c489c119a144d5335b840028c5a3385221cea7c",
      "f75e34c1080a4d7e73d7fefa03edf2c4754dc2547e34fd8c2e667054c23d3887aac15a2e5dbc320c0140bea902a1b6ca96fd3060b04f216947fac11d",
      "ff1bfb0c623308aa8590647943981f3cbec21f633f25769a151ea5ae72f48a097c61cbeba7f42890b92830c159fc1aeb072758efec6d2892753258a8",
      "6576900c7d8784343276a2d41a0c12377cd48fd951a207174e166d308a662284c5ae2019780a26fc6226191202188046b351bd67f11055a9c51bc899",
      "0a0c152d34b85ed437edba7fbf7f7cf7d565eb927fffe842a89795463aa0a1fb0b769db30b5c31e809cb5da4ee72d715ce5722dc11f4bb36705486ab",
      "5836f5834d5d04eacda0ac9d59f4f711b746ce0f8f998ac451d9b73cc4bd38c06fa474a907b1ad4e97d5124c38cc90a4ad485f2b15ad77c54dad3534",
      "a4fea7385ae9356d8527ce94cc471e614263a92a4d50840a0f7853ad0c2a1a34c4b99c9206d2c0b1ba042b9386b690f82993374da86454f022e63824",
      "e638382a28ed5ac2b02d9c92c42f6963d9e322824f8843c6b4101e1353b43b31ebe26f5c4c3ca1ed1b77865a6f49d6abc41cb542cd97748e3cdca4e6",
      "2e5a5df9eb48e96b51eec87b472b5ced9a9bd3c3f1df77df5d6688ffeebcb5118db0ad0d57f2b83d19cb8b463b05734e59dd567ce794e8b9857c4eb9",
      "9f5314e8940e3a05864e19e275b12282b12c4b00e997d949aac78436b2b4f3a98d744ddd8d48b37f3680ec8598afc3d04eb0da0b69ef816f273cee05",
      "d1f750fb7540de09db7bc17d3f05604f14d892099c8403272d614b5e70331caed320dc64092fa5c249bc38b6c30bd2b33ebe7bf7f8e6e132f5e4771f",
      "af4d9be53c23dbb4f0daf4dbfd7cae9cdc84a9277177a9ec886e47c07b6ae05a5538ea64d7399b6272d497abe41c55e828cc4dad3acad755d1bb1ebf",
      "d6f509cca0b01809ebad0a3a8da3567a04e3ae35162ebb8eda0d7b706b08464c7de0630dbe53ee468c27c71bc12a95ef404b057ea77942b2c8a8a81e",
      "166480c82f1e41210e6e4b0c57e798468526c29ae6752a94e0932028aadf07749d539590cb368ba38f6f2be27d43dbe0f1e1221af2c9af4edd23e6d1",
      "861252f34529586ea58a53d1e9d47d6ed5a14e0da95369ead4a3ee45ab4e65ab57feead4c83a95b45ebded754dee5eb6ebd4f63a15c0d775c2ac6442",
      "bc80683d757903d2c866af9092947b69c6e105d4ccdf51f8a02f40d546853c257f528f8503a4f19e1ca6d9ea314bf5405c83587617961f76f704fdc1",
      "ee5e55f48874c77605d5241d242ea0c8dfcb3ada30245fe702eaba8c746b031a5ca59467b2dfd30a9740c5b85eab72b953303e8946ea3d4a89d234cf",
      "03d7c821b05fd7add13d27e380a9e5d56848e54585d2e1c318e8068e72a4b9484a88836ba248ea2e265371211f04b54621cac0d547d50810a2eb68fd",
      "6506fc40a956e7d684312e0a546454d11404509e2e554be65240bd54e59bc7dcd692e1366a7db5882468206114a34ab37cbf2037272669482cf20a51",
      "23997e1d9353ede4d4443995535e7dd55684b5176a5d5773ed255f7e5d98533de6d4983995684ebd9a53d5e6d4be3915724e1d9d536db7d5e439957b",
      "5e799f5303e8550a9ea5f46d6ae4ed97d799bc9f7c76d222631c99bca3a1e0b3bdc467ed6d99eb8de56c3f77933a5bd9d9ef8e50b8161c8e6cd90590",
      "2ba61c61e688bc4d303ae27313b28d93d8aea091b442cb81eae10ad31aba5ee3b299f9e486d328b6bf58d2b7de97bf952338b46eb58dd6600634b92a",
      "d70c3834cbcdbcd0a3b1ba82ec375ed5768c0e417dca9e1c5556f5d0132af079465c0184a85fcc786485c5c81fbfcd358a3e55e2511a9803942adfbe",
      "66bb8e44215f57964bac71392596910203b13b0075d8f449a61c9fef3715a2afcaef446cd23c89398904882a2706a75a36b970da852d8aac22e6dcd7",
      "d711214794d6d66840695cee9631c85a3896735f5a6874d469e0ba8534d448e0b2d46c21d304c9a0b036c04c7a2b4d9cbcde72b7755fb87ffbfefb2f",
      "bef9e1a2ffc2eaab05898f03144d28a0bd5f9a2fa99d77ca5ef7da58a780d62bb3758a71775577ad0f1da5e9aa564701ef5ada51e58ec2f768c1c61d",
      "3669e97c607719388b655f51ceb2bb5e9bce027697b9b3199c2de36cac6dfb399bd4d9c8fb5edfe58123341cd1b20920474cf9c26c17798e602449cc",
      "8a24984748cf05007518cbff1d70b61bd3108b8c0cb4a4e4370d06356d4bddf4cc660419b2dcaa99798b022db9553d19e6c5469553b086e7a6d5e51d",
      "e85d08132d390bd6e8ad4aac36a8f21be4148e21a19dc9184c754e28e7ac6e1b5562663a6809d7f03c700308e8f5badceb0ad5ba1c2bc4faa48977b6",
      "d01ef12dd95d36ed30bd122fbe94cd61d1d8db84e6e631de04f53079c3ad801c52b7aa14aa66f53cc26391593a8d9496136348cf0e1abd46b1fbb521",
      "25c38297b10f86ba9e250d17094eb78558ebc9821f7914117e75794d1a7b4891963c57cca44a316e6d167e493932d234df7770eebf5cd66da786290d",
      "41682b9a0a2402c0a35a5e222571737b40c688d871cc0279713411beeb2d21e6a4ccb7d562e231b4a9450dcb10b890f62f38aed2f3fc7df6f1a974ba",
      "c5a376ba83f9bd80a9ee6ea92d4a5570be3a7b97e21c2b7c39a4cb5f6a7141ecba4674cc6259a3b39baa689af8fabef85eadafbd3d046a6b61244eef",
      "a8afb55a7e70fe35dfbbeb1997709e25e94648c46d055d53e5078e316daf81334bc788bab6f2cceccf9fc7bd52ca53fad5ac4d1ab8fc1abd025a5cf1",
      "1f51d6a967736336f6c4a3aca2e715389276cfa505f351b63e65258c15b90f4d1a40d7be8250d7a33a9aec2681b4cf7ee76320a51d4e9672ae8e5a10",
      "a973afda5eb3a351ab140513896f0a8d36641f69c34d1fea3884145f430f75a3dbe7968b3804e5eb77f0a1dea4363da98b9093dc511724f4c0771a5e",
      "bb161d07a4eba5dc5c99bbbfd3f3896e8e53c7bdea38612fb7c46ddd0ebe7eb8e83afeeb23211791fe749048acccf80212e95b258eede25838bb1d64",
      "b9f46dd4b4e2905700171101d2fe08d88ddabfd6541a6aa3a3b4832e734d224ac4d8f3984c5d918452e32c990b5c8e8fe33a9cb8be91a46b9ea004bd",
      "ca49dd2b8585f6a23836c78a2990b1282fa13d932c8470f2750e3e05527c9bc1ae0bd2a23ae9018008a28cc49725ed648d04164e1624256c418520fe",
      "2081e60e1da3067f8956e2ba7d1587469dc1125b5896a36752cb0b6a917b78a49cd7d60053c7a1ab4a42e04e8fdc858a3ec5f2b0f7c150d3c37f11a5",
      "a981a19296304b388c0a967fef2b4c1304292b705371662a0f5a11985ef815e29c0b19956f1ee2b4e78d28078b68c90b3e5710da3a46b23720b35ff8",
      "e15c7c2654a19a7198b08e5ce806bb65b1b51443bf607e38212bf0a852c27896d339c46fa3870e8974a9a643481ddaba91db9d01bb34d921d317bbfe",
      "c6ead2b76fef9f9eee2f2a053effecc4235031787279659ab3db23276e53dfbdf3ef557760a785b0df68d86947bcb52c76da1aef2da49c3e539b59ee",
      "d03987f439d4d0a18fd70cd3a1a10e59dd19ad437b776eec10689766bf908c276ea82f0fbc81c83c438a76eab4112c8786b964cda1740ef173e8e146",
      "221daae9f15187b43ad4f69affee1cd925d20eddde28f9965ce6a4a03979d27b9eca9ecc72bd956fab607878f7d5c3d5c967bfffe8e44aeb25e493d1",
      "42bb31bef058737f1d3ba6a763a07a66ecb5b1eb98c4aee1fc62237cb7e737a37ff70ca0686df02b10955d96491eecce4a7325ec8d28b74a7d6d4262",
      "5f7caba86dad38e722cba871a4c715bdbd4dcd8f8538f18d6e558ed4e058b30e5ae9787cc05ce7f20f83d46f99c7a23e7c8c0be6210ebbb0271d63ec",
      "08e71df3541d39c67a264331dfff8f858a38038b368d2384b6601247a631bbc42c871033269036d5e4b0a318daba8edf9b5445bc0189558e709b9ab0",
      "0ccbaf15e941aa65f380507c051794d62a1334a2b08964153338ac3c8bef4abd5230bb38b510fc450b947e2c54f58945cf7de938e3bcf17b9564366a",
      "4a127fd45e087dca71d2b893f64242ba4897cb9afd98c0ae69f0202d0f446db5f0332b6affd1d09869881f5381ae1c8ec49fd9a72d89639384633703",
      "b5084f0e9a130e286666b7e178ab9294f06a7c95a08ca27121ca767b228b12da0ccbfefd31103ba099744735dd5178afcf2b666d77de59137dfe8a42",
      "6d0ae7ceb19963a0a2081a31d8b1ccf4c2b951e81d93cb4d5045791c5e224d7b51747e15fe693e80a0b4bfe88c76046835a45c9e45426e1263d6bc5a",
      "40ecbfc0f75b4ff320ce261c8838da2b842911f236b48e0de14579608bc59e98f9f7d1ec6b1f2e8c8ae2cd250616d713e55e595b73010ab2b453135b",
      "07ae96a101f8ac65930d394e1212885a27d7603030b5994dc22770c874697f3cb47fe08f76c838d0e6dc715c408ea3687727fd48d7d44efb5ec433ff",
      "4ce09ff807646e2f8288785df3236a28369bc390252494f499ec9b920c09add03e035111e2c61263d7bf3b4fba806e3c7ee7e1e98bc7b74f0f97a78c",
      "7f722ad0ccc40f4e496875e41725a1796740a0ec820b2dd7c924b573c11f0a320ce2fca8bc7830e7d2763df3abda8583cfb02eb88141934b1cb31e6a",
      "8d26632d25e99f27821e5db7b51092ae9a767805d759e2c4c06c10dd9b9fd7bb5dd7115f2ca09c2d2ea8f21941649619d4c6e423a6aaf236749be626",
      "a63003abbd56af7a36d55caf50e4b82a13ec135e9bc6677487b9ce2f49ad0c8134f504b96d2850e7e5ba0e6d8991a11ed7010f532fec691a340a1fec",
      "defb4afaed511e38d2cae0ab7c4c5be37475cb5969830f2a9b2b7f38f1d91b68dd17f33a8e838fff8e2b8790d37b0bb6465a3939550e0457470392f3",
      "8600791e594193af9a87291ee53e96270733bbf3e3937a1e908514e556e6de01c4673ce195d6a12e110c94df7b3d30f2a147a0ca6df580e70e7c0d92",
      "7ea5344d995f0d82e2456be5a92379122fb218d986ee2b1d31f028a4f5af518def556d8203fa423032d7cf211a2983e64a758c7a825c3dde01eb0490",
      "f685814994227fe31492d94dd23811f664328f00fdd42ed045ad2390d02d941e794d55045b0fb781acfeae39a003a1693dea6dc5bfb91e04fb410f98",
      "033473912e976d557bf2a12f048d557148f48d8f794ba91d65a27c5db4f46b2089ef1eb5ad20d73c4e29906e2b5fa344d98071e5f824302bdc7cac57",
      "8f30920a37f2cd667696740d7598e3975041223ededd92b9e1a6d1dbdbbbc34f2367efc5bc46553e29119d22a33950688d7249f628adac80be1c9458",
      "f554544934ea22e7560656c889ebb4c35c959ef0bab15c9dab9a2e4a4b8556fb5c48d552f8e5e9897026a2cb426a2bf3400fee099a6387bbf7c97ee4",
      "6e8e503785d5497475d261b7a45927b5764fc0f58ff3d90efd718e06720f10da0f19da0f22da4f2b728e34ba507ab7e9e5fba7f71f789e873ffcfe74",
      "fe051f617b144af4f89272612f54bfc5f3dda8bf931be064103879064e2ec2b56bdbf17f3b5e72c7977eed7077bcf2aeef7e77f03b51803d54e0c413",
      "f6a0831399d8e2174e94c38d85ec011327aae2c55ef6088d13c761884ff208560451985d5e42a35d985768329878500ec9ece1d8f4e810fd126ce966",
      "9eaca25207966e923cdbdad2aacd10e3ad9a6d088b956f4e6ad0ecce04219739febac2c849af9bfa0a589862e29532e7f396e76e9e5e9bb02fb6865d",
      "cbdab1bf1d2bddb1e51d8bdff10b6cbe83177a211c8f86e3f5703c23bbfb6477b1387e9817fb741cff90e345727c4dd71e29c76ff5620798e74973fc",
      "6d9b536e77dcbdd807e8f813af5d8e8e57d2f35dee1ececd0fba3b4bff09bfabe3c3753cbd8e3fd8f11a3bbee5b3bebb492bffecfecddddfee2e3b94",
      "fdecdf4e0af95c6edcd14f34bea0dd1368a0244d9622ec028e81281b9d4499795b3a07632098b5c3114689976b8464ce019cd301a47443483ac975b3",
      "da7515bdd7214ee772c9d4248741150dd9c0cf1055fe681b49bc293b81206cecba8b77bf6d761f9edebf7bf8e2fd078f7ffde0e3c737df7ffb974bea",
      "f3b38fcf333d0fea03c716b18bdb67fa473bea760f9fe706549a315425819d478648b3bb77bfb181e2bb1fae7b8d7dfab3539d473ed579f0936a79c9",
      "e2fbdf1c4cbe5b6a775c5df9b6704e5ae2e593954fc1d3288ad0f2b288f8b024ef56bb895ce4c9da24e831229ec7d2f16bbade4fc747ea7952777fab",
      "e395f57cb7d70bc75f00d7cbe4c739991ddff4e6bf767cdc2ff797efbe77c743efb1a99d736dcccce16f1ec9db89a0c316374ee9304f979fee247663",
      "ba606153f8709ca757578a1cd764a59cf5fc3b5becad49d1d7bad733f27a93ea8eec0f55926e8a1ea2f5a3bcc03776dffc8ffb774ff717e5229ffdc7",
      "39d38fde70c916748eabf925b2c5512a8eea7114540a622e8d619e72b2751999dac7a8a3305d0ca13094c216a2ae990da1a8560fab2c1e54ecba2a6d",
      "2db83e6c4d263740443df072d593eee33b11ef5cd78dc4257ba3cdf556b94aaf5d0d2274e43354399f5d1d5a78d1201563598b4a3afa01f2a05cab75",
      "2b8443826fa5092680aa3cd01af77478bc9bfc9cbeba15121394c30ad5774ba3621fbc957ad07546af1554c21c08ceaa61ebb69da65daa20cb281e2f",
      "b87999fdf1f1dd37178bec8fa7feae115922077b221191fe259d9afd4ec67bbf63a72bb2d73bd9e9b0bcf761ee72605a8ed61736413bf0b16ad3fac9",
      "1249952f9c57df475e87cdfececd6200a4d5b0725b02ce42d95713dc88dad233ad005e405a1077f9ac7ded8d2ac9b2259ef8dc3c572462914bd6ed50",
      "1b059f8a0329dcaf59e51f1965ec8f0b5aa805c119023bf282d96d4dbc78d77fafda3b1e50978383f4b41f410c32193df96024e4e2e6a514a254ff4f",
      "6bf10b2f0ebbebfab25bc5cbb77aca4e2111dc94b22f5b5a1bc656631fa8490a3caae959a71d7e59015460a3d0ad717bcbae2df23a9fcbc42f3053b1",
      "cbaaf40f1d2b2bbd0ba2311e26dcda1ad4342e37f699d2f43394e7ba85ba3d45f7bea37b6f52a783a9d7e7d4e986eaf44cdd3bab3afd57f726ad2fef",
      "30fe9bfbb76f9f7e78f35f77572ed3df7c745264038bec20c94463cb4b7acbb8ae3327517a4ba7e6124340ada5e717b9b3132e36cbf56672f6dbbe29",
      "bd8debec6e4f065c4b0a479eb84267d7739b367474a6af5977fdbb2b6947935fab7b871338ccc1e1171e097198ca6e6f3b46f946051dc2e8d04a877c",
      "3a147527b28e63f77a0ddf688a3ebcbd70827c7a9c3b5b71fa8aedb00636d1f34bbab43bbdbaafdb793b1dbf9dbee05bf770b7c7b89cca66c9bc9de3",
      "cd1c848fc51e27c78d42ba8b57aae1a4ecce61e352c4dbdae089d578b3b6aa68f86808b4a2d058a471430d70e010745031c090c4347b48362a66097c",
      "d2ffa5cf3f12704a17c92ebfb9c17d9924d741f3821bfcc8ac6beccd69fb898a52817d35bb3776ddfbe2eb87afee2e12463efde4c4e2b86bb305a5d0",
      "76aebfa8b2fc3a948f1e6d383b04b3dd8a05fcbb206d9d32d6d0ceb570ddf53ac02c2586ac6734e7c7bf2ab094f24a9847ca3e43ab6f1cb10085e6d1",
      "b70ea1fdb11abb2162dea10c71da4c5a2d7d68c5c9f3a2454d5390f482a13c1a91db24d90babe105e714f1a8d52db571b41450d57e7e8070d62a206b",
      "fc0acf8a3cb1d455c31d70144de1336eac5e7aa03d0a27552c88e8115f982c148770f6e00ba3dd1e392db256eb2a4e0d3d70eac05c0da67be3bd31c3",
      "aa0b19a9696a44b314efce0de95127d06d9424bb2002be12c1094a57500b92d9522c9ba0a0e7255f48bcc65e2296c42f5163b2ecf4982443a558a13c",
      "0eb4aa3c0a894dab73ad64c924e52a8343e59cea11e3aafd14bd81c47d439afc427adc4a998f320fa7cb0a0a2f2049b499fec0411692ab308e514c8d",
      "e0345bf5a17d4649a0482bb53f5a42d1eab746b440f212542b0e39529385db09d28486a40e04046a4765d162cd4e12e2d62c59663c02bc412869d7f3",
      "a31097e5941fa44cad1ab89c241b2a5bc201193b837964b41e2c307f847faed30d11bc6d7ce4e850ff0ea02c127de5202439800b62a5971564d600bd",
      "15582665a9e9f56a723cf918e8c4a710af2a08b109565be52c3a1e501be3d94240a75af0aaa2d0a93a7ca647cb552717a7df8bd315c6ed1db3759871",
      "9b7d392dc1f6be614e7331af05d99ee5e1e4827819234e5ec9967ce264a86c792c4eb68b9b13e364ce38f9354e168e97abb367f438793f5b76909343",
      "e4661aede9484ece9293d9b4a73f7939527b2695936fe5646539a95bd7f95d4e1698972ae6e4936d59677b6eda9ec0e667b939b9704ec2dc9e55e7e6",
      "de6d197a4e1edf9eece764043a79834e76a19383789da9784d826e3dc8e2e99bcb732cd619f50dd55a47562fe94328d817d62c116142e8b7bfaadc43",
      "86cf74aba47290655a39905e0c8b646800abdac3059792199b81458dfec68a53114962541c8b2d1a131899d205d85489040c070b56d041315471bbde",
      "123f96b66d33ecfc769ab14dff4a979165d8b82fb15c4ed3d23a11a9c6f1fc66580569a9482010a441ad0d205deb26812134c8a38afeded6e5206cdc",
      "2de675edc553ed6d70a04fa69dd2e032925c5a6031e26d982ec81bf2f96034870dbe4a19d626379c0554d4d948563e728ff9d2a40da481751c22d0d0",
      "c77534bbb6268c4babe80e8f40eb361ea74174600dcb9530324486bd1efd3a7e157af385757dbd128e71173fcd7e32fb54e931b05be5c7f9c8c46f6b",
      "e2e61040e74f5b6a038397a171ca853c53b26600c96588ef155aef95034732479217515f55586c324528e6249312c38ad6b646d43a42875669072d18",
      "8c9a0a2875bb942800badc55a400e857c8995d21c048e0c8d9ae68068c630c2aa269491f815a52fa5115e50ff676814ff579556180ea4ce2b00546d2",
      "a83688f55a05b348c36e16906652d12362ac0790b84fc0cab4df802302f048919c7859b22712efc4751d1f60ce3f2a6afd71e4061564200053b512b9",
      "ca57c62565de98de01a0da071958f7f80059dfa0f0e977bc677ad2450a9797ee402225fa7de1f31efcf1a6564a47b06954beb6c0070209d675bfe030",
      "68951b6497251447372ce06c500075c2edfa5858ed61c826576150f8b5f8b1f49155e4f061292c0c8806abcce15cf9028c34fe306cd4c99836040354",
      "50b9da908363c3b86325ad53f4fe50798347401b36c42a7406703b845ef95a3d851312965e9f674ac9104bce28b22a575d70107fa4bc65611ed23465",
      "59fa512912bd3ace4badf8fef40c7d3f54020d96ce454f158cbc481b4b673d221477c1213115b2cd7e05faade5cccb570f008f15fdda05ab6a11466c",
      "da811d87426493b8dc491f50b6176e38c74476434826ac33aa7a0a167fd3e267fa0cc8d96cc04ab071d8ac2901ab1a6bc6b522aeb198f594d4863337",
      "f9a9e2b6c17b0cb421236898aca50d473f113f35043d71151fa007c68a9e0b889fdf5039549197a23f1553d244fd594950ac1a40874aec8654740aa8",
      "c88ed5d55939b3a0b3a48aba3531e9644e026b4b71228db0f3d65ca2aaa2c91909d5ca5987eb43a01a84ef974d3943b6447e2ecdd2fab5b9ca060ddd",
      "340bbb77f9b13dac799a0dfb13585c1f0cdd7320219bfe564c27883396812df606cb54f7e2fa86381d58542971cd68dfb0a17807cba5993ec3976abc",
      "dd6b386148abc586d27a7b56d7e8bb47d81ad411676b5c9bbcd608360e41f46eeb99b5895e5d9b043a0c075834763996631dd2ace3ffd5d6099f79c6",
      "3ab496b87e7fcda2ce8bc989c6cd6cf8195db3d380d15ea767200a93d614f33c35aceb74ba16aa106b7eacad93e1c900a6dd7bf12ea5f2ed6a3d967a",
      "af7269ce76699683e264dcb1755ae071a5b4359d25c9fddab05d0cf643bf02ce045d9f0d420133004a65833afafa42970f5dd9b831d92dc0889dafef",
      "dfaae8fc3386dfc8d7aedfc596f50016a3cd1ddc1a93efd7f56873be3602837ab0ef8822b008a8c596d7d299b4f1e04e2d36c560faf864482d1debf7",
      "23600cecfc2b7090027cbf9a0d89a7d21c3099c96befe008a5c01869e8b59f26d970af3a1b3379ad001242c04a51155b99d1246043339544f0d22fe3",
      "63644db4c3a8a177c1498fd564008cc9ce582c7d91ec0cb73898dd3415408a561e601a2645fe7b5e7fe7872d45d50bbf633ef425fe5eb3a963040381",
      "68a8025abbf1138cded1fa4a980e3d949ae3ac932fa1a72d4610788b4c3ea2de218a6c16fdf91fff0fbba13bb217cc0000"
    ), collapse = ""),
    stop("Unknown map asset: ", name))
  starts <- seq.int(1L, nchar(hex), by = 2L)
  rawToChar(memDecompress(as.raw(strtoi(substring(hex, starts, starts + 1L), 16L)), type = "gzip"))
}

# ---- RUN: starts automatically when the whole script is executed ----
analyse_demographics(root = PHMPT_ROOT, adsl_file = PHMPT_ADSL_FILE,
  output_dir = PHMPT_OUTPUT_DIR, excel_fr = PHMPT_EXCEL_FR, strict = PHMPT_STRICT,
  show_diagram = PHMPT_SHOW_DIAGRAM, png_dpi = PHMPT_PNG_DPI,
  sites_file = PHMPT_SITES_FILE, template_file = PHMPT_MAP_TEMPLATE,
  p1_file = PHMPT_P1_FILE, p2_file = PHMPT_P2_FILE,
  audit_asset_dir = PHMPT_AUDIT_ASSET_DIR, audit_dpi = PHMPT_AUDIT_DPI)
