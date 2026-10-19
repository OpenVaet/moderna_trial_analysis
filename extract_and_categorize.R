# PHMPT Moderna extraction and categorization v1.3.0 -- R 4.3+.
# Stand-alone: open extract_and_categorize.R in RStudio, adjust SETTINGS below,
# select all, and Run. The Source button and Rscript extract_and_categorize.R
# also work. The final line starts the job; no options() or manual call is needed.
# Requires package zip. Unsupported ZIP features fall back to an installed 7-Zip.
# First install: install.packages("zip")
# Earlier phmpt.autorun, phmpt.root and phmpt.sevenzip options, and the PHMPT_7ZIP
# environment variable, are not used; SETTINGS below control this script.
# ZIP archives are discovered in the local folders on every run. No archive
# filenames or remote download URLs are predefined; this script works offline.
# eua_zip/*.zip  -> eua_data/drops_data/<archive>/ + <extension>_data/
# prod_zip/*.zip -> prod_data/drops_data/<archive>/ + <extension>_data/
# Existing PDF/XPT paths are preserved. All extensions are now categorized.
# ZIP files inside a release are placed in zip_data; they are not recursively unpacked.
# Repeated runs compare all expected file paths/sizes, not only the first ZIP entry.
# Flat-folder collisions: identical bytes are reused; differing bytes get a stable
# __<MD5> suffix. Both source copies remain in drops_data, with a CSV provenance map.
# One combined report: analysis/file_types.html, with Release = EUA / Prod.
# Earlier separate HTML reports are preserved in the current run's previous_reports/.

# ---- SETTINGS: edit these values, then run the whole script ----
PHMPT_ROOT <- getwd()            # Or: "D:/Projects/moderna_trial_analysis"
PHMPT_GROUPS <- c("eua", "prod") # Or just "eua" or "prod"
PHMPT_BACKEND <- "auto"          # "auto", "zip", or "7zip"
PHMPT_SEVENZIP <- ""             # Empty = detect an installed 7-Zip automatically
# If needed: PHMPT_SEVENZIP <- "C:/Program Files/7-Zip/7z.exe"

# ---- FUNCTIONS ----

phmpt_extract_mkdir <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE))
    stop("Cannot create directory: ", path)
}

phmpt_extract_md5 <- function(path) {
  hash <- unname(tools::md5sum(path))
  if (length(hash) != 1L || is.na(hash)) stop("Cannot hash: ", path)
  hash
}

phmpt_extract_backup <- function(path) {
  target <- paste0(path, ".previous_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  while (file.exists(target) || dir.exists(target)) target <- paste0(target, "_")
  if (!file.rename(path, target)) stop("Cannot preserve prior extraction: ", path)
  message("  Preserved prior extraction: ", target)
  invisible(target)
}

phmpt_extract_members <- function(archive) {
  z <- zip::zip_list(archive)
  if (!nrow(z)) stop("Empty archive: ", archive)
  # Current zip versions expose entry type; older versions still protect path traversal.
  if ("type" %in% names(z) && any(!z$type %in% c("file", "directory")))
    stop("Archive contains non-regular entries: ", archive)
  path <- gsub("\\", "/", z$filename, fixed = TRUE)
  path <- sub("^(\\./)+", "", path, perl = TRUE)
  path <- gsub("/+", "/", path)
  unsafe <- grepl("^/|(^|/)\\.\\.(/|$)|:|[[:cntrl:]]", path, perl = TRUE) | !nzchar(path)
  if (any(unsafe)) stop("Unsafe archive member: ", z$filename[which(unsafe)[1L]])
  is_dir <- grepl("/$", path)
  if ("type" %in% names(z)) is_dir <- is_dir | z$type == "directory"
  keys <- tolower(sub("/$", "", path))
  if (anyDuplicated(keys)) stop("Duplicate/case-colliding paths within archive: ", archive)
  z$relative_path <- path; z$is_dir <- is_dir
  z
}

phmpt_extract_complete <- function(folder, entries) {
  entries <- entries[!entries$is_dir, , drop = FALSE]
  if (!dir.exists(folder)) return(FALSE)
  if (!nrow(entries)) return(TRUE)
  info <- file.info(file.path(folder, entries$relative_path))
  all(!is.na(info$size) & !is.na(info$isdir) & !info$isdir & info$size == entries$uncompressed_size)
}

phmpt_extract_find_7zip <- function(configured = PHMPT_SEVENZIP) {
  if (length(configured) != 1L || is.na(configured) || !is.character(configured))
    stop("PHMPT_SEVENZIP must be one executable path, or an empty string for automatic detection")
  if (nzchar(configured)) {
    found <- if (file.exists(configured) && !dir.exists(configured)) configured else unname(Sys.which(configured))
    if (!nzchar(found)) stop("Configured 7-Zip executable not found: ", configured,
      ". Set PHMPT_SEVENZIP at the top of the script to the installed executable.")
    return(normalizePath(found, winslash = "/", mustWork = TRUE))
  }
  candidates <- unname(Sys.which(c("7z", "7zz", "7za")))
  if (.Platform$OS.type == "windows") {
    bases <- unique(c(Sys.getenv(c("ProgramFiles", "ProgramFiles(x86)", "ProgramW6432")),
                      "C:/Program Files", "C:/Program Files (x86)"))
    bases <- bases[nzchar(bases)]
    candidates <- c(candidates, file.path(bases, "7-Zip", "7z.exe"))
    local <- Sys.getenv("LOCALAPPDATA", "")
    if (nzchar(local)) candidates <- c(candidates, file.path(local, "Programs", "7-Zip", "7z.exe"))
  }
  candidates <- unique(candidates[nzchar(candidates)])
  candidates <- candidates[file.exists(candidates) & !dir.exists(candidates)]
  if (!length(candidates)) return("")
  normalizePath(candidates[1L], winslash = "/", mustWork = TRUE)
}

phmpt_extract_7zip_command <- function(executable, args) {
  # Quote each argument, including Windows paths containing spaces.
  # An empty stdin prevents an encrypted archive from waiting for user input.
  empty_input <- tempfile("sevenzip_stdin_")
  file.create(empty_input)
  on.exit(unlink(empty_input), add = TRUE)
  output <- suppressWarnings(system2(executable, args = vapply(args, shQuote, ""),
    stdout = TRUE, stderr = TRUE, stdin = empty_input))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("7-Zip exited with status ", status, ": ",
    paste(tail(output, 20L), collapse = " | "))
  output
}

phmpt_extract_7zip <- function(archive, stage, executable) {
  # Keep the existing ZIP member/path checks; also reject link/special entries
  # before invoking an external extractor, including with older zip packages.
  listing <- phmpt_extract_7zip_command(executable, c("l", "-slt", "-ba", "--", archive))
  attrs <- sub("^Attributes = *", "", listing[grepl("^Attributes =", listing)])
  if (any(grepl("(^|[[:space:]])[lbpcs][-rwxstST]{9}($|[[:space:]])", attrs)) ||
      any(grepl("^(Symbolic Link|Hard Link) = .+", listing)))
    stop("7-Zip listing contains a link or other non-regular member: ", archive)
  if (any(grepl("^Encrypted = \\+", listing)))
    stop("Password-protected archive cannot be extracted unattended: ", archive)
  methods <- unique(sub("^Method = *", "", listing[grepl("^Method =", listing)]))
  phmpt_extract_7zip_command(executable,
    c("x", "-y", "-bd", "-bb0", paste0("-o", normalizePath(stage, winslash = "/", mustWork = TRUE)), "--", archive))
  # 7-Zip checks entry CRCs while extracting. Any nonzero status was rejected above.
  paste(methods, collapse = ", ")
}

phmpt_extract_validate_stage <- function(stage, entries) {
  if (!phmpt_extract_complete(stage, entries)) stop("Extracted paths/sizes differ from archive")
  paths <- list.files(stage, recursive = TRUE, all.files = TRUE, full.names = TRUE,
                      include.dirs = TRUE, no.. = TRUE)
  if (length(paths)) {
    links <- Sys.readlink(paths)
    if (any(!is.na(links) & nzchar(links))) stop("Extracted staging directory contains a symbolic link")
    regular <- paths[!file.info(paths)$isdir]
    relative <- substring(regular, nchar(stage) + 2L)
    expected <- entries$relative_path[!entries$is_dir]
    if (length(setdiff(relative, expected))) stop("Extracted staging directory contains unexpected files")
  }
  invisible(TRUE)
}

phmpt_extract_unpack <- function(archive, stage, entries, backend = c("auto", "zip", "7zip")) {
  backend <- match.arg(backend)
  attempts <- character()
  reset_stage <- function() {
    # stage is the private temporary directory created for this archive only.
    if (dir.exists(stage) && unlink(stage, recursive = TRUE) != 0L)
      stop("Cannot clear extraction staging directory: ", stage)
    phmpt_extract_mkdir(stage)
  }
  if (backend != "7zip") {
    message("  Trying zip::unzip")
    result <- tryCatch({
      zip::unzip(archive, exdir = stage, overwrite = FALSE)
      phmpt_extract_validate_stage(stage, entries)
      NULL
    }, error = function(e) conditionMessage(e), warning = function(w) conditionMessage(w))
    if (is.null(result)) return(list(backend = "zip", note = "zip::unzip; all member paths/sizes checked"))
    attempts <- paste0("zip::unzip: ", result)
    if (backend == "zip") stop(attempts)
    message("  zip::unzip failed; checking for 7-Zip fallback.")
  }
  executable <- phmpt_extract_find_7zip()
  if (!nzchar(executable)) stop(paste(c(attempts,
    "7-Zip was not found. Install 7-Zip from https://www.7-zip.org/ and rerun this script.",
    "Or set PHMPT_SEVENZIP at the top of the script to an existing 7-Zip executable.",
    "Completed extractions are reusable; do not delete the ZIP files or extracted data."), collapse = "\n"))
  reset_stage()
  message("  Trying 7-Zip: ", executable)
  methods <- tryCatch({
    m <- phmpt_extract_7zip(archive, stage, executable)
    phmpt_extract_validate_stage(stage, entries)
    m
  }, error = function(e) stop(paste(c(attempts, paste0("7-Zip: ", conditionMessage(e))), collapse = "\n"), call. = FALSE))
  list(backend = "7zip", note = paste0("7-Zip (", methods, "); CRC checks and all member paths/sizes checked",
    if (length(attempts)) paste0("; fallback after ", attempts) else ""))
}

phmpt_extract_archive <- function(archive, drops, state_dir, backup_dir, backend = "auto") {
  name <- tools::file_path_sans_ext(basename(archive))
  folder <- file.path(drops, name)
  state_path <- file.path(state_dir, paste0(name, ".rds"))
  entries <- phmpt_extract_members(archive)
  stamp <- list(bytes = unname(file.info(archive)$size), md5 = phmpt_extract_md5(archive))
  state <- if (file.exists(state_path)) tryCatch(readRDS(state_path), error = function(e) NULL) else NULL
  complete <- phmpt_extract_complete(folder, entries)
  if (complete && (is.null(state) || identical(state$archive, stamp))) {
    used <- if (!is.null(state$backend)) state$backend else "previous_extraction"
    saveRDS(list(archive = stamp, checked_at = Sys.time(), backend = used), state_path)
    return(list(status = "skipped_complete", backend = used, entries = entries, folder = folder,
                note = if (is.null(state)) "Adopted existing extraction after checking all paths/sizes" else "All paths/sizes checked; archive MD5 unchanged"))
  }
  # Extract into a sibling staging directory. Commit only when extraction succeeds.
  stage <- tempfile(paste0(name, "_"), tmpdir = file.path(dirname(drops), "extraction_staging"))
  phmpt_extract_mkdir(stage)
  on.exit(if (dir.exists(stage)) unlink(stage, recursive = TRUE), add = TRUE)
  unpacked <- phmpt_extract_unpack(archive, stage, entries, backend)
  if (dir.exists(folder)) {
    phmpt_extract_mkdir(backup_dir)
    backup <- file.path(backup_dir, name)
    while (file.exists(backup) || dir.exists(backup)) backup <- paste0(backup, "_")
    if (!file.rename(folder, backup)) stop("Cannot preserve existing extraction: ", folder)
    if (!file.rename(stage, folder)) {
      file.rename(backup, folder)
      stop("Cannot commit extraction: ", folder)
    }
  } else if (!file.rename(stage, folder)) stop("Cannot commit extraction: ", folder)
  saveRDS(list(archive = stamp, checked_at = Sys.time(), backend = unpacked$backend), state_path)
  list(status = "extracted", backend = unpacked$backend, entries = entries, folder = folder,
       note = paste0("Extracted successfully; ", unpacked$note))
}

phmpt_extract_copy <- function(source, destdir, hash_cache) {
  phmpt_extract_mkdir(destdir)
  target <- file.path(destdir, basename(source))
  hash <- function(path) {
    key <- normalizePath(path, winslash = "/", mustWork = TRUE)
    if (!exists(key, envir = hash_cache, inherits = FALSE))
      assign(key, phmpt_extract_md5(path), envir = hash_cache)
    get(key, envir = hash_cache, inherits = FALSE)
  }
  same <- function(a, b) file.exists(b) && !isTRUE(file.info(b)$isdir) &&
      file.info(a)$size == file.info(b)$size && identical(hash(a), hash(b))
  status <- "copied"
  if (file.exists(target)) {
    if (same(source, target)) return(list(path = target, status = "identical_existing", md5 = hash(source)))
    extension <- tools::file_ext(basename(source))
    stem <- tools::file_path_sans_ext(basename(source))
    target <- file.path(destdir, paste0(stem, "__", hash(source),
                        if (nzchar(extension)) paste0(".", extension) else ""))
    status <- "copied_with_hash_suffix"
    if (file.exists(target)) {
      if (same(source, target)) return(list(path = target, status = "identical_existing_variant", md5 = hash(source)))
      stop("Existing hash-suffixed file differs; inspect manually: ", target)
    }
  }
  temporary <- tempfile("copy_", tmpdir = destdir)
  on.exit(unlink(temporary), add = TRUE)
  if (!file.copy(source, temporary, overwrite = FALSE, copy.date = TRUE)) stop("Copy failed: ", source)
  if (!same(source, temporary)) stop("Copied bytes do not match: ", source)
  if (!file.rename(temporary, target)) stop("Cannot finalize copy: ", target)
  assign(normalizePath(target, winslash = "/", mustWork = TRUE), hash(source), envir = hash_cache)
  list(path = target, status = status, md5 = hash(source))
}

phmpt_extract_html_escape <- function(x) {
  x <- enc2utf8(as.character(x))
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

phmpt_extract_report <- function(records, archives, title, location, output) {
  esc <- phmpt_extract_html_escape
  bytes <- function(x) {
    pow <- pmax(0, pmin(4, floor(log(pmax(x, 1), 1024))))
    paste(format(round(x / 1024^pow, 2), trim = TRUE, big.mark = ","), c("B", "KB", "MB", "GB", "TB")[pow + 1L])
  }
  html_table <- function(x) {
    if (!nrow(x)) return("<p>No records.</p>")
    paste0("<table><thead><tr>", paste0("<th>", esc(names(x)), "</th>", collapse = ""),
      "</tr></thead><tbody>", paste(apply(x, 1L, function(r) paste0("<tr>",
        paste0("<td>", esc(r), "</td>", collapse = ""), "</tr>")), collapse = "\n"), "</tbody></table>")
  }
  summary <- data.frame()
  if (nrow(records)) {
    groups <- split(seq_len(nrow(records)),
      interaction(records$Release, records$extension, drop = TRUE, lex.order = TRUE))
    summary <- do.call(rbind, lapply(groups, function(i) {
      data.frame(Release = records$Release[i[1L]], Extension = records$extension[i[1L]],
        Count = length(i), total_bytes = sum(records$bytes[i]))
    }))
    summary <- summary[order(match(summary$Release, c("EUA", "Prod")), -summary$total_bytes), , drop = FALSE]
    total <- sum(summary$total_bytes)
    summary$Share <- sprintf("%.1f%%", if (total > 0) 100 * summary$total_bytes / total else 0)
    summary$Total <- bytes(summary$total_bytes)
    summary <- summary[c("Release", "Extension", "Count", "Total", "Share")]
  }
  largest <- records[head(order(records$bytes, decreasing = TRUE), 10L), c("Release", "source", "extension", "bytes"), drop = FALSE]
  if (nrow(largest)) largest$bytes <- bytes(largest$bytes)
  failed <- archives$status == "failed"
  text <- paste0("<!doctype html><html lang='en'><head><meta charset='utf-8'>",
     "<meta name='viewport' content='width=device-width,initial-scale=1'><title>", esc(title), "</title>",
     "<style>body{font:16px/1.5 system-ui,sans-serif;max-width:1180px;margin:40px auto;padding:0 20px;color:#223}",
     "table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:8px;border-bottom:1px solid #ddd;text-align:left;overflow-wrap:anywhere}",
     "th{background:#eef3f8}h2{margin-top:32px}.note{color:#596679}</style></head><body><h1>", esc(title), "</h1>",
     "<p>", esc(location), "<br>Generated ", esc(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "</p>",
     "<p><b>", nrow(records), " source files; ", bytes(sum(records$bytes)), "; ",
     sum(failed), " failed archive(s).</b></p>",
     "<p class='note'>Counts describe source files in the current ZIP inventory, before identical-file deduplication in categorized folders. ",
     "Share is the percentage of total source bytes across the releases included in this run. ",
     "Archives that failed extraction are excluded from file-type counts and shown below. Historical categorized files are retained.</p>",
     "<h2>By file type</h2>", html_table(summary), "<h2>Top 10 largest files</h2>", html_table(largest),
     "<h2>Archive status</h2>", html_table(archives), "</body></html>")
  writeLines(enc2utf8(text), output, useBytes = TRUE)
}

extract_phmpt_archives <- function(root = getwd(), groups = c("eua", "prod"),
                                  backend = c("auto", "zip", "7zip")) {
  backend <- match.arg(backend)
  if (!requireNamespace("zip", quietly = TRUE)) stop('Install package zip: install.packages("zip")', call. = FALSE)
  if (!length(groups) || any(!groups %in% c("eua", "prod"))) stop("groups must contain eua and/or prod")
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  analysis <- file.path(root, "analysis"); phmpt_extract_mkdir(analysis)
  # Each run keeps separate provenance logs and updates one combined HTML report.
  run <- tempfile(paste0("phmpt_extract_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = analysis)
  phmpt_extract_mkdir(run)
  failures <- 0L; processed <- 0L
  report_records <- list(); report_archives <- list()
  hash_cache <- new.env(parent = emptyenv())
  for (group in unique(groups)) {
    source_dir <- file.path(root, paste0(group, "_zip"))
    data_dir <- file.path(root, paste0(group, "_data"))
    drops <- file.path(data_dir, "drops_data")
    state <- file.path(data_dir, "extraction_state")
    backup <- file.path(data_dir, "extraction_backups", basename(run))
    for (path in c(drops, state, file.path(data_dir, "extraction_staging"),
                   file.path(data_dir, "pdf_data"), file.path(data_dir, "xpt_data"))) phmpt_extract_mkdir(path)
    archives <- sort(list.files(source_dir, pattern = "\\.zip$", ignore.case = TRUE, full.names = TRUE))
    if (anyDuplicated(tolower(basename(archives)))) stop("Archive names collide on Windows in ", source_dir)
    archive_log <- data.frame(archive = character(), status = character(), backend = character(), note = character())
    records <- data.frame(archive = character(), member = character(), source = character(),
        extension = character(), bytes = numeric(), destination = character(), copy_status = character(), md5 = character())
    manifest_file <- file.path(run, paste0(group, "_files.csv"))
    archive_file <- file.path(run, paste0(group, "_archives.csv"))
    for (i in seq_along(archives)) {
      archive <- normalizePath(archives[i], winslash = "/", mustWork = TRUE)
      message("[", toupper(group), " ", i, "/", length(archives), "] ", basename(archive))
      result <- tryCatch({
        extracted <- phmpt_extract_archive(archive, drops, state, backup, backend)
        entries <- extracted$entries[!extracted$entries$is_dir, , drop = FALSE]
        rows <- vector("list", nrow(entries))
        for (j in seq_len(nrow(entries))) {
          source <- file.path(extracted$folder, entries$relative_path[j])
          extension <- tolower(tools::file_ext(source))
          if (!nzchar(extension)) extension <- "no_extension"
          # Unusual extensions are retained in other_data; the exact extension is logged.
          category <- if (grepl("^[a-z0-9_]{1,24}$", extension)) paste0(extension, "_data") else "other_data"
          dest <- phmpt_extract_copy(source, file.path(data_dir, category), hash_cache)
          rows[[j]] <- data.frame(archive = basename(archive), member = entries$relative_path[j],
            source = source, extension = extension, bytes = entries$uncompressed_size[j],
            destination = dest$path, copy_status = dest$status, md5 = dest$md5)
          if (j %% 250L == 0L) message("  Categorized ", j, "/", nrow(entries), " files")
        }
        list(status = extracted$status, backend = extracted$backend, note = extracted$note,
             rows = if (length(rows)) do.call(rbind, rows) else records[0, ])
      }, error = function(e) list(status = "failed", backend = "", note = conditionMessage(e), rows = records[0, ]))
      archive_log <- rbind(archive_log, data.frame(archive = basename(archive), status = result$status,
                                                backend = result$backend, note = result$note))
      records <- rbind(records, result$rows)
      utils::write.csv(archive_log, archive_file, row.names = FALSE, fileEncoding = "UTF-8")
      utils::write.csv(records, manifest_file, row.names = FALSE, fileEncoding = "UTF-8")
      message("  ", result$status, ": ", result$note)
      failures <- failures + as.integer(result$status == "failed")
      processed <- processed + 1L
    }
    if (!length(archives)) message("No ZIP files in ", source_dir)
    utils::write.csv(archive_log, archive_file, row.names = FALSE, fileEncoding = "UTF-8")
    utils::write.csv(records, manifest_file, row.names = FALSE, fileEncoding = "UTF-8")
    release <- if (group == "eua") "EUA" else "Prod"
    report_records[[group]] <- data.frame(Release = rep(release, nrow(records)), records,
                                          check.names = FALSE)
    report_archives[[group]] <- data.frame(Release = rep(release, nrow(archive_log)), archive_log,
                                           check.names = FALSE)
  }
  report <- file.path(analysis, "file_types.html")
  releases <- ifelse(unique(groups) == "eua", "EUA", "Prod")
  phmpt_extract_report(do.call(rbind, report_records), do.call(rbind, report_archives),
    "Moderna releases - File Type Summary",
    paste0(root, " | Releases included: ", paste(releases, collapse = " / ")), report)
  # Retire the old per-release reports without deleting their contents.
  legacy <- file.path(analysis, c("eua_file_types.html", "prod_file_types.html"))
  legacy <- legacy[file.exists(legacy) & !dir.exists(legacy)]
  if (length(legacy)) {
    previous <- file.path(run, "previous_reports"); phmpt_extract_mkdir(previous)
    for (old in legacy) if (!file.rename(old, file.path(previous, basename(old))))
      warning("Could not move the previous separate report: ", old, call. = FALSE)
  }
  message("Report: ", report)
  message("Provenance logs: ", run)
  if (failures) stop(failures, " archive(s) failed. See the archive logs; rerun after correcting them.", call. = FALSE)
  if (!processed) warning("No archives processed. Check the project root and download folders.", call. = FALSE)
  invisible(run)
}

# ---- RUN: starts automatically when the whole script is executed ----
extract_phmpt_archives(root = PHMPT_ROOT, groups = PHMPT_GROUPS,
                      backend = PHMPT_BACKEND)
