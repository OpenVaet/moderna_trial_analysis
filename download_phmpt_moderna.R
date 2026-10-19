# PHMPT Moderna archive downloader v1.3 -- R 4.3+, Windows/Linux/macOS.
# Stand-alone: open in RStudio, adjust SETTINGS below, select all, and Run.
# RStudio's Source button and Rscript download_phmpt_moderna.R also work.
# First install: install.packages(c("curl", "xml2", "zip"))
# The final line starts the job. No options() or manual function call is needed.
# Earlier phmpt.autorun, phmpt.root, phmpt.dry_run, phmpt.discovery, and
# phmpt.useragent options are not used; SETTINGS below control this script.
# Every run discovers archive URLs from the current PHMPT listing. There is
# no bundled archive list, fixed storage host, or fallback to saved URLs/HTML.
# If discovery fails, the script stops before requesting any archive.
# Only filenames starting md- are selected. EUA takes priority over production.
# Downloads are sequential; interrupted .part files are resumed when the remote
# ETag (or Last-Modified + length) still matches. TLS verification stays enabled.
# A ZIP directory/size check is not a complete CRC check of every member.
# The listing uses the browser headers from the working Perl downloader,
# HTTP/1.1, and automatic gzip/deflate decoding. Archive requests use the same
# user agent but retain identity encoding for byte-exact download/resume checks.

# ---- SETTINGS: edit these values, then run the whole script ----
PHMPT_ROOT <- getwd()  # Or: "D:/Projects/moderna_trial_analysis"
PHMPT_DRY_RUN <- FALSE # TRUE = list archives and write manifest only
PHMPT_RETRIES <- 3L
PHMPT_PROGRESS <- TRUE
PHMPT_PAGE_URL <- "https://phmpt.org/multiple-file-downloads/"
PHMPT_USER_AGENT <- paste0(
  "Mozilla/5.0 (Windows NT 10.0; Win64; x64) ",
  "AppleWebKit/537.36 (KHTML, like Gecko) ",
  "Chrome/125.0.0.0 Safari/537.36")

# ---- FUNCTIONS ----

phmpt_require_download <- function() {
  pkgs <- c("curl", "xml2", "zip")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: install.packages(c(",
                           paste(sprintf('"%s"', missing), collapse = ", "), "))", call. = FALSE)
}

phmpt_discover <- function(page, base_url = "https://phmpt.org/multiple-file-downloads/") {
  doc <- xml2::read_html(page, options = c("RECOVER", "NOERROR", "NOWARNING", "NONET"))
  href <- xml2::xml_attr(xml2::xml_find_all(doc, "//a[@href]"), "href")
  href <- unique(xml2::url_absolute(href[!is.na(href)], base_url))
  href <- href[grepl("^https?://", href, ignore.case = TRUE)]
  name <- vapply(sub("[?#].*$", "", href), function(u) utils::URLdecode(basename(u)), character(1))
  keep <- grepl("^md-.*\\.zip$", name, ignore.case = TRUE)
  href <- href[keep]; name <- name[keep]
  if (!length(name)) stop("No md- ZIP links found. Check the page or its HTML structure.", call. = FALSE)
  # Keep literal filenames from the URL; refuse anything unsafe as a local filename.
  if (any(grepl('[<>:"/\\\\|?*]|[[:cntrl:]]', name))) stop("Unsafe archive filename on page.")
  group <- ifelse(grepl("(^|-)eua-", name, ignore.case = TRUE), "eua",
                  ifelse(grepl("(^|-)production-", name, ignore.case = TRUE), "prod", "unclassified"))
  out <- unique(data.frame(filename = name, url = href, group = group, stringsAsFactors = FALSE))
  if (anyDuplicated(tolower(out$filename))) {
    stop("Different URLs map to the same filename; resolve this before downloading.", call. = FALSE)
  }
  out$directory <- ifelse(group == "eua", "eua_zip", ifelse(group == "prod", "prod_zip", ""))
  out[order(out$group, out$filename), , drop = FALSE]
}

phmpt_user_agent <- function() {
  ua <- PHMPT_USER_AGENT
  if (!is.character(ua) || length(ua) != 1L || is.na(ua) ||
      !nzchar(trimws(ua)) || grepl("[[:cntrl:]]", ua)) {
    stop("PHMPT_USER_AGENT must be a non-empty, single-line string.", call. = FALSE)
  }
  ua
}

phmpt_fetch_listing <- function(page_url) {
  # Match the working Mojo::UserAgent request, including its HTTP/1.1 transport.
  # libcurl's enum value 2 means HTTP/1.1 (it does NOT mean HTTP/2).
  # Use accept_encoding rather than just an Accept-Encoding header so curl
  # also decompresses the response before it is parsed and saved as HTML.
  handle <- curl::new_handle(
      followlocation = TRUE, maxredirs = 10L, failonerror = TRUE,
      connecttimeout = 30, timeout = 120, http_version = 2L,
      accept_encoding = "gzip, deflate", cookiefile = "",
      useragent = phmpt_user_agent())
  curl::handle_setheaders(handle, .list = list(
      Accept = paste0("text/html,application/xhtml+xml,application/xml;q=0.9,",
                      "image/avif,image/webp,image/apng,*/*;q=0.8"),
      "Accept-Language" = "en-US,en;q=0.9,fr;q=0.8",
      Connection = "keep-alive",
      "Upgrade-Insecure-Requests" = "1",
      "Sec-Fetch-Dest" = "document",
      "Sec-Fetch-Mode" = "navigate",
      "Sec-Fetch-Site" = "none",
      "Sec-Fetch-User" = "?1",
      "Cache-Control" = "max-age=0"))
  response <- curl::curl_fetch_memory(page_url, handle)
  if (response$status_code != 200L) stop("Page returned HTTP ", response$status_code)
  response$content
}

phmpt_listing_index <- function(page_url, analysis) {
  message("Reading ", page_url)
  result <- tryCatch({
    page <- phmpt_fetch_listing(page_url)
    # Parse before saving: a 200 response can still contain an error page.
    index <- phmpt_discover(page, base_url = page_url)
    list(index = index, page = page)
  }, error = function(e) e)
  if (inherits(result, "error")) {
    stop("Unable to discover archives from the current PHMPT listing: ",
         conditionMessage(result),
         "\nNo archives were requested. Correct the problem and rerun the script.",
         call. = FALSE)
  }
  index <- result$index
  # Evidence of this discovery only; this file is never used as a fallback.
  writeBin(result$page, file.path(analysis, "phmpt_download_page.html"))
  index$listing_source <- "live_page"
  index$listing_page <- page_url
  index$listing_fetched_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  index
}

phmpt_zip_directory_ok <- function(path, expected_bytes = NA_real_) {
  if (!file.exists(path) || isTRUE(file.info(path)$isdir)) return(FALSE)
  bytes <- file.info(path)$size
  if (!is.finite(bytes) || bytes < 22 || (is.finite(expected_bytes) && bytes != expected_bytes)) return(FALSE)
  tryCatch({
    con <- file(path, "rb"); on.exit(close(con))
    sig <- readBin(con, "raw", n = 4L)
    if (length(sig) != 4L || !identical(sig[1:2], charToRaw("PK"))) return(FALSE)
    z <- zip::zip_list(path)
    nrow(z) > 0L && all(is.finite(z$uncompressed_size))
  }, error = function(e) FALSE)
}

phmpt_preserve_download <- function(path, reason) {
  if (!file.exists(path)) return(invisible(NULL))
  target <- paste0(path, ".", reason, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  while (file.exists(target)) target <- paste0(target, "_")
  if (!file.rename(path, target)) stop("Cannot preserve existing file: ", path)
  message("  Preserved previous bytes: ", basename(target))
  invisible(target)
}

phmpt_remote_info <- function(url) {
  x <- list(url = url, bytes = NA_real_, etag = "", modified = "", head_note = "")
  tryCatch({
    h <- curl::new_handle(nobody = TRUE, followlocation = TRUE, failonerror = TRUE,
                         connecttimeout = 30, timeout = 60, accept_encoding = "identity",
                         useragent = phmpt_user_agent())
    r <- curl::curl_fetch_memory(url, h)
    if (r$status_code < 200L || r$status_code >= 300L) stop("HTTP ", r$status_code)
    headers <- curl::parse_headers(r$headers)
    get_header <- function(key) {
      line <- headers[startsWith(tolower(headers), paste0(key, ":"))]
      if (!length(line)) "" else trimws(sub("^[^:]+:", "", tail(line, 1L)))
    }
    x$bytes <- suppressWarnings(as.numeric(get_header("content-length")))
    if (!length(x$bytes)) x$bytes <- NA_real_
    x$etag <- get_header("etag"); x$modified <- get_header("last-modified")
    x
  }, error = function(e) { x$head_note <- conditionMessage(e); x })
}

phmpt_same_remote <- function(old, new) {
  if (!is.list(old) || !identical(old$url, new$url)) return(FALSE)
  if (nzchar(new$etag)) return(identical(old$etag, new$etag))
  nzchar(new$modified) && identical(old$modified, new$modified) &&
    is.finite(new$bytes) && identical(old$bytes, new$bytes)
}

phmpt_download_one <- function(url, destination, retries = 3L, progress = TRUE) {
  remote <- phmpt_remote_info(url)
  part <- paste0(destination, ".part")
  meta <- paste0(part, ".rds")
  done_meta <- paste0(destination, ".download.rds")
  old_done <- if (file.exists(done_meta)) tryCatch(readRDS(done_meta), error = function(e) NULL) else NULL
  unchanged <- is.null(old_done) || !nzchar(remote$etag) || identical(old_done$etag, remote$etag)
  result <- function(status, note = "") list(status = status, bytes = file.info(destination)$size,
       remote_bytes = remote$bytes, etag = remote$etag, note = paste(c(note, remote$head_note)[nzchar(c(note, remote$head_note))], collapse = "; "))
  if (unchanged && phmpt_zip_directory_ok(destination, remote$bytes)) {
    saveRDS(remote, done_meta)
    return(result("skipped_existing", "ZIP directory and available remote length checked"))
  }
  if (file.exists(destination)) phmpt_preserve_download(destination, "previous")
  if (file.exists(part)) {
    old <- if (file.exists(meta)) tryCatch(readRDS(meta), error = function(e) NULL) else NULL
    if (!phmpt_same_remote(old, remote) ||
        (is.finite(remote$bytes) && file.info(part)$size > remote$bytes)) {
      phmpt_preserve_download(part, "unmatched_partial")
    }
  }
  saveRDS(remote, meta)
  last_error <- "Download did not complete"
  for (attempt in seq_len(retries)) {
    # A previous run can have finished receiving bytes before it renamed the file.
    if (is.finite(remote$bytes) && phmpt_zip_directory_ok(part, remote$bytes)) {
      if (!file.rename(part, destination)) stop("Cannot finalize ", destination)
      saveRDS(remote, done_meta); unlink(meta)
      return(result("completed_from_partial"))
    }
    resume <- file.exists(part) && file.info(part)$size > 0
    message("  Attempt ", attempt, "/", retries, if (resume) " (resuming)" else "")
    r <- tryCatch(curl::multi_download(url, destfiles = part, resume = resume,
          progress = progress, multiplex = FALSE, followlocation = TRUE,
          failonerror = TRUE, connecttimeout = 30, low_speed_limit = 1024,
          low_speed_time = 120, accept_encoding = "identity",
          useragent = phmpt_user_agent()), error = function(e) e)
    if (inherits(r, "error")) {
      last_error <- conditionMessage(r)
    } else {
      if (is.na(r$success[1L])) stop("Download interrupted; rerun to resume: ", basename(destination), call. = FALSE)
      code <- r$status_code[1L]
      if (isTRUE(r$success[1L]) && code %in% c(200L, 206L) &&
          phmpt_zip_directory_ok(part, remote$bytes)) {
        if (!file.rename(part, destination)) stop("Cannot finalize ", destination)
        saveRDS(remote, done_meta); unlink(meta)
        return(result("downloaded"))
      }
      last_error <- paste("HTTP", code, ifelse(is.na(r$error[1L]), "", r$error[1L]),
                          "ZIP incomplete, invalid, or transfer unsuccessful")
      # A complete but invalid body, or unsupported range request, needs a fresh GET.
      if ((isTRUE(r$success[1L]) && code %in% c(200L, 206L)) ||
          (resume && code %in% c(200L, 416L))) {
        phmpt_preserve_download(part, "restart")
      }
    }
    message("  ", last_error)
    if (attempt < retries) Sys.sleep(min(2^attempt, 10))
  }
  result("failed", last_error)
}

download_phmpt_moderna <- function(root = getwd(), dry_run = FALSE, retries = 3L,
    page_url = "https://phmpt.org/multiple-file-downloads/", progress = TRUE,
    discovery = "live") {
  phmpt_require_download()
  # Retain discovery="live" for existing manual calls; only live is supported.
  if (!identical(discovery, "live")) {
    stop("Only discovery = 'live' is supported. Archive URLs must come from the current listing.",
         call. = FALSE)
  }
  stopifnot(length(retries) == 1L, is.finite(retries), retries >= 1L)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  analysis <- file.path(root, "analysis"); dir.create(analysis, showWarnings = FALSE)
  manifest <- phmpt_listing_index(page_url, analysis)
  manifest$destination <- file.path(root, manifest$directory, manifest$filename)
  manifest$status <- "planned"; manifest$bytes <- NA_real_
  manifest$remote_bytes <- NA_real_; manifest$etag <- ""; manifest$note <- ""
  manifest$updated_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  manifest_file <- file.path(analysis, "phmpt_download_manifest.csv")
  flush <- function() utils::write.csv(manifest, manifest_file, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  flush()
  print(table(manifest$group))
  if (any(manifest$group == "unclassified")) {
    stop("Unclassified md- archives; inspect ", manifest_file, ". Nothing downloaded.", call. = FALSE)
  }
  if (isTRUE(dry_run)) {
    message("Discovery only: ", nrow(manifest), " archives. Manifest: ", manifest_file)
    return(invisible(manifest))
  }
  for (folder in unique(manifest$directory)) dir.create(file.path(root, folder), showWarnings = FALSE)
  for (i in seq_len(nrow(manifest))) {
    message("[", i, "/", nrow(manifest), "] ", manifest$filename[i], " -> ", manifest$directory[i])
    manifest$status[i] <- "in_progress"; flush()
    result <- tryCatch(phmpt_download_one(manifest$url[i], manifest$destination[i], retries, progress),
        error = function(e) {
          if (grepl("^Download interrupted", conditionMessage(e))) stop(e)
          list(status = "failed", bytes = NA_real_, remote_bytes = NA_real_, etag = "", note = conditionMessage(e))
        })
    for (key in names(result)) manifest[i, key] <- result[[key]]
    manifest$updated_at[i] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    flush(); message("  ", result$status)
  }
  message("Manifest: ", manifest_file)
  if (any(manifest$status == "failed")) stop(sum(manifest$status == "failed"),
      " download(s) failed. Successful ZIPs are retained; rerun to retry/resume.", call. = FALSE)
  invisible(manifest)
}

# ---- RUN: starts automatically when the whole script is executed ----
download_phmpt_moderna(root = PHMPT_ROOT, dry_run = PHMPT_DRY_RUN,
                      retries = PHMPT_RETRIES, progress = PHMPT_PROGRESS,
                      page_url = PHMPT_PAGE_URL)
