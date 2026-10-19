# PHMPT Moderna XPT-to-CSV conversion v1.0.0 -- R 4.3+.
# Stand-alone: open convert_xpt_to_csv.R in RStudio, adjust SETTINGS below,
# select all, and Run. Source and Rscript also work; the final call starts the job.
# First install: install.packages("haven")
# eua_data/xpt_data/*.xpt  -> eua_data/csv_data/*.csv
# prod_data/xpt_data/*.xpt -> prod_data/csv_data/*.csv
# Uses the categorized folders produced by extract_and_categorize.R.
# Finds both .xpt and .XPT files; preserves each source basename.
# New CSVs are comma-delimited UTF-8, with column names and no row-number column.
# Missing values remain NA, as in the earlier script. Keep the original XPTs:
# CSV does not retain SAS variable labels, formats or tagged-missing metadata.
# Existing CSVs are skipped by default, without checking their contents.
# Set PHMPT_OVERWRITE = TRUE to regenerate them, including older row-number CSVs.

# ---- SETTINGS: edit these values, then run the whole script ----
PHMPT_ROOT <- getwd()             # Or: "D:/Projects/moderna_trial_analysis"
PHMPT_GROUPS <- c("eua", "prod")  # Or just "eua" or "prod"
PHMPT_OVERWRITE <- FALSE          # TRUE = regenerate existing CSV files

# ---- FUNCTIONS ----

phmpt_csv_mkdir <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE))
    stop("Cannot create directory: ", path)
}

phmpt_csv_write <- function(data, output, overwrite = FALSE) {
  # Finish writing a temporary file before replacing any existing CSV.
  temporary <- tempfile("xpt_csv_", tmpdir = dirname(output), fileext = ".tmp")
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::write.csv(data, temporary, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  backup <- ""
  if (file.exists(output)) {
    if (!overwrite) stop("Output appeared during conversion; retained existing file: ", output)
    backup <- tempfile("previous_csv_", tmpdir = dirname(output), fileext = ".bak")
    if (!file.rename(output, backup)) stop("Cannot preserve existing CSV before replacement: ", output)
  }
  if (!file.rename(temporary, output)) {
    if (nzchar(backup) && !file.rename(backup, output))
      stop("Cannot commit CSV or restore its predecessor. Prior CSV preserved at: ", backup)
    stop("Cannot commit converted CSV: ", output)
  }
  if (nzchar(backup) && unlink(backup) != 0L)
    warning("CSV written; previous copy could not be removed: ", backup, call. = FALSE)
  invisible(output)
}

convert_phmpt_xpt_to_csv <- function(root = getwd(), groups = c("eua", "prod"),
                                    overwrite = FALSE) {
  if (!requireNamespace("haven", quietly = TRUE))
    stop('Install package haven first: install.packages("haven")', call. = FALSE)
  if (!length(groups) || any(!groups %in% c("eua", "prod")))
    stop("groups must contain eua and/or prod")
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite))
    stop("overwrite must be TRUE or FALSE")
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (!dir.exists(root)) stop("Project root must be a directory: ", root)
  analysis <- file.path(root, "analysis"); phmpt_csv_mkdir(analysis)
  run <- tempfile(paste0("phmpt_xpt_to_csv_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
                  tmpdir = analysis)
  phmpt_csv_mkdir(run)
  log_file <- file.path(run, "conversion_log.csv")
  log <- data.frame(Release = character(), source = character(), destination = character(),
    status = character(), rows = integer(), columns = integer(), note = character())
  record <- function(release, source, destination, status, rows = NA_integer_,
                     columns = NA_integer_, note = "") {
    log <<- rbind(log, data.frame(Release = release, source = source, destination = destination,
      status = status, rows = rows, columns = columns, note = note))
    utils::write.csv(log, log_file, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  }
  files_found <- 0L
  for (group in unique(groups)) {
    release <- if (group == "eua") "EUA" else "Prod"
    input_dir <- file.path(root, paste0(group, "_data"), "xpt_data")
    output_dir <- file.path(root, paste0(group, "_data"), "csv_data")
    if (!dir.exists(input_dir)) {
      record(release, input_dir, output_dir, "missing_input_directory",
             note = "Run extract_and_categorize.R first, or check PHMPT_ROOT.")
      warning("Missing XPT directory: ", input_dir, call. = FALSE)
      next
    }
    files <- sort(list.files(input_dir, pattern = "\\.xpt$", ignore.case = TRUE,
                             full.names = TRUE))
    files <- files[!dir.exists(files)]
    if (!length(files)) {
      record(release, input_dir, output_dir, "no_xpt_files")
      message("[", release, "] No XPT files in ", input_dir)
      next
    }
    outputs <- file.path(output_dir, sub("\\.xpt$", ".csv", basename(files), ignore.case = TRUE))
    if (anyDuplicated(tolower(outputs)))
      stop("Source filenames would produce colliding CSV names in ", output_dir)
    phmpt_csv_mkdir(output_dir)
    files_found <- files_found + length(files)
    for (i in seq_along(files)) {
      source <- files[i]; output <- outputs[i]
      message("[", release, " ", i, "/", length(files), "] ", basename(source))
      result <- tryCatch({
        if (dir.exists(output)) stop("CSV destination is a directory: ", output)
        existed <- file.exists(output)
        if (existed && !overwrite) {
          list(status = "skipped_existing", rows = NA_integer_, columns = NA_integer_,
               note = "Existing CSV retained without content verification; set PHMPT_OVERWRITE=TRUE to regenerate.")
        } else {
          data <- haven::read_xpt(source, .name_repair = "check_unique")
          phmpt_csv_write(data, output, overwrite = overwrite)
          result <- list(status = if (existed) "replaced" else "converted",
                         rows = nrow(data), columns = ncol(data), note = "")
          rm(data)
          result
        }
      }, error = function(e) list(status = "failed", rows = NA_integer_,
                                  columns = NA_integer_, note = conditionMessage(e)))
      record(release, source, output, result$status, result$rows, result$columns, result$note)
      message("  ", result$status, if (nzchar(result$note)) paste0(": ", result$note) else "")
    }
  }
  failures <- sum(log$status == "failed")
  message("Finished: ", sum(log$status == "converted"), " converted; ",
    sum(log$status == "replaced"), " replaced; ", sum(log$status == "skipped_existing"),
    " skipped; ", failures, " failed.\nLog: ", log_file)
  if (failures) stop(failures, " file(s) failed. Other files were processed; see conversion_log.csv.",
                     call. = FALSE)
  if (!files_found) warning("No XPT files found. Check PHMPT_ROOT and the selected release folders.",
                            call. = FALSE)
  invisible(list(log = log, output_dir = run))
}

# ---- RUN: starts automatically when the whole script is executed ----
convert_phmpt_xpt_to_csv(root = PHMPT_ROOT, groups = PHMPT_GROUPS,
                       overwrite = PHMPT_OVERWRITE)
