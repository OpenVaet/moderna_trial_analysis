# Standalone print regeneration from the accompanying, unchanged table CSVs.
# Set this extracted folder as the RStudio working directory and run all code.
# Alternatively: options(phmpt.figure_root = 'D:/Projects/your/extracted/folder')
# The two analyse_*.R scripts perform the complete analyses using the preceding
# audit-follow-up bundle's inputs; this renderer needs only these three CSVs.

render_V47_V49 <- function(root = getOption('phmpt.figure_root', getwd())) {
  root <- normalizePath(root, winslash = '/', mustWork = TRUE)
  helper <- file.path(root, 'audit_followup_helpers.R')
  if (!file.exists(helper)) stop('Set the working directory to this extracted folder.')
  env <- new.env(parent = globalenv())
  sys.source(helper, envir = env)
  figures <- data.frame(
    stem = c('V47_complete_episode_comparison', 'V48_common_support',
             'V49_AR_flag_analysis_trace'),
    title = c(
      'V47 All qualifying episodes and positivity-restricted sensitivity',
      'V48 Checklist coverage inside and outside the matched subset',
      'V49 Reaction records retained in the S3 safety analysis'),
    subtitle = c(
      'Day 3 is a calendar-date proxy; rows 3 and 4 share the same follow-up rule.',
      'Standardization gives both arms the same distribution across 365 shared clinical/site/month strata.',
      'These counts establish S3 treatment; later S27 table treatment remains unresolved.'),
    stringsAsFactors = FALSE
  )
  assets <- file.path(root, 'audit_assets')
  for (i in seq_len(nrow(figures))) {
    csv <- file.path(assets, paste0(figures$stem[i], '.csv'))
    if (!file.exists(csv)) stop('Missing figure data: ', csv)
    values <- read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE,
                       fileEncoding = 'UTF-8')
    png <- file.path(assets, paste0(figures$stem[i], '.png'))
    env$print_table_png(values, png, figures$title[i], figures$subtitle[i])
    message('Created ', png)
  }
  invisible(file.path(assets, paste0(figures$stem, '.png')))
}

if (isTRUE(getOption('phmpt.figure_autorun', TRUE))) render_V47_V49()
