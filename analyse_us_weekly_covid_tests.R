# =============================================================================
# US background COVID-19 testing: weekly reported PCR test results
# Version 1.1.0 | 2026-10-22 | Standalone: select all and run in RStudio.
# Base R only. No PHMPT trial files, site file or other scripts are required.
#
# Source: HHS/CDC j8mb-icvb, also used by analyse_background_covid.R.
# https://healthdata.gov/d/j8mb-icvb
#
# COUNTING RULE: sum new_results_reported for Positive + Negative + Inconclusive.
# Never sum the cumulative total_results_reported column across dates.
# These are reported specimen/test results, NOT distinct people, all test types,
# or an exhaustive census of all tests administered. Repeat tests are retained.
# Serology is excluded where possible in the source; this is a PCR series, not
# a combined PCR + antigen series. Dates are source test/result dates, sometimes
# substituted with collection, receipt or order dates. Reporting batches occur.
#
# The available series starts on 2020-03-01. Tests occurred before that date;
# January-February are NOT represented by zeroes or reconstructed here.
# Monday-Sunday weeks; retain partial boundary weeks. Empty END_DATE discovers
# the last observed source date automatically (currently 2024-06-01).
# Reporting coverage is plotted separately. After 2023-05-11 reporting changed;
# later counts are not comparable as a complete national testing census.
# Default geography = 50 states + DC (including Alaska and Hawaii).
# Missing state/date/outcome cells remain missing; bars sum available records.
# Negative source corrections are retained, not discarded or replaced by zero.
# A complete grid is a check of available rows, NOT proof of complete reporting.
#
# Outputs in analysis/us_weekly_covid_tests/: PDF, 300-dpi PNG, HTML report,
# daily/weekly CSV, state coverage, missing cells, corrections, source manifest.
# Reuse verified cache by default. REFRESH=TRUE downloads a new snapshot.
# To load functions only: options(phmpt.autorun=FALSE) before source().
# =============================================================================

PHMPT_ROOT <- getOption("phmpt.root", getwd())
PHMPT_US_END_DATE <- "" # empty = all available source dates; or use YYYY-MM-DD
PHMPT_US_OUTPUT_DIR <- file.path(PHMPT_ROOT, "analysis", "us_weekly_covid_tests")
PHMPT_US_CACHE_DIR <- file.path(PHMPT_ROOT, "cache", "us_weekly_covid_tests")
PHMPT_US_INCLUDE_TERRITORIES <- FALSE # TRUE adds AS, GU, MP, PR, VI; missing remains missing
PHMPT_US_REFRESH <- FALSE
PHMPT_US_EXCEL_FR <- FALSE            # also write semicolon/decimal-comma CSV copies
PHMPT_US_SHOW_REPORT <- interactive()
PHMPT_US_RAW_FILE <- ""               # optional previously exported source_rows.csv

ust_date <- function(x) as.Date(substr(as.character(x), 1L, 10L))
ust_monday <- function(x) { x <- as.Date(x); x - ((as.integer(x) + 3L) %% 7L) }
ust_fmt <- function(x) format(x, big.mark=",", scientific=FALSE, trim=TRUE)
ust_html <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed=TRUE)
  x <- gsub("<", "&lt;", x, fixed=TRUE); gsub(">", "&gt;", x, fixed=TRUE)
}
ust_query <- function(base, args) paste0(base, "?", paste(vapply(names(args),
  function(k) paste0(utils::URLencode(k, reserved=TRUE), "=",
    utils::URLencode(as.character(args[[k]]), reserved=TRUE)), ""), collapse="&"))
ust_read <- function(path) utils::read.csv(path, colClasses="character",
  check.names=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8-BOM", na.strings="")
ust_csv <- function(d, filename, ctx) {
  p <- file.path(ctx$out, filename)
  utils::write.csv(d, p, row.names=FALSE, na="", fileEncoding="UTF-8")
  if (ctx$excel_fr) utils::write.csv2(d, sub("[.]csv$", "_excel_fr.csv", p),
    row.names=FALSE, na="", fileEncoding="UTF-8")
}

ust_fetch <- function(url, filename, fields, ctx) {
  path <- file.path(ctx$cache, filename); meta <- paste0(path, ".meta.rds")
  valid <- function(p) isTRUE(tryCatch({
    d <- utils::read.csv(p, nrows=1L, check.names=FALSE)
    all(fields %in% names(d))
  }, error=function(e) FALSE))
  m <- if (!ctx$refresh && file.exists(path) && file.exists(meta))
    tryCatch(readRDS(meta), error=function(e) NULL) else NULL
  if (!is.null(m) && (!identical(m$url, url) ||
      !identical(m$md5, unname(tools::md5sum(path))) || !valid(path))) m <- NULL
  cached <- !is.null(m)
  if (!cached) {
    ok <- FALSE; last <- "Download failed"
    for (attempt in 1:3) {
      tmp <- tempfile("hhs_", tmpdir=ctx$cache)
      ok <- tryCatch({
        status <- utils::download.file(url, tmp, method="libcurl", mode="wb", quiet=TRUE,
          headers=c("User-Agent"="Mozilla/5.0 (compatible; US-background-analysis/1.0)",
                    "Accept"="text/csv"))
        if (status != 0L || !valid(tmp)) stop("Unexpected/empty download or missing columns")
        if (!file.copy(tmp, path, overwrite=TRUE)) stop("Cannot save cache file")
        TRUE
      }, error=function(e) { last <<- conditionMessage(e); FALSE })
      unlink(tmp)
      if (ok) break
      if (attempt < 3L) { message("Retry ", attempt+1L, ": ", filename); Sys.sleep(attempt) }
    }
    if (!ok) stop(last, "\nURL: ", url, "\nNo partial analysis was substituted.")
    m <- list(url=url, retrieved_utc=format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC"),
              md5=unname(tools::md5sum(path)))
    saveRDS(m, meta)
  }
  ctx$manifest[[length(ctx$manifest)+1L]] <- data.frame(
    file=normalizePath(path, winslash="/"), url=m$url,
    retrieved_utc=m$retrieved_utc, bytes=file.info(path)$size, md5=m$md5,
    reused_cache=cached, stringsAsFactors=FALSE)
  path
}

ust_download <- function(end, ctx) {
  base <- "https://healthdata.gov/resource/j8mb-icvb.csv"
  where <- if(is.na(end)) "date is not null" else paste0("date <= '", end, "T23:59:59'")
  tag <- if(is.na(end)) "all_available" else as.character(end)
  fields <- c("state", "state_name", "date", "overall_outcome",
              "new_results_reported", "total_results_reported")
  q <- ust_query(base, list("$select"="min(date) as first_date,max(date) as last_date,count(*) as n_rows",
                           "$where"=where))
  summary <- ust_read(ust_fetch(q, paste0("hhs_summary_", tag, ".csv"),
                               c("first_date", "last_date", "n_rows"), ctx))
  n_expected <- suppressWarnings(as.numeric(summary$n_rows))
  if (length(n_expected)!=1L || !is.finite(n_expected) || n_expected<1)
    stop("The source returned no records through the requested end date.")
  message("Reading HHS PCR results ",summary$first_date," to ",summary$last_date,"; ", n_expected, " source rows.")
  parts <- list(); offset <- 0L; limit <- 50000L
  repeat {
    url <- ust_query(base, list("$select"=paste(fields, collapse=","), "$where"=where,
      "$order"="state,date,overall_outcome", "$limit"=limit, "$offset"=offset))
    path <- ust_fetch(url, paste0("hhs_results_", tag, "_", offset, ".csv"), fields, ctx)
    d <- ust_read(path); parts[[length(parts)+1L]] <- d
    if (nrow(d)<limit) break
    offset <- offset+limit
    if (offset>5000000L) stop("Unexpectedly large result.")
  }
  d <- do.call(rbind, parts)
  if (nrow(d)!=n_expected || min(ust_date(d$date))!=ust_date(summary$first_date) ||
      max(ust_date(d$date))!=ust_date(summary$last_date))
    stop("Download does not reconcile with source row count/date range. Rerun with PHMPT_US_REFRESH=TRUE.")
  d
}

# Missing rows are explicit in the expanded grid. reported_tests is the sum of
# available values, not an imputation of zero for absent state/outcome cells.
# complete_grid_tests is available only when every expected cell is present.
ust_aggregate <- function(grid, period) {
  blocks <- split(seq_len(nrow(grid)), as.character(period))
  rows <- lapply(names(blocks), function(key) {
    z <- grid[blocks[[key]], , drop=FALSE]; v <- z$new_results_reported
    seen <- is.finite(v); n_seen <- sum(seen)
    count <- function(outcome) {
      a <- v[z$overall_outcome==outcome]
      if (any(is.finite(a))) sum(a, na.rm=TRUE) else NA_real_
    }
    p <- count("positive"); n <- count("negative"); inc <- count("inconclusive")
    determinate <- z$overall_outcome %in% c("positive", "negative")
    positivity_ok <- all(seen[determinate]) && all(v[determinate]>=0) && p+n>0
    data.frame(period_start=key, covered_start=as.character(min(z$date)),
      covered_end=as.character(max(z$date)), calendar_days=length(unique(z$date)),
      reported_tests=if(n_seen) sum(v[seen]) else NA_real_,
      positive_tests=p, negative_tests=n, inconclusive_tests=inc,
      complete_grid_tests=if(all(seen)) sum(v) else NA_real_,
      positivity_pct=if(isTRUE(positivity_ok)) 100*p/(p+n) else NA_real_,
      expected_cells=nrow(z), reported_cells=n_seen, missing_cells=sum(!seen),
      row_coverage_pct=100*n_seen/nrow(z),
      jurisdictions_with_any_record=length(unique(z$state[seen])),
      jurisdiction_days_with_any_record=nrow(unique(z[seen,c("state","date")])),
      negative_correction_rows=sum(v<0, na.rm=TRUE),
      net_negative_corrections=sum(v[v<0 & seen]), stringsAsFactors=FALSE)
  })
  d <- do.call(rbind, rows)
  for (nm in c("period_start", "covered_start", "covered_end")) d[[nm]] <- as.Date(d[[nm]])
  rownames(d) <- NULL; d
}

ust_plot <- function(w, start, end, scope) {
  old <- graphics::par(no.readonly=TRUE); on.exit(graphics::par(old))
  graphics::par(family="sans", col.axis="#38434A", col.lab="#26343C", las=1,
                xaxs="i", yaxs="i", oma=c(3.1,0,0,0))
  graphics::layout(matrix(c(1,2),2,1),heights=c(4.4,1.2))
  on.exit(graphics::layout(1),add=TRUE)
  x <- as.numeric(w$week_start+3L); y <- w$reported_tests/1e6
  if (!any(is.finite(y))) stop("No finite weekly totals to plot.")
  xr <- c(min(x)-5,max(x)+5); ylim <- range(c(0,y),na.rm=TRUE)
  ylim[2] <- max(ylim[2]*1.19,0.001); if(ylim[1]<0)ylim[1]<-ylim[1]*1.15
  step <- if(as.integer(end-start)>800) "3 months" else "month"
  months <- seq(as.Date(format(start,"%Y-%m-01")),as.Date(format(end,"%Y-%m-01")),by=step)
  labels <- paste(month.abb[as.integer(format(months,"%m"))],format(months,"%Y"),sep="\n")
  change <- as.Date("2023-05-11"); show_change <- start<=change && end>=change
  graphics::par(mar=c(.5,5.1,4.3,1))
  graphics::plot(x,y,type="n",xlim=xr,ylim=ylim,axes=FALSE,xlab="",ylab="Reported PCR test results (millions)")
  if(show_change)graphics::rect(as.numeric(change),ylim[1],xr[2],ylim[2],col="#F8F3EE",border=NA)
  ticks<-pretty(ylim);ticks<-ticks[ticks>=ylim[1]&ticks<=ylim[2]]
  graphics::abline(h=ticks,col="#E4E9ED",lwd=.7)
  graphics::rect(x-2.7,0,x+2.7,y,col=ifelse(w$missing_cells>0,"#BCCAD4","#245E82"),border=NA)
  graphics::abline(h=0,col="#A6B0B8",lwd=.7)
  graphics::axis(2,at=ticks,labels=format(ticks,trim=TRUE),tick=FALSE,line=-.4)
  graphics::title(main="Weekly COVID-19 PCR testing in the United States",adj=0,cex.main=1.3,col.main="#172D3C",line=2.4)
  graphics::mtext(paste0(scope," | ",format(start,"%d %b %Y")," to ",format(end,"%d %b %Y"),
                        " | Monday-Sunday weeks"),side=3,line=.8,adj=0,cex=.88)
  graphics::legend("topleft",bty="n",horiz=FALSE,cex=.74,
    legend=c("All state/date/outcome cells present","Some state/date/outcome cells absent"),
    fill=c("#245E82","#BCCAD4"),border=NA,inset=c(0,.01))
  if(show_change){
    graphics::abline(v=as.numeric(change),col="#B34D42",lty=2,lwd=1)
    graphics::text(as.numeric(change)+8,ylim[2]*.93,"11 May 2023\nReporting requirements changed",
      adj=c(0,1),cex=.72,col="#904335")
  }
  graphics::par(mar=c(3.8,5.1,.8,1))
  graphics::plot(x,w$row_coverage_pct,type="n",xlim=xr,ylim=c(0,105),axes=FALSE,xlab="",ylab="")
  if(show_change)graphics::rect(as.numeric(change),0,xr[2],105,col="#F8F3EE",border=NA)
  graphics::abline(h=c(0,50,100),col="#E4E9ED",lwd=.7)
  graphics::lines(x,w$row_coverage_pct,col="#245E82",lwd=1.5)
  if(show_change)graphics::abline(v=as.numeric(change),col="#B34D42",lty=2,lwd=1)
  graphics::axis(2,at=c(0,50,100),labels=c("0%","50%","100%"),tick=FALSE,line=-.4,cex.axis=.78)
  graphics::mtext("Source rows\npresent",side=2,line=3.4,las=0,cex=.78)
  graphics::axis(1,at=as.numeric(months),labels=labels,tick=FALSE,cex.axis=.76,padj=.3)
  notes<-c(
    "Source: HHS/CDC j8mb-icvb. Positive + negative + inconclusive PCR results; negative source corrections retained.",
    "Reporting coverage is the share of expected state/date/outcome cells present, not the share of all tests captured. Missing cells are not zero.",
    paste0("Boundary weeks: first ",w$calendar_days[1]," day(s); last ",tail(w$calendar_days,1),
      " day(s). Results after May 2023 have reduced reporting coverage; this is not a complete national census of testing."))
  for(i in seq_along(notes))graphics::mtext(notes[i],side=1,outer=TRUE,line=.1+(i-1)*.8,adj=.02,cex=.66,col="#48555E")
}

analyse_us_weekly_covid_tests <- function(root=PHMPT_ROOT, end_date=PHMPT_US_END_DATE,
    output_dir=PHMPT_US_OUTPUT_DIR, cache_dir=PHMPT_US_CACHE_DIR,
    include_territories=PHMPT_US_INCLUDE_TERRITORIES, refresh=PHMPT_US_REFRESH,
    excel_fr=PHMPT_US_EXCEL_FR, show_report=PHMPT_US_SHOW_REPORT, raw_file=PHMPT_US_RAW_FILE) {
  auto_end <- length(end_date)==1L && (is.na(end_date) || !nzchar(trimws(end_date)))
  end <- if(auto_end) as.Date(NA) else as.Date(end_date)
  if (!auto_end && (length(end)!=1L || is.na(end))) stop("Invalid end_date; use YYYY-MM-DD or empty for full scope.")
  states <- c(state.abb,"DC")
  if (include_territories) states <- c(states,"AS","GU","MP","PR","VI")
  old_options <- options(timeout=max(180, getOption("timeout",60)))
  on.exit(options(old_options), add=TRUE)
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  dir.create(cache_dir, recursive=TRUE, showWarnings=FALSE)
  ctx <- new.env(parent=emptyenv()); ctx$out <- output_dir; ctx$cache <- cache_dir
  ctx$refresh <- refresh; ctx$excel_fr <- excel_fr; ctx$manifest <- list()
  if (nzchar(raw_file)) {
    if (!file.exists(raw_file)) stop("Cannot read PHMPT_US_RAW_FILE: ",raw_file)
    raw <- ust_read(raw_file)
    ctx$manifest[[1]] <- data.frame(file=normalizePath(raw_file,winslash="/"),
      url="Local override: see original download manifest for provenance", retrieved_utc=NA_character_,
      bytes=file.info(raw_file)$size, md5=unname(tools::md5sum(raw_file)), reused_cache=TRUE)
  } else raw <- ust_download(end, ctx)
  ust_csv(do.call(rbind,ctx$manifest), "source_manifest.csv", ctx)
  fields <- c("state","state_name","date","overall_outcome","new_results_reported")
  if (!all(fields %in% names(raw))) stop("Missing source columns: ",paste(setdiff(fields,names(raw)),collapse=", "))
  ust_csv(raw, "source_rows.csv", ctx)
  d <- raw; d$date <- ust_date(d$date); d$state <- toupper(trimws(d$state))
  d$overall_outcome <- tolower(trimws(d$overall_outcome))
  d$new_results_reported <- suppressWarnings(as.numeric(d$new_results_reported))
  if (anyNA(d$date) || anyNA(d$state) || any(!nzchar(d$state))) stop("Missing/invalid source date or jurisdiction.")
  if(auto_end){
    if(!any(d$state%in%states))stop("No source rows in selected geography")
    end<-max(d$date[d$state%in%states])
  }
  d <- d[d$date<=end,,drop=FALSE]
  if (!nrow(d)) stop("No source rows through requested end date.")
  if (any(!d$overall_outcome %in% c("positive","negative","inconclusive")))
    stop("Unexpected outcome category; review source_rows.csv before summing.")
  bad <- !is.finite(d$new_results_reported) | d$new_results_reported!=floor(d$new_results_reported)
  if (any(bad)) { ust_csv(d[bad,,drop=FALSE], "invalid_counts.csv", ctx); stop("Invalid test counts; see invalid_counts.csv.") }
  dup <- duplicated(d[c("state","date","overall_outcome")]) | duplicated(d[c("state","date","overall_outcome")], fromLast=TRUE)
  if (any(dup)) { ust_csv(d[dup,,drop=FALSE], "duplicate_keys.csv", ctx); stop("Duplicate source keys; not summed. See duplicate_keys.csv.") }

  scope <- if (include_territories) "50 states + DC + 5 US territories" else "50 states + District of Columbia"
  ust_csv(d[!d$state %in% states,,drop=FALSE], "excluded_jurisdictions.csv", ctx)
  d <- d[d$state %in% states,,drop=FALSE]
  if (!nrow(d)) stop("No records for the requested geographic scope.")
  start <- min(d$date); last <- max(d$date)
  message("Available dates: ",start," to ",last,". Scope: ",scope)
  message("No earlier testing is inferred; this series does not cover the first January-February 2020 tests.")
  grid <- expand.grid(state=sort(states), date=seq(start,end,by="day"),
    overall_outcome=c("positive","negative","inconclusive"), stringsAsFactors=FALSE)
  grid <- merge(grid, d[c("state","date","overall_outcome","new_results_reported")],
                by=c("state","date","overall_outcome"), all.x=TRUE, sort=TRUE)
  ust_csv(grid[is.na(grid$new_results_reported),c("state","date","overall_outcome")], "missing_cells.csv", ctx)
  ust_csv(d[d$new_results_reported<0,,drop=FALSE], "negative_corrections.csv", ctx)
  coverage <- do.call(rbind,lapply(sort(states), function(s) {
    z <- d[d$state==s,,drop=FALSE]
    data.frame(state=s, first_record=if(nrow(z)) as.character(min(z$date)) else NA_character_,
      last_record=if(nrow(z)) as.character(max(z$date)) else NA_character_,
      days_with_any_record=length(unique(z$date)), reported_cells=nrow(z),
      expected_cells=as.integer(end-start+1L)*3,
      reported_tests=if(nrow(z)) sum(z$new_results_reported) else NA_real_)
  }))
  ust_csv(coverage, "jurisdiction_coverage.csv", ctx)
  daily <- ust_aggregate(grid, grid$date); names(daily)[1] <- "date"
  weekly <- ust_aggregate(grid, ust_monday(grid$date)); names(weekly)[1] <- "week_start"
  weekly$week_end <- weekly$week_start+6L
  weekly$partial_calendar_week <- weekly$calendar_days<7L
  weekly <- weekly[c("week_start","week_end",setdiff(names(weekly),c("week_start","week_end")))]
  # Each valid source record must contribute exactly once at both time scales.
  total <- sum(d$new_results_reported)
  stopifnot(isTRUE(all.equal(sum(weekly$reported_tests,na.rm=TRUE),total)),
            isTRUE(all.equal(sum(daily$reported_tests,na.rm=TRUE),total)))
  ust_csv(daily, "us_daily_tests.csv", ctx); ust_csv(weekly, "us_weekly_tests.csv", ctx)
  summary <- data.frame(geographic_scope=scope, first_available_date=start,
    last_available_date=last, requested_end_date=end, automatic_end_date=auto_end, reported_tests=total,
    last_date_jurisdictions=length(unique(d$state[d$date==last])),
    weeks=nrow(weekly), weeks_with_missing_cells=sum(weekly$missing_cells>0),
    negative_correction_rows=sum(d$new_results_reported<0),
    net_negative_corrections=sum(d$new_results_reported[d$new_results_reported<0]),
    stringsAsFactors=FALSE)
  ust_csv(summary,"summary.csv",ctx)
  pdf_path <- file.path(output_dir,"us_weekly_covid_tests.pdf")
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(pdf_path, width=12, height=8, family="sans", fallback_resolution=300)
  } else {
    grDevices::pdf(pdf_path, width=12, height=8, family="Helvetica", useDingbats=FALSE)
  }
  tryCatch(ust_plot(weekly,start,end,scope), finally=grDevices::dev.off())
  png_path <- file.path(output_dir,"us_weekly_covid_tests.png")
  grDevices::png(png_path, width=3600, height=2400, res=300,
                type=if(capabilities("cairo")) "cairo" else getOption("bitmapType"))
  tryCatch(ust_plot(weekly,start,end,scope), finally=grDevices::dev.off())
  rows <- paste(vapply(seq_len(nrow(weekly)),function(i) paste0("<tr><td>",weekly$week_start[i],
    "</td><td>",weekly$week_end[i],"</td><td>",ust_fmt(weekly$reported_tests[i]),
    "</td><td>",sprintf("%.1f%%",weekly$row_coverage_pct[i]),"</td><td>",
    weekly$calendar_days[i],"</td></tr>"),""), collapse="\n")
  html <- paste0('<!doctype html><html lang="en"><meta charset="utf-8"><title>US weekly PCR testing</title>',
    '<style>body{font:16px/1.55 Arial,sans-serif;color:#243645;max-width:1150px;margin:36px auto;padding:0 24px}',
    'img{width:100%}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:7px 12px;border-bottom:1px solid #d8dfe4;text-align:right}',
    'th:first-child,td:first-child{text-align:left}.note{background:#f1f5f7;padding:16px}@media print{body{margin:0}tr{break-inside:avoid}}</style>',
    '<h1>Weekly COVID-19 PCR testing in the United States</h1><p>',scope,'; ',start,' to ',end,
    '. Sum of reported results: <strong>',ust_fmt(total),'</strong>.</p><img src="us_weekly_covid_tests.png" alt="Weekly reported PCR tests">',
    '<p><a href="us_weekly_covid_tests.pdf">Print-ready PDF</a> | <a href="us_weekly_tests.csv">Weekly CSV</a> | ',
    '<a href="us_daily_tests.csv">Daily CSV</a> | <a href="source_manifest.csv">Source manifest</a></p>',
    '<div class="note"><strong>Coverage.</strong> The downloaded series begins ',start,
    '; it does not cover the earlier January-February 2020 tests. The latest available record is ',last,
    '. The first and last calendar weeks retain only days inside the observation window. Pale columns mark missing state/date/outcome cells; ',
    'the sum shown is of available records. Missing cells are not interpreted as zero. ',
    'A fully populated grid does not establish that every laboratory reported every test. The final source date has records from ',summary$last_date_jurisdictions,' jurisdictions.</div>',
    '<p><strong>Reporting change.</strong> On 11 May 2023, the US public health emergency ended and HHS could no longer require negative SARS-CoV-2 laboratory results through CELR. ',
    'The shaded section and dashed marker identify this change. Later archived counts are retained, but a decline may reflect loss of reporting as well as changes in testing. ',
    'No NREVSS sentinel counts are appended to this series. <a href="https://www.cdc.gov/mmwr/volumes/72/wr/mm7219e1.htm">CDC explanation</a>.</p>',
    '<p><strong>Measure.</strong> Positive + negative + inconclusive daily PCR results, including repeat tests. ',
    'These are reported specimens/tests, not unique people, cases, or an exhaustive total of all COVID tests administered. ',
    'The source excludes serology where possible and is not a complete antigen-testing series. ',
    'Its date generally means completion or result reporting; collection, receipt or order dates may be substituted. ',
    'Large spikes can reflect reporting batches. This is the archived data as retrieved, not an unrevised contemporary snapshot.</p>',
    '<p><strong>Aggregation.</strong> The daily new_results_reported field is summed once per state/date/outcome; ',
    'the cumulative field is not summed. National and territory aggregates are not mixed. ',
    'Negative source corrections are retained: ',summary$negative_correction_rows,' rows, net ',
    ust_fmt(summary$net_negative_corrections),' tests. See negative_corrections.csv, missing_cells.csv and jurisdiction_coverage.csv. ',
    'The weekly CSV also supplies a complete-grid total only when every expected cell exists. Positivity is ',
    'positive/(positive + negative), excluding inconclusive; it is left missing if determinate cells are absent or negative.</p>',
    '<p><a href="https://healthdata.gov/d/j8mb-icvb">HHS/CDC source</a> | ',
    '<a href="https://catalog.data.gov/dataset/covid-19-diagnostic-laboratory-testing-pcr-testing-time-series">Official source description</a></p>',
    '<h2>Weekly reported totals</h2><table><thead><tr><th>Week start</th><th>Week end</th><th>Reported tests</th>',
    '<th>Expected cells present</th><th>Calendar days included</th></tr></thead><tbody>',rows,'</tbody></table></html>')
  report <- file.path(output_dir,"us_weekly_covid_tests.html")
  writeLines(html,report,useBytes=TRUE)
  capture.output(sessionInfo(),file=file.path(output_dir,"session_info.txt"))
  message("Reported PCR results: ",ust_fmt(total),"; weeks: ",nrow(weekly),
          "; weeks with missing source cells: ",sum(weekly$missing_cells>0))
  message("Outputs: ",normalizePath(output_dir,winslash="/"))
  if (show_report) utils::browseURL(normalizePath(report,winslash="/"))
  invisible(list(weekly=weekly,daily=daily,summary=summary,report=report))
}

# ---- RUN: executes when the whole script is selected and run in RStudio ----
if (isTRUE(getOption("phmpt.autorun", TRUE))) analyse_us_weekly_covid_tests()
