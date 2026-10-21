# =============================================================================
# Moderna P301: background weekly COVID-19 cases and test positivity
# Version 1.0.0 | 2026-10-21 | Stand-alone: select all and run in RStudio.
# Place this script and clinical_trial_sites.csv in the project root, or set ROOT.
# Dependencies: jsonlite (installed from CRAN if absent); otherwise base R.
#
# County positivity: CDC nra9-vzzn, historical 7-day NAAT percent positivity.
# Cases: JHU CSSE US cumulative cases, differenced at weekly boundaries.
# State comparison/fallback: HHS j8mb-icvb, daily PCR test-result counts.
# Coordinates -> counties/states: US Census geocoder, 2020 geography vintage.
#
# IMPORTANT: cases are reported cases, not positive-test counts or infections.
# CDC NAAT includes PCR and other molecular amplification methods. State PCR
# positivity below excludes inconclusive results. These definitions differ.
# Public testing includes diagnostic and screening tests; positivity is NOT
# population infection prevalence, nor positivity among symptomatic trial
# participants. Site county is a geographic proxy, not participant residence.
# The script never averages county percentages to create a pooled percentage.
#
# Weeks are Monday-Sunday. Only COMPLETE weeks within START_DATE..END_DATE are
# used. The default includes the complete week containing the 26 March 2021 BLA
# cutoff (i.e. two further days of community data); set END_DATE="2021-03-26"
# for a strict cutoff, whose last complete week ends 21 March.
# To load functions only: options(phmpt.autorun=FALSE) before source().
# =============================================================================

PHMPT_ROOT <- getOption("phmpt.root", getwd())
PHMPT_SITES_FILE <- file.path(PHMPT_ROOT, "clinical_trial_sites.csv")
PHMPT_START_DATE <- "2020-07-27"
PHMPT_END_DATE <- "2021-03-28"
PHMPT_OUTPUT_DIR <- file.path(PHMPT_ROOT, "analysis", "background_covid")
PHMPT_CACHE_DIR <- file.path(PHMPT_ROOT, "cache", "background_covid")
PHMPT_INCLUDE_STATES <- TRUE
PHMPT_REFRESH <- FALSE             # FALSE reuses verified cached downloads
PHMPT_SITE_IDS <- character(0)      # optional: c("US300", "US365")
PHMPT_WRITE_PNG <- TRUE             # 300 dpi individual charts + HTML gallery
PHMPT_EXCEL_FR <- FALSE             # TRUE also writes *_excel_fr.csv (; / ,)
PHMPT_SHOW_REPORT <- interactive()

bg_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
bg_date <- function(x) as.Date(substr(as.character(x), 1L, 10L))
bg_fips <- function(x, width=5L) {
  z <- bg_num(x); ifelse(is.finite(z) & z == floor(z),
                         sprintf(paste0("%0", width, "d"), as.integer(z)), NA_character_)
}
bg_monday <- function(x) {x <- as.Date(x); x - ((as.integer(x) + 3L) %% 7L)}
bg_query <- function(base, args) paste0(base, "?", paste(vapply(names(args),
  function(k) paste0(utils::URLencode(k, reserved=TRUE), "=",
                     utils::URLencode(as.character(args[[k]]), reserved=TRUE)), ""), collapse="&"))
bg_html <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed=TRUE)
  x <- gsub("<", "&lt;", x, fixed=TRUE); gsub(">", "&gt;", x, fixed=TRUE)
}
bg_csv <- function(x, path, fr=FALSE) {
  utils::write.csv(x, path, row.names=FALSE, na="", fileEncoding="UTF-8")
  if (fr) utils::write.csv2(x, sub("[.]csv$", "_excel_fr.csv", path),
                           row.names=FALSE, na="", fileEncoding="UTF-8")
}
bg_read_csv <- function(path) utils::read.csv(path, colClasses="character",
  check.names=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8-BOM", na.strings="")
bg_problem <- function(ctx, stage, item, message) {
  ctx$problems[[length(ctx$problems)+1L]] <- data.frame(stage=stage, item=item,
    message=as.character(message), stringsAsFactors=FALSE)
  warning(stage, " / ", item, ": ", message, call.=FALSE, immediate.=TRUE)
}
bg_require <- function(d, fields, what) {
  missing <- setdiff(fields, names(d))
  if (length(missing)) stop(what, ": missing columns: ", paste(missing, collapse=", "))
}
bg_unique <- function(d, keys, what) {
  d <- unique(d)
  if (nrow(d) && anyDuplicated(d[keys]))
    stop(what, ": conflicting duplicate keys (", paste(keys, collapse=" + "), ").")
  d
}

# Download to a temporary file, validate, then cache with URL/time/MD5 metadata.
# A failed refresh never silently substitutes an old file. No TLS checks disabled.
bg_fetch <- function(url, filename, ctx, check=function(p) TRUE) {
  dest <- file.path(ctx$cache, filename); dir.create(dirname(dest), recursive=TRUE, showWarnings=FALSE)
  meta <- paste0(dest, ".meta.rds"); m <- NULL
  if (!ctx$refresh && file.exists(dest) && file.exists(meta)) {
    m <- tryCatch(readRDS(meta), error=function(e) NULL)
    if (is.null(m) || !identical(m$url, url) ||
        !identical(m$md5, unname(tools::md5sum(dest))) ||
        !isTRUE(tryCatch(check(dest), error=function(e) FALSE))) m <- NULL
  }
  cached <- !is.null(m)
  if (!cached) {
    last <- "Download failed"
    for (attempt in 1:3) {
      tmp <- tempfile("download_", tmpdir=dirname(dest))
      ok <- tryCatch({
        status <- suppressWarnings(utils::download.file(url, tmp, method="libcurl", mode="wb", quiet=TRUE,
          headers=c("User-Agent"="Mozilla/5.0 (compatible; PHMPT-background-analysis/1.0)",
                    "Accept"="text/csv,application/json,*/*")))
        if (status != 0L || !file.exists(tmp) || file.info(tmp)$size == 0) stop("Empty/failed download")
        if (!isTRUE(check(tmp))) stop("Unexpected response format")
        if (!file.copy(tmp, dest, overwrite=TRUE)) stop("Cannot save downloaded file")
        TRUE
      }, error=function(e) {last <<- conditionMessage(e); FALSE})
      unlink(tmp)
      if (ok) break
      if (attempt < 3L) {message("  Retry ", attempt+1L, "/3: ", filename); Sys.sleep(attempt*2)}
    }
    if (!ok) stop(last, " [", url, "]")
    m <- list(url=url, retrieved_utc=format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC"),
              md5=unname(tools::md5sum(dest)))
    saveRDS(m, meta)
  }
  ctx$manifest[[length(ctx$manifest)+1L]] <- data.frame(file=normalizePath(dest, winslash="/"),
    url=m$url, retrieved_utc=m$retrieved_utc, bytes=file.info(dest)$size, md5=m$md5,
    reused_cache=cached, stringsAsFactors=FALSE)
  dest
}
bg_json_check <- function(path) {jsonlite::fromJSON(path, simplifyVector=FALSE); TRUE}
bg_csv_check <- function(fields) function(path) {
  d <- utils::read.csv(path, nrows=1L, check.names=FALSE)
  all(fields %in% names(d))
}

# Pagination prevents the Socrata default 1,000-row limit from truncating data.
bg_socrata <- function(base, where, order, fields, tag, ctx) {
  pages <- list(); offset <- 0L; limit <- 50000L
  repeat {
    u <- bg_query(paste0(base, ".csv"), list("$select"=paste(fields, collapse=","),
      "$where"=where, "$order"=order, "$limit"=limit, "$offset"=offset))
    p <- bg_fetch(u, paste0(tag, "_", offset, ".csv"), ctx, bg_csv_check(fields))
    d <- bg_read_csv(p); bg_require(d, fields, tag)
    pages[[length(pages)+1L]] <- d
    if (nrow(d) < limit) break
    offset <- offset+limit
    if (offset > 5000000L) stop("Unexpectedly large API result: ", tag)
  }
  do.call(rbind, pages)
}

bg_geocode <- function(sites, ctx) {
  extra <- c("county_fips","county_name","state_fips","state_abbr","state_name","geography_status")
  sites[extra] <- NA_character_
  keys <- paste(sprintf("%.6f", sites$latitude), sprintf("%.6f", sites$longitude), sep="_")
  for (i in seq_along(unique(keys))) {
    key <- unique(keys)[i]; ix <- which(keys == key); s <- sites[ix[1L], ]
    message(sprintf("[Geography %d/%d] %s", i, length(unique(keys)), paste(sites$SITEID[ix],collapse=", ")))
    u <- bg_query("https://geocoding.geo.census.gov/geocoder/geographies/coordinates",
      list(x=sprintf("%.6f", s$longitude), y=sprintf("%.6f", s$latitude),
           benchmark="Public_AR_Census2020", vintage="Census2020_Census2020",
           format="json", layers="Counties,States"))
    tryCatch({
      p <- bg_fetch(u, paste0("geography/census_",key,".json"), ctx, bg_json_check)
      g <- jsonlite::fromJSON(p, simplifyVector=FALSE)$result$geographies
      if (length(g$Counties) != 1L || length(g$States) != 1L) stop("No unique county/state match")
      c <- g$Counties[[1L]]; st <- g$States[[1L]]
      if (substr(c$GEOID,1L,2L) != st$GEOID) stop("County/state FIPS mismatch")
      sites[ix,"county_fips"] <- c$GEOID; sites[ix,"county_name"] <- c$NAME
      sites[ix,"state_fips"] <- st$GEOID; sites[ix,"state_abbr"] <- st$STUSAB
      sites[ix,"state_name"] <- st$NAME
      sites[ix,"geography_status"] <- "Census2020 coordinate match"
    }, error=function(e) {
      sites[ix,"geography_status"] <<- paste("FAILED:",conditionMessage(e))
      bg_problem(ctx, "geocoding", paste(sites$SITEID[ix],collapse=";"), conditionMessage(e))
    })
  }
  sites
}

bg_county_data <- function(fips, start, end, ctx) {
  fields <- c("state_name","county_name","fips_code","date",
              "cases_per_100k_7_day_count","percent_test_results_reported")
  out <- list(); groups <- split(sort(fips), ceiling(seq_along(fips)/15))
  for (i in seq_along(groups)) {
    message("[CDC county data ", i, "/", length(groups), "]")
    where <- paste0("fips_code in ('",paste(groups[[i]],collapse="','"),"') and date between '",
                    start,"T00:00:00' and '",end,"T23:59:59'")
    tag <- paste0("cdc_county_",i,"_",start,"_",end)
    d <- tryCatch(bg_socrata("https://data.cdc.gov/resource/nra9-vzzn",where,
      "fips_code,date",fields,tag,ctx), error=function(e) {
        bg_problem(ctx,"CDC download",paste(groups[[i]],collapse=";"),conditionMessage(e)); NULL
      })
    if (!is.null(d)) out[[length(out)+1L]] <- d
  }
  if (!length(out)) return(as.data.frame(setNames(rep(list(character()),length(fields)),fields)))
  d <- do.call(rbind,out); d$fips_code <- bg_fips(d$fips_code); d$date <- bg_date(d$date)
  bg_unique(d,c("fips_code","date"),"CDC county data")
}
bg_state_data <- function(states, start, end, ctx) {
  fields <- c("state","state_name","date","overall_outcome","new_results_reported")
  where <- paste0("state in ('",paste(sort(states),collapse="','"),"') and date between '",
                  start,"T00:00:00' and '",end,"T23:59:59'")
  d <- tryCatch(bg_socrata("https://healthdata.gov/resource/j8mb-icvb",where,
    "state,date,overall_outcome",fields,paste0("hhs_state_",start,"_",end),ctx),
    error=function(e) {bg_problem(ctx,"HHS download","state PCR tests",conditionMessage(e)); NULL})
  if (is.null(d)) return(as.data.frame(setNames(rep(list(character()),length(fields)),fields)))
  d$date <- bg_date(d$date); d$new_results_reported <- bg_num(d$new_results_reported)
  d$overall_outcome <- tolower(trimws(d$overall_outcome))
  bg_unique(d,c("state","date","overall_outcome"),"HHS PCR data")
}

bg_cases <- function(jhu, level, id, weeks) {
  ix <- if (level == "county") which(jhu$FIPS == id) else which(jhu$Province_State == id & jhu$Country_Region == "US")
  if (level == "county" && length(ix)>1L) stop("Multiple JHU rows for county ", id)
  take <- function(dt) {
    col <- paste(as.integer(format(dt,"%m")),as.integer(format(dt,"%d")),format(dt,"%y"),sep="/")
    if (!length(ix) || !col %in% names(jhu)) return(NA_real_)
    z <- bg_num(jhu[ix,col]); if (any(!is.finite(z))) return(NA_real_)
    sum(z)
  }
  vapply(seq_len(nrow(weeks)),function(k) take(weeks$week_end[k])-take(weeks$week_start[k]-1L),0.0)
}

bg_hhs_week <- function(d, state, start, end) {
  z <- d[d$state == state & d$date >= start & d$date <= end,,drop=FALSE]
  one <- function(outcome) {
    v <- z[z$overall_outcome == outcome,,drop=FALSE]
    complete <- nrow(v)==7L && length(unique(v$date))==7L && all(is.finite(v$new_results_reported))
    valid <- complete && all(v$new_results_reported>=0)
    c(n_days=length(unique(v$date)), total=if(valid) sum(v$new_results_reported) else NA_real_,
      corrected=any(v$new_results_reported<0,na.rm=TRUE))
  }
  p <- one("positive"); n <- one("negative"); inc <- one("inconclusive")
  denominator <- p[["total"]]+n[["total"]]
  pct <- if (is.finite(denominator) && denominator>0) 100*p[["total"]]/denominator else NA_real_
  status <- if (is.finite(pct)) "reported" else if(p[["corrected"]] || n[["corrected"]])
    "negative daily test-count correction; not calculated" else if(is.finite(denominator) && denominator==0)
      "zero determinate tests" else "missing/incomplete positive or negative daily counts"
  list(positive_tests=p[["total"]], negative_tests=n[["total"]], inconclusive_tests=inc[["total"]],
       determinate_tests=denominator, positive_days=p[["n_days"]], negative_days=n[["n_days"]],
       positivity_pct=pct, positivity_status=status)
}

bg_weekly <- function(sites, weeks, jhu, county_daily, state_daily, include_states) {
  cc <- unique(sites[!is.na(sites$county_fips),c("county_fips","county_name","state_abbr","state_name")])
  ss <- unique(sites[!is.na(sites$state_abbr),c("state_abbr","state_name")])
  cc <- cc[order(cc$state_name,cc$county_name),]; ss <- ss[order(ss$state_name),]
  make <- function(level) {
    areas <- if(level=="county") cc else ss
    if (!nrow(areas)) return(data.frame())
    out <- vector("list",nrow(areas))
    for (i in seq_len(nrow(areas))) {
      a <- areas[i,]; id <- if(level=="county") a$county_fips else a$state_abbr
      selected <- if(level=="county") sites$county_fips %in% id else sites$state_abbr %in% id
      ids <- paste(sites$SITEID[selected],collapse="; ")
      x <- data.frame(geography_level=level,geography_id=id,
        geography_name=if(level=="county") paste(a$county_name,a$state_name,sep=", ") else a$state_name,
        state_abbr=a$state_abbr,site_ids=ids,n_sites=sum(selected),weeks,
        stringsAsFactors=FALSE)
      x$new_reported_cases <- bg_cases(jhu,level,if(level=="county") id else a$state_name,weeks)
      x$cases_status <- ifelse(is.na(x$new_reported_cases),"missing",
        ifelse(x$new_reported_cases<0,"negative net revision","reported net change"))
      x$case_source <- "JHU CSSE US time series; weekly change in cumulative reported cases"
      x$positivity_pct <- NA_real_; x$positivity_status <- "not reported"
      x$positive_tests <- x$negative_tests <- x$inconclusive_tests <- x$determinate_tests <- NA_real_
      x$positive_days <- x$negative_days <- NA_integer_
      x$cdc_case_rate_per_100k_raw <- NA_character_
      if(level=="county") {
        d <- county_daily[county_daily$fips_code %in% id,,drop=FALSE]
        hit <- match(as.character(x$week_end),as.character(d$date))
        raw <- d$percent_test_results_reported[hit]; val <- bg_num(raw)
        valid <- is.finite(val) & val>=0 & val<=100
        x$positivity_pct[valid] <- val[valid]
        x$positivity_status <- ifelse(valid,"reported",ifelse(is.na(hit),"no county record for week ending",
          ifelse(is.na(raw) | !nzchar(raw),"missing/suppressed","invalid/suppressed source value")))
        x$cdc_case_rate_per_100k_raw <- d$cases_per_100k_7_day_count[hit]
        x$positivity_source <- "CDC nra9-vzzn; 7-day NAAT positivity at Sunday endpoint"
        x$positivity_denominator <- "CDC published percentage; constituent test counts not supplied in this dataset"
      } else {
        for (k in seq_len(nrow(x))) {
          vals <- bg_hhs_week(state_daily,id,x$week_start[k],x$week_end[k])
          for (n in names(vals)) x[k,n] <- vals[[n]]
        }
        x$positivity_source <- "HHS j8mb-icvb; weekly positive/(positive+negative) PCR results"
        x$positivity_denominator <- "positive + negative PCR results; inconclusive excluded"
      }
      out[[i]] <- x
    }
    do.call(rbind,out)
  }
  list(county=make("county"),state=if(include_states) make("state") else data.frame())
}

bg_coverage <- function(d) {
  if (!nrow(d)) return(data.frame())
  do.call(rbind,lapply(split(d,d$geography_id),function(x) data.frame(
    geography_level=x$geography_level[1],geography_id=x$geography_id[1],geography_name=x$geography_name[1],
    site_ids=x$site_ids[1],n_sites=x$n_sites[1],weeks=nrow(x),
    weeks_with_cases=sum(is.finite(x$new_reported_cases)),weeks_with_positivity=sum(is.finite(x$positivity_pct)),
    positivity_coverage_pct=100*mean(is.finite(x$positivity_pct)),
    negative_case_revision_weeks=sum(x$new_reported_cases<0,na.rm=TRUE),stringsAsFactors=FALSE)))
}

bg_site_weekly <- function(sites, weeks, county, state) {
  d <- sites[rep(seq_len(nrow(sites)),each=nrow(weeks)),,drop=FALSE]
  d$week_start <- rep(weeks$week_start,nrow(sites)); d$week_end <- rep(weeks$week_end,nrow(sites))
  ck <- match(paste(d$county_fips,d$week_end),paste(county$geography_id,county$week_end))
  sk <- match(paste(d$state_abbr,d$week_end),paste(state$geography_id,state$week_end))
  get <- function(x, column, ix) if(nrow(x)) x[[column]][ix] else rep(NA_real_,length(ix))
  d$county_cases <- get(county,"new_reported_cases",ck); d$county_positivity_pct <- get(county,"positivity_pct",ck)
  d$state_cases <- get(state,"new_reported_cases",sk); d$state_positivity_pct <- get(state,"positivity_pct",sk)
  cgood <- is.finite(d$county_positivity_pct); sgood <- is.finite(d$state_positivity_pct)
  d$best_available_positivity_pct <- ifelse(cgood,d$county_positivity_pct,d$state_positivity_pct)
  d$positivity_geography_used <- ifelse(cgood,"county",ifelse(sgood,"state","missing"))
  d$positivity_geography_id <- ifelse(cgood,d$county_fips,ifelse(sgood,d$state_abbr,NA_character_))
  d$cases_at_selected_geography <- ifelse(cgood,d$county_cases,ifelse(sgood,d$state_cases,NA_real_))
  rownames(d) <- NULL; d
}

bg_plot <- function(d) {
  d <- d[order(d$week_start),]; x <- as.numeric(d$week_start)
  cases <- d$new_reported_cases; pct <- d$positivity_pct
  cp <- pretty(range(c(0,cases[is.finite(cases)]))); cy <- range(cp)
  if(diff(cy)==0) cy <- c(0,1)
  py <- c(0,max(5,ceiling(max(c(0,pct[is.finite(pct)]))/5)*5))
  old <- par(mar=c(7.7,5.5,6.6,5.4), mgp=c(3.2,0.7,0), tcl=-0.25, las=1,
             family="sans", cex=0.95, xaxs="i", yaxs="i")
  on.exit(par(old))
  plot(x,cases,type="n",ylim=cy,xlim=range(x)+c(-4,4),axes=FALSE,xlab="",ylab="")
  abline(h=cp,col="#E5E9EE",lwd=0.7)
  rect(x-2.8,pmin(0,cases),x+2.8,pmax(0,cases),col="#B7C6D7",border=NA)
  abline(h=0,col="#697785",lwd=0.7)
  axis(2,at=cp,labels=format(cp,trim=TRUE,scientific=FALSE,big.mark=","),cex.axis=0.85)
  ticks <- seq(as.Date(format(min(d$week_start),"%Y-%m-01")),
               as.Date(format(max(d$week_end),"%Y-%m-01")),by="month")
  ticks <- ticks[as.numeric(ticks)>=min(x)-4 & as.numeric(ticks)<=max(x)+4]
  if(!length(ticks)) ticks <- d$week_start[1L]
  ticks <- ticks[seq(1,length(ticks),by=max(1,ceiling(length(ticks)/12)))]
  labs <- paste(month.abb[as.integer(format(ticks,"%m"))],format(ticks,"%Y"),sep="\n")
  axis(1,at=as.numeric(ticks),labels=labs,cex.axis=0.85)
  mtext("New reported cases (weekly net change)",side=2,line=4.2,las=0,cex=0.9)
  par(new=TRUE)
  plot(x,pct,type="n",ylim=py,xlim=range(x)+c(-4,4),axes=FALSE,xlab="",ylab="")
  lines(x,pct,col="#AE332B",lwd=2.2)
  points(x,pct,col="#AE332B",pch=16,cex=0.45)
  axis(4,at=pretty(py),labels=paste0(pretty(py),"%"),col.axis="#AE332B",cex.axis=0.85)
  mtext("Test positivity (%)",side=4,line=3.6,las=0,col="#AE332B",cex=0.9)
  missing <- !is.finite(pct)
  if(any(missing)) points(x[missing],rep(py[2]*0.025,sum(missing)),pch=4,col="#AE332B",cex=0.65)
  mtext(d$geography_name[1],side=3,line=4.7,adj=0,font=2,cex=1.2)
  mtext(paste0(toupper(substr(d$geography_level[1],1,1)),substring(d$geography_level[1],2),
               " background | ",min(d$week_start)," to ",max(d$week_end)),side=3,line=3.3,adj=0,cex=0.93)
  sid <- strwrap(paste("Trial sites:",d$site_ids[1]),width=110)
  mtext(sid,side=3,line=seq(2.0,by=-0.8,length.out=length(sid)),adj=0,cex=0.72,col="#45515E")
  legend("topleft",legend=c("Reported cases (left axis)","Test positivity (right axis)"),
    pch=c(15,16),pt.cex=c(1.7,0.6),lty=c(0,1),lwd=c(0,2.2),
    col=c("#B7C6D7","#AE332B"),bty="n",cex=0.85,bg="white")
  mtext("Monday-Sunday weeks; each column and positivity point is positioned at its week's Monday. Month labels mark the first of the month.",
    side=1,line=3.1,adj=0,cex=0.7)
  mtext(paste0("Positivity available: ",sum(!missing),"/",nrow(d)," weeks. Missing values are gaps (x), not zero. Axes vary by chart."),
    side=1,line=4.1,adj=0,cex=0.7)
  note <- if(d$geography_level[1]=="county")
    "Sources: JHU CSSE reported cases; CDC nra9-vzzn, 7-day NAAT positivity. Constituent county test counts unavailable." else
    "Sources: JHU CSSE state case totals; HHS j8mb-icvb PCR tests. Positivity = positive / (positive + negative); inconclusive excluded."
  mtext(note,side=1,line=5.1,adj=0,cex=0.66)
  mtext("Retrospective surveillance: reporting/revision effects apply. Site geography does not establish participants' exposure or residence.",
    side=1,line=6.1,adj=0,cex=0.66,col="#45515E")
  invisible(NULL)
}

bg_charts <- function(d, output, level, pngs) {
  if (!nrow(d)) return(character())
  areas <- unique(d$geography_id); paths <- character(length(areas))
  pdf_file <- file.path(output,paste0(level,"_charts.pdf"))
  if(capabilities("cairo")) {
    grDevices::cairo_pdf(pdf_file,width=11.7,height=8.3,onefile=TRUE,family="sans",fallback_resolution=300)
  } else {
    grDevices::pdf(pdf_file,width=11.7,height=8.3,onefile=TRUE,useDingbats=FALSE,
      title=paste("Community COVID-19 activity:",level),paper="special")
  }
  tryCatch(for(id in areas) bg_plot(d[d$geography_id==id,]),finally=grDevices::dev.off())
  if(pngs) {
    dir.create(file.path(output,"charts"),showWarnings=FALSE)
    for(i in seq_along(areas)) {
      paths[i] <- paste0("charts/",level,"_",areas[i],".png")
      args <- list(filename=file.path(output,paths[i]),width=11.7,height=8.3,units="in",res=300)
      if(capabilities("cairo")) args$type <- "cairo"
      do.call(grDevices::png,args)
      tryCatch(bg_plot(d[d$geography_id==areas[i],]),finally=grDevices::dev.off())
    }
  }
  setNames(paths,areas)
}

bg_report <- function(output, sites, weeks, coverage, png_paths, problems) {
  n_good <- sum(!is.na(sites$county_fips))
  lines <- c('<!doctype html><html lang="en"><meta charset="utf-8"><title>Background COVID-19</title>',
    '<style>body{font:16px/1.5 Arial,sans-serif;color:#253344;max-width:1100px;margin:32px auto;padding:0 24px}h1{font-size:30px}h2{margin-top:38px}table{border-collapse:collapse;width:100%;font-size:14px}th,td{border-bottom:1px solid #dce3e9;padding:8px;text-align:left}th{background:#edf2f6}img{max-width:100%}.note{background:#f2f5f8;padding:14px}a{color:#1b5682}@media print{body{margin:0}article{break-before:page}a{color:inherit;text-decoration:none}.navigation{display:none}}</style>',
    '<h1>Community COVID-19 activity around the trial sites</h1>',
    paste0('<p>',nrow(sites),' sites; ',n_good,' matched to a county; ',length(unique(na.omit(sites$county_fips))),
      ' counties; ',length(unique(na.omit(sites$state_abbr))),' states. ',nrow(weeks),' complete Monday-Sunday weeks: ',
      min(weeks$week_start),' to ',max(weeks$week_end),'.</p>'),
    '<p class="navigation"><a href="county_charts.pdf">County charts (print PDF)</a> | <a href="state_charts.pdf">State charts (print PDF, if enabled)</a> | <a href="site_weekly.csv">Site-by-week CSV</a> | <a href="coverage.csv">Coverage CSV</a></p>',
    '<div class="note"><strong>Interpretation.</strong> Bars show weekly changes in reported case totals. Lines show test positivity, not infection prevalence. County and state charts are separate; their positivity definitions and geographic units differ. Geographic background is descriptive and cannot establish trial testing bias or an individual participant\'s infection.</div>',
    '<h2>Definitions and reproducibility</h2><ul>',
    '<li>Site latitude/longitude are matched to 2020 Census county/state boundaries. This is the county of the site address, not a verified recruitment catchment or participants\' home county. Coordinate errors in the input asset propagate into the assignment.</li>',
    '<li>County positivity: CDC historical dataset nra9-vzzn, as revised in the archive, not an unrevised contemporaneous snapshot. Use the published 7-day percentage on the Sunday endpoint. Daily rolling percentages are neither summed nor averaged. Positive and total test counts are not supplied in this county dataset; pooled county positivity and sampling confidence intervals cannot be calculated from it alone.</li>',
    '<li>Case counts: JHU CSSE archived US time series, cumulative Sunday value minus the preceding Sunday value. These are net changes in reported cases, not symptom-onset dates. Reporting backlogs, corrections and jurisdiction-specific inclusion of probable cases apply. Negative net revisions remain negative; they are not replaced with zero. State totals include all JHU rows for that state, including geographically unassigned/out-of-state categories.</li>',
    '<li>State positivity: sum positive PCR results / sum positive plus negative PCR results across the seven days, times 100. Inconclusive tests are retained separately and excluded. All seven dates for both positive and negative results are required; missing categories are not treated as zero. Negative daily test-count corrections are flagged and the affected percentage is left missing.</li>',
    '<li>The county CDC series refers to diagnostic/screening NAAT tests, while the state HHS series is labelled PCR testing. HHS dates generally refer to test completion/result reporting, with other dates substituted when unavailable. Cases and tests therefore need not have the same event-date basis or the same numerator.</li>',
    '<li>No interpolation or substitution is made in the plotted county series. State data appear only on separate charts. In site_weekly.csv the best_available_positivity_pct field selects county then state, with explicit geography and a case count at that same selected level. Use these provenance fields when analysing it.</li>',
    '<li>Several sites can share a county. The repeated area counts in site_weekly.csv must NOT be summed across sites. No trial-wide positivity is created by an unweighted average of counties or by counting shared counties more than once.</li>',
    '<li>Public testing includes people tested for different reasons, including screening. It is not directly comparable to positivity among symptomatic trial participants. Dual chart axes have independent scales; bar/line heights and crossings have no numerical equivalence.</li>',
    '<li>Only complete weeks are shown. Source URLs, retrieval timestamps and file checksums are in source_manifest.csv. Cached files and their metadata support reruns; REFRESH=TRUE requests new snapshots. API/download problems are listed in download_issues.csv.</li></ul>',
    '<p>Sources: <a href="https://data.cdc.gov/d/nra9-vzzn">CDC county archive</a>; <a href="https://healthdata.gov/d/j8mb-icvb">HHS state PCR results</a>; <a href="https://github.com/CSSEGISandData/COVID-19">JHU CSSE</a>; <a href="https://geocoding.geo.census.gov/geocoder/Geocoding_Services_API.html">US Census geocoder</a>.</p>',
    '<h2>Geographic coverage</h2><table><tr><th>Level</th><th>Area</th><th>Sites</th><th>Cases: weeks</th><th>Positivity: weeks</th></tr>')
  for(i in seq_len(nrow(coverage))) {
    z <- coverage[i,]; key <- paste(z$geography_level,z$geography_id,sep="_")
    lines <- c(lines,paste0('<tr><td>',z$geography_level,'</td><td><a href="#',key,'">',bg_html(z$geography_name),
      '</a></td><td>',z$n_sites,'</td><td>',z$weeks_with_cases,'/',z$weeks,'</td><td>',z$weeks_with_positivity,'/',z$weeks,'</td></tr>'))
  }
  lines <- c(lines,'</table>')
  if(nrow(problems)) lines <- c(lines,'<h2>Download/processing issues</h2><pre>',
    bg_html(paste(apply(problems,1,paste,collapse=" | "),collapse="\n")),'</pre>')
  for(i in seq_len(nrow(coverage))) {
    z <- coverage[i,]; key <- paste(z$geography_level,z$geography_id,sep="_")
    if(key %in% names(png_paths) && nzchar(png_paths[[key]])) lines <- c(lines,
      paste0('<article id="',key,'"><h2>',bg_html(z$geography_name),'</h2><img loading="lazy" src="',
             png_paths[[key]],'" alt="Weekly cases and positivity: ',bg_html(z$geography_name),'"></article>'))
  }
  writeLines(c(lines,'</html>'),file.path(output,"background_covid.html"),useBytes=TRUE)
}

analyse_background_covid <- function(root=PHMPT_ROOT,sites_file=file.path(root,"clinical_trial_sites.csv"),
    start_date=PHMPT_START_DATE,end_date=PHMPT_END_DATE,output_dir=file.path(root,"analysis","background_covid"),
    cache_dir=file.path(root,"cache","background_covid"),include_states=PHMPT_INCLUDE_STATES,refresh=PHMPT_REFRESH,
    site_ids=PHMPT_SITE_IDS,write_png=PHMPT_WRITE_PNG,excel_fr=PHMPT_EXCEL_FR,
    show_report=PHMPT_SHOW_REPORT) {
  if(!requireNamespace("jsonlite",quietly=TRUE)) {
    message("Installing the JSON reader package jsonlite...")
    utils::install.packages("jsonlite",repos="https://cloud.r-project.org")
    if(!requireNamespace("jsonlite",quietly=TRUE)) stop("Install package jsonlite, then rerun this script.")
  }
  old_timeout <- getOption("timeout"); options(timeout=max(120,old_timeout))
  on.exit(options(timeout=old_timeout),add=TRUE)
  start <- as.Date(start_date); end <- as.Date(end_date)
  if(anyNA(c(start,end)) || end<start) stop("Invalid start/end dates.")
  first <- bg_monday(start); if(first<start) first <- first+7L
  last <- bg_monday(end)-1L; if(end==bg_monday(end)+6L) last <- end
  if(last<first+6L) stop("The selected date range contains no complete Monday-Sunday week.")
  weeks <- data.frame(week_start=seq(first,last-6L,by="7 days"))
  weeks$week_end <- weeks$week_start+6L
  if(first!=start || last!=end) message("Partial boundary weeks omitted; chart window: ",first," to ",last)
  if(first<as.Date("2020-01-27") || last>as.Date("2022-10-16"))
    stop("This historical CDC source supports complete weeks from 2020-01-27 through 2022-10-16. Choose a range within it.")
  if(!file.exists(sites_file)) stop("Cannot find ",sites_file,". Set PHMPT_ROOT or PHMPT_SITES_FILE to your local project.")
  sites <- bg_read_csv(sites_file)
  if(ncol(sites)==1L) sites <- utils::read.csv2(sites_file,colClasses="character",check.names=FALSE,
    stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
  bg_require(sites,c("SITEID","latitude","longitude"),"clinical_trial_sites.csv")
  sites$SITEID <- trimws(sites$SITEID)
  if(anyNA(sites$SITEID) || any(!nzchar(sites$SITEID)) || anyDuplicated(sites$SITEID)) stop("SITEID must be nonempty and unique.")
  if(length(site_ids)) {
    if(length(setdiff(site_ids,sites$SITEID))) stop("Requested SITEID absent from clinical_trial_sites.csv.")
    sites <- sites[sites$SITEID %in% site_ids,,drop=FALSE]
  }
  sites$latitude <- bg_num(gsub(",",".",sites$latitude,fixed=TRUE))
  sites$longitude <- bg_num(gsub(",",".",sites$longitude,fixed=TRUE))
  if(!nrow(sites) || any(!is.finite(sites$latitude) | abs(sites$latitude)>90 |
                       !is.finite(sites$longitude) | abs(sites$longitude)>180)) stop("Missing/invalid site coordinates.")
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
  ctx <- new.env(parent=emptyenv()); ctx$cache <- cache_dir;ctx$refresh <- refresh
  ctx$manifest <- list();ctx$problems <- list()
  persist_logs <- function() {
    if(length(ctx$manifest)) bg_csv(unique(do.call(rbind,ctx$manifest)),file.path(output_dir,"source_manifest.csv"),excel_fr)
    p <- if(length(ctx$problems)) do.call(rbind,ctx$problems) else data.frame(stage=character(),item=character(),message=character())
    bg_csv(p,file.path(output_dir,"download_issues.csv"),excel_fr); p
  }
  on.exit(persist_logs(),add=TRUE)
  sites <- bg_geocode(sites,ctx)
  bg_csv(sites,file.path(output_dir,"sites_geography.csv"),excel_fr)
  if(!any(!is.na(sites$county_fips))) stop("No sites could be assigned to counties. See download_issues.csv.")
  county <- bg_county_data(unique(na.omit(sites$county_fips)),first,last,ctx)
  state <- if(include_states) bg_state_data(unique(na.omit(sites$state_abbr)),first,last,ctx) else data.frame()
  message("[Cases] Reading the archived JHU CSSE US series")
  jurl <- paste0("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/",
    "csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv")
  jp <- bg_fetch(jurl,"jhu_confirmed_US.csv",ctx,bg_csv_check(c("FIPS","Province_State","Country_Region")))
  jhu <- bg_read_csv(jp);jhu$FIPS <- bg_fips(jhu$FIPS)
  weekly <- bg_weekly(sites,weeks,jhu,county,state,include_states)
  coverage <- rbind(bg_coverage(weekly$county),bg_coverage(weekly$state))
  site_weekly <- bg_site_weekly(sites,weeks,weekly$county,weekly$state)
  bg_csv(weekly$county,file.path(output_dir,"county_weekly.csv"),excel_fr)
  if(include_states) bg_csv(weekly$state,file.path(output_dir,"state_weekly.csv"),excel_fr)
  bg_csv(site_weekly,file.path(output_dir,"site_weekly.csv"),excel_fr)
  bg_csv(coverage,file.path(output_dir,"coverage.csv"),excel_fr)
  params <- data.frame(parameter=c("script_version","run_utc","requested_start","requested_end","first_week_start",
    "last_week_end","sites_input","sites_input_md5","include_states","census_benchmark","census_vintage"),
    value=c("1.0.0",format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC"),as.character(start),as.character(end),
    as.character(first),as.character(last),normalizePath(sites_file,winslash="/"),unname(tools::md5sum(sites_file)),
    as.character(include_states),"Public_AR_Census2020","Census2020_Census2020"))
  bg_csv(params,file.path(output_dir,"run_parameters.csv"))
  message("[Charts] Writing county and state print PDFs",if(write_png) " and 300 dpi PNGs" else "")
  cp <- bg_charts(weekly$county,output_dir,"county",write_png)
  sp <- bg_charts(weekly$state,output_dir,"state",write_png)
  if(length(cp)) names(cp) <- paste0("county_",names(cp))
  if(length(sp)) names(sp) <- paste0("state_",names(sp))
  problems <- persist_logs()
  bg_report(output_dir,sites,weeks,coverage,c(cp,sp),problems)
  message("Done: ",nrow(sites)," sites; ",length(unique(na.omit(sites$county_fips)))," counties; ",
          nrow(weeks)," weeks. County positivity: ",sum(is.finite(weekly$county$positivity_pct)),"/",nrow(weekly$county)," area-weeks.")
  message("Report: ",normalizePath(file.path(output_dir,"background_covid.html"),winslash="/"))
  if(any(coverage$weeks_with_positivity<coverage$weeks)) message("Some positivity values are unavailable; inspect coverage.csv.")
  if(show_report) utils::browseURL(normalizePath(file.path(output_dir,"background_covid.html"),winslash="/"))
  invisible(list(sites=sites,weeks=weeks,county=weekly$county,state=weekly$state,
                 site_weekly=site_weekly,coverage=coverage,problems=problems))
}

# ---- RUN: starts automatically when the whole script is executed ----
if(isTRUE(getOption("phmpt.autorun",TRUE))) {
  PHMPT_BACKGROUND_RESULT <- analyse_background_covid(root=PHMPT_ROOT,sites_file=PHMPT_SITES_FILE,
    start_date=PHMPT_START_DATE,end_date=PHMPT_END_DATE,output_dir=PHMPT_OUTPUT_DIR,
    cache_dir=PHMPT_CACHE_DIR,include_states=PHMPT_INCLUDE_STATES,refresh=PHMPT_REFRESH,
    site_ids=PHMPT_SITE_IDS,write_png=PHMPT_WRITE_PNG,excel_fr=PHMPT_EXCEL_FR,
    show_report=PHMPT_SHOW_REPORT)
}
