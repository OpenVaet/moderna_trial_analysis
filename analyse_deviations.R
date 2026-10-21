# Moderna P301 protocol deviations -- v1.3.0 -- R 4.3+
# Standalone: save in the project root, Select All / Run or Source in RStudio.
# Reads your local PHMPT files. No downloads and no earlier analysis scripts.
# Base R + grid suffice for CSV/JMP and PDF/HTML output; haven reads XPT/SAS7BDAT.
# Tests compare participants with >=1 recorded deviation, not repeated records.
# Primary grouping: exact DVTERM + DVTERM1 (whitespace normalized only).
# Weekly charts count records for seven staff-related types, across ALL sites.

# ---- SETTINGS ----
PHMPT_ROOT <- ""                         # Auto-detect, or D:/Projects/moderna_trial_analysis
PHMPT_RELEASES <- c("EUA", "Prod")         # Missing releases are documented and skipped
PHMPT_ADSL_FILES <- c(EUA = "", Prod = "") # Optional explicit paths, relative to root or absolute
PHMPT_DV_FILES <- c(EUA = "", Prod = "")
PHMPT_ADDV_FILES <- c(EUA = "", Prod = "")
PHMPT_OUTPUT_DIR <- ""                   # New dated folder under analysis/
PHMPT_EUA_CUTOFF <- "2020-11-25"          # Data-analysis cutoff, not authorization date
PHMPT_BLA_CUTOFF <- "2021-03-26"
PHMPT_START <- "randomization"           # Or "first_dose" for a post-dose analysis
PHMPT_ALPHA <- 0.05
PHMPT_EXPECTED_MIN <- 5                  # Pearson only if ALL four expected cells >=5
PHMPT_MIN_RANDOMIZED <- 30000L           # Guard against accidentally using a subject subset
PHMPT_CSV_ENCODING <- "auto"             # UTF-8-BOM, then windows-1252 on decoding failure
PHMPT_NUMERIC_DATE_ORIGIN <- ""          # ADSL CSV numeric dates: set only if origin is documented
PHMPT_EXCEL_FR <- TRUE                   # Additional semicolon-delimited summary CSV
PHMPT_SITE_MIN_RECORDS <- 20L            # Combined arms, within each site/window/comparison
PHMPT_WEEKLY_PNG_DPI <- 300L             # A4 landscape PDF + one PNG per chart page
PHMPT_UNBLINDING_START <- "2020-12-14"   # EUA-era calendar reference, not first individual unblinding
# Earlier UBLDDT dates are retained and audited. Scheduled open-label visits began
# 2020-12-24; see weekly_staff_deviations_methods.txt for these distinct milestones.

# ---- INPUT AND VALIDATION ----
dv_root <- function(root = "") {
  if (nzchar(root)) return(normalizePath(root, winslash = "/", mustWork = TRUE))
  candidates <- getwd()
  for (f in sys.frames()) if (!is.null(f$ofile)) candidates <- c(candidates, dirname(f$ofile))
  a <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(a)) candidates <- c(candidates, dirname(sub("^--file=", "", a[1])))
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p)) candidates <- c(candidates, dirname(p))
  }
  for (p in unique(candidates[dir.exists(candidates)])) {
    p <- normalizePath(p, winslash = "/", mustWork = TRUE)
    for (i in seq_len(8L)) {
      if (dir.exists(file.path(p,"prod_data")) || dir.exists(file.path(p,"eua_data"))) return(p)
      parent <- dirname(p); if (parent == p) break; p <- parent
    }
  }
  stop("Cannot find prod_data/eua_data. Set PHMPT_ROOT to your extracted project root.", call. = FALSE)
}
dv_path <- function(root,p) {
  if (grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)",p)) p else file.path(root,p)
}
dv_csv <- function(x,path,sep=",") {
  con <- file(path,"wb"); on.exit(close(con)); writeBin(as.raw(c(239,187,191)),con)
  write.table(x,con,sep=sep,dec=if(sep==";")"," else ".",row.names=FALSE,col.names=TRUE,quote=TRUE,na="",qmethod="double")
}
dv_text <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""
  bad <- !validUTF8(x)
  if (any(bad)) x[bad] <- iconv(x[bad],from="windows-1252",to="UTF-8",sub=NA)
  if (anyNA(x)) stop("Text cannot be decoded losslessly. Specify the source encoding.")
  enc2utf8(x)
}
dv_squish <- function(x) trimws(gsub("[[:space:]\\p{Z}]+"," ",dv_text(x),perl=TRUE))
dv_need <- function(x,names) {
  absent <- setdiff(names,colnames(x))
  if (length(absent)) stop("Missing columns: ",paste(absent,collapse=", "))
}
dv_fields <- function(kind) switch(kind,
  ADSL=c("STUDYID","USUBJID","SITEID","RANDFL","TRT01P","RANDDT","TR01SDT","CUTOFFDT",
    "AP01EDT","AP02SDT","UBLDDT","EOSDT","DTHDT","SEX","CHILDBP"),
  DV=c("STUDYID","USUBJID","DVSEQ","DVTERM","DVSTDTC","DVCAT","DVSCAT","DVSPID","EPOCH"),
  ADDV=c("STUDYID","USUBJID","DVSEQ","DVTERM","DVTERM1","DVSTDTC","DVCAT","DVSCAT","DVSIG","TRT01P","RANDFL"))
dv_read <- function(path,kind,encoding="auto") {
  ext <- tolower(tools::file_ext(path))
  if (ext == "jmp") {
    if (kind != "ADSL") stop("This embedded JMP reader supports P301 ADSL only; export DV/ADDV JMP to CSV first.")
    x <- dv_read_jmp(path, dv_fields("ADSL"))
  } else if (ext == "csv") {
    encodings <- if (encoding == "auto") c("UTF-8-BOM","windows-1252") else encoding
    x <- NULL; errors <- character()
    for (enc in encodings) {
      y <- tryCatch(withCallingHandlers(read.csv(path,colClasses="character",na.strings=c("","NA"),
        check.names=FALSE,fileEncoding=enc),warning=function(w)stop(conditionMessage(w))),error=function(e)e)
      if (!inherits(y,"error")) {x <- y; break}
      errors <- c(errors,conditionMessage(y))
    }
    if (is.null(x)) stop("CSV decoding failed: ",paste(errors,collapse="; "))
  } else if (ext %in% c("xpt","sas7bdat")) {
    if (!requireNamespace("haven",quietly=TRUE)) stop('Install haven for XPT/SAS7BDAT: install.packages("haven")')
    x <- as.data.frame(if(ext=="xpt") haven::read_xpt(path,.name_repair="minimal") else haven::read_sas(path,.name_repair="minimal"))
    x[] <- lapply(x,function(z) if(inherits(z,"Date"))format(z,"%Y-%m-%d") else
      if(inherits(z,"POSIXt"))format(z,"%Y-%m-%d",tz="UTC") else dv_text(z))
  } else stop("Unsupported input format: ",ext)
  names(x) <- sub("^\ufeff","",names(x))
  if(anyDuplicated(names(x))) stop("Duplicate input column names.")
  x[] <- lapply(x,dv_text)
  required <- if(kind=="ADSL") c("STUDYID","USUBJID","RANDFL","TRT01P","RANDDT","CUTOFFDT") else
    c("STUDYID","USUBJID","DVSEQ","DVTERM","DVSTDTC")
  dv_need(x,required)
  x <- x[intersect(dv_fields(kind),names(x))]
  for(n in setdiff(dv_fields(kind),names(x)))x[[n]] <- ""
  x <- x[dv_fields(kind)]
  for(n in c("STUDYID","USUBJID"))x[[n]] <- trimws(x[[n]])
  if(!nrow(x) || any(!nzchar(x$USUBJID)) || any(toupper(x$STUDYID)!="MRNA-1273-P301"))
    stop("Empty input, missing subject IDs, or an unexpected study (requires mRNA-1273-P301).")
  if(kind!="ADSL") {
    seq <- suppressWarnings(as.numeric(x$DVSEQ))
    if(any(!is.finite(seq)|seq<0|seq!=floor(seq)))stop("Invalid DVSEQ.")
    x$DVSEQ <- format(seq,scientific=FALSE,trim=TRUE)
  }
  key <- if(kind=="ADSL")x$USUBJID else paste(x$USUBJID,x$DVSEQ,sep="|")
  if(anyDuplicated(key))stop("Duplicate subject/event keys; source must be reconciled before analysis.")
  rownames(x) <- NULL; x
}
dv_date <- function(x,field,origin="",strict=FALSE) {
  raw <- trimws(dv_text(x)); d <- as.Date(rep(NA_character_,length(raw)))
  iso <- grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}($|T)",raw)
  d[iso] <- suppressWarnings(as.Date(substr(raw[iso],1L,10L),"%Y-%m-%d"))
  invalid <- iso & (is.na(d) | format(d,"%Y-%m-%d") != substr(raw,1L,10L))
  invalid[is.na(invalid)] <- FALSE; d[invalid] <- as.Date(NA)
  num <- nzchar(raw) & grepl("^-?[0-9]+([.][0-9]+)?$",raw)
  if(any(num)) {
    if(!nzchar(origin))stop(field,": numeric dates found. Set PHMPT_NUMERIC_DATE_ORIGIN only from documented metadata; no origin is guessed.")
    z <- as.numeric(raw[num]); if(any(z!=floor(z)))stop(field,": numeric date includes a fractional day.")
    d[num] <- as.Date(z,origin=origin)
  }
  if(strict && any(nzchar(raw)&is.na(d)))stop(field,": incomplete or invalid date.")
  d
}
dv_normalize <- function(x,kind,release,cutoff,min_randomized,origin) {
  if(kind=="ADSL") {
    x$RANDFL <- toupper(trimws(x$RANDFL)); x$TRT01P <- trimws(x$TRT01P)
    x$TRT01P[toupper(x$TRT01P)=="PLACEBO"] <- "Placebo"
    x$TRT01P[toupper(x$TRT01P)=="MRNA-1273"] <- "mRNA-1273"
    if(any(!x$RANDFL %in% c("Y","N")))stop("Unexpected RANDFL coding.")
    if(any(x$RANDFL=="Y" & !x$TRT01P %in% c("mRNA-1273","Placebo")))stop("Unexpected randomized arm in TRT01P.")
    if(sum(x$RANDFL=="Y") < min_randomized)stop("ADSL appears to be a subject subset: ",sum(x$RANDFL=="Y")," randomized records.")
    for(n in c("RANDDT","TR01SDT","CUTOFFDT","AP01EDT","AP02SDT","UBLDDT","EOSDT","DTHDT"))
      x[[n]] <- format(dv_date(x[[n]],n,origin,strict=TRUE),"%Y-%m-%d")
    if(anyNA(x$CUTOFFDT) || any(x$CUTOFFDT!=cutoff))stop("ADSL CUTOFFDT must equal ",cutoff," for ",release,".")
    if(any(x$RANDFL=="Y" & is.na(x$RANDDT)))stop("Randomized participants without RANDDT.")
    if(release=="Prod" && all(is.na(x$AP01EDT)) && all(is.na(x$AP02SDT)) && all(is.na(x$UBLDDT)))
      stop("Production ADSL lacks blinded-period endpoints. Use the full BLA ADSL, not an earlier selected-column export.")
  }
  x[] <- lapply(x,function(z){z[is.na(z)]<-"";z})
  x
}
dv_candidates <- function(root,release,kind) {
  dir <- file.path(root,if(release=="Prod")"prod_data" else "eua_data")
  files <- list.files(dir,recursive=TRUE,full.names=TRUE,ignore.case=TRUE,
    pattern="\\.(csv|xpt|sas7bdat|jmp)$")
  files <- unique(files[grepl(paste0("(^|[_ -])",tolower(kind),"([_. -]|$)"),basename(files),ignore.case=TRUE)])
  files <- files[!grepl("p201|p101|dmid",basename(files),ignore.case=TRUE)]
  pr <- match(tolower(tools::file_ext(files)),c("jmp","xpt","sas7bdat","csv"))
  files[order(pr,grepl("[/\\\\]drops_data[/\\\\]",files),files)]
}
dv_select <- function(root,release,kind,explicit,out,cutoff,min_randomized,encoding,origin,optional=FALSE) {
  paths <- if(nzchar(explicit))dv_path(root,explicit) else dv_candidates(root,release,kind)
  audit <- data.frame(release=character(),kind=character(),path=character(),md5=character(),rows=integer(),
    status=character(),detail=character(),stringsAsFactors=FALSE)
  hashes <- character(); valid <- list(); signatures <- character()
  log <- file.path(out,paste0("input_discovery_",release,"_",kind,".csv"))
  if(!length(paths) && optional) {dv_csv(audit,log);return(NULL)}
  for(p in paths) {
    message("Checking ",release," ",kind,": ",basename(p))
    hash <- if(file.exists(p))unname(tools::md5sum(p)) else ""
    same <- if(nzchar(hash))match(hash,hashes) else NA_integer_
    if(!is.na(same)) {
      audit <- rbind(audit,data.frame(release=release,kind=kind,path=p,md5=hash,rows=audit$rows[same],
        status="IDENTICAL_COPY",detail=paste("Same bytes as",audit$path[same])));hashes<-c(hashes,hash);next
    }
    result <- tryCatch({
      x <- dv_normalize(dv_read(p,kind,encoding),kind,release,cutoff,min_randomized,origin)
      tmp <- tempfile(); y <- x[do.call(order,x[intersect(c("USUBJID","DVSEQ"),names(x))]),,drop=FALSE];rownames(y)<-NULL
      saveRDS(y,tmp,version=2,compress=FALSE);signature<-unname(tools::md5sum(tmp));unlink(tmp)
      list(path=normalizePath(p,winslash="/"),data=x,md5=hash,signature=signature)
    },error=function(e)e)
    ok <- !inherits(result,"error")
    audit <- rbind(audit,data.frame(release=release,kind=kind,path=p,md5=hash,rows=if(ok)nrow(result$data) else NA_integer_,
      status=if(ok)"VALID" else "REJECTED",detail=if(ok)"" else conditionMessage(result)))
    hashes <- c(hashes,hash)
    if(ok){valid[[length(valid)+1L]]<-result;signatures<-c(signatures,result$signature)}
  }
  dv_csv(audit,log)
  if(!length(valid))stop("No usable ",release," ",kind," input. See ",log,". ",paste(unique(audit$detail),collapse=" | "))
  if(length(unique(signatures))!=1L)stop("Conflicting ",release," ",kind," inputs. Set PHMPT_",kind,"_FILES explicitly; see ",log)
  selected <- valid[[1L]]
  audit$status[audit$status=="VALID"] <- "EQUIVALENT_FORMAT"
  audit$status[audit$path==selected$path] <- "SELECTED";dv_csv(audit,log)
  selected
}

# ---- COHORT, SOURCE RECONCILIATION AND WINDOWS ----
dv_context <- function(adsl,start_basis,origin) {
  x <- adsl
  x$arm <- x$TRT01P
  for(k in c("RANDDT","TR01SDT","CUTOFFDT","AP01EDT","AP02SDT","UBLDDT","EOSDT","DTHDT"))
    x[[k]] <- dv_date(x[[k]],k,origin,strict=TRUE)
  x$start <- if(start_basis=="first_dose")x$TR01SDT else x$RANDDT
  x$end <- x$CUTOFFDT
  x$end_basis <- "CUTOFFDT"
  for(k in c("AP01EDT","AP02SDT","UBLDDT","EOSDT","DTHDT")) {
    cap <- x[[k]]; if(k %in% c("AP02SDT","UBLDDT"))cap<-cap-1L
    take <- !is.na(cap) & cap<x$end
    x$end[take]<-cap[take];x$end_basis[take]<-if(k %in% c("AP02SDT","UBLDDT"))paste0(k," minus 1 day") else k
  }
  x$period_endpoint_missing <- is.na(x$AP01EDT)&is.na(x$AP02SDT)&is.na(x$UBLDDT)
  x$eligible <- x$RANDFL=="Y" & !is.na(x$start) & x$end>=x$start
  x$eligible[is.na(x$eligible)]<-FALSE
  x$exclusion <- ifelse(x$RANDFL!="Y","NOT_RANDOMIZED",ifelse(is.na(x$start),"START_DATE_MISSING",
    ifelse(x$end<x$start,"END_BEFORE_START","")))
  x
}
dv_records <- function(raw,addv,adsl,out,release) {
  if(is.null(raw)&&is.null(addv))stop("Neither DV nor ADDV is available.")
  ad <- if(is.null(addv))NULL else addv$data
  x <- if(is.null(raw))ad else raw$data
  x$source_row <- seq_len(nrow(x));x$source_dataset <- if(is.null(raw))"ADDV" else "DV"
  key <- function(z)paste(z$STUDYID,z$USUBJID,z$DVSEQ,sep="|")
  x$event_key <- key(x)
  if(!is.null(raw) && !is.null(ad)) {
    if(any(!key(ad)%in%x$event_key))stop("ADDV contains events absent from DV. Select matching snapshots; no sources were combined.")
    j <- match(x$event_key,key(ad));hit <- !is.na(j)
    conflicts <- list()
    for(k in c("DVTERM","DVSTDTC","DVCAT","DVSCAT")) {
      bad <- which(hit & dv_squish(x[[k]])!=dv_squish(ad[[k]][j]))
      if(length(bad))conflicts[[k]]<-data.frame(event_key=x$event_key[bad],field=k,DV=x[[k]][bad],ADDV=ad[[k]][j[bad]])
    }
    if(length(conflicts)) {
      dv_csv(do.call(rbind,conflicts),file.path(out,paste0(release,"_source_conflicts.csv")))
      stop("DV/ADDV content conflicts; see source_conflicts.csv.")
    }
    x$DVTERM1 <- ad$DVTERM1[j];x$DVSIG <- ad$DVSIG[j];x$in_addv <- hit
  } else {
    if(!"DVTERM1"%in%names(x))x$DVTERM1<-""
    if(!"DVSIG"%in%names(x))x$DVSIG<-""
    x$in_addv <- !is.null(ad)
  }
  x$DVTERM1 <- dv_text(x$DVTERM1);x$DVSIG<-dv_text(x$DVSIG)
  x$term <- dv_squish(paste(x$DVTERM,x$DVTERM1))
  if(any(!nzchar(x$term)))stop("Deviation records without any description.")
  x$sponsor_category <- dv_squish(x$DVSCAT)
  x$sponsor_category[x$sponsor_category==""] <- "Unspecified sponsor category"
  x$deviation_date <- dv_date(x$DVSTDTC,"DVSTDTC",origin="")
  x$date_status <- ifelse(is.na(x$deviation_date),ifelse(nzchar(x$DVSTDTC),"INCOMPLETE_OR_INVALID","MISSING"),"COMPLETE")
  j <- match(x$USUBJID,adsl$USUBJID)
  x$in_adsl <- !is.na(j)
  x$arm <- adsl$TRT01P[j];x$SITEID <- adsl$SITEID[j]
  if(!is.null(ad)) {
    h <- match(ad$USUBJID,adsl$USUBJID)
    bad <- which(!is.na(h) & nzchar(ad$TRT01P) & ad$TRT01P!=adsl$TRT01P[h])
    if(length(bad))stop("ADSL and ADDV disagree on original randomized treatment.")
  }
  x
}

# ---- TWO-ARM STATISTICS ----
dv_compare <- function(a,nv,c,np,expected_min=5,force_fisher=FALSE) {
  if(any(!is.finite(c(a,nv,c,np)))||nv<=0||np<=0||a<0||c<0||a>nv||c>np)stop("Invalid 2x2 counts.")
  tab <- matrix(c(a,nv-a,c,np-c),nrow=2L,byrow=TRUE,
    dimnames=list(arm=c("mRNA-1273","Placebo"),outcome=c("Affected","Not affected")))
  expected <- outer(rowSums(tab),colSums(tab))/sum(tab)
  if(any(colSums(tab)==0)) {p<-1;method<-"No variation";statistic<-NA_real_} else
  if(!force_fisher && all(expected>=expected_min)) {
    test<-stats::chisq.test(tab,correct=FALSE);p<-test$p.value;method<-"Chi-square";statistic<-unname(test$statistic)
  } else {test<-stats::fisher.test(tab,alternative="two.sided");p<-test$p.value;method<-"Fisher exact";statistic<-NA_real_}
  wilson <- function(x,n) {
    z<-stats::qnorm(.975);p<-x/n;center<-(p+z*z/(2*n))/(1+z*z/n)
    half<-z*sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    c(max(0,center-half),min(1,center+half))
  }
  pv<-a/nv;pp<-c/np;v<-wilson(a,nv);p0<-wilson(c,np);d<-pv-pp
  lower<-max(-1,d-sqrt((pv-v[1])^2+(p0[2]-pp)^2))
  upper<-min(1,d+sqrt((v[2]-pv)^2+(pp-p0[1])^2))
  data.frame(vaccine_n=a,vaccine_N=nv,vaccine_percent=100*pv,placebo_n=c,placebo_N=np,placebo_percent=100*pp,
    difference_pp=100*d,lower95_pp=100*lower,upper95_pp=100*upper,
    prevalence_ratio=if(pp==0)if(pv==0)NA_real_ else Inf else pv/pp,
    test=method,min_expected=min(expected),chi_square=statistic,df=if(method=="Chi-square")1L else NA_integer_,
    p_raw=p,higher_arm=if(pv>pp)"mRNA-1273" else if(pp>pv)"Placebo" else "Equal")
}
dv_window <- function(events,context,release,window,cap,expected_min,out,unblinding_start="2020-12-14") {
  ctx <- context;ctx$window_end<-pmin(ctx$end,as.Date(cap))
  pop <- ctx[ctx$eligible & ctx$window_end>=ctx$start,,drop=FALSE]
  den <- table(factor(pop$arm,levels=c("mRNA-1273","Placebo")))
  if(any(den==0))stop("An arm has no eligible participants in ",window)
  j<-match(events$USUBJID,pop$USUBJID)
  reason<-rep("INCLUDED",nrow(events))
  reason[is.na(j)]<-"OUTSIDE_ELIGIBLE_COHORT"
  reason[!is.na(j)&is.na(events$deviation_date)]<-"DATE_UNRESOLVED"
  before<-!is.na(j)&!is.na(events$deviation_date)&events$deviation_date<pop$start[j]
  after<-!is.na(j)&!is.na(events$deviation_date)&events$deviation_date>pop$window_end[j]
  reason[which(before)]<-"BEFORE_START";reason[which(after)]<-"AFTER_WINDOW"
  records<-events;records$window<-window;records$release<-release;records$window_status<-reason
  records$window_start<-pop$start[j];records$window_end<-pop$window_end[j]
  dv_csv(records,file.path(out,paste0("records_",window,".csv")))
  ev<-events[reason=="INCLUDED",,drop=FALSE]
  pop$followup_days<-as.integer(pop$window_end-pop$start)+1L;pop$window<-window
  dv_csv(pop,file.path(out,paste0("population_",window,".csv")))
  result<-list();sites<-list()
  for(level in c("term","sponsor_category")) {
    for(term in sort(unique(events[[level]]))) {
      z<-ev[ev[[level]]==term,,drop=FALSE]
      a<-length(unique(z$USUBJID[z$arm=="mRNA-1273"]));c<-length(unique(z$USUBJID[z$arm=="Placebo"]))
      r<-dv_compare(a,unname(den[1]),c,unname(den[2]),expected_min)
      r<-cbind(data.frame(release=release,window=window,level=level,deviation=term),r)
      r$vaccine_records<-sum(z$arm=="mRNA-1273");r$placebo_records<-sum(z$arm=="Placebo")
      r$affected_sites<-length(unique(z$SITEID[nzchar(z$SITEID)]))
      result[[length(result)+1L]]<-r
      if(nrow(z) && level=="term") {
        u<-unique(z[c("USUBJID","SITEID","arm")]);counts<-aggregate(rep(1L,nrow(u)),u[c("SITEID","arm")],sum)
        names(counts)[3]<-"affected_participants";counts$release<-release;counts$window<-window;counts$deviation<-term
        counts$site_arm_N<-vapply(seq_len(nrow(counts)),function(i)sum(pop$SITEID==counts$SITEID[i]&pop$arm==counts$arm[i]),integer(1))
        sites[[length(sites)+1L]]<-counts
      }
    }
  }
  summary<-do.call(rbind,lapply(c("mRNA-1273","Placebo"),function(arm){
    p<-pop[pop$arm==arm,,drop=FALSE];e<-ev[ev$arm==arm,,drop=FALSE]
    data.frame(release=release,window=window,arm=arm,N=nrow(p),any_deviation=length(unique(e$USUBJID)),records=nrow(e),
      person_days=sum(p$followup_days),median_followup_days=median(p$followup_days),
      observed_first_date=as.character(min(p$start)),observed_last_date=as.character(max(p$window_end)),
      missing_blind_endpoint=sum(p$period_endpoint_missing))
  }))
  site_fisher <- dv_site_window(ev,pop,release,window)
  list(results=do.call(rbind,result),sites=if(length(sites))do.call(rbind,sites) else NULL,
    summary=summary,site_fisher=site_fisher,weekly=dv_weekly_window(records,pop,release,window,unblinding_start))
}

# ---- SITE FISHER TESTS: THE SEVEN DEVIATION TYPES FROM THE EARLIER AUDIT ----
dv_site_rules <- function() data.frame(
  category=c("OBSERVATION_UNDER_30","DIARY_30_MISSING","NP_SWAB_MISSING","PREGNANCY_TEST",
    "EXAM_VITALS_INCOMPLETE","SAE_AESI_REPORT_DELAY","UNBLINDED_BLINDED_TASKS"),
  label=c("Observation shorter than 30 minutes","Missing 30-minute diary","Required nasopharyngeal swab missing",
    "Pregnancy testing not per protocol","Incomplete examination/vital signs","SAE/AESI reporting beyond 24 hours",
    "Unblinded staff performing blinded tasks"),
  pattern=c("^P08[.]01:","^P07[.]10:","^P07[.]06:","^P08[.]05:","^P08: INCOMPLETE ASSESSMENT",
    "^I12[.]01:","^P08: UNBLINDED PERSONNEL|^I06[.]02:"),stringsAsFactors=FALSE)
dv_site_category <- function(term) {
  rules<-dv_site_rules();out<-rep("",length(term));t<-toupper(dv_squish(term))
  for(i in seq_len(nrow(rules))) {
    hit<-grepl(rules$pattern[i],t,perl=TRUE)
    if(any(nzchar(out[hit])))stop("Overlapping site deviation mapping rules.")
    out[hit]<-rules$category[i]
  };out
}
dv_site_window <- function(events,pop,release,window) {
  rules<-dv_site_rules();events$site_category<-dv_site_category(events$term)
  focus<-events[nzchar(events$site_category),,drop=FALSE]
  ids<-sort(unique(pop$SITEID[nzchar(pop$SITEID)]));rows<-list()
  for(site in ids) {
    people<-pop[pop$SITEID==site,,drop=FALSE]
    nv<-sum(people$arm=="mRNA-1273");np<-sum(people$arm=="Placebo")
    z<-focus[focus$SITEID==site,,drop=FALSE]
    for(cat in c("ANY_OF_SEVEN",rules$category)) {
      e<-if(cat=="ANY_OF_SEVEN")z else z[z$site_category==cat,,drop=FALSE]
      a<-length(unique(e$USUBJID[e$arm=="mRNA-1273"]));c<-length(unique(e$USUBJID[e$arm=="Placebo"]))
      r<-dv_compare(a,max(1L,nv),c,max(1L,np),force_fisher=TRUE)
      # A missing arm cannot support a within-site comparison; never infer its denominator.
      if(nv==0L || np==0L) {
        r[1,]<-NA;r$vaccine_n<-a;r$vaccine_N<-nv;r$placebo_n<-c;r$placebo_N<-np
        r$vaccine_percent<-if(nv)100*a/nv else NA_real_
        r$placebo_percent<-if(np)100*c/np else NA_real_
        r$test<-"Not testable (one arm absent)";r$higher_arm<-NA_character_
      } else if(r$test=="No variation")r$test<-"Fisher exact (no variation)"
      label<-if(cat=="ANY_OF_SEVEN")"Any of the seven deviation types" else rules$label[match(cat,rules$category)]
      r<-cbind(data.frame(release=release,window=window,SITEID=site,
        level=if(cat=="ANY_OF_SEVEN")"site_combined" else "site_category",category=cat,deviation=label),r)
      r$vaccine_records<-sum(e$arm=="mRNA-1273");r$placebo_records<-sum(e$arm=="Placebo")
      r$total_records<-nrow(e);r$total_affected<-a+c
      r$vaccine_median_followup<-if(nv)median(people$followup_days[people$arm=="mRNA-1273"]) else NA_real_
      r$placebo_median_followup<-if(np)median(people$followup_days[people$arm=="Placebo"]) else NA_real_
      r$vaccine_childbearing_N<-sum(people$arm=="mRNA-1273"&people$CHILDBP=="Y")
      r$placebo_childbearing_N<-sum(people$arm=="Placebo"&people$CHILDBP=="Y")
      rows[[length(rows)+1L]]<-r
    }
  }
  list(results=if(length(rows))do.call(rbind,rows) else NULL,coverage=data.frame(release=release,window=window,
    sites=length(ids),population_missing_site=sum(!nzchar(pop$SITEID)),
    selected_records_missing_site=sum(!nzchar(focus$SITEID)),selected_records=nrow(focus)),
    mapping=unique(data.frame(release=rep(release,nrow(focus)),term=focus$term,category=focus$site_category)))
}
dv_site_finish <- function(runs,out,status,alpha,min_records,excel_fr,start_basis) {
  dv_csv(do.call(rbind,lapply(runs,function(z)z$coverage)),file.path(out,"site_fisher_coverage.csv"))
  x<-do.call(rbind,lapply(runs,function(z)z$results))
  if(is.null(x)) {
    note<-"Site analysis unavailable: all eligible participants lack SITEID. Overall comparisons remain available."
    writeLines(note,file.path(out,"site_fisher_methods.txt"));warning(note,call.=FALSE)
    return(list(results=data.frame(),selected=data.frame()))
  }
  x$p_holm<-NA_real_
  for(family in unique(x$level)) {
    j<-which(x$level==family & !is.na(x$p_raw))
    x$p_holm[j]<-p.adjust(x$p_raw[j],method="holm")
  }
  x$volume_qualified<-x$total_records>=min_records
  x$signal<-ifelse(is.na(x$p_raw),"NOT_TESTABLE",ifelse(x$p_holm<alpha,"HOLM_ADJUSTED",
    ifelse(x$p_raw<alpha,"NOMINAL_ONLY","NOT_SELECTED")))
  x$selected<-x$volume_qualified & !is.na(x$p_raw) & x$p_raw<alpha
  x<-x[order(x$level,x$window,x$p_raw,x$SITEID),];rownames(x)<-NULL
  selected<-x[x$selected,,drop=FALSE]
  dv_csv(x,file.path(out,"site_fisher_all_tests.csv"))
  dv_csv(x[x$volume_qualified,,drop=FALSE],file.path(out,"site_fisher_volume_qualified.csv"))
  dv_csv(selected,file.path(out,"site_fisher_imbalanced.csv"))
  if(excel_fr)dv_csv(selected,file.path(out,"site_fisher_imbalanced_Excel_FR.csv"),sep=";")
  dv_csv(unique(do.call(rbind,lapply(runs,function(z)z$mapping))),file.path(out,"site_deviation_mapping.csv"))
  dv_csv(dv_site_rules(),file.path(out,"site_deviation_rules.csv"))
  fam<-table(factor(x$level[!is.na(x$p_raw)],levels=c("site_combined","site_category")))
  notes<-c(
    paste("Release coverage:",paste(paste(status$Release,status$Status),collapse="; ")),
    paste0("Selection: at least ",min_records," deviation records across both arms within the same site, window and comparison, and two-sided Fisher exact p < ",alpha,". All volume-qualified rows and all tests are exported, including nonsignificant results."),
    "The combined table pools the seven specified deviation types; each person counts once if affected by any of them. The type-specific table applies the same volume threshold separately to each type. No global-significance prefilter is applied.",
    paste0("Within each site, the 2x2 table contains affected and unaffected participants by original randomized arm. Denominators include all eligible randomized participants at that site, starting at ",if(start_basis=="randomization")"randomization" else "first dose",". Repeated records are used only for the volume filter and the Records V/P column."),
    "Windows and individual follow-up caps match the overall analysis. Prod_EUA is production data retrospectively dated through the configured EUA cutoff; Prod_BLINDED ends at individual blinded-period endpoints or the administrative cutoff. Missing blinded endpoints use available administrative/study-end/death caps, as flagged in population exports.",
    paste0("All site comparisons use two-sided Fisher exact tests, irrespective of volume. Holm correction is separate for the combined family (",fam[["site_combined"]]," site/window comparisons) and the type-specific family (",fam[["site_category"]]," site/type/window comparisons), across analyzed releases and before the volume/significance filter. Overlapping windows are not independent replications."),
    "The report includes nominal signals even if Holm p is not significant. The two p columns use strict display thresholds through < 0.000001; CSVs retain numerical p-values. Zero-variation tables receive p=1; one-arm sites are not testable.",
    "Difference is vaccine minus placebo in percentage points, with pointwise Newcombe/Wilson 95% confidence intervals. These intervals are not adjusted for multiplicity and need not agree with an exact-test decision.",
    "These are participant-level burdens, not failure rates per required procedure. Pregnancy testing uses the whole-site arm denominator for consistency; testing-eligible participants and required testing occasions would be needed for a procedure-specific compliance analysis. Follow-up, illness burden and procedure opportunities may differ.",
    "Site attribution uses ADSL SITEID. Participants/records without a site are excluded only from site comparisons and counted in site_fisher_coverage.csv. A site imbalance does not by itself establish differential conduct or failed blinding.")
  # Reuse the print renderer; site-mode selection also requires volume_qualified.
  shown<-x
  shown$window<-paste(shown$window,ifelse(shown$level=="site_combined","All seven combined","By deviation type"),sep=" | ")
  shown$deviation<-paste(shown$SITEID,shown$deviation,sep=" | ");shown$level<-"term"
  dv_html(shown,data.frame(),status,notes,out,alpha,site_mode=TRUE)
  dv_pdf(shown,data.frame(),notes,out,alpha,site_mode=TRUE)
  writeLines(notes,file.path(out,"site_fisher_methods.txt"))
  message("Site analysis: ",nrow(selected)," nominal rows meet the ",min_records,"-record threshold; ",
    sum(selected$p_holm<alpha)," survive Holm correction.")
  list(results=x,selected=selected)
}

# ---- WEEKLY STAFF-RELATED DEVIATIONS: ALL SITES, NO SIGNIFICANCE FILTER ----
dv_monday <- function(x) {
  x<-as.Date(x);x-((as.integer(x)+3L)%%7L) # 1970-01-01 was a Thursday; locale independent
}
dv_unblind_breakdown <- function(records,pop,release,window,marker) {
  marker<-as.Date(marker);rules<-dv_site_rules();rows<-list()
  ev<-records[records$window_status=="INCLUDED",,drop=FALSE]
  ev$category<-dv_site_category(ev$term);ev<-ev[nzchar(ev$category),,drop=FALSE]
  if(!"UBLDDT"%in%names(pop))pop$UBLDDT<-as.Date(NA)
  valid<-!is.na(pop$UBLDDT) & pop$UBLDDT>=pop$start
  j<-match(ev$USUBJID,pop$USUBJID)
  phase<-rep("UNRESOLVED",nrow(ev))
  known<-which(valid[j]);phase[known]<-ifelse(ev$deviation_date[known]<pop$UBLDDT[j[known]],"BEFORE","ON_AFTER")
  for(basis in c("CALENDAR_MARKER","PARTICIPANT_UBLDDT")) {
    event_phase<-if(basis=="CALENDAR_MARKER")ifelse(ev$deviation_date<marker,"BEFORE","ON_AFTER") else phase
    for(period in if(basis=="CALENDAR_MARKER")c("BEFORE","ON_AFTER") else c("BEFORE","ON_AFTER","UNRESOLVED")) {
      start<-pop$start;end<-pop$window_end;keep<-rep(TRUE,nrow(pop))
      if(basis=="CALENDAR_MARKER") {
        if(period=="BEFORE")end<-pmin(end,marker-1L) else start<-pmax(start,marker)
      } else if(period=="UNRESOLVED")keep<-!valid else {
        keep<-valid
        if(period=="BEFORE")end<-pmin(end,pop$UBLDDT-1L) else start<-pmax(start,pop$UBLDDT)
      }
      days<-pmax(0L,as.integer(end-start)+1L);days[!keep | is.na(days)]<-0L
      for(arm in c("mRNA-1273","Placebo"))for(cat in c("ANY_OF_SEVEN",rules$category)) {
        z<-ev[ev$arm==arm & event_phase==period & (cat=="ANY_OF_SEVEN" | ev$category==cat),,drop=FALSE]
        people<-which(pop$arm==arm & days>0L);pd<-sum(days[people]);affected<-length(unique(z$USUBJID))
        label<-if(cat=="ANY_OF_SEVEN")"All seven deviation types combined" else rules$label[match(cat,rules$category)]
        rows[[length(rows)+1L]]<-data.frame(release=release,window=window,basis=basis,marker_date=marker,
          period=period,category=cat,label=label,arm=arm,records=nrow(z),affected_participants=affected,
          participants_observed=length(people),percent_affected=if(length(people))100*affected/length(people) else NA_real_,
          person_days=pd,person_weeks=pd/7,records_per_1000_person_weeks=if(pd>0)nrow(z)*7000/pd else NA_real_,
          records_on_unblinding_date=sum(z$deviation_date==pop$UBLDDT[match(z$USUBJID,pop$USUBJID)],na.rm=TRUE))
        stopifnot(affected<=length(people),pd>0 || nrow(z)==0L)
      }
    }
  }
  do.call(rbind,rows)
}
dv_weekly_full_followup <- function(events,context,release,cap,out,unblinding_start) {
  # Weekly/descriptive output only. Existing blinded-window tests are unchanged.
  ctx<-context;ctx$end<-pmin(ctx$CUTOFFDT,as.Date(cap));ctx$end_basis<-"Administrative cutoff"
  for(k in c("EOSDT","DTHDT")) {
    take<-!is.na(ctx[[k]]) & ctx[[k]]<ctx$end
    ctx$end[take]<-ctx[[k]][take];ctx$end_basis[take]<-k
  }
  ctx$eligible<-ctx$RANDFL=="Y" & !is.na(ctx$start) & ctx$end>=ctx$start
  ctx$eligible[is.na(ctx$eligible)]<-FALSE
  ctx$exclusion<-ifelse(ctx$RANDFL!="Y","NOT_RANDOMIZED",ifelse(is.na(ctx$start),"START_DATE_MISSING",
    ifelse(ctx$end<ctx$start,"END_BEFORE_START","")))
  pop<-ctx[ctx$eligible,,drop=FALSE];pop$window_end<-pop$end
  pop$followup_days<-as.integer(pop$window_end-pop$start)+1L;pop$window<-paste0(release,"_ALL_FOLLOWUP")
  j<-match(events$USUBJID,pop$USUBJID);reason<-rep("INCLUDED",nrow(events))
  reason[is.na(j)]<-"OUTSIDE_ELIGIBLE_COHORT"
  reason[!is.na(j)&is.na(events$deviation_date)]<-"DATE_UNRESOLVED"
  reason[which(!is.na(j)&!is.na(events$deviation_date)&events$deviation_date<pop$start[j])]<-"BEFORE_START"
  reason[which(!is.na(j)&!is.na(events$deviation_date)&events$deviation_date>pop$window_end[j])]<-"AFTER_WINDOW"
  records<-events;records$window<-paste0(release,"_ALL_FOLLOWUP");records$release<-release;records$window_status<-reason
  records$window_start<-pop$start[j];records$window_end<-pop$window_end[j]
  records$subject_unblinding_date<-pop$UBLDDT[j]
  records$calendar_unblinding_period<-ifelse(is.na(records$deviation_date),"UNRESOLVED_DATE",
    ifelse(records$deviation_date<as.Date(unblinding_start),"BEFORE","ON_AFTER"))
  known<-!is.na(j)&!is.na(records$subject_unblinding_date)&records$subject_unblinding_date>=pop$start[j]
  records$participant_unblinding_period<-"UNRESOLVED"
  use<-which(known & !is.na(records$deviation_date))
  records$participant_unblinding_period[use]<-ifelse(records$deviation_date[use]<records$subject_unblinding_date[use],"BEFORE","ON_AFTER")
  records$participant_unblinding_period[is.na(records$deviation_date)]<-"UNRESOLVED_DATE"
  dv_csv(records,file.path(out,paste0("records_",release,"_ALL_FOLLOWUP.csv")))
  dv_csv(pop,file.path(out,paste0("population_",release,"_ALL_FOLLOWUP.csv")))
  dv_weekly_window(records,pop,release,paste0(release,"_ALL_FOLLOWUP"),unblinding_start)
}
dv_unblind_inventory <- function(context,release,marker,out) {
  x<-context[context$RANDFL=="Y",c("USUBJID","SITEID","arm","RANDDT","start","UBLDDT","AP01EDT","AP02SDT","EOSDT","DTHDT","CUTOFFDT")]
  x$release<-release;x$before_calendar_reference<-!is.na(x$UBLDDT)&x$UBLDDT<as.Date(marker)
  x$date_status<-ifelse(is.na(x$UBLDDT),"MISSING",ifelse(x$UBLDDT<x$RANDDT,"BEFORE_RANDOMIZATION",
    ifelse(x$UBLDDT>x$CUTOFFDT,"AFTER_DATA_CUTOFF","RECORDED")))
  dv_csv(x,file.path(out,paste0("unblinding_subject_dates_",release,".csv")))
  daily<-x[!is.na(x$UBLDDT),,drop=FALSE]
  if(nrow(daily)) {
    d<-aggregate(rep(1L,nrow(daily)),daily[c("release","UBLDDT","arm")],sum)
    names(d)[4]<-"participants";d<-d[order(d$UBLDDT,d$arm),]
  } else d<-data.frame(release=character(),UBLDDT=as.Date(character()),arm=character(),participants=integer())
  dv_csv(d,file.path(out,paste0("unblinding_daily_",release,".csv")))
  data.frame(release=release,calendar_reference=as.Date(marker),randomized=nrow(x),
    first_recorded_unblinding=if(nrow(daily))min(daily$UBLDDT) else as.Date(NA),
    earlier_unblindings=sum(x$before_calendar_reference),missing_UBLDDT=sum(is.na(x$UBLDDT)),
    before_randomization=sum(x$date_status=="BEFORE_RANDOMIZATION"),after_cutoff=sum(x$date_status=="AFTER_DATA_CUTOFF"))
}
dv_weekly_window <- function(records,pop,release,window,unblinding_start="2020-12-14") {
  rules<-dv_site_rules();arms<-c("mRNA-1273","Placebo")
  records$weekly_category<-dv_site_category(records$term)
  selected<-records[nzchar(records$weekly_category),,drop=FALSE]
  ev<-selected[selected$window_status=="INCLUDED",,drop=FALSE]
  ev$week_start<-dv_monday(ev$deviation_date)
  first<-min(pop$start);last<-max(pop$window_end)
  weeks<-seq(dv_monday(first),dv_monday(last),by="7 days")
  rows<-list();coverage<-list()
  for(arm in arms) {
    p<-pop[pop$arm==arm,,drop=FALSE]
    person_days<-vapply(seq_along(weeks),function(i)
      sum(pmax(0L,as.integer(pmin(p$window_end,weeks[i]+6L)-pmax(p$start,weeks[i]))+1L)),numeric(1))
    n_observed<-vapply(seq_along(weeks),function(i)
      sum(p$start<=weeks[i]+6L & p$window_end>=weeks[i]),integer(1))
    for(cat in c("ANY_OF_SEVEN",rules$category)) {
      z<-ev[ev$arm==arm & (cat=="ANY_OF_SEVEN" | ev$weekly_category==cat),,drop=FALSE]
      src<-selected[which(selected$arm==arm & (cat=="ANY_OF_SEVEN" | selected$weekly_category==cat)),,drop=FALSE]
      label<-if(cat=="ANY_OF_SEVEN")"All seven deviation types combined" else rules$label[match(cat,rules$category)]
      counts<-tabulate(match(z$week_start,weeks),nbins=length(weeks))
      u<-unique(z[c("week_start","USUBJID")])
      affected<-tabulate(match(u$week_start,weeks),nbins=length(weeks))
      s<-unique(z[nzchar(z$SITEID),c("week_start","SITEID")])
      n_sites<-tabulate(match(s$week_start,weeks),nbins=length(weeks))
      rows[[length(rows)+1L]]<-data.frame(release=release,window=window,category=cat,label=label,
        week_start=weeks,week_end=weeks+6L,observed_start=pmax(weeks,first),observed_end=pmin(weeks+6L,last),
        partial_boundary_week=weeks<first | weeks+6L>last,arm=arm,records=counts,
        affected_participants=affected,participants_observed=n_observed,person_days=person_days,
        person_weeks=person_days/7,records_per_1000_person_weeks=ifelse(person_days>0,counts*7000/person_days,NA_real_),
        affected_sites=n_sites)
      coverage[[length(coverage)+1L]]<-data.frame(release=release,window=window,category=cat,label=label,arm=arm,
        N=nrow(p),sites=length(unique(p$SITEID[nzchar(p$SITEID)])),window_start=first,window_end=last,
        records=nrow(z),affected_participants=length(unique(z$USUBJID)),person_days=sum(person_days),
        records_missing_site=sum(!nzchar(z$SITEID)),participants_missing_site=sum(!nzchar(p$SITEID)),
        missing_blind_endpoint=sum(p$period_endpoint_missing),
        undated_records=sum(src$window_status=="DATE_UNRESOLVED"),
        before_start_records=sum(src$window_status=="BEFORE_START"),after_window_records=sum(src$window_status=="AFTER_WINDOW"))
      stopifnot(sum(counts)==nrow(z),all(affected<=counts),all(affected<=n_observed),
        sum(person_days)==sum(p$followup_days),all(counts[person_days==0]==0))
    }
  }
  list(data=do.call(rbind,rows),coverage=do.call(rbind,coverage),
    breakdown=dv_unblind_breakdown(records,pop,release,window,unblinding_start))
}
dv_weekly_page <- function(x,coverage,window,categories,eua_cutoff,page,total,font="sans",unblinding_start="2020-12-14") {
  # Explicit millimetre positions keep vector-PDF and 300-dpi PNG layouts identical.
  grid::grid.newpage()
  txt<-function(label,x,y,size=9,bold=FALSE,color="#182631",just="left")
    grid::grid.text(label,grid::unit(x,"mm"),grid::unit(y,"mm"),just=just,
      gp=grid::gpar(fontfamily=font,fontsize=size,fontface=if(bold)"bold" else "plain",col=color))
  line<-function(x1,y1,x2,y2,color="#ccd3d9",width=.5,lty=1)
    grid::grid.lines(grid::unit(c(x1,x2),"mm"),grid::unit(c(y1,y2),"mm"),gp=grid::gpar(col=color,lwd=width,lty=lty))
  rect<-function(x,y,w,h,fill,border=NA)
    grid::grid.rect(grid::unit(x,"mm"),grid::unit(y,"mm"),grid::unit(w,"mm"),grid::unit(h,"mm"),
      just=c("left","bottom"),gp=grid::gpar(fill=fill,col=border,lwd=.4))
  arms<-c("mRNA-1273","Placebo");colors<-c("#165A82","#E6A143");borders<-c("#153F5A","#85561F")
  z<-x[x$window==window,,drop=FALSE];meta<-coverage[coverage$window==window & coverage$category=="ANY_OF_SEVEN",]
  meta<-meta[match(arms,meta$arm),];weeks<-sort(unique(z$week_start));nweeks<-length(weeks)
  first<-min(meta$window_start);last<-max(meta$window_end);single<-length(categories)==1L
  title<-switch(window,Prod_EUA="Production data | through the EUA cutoff",
    Prod_BLINDED="Production data | through individual blinded-period ends",
    Prod_ALL_FOLLOWUP="Production data | full available follow-up, including after unblinding",
    EUA_SNAPSHOT="EUA release | through the EUA cutoff",window)
  txt("MODERNA P301  /  PROTOCOL DEVIATIONS",12,200,8,TRUE,color="#52677a")
  txt("Weekly staff-related deviations",12,190,18,TRUE)
  txt(title,12,182,10,TRUE)
  txt(paste0(first," to ",last," | All sites pooled | Original randomized arm | ",
    if(single)"Seven specified types combined" else "Breakdown by recorded deviation type"),12,175,8.5)
  for(a in seq_along(arms)) {
    lx<-12+(a-1)*81;rect(lx,164.4,4,3.3,colors[a],borders[a])
    txt(paste0(arms[a]," (N = ",format(meta$N[a],big.mark=",",trim=TRUE),")"),lx+6,166,9)
  }
  txt("Y: records per week",206,166,8.5,TRUE)
  line(12,160,285,160)
  panel<-function(cat,x0,y0,w,h) {
    d<-z[z$category==cat,,drop=FALSE];v<-d[d$arm==arms[1],];p<-d[d$arm==arms[2],]
    v<-v[match(weeks,v$week_start),];p<-p[match(weeks,p$week_start),]
    counts<-cbind(v$records,p$records)
    cm<-coverage[coverage$window==window & coverage$category==cat,,drop=FALSE];cm<-cm[match(arms,cm$arm),]
    label<-unique(d$label)
    if(!single)label<-paste(strwrap(label,width=46),collapse="\n")
    txt(label,x0,y0+h,if(single)12 else 10,TRUE,just=c("left","top"))
    txt(paste0("Records: ",cm$records[1]," / ",cm$records[2],"   |   Participants: ",
      cm$affected_participants[1]," / ",cm$affected_participants[2],"   (vaccine / placebo)"),
      x0,y0+h-if(single)8 else 10,if(single)9 else 7.4)
    left<-x0+13;right<-x0+w-1;bottom<-y0+13;top<-y0+h-if(single)18 else 20
    ymax<-max(1,counts);ticks<-pretty(c(0,ymax),n=if(single)5 else 3)
    ticks<-sort(unique(c(0,ceiling(ticks[ticks>=0]))));ymax<-max(ticks)
    sy<-function(y)bottom+y/ymax*(top-bottom);step<-(right-left)/nweeks
    for(t in ticks) {line(left,sy(t),right,sy(t),if(t==0)"#778794" else "#e0e5e9");txt(format(t,trim=TRUE),left-2,sy(t),8,just="right")}
    for(i in seq_len(nweeks))for(a in seq_along(arms))if(counts[i,a]>0)
      rect(left+(i-1)*step+if(a==1).10*step else .52*step,bottom,.36*step,sy(counts[i,a])-bottom,colors[a],borders[a])
    # Dashed line denotes the exact cutoff date within its week, not a whole-week boundary.
    cutoff<-as.Date(eua_cutoff)
    if(first<cutoff && last>cutoff) {
      cx<-left+(as.numeric(cutoff-weeks[1])+.5)/7*step
      line(cx,bottom,cx,top,"#454f59",.7,2)
      txt("EUA cutoff",cx-1,top+1,7,color="#454f59",just=c("right","bottom"))
    }
    unblind<-as.Date(unblinding_start)
    if(first<=unblind && last>=unblind) {
      ux<-left+as.numeric(unblind-weeks[1])/7*step
      line(ux,bottom,ux,top,"#C82032",1.2,2)
      lab<-if(single)paste0("EUA-era unblinding: ",unblind) else "Unblinding reference"
      align_left<-ux<left+.78*(right-left)
      txt(lab,if(align_left)ux+1 else ux-1,top+1,7.2,TRUE,"#C82032",
        just=c(if(align_left)"left" else "right","bottom"))
    }
    ix<-unique(round(seq(1,nweeks,length.out=min(nweeks,if(single)9L else 5L))))
    for(i in ix) {
      cx<-left+(i-.5)*step;line(cx,bottom,cx,bottom-1,"#778794")
      txt(format(weeks[i],"%d/%m"),cx,bottom-4,8,just="center")
    }
    txt("Week starting Monday (DD/MM)",(left+right)/2,y0+1,8,just="center")
    if(!any(counts))txt("No recorded deviations in this window",(left+right)/2,(bottom+top)/2,9,just="center",color="#52677a")
  }
  if(single)panel(categories,12,45,273,108) else {
    for(i in seq_along(categories))panel(categories[i],12+((i-1)%%2)*141,100-floor((i-1)/2)*66,132,55)
    txt("Each panel has its own y-axis scale. Compare arm heights within a panel.",12,28,8,color="#52677a")
  }
  note1<-paste0("Monday-Sunday bins by DVSTDTC. Red line: EUA-era unblinding reference, ",unblinding_start,
    if(as.Date(unblinding_start)>last)" (after this window)." else if(as.Date(unblinding_start)<first)" (before this window)." else ".")
  note2<-"All sites and all seven specified types are retained, irrespective of volume or significance. Repeated records are counted."
  note3<-if(window=="Prod_ALL_FOLLOWUP")
    "Arms retain original assignment after crossover; original placebo participants may have received vaccine. Counts are descriptive." else
    "Individual blinded-period caps still apply: after the red line does not mean after that participant's own unblinding."
  if(single) {
    txt("The seven-type total is a descriptive sum; it combines deviations with different clinical implications.",12,34,8,color="#52677a")
    txt("Weekly participant counts and person-days are exported in weekly_staff_deviations.csv.",12,29,8,color="#52677a")
  }
  for(i in seq_along(c(note1,note2,note3)))txt(c(note1,note2,note3)[i],12,22-(i-1)*4,7.4,color="#52677a")
  if(sum(meta$missing_blind_endpoint)>0 && window!="Prod_ALL_FOLLOWUP")txt(paste0("Blind endpoint missing for ",sum(meta$missing_blind_endpoint),
    " participants: available administrative/study-end/death caps used."),12,10,7,color="#52677a")
  txt(paste0("analyse_deviations.R v1.3.0 | ",window," | ",page," / ",total),285,5,7,color="#52677a",just="right")
}
dv_unblind_page <- function(breakdown,window,basis,page,total,font="sans") {
  grid::grid.newpage()
  txt<-function(label,x,y,size=9,bold=FALSE,color="#182631",just="left")
    grid::grid.text(label,grid::unit(x,"mm"),grid::unit(y,"mm"),just=just,
      gp=grid::gpar(fontfamily=font,fontsize=size,fontface=if(bold)"bold" else "plain",col=color))
  line<-function(y)grid::grid.lines(grid::unit(c(12,285),"mm"),grid::unit(c(y,y),"mm"),gp=grid::gpar(col="#dce2e7",lwd=.5))
  d<-breakdown[breakdown$window==window & breakdown$basis==basis,,drop=FALSE]
  calendar<-basis=="CALENDAR_MARKER";marker<-unique(d$marker_date)
  txt("MODERNA P301  /  PROTOCOL DEVIATIONS",12,200,8,TRUE,color="#52677a")
  txt(if(calendar)"Before / after the calendar marker" else "Before / after each participant's unblinding",12,190,17,TRUE)
  txt(paste0(window," | All sites | Original randomized arms | ",
    if(calendar)paste0("Reference: ",marker) else "Split at each recorded UBLDDT"),12,181,9)
  colors<-c("#165A82","#85561F");xs<-c(108,139,173,203,239,272)
  fmt<-function(n)ifelse(is.na(n),"-",format(n,big.mark=",",trim=TRUE,scientific=FALSE))
  rate<-function(n)if(is.na(n))"-" else if(n>0 && n<.01)"<0.01" else sprintf("%.2f",n)
  for(g in seq_along(c("BEFORE","ON_AFTER"))) {
    period<-c("BEFORE","ON_AFTER")[g];y<-171-(g-1)*69
    q<-d[d$period==period,,drop=FALSE];m<-q[q$category=="ANY_OF_SEVEN",];m<-m[match(c("mRNA-1273","Placebo"),m$arm),]
    heading<-if(calendar)paste(if(g==1)"Before" else "On/after",marker) else if(g==1)"Before recorded individual unblinding" else "On/after recorded individual unblinding"
    txt(heading,12,y,10,TRUE,color=if(calendar)"#C82032" else "#182631")
    txt(paste0("Participants observed: ",fmt(m$participants_observed[1])," / ",fmt(m$participants_observed[2]),
      " | Person-weeks: ",fmt(round(m$person_weeks[1],1))," / ",fmt(round(m$person_weeks[2],1))," (vaccine / placebo)"),12,y-5.5,8)
    txt("Deviation category",12,y-12,8,TRUE)
    labels<-c("Records V","Records P","People V","People P","Rate V*","Rate P*")
    for(i in seq_along(xs))txt(labels[i],xs[i],y-12,8,TRUE,colors[(i-1)%%2+1],just="right")
    line(y-15)
    cats<-c("ANY_OF_SEVEN",dv_site_rules()$category)
    for(i in seq_along(cats)) {
      v<-q[q$category==cats[i] & q$arm=="mRNA-1273",];p<-q[q$category==cats[i] & q$arm=="Placebo",]
      yy<-y-20-(i-1)*5.6
      txt(v$label,12,yy,8.2,i==1)
      cells<-c(if(v$participants_observed)fmt(v$records) else "-",if(p$participants_observed)fmt(p$records) else "-",
        if(v$participants_observed)fmt(v$affected_participants) else "-",if(p$participants_observed)fmt(p$affected_participants) else "-",
        rate(v$records_per_1000_person_weeks),rate(p$records_per_1000_person_weeks))
      for(k in seq_along(xs))txt(cells[k],xs[k],yy,8.3,i==1,just="right")
      line(yy-2.7)
    }
  }
  unknown<-d[d$period=="UNRESOLVED" & d$category=="ANY_OF_SEVEN",]
  if(nrow(unknown)) {
    unknown<-unknown[match(c("mRNA-1273","Placebo"),unknown$arm),]
    txt(paste0("Unresolved UBLDDT, reported separately: ",fmt(unknown$records[1])," / ",fmt(unknown$records[2]),
      " records; ",fmt(unknown$affected_participants[1])," / ",fmt(unknown$affected_participants[2]),
      " affected participants. Details by category are in the CSV."),12,31,8,color="#52677a")
  } else txt("Calendar periods do not identify each participant's actual blinding status.",12,31,8,color="#52677a")
  txt("* Records per 1,000 observed person-weeks, using the whole cohort; this is not a failure rate per required procedure.",12,23,8,color="#52677a")
  txt("People are deduplicated within each period/category and can appear in both periods. The combined row deduplicates across categories.",12,18,7.6,color="#52677a")
  txt(if(window=="Prod_ALL_FOLLOWUP")"Follow-up includes crossover: original placebo recipients may have received vaccine. Same-day events are assigned to on/after; ordering is unknown." else
    "Individual blinded-period limits remain in force. After the calendar marker, this table retains only follow-up allowed by those limits.",12,13,7.4,color="#52677a")
  txt(paste0("analyse_deviations.R v1.3.0 | ",basis," | ",page," / ",total),285,5,7,color="#52677a",just="right")
}
dv_weekly_finish <- function(runs,out,status,eua_cutoff,png_dpi,excel_fr,start_basis,unblinding_start="2020-12-14",unblinding_inventory=NULL) {
  x<-do.call(rbind,lapply(runs,function(z)z$data));coverage<-do.call(rbind,lapply(runs,function(z)z$coverage))
  breakdown<-do.call(rbind,lapply(runs,function(z)z$breakdown))
  x<-x[order(x$window,x$category,x$week_start,x$arm),];rownames(x)<-NULL
  dv_csv(x,file.path(out,"weekly_staff_deviations.csv"))
  dv_csv(coverage,file.path(out,"weekly_staff_deviations_coverage.csv"))
  dv_csv(breakdown,file.path(out,"deviations_before_after_unblinding.csv"))
  if(!is.null(unblinding_inventory))dv_csv(unblinding_inventory,file.path(out,"unblinding_reference_audit.csv"))
  if(excel_fr){dv_csv(x,file.path(out,"weekly_staff_deviations_Excel_FR.csv"),sep=";")
    dv_csv(breakdown,file.path(out,"deviations_before_after_unblinding_Excel_FR.csv"),sep=";")}
  notes<-c(
    paste("Release coverage:",paste(paste(status$Release,status$Status),collapse="; ")),
    "Descriptive weekly grouped column charts, mRNA-1273 versus placebo, by original randomized assignment. All sites are pooled, including included records whose SITEID is missing. No site-volume, significance or sponsor-significance filter is applied.",
    "The seven categories use exactly dv_site_rules(), exported as site_deviation_rules.csv: short post-dose observation; missing 30-minute diary; required nasopharyngeal swab missing; pregnancy testing not per protocol; incomplete examination/vital signs; delayed SAE/AESI reporting; unblinded staff performing blinded tasks. Unmapped terms are outside this predefined chart scope, not presumed unrelated to staff.",
    "Date is the resolved DVSTDTC/deviation_date, never the data-entry, awareness or review date. Weeks run Monday-Sunday, labelled by Monday. Undated/partial-date records are not assigned to a week or imputed; eligible-cohort unresolved records are counted separately by category and arm in the coverage export.",
    paste0("Cohort starts at ",if(start_basis=="randomization")"randomization" else "first dose",". Prior EUA/blinded charts retain their original individual caps. The new Prod_ALL_FOLLOWUP descriptive window ends at earliest of administrative cutoff, study end or death; it deliberately continues through unblinding and crossover. No new global/site hypothesis tests are added for this window. Last/first weeks may be partial. Zero-record weeks are explicitly included. A partial_boundary_week flag describes the overall date boundary, not individual follow-up completeness."),
    "Primary y-axis is deviation RECORDS, so a person may contribute several times. The combined series sums the seven mutually exclusive record categories. affected_participants counts distinct subjects per week/category/arm, deduplicated again across categories for ANY_OF_SEVEN. Summing weekly affected counts does not give the number of unique participants over the window; use the coverage export for that total.",
    "Weekly person_days sums days actually within each participant's analysis window, inclusively. participants_observed is the number with at least one day in that week; person_weeks=person_days/7. records_per_1000_person_weeks is an optional descriptive exposure normalization, not a per-procedure failure rate. Denominators include the full analysis cohort, not only those eligible for pregnancy testing or another particular procedure.",
    "Separate releases/windows are plotted separately, never summed. Prod_EUA uses production records dated through the EUA cutoff; it is not evidence of what was present in the original EUA submission. Overlapping windows are not independent observations. Prod_BLINDED uses individual blind-period caps where available; missing endpoints are flagged and use the documented fallback caps.",
    paste0("The dashed red calendar marker is ",unblinding_start," (PHMPT_UNBLINDING_START). The default 2020-12-14 denotes the start of the EUA-era unblinding reference period, not the earliest individual unblinding or the scheduled open-label program. Earlier individual UBLDDT dates are preserved in unblinding_subject_dates_*.csv, daily counts in unblinding_daily_*.csv, and summarized in unblinding_reference_audit.csv. Dates outside a plot's window are stated in its footer; the axis is not artificially extended."),
    "Milestone sources: Mena Lora et al., JAMA Network Open 2023, doi:10.1001/jamanetworkopen.2022.51974, describes EUA-related unblinding eligibility after December 14, 2020 (https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2800703). The COVE protocol amendment was December 23; scheduled participant decision/unblinding visits started December 24, as described by Tsiatis and Davidian, Biometrics, doi:10.1111/biom.13509 (https://pmc.ncbi.nlm.nih.gov/articles/PMC8444907/). These dates describe different mechanisms and must not be conflated.",
    "deviations_before_after_unblinding.csv contains two separate breakdowns. CALENDAR_MARKER splits all included dates/person-days strictly before versus on/after the marker. PARTICIPANT_UBLDDT splits at each participant's own recorded date; missing dates or dates before analysis start form an UNRESOLVED group. AP02SDT is not substituted for missing UBLDDT. Dates beyond the analysis cutoff remain recorded but yield no post-unblinding time within the window. Split counts use exact dates, not assignment of a whole calendar week to one side.",
    "Individual unblinding-day records and the entire date are assigned to ON_AFTER; time ordering within that day is unknown, and records_on_unblinding_date is exported. Period-specific denominators count participants with at least one observed day in that period. A participant may appear in both periods, and counts across types/periods should not be added to infer distinct people. A zero-exposure period has missing rates/percentages and is printed as a dash, not a measured zero incidence.",
    "The original randomized-arm labels are retained after crossover; the original placebo group may include vaccinated participants. Calendar changes, crossover, participant selection and procedure opportunities limit causal interpretation. The calendar comparison on Prod_BLINDED is still restricted to retained blinded follow-up after the marker; the full-follow-up view is needed to see individual post-unblinding records.",
    "One combined chart and two four-panel pages per available window, on A4 landscape vector PDF; PNGs reproduce each page at the configured DPI. Panels have separate y scales to keep low-volume categories visible, with a common calendar axis within each window. Color and luminance differ between arms, with the vaccine bar always left of placebo. All y axes start at zero.",
    "Raw temporal differences can reflect enrollment, follow-up, illness burden, required procedure opportunities, documentation or staff conduct. These charts cannot identify the mechanism or establish failed blinding. No per-week tests are performed.")
  writeLines(notes,file.path(out,"weekly_staff_deviations_methods.txt"))
  cats<-dv_site_rules()$category;pages<-list()
  for(w in unique(coverage$window)) {
    for(group in list("ANY_OF_SEVEN",cats[1:4],cats[5:7]))pages[[length(pages)+1L]]<-list(window=w,categories=group,kind="chart")
    m<-coverage[coverage$window==w,,drop=FALSE]
    if(min(m$window_start)<as.Date(unblinding_start) && max(m$window_end)>=as.Date(unblinding_start))
      pages[[length(pages)+1L]]<-list(window=w,kind="CALENDAR_MARKER")
    if(w=="Prod_ALL_FOLLOWUP")pages[[length(pages)+1L]]<-list(window=w,kind="PARTICIPANT_UBLDDT")
  }
  cairo<-isTRUE(capabilities("cairo"));font<-if(cairo)"sans" else "Helvetica"
  draw<-function(i) {
    pg<-pages[[i]]
    if(pg$kind=="chart")dv_weekly_page(x,coverage,pg$window,pg$categories,eua_cutoff,i,length(pages),font,unblinding_start)
    else dv_unblind_page(breakdown,pg$window,pg$kind,i,length(pages),font)
  }
  pdf_file<-file.path(out,"weekly_staff_deviations.pdf")
  if(cairo)grDevices::cairo_pdf(pdf_file,width=297/25.4,height=210/25.4,family=font,onefile=TRUE)
  else grDevices::pdf(pdf_file,width=297/25.4,height=210/25.4,family=font,onefile=TRUE)
  pdf_device<-grDevices::dev.cur()
  tryCatch({for(i in seq_along(pages))draw(i)},
    finally=grDevices::dev.off(pdf_device))
  png_status<-character(length(pages))
  for(i in seq_along(pages)) {
    path<-file.path(out,sprintf("weekly_staff_deviations_%02d_%s.png",i,pages[[i]]$window))
    png_status[i]<-tryCatch({
      if(cairo)grDevices::png(path,width=297,height=210,units="mm",res=png_dpi,type="cairo")
      else grDevices::png(path,width=297,height=210,units="mm",res=png_dpi)
      png_device<-grDevices::dev.cur()
      tryCatch(draw(i),
        finally=grDevices::dev.off(png_device))
      bytes<-readBin(path,"raw",n=file.info(path)$size)
      if(length(bytes)<20L || !identical(bytes[1:8],as.raw(c(137,80,78,71,13,10,26,10))) ||
          !identical(tail(bytes,12L),as.raw(c(0,0,0,0,73,69,78,68,174,66,96,130))))stop("PNG file incomplete after device close.")
      "WRITTEN"
    },error=function(e)paste("FAILED:",conditionMessage(e)))
  }
  dv_csv(data.frame(page=seq_along(pages),window=vapply(pages,function(z)z$window,character(1)),kind=vapply(pages,function(z)z$kind,character(1)),png_status=png_status),
    file.path(out,"weekly_staff_deviations_pages.csv"))
  if(any(png_status!="WRITTEN"))warning("Some weekly PNGs failed; see weekly_staff_deviations_pages.csv. The vector PDF is available.",call.=FALSE)
  message("Weekly charts: ",length(pages)," print pages; all sites pooled, no significance filter.")
  list(data=x,coverage=coverage,breakdown=breakdown,unblinding_inventory=unblinding_inventory,pdf=pdf_file)
}

# ---- PRINTABLE TABLES ----
dv_escape <- function(x) {
  x<-as.character(x);x[is.na(x)]<-"";x<-gsub("&","&amp;",x,fixed=TRUE);x<-gsub("<","&lt;",x,fixed=TRUE)
  x<-gsub(">","&gt;",x,fixed=TRUE);gsub('"',"&quot;",x,fixed=TRUE)
}
dv_p <- function(x) {
  # Display only: retain numeric p-values for testing, sorting and CSV exports.
  # Use the smallest strictly greater threshold, down to five zeros before 1.
  out <- rep("-",length(x))
  valid <- !is.na(x)
  out[valid] <- formatC(x[valid],format="f",digits=4)
  thresholds <- c(.01,.001,.0001,.00001,.000001)
  labels <- c("< 0.01","< 0.001","< 0.0001","< 0.00001","< 0.000001")
  for(i in seq_along(thresholds)) out[valid & x < thresholds[i]] <- labels[i]
  out
}
dv_display <- function(x,site_mode=FALSE) {
  data.frame(Deviation=x$deviation,
    Vaccine=if(site_mode)sprintf("%d/%d (%.2f%%)",x$vaccine_n,x$vaccine_N,x$vaccine_percent) else sprintf("%s (%.2f%%)",format(x$vaccine_n,big.mark=",",trim=TRUE),x$vaccine_percent),
    Placebo=if(site_mode)sprintf("%d/%d (%.2f%%)",x$placebo_n,x$placebo_N,x$placebo_percent) else sprintf("%s (%.2f%%)",format(x$placebo_n,big.mark=",",trim=TRUE),x$placebo_percent),
    Difference=sprintf("%+.2f [%+.2f, %+.2f]",x$difference_pp,x$lower95_pp,x$upper95_pp),
    Test=if(site_mode)paste(x$vaccine_records,x$placebo_records,sep=" / ") else x$test,
    P=dv_p(x$p_raw),Holm=dv_p(x$p_holm),stringsAsFactors=FALSE)
}
dv_notes <- function(start_basis,alpha,expected_min,tests) c(
  paste0("Participant-level comparison: each person is counted once per deviation term and window. Denominators include eligible randomized participants without deviations; original TRT01P defines the arm."),
  paste0("Start: ",if(start_basis=="randomization")"randomization date" else "first recorded blinded dose",". End: earliest of data cutoff, AP01EDT, AP02SDT minus one day, UBLDDT minus one day, study end or death. Date-only records cannot resolve within-day ordering."),
  "If all blinded-period endpoints are missing, the administrative cutoff and available study-end/death dates are used; blinding throughout that interval is not verified. These participants are flagged in population exports and counted by arm in denominators.csv.",
  "Prod_EUA is production data retrospectively dated through the EUA cutoff, not the original EUA submission. Releases and windows are never pooled.",
  paste0("Two-sided Fisher exact test if any expected cell is < ",expected_min,"; otherwise Pearson chi-square (1 df, no Yates correction). All-zero/all-affected outcomes are marked No variation (p=1)."),
  paste0("Holm correction covers all ",tests," term-by-window comparisons across successfully analyzed releases, including rows absent from the printed selection. Sponsor-category tests form a separate exploratory family."),
  paste0("Print selection: raw p < ",alpha," in either direction. Adjusted signals and nominal-only signals are separated. A small p-value does not measure the magnitude or practical importance of a difference."),
  "Difference is vaccine minus placebo in percentage points; 95% confidence intervals use Newcombe's combination of Wilson score intervals, without continuity correction. These intervals are pointwise, not multiplicity-adjusted, and need not reproduce the exact-test p-value.",
  "These are crude participant burdens, not rates per required procedure, dose, safety event or illness visit. Opportunity counts, follow-up duration and site clustering can differ. An imbalance alone does not establish faulty blinding, investigator intent or missing COVID cases.",
  "Full terms retain DVTERM1 and the sponsor wording (including not performed OR out of window). DVSIG is a sponsor classification, not the statistical significance of this comparison. No restriction to significant deviations or per-protocol subjects is made.",
  "Rows overlap across deviation terms and time windows; percentages must not be summed. Records with unresolved dates, pre-start records and records outside the cohort remain in the audit exports. No missing dates are imputed."
)
dv_html <- function(results,den,status,notes,out,alpha,site_mode=FALSE) {
  esc<-dv_escape;sections<-character()
  stem<-if(site_mode)"site_fisher_imbalances" else "deviation_imbalances"
  title<-if(site_mode)"Moderna P301: site-level deviation imbalances" else "Moderna P301: recorded deviation imbalances"
  all_csv<-if(site_mode)"site_fisher_all_tests.csv" else "deviations_all_tests.csv"
  selected_csv<-if(site_mode)"site_fisher_imbalanced.csv" else "deviations_imbalanced.csv"
  groups<-if(site_mode)paste0("Fisher p < ",alpha,"; Holm p shown") else c("Holm-adjusted signals","Nominal-only signals")
  for(window in unique(results$window)) {
    z<-results[results$window==window&results$level=="term",,drop=FALSE]
    n<-den[den$window==window,,drop=FALSE]
    heads<-if(site_mode)c("Site | Deviation","mRNA-1273<br>n/N (%)","Placebo<br>n/N (%)","Difference, pp<br>[95% CI]","Records<br>V / P","Fisher p","Holm p") else c("Recorded deviation","mRNA-1273<br>n (%)","Placebo<br>n (%)","Difference, pp<br>[95% CI]","Test","Raw p","Holm p")
    body<-character()
    for(group in groups) {
      hit<-!is.na(z$p_raw) & z$p_raw<alpha & if(site_mode)TRUE else if(group=="Holm-adjusted signals")z$p_holm<alpha else z$p_holm>=alpha
      if(site_mode)hit<-hit & z$volume_qualified
      d<-z[hit,,drop=FALSE];d<-d[order(d$p_holm,d$p_raw,-abs(d$difference_pp)),,drop=FALSE]
      rows<-if(nrow(d))apply(dv_display(d,site_mode),1,function(r)paste0('<tr>',paste0('<td>',esc(r),'</td>',collapse=''),'</tr>')) else
        '<tr><td colspan="7">No deviations meet this criterion.</td></tr>'
      col_widths<-if(site_mode)c(30,15,15,16,8,8,8) else c(38,11,11,18,8,7,7)
      body<-c(body,paste0('<h3>',esc(group),'</h3><table class="signals"><colgroup>',paste0('<col style="width:',col_widths,'%">',collapse=''),'</colgroup><thead><tr>',paste0('<th>',heads,'</th>',collapse=''),'</tr></thead><tbody>',paste(rows,collapse=''),'</tbody></table>'))
    }
    den_text<-if(site_mode)"n/N: affected / eligible randomized participants in that site and arm. Records V / P: vaccine / placebo." else paste(vapply(seq_len(nrow(n)),function(i)paste0(n$arm[i],": N=",format(n$N[i],big.mark=","),
      "; ",format(n$any_deviation[i],big.mark=",")," with any deviation; median follow-up ",n$median_followup_days[i]," days"),character(1)),collapse=" | ")
    sections<-c(sections,paste0('<section><h2>',esc(window),'</h2><p class="den">',esc(den_text),'</p>',paste(body,collapse=''),'</section>'))
  }
  status_rows<-apply(status,1,function(r)paste0('<tr>',paste0('<td>',esc(r),'</td>',collapse=''),'</tr>'))
  html<-paste0('<!doctype html><html lang="en"><meta charset="utf-8"><title>Moderna protocol deviation imbalances</title><style>',
    'body{font:14px/1.4 Arial,sans-serif;color:#182d40;max-width:1250px;margin:28px auto;padding:0 18px}h1{font-size:27px;margin-bottom:4px}h2{font-size:21px;border-bottom:2px solid #234c6b;padding-bottom:6px}h3{font-size:15px;margin:18px 0 6px}p{margin:6px 0 12px}.sub,.den{color:#4d6171}table{border-collapse:collapse;width:100%;margin-bottom:16px}th,td{padding:7px 8px;vertical-align:top;border-bottom:1px solid #d4dce3;text-align:right}th{background:#eaf0f5;font-weight:600}th:first-child,td:first-child{text-align:left}tbody tr:nth-child(even){background:#f6f8fa}.signals{table-layout:fixed}.signals td:first-child{overflow-wrap:anywhere}.signals td{font-size:12px}.signals th{font-size:12px}a{color:#20547d}.notes{font-size:12px}li{margin-bottom:6px}.toolbar{margin:15px 0}',
    '@page{size:A4 landscape;margin:12mm}@media print{body{max-width:none;margin:0;padding:0;font-size:10pt;color:#000}h1{font-size:17pt}h2{font-size:13pt}h3{font-size:11pt}section{break-before:page}thead{display:table-header-group}tr{break-inside:avoid}th,td{padding:4px 5px}.signals td,.signals th{font-size:8.5pt}.toolbar{display:none}.notes{font-size:9pt}a{color:inherit;text-decoration:none}}',
    '</style><h1>',title,'</h1><p class="sub">Exploratory comparisons by original randomized arm | ',format(Sys.Date()),' | analyse_deviations.R v1.3.0</p>',
    '<div class="toolbar"><button onclick="window.print()">Print / Save as PDF</button> &nbsp; <a href="',stem,'.pdf">Print-ready PDF</a> | <a href="',all_csv,'">All tests CSV</a> | <a href="',selected_csv,'">Selected rows CSV</a></div>',
    '<h2>Methods and interpretation</h2><ol class="notes">',paste0('<li>',esc(notes),'</li>',collapse=''),'</ol>',
    '<h3>Release availability</h3><table><thead><tr>',paste0('<th>',esc(names(status)),'</th>',collapse=''),'</tr></thead><tbody>',paste(status_rows,collapse=''),'</tbody></table>',
    paste(sections,collapse=''),'</html>')
  writeLines(enc2utf8(html),file.path(out,paste0(stem,".html")),useBytes=TRUE)
}
dv_pdf <- function(results,den,notes,out,alpha,site_mode=FALSE) {
  # Native, selectable vector/text PDF. No browser, LaTeX or extra package needed.
  file<-file.path(out,if(site_mode)"site_fisher_imbalances.pdf" else "deviation_imbalances.pdf")
  font <- if(capabilities("cairo"))"sans" else "Helvetica"
  if(capabilities("cairo")) grDevices::cairo_pdf(file,width=297/25.4,height=210/25.4,family=font,onefile=TRUE) else
    grDevices::pdf(file,width=297/25.4,height=210/25.4,family=font,paper="special",useDingbats=FALSE)
  on.exit(grDevices::dev.off())
  widths<-if(site_mode)c(81,41,41,44,22,22,22) else c(100,29,29,46,24,22,23)
  xs<-12+c(0,cumsum(widths))[1:7]
  page<-0L
  ascii<-function(x){z<-iconv(x,from="UTF-8",to="ASCII//TRANSLIT",sub="?");z[is.na(z)]<-"?";z}
  text<-function(label,x,y,size=8.5,bold=FALSE,just="left",color="#182d40")grid::grid.text(ascii(label),
    x=grid::unit(x,"mm"),y=grid::unit(210-y,"mm"),just=c(just,"top"),gp=grid::gpar(fontsize=size,fontfamily=font,fontface=if(bold)2 else 1,col=color))
  line<-function(y)grid::grid.lines(x=grid::unit(c(12,285),"mm"),y=grid::unit(c(210-y,210-y),"mm"),gp=grid::gpar(col="#cbd5df",lwd=.5))
  newpage<-function(title,subtitle="") {
    grid::grid.newpage();page<<-page+1L;text(if(site_mode)"MODERNA P301  /  SITE-LEVEL DEVIATIONS" else "MODERNA P301  /  PROTOCOL DEVIATIONS",12,9,9,TRUE,color="#476780")
    text(title,12,17,15,TRUE);text(subtitle,12,27,8.5);line(35)
    text("Participant-level exploratory comparison | analyse_deviations.R v1.3.0",12,199,7,color="#52677a")
    text(paste("Page",page),285,199,7,just="right",color="#52677a")
  }
  wrap<-function(s,width,size=8.5) {
    words<-strsplit(ascii(s),"[[:space:]]+")[[1]];lines<-character();current<-""
    for(word in words) {
      candidate<-if(nzchar(current))paste(current,word) else word
      g<-grid::textGrob(candidate,gp=grid::gpar(fontsize=size,fontfamily=font))
      if(nzchar(current)&&grid::convertWidth(grid::grobWidth(g),"mm",valueOnly=TRUE)>width){lines<-c(lines,current);current<-word}else current<-candidate
    };paste(c(lines,current),collapse="\n")
  }
  newpage("Methods and interpretation",paste("Generated",format(Sys.Date()),if(site_mode)"| Tests use participants; the volume filter uses deviation records." else "| Counts are participants, not repeated deviation records."))
  y<-41
  for(i in seq_along(notes)) {
    v<-wrap(paste0(i,". ",notes[i]),269,9);h<-length(strsplit(v,"\n")[[1]])*4.1+4
    if(y+h>191){newpage("Methods and interpretation (continued)");y<-41}
    text(v,12,y,9);y<-y+h
  }
  heads<-if(site_mode)c("Site | Deviation","mRNA-1273\nn/N (%)","Placebo\nn/N (%)","Difference, pp\n[95% CI]","Records\nV / P","Fisher p","Holm p") else c("Recorded deviation","mRNA-1273\nn (%)","Placebo\nn (%)","Difference, pp\n[95% CI]","Test","Raw p","Holm p")
  groups<-if(site_mode)paste0("Fisher p < ",alpha) else c("Holm-adjusted signals","Nominal-only signals")
  for(window in unique(results$window)) {
    z<-results[results$window==window&results$level=="term",,drop=FALSE];n<-den[den$window==window,,drop=FALSE]
    subtitle<-if(site_mode)"n/N = affected / eligible participants within site and arm. Records V / P = vaccine / placebo." else paste(vapply(seq_len(nrow(n)),function(i)paste0(n$arm[i],": N=",n$N[i],"; median follow-up ",n$median_followup_days[i]," d"),character(1)),collapse="   |   ")
    for(group in groups) {
      hit<-!is.na(z$p_raw) & z$p_raw<alpha & if(site_mode)TRUE else if(group=="Holm-adjusted signals")z$p_holm<alpha else z$p_holm>=alpha
      if(site_mode)hit<-hit & z$volume_qualified
      d<-z[hit,,drop=FALSE];d<-d[order(d$p_holm,d$p_raw,-abs(d$difference_pp)),,drop=FALSE];display<-dv_display(d,site_mode)
      header<-function(continued=FALSE){newpage(paste(window,"|",group,if(continued)"(continued)" else ""),subtitle)
        for(j in seq_along(heads))text(heads[j],xs[j]+1,39,8,TRUE)
        line(50)}
      header();y<-53
      if(!nrow(d)){text("No deviations meet this criterion.",13,y,10);next}
      for(i in seq_len(nrow(display))) {
        cells<-vapply(seq_len(ncol(display)),function(j)wrap(display[[j]][i],widths[j]-3,8.5),character(1))
        h<-max(vapply(strsplit(cells,"\n"),length,integer(1)))*3.6+5
        if(y+h>192){header(TRUE);y<-53}
        if(i%%2==0)grid::grid.rect(x=grid::unit(12,"mm"),y=grid::unit(210-y+1,"mm"),width=grid::unit(273,"mm"),height=grid::unit(h,"mm"),just=c("left","top"),gp=grid::gpar(fill="#f2f5f8",col=NA))
        for(j in seq_along(cells))text(cells[j],xs[j]+1,y,8.5)
        y<-y+h;line(y-2)
      }
    }
  }
  invisible(file)
}

# ---- MAIN ----
analyse_deviations <- function(root="",releases=c("EUA","Prod"),adsl_files=c(EUA="",Prod=""),
  dv_files=c(EUA="",Prod=""),addv_files=c(EUA="",Prod=""),output_dir="",
  eua_cutoff="2020-11-25",bla_cutoff="2021-03-26",start_basis="randomization",alpha=.05,
  expected_min=5,min_randomized=30000L,csv_encoding="auto",numeric_date_origin="",excel_fr=TRUE,
  site_min_records=20L,weekly_png_dpi=300L,unblinding_start="2020-12-14") {
  root<-dv_root(root)
  if(!all(releases%in%c("EUA","Prod"))||anyDuplicated(releases))stop("Choose distinct EUA/Prod releases.")
  if(!start_basis%in%c("randomization","first_dose"))stop("Invalid start_basis.")
  if(length(alpha)!=1L||!is.finite(alpha)||alpha<=0||alpha>=1)stop("Invalid alpha.")
  if(length(expected_min)!=1L||!is.finite(expected_min)||expected_min<5)stop("expected_min must be at least 5.")
  if(length(site_min_records)!=1L||!is.finite(site_min_records)||site_min_records<1||site_min_records!=floor(site_min_records))stop("site_min_records must be a positive integer.")
  if(length(weekly_png_dpi)!=1L||!is.finite(weekly_png_dpi)||weekly_png_dpi<72||weekly_png_dpi>1200)stop("weekly_png_dpi must be between 72 and 1200.")
  if(length(unblinding_start)!=1L || is.na(dv_date(unblinding_start,"unblinding_start",strict=TRUE)))stop("Set unblinding_start to a complete YYYY-MM-DD calendar reference.")
  for(cap in c(eua_cutoff,bla_cutoff))if(is.na(dv_date(cap,"cutoff",strict=TRUE)))stop("Missing cutoff.")
  if(as.Date(eua_cutoff)>as.Date(bla_cutoff))stop("EUA cutoff must not exceed BLA cutoff.")
  if(!nzchar(output_dir)) {
    dir.create(file.path(root,"analysis"),showWarnings=FALSE,recursive=TRUE)
    output_dir<-tempfile(paste0("deviations_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_"),tmpdir=file.path(root,"analysis"))
  } else output_dir<-dv_path(root,output_dir)
  if(dir.exists(output_dir)&&length(list.files(output_dir,all.files=TRUE,no..=TRUE)))stop("Output directory is not empty; choose a new directory.")
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  out<-output_dir;status<-list();all<-list();dens<-list();sites<-list();inputs<-list();site_runs<-list();weekly_runs<-list();unblinding_audits<-list()
  explicit<-function(x,r)if(r%in%names(x))unname(x[r]) else ""
  for(release in releases) {
    base<-file.path(root,if(release=="Prod")"prod_data" else "eua_data")
    any_explicit<-any(nzchar(c(explicit(adsl_files,release),explicit(dv_files,release),explicit(addv_files,release))))
    if(!any_explicit && !length(dv_candidates(root,release,"DV")) && !length(dv_candidates(root,release,"ADDV"))) {
      status[[release]]<-data.frame(Release=release,Status="UNAVAILABLE",Detail="No local DV/ADDV source found; not interpreted as zero deviations.")
      next
    }
    run<-tryCatch({
      cutoff<-if(release=="Prod")bla_cutoff else eua_cutoff
      select<-function(k,path,optional=FALSE)dv_select(root,release,k,path,out,cutoff,min_randomized,csv_encoding,numeric_date_origin,optional)
      adsl<-select("ADSL",explicit(adsl_files,release))
      raw<-select("DV",explicit(dv_files,release),TRUE)
      ad<-select("ADDV",explicit(addv_files,release),TRUE)
      context<-dv_context(adsl$data,start_basis,numeric_date_origin)
      events<-dv_records(raw,ad,adsl$data,out,release)
      dv_csv(context,file.path(out,paste0("subject_context_",release,".csv")))
      dv_csv(events,file.path(out,paste0("all_source_records_",release,".csv")))
      dictionary<-unique(events[c("term","DVTERM","DVTERM1","sponsor_category","DVSIG")])
      dv_csv(dictionary,file.path(out,paste0("term_dictionary_",release,".csv")))
      windows<-if(release=="Prod")c(Prod_EUA=eua_cutoff,Prod_BLINDED=bla_cutoff) else c(EUA_SNAPSHOT=eua_cutoff)
      completed<-lapply(names(windows),function(w)dv_window(events,context,release,w,windows[[w]],expected_min,out,unblinding_start))
      full_weekly<-if(release=="Prod")dv_weekly_full_followup(events,context,release,bla_cutoff,out,unblinding_start) else NULL
      unblinding_audit<-dv_unblind_inventory(context,release,unblinding_start,out)
      for(src in list(ADSL=adsl,DV=raw,ADDV=ad))if(!is.null(src))inputs[[length(inputs)+1L]]<-data.frame(release=release,path=src$path,md5=src$md5,rows=nrow(src$data))
      list(windows=completed,full_weekly=full_weekly,unblinding_audit=unblinding_audit,detail=paste0(nrow(adsl$data)," ADSL; ",nrow(events)," deviation records. ",
        if(is.null(raw))"ADDV only: raw DV reconciliation unavailable. " else "",
        if(is.null(ad))"DV only: DVTERM1/DVSIG not available; descriptions may be incomplete. " else "",
        sum(!events$in_adsl)," records outside the ADSL roster; ",sum(events$date_status!="COMPLETE")," unresolved dates."))
    },error=function(e)e)
    if(inherits(run,"error"))status[[release]]<-data.frame(Release=release,Status="FAILED",Detail=conditionMessage(run)) else {
      status[[release]]<-data.frame(Release=release,Status="ANALYZED",Detail=run$detail)
      for(w in run$windows){all[[length(all)+1L]]<-w$results;dens[[length(dens)+1L]]<-w$summary;sites[[length(sites)+1L]]<-w$sites;site_runs[[length(site_runs)+1L]]<-w$site_fisher;weekly_runs[[length(weekly_runs)+1L]]<-w$weekly}
      if(!is.null(run$full_weekly))weekly_runs[[length(weekly_runs)+1L]]<-run$full_weekly
      unblinding_audits[[length(unblinding_audits)+1L]]<-run$unblinding_audit
    }
  }
  status<-do.call(rbind,status);rownames(status)<-NULL;dv_csv(status,file.path(out,"release_status.csv"))
  if(!length(all))stop("No release was analyzable. See ",file.path(out,"release_status.csv"),". ",paste(status$Detail,collapse=" | "))
  results<-do.call(rbind,all);den<-do.call(rbind,dens);results$p_holm<-NA_real_
  for(level in unique(results$level)) {
    j<-which(results$level==level);results$p_holm[j]<-p.adjust(results$p_raw[j],method="holm")
  }
  results$signal<-ifelse(results$p_holm<alpha,"HOLM_ADJUSTED",ifelse(results$p_raw<alpha,"NOMINAL_ONLY","NOT_SELECTED"))
  results<-results[order(results$level,results$window,results$p_holm,results$p_raw,results$deviation),]
  selected<-results[results$level=="term"&results$p_raw<alpha,,drop=FALSE]
  dv_csv(results,file.path(out,"deviations_all_tests.csv"));dv_csv(selected,file.path(out,"deviations_imbalanced.csv"))
  if(excel_fr)dv_csv(selected,file.path(out,"deviations_imbalanced_Excel_FR.csv"),sep=";")
  dv_csv(den,file.path(out,"denominators.csv"))
  if(length(sites))dv_csv(do.call(rbind,sites),file.path(out,"site_term_counts.csv"))
  if(length(inputs))dv_csv(do.call(rbind,inputs),file.path(out,"source_manifest.csv"))
  ntests<-sum(results$level=="term")
  notes<-dv_notes(start_basis,alpha,expected_min,ntests)
  if(any(status$Status!="ANALYZED"))notes<-c(paste("Release coverage:",paste(paste(status$Release,status$Status),collapse="; "),
    ". Only analyzed releases enter multiplicity correction; see release_status.csv for reasons."),notes)
  dv_html(results,den,status,notes,out,alpha);dv_pdf(results,den,notes,out,alpha)
  site_results<-dv_site_finish(site_runs,out,status,alpha,site_min_records,excel_fr,start_basis)
  weekly_results<-dv_weekly_finish(weekly_runs,out,status,eua_cutoff,weekly_png_dpi,excel_fr,start_basis,unblinding_start,do.call(rbind,unblinding_audits))
  writeLines(c("MODERNA P301 DEVIATION ANALYSIS v1.3.0","",notes,"",
    "OUTPUTS: deviation_imbalances.pdf (A4 landscape); deviation_imbalances.html (print CSS); deviations_all_tests.csv (all results, including nonsignificant); deviations_imbalanced.csv (nominal p<alpha term rows); denominators.csv; release_status.csv; source_manifest.csv.",
    "SITE OUTPUTS: site_fisher_imbalances.pdf/html; site_fisher_all_tests.csv; site_fisher_volume_qualified.csv; site_fisher_imbalanced.csv; site_fisher_coverage.csv; site_deviation_mapping.csv; site_deviation_rules.csv; site_fisher_methods.txt. These cover the seven deviation types specified in the earlier audit, jointly and individually. The volume filter counts records; Fisher tests count affected participants against within-site randomized-arm denominators. No global-significance prefilter is applied. See site_fisher_methods.txt for multiplicity and interpretation.",
    "WEEKLY OUTPUTS: weekly_staff_deviations.pdf (A4 landscape vector charts); weekly_staff_deviations_*.png (one page per image, default 300 dpi); weekly_staff_deviations.csv (records, distinct affected participants and follow-up by arm/week/type); weekly_staff_deviations_coverage.csv; weekly_staff_deviations_methods.txt; weekly_staff_deviations_pages.csv. All sites and the seven specified types are retained, without a volume/significance filter. Combined and type-specific counts use Monday-Sunday bins by recorded deviation date, with original windows/caps. Separate windows/releases are never summed.",
    "UNBLINDING OUTPUTS: dashed red calendar reference at PHMPT_UNBLINDING_START; deviations_before_after_unblinding.csv (exact-date calendar and participant-UBLDDT splits); unblinding_reference_audit.csv; unblinding_daily_*.csv; unblinding_subject_dates_*.csv. New Prod_ALL_FOLLOWUP charts/records/population retain follow-up through unblinding/crossover up to administrative/study-end/death caps. Original global/site tests remain restricted to their original windows. See weekly_staff_deviations_methods.txt for distinct calendar versus participant-level interpretations and missing-date handling.",
    "Audit exports preserve all selected source events, source row/event keys, full descriptions, window inclusion reasons, source-only records, cohort exclusions, date resolution, sponsor significance and site counts. Separate releases are never concatenated into a single participant pool.",
    "Missing AP01EDT/AP02SDT/UBLDDT uses the administrative cutoff and available study-end/death caps; this limitation is counted in denominators.csv and flagged per subject. Off-study vaccination dates do not define protocol-deviation follow-up here.",
    "Input matching compares the columns used in this analysis; duplicate byte copies and equivalent formats are logged. Different valid snapshots require explicit PHMPT_*_FILES paths. ADSL randomized count below PHMPT_MIN_RANDOMIZED is rejected as a likely subset; this guard does not prove dataset completeness.",
    "KNOWN PRODUCTION SOURCES: 125752_S3_M5_mrna-1273-p301_A_D_adsl.jmp (md-production-050125.zip); 125752_S3_M5_mrna-1273-p301_S_dv.xpt (md-production-010225.zip); 125752_S3_M5_mrna-1273-p301_A_D_addv.xpt (md-production-020325.zip).",
    "This script intentionally does not reproduce the earlier audit's treated-subject population. Default is randomized participants from RANDDT; PHMPT_START='first_dose' selects a post-dose cohort from TR01SDT. Date caps and complete ADSL provenance are exported.",
    "METHOD REFERENCES: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/chisq.test.html ; https://stat.ethz.ch/R-manual/R-devel/library/stats/html/fisher.test.html ; https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html ; Newcombe RG (1998), Statistics in Medicine 17:873-890, https://pubmed.ncbi.nlm.nih.gov/9595617/",
    "The PDF transliterates non-ASCII text for portable base-R fonts. HTML/CSV retain original UTF-8 wording."),file.path(out,"README.txt"),useBytes=TRUE)
  dput(list(version="1.3.0",root=root,releases=releases,eua_cutoff=eua_cutoff,bla_cutoff=bla_cutoff,start_basis=start_basis,
    alpha=alpha,expected_min=expected_min,min_randomized=min_randomized,csv_encoding=csv_encoding,numeric_date_origin=numeric_date_origin,site_min_records=site_min_records,weekly_png_dpi=weekly_png_dpi,unblinding_start=unblinding_start),file.path(out,"run_settings.R"))
  capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
  fs<-list.files(out,full.names=TRUE);dv_csv(data.frame(file=basename(fs),md5=unname(tools::md5sum(fs))),file.path(out,"output_checksums.csv"))
  message("Wrote ",nrow(selected)," nominal term signals; ",sum(selected$p_holm<alpha)," survive Holm correction.\n",normalizePath(out,winslash="/"))
  if(any(status$Status=="FAILED"))warning("Some requested releases FAILED: review release_status.csv. Outputs cover analyzed releases only.",call.=FALSE)
  invisible(list(results=results,denominators=den,status=status,site_results=site_results,weekly_results=weekly_results,output_dir=out))
}

# The validated PHMPT P301 ADSL JMP reader is embedded below.

# Binary-layout references: JMPReader.jl; jmpio-python commit
# 5e9d49f07f052e266c183f9e19a8ce087367324a (fixed-width character / Float64 dates only).
dv_read_jmp <- function(path, fields) {
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

# ---- RUN: select all and Run, or Source ----
analyse_deviations(root=PHMPT_ROOT, releases=PHMPT_RELEASES,
  adsl_files=PHMPT_ADSL_FILES, dv_files=PHMPT_DV_FILES, addv_files=PHMPT_ADDV_FILES,
  output_dir=PHMPT_OUTPUT_DIR, eua_cutoff=PHMPT_EUA_CUTOFF, bla_cutoff=PHMPT_BLA_CUTOFF,
  start_basis=PHMPT_START, alpha=PHMPT_ALPHA, expected_min=PHMPT_EXPECTED_MIN,
  min_randomized=PHMPT_MIN_RANDOMIZED, csv_encoding=PHMPT_CSV_ENCODING,
  numeric_date_origin=PHMPT_NUMERIC_DATE_ORIGIN, excel_fr=PHMPT_EXCEL_FR,
  site_min_records=PHMPT_SITE_MIN_RECORDS,weekly_png_dpi=PHMPT_WEEKLY_PNG_DPI,
  unblinding_start=PHMPT_UNBLINDING_START)
