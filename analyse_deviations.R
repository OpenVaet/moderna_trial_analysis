# Moderna P301 protocol deviations -- v1.0.0 -- R 4.3+
# Standalone: save in the project root, Select All / Run or Source in RStudio.
# Reads your local PHMPT files. No downloads and no earlier analysis scripts.
# Base R + grid suffice for CSV/JMP and PDF/HTML output; haven reads XPT/SAS7BDAT.
# Tests compare participants with >=1 recorded deviation, not repeated records.
# Primary grouping: exact DVTERM + DVTERM1 (whitespace normalized only).

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
dv_compare <- function(a,nv,c,np,expected_min=5) {
  if(any(!is.finite(c(a,nv,c,np)))||nv<=0||np<=0||a<0||c<0||a>nv||c>np)stop("Invalid 2x2 counts.")
  tab <- matrix(c(a,nv-a,c,np-c),nrow=2L,byrow=TRUE,
    dimnames=list(arm=c("mRNA-1273","Placebo"),outcome=c("Affected","Not affected")))
  expected <- outer(rowSums(tab),colSums(tab))/sum(tab)
  if(any(colSums(tab)==0)) {p<-1;method<-"No variation";statistic<-NA_real_} else
  if(all(expected>=expected_min)) {
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
dv_window <- function(events,context,release,window,cap,expected_min,out) {
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
  list(results=do.call(rbind,result),sites=if(length(sites))do.call(rbind,sites) else NULL,summary=summary)
}

# ---- PRINTABLE TABLES ----
dv_escape <- function(x) {
  x<-as.character(x);x[is.na(x)]<-"";x<-gsub("&","&amp;",x,fixed=TRUE);x<-gsub("<","&lt;",x,fixed=TRUE)
  x<-gsub(">","&gt;",x,fixed=TRUE);gsub('"',"&quot;",x,fixed=TRUE)
}
dv_p <- function(x)ifelse(is.na(x),"-",ifelse(x<.0001,formatC(x,format="e",digits=2),formatC(x,format="f",digits=4)))
dv_display <- function(x) {
  data.frame(Deviation=x$deviation,
    Vaccine=sprintf("%s (%.2f%%)",format(x$vaccine_n,big.mark=",",trim=TRUE),x$vaccine_percent),
    Placebo=sprintf("%s (%.2f%%)",format(x$placebo_n,big.mark=",",trim=TRUE),x$placebo_percent),
    Difference=sprintf("%+.2f [%+.2f, %+.2f]",x$difference_pp,x$lower95_pp,x$upper95_pp),
    Test=x$test,P=dv_p(x$p_raw),Holm=dv_p(x$p_holm),stringsAsFactors=FALSE)
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
dv_html <- function(results,den,status,notes,out,alpha) {
  esc<-dv_escape;sections<-character()
  for(window in unique(results$window)) {
    z<-results[results$window==window&results$level=="term",,drop=FALSE]
    n<-den[den$window==window,,drop=FALSE]
    heads<-c("Recorded deviation","mRNA-1273<br>n (%)","Placebo<br>n (%)","Difference, pp<br>[95% CI]","Test","Raw p","Holm p")
    body<-character()
    for(group in c("Holm-adjusted signals","Nominal-only signals")) {
      hit<-z$p_raw<alpha & if(group=="Holm-adjusted signals")z$p_holm<alpha else z$p_holm>=alpha
      d<-z[hit,,drop=FALSE];d<-d[order(d$p_holm,d$p_raw,-abs(d$difference_pp)),,drop=FALSE]
      rows<-if(nrow(d))apply(dv_display(d),1,function(r)paste0('<tr>',paste0('<td>',esc(r),'</td>',collapse=''),'</tr>')) else
        '<tr><td colspan="7">No deviations meet this criterion.</td></tr>'
      body<-c(body,paste0('<h3>',group,'</h3><table class="signals"><colgroup><col style="width:38%"><col style="width:11%"><col style="width:11%"><col style="width:18%"><col style="width:8%"><col style="width:7%"><col style="width:7%"></colgroup><thead><tr>',paste0('<th>',heads,'</th>',collapse=''),'</tr></thead><tbody>',paste(rows,collapse=''),'</tbody></table>'))
    }
    den_text<-paste(vapply(seq_len(nrow(n)),function(i)paste0(n$arm[i],": N=",format(n$N[i],big.mark=","),
      "; ",format(n$any_deviation[i],big.mark=",")," with any deviation; median follow-up ",n$median_followup_days[i]," days"),character(1)),collapse=" | ")
    sections<-c(sections,paste0('<section><h2>',esc(window),'</h2><p class="den">',esc(den_text),'</p>',paste(body,collapse=''),'</section>'))
  }
  status_rows<-apply(status,1,function(r)paste0('<tr>',paste0('<td>',esc(r),'</td>',collapse=''),'</tr>'))
  html<-paste0('<!doctype html><html lang="en"><meta charset="utf-8"><title>Moderna protocol deviation imbalances</title><style>',
    'body{font:14px/1.4 Arial,sans-serif;color:#182d40;max-width:1250px;margin:28px auto;padding:0 18px}h1{font-size:27px;margin-bottom:4px}h2{font-size:21px;border-bottom:2px solid #234c6b;padding-bottom:6px}h3{font-size:15px;margin:18px 0 6px}p{margin:6px 0 12px}.sub,.den{color:#4d6171}table{border-collapse:collapse;width:100%;margin-bottom:16px}th,td{padding:7px 8px;vertical-align:top;border-bottom:1px solid #d4dce3;text-align:right}th{background:#eaf0f5;font-weight:600}th:first-child,td:first-child{text-align:left}tbody tr:nth-child(even){background:#f6f8fa}.signals{table-layout:fixed}.signals td:first-child{overflow-wrap:anywhere}.signals td{font-size:12px}.signals th{font-size:12px}a{color:#20547d}.notes{font-size:12px}li{margin-bottom:6px}.toolbar{margin:15px 0}',
    '@page{size:A4 landscape;margin:12mm}@media print{body{max-width:none;margin:0;padding:0;font-size:10pt;color:#000}h1{font-size:17pt}h2{font-size:13pt}h3{font-size:11pt}section{break-before:page}thead{display:table-header-group}tr{break-inside:avoid}th,td{padding:4px 5px}.signals td,.signals th{font-size:8.5pt}.toolbar{display:none}.notes{font-size:9pt}a{color:inherit;text-decoration:none}}',
    '</style><h1>Moderna P301: recorded deviation imbalances</h1><p class="sub">Exploratory comparisons by original randomized arm | ',format(Sys.Date()),' | analyse_deviations.R v1.0.0</p>',
    '<div class="toolbar"><button onclick="window.print()">Print / Save as PDF</button> &nbsp; <a href="deviation_imbalances.pdf">Print-ready PDF</a> | <a href="deviations_all_tests.csv">All tests CSV</a> | <a href="deviations_imbalanced.csv">Selected rows CSV</a></div>',
    '<h2>Methods and interpretation</h2><ol class="notes">',paste0('<li>',esc(notes),'</li>',collapse=''),'</ol>',
    '<h3>Release availability</h3><table><thead><tr>',paste0('<th>',esc(names(status)),'</th>',collapse=''),'</tr></thead><tbody>',paste(status_rows,collapse=''),'</tbody></table>',
    paste(sections,collapse=''),'</html>')
  writeLines(enc2utf8(html),file.path(out,"deviation_imbalances.html"),useBytes=TRUE)
}
dv_pdf <- function(results,den,notes,out,alpha) {
  # Native, selectable vector/text PDF. No browser, LaTeX or extra package needed.
  file<-file.path(out,"deviation_imbalances.pdf")
  font <- if(capabilities("cairo"))"sans" else "Helvetica"
  if(capabilities("cairo")) grDevices::cairo_pdf(file,width=297/25.4,height=210/25.4,family=font,onefile=TRUE) else
    grDevices::pdf(file,width=297/25.4,height=210/25.4,family=font,paper="special",useDingbats=FALSE)
  on.exit(grDevices::dev.off())
  widths<-c(100,29,29,46,24,22,23);xs<-12+c(0,cumsum(widths))[1:7]
  page<-0L
  ascii<-function(x){z<-iconv(x,from="UTF-8",to="ASCII//TRANSLIT",sub="?");z[is.na(z)]<-"?";z}
  text<-function(label,x,y,size=8.5,bold=FALSE,just="left",color="#182d40")grid::grid.text(ascii(label),
    x=grid::unit(x,"mm"),y=grid::unit(210-y,"mm"),just=c(just,"top"),gp=grid::gpar(fontsize=size,fontfamily=font,fontface=if(bold)2 else 1,col=color))
  line<-function(y)grid::grid.lines(x=grid::unit(c(12,285),"mm"),y=grid::unit(c(210-y,210-y),"mm"),gp=grid::gpar(col="#cbd5df",lwd=.5))
  newpage<-function(title,subtitle="") {
    grid::grid.newpage();page<<-page+1L;text("MODERNA P301  /  PROTOCOL DEVIATIONS",12,9,9,TRUE,color="#476780")
    text(title,12,17,15,TRUE);text(subtitle,12,27,8.5);line(35)
    text("Participant-level exploratory comparison | analyse_deviations.R v1.0.0",12,199,7,color="#52677a")
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
  newpage("Methods and interpretation",paste("Generated",format(Sys.Date()),"| Counts are participants, not repeated deviation records."))
  y<-41
  for(i in seq_along(notes)) {
    v<-wrap(paste0(i,". ",notes[i]),269,9);h<-length(strsplit(v,"\n")[[1]])*4.1+4
    if(y+h>191){newpage("Methods and interpretation (continued)");y<-41}
    text(v,12,y,9);y<-y+h
  }
  heads<-c("Recorded deviation","mRNA-1273\nn (%)","Placebo\nn (%)","Difference, pp\n[95% CI]","Test","Raw p","Holm p")
  for(window in unique(results$window)) {
    z<-results[results$window==window&results$level=="term",,drop=FALSE];n<-den[den$window==window,,drop=FALSE]
    subtitle<-paste(vapply(seq_len(nrow(n)),function(i)paste0(n$arm[i],": N=",n$N[i],"; median follow-up ",n$median_followup_days[i]," d"),character(1)),collapse="   |   ")
    for(group in c("Holm-adjusted signals","Nominal-only signals")) {
      hit<-z$p_raw<alpha & if(group=="Holm-adjusted signals")z$p_holm<alpha else z$p_holm>=alpha
      d<-z[hit,,drop=FALSE];d<-d[order(d$p_holm,d$p_raw,-abs(d$difference_pp)),,drop=FALSE];display<-dv_display(d)
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
  expected_min=5,min_randomized=30000L,csv_encoding="auto",numeric_date_origin="",excel_fr=TRUE) {
  root<-dv_root(root)
  if(!all(releases%in%c("EUA","Prod"))||anyDuplicated(releases))stop("Choose distinct EUA/Prod releases.")
  if(!start_basis%in%c("randomization","first_dose"))stop("Invalid start_basis.")
  if(length(alpha)!=1L||!is.finite(alpha)||alpha<=0||alpha>=1)stop("Invalid alpha.")
  if(length(expected_min)!=1L||!is.finite(expected_min)||expected_min<5)stop("expected_min must be at least 5.")
  for(cap in c(eua_cutoff,bla_cutoff))if(is.na(dv_date(cap,"cutoff",strict=TRUE)))stop("Missing cutoff.")
  if(as.Date(eua_cutoff)>as.Date(bla_cutoff))stop("EUA cutoff must not exceed BLA cutoff.")
  if(!nzchar(output_dir)) {
    dir.create(file.path(root,"analysis"),showWarnings=FALSE,recursive=TRUE)
    output_dir<-tempfile(paste0("deviations_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_"),tmpdir=file.path(root,"analysis"))
  } else output_dir<-dv_path(root,output_dir)
  if(dir.exists(output_dir)&&length(list.files(output_dir,all.files=TRUE,no..=TRUE)))stop("Output directory is not empty; choose a new directory.")
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  out<-output_dir;status<-list();all<-list();dens<-list();sites<-list();inputs<-list()
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
      completed<-lapply(names(windows),function(w)dv_window(events,context,release,w,windows[[w]],expected_min,out))
      for(src in list(ADSL=adsl,DV=raw,ADDV=ad))if(!is.null(src))inputs[[length(inputs)+1L]]<-data.frame(release=release,path=src$path,md5=src$md5,rows=nrow(src$data))
      list(windows=completed,detail=paste0(nrow(adsl$data)," ADSL; ",nrow(events)," deviation records. ",
        if(is.null(raw))"ADDV only: raw DV reconciliation unavailable. " else "",
        if(is.null(ad))"DV only: DVTERM1/DVSIG not available; descriptions may be incomplete. " else "",
        sum(!events$in_adsl)," records outside the ADSL roster; ",sum(events$date_status!="COMPLETE")," unresolved dates."))
    },error=function(e)e)
    if(inherits(run,"error"))status[[release]]<-data.frame(Release=release,Status="FAILED",Detail=conditionMessage(run)) else {
      status[[release]]<-data.frame(Release=release,Status="ANALYZED",Detail=run$detail)
      for(w in run$windows){all[[length(all)+1L]]<-w$results;dens[[length(dens)+1L]]<-w$summary;sites[[length(sites)+1L]]<-w$sites}
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
  writeLines(c("MODERNA P301 DEVIATION ANALYSIS v1.0.0","",notes,"",
    "OUTPUTS: deviation_imbalances.pdf (A4 landscape); deviation_imbalances.html (print CSS); deviations_all_tests.csv (all results, including nonsignificant); deviations_imbalanced.csv (nominal p<alpha term rows); denominators.csv; release_status.csv; source_manifest.csv.",
    "Audit exports preserve all selected source events, source row/event keys, full descriptions, window inclusion reasons, source-only records, cohort exclusions, date resolution, sponsor significance and site counts. Separate releases are never concatenated into a single participant pool.",
    "Missing AP01EDT/AP02SDT/UBLDDT uses the administrative cutoff and available study-end/death caps; this limitation is counted in denominators.csv and flagged per subject. Off-study vaccination dates do not define protocol-deviation follow-up here.",
    "Input matching compares the columns used in this analysis; duplicate byte copies and equivalent formats are logged. Different valid snapshots require explicit PHMPT_*_FILES paths. ADSL randomized count below PHMPT_MIN_RANDOMIZED is rejected as a likely subset; this guard does not prove dataset completeness.",
    "KNOWN PRODUCTION SOURCES: 125752_S3_M5_mrna-1273-p301_A_D_adsl.jmp (md-production-050125.zip); 125752_S3_M5_mrna-1273-p301_S_dv.xpt (md-production-010225.zip); 125752_S3_M5_mrna-1273-p301_A_D_addv.xpt (md-production-020325.zip).",
    "This script intentionally does not reproduce the earlier audit's treated-subject population. Default is randomized participants from RANDDT; PHMPT_START='first_dose' selects a post-dose cohort from TR01SDT. Date caps and complete ADSL provenance are exported.",
    "METHOD REFERENCES: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/chisq.test.html ; https://stat.ethz.ch/R-manual/R-devel/library/stats/html/fisher.test.html ; https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html ; Newcombe RG (1998), Statistics in Medicine 17:873-890, https://pubmed.ncbi.nlm.nih.gov/9595617/",
    "The PDF transliterates non-ASCII text for portable base-R fonts. HTML/CSV retain original UTF-8 wording."),file.path(out,"README.txt"),useBytes=TRUE)
  dput(list(version="1.0.0",root=root,releases=releases,eua_cutoff=eua_cutoff,bla_cutoff=bla_cutoff,start_basis=start_basis,
    alpha=alpha,expected_min=expected_min,min_randomized=min_randomized,csv_encoding=csv_encoding,numeric_date_origin=numeric_date_origin),file.path(out,"run_settings.R"))
  capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
  fs<-list.files(out,full.names=TRUE);dv_csv(data.frame(file=basename(fs),md5=unname(tools::md5sum(fs))),file.path(out,"output_checksums.csv"))
  message("Wrote ",nrow(selected)," nominal term signals; ",sum(selected$p_holm<alpha)," survive Holm correction.\n",normalizePath(out,winslash="/"))
  if(any(status$Status=="FAILED"))warning("Some requested releases FAILED: review release_status.csv. Outputs cover analyzed releases only.",call.=FALSE)
  invisible(list(results=results,denominators=den,status=status,output_dir=out))
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
  numeric_date_origin=PHMPT_NUMERIC_DATE_ORIGIN, excel_fr=PHMPT_EXCEL_FR)
