# =============================================================================
# Moderna P301 | Day 57 anti-N / PCR ascertainment audit | v1.0.0 | 2026-10-21
# Standalone R >= 4.3. Select All / Run in RStudio from your project.
# CSV/JMP and reports: base R only. XPT/SAS7BDAT require installed 'haven'.
# No downloads, installations, efficacy estimates or automatic infection imputation.
# Prior complete reconciled export OR production ADSL + IS (+ADIS) + MB (+ADMB).
# Explicit export can be its extracted folder or the original audit ZIP.
# options(phmpt.autorun=FALSE) loads functions without executing the analysis.
# =============================================================================
PHMPT_ROOT <- getOption('phmpt.root', '')
PHMPT_EXPORT <- getOption('phmpt.export', '')
PHMPT_ADSL_FILE <- ''
PHMPT_IS_FILE <- ''
PHMPT_ADIS_FILE <- ''
PHMPT_MB_FILE <- ''
PHMPT_ADMB_FILE <- ''
PHMPT_OUTPUT_DIR <- getOption('phmpt.output', '')
PHMPT_EUA_CUTOFF <- '2020-11-25' # calendar restriction of BLA export, not EUA snapshot
PHMPT_BLA_CUTOFF <- '2021-03-26'
PHMPT_DAY57_EARLY <- 3L         # analyst visit window, not an asserted protocol rule
PHMPT_DAY57_LATE <- 7L
PHMPT_BASELINE_LOOKBACK <- 28L
PHMPT_MIN_RANDOMIZED <- 30000L
PHMPT_CSV_ENCODING <- 'auto'
PHMPT_NUMERIC_DATE_ORIGIN <- ''
PHMPT_Q_GRID <- seq(0.05, 1, by=0.05) # illustrative detection assumptions, not a CI
PHMPT_SPECIFICITY_GRID <- c(1, 0.999, 0.995) # illustrative, not assay validation estimates
PHMPT_SHOW_REPORT <- interactive()


nt_bool <- function(x) ws_norm(x)%in%c("TRUE","Y","YES","1")
nt_result <- function(x) {
  x<-ws_norm(x);r<-rep("UNKNOWN",length(x))
  r[x%in%c("NEGATIVE","NOT DETECTED","NONREACTIVE","NON-REACTIVE","RNA NOT DETECTED")]<-"NEGATIVE"
  r[x%in%c("POSITIVE","DETECTED","REACTIVE","RNA DETECTED","PRESUMPTIVE POSITIVE")]<-"POSITIVE"
  r[x%in%c("INDETERMINATE","EQUIVOCAL","INVALID","INCONCLUSIVE")]<-"INDETERMINATE"
  r
}
nt_dates_out <- function(x) {
  cols<-intersect(c("date","baseline_date","target_date","window_start","window_end","week","TR01SDT","DOSE2DT",
    "RANDDT","CUTOFFDT","AP01EDT","UBLDDT","AP02SDT","TR02SDT","EUAVACDT","EOSDT","DTHDT","PDVDT",
    "end_blinded","end_full","date_min","date_max"),names(x))
  for(k in cols)if(is.numeric(x[[k]]))x[[k]]<-ws_iso(x[[k]])
  x
}
nt_optional <- function(root,domain,path,out,encoding,origin,min_n,submission) {
  if(!nzchar(path)) {
    p<-ws_candidates(root,"Prod",domain)
    if(nzchar(submission))p<-p[ws_submission(p)==submission]
    if(!length(p))return(NULL)
  }
  ws_pick(root,"Prod",domain,path,out,encoding,origin,min_n,submission)
}
nt_read_export <- function(path,kind,out) {
  if(!file.exists(path))stop("Reconciled export not found: ",path)
  x<-read.csv(path,colClasses="character",check.names=FALSE,fileEncoding="UTF-8-BOM",na.strings=c("","NA"));x[]<-lapply(x,ws_text)
  need<-c("record_key","subject","sequence","raw_present","analysis_present","source_refs","positive","negative","not_done",
    "derived","result_issue","date_issue","baseline_flag","visits","methods","result","date","method_ok","raw_values")
  if(!all(need%in%names(x)))stop("Expected the validated ",kind,"_records.csv export; missing: ",paste(setdiff(need,names(x)),collapse=", "))
  if(anyDuplicated(x$record_key))stop("Duplicate keys in reconciled export.")
  raw_visit<-vapply(strsplit(x$visits," | ",fixed=TRUE),function(z)if(length(z))z[1] else "",character(1))
  p<-data.frame(record_id=x$record_key,subject=ws_sid(x$subject),sequence=x$sequence,
    date=ws_dates(x$date,"export.date",strict=TRUE),date_source="Previously reconciled collection date",
    raw_visit=raw_visit,analysis_visit=x$visits,raw_result=x$raw_values,analysis_result="",
    positive=nt_bool(x$positive),negative=nt_bool(x$negative),indeterminate=x$result=="INDETERMINATE",
    not_done=nt_bool(x$not_done),result_conflict=nt_bool(x$result_issue),date_conflict=nt_bool(x$date_issue),
    derived=nt_bool(x$derived),baseline_flag=nt_bool(x$baseline_flag),method=x$methods,
    raw_present=nt_bool(x$raw_present),analysis_present=nt_bool(x$analysis_present),source_refs=x$source_refs,
    method_ok=if(kind=="serology")TRUE else nt_bool(x$method_ok),stringsAsFactors=FALSE)
  nt_finish_records(p,kind,out)
}
nt_reconcile <- function(raw,analysis,kind,origin,out) {
  sero<-kind=="serology";pref<-if(sero)"IS" else "MB"
  code<-paste0(pref,"TESTCD");sq<-paste0(pref,"SEQ");dt<-paste0(pref,"DTC")
  raw$.source_row<-seq_len(nrow(raw));r<-raw[ws_norm(raw[[code]])=="SARSCOV2",,drop=FALSE]
  if(!nrow(r))stop("No SARSCOV2 records in ",pref)
  if(sero) {
    label<-ws_norm(ws_col(r,"ISTEST"))
    if(any(grepl("SPIKE|ANTI[- ]S($| )",label)&!grepl("NUCLEO",label)))stop("SARSCOV2 label appears spike-specific; inspect IS assay identity.")
  }
  raw_key<-paste(r$subject,r[[sq]],sep="|")
  raw_key[!nzchar(r[[sq]])]<-paste(r$subject[!nzchar(r[[sq]])],"RAW_MISSINGSEQ",r$.source_row[!nzchar(r[[sq]])],sep="|")
  if(anyDuplicated(raw_key))stop("Repeated raw ",pref," subject/sequence; resolve the source snapshot.")
  a<-NULL
  if(!is.null(analysis)) {
    analysis$.source_row<-seq_len(nrow(analysis));a<-analysis[ws_norm(analysis$PARAMCD)=="SARSCOV2",,drop=FALSE]
    if(sero&&any(!grepl("ELECSYS",ws_norm(a$PARAM))))stop("Unexpected ADIS SARSCOV2 label: expected the Elecsys anti-N assay.")
    ws_csv(a[nzchar(a$DTYPE),,drop=FALSE],file.path(out,paste0(kind,"_excluded_derived_analysis.csv")))
    a<-a[!nzchar(a$DTYPE),,drop=FALSE]
    a$key<-paste(a$subject,a[[sq]],sep="|")
    a$key[!nzchar(a[[sq]])]<-paste(a$subject[!nzchar(a[[sq]])],"ANALYSIS_MISSINGSEQ",a$.source_row[!nzchar(a[[sq]])],sep="|")
    if(anyDuplicated(a$key))stop("Repeated analysis ",pref," subject/sequence; resolve the source snapshot.")
  }
  keys<-unique(c(raw_key,if(!is.null(a))a$key));ri<-match(keys,raw_key);ai<-if(is.null(a))rep(NA_integer_,length(keys)) else match(keys,a$key)
  val<-function(x,j,k){if(is.null(x)||!k%in%names(x))return(rep("",length(keys)));z<-x[[k]][j];z[is.na(z)]<-"";z}
  rv<-function(k)val(r,ri,k);av<-function(k)val(a,ai,k)
  first<-function(x,y)ifelse(nzchar(x),x,y)
  d1<-ws_dates(rv(dt),paste0(pref,".",dt),origin);d2<-ws_dates(av(dt),paste0("AD",pref,".",dt),origin)
  d3<-ws_dates(av("ADT"),paste0("AD",pref,".ADT"),origin)
  d<-ifelse(is.finite(d1),d1,ifelse(is.finite(d2),d2,d3))
  bad_date<-(is.finite(d1)&is.finite(d2)&d1!=d2)|(is.finite(d1)&is.finite(d3)&d1!=d3)|(is.finite(d2)&is.finite(d3)&d2!=d3)
  partial<-(nzchar(rv(dt))&!is.finite(d1))|(nzchar(av(dt))&!is.finite(d2))
  # Do not overwrite an explicitly partial/malformed collection date with ADT.
  bad_date<-bad_date|partial
  c1<-nt_result(rv(paste0(pref,"STRESC")));c2<-nt_result(rv(paste0(pref,"ORRES")))
  c3<-nt_result(av("AVALC"));c4<-if(sero)nt_result(av("ISSTRESC")) else rep("UNKNOWN",length(keys))
  pos<-c1=="POSITIVE"|c2=="POSITIVE"|c3=="POSITIVE"|c4=="POSITIVE"
  neg<-c1=="NEGATIVE"|c2=="NEGATIVE"|c3=="NEGATIVE"|c4=="NEGATIVE"
  ind<-c1=="INDETERMINATE"|c2=="INDETERMINATE"|c3=="INDETERMINATE"|c4=="INDETERMINATE"
  nd<-ws_norm(rv(paste0(pref,"STAT")))%in%c("NOT DONE","NOT PERFORMED","CANCELLED","CANCELED")
  meth<-rv(paste0(pref,"METHOD"))
  p<-data.frame(record_id=keys,subject=first(rv("subject"),av("subject")),sequence=first(rv(sq),av(sq)),date=d,
    date_source=ifelse(is.finite(d1),paste0(pref,".",dt),ifelse(is.finite(d2),paste0("AD",pref,".",dt),paste0("AD",pref,".ADT"))),
    raw_visit=first(rv("VISIT"),av("VISIT")),analysis_visit=av("AVISIT"),
    raw_result=paste(rv(paste0(pref,"STRESC")),rv(paste0(pref,"ORRES")),sep=" / "),analysis_result=av("AVALC"),
    positive=pos,negative=neg,indeterminate=ind,not_done=nd,result_conflict=(pos&neg)|(nd&(pos|neg|ind)),
    date_conflict=bad_date,derived=FALSE,baseline_flag=nt_bool(av("ABLFL")),method=meth,
    raw_present=!is.na(ri),analysis_present=!is.na(ai),
    source_refs=paste0(ifelse(!is.na(ri),paste0(pref,":",rv(".source_row")),""),";",ifelse(!is.na(ai),paste0("AD",pref,":",av(".source_row")),"")),
    method_ok=if(sero)TRUE else grepl("PCR|POLYMERASE",ws_norm(meth))&!grepl("ANTIGEN|ANTIBODY|SEROLOG",ws_norm(meth)),stringsAsFactors=FALSE)
  nt_finish_records(p,kind,out)
}
nt_finish_records <- function(p,kind,out) {
  v<-ws_norm(p$raw_visit)
  p$open_label<-grepl("^OL-|OPEN LABEL",v)
  p$pdv<-grepl("PARTICIPANT DECISION|OL-D1($| )|OPEN LABEL DAY 1($| )",v)
  p$day1_label<-grepl("^VISIT 1 DAY 1($| )|^DAY 1$|^BASELINE$",v)&!p$open_label&!p$pdv
  p$day57_label<-grepl("^VISIT 3 DAY 57($| )|^DAY 57$",v)&!p$open_label&!p$pdv
  p$dated<-is.finite(p$date)&!p$date_conflict
  p$tested<-p$dated&!p$not_done&!p$derived&p$method_ok&(p$positive|p$negative|p$indeterminate)
  # A positive/negative discrepancy still documents an assay, but cannot establish
  # a negative baseline or an interpretable binary result.
  p$binary_result<-p$tested&!p$result_conflict&(p$positive|p$negative)
  if(any(p$binary_result)) {
    key<-paste(p$subject,p$date,sep="|")
    conflict<-intersect(key[p$binary_result&p$positive],key[p$binary_result&p$negative])
    if(length(conflict)){bad<-key%in%conflict;p$result_conflict[bad]<-TRUE;p$binary_result[bad]<-FALSE}
  }
  ws_csv(nt_dates_out(p),file.path(out,paste0(kind,"_records_reconciled.csv")))
  ws_csv(nt_dates_out(p[!p$tested|p$result_conflict,,drop=FALSE]),file.path(out,paste0(kind,"_quality_review.csv")))
  p
}

ws_norm <- function(x) toupper(trimws(gsub("[[:space:]\\p{Z}]+"," ",ws_text(x),perl=TRUE)))
ws_col <- function(x,n) if(n %in% names(x))ws_text(x[[n]]) else rep("",nrow(x))
ws_week <- function(day) floor((as.numeric(day)+3)/7)*7-3
ws_iso <- function(x) {z<-format(as.Date(x,origin="1970-01-01"),"%Y-%m-%d");z[is.na(z)]<-"";z}
ws_n <- function(x)format(x,big.mark=",",scientific=FALSE,trim=TRUE)
ws_sid <- function(x) {
  z<-ws_norm(x);z<-sub("^.*(US[0-9]{3}-?[0-9]{4})$","\\1",z)
  z<-gsub("-","",z,fixed=TRUE)
  if(any(!grepl("^US[0-9]{7}$",z)))stop("Invalid P301 subject ID: ",paste(head(unique(z[!grepl("^US[0-9]{7}$",z)])),collapse=", "))
  z
}
ws_dates <- function(x,field,origin="",strict=FALSE) {
  s<-trimws(ws_text(x));s[s %in% c("NA",".")]<-"";out<-rep(NA_real_,length(s))
  iso<-grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}($|T| )",s)
  v<-substr(s[iso],1,10);d<-suppressWarnings(as.Date(v,"%Y-%m-%d"))
  d[!is.na(d)&format(d,"%Y-%m-%d")!=v]<-as.Date(NA);out[iso]<-as.numeric(d)
  sas<-grepl("^[0-9]{1,2}[A-Za-z]{3}[0-9]{4}$",s)
  if(any(sas)) {
    a<-toupper(s[sas]);dd<-as.integer(sub("[A-Z].*$","",a));mm<-match(sub("^[0-9]+([A-Z]{3}).*$","\\1",a),toupper(month.abb))
    yy<-as.integer(substr(a,nchar(a)-3,nchar(a)));v<-sprintf("%04d-%02d-%02d",yy,mm,dd)
    d<-suppressWarnings(as.Date(v,"%Y-%m-%d"));d[!is.na(d)&format(d,"%Y-%m-%d")!=v]<-as.Date(NA);out[sas]<-as.numeric(d)
  }
  num<-grepl("^-?[0-9]+([.][0-9]+)?$",s)&nzchar(s)
  if(any(num)) {
    if(!nzchar(origin))stop(field,": unformatted numeric dates. Set PHMPT_NUMERIC_DATE_ORIGIN from metadata, or use the original XPT/JMP.")
    v<-as.numeric(s[num]);if(any(v!=floor(v)))stop(field,": fractional numeric dates are unsupported.")
    out[num]<-as.numeric(as.Date(v,origin=origin))
  }
  if(strict&&any(nzchar(s)&!is.finite(out)))stop(field,": invalid/partial dates: ",paste(head(unique(s[nzchar(s)&!is.finite(out)]),4),collapse=", "))
  out
}
ws_specs <- function(domain) switch(domain,
  ADSL=list(required=c("STUDYID","USUBJID","RANDFL","TRT01P","RANDDT","DOSE1FL","DOSE2FL","TR01SDT","DOSE2DT","CUTOFFDT","ELECSBL","RTPCRBL"),
    optional=c("SUBJID","SITEID","SCOV2BL","AP01EDT","UBLDDT","AP02SDT","TR02SDT","EUAVACDT","EOSDT","DTHDT","PDVDT")),
  IS=list(required=c("STUDYID","USUBJID","ISSEQ","ISTESTCD","ISSTRESC","ISDTC"),optional=c("ISTEST","ISORRES","ISSTAT","VISIT","VISITNUM","ISMETHOD")),
  ADIS=list(required=c("STUDYID","USUBJID","ISSEQ","PARAMCD","PARAM","AVALC"),optional=c("ISDTC","ADT","ISSTRESC","VISIT","AVISIT","ABLFL","DTYPE")),
  MB=list(required=c("STUDYID","USUBJID","MBSEQ","MBTESTCD","MBSTRESC","MBDTC"),optional=c("MBORRES","MBSTAT","VISIT","MBMETHOD","MBTEST")),
  ADMB=list(required=c("STUDYID","USUBJID","MBSEQ","PARAMCD","AVALC"),optional=c("MBDTC","ADT","VISIT","AVISIT","ABLFL","DTYPE")),
  stop("Unsupported source domain: ",domain))
ws_read <- function(path,domain,encoding="auto") {
  spec<-ws_specs(domain);ext<-tolower(tools::file_ext(path))
  if(ext=="jmp")x<-ws_read_jmp(path,spec$required,spec$optional)
  else if(ext=="csv") {
    x<-NULL;errors<-character()
    for(enc in if(encoding=="auto")c("UTF-8-BOM","windows-1252") else encoding) {
      z<-tryCatch(withCallingHandlers({
        header<-read.csv(path,nrows=0,check.names=FALSE,fileEncoding=enc)
        upper<-toupper(sub("^\ufeff","",names(header)))
        if(anyDuplicated(upper))stop("Duplicate CSV column names.")
        classes<-ifelse(upper %in% c(spec$required,spec$optional),"character","NULL")
        read.csv(path,colClasses=classes,check.names=FALSE,na.strings=c("","NA"),fileEncoding=enc)
      },warning=function(w)stop(conditionMessage(w))),error=function(e)e)
      if(!inherits(z,"error")){x<-z;break};errors<-c(errors,conditionMessage(z))
    }
    if(is.null(x))stop("CSV read failed: ",paste(errors,collapse="; "))
  } else if(ext %in% c("xpt","sas7bdat")) {
    if(!requireNamespace("haven",quietly=TRUE))stop('Install haven: install.packages("haven")')
    readfun<-if(ext=="xpt")haven::read_xpt else haven::read_sas
    header<-readfun(path,n_max=0,.name_repair="minimal")
    pick<-names(header)[toupper(names(header)) %in% c(spec$required,spec$optional)]
    x<-as.data.frame(readfun(path,col_select=tidyselect::all_of(pick),.name_repair="minimal"))
  } else stop("Unsupported source type: ",ext)
  names(x)<-toupper(sub("^\ufeff","",names(x)))
  if(anyDuplicated(names(x)))stop("Duplicate column names.")
  missing<-setdiff(spec$required,names(x));if(length(missing))stop("Missing ",domain," fields: ",paste(missing,collapse=", "))
  x[]<-lapply(x,ws_text)
  x<-x[ws_norm(x$STUDYID)=="MRNA-1273-P301",,drop=FALSE]
  if(!nrow(x))stop("No P301 rows in this table.")
  for(n in setdiff(spec$optional,names(x)))x[[n]]<-rep("",nrow(x))
  x<-x[c(spec$required,spec$optional)];x$subject<-ws_sid(x$USUBJID)
  for(n in intersect(c("ISSEQ","MBSEQ"),names(x))) {
    ok<-grepl("^[0-9]+([.]0+)?$",x[[n]]);x[[n]][ok]<-as.character(as.numeric(x[[n]][ok]))
  }
  if(domain=="ADSL"&&anyDuplicated(x$subject))stop("Duplicate subject in ADSL.")
  if(domain=="AE") {
    if(any(!nzchar(x$AESEQ)))stop("Missing AESEQ.")
    x<-unique(x)
    if(anyDuplicated(paste(x$subject,x$AESEQ,sep="|")))stop("Conflicting rows share subject/AESEQ.")
  }
  x
}
ws_candidates <- function(root,release,domain) {
  p<-list.files(file.path(root,if(release=="Prod")"prod_data" else "eua_data"),recursive=TRUE,full.names=TRUE,
    pattern="\\.(csv|xpt|sas7bdat|jmp)$",ignore.case=TRUE)
  token<-if(domain=="ADSYMP")"(adsymp|adsy)" else tolower(domain)
  p<-p[grepl(paste0("(^|[_ -])",token,"([_. -]|$)"),basename(p),ignore.case=TRUE)]
  p<-p[!grepl("p101|p201|dmid",basename(p),ignore.case=TRUE)]
  if(release=="Prod") {s3<-grepl("_s3_m[0-9]+_",basename(p),ignore.case=TRUE);if(any(s3))p<-p[s3]}
  unique(p[order(match(tolower(tools::file_ext(p)),c("csv","jmp","xpt","sas7bdat")),p)])
}
ws_submission <- function(p) {
  z<-tolower(basename(p));m<-regexpr("[0-9]+_s[0-9]+_m[0-9]+_mrna-1273-p301",z)
  out<-rep("",length(p));ok<-m>0
  out[ok]<-substring(z[ok],m[ok],m[ok]+attr(m,"match.length")[ok]-1L);out
}
ws_signature <- function(x,origin="") {
  dates<-intersect(c("RANDDT","TR01SDT","DOSE2DT","CUTOFFDT","AP01EDT","UBLDDT","AP02SDT","TR02SDT","EUAVACDT","EOSDT","DTHDT","AESTDTC","AEENDTC","ADT","ASTDT","AENDT"),names(x))
  for(n in dates) {
    d<-ws_dates(x[[n]],n,origin)
    x[[n]]<-ifelse(is.finite(d),ws_iso(d),trimws(x[[n]]))
  }
  x<-unique(x);x<-x[do.call(order,x),,drop=FALSE];rownames(x)<-NULL
  p<-tempfile();on.exit(unlink(p));saveRDS(x,p,compress=FALSE,version=2);unname(tools::md5sum(p))
}
ws_pick <- function(root,release,domain,explicit,out,encoding,origin,min_randomized,submission="") {
  if(nzchar(explicit))paths<-ws_path(root,explicit) else paths<-ws_candidates(root,release,domain)
  if(!length(paths))stop("No ",release," ",domain," found. Set PHMPT_",switch(domain,ADSYMP="SYMPTOMS",domain),"_FILE.")
  logs<-list();valid<-list();seen<-character()
  for(i in seq_along(paths)) {
    p<-paths[i];message("Checking ",release," ",domain,": ",basename(p))
    h<-if(file.exists(p))unname(tools::md5sum(p)) else ""
    other<-!nzchar(explicit)&&nzchar(submission)&&ws_submission(p)!=submission
    duplicate<-nzchar(h)&&h %in% seen
    x<-if(other||duplicate)NULL else tryCatch(ws_read(p,domain,encoding),error=function(e)e)
    cutoff<-NA_real_
    if(!other&&!duplicate&&!inherits(x,"error")&&domain=="ADSL") {
      x<-tryCatch({
        if(sum(ws_norm(x$RANDFL)=="Y")<min_randomized)stop("Incomplete randomized ADSL population.")
        dc<-ws_dates(x$CUTOFFDT,"CUTOFFDT",origin,TRUE)
        if(anyNA(dc)||length(unique(dc))!=1L)stop("ADSL has no unique complete cutoff.")
        cutoff<-dc[1];x
      },error=function(e)e)
    }
    good<-!other&&!duplicate&&!inherits(x,"error")
    logs[[i]]<-data.frame(domain=domain,release=release,path=normalizePath(p,winslash="/",mustWork=FALSE),
      md5=h,rows=if(good)nrow(x) else NA_integer_,cutoff=ws_iso(cutoff),
      status=if(other)"OTHER_SUBMISSION" else if(duplicate)"IDENTICAL_BYTES" else if(good)"VALID" else "REJECTED",
      detail=if(!other&&!duplicate&&!good)conditionMessage(x) else "",stringsAsFactors=FALSE)
    if(good)valid[[length(valid)+1L]]<-list(x=x,index=i,cutoff=cutoff,signature=ws_signature(x,origin))
    if(!other)seen<-c(seen,h)
  }
  log<-do.call(rbind,logs);logfile<-file.path(out,paste0("discovery_",domain,".csv"));ws_csv(log,logfile)
  if(!length(valid))stop("No usable ",domain,". See ",logfile,". ",paste(unique(log$detail[nzchar(log$detail)]),collapse=" | "))
  if(domain=="ADSL") {
    latest<-max(vapply(valid,function(a)a$cutoff,numeric(1)))
    valid<-Filter(function(a)a$cutoff==latest,valid)
  }
  if(length(unique(vapply(valid,function(a)a$signature,character(1))))>1L)
    stop("Different ",domain," snapshots remain. Select one with PHMPT_",switch(domain,ADSYMP="SYMPTOMS",domain),"_FILE; see ",logfile,". Tables are never added together.")
  # Unreadable alternative formats of a readable file are harmless. Other errors
  # may hide a different snapshot and must be resolved explicitly.
  stems<-tolower(tools::file_path_sans_ext(log$path))
  unresolved<-log$status=="REJECTED" & !stems %in% stems[log$status %in% c("VALID","IDENTICAL_BYTES")]
  if(any(unresolved))stop("Unresolved ",domain," candidates; see ",logfile,". Set the explicit input path to choose a reviewed source.")
  chosen<-valid[[1]];log$status[log$status=="VALID"]<-"NOT_SELECTED";log$status[chosen$index]<-"SELECTED"
  ws_csv(log,logfile);list(x=chosen$x,manifest=log[chosen$index,,drop=FALSE])
}


ws_context <- function(x,origin="") {
  for(n in c("RANDFL","DOSE1FL","DOSE2FL")) {
    x[[n]]<-ws_norm(x[[n]]);if(any(!x[[n]] %in% c("Y","N")))stop(n,": missing/unexpected flag.")
  }
  dates<-c("RANDDT","TR01SDT","DOSE2DT","CUTOFFDT","AP01EDT","UBLDDT","AP02SDT","TR02SDT","EUAVACDT","EOSDT","DTHDT")
  for(n in dates)x[[n]]<-ws_dates(x[[n]],n,origin,TRUE)
  for(p in list(c("DOSE1FL","TR01SDT"),c("DOSE2FL","DOSE2DT")))
    if(any((x[[p[1]]]=="Y")!=is.finite(x[[p[2]]])))stop(p[1]," and ",p[2]," disagree.")
  if(any(x$DOSE2FL=="Y"&x$DOSE1FL!="Y")||any(x$DOSE2DT<x$TR01SDT,na.rm=TRUE))stop("Inconsistent dose chronology.")
  x<-x[x$RANDFL=="Y",,drop=FALSE]
  if(any(!is.finite(x$RANDDT))||any(x$RANDDT>x$CUTOFFDT))stop("Invalid randomized start dates.")
  x$arm<-ifelse(ws_norm(x$TRT01P)=="MRNA-1273","mRNA-1273",ifelse(ws_norm(x$TRT01P)=="PLACEBO","Placebo","UNKNOWN"))
  if(any(x$arm=="UNKNOWN"))stop("Unexpected randomized arm.")
  if(any(x$TR01SDT<x$RANDDT|x$TR01SDT>x$CUTOFFDT|x$DOSE2DT>x$CUTOFFDT,na.rm=TRUE))stop("Doses outside randomization/cutoff dates.")
  x$end_full<-x$CUTOFFDT
  for(n in c("EOSDT","DTHDT"))x$end_full<-pmin(x$end_full,ifelse(is.finite(x[[n]]),x[[n]],Inf))
  # Exclude the unblinding/crossover/outside-vaccination day when time is unavailable.
  x$end_blinded<-x$end_full
  x$has_blind_boundary<-is.finite(x$AP01EDT)|is.finite(x$UBLDDT)|is.finite(x$AP02SDT)|is.finite(x$TR02SDT)|is.finite(x$EUAVACDT)
  x$end_blinded<-pmin(x$end_blinded,ifelse(is.finite(x$AP01EDT),x$AP01EDT,Inf))
  for(n in c("UBLDDT","AP02SDT","TR02SDT","EUAVACDT"))x$end_blinded<-pmin(x$end_blinded,ifelse(is.finite(x[[n]]),x[[n]]-1,Inf))
  x
}


ws_text <- function(x) {
  if(inherits(x,"Date")) x<-format(x,"%Y-%m-%d")
  else if(inherits(x,"POSIXt")) x<-format(x,"%Y-%m-%d",tz="UTC")
  x<-as.character(x);x[is.na(x)]<-""
  bad<-!validUTF8(x)
  if(any(bad))x[bad]<-iconv(x[bad],from="windows-1252",to="UTF-8",sub=NA)
  if(anyNA(x))stop("Undecodable input text. Set PHMPT_CSV_ENCODING.")
  enc2utf8(x)
}
ws_csv <- function(x,path) {
  con<-file(path,"wb");on.exit(close(con));writeBin(as.raw(c(239,187,191)),con)
  write.table(x,con,sep=",",row.names=FALSE,quote=TRUE,na="",qmethod="double")
}
ws_path <- function(root,p) if(grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)",p))p else file.path(root,p)
ws_root <- function(root="") {
  if(nzchar(root))return(normalizePath(root,winslash="/",mustWork=TRUE))
  cand<-getwd()
  for(f in sys.frames())if(!is.null(f$ofile))cand<-c(cand,dirname(f$ofile))
  a<-grep("^--file=",commandArgs(),value=TRUE)
  if(length(a))cand<-c(cand,dirname(sub("^--file=","",a[1])))
  if(interactive()&&requireNamespace("rstudioapi",quietly=TRUE)&&rstudioapi::isAvailable()) {
    p<-tryCatch(rstudioapi::getActiveDocumentContext()$path,error=function(e)"")
    if(nzchar(p))cand<-c(cand,dirname(p))
  }
  for(p in unique(cand[dir.exists(cand)]))for(i in seq_len(8L)) {
    if(any(dir.exists(file.path(p,c("prod_data","eua_data")))))return(normalizePath(p,winslash="/",mustWork=TRUE))
    parent<-dirname(p);if(parent==p)break;p<-parent
  }
  stop("Cannot locate eua_data/prod_data. Set PHMPT_ROOT to the project root.")
}

ws_read_jmp <- function(path, fields, optional=character()) {
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
  selected <- unique(c(fields, intersect(optional, column_names)))
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






d57_read <- function(p) {
  x<-read.csv(p,colClasses='character',check.names=FALSE,fileEncoding='UTF-8-BOM',na.strings=c('','NA'))
  x[]<-lapply(x,ws_text);x
}
d57_root <- function(root) {
  if(nzchar(root))return(normalizePath(root,winslash='/',mustWork=TRUE))
  z<-tryCatch(ws_root(''),error=function(e)NULL)
  if(!is.null(z))return(z)
  normalizePath(getwd(),winslash='/',mustWork=TRUE)
}
d57_inputs <- function(root,export,out,cutoff) {
  if(!nzchar(export)&&!length(ws_candidates(root,'Prod','ADSL'))) {
    dirs<-unique(c(root,dirname(list.files(root,pattern='^subject_context[.]csv$',recursive=TRUE,full.names=TRUE))))
    good<-vapply(dirs,function(p)all(file.exists(file.path(p,c('subject_context.csv','subjects.csv','serology_records.csv','pcr_records.csv')))),logical(1))
    candidates<-dirs[good]
    if(!length(candidates))candidates<-list.files(root,pattern='^prod_anti_n_pcr_v.*[.]zip$',recursive=TRUE,full.names=TRUE)
    if(length(candidates)!=1L)stop('Cannot select a single complete export. Set PHMPT_EXPORT to the prior audit folder/ZIP, or supply production tables. Candidates: ',paste(candidates,collapse='; '))
    export<-candidates[1]
  }
  if(nzchar(export)) {
    p<-ws_path(root,export);original<-p
    if(!dir.exists(p)) {
      if(!file.exists(p)||tolower(tools::file_ext(p))!='zip')stop('PHMPT_EXPORT must be an existing folder or ZIP.')
      members<-unzip(p,list=TRUE)$Name
      need<-c('subject_context.csv','subjects.csv','serology_records.csv','pcr_records.csv')
      chosen<-members[basename(members)%in%c(need,'source_manifest.csv')]
      if(anyDuplicated(basename(chosen))||!all(need%in%basename(chosen)))stop('ZIP does not contain one complete validated audit export.')
      if(any(grepl('(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))',chosen)))stop('Unsafe archive member path.')
      p<-file.path(out,'input_export');dir.create(p)
      unzip(original,files=chosen,exdir=p,junkpaths=TRUE)
    }
    ctx<-d57_read(file.path(p,'subject_context.csv'));old<-d57_read(file.path(p,'subjects.csv'))
    req<-c('subject','USUBJID','RANDFL','TRT01P','ELECSBL','RTPCRBL','SCOV2BL','TR01SDT','DOSE2DT','context_conflict')
    if(!all(req%in%names(ctx)))stop('Incomplete subject_context.csv schema.')
    if(anyDuplicated(ctx$subject))stop('Duplicate context subjects.')
    ctx<-ctx[ws_norm(ctx$RANDFL)=='Y',,drop=FALSE]
    if(nrow(ctx)<PHMPT_MIN_RANDOMIZED)stop('Incomplete randomized context export.')
    if(any(nt_bool(ctx$context_conflict)))stop('Resolve randomized context conflicts first.')
    old<-old[old$scope=='STRICT_BLINDED',,drop=FALSE]
    if(anyDuplicated(old$subject)||!all(ctx$subject%in%old$subject))stop('Incomplete strict-window subject export.')
    old<-old[match(ctx$subject,old$subject),,drop=FALSE]
    ctx$context_usable<-nt_bool(old$context_usable)
    ctx$arm<-ifelse(ws_norm(ctx$TRT01P)=='MRNA-1273','mRNA-1273',ifelse(ws_norm(ctx$TRT01P)=='PLACEBO','Placebo','UNKNOWN'))
    for(k in c('TR01SDT','DOSE2DT'))ctx[[k]]<-ws_dates(ctx[[k]],k,strict=TRUE)
    ctx$end_blinded<-pmin(ws_dates(old$scope_end,'scope_end',strict=TRUE),cutoff)
    # Reuse validated strict boundaries from the full previous export, not PDV-inclusive dates.
    ctx$end_full<-cutoff;ctx$CUTOFFDT<-cutoff
    ctx$export_prior_eligible<-nt_bool(old$eligible)
    ctx$export_prior_baseline_conflict<-nt_bool(old$baseline_conflict)
    if(any(ctx$arm=='UNKNOWN'))stop('Invalid randomized arm in export.')
    if('EFFCODT'%in%names(ctx)&&max(ws_dates(ctx$EFFCODT,'EFFCODT',strict=TRUE),na.rm=TRUE)!=cutoff)stop('Export does not match BLA cutoff.')
    sero<-nt_read_export(file.path(p,'serology_records.csv'),'serology',out)
    pcr<-nt_read_export(file.path(p,'pcr_records.csv'),'pcr',out)
    files<-file.path(p,c('subject_context.csv','subjects.csv','serology_records.csv','pcr_records.csv','source_manifest.csv'))
    files<-files[file.exists(files)]
    manifest<-data.frame(path=normalizePath(files,winslash='/'),md5=unname(tools::md5sum(files)),input_mode='COMPLETE_PRIOR_EXPORT')
    if(file.exists(file.path(p,'source_manifest.csv')))file.copy(file.path(p,'source_manifest.csv'),file.path(out,'upstream_source_manifest.csv'))
    writeLines(c('Input: complete prior audit export.',paste('Original location:',original),
      'Strict boundaries and reconciled assay/date flags are inherited from the prior validated export.',
      'No viral-load or complete symptom data are available in this export. No symptom/viral-load adjustment is claimed.'),file.path(out,'input_notes.txt'))
  } else {
    adsl<-ws_pick(root,'Prod','ADSL',PHMPT_ADSL_FILE,out,PHMPT_CSV_ENCODING,PHMPT_NUMERIC_DATE_ORIGIN,PHMPT_MIN_RANDOMIZED)
    ctx<-ws_context(adsl$x,PHMPT_NUMERIC_DATE_ORIGIN)
    if(any(ctx$CUTOFFDT!=cutoff))stop('ADSL cutoff differs from requested BLA snapshot.')
    manifest<-adsl$manifest;sub<-ws_submission(manifest$path)
    get<-function(domain,path,optional=FALSE) {
      z<-if(optional)nt_optional(root,domain,path,out,PHMPT_CSV_ENCODING,PHMPT_NUMERIC_DATE_ORIGIN,PHMPT_MIN_RANDOMIZED,sub) else
        ws_pick(root,'Prod',domain,path,out,PHMPT_CSV_ENCODING,PHMPT_NUMERIC_DATE_ORIGIN,PHMPT_MIN_RANDOMIZED,sub)
      if(!is.null(z))manifest<<-rbind(manifest,z$manifest)
      if(is.null(z))NULL else z$x
    }
    sero<-nt_reconcile(get('IS',PHMPT_IS_FILE),get('ADIS',PHMPT_ADIS_FILE,TRUE),'serology',PHMPT_NUMERIC_DATE_ORIGIN,out)
    pcr<-nt_reconcile(get('MB',PHMPT_MB_FILE),get('ADMB',PHMPT_ADMB_FILE,TRUE),'pcr',PHMPT_NUMERIC_DATE_ORIGIN,out)
    writeLines('Input: production tables. Only laboratory records and randomized subject context are used; no symptom/viral-load adjustment is claimed.',file.path(out,'input_notes.txt'))
  }
  ws_csv(manifest,file.path(out,'source_manifest.csv'))
  list(ctx=ctx,sero=sero,pcr=pcr)
}
d57_baseline <- function(ctx,sero,pcr) {
  ctx$baseline_date<-ctx$TR01SDT
  ctx$baseline_flags_negative<-nt_result(ctx$ELECSBL)=='NEGATIVE'&nt_result(ctx$RTPCRBL)=='NEGATIVE'&nt_result(ctx$SCOV2BL)!='POSITIVE'
  for(kind in c('n','pcr')) {
    r<-if(kind=='n')sero else pcr;j<-match(r$subject,ctx$subject)
    before<-!is.na(j)&r$dated&!r$open_label&!r$pdv&r$date>=ctx$baseline_date[j]-PHMPT_BASELINE_LOOKBACK&r$date<=ctx$baseline_date[j]
    count<-function(flag)tabulate(j[which(before&flag)],nrow(ctx))>0
    ctx[[paste0(kind,'_baseline_negative')]]<-count(r$binary_result&r$negative)
    unplaced<-!is.na(j)&!r$dated&!r$derived&!r$open_label&!r$pdv&
      (r$positive|r$result_conflict)&(r$baseline_flag|r$day1_label|!nzchar(ws_norm(r$raw_visit)))
    ctx[[paste0(kind,'_unplaced_positive_or_conflict')]]<-tabulate(j[which(unplaced)],nrow(ctx))>0
    ctx[[paste0(kind,'_baseline_issue')]]<-count(r$positive|r$result_conflict)|ctx[[paste0(kind,'_unplaced_positive_or_conflict')]]
  }
  if(!'context_usable'%in%names(ctx))ctx$context_usable<-is.finite(ctx$end_blinded)
  ctx$usable_boundary<-ctx$context_usable&is.finite(ctx$end_blinded)&is.finite(ctx$baseline_date)
  ctx$eligible_flags<-ctx$baseline_flags_negative&ctx$usable_boundary&!ctx$n_baseline_issue&!ctx$pcr_baseline_issue
  ctx$eligible<-ctx$eligible_flags&ctx$n_baseline_negative&ctx$pcr_baseline_negative
  ctx$target_dose1<-ctx$baseline_date+56
  ctx$target_dose2<-ifelse(is.finite(ctx$DOSE2DT),ctx$DOSE2DT+28,ctx$target_dose1)
  ctx
}
d57_empty <- function(r) r[FALSE,,drop=FALSE]
d57_min <- function(x) if(length(x))min(x) else NA_real_
d57_max <- function(x) if(length(x))max(x) else NA_real_
d57_select <- function(ctx,sero,pcr,cutoff,mode='PRIMARY',baseline='VERIFIED') {
  si<-split(seq_len(nrow(sero)),sero$subject);pi<-split(seq_len(nrow(pcr)),pcr$subject)
  sim<-match(ctx$subject,names(si));pim<-match(ctx$subject,names(pi))
  subrec<-function(x,k)structure(lapply(x,function(v)v[k]),class='data.frame',row.names=.set_row_names(length(k)))
  target<-if(mode=='NOMINAL_DOSE1')ctx$target_dose1 else ctx$target_dose2
  due<-(if(baseline=='FLAGS')ctx$eligible_flags else ctx$eligible)&is.finite(target)&target+PHMPT_DAY57_LATE<=cutoff
  rows<-vector('list',sum(due));j<-0L;sites<-ws_col(ctx,'SITEID')
  for(i in which(due)) {
    j<-j+1L;s<-subrec(sero,if(is.na(sim[i]))integer() else si[[sim[i]]]);p<-subrec(pcr,if(is.na(pim[i]))integer() else pi[[pim[i]]])
    end<-min(cutoff,ctx$end_blinded[i]);base<-ctx$baseline_date[i];t<-target[i]
    initial<-!s$open_label&!s$pdv&!s$derived
    bounded<-s$dated&s$date>base&s$date<=end
    near<-s$date>=t-PHMPT_DAY57_EARLY&s$date<=t+PHMPT_DAY57_LATE
    # Never relabel a Day 1/29 or illness-labelled specimen as a routine Day 57 visit.
    other_visit<-grepl('DAY 1($| )|DAY 29($| )|ILLNESS',ws_norm(s$raw_visit))&!s$day57_label
    candidate<-initial&bounded&((near&!other_visit)|(mode=='LABELLED_EXTENDED'&s$day57_label))
    candidate[is.na(candidate)]<-FALSE
    pick<-s[candidate,,drop=FALSE];draw<-NA_real_;status<-'NO_ELIGIBLE_SPECIMEN';ids<-refs<-visit<-'';positive<-negative<-FALSE
    if(nrow(pick)) {
      # Selection is fixed before looking at results. Prefer the original labelled visit,
      # then closest date, earlier date on a tie. Keep all records on that date.
      k<-order(!pick$day57_label,abs(pick$date-t),pick$date,pick$record_id)[1]
      draw<-pick$date[k];pick<-pick[pick$date==draw,,drop=FALSE]
      ids<-paste(pick$record_id,collapse=';');refs<-paste(pick$source_refs,collapse=';');visit<-paste(unique(pick$raw_visit),collapse=' | ')
      positive<-any(pick$binary_result&pick$positive);negative<-any(pick$binary_result&pick$negative)
      status<-if(any(pick$result_conflict)||(positive&&negative))'CONFLICT' else if(positive)'POSITIVE' else if(negative)'NEGATIVE' else if(any(pick$tested))'INDETERMINATE' else 'NO_RESULT'
    }
    # When no draw exists, PCR through the administrative visit-window end is exported
    # for missingness review only. It never enters the matched-sample calibration.
    pend<-if(is.finite(draw))draw else min(t+PHMPT_DAY57_LATE,end)
    in_p<-p$dated&p$date>base&p$date<=pend&!p$derived&!p$open_label&!p$pdv
    validpos<-in_p&p$binary_result&p$positive
    validneg<-in_p&p$binary_result&p$negative
    unresolved<-!p$derived&!p$open_label&!p$pdv&((!p$dated&(p$positive|p$result_conflict))|
      (in_p&(p$result_conflict|(p$positive&!p$method_ok))))
    pd<-d57_min(p$date[which(validpos)])
    pc<-if(is.finite(pd))'POSITIVE' else if(any(unresolved,na.rm=TRUE))'UNRESOLVED' else if(any(validneg,na.rm=TRUE))'NEGATIVE_ONLY' else 'NO_COMPLETED_PCR'
    before_n<-s$dated&s$binary_result&s$date>base&s$date<=pend&initial
    prevpos<-s$date[which(before_n&s$positive)]
    lastneg<-s$date[which(s$dated&s$binary_result&s$negative&initial&s$date>=base&s$date<pend)]
    later<-p$dated&p$binary_result&p$positive&p$date>pend&p$date<=min(pend+14,end)&!p$derived&!p$open_label&!p$pdv
    rows[[j]]<-list(subject=ctx$subject[i],arm=ctx$arm[i],SITEID=sites[i],mode=mode,baseline=baseline,
      baseline_date=base,dose2_date=ctx$DOSE2DT[i],target_date=t,scope_end=end,draw_date=draw,
      draw_offset=draw-t,day57_label=any(pick$day57_label),draw_status=status,anti_n_positive=status=='POSITIVE',
      binary=status%in%c('POSITIVE','NEGATIVE'),pcr_status=pc,first_positive_pcr=pd,pcr_lag=draw-pd,
      negative_pcr_records=sum(validneg,na.rm=TRUE),pcr_positive_record_ids=paste(p$record_id[which(validpos)],collapse=';'),
      pcr_unresolved_record_ids=paste(p$record_id[which(unresolved)],collapse=';'),
      later_positive_pcr_14d=any(later,na.rm=TRUE),first_positive_anti_n=d57_min(prevpos),last_negative_anti_n=d57_max(lastneg),
      ever_anti_n_positive=length(prevpos)>0,prior_positive_current_negative=status=='NEGATIVE'&&length(prevpos)>0,
      serology_record_ids=ids,serology_source_refs=refs,raw_visit=visit)
  }
  if(!length(rows))stop('No Day 57-eligible participants at this cutoff.')
  as.data.frame(setNames(lapply(names(rows[[1]]),function(n)unlist(lapply(rows,`[[`,n),use.names=FALSE)),names(rows[[1]])),stringsAsFactors=FALSE)
}
d57_ci <- function(k,n) if(n)unname(binom.test(k,n)$conf.int) else c(NA_real_,NA_real_)
d57_summary <- function(d,scope) {
  do.call(rbind,lapply(c('mRNA-1273','Placebo'),function(arm) {
    z<-d[d$arm==arm,,drop=FALSE];b<-z[z$binary,,drop=FALSE];clear<-b$pcr_status!='UNRESOLVED';npos<-b$pcr_status=='POSITIVE'
    n0<-sum(clear&!npos);k<-sum(b$anti_n_positive&clear&!npos);n<-sum(clear);ci<-d57_ci(k,n)
    data.frame(scope=scope,mode=unique(d$mode),baseline=unique(d$baseline),arm=arm,due=nrow(z),binary=nrow(b),
      binary_coverage_pct=100*nrow(b)/nrow(z),missing_or_uninterpretable=nrow(z)-nrow(b),
      unresolved_pcr=sum(!clear),evaluable=n,pcr_positive=sum(npos),pcr_positive_anti_n_positive=sum(npos&b$anti_n_positive),
      no_positive_pcr=n0,anti_n_positive_no_pcr=k,percent_evaluable=100*k/n,ci_low_pct=100*ci[1],ci_high_pct=100*ci[2],
      percent_all_due_documented=100*k/nrow(z),anti_n_positive_negative_pcr=sum(b$anti_n_positive&b$pcr_status=='NEGATIVE_ONLY'),
      anti_n_positive_no_completed_pcr=sum(b$anti_n_positive&b$pcr_status=='NO_COMPLETED_PCR'),
      later_pcr_14d_among_anti_n_only=sum(b$anti_n_positive&clear&!npos&b$later_positive_pcr_14d),
      prior_positive_current_negative=sum(b$prior_positive_current_negative),
      serology_only_infection_lower_bound_pct=100*k/n,
      serology_only_infection_upper_bound_pct=100*n0/n,stringsAsFactors=FALSE)
  }))
}
d57_calibration <- function(d,scope) {
  lag<-d$pcr_lag
  masks<-list(ALL_LAGS=is.finite(lag)&lag>=0,LAG_0_13=is.finite(lag)&lag>=0&lag<=13,
    LAG_14_27=is.finite(lag)&lag>=14&lag<=27,LAG_28_PLUS=is.finite(lag)&lag>=28,
    LAG_14_PLUS=is.finite(lag)&lag>=14,
    PCR_BEFORE_DOSE2_PLUS14=is.finite(lag)&(!is.finite(d$dose2_date)|d$first_positive_pcr<d$dose2_date+14),
    PCR_AFTER_DOSE2_PLUS14=is.finite(lag)&is.finite(d$dose2_date)&d$first_positive_pcr>=d$dose2_date+14)
  do.call(rbind,lapply(names(masks),function(stratum)do.call(rbind,lapply(c('mRNA-1273','Placebo'),function(arm) {
    ok<-d$arm==arm&d$binary&d$pcr_status=='POSITIVE'&masks[[stratum]];n<-sum(ok);k<-sum(d$anti_n_positive&ok);ci<-d57_ci(k,n)
    data.frame(scope=scope,stratum=stratum,arm=arm,positive=k,negative=n-k,total=n,response=if(n)k/n else NA_real_,
      ci_low=ci[1],ci_high=ci[2],median_lag=if(n)median(lag[ok]) else NA_real_,
      min_lag=if(n)min(lag[ok]) else NA_real_,max_lag=if(n)max(lag[ok]) else NA_real_)
  }))))
}
d57_standardize <- function(cal) {
  a<-cal[cal$stratum%in%c('LAG_0_13','LAG_14_27','LAG_28_PLUS'),]
  out<-list()
  for(sc in unique(a$scope)) {
    x<-a[a$scope==sc,];v<-x[x$arm=='mRNA-1273',];p<-x[x$arm=='Placebo',];p<-p[match(v$stratum,p$stratum),]
    both<-v$total>0&p$total>0;pool<-v$total+p$total;den<-sum(pool[both]);w<-if(den)pool/den else rep(NA_real_,nrow(v))
    for(arm in c('mRNA-1273','Placebo')) {
      r<-if(arm=='mRNA-1273')v$response else p$response
      out[[length(out)+1L]]<-data.frame(scope=sc,arm=arm,standardized_response=if(den)sum(w[both]*r[both]) else NA_real_,
        pooled_common_support=den,pooled_excluded=sum(pool[!both]),common_strata=paste(v$stratum[both],collapse=';'),
        note='Same pooled lag weights in both arms; PCR-detected/evaluable subgroup only. Not the response probability for PCR-undetected infections.')
    }
  };do.call(rbind,out)
}
# Conditional expected-count model. q = P(N+ | infected, no positive PCR, measured).
# f = 1-specificity. r = K/N0 = q*pi + f*(1-pi), so U=(K-f*N0)/(q-f).
# With f>0, double-missed is (1-q)*U, NOT U-K. Model incompatibility is flagged.
d57_model <- function(k,n0,q,specificity=1) {
  f<-1-specificity
  valid<-is.finite(q)&q>f&q<=1&n0>0
  u<-rep(NA_real_,length(q));u[valid]<-(k-f*n0)/(q[valid]-f)
  valid<-valid&is.finite(u)&u>= -1e-8&u<=n0+1e-8
  u[!valid]<-NA_real_
  data.frame(q=q,specificity=specificity,compatible=valid,implied_no_pcr_infections=u,
    expected_true_anti_n_positive=q*u,expected_false_anti_n_positive=f*(n0-u),
    implied_missed_both=(1-q)*u)
}
d57_scenarios <- function(summary,cal) {
  out<-list();trans<-list();parity<-list()
  for(sc in unique(summary$scope)) {
    x<-summary[summary$scope==sc&summary$mode=='PRIMARY'&summary$baseline=='VERIFIED',]
    if(nrow(x)!=2)next
    v<-x[x$arm=='mRNA-1273',];p<-x[x$arm=='Placebo',]
    parity[[length(parity)+1L]]<-data.frame(scope=sc,
      q_v_over_q_p_equal_counts=if(p$anti_n_positive_no_pcr>0)v$anti_n_positive_no_pcr/p$anti_n_positive_no_pcr else NA_real_,
      q_v_over_q_p_equal_rates=if(p$anti_n_positive_no_pcr>0)(v$anti_n_positive_no_pcr/v$evaluable)/(p$anti_n_positive_no_pcr/p$evaluable) else NA_real_,
      specificity=1,note='Threshold for PCR-undetected infection only; observed-serology cohort. Not total infections or efficacy.')
    for(i in 1:2) {
      r<-x[i,]
      for(sp in PHMPT_SPECIFICITY_GRID) {
        z<-d57_model(r$anti_n_positive_no_pcr,r$no_positive_pcr,PHMPT_Q_GRID,sp)
        z$scope<-sc;z$arm<-r$arm;z$evaluable<-r$evaluable;z$observed_anti_n_only<-r$anti_n_positive_no_pcr
        z$implied_no_pcr_percent<-100*z$implied_no_pcr_infections/r$evaluable
        z$implied_missed_both_percent<-100*z$implied_missed_both/r$evaluable
        out[[length(out)+1L]]<-z
      }
      for(st in c('ALL_LAGS','LAG_14_PLUS','LAG_28_PLUS','PCR_AFTER_DOSE2_PLUS14')) {
        b<-cal[cal$scope==sc&cal$arm==r$arm&cal$stratum==st,]
        for(mult in c(.5, .75, 1, 1.25)) {
          q<-b$response*mult
          z<-d57_model(r$anti_n_positive_no_pcr,r$no_positive_pcr,q,1)
          z$scope<-sc;z$arm<-r$arm;z$calibration<-st;z$calibration_positive<-b$positive;z$calibration_total<-b$total
          z$transfer_multiplier<-mult;z$observed_response<-b$response
          z$note<-'Illustrative transfer assumption; not validated sensitivity or an estimate of actual infections.'
          trans[[length(trans)+1L]]<-z
        }
      }
    }
  };list(grid=do.call(rbind,out),transfer=do.call(rbind,trans),parity=do.call(rbind,parity))
}
d57_tests <- function(d,scope) {
  measures<-list(DAY57_BINARY_COVERAGE=list(eligible=rep(TRUE,nrow(d)),yes=d$binary),
    ANTI_N_ONLY_OBSERVED=list(eligible=d$binary&d$pcr_status!='UNRESOLVED',yes=d$anti_n_positive&d$pcr_status%in%c('NEGATIVE_ONLY','NO_COMPLETED_PCR')),
    ANTI_N_RESPONSE_PCR_LAG14=list(eligible=d$binary&d$pcr_status=='POSITIVE'&is.finite(d$pcr_lag)&d$pcr_lag>=14,yes=d$anti_n_positive))
  do.call(rbind,lapply(names(measures),function(nm) {
    a<-measures[[nm]];z<-d[a$eligible,,drop=FALSE];y<-a$yes[a$eligible]
    n<-sapply(c('mRNA-1273','Placebo'),function(ar)sum(z$arm==ar));k<-sapply(c('mRNA-1273','Placebo'),function(ar)sum(y[z$arm==ar]))
    tab<-cbind(k,n-k);pv<-if(all(n>0)&&all(colSums(tab)>0))fisher.test(tab)$p.value else NA_real_
    data.frame(scope=scope,comparison=nm,v_yes=k[1],v_total=n[1],p_yes=k[2],p_total=n[2],fisher_p=pv,
      interpretation=if(nm=='DAY57_BINARY_COVERAGE')'Documentation/selection comparison.' else 'Descriptive post-randomization comparison; does not establish investigator bias.')
  }))
}
d57_write <- function(x,path) {
  for(k in intersect(names(x),c('baseline_date','dose2_date','target_date','scope_end','draw_date','first_positive_pcr','first_positive_anti_n','last_negative_anti_n')))
    x[[k]]<-ws_iso(x[[k]])
  ws_csv(x,path)
}

d57_fmt <- function(k,n,digits=2) ifelse(n>0,sprintf(paste0('%s / %s (%.',digits,'f%%)'),ws_n(k),ws_n(n),100*k/n),'Not estimable')
d57_pct <- function(x,digits=1) ifelse(is.finite(x),sprintf(paste0('%.',digits,'f%%'),100*x),'Not estimable')
d57_methods <- function() c(
  'Population: randomized original-arm participants with negative ADSL baseline PCR/anti-N flags and both negative laboratory results in the 28 days ending on dose 1. Any baseline positive/conflicting result excludes the participant, as does an unplaceable positive/conflict with no date and no usable visit label. Undated open-label results are not treated as baseline records. Participants lacking usable start/strict-boundary context are excluded and counted in population_flow.csv; missing boundaries are not imputed. A flags-only sensitivity cohort is exported.',
  paste0('Primary target: dose 2 + 28 days, otherwise dose 1 + 56. Analyst window: ',PHMPT_DAY57_EARLY,' days early to ',PHMPT_DAY57_LATE,' days late. The entire target window must fall before the cutoff. Participants who discontinued or unblinded before a specimen remain in the coverage denominator.'),
  'Sample selection precedes examination of outcomes: original Day 57 label preferred within the window, then closest collection date, then earlier date. Same-day records are reconciled. Uninterpretable/NOT DONE samples are not replaced by a later favorable result. Open-label Day 57 and decision-visit samples are excluded.',
  'PCR classification ends on the selected blood-collection date, not the BLA cutoff. Positive-PCR records after that draw are separately flagged for the following 14 days. No recorded positive PCR is split into negative-results-only versus no completed PCR. Unresolved positive/date/result records form a separate exclusion category.',
  'Baseline-to-Day 57 anti-N positivity is cumulative evidence of infection; it does not date infection. Day 29 positives are not recast as new infections after dose 2. Prior positivity with a negative Day 57 result is separately exported. No reinfection counting or symptomatic-case adjudication is attempted.',
  'Calibration uses the same Day 57 specimen in participants with prior positive PCR. PCR-to-blood lag strata and timing relative to dose 2 are reported. At least 14 days is an analyst response-lag restriction, not a guarantee of seroconversion. First positive PCR is not the infection date.',
  'Exact binomial 95% intervals describe observed response proportions. Fisher tests are descriptive, two-sided, and Holm-adjusted across the exported comparisons. Conditioning on infection/PCR detection is post-randomization selection; these contrasts do not identify a causal vaccine or investigator effect.',
  'Common-lag standardization uses pooled stratum weights only where both arms have PCR-positive measured participants. It equalizes that measured lag distribution, not infection timing, viral burden, symptom severity, or the unknown distribution among PCR-undetected infections.',
  'Sensitivity parameter q is the unobserved probability of Day 57 anti-N positivity among infected, PCR-undetected participants with measured serology. This is not the empirically observed response s among PCR-positive participants. Transfer scenarios q = multiplier * s are explicitly unvalidated assumptions.',
  'For N0 participants without a recorded positive PCR and K anti-N positives, f=1-specificity and U=(K-f*N0)/(q-f). The additional infections missed by both methods are (1-q)*U. Incompatible expected-count scenarios are marked, not truncated into feasible estimates. The q/specificity grids are illustrative, not confidence intervals or endorsed plausible ranges.',
  'Under specificity 1, U=K/q and equal arm-specific U/evaluable rates require q_v/q_p=(K_v/evaluable_v)/(K_p/evaluable_p). This is a threshold for PCR-undetected infection in the observed cohort, not equality of all infections, and not efficacy.',
  'Missing or uninterpretable serology is never treated as negative. Coverage and worst-case bounds on Day 57 positivity are reported; no missing-at-random assumption is silently imposed. Sampling intervals and test dependence remain limits to inference.',
  'The production data are retrospectively restricted to EUA (25 November 2020) and BLA (26 March 2021) collection dates. These are not separate original database snapshots. The earlier 30/56 and 69/226 whole-blinded-period results and 46.9%/92.1% response fractions are not inputs to this Day 57 calculation.',
  'Serology-only cases were detected through serology. This analysis does not establish their absence from derived infection datasets or a published endpoint. No alternative efficacy, number of symptomatic COVID-19 cases, or proven investigator bias is calculated.',
  'Source for the biological concern: Follmann et al., Ann Intern Med 2022, doi:10.7326/M22-1300, https://pmc.ncbi.nlm.nih.gov/articles/PMC9258784/ . The study associates anti-N response with vaccination and viral burden; its PDV response proportions are not Day 57 calibration values.'
)
d57_html_escape <- function(x) {
  x<-gsub('&','&amp;',as.character(x),fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE);gsub('>','&gt;',x,fixed=TRUE)
}
d57_html_table <- function(x) paste0('<table><thead><tr>',paste0('<th>',d57_html_escape(names(x)),'</th>',collapse=''),'</tr></thead><tbody>',
  paste(vapply(seq_len(nrow(x)),function(i)paste0('<tr>',paste0('<td>',d57_html_escape(unlist(x[i,])),'</td>',collapse=''),'</tr>'),character(1)),collapse=''),'</tbody></table>')
d57_display <- function(summary,scope='BLA') {
  z<-summary[summary$scope==scope&summary$mode=='PRIMARY'&summary$baseline=='VERIFIED',];z<-z[match(c('mRNA-1273','Placebo'),z$arm),]
  rows<-c('Baseline-negative, Day 57 window due','Interpretable Day 57 anti-N result',
    'Missing/uninterpretable Day 57 result','Unresolved PCR among measured subjects','Evaluable paired classification',
    'Anti-N positive; no positive PCR by draw','  With completed negative PCR only','  With no completed follow-up PCR',
    'PCR positive by draw; anti-N positive','PCR positive by draw; anti-N negative',
    'Later positive PCR within 14 days (anti-N-only group)','Prior anti-N positive, Day 57 negative')
  col<-function(r)c(ws_n(r$due),d57_fmt(r$binary,r$due),ws_n(r$missing_or_uninterpretable),ws_n(r$unresolved_pcr),
    ws_n(r$evaluable),d57_fmt(r$anti_n_positive_no_pcr,r$evaluable),ws_n(r$anti_n_positive_negative_pcr),ws_n(r$anti_n_positive_no_completed_pcr),
    ws_n(r$pcr_positive_anti_n_positive),ws_n(r$pcr_positive-r$pcr_positive_anti_n_positive),ws_n(r$later_pcr_14d_among_anti_n_only),ws_n(r$prior_positive_current_negative))
  data.frame(Finding=rows,'mRNA-1273'=col(z[1,]),Placebo=col(z[2,]),check.names=FALSE)
}
d57_cal_display <- function(cal,scope='BLA') {
  x<-cal[cal$scope==scope,];strata<-unique(x$stratum)
  labels<-c(ALL_LAGS='All PCR-to-draw lags',LAG_0_13='PCR 0-13 days before draw',LAG_14_27='PCR 14-27 days before draw',LAG_28_PLUS='PCR >=28 days before draw',LAG_14_PLUS='PCR >=14 days before draw',
    PCR_BEFORE_DOSE2_PLUS14='PCR before dose 2 +14 / no dose 2',PCR_AFTER_DOSE2_PLUS14='PCR on/after dose 2 +14')
  get<-function(ar)vapply(strata,function(st){r<-x[x$arm==ar&x$stratum==st,];paste0(d57_fmt(r$positive,r$total,1),'\n95% CI ',d57_pct(r$ci_low),' to ',d57_pct(r$ci_high))},character(1))
  data.frame(Stratum=unname(labels[strata]),'mRNA-1273'=get('mRNA-1273'),Placebo=get('Placebo'),check.names=FALSE)
}
d57_window_display <- function(summary) {
  keys<-unique(summary[c('scope','mode','baseline')]);rows<-lapply(seq_len(nrow(keys)),function(i){
    k<-keys[i,];a<-summary[summary$scope==k$scope&summary$mode==k$mode&summary$baseline==k$baseline,]
    lab<-if(k$baseline=='FLAGS')'Primary window, baseline flags' else switch(k$mode,PRIMARY='Primary window, laboratory verified',NOMINAL_DOSE1='Dose 1 +56 target',LABELLED_EXTENDED='Original-labelled visit, extended dates')
    v<-a[a$arm=='mRNA-1273',];p<-a[a$arm=='Placebo',]
    data.frame(Definition=paste(k$scope,lab,sep=': '),'mRNA-1273'=d57_fmt(v$anti_n_positive_no_pcr,v$evaluable),Placebo=d57_fmt(p$anti_n_positive_no_pcr,p$evaluable),check.names=FALSE)
  });do.call(rbind,rows)
}
d57_pdf_page <- function(title,subtitle,page) {
  grid::grid.newpage();grid::grid.rect(gp=grid::gpar(fill='white',col=NA))
  grid::grid.text(title,x=.055,y=.95,just='left',gp=grid::gpar(fontsize=19,fontface='bold',col='#17324d'))
  grid::grid.text(subtitle,x=.055,y=.905,just=c('left','top'),gp=grid::gpar(fontsize=10,col='#4c5b68'))
  grid::grid.lines(x=c(.055,.945),y=c(.045,.045),gp=grid::gpar(col='#ccd5dd'))
  grid::grid.text('Moderna P301 | Day 57 ascertainment audit | v1.0.0',x=.055,y=.025,just='left',gp=grid::gpar(fontsize=8,col='#596875'))
  grid::grid.text(page,x=.945,y=.025,just='right',gp=grid::gpar(fontsize=8,col='#596875'))
}
d57_pdf_table <- function(x,top=.84,row_height=.043,widths=c(.5,.195,.195),font=10) {
  xs<-.055+c(0,cumsum(widths))[seq_along(widths)]
  for(i in 0:nrow(x)) {
    y<-top-i*row_height
    grid::grid.rect(x=.5,y=y-row_height/2,width=.89,height=row_height,gp=grid::gpar(fill=if(i==0)'#17324d' else if(i%%2)'#eef3f7' else 'white',col=NA))
    values<-if(i==0)names(x) else as.character(unlist(x[i,]))
    for(j in seq_along(values))grid::grid.text(values[j],x=xs[j]+.009,y=y-row_height/2,just='left',gp=grid::gpar(fontsize=font,col=if(i==0)'white' else '#1d2d3b',fontface=if(i==0)'bold' else 'plain',lineheight=1.15))
  }
  invisible(top-(nrow(x)+1)*row_height)
}
d57_pdf_text <- function(txt,y,width=114,size=10) {
  t<-paste(strwrap(txt,width=width),collapse='\n')
  grid::grid.text(t,x=.055,y=y,just=c('left','top'),gp=grid::gpar(fontsize=size,col='#263848',lineheight=1.3))
}
d57_report <- function(summary,cal,std,scenario,tests,missing,out) {
  z<-summary[summary$scope=='BLA'&summary$mode=='PRIMARY'&summary$baseline=='VERIFIED',];v<-z[z$arm=='mRNA-1273',];p<-z[z$arm=='Placebo',]
  b<-cal[cal$scope=='BLA'&cal$stratum=='LAG_14_PLUS',];bv<-b[b$arm=='mRNA-1273',];bp<-b[b$arm=='Placebo',]
  par<-scenario$parity[scenario$parity$scope=='BLA',]
  audit<-sprintf(paste0('Among baseline-PCR/anti-N-negative participants with interpretable Day 57 serology and a classifiable PCR history, ',
    '%s of %s mRNA-1273 recipients (%.3f%%) and %s of %s placebo recipients (%.3f%%) had anti-N positivity without a recorded positive PCR through the blood draw. ',
    'These participants were identified by serology, not missed by both methods. In participants whose first positive PCR preceded that same draw by at least 14 days, anti-N positivity was %s versus %s. ',
    'Transferring those response probabilities to PCR-undetected infections is an unverified assumption because infection timing, viral burden and ascertainment differ. ',
    'Under perfect anti-N specificity, equal rates of PCR-undetected infection in the measured cohort would require the vaccine/placebo anti-N detection-probability ratio to be %.3f. ',
    'The analysis therefore reports observable discordance and assumption-dependent scenarios; it does not establish the number missed by both methods, equal total infections, or an alternative efficacy.'),
    ws_n(v$anti_n_positive_no_pcr),ws_n(v$evaluable),v$percent_evaluable,ws_n(p$anti_n_positive_no_pcr),ws_n(p$evaluable),p$percent_evaluable,
    d57_fmt(bv$positive,bv$total,1),d57_fmt(bp$positive,bp$total,1),par$q_v_over_q_p_equal_rates)
  writeLines(audit,file.path(out,'audit_paragraph.txt'))
  methods<-d57_methods();writeLines(paste0(seq_along(methods),'. ',methods),file.path(out,'METHODS.txt'))
  obs<-d57_display(summary);cd<-d57_cal_display(cal)
  tr<-scenario$transfer;tr<-tr[tr$scope=='BLA'&tr$transfer_multiplier==1,]
  td<-do.call(rbind,lapply(unique(tr$calibration),function(st){a<-tr[tr$calibration==st&tr$arm=='mRNA-1273',];b<-tr[tr$calibration==st&tr$arm=='Placebo',]
    fmt<-function(r)if(r$compatible)sprintf('%.1f / %.1f',r$implied_no_pcr_infections,r$implied_missed_both) else 'Not estimable'
    data.frame(Calibration=switch(st,ALL_LAGS='All PCR-to-draw lags',LAG_14_PLUS='PCR >=14 days before draw',LAG_28_PLUS='PCR >=28 days before draw',PCR_AFTER_DOSE2_PLUS14='PCR on/after dose 2 +14'), 'mRNA-1273'=fmt(a),Placebo=fmt(b),check.names=FALSE)}))
  names(td)[2:3]<-c('mRNA-1273: U / both','Placebo: U / both')
  html<-paste0('<!doctype html><html><meta charset="utf-8"><title>Day 57 anti-N ascertainment audit</title><style>',
    'body{font:15px/1.5 Arial,sans-serif;max-width:1120px;margin:40px auto;color:#1d2d3b;padding:0 24px}h1,h2{color:#17324d}table{border-collapse:collapse;width:100%;margin:18px 0}th{background:#17324d;color:white;text-align:left}th,td{padding:9px;border-bottom:1px solid #d7dfe6;white-space:pre-line}tr:nth-child(even){background:#eef3f7}.note{background:#edf4f8;padding:18px}li{margin:10px 0}@media print{body{font-size:10pt;margin:0}h2{break-after:avoid}tr{break-inside:avoid}}</style>',
    '<h1>Day 57: PCR and anti-N ascertainment</h1><p class="note">Observed findings and sensitivity analysis. No alternative efficacy is calculated.</p>',
    '<p>',d57_html_escape(audit),'</p><h2>BLA restriction: observed comparison</h2>',d57_html_table(obs),
    '<h2>EUA calendar restriction</h2>',d57_html_table(d57_display(summary,'EUA')),
    '<h2>Visit-definition and baseline sensitivity</h2>',d57_html_table(d57_window_display(summary)),
    '<h2>Response in PCR-positive subjects, using the same Day 57 sample</h2>',d57_html_table(cd),
    '<h2>Transfer assumption: q equals the observed PCR-positive response</h2><p>U = implied infections without positive PCR; both = additional infections missed by both methods. These are expected-count scenarios, not measured people or validated estimates. Specificity = 100%.</p>',d57_html_table(td),
    '<h2>Equal-rate threshold</h2>',d57_html_table(scenario$parity),
    '<h2>Common-lag descriptive standardization</h2>',d57_html_table(std),
    '<h2>Missingness bounds</h2>',d57_html_table(missing),
    '<h2>Descriptive tests</h2>',d57_html_table(tests),
    '<h2>Methods and remaining limitations</h2><ol>',paste0('<li>',d57_html_escape(methods),'</li>',collapse=''),'</ol></html>')
  writeLines(html,file.path(out,'day57_anti_n_report.html'),useBytes=TRUE)
  pdf<-file.path(out,'day57_anti_n_audit.pdf')
  if(capabilities('cairo'))grDevices::cairo_pdf(pdf,width=11.69,height=8.27,family='sans',onefile=TRUE) else grDevices::pdf(pdf,width=11.69,height=8.27,onefile=TRUE)
  on.exit(grDevices::dev.off())
  d57_pdf_page('Day 57: the observed comparison','Original randomized arms | Baseline PCR/anti-N negative | Strict blinded samples | BLA cutoff',1)
  y<-d57_pdf_table(obs,row_height=.044,font=10)
  d57_pdf_text('The anti-N-only group has no recorded positive PCR by its blood draw. It was detected by serology. Percentages use the evaluable paired cohort; the coverage row uses all participants whose Day 57 window was due.',y-.035)
  d57_pdf_page('Calibration uses the same blood draw','Anti-N positivity among PCR-positive participants | Exact binomial 95% confidence intervals',2)
  y<-d57_pdf_table(cd,row_height=.077,font=10)
  d57_pdf_text('Lag strata expose differences in time available for antibody response. Early infections may precede full vaccination. These measured response rates cannot automatically be transferred to infections that escaped PCR.',y-.035)
  d57_pdf_page('What the adjustment assumes','Illustrative q = observed response among PCR-positive participants | Specificity assumed 100%',3)
  y<-d57_pdf_table(td,row_height=.065,font=10)
  d57_pdf_text('Each cell is U / additional missed by both. U = observed anti-N-only positives divided by assumed detection probability q. These figures apply only to the measured, classifiable Day 57 cohort. They do not include PCR-detected infections.',y-.04)
  d57_pdf_text(sprintf('Equal PCR-undetected counts require q_v / q_p = %.3f. Equal rates per evaluable participant require q_v / q_p = %.3f. This threshold is an assumption boundary, not evidence that parity occurred.',par$q_v_over_q_p_equal_counts,par$q_v_over_q_p_equal_rates),y-.17)
  d57_pdf_text('Missing serology, unknown infection dates, viral burden and correlated PCR/antibody detection remain unresolved. The supplied CSVs vary q independently in each arm and include false-positive scenarios. Incompatible expected-count assumptions are explicitly flagged.',y-.30)
  d57_pdf_page('Sensitivity: where equal rates would occur','PCR-undetected infections per evaluable participant | Perfect specificity | No transfer assumption imposed',4)
  grid::pushViewport(grid::viewport(x=.50,y=.50,width=.72,height=.55))
  q<-PHMPT_Q_GRID;z<-outer(q,q,function(qv,qp)log2((v$anti_n_positive_no_pcr/v$evaluable/qv)/(p$anti_n_positive_no_pcr/p$evaluable/qp)))
  feasible_v<-d57_model(v$anti_n_positive_no_pcr,v$no_positive_pcr,q,1)$compatible
  feasible_p<-d57_model(p$anti_n_positive_no_pcr,p$no_positive_pcr,q,1)$compatible
  z[!outer(feasible_v,feasible_p,`&`)]<-NA_real_
  cols<-colorRampPalette(c('#22618c','#ffffff','#aa453d'))(101)
  breaks<-seq(-4,4,length.out=102)
  for(i in seq_along(q))for(j in seq_along(q)) {
    col<-if(!is.finite(z[i,j]))'#dddddd' else cols[max(1,min(101,findInterval(z[i,j],breaks)))]
    grid::grid.rect(x=(i-.5)/length(q),y=(j-.5)/length(q),width=1/length(q),height=1/length(q),gp=grid::gpar(fill=col,col=NA))
  }
  # Values are located at cell centers. Clip parity line to grid boundaries.
  qp<-seq(min(q),max(q),length.out=500);qv<-par$q_v_over_q_p_equal_rates*qp;ok<-qv>=min(q)&qv<=max(q)
  rescale<-function(a)(a-min(q))/(max(q)-min(q))*(1-1/length(q))+.5/length(q)
  if(any(ok))grid::grid.lines(x=rescale(qv[ok]),y=rescale(qp[ok]),gp=grid::gpar(col='#111111',lwd=2,lty=2))
  ticks<-seq(.2,1,.2);grid::grid.xaxis(at=rescale(ticks),label=paste0(100*ticks,'%'),gp=grid::gpar(fontsize=10));grid::grid.yaxis(at=rescale(ticks),label=paste0(100*ticks,'%'),gp=grid::gpar(fontsize=10))
  grid::grid.text('Assumed anti-N detection q: mRNA-1273',x=.5,y=-.12,gp=grid::gpar(fontsize=11))
  grid::grid.text('Assumed anti-N detection q: placebo',x=-.14,y=.5,rot=90,gp=grid::gpar(fontsize=11))
  grid::popViewport()
  d57_pdf_text('Blue: lower implied PCR-undetected infection rate in mRNA-1273; red: higher. Dashed: equality. Colors reflect assumptions, not significance. Grey: incompatible values. See CSV for counts.',.115,size=9)
  d57_pdf_page('Scope and interpretation','What is addressed, and what remains unidentified',5)
  lines<-c('Addressed: consistent baseline, original randomized arms, one original-series Day 57 specimen, matched PCR/serology endpoints, strict blinded boundaries, missingness, sample timing and duplicate/conflicting records.',
    'Exposed rather than assumed away: anti-N detection among PCR-undetected infections. Its probability cannot be measured from infections that neither test detected. PCR-positive response is an imperfect calibration group, even after lag standardization.',
    'Incomplete detection is not equivalent to biased efficacy. This audit concerns infection ascertainment through Day 57. It neither reconstructs symptomatic endpoint eligibility nor identifies the reason that a participant was or was not PCR-tested.',
    'All serology-only positives are already detected through an alternative test. Whether a given published endpoint included them is a separate tracing task. Negative anti-N results cannot by themselves count or rule out prior infections.',
    'The 30/56 strict and 69/226 boundary-inclusive counts quoted in the previous audit use longer observation windows. The corresponding later-serology response rates are not inputs to this Day 57 analysis.',
    'Biological context: Anti-N response varies with vaccination and viral burden. The participant decision-visit response is not a Day 57 sensitivity estimate.')
  for(i in seq_along(lines))d57_pdf_text(lines[i],.84-(i-1)*.12,width=112,size=10)
  d57_pdf_page('How the visit definition changes the result','Anti-N positive without prior positive PCR / evaluable participants | Each row has its own matched draw endpoint',6)
  yy<-d57_pdf_table(d57_window_display(summary),row_height=.061,font=9.5)
  d57_pdf_text('Extended dates admit the original labelled Day 57 visit beyond the primary analyst window. Open-label and decision-visit specimens remain excluded. A changed count is a change of ascertainment window, not evidence that infection incidence changed.',yy-.035,size=10)
  d57_pdf_text('The baseline-flags sensitivity still excludes baseline contradictions and unplaceable positive/conflict records; it relaxes the requirement for two independently documented negative baseline laboratory results. The EUA restriction is applied to production collection dates, not an original EUA database snapshot.',yy-.15,size=10)
  invisible(pdf)
}

analyse_day57_anti_n_sensitivity <- function(root=PHMPT_ROOT,export=PHMPT_EXPORT,output_dir=PHMPT_OUTPUT_DIR,show_report=PHMPT_SHOW_REPORT) {
  for(v in c(PHMPT_DAY57_EARLY,PHMPT_DAY57_LATE,PHMPT_BASELINE_LOOKBACK))if(length(v)!=1||!is.finite(v)||v<0||v!=floor(v))stop('Invalid day-window setting.')
  if(length(PHMPT_Q_GRID)<2||any(!is.finite(PHMPT_Q_GRID)|PHMPT_Q_GRID<=0|PHMPT_Q_GRID>1)||is.unsorted(PHMPT_Q_GRID,strictly=TRUE))stop('Q grid must be increasing probabilities in (0,1].')
  if(length(PHMPT_Q_GRID)>2&&max(abs(diff(PHMPT_Q_GRID)-diff(PHMPT_Q_GRID)[1]))>1e-8)stop('Use an equally spaced Q grid for the printed probability surface.')
  if(any(!is.finite(PHMPT_SPECIFICITY_GRID)|PHMPT_SPECIFICITY_GRID<=0|PHMPT_SPECIFICITY_GRID>1))stop('Invalid specificity grid.')
  root<-d57_root(root);out<-if(nzchar(output_dir))ws_path(root,output_dir) else file.path(root,'analysis',paste0('day57_anti_n_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_',Sys.getpid()))
  if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop('Output folder is not empty: ',out)
  dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out,winslash='/',mustWork=TRUE)
  complete<-FALSE;writeLines('RUNNING',file.path(out,'RUN_STATUS.txt'))
  on.exit(if(!complete)writeLines('FAILED - outputs incomplete; inspect R console and audit CSVs.',file.path(out,'RUN_STATUS.txt')),add=TRUE)
  cuts<-c(EUA=ws_dates(PHMPT_EUA_CUTOFF,'EUA cutoff',strict=TRUE),BLA=ws_dates(PHMPT_BLA_CUTOFF,'BLA cutoff',strict=TRUE))
  if(any(!is.finite(cuts))||cuts['EUA']>cuts['BLA'])stop('Invalid cutoffs.')
  message('Day 57 anti-N ascertainment audit v1.0.0 | ',out)
  inp<-d57_inputs(root,export,out,cuts['BLA']);ctx<-d57_baseline(inp$ctx,inp$sero,inp$pcr)
  d57_write(ctx,file.path(out,'baseline_roster.csv'))
  if('export_prior_eligible'%in%names(ctx))d57_write(ctx[ctx$eligible!=ctx$export_prior_eligible,,drop=FALSE],file.path(out,'baseline_changes_from_prior_audit.csv'))
  flow<-do.call(rbind,lapply(c('mRNA-1273','Placebo'),function(ar){x<-ctx[ctx$arm==ar,];data.frame(arm=ar,randomized=nrow(x),negative_adsl_flags=sum(x$baseline_flags_negative),
    baseline_negative_unusable_boundary=sum(x$baseline_flags_negative&!x$usable_boundary),negative_flags_without_lab_conflict=sum(x$eligible_flags),laboratory_verified=sum(x$eligible),no_dose2_among_verified=sum(x$eligible&!is.finite(x$DOSE2DT)))}))
  ws_csv(flow,file.path(out,'population_flow.csv'))
  summaries<-cal<-tests<-missing<-list();primary<-list()
  for(sc in names(cuts)) {
    for(mode in c('PRIMARY','NOMINAL_DOSE1','LABELLED_EXTENDED')) {
      message(sc,': ',mode)
      d<-d57_select(ctx,inp$sero,inp$pcr,cuts[sc],mode)
      d57_write(d,file.path(out,paste0('subjects_',sc,'_',mode,'.csv')))
      summaries[[length(summaries)+1L]]<-d57_summary(d,sc)
      if(mode=='PRIMARY') {
        primary[[sc]]<-d;cal[[sc]]<-d57_calibration(d,sc);tests[[sc]]<-d57_tests(d,sc)
        missing[[sc]]<-do.call(rbind,lapply(c('mRNA-1273','Placebo'),function(ar){z<-d[d$arm==ar,];pos<-sum(z$anti_n_positive);miss<-sum(!z$binary)
          data.frame(scope=sc,arm=ar,due=nrow(z),positive_day57=pos,unmeasured_or_uninterpretable=miss,
            positivity_lower_pct=100*pos/nrow(z),positivity_upper_pct=100*(pos+miss)/nrow(z),
            pcr_positive_without_binary_serology=sum(!z$binary&z$pcr_status=='POSITIVE'),
            note='Worst-case positivity bounds, not confidence intervals or infection bounds. Missing samples are not assumed negative.')}))
      }
    }
    d<-d57_select(ctx,inp$sero,inp$pcr,cuts[sc],'PRIMARY','FLAGS')
    d57_write(d,file.path(out,paste0('subjects_',sc,'_FLAGS_SENSITIVITY.csv')))
    summaries[[length(summaries)+1L]]<-d57_summary(d,sc)
  }
  summary<-do.call(rbind,summaries);cal<-do.call(rbind,cal);tests<-do.call(rbind,tests);missing<-do.call(rbind,missing)
  tests$fisher_p_holm<-p.adjust(tests$fisher_p,'holm');std<-d57_standardize(cal);scenario<-d57_scenarios(summary,cal)
  ws_csv(summary,file.path(out,'observed_summary.csv'));ws_csv(cal,file.path(out,'pcr_positive_day57_response.csv'))
  ws_csv(std,file.path(out,'lag_standardized_response.csv'));ws_csv(tests,file.path(out,'descriptive_tests.csv'));ws_csv(missing,file.path(out,'missing_serology_bounds.csv'))
  ws_csv(scenario$grid,file.path(out,'independent_q_scenarios.csv'));ws_csv(scenario$transfer,file.path(out,'calibration_transfer_scenarios.csv'));ws_csv(scenario$parity,file.path(out,'parity_thresholds.csv'))
  # Explicit two-arm Cartesian scenario surface; no vaccine efficacy column.
  grids<-list()
  for(sc in names(cuts))for(sp in PHMPT_SPECIFICITY_GRID) {
    g<-scenario$grid;v<-g[g$scope==sc&g$specificity==sp&g$arm=='mRNA-1273',];p<-g[g$scope==sc&g$specificity==sp&g$arm=='Placebo',]
    a<-expand.grid(v=seq_len(nrow(v)),p=seq_len(nrow(p)));vv<-v[a$v,];pp<-p[a$p,]
    grids[[length(grids)+1L]]<-data.frame(scope=sc,specificity=sp,q_v=vv$q,q_p=pp$q,compatible=vv$compatible&pp$compatible,
      implied_v_no_pcr=vv$implied_no_pcr_infections,implied_p_no_pcr=pp$implied_no_pcr_infections,
      implied_v_missed_both=vv$implied_missed_both,implied_p_missed_both=pp$implied_missed_both,
      difference_no_pcr_per_10000=10000*(vv$implied_no_pcr_infections/vv$evaluable-pp$implied_no_pcr_infections/pp$evaluable),
      ratio_no_pcr_rates=(vv$implied_no_pcr_infections/vv$evaluable)/(pp$implied_no_pcr_infections/pp$evaluable))
  }
  ws_csv(do.call(rbind,grids),file.path(out,'two_arm_sensitivity_surface.csv'))
  params<-list(version='1.0.0',EUA=PHMPT_EUA_CUTOFF,BLA=PHMPT_BLA_CUTOFF,early=PHMPT_DAY57_EARLY,late=PHMPT_DAY57_LATE,
    baseline_lookback=PHMPT_BASELINE_LOOKBACK,q_grid=paste(PHMPT_Q_GRID,collapse=';'),specificity_grid=paste(PHMPT_SPECIFICITY_GRID,collapse=';'),
    built_utc=format(Sys.time(),tz='UTC',usetz=TRUE),baseline='ADSL negative + both laboratory negative + no baseline conflict',target='dose2+28; fallback dose1+56')
  ws_csv(data.frame(parameter=names(params),value=as.character(unlist(params))),file.path(out,'parameters.csv'))
  capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
  stopifnot(all(summary$evaluable+summary$unresolved_pcr==summary$binary),all(summary$binary+summary$missing_or_uninterpretable==summary$due),
    all(summary$anti_n_positive_no_pcr<=summary$no_positive_pcr),all(summary$no_positive_pcr+summary$pcr_positive==summary$evaluable))
  d57_report(summary,cal,std,scenario,tests,missing,out)
  d57_export_audit_assets(root,out,summary,cal,std,scenario)
  writeLines('COMPLETE',file.path(out,'RUN_STATUS.txt'));complete<-TRUE
  message('Complete: ',file.path(out,'day57_anti_n_report.html'))
  print(d57_display(summary),row.names=FALSE,right=FALSE)
  if(isTRUE(show_report))try(utils::browseURL(file.path(out,'day57_anti_n_report.html')),silent=TRUE)
  invisible(list(output_dir=out,summary=summary,calibration=cal,primary=primary,scenarios=scenario))
}
# ---- RUN: standalone, no command-line arguments needed ----
# Audit PNG export extension v1.0.0 | 2026-10-23
# Standalone; stable Vxx filenames in <project>/audit_assets plus this run's copy.
# Long tables continue as _02.png, _03.png, etc. References remain in HTML/provenance.
PHMPT_AUDIT_ASSET_DIR <- getOption("phmpt.audit_asset_dir", "audit_assets")
PHMPT_AUDIT_DPI <- getOption("phmpt.audit_dpi", 300L)

# ---- Standalone audit-asset renderer (V08-V35), version 1.0.0 ----
# Embedded in each script: no shared R file or browser is required.
ax_clean <- function(x) {
  x<-as.character(x);x[is.na(x)]<-"Unavailable"
  x<-gsub("https?://[^[:space:]<>]+","",x,perl=TRUE)
  x<-gsub("(?i)doi:[[:space:]]*[^[:space:]]+","",x,perl=TRUE)
  x<-gsub("\\[[0-9]+([,; -]+[0-9]+)*\\]","",x,perl=TRUE)
  trimws(x)
}
ax_path <- function(root,p)if(grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)",p))p else file.path(root,p)
ax_root <- function(root="",out=getwd()) {
  if(nzchar(root))return(normalizePath(root,winslash="/",mustWork=TRUE))
  for(start in unique(c(getwd(),out))) {
    p<-normalizePath(start,winslash="/",mustWork=TRUE)
    repeat {
      if(any(dir.exists(file.path(p,c("prod_data","eua_data")))))return(p)
      if(basename(p)=="analysis")return(dirname(p))
      up<-dirname(p);if(up==p)break;p<-up
    }
  }
  stop("Set PHMPT_ROOT or phmpt.root to the project folder for audit_assets exports.")
}
ax_read <- function(out,file,required=TRUE) {
  p<-file.path(out,file)
  if(!file.exists(p)){if(required)stop("Missing export: ",p);return(data.frame())}
  read.csv(p,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
}
ax_csv <- function(x,path) {
  con<-file(path,"wb");on.exit(close(con));writeBin(as.raw(c(239,187,191)),con)
  write.table(x,con,sep=",",row.names=FALSE,quote=TRUE,na="",qmethod="double")
}
ax_num <- function(x)ifelse(is.na(x),"Unavailable",format(x,big.mark=",",scientific=FALSE,trim=TRUE))
ax_rate <- function(k,n,digits=1)ifelse(is.na(k)|is.na(n),"Unavailable",ifelse(n>0,
  paste0(ax_num(k),"/",ax_num(n),sprintf(paste0(" (%.",digits,"f%%)"),100*k/n)),"0/0 (not estimable)"))
ax_p <- function(p) {
  z<-ifelse(is.finite(p),formatC(p,digits=4,format="f"),"Not estimable")
  for(i in 2:6)z[is.finite(p)&p<10^-i]<-paste0("< ",format(10^-i,scientific=FALSE,trim=TRUE))
  z
}
ax_text <- function(s,x,y,size=11.5,bold=FALSE,col="#25343D",just=c("left","top"))
  grid::grid.text(ax_clean(s),x=grid::unit(x,"in"),y=grid::unit(1,"npc")-grid::unit(y,"in"),just=just,
    gp=grid::gpar(fontfamily="sans",fontsize=size,fontface=if(bold)"bold" else "plain",col=col,lineheight=1.12))
ax_height <- function(s,size=11.5)grid::convertHeight(grid::grobHeight(grid::textGrob(s,
  gp=grid::gpar(fontfamily="sans",fontsize=size,lineheight=1.12))),"in",valueOnly=TRUE)
ax_wrap <- function(s,width,size=11.5,bold=FALSE) {
  measure<-function(t)grid::convertWidth(grid::grobWidth(grid::textGrob(t,
    gp=grid::gpar(fontfamily="sans",fontsize=size,fontface=if(bold)"bold" else "plain"))),"in",valueOnly=TRUE)
  paragraphs<-strsplit(ax_clean(s),"\n",fixed=TRUE)[[1]];all<-character()
  for(p in paragraphs) {
    words<-strsplit(gsub("[[:space:]]+"," ",trimws(p))," ",fixed=TRUE)[[1]];lines<-character();line<-""
    for(w in words) {
      z<-paste(c(if(nzchar(line))line,w),collapse=" ")
      if(measure(z)<=width){line<-z;next}
      if(nzchar(line)){lines<-c(lines,line);line<-""}
      while(nchar(w)>1L&&measure(w)>width) {
        n<-nchar(w)-1L;while(n>1L&&measure(substr(w,1,n))>width)n<-n-1L
        lines<-c(lines,substr(w,1,n));w<-substring(w,n+1L)
      };line<-w
    };all<-c(all,lines,line)
  };paste(all,collapse="\n")
}
ax_begin <- function(root,out,owner,dpi=PHMPT_AUDIT_DPI,asset_dir=PHMPT_AUDIT_ASSET_DIR) {
  if(length(dpi)!=1||!is.finite(dpi)||dpi<150||dpi>600)stop("Audit DPI must be 150-600.")
  if(length(asset_dir)!=1||is.na(asset_dir)||!nzchar(asset_dir))stop("Audit asset directory is empty.")
  a<-new.env(parent=emptyenv());a$root<-ax_root(root,out);a$out<-normalizePath(out,winslash="/",mustWork=TRUE)
  a$owner<-owner;a$dpi<-dpi;a$dest<-ax_path(a$root,asset_dir)
  a$stage<-tempfile(paste0("audit_",owner,"_"),tmpdir=out);dir.create(a$stage)
  a$rows<-list();a$files<-character();a
}
ax_png <- function(a,id,slug,draw,height=7,page=1L,pages=1L,status="AVAILABLE",rows=NA_integer_) {
  file<-paste0(id,"_",slug,if(page>1)sprintf("_%02d",page) else "",".png");path<-file.path(a$stage,file)
  grDevices::png(path,width=8.5,height=height,units="in",res=a$dpi,pointsize=13,bg="white",
    type=if(capabilities("cairo"))"cairo" else if(.Platform$OS.type=="windows")"windows" else getOption("bitmapType","Xlib"))
  device<-grDevices::dev.cur();ok<-FALSE
  on.exit({if(device%in%grDevices::dev.list())grDevices::dev.off(device);if(!ok)unlink(path)})
  draw();grDevices::dev.off(device)
  bytes<-readBin(path,"raw",n=file.info(path)$size)
  if(length(bytes)<32||!identical(bytes[1:8],as.raw(c(137,80,78,71,13,10,26,10)))||
    !identical(tail(bytes,12),as.raw(c(0,0,0,0,73,69,78,68,174,66,96,130))))stop("Incomplete PNG: ",path)
  dims<-vapply(list(bytes[17:20],bytes[21:24]),function(b)sum(as.integer(b)*256^(3:0)),numeric(1))
  if(any(abs(dims-round(c(8.5,height)*a$dpi))>1))stop("Invalid PNG dimensions.")
  a$rows[[length(a$rows)+1L]]<-data.frame(visual=id,file=file,page=page,pages=pages,status=status,rows=rows,width_px=dims[1],height_px=dims[2],dpi=a$dpi)
  a$files<-c(a$files,file);ok<-TRUE;invisible(path)
}
ax_data <- function(a,id,x,suffix="data") {
  file<-paste0(id,"_",suffix,".csv");ax_csv(x,file.path(a$stage,file));a$files<-c(a$files,file)
}
ax_table <- function(a,id,slug,title,data,subtitle="",note="",widths=NULL,status="AVAILABLE",body_size=11.5,highlight=integer()) {
  data<-as.data.frame(data,check.names=FALSE);ax_data(a,id,data)
  if(!ncol(data))data<-data.frame(Result="No matching rows",check.names=FALSE)
  if(!nrow(data)){data[1,]<-rep("No matching rows",ncol(data));status<-"NO_MATCHING_ROWS"}
  if(is.null(widths))widths<-if(ncol(data)==3)c(.5,.25,.25) else rep(1/ncol(data),ncol(data))
  widths<-widths/sum(widths)*7.8;left<-.35+c(0,head(cumsum(widths),-1))
  tmp<-tempfile(fileext=".pdf");if(capabilities("cairo"))grDevices::cairo_pdf(tmp,width=8.5,height=13,family="sans") else grDevices::pdf(tmp,width=8.5,height=13)
  dev<-grDevices::dev.cur();on.exit({if(dev%in%grDevices::dev.list())grDevices::dev.off(dev);unlink(tmp)})
  wraprow<-function(z,size=body_size,bold=FALSE)vapply(seq_along(z),function(j)ax_wrap(z[j],widths[j]-.16,size,bold),character(1))
  heads<-wraprow(names(data),11.2,TRUE);head_h<-max(vapply(heads,ax_height,numeric(1),size=11.2))+.22
  title<-ax_wrap(paste(id,title),7.8,17,TRUE);sub<-ax_wrap(subtitle,7.8,11.2);foot<-ax_wrap(note,7.8,10.6)
  top<-.25+ax_height(title,17)+.18+if(nzchar(sub))ax_height(sub,11.2)+.18 else 0
  bottom<-.42+if(nzchar(foot))ax_height(foot,10.6) else 0
  budget<-12.8-top-head_h-bottom
  cells<-list();heights<-numeric();origin<-integer()
  # Split a very long source row into labelled continuation rows instead of truncating it.
  for(i in seq_len(nrow(data))) {
    line<-lapply(wraprow(as.character(data[i,])),strsplit,split="\n",fixed=TRUE);line<-lapply(line,`[[`,1)
    cap<-max(1L,floor((min(budget,5.5)-.25)/(body_size/72*1.4)))
    nchunks<-max(1L,ceiling(max(lengths(line))/cap))
    for(k in seq_len(nchunks)) {
      row<-vapply(line,function(z)paste(z[seq_along(z)>=(k-1)*cap+1L & seq_along(z)<=k*cap],collapse="\n"),character(1))
      if(k>1&&!nzchar(row[1]))row[1]<-paste0(as.character(data[i,1])," (continued)")
      row[1]<-ax_wrap(row[1],widths[1]-.16,body_size)
      h<-max(vapply(row,ax_height,numeric(1),size=body_size))+.24
      if(h>budget)stop("Audit row is too tall even after wrapping: ",id)
      cells[[length(cells)+1L]]<-row;heights<-c(heights,h);origin<-c(origin,i)
    }
  }
  pages<-list();page<-integer();used<-0
  for(i in seq_along(cells)) {
    if(length(page)&&used+heights[i]>budget){pages[[length(pages)+1L]]<-page;page<-integer();used<-0}
    page<-c(page,i);used<-used+heights[i]
  };if(length(page))pages[[length(pages)+1L]]<-page
  grDevices::dev.off(dev)
  for(p in seq_along(pages)) {
    ix<-pages[[p]];height<-max(3.6,top+head_h+sum(heights[ix])+bottom)
    draw<-function() {
      grid::grid.newpage();ax_text(title,.35,.25,17,TRUE,col="#17384B")
      sy<-.25+ax_height(title,17)+.18;if(nzchar(sub))ax_text(sub,.35,sy,11.2)
      y<-top
      for(r in c(0L,ix)) {
        h<-if(r==0)head_h else heights[r]
        grid::grid.rect(x=grid::unit(.35,"in"),y=grid::unit(1,"npc")-grid::unit(y,"in"),
          width=grid::unit(7.8,"in"),height=grid::unit(h,"in"),just=c("left","top"),
          gp=grid::gpar(fill=if(r==0)"#E8EFF3" else if(origin[r]%in%highlight)"#FFF0D9" else if(origin[r]%%2==0)"#F3F6F8" else "white",col=NA))
        values<-if(r==0)heads else cells[[r]]
        for(j in seq_along(values))ax_text(values[j],left[j]+.08,y+.11,if(r==0)11.2 else body_size,r==0||(r>0&&j==1&&origin[r]%in%highlight),
          col=if(r!=0&&grepl("mRNA|Vaccine",names(data)[j],ignore.case=TRUE))"#1E667D" else if(r!=0&&grepl("Placebo",names(data)[j]))"#A65B22" else "#25343D")
        y<-y+h
      }
      if(nzchar(foot))ax_text(foot,.35,y+.16,10.6,col="#52616C")
      ax_text(paste0("Page ",p," of ",length(pages)),8.15,height-.17,9.5,col="#64737D",just=c("right","bottom"))
    }
    ax_png(a,id,slug,draw,height,p,length(pages),status,length(unique(origin[ix])))
  };invisible(NULL)
}
ax_unavailable <- function(a,id,slug,title,reason)
  ax_table(a,id,slug,title,data.frame(Status="Not available in this run",Reason=reason,check.names=FALSE),widths=c(.3,.7),status="UNAVAILABLE")
ax_finish <- function(a,notes=character()) {
  if(!length(a$rows))stop("No audit assets were rendered.")
  manifest<-do.call(rbind,a$rows);manifest$renderer_version<-"1.0.0";manifest$generated_at<-format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z")
  manifest$run_directory<-a$out;manifest$md5<-unname(tools::md5sum(file.path(a$stage,manifest$file)))
  mf<-paste0(a$owner,"_manifest.csv");ax_csv(manifest,file.path(a$stage,mf))
  readme<-paste0(a$owner,"_README.txt")
  writeLines(c("AUDIT PNG ASSETS", "8.5 inches wide; 300 dpi by default (2550 pixels). White background. Upload PNG originals to Substack.",
    "Tables wrap and continue as _02, _03, etc.; no rows are silently omitted. Printed figures have no bibliography, footnote markers, DOI or source URLs.",
    "Statistical definitions, endpoint labels, observation windows and denominators remain on figures. Source references stay in the main report and original provenance/method exports.",
    "A status of UNAVAILABLE is an explicit missing-input/model notice, not a zero finding. Data CSVs preserve the displayed inputs.",notes),file.path(a$stage,readme))
  a$files<-unique(c(a$files,mf,readme));destinations<-unique(c(file.path(a$out,"audit_assets"),a$dest))
  # Complete rendering precedes publishing. Only files listed in this owner's
  # previous manifest and this batch may be replaced or removed.
  for(dest in destinations) {
    dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    oldmf<-file.path(dest,mf);old<-if(file.exists(oldmf))ax_read(dest,mf)$file else character()
    owned_old<-old[grepl("^V[0-9]{2}_[A-Za-z0-9_]+[.]png$",old)]
    existing<-intersect(unique(c(a$files,owned_old)),list.files(dest))
    backup<-tempfile("backup_",tmpdir=a$stage);dir.create(backup)
    if(length(existing)&&!all(file.copy(file.path(dest,existing),backup,overwrite=TRUE)))stop("Cannot back up previous audit assets.")
    success<-FALSE
    tryCatch({
      for(f in a$files) {
        if(!file.copy(file.path(a$stage,f),file.path(dest,f),overwrite=TRUE))stop("Cannot publish ",f)
        if(unname(tools::md5sum(file.path(a$stage,f)))!=unname(tools::md5sum(file.path(dest,f))))stop("Audit copy mismatch: ",f)
      }
      stale<-setdiff(owned_old,a$files)
      if(length(stale)&&unlink(file.path(dest,stale))!=0L)stop("Cannot remove obsolete continuation pages.")
      success<-TRUE
    },finally={if(!success){unlink(file.path(dest,setdiff(a$files,existing)));file.copy(file.path(backup,existing),dest,overwrite=TRUE)}})
    unlink(backup,recursive=TRUE)
  }
  unlink(a$stage,recursive=TRUE)
  message("Audit PNGs: ",a$dest," | ",paste(unique(manifest$visual),collapse=", "))
  invisible(list(directory=a$dest,manifest=manifest))
}

d57_export_audit_assets <- function(root,out,summary,cal,std,scenario) {
  a<-ax_begin(root,out,"V25_V27_day57")
  subjects<-ax_read(out,"subjects_BLA_PRIMARY.csv")
  arms<-c("mRNA-1273","Placebo")
  prior<-vapply(arms,function(ar)sum(subjects$arm==ar & subjects$binary &
    subjects$pcr_status%in%c("NEGATIVE_ONLY","NO_COMPLETED_PCR") & subjects$prior_positive_current_negative,na.rm=TRUE),numeric(1))
  z<-summary[summary$scope=="BLA"&summary$mode=="PRIMARY"&summary$baseline=="VERIFIED",]
  z<-z[match(arms,z$arm),];if(anyNA(z$arm))stop("Day 57 BLA primary paired cohort unavailable.")
  d<-d57_display(summary)
  d<-rbind(d,data.frame(Finding=c("Earlier positive, Day 57 negative, no positive PCR","Ever anti-N positive through Day 57, no positive PCR"),
    `mRNA-1273`=ax_num(c(prior[1],z$anti_n_positive_no_pcr[1]+prior[1])),Placebo=ax_num(c(prior[2],z$anti_n_positive_no_pcr[2]+prior[2])),check.names=FALSE))
  ax_table(a,"V25","day57_observed_results","Day 57 anti-N results and PCR history",d,
    "Baseline PCR/anti-N negative | Strict blinded samples | BLA cutoff | Original randomized arms",
    "PCR history ends at the selected blood draw. Coverage uses participants whose Day 57 window was due; anti-N-only percentages use the evaluable paired cohort. Earlier positive serology is retained separately. The generic prior-positive row may include participants with positive PCR; the final two rows explicitly exclude them.")
  c<-d57_cal_display(cal);s<-std[std$scope=="BLA",];s<-s[match(arms,s$arm),]
  c<-rbind(c,data.frame(Stratum="Standardized to pooled PCR-to-draw interval weights",`mRNA-1273`=d57_pct(s$standardized_response[1]),Placebo=d57_pct(s$standardized_response[2]),check.names=FALSE))
  ax_table(a,"V26","day57_pcr_positive_calibration","Day 57 response among PCR-positive participants",c,
    "Same Day 57 blood draw | BLA restriction | Exact binomial 95% confidence intervals",
    "Interval standardization uses common pooled weights for 0-13, 14-27 and >=28 days. These are selected PCR-detected infections with available serology; the fractions do not identify detection probabilities among infections that escaped PCR.")
  t<-scenario$transfer;t<-t[t$scope=="BLA"&t$calibration=="LAG_14_PLUS"&t$transfer_multiplier==1,];t<-t[match(arms,t$arm),]
  if(nrow(t)!=2||anyNA(t$arm))stop("Day 57 >=14-day transfer scenario unavailable.")
  K<-z$anti_n_positive_no_pcr;U<-t$implied_no_pcr_infections;remaining<-U-K-prior
  fmt<-function(x)ifelse(is.finite(x),sprintf("%.1f",x),"Not estimable")
  values<-list("Observed Day 57 anti-N positive; no positive PCR"=ax_num(K),
    "Assumed detection q (transferred PCR-positive fraction)"=paste0(ax_rate(t$calibration_positive,t$calibration_total)," = ",sprintf("%.3f",t$q)),
    "Scenario-implied PCR-undetected infections: K / q"=fmt(U),
    "Less Day 57 positives: K / q - K"=fmt(U-K),
    "Earlier positive serology, Day 57 negative"=ax_num(prior),
    "Remaining after both current and earlier serology"=ifelse(is.finite(remaining)&remaining>=0,fmt(remaining),"Incompatible with known positives"),
    "Minimum infections already identified through serology"=ax_num(K+prior))
  d<-data.frame(Calculation=names(values),`mRNA-1273`=vapply(values,`[`,character(1),1),Placebo=vapply(values,`[`,character(1),2),check.names=FALSE)
  ax_table(a,"V27","day57_sensitivity_scenario","Illustrative anti-N detection scenario",d,
    "Perfect specificity assumed | Transfer of >=14-day PCR-positive response fractions | Measured Day 57 cohort",
    "Fractional values are scenario expectations, not identified participants or an alternative efficacy. The scenario must accommodate every earlier positive record. PCR-undetected infections are not total infections; PCR-detected cases are separate.")
  ax_data(a,"V27",data.frame(arm=arms,K=K,q=t$q,implied=U,earlier_positive=prior,remaining=remaining,compatible_with_prior=is.finite(remaining)&remaining>=0),"numeric_scenario")
  ax_finish(a)
}

if(isTRUE(getOption('phmpt.autorun',TRUE)))analyse_day57_anti_n_sensitivity()
