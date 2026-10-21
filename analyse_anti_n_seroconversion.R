# =============================================================================
# Moderna P301 - anti-N seroconversion after confirmed COVID-19 | v1.0.0
# 2026-10-21 | Stand-alone R >=4.3: select all and run in RStudio.
# Source: Follmann et al. Ann Intern Med. 2022. DOI: 10.7326/M22-1300.
# Always recalculates the published 21/52 versus 605/648 comparison.
# Also attempts a PHMPT participant-level reconstruction using full Prod tables.
# An unavailable reconstruction or viral-load model is explicitly reported; the
# published benchmark is NEVER relabelled as a result calculated from local data.
# Base R for CSV/JMP and graphics; haven is needed only for XPT/SAS7BDAT.
# No downloads, additional scripts, or installation are triggered automatically.
# Load functions only: options(phmpt.autorun=FALSE) before source().
# =============================================================================
PHMPT_ROOT <- getOption('phmpt.root','')
PHMPT_OUTPUT_DIR <- getOption('phmpt.output_dir','')
PHMPT_ADSL_FILE <- ''
PHMPT_ADEFF_FILE <- ''
PHMPT_IS_FILE <- ''
PHMPT_ADIS_FILE <- ''
PHMPT_MB_FILE <- ''
PHMPT_ADMB_FILE <- ''
PHMPT_CASE_PARAMCD <- ''   # blank: unique adjudicated, non-severe primary COVID label
PHMPT_VIRAL_PARAMCD <- ''  # blank: unique explicitly log10 copies/mL parameter in ADMB
# Optional reviewed viral-load asset: subject,PCRDT,LOG10_COPIES_ML,SOURCE_REF.
# Exact diagnosis-date match; Ct is never treated as log10 viral copies.
PHMPT_VIRAL_FILE <- ''
PHMPT_BLA_CUTOFF <- '2021-03-26'
PHMPT_MIN_RANDOMIZED <- 30000L
PHMPT_CSV_ENCODING <- 'auto'
PHMPT_NUMERIC_DATE_ORIGIN <- ''
PHMPT_PNG_DPI <- 300L
PHMPT_RECONSTRUCT <- TRUE
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
  ADSL=list(required=c("STUDYID","USUBJID","FASFL","RANDFL","TRT01P","RANDDT","DOSE1FL","DOSE2FL","TR01SDT","DOSE2DT","CUTOFFDT","ELECSBL","RTPCRBL"),
    optional=c("SUBJID","SITEID","SCOV2BL","AP01EDT","UBLDDT","AP02SDT","TR02SDT","EUAVACDT","EOSDT","DTHDT","PDVDT")),
  IS=list(required=c("STUDYID","USUBJID","ISSEQ","ISTESTCD","ISSTRESC","ISDTC"),optional=c("ISTEST","ISORRES","ISSTAT","VISIT","VISITNUM","ISMETHOD")),
  ADIS=list(required=c("STUDYID","USUBJID","ISSEQ","PARAMCD","PARAM","AVALC"),optional=c("ISDTC","ADT","ISSTRESC","VISIT","AVISIT","ABLFL","DTYPE")),
  MB=list(required=c("STUDYID","USUBJID","MBSEQ","MBTESTCD","MBSTRESC","MBDTC"),optional=c("MBORRES","MBSTAT","VISIT","MBMETHOD","MBTEST")),
  ADMB=list(required=c("STUDYID","USUBJID","MBSEQ","PARAMCD","AVALC"),optional=c("MBDTC","ADT","VISIT","AVISIT","ABLFL","DTYPE","PARAM","AVAL","AVALU","MBSTRESU","MBORRESU","MBLLOQ")),
  ADEFF=list(required=c("STUDYID","USUBJID","PARAMCD","PARAM","ADT","PCRDT","ANL01FL"),optional=c("DTYPE","AVALC","AVAL","CNSR")),
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






# ---- Published benchmark and transparent, recalculated uncertainty ----
sc_p <- function(p) {
  if(!is.finite(p))return('Not estimable')
  for(t in c(1e-6,1e-5,1e-4,1e-3,1e-2))if(p<t)return(paste('<',formatC(t,format='f',digits=round(-log10(t)))))
  sprintf('%.4f',p)
}
sc_pct <- function(k,n)if(n>0)sprintf('%d/%d (%.1f%%)',k,n,100*k/n) else '0/0 (not estimable)'
sc_ci <- function(k,n)if(n>0)as.numeric(stats::binom.test(k,n)$conf.int) else c(NA_real_,NA_real_)
sc_stats <- function(k,n,label) {
  ci<-t(mapply(sc_ci,k,n));d<-data.frame(source=label,arm=c('mRNA-1273','Placebo'),positive=k,negative=n-k,evaluable=n,
    proportion=ifelse(n>0,k/n,NA_real_),ci_lower=ci[,1],ci_upper=ci[,2],ci_method='Clopper-Pearson exact binomial')
  mat<-matrix(c(k[1],n[1]-k[1],k[2],n[2]-k[2]),nrow=2,byrow=TRUE,
    dimnames=list(c('mRNA-1273','Placebo'),c('Anti-N positive','Anti-N negative')))
  fisher<-if(all(n>0)&&all(colSums(mat)>0))stats::fisher.test(mat) else NULL
  contrast<-data.frame(source=label,difference_pp=100*(k[1]/n[1]-k[2]/n[2]),
    positivity_ratio=(k[1]/n[1])/(k[2]/n[2]),fisher_p=if(is.null(fisher))NA_real_ else fisher$p.value)
  list(arms=d,contrast=contrast,table=mat)
}
sc_benchmark <- function(out) {
  d<-data.frame(arm=c('mRNA-1273','Placebo'),confirmed_cases=c(58L,754L),no_PDV=c(0L,51L),PDV_without_serology=c(6L,55L),
    evaluable=c(52L,648L),anti_N_positive=c(21L,605L),abstract_ci_lower=c(.27,.92),abstract_ci_upper=c(.54,.95),
    table1_ci_lower=c(.271,.915),table1_ci_upper=c(.537,.953))
  stopifnot(all(d$confirmed_cases-d$no_PDV-d$PDV_without_serology==d$evaluable))
  z<-sc_stats(d$anti_N_positive,d$evaluable,'Recalculated from published counts')
  ws_csv(d,file.path(out,'published_reference_counts.csv'))
  ws_csv(z$arms,file.path(out,'published_counts_recalculated.csv'))
  ws_csv(z$contrast,file.path(out,'published_counts_contrasts.csv'))
  ws_csv(data.frame(predictor='Illness-visit viral level, per 1 log10 copies/mL',OR=1.90,ci_lower=1.59,ci_upper=2.28,
    adjustment='Randomized arm',status='Published reference only; cannot be refitted from aggregate counts'),file.path(out,'published_viral_load_reference.csv'))
  z$reference<-d;z
}
sc_inventory <- function(x,path) {
  cols<-intersect(c('PARAMCD','PARAM','AVALU','MBSTRESU','MBORRESU'),names(x))
  z<-unique(x[cols]);z<-z[order(z$PARAMCD),,drop=FALSE];ws_csv(z,path)
}
sc_case_param <- function(e,explicit,out) {
  sc_inventory(e,file.path(out,'adeff_parameter_inventory.csv'))
  if(nzchar(explicit)) {
    if(!explicit%in%e$PARAMCD)stop('Selected PHMPT_CASE_PARAMCD is absent from ADEFF.')
    return(explicit)
  }
  label<-ws_norm(e$PARAM)
  ok<-grepl('^COVID[- ]?19( DISEASE| CASE)? (PER ADJUDICATION|AS ADJUDICATED)$',label)
  codes<-unique(e$PARAMCD[ok])
  if(length(codes)!=1L)stop('Select PHMPT_CASE_PARAMCD from adeff_parameter_inventory.csv. A unique adjudicated primary COVID endpoint was not identifiable; no endpoint was guessed.')
  codes
}
# Each participant's recorded PDV sample is retained, including the boundary date.
# No arbitrary "nearest day 53" sample or later post-crossover sample is substituted.
sc_pdv <- function(ctx,s,cutoff,out) {
  label<-ws_norm(s$raw_visit)
  exact<-grepl('^(PARTICIPANT DECISION( VISIT)?|OL-D1|OPEN LABEL DAY 1)$',label)
  raw<-s[exact&s$dated&s$date<=cutoff,,drop=FALSE]
  bys<-split(raw,raw$subject);q<-ctx[c('subject','arm')]
  q$pdv_date<-NA_real_;q$positive<-NA_integer_;q$pdv_status<-'NO_PDV_SEROLOGY_RECORD';q$record_ids<-''
  declared<-ws_dates(ws_col(ctx,'PDVDT'),'PDVDT',strict=TRUE)
  for(i in seq_len(nrow(q))) {
    z<-bys[[q$subject[i]]]
    if(is.null(z)||!nrow(z))next
    dates<-unique(z$date)
    if(length(dates)!=1L){q$pdv_status[i]<-'MULTIPLE_PDV_DATES_REVIEW';next}
    q$pdv_date[i]<-dates[1];q$record_ids[i]<-paste(z$record_id,collapse='; ')
    if(is.finite(declared[i])&&declared[i]!=dates[1]){q$pdv_status[i]<-'ADSL_PDV_DATE_CONFLICT';next}
    if(is.finite(ctx$AP02SDT[i])&&dates[1]>ctx$AP02SDT[i]){q$pdv_status[i]<-'SAMPLE_AFTER_OPEN_LABEL_START';next}
    if(any(z$result_conflict)){q$pdv_status[i]<-'PDV_RESULT_CONFLICT';next}
    ans<-unique(as.integer(z$positive[z$binary_result]))
    if(length(ans)==1L){q$positive[i]<-ans;q$pdv_status[i]<-'EVALUABLE'}
    else q$pdv_status[i]<-'NO_INTERPRETABLE_PDV_RESULT'
  }
  q
}
sc_cohort <- function(ctx,e,s,p,case_code,cutoff,origin,out) {
  q<-sc_pdv(ctx,s,cutoff,out)
  e$.ADEFF_row<-seq_len(nrow(e));z<-e[e$PARAMCD==case_code,,drop=FALSE]
  z$illness_date<-ws_dates(z$ADT,'ADEFF.ADT',origin)
  z$pcr_date<-ws_dates(z$PCRDT,'ADEFF.PCRDT',origin)
  z$selected<-ws_norm(z$ANL01FL)=='Y' & !nzchar(z$DTYPE)
  # A censoring row is not an observed event, if this variable is supplied.
  z$selected<-z$selected & (!nzchar(z$CNSR)|z$CNSR=='0')
  z<-z[order(z$subject,z$illness_date,z$pcr_date,z$.ADEFF_row),,drop=FALSE]
  ix<-match(z$subject,ctx$subject);j<-match(z$subject,q$subject)
  z$arm<-ctx$arm[ix];z$reason<-''
  reason<-function(flag,txt){k<-which(flag&!nzchar(z$reason));z$reason[k]<<-txt}
  reason(is.na(ix),'NOT_ON_RANDOMIZED_ROSTER')
  reason(!z$selected,'NOT_SELECTED_EVENT_RECORD')
  reason(!is.finite(z$illness_date)|!is.finite(z$pcr_date),'MISSING_OR_PARTIAL_CASE_DATES')
  reason(ws_norm(ctx$FASFL[ix])!='Y','OUTSIDE_FULL_ANALYSIS_SET')
  reason(nt_result(ctx$ELECSBL[ix])!='NEGATIVE'|nt_result(ctx$RTPCRBL[ix])!='NEGATIVE','NOT_BASELINE_DOUBLE_NEGATIVE')
  reason(nt_result(ctx$SCOV2BL[ix])=='POSITIVE','BASELINE_FLAGS_CONFLICT')
  reason(!is.finite(ctx$DOSE2DT[ix]),'NO_SECOND_DOSE_DATE')
  reason(z$illness_date<ctx$DOSE2DT[ix]+14,'ILLNESS_BEFORE_DOSE2_PLUS14')
  reason(z$illness_date>cutoff|z$pcr_date>cutoff,'AFTER_BLA_CUTOFF')
  reason(z$illness_date>ctx$end_blinded[ix],'ILLNESS_AFTER_BLINDED_FOLLOWUP')
  reason(is.finite(q$pdv_date[j])&z$illness_date>=q$pdv_date[j],'ILLNESS_NOT_BEFORE_PDV')
  reason(is.finite(q$pdv_date[j])&z$pcr_date>q$pdv_date[j],'DIAGNOSTIC_PCR_AFTER_PDV')
  # Verify the positive result on the ADEFF PCR date. Absence is kept as a gap.
  pk<-paste(p$subject,p$date,sep='|');positive_keys<-pk[p$binary_result&p$positive]
  z$pcr_verified<-paste(z$subject,z$pcr_date,sep='|')%in%positive_keys
  reason(!z$pcr_verified,'NO_VERIFIED_POSITIVE_PCR_ON_ADEFF_DATE')
  accepted<-which(!nzchar(z$reason));dup<-duplicated(z$subject[accepted])
  z$reason[accepted[dup]]<-'LATER_OR_DUPLICATE_CASE_RECORD'
  z$included_case<-!nzchar(z$reason)
  ws_csv(transform(z,illness_date=ws_iso(illness_date),pcr_date=ws_iso(pcr_date)),file.path(out,'case_selection_audit.csv'))
  z<-z[z$included_case,,drop=FALSE]
  if(!nrow(z))stop('No eligible confirmed cases remain; inspect case_selection_audit.csv and the selected ADEFF parameter.')
  j<-match(z$subject,q$subject)
  z$pdv_date<-q$pdv_date[j];z$anti_n_positive<-q$positive[j];z$pdv_status<-q$pdv_status[j];z$serology_record_ids<-q$record_ids[j]
  z$days_illness_to_pdv<-z$pdv_date-z$illness_date
  z$evaluable<-z$pdv_status=='EVALUABLE'&is.finite(z$days_illness_to_pdv)&z$days_illness_to_pdv>0
  # Prior positive tests are exported for adjudication review; no extra ad hoc
  # exclusion is invented beyond the selected primary endpoint and baseline rules.
  previous<-function(records) {
    v<-split(records[records$binary_result&records$positive,,drop=FALSE],records$subject[records$binary_result&records$positive])
    vapply(seq_len(nrow(z)),function(i){w<-v[[z$subject[i]]];!is.null(w)&&any(w$date<z$illness_date[i])},logical(1))
  }
  z$earlier_positive_anti_n_for_review<-previous(s);z$earlier_positive_pcr_for_review<-previous(p)
  k<-n<-integer(2)
  for(i in 1:2){v<-z[z$evaluable&z$arm==c('mRNA-1273','Placebo')[i],];n[i]<-nrow(v);k[i]<-sum(v$anti_n_positive)}
  st<-sc_stats(k,n,'Reconstructed from local PHMPT records')
  flow<-data.frame(arm=c('mRNA-1273','Placebo'),selected_confirmed_cases=sapply(c('mRNA-1273','Placebo'),function(a)sum(z$arm==a)),
    evaluable_PDV=n,anti_N_positive=k,not_evaluable=sapply(c('mRNA-1273','Placebo'),function(a)sum(z$arm==a&!z$evaluable)))
  ws_csv(flow,file.path(out,'local_cohort_flow.csv'));ws_csv(st$arms,file.path(out,'local_seroconversion_by_arm.csv'))
  ws_csv(st$contrast,file.path(out,'local_seroconversion_contrasts.csv'))
  ws_csv(transform(z,illness_date=ws_iso(illness_date),pcr_date=ws_iso(pcr_date),pdv_date=ws_iso(pdv_date)),file.path(out,'local_subject_level.csv'))
  status<-data.frame(metric=c('Confirmed cases: mRNA','Confirmed cases: placebo','PDV evaluable: mRNA','PDV evaluable: placebo','Anti-N positive: mRNA','Anti-N positive: placebo'),
    published=c(58,754,52,648,21,605),local=c(flow$selected_confirmed_cases,n,k))
  status$matches<-status$published==status$local;ws_csv(status,file.path(out,'published_vs_local_counts.csv'))
  timing<-do.call(rbind,lapply(c('All evaluable','mRNA-1273','Placebo'),function(a) {
    v<-z$days_illness_to_pdv[z$evaluable&(a=='All evaluable'|z$arm==a)]
    data.frame(group=a,n=length(v),median_days=if(length(v))median(v) else NA_real_,min_days=if(length(v))min(v) else NA_real_,max_days=if(length(v))max(v) else NA_real_)
  }))
  ws_csv(timing,file.path(out,'local_illness_to_pdv_timing.csv'))
  ws_csv(data.frame(metric=c('Median days','Minimum days','Maximum days'),published=c(53,5,150),local=unlist(timing[1,c('median_days','min_days','max_days')],use.names=FALSE)),file.path(out,'published_vs_local_timing.csv'))
  list(subjects=z,stats=st,flow=flow,agreement=status,timing=timing,case_code=case_code)
}
sc_viral <- function(admb,cohort,viral_file,paramcd,root,origin,out) {
  z<-cohort$subjects[cohort$subjects$evaluable,,drop=FALSE]
  if(nzchar(viral_file)) {
    path<-ws_path(root,viral_file);v<-read.csv(path,stringsAsFactors=FALSE,check.names=FALSE,fileEncoding='UTF-8-BOM')
    needed<-c('subject','PCRDT','LOG10_COPIES_ML','SOURCE_REF');if(!all(needed%in%names(v)))stop('Viral CSV requires: ',paste(needed,collapse=', '))
    v$subject<-ws_sid(v$subject);v$pcr_date<-ws_dates(v$PCRDT,'viral.PCRDT',origin,TRUE)
    v$log10_copies_ml<-suppressWarnings(as.numeric(v$LOG10_COPIES_ML));v$source_ref<-ws_text(v$SOURCE_REF)
    if(any(!is.finite(v$log10_copies_ml))||any(!nzchar(v$source_ref)))stop('Viral CSV contains missing/non-numeric exact log10 values or missing source references.')
    source_label<-basename(path);ws_csv(data.frame(path=normalizePath(path,winslash='/'),md5=unname(tools::md5sum(path))),file.path(out,'viral_input_checksum.csv'))
  } else {
    if(is.null(admb))stop('No quantitative viral-load source supplied. Qualitative PCR results cannot reproduce the viral-load regression.')
    sc_inventory(admb,file.path(out,'admb_parameter_inventory.csv'))
    lab<-ws_norm(admb$PARAM);units<-ws_norm(paste(admb$AVALU,admb$MBSTRESU,admb$MBORRESU))
    auto<-grepl('LOG[ _-]*10',paste(lab,units))&grepl('COP',paste(lab,units))&grepl('ML',paste(lab,units))&!grepl('CYCLE|THRESHOLD',lab)
    code<-if(nzchar(paramcd))paramcd else unique(admb$PARAMCD[auto])
    if(length(code)!=1L||!any(admb$PARAMCD==code))stop('No unique log10 copies/mL parameter identified. Review admb_parameter_inventory.csv or supply PHMPT_VIRAL_FILE. Ct-to-copy conversion is not inferred.')
    a<-admb[admb$PARAMCD==code,,drop=FALSE]
    # Even an explicit parameter must state the scale and units; otherwise use
    # the reviewed CSV interface with documented calibration/provenance.
    if(!all(auto[admb$PARAMCD==code]))stop('Selected viral parameter does not explicitly identify log10 copies/mL. Supply a reviewed viral CSV after validating units.')
    d<-ws_dates(a$MBDTC,'quantitative MBDTC',origin);other<-ws_dates(a$ADT,'quantitative ADT',origin)
    bad<-nzchar(a$MBDTC)&!is.finite(d);conflict<-is.finite(d)&is.finite(other)&d!=other
    d[!nzchar(a$MBDTC)]<-other[!nzchar(a$MBDTC)]
    vals<-suppressWarnings(as.numeric(a$AVAL))
    day1<-grepl('^ILLNESS VISIT DAY 1($| )',ws_norm(a$VISIT))
    keep<-is.finite(d)&is.finite(vals)&!bad&!conflict&!nzchar(a$DTYPE)&!grepl('[<>]',a$AVALC)&day1
    v<-data.frame(subject=a$subject[keep],pcr_date=d[keep],log10_copies_ml=vals[keep],source_ref=paste0('ADMB:P301-row:',which(admb$PARAMCD==code)[keep],':',code))
    source_label<-paste0('ADMB parameter ',code)
    ws_csv(a[!keep,,drop=FALSE],file.path(out,'viral_values_excluded.csv'))
  }
  v<-unique(v[c('subject','pcr_date','log10_copies_ml','source_ref')]);key<-paste(v$subject,v$pcr_date,sep='|')
  if(anyDuplicated(key))stop('Several viral values match one subject/diagnosis date. Resolve sample identity before fitting; no averaging or maximum is substituted.')
  j<-match(paste(z$subject,z$pcr_date,sep='|'),key)
  z$log10_copies_ml<-v$log10_copies_ml[j];z$viral_source_ref<-v$source_ref[j]
  ws_csv(transform(z,pcr_date=ws_iso(pcr_date),pdv_date=ws_iso(pdv_date),illness_date=ws_iso(illness_date)),file.path(out,'viral_model_subjects.csv'))
  d<-z[is.finite(z$log10_copies_ml),];d$arm<-factor(d$arm,levels=c('mRNA-1273','Placebo'))
  if(nrow(d)<5||length(unique(d$arm))!=2L||length(unique(d$anti_n_positive))!=2L)stop('Insufficient complete participant-level data to fit the viral-load model.')
  if(length(unique(d$log10_copies_ml))<2)stop('No variation in viral load.')
  warnings_seen<-character()
  fit<-withCallingHandlers(glm(anti_n_positive~arm+log10_copies_ml,family=binomial(),data=d),warning=function(w){warnings_seen<<-c(warnings_seen,conditionMessage(w));invokeRestart('muffleWarning')})
  if(!fit$converged||length(warnings_seen))stop('Viral-load model requires review: ',paste(warnings_seen,collapse='; '))
  co<-summary(fit)$coefficients
  tab<-data.frame(term=rownames(co),log_odds=co[,1],se=co[,2],OR=exp(co[,1]),ci_lower=exp(co[,1]-qnorm(.975)*co[,2]),
    ci_upper=exp(co[,1]+qnorm(.975)*co[,2]),p=co[,4],n=nrow(d),omitted_missing_load=nrow(z)-nrow(d),source=source_label,CI_method='Wald on log-odds scale',row.names=NULL)
  if(any(!is.finite(tab$OR)|!is.finite(tab$ci_lower)|!is.finite(tab$ci_upper)))stop('Non-finite regression estimate; inspect separation and quantitative values.')
  ws_csv(tab,file.path(out,'local_viral_load_model.csv'));tab
}
sc_reconstruct <- function(root,out,paths,case_param,viral_param,viral_file,cutoff,encoding,origin,min_n) {
  adsl<-ws_pick(root,'Prod','ADSL',paths$ADSL,out,encoding,origin,min_n)
  ctx<-ws_context(adsl$x,origin);ctx$PDVDT<-ws_dates(ctx$PDVDT,'PDVDT',origin,TRUE)
  # sc_pdv parses string dates; restore ISO for that optional field.
  ctx$PDVDT<-ws_iso(ctx$PDVDT)
  if(unique(ctx$CUTOFFDT)!=cutoff)stop('Choose the complete 26 March 2021 BLA snapshot. Selected ADSL cutoff: ',ws_iso(unique(ctx$CUTOFFDT)))
  manifest<-adsl$manifest;submission<-ws_submission(manifest$path)
  on.exit(ws_csv(manifest,file.path(out,'source_manifest.csv')),add=TRUE)
  get<-function(domain,optional=FALSE) {
    a<-if(optional)nt_optional(root,domain,paths[[domain]],out,encoding,origin,min_n,submission) else ws_pick(root,'Prod',domain,paths[[domain]],out,encoding,origin,min_n,submission)
    if(is.null(a))return(NULL);manifest<<-rbind(manifest,a$manifest);a$x
  }
  e<-get('ADEFF');code<-sc_case_param(e,case_param,out)
  raw_s<-get('IS');as<-get('ADIS',TRUE);s<-nt_reconcile(raw_s,as,'serology',origin,out)
  raw_p<-get('MB');ap<-get('ADMB',TRUE);p<-nt_reconcile(raw_p,ap,'pcr',origin,out)
  if(!is.null(ap))sc_inventory(ap,file.path(out,'admb_parameter_inventory.csv'))
  result<-sc_cohort(ctx,e,s,p,code,cutoff,origin,out)
  v<-tryCatch(sc_viral(ap,result,viral_file,viral_param,root,origin,out),error=function(e)e)
  result['viral']<-list(if(inherits(v,'error'))NULL else v)
  result$viral_status<-if(inherits(v,'error'))conditionMessage(v) else 'Fitted from participant-level records; compare estimates and cohort with the paper.'
  result
}

sc_html <- function(x) {
  x<-gsub('&','&amp;',as.character(x),fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE);x<-gsub('>','&gt;',x,fixed=TRUE);gsub('"','&quot;',x,fixed=TRUE)
}
sc_html_table <- function(d)paste0('<table><thead><tr>',paste0('<th>',sc_html(names(d)),'</th>',collapse=''),'</tr></thead><tbody>',
  paste(apply(d,1,function(r)paste0('<tr>',paste0('<td>',sc_html(r),'</td>',collapse=''),'</tr>')),collapse='\n'),'</tbody></table>')
sc_display <- function(z,published=FALSE) {
  a<-z$arms
  rows<-list('Evaluable anti-N results at PDV'=as.character(a$evaluable),
    'Anti-N positive'=mapply(sc_pct,a$positive,a$evaluable),
    'Anti-N negative'=mapply(sc_pct,a$negative,a$evaluable),
    'Recalculated 95% CI (exact binomial)'=sprintf('%.1f%% to %.1f%%',100*a$ci_lower,100*a$ci_upper))
  if(published)rows<-c(list('PCR-confirmed cases before PDV (published)'=c('58','754'),
    'No PDV (published)'=c('0','51'),'PDV but no serology (published)'=c('6','55')),rows,
    list('95% CI reported in the abstract'=c('27% to 54%','92% to 95%')))
  data.frame(Measure=names(rows),`mRNA-1273`=vapply(rows,`[`,character(1),1),Placebo=vapply(rows,`[`,character(1),2),row.names=NULL,check.names=FALSE)
}
sc_pdf_page <- function(table,title,subtitle,contrast,viral_line,notes) {
  grid::grid.newpage()
  txt<-function(label,x,y,size=11,bold=FALSE,col='#111111',just='left')grid::grid.text(label,x=x,y=y,just=just,
    gp=grid::gpar(fontsize=size,fontface=if(bold)'bold' else 'plain',col=col))
  txt('MODERNA P301 / ANTI-N SEROLOGY',.055,.95,10,TRUE,'#555555')
  txt(title,.055,.902,18,TRUE);txt(subtitle,.055,.86,10)
  x<-c(.055,.57,.79);headers<-names(table);top<-.805;rowh<-.046
  grid::grid.rect(x=.5,y=top,width=.90,height=.047,gp=grid::gpar(fill='#E7ECEF',col=NA))
  for(j in 1:3)txt(headers[j],x[j],top,10,TRUE)
  for(i in seq_len(nrow(table))) {
    y<-top-i*rowh
    if(i%%2==0)grid::grid.rect(x=.5,y=y,width=.90,height=rowh,gp=grid::gpar(fill='#F4F6F7',col=NA))
    for(j in 1:3)txt(as.character(table[i,j]),x[j],y,if(j==1)9.5 else 10,FALSE)
  }
  bottom<-top-(nrow(table)+1)*rowh
  txt(paste0('Fisher exact test (two-sided): p ',if(startsWith(sc_p(contrast$fisher_p),'<'))'' else '= ',sc_p(contrast$fisher_p),
    '     |     Difference: ',sprintf('%+.1f',contrast$difference_pp),' percentage points'),.055,bottom,10,TRUE)
  for(i in seq_along(viral_line))txt(viral_line[i],.055,bottom-.055-.026*(i-1),10,i==1)
  wrapped<-unlist(lapply(notes,function(n)strwrap(n,width=143)))
  y<-bottom-.15
  for(line in wrapped){txt(line,.055,y,8.4,FALSE,'#444444');y<-y-.022}
  txt('Source: Follmann et al. Ann Intern Med. 2022. doi:10.7326/M22-1300 | Script v1.0.0',.055,.035,8,FALSE,'#555555')
}
sc_report <- function(b,local,status,out,dpi) {
  table<-sc_display(b,TRUE);ws_csv(table,file.path(out,'anti_n_seroconversion_print_table.csv'))
  notes<-c('PDV = participant decision visit. The 700 evaluable participants are a subset of 812 confirmed cases, not the entire randomized population.',
    'Percentages, exact binomial confidence intervals and Fisher p are recalculated from the published 2 x 2 counts. The abstract intervals are quoted separately; they are not the exact-binomial intervals.',
    'The reported median illness-to-PDV interval is 53 days (range 5-150); this is not a sample-selection window. No participant dates are reconstructed from that median.',
    'The viral-load OR is adjusted for randomized arm. It cannot be recalculated from arm-level counts. Negative anti-N at PDV does not establish absence of earlier infection or absence of all previous antibody responses.',
    'These results describe anti-N detection among confirmed illnesses. They do not count PCR-missed infections or supply a correction factor for vaccine efficacy.')
  viral_line<-c('Published viral-load association: OR 1.90 (95% CI 1.59-2.28)',
    'Per 1 log10 increase in copies/mL at illness presentation; not refitted from aggregate counts.')
  pdf<-file.path(out,'anti_n_seroconversion_table.pdf')
  if(capabilities('cairo'))grDevices::cairo_pdf(pdf,width=11.69,height=8.27,family='sans',onefile=TRUE) else grDevices::pdf(pdf,width=11.69,height=8.27,family='Helvetica')
  tryCatch({
    sc_pdf_page(table,'Anti-N detection after PCR-confirmed COVID-19','Published-count recalculation | Elecsys assay | Blinded phase through 26 March 2021',b$contrast,viral_line,notes)
    if(!is.null(local)) {
      vl<-if(is.null(local$viral))c('Viral-load regression not reproduced.',local$viral_status) else {
        z<-local$viral[local$viral$term=='log10_copies_ml',];c(sprintf('Local viral-load OR %.3f (95%% CI %.3f-%.3f)',z$OR,z$ci_lower,z$ci_upper),paste('Complete cases:',z$n,'; missing viral load:',z$omitted_missing_load))
      }
      sc_pdf_page(sc_display(local$stats),'PHMPT participant-level reconstruction',paste('Selected ADEFF parameter:',local$case_code,'| Review agreement and selection CSVs'),local$stats$contrast,
        unlist(lapply(vl,strwrap,width=125)),c('Aggregate agreement with the paper does not independently establish identical participant membership.',
          'PDV serology is retained on the unblinding boundary. Missing, conflicting and later samples are not counted as negative.',
          'Exact-binomial intervals are analyst recalculations; the optional logistic model uses Wald intervals.',
          'Earlier positive records are flagged for review in local_subject_level.csv. No estimate is forced to match the published targets.'))
    }
  },finally=grDevices::dev.off())
  grDevices::png(file.path(out,'anti_n_seroconversion_table.png'),width=11.69,height=8.27,units='in',res=dpi,
    type=if(capabilities('cairo'))'cairo' else getOption('bitmapType','Xlib'))
  tryCatch(sc_pdf_page(table,'Anti-N detection after PCR-confirmed COVID-19','Published-count recalculation | Elecsys assay | Blinded phase through 26 March 2021',b$contrast,viral_line,notes),finally=grDevices::dev.off())
  paragraph<-paste0('However, seroconversion alone did not satisfy the primary symptomatic COVID-19 endpoint, which required the specified symptoms together with a positive PCR result; conversely, a negative anti-N result could not reliably exclude previous infection. Follmann et al. illustrated this limitation in 700 participants with PCR-confirmed COVID-19 and evaluable serology during the blinded phase: at the participant decision visit, a median of 53 days after diagnosis, anti-N antibodies were detected in 21 of 52 mRNA-1273 recipients (40%; reported 95% CI, 27%-54%) versus 605 of 648 placebo recipients (93%; 92%-95%). Thus, 31 vaccine recipients (59.6%) and 43 placebo recipients (6.6%) were anti-N-negative at that assessment despite confirmed illness. A tenfold increase in viral copies at diagnosis was associated with 90% higher odds of anti-N positivity (arm-adjusted OR, 1.90; 95% CI, 1.59-2.28). These results show why anti-N testing could undercount previous infection differently between the arms, particularly in vaccine recipients. They do not establish how many infections escaped PCR detection, or justify applying these proportions directly to asymptomatic infections or recalculating vaccine efficacy.')
  introduction<-paste0('The following table checks the published comparison of anti-N positivity after PCR-confirmed COVID-19. The script recalculates each arm\'s proportion, the number remaining anti-N-negative, exact binomial 95% confidence intervals and a two-sided Fisher exact test from the reported counts. It also attempts a separate reconstruction from the full PHMPT production tables, retaining serology collected at the participant decision visit. The viral-load association is shown as a published reference unless participant-level quantitative data permit the arm-adjusted logistic model to be refitted. Missing data and differences from the published cohort are reported explicitly.')
  writeLines(c(paragraph,'','Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC9258784/'),file.path(out,'audit_paragraph.txt'))
  writeLines(introduction,file.path(out,'table_introduction.txt'))
  html<-c('<!doctype html><html lang="en"><meta charset="utf-8"><title>Anti-N seroconversion audit</title>',
    '<style>body{font:17px/1.55 Georgia,serif;max-width:1000px;margin:30px auto;padding:0 24px;color:#171717}h1,h2{font-family:Arial,sans-serif}table{border-collapse:collapse;width:100%;font:14px/1.4 Arial;margin:22px 0}th,td{padding:10px;border-bottom:1px solid #ccc;text-align:left}th{background:#e7ecef}tr:nth-child(even){background:#f5f6f7}.status{font:14px/1.5 Arial;background:#f1f3f5;padding:16px}.notes{font:13px/1.5 Arial}button{font:14px Arial;padding:9px;margin:6px}a{color:#165c87}@media print{.controls{display:none}body{font-size:11pt;margin:0;padding:0}table{font-size:9pt}tr{break-inside:avoid}}</style>',
    '<div class="controls"><button onclick="selectBlock(\'copy\')">Select audit paragraph and table</button><a href="anti_n_seroconversion_table.pdf">Print-ready PDF</a></div>',
    paste0('<p class="status">Participant-level status: ',sc_html(status),'</p>'),'<article id="copy"><h1>Anti-N detection after confirmed illness</h1>',
    paste0('<p>',sc_html(paragraph),'</p><p><a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC9258784/">Follmann et al., 2022</a></p>'),
    paste0('<p>',sc_html(introduction),'</p><h2>Published-count recalculation</h2>'),sc_html_table(table),
    paste0('<p>Fisher p ',if(startsWith(sc_p(b$contrast$fisher_p),'<'))'' else '= ',sc_html(sc_p(b$contrast$fisher_p)),
      '; difference: ',sprintf('%+.1f',b$contrast$difference_pp),' percentage points (mRNA-1273 minus placebo).</p>'),
    '<p><b>Published viral-load OR: 1.90 (95% CI 1.59-2.28).</b> This reference value has not been refitted from aggregate counts.</p>',
    '<div class="notes">',paste0('<p>',sc_html(notes),'</p>'),'</div>')
  if(!is.null(local)) {
    html<-c(html,'<h2>PHMPT reconstruction</h2>',sc_html_table(sc_display(local$stats)),sc_html_table(local$agreement),sc_html_table(local$timing),
      paste0('<p>',sc_html(local$viral_status),'</p>'))
    if(!is.null(local$viral))html<-c(html,sc_html_table(local$viral[local$viral$term!='(Intercept)',c('term','OR','ci_lower','ci_upper','p','n')]))
  }
  html<-c(html,'</article><script>function selectBlock(id){const r=document.createRange();r.selectNodeContents(document.getElementById(id));const s=getSelection();s.removeAllRanges();s.addRange(r);}</script></html>')
  writeLines(html,file.path(out,'anti_n_seroconversion_report.html'),useBytes=TRUE)
}

analyse_anti_n_seroconversion <- function(root=PHMPT_ROOT,output_dir=PHMPT_OUTPUT_DIR,
  reconstruct=PHMPT_RECONSTRUCT,adsl_file=PHMPT_ADSL_FILE,adeff_file=PHMPT_ADEFF_FILE,
  is_file=PHMPT_IS_FILE,adis_file=PHMPT_ADIS_FILE,mb_file=PHMPT_MB_FILE,admb_file=PHMPT_ADMB_FILE,
  case_param=PHMPT_CASE_PARAMCD,viral_param=PHMPT_VIRAL_PARAMCD,viral_file=PHMPT_VIRAL_FILE,
  cutoff=PHMPT_BLA_CUTOFF,encoding=PHMPT_CSV_ENCODING,origin=PHMPT_NUMERIC_DATE_ORIGIN,
  min_randomized=PHMPT_MIN_RANDOMIZED,png_dpi=PHMPT_PNG_DPI,show_report=PHMPT_SHOW_REPORT) {
  if(length(png_dpi)!=1L||!is.finite(png_dpi)||png_dpi<72)stop('Invalid PNG resolution.')
  if(!nzchar(root))root<-tryCatch(ws_root(''),error=function(e)getwd())
  root<-normalizePath(root,winslash='/',mustWork=TRUE)
  out<-if(nzchar(output_dir))ws_path(root,output_dir) else file.path(root,'analysis',paste0('anti_n_seroconversion_',format(Sys.time(),'%Y%m%d_%H%M%S')))
  if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop('Use a new output folder; existing results are not overwritten: ',out)
  dir.create(out,recursive=TRUE,showWarnings=FALSE);writeLines('RUNNING',file.path(out,'RUN_STATUS.txt'))
  complete<-FALSE;on.exit(if(!complete)writeLines('FAILED',file.path(out,'RUN_STATUS.txt')),add=TRUE)
  b<-sc_benchmark(out);local<-NULL
  cutoff<-ws_dates(cutoff,'BLA cutoff',strict=TRUE)
  if(length(cutoff)!=1L||!is.finite(cutoff))stop('Invalid cutoff.')
  status<-'Published-count calculation only; participant-level reconstruction was disabled.'
  if(reconstruct) {
    paths<-list(ADSL=adsl_file,ADEFF=adeff_file,IS=is_file,ADIS=adis_file,MB=mb_file,ADMB=admb_file)
    result<-tryCatch(sc_reconstruct(root,out,paths,case_param,viral_param,viral_file,cutoff,encoding,origin,min_randomized),error=function(e)e)
    if(inherits(result,'error')) {
      status<-paste('PHMPT reconstruction unavailable/incomplete:',conditionMessage(result))
      writeLines(status,file.path(out,'RECONSTRUCTION_GAP.txt'));message(status)
    } else {
      local<-result;status<-paste(sum(local$agreement$matches),'of',nrow(local$agreement),'published aggregate counts match. Subject identity is not independently verified against the authors\' analysis file.',local$viral_status)
    }
  }
  sc_report(b,local,status,out,png_dpi)
  ws_csv(data.frame(subject=character(),PCRDT=character(),LOG10_COPIES_ML=numeric(),SOURCE_REF=character()),file.path(out,'viral_load_input_template.csv'))
  writeLines(c('analyse_anti_n_seroconversion.R v1.0.0',status,
    'Published source: https://pmc.ncbi.nlm.nih.gov/articles/PMC9258784/; DOI 10.7326/M22-1300.',
    'Recalculated confidence intervals: Clopper-Pearson exact binomial. Published intervals remain separate; no rounding is tuned to force agreement.',
    'Case source: full ADEFF, primary adjudicated COVID parameter, ANL01FL=Y, non-derived/non-censored event; ADT >= dose2+14; before PDV and within blinded follow-up.',
    'Population: randomized FAS, baseline PCR and Elecsys negative. FAS is used to investigate the paper\'s 58/754 case totals; an extra PPROTFL restriction is not silently imposed.',
    'A positive PCR on ADEFF.PCRDT is required for local verification. Missing/conflicting lab evidence is recorded in the selection audit.',
    'Serology: original IS SARSCOV2, reconciled with ADIS Elecsys. Explicit PDV/OL-D1 sample at the boundary is included; later open-label samples are excluded.',
    'The PDV result is used, not any eventual seroconversion. Multiple PDV dates, inconsistent results or missing samples require review and are never counted as negative.',
    '53 days is a reported median, not a required collection day or an invented +/- window. Cases outside 5-150 days are retained if valid and exposed in the subject-level audit.',
    'Viral model: anti_n_positive ~ randomized arm + log10 copies/mL on the ADEFF diagnostic PCR date. Complete cases; Wald CIs. No Ct-to-copy calibration, below-limit substitution or time-window substitution is invented.',
    'Supply a reviewed viral CSV only with exact LOG10_COPIES_ML values and SOURCE_REF. Do not fill the template using the published odds ratio.',
    'The previous prod_anti_n_pcr export alone is insufficient: candidate_adeff.csv is a selected subset and the export lacks quantitative viral copies.',
    'The script inventories ADEFF/ADMB parameters and never automatically changes endpoint definitions to obtain target counts.'),file.path(out,'README.txt'))
  capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
  ws_csv(data.frame(parameter=c('script_version','BLA_cutoff','reconstruct','CI_method'),value=c('1.0.0',ws_iso(cutoff),reconstruct,'Exact binomial')),file.path(out,'run_parameters.csv'))
  complete<-TRUE;writeLines(c('COMPLETE_PUBLISHED_COUNT_RECALCULATION',if(is.null(local))'LOCAL_RECONSTRUCTION_NOT_COMPLETED' else 'LOCAL_RECONSTRUCTION_COMPLETED_REVIEW_AGREEMENT',status),file.path(out,'RUN_STATUS.txt'))
  message('Report: ',normalizePath(file.path(out,'anti_n_seroconversion_report.html'),winslash='/'))
  if(show_report)utils::browseURL(normalizePath(file.path(out,'anti_n_seroconversion_report.html'),winslash='/'))
  invisible(list(output_dir=out,benchmark=b,local=local,status=status))
}
# ---- RUN: whole-script execution starts the analysis ----
if(isTRUE(getOption('phmpt.autorun',TRUE)))PHMPT_ANTI_N_RESULT<-analyse_anti_n_seroconversion()
