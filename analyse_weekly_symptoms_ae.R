# Moderna P301 - weekly doses, symptom groups and AE recording share - v1.1.0
# Stand-alone, R >= 4.3: save in the project root, then Select All / Run in RStudio.
# CSV/JMP and charts: base R. XPT/SAS7BDAT: install.packages("haven").
# No positivity filter. This measures documented symptoms, not COVID diagnoses.

# ---- SETTINGS ----
PHMPT_ROOT <- ""                 # Auto-detect, or D:/Projects/moderna_trial_analysis
PHMPT_RELEASE <- "auto"          # Prod preferred; otherwise EUA. Never pools releases.
PHMPT_ADSL_FILE <- ""            # Optional explicit path; CSV, XPT, SAS7BDAT or JMP
PHMPT_AE_FILE <- ""              # Raw SDTM AE (not ADAE)
PHMPT_SYMPTOMS_FILE <- ""        # ADaM ADSYMP: includes the COVID symptom checklist
PHMPT_MAPPING_FILE <- ""         # Defaults to assets/ae_covid_symptom_mapping.csv
# If no external mapping is present, the embedded, previously supplied v1.0.0
# field-specific mapping is used and exported. No other R script is sourced.
PHMPT_WINDOWS <- c("BLINDED", "EUA", "FULL")
# Audit V09/V10 PNGs default to FULL so post-unblinding doses are visible.
# Set options(phmpt.audit_window="BLINDED") before running for blinded PNGs.
PHMPT_EUA_CUTOFF <- "2020-11-25" # EUA data cutoff, not authorization date
PHMPT_AE_TIMING <- "interval"    # "interval": known active AE weeks; "onset": onset week only
PHMPT_OUTPUT_DIR <- ""           # Default: analysis/weekly_symptoms_ae_<timestamp>
PHMPT_CSV_ENCODING <- "auto"     # UTF-8-BOM then Windows-1252
PHMPT_NUMERIC_DATE_ORIGIN <- "" # Unformatted numeric CSV dates: set documented origin
PHMPT_MIN_RANDOMIZED <- 30000L   # Reject incomplete ADSL subsets
PHMPT_PNG_DPI <- 300L
PHMPT_SHOW_CHART <- interactive()
# To load functions only: options(phmpt.autorun=FALSE) before source().

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
  ADSL=list(required=c("STUDYID","USUBJID","RANDFL","TRT01P","RANDDT","DOSE1FL","DOSE2FL","TR01SDT","DOSE2DT","CUTOFFDT"),
    optional=c("SUBJID","OLDOS1FL","OLDOS2FL","OLDOS2DT","TRT02A","AP01EDT","UBLDDT","AP02SDT","TR02SDT","EUAVACDT","EOSDT","DTHDT")),
  AE=list(required=c("STUDYID","USUBJID","AESEQ","AETERM","AEDECOD","AESTDTC"),
    optional=c("AELLT","AEENDTC","AEENRF","AEOUT","AESER","AEREL")),
  ADSYMP=list(required=c("STUDYID","USUBJID","PARAMCD","AVALC","ADT"),
    optional=c("ASEQ","PARAM","ASTDT","AENDT","SRCDOM","SRCSEQ","SRCVAR","LNKID","SEV","CUTOFFDT")))
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
  source_fields<-names(x)
  missing<-setdiff(spec$required,names(x));if(length(missing))stop("Missing ",domain," fields: ",paste(missing,collapse=", "))
  x[]<-lapply(x,ws_text)
  x<-x[ws_norm(x$STUDYID)=="MRNA-1273-P301",,drop=FALSE]
  if(!nrow(x))stop("No P301 rows in this table.")
  for(n in setdiff(spec$optional,names(x)))x[[n]]<-rep("",nrow(x))
  x<-x[c(spec$required,spec$optional)];x$subject<-ws_sid(x$USUBJID)
  for(n in intersect(c("AESEQ","ASEQ","SRCSEQ"),names(x))) {
    ok<-grepl("^[0-9]+([.]0+)?$",x[[n]]);x[[n]][ok]<-as.character(as.numeric(x[[n]][ok]))
  }
  if(domain=="ADSL"&&anyDuplicated(x$subject))stop("Duplicate subject in ADSL.")
  if(domain=="AE") {
    if(any(!nzchar(x$AESEQ)))stop("Missing AESEQ.")
    x<-unique(x)
    if(anyDuplicated(paste(x$subject,x$AESEQ,sep="|")))stop("Conflicting rows share subject/AESEQ.")
  }
  attr(x,"source_fields")<-source_fields
  x
}
ws_candidates <- function(root,release,domain) {
  p<-list.files(file.path(root,if(release=="Prod")"prod_data" else "eua_data"),recursive=TRUE,full.names=TRUE,
    pattern="\\.(csv|xpt|sas7bdat|jmp)$",ignore.case=TRUE)
  token<-if(domain=="ADSYMP")"(adsymp|adsy)" else tolower(domain)
  p<-p[grepl(paste0("(^|[_ -])",token,"([_. -]|$)"),basename(p),ignore.case=TRUE)]
  p<-p[!grepl("p101|p201|dmid",basename(p),ignore.case=TRUE)]
  p[order(match(tolower(tools::file_ext(p)),c("csv","jmp","xpt","sas7bdat")),p)]
}
ws_submission <- function(p) {
  z<-regmatches(tolower(basename(p)),regexpr("[0-9]+_s[0-9]+_m[0-9]+_mrna-1273-p301",tolower(basename(p))))
  if(length(z))z else ""
}
ws_signature <- function(x,origin="") {
  attr(x,"source_fields")<-NULL
  dates<-intersect(c("RANDDT","TR01SDT","DOSE2DT","OLDOS2DT","CUTOFFDT","AP01EDT","UBLDDT","AP02SDT","TR02SDT","EUAVACDT","EOSDT","DTHDT","AESTDTC","AEENDTC","ADT","ASTDT","AENDT"),names(x))
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
  if(domain=="ADSL"&&length(valid)>1L) {
    full<-vapply(valid,function(v)all(c("OLDOS1FL","OLDOS2FL","OLDOS2DT","TR02SDT")%in%attr(v$x,"source_fields")),logical(1))
    drop<-rep(FALSE,length(valid))
    if(any(full)&&any(!full))for(j in which(!full))for(k in which(full)) {
      common<-intersect(attr(valid[[j]]$x,"source_fields"),attr(valid[[k]]$x,"source_fields"))
      common<-sort(intersect(c(common,"subject"),intersect(names(valid[[j]]$x),names(valid[[k]]$x))))
      if(identical(ws_signature(valid[[j]]$x[common],origin),ws_signature(valid[[k]]$x[common],origin))) {
        drop[j]<-TRUE;log$status[valid[[j]]$index]<-"INCOMPLETE_MIRROR_NOT_USED"
        log$detail[valid[[j]]$index]<-"Common fields agree with a complete ADSL containing open-label dose fields.";break
      }
    }
    valid<-valid[!drop]
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
ws_mapping <- function(root,path,out) {
  if(!nzchar(path)) {
    p<-file.path(root,c("assets/ae_covid_symptom_mapping.csv","ae_covid_symptom_mapping.csv"))
    path<-if(any(file.exists(p)))p[which(file.exists(p))[1]] else ""
  } else path<-ws_path(root,path)
  if(nzchar(path)) {
    if(!file.exists(path))stop("Mapping asset not found: ",path)
    x<-read.csv(path,colClasses="character",check.names=FALSE,fileEncoding="UTF-8-BOM",na.strings=c("","NA"))
    source<-normalizePath(path,winslash="/");hash<-unname(tools::md5sum(path))
  } else {x<-ws_embedded_mapping();source<-"Embedded field-specific mapping v1.0.0";hash<-"";message(source)}
  need<-c("match_field","normalized_term","symptom_code","symptom_group","match_class","include_default")
  if(!all(need %in% names(x)))stop("Mapping lacks required fields: ",paste(setdiff(need,names(x)),collapse=", "))
  x[]<-lapply(x,ws_text);x$normalized_term<-ws_norm(x$normalized_term)
  x$match_field<-ws_norm(x$match_field);x$include_default<-ws_norm(x$include_default)
  x$symptom_code<-ws_norm(x$symptom_code);x$symptom_group<-ws_norm(x$symptom_group)
  x$match_class<-ws_norm(x$match_class)
  if(any(!nzchar(x$normalized_term)))stop("Mapping has an empty normalized_term.")
  if(any(!x$include_default %in% c("TRUE","FALSE")))stop("Mapping include_default must be TRUE/FALSE.")
  x$include_default<-x$include_default=="TRUE"
  if(any(!x$match_field %in% c("AETERM","AELLT","AEDECOD")))stop("Unknown mapping match_field.")
  if(any(!x$symptom_code %in% names(ws_param_map)))stop("Unknown mapping symptom_code.")
  if(any(x$symptom_group!=unname(ws_param_map[x$symptom_code])))stop("Mapping symptom groups differ from the agreed 12-group classification.")
  if(any(x$include_default & x$match_class!="DIRECT"))stop("Review-only mapping rows cannot be included without explicit reclassification.")
  # Several verbatim spellings can define the SAME normalized matching rule.
  # Keep the complete source rows; deduplicate by decisions, not by spelling,
  # labels, notes or dictionary provenance. Never silently resolve disagreement.
  key<-paste(x$match_field,x$normalized_term,x$symptom_code,sep="\r")
  group<-match(key,unique(key));groups<-split(seq_len(nrow(x)),group)
  x$mapping_source_row<-seq_len(nrow(x))
  x$mapping_rule_id<-sprintf("MAP%04d",group)
  ws_csv(x,file.path(out,"ae_covid_symptom_mapping_input.csv"))
  # Duration/qualification metadata are retained and checked for agreement even
  # though this burden analysis does not impose protocol-trigger eligibility.
  decisions<-intersect(c("symptom_group","match_class","include_default",
                        "min_duration_hours","qualification_check"),names(x))
  conflicts<-vapply(groups,function(ix) {
    z<-x[ix,decisions,drop=FALSE]
    z[]<-lapply(z,function(v)if(is.logical(v))v else ws_norm(v))
    nrow(unique(z))>1L
  },logical(1))
  bad<-unlist(groups[conflicts],use.names=FALSE)
  conflict_file<-file.path(out,"ae_covid_symptom_mapping_conflicts.csv")
  ws_csv(x[bad,,drop=FALSE],conflict_file)
  if(length(bad))stop("Conflicting mapping decisions for ",sum(conflicts),
    " normalized rule(s). Review ",normalizePath(conflict_file,winslash="/",mustWork=FALSE),
    ". No conflicting rule was selected automatically.")
  keep<-!duplicated(key);used<-x[keep,,drop=FALSE]
  used$mapping_source_rows<-vapply(group[keep],function(g)paste(which(group==g),collapse=";"),character(1))
  if("ae_term" %in% names(x))used$mapping_term_variants<-vapply(group[keep],function(g)
    paste(unique(x$ae_term[group==g]),collapse=" | "),character(1))
  collapsed<-x[duplicated(key)|duplicated(key,fromLast=TRUE),,drop=FALSE]
  ws_csv(collapsed,file.path(out,"ae_covid_symptom_mapping_duplicates.csv"))
  ws_csv(used,file.path(out,"ae_covid_symptom_mapping_used.csv"))
  ws_csv(data.frame(path=source,md5=hash,input_rows=nrow(x),effective_rules=nrow(used),
    equivalent_rows_collapsed=nrow(x)-nrow(used)),file.path(out,"mapping_source.csv"))
  if(nrow(x)>nrow(used))message("Mapping: ",nrow(x)," source rows -> ",nrow(used),
    " effective rules; ",nrow(x)-nrow(used)," equivalent duplicate row(s) collapsed. Original rows retained in mapping audit CSVs.")
  used
}
ws_param_map <- c(FEVER="FEVER",CHILLS="CHILLS",COUGH="COUGH",DIBREATH="DYSPNEA",SHBREATH="DYSPNEA",
  FATIGUE="FATIGUE",BODYACHE="MYALGIA",MYALGIA="MYALGIA",HEADACHE="HEADACHE",NLSMELL="SMELL_TASTE",
  NLTASTE="SMELL_TASTE",STHROAT="SORE_THROAT",NASALCON="NASAL",RUNNOSE="NASAL",NAUSEA="NAUSEA_VOMITING",
  VOMIT="NAUSEA_VOMITING",DIARRHEA="DIARRHEA")
ws_context <- function(x,origin="") {
  for(n in c("RANDFL","DOSE1FL","DOSE2FL")) {
    x[[n]]<-ws_norm(x[[n]]);if(any(!x[[n]] %in% c("Y","N")))stop(n,": missing/unexpected flag.")
  }
  dates<-c("RANDDT","TR01SDT","DOSE2DT","OLDOS2DT","CUTOFFDT","AP01EDT","UBLDDT","AP02SDT","TR02SDT","EUAVACDT","EOSDT","DTHDT")
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
  ws_crossover_context(x)
}
# Doses 3/4 are the two documented open-label injections in the ORIGINAL placebo
# arm, not boosters inferred from EUAVACDT or a change of treatment assignment.
# Source date/flag conflicts are exported and excluded, never silently repaired.
ws_crossover_context <- function(x) {
  source_fields<-attr(x,"source_fields");if(is.null(source_fields))source_fields<-names(x)
  for(n in c("OLDOS1FL","OLDOS2FL","TRT02A"))if(!n%in%names(x))x[[n]]<-rep("",nrow(x))
  if(!"OLDOS2DT"%in%names(x))x$OLDOS2DT<-rep(NA_real_,nrow(x))
  for(n in c("OLDOS1FL","OLDOS2FL"))x[[n]]<-ws_norm(x[[n]])
  for(d in 3:4) {
    flag<-if(d==3)"OLDOS1FL" else "OLDOS2FL";date<-if(d==3)"TR02SDT" else "OLDOS2DT"
    f<-x[[flag]];dt<-x[[date]];status<-rep("NO_RECORDED_DOSE",nrow(x))
    status[!all(c(flag,date)%in%source_fields)]<-"SOURCE_FIELDS_UNAVAILABLE"
    status[nzchar(f)&!f%in%c("Y","N")]<-"UNEXPECTED_DOSE_FLAG"
    status[(f=="Y")!=is.finite(dt)]<-"DOSE_FLAG_DATE_CONFLICT"
    ok<-f=="Y"&is.finite(dt);status[ok]<-"INCLUDED"
    status[ok&dt<x$RANDDT]<-"BEFORE_RANDOMIZATION"
    status[ok&((is.finite(x$DOSE2DT)&dt<=x$DOSE2DT)|(is.finite(x$TR01SDT)&dt<=x$TR01SDT))]<-"ORIGINAL_DOSE_CHRONOLOGY_CONFLICT"
    status[ok&is.finite(x$UBLDDT)&dt<x$UBLDDT]<-"BEFORE_RECORDED_UNBLINDING"
    status[ok&is.finite(x$AP02SDT)&dt<x$AP02SDT]<-"BEFORE_OPEN_LABEL_PERIOD"
    status[ok&dt>x$end_full]<-"OUTSIDE_FULL_FOLLOWUP"
    tr<-ws_norm(x$TRT02A)
    status[ok&nzchar(tr)&tr!="MRNA-1273"]<-"OPEN_LABEL_TREATMENT_CONFLICT"
    if(d==4)status[ok&(!is.finite(x$TR02SDT)|dt<=x$TR02SDT|x$OLDOS1FL!="Y")]<-"OPEN_LABEL_DOSE_1_CHRONOLOGY_UNRESOLVED"
    status[x$arm!="Placebo"]<-"NOT_ORIGINAL_PLACEBO_ARM"
    x[[paste0("dose_",d,"_source_available")]]<-rep(all(c(flag,date)%in%source_fields),nrow(x))
    x[[paste0("dose_",d,"_status")]]<-status
    x[[paste0("dose_",d,"_eligible")]]<-status=="INCLUDED"
    # A dedicated open-label dose flag supports period classification even if
    # UBLDDT is missing; keep the missing date visible in the evidence export.
  }
  x
}
ws_crossover_evidence <- function(ctx) {
  do.call(rbind,lapply(3:4,function(d){
    flag<-if(d==3)"OLDOS1FL" else "OLDOS2FL";date<-if(d==3)"TR02SDT" else "OLDOS2DT"
    z<-ctx[ctx$arm=="Placebo",,drop=FALSE]
    data.frame(subject=z$subject,original_arm=z$arm,dose=d,source_flag=flag,flag_value=z[[flag]],source_date_field=date,
      dose_date=ws_iso(z[[date]]),original_dose_2_date=ws_iso(z$DOSE2DT),original_dose_2_available=is.finite(z$DOSE2DT),unblinding_date=ws_iso(z$UBLDDT),unblinding_date_available=is.finite(z$UBLDDT),
      open_label_period_start=ws_iso(z$AP02SDT),open_label_actual_treatment=z$TRT02A,
      eligible_full_followup=z[[paste0("dose_",d,"_eligible")]],status=z[[paste0("dose_",d,"_status")]],stringsAsFactors=FALSE)
  }))
}
ws_dose_fills <- c("#E5EBEF","#EFE8DA","#DBE9F2","#E8E1F0")
ws_dose_borders <- c("#A7B4BD","#B3A48F","#8CAABF","#AE9BC1")
ws_background_doses <- function(q,xx,show_crossover=FALSE) {
  if(show_crossover) {
    left<-c(-.46,-.22,.02,.26);right<-c(-.26,-.02,.22,.46)
    for(i in 1:4)rect(xx+left[i],0,xx+right[i],q[[paste0("dose_",i)]],col=ws_dose_fills[i],border=ws_dose_borders[i],lwd=.65)
  }else {
    rect(xx-.46,0,xx-.02,q$dose_1,col=ws_dose_fills[1],border=ws_dose_borders[1],lwd=.6)
    rect(xx+.02,0,xx+.46,q$dose_2,col=ws_dose_fills[2],border=ws_dose_borders[2],lwd=.6)
  }
}

ws_events <- function(ae,sy,mapping,origin,out,ae_timing) {
  evidence<-list();unmatched<-list()
  for(field in c("AETERM","AELLT","AEDECOD")) {
    terms<-ws_norm(ae[[field]]);rules<-mapping[mapping$match_field==field,,drop=FALSE]
    ix<-data.frame(row_id=seq_len(nrow(ae)),normalized_term=terms,stringsAsFactors=FALSE)
    hit<-merge(ix,rules,by="normalized_term",sort=FALSE)
    if(nrow(hit)) {hit$priority<-match(field,c("AETERM","AELLT","AEDECOD"));evidence[[length(evidence)+1L]]<-hit}
    bad<-nzchar(terms)&!terms %in% rules$normalized_term
    if(any(bad)) {
      u<-as.data.frame(table(terms[bad]),stringsAsFactors=FALSE);names(u)<-c("normalized_term","records")
      u$match_field<-field;unmatched[[length(unmatched)+1L]]<-u
    }
  }
  if(length(unmatched))ws_csv(do.call(rbind,unmatched),file.path(out,"unmapped_ae_terms.csv"))
  ev<-data.frame(subject=character(),source=character(),source_id=character(),source_domain=character(),source_sequence=character(),symptom_code=character(),symptom_group=character(),
    lo=numeric(),hi=numeric(),date_quality=character(),term=character(),stringsAsFactors=FALSE)
  if(length(evidence)) {
    h<-do.call(rbind,evidence);h<-h[order(h$row_id,h$symptom_code,h$priority),,drop=FALSE]
    h$selected_evidence<-!duplicated(paste(h$row_id,h$symptom_code,sep="|"))
    h$USUBJID<-ae$USUBJID[h$row_id];h$AESEQ<-ae$AESEQ[h$row_id]
    ws_csv(h,file.path(out,"ae_matching_evidence.csv"))
    # Apply evidence priority BEFORE filtering include_default: a specific
    # review-only term must not be overridden by a generic direct preferred term.
    h<-h[h$selected_evidence & h$include_default,,drop=FALSE];j<-h$row_id
    lo<-ws_dates(ae$AESTDTC,"AESTDTC",origin);ends<-ws_dates(ae$AEENDTC,"AEENDTC",origin)
    quality<-ifelse(!is.finite(lo),"UNKNOWN_START",ifelse(is.finite(ends)&ends<lo,"END_BEFORE_START",
      ifelse(is.finite(ends),"KNOWN_INTERVAL","START_ONLY_END_UNKNOWN")))
    hi<-if(ae_timing=="onset")lo else ifelse(is.finite(ends),ends,lo)
    hi[quality=="END_BEFORE_START"]<-NA_real_
    ev<-data.frame(subject=ae$subject[j],source=rep("AE",length(j)),source_id=if(length(j))paste0("AE:",ae$AESEQ[j]) else character(),
      source_domain=rep("AE",length(j)),source_sequence=ae$AESEQ[j],symptom_code=h$symptom_code,
      symptom_group=h$symptom_group,lo=lo[j],hi=hi[j],date_quality=quality[j],term=ae$AETERM[j],stringsAsFactors=FALSE)
  }
  code<-ws_norm(sy$PARAMCD);value<-ws_norm(sy$AVALC);core<-code %in% names(ws_param_map)
  inventory<-as.data.frame(table(parameter=code,result=value,source_domain=sy$SRCDOM),stringsAsFactors=FALSE);inventory<-inventory[inventory$Freq>0,]
  inventory$in_mapping<-inventory$parameter %in% names(ws_param_map)
  ws_csv(inventory,file.path(out,"symptom_parameter_inventory.csv"))
  positive<-core & value %in% c("Y","YES");j<-which(positive)
  date<-ws_dates(sy$ADT,"ADSYMP.ADT",origin)
  # The same symptom/date can have disagreeing assessments. Retain positive
  # evidence, flag the conflict, and never interpret missing/unknown as negative.
  keys<-paste(sy$subject,code,ws_iso(date),sep="|")
  conflict<-keys %in% keys[core & value %in% c("N","NO")]
  if(any(positive&conflict))ws_csv(sy[positive&conflict,,drop=FALSE],file.path(out,"contradictory_symptom_assessments.csv"))
  se<-data.frame(subject=sy$subject[j],source=rep("COVID_FORM",length(j)),source_id=if(length(j))paste0("ADSYMP:",ifelse(nzchar(sy$ASEQ[j]),sy$ASEQ[j],paste0("ROW_",j))) else character(),
    source_domain=sy$SRCDOM[j],source_sequence=sy$SRCSEQ[j],
    symptom_code=code[j],symptom_group=unname(ws_param_map[code[j]]),lo=date[j],hi=date[j],
    date_quality=ifelse(is.finite(date[j]),ifelse(conflict[j],"POSITIVE_WITH_NEGATIVE_SAME_DAY","POSITIVE_ASSESSMENT_DAY"),"UNKNOWN_ASSESSMENT_DATE"),
    term=sy$PARAM[j],stringsAsFactors=FALSE)
  ev<-unique(rbind(ev,se));ev$event_id<-seq_len(nrow(ev));ev
}
ws_expand <- function(e) {
  if(!nrow(e))return(data.frame(event_id=integer(),subject=character(),arm=character(),source=character(),symptom_group=character(),week=numeric()))
  first<-ws_week(e$lo);n<-as.integer((ws_week(e$hi)-first)/7+1)
  if(any(n<1)||sum(n)>1e7)stop("Implausible/oversized weekly expansion; inspect dates.")
  j<-rep.int(seq_len(nrow(e)),n)
  data.frame(event_id=e$event_id[j],subject=e$subject[j],arm=e$arm[j],source=e$source[j],symptom_group=e$symptom_group[j],
    week=rep.int(first,n)+7*(sequence(n)-1),stringsAsFactors=FALSE)
}
ws_window <- function(ctx,ev,window,eua_cutoff,release,out) {
  end<-switch(window,FULL=ctx$end_full,BLINDED=ctx$end_blinded,EUA=pmin(ctx$end_blinded,eua_cutoff))
  if(window!="FULL"&&max(ctx$CUTOFFDT)>as.numeric(as.Date("2020-12-18"))&&!any(ctx$has_blind_boundary))
    stop("Blinded follow-up requested but ADSL has no blinding/crossover boundary information. Supply complete ADSL or choose FULL only.")
  pop<-data.frame(subject=ctx$subject,arm=ctx$arm,lo=ctx$RANDDT,hi=end,stringsAsFactors=FALSE)
  pop<-pop[is.finite(pop$lo)&is.finite(pop$hi)&pop$hi>=pop$lo,,drop=FALSE]
  if(!nrow(pop))stop("Empty follow-up window: ",window)
  weeks<-seq(ws_week(min(pop$lo)),ws_week(max(pop$hi)),by=7)
  d<-merge(ev,pop,by="subject",suffixes=c("","_followup"),sort=FALSE)
  eligible<-is.finite(d$lo)&is.finite(d$hi)&d$hi>=d$lo&d$hi>=d$lo_followup&d$lo<=d$hi_followup
  exclusions<-d[!eligible,,drop=FALSE]
  exclusions$reason<-ifelse(!is.finite(exclusions$lo)|!is.finite(exclusions$hi),"UNKNOWN_OR_INVALID_DATE","OUTSIDE_WINDOW")
  ws_csv(exclusions,file.path(out,paste0(window,"_excluded_events.csv")))
  d<-d[eligible,,drop=FALSE];d$lo<-pmax(d$lo,d$lo_followup);d$hi<-pmin(d$hi,d$hi_followup)
  links<-ws_expand(d);links$week_start<-ws_iso(links$week)
  ws_csv(links,file.path(out,paste0(window,"_source_event_week_links.csv")))
  k<-c("subject","arm","symptom_group","week")
  if(nrow(links)) {
    a<-unique(links[links$source=="AE",k,drop=FALSE]);a$ae<-rep(TRUE,nrow(a))
    c<-unique(links[links$source=="COVID_FORM",k,drop=FALSE]);c$covid_form<-rep(TRUE,nrow(c))
    units<-merge(a,c,by=k,all=TRUE,sort=FALSE);units$ae[is.na(units$ae)]<-FALSE;units$covid_form[is.na(units$covid_form)]<-FALSE
    units$route<-ifelse(units$ae&units$covid_form,"BOTH",ifelse(units$ae,"AE_ONLY","COVID_ONLY"))
  } else units<-data.frame(subject=character(),arm=character(),symptom_group=character(),week=numeric(),ae=logical(),covid_form=logical(),route=character())
  units$week_start<-ws_iso(units$week)
  ws_csv(units,file.path(out,paste0(window,"_subject_symptom_weeks.csv")))
  denom<-ws_expand(transform(pop,event_id=seq_len(nrow(pop)),source="FOLLOWUP",symptom_group=""))
  dx<-match(denom$subject,pop$subject)
  denom$person_days<-pmin(pop$hi[dx],denom$week+6)-pmax(pop$lo[dx],denom$week)+1
  rows<-lapply(c("mRNA-1273","Placebo"),function(arm) {
    u<-units[units$arm==arm,,drop=FALSE];ct<-function(route)tabulate(match(u$week[u$route %in% route],weeks),length(weeks))
    p<-pop[pop$arm==arm,,drop=FALSE];z<-ctx[ctx$arm==arm,,drop=FALSE];endz<-end[ctx$arm==arm]
    dose<-function(flag,date) {keep<-z[[flag]]=="Y"&is.finite(z[[date]])&z[[date]]>=z$RANDDT&z[[date]]<=endz;tabulate(match(ws_week(z[[date]][keep]),weeks),length(weeks))}
    crossover<-function(d) {
      if(window!="FULL"||arm!="Placebo")return(integer(length(weeks)))
      if(!all(z[[paste0("dose_",d,"_source_available")]]))return(rep(NA_integer_,length(weeks)))
      date<-if(d==3)"TR02SDT" else "OLDOS2DT"
      keep<-z[[paste0("dose_",d,"_eligible")]]&is.finite(z[[date]])&z[[date]]<=endz
      tabulate(match(ws_week(z[[date]][keep]),weeks),length(weeks))
    }
    de<-denom[denom$arm==arm,,drop=FALSE]
    pd<-tapply(de$person_days,de$week,sum);days<-as.numeric(pd[match(as.character(weeks),names(pd))]);days[is.na(days)]<-0
    subjects<-unique(u[c("subject","week")])
    out<-data.frame(window=window,release=release,arm=arm,week_start=ws_iso(weeks),week_end=ws_iso(weeks+6),
      dose_1=dose("DOSE1FL","TR01SDT"),dose_2=dose("DOSE2FL","DOSE2DT"),dose_3=crossover(3),dose_4=crossover(4),
      covid_only=ct("COVID_ONLY"),ae_only=ct("AE_ONLY"),both=ct("BOTH"),
      subjects_with_symptoms=tabulate(match(subjects$week,weeks),length(weeks)),
      subjects_in_followup=tabulate(match(de$week,weeks),length(weeks)),person_days=days,stringsAsFactors=FALSE)
    out$total_symptom_groups<-out$covid_only+out$ae_only+out$both
    out$ae_symptom_groups<-out$ae_only+out$both
    out$pct_covid_only<-ifelse(out$total_symptom_groups>0,100*out$covid_only/out$total_symptom_groups,NA_real_)
    out$pct_ae_any<-ifelse(out$total_symptom_groups>0,100*out$ae_symptom_groups/out$total_symptom_groups,NA_real_)
    out$pct_ae_only<-ifelse(out$total_symptom_groups>0,100*out$ae_only/out$total_symptom_groups,NA_real_)
    out$pct_both<-ifelse(out$total_symptom_groups>0,100*out$both/out$total_symptom_groups,NA_real_)
    out
  })
  weekly<-do.call(rbind,rows);rownames(weekly)<-NULL
  if(sum(weekly$total_symptom_groups)!=nrow(units)||any(weekly$ae_symptom_groups>weekly$total_symptom_groups))stop("Weekly reconciliation failed.")
  ws_csv(weekly,file.path(out,paste0(window,"_weekly_counts.csv")))
  list(weekly=weekly,units=units,population=pop)
}
# Keep every percentage readable, including short/zero-height stack segments.
ws_pct_label <- function(x) ifelse(is.finite(x),sprintf("%.1f%%",x),"")
ws_label_positions <- function(target,height,lower,upper,gap) {
  if(!length(target))return(numeric())
  ord<-order(target);t<-target[ord];h<-height[ord]/2;y<-t
  if(sum(height)+gap*(length(y)-1)>upper-lower)stop("Not enough vertical space for percentage labels; enlarge the output page.")
  y[1]<-max(y[1],lower+h[1])
  if(length(y)>1)for(i in 2:length(y))y[i]<-max(y[i],y[i-1]+h[i-1]+h[i]+gap)
  y[length(y)]<-min(y[length(y)],upper-h[length(y)])
  if(length(y)>1)for(i in (length(y)-1):1)y[i]<-min(y[i],y[i+1]-h[i+1]-h[i]-gap)
  if(y[1]<lower+h[1])y<-y+(lower+h[1]-y[1])
  ans<-numeric(length(y));ans[ord]<-y;ans
}
ws_plot_arm <- function(w,ctx,window,release,ae_timing,arm) {
  old<-par(no.readonly=TRUE);on.exit({layout(1);par(old)})
  layout(matrix(1:3,ncol=1),heights=c(.92,6.0,.83))
  par(mar=c(0,.8,0,.8),family="sans");plot.new();plot.window(c(0,1),c(0,1))
  text(.02,.85,paste0("Moderna P301 | ",arm," | Weekly symptom documentation"),adj=0,font=2,cex=1.17,col="#183746")
  label<-switch(window,BLINDED="Individual blinded follow-up",EUA="Through the EUA data cutoff",FULL="Full follow-up by original randomized arm")
  text(.02,.57,paste0(label," | ",release," | Source cutoff ",ws_iso(ctx$CUTOFFDT[1])),adj=0,cex=.83,col="#405865")
  fill<-c(ws_dose_fills[1:2],"#367E8A","#A65B31","#785E91",NA)
  crossover<-window=="FULL"&&arm=="Placebo"
  nd<-if(crossover)4L else 2L
  legend(.02,.37,legend=if(crossover)c("Dose 1","Dose 2","Dose 3 (open-label 1)","Dose 4 (open-label 2)") else c("Dose 1","Dose 2"),
    fill=ws_dose_fills[seq_len(nd)],border=ws_dose_borders[seq_len(nd)],ncol=nd,bty="n",cex=.73,x.intersp=.65,y.intersp=.8)
  legend(.02,.11,legend=c("Symptom table only","AE only","Both sources","AE share (right axis)"),
    fill=c(fill[3:5],NA),border=NA,col=c(rep(NA,3),"#222838"),lty=c(rep(NA,3),1),lwd=c(rep(NA,3),1.8),
    pch=c(rep(NA,3),16),pt.cex=.6,ncol=4,bty="n",cex=.73,x.intersp=.65,y.intersp=.8)
  maxn<-max(c(1,w$dose_1,w$dose_2,w$dose_3,w$dose_4,w$total_symptom_groups),na.rm=TRUE);top<-max(pretty(c(0,maxn*1.08),n=5))
  q<-w[w$arm==arm,,drop=FALSE];q<-q[order(q$week_start),,drop=FALSE];n<-nrow(q);xx<-seq_len(n)
  if(!n)stop("No weekly rows for arm ",arm)
  par(mar=c(3.3,4.8,2.1,4.8),mgp=c(2.7,.65,0),las=1,bty="n",xaxs="i",yaxs="i")
  plot(c(.35,n+.7),c(0,top),type="n",axes=FALSE,xlab="",ylab="Weekly count",cex.lab=.88)
  ticks<-pretty(c(0,top),n=5);abline(h=ticks,col="#DFE5E8",lwd=.55)
  ws_background_doses(q,xx,crossover)
  base<-rep(0,n);centers<-matrix(0,nrow=n,ncol=3)
  for(i in 1:3) {
    vals<-q[[c("covid_only","ae_only","both")[i]]]
    rect(xx-.40,base,xx+.18,base+vals,col=fill[i+2],border=NA)
    centers[,i]<-base+vals/2;base<-base+vals
  }
  axis(2,at=ticks,labels=ws_n(ticks),tick=FALSE,cex.axis=.8)
  idx<-seq_len(n);dates<-as.Date(q$week_start[idx])
  labs<-paste(format(dates,"%d"),month.abb[as.integer(format(dates,"%m"))],format(dates,"%y"))
  axis(1,at=xx[idx],labels=FALSE,tck=-.012,col="#73818A")
  text(xx[idx],-.025*top,labels=labs,srt=55,adj=1,xpd=NA,cex=.67,col="#3C515E")
  lines(xx,q$pct_ae_any/100*top,col="#222838",lwd=1.6)
  points(xx,q$pct_ae_any/100*top,pch=21,bg="white",col="#222838",cex=.52,lwd=.85)
  axis(4,at=(0:4)*top/4,labels=paste0((0:4)*25,"%"),tick=FALSE,cex.axis=.8,col.axis="#222838")
  mtext("Symptom groups recorded as AEs (%)",side=4,line=2.9,las=0,cex=.8)
  title(main=paste0("Randomized N = ",ws_n(sum(ctx$arm==arm))," | ",ws_n(sum(q$total_symptom_groups))," symptom-group weeks"),
    adj=0,line=.55,cex.main=.91,col.main="#183746")
  # Rotated labels fit one lane for stack shares and a second for the AE line.
  # Percentages always use the weekly union as denominator, not arm population.
  cex_label<-1.05 # about 8.3 pt after base R layout scaling
  inch_x<-diff(par("usr")[1:2])/par("pin")[1];inch_y<-diff(par("usr")[3:4])/par("pin")[2]
  box_size<-function(labels)list(width=(strheight(labels,units="inches",cex=cex_label)+.025)*inch_x,
    height=(strwidth(labels,units="inches",cex=cex_label)+.035)*inch_y)
  draw_badge<-function(x,y,lab,col,text_col,border) {
    sz<-box_size(lab);rect(x-sz$width/2,y-sz$height/2,x+sz$width/2,y+sz$height/2,col=col,border=border,lwd=.45)
    text(x,y,lab,srt=90,cex=cex_label,col=text_col)
  }
  for(j in seq_len(n))if(q$total_symptom_groups[j]>0) {
    counts<-unlist(q[j,c("covid_only","ae_only","both")],use.names=FALSE)
    labs<-ws_pct_label(100*counts/q$total_symptom_groups[j]);sz<-box_size(labs)
    yy<-ws_label_positions(centers[j,],sz$height,0,top,.022*inch_y)
    for(i in 1:3) {
      if(abs(yy[i]-centers[j,i])>.004*top) {
        # Bent colored leaders link moved labels to the original segment center.
        segments(xx[j]-.11,centers[j,i],xx[j]-.43,centers[j,i],col=fill[i+2],lwd=.7)
        segments(xx[j]-.43,centers[j,i],xx[j]-.43,yy[i],col=fill[i+2],lwd=.7)
        segments(xx[j]-.43,yy[i],xx[j]-.11,yy[i],col=fill[i+2],lwd=.7)
      }
      inside<-yy[i]-sz$height[i]/2>=centers[j,i]-counts[i]/2 && yy[i]+sz$height[i]/2<=centers[j,i]+counts[i]/2
      # An offset badge stays white: its background must not visually extend
      # the colored count bar beyond the actual observed count.
      draw_badge(xx[j]-.11,yy[i],labs[i],if(inside)fill[i+2] else "white",
        if(inside)"white" else fill[i+2],if(inside)"white" else fill[i+2])
    }
  }
  for(j in which(is.finite(q$pct_ae_any))) {
    lab<-ws_pct_label(q$pct_ae_any[j]);sz<-box_size(lab);actual<-q$pct_ae_any[j]/100*top
    yy<-min(top-sz$height/2,max(sz$height/2,actual));label_x<-xx[j]+.37
    segments(xx[j],actual,label_x,yy,col="#222838",lwd=.6)
    draw_badge(label_x,yy,lab,"white","#171D28","#A9B0B7")
  }
  par(mar=c(0,.8,0,.8));plot.new();plot.window(c(0,1),c(0,1))
  notes<-c("Colored labels: % of that week's symptom-group total. Black labels: AE share = (AE only + both) / all three sources.",
    "Small segments use offset labels with connecting lines. Values rounded to 0.1%; empty weeks have no percentage. One arm per page.",
    paste0("One participant x symptom group x week; Monday start. AE timing: ",if(ae_timing=="interval")"known active intervals" else "onset week only",". Symptom table: positive assessments."),
    if(window=="FULL")if(anyNA(q[c("dose_3","dose_4")]))"Original arms retained. Open-label dose counts unavailable in this ADSL; see dose evidence CSV." else "Original arms retained. Placebo chart includes open-label doses 3/4; symptoms are not confirmed infections."
      else "Original doses 1/2 only. Counts depend on documentation and follow-up; these symptoms are not confirmed infections.")
  text(.02,seq(.87,.12,length.out=4),notes,adj=0,cex=.70,col="#405865")
}
ws_plot <- function(w,ctx,window,release,ae_timing,arm=NULL) {
  arms<-if(is.null(arm))c("mRNA-1273","Placebo") else match.arg(arm,c("mRNA-1273","Placebo"))
  for(a in arms)ws_plot_arm(w,ctx,window,release,ae_timing,a)
}
ws_export_plot <- function(draw,path,dpi) {
  if(capabilities("cairo"))grDevices::cairo_pdf(paste0(path,".pdf"),width=297/25.4,height=210/25.4,family="sans",onefile=TRUE)
  else grDevices::pdf(paste0(path,".pdf"),width=297/25.4,height=210/25.4,paper="special",family="Helvetica",onefile=TRUE)
  tryCatch(draw(),finally=grDevices::dev.off())
  for(arm in c("mRNA-1273","Placebo")) {
    name<-paste0(path,if(arm=="mRNA-1273")"_mrna1273" else "_placebo")
    grDevices::png(paste0(name,".png"),width=297/25.4,height=210/25.4,units="in",res=dpi,type=if(capabilities("cairo"))"cairo" else getOption("bitmapType"))
    tryCatch(draw(arm),finally=grDevices::dev.off())
    if(capabilities("cairo")) {grDevices::svg(paste0(name,".svg"),width=297/25.4,height=210/25.4);tryCatch(draw(arm),finally=grDevices::dev.off())}
  }
}

analyse_weekly_symptoms_ae <- function(root="",release="auto",adsl_file="",ae_file="",symptoms_file="",mapping_file="",
  windows=c("BLINDED","EUA","FULL"),eua_cutoff="2020-11-25",ae_timing="interval",output_dir="",encoding="auto",
  numeric_date_origin="",min_randomized=30000L,png_dpi=300L,show_chart=FALSE) {
  if(!release %in% c("auto","EUA","Prod"))stop("Release must be auto, EUA or Prod.")
  if(!length(windows)||any(!windows %in% c("BLINDED","EUA","FULL")))stop("Invalid windows.")
  windows<-unique(windows)
  if(!ae_timing %in% c("interval","onset"))stop("AE timing must be interval or onset.")
  if(length(png_dpi)!=1L||!is.finite(png_dpi)||png_dpi<72||png_dpi>1200)stop("PNG DPI must be between 72 and 1200.")
  cutoff<-ws_dates(eua_cutoff,"EUA cutoff","",TRUE);if(length(cutoff)!=1L||!is.finite(cutoff))stop("Complete EUA cutoff required.")
  root<-ws_root(root)
  out<-if(nzchar(output_dir))ws_path(root,output_dir) else file.path(root,"analysis",paste0("weekly_symptoms_ae_",format(Sys.time(),"%Y%m%d_%H%M%S")))
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  if(release=="auto") {
    release<-if(length(ws_candidates(root,"Prod","ADSL")))"Prod" else "EUA"
    if(nzchar(adsl_file)) {
      p<-tolower(gsub("\\\\","/",ws_path(root,adsl_file)))
      if(grepl("/prod_data/",p,fixed=TRUE))release<-"Prod" else if(grepl("/eua_data/",p,fixed=TRUE))release<-"EUA"
      else stop("For an explicit ADSL outside eua_data/prod_data, set PHMPT_RELEASE to EUA or Prod.")
    }
  }
  a<-ws_pick(root,release,"ADSL",adsl_file,out,encoding,numeric_date_origin,min_randomized)
  ctx<-ws_context(a$x,numeric_date_origin);submission<-ws_submission(a$manifest$path)
  dose_evidence<-ws_crossover_evidence(ctx)
  ws_csv(dose_evidence,file.path(out,"placebo_open_label_dose_evidence.csv"))
  ws_csv(dose_evidence[!dose_evidence$status%in%c("INCLUDED","NO_RECORDED_DOSE"),,drop=FALSE],file.path(out,"placebo_open_label_dose_issues.csv"))
  ds<-as.data.frame(table(dose=dose_evidence$dose,status=dose_evidence$status),stringsAsFactors=FALSE);ds<-ds[ds$Freq>0,,drop=FALSE]
  ws_csv(ds,file.path(out,"placebo_open_label_dose_status_counts.csv"))
  b<-ws_pick(root,release,"AE",ae_file,out,encoding,numeric_date_origin,min_randomized,submission)
  c<-ws_pick(root,release,"ADSYMP",symptoms_file,out,encoding,numeric_date_origin,min_randomized,submission)
  manifest<-rbind(a$manifest,b$manifest,c$manifest);ws_csv(manifest,file.path(out,"selected_inputs.csv"))
  # If ADSYMP carries a cutoff, require agreement with the selected ADSL.
  cd<-ws_dates(c$x$CUTOFFDT,"ADSYMP.CUTOFFDT",numeric_date_origin,TRUE)
  if(any(is.finite(cd)&cd!=ctx$CUTOFFDT[1]))stop("ADSYMP and ADSL cutoffs disagree.")
  mapping<-ws_mapping(root,mapping_file,out)
  ev<-ws_events(b$x,c$x,mapping,numeric_date_origin,out,ae_timing)
  ev$start_date<-ws_iso(ev$lo);ev$end_date<-ws_iso(ev$hi)
  ev$in_randomized_population<-ev$subject %in% ctx$subject
  ev$arm<-ctx$arm[match(ev$subject,ctx$subject)]
  ws_csv(ev,file.path(out,"symptom_event_evidence.csv"))
  ws_csv(ev[!is.finite(ev$lo)|!is.finite(ev$hi)|ev$hi<ev$lo,,drop=FALSE],file.path(out,"undated_or_invalid_symptom_events.csv"))
  ws_csv(ev[!ev$in_randomized_population,,drop=FALSE],file.path(out,"events_outside_randomized_population.csv"))
  # Arm is attached inside each window to avoid duplicate arm columns on merge.
  ev$arm<-NULL
  context<-ctx;for(n in c("RANDDT","TR01SDT","DOSE2DT","OLDOS2DT","CUTOFFDT","AP01EDT","UBLDDT","AP02SDT","TR02SDT","EUAVACDT","EOSDT","DTHDT","end_full","end_blinded"))context[[n]]<-ws_iso(context[[n]])
  ws_csv(context,file.path(out,"subject_followup_and_doses.csv"))
  results<-list()
  for(window in windows) {
    message("Aggregating ",window," ...")
    r<-ws_window(ctx,ev,window,cutoff,release,out);results[[window]]<-r$weekly
    draw<-local({w<-r$weekly;win<-window;function(arm=NULL)ws_plot(w,ctx,win,release,ae_timing,arm)})
    ws_export_plot(draw,file.path(out,paste0("weekly_symptoms_ae_",tolower(window))),png_dpi)
  }
  all<-do.call(rbind,results);rownames(all)<-NULL;ws_csv(all,file.path(out,"weekly_counts_all_windows.csv"))
  totals<-do.call(rbind,lapply(split(all,paste(all$window,all$arm)),function(z) {
    den<-sum(z$total_symptom_groups);data.frame(window=z$window[1],arm=z$arm[1],dose_1=sum(z$dose_1),dose_2=sum(z$dose_2),dose_3=sum(z$dose_3),dose_4=sum(z$dose_4),
      covid_only=sum(z$covid_only),ae_only=sum(z$ae_only),both=sum(z$both),symptom_group_weeks=den,
      pct_ae_any=if(den>0)100*sum(z$ae_symptom_groups)/den else NA_real_)
  }));rownames(totals)<-NULL;ws_csv(totals,file.path(out,"window_totals.csv"))
  notes<-c("MODERNA P301 - WEEKLY DOSES, SYMPTOM GROUPS AND AE RECORDING SHARE - v1.1.0", "",
    paste("Release:",release,"; source cutoff:",ws_iso(ctx$CUTOFFDT[1]),"; randomized population:",nrow(ctx)),
    paste("AE timing:",ae_timing,"; windows:",paste(windows,collapse=", ")),"",
    "COUNTING UNIT: a distinct participant x symptom group x Monday-Sunday week. Repeated records, days and related parameters within that group/week count once. One participant with three groups contributes three units. This is not a count of visits, raw AE records, clinical episodes or infections. Units in different weeks are not unique participants.",
    "The 17 symptom parameters collapse to the same 12 groups used in the agreed AE mapping: breathing; muscle/body aches; smell/taste; nasal symptoms; and nausea/vomiting each share one group. Other groups remain separate. AE labels map by exact field plus normalized term, without automatic expansion of diagnoses.",
    "COLUMNS: symptom table only (COVID_FORM in CSV evidence), AE only and both sources are disjoint and sum to the union. ADSYMP includes COVID-form symptoms and qualifying vital-sign evidence, such as fever. BOTH means evidence in both sources for the same participant/group/week; it does not prove that the records refer to the same illness. Repeated separate episodes within a week are deliberately collapsed.",
    "LABELS: all three weekly source shares and every defined AE-share point are displayed to one decimal place. Small stack segments have offset labels connected to the original segment. Independently rounded source shares may sum to 99.9% or 100.1%. Zero source shares are shown when the weekly total is positive; empty weeks have no percentage.",
    "LINE: 100 * (AE_ONLY + BOTH) / (COVID_ONLY + AE_ONLY + BOTH). An empty denominator is NA, so the line breaks. This is a documentation-route proportion, not the probability that a reported symptom was reclassified, not an adverse-event incidence rate and not a treatment effect. pct_ae_only is also exported separately.",
    "DATES: COVID_FORM evidence uses ADSYMP.ADT and AVALC Y/YES. ADT is the assessment/analysis date, not necessarily onset. Only mapped checklist parameters are counted. Contradictory positive/negative same-day records are logged; positive evidence is retained. Negative, unknown, unperformed or missing assessments do not count as positive.",
    "AE interval mode uses AESTDTC through a complete AEENDTC, inclusive, clipped to follow-up. Missing/partial end dates contribute only the onset day; ongoing duration is not imputed. Onset mode uses AESTDTC only. Missing/partial starts and end-before-start records are excluded and logged. Different collection frequencies and incomplete durations can affect the line.",
    "FOLLOW-UP: RANDFL Y, planned original arm TRT01P, from RANDDT. FULL ends at the earliest CUTOFFDT/EOSDT/DTHDT. BLINDED additionally ends at AP01EDT or the day before UBLDDT/AP02SDT/TR02SDT/EUAVACDT, whichever comes first. Missing individual boundaries fall back to the full endpoint and are exposed in the context CSV. EUA uses the blinded endpoint capped at the configured EUA date.",
    "FULL retains original randomization labels after unblinding/crossover; placebo-assigned participants may subsequently receive vaccine. It is a descriptive follow-up view. The blinded comparison also has changing follow-up: weekly subjects_in_followup and person_days are exported. No causal conclusion is established by these curves.",
    "DOSES: background columns show original DOSE1FL/TR01SDT and DOSE2FL/DOSE2DT administrations in the same population/window. Placebo injections are included. FULL also shows doses 3/4 for the original placebo arm, using OLDOS1FL/TR02SDT and OLDOS2FL/OLDOS2DT (open-label injections 1/2). Recorded flags and complete dates must agree and fall within follow-up, after any recorded original doses and on/after any recorded unblinding/open-label-period start. Dose 3/4 labels refer to the scheduled open-label injections 1/2; an incomplete original series is exposed in the evidence CSV, not grounds to drop an otherwise documented open-label dose. These doses are zero in BLINDED/EUA windows. In FULL, unavailable source fields produce NA, not a zero-dose count. Missing source fields, inconsistent records and exclusions are exported in the open-label dose evidence/status/issue CSVs. Missing UBLDDT remains explicit; a dedicated open-label dose flag is not replaced with an imputed unblinding date. EUAVACDT/outside-study vaccination dates are not counted as administered study doses. mRNA-1273-assigned open-label doses are outside this requested overlay. The left axis is unscaled weekly counts shared between arms; the right axis is a fixed 0-100% for the AE share. Bar and line heights use different units.",
    "MAPPING: an external assets/ae_covid_symptom_mapping.csv takes precedence; otherwise the embedded field-specific v1.0.0 asset is used. Specific AETERM then AELLT evidence takes precedence over AEDECOD for each symptom parameter, before include_default is applied. Review-only matches are excluded. Equivalent normalized mapping rows are consolidated only when their classification, inclusion and qualification decisions agree. Original mapping rows, duplicate aliases and genuine conflicts are exported separately. All matching evidence, excluded dates and unmapped terms are exported for review.",
    "SCOPE: COVID forms plus mapped raw AE records, regardless of PCR/serology, severity or investigator causal attribution. Solicited reactogenicity diaries are a separate collection source and are not included. No duration/fever threshold for protocol-trigger eligibility is imposed. Compatible symptoms do not establish COVID-19 or a missed PCR test.",
    "SOURCE CONTROL: one selected snapshot per domain/release; duplicate files are not pooled. Matching submission identifiers are preferred when available. An older ADSL CSV projection lacking open-label dose fields is superseded only when its common fields agree with a complete same-cutoff ADSL. Conflicting snapshots require an explicit path. Selected file paths and MD5 checksums are exported; inspect them before publication. Raw AE and ADSYMP often have no cutoff field, so submission/date coverage should be reviewed using the source documents.",
    "OUTPUTS: a two-page print-ready PDF per window (mRNA-1273 then placebo, matching scales), plus separate _mrna1273 and _placebo PNG/SVG images (SVG when Cairo is available); weekly counts; participant-group-week membership; event-week links; original event evidence; mapping evidence; follow-up context; source discovery and methods.","",
    "Source metadata: production P301 ADaM define.xml ADSYMP.ADT, AVALC and ADSL dosing/window variables; the previously reviewed analyst AE-to-symptom mapping.")
  writeLines(notes,file.path(out,"README.txt"));capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
  ws_export_audit_assets(root,out,all,ctx,release,ae_timing)
  if(show_chart)ws_plot(results[[1]],ctx,windows[1],release,ae_timing)
  message("Reports and CSV assets: ",normalizePath(out,winslash="/"));print(totals,row.names=FALSE)
  invisible(list(output_dir=out,weekly=all,totals=totals,context=ctx))
}

# ---- EMBEDDED STAND-ALONE UTILITIES AND REFERENCE MAPPING ----
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


ws_embedded_mapping <- function() read.csv(text=paste(c(
  "dictionary_version,match_field,normalized_term,symptom_code,symptom_label,symptom_group,match_class,include_default",
  "1.0.0,AEDECOD,POST PROCEDURAL FEVER,FEVER,Fever,FEVER,DIRECT,TRUE",
  "1.0.0,AEDECOD,PYREXIA,FEVER,Fever,FEVER,DIRECT,TRUE",
  "1.0.0,AELLT,FEBRILE REACTION,FEVER,Fever,FEVER,DIRECT,TRUE",
  "1.0.0,AELLT,FEVER,FEVER,Fever,FEVER,DIRECT,TRUE",
  "1.0.0,AELLT,FEVER CHILLS,FEVER,Fever,FEVER,DIRECT,TRUE",
  "1.0.0,AELLT,FEVER OF UNKNOWN ORIGIN,FEVER,Fever,FEVER,DIRECT,TRUE",
  "1.0.0,AELLT,INTERMITTENT FEVER,FEVER,Fever,FEVER,DIRECT,TRUE",
  "1.0.0,AELLT,POST PROCEDURAL FEVER,FEVER,Fever,FEVER,DIRECT,TRUE",
  "1.0.0,AELLT,POSTOPERATIVE FEVER,FEVER,Fever,FEVER,DIRECT,TRUE",
  "1.0.0,AELLT,PYREXIA,FEVER,Fever,FEVER,DIRECT,TRUE",
  "1.0.0,AETERM,FEVER,FEVER,Fever,FEVER,DIRECT,TRUE",
  "1.0.0,AETERM,FEVER OF UNKNOWN ORIGIN,FEVER,Fever,FEVER,DIRECT,TRUE",
  "1.0.0,AETERM,INTERMITTENT FEVER,FEVER,Fever,FEVER,DIRECT,TRUE",
  "1.0.0,AETERM,PYREXIA,FEVER,Fever,FEVER,DIRECT,TRUE",
  "1.0.0,AEDECOD,BODY TEMPERATURE INCREASED,FEVER,Fever,FEVER,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,BODY TEMPERATURE INCREASED,FEVER,Fever,FEVER,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,FEVERISH,FEVER,Fever,FEVER,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,HIGH TEMPERATURE,FEVER,Fever,FEVER,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,LOW GRADE FEVER,FEVER,Fever,FEVER,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,SLIGHT TEMPERATURE,FEVER,Fever,FEVER,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,FEVERISH,FEVER,Fever,FEVER,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,LOW GRADE FEVER,FEVER,Fever,FEVER,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,CHILLS,CHILLS,Chills,CHILLS,DIRECT,TRUE",
  "1.0.0,AELLT,CHILLS,CHILLS,Chills,CHILLS,DIRECT,TRUE",
  "1.0.0,AELLT,FEVER CHILLS,CHILLS,Chills,CHILLS,DIRECT,TRUE",
  "1.0.0,AELLT,RIGORS,CHILLS,Chills,CHILLS,DIRECT,TRUE",
  "1.0.0,AETERM,CHILL,CHILLS,Chills,CHILLS,DIRECT,TRUE",
  "1.0.0,AETERM,CHILLS,CHILLS,Chills,CHILLS,DIRECT,TRUE",
  "1.0.0,AETERM,RIGORS,CHILLS,Chills,CHILLS,DIRECT,TRUE",
  "1.0.0,AEDECOD,ALLERGIC COUGH,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AEDECOD,COUGH,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AEDECOD,PRODUCTIVE COUGH,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AELLT,ALLERGIC COUGH,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AELLT,CHRONIC COUGH,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AELLT,COUGH,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AELLT,COUGH AGGRAVATED,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AELLT,COUGH INCREASED,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AELLT,COUGH NONPRODUCTIVE,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AELLT,COUGHING,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AELLT,DRUG-INDUCED COUGH,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AELLT,DRY COUGH,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AELLT,NOCTURNAL COUGH,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AELLT,PERSISTENT COUGH,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AELLT,PRODUCTIVE COUGH,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AELLT,SMOKER'S COUGH,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AETERM,CHRONIC COUGH,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AETERM,COUGH,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AETERM,COUGHING,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AETERM,DRY COUGH,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AETERM,PRODUCTIVE COUGH,COUGH,Cough,COUGH,DIRECT,TRUE",
  "1.0.0,AEDECOD,DYSPNOEA,SHBREATH,Shortness of breath,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AEDECOD,DYSPNOEA EXERTIONAL,SHBREATH,Shortness of breath,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AELLT,DYSPNEA,SHBREATH,Shortness of breath,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AELLT,DYSPNEA EXACERBATED,SHBREATH,Shortness of breath,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AELLT,DYSPNEA EXERTIONAL,SHBREATH,Shortness of breath,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AELLT,EXERTIONAL DYSPNEA,SHBREATH,Shortness of breath,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AELLT,INCREASED SHORTNESS OF BREATH,SHBREATH,Shortness of breath,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AELLT,SHORT OF BREATH,SHBREATH,Shortness of breath,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AELLT,SHORTNESS OF BREATH,SHBREATH,Shortness of breath,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AETERM,DYSPNEA,SHBREATH,Shortness of breath,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AETERM,EXERTIONAL DYSPNEA,SHBREATH,Shortness of breath,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AETERM,INCREASED SHORTNESS OF BREATH,SHBREATH,Shortness of breath,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AETERM,SHORT OF BREATH,SHBREATH,Shortness of breath,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AETERM,SHORTNESS OF BREATH,SHBREATH,Shortness of breath,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AELLT,BREATHING DIFFICULT,DIBREATH,Difficulty breathing,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AELLT,DIFFICULTY BREATHING,DIBREATH,Difficulty breathing,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AELLT,LABORED BREATHING,DIBREATH,Difficulty breathing,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AETERM,DIFFICULT BREATHING,DIBREATH,Difficulty breathing,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AETERM,DIFFICULTY BREATHING,DIBREATH,Difficulty breathing,DYSPNEA,DIRECT,TRUE",
  "1.0.0,AEDECOD,RESPIRATORY DISTRESS,DIBREATH,Difficulty breathing,DYSPNEA,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,RESPIRATORY DISTRESS,DIBREATH,Difficulty breathing,DYSPNEA,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,RESPIRATORY DISTRESS,DIBREATH,Difficulty breathing,DYSPNEA,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,FATIGUE,FATIGUE,Fatigue,FATIGUE,DIRECT,TRUE",
  "1.0.0,AELLT,FATIGUE,FATIGUE,Fatigue,FATIGUE,DIRECT,TRUE",
  "1.0.0,AELLT,FATIGUE AGGRAVATED,FATIGUE,Fatigue,FATIGUE,DIRECT,TRUE",
  "1.0.0,AELLT,FATIGUE EXTREME,FATIGUE,Fatigue,FATIGUE,DIRECT,TRUE",
  "1.0.0,AELLT,TIREDNESS,FATIGUE,Fatigue,FATIGUE,DIRECT,TRUE",
  "1.0.0,AETERM,FATIGUE,FATIGUE,Fatigue,FATIGUE,DIRECT,TRUE",
  "1.0.0,AETERM,TIREDNESS,FATIGUE,Fatigue,FATIGUE,DIRECT,TRUE",
  "1.0.0,AEDECOD,ASTHENIA,FATIGUE,Fatigue,FATIGUE,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,MALAISE,FATIGUE,Fatigue,FATIGUE,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,MENTAL FATIGUE,FATIGUE,Fatigue,FATIGUE,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,MUSCLE FATIGUE,FATIGUE,Fatigue,FATIGUE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,ASTHENIA,FATIGUE,Fatigue,FATIGUE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,MALAISE,FATIGUE,Fatigue,FATIGUE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,MENTAL FATIGUE,FATIGUE,Fatigue,FATIGUE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,MUSCLE FATIGUE,FATIGUE,Fatigue,FATIGUE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,WEAKNESS,FATIGUE,Fatigue,FATIGUE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,WEAKNESS GENERALIZED,FATIGUE,Fatigue,FATIGUE,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,GENERALIZED WEAKNESS,FATIGUE,Fatigue,FATIGUE,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,MALAISE,FATIGUE,Fatigue,FATIGUE,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,MENTAL FATIGUE,FATIGUE,Fatigue,FATIGUE,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,MUSCLE FATIGUE,FATIGUE,Fatigue,FATIGUE,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,WEAKNESS,FATIGUE,Fatigue,FATIGUE,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,MYALGIA,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,GENERALISED MUSCLE ACHES,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,GENERALIZED MUSCLE ACHES,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,LOCALISED MUSCLE PAIN,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,LOCALIZED MUSCLE PAIN,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,MUSCLE ACHE,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,MUSCLE PAIN,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,MUSCLE SORENESS,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,MUSCULAR PAIN,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,MUSCULAR PAINS,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,MYALGIA,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,MYALGIA AGGRAVATED,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,MYALGIA OF LOWER EXTREMITIES,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,MYALGIA UPPER EXTREMITIES,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,PAIN MUSCLE,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,POLYMYALGIA,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,TENDERNESS MUSCLE,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AETERM,GENERALIZED MUSCLE ACHES,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AETERM,MUSCLE ACHES,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AETERM,MUSCLE ACHE,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AETERM,MUSCLE ACHES,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AETERM,MUSCLE PAIN,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AETERM,MUSCLE SORENESS,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AETERM,MYALGIA,MYALGIA,Muscle aches (myalgia),MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,GENERAL BODY PAIN,BODYACHE,Body aches,MYALGIA,DIRECT,TRUE",
  "1.0.0,AELLT,GENERALIZED ACHING,BODYACHE,Body aches,MYALGIA,DIRECT,TRUE",
  "1.0.0,AETERM,(WHOLE BODY) BODY ACHES,BODYACHE,Body aches,MYALGIA,DIRECT,TRUE",
  "1.0.0,AETERM,ALL OVER BODY ACHES,BODYACHE,Body aches,MYALGIA,DIRECT,TRUE",
  "1.0.0,AETERM,BODY ACHE,BODYACHE,Body aches,MYALGIA,DIRECT,TRUE",
  "1.0.0,AETERM,BODY ACHES,BODYACHE,Body aches,MYALGIA,DIRECT,TRUE",
  "1.0.0,AETERM,GENERAL BODY ACHES,BODYACHE,Body aches,MYALGIA,DIRECT,TRUE",
  "1.0.0,AETERM,GENERALIZED PAIN,BODYACHE,Body aches,MYALGIA,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,GENERALIZED PAIN,BODYACHE,Body aches,MYALGIA,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,CLUSTER HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AEDECOD,EXERTIONAL HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AEDECOD,HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AEDECOD,NEW DAILY PERSISTENT HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AEDECOD,POST-TRAUMATIC HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AEDECOD,PRIMARY COUGH HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AEDECOD,PRIMARY HEADACHE ASSOCIATED WITH SEXUAL ACTIVITY,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AEDECOD,PROCEDURAL HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AEDECOD,SINUS HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AEDECOD,TENSION HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,CLUSTER HEADACHES,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,COITAL HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,DRUG-INDUCED HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,EXERTIONAL HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,FOREHEAD HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,FREQUENT HEADACHES,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,FRONTAL HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,HEAD PAIN,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,HEADACHE (EXCL MIGRAINE) AGGRAVATED,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,HEADACHE AGGRAVATED,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,HEADACHE DULL,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,HEADACHE POST-TRAUMATIC,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,HEADACHE POSTOPERATIVE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,HEADACHE RECURRENT,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,HEADACHE SINUS,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,HEADACHE TEMPORAL,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,HEADACHE TENSION,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,HEADACHE UNILATERAL,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,INTERMITTENT HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,MIGRAINE HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,NEW DAILY PERSISTENT HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,OCCIPITAL HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,PAIN HEAD,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,POST PROCEDURAL HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,POST-TRAUMATIC HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,PRIMARY COUGH HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,SINUS HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,TEMPORAL HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,TENSION HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AELLT,TENSION HEADACHES,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AETERM,CLUSTER HEADACHES,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AETERM,FRONTAL HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AETERM,HEAD PAIN,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AETERM,HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AETERM,INTERMITTENT HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AETERM,MIGRAINE HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AETERM,OCCIPITAL HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AETERM,SINUS HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AETERM,TENSION HEADACHE,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AETERM,TENSION HEADACHES,HEADACHE,Headache,HEADACHE,DIRECT,TRUE",
  "1.0.0,AEDECOD,MIGRAINE,HEADACHE,Headache,HEADACHE,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,MIGRAINE WITH AURA,HEADACHE,Headache,HEADACHE,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,MIGRAINE WITHOUT AURA,HEADACHE,Headache,HEADACHE,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,VESTIBULAR MIGRAINE,HEADACHE,Headache,HEADACHE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,COMPLICATED MIGRAINE,HEADACHE,Headache,HEADACHE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,MIGRAINE,HEADACHE,Headache,HEADACHE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,MIGRAINE AGGRAVATED,HEADACHE,Headache,HEADACHE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,MIGRAINE AURA,HEADACHE,Headache,HEADACHE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,MIGRAINE WITH AURA,HEADACHE,Headache,HEADACHE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,MIGRAINE WITHOUT AURA,HEADACHE,Headache,HEADACHE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,OPHTHALMIC MIGRAINE,HEADACHE,Headache,HEADACHE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,VESTIBULAR MIGRAINE,HEADACHE,Headache,HEADACHE,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,COMPLICATED MIGRAINE,HEADACHE,Headache,HEADACHE,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,MIGRAINE,HEADACHE,Headache,HEADACHE,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,MIGRAINE AURA,HEADACHE,Headache,HEADACHE,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,ANOSMIA,NLSMELL,New loss of smell,SMELL_TASTE,DIRECT,TRUE",
  "1.0.0,AELLT,ANOSMIA,NLSMELL,New loss of smell,SMELL_TASTE,DIRECT,TRUE",
  "1.0.0,AELLT,LOSS OF SMELL,NLSMELL,New loss of smell,SMELL_TASTE,DIRECT,TRUE",
  "1.0.0,AETERM,ANOSMIA,NLSMELL,New loss of smell,SMELL_TASTE,DIRECT,TRUE",
  "1.0.0,AETERM,LOSS OF SMELL,NLSMELL,New loss of smell,SMELL_TASTE,DIRECT,TRUE",
  "1.0.0,AETERM,NEW LOSS OF SMELL,NLSMELL,New loss of smell,SMELL_TASTE,DIRECT,TRUE",
  "1.0.0,AEDECOD,HYPOSMIA,NLSMELL,New loss of smell,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,PAROSMIA,NLSMELL,New loss of smell,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,ALTERED SMELL SENSATION,NLSMELL,New loss of smell,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,DIMINISHED SENSE OF SMELL,NLSMELL,New loss of smell,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,HYPOSMIA,NLSMELL,New loss of smell,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,PHANTOSMIA,NLSMELL,New loss of smell,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,SMELL ALTERATION,NLSMELL,New loss of smell,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,SMELL CHANGE,NLSMELL,New loss of smell,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,DIMINISHED SENSE OF SMELL,NLSMELL,New loss of smell,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,HYPOSMIA,NLSMELL,New loss of smell,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,PHANTOSMIA,NLSMELL,New loss of smell,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,AGEUSIA,NLTASTE,New loss of taste,SMELL_TASTE,DIRECT,TRUE",
  "1.0.0,AELLT,AGEUSIA,NLTASTE,New loss of taste,SMELL_TASTE,DIRECT,TRUE",
  "1.0.0,AELLT,LOSS OF TASTE,NLTASTE,New loss of taste,SMELL_TASTE,DIRECT,TRUE",
  "1.0.0,AETERM,AGEUSIA,NLTASTE,New loss of taste,SMELL_TASTE,DIRECT,TRUE",
  "1.0.0,AETERM,LOSS OF TASTE,NLTASTE,New loss of taste,SMELL_TASTE,DIRECT,TRUE",
  "1.0.0,AETERM,NEW LOSS OF TASTE,NLTASTE,New loss of taste,SMELL_TASTE,DIRECT,TRUE",
  "1.0.0,AEDECOD,DYSGEUSIA,NLTASTE,New loss of taste,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,HYPOGEUSIA,NLTASTE,New loss of taste,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,TASTE DISORDER,NLTASTE,New loss of taste,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,DYSGEUSIA,NLTASTE,New loss of taste,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,HYPOGEUSIA,NLTASTE,New loss of taste,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,IMPAIRED TASTE,NLTASTE,New loss of taste,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,TASTE ALTERATION,NLTASTE,New loss of taste,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,TASTE ALTERED,NLTASTE,New loss of taste,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,TASTE BITTER,NLTASTE,New loss of taste,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,TASTE CHANGED,NLTASTE,New loss of taste,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,TASTE DIMINISHED,NLTASTE,New loss of taste,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,TASTE METALLIC,NLTASTE,New loss of taste,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,TASTE SALTY,NLTASTE,New loss of taste,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,DYSGEUSIA,NLTASTE,New loss of taste,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,HYPOGEUSIA,NLTASTE,New loss of taste,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,TASTE ALTERATION,NLTASTE,New loss of taste,SMELL_TASTE,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,OROPHARYNGEAL PAIN,STHROAT,Sore throat,SORE_THROAT,DIRECT,TRUE",
  "1.0.0,AELLT,CHRONIC SORE THROAT,STHROAT,Sore throat,SORE_THROAT,DIRECT,TRUE",
  "1.0.0,AELLT,OROPHARYNGEAL PAIN,STHROAT,Sore throat,SORE_THROAT,DIRECT,TRUE",
  "1.0.0,AELLT,PAIN THROAT,STHROAT,Sore throat,SORE_THROAT,DIRECT,TRUE",
  "1.0.0,AELLT,SORE THROAT,STHROAT,Sore throat,SORE_THROAT,DIRECT,TRUE",
  "1.0.0,AELLT,STREPTOCOCCAL SORE THROAT,STHROAT,Sore throat,SORE_THROAT,DIRECT,TRUE",
  "1.0.0,AELLT,THROAT PAIN,STHROAT,Sore throat,SORE_THROAT,DIRECT,TRUE",
  "1.0.0,AELLT,THROAT SORE,STHROAT,Sore throat,SORE_THROAT,DIRECT,TRUE",
  "1.0.0,AETERM,SORE THROAT,STHROAT,Sore throat,SORE_THROAT,DIRECT,TRUE",
  "1.0.0,AETERM,SORE THROAT,STHROAT,Sore throat,SORE_THROAT,DIRECT,TRUE",
  "1.0.0,AETERM,THROAT PAIN,STHROAT,Sore throat,SORE_THROAT,DIRECT,TRUE",
  "1.0.0,AEDECOD,DRY THROAT,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,OROPHARYNGEAL DISCOMFORT,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,PHARYNGEAL PARAESTHESIA,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,THROAT IRRITATION,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,THROAT TIGHTNESS,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,DRY THROAT,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,ITCHY THROAT,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,OROPHARYNGEAL DISCOMFORT,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,PHARYNX IRRITATED SENSATION OF,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,THROAT DISCOMFORT,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,THROAT IRRITATION,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,THROAT TIGHTNESS,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,TINGLING THROAT,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,DRY THROAT,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,ITCHY THROAT,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,THROAT DISCOMFORT,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,THROAT IRRITATION,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,THROAT TIGHTNESS,STHROAT,Sore throat,SORE_THROAT,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,NASAL CONGESTION,NASALCON,Nasal congestion,NASAL,DIRECT,TRUE",
  "1.0.0,AELLT,CHRONIC NASAL CONGESTION,NASALCON,Nasal congestion,NASAL,DIRECT,TRUE",
  "1.0.0,AELLT,CONGESTION NASAL,NASALCON,Nasal congestion,NASAL,DIRECT,TRUE",
  "1.0.0,AELLT,NASAL CONGESTION,NASALCON,Nasal congestion,NASAL,DIRECT,TRUE",
  "1.0.0,AELLT,NASAL SINUS CONGESTION,NASALCON,Nasal congestion,NASAL,DIRECT,TRUE",
  "1.0.0,AELLT,NASAL STUFFINESS,NASALCON,Nasal congestion,NASAL,DIRECT,TRUE",
  "1.0.0,AELLT,NOSE CONGESTION,NASALCON,Nasal congestion,NASAL,DIRECT,TRUE",
  "1.0.0,AETERM,CONGESTION NASAL,NASALCON,Nasal congestion,NASAL,DIRECT,TRUE",
  "1.0.0,AETERM,NASAL CONGESTION,NASALCON,Nasal congestion,NASAL,DIRECT,TRUE",
  "1.0.0,AETERM,NASAL CONGESTION,NASALCON,Nasal congestion,NASAL,DIRECT,TRUE",
  "1.0.0,AETERM,NOSE CONGESTION,NASALCON,Nasal congestion,NASAL,DIRECT,TRUE",
  "1.0.0,AEDECOD,RESPIRATORY TRACT CONGESTION,NASALCON,Nasal congestion,NASAL,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,SINUS CONGESTION,NASALCON,Nasal congestion,NASAL,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,UPPER RESPIRATORY TRACT CONGESTION,NASALCON,Nasal congestion,NASAL,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,RESPIRATORY TRACT CONGESTION,NASALCON,Nasal congestion,NASAL,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,SINUS CONGESTION,NASALCON,Nasal congestion,NASAL,RELATED_REVIEW,FALSE",
  "1.0.0,AELLT,UPPER RESPIRATORY TRACT CONGESTION,NASALCON,Nasal congestion,NASAL,RELATED_REVIEW,FALSE",
  "1.0.0,AETERM,SINUS CONGESTION,NASALCON,Nasal congestion,NASAL,RELATED_REVIEW,FALSE",
  "1.0.0,AEDECOD,RHINORRHOEA,RUNNOSE,Runny nose (rhinorrhea),NASAL,DIRECT,TRUE",
  "1.0.0,AELLT,NASAL DISCHARGE,RUNNOSE,Runny nose (rhinorrhea),NASAL,DIRECT,TRUE",
  "1.0.0,AELLT,NASAL MUCUS INCREASED,RUNNOSE,Runny nose (rhinorrhea),NASAL,DIRECT,TRUE",
  "1.0.0,AELLT,RHINORRHEA,RUNNOSE,Runny nose (rhinorrhea),NASAL,DIRECT,TRUE",
  "1.0.0,AELLT,RUNNY NOSE,RUNNOSE,Runny nose (rhinorrhea),NASAL,DIRECT,TRUE",
  "1.0.0,AETERM,RHINORRHEA,RUNNOSE,Runny nose (rhinorrhea),NASAL,DIRECT,TRUE",
  "1.0.0,AETERM,RUNNY NOSE,RUNNOSE,Runny nose (rhinorrhea),NASAL,DIRECT,TRUE",
  "1.0.0,AETERM,RUNNY NOSE,RUNNOSE,Runny nose (rhinorrhea),NASAL,DIRECT,TRUE",
  "1.0.0,AEDECOD,NAUSEA,NAUSEA,Nausea,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AEDECOD,PROCEDURAL NAUSEA,NAUSEA,Nausea,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AELLT,NAUSEA,NAUSEA,Nausea,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AELLT,NAUSEA AGGRAVATED,NAUSEA,Nausea,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AELLT,POST PROCEDURAL NAUSEA,NAUSEA,Nausea,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AELLT,POSTOPERATIVE NAUSEA,NAUSEA,Nausea,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AELLT,PROCEDURAL NAUSEA,NAUSEA,Nausea,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AETERM,NAUSEA,NAUSEA,Nausea,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AETERM,NAUSEA/VOMITING,NAUSEA,Nausea,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AETERM,POSTOPERATIVE NAUSEA,NAUSEA,Nausea,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AEDECOD,PROCEDURAL VOMITING,VOMIT,Vomiting,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AEDECOD,VOMITING,VOMIT,Vomiting,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AELLT,EMESIS,VOMIT,Vomiting,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AELLT,HYPEREMESIS,VOMIT,Vomiting,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AELLT,PROCEDURAL VOMITING,VOMIT,Vomiting,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AELLT,VOMITED,VOMIT,Vomiting,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AELLT,VOMITING,VOMIT,Vomiting,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AETERM,EMESIS,VOMIT,Vomiting,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AETERM,HYPEREMESIS,VOMIT,Vomiting,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AETERM,NAUSEA/VOMITING,VOMIT,Vomiting,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AETERM,VOMITED,VOMIT,Vomiting,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AETERM,VOMITING,VOMIT,Vomiting,NAUSEA_VOMITING,DIRECT,TRUE",
  "1.0.0,AEDECOD,BACTERIAL DIARRHOEA,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AEDECOD,DIARRHOEA,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AEDECOD,DIARRHOEA INFECTIOUS,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AEDECOD,VIRAL DIARRHOEA,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AELLT,ACUTE DIARRHEA,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AELLT,BACTERIAL DIARRHEA,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AELLT,DIARRHEA,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AELLT,DIARRHEA AGGRAVATED,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AELLT,DIARRHOEA,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AELLT,LOOSE BOWEL,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AELLT,LOOSE BOWELS,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AELLT,LOOSE STOOLS,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AELLT,STOOLS WATERY,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AELLT,TRAVELLER'S DIARRHEA,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AELLT,VIRAL DIARRHEA,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AELLT,WATERY DIARRHEA,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AETERM,ACUTE DIARRHEA,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AETERM,DIARRHEA,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE",
  "1.0.0,AETERM,LOOSE STOOLS,DIARRHEA,Diarrhea,DIARRHEA,DIRECT,TRUE"
),collapse="\n"),colClasses="character",check.names=FALSE)
# ---- RUN: select and run the whole script in RStudio ----
# Audit PNG export extension v1.1.0 | 2026-10-23
# Standalone; stable Vxx filenames in <project>/audit_assets plus this run's copy.
# Long tables continue as _02.png, _03.png, etc. References remain in HTML/provenance.
PHMPT_AUDIT_ASSET_DIR <- getOption("phmpt.audit_asset_dir", "audit_assets")
PHMPT_AUDIT_DPI <- getOption("phmpt.audit_dpi", 300L)
PHMPT_AUDIT_WINDOW <- getOption("phmpt.audit_window", "FULL")

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

ws_export_audit_assets <- function(root,out,weekly,ctx,release,ae_timing,window=PHMPT_AUDIT_WINDOW) {
  a<-ax_begin(root,out,"V09_V10_weekly_symptoms")
  w<-weekly[weekly$window==window,];if(!nrow(w))stop("Requested audit symptom window unavailable: ",window)
  dates<-sort(unique(w$week_start));blocks<-split(dates,ceiling(seq_along(dates)/18))
  height<-1.5+5.15*length(blocks)+.85
  top<-max(pretty(c(0,max(c(1,w$dose_1,w$dose_2,w$dose_3,w$dose_4,w$total_symptom_groups),na.rm=TRUE)*1.1),n=5))
  fills<-c("#367E8A","#A65B31","#785E91");arms<-c("mRNA-1273","Placebo")
  for(ar in arms) {
    id<-if(ar=="mRNA-1273")"V09" else "V10"
    draw<-function() {
      old<-par(no.readonly=TRUE);on.exit(par(old))
      for(k in seq_along(blocks)) {
        panel_top<-1.5+(k-1)*5.15;panel_bottom<-panel_top+5.15
        par(fig=c(0,1,1-panel_bottom/height,1-panel_top/height),mai=c(.77,.86,.34,.86),
          new=k>1,family="sans",las=1,xaxs="i",yaxs="i",bty="n",mgp=c(2.5,.5,0),cex=1)
        q<-w[w$arm==ar&w$week_start%in%blocks[[k]],];q<-q[order(q$week_start),];xx<-seq_len(nrow(q));n<-length(xx)
        plot(c(.45,n+.55),c(0,top),type="n",axes=FALSE,xlab="",ylab="")
        ticks<-pretty(c(0,top),n=5);ticks<-ticks[ticks>=0&ticks<=top];abline(h=ticks,col="#E2E7EA",lwd=.6)
        ws_background_doses(q,xx,window=="FULL"&&ar=="Placebo")
        base<-rep(0,n);centres<-matrix(0,n,3)
        for(i in 1:3){v<-q[[c("covid_only","ae_only","both")[i]]];rect(xx-.40,base,xx+.12,base+v,col=fills[i],border=NA);centres[,i]<-base+v/2;base<-base+v}
        axis(2,at=ticks,labels=ax_num(ticks),tick=FALSE,cex.axis=.82)
        axis(4,at=(0:4)*top/4,labels=paste0((0:4)*25,"%"),tick=FALSE,cex.axis=.82)
        mtext("Weekly count",2,line=3.05,las=0,cex=.86);mtext("AE share (%)",4,line=2.75,las=0,cex=.86)
        axis(1,at=xx,labels=FALSE,tck=-.012,col="#87949C")
        dates<-as.Date(q$week_start);labs<-paste(format(dates,"%d"),month.abb[as.integer(format(dates,"%m"))],format(dates,"%y"))
        text(xx,-.026*top,labs,srt=60,adj=1,xpd=NA,cex=.76)
        lines(xx,q$pct_ae_any/100*top,col="#202936",lwd=1.5);points(xx,q$pct_ae_any/100*top,pch=21,bg="white",col="#202936",cex=.55)
        label_cex<-.82;ix<-diff(par("usr")[1:2])/par("pin")[1];iy<-diff(par("usr")[3:4])/par("pin")[2]
        size<-function(l)list(w=(strheight(l,units="inches",cex=label_cex)+.02)*ix,h=(strwidth(l,units="inches",cex=label_cex)+.035)*iy)
        badge<-function(x,y,l,col,inside=FALSE) {
          sz<-size(l);rect(x-sz$w/2,y-sz$h/2,x+sz$w/2,y+sz$h/2,col=if(inside)col else "white",border=col,lwd=.45)
          text(x,y,l,srt=90,cex=label_cex,col=if(inside)"white" else col)
        }
        for(j in seq_len(n))if(q$total_symptom_groups[j]>0) {
          counts<-as.numeric(unlist(q[j,c("covid_only","ae_only","both")],use.names=FALSE));labs<-sprintf("%.1f%%",100*counts/q$total_symptom_groups[j]);sz<-size(labs)
          yy<-ws_label_positions(centres[j,],sz$h,0,top,.025*iy)
          for(i in 1:3){inside<-yy[i]-sz$h[i]/2>=centres[j,i]-counts[i]/2 && yy[i]+sz$h[i]/2<=centres[j,i]+counts[i]/2
            if(!inside)segments(xx[j]-.14,centres[j,i],xx[j]-.14,yy[i],col=fills[i],lwd=.8)
            badge(xx[j]-.14,yy[i],labs[i],fills[i],inside)}
        }
        for(j in which(is.finite(q$pct_ae_any))) {
          lab<-sprintf("%.1f%%",q$pct_ae_any[j]);sz<-size(lab);actual<-q$pct_ae_any[j]/100*top
          yy<-min(top-sz$h/2,max(sz$h/2,actual));segments(xx[j],actual,xx[j]+.31,yy,col="#202936",lwd=.7);badge(xx[j]+.31,yy,lab,"#202936")
        }
        mtext(paste0(min(q$week_start)," to ",max(q$week_start)," | Part ",k," of ",length(blocks)),3,line=.4,adj=0,cex=.86)
      }
      grid::pushViewport(grid::viewport())
      ax_text(paste(id,"Weekly symptom recording:",ar),.35,.23,17,TRUE,col="#17384B")
      ax_text(paste0(release," | ",window," follow-up | Randomized N = ",ax_num(sum(ctx$arm==ar))," | AE timing: ",ae_timing),.35,.67,11.5)
      crossover<-window=="FULL"&&ar=="Placebo"
      labels<-c("Dose 1","Dose 2",if(crossover)c("Dose 3 (open-label 1)","Dose 4 (open-label 2)"),"Symptom table only","AE only","Both sources","AE share")
      colors<-c(ws_dose_fills[seq_len(if(crossover)4L else 2L)],fills,"#202936")
      xs<-if(crossover)c(.35,1.65,3.0,5.45,.35,2.85,4.35,6.0) else c(.35,1.65,.35,2.85,4.35,6.0)
      ys<-if(crossover)c(rep(1.04,4),rep(1.30,4)) else c(rep(1.04,2),rep(1.30,4))
      for(i in seq_along(labels)){grid::grid.rect(x=grid::unit(xs[i]+.06,"in"),y=grid::unit(1,"npc")-grid::unit(ys[i]+.06,"in"),width=grid::unit(.12,"in"),height=grid::unit(.12,"in"),gp=grid::gpar(fill=colors[i],col=if(i<=(if(crossover)4 else 2))ws_dose_borders[i] else NA));ax_text(labels[i],xs[i]+.17,ys[i],10.5)}
      dose_note<-if(window=="FULL"&&anyNA(w[w$arm==ar,c("dose_3","dose_4")]))"Open-label dose counts unavailable in this ADSL; see dose evidence CSV." else if(window=="FULL")"Doses 3/4 are open-label injections in the original placebo arm." else "Doses 1/2 only within this follow-up window."
      ax_text(ax_wrap(paste("One participant x symptom group x week. Colored badges show source shares; black badges show (AE only + both) / all sources. Zero-denominator weeks have no percentage. Offset labels retain the actual bar heights.",dose_note,"Both arms and all panels share the same scales."),7.8,10.6),.35,height-.72,10.6,col="#52616C")
      grid::popViewport()
    }
    ax_png(a,id,if(ar=="mRNA-1273")"weekly_symptoms_mrna1273" else "weekly_symptoms_placebo",draw,height)
    ax_data(a,id,w[w$arm==ar,])
  };ax_finish(a,paste("Selected audit window:",window))
}

if(isTRUE(getOption("phmpt.autorun",TRUE))) {
  analyse_weekly_symptoms_ae(root=PHMPT_ROOT,release=PHMPT_RELEASE,
    adsl_file=PHMPT_ADSL_FILE,ae_file=PHMPT_AE_FILE,symptoms_file=PHMPT_SYMPTOMS_FILE,
    mapping_file=PHMPT_MAPPING_FILE,windows=PHMPT_WINDOWS,eua_cutoff=PHMPT_EUA_CUTOFF,
    ae_timing=PHMPT_AE_TIMING,output_dir=PHMPT_OUTPUT_DIR,encoding=PHMPT_CSV_ENCODING,
    numeric_date_origin=PHMPT_NUMERIC_DATE_ORIGIN,min_randomized=PHMPT_MIN_RANDOMIZED,
    png_dpi=PHMPT_PNG_DPI,show_chart=PHMPT_SHOW_CHART)
}
