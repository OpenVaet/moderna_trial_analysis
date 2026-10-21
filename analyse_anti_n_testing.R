# =============================================================================
# Moderna P301 - anti-N testing coverage | v1.0.0 | 2026-10-21
# Stand-alone R >=4.3: save in project root, Select All / Run in RStudio.
# Base R for CSV/JMP and graphics; install.packages("haven") for XPT/SAS7BDAT.
# Main cohort: randomized participants with ADSL RTPCRBL and ELECSBL negative.
# Day 1 in that cohort is conditional on baseline testing. An all-randomized
# Day 1 comparison is also provided, plus a laboratory-confirmed sensitivity set.
# Production snapshot only, retrospectively restricted to each collection-date
# cutoff. An EUA-date row is NOT the original November 2020 database snapshot.
# Counts measure documented assay results, not orders, blood draws or positivity.
# =============================================================================
PHMPT_ROOT <- ""
PHMPT_ADSL_FILE <- ""
PHMPT_IS_FILE <- ""                 # SDTM IS: Elecsys SARSCOV2
PHMPT_ADIS_FILE <- ""               # optional analysis reconciliation
PHMPT_MB_FILE <- ""                 # SDTM MB: baseline PCR verification
PHMPT_ADMB_FILE <- ""               # optional analysis reconciliation
# Optional alternative: the exact CSV exports from prod_anti_n_pcr_audit.R.
# Leave blank for normal direct discovery of local production source tables.
PHMPT_SEROLOGY_RECORDS_FILE <- ""
PHMPT_PCR_RECORDS_FILE <- ""
PHMPT_EUA_CUTOFF <- "2020-11-25" # final EUA safety cutoff; efficacy cutoff was 2020-11-21
PHMPT_BLA_CUTOFF <- "2021-03-26"
PHMPT_DAY57_ANCHOR <- "dose2"       # dose2 +28 days; dose1 +56 if no dose 2
# "dose1" uses the nominal first-dose +56 calendar for every participant.
PHMPT_DAY57_EARLY <- 3L
PHMPT_DAY57_LATE <- 7L
PHMPT_BASELINE_LOOKBACK <- 28L
PHMPT_MIN_RANDOMIZED <- 30000L
PHMPT_CSV_ENCODING <- "auto"
PHMPT_NUMERIC_DATE_ORIGIN <- ""
PHMPT_OUTPUT_DIR <- ""
PHMPT_PNG_DPI <- 300L
PHMPT_SHOW_REPORT <- interactive()
# To load functions only: options(phmpt.autorun=FALSE) before source().

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

nt_subjects <- function(ctx,sero,pcr,anchor,early,late,lookback,cutoffs,out) {
  ctx$baseline_date<-ifelse(is.finite(ctx$TR01SDT),ctx$TR01SDT,ctx$RANDDT)
  ctx$baseline_anchor<-ifelse(is.finite(ctx$TR01SDT),"DOSE1","RANDOMIZATION_NO_DOSE1")
  ctx$baseline_flags_negative<-nt_result(ctx$ELECSBL)=="NEGATIVE"&nt_result(ctx$RTPCRBL)=="NEGATIVE"
  ctx$baseline_composite_conflict<-ctx$baseline_flags_negative&nt_result(ctx$SCOV2BL)=="POSITIVE"
  ctx$cohort<-ctx$baseline_flags_negative&!ctx$baseline_composite_conflict
  ctx$target_date<-ctx$TR01SDT+56L
  second<-anchor=="dose2"&is.finite(ctx$DOSE2DT)
  ctx$target_date[second]<-ctx$DOSE2DT[second]+28L
  ctx$target_anchor<-ifelse(second,"DOSE2_PLUS_28",ifelse(is.finite(ctx$TR01SDT),"DOSE1_PLUS_56","NO_FIRST_DOSE"))
  ctx$window_start<-ctx$target_date-early;ctx$window_end<-ctx$target_date+late;ctx$week<-ws_week(ctx$target_date)
  ctx$left_before_window_end<-is.finite(ctx$window_end)&ctx$end_full<ctx$window_end
  ctx$unblinded_before_window_end<-is.finite(ctx$window_end)&ctx$end_blinded<ctx$window_end
  n<-nrow(ctx)
  aggregate_records<-function(p,kind) {
    j<-match(p$subject,ctx$subject);valid<-!is.na(j)
    ws_csv(nt_dates_out(p[!valid,,drop=FALSE]),file.path(out,paste0(kind,"_outside_randomized_roster.csv")))
    p<-p[valid,,drop=FALSE];j<-j[valid]
    count<-function(flag)tabulate(j[which(flag)],n)>0L
    before<-p$dated&p$date>=ctx$baseline_date[j]-lookback&p$date<=ctx$baseline_date[j]&!p$open_label&!p$pdv
    list(p=p,j=j,count=count,baseline_negative=count(before&p$binary_result&p$negative),
      baseline_positive=count(before&p$positive),baseline_conflict=count(before&p$result_conflict),
      baseline_test=count(before&p$tested),day1_test=count(p$tested&p$date==ctx$baseline_date[j]&!p$open_label&!p$pdv),
      day1_label_test=count(p$tested&p$day1_label&before))
  }
  ns<-aggregate_records(sero,"serology");ps<-aggregate_records(pcr,"pcr")
  for(kind in c("n","pcr")) {
    a<-if(kind=="n")ns else ps
    for(field in c("baseline_negative","baseline_positive","baseline_conflict","baseline_test","day1_test","day1_label_test"))ctx[[paste0(kind,"_",field)]]<-a[[field]]
  }
  ctx$baseline_lab_conflict<-ctx$n_baseline_positive|ctx$pcr_baseline_positive|ctx$n_baseline_conflict|ctx$pcr_baseline_conflict
  ctx$cohort_verified<-ctx$cohort&ctx$n_baseline_negative&ctx$pcr_baseline_negative&!ctx$baseline_lab_conflict
  p<-ns$p;j<-ns$j
  p$cohort<-ctx$cohort[j];p$target_date<-ctx$target_date[j]
  p$window_start<-ctx$window_start[j];p$window_end<-ctx$window_end[j]
  p$initial_phase<-!p$open_label&!p$pdv
  p$near_day57<-p$tested&p$initial_phase&is.finite(p$target_date)&p$date>=p$window_start&p$date<=p$window_end
  p$labelled_day57<-p$tested&p$day57_label&p$date>ctx$baseline_date[j]
  p$day57_coverage<-p$near_day57|p$labelled_day57
  p$after_baseline<-p$tested&p$date>ctx$baseline_date[j]&!p$baseline_flag&!p$day1_label
  p$after_day57_window<-p$after_baseline&is.finite(p$window_end)&p$date>p$window_end
  for(scope in names(cutoffs)) {
    cap<-cutoffs[[scope]];bycut<-p$date<=cap
    eligible<-ctx$cohort&ctx$baseline_date<=cap
    due<-eligible&ctx$DOSE1FL=="Y"&is.finite(ctx$window_end)&ctx$window_end<=cap
    ctx[[paste0(scope,"_eligible")]]<-eligible
    ctx[[paste0(scope,"_due")]]<-due
    ctx[[paste0(scope,"_not_yet_mature")]]<-eligible&is.finite(ctx$window_end)&ctx$window_end>cap
    for(field in c("day57_coverage","near_day57","labelled_day57","after_baseline","after_day57_window")) {
      ctx[[paste0(scope,"_",field)]]<-ns$count(p[[field]]&bycut)
    }
    ctx[[paste0(scope,"_day57_binary")]]<-ns$count(p$day57_coverage&bycut&p$binary_result)
    ctx[[paste0(scope,"_day57_openlabel_only")]]<-ns$count(p$tested&p$open_label&is.finite(p$target_date)&p$date>=p$window_start&p$date<=p$window_end&bycut)
    ctx[[paste0(scope,"_labelled_late")]]<-ns$count(p$labelled_day57&!p$near_day57&bycut)
  }
  ws_csv(nt_dates_out(p),file.path(out,"serology_subject_visit_links.csv"))
  ws_csv(nt_dates_out(ps$p),file.path(out,"pcr_baseline_audit_records.csv"))
  ws_csv(nt_dates_out(ctx),file.path(out,"subject_testing_status.csv"))
  ws_csv(nt_dates_out(ctx[ctx$cohort,,drop=FALSE]),file.path(out,"baseline_negative_subjects.csv"))
  ws_csv(nt_dates_out(ctx[ctx$cohort&!ctx$cohort_verified,,drop=FALSE]),file.path(out,"baseline_verification_review.csv"))
  ctx
}

nt_p <- function(p) {
  vapply(p,function(v) {
    if(!is.finite(v))return("NA")
    for(k in 6:2)if(v<10^-k)return(paste0("< ",format(10^-k,scientific=FALSE,trim=TRUE)))
    sprintf("%.3f",v)
  },character(1))
}
nt_wilson <- function(k,n) {
  if(!n)return(c(lower=NA_real_,upper=NA_real_))
  z<-qnorm(.975);p<-k/n;den<-1+z*z/n;c(lower=(p+z*z/(2*n)-z*sqrt(p*(1-p)/n+z*z/(4*n*n)))/den,
    upper=(p+z*z/(2*n)+z*sqrt(p*(1-p)/n+z*z/(4*n*n)))/den)
}
nt_compare <- function(ctx,eligible,tested,endpoint,scope,cutoff,cohort_label,conditional=FALSE) {
  counts<-vapply(c("mRNA-1273","Placebo"),function(a){q<-eligible&ctx$arm==a;c(tested=sum(q&tested),total=sum(q))},numeric(2))
  v<-counts[1,1];vn<-counts[2,1];p<-counts[1,2];pn<-counts[2,2]
  mat<-rbind(c(v,vn-v),c(p,pn-p));stat<-pv<-pf<-expected<-NA_real_;method<-"Not estimable"
  status<-if(!vn||!pn)"EMPTY_ARM" else if(any(colSums(mat)==0))"NO_OUTCOME_VARIATION" else "OK"
  if(status=="OK"&&!conditional) {
    chi<-suppressWarnings(chisq.test(mat,correct=FALSE));stat<-unname(chi$statistic);pv<-chi$p.value;expected<-min(chi$expected)
    if(expected<5){pf<-fisher.test(mat,alternative="two.sided")$p.value;method<-"Fisher (small expected count)"}
    else method<-"Pearson chi-square"
  }
  if(conditional){status<-"CONDITIONAL_ON_BASELINE_TESTING";method<-"Descriptive: selected on baseline results"}
  vc<-nt_wilson(v,vn);pc<-nt_wilson(p,pn)
  data.frame(cohort=cohort_label,endpoint=endpoint,scope=scope,cutoff=ws_iso(cutoff),
    vaccine_tested=v,vaccine_total=vn,vaccine_untested=vn-v,vaccine_pct=if(vn)100*v/vn else NA,
    placebo_tested=p,placebo_total=pn,placebo_untested=pn-p,placebo_pct=if(pn)100*p/pn else NA,
    difference_pp=if(vn&&pn)100*(v/vn-p/pn) else NA,
    vaccine_ci_lower=100*vc[1],vaccine_ci_upper=100*vc[2],placebo_ci_lower=100*pc[1],placebo_ci_upper=100*pc[2],
    pearson_chi_square=stat,df=if(is.finite(stat))1L else NA_integer_,pearson_p=pv,min_expected=expected,
    fisher_p=pf,selected_p=if(is.finite(pf))pf else pv,selected_test=method,status=status,stringsAsFactors=FALSE,row.names=NULL)
}
nt_summaries <- function(ctx,cutoffs,verified=FALSE) {
  cohort<-if(verified)ctx$cohort_verified else ctx$cohort
  label<-if(verified)"Baseline negative, laboratory-confirmed subset" else "Baseline negative, ADSL-defined"
  ans<-list();add<-function(e,t,name,scope,conditional=FALSE) {
    ans[[length(ans)+1L]]<<-nt_compare(ctx,e,t,name,scope,cutoffs[[scope]],label,conditional)
  }
  # A fixed baseline visit comparison; selecting negative anti-N at baseline
  # makes this cohort unsuitable for inferring initial testing access.
  for(sc in names(cutoffs)) {
    add(cohort&ctx$baseline_date<=cutoffs[[sc]],ctx$n_day1_label_test|ctx$n_day1_test,"Day 1 visit/date anti-N result",sc,TRUE)
    e<-ctx[[paste0(sc,"_eligible")]]&cohort;d<-ctx[[paste0(sc,"_due")]]&cohort
    add(d,ctx[[paste0(sc,"_day57_coverage")]],"Day 57 labelled visit or nearby test",sc)
    add(d,ctx[[paste0(sc,"_near_day57")]],"Anti-N result within Day 57 window",sc)
    add(d,ctx[[paste0(sc,"_labelled_day57")]],"Labelled Day 57 visit (any timing)",sc)
    add(e,ctx[[paste0(sc,"_after_baseline")]],"At least one post-baseline anti-N test",sc)
    add(d&ctx$window_end<cutoffs[[sc]],ctx[[paste0(sc,"_after_day57_window")]],"At least one test after Day 57 window",sc)
  }
  z<-do.call(rbind,ans);z$selected_p_holm<-p.adjust(z$selected_p,"holm")
  z$p_display<-nt_p(z$selected_p);z$p_holm_display<-nt_p(z$selected_p_holm);z
}
nt_weekly <- function(ctx,cap,verified=FALSE) {
  cohort<-if(verified)ctx$cohort_verified else ctx$cohort
  good<-cohort&ctx$BLA_due
  if(!any(good))stop("No completed Day 57 observation windows by BLA cutoff.")
  weeks<-seq(min(ctx$week[good]),ws_week(cap),7L)
  rows<-list()
  for(w in weeks)for(a in c("mRNA-1273","Placebo")) {
    i<-which(good&ctx$arm==a&ctx$week==w);n<-length(i);k<-sum(ctx$BLA_day57_coverage[i]);ci<-nt_wilson(k,n)
    rows[[length(rows)+1L]]<-data.frame(week=w,week_start=ws_iso(w),week_end=ws_iso(w+6),arm=a,
      due=n,tested=k,not_tested=n-k,testing_pct=if(n)100*k/n else NA_real_,ci_lower=100*ci[1],ci_upper=100*ci[2],
      within_window=sum(ctx$BLA_near_day57[i]),within_window_pct=if(n)100*sum(ctx$BLA_near_day57[i])/n else NA_real_,
      labelled_visit=sum(ctx$BLA_labelled_day57[i]),
      labelled_outside_window=sum(ctx$BLA_labelled_late[i]),
      left_before_window_end=sum(ctx$left_before_window_end[i]),unblinded_before_window_end=sum(ctx$unblinded_before_window_end[i]),
      no_second_dose=sum(!is.finite(ctx$DOSE2DT[i])),fallback_dose1=sum(ctx$target_anchor[i]=="DOSE1_PLUS_56"),
      cohort=if(verified)"LAB_VERIFIED" else "ADSL_BASELINE_NEGATIVE",row.names=NULL)
  }
  z<-do.call(rbind,rows);stopifnot(sum(z$due)==sum(good),sum(z$tested)==sum(good&ctx$BLA_day57_coverage),all(z$tested<=z$due));z
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





nt_pdf <- function(path,width,height) {
  if(capabilities("cairo")) grDevices::cairo_pdf(path,width=width,height=height,onefile=TRUE,family="sans")
  else grDevices::pdf(path,width=width,height=height,onefile=TRUE)
}
nt_display <- function(z) {
  pct<-function(k,n,p)ifelse(n>0,paste0(ws_n(k)," / ",ws_n(n)," (",sprintf("%.2f",p),"%)"),"No eligible subjects")
  data.frame(Endpoint=z$endpoint,Cutoff=z$scope,
    `mRNA-1273: tested / eligible`=pct(z$vaccine_tested,z$vaccine_total,z$vaccine_pct),
    `Placebo: tested / eligible`=pct(z$placebo_tested,z$placebo_total,z$placebo_pct),
    `Difference (pp)`=ifelse(is.finite(z$difference_pp),sprintf("%+.2f",z$difference_pp),"NA"),
    `Chi-square`=ifelse(is.finite(z$pearson_chi_square),sprintf("%.3f",z$pearson_chi_square),"NA"),
    p=z$p_display,`p (Holm)`=z$p_holm_display,check.names=FALSE)
}
nt_html <- function(x) {
  x<-gsub("&","&amp;",as.character(x),fixed=TRUE);x<-gsub("<","&lt;",x,fixed=TRUE);gsub(">","&gt;",x,fixed=TRUE)
}
nt_table_html <- function(x) paste0("<table><thead><tr>",paste0("<th>",nt_html(names(x)),"</th>",collapse=""),
  "</tr></thead><tbody>",paste(apply(x,1,function(r)paste0("<tr>",paste0("<td>",nt_html(r),"</td>",collapse=""),"</tr>")),collapse="\n"),"</tbody></table>")
nt_plot <- function(w,out,dpi,cutoff,anchor) {
  draw<-function() {
    old<-par(no.readonly=TRUE);on.exit(par(old));layout(matrix(1:3,ncol=1),heights=c(6,2.3,1.25))
    colors<-c("mRNA-1273"="#226A8A","Placebo"="#AA5A21");weeks<-sort(unique(w$week));lim<-range(weeks)+c(-4,4)
    ticks<-weeks[seq(1,length(weeks),by=if(length(weeks)>18)2L else 1L)]
    par(mar=c(.7,4.8,4.3,1.2),mgp=c(2.8,.7,0),las=1)
    plot(NA,xlim=lim,ylim=c(0,103),xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="Subjects tested (%)",bty="l")
    abline(h=seq(0,100,20),col="#E4E8EA");axis(2,at=seq(0,100,20),labels=paste0(seq(0,100,20),"%"))
    for(a in names(colors)) {
      z<-w[w$arm==a,];xx<-z$week+if(a=="mRNA-1273")-.65 else .65;ok<-z$due>0
      arrows(xx[ok],z$ci_lower[ok],xx[ok],z$ci_upper[ok],angle=90,code=3,length=.025,col=adjustcolor(colors[a],.5),lwd=.8)
      lines(xx,z$testing_pct,col=colors[a],lwd=1.5);points(xx,z$testing_pct,pch=if(a=="mRNA-1273")16 else 17,col=colors[a],cex=.65)
    }
    title("Anti-N testing coverage at the original Day 57 visit",adj=0,cex.main=1.18,line=2.6)
    mtext(paste0("Baseline PCR- and anti-N-negative population | Results through ",ws_iso(cutoff)),side=3,line=1.15,adj=0,cex=.85)
    legend("bottomleft",legend=names(colors),col=colors,lty=1,pch=c(16,17),bty="n",horiz=TRUE,inset=c(.01,.03),cex=.9)
    par(mar=c(3.4,4.8,.5,1.2),mgp=c(2.8,.7,0),las=1)
    plot(NA,xlim=lim,ylim=c(0,max(w$due)*1.12),xaxs="i",yaxs="i",xaxt="n",xlab="",ylab="Subjects due",bty="l")
    for(a in names(colors)) {
      z<-w[w$arm==a,];xx<-z$week+if(a=="mRNA-1273")-1.5 else 1.5
      rect(xx-1.35,0,xx+1.35,z$due,col=adjustcolor(colors[a],.55),border=NA)
    }
    axis(1,at=ticks,labels=format(as.Date(ticks,origin="1970-01-01"),"%d %b\n%Y"),cex.axis=.65,tck=-.025)
    mtext("Scheduled Day 57 week (Monday start)",side=1,line=2.4,cex=.8)
    par(mar=c(.2,.6,.1,.6));plot.new()
    notes<-c("Points: tested / subjects due; bars: denominator. Whiskers: 95% Wilson confidence intervals.",
      "Coverage: original labelled Day 57 result, or a result in its nearby window; completed windows only.",
      paste0("Target: ",if(anchor=="dose2")"dose 2 +28 days (dose 1 +56 if no dose 2)" else "dose 1 +56 days", ". Open-label Day 57 excluded."),
      "Missing weeks are gaps. Late labelled visits count as coverage; within-window rates are reported separately.")
    text(.01,seq(.93,.16,length.out=length(notes)),labels=notes,adj=c(0,1),cex=.73,col="#424D55")
  }
  nt_pdf(file.path(out,"weekly_day57_testing.pdf"),width=11.7,height=8.3);tryCatch(draw(),finally=dev.off())
  png(file.path(out,"weekly_day57_testing.png"),width=11.7,height=8.3,units="in",res=dpi);tryCatch(draw(),finally=dev.off())
}
nt_summary_pdf <- function(z,allrand,out,methods) {
  nt_pdf(file.path(out,"anti_n_testing_summary.pdf"),width=11.7,height=8.3);on.exit(dev.off())
  draw<-function(d,title,subtitle) {
    par(mar=rep(.4,4));plot.new();plot.window(c(0,1),c(0,1),xaxs="i",yaxs="i")
    text(.02,.96,title,adj=0,font=2,cex=1.25);text(.02,.915,subtitle,adj=0,cex=.85)
    labels<-c("Visit / endpoint","Cutoff","mRNA-1273\ntested / eligible (%)","Placebo\ntested / eligible (%)","Diff.\npp","Chi-sq.","p","p Holm")
    x<-c(.02,.365,.425,.59,.765,.824,.88,.945);height<-.048;top<-.82
    text(x,rep(.855,8),labels,adj=c(0,.5),cex=.7,font=2)
    for(i in seq_len(nrow(d))) {
      y<-top-(i-.5)*height
      if(i%%2L==1L)rect(.014,y-height/2,.99,y+height/2,col="#EEF3F5",border=NA)
      vals<-as.character(d[i,]);text(x,rep(y,8),vals,adj=c(0,.5),cex=.64)
    }
    foot<-c("Counts are unique subjects with documented anti-N assay results. Difference = vaccine minus placebo, percentage points.",
      "Day 57 denominator: baseline-negative dosed subjects whose whole visit window had elapsed by the cutoff; withdrawals remain included.",
      "Day 1 in the selected cohort is descriptive. Pearson chi-square has df=1, no continuity correction; Fisher p is selected if expected <5.",
      "Holm p adjusts the displayed family of exploratory comparisons. Full precision, selected test and source records are in the CSV files.",
      "March 2021 production snapshot, retrospectively date-restricted. Missing documentation is not proof that an assay was never performed.")
    text(.02,seq(.19,.055,length.out=length(foot)),foot,adj=0,cex=.7,col="#424D55")
  }
  draw(nt_display(z),"Moderna P301 | Anti-N testing coverage","Main cohort: randomized participants with negative baseline ADSL PCR and Elecsys anti-N flags")
  draw(nt_display(allrand),"Day 1 coverage | All randomized participants","Supplement: avoids conditioning the Day 1 comparison on an already documented negative baseline result")
  par(mar=rep(.6,4));plot.new();text(0,.98,"Definitions and interpretation",adj=0,font=2,cex=1.2)
  y<-.91
  for(m in methods) {
    ln<-strwrap(m,width=135);h<-length(ln)*.019+.018
    if(y-h<.04){plot.new();text(0,.98,"Definitions and interpretation (continued)",adj=0,font=2,cex=1.2);y<-.91}
    text(0,y,paste(ln,collapse="\n"),adj=c(0,1),cex=.7);y<-y-h
  }
}
nt_report <- function(summary,verified,allrand,ctx,w,out,cutoffs,anchor,early,late,lookback) {
  methods<-c(
    "Population: randomized P301 participants, original TRT01P arms, with ADSL RTPCRBL=Negative and ELECSBL=Negative. A contradictory positive composite SCOV2BL flag is excluded. Later infections do not remove subjects from this cohort.",
    paste0("Baseline verification: negative PCR and anti-N records on or before dose 1 (randomization if undosed), within ",lookback," days, with no conflicting positive result. The separate laboratory-confirmed sensitivity cohort requires both. Main-cohort records needing review are exported, not silently reassigned."),
    "Day 1 coverage: a completed original Day 1-labelled assay in the baseline window, or an assay on the baseline anchor date. Because the main cohort is selected on negative baseline results, this row is descriptive. All-randomized Day 1 comparisons are supplied separately.",
    paste0("Day 57 target: ",if(anchor=="dose2")"second-dose date +28 days, with first-dose date +56 days for subjects without dose 2" else "first-dose date +56 days", ". Nearby window: target minus ",early," through plus ",late," days, inclusive. This audit rule does not itself classify protocol deviations."),
    "Day 57 coverage: at least one original Day 57-labelled result by cutoff, OR an initial-series anti-N result within the nearby window. Late labelled samples count as coverage; within-window and labelled-only results are also reported. A nearby illness sample documents antibody testing but does not prove a scheduled visit occurred. PDV and open-label samples do not qualify as nearby original Day 57 tests.",
    "Day 57 denominator: baseline-negative dosed subjects whose complete nearby window ends on or before the cutoff. Subjects without dose 1 are excluded. Withdrawal, death and unblinding do not silently remove subjects; their counts are exported. Comparisons measure documented coverage, not whether staff withheld a test.",
    "Cumulative coverage: at least one completed assay after baseline by the cutoff. Baseline-labelled and baseline-flagged records are excluded. A further row counts testing after the Day 57 window. BLA cumulative coverage includes PDV/open-label samples, grouped by original randomized assignment; it is not an efficacy comparison.",
    "Assay: IS ISTESTCD=SARSCOV2, reconciled with ADIS PARAMCD=SARSCOV2 / Elecsys when available. Other antibody assays are excluded. Dated positive, negative, or explicit indeterminate/invalid results count as assay evidence. Not-done, derived, unknown-result and conflicting/partial-date records do not. Positive/negative conflicts indicate an assay but cannot verify negative baseline status.",
    "Reconciliation: raw and analysis rows link on subject plus ISSEQ (MBSEQ for PCR), avoiding double counting. These keys identify records, not necessarily separate physical specimens. Raw collection date is preferred; contradictory or explicitly partial dates are flagged. Analysis-only assay rows are retained and identified. Counts are unique subjects.",
    paste0("Cutoffs: EUA ",ws_iso(cutoffs$EUA),"; BLA ",ws_iso(cutoffs$BLA),". The default EUA date is the final safety cutoff (the final efficacy cutoff was 21 November 2020); change it at the top if needed. Both are retrospective collection-date restrictions of the selected production snapshot, not a reconstruction of what was available in the original EUA database."),
    "Inference: Pearson 2x2 tested/untested chi-square, df=1, without continuity correction. If any expected count is below 5, Pearson values remain in CSV and two-sided Fisher becomes the selected p. Empty arms or constant outcomes are not estimable. Holm adjustment is within each summary family; endpoints overlap and are exploratory. No site, calendar-time or follow-up adjustment is claimed.",
    "Weekly figure: target-week cohorts (Monday start), numerator Day 57 coverage by BLA, denominator completed windows. Each subject appears in one due-week only. Missing denominators remain gaps; 95% Wilson intervals and denominator columns are shown. Within-window rates and subject identifiers are available in CSV.")
  writeLines(c("Moderna P301 anti-N testing coverage v1.0.0",methods,"",
    "Protocol reference: Amendment 3, 20 August 2020, Table 14 and section 8.1.5:",
    "https://docs.house.gov/meetings/IF/IF02/20200930/111063/HHRG-116-IF02-20200930-SD008.pdf",
    "FDA EUA review (cutoff dates, p.19): https://www.fda.gov/media/144673/download"),file.path(out,"METHODS.txt"),useBytes=TRUE)
  flow<-do.call(rbind,lapply(c("mRNA-1273","Placebo"),function(a) {
    j<-ctx$arm==a;data.frame(arm=a,randomized=sum(j),baseline_negative=sum(j&ctx$cohort),
      lab_verified=sum(j&ctx$cohort_verified),baseline_review=sum(j&ctx$cohort&!ctx$cohort_verified),
      no_first_dose=sum(j&ctx$cohort&ctx$DOSE1FL!="Y"),target_fallback_dose1=sum(j&ctx$cohort&ctx$target_anchor=="DOSE1_PLUS_56"),
      EUA_mature_windows=sum(j&ctx$EUA_due),BLA_mature_windows=sum(j&ctx$BLA_due),
      BLA_immature_windows=sum(j&ctx$BLA_not_yet_mature),
      left_before_window_end=sum(j&ctx$cohort&ctx$left_before_window_end),
      unblinded_before_window_end=sum(j&ctx$cohort&ctx$unblinded_before_window_end),row.names=NULL)
  }))
  ws_csv(flow,file.path(out,"population_flow.csv"))
  html<-c('<!doctype html><html lang="en"><meta charset="utf-8"><title>Moderna anti-N testing audit</title>',
    '<style>body{font:15px Arial,sans-serif;max-width:1400px;margin:35px auto;color:#202f38;padding:0 18px}h1{font-size:29px}h2{margin-top:32px}table{border-collapse:collapse;width:100%;font-size:12px}th,td{padding:9px 6px;border-bottom:1px solid #d6dfe4;text-align:right}th:first-child,td:first-child{text-align:left}th{background:#e5eef2}tr:nth-child(even){background:#f5f7f8}img{width:100%}.note{border-left:4px solid #226a8a;padding:12px;background:#f0f5f7}li{margin:10px 0;line-height:1.45}@page{size:A4 landscape;margin:12mm}@media print{body{margin:0;padding:0}h2{break-after:avoid}table,img{break-inside:avoid}details{display:block}a{color:inherit}}</style>',
    '<h1>Moderna P301 | Anti-N testing coverage</h1>',
    paste0('<p>Production data through ',ws_iso(cutoffs$BLA),'; EUA date restriction ',ws_iso(cutoffs$EUA),'. Original randomized arms.</p>'),
    '<p class="note">Day 1 testing in a cohort selected for negative baseline anti-N is conditional on testing. Use the all-randomized supplement to assess Day 1 testing imbalance. Missing records indicate missing documentation, not necessarily a missed blood draw.</p>',
    '<h2>Baseline-negative cohort</h2>',nt_table_html(nt_display(summary)),
    '<p>pp: percentage points. p: Pearson, or Fisher when expected counts are small. The CSV retains both tests and exact p-values.</p>',
    '<h2>Weekly Day 57 coverage</h2><img src="weekly_day57_testing.png" alt="Day 57 anti-N testing percentages by arm, with confidence intervals and numbers due">',
    '<p><a href="weekly_day57_testing.pdf">Print chart (PDF)</a> | <a href="anti_n_testing_summary.pdf">Print summary (PDF)</a> | <a href="weekly_day57_testing.csv">Weekly CSV</a> | <a href="baseline_negative_subjects.csv">Subject-level CSV</a></p>',
    '<h2>All-randomized baseline supplement</h2>',nt_table_html(nt_display(allrand)),
    '<h2>Laboratory-confirmed sensitivity cohort</h2>',nt_table_html(nt_display(verified)),
    '<h2>Population accounting</h2>',nt_table_html(flow),
    '<h2>Definitions and interpretation</h2><ol>',paste0('<li>',nt_html(methods),'</li>'),'</ol>',
    '<p>Protocol: Amendment 3, 20 August 2020, Table 14 and section 8.1.5. <a href="https://docs.house.gov/meetings/IF/IF02/20200930/111063/HHRG-116-IF02-20200930-SD008.pdf">Protocol PDF</a>.</p></html>')
  writeLines(html,file.path(out,"anti_n_testing_report.html"),useBytes=TRUE)
  nt_summary_pdf(summary,allrand,out,methods)
}

analyse_anti_n_testing <- function(root=PHMPT_ROOT,adsl_file=PHMPT_ADSL_FILE,is_file=PHMPT_IS_FILE,adis_file=PHMPT_ADIS_FILE,
  mb_file=PHMPT_MB_FILE,admb_file=PHMPT_ADMB_FILE,serology_records_file=PHMPT_SEROLOGY_RECORDS_FILE,pcr_records_file=PHMPT_PCR_RECORDS_FILE,
  eua_cutoff=PHMPT_EUA_CUTOFF,bla_cutoff=PHMPT_BLA_CUTOFF,day57_anchor=PHMPT_DAY57_ANCHOR,
  day57_early=PHMPT_DAY57_EARLY,day57_late=PHMPT_DAY57_LATE,baseline_lookback=PHMPT_BASELINE_LOOKBACK,
  min_randomized=PHMPT_MIN_RANDOMIZED,encoding=PHMPT_CSV_ENCODING,numeric_date_origin=PHMPT_NUMERIC_DATE_ORIGIN,
  output_dir=PHMPT_OUTPUT_DIR,png_dpi=PHMPT_PNG_DPI,show_report=PHMPT_SHOW_REPORT) {
  root<-ws_root(root);day57_anchor<-match.arg(day57_anchor,c("dose2","dose1"))
  for(v in c(day57_early,day57_late,baseline_lookback))if(!is.finite(v)||v<0||v!=floor(v))stop("Visit window settings must be non-negative whole days.")
  cutoffs<-list(EUA=ws_dates(eua_cutoff,"EUA cutoff",strict=TRUE),BLA=ws_dates(bla_cutoff,"BLA cutoff",strict=TRUE))
  if(any(lengths(cutoffs)!=1L)||any(!is.finite(unlist(cutoffs)))||cutoffs$EUA>cutoffs$BLA)stop("Invalid cutoff dates.")
  out<-if(nzchar(output_dir))ws_path(root,output_dir) else file.path(root,"analysis",paste0("anti_n_testing_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid()))
  if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop("Output folder is not empty: ",out)
  dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out,winslash="/",mustWork=TRUE)
  writeLines("RUNNING",file.path(out,"RUN_STATUS.txt"));complete<-FALSE
  on.exit(if(!complete)writeLines("FAILED - incomplete results; see R console error and discovery/quality CSVs.",file.path(out,"RUN_STATUS.txt")),add=TRUE)
  message("Anti-N testing audit v1.0.0 | output: ",out)
  adsl<-ws_pick(root,"Prod","ADSL",adsl_file,out,encoding,numeric_date_origin,min_randomized)
  ctx<-ws_context(adsl$x,numeric_date_origin)
  if(unique(ctx$CUTOFFDT)!=cutoffs$BLA)stop("ADSL cutoff is ",ws_iso(unique(ctx$CUTOFFDT)),"; expected BLA cutoff ",ws_iso(cutoffs$BLA),". Select the intended production snapshot.")
  manifest<-adsl$manifest;submission<-ws_submission(manifest$path)
  if(!nzchar(submission)&&!nzchar(is_file)&&!nzchar(serology_records_file))warning("ADSL filename has no submission tag. IS/ADIS consistency still checked; review source_manifest.csv.")
  get_domain<-function(domain,explicit,optional=FALSE) {
    z<-if(optional)nt_optional(root,domain,explicit,out,encoding,numeric_date_origin,min_randomized,submission) else ws_pick(root,"Prod",domain,explicit,out,encoding,numeric_date_origin,min_randomized,submission)
    if(!is.null(z))manifest<<-rbind(manifest,z$manifest)
    if(is.null(z))NULL else z$x
  }
  get_export<-function(file,kind) {
    p<-ws_path(root,file);z<-nt_read_export(p,kind,out)
    manifest<<-rbind(manifest,data.frame(domain=paste0(kind,"_RECONCILED_EXPORT"),release="Prod",path=normalizePath(p,winslash="/"),md5=unname(tools::md5sum(p)),rows=nrow(z),cutoff="",status="SELECTED",detail="User-selected prior audit export; verify its provenance matches selected ADSL."))
    z
  }
  if(nzchar(serology_records_file))sero<-get_export(serology_records_file,"serology") else {
    raw<-get_domain("IS",is_file);a<-get_domain("ADIS",adis_file,TRUE)
    sero<-nt_reconcile(raw,a,"serology",numeric_date_origin,out)
  }
  if(nzchar(pcr_records_file))pcr<-get_export(pcr_records_file,"pcr") else {
    raw<-get_domain("MB",mb_file);a<-get_domain("ADMB",admb_file,TRUE)
    pcr<-nt_reconcile(raw,a,"pcr",numeric_date_origin,out)
  }
  ws_csv(manifest,file.path(out,"source_manifest.csv"))
  ctx<-nt_subjects(ctx,sero,pcr,day57_anchor,day57_early,day57_late,baseline_lookback,cutoffs,out)
  if(!any(ctx$cohort))stop("No baseline-negative randomized subjects.")
  if(anyNA(ctx$cohort)||anyNA(ctx$EUA_eligible)||anyNA(ctx$BLA_due))stop("Incomplete subject eligibility logic; inspect baseline dates.")
  summary<-nt_summaries(ctx,cutoffs);verified<-nt_summaries(ctx,cutoffs,TRUE)
  allrand<-do.call(rbind,lapply(names(cutoffs),function(sc)rbind(
    nt_compare(ctx,ctx$baseline_date<=cutoffs[[sc]],ctx$n_day1_label_test|ctx$n_day1_test,"Day 1 visit/date anti-N result",sc,cutoffs[[sc]],"All randomized"),
    nt_compare(ctx,ctx$baseline_date<=cutoffs[[sc]],ctx$n_baseline_test,"Any anti-N result in baseline window",sc,cutoffs[[sc]],"All randomized"))))
  allrand$selected_p_holm<-p.adjust(allrand$selected_p,"holm");allrand$p_display<-nt_p(allrand$selected_p);allrand$p_holm_display<-nt_p(allrand$selected_p_holm)
  w<-nt_weekly(ctx,cutoffs$BLA)
  ws_csv(summary,file.path(out,"testing_comparisons.csv"));ws_csv(verified,file.path(out,"testing_comparisons_lab_verified.csv"))
  ws_csv(allrand,file.path(out,"baseline_all_randomized.csv"));ws_csv(w,file.path(out,"weekly_day57_testing.csv"))
  if(any(ctx$cohort_verified&ctx$BLA_due))ws_csv(nt_weekly(ctx,cutoffs$BLA,TRUE),file.path(out,"weekly_day57_lab_verified.csv"))
  params<-list(version="1.0.0",root=root,EUA=eua_cutoff,BLA=bla_cutoff,day57_anchor=day57_anchor,
    day57_early=day57_early,day57_late=day57_late,baseline_lookback=baseline_lookback,min_randomized=min_randomized,
    numeric_date_origin=numeric_date_origin,created_utc=format(Sys.time(),tz="UTC",usetz=TRUE))
  ws_csv(data.frame(parameter=names(params),value=as.character(unlist(params))),file.path(out,"parameters.csv"))
  capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
  nt_plot(w,out,png_dpi,cutoffs$BLA,day57_anchor)
  nt_report(summary,verified,allrand,ctx,w,out,cutoffs,day57_anchor,day57_early,day57_late,baseline_lookback)
  writeLines("COMPLETE",file.path(out,"RUN_STATUS.txt"));complete<-TRUE
  message("Complete. Randomized: ",ws_n(nrow(ctx)),"; baseline negative: ",ws_n(sum(ctx$cohort)),"; laboratory-confirmed subset: ",ws_n(sum(ctx$cohort_verified)))
  message("Report: ",file.path(out,"anti_n_testing_report.html"))
  print(nt_display(summary),row.names=FALSE,right=FALSE)
  if(isTRUE(show_report))try(utils::browseURL(file.path(out,"anti_n_testing_report.html")),silent=TRUE)
  invisible(list(output_dir=out,summary=summary,weekly=w,subjects=ctx))
}

# ---- RUN: executes when the whole script is selected and run, or sourced ----
if (isTRUE(getOption("phmpt.autorun", TRUE))) {
  analyse_anti_n_testing()
}
