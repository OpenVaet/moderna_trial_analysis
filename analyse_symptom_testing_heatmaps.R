# =============================================================================
# Moderna P301 - weekly symptom-to-visit-to-PCR heatmaps | v1.0.4 | 2026-10-21
# Stand-alone: save in the project root, select all and run in RStudio.
# CSV/JMP: base R. For XPT/SAS7BDAT input only: install.packages("haven").
# First run analyse_background_covid.R to create analysis/background_covid/
# site_weekly.csv. This script reads that asset; it does not download it again.
#
# PRIMARY UNIT: a reconstructed symptomatic EPISODE, not each symptom entry.
# AE symptoms and positive COVID-form symptoms are combined before testing is
# inspected. Overlapping intervals / gaps <= 7 unrecorded days are one episode.
# Episodes are assigned to their FIRST recorded symptom week, once only.
# Default denominator: incident episodes with follow-up through onset day + 7.
# Numerator: a documented illness-visit association and a diagnostic PCR result
# collected on/after that visit, <=3 days later, and <=7 days after episode onset.
# Day 0 and day 7 are included. These are analyst calendar-day linkage rules,
# not a verified application of protocol trigger/duration or 72-hour rules.
#
# RAW visit evidence: SV actual illness-visit date, FAOT physical-examination
# assessment at an illness visit, or an ORIGINAL MB illness-visit label/date.
# ADMB AVISIT and ADMB VISIT are NOT visit evidence: analysis programming can
# derive illness labels for local tests with blank original visit fields.
# INDEPENDENT mode requires SV/FAOT, not a visit label on the PCR record itself.
# RAW mode is a documentation-association measure, not proof of a physical visit
# or that symptoms caused a test. Both endpoints are exported for review.
# Diagnostic PCR = original MB illness label or MBGRPID=003 / ADMB LOCALFL=Y;
# a compatible PCR method and documented positive/negative/indeterminate result
# are required. Routine scheduled swabs alone do not satisfy the numerator.
#
# Use original randomized arm. BLINDED is the main comparison; FULL includes
# post-unblinding/crossover data and is descriptive only. The illustrative PDF
# adds exploratory site tests; these are not evidence of investigator intent.
# To load functions only: options(phmpt.autorun=FALSE) before source().
# =============================================================================

PHMPT_ROOT <- ""                   # auto-detect; or D:/Projects/moderna_trial_analysis
PHMPT_RELEASE <- "auto"            # Prod preferred, otherwise EUA; never pools them
PHMPT_WINDOW <- "BLINDED"           # BLINDED, EUA (through 2020-11-25), or FULL
PHMPT_ADSL_FILE <- ""
PHMPT_AE_FILE <- ""
PHMPT_SYMPTOMS_FILE <- ""           # ADSYMP / ADSY
PHMPT_MB_FILE <- ""                 # required raw SDTM MB
PHMPT_ADMB_FILE <- ""               # optional result reconciliation/baseline flags
PHMPT_SV_FILE <- ""                 # optional independent visit evidence
PHMPT_FAOT_FILE <- ""               # optional independent visit evidence
PHMPT_MAPPING_FILE <- ""            # existing mapping asset, otherwise embedded mapping
PHMPT_BACKGROUND_FILE <- "analysis/background_covid/site_weekly.csv"
PHMPT_EXPECTED_PROD_CUTOFF <- "2021-03-26"
PHMPT_EUA_CUTOFF <- "2020-11-25"
PHMPT_START_DATE <- ""              # optional displayed onset-date limit
PHMPT_END_DATE <- ""                # optional analysis/follow-up cutoff
PHMPT_VISIT_EVIDENCE <- "RAW"       # RAW or INDEPENDENT (SV/FAOT required)
PHMPT_GAP_DAYS <- 7L
PHMPT_ONSET_TO_PCR_DAYS <- 7L
PHMPT_VISIT_TO_PCR_DAYS <- 3L
PHMPT_AE_TIMING <- "interval"       # interval; onset for a sensitivity analysis
PHMPT_SMALL_CELL_N <- 5L            # mark, rather than hide, small denominators
PHMPT_CSV_ENCODING <- "auto"
PHMPT_NUMERIC_DATE_ORIGIN <- ""      # set only from documented date metadata
PHMPT_MIN_RANDOMIZED <- 30000L
PHMPT_OUTPUT_DIR <- ""              # default analysis/symptom_testing_heatmaps_<window>
PHMPT_WRITE_ANOMALY_EXAMPLES <- TRUE # five-page PDF + separate PDF/PNG/SVG figures
PHMPT_WRITE_REPORT_INSERT <- TRUE   # numbered, copy-ready article section
PHMPT_REPORT_START <- 25L           # follows paragraph 24 in the current draft
PHMPT_REFERENCE_START <- 32L        # follows reference 31 in the current draft
PHMPT_ANOMALY_HIGH_POSITIVITY <- 10 # percent; low <10, high >=10
PHMPT_WRITE_PNG <- TRUE
PHMPT_PNG_DPI <- 300L
PHMPT_SHOW_REPORT <- interactive()

hm_join <- function(x) paste(sort(unique(x[!is.na(x)&nzchar(x)])),collapse="; ")
hm_bool <- function(x) ws_norm(x) %in% c("Y","YES","TRUE","1")
hm_empty_manifest <- function() data.frame(domain=character(),release=character(),path=character(),md5=character(),rows=integer(),cutoff=character(),status=character(),detail=character())
hm_read_csv <- function(path) {
  x<-read.csv(path,colClasses="character",check.names=FALSE,fileEncoding="UTF-8-BOM",na.strings=c("","NA"))
  x[]<-lapply(x,ws_text);x
}
hm_require <- function(x,cols,label) {
  n<-setdiff(cols,names(x));if(length(n))stop(label,": missing fields ",paste(n,collapse=", "))
}
hm_dates_out <- function(x) {
  for(k in intersect(c("start","last_symptom","follow_start","follow_end","visit_date","test_date",
      "lo","hi","first_dose","dose2","week","onset_deadline",
      "RANDDT","TR01SDT","DOSE2DT","CUTOFFDT","AP01EDT","UBLDDT","AP02SDT","TR02SDT",
      "EUAVACDT","EOSDT","DTHDT","end_full","end_blinded"),names(x)))x[[k]]<-ws_iso(x[[k]])
  x
}
hm_metadata <- function(out,name,x) ws_csv(x,file.path(out,name))

# Extend the verified source readers below. Preserve raw visit fields separately.
hm_install_specs <- function() {
  base_specs<-ws_specs
  specs<-function(domain) {
    s<-switch(domain,
      MB=list(required=c("STUDYID","USUBJID","MBSEQ","MBTESTCD","MBSTRESC","MBDTC"),
        optional=c("MBTEST","MBORRES","MBSTAT","MBMETHOD","MBGRPID","VISIT","VISITNUM","MBREASND","MBSPEC","MBNAM")),
      ADMB=list(required=c("STUDYID","USUBJID","MBSEQ","PARAMCD","AVALC"),
        optional=c("ADT","MBDTC","ABLFL","LOCALFL","VISIT","AVISIT")),
      SV=list(required=c("STUDYID","USUBJID","VISIT","SVSTDTC"),
        optional=c("SVSEQ","VISITNUM","SVENDTC","SVOCCUR","SVTERM","SVUPDES")),
      FAOT=list(required=c("STUDYID","USUBJID","FASEQ","FATESTCD","FACAT","FADTC","VISIT"),
        optional=c("FAOBJ","FAORRES","FASTRESC","VISITNUM")),base_specs(domain))
    if(domain=="ADSL")s$optional<-unique(c(s$optional,"SITEID","TR02SDT"))
    s
  }
  specs
}

hm_pick_optional <- function(root,release,domain,path,out,encoding,origin,min_n,submission) {
  if(!nzchar(path)) {
    p<-ws_candidates(root,release,domain)
    if(nzchar(submission))p<-p[ws_submission(p)==submission]
    if(!length(p))return(list(x=NULL,manifest=hm_empty_manifest()))
  }
  ws_pick(root,release,domain,path,out,encoding,origin,min_n,submission)
}

hm_tests <- function(mb,admb,origin,out) {
  m<-mb[ws_norm(mb$MBTESTCD)=="SARSCOV2",,drop=FALSE]
  if(!nrow(m))stop("No MBTESTCD=SARSCOV2 records in MB.")
  m$key<-paste(m$subject,m$MBSEQ,sep=":")
  if(any(!nzchar(m$MBSEQ)))stop("SARS-CoV-2 MB records have missing MBSEQ.")
  if(anyDuplicated(m$key))stop("Conflicting/repeated MB subject+MBSEQ. Resolve before analysis.")
  a<-if(is.null(admb))NULL else admb[ws_norm(admb$PARAMCD)=="SARSCOV2",,drop=FALSE]
  if(!is.null(a)) {
    a$key<-paste(a$subject,a$MBSEQ,sep=":")
    # Identical keys with divergent source values are not silently collapsed.
    a<-unique(a)
    if(anyDuplicated(a$key))stop("Repeated ADMB subject+MBSEQ; choose a reconciled source snapshot.")
    extra<-a[!a$key %in% m$key,,drop=FALSE]
    ws_csv(extra,file.path(out,"admb_records_without_raw_mb.csv"))
    j<-match(m$key,a$key)
  } else j<-rep(NA_integer_,nrow(m))
  av<-function(k) if(is.null(a))rep("",nrow(m)) else {v<-a[[k]][j];v[is.na(v)]<-"";v}
  result_class<-function(x) {
    x<-ws_norm(x);out<-rep("UNRESOLVED",length(x))
    out[x %in% c("POSITIVE","DETECTED","RNA DETECTED","PRESUMPTIVE POSITIVE","Y","YES","POS","DETECTABLE")]<-"POSITIVE"
    out[x %in% c("NEGATIVE","NOT DETECTED","RNA NOT DETECTED","NONREACTIVE","NON-REACTIVE","N","NO","NEG")]<-"NEGATIVE"
    out[x %in% c("INDETERMINATE","INCONCLUSIVE","INVALID","EQUIVOCAL")]<-"INDETERMINATE"
    out
  }
  mr<-result_class(m$MBSTRESC);ar<-result_class(av("AVALC"));result<-ifelse(mr!="UNRESOLVED",mr,ar)
  raw_date<-ws_dates(m$MBDTC,"MB.MBDTC",origin)
  ad_raw<-ws_dates(av("MBDTC"),"ADMB.MBDTC",origin)
  ad_date<-ws_dates(av("ADT"),"ADMB.ADT",origin)
  conflict_date<-(is.finite(raw_date)&is.finite(ad_raw)&raw_date!=ad_raw)|
    (is.finite(raw_date)&is.finite(ad_date)&raw_date!=ad_date)
  # No ADT fallback: partial/raw missing collection dates are not imputed.
  method<-ws_norm(m$MBMETHOD)
  pcr_method<-grepl("PCR|POLYMERASE",method)
  contradictory_method<-grepl("ANTIGEN|SEROLOG|ANTIBODY",method)
  not_done<-ws_norm(m$MBSTAT) %in% c("NOT DONE","NOT PERFORMED","CANCELLED","CANCELED")
  illness<-grepl("^ILLNESS VISIT",ws_norm(m$VISIT))
  group<-sub("^0+","",trimws(m$MBGRPID));local<-group=="3"|hm_bool(av("LOCALFL"))
  baseline<-hm_bool(av("ABLFL"))|grepl("^VISIT 1 DAY 1($| )",ws_norm(m$VISIT))
  usable<-is.finite(raw_date)&!conflict_date&!not_done&pcr_method&!contradictory_method&result!="UNRESOLVED"
  p<-data.frame(test_id=m$key,subject=m$subject,MBSEQ=m$MBSEQ,test_date=raw_date,
    collection_text=m$MBDTC,raw_visit=m$VISIT,raw_visit_number=m$VISITNUM,
    MBGRPID=m$MBGRPID,local_coded=local,original_illness_label=illness,
    baseline=baseline,method=m$MBMETHOD,pcr_method=pcr_method,method_conflict=contradictory_method,
    raw_result=m$MBSTRESC,analysis_result=av("AVALC"),result=result,
    result_conflict=mr!="UNRESOLVED"&ar!="UNRESOLVED"&mr!=ar,
    not_done=not_done,date_conflict=conflict_date,usable=usable,
    diagnostic=illness|local,indeterminate=result=="INDETERMINATE",stringsAsFactors=FALSE)
  p$usable_diagnostic<-p$usable&p$diagnostic&!p$baseline
  # A positive/negative disagreement does not erase evidence that an assay ran;
  # its result remains flagged. This task estimates testing, not positivity.
  ws_csv(hm_dates_out(p),file.path(out,"pcr_records_reconciled.csv"))
  inv<-as.data.frame(table(raw_result=m$MBSTRESC,analysis_result=av("AVALC"),status=m$MBSTAT),stringsAsFactors=FALSE)
  ws_csv(inv[inv$Freq>0,,drop=FALSE],file.path(out,"pcr_result_inventory.csv"))
  p
}

hm_visits <- function(sv,faot,pcr,origin,out) {
  blank<-data.frame(subject=character(),visit_date=numeric(),source=character(),source_ref=character(),label=character(),assessment=character(),independent=logical())
  rows<-list(blank)
  if(!is.null(sv)) {
    keep<-grepl("^ILLNESS VISIT",ws_norm(sv$VISIT))&!ws_norm(sv$SVOCCUR)%in%c("N","NO")
    rows[[length(rows)+1L]]<-data.frame(subject=sv$subject[keep],visit_date=ws_dates(sv$SVSTDTC[keep],"SV.SVSTDTC",origin),
      source=rep("SV",sum(keep)),source_ref=if(any(keep))paste0("SV:",ifelse(nzchar(sv$SVSEQ[keep]),sv$SVSEQ[keep],which(keep))) else character(),
      label=sv$VISIT[keep],assessment=sv$SVOCCUR[keep],independent=rep(TRUE,sum(keep)),stringsAsFactors=FALSE)
  }
  if(!is.null(faot)) {
    keep<-grepl("^ILLNESS VISIT",ws_norm(faot$VISIT))&ws_norm(faot$FACAT)=="PHYSICAL EXAMINATION"&ws_norm(faot$FATESTCD)=="ASSESS"
    rows[[length(rows)+1L]]<-data.frame(subject=faot$subject[keep],visit_date=ws_dates(faot$FADTC[keep],"FAOT.FADTC",origin),
      source=rep("FAOT",sum(keep)),source_ref=if(any(keep))paste0("FAOT:",faot$FASEQ[keep]) else character(),label=faot$VISIT[keep],
      assessment=faot$FASTRESC[keep],independent=rep(TRUE,sum(keep)),stringsAsFactors=FALSE)
  }
  q<-pcr[pcr$original_illness_label,,drop=FALSE]
  if(nrow(q))rows[[length(rows)+1L]]<-data.frame(subject=q$subject,visit_date=q$test_date,source="MB_ORIGINAL_VISIT",
    source_ref=paste0("MB:",q$MBSEQ),label=q$raw_visit,assessment=ifelse(q$not_done,"PCR NOT DONE",q$result),
    independent=FALSE,stringsAsFactors=FALSE)
  raw<-do.call(rbind,rows)
  ws_csv(hm_dates_out(raw),file.path(out,"illness_visit_source_evidence.csv"))
  raw<-raw[is.finite(raw$visit_date),,drop=FALSE]
  if(!nrow(raw))return(data.frame(visit_id=character(),subject=character(),visit_date=numeric(),sources=character(),source_refs=character(),labels=character(),independent=logical()))
  groups<-split(seq_len(nrow(raw)),paste(raw$subject,raw$visit_date,sep="|"))
  v<-do.call(rbind,lapply(groups,function(ix) {
    z<-raw[ix,];data.frame(subject=z$subject[1],visit_date=z$visit_date[1],sources=hm_join(z$source),
      source_refs=hm_join(paste0(z$source,"/",z$source_ref)),labels=hm_join(z$label),independent=any(z$independent))
  }))
  rownames(v)<-NULL;v$visit_id<-sprintf("V%07d",seq_len(nrow(v)))
  v<-v[c("visit_id","subject","visit_date","sources","source_refs","labels","independent")]
  ws_csv(hm_dates_out(v),file.path(out,"illness_visits.csv"));v
}

hm_episodes <- function(events,ctx,gap,horizon,display_start,display_end,out) {
  p<-ctx[ctx$DOSE1FL=="Y"&is.finite(ctx$TR01SDT)&ctx$follow_end>=ctx$TR01SDT,,drop=FALSE]
  j<-match(events$subject,p$subject)
  e<-events[!is.na(j),,drop=FALSE];j<-j[!is.na(j)]
  good<-is.finite(e$lo)&is.finite(e$hi)&e$hi>=e$lo&e$lo<=p$follow_end[j]
  bad<-e[!good,,drop=FALSE];bad$reason<-ifelse(!is.finite(bad$lo)|!is.finite(bad$hi),"UNKNOWN_DATE","INVALID_OR_AFTER_FOLLOWUP")
  ws_csv(hm_dates_out(bad),file.path(out,"excluded_symptom_records.csv"))
  e<-e[good,,drop=FALSE];j<-j[good];e$hi<-pmin(e$hi,p$follow_end[j])
  e<-e[order(e$subject,e$lo,e$hi,e$source,e$source_id,e$symptom_group),,drop=FALSE]
  if(!nrow(e))stop("No dated positive symptoms/AE matches within selected subject follow-up.")
  bysub<-split(seq_len(nrow(e)),e$subject);block<-integer(nrow(e));nblock<-0L
  for(ix in bysub) {
    last<-c(-Inf,head(cummax(e$hi[ix]),-1L));g<-cumsum(e$lo[ix]>last+gap+1L)
    block[ix]<-g+nblock;nblock<-nblock+max(g)
  }
  e$episode_id<-sprintf("E%07d",block)
  groups<-split(seq_len(nrow(e)),e$episode_id)
  ep<-do.call(rbind,lapply(groups,function(ix) {
    z<-e[ix,];onset<-min(z$lo);on<-z$lo==onset
    data.frame(episode_id=z$episode_id[1],subject=z$subject[1],start=onset,last_symptom=max(z$hi),
      source_class=hm_join(z$source),onset_source=hm_join(z$source[on]),
      symptom_groups=hm_join(z$symptom_group),onset_groups=hm_join(z$symptom_group[on]),
      n_groups=length(unique(z$symptom_group)),n_onset_groups=length(unique(z$symptom_group[on])),
      has_ae=any(z$source=="AE"),has_covid_form=any(z$source=="COVID_FORM"),
      unknown_ae_end=any(z$date_quality=="START_ONLY_END_UNKNOWN"),source_records=nrow(z),stringsAsFactors=FALSE)
  }))
  rownames(ep)<-NULL;ix<-match(ep$subject,p$subject)
  ep$SITEID<-p$SITEID[ix];ep$arm<-p$arm[ix];ep$first_dose<-p$TR01SDT[ix];ep$dose2<-p$DOSE2DT[ix]
  ep$follow_start<-ep$first_dose;ep$follow_end<-p$follow_end[ix];ep$week<-ws_week(ep$start)
  ep$span_days<-ep$last_symptom-ep$start+1L;ep$preexisting<-ep$start<ep$follow_start
  ep$onset_in_display<-ep$start>=display_start&ep$start<=display_end
  ep$early_postdose<-ep$start>=ep$first_dose&ep$start<=ep$first_dose+6L |
    is.finite(ep$dose2)&ep$start>=ep$dose2&ep$start<=ep$dose2+6L
  ep$complete_followup<-ep$start+horizon<=ep$follow_end
  # Short-gap sensitivity must not allow the same test to satisfy two episodes.
  ep<-ep[order(ep$subject,ep$start),,drop=FALSE]
  nxt<-c(ep$start[-1L],Inf);same<-c(ep$subject[-1L]==head(ep$subject,-1L),FALSE)
  ep$next_episode_within_horizon<-same & nxt<=ep$start+horizon
  ep$eligible<-!ep$preexisting&ep$onset_in_display&ep$complete_followup&!ep$next_episode_within_horizon
  ep$exclusion_reason<-ifelse(ep$preexisting,"PREEXISTING",ifelse(!ep$onset_in_display,"OUTSIDE_DISPLAY",
    ifelse(!ep$complete_followup,"INCOMPLETE_FOLLOWUP",ifelse(ep$next_episode_within_horizon,"OVERLAPPING_TEST_WINDOW",""))))
  ws_csv(hm_dates_out(e),file.path(out,"episode_symptom_source_links.csv"))
  ep
}

hm_link <- function(ep,pcr,visits,horizon,visit_delay,visit_mode,independent_available,out) {
  pi<-split(seq_len(nrow(pcr)),pcr$subject);vi<-split(seq_len(nrow(visits)),visits$subject)
  logicals<-c("visit_by_horizon","independent_visit_by_horizon","diagnostic_pcr_by_horizon",
    "any_nonbaseline_pcr_by_horizon","raw_visit_pcr","independent_visit_pcr","tested",
    "indeterminate_test_in_chain","not_done_in_window")
  for(k in logicals)ep[[k]]<-FALSE
  ep$first_pcr_delay<-ep$first_visit_delay<-NA_real_;ep$n_matching_tests<-0L
  ep$linked_visit_ids<-ep$linked_test_ids<-ep$link_evidence<-""
  links<-list();vl<-list()
  for(i in seq_len(nrow(ep))) {
    if(ep$preexisting[i]||!ep$onset_in_display[i])next
    lower<-ep$start[i];upper<-min(lower+horizon,ep$follow_end[i]);s<-ep$subject[i]
    v<-visits[vi[[s]],,drop=FALSE];t<-pcr[pi[[s]],,drop=FALSE]
    v<-v[is.finite(v$visit_date)&v$visit_date>=lower&v$visit_date<=upper,,drop=FALSE]
    t<-t[is.finite(t$test_date)&t$test_date>=lower&t$test_date<=upper,,drop=FALSE]
    ep$visit_by_horizon[i]<-nrow(v)>0;ep$independent_visit_by_horizon[i]<-any(v$independent)
    ep$diagnostic_pcr_by_horizon[i]<-any(t$usable_diagnostic)
    ep$any_nonbaseline_pcr_by_horizon[i]<-any(t$usable&!t$baseline)
    ep$not_done_in_window[i]<-any(t$not_done)
    if(nrow(v)) {v$episode_id<-ep$episode_id[i];v$delay_from_onset<-v$visit_date-lower;vl[[length(vl)+1L]]<-v}
    hit<-list()
    for(k in which(t$usable_diagnostic)) {
      vv<-v[v$visit_date<=t$test_date[k]&t$test_date[k]-v$visit_date<=visit_delay,,drop=FALSE]
      if(!nrow(vv))next
      hit[[length(hit)+1L]]<-data.frame(episode_id=ep$episode_id[i],subject=s,
        visit_id=vv$visit_id,visit_date=vv$visit_date,visit_sources=vv$sources,independent_visit=vv$independent,
        test_id=t$test_id[k],test_date=t$test_date[k],test_result=t$result[k],
        result_conflict=t$result_conflict[k],test_raw_visit=t$raw_visit[k],MBGRPID=t$MBGRPID[k],
        onset_to_visit=vv$visit_date-lower,onset_to_test=t$test_date[k]-lower,
        visit_to_test=t$test_date[k]-vv$visit_date,stringsAsFactors=FALSE)
    }
    if(length(hit)) {
      h<-do.call(rbind,hit);h$episode_eligible<-ep$eligible[i];links[[length(links)+1L]]<-h
      ep$raw_visit_pcr[i]<-TRUE;ep$independent_visit_pcr[i]<-any(h$independent_visit)
      selected<-if(visit_mode=="INDEPENDENT")h[h$independent_visit,,drop=FALSE] else h
      if(nrow(selected)) {
        ep$tested[i]<-TRUE;ep$first_pcr_delay[i]<-min(selected$onset_to_test)
        ep$first_visit_delay[i]<-min(selected$onset_to_visit);ep$n_matching_tests[i]<-length(unique(selected$test_id))
        ep$linked_visit_ids[i]<-hm_join(selected$visit_id);ep$linked_test_ids[i]<-hm_join(selected$test_id)
        ep$link_evidence[i]<-hm_join(selected$visit_sources)
        ep$indeterminate_test_in_chain[i]<-any(selected$test_result=="INDETERMINATE")
      }
    }
  }
  blank<-data.frame(episode_id=character(),subject=character(),visit_id=character(),visit_date=numeric(),visit_sources=character(),
    independent_visit=logical(),test_id=character(),test_date=numeric(),test_result=character(),result_conflict=logical(),
    test_raw_visit=character(),MBGRPID=character(),onset_to_visit=numeric(),onset_to_test=numeric(),visit_to_test=numeric(),episode_eligible=logical())
  h<-if(length(links))do.call(rbind,links) else blank
  ws_csv(hm_dates_out(h),file.path(out,"episode_visit_pcr_links.csv"))
  if(length(vl))ws_csv(hm_dates_out(do.call(rbind,vl)),file.path(out,"episode_visit_links.csv"))
  if(!independent_available)ep$independent_visit_pcr<-ep$independent_visit_by_horizon<-NA
  if(any(ep$tested&!ep$diagnostic_pcr_by_horizon)||any(ep$tested&!ep$visit_by_horizon))stop("Visit/PCR linkage invariant failed.")
  ep
}

hm_background <- function(path,sites,weeks,out) {
  if(!file.exists(path))stop("Background asset missing: ",path,". Run analyse_background_covid.R first, or set PHMPT_BACKGROUND_FILE.")
  b<-hm_read_csv(path)
  required<-c("SITEID","week_start","county_fips","county_name","state_abbr","state_name",
              "county_cases","county_positivity_pct","state_cases","state_positivity_pct")
  hm_require(b,required,"Background site_weekly.csv")
  b<-unique(b);if(anyDuplicated(paste(b$SITEID,b$week_start)))stop("Background asset has conflicting site/week duplicates.")
  b$week<-ws_dates(b$week_start,"background.week_start",strict=TRUE)
  if(any(b$week!=ws_week(b$week)))stop("Background data must use Monday week starts.")
  for(k in c("county_cases","county_positivity_pct","state_cases","state_positivity_pct"))b[[k]]<-suppressWarnings(as.numeric(b[[k]]))
  if(any(b$county_positivity_pct<0|b$county_positivity_pct>100|b$state_positivity_pct<0|b$state_positivity_pct>100,na.rm=TRUE))stop("Background positivity outside 0-100%.")
  outrows<-list();cross<-list()
  for(i in seq_along(sites)) {
    s<-sites[i];z<-b[b$SITEID==s,,drop=FALSE]
    geo<-unique(z[c("county_fips","county_name","state_abbr","state_name")])
    if(nrow(geo)>1L)stop("Changing/conflicting background geography for site ",s)
    county<-nrow(geo)==1L && nzchar(geo$county_fips) && any(is.finite(z$county_positivity_pct[z$week %in% weeks]))
    state<-!county && nrow(geo)==1L && nzchar(geo$state_abbr) && any(is.finite(z$state_positivity_pct[z$week %in% weeks]))
    if(county) {level<-"county";id<-geo$county_fips;name<-paste(geo$county_name,geo$state_name,sep=", ")}
    else if(state) {level<-"state";id<-geo$state_abbr;name<-paste0(geo$state_name," (state fallback)")}
    else if(nrow(geo)==1L && nzchar(geo$county_fips)) {level<-"county";id<-geo$county_fips;name<-paste(geo$county_name,geo$state_name,sep=", ")}
    else {level<-"missing";id<-s;name<-paste(s,"- background unavailable")}
    m<-match(weeks,z$week);cases<-pct<-rep(NA_real_,length(weeks))
    if(level!="missing") {
      cases<-z[[paste0(level,"_cases")]][m];pct<-z[[paste0(level,"_positivity_pct")]][m]
    }
    cross[[i]]<-data.frame(SITEID=s,region_key=paste(level,id,sep="_"),region_level=level,region_id=id,region_name=name,
      weeks_with_cases=sum(is.finite(cases)),weeks_with_positivity=sum(is.finite(pct)),weeks=length(weeks))
    outrows[[i]]<-data.frame(SITEID=s,region_key=paste(level,id,sep="_"),week=weeks,cases=cases,positivity=pct)
  }
  map<-do.call(rbind,cross);d<-do.call(rbind,outrows)
  # Shared-county data are reused, never summed/averaged across sites.
  area<-unique(d[c("region_key","week","cases","positivity")])
  if(anyDuplicated(area[c("region_key","week")]))stop("Sites sharing a region have inconsistent background values.")
  ws_csv(map,file.path(out,"site_region_coverage.csv"));ws_csv(hm_dates_out(d),file.path(out,"site_background_used.csv"))
  list(map=map,area=area)
}

hm_cells <- function(ep,ctx,sites,weeks,horizon,small_n,mode) {
  x<-expand.grid(SITEID=sites,arm=c("mRNA-1273","Placebo"),week=weeks,KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
  count_fields<-c("symptom_episodes","tested_episodes","visit_episodes","diagnostic_pcr_episodes","any_pcr_episodes",
    "raw_visit_pcr_episodes","independent_visit_pcr_episodes","ae_only_episodes","covid_only_episodes","both_sources_episodes",
    "early_postdose_episodes","early_postdose_tested","outside_postdose_episodes","outside_postdose_tested",
    "incomplete_episodes","preexisting_episodes","overlapping_window_episodes","subjects_with_opportunity","subjects_observed")
  for(k in count_fields)x[[k]]<-0L
  key<-paste(x$SITEID,x$arm,x$week,sep="|")
  ix<-match(paste(ep$SITEID,ep$arm,ep$week,sep="|"),key)
  count<-function(flag)tabulate(ix[which(flag)],nrow(x))
  good<-ep$eligible
  x$symptom_episodes<-count(good);x$tested_episodes<-count(good&ep$tested)
  x$visit_episodes<-count(good&ep$visit_by_horizon)
  x$diagnostic_pcr_episodes<-count(good&ep$diagnostic_pcr_by_horizon)
  x$any_pcr_episodes<-count(good&ep$any_nonbaseline_pcr_by_horizon)
  x$raw_visit_pcr_episodes<-count(good&ep$raw_visit_pcr)
  x$independent_visit_pcr_episodes<-if(all(is.na(ep$independent_visit_pcr)))NA_integer_ else count(good&ep$independent_visit_pcr)
  x$ae_only_episodes<-count(good&ep$has_ae&!ep$has_covid_form)
  x$covid_only_episodes<-count(good&!ep$has_ae&ep$has_covid_form)
  x$both_sources_episodes<-count(good&ep$has_ae&ep$has_covid_form)
  x$early_postdose_episodes<-count(good&ep$early_postdose);x$early_postdose_tested<-count(good&ep$early_postdose&ep$tested)
  x$outside_postdose_episodes<-count(good&!ep$early_postdose);x$outside_postdose_tested<-count(good&!ep$early_postdose&ep$tested)
  x$incomplete_episodes<-count(!ep$preexisting&ep$onset_in_display&!ep$complete_followup)
  x$preexisting_episodes<-count(ep$preexisting&ep$onset_in_display)
  x$overlapping_window_episodes<-count(!ep$preexisting&ep$onset_in_display&ep$complete_followup&ep$next_episode_within_horizon)
  for(s in sites)for(a in c("mRNA-1273","Placebo")) {
    rows<-which(x$SITEID==s&x$arm==a);p<-ctx[ctx$SITEID==s&ctx$arm==a&ctx$DOSE1FL=="Y",]
    for(r in rows) {
      x$subjects_observed[r]<-sum(is.finite(p$TR01SDT)&p$TR01SDT<=x$week[r]+6&p$follow_end>=x$week[r])
      x$subjects_with_opportunity[r]<-sum(is.finite(p$TR01SDT)&p$TR01SDT<=x$week[r]+6&
        p$follow_end-horizon>=x$week[r]&p$follow_end-horizon>=p$TR01SDT)
    }
  }
  x$testing_pct<-ifelse(x$symptom_episodes>0,100*x$tested_episodes/x$symptom_episodes,NA_real_)
  x$independent_visit_pcr_pct<-ifelse(x$symptom_episodes>0,100*x$independent_visit_pcr_episodes/x$symptom_episodes,NA_real_)
  x$outside_postdose_pct<-ifelse(x$outside_postdose_episodes>0,100*x$outside_postdose_tested/x$outside_postdose_episodes,NA_real_)
  x$small_denominator<-x$symptom_episodes>0&x$symptom_episodes<small_n
  x$cell_status<-ifelse(x$symptom_episodes>0,"ESTIMABLE",ifelse(x$incomplete_episodes>0,"CENSORED_ONLY",
    ifelse(x$subjects_with_opportunity==0,"NO_OBSERVATION_OPPORTUNITY","NO_ELIGIBLE_SYMPTOM_EPISODE")))
  x$window<-mode;x$week_start<-ws_iso(x$week);x$week_end<-ws_iso(x$week+6)
  stopifnot(sum(x$symptom_episodes)==sum(ep$eligible),sum(x$tested_episodes)==sum(ep$eligible&ep$tested),
    all(x$tested_episodes<=x$symptom_episodes),all(x$ae_only_episodes+x$covid_only_episodes+x$both_sources_episodes==x$symptom_episodes))
  x
}

# Fixed red-to-green scale for BOTH arms; missing/zero-denominator cells are gray.
hm_palette <- grDevices::colorRampPalette(c("#A50026","#D73027","#F46D43","#FEE08B","#D9EF8B","#66BD63","#1A9850","#006837"))(101)
hm_ticks <- function(weeks) {
  dates<-as.Date(weeks,origin="1970-01-01");pick<-which(!duplicated(format(dates,"%Y-%m")))
  if(length(pick)>1L&&pick[2]-pick[1]<3L)pick<-pick[-1L]
  pick<-pick[seq(1,length(pick),by=max(1,ceiling(length(pick)/12)))]
  list(at=pick,label=paste(month.abb[as.integer(format(dates[pick],"%m"))],format(dates[pick],"%Y"),sep="\n"))
}
hm_heat <- function(cells,sites,weeks,arm,labels=TRUE,main="",overview=FALSE) {
  z<-cells[cells$arm==arm&cells$SITEID %in% sites,,drop=FALSE];n<-length(sites);nw<-length(weeks)
  par(mar=if(overview)c(4.7,5.2,3.0,4.3) else c(2.6,5.2,0.3,4.3),mgp=c(2.5,.4,0),tcl=-.18,
      family="sans",las=1,xaxs="i",yaxs="i")
  plot.new();plot.window(xlim=c(.5,nw+.5),ylim=c(.5,n+.5))
  xx<-match(z$week,weeks);yy<-n+1L-match(z$SITEID,sites)
  fill<-rep("#ECEFF2",nrow(z));fill[z$cell_status=="NO_OBSERVATION_OPPORTUNITY"]<-"#B6BEC7"
  ok<-is.finite(z$testing_pct);fill[ok]<-hm_palette[pmax(1L,pmin(101L,round(z$testing_pct[ok])+1L))]
  rect(xx-.49,yy-.47,xx+.49,yy+.47,col=fill,border="white",lwd=.3)
  cens<-z$cell_status=="CENSORED_ONLY"
  if(any(cens))rect(xx[cens]-.49,yy[cens]-.47,xx[cens]+.49,yy[cens]+.47,col="#77828D",density=18,angle=45,border=NA)
  dark<-ok&(z$testing_pct<16|z$testing_pct>85);ink<-ifelse(dark,"white","#243444")
  if(labels) {
    lab<-ifelse(ok,paste0(z$tested_episodes,"/",z$symptom_episodes,ifelse(z$small_denominator,"*","")),
      ifelse(cens,"C",ifelse(z$cell_status=="NO_ELIGIBLE_SYMPTOM_EPISODE","-","")))
    text(xx,yy,labels=lab,col=ink,cex=max(.38,min(.65,24/nw)),font=ifelse(dark,2,1))
  } else if(any(z$small_denominator))points(xx[z$small_denominator],yy[z$small_denominator],pch=16,cex=.24,col=ink[z$small_denominator])
  axis(2,at=seq_len(n),labels=rev(sites),tick=FALSE,cex.axis=if(overview).65 else .77,line=.1)
  ticks<-hm_ticks(weeks);axis(1,at=ticks$at,labels=ticks$label,cex.axis=if(overview).75 else .67)
  if(nzchar(main))mtext(main,side=3,line=1.2,adj=0,font=2,cex=1.08)
}
hm_bg_panel <- function(background,weeks,arm,level) {
  b<-background[match(weeks,background$week),];x<-seq_along(weeks)
  cy<-range(pretty(range(c(0,b$cases[is.finite(b$cases)]))));if(diff(cy)==0)cy<-c(0,1)
  py<-c(0,max(5,ceiling(max(c(0,b$positivity[is.finite(b$positivity)]))/5)*5))
  par(mar=c(.25,5.2,1.9,4.3),mgp=c(2.7,.4,0),tcl=-.2,family="sans",las=1,xaxs="i",yaxs="i")
  plot(x,b$cases,type="n",xlim=c(.5,length(weeks)+.5),ylim=cy,axes=FALSE,xlab="",ylab="")
  abline(h=pretty(cy),col="#E4E9EE",lwd=.5)
  rect(x-.4,pmin(0,b$cases),x+.4,pmax(0,b$cases),col="#BAC8D8",border=NA)
  axis(2,at=pretty(cy),labels=format(pretty(cy),trim=TRUE,big.mark=",",scientific=FALSE),cex.axis=.68)
  mtext("Reported cases",2,line=3.7,las=0,cex=.7)
  par(new=TRUE);plot(x,b$positivity,type="n",xlim=c(.5,length(weeks)+.5),ylim=py,axes=FALSE,xlab="",ylab="")
  lines(x,b$positivity,col="#AC302C",lwd=1.6)
  axis(4,at=pretty(py),labels=paste0(pretty(py),"%"),col.axis="#AC302C",cex.axis=.68)
  mtext("Community positivity",4,line=2.8,las=0,col="#AC302C",cex=.7)
  mtext(paste0(arm,"  |  Local background above; symptom-to-visit-to-PCR rate below"),side=3,line=.65,adj=0,font=2,cex=.83)
  if(!any(is.finite(b$positivity)))text(mean(x),py[2]*.55,"Background positivity unavailable",col="#67737E",cex=.78)
  if(!any(is.finite(b$cases)))mtext("Reported cases unavailable",3,line=-1.2,adj=.02,cex=.72,col="#67737E")
}
hm_legend <- function(small_n) {
  par(mar=c(.2,5.2,.25,4.3),family="sans");plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
  for(i in 0:100)rect(.01+.38*i/101,.45,.01+.38*(i+1)/101,.75,col=hm_palette[i+1],border=NA)
  text(.01,.23,"0%",adj=c(0,.5),cex=.6)
  text(seq(.105,.295,length.out=3),.23,c("25%","50%","75%"),cex=.6)
  text(.39,.23,"100%",adj=c(1,.5),cex=.6)
  text(.20,.98,"Documented visit + PCR / symptomatic episodes",adj=c(.5,1),cex=.63)
  rect(.43,.49,.45,.75,col="#ECEFF2",border="white");text(.46,.62,"No eligible episode",adj=0,cex=.61)
  rect(.68,.49,.70,.75,col="#B6BEC7",border="white");text(.71,.62,"No observation opportunity",adj=0,cex=.61)
  text(.43,.16,paste0("n/N = tested / eligible; * (or dot) means N < ",small_n,". C = censored-only episodes."),adj=0,cex=.61)
}
hm_region_page <- function(cells,map,bg,weeks,region,window,visit_mode,small_n,horizon,visit_delay) {
  s<-sort(map$SITEID[map$region_key==region]);name<-map$region_name[match(region,map$region_key)]
  b<-bg[bg$region_key==region,,drop=FALSE];n<-length(s)
  layout(matrix(1:7,ncol=1),heights=c(.5,1.65,max(.9,.30*n+.55),1.65,max(.9,.30*n+.55),.75,.65))
  par(cex=1)
  par(mar=c(0,5.2,.15,4.3),family="sans");plot.new()
  text(0,.70,name,adj=0,cex=1.05,font=2)
  text(0,.18,paste0(window," | ",ws_iso(min(weeks))," to ",ws_iso(max(weeks)+6),
    " | ",if(visit_mode=="RAW")"Original visit-record association" else "Independent SV/FAOT visit evidence"),adj=0,cex=.72,col="#4D5966")
  for(a in c("mRNA-1273","Placebo")) {
    hm_bg_panel(b,weeks,a,map$region_level[match(region,map$region_key)])
    hm_heat(cells,s,weeks,a,labels=TRUE)
  }
  hm_legend(small_n)
  par(mar=c(0,5.2,0,4.3));plot.new()
  note<-paste0("Episodes assigned to onset week. PCR: days 0-",horizon," after onset; same day to +",visit_delay,
    " after recorded visit. Complete follow-up required.")
  text(0,.86,note,adj=0,cex=.6)
  text(0,.55,"AE and COVID-form symptoms included. Fixed 0-100% colors for both arms; unadjusted descriptive rates. Gray is not 0%.",adj=0,cex=.6)
  text(0,.24,"Community series are identical for both arms: case counts (columns), test positivity (line). Sources/limitations and all cell counts are in the report.",adj=0,cex=.59)
}
hm_overview_page <- function(cells,sites,weeks,arm,window,small_n) {
  layout(matrix(1:3,ncol=1),heights=c(max(8,length(sites)*.19),.75,.65))
  par(cex=1)
  hm_heat(cells,sites,weeks,arm,labels=FALSE,main=paste(arm,"-",window,"- weekly symptom-to-visit-to-PCR rate"),overview=TRUE)
  hm_legend(small_n);par(mar=c(0,5.2,0,4.3));plot.new()
  text(0,.72,"Rows are trial sites; columns are onset weeks. Both arms use the same site order and fixed color scale. No regional percentages are pooled.",adj=0,cex=.68)
  text(0,.28,"Use the county/state pages for local background waves and n/N cell labels. Missing observation opportunity and zero symptoms are gray, not red.",adj=0,cex=.68)
}
hm_pdf <- function(path,width,height) {
  if(capabilities("cairo"))grDevices::cairo_pdf(path,width=width,height=height,onefile=TRUE,family="sans",fallback_resolution=300)
  else grDevices::pdf(path,width=width,height=height,onefile=TRUE,useDingbats=FALSE)
}
hm_png <- function(path,width,height,dpi,draw) {
  args<-list(filename=path,width=width,height=height,units="in",res=dpi)
  if(capabilities("cairo"))args$type<-"cairo"
  do.call(grDevices::png,args);tryCatch(draw(),finally=grDevices::dev.off())
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
  ADSL=list(required=c("STUDYID","USUBJID","RANDFL","TRT01P","RANDDT","DOSE1FL","DOSE2FL","TR01SDT","DOSE2DT","CUTOFFDT"),
    optional=c("SUBJID","AP01EDT","UBLDDT","AP02SDT","TR02SDT","EUAVACDT","EOSDT","DTHDT")),
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
  missing<-setdiff(spec$required,names(x));if(length(missing))stop("Missing ",domain," fields: ",paste(missing,collapse=", "))
  x[]<-lapply(x,ws_text)
  x<-x[ws_norm(x$STUDYID)=="MRNA-1273-P301",,drop=FALSE]
  if(!nrow(x))stop("No P301 rows in this table.")
  for(n in setdiff(spec$optional,names(x)))x[[n]]<-rep("",nrow(x))
  x<-x[c(spec$required,spec$optional)];x$subject<-ws_sid(x$USUBJID)
  for(n in intersect(c("AESEQ","ASEQ","SRCSEQ","MBSEQ","SVSEQ","FASEQ"),names(x))) {
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
  p[order(match(tolower(tools::file_ext(p)),c("csv","jmp","xpt","sas7bdat")),p)]
}
ws_submission <- function(p) {
  z<-tolower(basename(p));m<-regexpr("[0-9]+_s[0-9]+_m[0-9]+_mrna-1273-p301",z)
  out<-rep("",length(p));ok<-m>0
  out[ok]<-substring(z[ok],m[ok],m[ok]+attr(m,"match.length")[ok]-1L)
  out
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


hm_html <- function(x) {
  x<-gsub("&","&amp;",as.character(x),fixed=TRUE);x<-gsub("<","&lt;",x,fixed=TRUE)
  x<-gsub(">","&gt;",x,fixed=TRUE);gsub('"',"&quot;",x,fixed=TRUE)
}
hm_table <- function(d) {
  if(!nrow(d))return("<p>No records.</p>")
  paste0("<table><thead><tr>",paste0("<th>",hm_html(names(d)),"</th>",collapse=""),"</tr></thead><tbody>",
    paste(apply(d,1,function(r)paste0("<tr>",paste0("<td>",hm_html(r),"</td>",collapse=""),"</tr>")),collapse="\n"),"</tbody></table>")
}
# Embedded, package-free figure example renderer; shared verbatim by both scripts.
# Examples are selected for audit illustration, not a confirmatory test family.
phmpt_anomaly_examples <- function(data_dir, output_dir=data_dir, high_cutoff=10,
  example_sites=c(coverage='US372',pcr='US325',source='US352',timing='US372',waves='US302')) {
  if(!is.finite(high_cutoff)||high_cutoff<=0||high_cutoff>=100)stop('Invalid positivity threshold.')
  required<-c('coverage','pcr','source','timing','waves')
  if(!all(required%in%names(example_sites)))stop('example_sites must name coverage, pcr, source, timing, waves.')
  rd<-function(n)read.csv(file.path(data_dir,n),stringsAsFactors=FALSE,check.names=FALSE,fileEncoding='UTF-8-BOM',na.strings=c('','NA'))
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  figdir<-file.path(output_dir,'anomaly_examples');dir.create(figdir,showWarnings=FALSE)
  csv<-function(z,n)write.csv(z,file.path(figdir,n),row.names=FALSE,na='')
  e<-rd('episodes.csv');bg<-rd('site_background_used.csv');map<-rd('site_region_coverage.csv');roster<-rd('subject_followup.csv')
  stopifnot(!anyDuplicated(e$episode_id),!anyDuplicated(bg[c('SITEID','week')]),!anyDuplicated(map$SITEID),!anyDuplicated(roster$subject))
  x<-merge(e[e$eligible,],bg[c('SITEID','week','positivity')],by=c('SITEID','week'),all.x=TRUE,sort=FALSE)
  stopifnot(nrow(x)==sum(e$eligible),!anyNA(x$tested),all(na.omit(x$positivity)>=0&na.omit(x$positivity)<=100))
  missing_background<-sum(is.na(x$positivity));x<-x[!is.na(x$positivity),]
  x<-x[order(x$subject,x$start,x$episode_id),];x$vaccine<-as.integer(x$arm=='mRNA-1273');x$high<-as.integer(x$positivity>=high_cutoff)
  sites<-sort(unique(roster$SITEID));M<-length(sites)
  pars<-if(file.exists(file.path(data_dir,'run_parameters.csv')))rd('run_parameters.csv') else data.frame(parameter=character(),value=character())
  parameter<-function(k,default){v<-pars$value[match(k,pars$parameter)];if(length(v)&&!is.na(v)&&nzchar(v))as.character(v) else default}
  horizon<-as.integer(parameter('onset_to_pcr_days','7'));visit_days<-as.integer(parameter('visit_to_pcr_days','3'))
  window<-parameter('window','selected window');visit_mode<-parameter('visit_evidence','RAW')
  count<-function(d,a){z<-d[d$arm==a,];c(y=sum(z$tested),n=nrow(z),pct=if(nrow(z))100*mean(z$tested) else NA_real_)}
  pair<-function(d){v<-count(d,'mRNA-1273');p<-count(d,'Placebo');c(v_y=v[1],v_n=v[2],v_pct=v[3],p_y=p[1],p_n=p[2],p_pct=p[3])}
  pair<-local({fun<-pair;function(d)setNames(fun(d),c('v_y','v_n','v_pct','p_y','p_n','p_pct'))})
  first<-function(d)d[!duplicated(d$subject),,drop=FALSE]
  pfmt<-function(p){if(!is.finite(p))'not estimable' else if(p<1e-6)'<0.000001' else formatC(p,format='f',digits=5)}
  pct<-function(p)ifelse(is.finite(p),sprintf('%.1f%%',p),'NA')
  lab<-function(y,n)if(n>0)sprintf('%d/%d (%.1f%%)',y,n,100*y/n) else '0/0 (no episodes)'
  region<-function(s){v<-map$region_name[match(s,map$SITEID)];ifelse(is.na(v),s,paste(s,v,sep=' | '))}
  high<-x[x$high==1,];fh<-first(high)
  screen<-do.call(rbind,lapply(sites,function(s){z<-pair(fh[fh$SITEID==s,]);p<-NA_real_
    tab<-rbind(c(z['v_y'],z['v_n']-z['v_y']),c(z['p_y'],z['p_n']-z['p_y']))
    if(all(rowSums(tab)>0)&&all(colSums(tab)>0))p<-fisher.test(tab)$p.value
    data.frame(SITEID=s,as.list(z),fisher_p=p,stringsAsFactors=FALSE,row.names=NULL)}))
  screen$fisher_q_BH<-p.adjust(ifelse(is.finite(screen$fisher_p),screen$fisher_p,1),'BH');csv(screen,'first_high_episode_site_screen.csv')
  # Site-specific arm x high-positivity interaction. Saturated logistic contrast;
  # participant-clustered HC0 covariance, G/(G-1), t(G-1) reference.
  interaction_test<-function(d){
    d$cell<-paste0(d$vaccine,d$high);keys<-c('11','10','01','00');sign<-c(1,-1,-1,1)
    n<-vapply(keys,function(k)sum(d$cell==k),integer(1));y<-vapply(keys,function(k)sum(d$tested[d$cell==k]),integer(1));p<-y/n
    G<-length(unique(d$subject));beta<-se<-pv<-NA_real_;status<-'EMPTY_OR_BOUNDARY_CELL'
    if(G>1&&all(n>0)&&all(y>0&y<n)){
      beta<-sum(sign*qlogis(p));j<-match(d$cell,keys)
      score<-rowsum(sign[j]*(d$tested-p[j])/(n[j]*p[j]*(1-p[j])),d$subject,reorder=FALSE)
      se<-sqrt(sum(score^2)*G/(G-1));pv<-2*pt(-abs(beta/se),G-1);status<-'ESTIMABLE'
    }
    crit<-if(G>1)qt(.975,G-1) else NA_real_
    data.frame(v_high_y=y[1],v_high_n=n[1],v_low_y=y[2],v_low_n=n[2],p_high_y=y[3],p_high_n=n[3],p_low_y=y[4],p_low_n=n[4],
      vaccine_change_pp=100*(p[1]-p[2]),placebo_change_pp=100*(p[3]-p[4]),difference_in_changes_pp=100*sum(sign*p),
      interaction_OR=exp(beta),ci_lower=exp(beta-crit*se),ci_upper=exp(beta+crit*se),interaction_p=pv,
      subjects=G,status=status,row.names=NULL)
  }
  subsets<-list(ALL_EPISODES=x,OUTSIDE_POSTDOSE=x[!x$early_postdose,],FIRST_EPISODE_OVERALL=first(x))
  interactions<-do.call(rbind,lapply(names(subsets),function(nm){d<-subsets[[nm]]
    z<-do.call(rbind,lapply(sites,function(s)cbind(subset=nm,SITEID=s,interaction_test(d[d$SITEID==s,]))))
    z$interaction_q_BH<-p.adjust(ifelse(is.finite(z$interaction_p),z$interaction_p,1),'BH');z}))
  csv(interactions,'site_high_low_interactions.csv')
  # All selected episodes and exact figure counts remain reviewable.
  selected<-x[x$SITEID%in%example_sites,];csv(selected,'selected_site_episodes.csv')
  source_site<-example_sites['source'];sz<-high[high$SITEID==source_site,]
  src_levels<-c('AE','COVID_FORM','AE; COVID_FORM');src_names<-c('AE only','COVID form only','Both at onset')
  source_counts<-do.call(rbind,lapply(c('mRNA-1273','Placebo'),function(a)do.call(rbind,lapply(seq_along(src_levels),function(i){
    z<-sz[sz$arm==a&sz$onset_source==src_levels[i],];N<-sum(sz$arm==a)
    data.frame(arm=a,source=src_levels[i],episodes=nrow(z),tested=sum(z$tested),share_pct=if(N)100*nrow(z)/N else NA_real_)}))))
  csv(source_counts,'recording_source_counts.csv')
  tz<-high[high$SITEID==example_sites['timing'],]
  timing_counts<-do.call(rbind,lapply(c('All episodes','Days 0-6','Outside days 0-6'),function(k){d<-switch(k,'All episodes'=tz,'Days 0-6'=tz[tz$early_postdose,],'Outside days 0-6'=tz[!tz$early_postdose,]);data.frame(group=k,as.list(pair(d)))}))
  csv(timing_counts,'postdose_timing_counts.csv')
  case<-high[high$SITEID==example_sites['pcr']&high$vaccine==1&high$onset_source=='COVID_FORM'&!high$tested,]
  case<-case[order(case$start,case$subject),c('episode_id','subject','start','last_symptom','early_postdose','onset_groups','positivity')]
  case$first_later_diagnostic_pcr<-NA_character_;case$days_to_later_pcr<-NA_integer_;case$later_result<-NA_character_
  case$within_same_episode_span<-FALSE;case$not_done_in_horizon<-FALSE;case$first_not_done_day<-NA_integer_
  pcr_available<-file.exists(file.path(data_dir,'pcr_records_reconciled.csv'))
  if(pcr_available){pcr<-rd('pcr_records_reconciled.csv')
    for(i in seq_len(nrow(case))){z<-pcr[pcr$subject==case$subject[i]&!is.na(pcr$test_date),];days<-as.integer(as.Date(z$test_date)-as.Date(case$start[i]))
      case$not_done_in_horizon[i]<-any(z$not_done&days>=0&days<=horizon,na.rm=TRUE)
      if(case$not_done_in_horizon[i])case$first_not_done_day[i]<-min(days[which(z$not_done&days>=0&days<=horizon)])
      zz<-z[which(z$usable_diagnostic&days>=0),];zz<-zz[order(zz$test_date),]
      if(nrow(zz)){case$first_later_diagnostic_pcr[i]<-zz$test_date[1];case$days_to_later_pcr[i]<-as.integer(as.Date(zz$test_date[1])-as.Date(case$start[i]));case$later_result[i]<-zz$result[1]
        case$within_same_episode_span[i]<-zz$test_date[1]<=case$last_symptom[i]}
    }
  }
  csv(case,'pcr_episode_review.csv')
  # Page coordinates are normalized; output is vector PDF, landscape A4.
  ink<-'#203148';muted<-'#536579';vcol<-'#21679C';pcol<-'#A86116';linecol<-'#D7DFE7'
  txt<-function(x,y,label,size=11,col=ink,just='left',bold=FALSE){grid::grid.text(label,x=x,y=y,just=just,gp=grid::gpar(fontfamily='Helvetica',fontsize=size,col=col,fontface=if(bold)'bold' else 'plain'))}
  rect<-function(x,y,w,h,fill,border=NA){grid::grid.rect(x=x,y=y,width=w,height=h,just=c('left','bottom'),gp=grid::gpar(fill=fill,col=border))}
  ln<-function(x,y,col=linecol,lwd=1,lty=1){grid::grid.lines(x=x,y=y,gp=grid::gpar(col=col,lwd=lwd,lty=lty))}
  point<-function(x,y,col,pch=16,size=3){grid::grid.points(x=x,y=y,default.units='npc',pch=pch,size=grid::unit(size,'mm'),gp=grid::gpar(col=col,fill=if(pch==1)NA else col))}
  wrap<-function(s,width=140)paste(strwrap(s,width=width),collapse='\n')
  footer<-function(number,notes){ln(c(.065,.94),c(.155,.155));txt(.065,.12,wrap(notes,142),9.3,muted,just=c('left','top'))
    txt(.065,.035,paste('P301 |',window,'|',visit_mode,'visit evidence | Exploratory audit examples'),8.5,muted);txt(.94,.035,paste0(number,' / 5'),8.5,muted,just='right')}
  header<-function(number,title,subtitle,unit){grid::grid.newpage();rect(0,.985,1,.015,vcol)
    txt(.065,.945,paste0('MODERNA P301  /  AUDIT EXAMPLE ',number),10,muted,bold=TRUE)
    txt(.065,.895,title,21,bold=TRUE);txt(.065,.845,subtitle,12,muted);txt(.065,.800,unit,10.5,muted)}
  axes<-function(y0=.28,h=.43,xleft=.18,xright=.91){for(a in seq(0,100,20)){y<-y0+h*a/100;ln(c(xleft,xright),c(y,y));txt(xleft-.02,y,paste0(a,'%'),10,muted,just='right')}}
  bar<-function(x,value,label,bottom=.28,height=.43,width=.14,col=vcol,cex=12){if(is.finite(value)){
    rect(x-width/2,bottom,width,height*value/100,col);txt(x,bottom+height*value/100+.025,label,cex,col,just='centre',bold=TRUE)}}
  empty<-function(number,title,s){header(number,title,region(s),'Requested example is unavailable in the selected data.');txt(.065,.60,'No eligible observations for this example.',14);footer(number,'No values were inferred or substituted.')}
  coverage_page<-function(){s<-example_sites['coverage'];z<-screen[screen$SITEID==s,]
    if(!nrow(z)||!z$v_n||!z$p_n)return(empty(1,'Unequal testing coverage during high positivity',s))
    header(1,'Unequal testing coverage during high positivity',region(s),paste0('County positivity >=',high_cutoff,'%; first eligible high-period episode per participant.'))
    axes();bar(.40,z$v_pct,lab(z$v_y,z$v_n),col=vcol);bar(.72,z$p_pct,lab(z$p_y,z$p_n),col=pcol)
    txt(.40,.235,'mRNA-1273',12,vcol,just='centre',bold=TRUE);txt(.72,.235,'Placebo',12,pcol,just='centre',bold=TRUE)
    txt(.065,.19,sprintf('Difference: %+.1f percentage points  |  Fisher p=%s  |  BH q=%s (%d sites)',z$v_pct-z$p_pct,pfmt(z$fisher_p),formatC(z$fisher_q_BH,format='f',digits=3),M),11,bold=TRUE)
    footer(1,paste0('Tested = recorded illness association plus diagnostic PCR by onset +',horizon,' days (PCR within ',visit_days,' days of visit). ',
      'This is a nominal site signal, not proof of differential handling. First-episode selection avoids repeated participants in Fisher tests.'))
  }
  pcr_page<-function(){s<-example_sites['pcr'];cc<-head(case,6)
    if(!nrow(cc)||!pcr_available)return(empty(2,'Missing timely endpoint: PCR record review',s))
    header(2,'Missing timely endpoint: PCR record review',region(s),paste0('Vaccine-arm COVID-form-onset episodes during positivity >=',high_cutoff,'%; selected episodes did not meet the visit + PCR endpoint.'))
    maxday<-max(c(horizon+3,cc$days_to_later_pcr),na.rm=TRUE)+3;xx<-function(d).30+.61*d/maxday
    rect(xx(0),.245,xx(horizon)-xx(0),.49,'#E8F1F5');ln(rep(xx(horizon),2),c(.245,.735),vcol,1,2)
    txt((xx(0)+xx(horizon))/2,.755,paste0('Days 0-',horizon),10,vcol,just='centre',bold=TRUE)
    yy<-seq(.68,.32,length.out=nrow(cc));if(nrow(cc)==1)yy<-.5
    for(i in seq_len(nrow(cc))){c<-cc[i,];y<-yy[i];txt(.065,y+.012,c$subject,11,bold=TRUE);txt(.065,y-.018,paste(c$episode_id,'|',c$start),8.6,muted)
      ln(c(xx(0),xx(maxday)),c(y,y),linecol);span<-max(0,as.integer(as.Date(c$last_symptom)-as.Date(c$start)))
      ln(c(xx(0),xx(min(span,maxday))),c(y,y),'#8398A8',5);point(xx(0),y,vcol,18,2.5)
      if(is.finite(c$days_to_later_pcr)){d<-c$days_to_later_pcr;co<-if(c$within_same_episode_span)pcol else muted
        point(xx(d),y,co,if(c$within_same_episode_span)16 else 1,3)
        label<-paste0('+',d,' days: ',tolower(c$later_result));txt(xx(d),y+.032,label,9.7,co,just='centre',bold=TRUE)
        txt(xx(d),y-.032,if(c$within_same_episode_span)'within episode span' else 'outside this episode span',8.7,co,just='centre')
      }else txt(xx(horizon+1),y+.025,'No later diagnostic PCR identified',10,muted)
      if(c$not_done_in_horizon){point(xx(c$first_not_done_day),y,'#A33E35',4,4);txt(xx(horizon+1),y-.018,paste0('NOT DONE record within days 0-',horizon),9.7,'#A33E35')}
    }
    for(d in sort(unique(c(0,horizon,seq(0,floor(maxday/5)*5,5))))){ln(c(xx(d),xx(d)),c(.245,.235),muted);txt(xx(d),.21,as.character(d),9,muted,just='centre')}
    txt(.60,.178,'Days since first recorded symptom',10,muted,just='centre')
    footer(2,'Blue shading: endpoint window. Gray line: reconstructed symptom span. Diamond: first symptom record. Cross: NOT DONE. A filled PCR dot falls within the span; a hollow dot falls outside it. A later PCR need not belong to the same illness; missing timely testing does not establish infection.')
  }
  source_page<-function(){s<-example_sites['source'];if(!nrow(sz))return(empty(3,'Different recording sources at symptom onset',s))
    header(3,'Different recording sources at symptom onset',region(s),paste0('All eligible episodes during positivity >=',high_cutoff,'%; each episode classified by its first recorded symptom date.'))
    cols<-c('#58778D','#3C8D7B','#C18736');N<-vapply(c('mRNA-1273','Placebo'),function(a)sum(sz$arm==a),integer(1))
    for(j in 1:2){a<-c('mRNA-1273','Placebo')[j];z<-source_counts[source_counts$arm==a,];y<-c(.65,.50)[j];txt(.065,y+.035,paste0(a,'\nN=',N[j]),11,if(j==1)vcol else pcol,bold=TRUE)
      left<-.25
      for(k in 1:3){fraction<-if(N[j])z$episodes[k]/N[j] else 0;w<-.65*fraction;rect(left,y,w,.07,cols[k])
        if(fraction>=.15)txt(left+w/2,y+.035,paste0(src_names[k],'  ',pct(100*fraction)),11,'white',just='centre',bold=TRUE)
        left<-left+w}
      small<-which(z$share_pct<15 & z$episodes>0);extra<-if(length(small))paste(paste0(src_names[small],': ',z$episodes[small],'/',N[j],' (',pct(z$share_pct[small]),')'),collapse='; ') else ''
      if(any(z$episodes==0))extra<-paste(extra,paste(paste0(src_names[z$episodes==0],': 0'),collapse='; '),sep='   ')
      txt(.25,y-.028,trimws(extra),9.5,muted)
    }
    txt(.25,.395,'Onset source',11,bold=TRUE);txt(.62,.395,'mRNA-1273 tested',11,vcol,just='centre',bold=TRUE);txt(.82,.395,'Placebo tested',11,pcol,just='centre',bold=TRUE)
    for(k in 1:3){y<-.35-(k-1)*.055;txt(.25,y,src_names[k],10.5);for(j in 1:2){z<-source_counts[source_counts$arm==c('mRNA-1273','Placebo')[j]&source_counts$source==src_levels[k],];txt(c(.62,.82)[j],y,lab(z$tested,z$episodes),10.5,just='centre')}}
    footer(3,'Bars show episode composition, not testing percentages. The lower table shows tested/eligible episodes within each category. Recording source may itself reflect the decision being audited; adjusting for it does not establish that the arm difference is clinically justified.')
  }
  timing_page<-function(){s<-example_sites['timing'];if(!nrow(tz))return(empty(4,'Early post-dose symptoms contribute to the gap',s))
    header(4,'Early post-dose symptoms contribute to the gap',region(s),paste0('Positivity >=',high_cutoff,'%; all eligible episodes retained, then partitioned by onset timing.'))
    axes(y0=.30,h=.39,xleft=.15,xright=.94);centres<-c(.28,.54,.80)
    for(i in 1:3){z<-timing_counts[i,];bar(centres[i]-.052,z$v_pct,sub(' (','\n(',lab(z$v_y,z$v_n),fixed=TRUE),bottom=.30,height=.39,width=.085,col=vcol,cex=10)
      bar(centres[i]+.052,z$p_pct,sub(' (','\n(',lab(z$p_y,z$p_n),fixed=TRUE),bottom=.30,height=.39,width=.085,col=pcol,cex=10)
      txt(centres[i],.257,as.character(z$group),10.5,just='centre',bold=TRUE)}
    txt(.36,.205,'mRNA-1273',10.5,vcol,just='centre',bold=TRUE);txt(.65,.205,'Placebo',10.5,pcol,just='centre',bold=TRUE)
    footer(4,'Days 0-6 means onset during the first seven calendar days after either original dose. The two timing subsets partition all episodes. A change after excluding early episodes localizes the signal; it does not explain its cause or rule out differential handling. Unlike page 1, participants may contribute repeated episodes.')
  }
  waves_page<-function(){s<-example_sites['waves'];z<-interactions[interactions$SITEID==s&interactions$subset=='ALL_EPISODES',]
    if(!nrow(z)||any(z[c('v_high_n','v_low_n','p_high_n','p_low_n')]==0))return(empty(5,'Opposing changes across positivity periods',s))
    header(5,'Opposing changes across positivity periods',region(s),paste0('All eligible episodes; early post-dose symptoms included. Low <',high_cutoff,'%; high >=',high_cutoff,'%.'))
    axes(y0=.30,h=.42,xleft=.23,xright=.80);xx<-c(.33,.71)
    vals<-list(c(100*z$v_low_y/z$v_low_n,100*z$v_high_y/z$v_high_n),c(100*z$p_low_y/z$p_low_n,100*z$p_high_y/z$p_high_n))
    labels<-list(c(lab(z$v_low_y,z$v_low_n),lab(z$v_high_y,z$v_high_n)),c(lab(z$p_low_y,z$p_low_n),lab(z$p_high_y,z$p_high_n)))
    for(j in 1:2){col<-c(vcol,pcol)[j];yy<-.30+.42*vals[[j]]/100;ln(xx,yy,col,2.6,if(j==1)1 else 2);point(xx,yy,col,if(j==1)16 else 15,3.3)
      txt(xx[1]-.018,yy[1]+.025,labels[[j]][1],11,col,just='right',bold=TRUE);txt(xx[2]+.018,yy[2]+.025,labels[[j]][2],11,col,just='left',bold=TRUE)}
    txt(.33,.25,'Low positivity',11,just='centre',bold=TRUE);txt(.71,.25,'High positivity',11,just='centre',bold=TRUE)
    txt(.065,.205,paste0('mRNA-1273: ',sprintf('%+.1f pp',z$vaccine_change_pp),'   |   Placebo: ',sprintf('%+.1f pp',z$placebo_change_pp),'   |   Difference between changes: ',sprintf('%+.1f pp',z$difference_in_changes_pp)),11,bold=TRUE)
    footer(5,paste0('Blue solid/circles = mRNA-1273; ochre dashed/squares = placebo. Site-specific logistic interaction OR=',sprintf('%.3f',z$interaction_OR),
      ', 95% CI ',sprintf('%.3f-%.3f',z$ci_lower,z$ci_upper),'; participant-clustered p=',pfmt(z$interaction_p),'. These exploratory patterns do not establish investigator intent.'))
  }
  plotters<-list(coverage_page,pcr_page,source_page,timing_page,waves_page)
  filenames<-c('01_coverage_US372.pdf','02_pcr_timing_US325.pdf','03_recording_source_US352.pdf','04_postdose_timing_US372.pdf','05_high_low_US302.pdf')
  filenames<-paste0(sprintf('%02d',1:5),'_',c('coverage','pcr_timing','recording_source','postdose_timing','high_low'),'_',unname(example_sites[required]),'.pdf')
  make_pdf<-function(path,funs){
    if(capabilities('cairo'))grDevices::cairo_pdf(path,width=11.69,height=8.27,family='sans',onefile=TRUE)
    else grDevices::pdf(path,width=11.69,height=8.27,family='Helvetica',onefile=TRUE,useDingbats=FALSE,title='Moderna P301 - site anomaly examples')
    tryCatch(for(f in funs)f(),finally=grDevices::dev.off())}
  combined<-file.path(output_dir,'site_anomaly_examples.pdf');make_pdf(combined,plotters)
  dpi<-as.integer(getOption('phmpt.report_png_dpi',300L))
  if(length(dpi)!=1L || is.na(dpi) || dpi<72L)stop('Invalid report PNG resolution.')
  for(i in seq_along(plotters)) {
    make_pdf(file.path(figdir,filenames[i]),plotters[i])
    grDevices::png(file.path(figdir,sub('[.]pdf$','.png',filenames[i])),width=11.69,height=8.27,
      units='in',res=dpi,type=if(capabilities('cairo'))'cairo' else getOption('bitmapType','Xlib'))
    tryCatch(plotters[[i]](),finally=grDevices::dev.off())
    if(capabilities('cairo')) {
      grDevices::svg(file.path(figdir,sub('[.]pdf$','.svg',filenames[i])),width=11.69,height=8.27,family='sans')
      tryCatch(plotters[[i]](),finally=grDevices::dev.off())
    }
  }
  selection<-data.frame(page=1:5,anomaly=required,SITEID=unname(example_sites[required]),file=filenames,
    rationale=c('Strongest nominal high-period Fisher signal with substantive denominators; not largest sparse percentage gap.',
      'Vaccine COVID-form episodes without a qualifying test illustrate record-level PCR timing.',
      'Strong contrast in onset recording source; source may itself be part of the mechanism.',
      'Early post-dose episodes contribute strongly; retain and partition rather than dismiss.',
      'Largest opposing absolute changes among supplied sites with at least ten episodes in all four cells.'),stringsAsFactors=FALSE)
  csv(selection,'example_selection.csv')
  csv(data.frame(parameter=c('module_version','high_positivity_threshold','onset_to_pcr_days','visit_to_pcr_days','window','visit_evidence','sites_in_multiplicity_family','episodes_missing_background','selection_type'),
    value=c('1.2.0',high_cutoff,horizon,visit_days,window,visit_mode,M,missing_background,'Fixed illustrative sites from the reviewed export; not automatically reselected')),'figure_parameters.csv')
  csv(data.frame(file=c('episodes.csv','site_background_used.csv','site_region_coverage.csv','subject_followup.csv'),
    md5=unname(tools::md5sum(file.path(data_dir,c('episodes.csv','site_background_used.csv','site_region_coverage.csv','subject_followup.csv'))))),'figure_input_checksums.csv')
  writeLines(c('Five illustrative audit patterns; one vector plot per PDF page, plus separate single-page PDFs and PNGs; SVGs are added when Cairo is available.',
    'The examples are post hoc selections from the reviewed export, not a representative or confirmatory sample.',
    'Page 1 uses the first eligible high-period episode per participant. Pages 3-5 use all episodes; do not interchange the denominators.',
    'The main results retain post-dose days 0-6; the timing panel partitions them, without attributing their cause.',
    'Standard site interaction inference is unavailable for empty or all-tested/all-untested cells. Such cells are labelled, not treated as evidence of no interaction.',
    'BH corrections retain all roster sites; untestable p-values contribute 1. Raw all-episode Fisher p-values are not used.',
    'All figures depict documented testing, not confirmed infections, missed infection counts, or vaccine efficacy.',
    'Method reference: https://sandwich.r-forge.r-project.org/reference/vcovCL.html'),file.path(figdir,'README.txt'))
  message('Anomaly examples: ',normalizePath(combined,winslash='/'))
  invisible(list(pdf=combined,selection=selection,screen=screen,interactions=interactions,cases=case,source=source_counts,timing=timing_counts))
}

# Shared report-insertion writer, embedded in both stand-alone scripts.
# This formats existing calculations; no symptom, visit or PCR linkage is changed.
phmpt_report_insert <- function(data_dir, output_dir, start_paragraph=25L,
  reference_start=32L, audit_results_dir=output_dir,
  repo_url='https://github.com/OpenVaet/moderna_trial_analysis/blob/main', audit_root=PHMPT_ROOT) {
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  if(length(start_paragraph)!=1L || is.na(start_paragraph) || start_paragraph<1)stop('Invalid starting paragraph')
  if(length(reference_start)!=1L || is.na(reference_start) || reference_start<1)stop('Invalid starting footnote')
  rd<-function(p)read.csv(p,stringsAsFactors=FALSE,check.names=FALSE,fileEncoding='UTF-8-BOM',na.strings=c('','NA'))
  esc<-function(s){s<-gsub('&','&amp;',as.character(s),fixed=TRUE);s<-gsub('<','&lt;',s,fixed=TRUE);s<-gsub('>','&gt;',s,fixed=TRUE);gsub('"','&quot;',s,fixed=TRUE)}
  nfmt<-function(n)format(n,big.mark=',',trim=TRUE,scientific=FALSE)
  num<-function(n,d=1)ifelse(is.finite(n),formatC(n,format='f',digits=d,decimal.mark='.'),'not estimable')
  rate<-function(k,n)if(n>0)paste0(nfmt(k),'/',nfmt(n),' (',num(100*k/n),'%)') else '0/0 (no eligible episodes)'
  pf<-function(p){if(!is.finite(p))return('not estimable');for(v in c(1e-6,1e-5,1e-4,1e-3,1e-2))if(p<v)return(paste('<',formatC(v,format='f',digits=ceiling(-log10(v)))));formatC(p,format='f',digits=4)}
  ptext<-function(p){s<-pf(p);paste0('p ',if(startsWith(s,'<'))'' else '= ',s)}
  figdir<-file.path(output_dir,'anomaly_examples')
  files<-c('first_high_episode_site_screen.csv','site_high_low_interactions.csv','pcr_episode_review.csv',
    'recording_source_counts.csv','postdose_timing_counts.csv','example_selection.csv','figure_parameters.csv')
  if(!all(file.exists(file.path(figdir,files))))stop('Generate phmpt_anomaly_examples() before the report insertion.')
  e<-rd(file.path(data_dir,'episodes.csv'));b<-rd(file.path(data_dir,'site_background_used.csv'))
  map<-rd(file.path(data_dir,'site_region_coverage.csv'));pars<-rd(file.path(data_dir,'run_parameters.csv'))
  fp<-rd(file.path(figdir,'figure_parameters.csv'));sel<-rd(file.path(figdir,'example_selection.csv'))
  getp<-function(p,k,default=''){v<-p$value[match(k,p$parameter)];if(!length(v)||is.na(v)||!nzchar(v))default else as.character(v)}
  horizon<-as.integer(getp(pars,'onset_to_pcr_days','7'));visit_days<-as.integer(getp(pars,'visit_to_pcr_days','3'))
  gap<-as.integer(getp(pars,'gap_days','7'));cut<-as.numeric(getp(fp,'high_positivity_threshold','10'))
  window<-getp(pars,'window','selected window');mode<-getp(pars,'visit_evidence','RAW')
  screen<-rd(file.path(figdir,files[1]));waves<-rd(file.path(figdir,files[2]));cases<-rd(file.path(figdir,files[3]))
  sources<-rd(file.path(figdir,files[4]));timing<-rd(file.path(figdir,files[5]))
  x<-merge(e[e$eligible,],b[c('SITEID','week','positivity')],by=c('SITEID','week'),all.x=TRUE,sort=FALSE)
  if(nrow(x)!=sum(e$eligible)||anyDuplicated(b[c('SITEID','week')]))stop('Invalid background join')
  missing_bg<-sum(!is.finite(x$positivity));M<-nrow(screen)
  region<-function(s){v<-map$region_name[match(s,map$SITEID)];ifelse(is.na(v),s,paste(s,v,sep=' - '))}
  pair<-function(d) {
    v<-d[d$arm=='mRNA-1273',];p<-d[d$arm=='Placebo',]
    c(vk=sum(v$tested),vn=nrow(v),pk=sum(p$tested),pn=nrow(p))
  }
  # Reuse adjusted models only if their input hashes match THIS source export.
  model_file<-file.path(audit_results_dir,'adjusted_model_key_coefficients.csv')
  check_file<-file.path(audit_results_dir,'input_checksums.csv');models<-NULL
  checks_ok<-FALSE
  if(all(file.exists(c(model_file,check_file)))) {
    ch<-rd(check_file);required_hashes<-c('episodes.csv','site_background_used.csv','site_region_coverage.csv','subject_followup.csv','site_arm_weekly_testing.csv','pcr_records_reconciled.csv','illness_visits.csv')
    if(all(required_hashes%in%ch$file) && all(file.exists(file.path(data_dir,required_hashes)))) {
      checks_ok<-all(tolower(ch$md5[match(required_hashes,ch$file)])==tolower(unname(tools::md5sum(file.path(data_dir,required_hashes)))))
    }
    if(checks_ok)models<-rd(model_file)
  }
  modelrow<-function(model,term='vaccine') {
    if(is.null(models))return(NULL)
    q<-models[models$model==model & models$term==term & models$cluster_unit=='SITEID',,drop=FALSE]
    if(nrow(q)!=1L || !is.finite(q$OR)||!isTRUE(q$converged))NULL else q
  }
  body<-character();plain<-character();number<-as.integer(start_paragraph);tables<-list();audit_tables<-list();captions<-list();fig_index<-0L
  addh<-function(s){body<<-c(body,paste0('<h2>',esc(s),'</h2>'));plain<<-c(plain,s,'')}
  addp<-function(s,refs=integer()) {
    marker<-if(length(refs))paste0(' <sup>',paste0('<a href="#ref-',refs,'">',refs,'</a>',collapse=', '),'</sup>') else ''
    body<<-c(body,paste0('<p class="numbered"><span class="number">',number,'.</span> ',esc(s),marker,'</p>'))
    plain<<-c(plain,paste0(number,'. ',s,if(length(refs))paste0(' [',paste(refs,collapse=', '),']') else ''),'')
    number<<-number+1L
  }
  addtable<-function(title,headers,rows,note='') {
    rows<-as.data.frame(rows,stringsAsFactors=FALSE);names(rows)<-headers
    audit_tables[[title]]<<-rows
    ht<-paste0('<table><caption>',esc(title),'</caption><thead><tr>',paste0('<th>',esc(headers),'</th>',collapse=''),'</tr></thead><tbody>',
      paste(apply(rows,1,function(z)paste0('<tr>',paste0('<td>',esc(z),'</td>',collapse=''),'</tr>')),collapse=''),'</tbody></table>')
    body<<-c(body,ht,if(nzchar(note))paste0('<p class="table-note">',esc(note),'</p>') else '')
    plain<<-c(plain,title,paste(headers,collapse='\t'),apply(rows,1,paste,collapse='\t'),note,'')
    name<-sprintf('report_table_%02d.csv',length(tables)+1L);write.csv(rows,file.path(output_dir,name),row.names=FALSE,na='');tables[[name]]<<-rows
  }
  b64<-function(path) {
    raw<-readBin(path,'raw',n=file.info(path)$size);n<-length(raw);pad<-(3L-n%%3L)%%3L
    v<-matrix(c(as.integer(raw),rep(0L,pad)),nrow=3L)
    ix<-rbind(v[1,]%/%4L,(v[1,]%%4L)*16L+v[2,]%/%16L,(v[2,]%%16L)*4L+v[3,]%/%64L,v[3,]%%64L)
    alphabet<-strsplit('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/','',fixed=TRUE)[[1]]
    out<-alphabet[as.vector(ix)+1L];if(pad)out[(length(out)-pad+1L):length(out)]<-'=';paste0(out,collapse='')
  }
  figure<-function(anomaly,caption) {
    s<-sel[sel$anomaly==anomaly,,drop=FALSE];if(nrow(s)!=1L)return(invisible(NULL))
    fig_index<<-fig_index+1L;png<-sub('[.]pdf$','.png',s$file);path<-file.path(figdir,png)
    label<-paste0('Figure ',fig_index,'. ',caption)
    if(file.exists(path))body<<-c(body,paste0('<figure><img src="data:image/png;base64,',b64(path),'" alt="',esc(label),'"/><figcaption>',esc(label),'</figcaption></figure>'))
    else body<<-c(body,paste0('<p class="figure-caption">',esc(label),' (Use ',esc(s$file),'.)</p>'))
    plain<<-c(plain,label,paste0('[Insert image: anomaly_examples/',png,']'),'')
    captions[[length(captions)+1L]]<<-data.frame(figure=fig_index,anomaly=anomaly,SITEID=s$SITEID,png=png,pdf=s$file,caption=label)
  }
  addh('Symptomatic episodes and COVID testing')
  addp(sprintf('The next analysis examines whether COVID-compatible symptoms documented in either reporting channel were followed by an illness-visit association and a diagnostic PCR result. Symptoms recorded on COVID forms and the matched AE terms described above were combined into reconstructed episodes: overlapping intervals, or intervals separated by no more than %d unrecorded days, were grouped together. Each episode contributes once, in the week of its first recorded symptom. For COVID forms, the assessment date may differ from actual symptom onset; an AE with an unknown end date contributes its onset only. This differs from the participant-symptom-group-week observations used in the preceding charts. All eligible episodes, including those beginning during post-dose days 0-6, remain in the main analysis.',gap),c(reference_start,reference_start+1L))
  addp(sprintf('For this audit, a tested episode requires a documented illness-visit association followed by a diagnostic PCR result on the same date or within %d calendar days, and no later than %d days after the first recorded symptom. Complete follow-up through that horizon is required for both arms. These are reconstruction and linkage rules, rather than a direct test of compliance with every protocol trigger, duration requirement or 72-hour deadline. In the %s visit-evidence setting, %s. The analysis uses the %s observation window, from %s to %s, subject to participant-specific follow-up boundaries.',
    visit_days,horizon,mode,if(mode=='RAW')'an original illness-visit label on a laboratory record can establish the association; the measure does not independently prove physical attendance' else 'the association requires independent SV or FAOT visit documentation',
    tolower(window),getp(pars,'display_start'),getp(pars,'display_end')))
  if(window=='FULL')addp('The full observation window includes follow-up after unblinding and crossover. Original randomized-arm labels are retained, so comparisons over that window are descriptive and do not represent two consistently distinct treatments.')
  rates<-lapply(list('All eligible episodes'=rep(TRUE,nrow(x)),'Lower positivity'=is.finite(x$positivity)&x$positivity<cut,'Higher positivity'=is.finite(x$positivity)&x$positivity>=cut),function(k)pair(x[k,]))
  allp<-rates[[1]];hp<-rates[[3]]
  addp(paste0('Across all sites, the documented visit-plus-PCR endpoint was met for ',rate(allp['vk'],allp['vn']),' of vaccine-arm episodes and ',rate(allp['pk'],allp['pn']),' of placebo-arm episodes. Local background positivity was then divided into lower (<',num(cut,0),'%) and higher (>=',num(cut,0),'%) periods. This threshold is an exploratory analytical choice, not a protocol criterion. County or state background testing is a contextual proxy for the site location, not a measure of infection in trial participants. During higher-positivity periods, the corresponding proportions were ',rate(hp['vk'],hp['vn']),' and ',rate(hp['pk'],hp['pn']),'.',if(missing_bg)paste0(' Background positivity was unavailable for ',missing_bg,' episodes, which remain in the overall totals but not in the high/low comparisons.') else ''))
  rows<-do.call(rbind,lapply(seq_along(rates),function(i){z<-rates[[i]];c(names(rates)[i],rate(z['vk'],z['vn']),rate(z['pk'],z['pn']),paste0(num(100*(z['vk']/z['vn']-z['pk']/z['pn'])),' pp'))}))
  addtable('Documented testing after symptomatic episodes',c('Period','mRNA-1273','Placebo','Difference'),rows,'Numerators are episodes meeting the documented visit-plus-PCR endpoint. Repeated episodes from one participant can contribute to this descriptive table. Differences are mRNA-1273 minus placebo.')
  ma<-modelrow('HIGH10_SYMPTOMS_TIMING');ms<-modelrow('HIGH10_SOURCE_SENSITIVITY');mi<-modelrow('HIGH_LOW_INTERACTION','vaccine:high10')
  # The existing model names specifically denote the 10% threshold.
  if(cut==10 && !is.null(ma)) {
    ptxt<-paste0('During higher-positivity periods, the odds of documented testing remained lower in the vaccine arm after adjustment for site, calendar week, symptom groups at onset and post-dose timing (odds ratio ',num(ma$OR,3),', 95% CI ',num(ma$ci_lower,3),' to ',num(ma$ci_upper,3),'; ',ptext(ma$p),'). Uncertainty estimates account for clustering within sites, including repeated episodes from the same participant.')
    if(!is.null(ms))ptxt<-paste0(ptxt,' Including the initial recording source as an additional sensitivity adjustment gave an odds ratio of ',num(ms$OR,3),' (95% CI ',num(ms$ci_lower,3),' to ',num(ms$ci_upper,3),'; ',ptext(ms$p),'). Recording source may itself reflect the decision under review; this adjustment does not establish that the original difference was clinically justified. Severity and complete clinical presentation were not matched by these models.')
    addp(ptxt,reference_start+1L)
    specs<-c(HIGH10_SITE_WEEK='Site and calendar week',HIGH10_SYMPTOMS_TIMING='Plus onset symptoms and post-dose timing',HIGH10_SOURCE_SENSITIVITY='Plus onset recording source (sensitivity)',HIGH10_OUTSIDE_POSTDOSE_ADJUSTED='Outside days 0-6: site, week and symptoms')
    mr<-lapply(names(specs),function(k){q<-modelrow(k);if(is.null(q))return(NULL);c(specs[[k]],num(q$OR,3),paste(num(q$ci_lower,3),num(q$ci_upper,3),sep=' to '),pf(q$p))})
    addtable('Adjusted comparisons during higher background positivity',c('Comparison','Odds ratio','95% CI','p'),do.call(rbind,mr),'Odds ratios compare documented testing in vaccine versus placebo episodes. They are not vaccine-efficacy estimates. The post-dose restriction is a sensitivity analysis; early episodes remain in the main comparison.')
  }
  addh('Site-level imbalances and record review')
  nominal<-screen[is.finite(screen$fisher_p)&screen$fisher_p<.05,,drop=FALSE];nominal<-nominal[order(nominal$fisher_p),]
  corrected<-sum(screen$fisher_q_BH<.05,na.rm=TRUE)
  addp(paste0('For the site-level screen, each participant contributes only their first eligible episode during higher-positivity periods. This avoids treating repeated episodes from the same person as independent observations in the Fisher exact tests. ',nrow(nominal),' of ',M,' sites had a nominal two-sided p-value below 0.05. After Benjamini-Hochberg correction across all ',M,' sites, ',if(corrected==0)'none remained below 0.05' else paste(corrected,'remained below 0.05'),'. The individual signals therefore identify records for further review rather than establish site-specific differential treatment.'))
  if(nrow(nominal)) {
    rows<-cbind(region(nominal$SITEID),mapply(rate,nominal$v_y,nominal$v_n),mapply(rate,nominal$p_y,nominal$p_n),vapply(nominal$fisher_p,pf,character(1)),num(nominal$fisher_q_BH,3))
    addtable('Sites with nominal testing imbalances during higher positivity',c('Site / county','mRNA-1273','Placebo','Fisher p','BH q'),rows,paste0('First eligible higher-positivity episode per participant. Multiplicity correction includes all ',M,' sites. Small denominators, including US313, limit interpretation.'))
  }
  figure('coverage',paste0('Documented testing at ',region(sel$SITEID[sel$anomaly=='coverage']),', using each participant\'s first eligible episode during higher-positivity periods.'))
  cs<-sel$SITEID[sel$anomaly=='pcr'];within<-cases[!is.na(cases$within_same_episode_span)&cases$within_same_episode_span & is.finite(cases$days_to_later_pcr),]
  days<-sort(unique(within$days_to_later_pcr));nd<-sum(cases$not_done_in_horizon,na.rm=TRUE)
  cp<-paste0('At ',region(cs),', ',nrow(cases),' vaccine-arm episodes beginning on COVID symptom forms during higher-positivity periods did not meet the visit-plus-PCR endpoint within ',horizon,' days.')
  if(nrow(within))cp<-paste0(cp,' For ',nrow(within),' of these episodes, a later diagnostic PCR was documented within the reconstructed symptom interval, at ',paste(days,collapse=' and '),' days after the first recorded symptom; the recorded results were ',paste(tolower(unique(within$later_result)),collapse=' and '),'.')
  if(nd)cp<-paste0(cp,' A NOT DONE record was present within the endpoint window for ',nd,' episode',if(nd==1)'' else 's','.')
  addp(paste0(cp,' These observations distinguish delayed testing, testing recorded outside the reconstructed illness, and absence of a qualifying result. They do not establish that an untested episode was COVID-19, or that a later negative result excludes infection at the initial presentation.'))
  figure('pcr',paste0('Recorded PCR timing for the selected vaccine-arm COVID-form episodes at ',cs,'. The shaded interval is the audit linkage window, not a restatement of the protocol deadline.'))
  ss<-sel$SITEID[sel$anomaly=='source'];vn<-sum(sources$episodes[sources$arm=='mRNA-1273']);pn<-sum(sources$episodes[sources$arm=='Placebo'])
  va<-sum(sources$episodes[sources$arm=='mRNA-1273'&sources$source=='AE']);pa<-sum(sources$episodes[sources$arm=='Placebo'&sources$source=='AE'])
  addp(paste0('At ',region(ss),', the higher-positivity episodes differed substantially in their initial recording source: ',rate(va,vn),' in the vaccine arm began only in AE records, compared with ',rate(pa,pn),' in the placebo arm. The figure separates the composition of episodes from the proportion tested within each source category. Here, source means the record present at the first recorded symptom date; it differs from the preceding weekly both-source measure. This concentration of the testing difference warrants review of how complaints entered the COVID assessment pathway.'))
  figure('source',paste0('Initial recording source and subsequent documented testing at ',ss,' during higher-positivity periods. The bars describe source composition; the lower table describes testing within each category.'))
  ts<-sel$SITEID[sel$anomaly=='timing'];early<-timing[timing$group=='Days 0-6',];late<-timing[timing$group=='Outside days 0-6',]
  tp<-paste0('The role of post-dose timing was examined without removing early symptoms from the main analysis. At ',region(ts),', among all higher-positivity episodes, ',rate(early$v_y,early$v_n),' of vaccine-arm and ',rate(early$p_y,early$p_n),' of placebo-arm episodes beginning on days 0-6 met the endpoint. Outside this interval, the proportions were ',rate(late$v_y,late$v_n),' and ',rate(late$p_y,late$p_n),'. This partition localizes the imbalance; it does not determine whether clinical differences, discretionary testing decisions or documentation explain it. The protocol\'s special provisions for early post-dose symptoms therefore remain relevant to record review.')
  addp(tp);figure('timing',paste0('All eligible higher-positivity episodes at ',ts,', partitioned by onset during or outside post-dose days 0-6. Participants may contribute repeated episodes; these denominators differ from the first-episode Fisher screen.'))
  addh('Testing across lower and higher local positivity')
  w<-waves[waves$subset=='ALL_EPISODES',];chosen<-w[w$v_low_n>=10&w$v_high_n>=10&w$p_low_n>=10&w$p_high_n>=10&w$vaccine_change_pp<0&w$placebo_change_pp>0,]
  chosen<-chosen[order(chosen$difference_in_changes_pp),]
  ws<-sel$SITEID[sel$anomaly=='waves'];z<-w[w$SITEID==ws,]
  addp(paste0('Within-site comparisons identified ',nrow(chosen),' sites with at least ten episodes in every arm-period cell where vaccine-arm testing was lower during higher background positivity while placebo-arm testing was higher. At ',region(ws),', the vaccine proportion changed from ',rate(z$v_low_y,z$v_low_n),' to ',rate(z$v_high_y,z$v_high_n),', while placebo changed from ',rate(z$p_low_y,z$p_low_n),' to ',rate(z$p_high_y,z$p_high_n),'. The vaccine change minus the placebo change was ',num(z$difference_in_changes_pp),' percentage points. These are comparisons between sets of weeks, not necessarily one chronological transition. The site-specific arm-by-period interaction had participant-clustered ',ptext(z$interaction_p),'; the descriptive pattern does not by itself establish an interaction.'))
  if(nrow(chosen)) {
    rows<-cbind(chosen$SITEID,mapply(rate,chosen$v_low_y,chosen$v_low_n),mapply(rate,chosen$v_high_y,chosen$v_high_n),mapply(rate,chosen$p_low_y,chosen$p_low_n),mapply(rate,chosen$p_high_y,chosen$p_high_n),vapply(chosen$interaction_p,pf,character(1)))
    addtable('Sites with opposing changes across positivity periods',c('Site','mRNA low','mRNA high','Placebo low','Placebo high','Interaction p'),rows,'All eligible episodes, including post-dose days 0-6. Interaction p-values account for repeated participants but are nominal; the full site family and non-estimable cells are retained in the CSV audit.')
  }
  figure('waves',paste0('Testing at ',region(ws),' in lower and higher background positivity. Lines connect two pooled period categories, not successive calendar weeks.'))
  conclusion<-'Taken together, these findings identify differences in documented case ascertainment that warrant examination of the original symptom descriptions, their timing, clinical severity, visit records and reasons for performing or omitting a swab. The early post-dose findings and the recording-source differences remain part of that investigation. Neither the figures nor the site screens alone establish failed blinding, investigator intent, the number of missed infections or a resulting bias in vaccine efficacy.'
  if(cut==10 && !is.null(mi))conclusion<-paste0('Across sites, the adjusted arm-by-positivity interaction was ',num(mi$OR,3),' (95% CI ',num(mi$ci_lower,3),' to ',num(mi$ci_upper,3),'; ',ptext(mi$p),'). ',conclusion)
  addp(conclusion)
  notes<-data.frame(number=reference_start+0:1,title=c('Reconstruction of symptomatic episodes and weekly testing heatmaps','Site-level testing audit and reproducible comparisons'),
    url=paste0(sub('/$','',repo_url),'/',c('analyse_symptom_testing_heatmaps.R','analyse_site_testing_imbalance.R')))
  refs<-paste0('<h3>References to add</h3><div class="references">',paste0('<p id="ref-',notes$number,'">',notes$number,'. <a href="',esc(notes$url),'">',esc(notes$url),'</a></p>',collapse=''),'</div>')
  body<-c(body,refs);plain<-c(plain,'References to add',paste0(notes$number,'. ',notes$url))
  article<-paste(body,collapse='\n')
  css<-'body{margin:28px auto;padding:0 24px;max-width:980px;color:#171717;font:19px/1.6 Georgia,serif}h1,h2,h3{font-family:Arial,sans-serif;line-height:1.2}h1{font-size:28px}h2{font-size:27px;margin:36px 0 20px}h3{font-size:20px}.numbered{padding-left:2.6em;position:relative}.number{position:absolute;left:0}table{border-collapse:collapse;width:100%;font:14px/1.4 Arial,sans-serif;margin:24px 0 8px}caption{text-align:left;font-weight:bold;margin-bottom:10px}th,td{padding:9px 7px;border-bottom:1px solid #ccc;text-align:left;vertical-align:top}th{background:#eee}tbody tr:nth-child(even){background:#fafafa}.table-note,figcaption{font:14px/1.45 Arial,sans-serif;color:#444}.references{font:14px/1.5 Arial,sans-serif;overflow-wrap:anywhere}figure{margin:26px 0}img{max-width:100%;height:auto}.editor{font:15px/1.45 Arial,sans-serif;background:#f2f5f7;padding:18px;border:1px solid #d5dce1;margin-bottom:30px}.editor p{margin:8px 0}button{font:inherit;padding:8px 12px;margin:5px 8px 5px 0;cursor:pointer}a{color:#155fa0}.status{font-size:14px}sup{font-size:12px}@media print{.editor{display:none}body{margin:0;padding:0;font-size:12pt;max-width:none}table{font-size:9pt}h2{font-size:18pt}h2,h3{break-after:avoid}tr,figure{break-inside:avoid}a{color:inherit}}'
  status<-if(!is.null(ma)&&cut==10)'Adjusted comparisons are included and their input hashes match this export.' else 'This insert contains descriptive and site-level comparisons. Run analyse_site_testing_imbalance.R to add the adjusted models; no earlier-run estimates have been substituted.'
  guide<-paste0('<aside class="editor"><h1>Insert into the current draft</h1><p>Place after paragraph 24 and the US375 chart, before <b>Anti-N-Antibodies</b>. This section uses paragraphs <b>',start_paragraph,'-',number-1L,'</b>; renumber the following anti-N paragraph to <b>',number,'</b>. Add footnotes <b>',reference_start,' and ',reference_start+1L,'</b> after the existing reference 31. Figure numbers are local to this insertion.</p><p>',status,'</p><p>Copy the formatted text and tables, then insert the supplied PNG figures at their captions. You can also copy the whole section including embedded images; image-paste support depends on the destination editor. This HTML remains readable offline.</p><button onclick="copySection(false)">Copy text and tables</button><button onclick="copySection(true)">Copy including figures</button><button onclick="selectSection()">Select section</button><p id="copy-status" class="status" role="status"></p></aside>')
  js<-"<script>function selectSection(){const r=document.createRange();r.selectNodeContents(document.getElementById('report-insert'));const s=getSelection();s.removeAllRanges();s.addRange(r);}async function copySection(images){const node=document.getElementById('report-insert').cloneNode(true);if(!images)node.querySelectorAll('img').forEach(x=>x.remove());node.querySelectorAll('table').forEach(x=>x.style.cssText='border-collapse:collapse;width:100%;font-family:Arial;font-size:11pt');node.querySelectorAll('td,th').forEach(x=>x.style.cssText='border:1px solid #cccccc;padding:6px;text-align:left');node.querySelectorAll('.numbered').forEach(x=>x.style.cssText='margin:0 0 14px 0');node.querySelectorAll('.number').forEach(x=>x.style.cssText='position:static');const html=node.innerHTML;const text=document.getElementById('report-plain').value;try{if(!navigator.clipboard||!window.ClipboardItem)throw new Error('Clipboard unavailable');await navigator.clipboard.write([new ClipboardItem({'text/html':new Blob([html],{type:'text/html'}),'text/plain':new Blob([text],{type:'text/plain'})})]);document.getElementById('copy-status').textContent='Copied. Paste into your report editor.';}catch(e){selectSection();document.getElementById('copy-status').textContent='Section selected. Press Ctrl+C (or Cmd+C), then paste into your report editor.';}}</script>"
  writeLines(paste0('<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Moderna audit - section for the draft</title><style>',css,'</style></head><body>',guide,'<article id="report-insert">',article,'</article><textarea id="report-plain" hidden>',esc(paste(plain,collapse='\n')),'</textarea>',js,'</body></html>'),file.path(output_dir,'report_insert.html'),useBytes=TRUE)
  file.copy(file.path(output_dir,'report_insert.html'),file.path(output_dir,'site_testing_positivity_findings.html'),overwrite=TRUE)
  writeLines(plain,file.path(output_dir,'report_insert.txt'),useBytes=TRUE)
  write.csv(notes,file.path(output_dir,'report_references.csv'),row.names=FALSE)
  if(length(captions))write.csv(do.call(rbind,captions),file.path(output_dir,'report_figure_captions.csv'),row.names=FALSE)
  writeLines(c('Integration into Moderna mRNA-1273 Trial - Why more Transparency is Required - v0.1',
    paste0('Insert after the existing symptom-recording paragraph 24 and US375 illustration. New paragraphs ',start_paragraph,'-',number-1L,'. Renumber the following Anti-N paragraph to ',number,'.'),
    'The draft has two paragraphs numbered 13 and two numbered 24; this insertion does not rewrite unrelated draft passages.',
    'report_insert.html is self-contained. Use its copy controls; figures are also supplied as PNG (default 300 dpi), vector PDF and SVG when Cairo is available.',
    'report_insert.txt contains the same text and tab-separated tables for plain-text workflows.',
    'The heatmap script generates the descriptive insertion. The standalone site audit adds adjusted comparisons when their input checksums match.',
    'Values, tables, captions and statements about multiplicity are rebuilt from the selected export. Fixed example sites remain identified as illustrative.',
    'P-values retain full precision in the statistical CSVs; display thresholds run from <0.01 to <0.000001.',
    status),file.path(output_dir,'REPORT_INTEGRATION.txt'),useBytes=TRUE)
  phmpt_export_audit_assets(audit_root,output_dir,data_dir,audit_tables)
  message('Copy-ready report section: ',normalizePath(file.path(output_dir,'report_insert.html'),winslash='/'))
  invisible(list(html=file.path(output_dir,'report_insert.html'),paragraph_start=start_paragraph,paragraph_end=number-1L,adjusted_models=!is.null(ma)&&cut==10))
}


hm_region_cohort <- function(ctx,map) {
  # Header counts use the entire randomized roster, including subjects without
  # symptoms or dose 1; they are distinct from the episode/testing denominator.
  need<-c("subject","SITEID","arm")
  if(!all(need %in% names(ctx)))stop("Region cohort counts require subject, SITEID and arm in ADSL context.")
  if(anyDuplicated(ctx$subject))stop("Duplicate randomized subject in region cohort roster.")
  m<-unique(map[c("SITEID","region_key","region_name","region_level","region_id")])
  if(anyDuplicated(m$SITEID))stop("A site is assigned to more than one background region.")
  j<-match(ctx$SITEID,m$SITEID)
  if(anyNA(j))stop("Some randomized sites lack a background-region assignment.")
  regions<-unique(m[c("region_key","region_name","region_level","region_id")])
  regions<-regions[order(regions$region_name),,drop=FALSE]
  regions$site_count<-vapply(regions$region_key,function(k)sum(m$region_key==k),integer(1))
  regions$subjects<-vapply(regions$region_key,function(k)sum(m$region_key[j]==k),integer(1))
  regions$mrna1273_subjects<-vapply(regions$region_key,function(k)sum(m$region_key[j]==k & ctx$arm=="mRNA-1273"),integer(1))
  regions$placebo_subjects<-vapply(regions$region_key,function(k)sum(m$region_key[j]==k & ctx$arm=="Placebo"),integer(1))
  regions$trial_cohort_subjects<-nrow(ctx)
  regions$trial_cohort_pct<-100*regions$subjects/nrow(ctx)
  stopifnot(sum(regions$subjects)==nrow(ctx),all(regions$subjects==regions$mrna1273_subjects+regions$placebo_subjects))
  regions$display_name<-paste0(regions$region_name," (",ws_n(regions$subjects)," subjects; ",
    formatC(regions$trial_cohort_pct,format="f",digits=2,decimal.mark="."),"% of trial cohort)")
  rownames(regions)<-NULL;regions
}
hm_report <- function(out,cells,bg,summary,methods,pngs,manifest,quality,ctx) {
  regions<-hm_region_cohort(ctx,bg$map)
  ws_csv(regions,file.path(out,"region_cohort_counts.csv"))
  nav<-paste0('<option value="r_',regions$region_key,'">',hm_html(regions$display_name),'</option>')
  html<-c('<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Weekly symptom testing by site</title>',
    '<style>body{font:16px/1.5 Arial,sans-serif;color:#253448;max-width:1300px;margin:30px auto;padding:0 24px}h1{font-size:30px}h2{margin-top:32px}p,li{max-width:1150px}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:7px;border-bottom:1px solid #dbe3eb;text-align:left}th{background:#edf2f6}a{color:#176090}img{max-width:100%;height:auto}.note{background:#f2f5f8;padding:15px}.scroll{overflow:auto}select{font:inherit;max-width:95%;padding:8px}.region{border-top:1px solid #dbe3eb;margin-top:30px}@media print{.region{break-before:page}nav,select{display:none}body{padding:0;margin:0}a{color:inherit;text-decoration:none}}</style></head><body>',
    '<h1>Weekly symptom-to-visit-to-PCR rates by trial site</h1>',
    '<p class="note">Each heatmap cell is a documented visit + diagnostic PCR count divided by eligible reconstructed symptomatic episodes. Symptoms from AE records and COVID symptom assessments are combined. This is a documentation measure, not a determination that an investigator withheld testing or that an untested illness was COVID-19.</p>',
    '<nav><a href="all_sites_heatmaps.pdf">Two all-site heatmaps (PDF)</a> | <a href="regional_heatmaps.pdf">Local background + heatmaps (print PDF)</a> | <a href="site_arm_weekly_testing.csv">Every heatmap cell (CSV)</a> | <a href="episodes.csv">Episode audit (CSV)</a> | <a href="episode_visit_pcr_links.csv">Visit/PCR links (CSV)</a></nav>',
    if(file.exists(file.path(out,'report_insert.html')))'<p class="note"><a href="report_insert.html"><b>Copy-ready section for the current audit draft</b></a> — numbered paragraphs, tables, figures and references. <a href="report_insert.txt">Plain text</a>.</p>',
    '<h2>Counting and evidence rules</h2><ul>',paste0('<li>',hm_html(methods),'</li>'),'</ul>',
    if(file.exists(file.path(out,'site_anomaly_examples.pdf')))'<p><a href="site_anomaly_examples.pdf">Five illustrative site anomaly patterns (PDF, one figure per page)</a>. Exact counts and exploratory comparisons: <a href="anomaly_examples/README.txt">example methods</a> and <a href="anomaly_examples/example_selection.csv">selection record</a>. Early post-dose episodes remain in the main results.</p>',
    '<h2>Arm totals</h2><div class="scroll">',hm_table(summary),'</div>',
    '<h2>Input and coverage checks</h2>',hm_table(quality),
    '<p>Read <a href="source_manifest.csv">source_manifest.csv</a>, <a href="run_parameters.csv">run_parameters.csv</a>, <a href="site_region_coverage.csv">site_region_coverage.csv</a> and <a href="symptom_mapping/ae_covid_symptom_mapping_used.csv">the exact symptom mapping</a> before publication. Unknown dates, source conflicts and exclusions remain in the CSV audit files.</p>',
    '<h2>All-site overview</h2><p>The same site order and 0-100% color scale are used for both arms. Background series appear on the regional pages; no average across unrelated counties is substituted.</p>')
  if(pngs)html<-c(html,'<img loading="lazy" src="charts/all_sites_mrna1273.png" alt="All-site mRNA-1273 heatmap">',
    '<img loading="lazy" src="charts/all_sites_placebo.png" alt="All-site placebo heatmap">')
  html<-c(html,'<h2>Local COVID waves and testing rates</h2>',
    paste0('<p>County/state heading counts include all ',ws_n(nrow(ctx)),' randomized subjects, both arms combined, assigned to the sites in that region; the percentage uses the entire randomized trial roster. These are site locations, not participant residences, and are not the symptomatic-episode denominator. <a href="region_cohort_counts.csv">County/state cohort counts (CSV)</a>.</p>'),
    '<p>County-level background is preferred. If unavailable throughout the selected period, a labelled state series is used for that region. Within each page, both arms share identical case/positivity data and axis ranges. Background gaps remain gaps.</p>',
    paste0('<label for="region">Go to county/state: </label><select id="region" onchange="if(this.value)document.getElementById(this.value).scrollIntoView({behavior:\'smooth\'});"><option value="">Choose an area</option>',paste(nav,collapse=""),'</select>'))
  for(i in seq_len(nrow(regions))) {
    k<-regions$region_key[i]
    html<-c(html,paste0('<section class="region" id="r_',k,'"><h2>',hm_html(regions$display_name[i]),'</h2>'))
    if(pngs)html<-c(html,paste0('<img loading="lazy" src="charts/region_',k,'.png" alt="Community background and treatment-arm heatmaps">'))
    else html<-c(html,'<p>See regional_heatmaps.pdf for this chart.</p>')
    html<-c(html,'</section>')
  }
  html<-c(html,'</body></html>');writeLines(html,file.path(out,"index.html"),useBytes=TRUE)
}

analyse_symptom_testing_heatmaps <- function(root="",release="auto",window="BLINDED",
  adsl_file="",ae_file="",symptoms_file="",mb_file="",admb_file="",sv_file="",faot_file="",mapping_file="",
  background_file="analysis/background_covid/site_weekly.csv",expected_prod_cutoff="2021-03-26",eua_cutoff="2020-11-25",
  start_date="",end_date="",visit_evidence="RAW",gap_days=7L,onset_to_pcr_days=7L,visit_to_pcr_days=3L,
  ae_timing="interval",small_cell_n=5L,encoding="auto",numeric_date_origin="",min_randomized=30000L,
  output_dir="",write_png=TRUE,png_dpi=300L,show_report=interactive(),write_anomaly_examples=TRUE,anomaly_high_positivity=10,write_report_insert=TRUE,
  report_start=25L,reference_start=32L) {
  window<-match.arg(toupper(window),c("BLINDED","EUA","FULL"))
  visit_evidence<-match.arg(toupper(visit_evidence),c("RAW","INDEPENDENT"))
  ae_timing<-match.arg(ae_timing,c("interval","onset"))
  ints<-c(gap_days,onset_to_pcr_days,visit_to_pcr_days,small_cell_n,png_dpi)
  if(any(!is.finite(ints)|ints<0|ints!=floor(ints))||onset_to_pcr_days<1||small_cell_n<1||png_dpi<72)
    stop("Invalid day-window, small-cell or DPI setting.")
  root<-ws_root(root)
  background_path<-ws_path(root,background_file)
  if(!file.exists(background_path))stop("Run analyse_background_covid.R first: missing ",background_path)
  if(!nzchar(output_dir))output_dir<-file.path(root,"analysis",paste0("symptom_testing_heatmaps_",tolower(window),"_",format(Sys.time(),"%Y%m%d_%H%M%S")))
  out<-ws_path(root,output_dir)
  if(dir.exists(out)&&length(list.files(out,all.files=FALSE)))stop("Choose a new/empty PHMPT_OUTPUT_DIR; existing results are preserved: ",out)
  dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(file.path(out,"symptom_mapping"),showWarnings=FALSE)
  writeLines("RUNNING",file.path(out,"STATUS.txt"));complete<-FALSE
  on.exit(if(!complete)writeLines("FAILED - inspect source discovery logs and the R error",file.path(out,"STATUS.txt")),add=TRUE)
  if(tolower(release)=="auto")release<-if(length(ws_candidates(root,"Prod","ADSL")))"Prod" else "EUA"
  release<-if(tolower(release)=="prod")"Prod" else if(tolower(release)=="eua")"EUA" else stop("Release must be auto, Prod or EUA.")
  message("[1/7] Reading ",release," P301 source data")
  a<-ws_pick(root,release,"ADSL",adsl_file,out,encoding,numeric_date_origin,min_randomized)
  ctx<-ws_context(a$x,numeric_date_origin);submission<-ws_submission(a$manifest$path)
  if(release=="Prod"&&nzchar(expected_prod_cutoff)&&any(ctx$CUTOFFDT!=ws_dates(expected_prod_cutoff,"expected cutoff",strict=TRUE)))
    stop("Selected production ADSL is not the expected ",expected_prod_cutoff," snapshot. Set PHMPT_ADSL_FILE or review PHMPT_EXPECTED_PROD_CUTOFF.")
  from_id<-substr(ctx$subject,1,5);site<-ws_norm(ctx$SITEID)
  if(any(nzchar(site)&site!=from_id))stop("ADSL SITEID disagrees with the P301 subject-ID site prefix.")
  ctx$SITEID<-ifelse(nzchar(site),site,from_id)
  if(window!="FULL"&&max(ctx$CUTOFFDT)>as.numeric(as.Date("2020-12-18"))&&!any(ctx$has_blind_boundary))
    stop("Cannot determine participant-specific blinded follow-up. Supply complete ADSL.")
  ctx$follow_end<-switch(window,BLINDED=ctx$end_blinded,EUA=pmin(ctx$end_blinded,ws_dates(eua_cutoff,"EUA cutoff",strict=TRUE)),FULL=ctx$end_full)
  if(nzchar(end_date))ctx$follow_end<-pmin(ctx$follow_end,ws_dates(end_date,"end_date",strict=TRUE))
  ctx$follow_start<-ctx$TR01SDT
  if(!any(ctx$DOSE1FL=="Y"&ctx$follow_end>=ctx$TR01SDT,na.rm=TRUE))stop("No observed post-dose follow-up in the selected window.")
  first<-if(nzchar(start_date))ws_dates(start_date,"start_date",strict=TRUE) else min(ctx$TR01SDT,na.rm=TRUE)
  last<-max(ctx$follow_end,na.rm=TRUE)
  if(!is.finite(first)||!is.finite(last)||first>last)stop("No valid display date range.")
  weeks<-seq(ws_week(first),ws_week(last),by=7);sites<-sort(unique(ctx$SITEID))
  bg<-hm_background(background_path,sites,weeks,out)
  b<-ws_pick(root,release,"AE",ae_file,out,encoding,numeric_date_origin,min_randomized,submission)
  c<-ws_pick(root,release,"ADSYMP",symptoms_file,out,encoding,numeric_date_origin,min_randomized,submission)
  m<-ws_pick(root,release,"MB",mb_file,out,encoding,numeric_date_origin,min_randomized,submission)
  ad<-hm_pick_optional(root,release,"ADMB",admb_file,out,encoding,numeric_date_origin,min_randomized,submission)
  sv<-hm_pick_optional(root,release,"SV",sv_file,out,encoding,numeric_date_origin,min_randomized,submission)
  fa<-hm_pick_optional(root,release,"FAOT",faot_file,out,encoding,numeric_date_origin,min_randomized,submission)
  independent_available<-!is.null(sv$x)||!is.null(fa$x)
  if(visit_evidence=="INDEPENDENT"&&!independent_available)stop("INDEPENDENT mode requires raw SV and/or FAOT. Set PHMPT_SV_FILE/PHMPT_FAOT_FILE.")
  if(!independent_available)warning("SV/FAOT unavailable: primary RAW endpoint uses original MB illness labels only; independent-visit rates will be missing.",call.=FALSE,immediate.=TRUE)
  cd<-ws_dates(c$x$CUTOFFDT,"ADSYMP.CUTOFFDT",numeric_date_origin,TRUE)
  if(any(is.finite(cd)&cd!=ctx$CUTOFFDT[1]))stop("ADSYMP and ADSL cutoff dates disagree.")
  manifest<-do.call(rbind,lapply(list(a,b,c,m,ad,sv,fa),function(z)z$manifest))
  background_manifest<-data.frame(domain="BACKGROUND",release="Public surveillance",path=normalizePath(background_path,winslash="/"),
    md5=unname(tools::md5sum(background_path)),rows=NA_integer_,cutoff="",status="SELECTED",detail="Output of analyse_background_covid.R")
  manifest<-rbind(manifest,background_manifest);ws_csv(manifest,file.path(out,"source_manifest.csv"))
  message("[2/7] Mapping AE and COVID-form symptoms")
  mapping<-ws_mapping(root,mapping_file,file.path(out,"symptom_mapping"))
  events<-ws_events(b$x,c$x,mapping,numeric_date_origin,file.path(out,"symptom_mapping"),ae_timing)
  ws_csv(hm_dates_out(events),file.path(out,"symptom_event_evidence.csv"))
  ws_csv(hm_dates_out(ctx),file.path(out,"subject_followup.csv"))
  message("[3/7] Reconstructing symptom episodes before matching PCR/visit records")
  ep<-hm_episodes(events,ctx,gap_days,onset_to_pcr_days,first,last,out)
  message("[4/7] Linking original illness-visit evidence and diagnostic PCR records")
  pcr<-hm_tests(m$x,ad$x,numeric_date_origin,out)
  visits<-hm_visits(sv$x,fa$x,pcr,numeric_date_origin,out)
  ep<-hm_link(ep,pcr,visits,onset_to_pcr_days,visit_to_pcr_days,visit_evidence,independent_available,out)
  ws_csv(hm_dates_out(ep),file.path(out,"episodes.csv"))
  message("[5/7] Creating site/arm/week denominators and coverage")
  cells<-hm_cells(ep,ctx,sites,weeks,onset_to_pcr_days,small_cell_n,window)
  ix<-match(cells$SITEID,bg$map$SITEID)
  for(k in c("region_key","region_level","region_id","region_name"))cells[[k]]<-bg$map[[k]][ix]
  ws_csv(hm_dates_out(cells),file.path(out,"site_arm_weekly_testing.csv"))
  summary<-do.call(rbind,lapply(c("mRNA-1273","Placebo"),function(arm) {
    z<-ep[ep$arm==arm&ep$eligible,];data.frame(arm=arm,eligible_episodes=nrow(z),subjects=length(unique(z$subject)),
      visit_plus_pcr=sum(z$tested),testing_pct=if(nrow(z))round(100*mean(z$tested),2) else NA_real_,
      independent_visit_plus_pcr=if(independent_available)sum(z$independent_visit_pcr) else NA_integer_,
      diagnostic_pcr_regardless_of_visit=sum(z$diagnostic_pcr_by_horizon),
      incomplete_followup=sum(ep$arm==arm&ep$onset_in_display&!ep$preexisting&!ep$complete_followup))
  }));ws_csv(summary,file.path(out,"arm_totals.csv"))
  quality<-data.frame(check=c("Randomized subjects","Sites","Dated reconstructed episodes","Eligible episodes","Weeks",
    "Raw MB SARS-CoV-2 records","Usable diagnostic PCR records","Recorded illness-visit dates","Independent visit sources present",
    "Episodes spanning >28 days","Unknown AE-end episodes","Small nonempty site-arm-week denominators",
    "Sites missing some background positivity weeks"),value=as.character(c(nrow(ctx),length(sites),nrow(ep),sum(ep$eligible),length(weeks),
    nrow(pcr),sum(pcr$usable_diagnostic),nrow(visits),independent_available,sum(ep$span_days>28),sum(ep$unknown_ae_end),
    sum(cells$small_denominator),sum(bg$map$weeks_with_positivity<bg$map$weeks))))
  ws_csv(quality,file.path(out,"quality_checks.csv"))
  methods<-c(
    paste0("Population/window: ",release," / ",window,"; original randomized arm, dosed participants with observed follow-up. Displayed onset dates ",ws_iso(first)," to ",ws_iso(last),". EUA uses the later selected snapshot restricted by date, not an unrevised EUA database."),
    if(window=="FULL")"FULL includes post-unblinding/crossover follow-up in the original arm; those arms no longer represent consistently distinct treatments." else "Participant-specific unblinding, crossover, outside vaccination, end of study, death and data cutoff determine follow-up. Unblinding/crossover dates themselves are excluded when times are unavailable.",
    "Denominator: incident reconstructed symptomatic episodes from the union of mapped AE records and positive COVID-form/ADSYMP assessments. One episode can include several symptoms and either/both reporting sources; it counts once, in its first recorded symptom week. This differs from the earlier participant-by-symptom-group-by-week burden chart. Solicited diary-only reactions are not added in this script.",
    paste0("Episode rule: merge overlapping symptom intervals and intervals separated by at most ",gap_days," unrecorded calendar days. AE timing is ",ae_timing,"; an unknown end contributes onset only. COVID-form ADT is an assessment date, not necessarily actual illness onset. Pre-existing episodes are excluded. Absence of a record does not establish a symptom-free day; long episodes and unknown ends remain flagged."),
    paste0("Complete observation through onset day +",onset_to_pcr_days," is required for every denominator episode, including ones tested early. Days 0 through ",onset_to_pcr_days," are inclusive. Episodes with another onset within that horizon are excluded, preventing ambiguous reuse of a test under shorter-gap sensitivity settings. Partial final weeks and censored episodes are flagged, not counted as untested."),
    paste0("Numerator: diagnostic PCR with documented result and a recorded illness-visit association on or after episode onset; PCR collection on the visit date or up to ",visit_to_pcr_days," calendar days later, and no later than onset +",onset_to_pcr_days,". Same-calendar-day linkage is allowed because clock times are not uniformly available; ordering within that day cannot be established. These are analyst linkage rules, not the protocol's complete suspected-case criteria."),
    paste0("Visit mode: ",visit_evidence,". RAW uses original SV illness-visit dates, FAOT PHYSICAL EXAMINATION/ASSESS illness-visit records, or the original raw MB illness-visit label/date. INDEPENDENT uses SV/FAOT only. FAOT assessment N can still mark a recorded visit with an examination not performed; neither mode establishes independently verified attendance. ADMB VISIT/AVISIT and LOCALFL alone never create a visit."),
    "The RAW endpoint can use the visit label on the PCR record itself: this is an explicit documented association, not independent visit confirmation. The CSV also supplies independent_visit_pcr_pct and diagnostic-PCR-without-visit counts. If SV/FAOT are unavailable, independent rates are missing, not zero. Missing visit documentation can lower the main rate even when a PCR is documented.",
    "PCR identity and method come from raw MBTESTCD=SARSCOV2 and MBMETHOD containing PCR/POLYMERASE, without explicit antigen/serology/antibody method contradiction. ADMB is reconciled by subject + MBSEQ for results, baseline and local-coded flags. Raw collection dates are required; conflicting dates are excluded. Positive, negative and explicit indeterminate/inconclusive/invalid/equivocal assay results establish testing. NOT DONE, comment-only and unresolved values do not. A disagreement about the result is flagged but does not erase evidence an assay was performed. Free-text laboratory comments are not a full independent validation of laboratory method.",
    "Diagnostic eligibility: original raw MB illness label or group 003/LOCALFL=Y, excluding baseline records. A routine scheduled swab alone does not satisfy the plotted endpoint. All nonbaseline PCR evidence is exported separately. Multiple result representations/test records cannot increase an episode's binary numerator beyond one; distinct MBSEQ values do not by themselves establish distinct physical specimens.",
    "Red means low documented testing proportion; green means high. Both arms share a fixed 0-100% scale and site order. Zero symptoms is light gray, no eligible observation opportunity is darker gray, and censored-only episodes are hatched. n/N labels give tested/eligible episodes. An asterisk or dot flags small denominators; 1/1 is not strong evidence of consistent practice.",
    "The county-specific background from the previous script is plotted above each arm's regional heatmap: reported cases as columns, positivity as a line. Both arms see exactly the same community series. State fallback is selected for an entire region when county positivity is unavailable throughout the displayed period; county/state percentages are not silently spliced week by week or averaged across unrelated regions. Case and positivity axes are independent.",
    "Community sources retain the previous script's definitions: JHU archived net reported-case changes; CDC published county 7-day NAAT positivity; HHS state positive/(positive+negative) PCR results, inconclusive excluded. No county test denominator is available in that CDC series. Site county is a proxy, not participants' residence. Public testing includes screening and is not directly comparable to symptom-triggered trial testing.",
    "These are unadjusted descriptive rates. Symptoms, severity, vaccination timing, follow-up and reporting-source mix can vary by arm and calendar time. The CSV includes early post-dose and outside-post-dose counts to support follow-up analyses. A visual association with a COVID wave does not establish differential investigator treatment, missed infections or a causal vaccine effect.")
  writeLines(methods,file.path(out,"METHODS.txt"),useBytes=TRUE)
  pars<-data.frame(parameter=c("version","release","window","visit_evidence","gap_days","onset_to_pcr_days","visit_to_pcr_days",
    "ae_timing","small_cell_n","display_start","display_end","independent_visit_sources","background_file","generated_utc"),
    value=c("1.0.4",release,window,visit_evidence,gap_days,onset_to_pcr_days,visit_to_pcr_days,ae_timing,small_cell_n,
    ws_iso(first),ws_iso(last),paste(c(if(!is.null(sv$x))"SV",if(!is.null(fa$x))"FAOT"),collapse=";"),
    normalizePath(background_path,winslash="/"),format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC")))
  ws_csv(pars,file.path(out,"run_parameters.csv"))
  if(write_anomaly_examples || write_report_insert)phmpt_anomaly_examples(out,out,high_cutoff=anomaly_high_positivity)
  if(write_report_insert)phmpt_report_insert(out,out,report_start,reference_start,audit_root=root)
  message("[6/7] Drawing the two overview heatmaps and regional print pages")
  h<-max(11,3+length(sites)*.19)
  hm_pdf(file.path(out,"all_sites_heatmaps.pdf"),15.5,h)
  tryCatch(for(arm in c("mRNA-1273","Placebo"))hm_overview_page(cells,sites,weeks,arm,window,small_cell_n),finally=grDevices::dev.off())
  regs<-unique(bg$map[c("region_key","region_name")]);regs<-regs[order(regs$region_name),]
  hm_pdf(file.path(out,"regional_heatmaps.pdf"),11.7,8.3)
  tryCatch(for(r in regs$region_key)hm_region_page(cells,bg$map,bg$area,weeks,r,window,visit_evidence,small_cell_n,onset_to_pcr_days,visit_to_pcr_days),finally=grDevices::dev.off())
  if(write_png) {
    dir.create(file.path(out,"charts"),showWarnings=FALSE)
    for(arm in c("mRNA-1273","Placebo")) {
      file<-if(arm=="mRNA-1273")"all_sites_mrna1273.png" else "all_sites_placebo.png"
      hm_png(file.path(out,"charts",file),15.5,h,png_dpi,function()hm_overview_page(cells,sites,weeks,arm,window,small_cell_n))
    }
    for(r in regs$region_key)hm_png(file.path(out,"charts",paste0("region_",r,".png")),11.7,8.3,png_dpi,
      function()hm_region_page(cells,bg$map,bg$area,weeks,r,window,visit_evidence,small_cell_n,onset_to_pcr_days,visit_to_pcr_days))
  }
  message("[7/7] Writing report, audit links and completion marker")
  hm_report(out,cells,bg,summary,methods,write_png,manifest,quality,ctx)
  writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
  outputs<-list.files(out,recursive=TRUE,full.names=TRUE);outputs<-outputs[basename(outputs)!="STATUS.txt"]
  ws_csv(data.frame(file=substring(outputs,nchar(out)+2),bytes=file.info(outputs)$size,md5=unname(tools::md5sum(outputs))),file.path(out,"output_checksums.csv"))
  complete<-TRUE;writeLines("COMPLETE",file.path(out,"STATUS.txt"))
  print(summary,row.names=FALSE);message("Report: ",normalizePath(file.path(out,"index.html"),winslash="/"))
  if(show_report)utils::browseURL(normalizePath(file.path(out,"index.html"),winslash="/"))
  invisible(list(output_dir=out,cells=cells,episodes=ep,pcr=pcr,visits=visits,background=bg,summary=summary,context=ctx))
}

# Register extended readers and prefer the March 2021 S3 BLA submission when it
# is identifiable in filenames. Explicit source paths always remain available.
ws_specs <- hm_install_specs()
.hm_original_candidates <- ws_candidates
ws_candidates <- function(root,release,domain) {
  p<-.hm_original_candidates(root,release,domain)
  if(release=="Prod") {
    s3<-grepl("_s3_m[0-9]+_",basename(p),ignore.case=TRUE)
    if(any(s3))p<-p[s3]
  }
  p
}

# ---- RUN: selecting/running the whole script executes the analysis ----
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

phmpt_export_audit_assets <- function(root,out,data_dir,audit_tables) {
  a<-ax_begin(root,out,"V15_V23_episodes")
  fp<-ax_read(file.path(out,"anomaly_examples"),"figure_parameters.csv")
  setting<-function(key,default){z<-fp$value[fp$parameter==key];if(length(z)==1)as.character(z) else default}
  cut<-as.numeric(setting("high_positivity_threshold","10"));horizon<-as.numeric(setting("onset_to_pcr_days","7"))
  window<-setting("window","UNKNOWN")
  high_label<-paste0("High positivity >=",cut,"% | Low <",cut,"%")
  gettab<-function(name)if(name%in%names(audit_tables))audit_tables[[name]] else NULL
  specs<-list(
    V15=c("Documented testing after symptomatic episodes","episode_testing_overview"),
    V16=c("Adjusted comparisons during higher background positivity","adjusted_testing_models"),
    V17=c("Sites with nominal testing imbalances during higher positivity","site_testing_screen"),
    V22=c("Sites with opposing changes across positivity periods","site_positivity_interactions"))
  for(id in names(specs)) {
    s<-specs[[id]];d<-gettab(s[1])
    if(is.null(d)){ax_unavailable(a,id,s[2],s[1],if(id=="V16")"No converged, input-hash-matched adjusted models. Run analyse_site_testing_imbalance.R on this export." else "No matching comparison in this run.");next}
    if(id=="V22") {
      d<-do.call(rbind,lapply(seq_len(nrow(d)),function(i)data.frame(Site=d[i,1],Arm=c("mRNA-1273","Placebo"),
        `Low positivity`=c(d[i,2],d[i,4]),`High positivity`=c(d[i,3],d[i,5]),`Interaction p`=c(d[i,6],""),check.names=FALSE)))
    }
    note<-switch(id,V15="Endpoint: documented COVID visit plus qualifying PCR in the configured linkage window. Repeated episodes can contribute. Differences are vaccine minus placebo; missing background remains in overall totals.",
      V16="Site-clustered uncertainty. Odds ratios compare testing in vaccine versus placebo episodes. Onset recording source is a sensitivity adjustment. Early post-dose symptoms remain in the main comparison.",
      V17="First eligible high-positivity episode per participant. Fisher tests are two-sided; BH q uses the complete site family. Nominal signals select records for review.",
      V22="Low/high are pooled positivity categories, not consecutive weeks. Every arm-period cell has at least ten episodes. Interaction p values account for repeated participants and are nominal.")
    ax_table(a,id,s[2],s[1],d,paste0("COVID forms plus matched AEs | ",window," | ",high_label),note,
      widths=switch(id,V15=c(.30,.25,.25,.20),V16=c(.45,.16,.23,.16),V17=c(.30,.22,.22,.13,.13),V22=c(.15,.20,.25,.25,.15)))
  }
  ex<-file.path(out,"anomaly_examples");sel<-ax_read(ex,"example_selection.csv",FALSE)
  if(!nrow(sel))sel<-ax_read(ex,"selected_examples.csv",FALSE)
  site<-function(kind,fallback){z<-sel$SITEID[sel$anomaly==kind];if(length(z)==1&&nzchar(z))z else fallback}
  chart<-function(id,slug,title,subtitle,note,fun,height=6.6){
    draw<-function(){old<-par(no.readonly=TRUE);on.exit(par(old));par(fig=c(0,1,.18,.80),mai=c(.65,.85,.3,.35),family="sans",las=1,bty="n",cex=1)
      fun();grid::pushViewport(grid::viewport());ax_text(ax_wrap(paste(id,title),7.8,17,TRUE),.35,.24,17,TRUE,col="#17384B")
      ax_text(ax_wrap(subtitle,7.8,11.2),.35,.82,11.2)
      ax_text(ax_wrap(note,7.8,10.6),.35,height-.95,10.6,col="#52616C");grid::popViewport()}
    ax_png(a,id,slug,draw,height)
  }
  col<-c("#1E667D","#A65B22")
  bars<-function(k,n,labels,main="") {
    pct<-ifelse(n>0,100*k/n,NA_real_);x<-barplot(pct,ylim=c(0,116),col=rep(col,length.out=length(k)),border=NA,axes=FALSE,names.arg=labels,cex.names=.86,space=.55)
    axis(2,at=seq(0,100,20),labels=paste0(seq(0,100,20),"%"),tick=FALSE,cex.axis=.86)
    abline(h=0,col="#9EA9AF");text(x,pmin(108,ifelse(is.finite(pct),pct,0)+6),ax_rate(k,n),cex=.84,font=2,xpd=NA)
    if(nzchar(main))mtext(main,3,adj=0,line=.6,cex=.9);invisible(x)
  }
  z<-ax_read(ex,"first_high_episode_site_screen.csv");s<-site("coverage","US372");z<-z[z$SITEID==s,]
  if(nrow(z)==1){chart("V18","site_testing_example",paste(s,"Documented testing after symptoms"),paste0("First eligible high-positivity episode per participant | ",high_label),
    paste0("Visit-plus-PCR endpoint. Two-sided Fisher p ",ax_p(z$fisher_p),"; BH q ",ax_p(z$fisher_q_BH)," across the complete site family."),
    function()bars(c(z$v_y,z$p_y),c(z$v_n,z$p_n),c("mRNA-1273","Placebo")));ax_data(a,"V18",z)
  }else ax_unavailable(a,"V18","site_testing_example","Site testing example","Selected site not available.")
  z<-ax_read(ex,"pcr_episode_review.csv");s<-site("pcr","US325")
  if(nrow(z)) {
    chart("V19","delayed_pcr_example",paste(s,"PCR timing after COVID-form symptoms"),paste0("Vaccine-arm high-positivity episodes without the qualifying endpoint within ",horizon," days"),
      paste0("Blue spans: reconstructed symptom interval. Circle: later completed PCR (labelled result). Cross: NOT DONE within ",horizon," days. Shading: audit linkage interval. Missing later tests remain explicitly labelled."),function(){
      nn<-nrow(z);late<-z$days_to_later_pcr;duration<-as.numeric(as.Date(z$last_symptom)-as.Date(z$start));limit<-max(c(22,late,duration),na.rm=TRUE)+2
      par(mai=c(.65,1.65,.3,.45));plot(c(-1,limit),c(.3,nn+.8),type="n",axes=FALSE,xlab="",ylab="")
      rect(0,0,horizon,nn+1,col="#EDF1F3",border=NA);axis(1,at=pretty(c(0,limit)),cex.axis=.85);mtext("Days after first recorded symptom",1,line=2,cex=.86)
      yy<-rev(seq_len(nn));axis(2,at=yy,labels=z$subject,tick=FALSE,cex.axis=.79)
      for(i in seq_len(nn)) {
        segments(0,yy[i],duration[i],yy[i],col=col[1],lwd=5);points(0,yy[i],pch=16,col=col[1],cex=.6)
        if(is.finite(late[i])){points(late[i],yy[i],pch=21,bg="white",cex=1.1);text(late[i],yy[i]+.23,paste0("Day ",late[i],": ",tolower(z$later_result[i])),cex=.76,adj=c(.5,0),xpd=NA)}
        else text(8,yy[i]+.23,"No later completed PCR",adj=c(0,0),cex=.76)
        if(isTRUE(z$not_done_in_horizon[i])){nd<-suppressWarnings(as.numeric(z$first_not_done_day[i]));if(is.finite(nd))points(nd,yy[i],pch=4,lwd=1.8,cex=1.2)}
      }
    },height=max(6.6,3.5+.7*nrow(z)));ax_data(a,"V19",z)
  }else ax_unavailable(a,"V19","delayed_pcr_example","PCR timing review","No selected episodes.")
  z<-ax_read(ex,"recording_source_counts.csv");s<-site("source","US352")
  d<-data.frame(Arm=z$arm,`Onset source`=z$source,`Share of episodes`=paste0(z$episodes," (",sprintf("%.1f%%",z$share_pct),")"),
    `Tested within source`=ax_rate(z$tested,z$episodes),check.names=FALSE)
  # A compact table gives both the composition and the within-source denominator.
  ax_table(a,"V20","source_composition_example",paste(s,"Recording source and documented testing"),d,
    paste0("All eligible episodes | Source at first symptom date | ",high_label),
    "AE = AE records only; COVID_FORM = COVID forms only; AE; COVID_FORM = both on the first symptom date. The share column describes source composition; the final column describes visit-plus-PCR testing within that source.",widths=c(.19,.24,.27,.30))
  z<-ax_read(ex,"postdose_timing_counts.csv");s<-site("timing","US372")
  if(nrow(z)) {
    chart("V21","postdose_timing_example",paste(s,"Testing by post-dose symptom timing"),paste0("All eligible episodes | Repeated episodes retained | ",high_label),
      "Each pair shows mRNA-1273 (blue) and placebo (brown). Post-dose days 0-6 remain in the main comparison. The partition localizes the difference without identifying its cause.",function(){
      labs<-rep(c("Vaccine","Placebo"),nrow(z));k<-as.vector(rbind(z$v_y,z$p_y));n<-as.vector(rbind(z$v_n,z$p_n));x<-bars(k,n,labs)
      for(i in seq_len(nrow(z)))mtext(z$group[i],1,at=mean(x[(2*i-1):(2*i)]),line=2.2,cex=.82)
    });ax_data(a,"V21",z)
  }else ax_unavailable(a,"V21","postdose_timing_example","Testing by post-dose timing","No selected episodes.")
  z<-ax_read(ex,"site_high_low_interactions.csv");s<-site("waves","US302");z<-z[z$subset=="ALL_EPISODES"&z$SITEID==s,]
  if(nrow(z)==1) {
    chart("V23","positivity_testing_example",paste(s,"Testing in lower and higher positivity"),paste0("All eligible episodes | Pooled categories | ",high_label),
      paste0("Vaccine change minus placebo change: ",sprintf("%+.1f",z$difference_in_changes_pp)," percentage points. Participant-clustered interaction p ",ax_p(z$interaction_p)," (nominal). Lines connect pooled categories, not successive calendar weeks."),function(){
      plot(c(.65,2.35),c(0,106),type="n",axes=FALSE,xlab="",ylab="");abline(h=seq(0,100,20),col="#E1E7EB")
      axis(2,at=seq(0,100,20),labels=paste0(seq(0,100,20),"%"),tick=FALSE,cex.axis=.85);axis(1,at=1:2,labels=c("Lower positivity","Higher positivity"),tick=FALSE,cex.axis=.85)
      for(i in 1:2){p<-if(i==1)"v" else "p";k<-unlist(z[paste0(p,c("_low_y","_high_y"))],use.names=FALSE);n<-unlist(z[paste0(p,c("_low_n","_high_n"))],use.names=FALSE)
        y<-100*k/n;lines(1:2,y,col=col[i],lwd=2);points(1:2,y,pch=16,col=col[i]);other<-if(i==1)100*c(z$p_low_y/z$p_low_n,z$p_high_y/z$p_high_n) else 100*c(z$v_low_y/z$v_low_n,z$v_high_y/z$v_high_n);text(1:2,y+ifelse(y>=other,8,-8),ax_rate(k,n),col=col[i],cex=.87,font=2)}
      legend("topright",c("mRNA-1273","Placebo"),col=col,lwd=2,bty="n",cex=.86)
    });ax_data(a,"V23",z)
  }else ax_unavailable(a,"V23","positivity_testing_example","Testing across positivity periods","Selected site unavailable.")
  ax_finish(a)
}

if(isTRUE(getOption("phmpt.autorun",TRUE))) {
  PHMPT_HEATMAP_RESULT <- analyse_symptom_testing_heatmaps(root=PHMPT_ROOT,release=PHMPT_RELEASE,window=PHMPT_WINDOW,
    adsl_file=PHMPT_ADSL_FILE,ae_file=PHMPT_AE_FILE,symptoms_file=PHMPT_SYMPTOMS_FILE,mb_file=PHMPT_MB_FILE,
    admb_file=PHMPT_ADMB_FILE,sv_file=PHMPT_SV_FILE,faot_file=PHMPT_FAOT_FILE,mapping_file=PHMPT_MAPPING_FILE,
    background_file=PHMPT_BACKGROUND_FILE,expected_prod_cutoff=PHMPT_EXPECTED_PROD_CUTOFF,eua_cutoff=PHMPT_EUA_CUTOFF,
    start_date=PHMPT_START_DATE,end_date=PHMPT_END_DATE,visit_evidence=PHMPT_VISIT_EVIDENCE,gap_days=PHMPT_GAP_DAYS,
    onset_to_pcr_days=PHMPT_ONSET_TO_PCR_DAYS,visit_to_pcr_days=PHMPT_VISIT_TO_PCR_DAYS,ae_timing=PHMPT_AE_TIMING,
    small_cell_n=PHMPT_SMALL_CELL_N,encoding=PHMPT_CSV_ENCODING,numeric_date_origin=PHMPT_NUMERIC_DATE_ORIGIN,
    min_randomized=PHMPT_MIN_RANDOMIZED,output_dir=PHMPT_OUTPUT_DIR,write_png=PHMPT_WRITE_PNG,png_dpi=PHMPT_PNG_DPI,
    show_report=PHMPT_SHOW_REPORT,write_anomaly_examples=PHMPT_WRITE_ANOMALY_EXAMPLES,
    anomaly_high_positivity=PHMPT_ANOMALY_HIGH_POSITIVITY,write_report_insert=PHMPT_WRITE_REPORT_INSERT,
    report_start=PHMPT_REPORT_START,reference_start=PHMPT_REFERENCE_START)
}
