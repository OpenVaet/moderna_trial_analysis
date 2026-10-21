# =============================================================================
# Moderna P301 - weekly symptom-to-visit-to-PCR heatmaps | v1.0.1 | 2026-10-21
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
# post-unblinding/crossover data and is descriptive only. No significance tests.
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
hm_report <- function(out,cells,bg,summary,methods,pngs,manifest,quality) {
  regions<-unique(bg$map[c("region_key","region_name")]);regions<-regions[order(regions$region_name),]
  nav<-paste0('<option value="r_',regions$region_key,'">',hm_html(regions$region_name),'</option>')
  html<-c('<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Weekly symptom testing by site</title>',
    '<style>body{font:16px/1.5 Arial,sans-serif;color:#253448;max-width:1300px;margin:30px auto;padding:0 24px}h1{font-size:30px}h2{margin-top:32px}p,li{max-width:1150px}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:7px;border-bottom:1px solid #dbe3eb;text-align:left}th{background:#edf2f6}a{color:#176090}img{max-width:100%;height:auto}.note{background:#f2f5f8;padding:15px}.scroll{overflow:auto}select{font:inherit;max-width:95%;padding:8px}.region{border-top:1px solid #dbe3eb;margin-top:30px}@media print{.region{break-before:page}nav,select{display:none}body{padding:0;margin:0}a{color:inherit;text-decoration:none}}</style></head><body>',
    '<h1>Weekly symptom-to-visit-to-PCR rates by trial site</h1>',
    '<p class="note">Each heatmap cell is a documented visit + diagnostic PCR count divided by eligible reconstructed symptomatic episodes. Symptoms from AE records and COVID symptom assessments are combined. This is a documentation measure, not a determination that an investigator withheld testing or that an untested illness was COVID-19.</p>',
    '<nav><a href="all_sites_heatmaps.pdf">Two all-site heatmaps (PDF)</a> | <a href="regional_heatmaps.pdf">Local background + heatmaps (print PDF)</a> | <a href="site_arm_weekly_testing.csv">Every heatmap cell (CSV)</a> | <a href="episodes.csv">Episode audit (CSV)</a> | <a href="episode_visit_pcr_links.csv">Visit/PCR links (CSV)</a></nav>',
    '<h2>Counting and evidence rules</h2><ul>',paste0('<li>',hm_html(methods),'</li>'),'</ul>',
    '<h2>Arm totals</h2><div class="scroll">',hm_table(summary),'</div>',
    '<h2>Input and coverage checks</h2>',hm_table(quality),
    '<p>Read <a href="source_manifest.csv">source_manifest.csv</a>, <a href="run_parameters.csv">run_parameters.csv</a>, <a href="site_region_coverage.csv">site_region_coverage.csv</a> and <a href="symptom_mapping/ae_covid_symptom_mapping_used.csv">the exact symptom mapping</a> before publication. Unknown dates, source conflicts and exclusions remain in the CSV audit files.</p>',
    '<h2>All-site overview</h2><p>The same site order and 0-100% color scale are used for both arms. Background series appear on the regional pages; no average across unrelated counties is substituted.</p>')
  if(pngs)html<-c(html,'<img loading="lazy" src="charts/all_sites_mrna1273.png" alt="All-site mRNA-1273 heatmap">',
    '<img loading="lazy" src="charts/all_sites_placebo.png" alt="All-site placebo heatmap">')
  html<-c(html,'<h2>Local COVID waves and testing rates</h2>',
    '<p>County-level background is preferred. If unavailable throughout the selected period, a labelled state series is used for that region. Within each page, both arms share identical case/positivity data and axis ranges. Background gaps remain gaps.</p>',
    paste0('<label for="region">Go to county/state: </label><select id="region" onchange="if(this.value)document.getElementById(this.value).scrollIntoView({behavior:\'smooth\'});"><option value="">Choose an area</option>',paste(nav,collapse=""),'</select>'))
  for(i in seq_len(nrow(regions))) {
    k<-regions$region_key[i]
    html<-c(html,paste0('<section class="region" id="r_',k,'"><h2>',hm_html(regions$region_name[i]),'</h2>'))
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
  output_dir="",write_png=TRUE,png_dpi=300L,show_report=interactive()) {
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
    value=c("1.0.1",release,window,visit_evidence,gap_days,onset_to_pcr_days,visit_to_pcr_days,ae_timing,small_cell_n,
    ws_iso(first),ws_iso(last),paste(c(if(!is.null(sv$x))"SV",if(!is.null(fa$x))"FAOT"),collapse=";"),
    normalizePath(background_path,winslash="/"),format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC")))
  ws_csv(pars,file.path(out,"run_parameters.csv"))
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
  hm_report(out,cells,bg,summary,methods,write_png,manifest,quality)
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
if(isTRUE(getOption("phmpt.autorun",TRUE))) {
  PHMPT_HEATMAP_RESULT <- analyse_symptom_testing_heatmaps(root=PHMPT_ROOT,release=PHMPT_RELEASE,window=PHMPT_WINDOW,
    adsl_file=PHMPT_ADSL_FILE,ae_file=PHMPT_AE_FILE,symptoms_file=PHMPT_SYMPTOMS_FILE,mb_file=PHMPT_MB_FILE,
    admb_file=PHMPT_ADMB_FILE,sv_file=PHMPT_SV_FILE,faot_file=PHMPT_FAOT_FILE,mapping_file=PHMPT_MAPPING_FILE,
    background_file=PHMPT_BACKGROUND_FILE,expected_prod_cutoff=PHMPT_EXPECTED_PROD_CUTOFF,eua_cutoff=PHMPT_EUA_CUTOFF,
    start_date=PHMPT_START_DATE,end_date=PHMPT_END_DATE,visit_evidence=PHMPT_VISIT_EVIDENCE,gap_days=PHMPT_GAP_DAYS,
    onset_to_pcr_days=PHMPT_ONSET_TO_PCR_DAYS,visit_to_pcr_days=PHMPT_VISIT_TO_PCR_DAYS,ae_timing=PHMPT_AE_TIMING,
    small_cell_n=PHMPT_SMALL_CELL_N,encoding=PHMPT_CSV_ENCODING,numeric_date_origin=PHMPT_NUMERIC_DATE_ORIGIN,
    min_randomized=PHMPT_MIN_RANDOMIZED,output_dir=PHMPT_OUTPUT_DIR,write_png=PHMPT_WRITE_PNG,png_dpi=PHMPT_PNG_DPI,
    show_report=PHMPT_SHOW_REPORT)
}
