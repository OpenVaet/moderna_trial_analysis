# Moderna P301: protocol-triggered illnesses, investigation and diagnostic specimen collection.
# Standalone BASE R. Run whole file in RStudio. Defaults work in the supplied bundle's root.
# Otherwise set options(trigger.heat="...zip", trigger.burden="...zip", trigger.lab="...zip",
# trigger.adsl=".../BLA_ADSL_selected.csv", phmpt.root=getwd()). Directories also accepted.
# Optional verified source overlays and their schemas are described in README.md.
# No source timestamps are invented. Missing records do not establish a breach.
# Do not mix a legacy event export with a strict M3 result.
m3_valid_export <- function(path) {
 if(!length(path)||!nzchar(path)||!file.exists(path))return(FALSE)
 if(dir.exists(path)){f<-file.path(path,"M3_mapping_version.txt");return(file.exists(f)&&identical(trimws(readLines(f,n=1L,warn=FALSE)),"1.1.0-M3"))}
 if(!grepl("[.]zip$",path,ignore.case=TRUE))return(FALSE)
 n<-utils::unzip(path,list=TRUE)$Name;n<-n[grepl("^([^/]+/)?M3_mapping_version[.]txt$",n)]
 if(length(n)!=1L)return(FALSE)
 con<-unz(path,n);on.exit(close(con));identical(trimws(readLines(con,n=1L,warn=FALSE)),"1.1.0-M3")
}
m3_find_export <- function(root, signature) {
 d<-file.path(root,"analysis");if(!dir.exists(d))return("")
 f<-list.files(d,"^M3_mapping_version[.]txt$",recursive=TRUE,full.names=TRUE);v<-dirname(f)
 v<-v[file.exists(file.path(v,signature)) & vapply(v,m3_valid_export,logical(1))]
 if(identical(signature,"episode_symptom_source_links.csv")&&length(v)){
  blinded<-function(d){f<-file.path(d,"run_parameters.csv");if(!file.exists(f))return(FALSE)
    z<-read.csv(f,fileEncoding="UTF-8-BOM");any(z$parameter=="window"&z$value=="BLINDED")}
  v<-v[vapply(v,blinded,logical(1))]
 }
 if(!length(v))return("")
 normalizePath(v[order(file.info(file.path(v,signature))$mtime,v,decreasing=TRUE)[1L]],winslash="/")
}
m3_require_export <- function(x, producer) {
 if(!m3_valid_export(x))stop("Run ",producer," with the M3 correction first; this input has no verified M3 mapping marker: ",x,call.=FALSE)
 invisible(x)
}

tr_run <- function(root=getOption("phmpt.root",getwd()),out=file.path(root,"analysis","protocol_trigger_audit")){
 path<-function(x){if(file.exists(x))return(normalizePath(x));if(file.exists(file.path(root,x)))return(normalizePath(file.path(root,x)));stop("Input not found: ",x)}
 pick<-function(opt,stem,inc){u<-getOption(opt,"");if(nzchar(u))return(path(u));if(file.exists(file.path(root,"inputs",inc)))return(path(file.path(root,"inputs",inc)))
  f<-unique(unlist(lapply(c(root,file.path(root,"analysis"),file.path(root,"reference_v1_12")),function(d)if(dir.exists(d))list.files(d,pattern=stem,full.names=TRUE))))
  if(length(f)!=1L)stop("Set options(",opt,"='path') to one input. Candidates: ",paste(f,collapse=", "));path(f)}
 heat<-if(!nzchar(getOption("trigger.heat","")))m3_find_export(root,"episode_symptom_source_links.csv") else ""
 if(!nzchar(heat))heat<-pick("trigger.heat","^symptom_testing_heatmaps_blinded_","heat")
 m3_require_export(heat,"analyse_symptom_testing_heatmaps.R")
 burden<-pick("trigger.burden","^prod_symptom_burden_","burden");lab<-pick("trigger.lab","^prod_pcr_lab_metadata_","lab")
 adsl<-pick("trigger.adsl","^BLA_ADSL_selected[.]csv$","BLA_ADSL_selected.csv")
 dir.create(out,showWarnings=FALSE,recursive=TRUE);dir.create(file.path(out,"restricted"),showWarnings=FALSE)
 save<-function(z,f)write.csv(z,file.path(out,f),row.names=FALSE,na="",fileEncoding="UTF-8")
 readf<-function(f)read.csv(f,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM",na.strings=c("","NA"))
 get<-function(b,n){if(dir.exists(b)){p<-list.files(b,recursive=TRUE);m<-p[basename(p)==n];if(length(m)!=1L)stop("Ambiguous/missing ",n);return(readf(file.path(b,m)))}
  p<-utils::unzip(b,list=TRUE)$Name;m<-p[basename(p)==n];if(length(m)!=1L)stop("Ambiguous/missing ",n);d<-tempfile();dir.create(d);on.exit(unlink(d,recursive=TRUE))
  utils::unzip(b,files=m,exdir=d,junkpaths=TRUE);readf(file.path(d,n))}
 ck<-function(ok,msg){if(!isTRUE(all(ok)))stop(msg)}
 key<-function(...)paste(...,sep="|");dy<-function(x)as.integer(as.Date(as.character(x),format="%Y-%m-%d"))
 dayfmt<-function(x)as.character(as.Date(as.numeric(x),origin="1970-01-01"))
 iso<-function(x){x<-as.character(x);ok<-!is.na(x)&grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}(:[0-9.]+)?$",x);ans<-rep(NA_real_,length(x));ans[ok]<-as.numeric(as.POSIXct(sub("T"," ",x[ok],fixed=TRUE),tz="UTC"));ans}
 # Timestamp values are site-local wall times, not conversions from different zones.
 optfile<-function(opt){u<-getOption(opt,"");if(nzchar(u))readf(path(u)) else NULL}
 ep<-get(heat,"episodes.csv"); ev<-get(heat,"episode_symptom_source_links.csv"); m<-get(heat,"pcr_records_reconciled.csv")
 visit<-get(heat,"illness_visits.csv"); mb<-get(lab,"mb_pcr_metadata.csv"); grade<-get(burden,"positive_symptom_records.csv")
 ae<-get(burden,"ae_term_mapping.csv"); am<-get(heat,"ae_matching_evidence.csv"); a<-readf(adsl)
 ck(!anyDuplicated(ep$episode_id)&&!anyDuplicated(ev$event_id)&&!anyDuplicated(m$test_id),"Duplicate identifiers")
 ck(!anyDuplicated(a$SUBJID),"Duplicate ADSL IDs"); j<-match(ep$subject,a$SUBJID);ck(!anyNA(j),"ADSL mismatch")
 ep$pp<-a$PPROTFL[j]=="Y" & a$SCOV2BL[j]=="Negative"; ep$effend<-dy(a$EFFCODT[j]);ep$d1<-dy(ep$first_dose);ep$d2<-dy(ep$dose2)
 ep$s<-dy(ep$start);ep$fe<-dy(ep$follow_end);ep$fs<-dy(ep$follow_start); ep$last<-dy(ep$last_symptom)
 k<-match(ev$episode_id,ep$episode_id);ck(!anyNA(k)&ev$subject==ep$subject[k],"Symptom join mismatch")
 ev$l<-dy(ev$lo);ev$h<-dy(ev$hi);ev$d1<-ep$d1[k];ev$d2<-ep$d2[k];ev$site<-ep$SITEID[k]
 ev$h<-pmin(ev$h,ep$fe[k],na.rm=FALSE); ev$conflict<-grepl("NEGATIVE_SAME_DAY",ev$date_quality)
 # Verify AE row linkage by subject, sequence/row and verbatim term before enrichment.
 am<-unique(am[,c("USUBJID","AESEQ","row_id")]);am$subject<-sub("^.*-(US[0-9]{3})-([0-9]+)$","\\1\\2",am$USUBJID)
 ck(!anyDuplicated(key(am$subject,am$AESEQ)),"Nonunique AE sequence")
 ae$row<-as.integer(sub("^.*:","",ae$source_ref));ck(!anyDuplicated(key(ae$subject,ae$row)),"Nonunique AE raw source row")
 ii<-which(ev$source=="AE");mr<-am$row_id[match(key(ev$subject[ii],ev$source_sequence[ii]),key(am$subject,am$AESEQ))]; ar<-match(key(ev$subject[ii],mr),key(ae$subject,ae$row))
 ck(!anyNA(ar),"AE source-row join failed");ck(toupper(trimws(ev$term[ii]))==toupper(trimws(ae$AETERM[ar])),"AE term disagrees across exports")
 ev$st<-ev$et<-ev$temperature<-rep(NA_real_,nrow(ev));ev$st[ii]<-iso(ae$AESTDTC[ar]);ev$et[ii]<-iso(ae$AEENDTC[ar]);ev$related<-NA_character_;ev$related[ii]<-ae$AEREL[ar]
 sevnum<-function(x)match(toupper(as.character(x)),c("MILD","MODERATE","SEVERE","LIFE THREATENING"))
 ag<-grade[grade$source=="AE",]; aid<-match(key(ev$subject[ii],ae$source_ref[ar],ev$symptom_group[ii]),key(ag$subject,ag$source_ref,ag$concept))
 ev$sev<-NA_integer_;ev$sev[ii]<-sevnum(ag$severity[aid]); cg<-grade[grade$source=="COVID_FORMS",]
 cg$key<-key(cg$subject,cg$term,dy(cg$start_date));cg$val<-sevnum(cg$severity)
 # If same-day same-item grades disagree, leave unknown rather than borrowing another grade.
 gu<-unique(cg[,c("key","val")]);amb<-unique(gu$key[duplicated(gu$key)]);jj<-which(ev$source=="COVID_FORM"); gkey<-key(ev$subject[jj],ev$term[jj],ev$l[jj]);idx<-match(gkey,gu$key);ev$sev[jj]<-ifelse(gkey%in%amb,NA_integer_,gu$val[idx])
 # Numeric same-day diary temperatures corroborate an existing matched fever complaint only.
 # No new complaints from solicited data; no fever inferred from a label or subjective grade.
 tt<-grade[grade$source=="SOLICITED"&grade$concept=="FEVER"&grade$term=="Fever (C)",];tt$v<-suppressWarnings(as.numeric(tt$result))
 tt<-tt[is.finite(tt$v)&tt$v>=30&tt$v<=45,]; tm<-tapply(tt$v,key(tt$subject,dy(tt$start_date)),max)
 fi<-which(ev$symptom_group=="FEVER");ev$temperature[fi]<-unname(tm[key(ev$subject[fi],ev$l[fi])])
 # A single unambiguous temperature with explicit unit in the verbatim AE; exclude ranges, multiple values.
 for(i in fi[ev$source[fi]=="AE"]){s<-toupper(ev$term[i]);nums<-regmatches(s,gregexpr("[0-9]+([.][0-9]+)?",s))[[1]];if(length(nums)!=1)next
  n<-as.numeric(nums); val<-if(grepl("[0-9] *°? *F(\\b|AHRENHEIT)",s,perl=TRUE)&&n>=86&&n<=113)(n-32)*5/9 else if(grepl("[0-9] *°? *C(\\b|ELSIUS)",s,perl=TRUE)&&n>=30&&n<=45)n else NA_real_
  if(is.finite(val))ev$temperature[i]<-max(val,ev$temperature[i],na.rm=TRUE)}
 # Optional verified symptom-level source abstraction (e.g. VS or continuous form intervals).
 v<-optfile("trigger.clinical");ev$verified_interval<-FALSE
 if(!is.null(v)){ck(all(c("subject","source_id","start_iso","end_iso","temp_c","severity","source_ref")%in%names(v)),"clinical overlay schema")
  ck(!anyNA(v$source_ref)&!anyDuplicated(key(v$subject,v$source_id)),"clinical overlay duplicate or sourceless")
  ix<-match(key(ev$subject,ev$source_id),key(v$subject,v$source_id)); sel<-which(!is.na(ix))
  for(i in sel){t<-v[ix[i],];if(!is.na(t$temp_c)){ck(t$temp_c>=30&t$temp_c<=45,"Invalid temperature");ev$temperature[i]<-t$temp_c}
   if(!is.na(t$start_iso)&&!is.na(t$end_iso)){sd<-dy(t$start_iso);ed<-dy(t$end_iso);ck(sd<=ev$l[i]&ed>=ev$l[i],"Clinical interval does not contain source date")
    ev$l[i]<-sd;ev$h[i]<-min(ed,ep$fe[k[i]]);ev$st[i]<-iso(t$start_iso);ev$et[i]<-iso(t$end_iso);ev$verified_interval[i]<-TRUE}
   if(!is.na(t$severity))ev$sev[i]<-sevnum(t$severity)}}
 cuts<-dy(c("2020-06-26","2020-07-31","2020-08-20","2020-09-30","2020-11-11","2020-12-23","2021-02-10","2021-03-23"))
 ps<-optfile("trigger.site_versions")
 if(!is.null(ps)){ck(all(c("SITEID","amendment","effective_date","source_ref")%in%names(ps)),"Site version schema")
  ck(!anyNA(ps$source_ref)&!anyDuplicated(key(ps$SITEID,ps$amendment)),"Invalid site history")
  for(s in unique(ps$SITEID)){zz<-ps[ps$SITEID==s,];zz<-zz[order(zz$amendment),];ck(identical(as.integer(zz$amendment),1:8)&&all(is.finite(dy(zz$effective_date)))&&all(diff(dy(zz$effective_date))>=0),"Supply verified effective dates for amendments 1 through 8 at each site in the overlay")}}
 version<-function(date,site,lag=0){res<-findInterval(date,c(cuts[1],cuts[-1]+lag)); if(!is.null(ps))for(s in intersect(site,ps$SITEID)){
  ix<-which(site==s);z<-ps[ps$SITEID==s,];z<-z[order(z$amendment),];res[ix]<-findInterval(date[ix],dy(z$effective_date))}
  res}
 recent<-function(d,d1,d2){(is.finite(d1)&d>=d1&d<d1+7)|(is.finite(d2)&d>=d2&d<d2+7)}
 # Repeated assessments are not by themselves a continuously symptomatic interval.
 # Sensitivity only: four uninterrupted daily positives of the SAME item, no contradictory day.
 ev$repeat_day<-NA_integer_; kk<-key(ev$episode_id,ev$symptom_code)
 bad_by_key<-split(ev$l[ev$conflict],kk[ev$conflict]); for(ids in split(jj[!ev$conflict[jj]],kk[jj[!ev$conflict[jj]]])){ds<-sort(unique(ev$l[ids]));if(length(ds)<4)next
  ds<-setdiff(ds,bad_by_key[[kk[ids[1]]]]);if(!length(ds))next
  for(run in split(ds,cumsum(c(TRUE,diff(ds)!=1))))if(length(run)>=4){i<-ids[match(run[3],ev$l[ids])];ev$repeat_day[i]<-run[3]}}
 # In date-only AE intervals ending two dates later, elapsed time can straddle 48 h.
 # If boundary times exist, use the narrowest supported range, never inclusive-day counts.
 startlo<-ifelse(is.finite(ev$st),ev$st,ev$l*86400);starthi<-ifelse(is.finite(ev$st),ev$st,(ev$l+1)*86400-0.001)
 endlo<-ifelse(is.finite(ev$et),ev$et,ev$h*86400);endhi<-ifelse(is.finite(ev$et),ev$et,(ev$h+1)*86400-0.001)
 endlo<-pmin(endlo,(ep$fe[k]+1)*86400-.001);endhi<-pmin(endhi,(ep$fe[k]+1)*86400-.001)
 interval<-ev$verified_interval | (ev$source=="AE" & ev$date_quality=="KNOWN_INTERVAL")
 duration_yes<-interval & endlo-starthi>=48*3600;duration_no<-interval & endhi-startlo<48*3600
 classify<-function(repeated=FALSE,lag=0){
  va<-version(ev$l,ev$site,lag); group<-ev$symptom_group
  fever<-group=="FEVER";resp<-group%in%c("COUGH","DYSPNEA");upper<-group%in%c("NASAL","SORE_THROAT")
  req48<-va>=2 & !(resp|fever|group=="CHILLS");status<-rep("MANDATORY",nrow(ev));dd<-ev$l
  status[req48]<-"DURATION_INDETERMINATE";status[req48&duration_no]<-"DOES_NOT_MEET"
  yes<-req48&duration_yes;status[yes]<-"MANDATORY";dd[yes]<-ev$l[yes]+2
  if(repeated){id<-which(req48&is.finite(ev$repeat_day));status[id]<-"MANDATORY";dd[id]<-ev$repeat_day[id]}
  status[fever & !is.finite(ev$temperature)]<-"TEMPERATURE_INDETERMINATE"
  status[fever & is.finite(ev$temperature)&ev$temperature<38]<-"DOES_NOT_MEET"
  early<-va>=2 & recent(dd,ev$d1,ev$d2)
  # Upper-airway symptoms trigger evaluation in week 1. Do not equate evaluation
  # with a mandatory immediate swab for an upper-airway symptom of unknown duration.
  status[early & upper & status!="MANDATORY"]<-"EVALUATION_ONLY"
  overlap<-group%in%c("FEVER","CHILLS","FATIGUE","MYALGIA","HEADACHE","NAUSEA_VOMITING")
  latestdose<-pmax(ifelse(is.finite(ev$d1)&dd>=ev$d1,ev$d1,-Inf),ifelse(is.finite(ev$d2)&dd>=ev$d2,ev$d2,-Inf))
  disc<-which(early & overlap & ev$l>=latestdose & status=="MANDATORY");status[disc]<-"CLINICAL_JUDGMENT"
  # Continuous symptoms still present after the first week can qualify then.
  d<-pmax(ifelse(is.finite(ev$d1)&dd>=ev$d1,ev$d1,-Inf),ifelse(is.finite(ev$d2)&dd>=ev$d2,ev$d2,-Inf))+7
  regain<-disc[interval[disc]&ev$h[disc]>=d[disc]&!recent(d[disc],ev$d1[disc],ev$d2[disc])]
  status[regain]<-"MANDATORY";dd[regain]<-d[regain]
  status[va<1|ev$conflict|dd<ep$fs[k]|dd>ep$fe[k]]<-"INDETERMINATE_OTHER"
  data.frame(status=status,date=dd,version=va)}
 main<-classify(); alt<-classify(TRUE);lag14<-classify(FALSE,14)
 ev$rule<-main$status;ev$qual_day<-main$date;ev$version_at_symptom<-main$version
 ev$qual_repeat<-ifelse(alt$status=="MANDATORY",alt$date,NA_integer_);ev$qual_lag14<-ifelse(lag14$status=="MANDATORY",lag14$date,NA_integer_)
 # Classify recorded specimen; every MBSEQ used must cross-reference the same subject.
 j<-match(key(m$subject,m$MBSEQ),key(mb$subject,mb$MBSEQ));ck(!anyNA(j)&&!anyDuplicated(key(mb$subject,mb$MBSEQ)),"MB specimen linkage")
 m$spec<-ifelse(mb$MBSPEC[j]=="SWABBED MATERIAL","SWAB",ifelse(mb$MBSPEC[j]=="SALIVA","SALIVA","UNKNOWN"));m$spec[is.na(m$spec)]<-"UNKNOWN"
 m$day<-dy(m$test_date);m$clock<-iso(m$collection_text)
 m$valid<-m$diagnostic & !m$baseline & m$pcr_method & !m$method_conflict & !m$date_conflict & is.finite(m$day)
 # Dated specimen plus binary/indeterminate assay or explicit specimen material, excluding NOT DONE.
 m$sample<-m$valid & !m$not_done & (m$usable | m$spec!="UNKNOWN")
 refs<-optfile("trigger.reports");if(!is.null(refs)){ck(all(c("episode_id","report_iso","source_ref")%in%names(refs))&&
 !anyNA(refs$source_ref)&&!anyDuplicated(refs$episode_id)&&all(refs$episode_id%in%ep$episode_id),"report overlay schema");refs$time<-iso(refs$report_iso);ck(all(is.finite(refs$time)),"Supply full local report times")}
 form_file<-getOption("trigger.forms",file.path(root,"inputs","form_days.csv"))
 if(!file.exists(form_file))stop("Supply the checklist-day index (trigger.forms); see README.")
 forms<-readf(form_file)
 rawform<-optfile("trigger.raw_adsymp")
 if(!is.null(rawform)){ck(all(c("SUBJID","ADT","SRCDOM","PARAMCD")%in%names(rawform)),"Raw ADSYMP CSV must contain SUBJID, ADT, SRCDOM, PARAMCD.")
  core<-c("BODYACHE","CHILLS","COUGH","DIARRHEA","DIBREATH","FATIGUE","HEADACHE","MYALGIA","NASALCON","NAUSEA","NLSMELL","NLTASTE","RUNNOSE","SHBREATH","STHROAT","VOMIT")
  rawform<-rawform[rawform$SRCDOM=="FA"&rawform$PARAMCD%in%core & !is.na(rawform$ADT),]
  ck(nrow(rawform)>0 && all(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}",rawform$ADT)),"Require ISO date strings in raw ADSYMP CSV")
  forms<-unique(data.frame(subject=rawform$SUBJID,date=substr(rawform$ADT,1,10),verified_checklist_rows=1L,untyped_rows=0L))}
 ck(all(c("subject","date","verified_checklist_rows","untyped_rows")%in%names(forms)),"Form-day index missing fields")
 dates<-dy(forms$date);ck(all(is.finite(dates)),"Invalid form-index dates")
 knownform<-split(dates[forms$verified_checklist_rows>0],forms$subject[forms$verified_checklist_rows>0])
 possibleform<-split(dates[forms$untyped_rows>0],forms$subject[forms$untyped_rows>0])
 sm<-split(seq_len(nrow(m)),m$subject);se<-split(seq_len(nrow(ev)),ev$episode_id);vis<-split(dy(visit$visit_date[visit$independent]),visit$subject[visit$independent])
 # Ownership: for two episodes in one person never attribute a later episode's sample to the earlier one.
 nxt<-rep(Inf,nrow(ep));for(ix in split(seq_len(nrow(ep)),ep$subject)){ix<-ix[order(ep$s[ix])];if(length(ix)>1)nxt[head(ix,-1)]<-ep$s[tail(ix,-1)]}
 mn<-function(x){x<-x[is.finite(x)];if(length(x))min(x)else NA_real_}
 mx<-function(x){x<-x[is.finite(x)];if(length(x))max(x)else NA_real_}
 resclass<-function(z){if(!nrow(z))return("UNTESTED");ok<-!z$result_conflict;if(any(ok&z$result=="POSITIVE"))"POSITIVE" else if(any(ok&z$result=="NEGATIVE"))"NEGATIVE" else "INDETERMINATE"}
 # An early upper-airway symptom requires evaluation even if an explicit swab rule
 # is first met later in the same episode. The investigation clock is kept separate.
 upper_early<-ev$symptom_group%in%c("NASAL","SORE_THROAT") & main$version>=2 & recent(ev$l,ev$d1,ev$d2) & !ev$conflict & ev$l>=ep$fs[k]&ev$l<=ep$fe[k]
 early_eval_ids<-unique(ev$episode_id[upper_early])
 one<-function(i,evalmode=FALSE){
  e<-ep[i,];v<-ev[se[[e$episode_id]],,drop=FALSE];if(evalmode){
    u<-upper_early[se[[e$episode_id]]];v$rule[u]<-"EVALUATION_ONLY";v$qual_day[u]<-v$l[u]}
  req<-v$rule%in%if(evalmode)c("MANDATORY","EVALUATION_ONLY")else"MANDATORY";tt<-mn(v$qual_day[req])
  judge<-any(v$rule=="CLINICAL_JUDGMENT");eo<-any(v$rule=="EVALUATION_ONLY");unc<-any(grepl("INDETERMINATE",v$rule))
  cl<-if(is.finite(tt)&&evalmode)"REQUIRED_EVALUATION" else if(is.finite(tt))"MANDATORY" else if(eo)"EVALUATION_ONLY" else if(judge)"CLINICAL_JUDGMENT" else if(unc)"INDETERMINATE" else "DOES_NOT_MEET"
  short<-function(x)paste(sort(unique(x)),collapse="+")
  # Use a separate threshold for evaluation-only and discretionary outputs.
  if(!is.finite(tt))tt<-if(cl%in%c("EVALUATION_ONLY","CLINICAL_JUDGMENT"))mn(v$qual_day[v$rule==cl])else NA_real_
  ad<-if(is.finite(tt))tt else e$s;qual<-v[v$qual_day==ad & if(evalmode)req else v$rule==cl,,drop=FALSE]
  obs<-v[v$l<=ad,,drop=FALSE];active<-obs[obs$h>=ad,,drop=FALSE];sev<-mx(active$sev)
  gs<-qual$symptom_group;fam<-if("DYSPNEA"%in%gs)"DYSPNEA" else if("COUGH"%in%gs)"COUGH" else if("CHILLS"%in%gs)"CHILLS" else if("FEVER"%in%gs)"FEVER" else if(length(gs))short(gs) else "UNKNOWN"
  candidate_dates<-function(x)x[x>=e$s & x<=min(e$fe,nxt[i]-1, max(e$last,ad+7))]
  formdates<-candidate_dates(c(v$l[v$source_domain=="FA" & v$source=="COVID_FORM"],knownform[[e$subject]]))
  form<-mn(formdates);fd<-form-ad; posform<-mn(candidate_dates(possibleform[[e$subject]]))
  form_possible<-is.finite(posform)&&posform<=ad+3
  msub<-m[sm[[e$subject]],,drop=FALSE];owned<-msub$day>=e$s & msub$day<=min(e$fe,nxt[i]-1)&is.finite(msub$day)
  msub<-msub[owned & msub$valid,,drop=FALSE];u<-msub[msub$sample,,drop=FALSE]
  z<-u[u$day>=ad & u$day<=ad+7,,drop=FALSE];z3<-z[z$day<=ad+3,,drop=FALSE]
  sw<-z[z$spec=="SWAB",,drop=FALSE];delta<-mn(sw$day)-ad
  # Site awareness clock, if independently supplied; it is never filled from symptom onset/ADT.
  outcome<-"NO_VERIFIED_REPORT";refday<-NA_real_
  if(!is.null(refs)){at<-refs$time[match(e$episode_id,refs$episode_id)]
   if(is.finite(at)&&cl=="MANDATORY"){rd<-floor(at/86400);ck(rd>=ad & rd<=e$fe && rd<nxt[i],"Report outside episode/trigger range")
    refday<-rd; tend<-at+72*3600;swall<-u[u$spec=="SWAB",,drop=FALSE]
    low<-ifelse(is.finite(swall$clock),swall$clock,swall$day*86400);high<-ifelse(is.finite(swall$clock),swall$clock,(swall$day+1)*86400-.001)
    outcome<-tr_clock(low,high,at,(min(e$fe,nxt[i]-1)+1)*86400)
   }}
  # Match symptom duration as observed by the trigger. Future AE end and full
  # illness duration are never covariates. Count other supported symptom groups by that date.
  age<-if(nrow(qual))max(ad-qual$l,na.rm=TRUE)else 0
  ageband<-if(age==0)"same date" else if(age==1)"1 day" else if(age<=3)"2-3 days" else "4+ days"
  pattern<-if(length(unique(obs$symptom_group))<=1)"one"else if(length(unique(obs$symptom_group))==2)"two"else "three+"
  dur<-ad-e$s;band<-if(dur==0)"onset" else if(dur==1)"1 day" else if(dur<=3)"2-3 days" else "4+ days"
  dist<-ad-pmax(ifelse(is.finite(e$d1)&ad>=e$d1,e$d1,-Inf),ifelse(is.finite(e$d2)&ad>=e$d2,e$d2,-Inf))
  pd<-if(!is.finite(dist))"pre-dose" else if(dist<=6)"0-6" else if(dist<=13)"7-13" else if(dist<=27)"14-27" else "28+"
  indvisit<-vis[[e$subject]];indvisit<-indvisit[indvisit>=e$s&indvisit<=min(ad+3,e$fe,nxt[i]-1)]
  # A form already present at qualification counts as documented by the deadline.
  # Conditional swab requires collection on/after first form date and by trigger day +3.
  data.frame(episode_id=e$episode_id,subject=e$subject,site=e$SITEID,arm=e$arm,onset=dayfmt(e$s),classification=cl,
    onset_source=e$onset_source,trigger_date=dayfmt(tt),trigger_amendment=if(is.finite(tt))version(tt,e$SITEID)else NA,
    version_basis=if(!is.null(ps)&&e$SITEID%in%ps$SITEID)"SITE_LOG"else"SPONSOR_DATE_PROXY",
    trigger_family=fam,trigger_groups=short(gs),trigger_source=short(qual$source),severity=if(is.finite(sev))c("MILD","MODERATE","SEVERE","LIFE THREATENING")[sev]else"UNKNOWN",
    elapsed_onset_band=band,trigger_duration_band=ageband,recorded_group_count_band=pattern,dose_phase=if(!is.finite(e$d1)||ad<e$d1)"pre-dose" else if(is.finite(e$d2)&&ad>=e$d2)"dose2"else"dose1",
    days_since_dose_band=pd,form_first=dayfmt(form),form_already_before=is.finite(fd)&&fd<0,
    form_by_d3=is.finite(fd)&&fd<=3,form_upper_d3=(is.finite(fd)&&fd<=3)||form_possible,form_by_d7=is.finite(fd)&&fd<=7,form_in_episode=is.finite(form),independent_visit_by_d3=length(indvisit)>0,
    form_then_swab_by_d3=is.finite(form)&&any(sw$day>=form &sw$day<=ad+3),
    order_proven_swab_after_form_d3=is.finite(form)&&any(sw$day>form &sw$day<=ad+3),
    swab_first_delta=delta,swab_d0_3=is.finite(delta)&&delta<=3,swab_d0_2=is.finite(delta)&&delta<=2,
    swab_definitely_after_proxy_within72=any((sw$day-ad)%in%c(1,2)), # relative to the proxy date, not site awareness
    swab_boundary=if(!is.finite(delta)||delta>3)"NONE_BY_D3"else c("D0_ORDER_UNKNOWN","D1","D2","D3_BOUNDARY")[delta+1],
    any_sample_d0_3=nrow(z3)>0,saliva_d0_3=any(z3$spec=="SALIVA"),unknown_specimen_d0_3=any(z3$spec=="UNKNOWN"),
    swab_d0_7=nrow(sw)>0,any_sample_d0_7=nrow(z)>0,pcr_result_d0_3=resclass(z3),pcr_result_d0_7=resclass(z),
    positive_and_negative_d7=any(z$result=="POSITIVE")&&any(z$result=="NEGATIVE"),not_done_d7=any(msub$not_done & msub$day>=ad & msub$day<=ad+7),
    sampled_before_trigger=any(u$day<ad),any_sample_whole_observed_illness=any(u$day>=e$s & u$day<=max(e$last,ad+7)),first_swab_date=dayfmt(mn(sw$day)),verified_report_day=dayfmt(refday),exact72=outcome,
    complete3=is.finite(tt)&&ad+3<=min(e$fe,nxt[i]-1),complete7=is.finite(tt)&&ad+7<=min(e$fe,nxt[i]-1),
    orig_eligible=e$eligible,old_visit_plus_pcr=e$tested,pp14=isTRUE(e$pp)&&is.finite(e$d2)&&e$s>=e$d2+14&&ad<=e$effend,
    incident=!e$preexisting && e$onset_in_display,ae_independent_at_trigger=any(qual$source=="AE")&&(!is.finite(form)||form>=ad),
    possible_earlier_trigger=any(grepl("INDETERMINATE",v$rule)&v$l<ad),unknown_temp=any(v$rule=="TEMPERATURE_INDETERMINATE"),unknown_duration=any(v$rule=="DURATION_INDETERMINATE"),
    alternative_repeat_date=dayfmt(mn(v$qual_repeat)),alternative_lag14_date=dayfmt(mn(v$qual_lag14)),judgment_evidence=judge,
    respiratory_onset=any(v$symptom_group%in%c("COUGH","DYSPNEA")&v$l==e$s),stringsAsFactors=FALSE)
 }
 message("Classifying episodes and linking diagnostic specimens...")
 epout<-do.call(rbind,lapply(seq_len(nrow(ep)),one))
 evaluation<-epout; altids<-which(ep$episode_id%in%early_eval_ids)
 if(length(altids))evaluation[altids,]<-do.call(rbind,lapply(altids,function(i)one(i,TRUE)))
 evaluation$classification[evaluation$classification=="MANDATORY"]<-"REQUIRED_EVALUATION"
 evalset<-subset(evaluation,incident & classification=="REQUIRED_EVALUATION"&complete3)
 save(evalset,"restricted/required_evaluation_episodes.csv")
 ck(!anyDuplicated(epout$episode_id)&&nrow(epout)==nrow(ep),"Episode inventory mismatch")
 ev$qual_day<-dayfmt(ev$qual_day);ev$qual_repeat<-dayfmt(ev$qual_repeat);ev$qual_lag14<-dayfmt(ev$qual_lag14)
 save(ev[c("episode_id","subject","source","source_id","source_sequence","source_domain","symptom_code","symptom_group","date_quality","lo","hi","term","sev","related","temperature","rule","qual_day","qual_repeat","qual_lag14","version_at_symptom")],"restricted/symptom_decisions.csv")
 save(epout,"restricted/episodes.csv")
 flow<-as.data.frame(with(epout[epout$incident,],table(arm,classification)));save(flow,"eligibility_flow.csv")
 mainset<-subset(epout,incident & classification=="MANDATORY"&complete3)
 seven<-subset(mainset,complete7 & orig_eligible)
 measures<-list(Required_investigation=evalset,All=mainset,PP_D2_plus14=subset(mainset,pp14),Onset_cough_or_dyspnea=subset(mainset,respiratory_onset),
  Trigger_documented_in_AE_no_earlier_form=subset(mainset,ae_independent_at_trigger),Common_complete_7d=seven,No_linked_positive_7d=subset(seven,pcr_result_d0_7!="POSITIVE"),
  First_week_judgment=subset(epout,incident&classification=="CLINICAL_JUDGMENT"&complete3),
  First_week_evaluation_only=subset(epout,incident&classification=="EVALUATION_ONLY"&complete3))
 sums<-function(x,label){
  b<-lapply(c("mRNA-1273","Placebo"),function(arm){d<-x[x$arm==arm,];form<-d$form_by_d3
    c(n=nrow(d),participants=length(unique(d$subject)),form=sum(form),possible_form_upper=sum(d$form_upper_d3),swab=sum(d$swab_d0_3),
      swab_day1_2=sum(d$swab_definitely_after_proxy_within72),d0_2=sum(d$swab_d0_2),sample=sum(d$any_sample_d0_3),
      saliva=sum(d$saliva_d0_3),unknown_specimen=sum(d$unknown_specimen_d0_3),no_form_swab=sum(!form&d$swab_d0_3),
      form_then_swab=sum(d$form_then_swab_by_d3 & form),strict_after_form=sum(d$order_proven_swab_after_form_d3 & form),
      form_before=sum(d$form_already_before),verified_reports=sum(d$exact72!="NO_VERIFIED_REPORT"),
      seven_n=sum(d$complete7&d$orig_eligible),swab7=sum(d$swab_d0_7&d$complete7&d$orig_eligible),
      old7=sum(d$old_visit_plus_pcr & d$complete7 & d$orig_eligible))})
  data.frame(scope=label,arm=c("mRNA-1273","Placebo"),do.call(rbind,b),row.names=NULL)}
 summary<-do.call(rbind,Map(sums,measures,names(measures)));save(summary,"arm_summary.csv")
 results<-do.call(rbind,lapply(names(measures),function(nm){d<-measures[[nm]]
   do.call(rbind,lapply(c(3,7),function(t){if(t==7)d<-d[d$complete7,]
   r<-d[[paste0("pcr_result_d0_",t)]];data.frame(scope=nm,days=t,as.data.frame(table(arm=d$arm,result=factor(r,levels=c("POSITIVE","NEGATIVE","INDETERMINATE","UNTESTED"))))) }))}));save(results,"pcr_outcomes.csv")
 save(as.data.frame(with(mainset,table(arm,exact72))),"verified_report_clock.csv")
 sa<-subset(epout,incident)
 sensitivity<-do.call(rbind,lapply(c("trigger_date","alternative_repeat_date","alternative_lag14_date"),function(field){
  ok<-if(field=="trigger_date")sa$classification=="MANDATORY" else !is.na(sa[[field]]);data.frame(rule=field,as.data.frame(table(sa$arm[ok])))}))
 save(sensitivity,"trigger_sensitivities.csv")
 valid_sw<-split(m$day[m$sample&m$spec=="SWAB"],m$subject[m$sample&m$spec=="SWAB"])
 altcov<-do.call(rbind,lapply(c("trigger_date","alternative_repeat_date","alternative_lag14_date"),function(field){
  anchors<-dy(epout[[field]]);eligible<-epout$incident & is.finite(anchors)&anchors+3<=pmin(ep$fe,nxt-1)
  if(field=="trigger_date")eligible<-eligible&epout$classification=="MANDATORY"
  ix<-which(eligible)
  F<-vapply(ix,function(i){ds<-knownform[[ep$subject[i]]];any(ds>=ep$s[i]&ds<=anchors[i]+3)},logical(1))
  W<-vapply(ix,function(i){ds<-valid_sw[[ep$subject[i]]];any(ds>=anchors[i]&ds<=anchors[i]+3)},logical(1))
  z<-data.frame(arm=ep$arm[ix],n=1L,form=as.integer(F),swab=as.integer(W))
  data.frame(rule=field,aggregate(z[c("n","form","swab")],list(arm=z$arm),sum))}))
 save(altcov,"sensitivity_coverage.csv")
 xx<-altcov[altcov$rule=="trigger_date",];yy<-summary[summary$scope=="All",];j<-match(xx$arm,yy$arm);ck(all(xx$n==yy$n[j]&xx$form==yy$form[j]&xx$swab==yy$swab[j]),"Independent anchor summary disagrees")

 # Comparison over within-site common support. No outcome or source variables in the strata.
 # Overlap weights nV*nP/(nV+nP); same weights in both arms; entire sites resampled for CI.
 standard<-function(z,stage,period,scope,complete_grade=FALSE){
  if(stage=="swab_given_form")z<-z[z$form_by_d3,]
  if(complete_grade)z<-z[z$severity!="UNKNOWN",]
  t<-as.Date(z$trigger_date);p<-if(period=="month")format(t,"%Y-%m")else as.character(t-(as.integer(format(t,"%u"))-1))
  z$outcome<-if(stage=="form")z$form_by_d3 else if(stage=="swab_given_form")z$form_then_swab_by_d3 else z$swab_d0_3
  z$stratum<-key(z$site,p,z$trigger_family,z$severity,z$elapsed_onset_band,z$trigger_duration_band,z$recorded_group_count_band,z$dose_phase,z$days_since_dose_band)
  a1<-z[z$arm=="mRNA-1273",];a0<-z[z$arm=="Placebo",]
  common<-intersect(unique(a1$stratum),unique(a0$stratum));c1<-a1[a1$stratum%in%common,];c0<-a0[a0$stratum%in%common,]
  est<-lo<-hi<-r1<-r0<-NA_real_;sites<-0L
  if(length(common)){f<-function(d) {ix<-split(d,d$stratum);do.call(rbind,lapply(ix,function(x)data.frame(st=x$stratum[1],site=x$site[1],n=nrow(x),r=mean(x$outcome))))}
   x<-merge(f(c1),f(c0),by=c("st","site"));x$w<-x$n.x*x$n.y/(x$n.x+x$n.y);r1<-weighted.mean(x$r.x,x$w);r0<-weighted.mean(x$r.y,x$w);est<-100*(r1-r0)
   g<-aggregate(cbind(num=x$w*(x$r.x-x$r.y),den=x$w),list(site=x$site),sum);sites<-nrow(g)
   if(sites>=3){ww<-rmultinom(2000,sites,rep(1/sites,sites));q<-quantile(100*as.vector(crossprod(g$num,ww))/as.vector(crossprod(g$den,ww)),c(.025,.975),names=FALSE);lo<-q[1];hi<-q[2]}}
  data.frame(scope=scope,stage=stage,period=period,known_grade_only=complete_grade,vaccine_all=nrow(a1),placebo_all=nrow(a0),
   vaccine_supported=nrow(c1),placebo_supported=nrow(c0),common_strata=length(common),sites=sites,vaccine_standardized_pct=100*r1,
   placebo_standardized_pct=100*r0,gap_pp=est,lower95=lo,upper95=hi)}
 set.seed(30172);sts<-list()
 for(nm in c("All","PP_D2_plus14","Trigger_documented_in_AE_no_earlier_form","Onset_cough_or_dyspnea","Required_investigation","No_linked_positive_7d"))for(stage in c("form","swab_given_form","swab_all"))for(period in c("month","week"))for(g in if(period=="month")c(FALSE,TRUE) else FALSE)
 sts[[length(sts)+1L]]<-standard(measures[[nm]],stage,period,nm,g)
 std<-do.call(rbind,sts);save(std,"matched_stage_comparisons.csv")
 # Reviewer-facing packet: all mandatory cough/dyspnea triggers without a documented
 # diagnostic specimen in 7 days + randomly interleaved tested controls. No labels, IDs,
 # source names, arm, site, absolute date, dose number, outcomes, or crosswalk in the packet.
 pool<-subset(seven,grepl("COUGH|DYSPNEA",trigger_groups));un<-which(!pool$any_sample_d0_7)
 co<-which(pool$any_sample_d0_7);set.seed(74301)
 control<-if(length(co))co[sample.int(length(co),min(length(un),length(co)))] else integer()
 ids<-c(un,control);if(length(ids))ids<-ids[sample.int(length(ids))]
 pack<-pool[ids,,drop=FALSE];pack$review_id<-sprintf("R%05d",seq_len(nrow(pack)))
 save(pack[,c("review_id","episode_id","subject","arm","site","any_sample_d0_7","sampled_before_trigger","any_sample_whole_observed_illness")],"restricted/reviewer_key.csv")
 rdir<-file.path(out,"reviewers_only");dir.create(rdir,showWarnings=FALSE)
 blank<-data.frame(review_id=pack$review_id,days_since_dose_band=pack$days_since_dose_band,amendment=pack$trigger_amendment,
   evaluation_warranted="",swab_warranted="",earliest_qualifying_relative_day="",confidence="",reason_or_records_needed="")
 write.csv(blank,file.path(rdir,"independent_review_sheet.csv"),row.names=FALSE,na="")
 ex<-match(ev$episode_id,pack$episode_id);ok<-which(!is.na(ex))
 # Only evidence through the trigger enters this packet; later investigation and PCR
 # language in AE terms could otherwise disclose outcome. Withhold whole text if flagged.
 dt<-dy(pack$trigger_date[ex[ok]]);ok<-ok[ev$l[ok]<=dt];ii<-ex[ok]
 terms<-as.character(ev$term[ok]);mask<-grepl("VACC|PLACEBO|M.?RNA|1273|PCR|POSITIVE|NEGATIVE|NASAL SWAB|TESTED|COVID|SARS|DOSE|INJECT|US[0-9]{3}|UNBLIND|HOSPITALIZED|LABORATORY|\\b(19|20)[0-9]{2}\\b",terms,ignore.case=TRUE)
 terms[mask]<-"[description requires masked source abstraction]"
 details<-data.frame(review_id=pack$review_id[ii],day_from_trigger=ev$l[ok]-dy(pack$trigger_date[ii]),
   symptom=ev$symptom_group[ok],source_description=terms,grade=c("Mild","Moderate","Severe","Life-threatening")[ev$sev[ok]],temperature_C=ev$temperature[ok],
   continuous_duration_to_trigger_evidenced=ifelse(interval[ok] & duration_yes[ok] & ev$l[ok]+2<=dy(pack$trigger_date[ii]),"48h supported","not established"))
 details<-unique(details);details<-details[order(details$review_id,details$day_from_trigger,details$symptom),]
 write.csv(details,file.path(rdir,"clinical_excerpts.csv"),row.names=FALSE,na="")
 writeLines(c("Give ONLY this reviewers_only folder to reviewers. Each reviewer completes a separate copy; resolve disagreements after locking both.",
 "These excerpts are coded form terms or verbatim AE descriptions when available, not complete narratives. Day 0 is an algorithmic candidate date, not verified notification.",
 "Review the indication for investigation and for swabbing separately. Use YES / NO / DISCRETION / INDETERMINATE. Record the earliest defensible day relative to day 0, rationale and source still needed.",
 "Limited evidence is a valid finding. Request independently redacted CRFs, temperature records, contact logs and clinical notes through the custodian before final adjudication.",
 "Reviewers must not receive the analytical results or restricted key before submitting judgments."),file.path(rdir,"READ_FIRST.txt"))
 save(data.frame(amendment=seq_along(cuts),sponsor_date=dayfmt(cuts),rule_family=c("any listed symptom",rep("48h rule, immediate exceptions and week-one judgment",7))),"protocol_calendar.csv")
 provenance<-data.frame(role=c("heat","symptom_burden","PCR specimen metadata","ADSL","checklist index"),
 source=c(heat,burden,lab,adsl,form_file));save(provenance,"sources.csv")
 tr_add_outputs(out)
 tr_report(out)
 invisible(list(summary=summary,models=std,flow=flow))
}
tr_clock<-function(earliest,latest,report,observed_until){
 deadline<-report+72*3600
 if(any(earliest>=report & latest<=deadline))"WITHIN_72H"
 else if(any(latest>=report & earliest<=deadline))"ORDER_OR_BOUNDARY_UNCERTAIN"
 else if(observed_until<deadline)"CENSORED" else "NO_DOCUMENTED_SWAB_WITHIN_72H"
}
tr_add_outputs<-function(out){
 read<-function(n)read.csv(file.path(out,"restricted",n),stringsAsFactors=FALSE)
 e<-read("episodes.csv");d<-e[e$incident & e$complete3 & e$classification=="MANDATORY",];q<-read("required_evaluation_episodes.csv")
 z<-rbind(data.frame(site=q$site,arm=q$arm,stage="Checklist by D3, required evaluation",yes=q$form_by_d3),
 data.frame(site=d$site,arm=d$arm,stage="Swab D0-3, mandatory",yes=d$swab_d0_3))
 f<-d[d$form_by_d3,];z<-rbind(z,data.frame(site=f$site,arm=f$arm,stage="Swab on/after checklist, mandatory",yes=f$form_then_swab_by_d3))
 x<-aggregate(cbind(events=as.integer(z$yes),episodes=rep(1L,nrow(z))),z[c("site","arm","stage")],sum)
 x$pct<-100*x$events/x$episodes;write.csv(x,file.path(out,"site_stage_counts.csv"),row.names=FALSE)
 # Optional dated eCRF trail, to distinguish actual form opening from assessment dates.
 logfile<-getOption("trigger.form_log","")
 if(nzchar(logfile)){
  l<-read.csv(logfile,stringsAsFactors=FALSE);stopifnot(all(c("episode_id","opened_iso","source_ref")%in%names(l)),
   !anyDuplicated(l$episode_id),all(l$episode_id%in%e$episode_id),all(!is.na(l$source_ref)&nzchar(l$source_ref)),
   all(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}",l$opened_iso)))
  j<-match(l$episode_id,q$episode_id); l$proxy_evaluation_date<-q$trigger_date[j];l$recorded_checklist_date<-q$form_first[j]
  l$open_minus_trigger_days<-as.integer(as.Date(substr(l$opened_iso,1,10))-as.Date(l$proxy_evaluation_date))
  write.csv(l,file.path(out,"restricted","verified_form_opening_audit.csv"),row.names=FALSE,na="")
 }
}

tr_report<-function(out){
 get<-function(f)read.csv(file.path(out,f),stringsAsFactors=FALSE,check.names=FALSE)
 ar<-get("arm_summary.csv");st<-get("matched_stage_comparisons.csv");fl<-get("eligibility_flow.csv");rr<-get("pcr_outcomes.csv");sens<-get("sensitivity_coverage.csv")
 e<-get("restricted/episodes.csv");review<-get("restricted/reviewer_key.csv")
 esc<-function(x){x<-gsub("&","&amp;",as.character(x),fixed=TRUE);x<-gsub("<","&lt;",x,fixed=TRUE);gsub(">","&gt;",x,fixed=TRUE)}
 tablehtml<-function(d){paste0("<table><thead><tr>",paste0("<th>",esc(names(d)),"</th>",collapse=""),"</tr></thead><tbody>",
 paste(apply(d,1,function(x)paste0("<tr>",paste0("<td>",esc(x),"</td>",collapse=""),"</tr>")),collapse=""),"</tbody></table>")}
 rate<-function(n,d)ifelse(d>0,sprintf("%s/%s (%.1f%%)",n,d,100*n/d),"—")
 pick<-function(label){z<-ar[ar$scope==label,];z[match(c("mRNA-1273","Placebo"),z$arm),]}
 row<-function(label,z,num,den="n"){data.frame(Measure=label,'mRNA-1273'=rate(z[[num]][1],z[[den]][1]),Placebo=rate(z[[num]][2],z[[den]][2]),check.names=FALSE)}
 all<-pick("All");req<-pick("Required_investigation");seven<-pick("Common_complete_7d");ax<-pick("Trigger_documented_in_AE_no_earlier_form")
 primary<-rbind(row("Checklist documented by D3 | all investigation triggers",req,"form"),
 row("Checklist documented by D3 | swab-eligible triggers",all,"form"),
 row("Checklist already documented before trigger | swab-eligible",all,"form_before"),
 row("Swab from trigger date to D3 | all swab-eligible",all,"swab"),
 row("Swab on/after first checklist date, by D3 | checklist documented",all,"form_then_swab","form"),
 row("Swab strictly after first checklist date, by D3 | checklist documented",all,"strict_after_form","form"),
 row("Any diagnostic specimen, including saliva, D0–3",all,"sample"),row("Documented swab D0–2",all,"d0_2"),
 row("At least one swab on D1–2",all,"swab_day1_2"),
 row("Swab without documented checklist by D3 | all swab-eligible",all,"no_form_swab"),
 row("Swab D0–7 from qualifying date | common 7-day cohort",seven,"swab7"),
 row("Original visit + PCR within 7 days of onset | same cohort",seven,"old7"))
 mtable<-function(nms,per="month"){a<-st[st$scope%in%nms&st$period==per&!st$known_grade_only &
      ((st$scope=="Required_investigation"&st$stage=="form")|(st$scope=="All"&st$stage!="form")|!st$scope%in%c("All","Required_investigation")),]
 data.frame(Population=a$scope,Decision=a$stage,'Matched episodes V / P'=paste(a$vaccine_supported,a$placebo_supported,sep=" / "),
  Sites=a$sites,'mRNA-1273'=sprintf("%.1f%%",a$vaccine_standardized_pct),'Placebo'=sprintf("%.1f%%",a$placebo_standardized_pct),
  'Difference, pp (95% interval)'=sprintf("%+.1f (%.1f to %.1f)",a$gap_pp,a$lower95,a$upper95),check.names=FALSE)}
 pp<-pick("PP_D2_plus14"); ppR<-rbind(row("Checklist documented by D3",pp,"form"),row("Swab D0–3",pp,"swab"),row("Swab on/after checklist date, by D3 | checklist documented",pp,"form_then_swab","form"))
 flow<-reshape(fl,idvar="classification",timevar="arm",direction="wide");names(flow)<-sub("Freq[.]","",names(flow))
 outcomes<-reshape(rr[rr$scope=="All",c("days","result","arm","Freq")],idvar=c("days","result"),timevar="arm",direction="wide");names(outcomes)<-sub("Freq[.]","",names(outcomes))
 srt<-data.frame(Rule=sens$rule,Arm=sens$arm,Eligible=sens$n,Checklist=rate(sens$form,sens$n),Swab=rate(sens$swab,sens$n))
 i<-e$incident & e$complete3 & e$classification=="MANDATORY"
 unk<-tapply(e$classification[e$incident]=="INDETERMINATE",e$arm[e$incident],sum)
 ro<-table(review$any_sample_d0_7)
 blend<-function(scope,stage,period="month"){x<-st[st$scope==scope&st$stage==stage&st$period==period&!st$known_grade_only,]
 sprintf("%+.1f percentage points (site-bootstrap 95%% interval %.1f to %.1f; %d vaccine and %d placebo episodes in comparable strata)",x$gap_pp,x$lower95,x$upper95,x$vaccine_supported,x$placebo_supported)}
 headline<-paste0("An investigation checklist was documented by the third calendar day for ",rate(req$form[1],req$n[1])," vaccine-arm and ",rate(req$form[2],req$n[2])," placebo-arm episodes that warranted evaluation. Among episodes that warranted a swab and had a checklist, a swab was collected on or after the first checklist date and by the third day in ",rate(all$form_then_swab[1],all$form[1])," versus ",rate(all$form_then_swab[2],all$form[2]),".")
 html<-c('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Moderna P301 | Protocol-triggered illness audit</title><style>body{max-width:64rem;margin:2rem auto;padding:0 1rem;color:#202e38;font:17px/1.5 system-ui,sans-serif}h1,h2{line-height:1.18}h2{margin-top:1.8em}p{max-width:84ch}table{border-collapse:collapse;width:100%;font-size:.88em;margin:1em 0}th,td{padding:.42em;text-align:left;border-bottom:1px solid #c9d3db}th{background:#eff4f6}td:first-child{overflow-wrap:anywhere}.lead{font-family:Georgia,serif;font-size:1.12em}@media print{thead{display:table-header-group}tr{break-inside:avoid}} </style><body>',
 "<h1>Which protocol-triggered illnesses were investigated and swabbed?</h1>",
 paste0('<p class="lead">',headline,"</p>"),paste0("<p>Within comparable site-month strata, the vaccine-minus-placebo difference in checklist documentation was ",blend("Required_investigation","form"),". For swabbing after checklist documentation it was ",blend("All","swab_given_form"),". These comparisons concern the episodes represented in both arms.</p>"),
 "<p><strong>Clock:</strong> D0 is the first date on which the available clinical evidence supports a qualifying symptom, not necessarily the time staff were notified. D0–3 includes a same-day sample of uncertain order and third-day samples potentially beyond 72 hours. D1–2 is unambiguously within 72 hours of the proxy date for any time of day, subject to local-clock precision. None of these measures proves when the site knew.</p>",
 paste0("<p>Verified report-to-specimen clock statuses: ",paste(names(table(e$exact72[i])),as.integer(table(e$exact72[i])),collapse="; "),". Without a verified qualifying-report timestamp an episode cannot be classified as meeting or missing the exact deadline.</p>"),
 "<h2>From symptoms to a form, then to a specimen</h2>",tablehtml(primary),
 paste0("<p>A dated local- or illness-labelled diagnostic PCR with swabbed material supports a swab; saliva is counted separately in the any-specimen analysis. The physical location is not inferred. The original comparator also requires a visit and permits a PCR up to seven days from <em>initial onset</em>; the new seven-day comparison starts at <em>qualification</em>. Neither comparator is a literal extension of the other. Possible untyped checklist rows would raise the stage-one counts to at most ",paste(req$possible_form_upper,collapse=" versus "),".</p>"),
 "<h2>Review of the first decision independently of a positive symptom checklist</h2>",
 paste0("<p>When a qualifying symptom was documented in an AE and there was no earlier verified checklist, subsequent or same-day checklist documentation was ",rate(ax$form[1],ax$n[1])," versus ",rate(ax$form[2],ax$n[2]),". Same-day records do not establish which came first. This subset matters because a symptom first found on a checklist already provides evidence for the outcome under study.</p>"),
 "<h2>Comparable episodes</h2><p>Exact strata use site, calendar month (or week), trigger symptom family, recorded severity, time since episode onset, observed trigger duration, number of symptom groups observed by the trigger, dose number and time since dose. Rates use the same overlap weights in both arms. The intervals resample sites; they describe comparisons within common-support strata and do not address unmeasured severity or unreported complaints. Observed AE grade can summarize an entire event. The CSV includes denominators before matching, known-grade-only checks, AE-origin episodes and other sensitivity analyses.</p>",
 tablehtml(mtable(c("Required_investigation","All","Trigger_documented_in_AE_no_earlier_form"))),
 "<h3>Exact site and calendar week</h3>",tablehtml(mtable(c("Required_investigation","All"),"week")),
 "<h2>What did the diagnostic specimens show?</h2>",tablehtml(outcomes),
 "<p>Results refer to any documented diagnostic specimen after the trigger (including saliva), by episode, not to adjudicated efficacy cases. POSITIVE takes precedence, then a nonconflicting NEGATIVE; an unresolved or conflicting result without a valid binary result is INDETERMINATE. UNTESTED means no documented qualifying specimen in that window. NOT DONE rows and samples before the trigger remain visible in the episode file.</p>",
 "<h2>Blinded per-protocol sensitivity: onset at least 14 days after dose 2</h2>",tablehtml(ppR),"<p>This uses baseline-negative, per-protocol participants within their follow-up and efficacy cutoff. It does not reconstruct first-event censoring or estimate vaccine efficacy.</p>",
 "<h2>Which illnesses entered the denominator?</h2>",tablehtml(flow),
 "<p>The flow includes incident episodes before applying the three-day opportunity rule. The main form denominator includes required evaluations; the main swab denominator excludes evaluation-only upper-airway complaints and first-week discretionary systemic reactions. The latter two groups are reported separately. Undocumented fever temperature, uncertain persistence and source conflicts are not coded as ruled out. Multiple indications can occur within one illness; the first supported mandatory indication defines its swab analysis date.</p>",
 tablehtml(srt),"<p>The strict rule requires an evidenced continuous interval for symptoms subject to the 48-hour rule. A sequence of positive checklist entries alone does not prove continuous symptoms. The sensitivity treats four consecutive positive assessment dates for the same item as evidence of a continuous spell and uses day three as its proxy trigger. The delayed-adoption scenario shifts sponsor amendment dates by 14 days unless a verified site history is supplied.</p>",
 "<h2>Blinded clinical review</h2>",paste0("<p>The reviewers-only packet contains ",nrow(review)," randomly ordered episodes. It mixes ",sum(!review$any_sample_d0_7)," respiratory-trigger episodes without any linked diagnostic specimen through day seven and ",sum(review$any_sample_d0_7)," tested controls. Only the custodian receives the linkage key. Original AE terms and coded form items are included where available and maskable. They do not replace complete narratives or contact logs. No human adjudication is claimed.</p>"),
 "<h2>Records for verification</h2><p>Request the site-approved protocol effective dates, the first qualifying contact with staff, the CRF creation audit trail, specimen collection timestamps and types, and the original narratives for discordant or untested cases. A checklist assessment date establishes documentation, not the exact moment a form was opened. Earlier testing or a justified home collection should be reviewed before alleging a missed swab.</p>",
 '<p>Protocol sources: <a href="https://phmpt.org/wp-content/uploads/2025/04/125752_S3_M5_5351_mrna1273p301-16-1-1.pdf">P301 sponsor protocols and amendments</a>; <a href="https://docs.house.gov/meetings/IF/IF02/20200930/111063/HHRG-116-IF02-20200930-SD008.pdf">Amendment 3, §§8.1.2–8.1.3 and Amendment 2 history</a>; <a href="https://cdn.clinicaltrials.gov/large-docs/27/NCT04470427/Prot_000.pdf">version dates and later protocol</a>.</p></body></html>')
 writeLines(enc2utf8(html),file.path(out,"protocol_triggered_illness_report.html"),useBytes=TRUE)
}

if(isTRUE(getOption("phmpt.autorun",TRUE))) tr_run()
