# Moderna COVE | C1 descriptive sensitivity checks | Standalone base R
# Run from the project root in RStudio (select the entire script and run).
# Inputs: the blinded heatmap export, the PCR-by-lab export, and BLA_ADSL_selected.csv
# (from the demographics reference bundle, or an equivalent CSV with the 5 fields below).
# Exports can be unzipped directories or ZIP files. On ambiguous discovery, set:
# options(phmpt.c1_heat="path/to/symptom_testing_heatmaps_blinded_...zip",
#         phmpt.c1_lab="path/to/pcr_testing_by_lab_...zip",
#         phmpt.c1_adsl="path/to/BLA_ADSL_selected.csv", phmpt.root=getwd())
# Output: analysis/symptom_critique_c1. No packages required.
# V02 outputs: analysis/symptom_critique_c1/{coverage,headline_pcr_coverage,V02_table_data}.csv,
# table_inputs_intro/headline_pcr_coverage.csv (three rows), audit_assets/V02_table_data.csv (six rows).
# The last row uses protocol episodes from arm_summary.csv or the matching published V19 row.
# To select another protocol run: options(phmpt.c1_trigger_file=".../arm_summary.csv").
# Use the accompanying six-row analyse_intro_figures.R to redraw the PNG.
# Arm and vaccine-minus-placebo intervals use the SAME 4,000 resampled site vectors.
# Descriptive, not a re-fit of the efficacy model. A positive is a linked diagnostic
# result in the specified window, not necessarily an adjudicated primary case.
# No linked positive includes untested participants/episodes. An absent test is
# never recoded as a negative result. Episodes may repeat within a participant.
# The timing restriction does NOT replace the adjudicated first-event/censoring rows.

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

analyse_symptom_critique <- function(root=getOption("phmpt.root",getwd()),
   heat=getOption("phmpt.c1_heat",""), lab=getOption("phmpt.c1_lab",""),
   adsl=getOption("phmpt.c1_adsl",""), out=file.path(root,"analysis","symptom_critique_c1")) {
  root<-normalizePath(root,mustWork=TRUE)
  resolve<-function(x){if(file.exists(x)) return(normalizePath(x))
    f<-file.path(root,x);if(file.exists(f))return(normalizePath(f));stop("Not found: ",x)}
  locate<-function(supplied,pattern,option){
    if(nzchar(supplied))return(resolve(supplied))
    paths<-unique(unlist(lapply(unique(c(root,file.path(root,"analysis"),file.path(root,"reference_v1_12"))),function(d)
      if(dir.exists(d))list.files(d,pattern=pattern,full.names=TRUE,ignore.case=TRUE))))
    if(length(paths)!=1)stop("Set options(",option,"='...') to one exact file or export. Candidates: ",paste(paths,collapse=", "))
    normalizePath(paths)
  }
  if(!nzchar(heat))heat<-m3_find_export(root,"episode_symptom_source_links.csv")
  if(!nzchar(lab))lab<-m3_find_export(root,"symptom_summary_global.csv")
  heat<-locate(heat,"^symptom_testing_heatmaps_blinded_","phmpt.c1_heat")
  lab<-locate(lab,"^pcr_testing_by_lab_[0-9]","phmpt.c1_lab")
  m3_require_export(heat,"analyse_symptom_testing_heatmaps.R");m3_require_export(lab,"analyse_pcr_testing_by_lab.R")
  adsl<-locate(adsl,"^BLA_ADSL_selected[.]csv$","phmpt.c1_adsl")
  getcsv<-function(bundle,suffix) {
    packed<-!dir.exists(bundle)
    paths<-if(packed)utils::unzip(bundle,list=TRUE)$Name else list.files(bundle,recursive=TRUE)
    match<-paths[paths==suffix|endsWith(paths,paste0("/",suffix))]
    if(length(match)!=1L)stop("Expected one ",suffix," in ",bundle,", found ",length(match))
    if(packed){
      tmp<-tempfile("c1-");dir.create(tmp);on.exit(unlink(tmp,recursive=TRUE))
      utils::unzip(bundle,files=match,exdir=tmp,junkpaths=TRUE)
      con<-file.path(tmp,basename(match))
    }else con<-file.path(bundle,match)
    d<-read.csv(con,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
    names(d)<-sub("^\ufeff","",names(d));d
  }
  tf<-function(x){x<-toupper(as.character(x));if(anyNA(x)||!all(x%in%c("TRUE","FALSE")))stop("Unexpected boolean");x=="TRUE"}
  A<-read.csv(adsl,stringsAsFactors=FALSE,check.names=FALSE,fileEncoding="UTF-8-BOM")
  stopifnot(all(c("SUBJID","PPROTFL","DOSE2DT","SCOV2BL","EFFCODT")%in%names(A)),!anyDuplicated(A$SUBJID))
  adj<-function(subj){j<-match(subj,A$SUBJID);if(anyNA(j))stop("Unmatched ADSL ID");A[j,,drop=FALSE]}
  daywindow<-function(subj,date,last=date){
    d<-adj(subj); t<-as.Date(date); t2<-as.Date(last)
    ok<-d$PPROTFL=="Y" & toupper(d$SCOV2BL)=="NEGATIVE" & t>=as.Date(d$DOSE2DT)+14 & t2<=as.Date(d$EFFCODT)
    ok[is.na(ok)]<-FALSE;ok
  }
  ep<-getcsv(heat,"episodes.csv"); ep<-ep[tf(ep$eligible),]; ep$hit<-tf(ep$tested)
  stopifnot(!anyDuplicated(ep$episode_id))
  el<-getcsv(heat,"episode_visit_pcr_links.csv"); el<-el[tf(el$episode_eligible),]
  stopifnot(all(el$onset_to_test>=0 & el$onset_to_test<=7),all(el$episode_id%in%ep$episode_id))
  ep$pos<-ep$episode_id%in%el$episode_id[toupper(el$test_result)=="POSITIVE" & !tf(el$result_conflict)]
  stopifnot(all(ep$hit[ep$pos]),identical(sort(ep$episode_id[ep$hit]),sort(unique(el$episode_id))))
  ep$pp<-daywindow(ep$subject,ep$start,as.Date(ep$start)+7)
  ep$resp<-vapply(strsplit(ep$onset_groups,";"),function(x)any(trimws(x)%in%c("COUGH","DYSPNEA")),logical(1))
  pfx<-"symptom_links/Prod_BLA_SNAPSHOT_STRICT_BLINDED/"
  anchors<-getcsv(lab,paste0(pfx,"symptom_evidence_dates.csv"))
  links<-getcsv(lab,paste0(pfx,"symptom_date_to_pcr_links.csv"));links<-links[tf(links$diagnostic),]
  tests<-getcsv(lab,"Prod/pcr_record_audit.csv")
  ky<-function(s,r)paste(s,r,sep=":")
  testkey<-ky(tests$subject,tests$record_id);stopifnot(!anyDuplicated(testkey))
  j<-match(ky(links$subject,links$record_id),testkey)
  stopifnot(!anyNA(j),all(tf(tests$usable_dated[j])),!anyDuplicated(anchors$anchor_id),all(links$anchor_id%in%anchors$anchor_id))
  anchors$hit<-anchors$anchor_id%in%links$anchor_id
  anchors$pos<-anchors$anchor_id%in%links$anchor_id[tf(tests$positive[j])&!tf(tests$result_conflict[j])]
  anchors$form<-tf(anchors$has_covid_form)
  stopifnot("early_postdose"%in%names(anchors),!anyDuplicated(unique(anchors[c("subject","arm")])$subject))
  ai<-match(links$anchor_id,anchors$anchor_id)
  stopifnot(all(links$subject==anchors$subject[ai]),all(as.Date(links$test_date)>=as.Date(anchors$evidence_date[ai])),
            all(as.Date(links$test_date)<=as.Date(anchors$link_end[ai])),
            all(as.Date(links$test_date)<=as.Date(anchors$evidence_date[ai])+7))
  # The upstream BY_CUTOFF participant metric permits truncated 7-day follow-up.
  # Verify there are no diagnostic links beyond the individual's efficacy cutoff
  # in the PP subset; if there are, rebuild links within that boundary.
  anchors$pp<-daywindow(anchors$subject,anchors$evidence_date)
  ix<-match(links$anchor_id,anchors$anchor_id)
  stopifnot(all(as.Date(links$test_date[anchors$pp[ix]])<=as.Date(adj(links$subject[anchors$pp[ix]])$EFFCODT)))
  people<-function(d)aggregate(d[c("hit","pos","form")],d[c("subject","arm")],any)
  s<-people(anchors); sp<-people(anchors[anchors$pp,])
  old<-getcsv(lab,"symptom_subjects/Prod_BLA_SNAPSHOT_STRICT_BLINDED_SYMPTOMS__COVID_FORM_PLUS_AE__BY_CUTOFF__DIAGNOSTIC.csv")
  k<-match(s$subject,old$subject);stopifnot(nrow(old)==nrow(s),!anyNA(k),all(s$hit==tf(old$any_tested[k])),all(s$form==tf(old$has_covid_form_evidence[k])))
  # Paired site-cluster bootstrap: resample entire sites, both arms and all
  # repeated episodes together. Percentile intervals describe unadjusted gaps.
  sites<-sort(unique(as.character(ep$SITEID)));stopifnot(all(substr(ep$subject,1,5)==ep$SITEID))
  set.seed(30114); weights<-stats::rmultinom(4000L,length(sites),rep(1/length(sites),length(sites)))
  sumrow<-function(label,x){
    arms<-c("mRNA-1273","Placebo"); n<-vapply(arms,function(a)sum(x$arm==a),integer(1))
    y<-vapply(arms,function(a)sum(x$arm==a & x$hit),integer(1));pc<-100*y/n
    sj<-match(substr(x$subject,1,5),sites);stopifnot(!anyNA(sj))
    bo<-lapply(arms,function(a){
      ix<-x$arm==a;num<-tabulate(sj[ix&x$hit],length(sites));den<-tabulate(sj[ix],length(sites))
      as.vector(crossprod(num,weights))/as.vector(crossprod(den,weights))
    })
    ci<-unname(quantile(100*(bo[[1]]-bo[[2]]),c(.025,.975),na.rm=TRUE))
    ci_arm<-lapply(bo,function(v)unname(quantile(100*v,c(.025,.975),na.rm=TRUE)))
    data.frame(scope=label,vaccine_tested=y[1],vaccine_total=n[1],vaccine_pct=unname(pc[1]),
      placebo_tested=y[2],placebo_total=n[2],placebo_pct=unname(pc[2]),gap_pp=unname(pc[1]-pc[2]),low95=ci[1],high95=ci[2],
      vaccine_low95=ci_arm[[1]][1],vaccine_high95=ci_arm[[1]][2],
      placebo_low95=ci_arm[[2]][1],placebo_high95=ci_arm[[2]][2],row.names=NULL)
  }
  panel<-function(label,x) do.call(rbind,list(sumrow(label,x),sumrow(paste(label,"| linked positive"),x[x$pos,]),
    sumrow(paste(label,"| no linked positive"),x[!x$pos,])))
  results<-rbind(panel("Participants: blinded",s),panel("Participants: PP, dose 2 + 14 days",sp),
    panel("Episodes: blinded",ep),panel("Episodes: PP, dose 2 + 14 days",ep[ep$pp,]),
    panel("Episodes: PP, cough/dyspnea recorded at onset",ep[ep$pp & ep$resp,]))
  form<-as.data.frame(with(s[s$form,],table(arm,pos)));names(form)[2:3]<-c("linked_positive","participants_with_COVID_forms")
  u<-as.data.frame(with(sp[!sp$form & !sp$hit,],table(arm)));names(u)[2]<-"AE_only_without_linked_diagnostic_test"
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  write.csv(form,file.path(out,"COVID_form_participants_by_result.csv"),row.names=FALSE)
  write.csv(u,file.path(out,"PP_timing_AE_only_untested_people.csv"),row.names=FALSE)
  selected<-match(c("Participants: blinded | no linked positive","Participants: PP, dose 2 + 14 days | no linked positive"),results$scope)
  stopifnot(!anyNA(selected))
  # Fix membership by absence of ANY linked positive during full blinded follow-up.
  # Drop early symptom dates, not entire participants; link a test only to a retained date.
  remaining<-anchors[!(anchors$subject%in%s$subject[s$pos]) & !tf(anchors$early_postdose),]
  outside<-people(remaining);stopifnot(!any(outside$pos))
  scoped<-list(s[!s$pos,],outside,sp[!sp$pos,])
  sensitivity<-sumrow("Participants: blinded | no linked positive | excluding post-dose days 0-6",outside)
  headline<-rbind(results[selected[1],],sensitivity,results[selected[2],]);headline$interval_method<-"Paired site-cluster percentile bootstrap; both arms retained in each site"
  headline$bootstrap_draws<-ncol(weights);headline$n_sites_in_resampling_frame<-length(sites);headline$seed<-30114L
  for(i in seq_len(nrow(headline))){
    x<-scoped[[i]]
    headline$sites_with_eligible_participants[i]<-length(unique(substr(x$subject,1,5)))
  }
  # Row 6 counts illness EPISODES conditional on a checklist, not unique participants.
  # Prefer an explicitly selected protocol export. Also support projects that have
  # the numeric V19 table but not the full aggregate input directory.
  absolute<-function(x)if(grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)",x))x else file.path(root,x)
  intro<-absolute(getOption("phmpt.intro_input_dir","table_inputs_intro"))
  assets<-absolute(getOption("phmpt.audit_asset_dir","audit_assets"))
  override<-getOption("phmpt.c1_trigger_file",getOption("phmpt.intro_trigger_file",""))
  if(!is.null(override) && nzchar(override)){
    trigger_file<-absolute(override)
  }else{
    candidates<-c(file.path(intro,"arm_summary.csv"),file.path(root,"table_inputs_v022","protocol","arm_summary.csv"),
      file.path(root,"audit_assets","V19_table_data.csv"))
    found<-candidates[file.exists(candidates)]
    if(!length(found))stop("To export all six V02 rows supply options(phmpt.c1_trigger_file='.../arm_summary.csv') or the numeric audit_assets/V19_table_data.csv.")
    trigger_file<-found[1L]
  }
  protocol<-read.csv(trigger_file,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
  arms<-c("mRNA-1273","Placebo")
  if(all(c("decision","vaccine_y","vaccine_n","placebo_y","placebo_n")%in%names(protocol))){
    v19<-protocol[protocol$decision=="Swab collected on or after the first checklist date and by day 3, among swab-indicated episodes with a checklist by day 3",,drop=FALSE]
    if(nrow(v19)!=1L)stop("Expected exactly one V19 checklist-conditional swab row in ",trigger_file)
    protocol<-data.frame(scope="All",arm=arms,form_then_swab=c(v19$vaccine_y,v19$placebo_y),form=c(v19$vaccine_n,v19$placebo_n))
  }
  required<-c("scope","arm","form","form_then_swab")
  if(!all(required%in%names(protocol)))stop("Not a compatible protocol summary or V19 data table: ",trigger_file)
  tr<-protocol[protocol$scope=="All",]
  stopifnot(nrow(tr)==2L,setequal(tr$arm,arms));tr<-tr[match(arms,tr$arm),]
  lab_global<-getcsv(lab,"symptom_summary_global.csv")
  need<-c("release","period","scope","stratum","granularity","endpoint","arm","n_eligible","n_tested")
  stopifnot(all(need%in%names(lab_global)))
  g<-subset(lab_global,release=="Prod"&period=="BLA_SNAPSHOT"&scope=="STRICT_BLINDED"&granularity=="GLOBAL"&endpoint=="Any PCR")
  row<-function(label,y,n,unit="Participants",h=NULL){
    stopifnot(length(y)==2,length(n)==2,all(is.finite(c(y,n))),all(c(y,n)==floor(c(y,n))),all(n>0 & y>=0 & y<=n))
    rates<-100*y/n; lo<-hi<-rep(NA_real_,2); gap<-gl<-gh<-NA_real_
    if(!is.null(h)){
      lo<-c(h$vaccine_low95,h$placebo_low95);hi<-c(h$vaccine_high95,h$placebo_high95)
      gap<-h$gap_pp;gl<-h$low95;gh<-h$high95
      stopifnot(all(is.finite(c(lo,hi,gap,gl,gh))),all(lo<=hi & lo>=0 & hi<=100),gl<=gh,
        max(abs(rates-c(h$vaccine_pct,h$placebo_pct)))<1e-8,abs(gap-(rates[1]-rates[2]))<1e-8)
    }
    data.frame(label,unit,vaccine_tested=y[1],vaccine_total=n[1],vaccine_pct=rates[1],vaccine_low95=lo[1],vaccine_high95=hi[1],
      placebo_tested=y[2],placebo_total=n[2],placebo_pct=rates[2],placebo_low95=lo[2],placebo_high95=hi[2],
      gap_pp=gap,gap_low95=gl,gap_high95=gh,row.names=NULL)
  }
  labels<-c("Symptoms on either record","Symptoms on a COVID checklist","All blinded follow-up; no linked positive PCR",
    "Same full-window no-positive cohort; symptom dates outside days 0-6 after either dose","Per-protocol window; no linked positive PCR",
    "Once a checklist existed: swab by day 3")
  rows<-list()
  for(i in 1:2){
    pop<-c("COVID_FORM_PLUS_AE","COVID_FORM_ONLY")[i]
    z<-g[g$stratum==paste0("SYMPTOMS__",pop,"__BY_CUTOFF__DIAGNOSTIC"),];stopifnot(nrow(z)==2L,setequal(z$arm,arms))
    z<-z[match(arms,z$arm),]
    a<-if(i==1L)s else people(anchors[anchors$form,])
    stopifnot(all(z$n_eligible==vapply(arms,function(ar)sum(a$arm==ar),integer(1))),
              all(z$n_tested==vapply(arms,function(ar)sum(a$arm==ar & a$hit),integer(1))))
    rows[[i]]<-row(labels[i],z$n_tested,z$n_eligible)
  }
  for(i in seq_len(nrow(headline))){h<-headline[i,];rows[[i+2L]]<-row(labels[i+2L],c(h$vaccine_tested,h$placebo_tested),c(h$vaccine_total,h$placebo_total),h=h)}
  rows[[6]]<-row(labels[6],tr$form_then_swab,tr$form,"Illness episodes")
  v02<-do.call(rbind,rows);stopifnot(nrow(v02)==6L,nrow(headline)==3L)
  # Keep all previous scopes in the same order. Add the first-week sensitivity to the detailed export too.
  results<-rbind(results,sensitivity)
  for(d in c(intro,assets))dir.create(d,recursive=TRUE,showWarnings=FALSE)
  csv<-function(x,file)write.csv(x,file,row.names=FALSE,na="",fileEncoding="UTF-8")
  csv(results,file.path(out,"coverage.csv"))
  for(d in unique(c(out,intro)))csv(headline,file.path(d,"headline_pcr_coverage.csv"))
  for(d in unique(c(out,assets)))csv(v02,file.path(d,"V02_table_data.csv"))
  csv(lab_global,file.path(intro,"symptom_summary_global.csv"))
  if(!identical(normalizePath(trigger_file,winslash="/"),normalizePath(file.path(intro,"arm_summary.csv"),winslash="/",mustWork=FALSE)))
    csv(protocol,file.path(intro,"arm_summary.csv"))
  csv(data.frame(role=c("heatmap","lab","ADSL","protocol"),input=c(heat,lab,adsl,normalizePath(trigger_file,winslash="/"))),
      file.path(out,"V02_sources.csv"))
  print(headline[,c("scope","vaccine_tested","vaccine_total","placebo_tested","placebo_total","gap_pp","low95","high95")],row.names=FALSE,digits=4)
  message("Headline CSV (3 rows): ",file.path(intro,"headline_pcr_coverage.csv"),
    "\nV02 figure CSV (6 rows): ",file.path(assets,"V02_table_data.csv"),"\nDetailed results: ",normalizePath(out))
  invisible(list(coverage=results,headline=headline,v02=v02,form=form,untested=u))
}
if(isTRUE(getOption("phmpt.autorun",TRUE)))analyse_symptom_critique()
