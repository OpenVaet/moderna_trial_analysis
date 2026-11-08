# V20: two-panel RR funnels; pooled placebo prevalence by default. Set options(phmpt.v17_p0_randomized=0.10) for a fixed benchmark in panel A.
# Moderna P301 - exploratory site testing audit | v1.2.0 | 2026-10-21
# Stand-alone base R >=4.3: select all and run in RStudio from the project root.
# Generates report_insert.html/txt, numbered prose, tables, references and five
# illustrations as PDF/PNG (and SVG with Cairo). Does not modify the draft.
# Primary high-positivity rule: county weekly NAAT positivity >=10%.
# Command line: Rscript analyse_site_testing_imbalance.R <export_folder> <output_folder>
# To load functions without running: options(phmpt.autorun=FALSE).


audit_release_dir <- function(root, legacy) {
  modern<-paste0(legacy,"_moderna")
  file.path(root,ifelse(dir.exists(file.path(root,modern)),modern,legacy))
}

DATA_DIR <- getOption('phmpt.audit_data_dir','')   # blank: discover latest BLINDED heatmap export
OUTPUT_DIR <- getOption('phmpt.audit_output_dir','') # blank: <export>/positivity_audit
PHMPT_REPORT_START <- getOption('phmpt.report_start_paragraph',25L)
PHMPT_REFERENCE_START <- getOption('phmpt.report_reference_start',32L)
PHMPT_REPORT_ONLY <- getOption('phmpt.report_only',FALSE) # reuse existing, hash-matched model CSVs

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

v17_check_blinded_source <- function(data_dir) {
  m3_require_export(data_dir,"analyse_symptom_testing_heatmaps.R")
  pars<-read.csv(file.path(data_dir,"run_parameters.csv"),stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
  val<-function(k)pars$value[pars$parameter==k]
  if(!identical(val("window"),"BLINDED")||!identical(val("release"),"Prod"))
    stop("V20 requires the blinded production heatmap export.",call.=FALSE)
  ver<-"1.1.0-M3";reported<-val("dictionary_version")
  if(length(reported)&&!identical(reported,ver))stop("Heatmap run parameters disagree with the mapping marker.",call.=FALSE)
  for(n in c("symptom_event_evidence.csv","episode_symptom_source_links.csv")) {
    f<-file.path(data_dir,n)
    if(!file.exists(f))stop("Incomplete heatmap export: missing ",n,call.=FALSE)
    cols<-names(read.csv(f,nrows=0,check.names=FALSE,fileEncoding="UTF-8-BOM"))
    if(!"dictionary_version"%in%cols)stop("No mapping version on source rows in ",n,". Regenerate the blinded heatmap.",call.=FALSE)
    v<-read.csv(f,colClasses=ifelse(cols=="dictionary_version","character","NULL"),fileEncoding="UTF-8-BOM")[[1L]]
    if(!length(v)||anyNA(v)||any(v!=ver))stop("Unverified or mixed symptom mapping in ",n,".",call.=FALSE)
  }
  invisible(TRUE)
}

phmpt_find_export <- function(data_dir='',root=getwd()) {
  required<-c('episodes.csv','site_background_used.csv','site_region_coverage.csv',
    'subject_followup.csv','site_arm_weekly_testing.csv','pcr_records_reconciled.csv','illness_visits.csv','run_parameters.csv')
  valid<-function(p)all(file.exists(file.path(p,required)))
  if(nzchar(data_dir)) {
    if(!valid(data_dir))stop('DATA_DIR must contain a complete heatmap CSV export: ',data_dir)
    return(normalizePath(data_dir,winslash='/'))
  }
  if(valid(root))return(normalizePath(root,winslash='/'))
  paths<-list.files(file.path(root,'analysis'),pattern='^episodes[.]csv$',recursive=TRUE,full.names=TRUE)
  dirs<-unique(dirname(paths));dirs<-dirs[vapply(dirs,valid,logical(1))]
  if(!length(dirs))stop('No complete heatmap export found. Run analyse_symptom_testing_heatmaps.R first, or set DATA_DIR to its output folder.')
  blinded<-vapply(dirs,function(p){z<-read.csv(file.path(p,'run_parameters.csv'),stringsAsFactors=FALSE);any(z$parameter=='window'&z$value=='BLINDED')},logical(1))
  if(any(blinded))dirs<-dirs[blinded]
  chosen<-dirs[which.max(file.info(file.path(dirs,'episodes.csv'))$mtime)]
  message('Using heatmap export: ',chosen)
  normalizePath(chosen,winslash='/')
}

analyse_site_testing_imbalance <- function(data_dir='',output_dir='',report_only=FALSE,
  start_paragraph=25L,reference_start=32L,audit_results_dir='') {
  if(!nzchar(data_dir))data_dir<-m3_find_export(getOption("phmpt.root",getwd()),"episode_symptom_source_links.csv")
  DATA_DIR<-phmpt_find_export(data_dir)
  v17_check_blinded_source(DATA_DIR)
  OUTPUT_DIR<-if(nzchar(output_dir))output_dir else file.path(DATA_DIR,'positivity_audit')
  dir.create(OUTPUT_DIR,recursive=TRUE,showWarnings=FALSE)
  if(!isTRUE(report_only)) {
read<-function(n)read.csv(file.path(DATA_DIR,n),stringsAsFactors=FALSE,check.names=FALSE,fileEncoding='UTF-8-BOM',na.strings=c('','NA'))
save_csv<-function(z,n)write.csv(z,file.path(OUTPUT_DIR,n),row.names=FALSE,na='')
e<-read('episodes.csv');b<-read('site_background_used.csv');m<-read('site_region_coverage.csv');a<-read('subject_followup.csv');cells<-read('site_arm_weekly_testing.csv')
stopifnot(!anyDuplicated(e$episode_id),!anyDuplicated(a$subject),!anyDuplicated(m$SITEID),!anyDuplicated(b[c('SITEID','week')]))
x<-merge(e[e$eligible,],b,by=c('SITEID','week'),all.x=TRUE,sort=FALSE)
x<-merge(x,m[c('SITEID','region_name')],by='SITEID',all.x=TRUE,sort=FALSE)
stopifnot(nrow(x)==sum(e$eligible),!anyNA(x$positivity),all(x$positivity>=0&x$positivity<=100),!anyNA(x$region_name))
# Check every site's episode/test totals against the exported heatmap cells.
v1<-aggregate(cbind(episodes=rep(1L,nrow(x)),tested=as.integer(x$tested)),x[c('SITEID','arm','week')],sum)
v2<-merge(v1,cells[c('SITEID','arm','week','symptom_episodes','tested_episodes')],by=c('SITEID','arm','week'),all.x=TRUE)
stopifnot(all(v2$episodes==v2$symptom_episodes),all(v2$tested==v2$tested_episodes))
x<-x[order(x$subject,x$start,x$episode_id),];x$vaccine<-as.integer(x$arm=='mRNA-1273');x$high10<-as.integer(x$positivity>=10)
x$early<-as.integer(x$early_postdose);x$month<-substr(x$start,1,7)
patterns<-sort(unique(unlist(strsplit(x$onset_groups,'; ',fixed=TRUE))))
for(p in patterns)x[[paste0('onset_',p)]]<-as.integer(vapply(strsplit(x$onset_groups,'; ',fixed=TRUE),function(z)p%in%z,logical(1)))
# County-specific upper quartile, one county/week each; not weighted by site count.
br<-unique(b[c('region_key','week','positivity')]);q75<-tapply(br$positivity,br$region_key,quantile,probs=.75,names=FALSE)
x$county_q75<-unname(q75[x$region_key]);x$high_q75<-x$positivity>=x$county_q75
first<-function(z)z[!duplicated(z$subject),,drop=FALSE]
sites<-sort(m$SITEID)
wilson<-function(k,n){if(!n)return(c(NA_real_,NA_real_));z<-qnorm(.975);p<-k/n;d<-1+z*z/n;c((p+z*z/(2*n)-z*sqrt(p*(1-p)/n+z*z/(4*n*n)))/d,(p+z*z/(2*n)+z*sqrt(p*(1-p)/n+z*z/(4*n*n)))/d)}
compare<-function(z,outcome='tested') {
 v<-z[z$vaccine==1,];p<-z[z$vaccine==0,];vn<-nrow(v);pn<-nrow(p);vk<-sum(v[[outcome]]);pk<-sum(p[[outcome]])
 mat<-rbind(c(vk,vn-vk),c(pk,pn-pk));pv<-chi<-cp<-or<-lo<-hi<-NA_real_
 if(vn&&pn&&all(colSums(mat)>0)){
  ft<-fisher.test(mat);pv<-ft$p.value;or<-unname(ft$estimate);lo<-ft$conf.int[1];hi<-ft$conf.int[2]
  ct<-suppressWarnings(chisq.test(mat,correct=FALSE));chi<-unname(ct$statistic);cp<-ct$p.value
 }
 vc<-wilson(vk,vn);pc<-wilson(pk,pn)
 data.frame(vaccine_tested=vk,vaccine_n=vn,vaccine_pct=if(vn)100*vk/vn else NA_real_,
  placebo_tested=pk,placebo_n=pn,placebo_pct=if(pn)100*pk/pn else NA_real_,
  difference_pp=if(vn&&pn)100*(vk/vn-pk/pn) else NA_real_,
  odds_ratio=or,or_lower=lo,or_upper=hi,fisher_p=pv,pearson_chi_square=chi,pearson_p=cp,
  vaccine_ci_lower=100*vc[1],vaccine_ci_upper=100*vc[2],placebo_ci_lower=100*pc[1],placebo_ci_upper=100*pc[2],
  unique_subjects=length(unique(z$subject)),weeks=length(unique(z$week)),
  min_positivity=if(nrow(z))min(z$positivity) else NA_real_,max_positivity=if(nrow(z))max(z$positivity) else NA_real_,row.names=NULL)
}
week_mh<-function(z,outcome='tested',strata='week') {
 if(!nrow(z))return(data.frame(mh_or=NA_real_,mh_p=NA_real_,informative_strata=0,informative_episodes=0))
 groups<-split(z,interaction(z[strata],drop=TRUE));tabs<-list();ne<-0L
 for(g in groups){t<-table(factor(g$vaccine,levels=c(1,0)),factor(g[[outcome]],levels=c(TRUE,FALSE)))
  if(all(rowSums(t)>0)&&all(colSums(t)>0)){tabs[[length(tabs)+1L]]<-t;ne<-ne+nrow(g)}}
 if(!length(tabs))return(data.frame(mh_or=NA_real_,mh_p=NA_real_,informative_strata=0,informative_episodes=0))
 ar<-array(unlist(tabs),dim=c(2,2,length(tabs)))
 aa<-ar[1,1,];bb<-ar[1,2,];cc<-ar[2,1,];dd<-ar[2,2,];nn<-aa+bb+cc+dd
 expected<-(aa+bb)*(aa+cc)/nn
 variance<-(aa+bb)*(cc+dd)*(aa+cc)*(bb+dd)/(nn^2*(nn-1))
 statistic<-sum(aa-expected)^2/sum(variance)
 data.frame(mh_or=sum(aa*dd/nn)/sum(bb*cc/nn),mh_p=pchisq(statistic,1,lower.tail=FALSE),informative_strata=length(tabs),informative_episodes=ne)
}
sets<-list(HIGH10=x$positivity>=10,LOW10=x$positivity<10,HIGH10_OUTSIDE_POSTDOSE=x$positivity>=10&!x$early_postdose,
 HIGH15=x$positivity>=15,HIGH15_OUTSIDE_POSTDOSE=x$positivity>=15&!x$early_postdose,
 COUNTY_UPPER_QUARTILE=x$high_q75,ALL=rep(TRUE,nrow(x)))
all_results<-list();totals<-list();members<-list()
for(nm in names(sets)) {
 z<-x[sets[[nm]],,drop=FALSE]
 for(unit in c('ALL_EPISODES','FIRST_EPISODE_PER_SUBJECT')) {
  zz<-if(unit=='ALL_EPISODES')z else first(z)
  totals[[length(totals)+1L]]<-cbind(subset=nm,unit=unit,compare(zz),week_mh(zz,strata=c('SITEID','week')))
  # primary formal site comparisons use independent participants, with Holm/BH
  # across all 99 sites (untestable sites conservatively contribute p=1).
  rows<-lapply(sites,function(s){d<-zz[zz$SITEID==s,];cbind(SITEID=s,region_name=m$region_name[match(s,m$SITEID)],
   subset=nm,unit=unit,compare(d),week_mh(d))})
  tab<-do.call(rbind,rows);pval<-ifelse(is.finite(tab$fisher_p),tab$fisher_p,1)
  tab$fisher_q_BH_99<-p.adjust(pval,'BH');tab$fisher_p_Holm_99<-p.adjust(pval,'holm')
  tab$mh_q_BH_99<-p.adjust(ifelse(is.finite(tab$mh_p),tab$mh_p,1),'BH')
  tab$adequate_10_each<-tab$vaccine_n>=10&tab$placebo_n>=10
  all_results[[length(all_results)+1L]]<-tab
  if(unit=='FIRST_EPISODE_PER_SUBJECT'&&nm=='HIGH10')members[[nm]]<-zz
 }
}
site<-do.call(rbind,all_results);tot<-do.call(rbind,totals)
save_csv(site,'site_comparisons.csv');save_csv(tot,'overall_comparisons.csv')
save_csv(members[[1]],'primary_first_high_episode_subjects.csv')
# Within-high-positivity checks by onset collection route and symptom pattern.
by_source<-do.call(rbind,lapply(c('HIGH10','LOW10'),function(nm)do.call(rbind,lapply(unique(x$onset_source),function(src){z<-x[sets[[nm]]&x$onset_source==src,];cbind(subset=nm,onset_source=src,compare(z))}))))
save_csv(by_source,'onset_source_comparisons.csv')
by_pattern<-do.call(rbind,lapply(patterns,function(p){z<-x[x$high10==1&x[[paste0('onset_',p)]]==1,];cbind(symptom=p,compare(z))}))
save_csv(by_pattern,'onset_pattern_comparisons.csv')
# Global adjusted logistic models, robust covariance summed by subject and by
# site. Site clusters contain all repeated subjects within that site; the main
# report uses site-clustered estimates with t(G-1) reference. Subject-clustered
# estimates are supplementary.
# HC0 sandwich with G/(G-1) cluster correction. Fits modelled log odds,
# not vaccine efficacy or a causal estimate of investigator behavior.
cluster_glm<-function(formula,d,label) {
 warnings_seen<-character();fit<-withCallingHandlers(glm(formula,data=d,family=binomial(),x=TRUE,y=TRUE,control=glm.control(maxit=100)),warning=function(w){warnings_seen<<-c(warnings_seen,conditionMessage(w));invokeRestart('muffleWarning')})
 keep<-which(!is.na(coef(fit)));X<-fit$x[,keep,drop=FALSE];mu<-fitted(fit);bet<-coef(fit)[keep]
 bread<-solve(crossprod(X,X*as.numeric(mu*(1-mu))))
 tables<-list();covariances<-list()
 for(unit in c('subject','SITEID')) {
  score<-rowsum(X*as.numeric(fit$y-mu),d[[unit]],reorder=FALSE);G<-nrow(score)
  V<-bread%*%crossprod(score)%*%bread*G/(G-1);covariances[[unit]]<-V
  se<-sqrt(pmax(diag(V),0));reference<-if(unit=='SITEID')'t(G-1)' else 'normal'
  crit<-if(unit=='SITEID')qt(.975,G-1) else qnorm(.975)
  pv<-if(unit=='SITEID')2*pt(-abs(bet/se),G-1) else 2*pnorm(-abs(bet/se))
  tables[[unit]]<-data.frame(model=label,term=names(bet),log_odds=bet,se=se,OR=exp(bet),ci_lower=exp(bet-crit*se),ci_upper=exp(bet+crit*se),p=pv,
    episodes=nrow(d),subjects=length(unique(d$subject)),clusters=G,cluster_unit=unit,reference=reference,converged=fit$converged,warnings=paste(unique(warnings_seen),collapse=' | '),row.names=NULL)
 }
 list(fit=fit,V=covariances$subject,V_site=covariances$SITEID,table=do.call(rbind,tables))
}
symptoms<-paste0('onset_',patterns,collapse=' + ')
model_specs<-list(
 HIGH10_SITE_WEEK=list(d=x[x$high10==1,],f='tested ~ vaccine + factor(SITEID) + factor(week)'),
 HIGH10_SYMPTOMS_TIMING=list(d=x[x$high10==1,],f=paste('tested ~ vaccine + factor(SITEID) + factor(week) + early +',symptoms)),
 HIGH10_SOURCE_SENSITIVITY=list(d=x[x$high10==1,],f=paste('tested ~ vaccine + factor(SITEID) + factor(week) + early + factor(onset_source) +',symptoms)),
 HIGH10_OUTSIDE_POSTDOSE_ADJUSTED=list(d=x[x$high10==1&x$early==0,],f=paste('tested ~ vaccine + factor(SITEID) + factor(week) +',symptoms)),
 HIGH_LOW_INTERACTION=list(d=x,f=paste('tested ~ vaccine * high10 + factor(SITEID) + factor(week) + early +',symptoms)),
 HIGH_LOW_SOURCE_INTERACTION=list(d=x,f=paste('tested ~ vaccine * high10 + factor(SITEID) + factor(week) + early + factor(onset_source) +',symptoms)))
models<-list()
for(nm in names(model_specs)) {
 message('Model: ',nm)
 z<-model_specs[[nm]];result<-tryCatch(cluster_glm(as.formula(z$f),z$d,nm),error=function(e)e)
 if(inherits(result,'error')){writeLines(conditionMessage(result),file.path(OUTPUT_DIR,paste0(nm,'_MODEL_ERROR.txt')));next}
 models[[nm]]<-result;save_csv(result$table,paste0(nm,'_coefficients.csv'))
}
key<-do.call(rbind,lapply(models,function(z)z$table[z$table$term%in%c('vaccine','high10','vaccine:high10'),]))
save_csv(key,'adjusted_model_key_coefficients.csv')
# Optional checkpoint; CSV outputs are sufficient for reviewing the findings.
if(isTRUE(getOption('phmpt.save_model_objects',FALSE)))
 saveRDS(list(data=x,site=site,totals=tot,models=models),file.path(OUTPUT_DIR,'analysis_objects.rds'))
capture.output(sessionInfo(),file=file.path(OUTPUT_DIR,'sessionInfo.txt'))
cat('\nPrimary site screen, first high-period episode per participant:\n')
primary<-site[site$subset=='HIGH10'&site$unit=='FIRST_EPISODE_PER_SUBJECT',]
print(primary[order(primary$fisher_p),c('SITEID','region_name','vaccine_tested','vaccine_n','placebo_tested','placebo_n','difference_pp','fisher_p','fisher_q_BH_99','mh_p')][1:15,],row.names=FALSE)
print(key,row.names=FALSE)

# Case review: exploratory record-level follow-up of the nominal site signals.
high<-x[x$high10==1,]
pcr<-read('pcr_records_reconciled.csv')
visits<-read('illness_visits.csv')
episodes<-e
priority<-c('US325','US352','US313','US372','US359')
save_csv(high[high$SITEID%in%priority,],'priority_site_high10_episodes.csv')
route<-aggregate(cbind(episodes=rep(1,nrow(high)),tested=as.integer(high$tested),early=as.integer(high$early_postdose)),high[c('SITEID','arm','onset_source')],sum)
save_csv(route,'site_onset_source_high10.csv')
case<-high[high$SITEID=='US325'&high$vaccine==1&high$onset_source=='COVID_FORM'&!high$tested,]
case<-case[c('episode_id','subject','start','last_symptom','onset_groups','positivity','early_postdose','tested')]
case$first_later_diagnostic_pcr<-'';case$days_to_later_pcr<-NA_integer_;case$later_pcr_result<-'';case$within_same_episode_span<-FALSE
case$next_test_in_other_episode<-FALSE;case$not_done_in_first_7_days<-FALSE
for(i in seq_len(nrow(case))) {
 z<-pcr[pcr$subject==case$subject[i]&!is.na(pcr$test_date),];day<-as.Date(z$test_date)-as.Date(case$start[i])
 case$not_done_in_first_7_days[i]<-any(z$not_done&day>=0&day<=7)
 zz<-z[z$usable_diagnostic&day>=0,];zz<-zz[order(zz$test_date),]
 if(nrow(zz)){
  d<-zz$test_date[1];case$first_later_diagnostic_pcr[i]<-d;case$days_to_later_pcr[i]<-as.integer(as.Date(d)-as.Date(case$start[i]))
  case$later_pcr_result[i]<-zz$result[1];case$within_same_episode_span[i]<-d<=case$last_symptom[i]
  other<-episodes[episodes$subject==case$subject[i]&episodes$episode_id!=case$episode_id[i],]
  case$next_test_in_other_episode[i]<-any(other$start<=d&other$last_symptom>=d)
 }
}
save_csv(case,'US325_covid_form_untested_review.csv')
ends<-do.call(rbind,lapply(c('tested','independent_visit_pcr','diagnostic_pcr_by_horizon','any_nonbaseline_pcr_by_horizon'),function(end)do.call(rbind,lapply(c('mRNA-1273','Placebo'),function(a){z<-high[high$arm==a,];data.frame(endpoint=end,arm=a,tested=sum(z[[end]]),episodes=nrow(z),pct=100*mean(z[[end]]))}))))
save_csv(ends,'endpoint_sensitivity.csv')
print(case,row.names=FALSE)
selected_subjects<-case$subject
save_csv(pcr[pcr$subject%in%selected_subjects,],'US325_selected_subject_pcr_records.csv')
save_csv(visits[visits$subject%in%selected_subjects,],'US325_selected_subject_illness_visits.csv')

# Equalize onset-source composition descriptively (no causal interpretation).
src<-by_source[by_source$subset=='HIGH10',]
src$pooled_weight<-(src$vaccine_n+src$placebo_n)/sum(src$vaccine_n+src$placebo_n)
standardized<-data.frame(vaccine_pct=sum(src$pooled_weight*src$vaccine_pct),placebo_pct=sum(src$pooled_weight*src$placebo_pct))
standardized$difference_pp<-standardized$vaccine_pct-standardized$placebo_pct
save_csv(standardized,'source_standardized_rates.csv')
inputs<-c('episodes.csv','site_background_used.csv','site_region_coverage.csv','subject_followup.csv','site_arm_weekly_testing.csv','pcr_records_reconciled.csv','illness_visits.csv')
save_csv(data.frame(file=inputs,md5=unname(tools::md5sum(file.path(DATA_DIR,inputs)))),'input_checksums.csv')
message('Audit outputs written to: ',normalizePath(OUTPUT_DIR,winslash='/'))
  }
  if(!nzchar(audit_results_dir))audit_results_dir<-OUTPUT_DIR
  phmpt_anomaly_examples(DATA_DIR,OUTPUT_DIR,high_cutoff=10)
  result<-phmpt_report_insert(DATA_DIR,OUTPUT_DIR,start_paragraph,reference_start,audit_results_dir)
  writeLines(c('COMPLETE',paste0('script_version=1.2.0'),paste0('report_only=',report_only)),file.path(OUTPUT_DIR,'REPORT_STATUS.txt'))
  invisible(result)
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
  screen$fisher_q_BH<-p.adjust(ifelse(is.finite(screen$fisher_p),screen$fisher_p,1),'BH');screen$fisher_p_bonferroni<-p.adjust(ifelse(is.finite(screen$fisher_p),screen$fisher_p,1),'bonferroni');csv(screen,'first_high_episode_site_screen.csv')
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
    txt(.065,.19,sprintf('Difference: %+.1f percentage points  |  Fisher p=%s  |  Bonferroni p=%s (%d sites)',z$v_pct-z$p_pct,pfmt(z$fisher_p),formatC(z$fisher_p_bonferroni,format='f',digits=3),M),11,bold=TRUE)
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
    'Bonferroni and supplementary BH adjustments retain all roster sites; untestable p-values contribute 1. Raw all-episode Fisher p-values are not used.',
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
  screen$fisher_p_bonferroni<-p.adjust(ifelse(is.finite(screen$fisher_p),screen$fisher_p,1),method="bonferroni")
  nominal<-screen[is.finite(screen$fisher_p)&screen$fisher_p<.05,,drop=FALSE];nominal<-nominal[order(nominal$fisher_p),]
  corrected<-sum(screen$fisher_p_bonferroni<.05,na.rm=TRUE)
  addp(paste0('For the site-level screen, each participant contributes only their first eligible episode during higher-positivity periods. This avoids treating repeated episodes from the same person as independent observations in the Fisher exact tests. ',nrow(nominal),' of ',M,' sites had a nominal two-sided p-value below 0.05. After Bonferroni adjustment across all ',M,' sites, ',if(corrected==0)'none remained below 0.05' else paste(corrected,'remained below 0.05'),'. The individual signals therefore identify records for further review rather than establish site-specific differential treatment.'))
  body<-c(body,'__V20_SITE_PLOT__');plain<-c(plain,'[V20: Site recruitment and relative testing coverage]','')
  if(nrow(nominal)) {
    rows<-cbind(region(nominal$SITEID),mapply(rate,nominal$v_y,nominal$v_n),mapply(rate,nominal$p_y,nominal$p_n),vapply(nominal$fisher_p,pf,character(1)),num(nominal$fisher_p_bonferroni,3))
    addtable('Sites with nominal testing imbalances during higher positivity',c('Site / county','mRNA-1273','Placebo','Fisher p','Bonferroni p'),rows,paste0('First eligible higher-positivity episode per participant. Multiplicity correction includes all ',M,' sites. Small denominators, including US313, limit interpretation.'))
  }
  figure('coverage',paste0('Documented testing at ',region(sel$SITEID[sel$anomaly=='coverage']),', using each participant\'s first eligible episode during higher-positivity periods.'))
  cs<-sel$SITEID[sel$anomaly=='pcr'];within<-cases[!is.na(cases$within_same_episode_span)&cases$within_same_episode_span & is.finite(cases$days_to_later_pcr),]
  days<-sort(unique(within$days_to_later_pcr));nd<-sum(cases$not_done_in_horizon,na.rm=TRUE)
  cp<-paste0('At ',region(cs),', ',nrow(cases),' vaccine-arm episodes beginning on COVID symptom forms during higher-positivity periods did not meet the visit-plus-PCR endpoint within ',horizon,' days.')
  compare_ep<-rd(file.path(figdir,"selected_site_episodes.csv")); compare_ep<-compare_ep[compare_ep$SITEID==cs & compare_ep$high==1 & compare_ep$onset_source=="COVID_FORM",]
  covered<-function(a){u<-compare_ep[compare_ep$arm==a,];rate(sum(u$tested),nrow(u))}
  cp<-paste0(cp," Under the same selection the endpoint was met in ",covered("mRNA-1273")," vaccine episodes and ",covered("Placebo")," placebo episodes.")
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
  phmpt_export_audit_assets(audit_root,output_dir,data_dir,audit_tables)
  vp<-file.path(output_dir,"audit_assets","V20_site_recruitment_testing_ratio.png")
  body[body=="__V20_SITE_PLOT__"]<-paste0('<figure><img src="data:image/png;base64,',b64(vp),'" alt="Site recruitment and vaccine-to-placebo testing rate ratio"><figcaption>',v17_caption(),'</figcaption></figure>')
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
  message('Copy-ready report section: ',normalizePath(file.path(output_dir,'report_insert.html'),winslash='/'))
  invisible(list(html=file.path(output_dir,'report_insert.html'),paragraph_start=start_paragraph,paragraph_end=number-1L,adjusted_models=!is.null(ma)&&cut==10))
}


# ---- RUN: selecting/running the whole script executes the analysis ----
# Audit PNG export extension v1.0.0 | 2026-10-23
# Standalone; stable Vxx filenames in <project>/audit_assets plus this run's copy.
# Long tables continue as _02.png, _03.png, etc. References remain in HTML/provenance.
PHMPT_AUDIT_ASSET_DIR <- getOption("phmpt.audit_asset_dir", "audit_assets")
PHMPT_AUDIT_DPI <- getOption("phmpt.audit_dpi", 300L)
PHMPT_ROOT <- getOption("phmpt.root", "")

# ---- Standalone audit-asset renderer (V10-V39), version 1.0.0 ----
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
      if(any(dir.exists(audit_release_dir(p,c("prod_data","eua_data")))))return(p)
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
  manifest<-do.call(rbind,a$rows);manifest<-manifest[order(manifest$visual,manifest$page),,drop=FALSE];manifest$renderer_version<-"1.0.0";manifest$generated_at<-format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z")
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

# New V20: same participant-level high-positivity screen as the following table.
v17_site_data <- function(data_dir,out) {
  roster<-ax_read(data_dir,"subject_followup.csv")
  screen<-ax_read(file.path(out,"anomaly_examples"),"first_high_episode_site_screen.csv")
  map<-ax_read(data_dir,"site_region_coverage.csv")
  needed<-function(d,n)if(!all(n%in%names(d)))stop("V20 missing fields: ",paste(setdiff(n,names(d)),collapse=", "))
  needed(roster,c("subject","SITEID","RANDFL","arm"));needed(screen,c("SITEID","v_y","v_n","p_y","p_n"))
  clean<-function(x)trimws(as.character(x))
  roster<-roster[!is.na(roster$RANDFL)&toupper(clean(roster$RANDFL))=="Y",,drop=FALSE]
  roster$SITEID<-clean(roster$SITEID);roster$subject<-clean(roster$subject)
  if(!nrow(roster)||anyNA(roster[c("subject","SITEID","arm")])||any(!nzchar(roster$subject)|!nzchar(roster$SITEID))||
     anyDuplicated(roster$subject)||any(!roster$arm%in%c("mRNA-1273","Placebo")))stop("V20: invalid or duplicated randomized roster.")
  sites<-sort(unique(roster$SITEID));screen$SITEID<-clean(screen$SITEID)
  if(anyDuplicated(screen$SITEID)||!setequal(sites,screen$SITEID))stop("V20: site screen and randomized roster do not cover identical sites.")
  m<-screen[match(sites,screen$SITEID),c("SITEID","v_y","v_n","p_y","p_n")]
  if(anyNA(m)||any(vapply(m[-1],function(x)!is.numeric(x)||any(!is.finite(x)|x<0|x!=floor(x)),logical(1)))||
     any(m$v_y>m$v_n|m$p_y>m$p_n))stop("V20: invalid arm counts in site screen.")
  n<-function(arm)tabulate(match(roster$SITEID[roster$arm==arm],sites),length(sites))
  m$vaccine_randomized<-n("mRNA-1273");m$placebo_randomized<-n("Placebo")
  m$randomized_total<-m$vaccine_randomized+m$placebo_randomized
  if(any(m$v_n>m$vaccine_randomized|m$p_n>m$placebo_randomized) || sum(m$randomized_total)!=nrow(roster))stop("V20: participant counts exceed the roster.")
  m$region_name<-map$region_name[match(m$SITEID,map$SITEID)]
  m$vaccine_testing_rate<-ifelse(m$v_n>0,m$v_y/m$v_n,NA_real_)
  m$placebo_testing_rate<-ifelse(m$p_n>0,m$p_y/m$p_n,NA_real_)
  m$ratio_status<-ifelse(m$v_n==0&m$p_n==0,"NO_ELIGIBLE_EITHER_ARM",ifelse(m$v_n==0,"NO_ELIGIBLE_VACCINE",
    ifelse(m$p_n==0,"NO_ELIGIBLE_PLACEBO",ifelse(m$p_y==0 & m$v_y==0,"BOTH_RATES_ZERO",ifelse(m$p_y==0,"POSITIVE_OVER_ZERO","FINITE")))))
  m$testing_rate_ratio<-m$vaccine_testing_rate/m$placebo_testing_rate
  m$testing_rate_ratio[is.nan(m$testing_rate_ratio)]<-NA_real_
  m
}
# Fix each site's arm denominators and its pooled number with a qualifying test.
# Enumerate attainable tables, using precisely the two-sided test used for V21.
v17_fisher_ranges <- function(m,alpha=.05,adjustment="nominal") {
  adjustment<-match.arg(adjustment,c("nominal","bonferroni"))
  M<-nrow(m)
  if(!M||anyDuplicated(m$SITEID))stop("V20 requires one row for each roster site.")
  stopifnot(length(alpha)==1,is.finite(alpha),alpha>0,alpha<1)
  m$fisher_status<-rep("TESTABLE",M);m$fisher_p<-m$bonferroni_p<-NA_real_
  m$nominal_low<-m$nominal_high<-m$bonferroni_low<-m$bonferroni_high<-NA_real_
  tables<-vector("list",M)
  for(i in seq_len(M)) {
    nv<-m$v_n[i];np<-m$p_n[i];T<-m$v_y[i]+m$p_y[i]
    if(!nv||!np){m$fisher_status[i]<-"EMPTY_ARM";next}
    if(T==0L||T==nv+np){m$fisher_status[i]<-"NO_OUTCOME_VARIATION";next}
    k<-seq.int(max(0L,T-np),min(nv,T))
    f<-vapply(k,function(j)stats::fisher.test(matrix(c(j,nv-j,T-j,np-(T-j)),nrow=2L,byrow=TRUE),
                                                alternative="two.sided")$p.value,numeric(1))
    rr<-(k/nv)/((T-k)/np)
    o<-match(m$v_y[i],k);if(is.na(o))stop("V20 observed table outside support.")
    m$fisher_p[i]<-f[o];m$bonferroni_p[i]<-min(1,M*f[o])
    for(method in c("nominal","bonferroni")) {
      accepted<-if(method=="nominal")f>=alpha else pmin(1,M*f)>=alpha
      ii<-which(accepted)
      if(!length(ii)||length(ii)!=max(ii)-min(ii)+1L)stop("V20: noncontiguous Fisher acceptance set.")
      m[[paste0(method,"_low")]][i]<-rr[min(ii)]
      m[[paste0(method,"_high")]][i]<-rr[max(ii)]
    }
    tables[[i]]<-data.frame(SITEID=m$SITEID[i],vaccine_tested=k,placebo_tested=T-k,
      vaccine_eligible=nv,placebo_eligible=np,testing_rate_ratio=rr,fisher_p=f,
      bonferroni_p=pmin(1,M*f),in_nominal_range=f>=alpha,
      in_bonferroni_range=pmin(1,M*f)>=alpha,is_observed=k==m$v_y[i])
  }
  m$reference_method<-adjustment;m$alpha<-alpha;m$family_sites<-M
  m$reference_low<-m[[paste0(adjustment,"_low")]]
  m$reference_high<-m[[paste0(adjustment,"_high")]]
  m$display_p<-if(adjustment=="nominal")m$fisher_p else m$bonferroni_p
  m$outside_reference<-ifelse(is.finite(m$display_p),m$display_p<alpha,NA)
  # Test decision must match membership for the actual attainable ratio,
  # including +Inf in a positive/zero cell. NA is never a negative finding.
  j<-which(m$fisher_status=="TESTABLE");actual<-m$testing_rate_ratio[j]
  stopifnot(!anyNA(actual),all(m$outside_reference[j]==
    (actual<m$reference_low[j]|actual>m$reference_high[j])))
  all_tables<-if(length(j))do.call(rbind,tables[j]) else data.frame(SITEID=character(),vaccine_tested=integer(),placebo_tested=integer(),
    vaccine_eligible=integer(),placebo_eligible=integer(),testing_rate_ratio=numeric(),fisher_p=numeric(),bonferroni_p=numeric(),
    in_nominal_range=logical(),in_bonferroni_range=logical(),is_observed=logical())
  list(data=m,tables=all_tables)
}
# Descriptive local-linear fit to the observed (untransformed) ratio. Use the
# true randomized counts, never the slight x offsets used to draw tied sites.
v17_loess_fit <- function(m,span=getOption("phmpt.v17_loess_span",0.40)) {
  if(!is.numeric(span)||length(span)!=1L||!is.finite(span)||span<=0||span>1)
    stop("phmpt.v17_loess_span must be one number in (0,1]; default 0.40.")
  use<-is.finite(m$randomized_total)&is.finite(m$testing_rate_ratio)
  z<-data.frame(enrolled=m$randomized_total[use],ratio=m$testing_rate_ratio[use])
  empty<-data.frame(enrolled=numeric(),fitted_ratio=numeric())
  result<-function(curve=empty,status="UNAVAILABLE",reason="")list(curve=curve,
    info=data.frame(status=status,reason=reason,span=span,degree=1L,finite_sites=nrow(z),
      excluded_sites=nrow(m)-nrow(z),distinct_enrollments=length(unique(z$enrolled)),
      x_min=if(nrow(z))min(z$enrolled) else NA_real_,x_max=if(nrow(z))max(z$enrolled) else NA_real_,
      weight_basis="One observation per site; no enrollment or precision weighting",stringsAsFactors=FALSE))
  if(length(unique(z$enrolled))<3L||floor(span*nrow(z))<3L)
    return(result(reason="Too few finite site ratios or distinct recruitment counts for this span."))
  z<-z[order(z$enrolled,z$ratio),,drop=FALSE]
  # Include each observed enrollment on the grid and do not extrapolate.
  gx<-sort(unique(c(z$enrolled,seq(min(z$enrolled),max(z$enrolled),length.out=301))))
  pred<-tryCatch(withCallingHandlers({
    fit<-stats::loess(ratio~enrolled,data=z,span=span,degree=1L,family="gaussian",
      control=stats::loess.control(surface="direct",statistics="none"))
    as.numeric(stats::predict(fit,newdata=data.frame(enrolled=gx)))
  },warning=function(w)stop(conditionMessage(w))),error=function(e)e)
  if(inherits(pred,"error"))return(result(reason=conditionMessage(pred)))
  if(length(pred)!=length(gx)||any(!is.finite(pred)))return(result(reason="Non-finite fitted values."))
  result(data.frame(enrolled=gx,fitted_ratio=pred),status="AVAILABLE")
}

v17_plot_fisher <- function(a,data_dir,out,window,cut,horizon,adjustment="nominal",tables_only=FALSE,m0=NULL) {
  ran<-v17_fisher_ranges(if(is.null(m0))v17_site_data(data_dir,out) else m0,adjustment=adjustment)
  m<-ran$data;ok<-m$ratio_status=="FINITE";testable<-m$fisher_status=="TESTABLE"
  flagged<-!is.na(m$outside_reference)&m$outside_reference
  stopifnot(all(ok==is.finite(m$testing_rate_ratio)))
  m$analysis_window<-window;m$minimum_county_positivity_pct<-cut;m$pcr_horizon_days<-horizon
  # Distinct recruitment counts stay at their exact x coordinates. Only sites
  # with identical enrollment share a horizontal slot (very slight separation).
  x<-m$randomized_total;ux<-sort(unique(x));n<-length(ux);span<-max(1,max(x))
  site_x<-as.numeric(x);half<-numeric(length(x))
  for(j in seq_along(ux)){
    ii<-which(x==ux[j]);near<-min(c(Inf,if(j>1)ux[j]-ux[j-1],if(j<n)ux[j+1]-ux[j]))
    w<-min(.0045*span,.48*near);edges<-seq(ux[j]-w,ux[j]+w,length.out=length(ii)+1L)
    site_x[ii]<-(head(edges,-1)+tail(edges,-1))/2;half[ii]<-.93*w/length(ii)
  }
  m$plot_x<-site_x
  trend<-v17_loess_fit(m);curve<-trend$curve
  if(trend$info$status!="AVAILABLE")message("V20 LOESS omitted: ",trend$info$reason)
  ax_data(a,"V20",curve,"loess_curve");ax_data(a,"V20",trend$info,"loess_method")
  top<-max(pretty(c(0,1,m$testing_rate_ratio[ok],curve$fitted_ratio),n=5));cap<-1.055*top
  m$above_plot_limit<-ifelse(testable,m$reference_high>cap,NA)
  csv<-m;names(csv)[match(c("v_y","v_n","p_y","p_n"),names(csv))]<-c("vaccine_tested","vaccine_eligible","placebo_tested","placebo_eligible")
  ax_data(a,"V20",csv);ax_data(a,"V20",ran$tables,"attainable_tables")
  files<-c(file.path(data_dir,"subject_followup.csv"),file.path(out,"anomaly_examples","first_high_episode_site_screen.csv"))
  ax_data(a,"V20",data.frame(input=basename(files),md5=unname(tools::md5sum(files))),"inputs")
  if(tables_only)return(invisible(NULL))
  mode<-m$reference_method[1];method<-if(mode=="nominal")"Nominal two-sided Fisher test" else paste0("Bonferroni adjustment over ",nrow(m)," sites")
  fill<-"#D2E7EF";edge<-"#90B6C5";navy<-"#285867";red<-"#A13D25"
  trend_col<-"#123D78"
  draw<-function(){
    old<-par(no.readonly=TRUE);on.exit(par(old))
    par(fig=c(0,1,.26,.76),mai=c(.49,.82,.08,.24),family="sans",las=1,bty="n",mgp=c(2.45,.48,0))
    nd<--top*.14;tick<-pretty(c(0,top),n=5);tick<-tick[tick>=0&tick<=top]
    plot(range(c(0,x))*c(1,1.04),c(nd*1.6,cap*1.035),type="n",axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
    u<-par("usr");rect(u[1],nd*1.55,u[2],nd*.45,col="#F3F5F6",border=NA)
    for(i in which(testable)){
      lo<-m$reference_low[i];hi<-m$reference_high[i]
      if(lo<=cap){
        rect(site_x[i]-half[i],lo,site_x[i]+half[i],min(hi,cap),col=fill,border=NA)
        segments(site_x[i],lo,site_x[i],min(hi,cap),col=edge,lwd=.65)
        segments(site_x[i]-half[i],lo,site_x[i]+half[i],lo,col=edge,lwd=.65)
        if(hi<=cap)segments(site_x[i]-half[i],hi,site_x[i]+half[i],hi,col=edge,lwd=.65)
        else points(site_x[i],cap,pch=24,cex=.45,bg=edge,col=edge)
      }
    }
    abline(h=tick,col=grDevices::adjustcolor("#A8BBC4",alpha.f=.36),lwd=.6)
    abline(h=1,col="#7B6046",lwd=1.25,lty=2)
    axis(1,tck=-.012,col="#A8B4BA",col.axis="#354854",cex.axis=.78);axis(2,at=tick,labels=format(tick,trim=TRUE),tick=FALSE,cex.axis=.78)
    axis(2,at=nd,labels="N/A",tick=FALSE,cex.axis=.72,col.axis="#607078")
    mtext("Randomized participants at the site (both arms)",1,line=1.96,cex=.86)
    mtext("Testing rate ratio: mRNA-1273 / placebo",2,line=2.38,las=0,cex=.84)
    if(nrow(curve)){
      # The gray N/A strip has no numeric y values; keep the curve out of it.
      clip(u[1],u[2],0,u[4])
      lines(curve$enrolled,curve$fitted_ratio,col="white",lwd=3.6)
      lines(curve$enrolled,curve$fitted_ratio,col=trend_col,lwd=2.15)
      clip(u[1],u[2],u[3],u[4])
    }
    points(site_x[ok&testable&!flagged],m$testing_rate_ratio[ok&testable&!flagged],pch=21,bg="#477F8A",col=navy,cex=.94,lwd=.72)
    points(site_x[ok&!testable],m$testing_rate_ratio[ok&!testable],pch=21,bg="white",col="#65717A",cex=.94)
    points(site_x[ok&flagged],m$testing_rate_ratio[ok&flagged],pch=23,bg="#CA6242",col=red,cex=1.15,lwd=1.1)
    points(site_x[!ok],rep(nd,sum(!ok)),pch=5,col=ifelse(flagged[!ok],red,"#68767E"),cex=.82)
    # The five nominal signals are labelled when present; no labels for null sites.
    for(j in which(ok&flagged)){
      lab<-m$SITEID[j];xj<-site_x[j];yj<-m$testing_rate_ratio[j]
      right<-xj<.87*max(x);tx<-xj+if(right)span*.013 else -span*.013
      yl<-if(yj==0)yj+.035*top else yj
      w<-strwidth(lab,cex=.70);ht<-strheight(lab,cex=.70)
      rect(if(right)tx-1 else tx-w-1,yl-ht*.6,if(right)tx+w+1 else tx+1,yl+ht*.6,col="white",border=NA)
      text(tx,yl,lab,adj=c(if(right)0 else 1,.5),cex=.70,col=red)
    }
    grid::pushViewport(grid::viewport());on.exit(grid::popViewport(),add=TRUE)
    ax_text("V20  Site recruitment and relative testing coverage",.35,.20,17,TRUE,col="#17384B")
    ax_text(paste0(window," follow-up | County positivity >= ",cut,"% | ",sum(ok)," of ",nrow(m)," ratios finite"),.35,.61,11.1)
    sw<-function(x,y,colour,label){grid::grid.rect(x=grid::unit(x+.065,"in"),y=grid::unit(1,"npc")-grid::unit(y+.065,"in"),
      width=grid::unit(.13,"in"),height=grid::unit(.13,"in"),gp=grid::gpar(fill=colour,col=NA));ax_text(label,x+.20,y-.025,10.2)}
    sw(.35,1.04,fill,"Within-site range: p >= 0.05");sw(4.45,1.04,"#CA6242","Observed outside range: p < 0.05")
    ax_text(ax_wrap(paste0(method," | ",sum(flagged)," sites flagged"),3.95,10.3),.35,1.40,10.3,col="#435760")
    if(nrow(curve)){
      grid::grid.lines(x=grid::unit(c(4.45,4.75),"in"),y=grid::unit(1,"npc")-grid::unit(rep(1.47,2),"in"),
        gp=grid::gpar(col=trend_col,lwd=2.15))
      ax_text(sprintf("LOESS: span %.2f; equal site weights",trend$info$span),4.84,1.40,10.2,col=trend_col)
    }else ax_text("LOESS unavailable; see method CSV",4.45,1.40,10.2,col="#435760")
    ax_text(ax_wrap(paste0("Each blue strip spans the attainable testing ratios for one site that do not cross the stated test threshold. It holds fixed both arms' eligible counts and the site's total tested. Triangles mark ranges extending above the axis. Tied enrollment counts are separated slightly."),7.8,10.15),.35,5.82,10.15,col="#354854")
    ax_text(ax_wrap(paste0("Dots use the first eligible symptomatic episode per participant and a visit plus qualifying PCR within ",horizon," days. The dashed line is equal coverage. Hollow dots have no outcome variation; N/A diamonds have no finite ratio. ",
      sum(!is.na(m$bonferroni_p)&m$bonferroni_p<.05)," sites meet the Bonferroni criterion across ",nrow(m)," sites."),7.8,10.15),.35,6.57,10.15,col="#354854")
  }
  ax_png(a,"V20","site_recruitment_testing_ratio",draw,7.50,rows=nrow(m))
}

# Null risk in each arm is p0. Pool participants (sum(tested)/sum(denominators)),
# including sites without a finite RR. The log-scale SE is approximated under
# this null, not from log(observed RR); a zero count needs no correction.
v17_funnel <- function(m, population, p0=NULL) {
  population<-match.arg(population,c("randomized","eligible"))
  nv<-if(population=="randomized")m$vaccine_randomized else m$v_n
  np<-if(population=="randomized")m$placebo_randomized else m$p_n
  yv<-m$v_y;yp<-m$p_y; N<-nv+np
  if(!nrow(m)||anyNA(c(nv,np,yv,yp))||any(yv>nv|yp>np|yv<0|yp<0))stop("Invalid V20 counts.")
  pv<-if(sum(nv)>0)sum(yv)/sum(nv) else NA_real_
  pp<-if(sum(np)>0)sum(yp)/sum(np) else NA_real_
  source<-if(is.null(p0))"pooled_placebo" else "specified"
  if(is.null(p0))p0<-pp
  if(length(p0)!=1L||!is.numeric(p0)||!is.finite(p0)||p0<=0||p0>=1)
    stop("V20: a benchmark strictly between 0 and 1 is required for ",population,".")
  use<-nv>0 & np>0;ne<-rep(NA_real_,nrow(m));ne[use]<-4/(1/nv[use]+1/np[use])
  # Preserve the actual site recruitment on A. On B use balanced-equivalent
  # volume so a single funnel is exact for the actual (unequal) arm sizes.
  x<-if(population=="randomized")N else ne
  RR<-rep(NA_real_,nrow(m));RR[use]<-(yv[use]/nv[use])/(yp[use]/np[use]);RR[is.nan(RR)]<-NA_real_
  status<-ifelse(!use,"EMPTY_ARM",ifelse(yp==0 & yv==0,"BOTH_ZERO",ifelse(yp==0,"INFINITE","FINITE")))
  searm<-ifelse(use,sqrt((1-p0)/p0*(1/nv+1/np)),NA_real_)
  seplot<-ifelse(is.finite(x)&x>0,2*sqrt((1-p0)/p0/x),NA_real_);z<-1.96
  classify<-function(r,lo,hi)ifelse(is.na(r)|is.na(lo)|is.na(hi),"UNDEFINED",
    ifelse(r<lo,"BELOW",ifelse(r>hi,"ABOVE","INSIDE")))
  bounds<-function(se,k)exp(k*z*se)
  sd<-data.frame(SITEID=m$SITEID,region_name=m$region_name,population=population,
    vaccine_tested=yv,vaccine_denominator=nv,placebo_tested=yp,placebo_denominator=np,
    vaccine_randomized=m$vaccine_randomized,placebo_randomized=m$placebo_randomized,
    vaccine_eligible=m$v_n,placebo_eligible=m$p_n,denominator_total=N,effective_n=ne,
    plot_volume=x,vaccine_rate=ifelse(nv>0,yv/nv,NA_real_),placebo_rate=ifelse(np>0,yp/np,NA_real_),
    RR=RR,ratio_status=status,p0=p0,se_log_rr=searm,lower95=bounds(searm,-1),upper95=bounds(searm,1),
    plot_se_log_rr=seplot,plot_lower95=bounds(seplot,-1),plot_upper95=bounds(seplot,1))
  sd$position<-classify(RR,sd$lower95,sd$upper95)
  sd$plot_position<-classify(RR,sd$plot_lower95,sd$plot_upper95)
  sd$sparse_under_null<-ifelse(use,pmin(nv,np)*min(p0,1-p0)<5,NA)
  sd$fisher_p<-vapply(seq_along(nv),function(i){
    if(!use[i]||yv[i]+yp[i]==0||yv[i]+yp[i]==N[i])return(NA_real_)
    stats::fisher.test(matrix(c(yv[i],nv[i]-yv[i],yp[i],np[i]-yp[i]),2,byrow=TRUE))$p.value
  },numeric(1))
  sd$bonferroni_p<-pmin(1,nrow(m)*sd$fisher_p)
  xx<-x[is.finite(x)&x>0]
  gx<-if(length(xx))sort(unique(c(xx,seq(min(xx),max(xx),length.out=400)))) else numeric()
  gs<-2*sqrt((1-p0)/p0/gx)
  curve<-data.frame(population=rep(population,length(gx)),plot_volume=gx,
    p0=rep(p0,length(gx)),se_log_rr=gs,lower95=bounds(gs,-1),upper95=bounds(gs,1))
  fit<-v17_loess_fit(data.frame(randomized_total=x,testing_rate_ratio=RR))
  fit$curve<-data.frame(population=rep(population,nrow(fit$curve)),plot_volume=fit$curve$enrolled,
    fitted_ratio=fit$curve$fitted_ratio)
  info<-data.frame(population=population,benchmark_source=source,p0=p0,coefficient=2*sqrt((1-p0)/p0),
    pooled_vaccine_tested=sum(yv),pooled_vaccine_denominator=sum(nv),pooled_vaccine_rate=pv,
    pooled_placebo_tested=sum(yp),pooled_placebo_denominator=sum(np),pooled_placebo_rate=pp,
    sites=nrow(m),finite=sum(is.finite(RR)),infinite=sum(is.infinite(RR)),undefined=sum(is.na(RR)),
    below=sum(sd$position=="BELOW"),above=sum(sd$position=="ABOVE"),above_infinite=sum(is.infinite(RR)&sd$position=="ABOVE"),
    plot_below=sum(sd$plot_position=="BELOW"),plot_above=sum(sd$plot_position=="ABOVE"),
    plot_vs_actual_arm_disagreements=sum(sd$position!=sd$plot_position),
    sparse_defined=sum(sd$sparse_under_null & !is.na(RR),na.rm=TRUE),
    fisher_nominal=sum(sd$fisher_p<.05,na.rm=TRUE),fisher_bonferroni=sum(sd$bonferroni_p<.05,na.rm=TRUE),
    phi=1,z=z,loess_span=fit$info$span,loess_status=fit$info$status,loess_reason=fit$info$reason,
    x_definition=if(population=="randomized")"actual total randomized" else "4/(1/n_vaccine + 1/n_placebo)")
  list(sites=sd,summary=info,curve=curve,loess=fit$curve)
}

v17_caption <- function(){
  paste0("V20. Testing during higher county positivity. Panel A uses all randomized participants; panel B uses only participants with an eligible symptomatic episode. The numerator is identical: a documented COVID visit plus qualifying PCR on the first eligible episode. The shaded areas are approximate, pointwise 95% limits under equal arm rates; each panel identifies its placebo benchmark. The dark blue LOESS lines describe the finite site ratios. The exact conditional site comparisons follow in V21.")
}

v17_plot <- function(a,data_dir,out,window,cut,horizon) {
  # The numbered V20 is always the two-panel funnel. Exact tests are reported in V21.
  v17_check_blinded_source(data_dir)
  # Retain legacy conditional Fisher tables and the original enrollment-based
  # LOESS CSVs so rerunning the launcher preserves existing downstream links.
  m<-v17_site_data(data_dir,out)
  v17_plot_fisher(a,data_dir,out,window,cut,horizon,adjustment="nominal",tables_only=TRUE,m0=m)
  panels<-lapply(c("randomized","eligible"),function(pop) v17_funnel(m,pop,getOption(paste0("phmpt.v17_p0_",pop))))
  tab<-function(k)do.call(rbind,lapply(panels,'[[',k))
  s<-tab("summary");s$window<-window;s$minimum_county_positivity_pct<-cut;s$onset_to_pcr_days<-horizon
  ax_data(a,"V20",tab("sites"),"funnel_sites");ax_data(a,"V20",s,"funnel_summary")
  ax_data(a,"V20",tab("curve"),"funnel_limits");ax_data(a,"V20",tab("loess"),"funnel_loess")
  draw<-function(){
    op<-par(no.readonly=TRUE);on.exit(par(op))
    # Grid text is placed in inches, measured from the top; both panels use
    # identical RR ticks and a separate marker for positive/zero comparisons.
    put_text<-ax_text;text_queue<-list();ax_text<-function(...)text_queue[[length(text_queue)+1L]]<<-list(...)
    H<-9.70;cap<-max(2.5,ceiling(2*max(c(1,unlist(lapply(panels,function(z)
      z$sites$RR[is.finite(z$sites$RR)])))))/2)
    dark<-"#123D78";ink<-"#17384B";rust<-"#9E402A";dot<-"#286577";fill<-"#DCEBF0";edge<-"#7294A0"
    ax_text("V20  Relative testing during weeks of higher positivity",.40,.23,18,TRUE,col=ink)
    ax_text(paste0(window," follow-up   |   County positivity ",intToUtf8(0x2265)," ",cut,"%   |   One first eligible episode per participant"),.40,.65,11)
    ax_text("Shading: pointwise 95% limits       Dashed: equal testing (RR = 1)       Dark blue: LOESS",.40,.99,11,col=dark)
    for(j in seq_along(panels)){
      q<-panels[[j]];dd<-q$sites;su<-q$summary;cu<-q$curve;lo<-q$loess
      y<-if(j==1)1.48 else 4.92
      ax_text(if(j==1)"A   Among all randomized participants" else "B   Among participants with an eligible symptomatic episode",.40,y,14,TRUE,col=ink)
      pref<-if(su$benchmark_source=="specified")"Specified benchmark" else "Pooled placebo"
      label<-sprintf("%s: %.2f%%",pref,100*su$p0)
      if(su$benchmark_source=="pooled_placebo")label<-paste0(label,sprintf(" (%s/%s)",format(su$pooled_placebo_tested,big.mark=",",trim=TRUE),format(su$pooled_placebo_denominator,big.mark=",",trim=TRUE)))
      ax_text(sprintf("%s   |   %d below, %d above the funnel   |   %d undefined",label,su$plot_below,su$plot_above,su$undefined),.40,y+.34,10.3)
      xx<-cu$plot_volume
      xr<-if(length(xx))c(0,max(1,1.04*max(xx))) else c(0,1)
      par(fig=c(0,1,(H-(y+2.96))/H,(H-(y+.67))/H),mai=c(.45,.83,.05,.32),
        new=(j!=1),bty="n",las=1,mgp=c(2.1,.45,0),family="sans")
      plot(xr,c(-.09,cap*1.14),type="n",xaxs="i",yaxs="i",axes=FALSE,xlab="",ylab="")
      if(length(xx)){
        polygon(c(xx,rev(xx)),c(cu$lower95,rev(pmin(cap,cu$upper95))),col=fill,border=NA)
        lines(xx,cu$lower95,col=edge,lwd=.85)
        hi<-cu$upper95;hi[hi>cap]<-NA_real_;lines(xx,hi,col=edge,lwd=.85)
      }
      tk<-pretty(c(0,cap),n=5);tk<-tk[tk>=0&tk<=cap]
      abline(h=tk,col=grDevices::adjustcolor("#92A9B2",.23),lwd=.6)
      abline(h=1,col="#655644",lwd=1.15,lty=2)
      axis(2,at=tk,tick=FALSE,cex.axis=.76);axis(1,tck=-.024,cex.axis=.78,col="#9CB0B8")
      mtext("Testing ratio: mRNA-1273 / placebo",2,line=2.25,las=0,cex=.83)
      mtext(if(j==1)"Participants randomized at site (both arms)" else "Eligible participants at site (effective total for unequal arm sizes)",1,line=1.50,cex=.78)
      if(nrow(lo)){lines(lo$plot_volume,lo$fitted_ratio,col="white",lwd=3.5);lines(lo$plot_volume,lo$fitted_ratio,col=dark,lwd=2.2)}
      ok<-is.finite(dd$RR);inside<-dd$plot_position=="INSIDE"
      points(dd$plot_volume[ok&inside],dd$RR[ok&inside],pch=21,col=dot,bg="#70A0AA",cex=.80,lwd=.75)
      points(dd$plot_volume[ok&!inside],dd$RR[ok&!inside],pch=23,col=rust,bg="#CC7555",cex=.95,lwd=.85)
      inf<-is.infinite(dd$RR)
      points(dd$plot_volume[inf],rep(cap*1.065,sum(inf)),pch=24,col=rust,bg="#CC7555",cex=.95)
      ax_text(sprintf("Positive / zero (triangles): %d. Undefined ratios are omitted. LOESS span: %.2f%s.",
        su$infinite,su$loess_span,if(su$loess_status=="AVAILABLE")"" else " (not estimable)"),.41,y+3.20,9.3,col="#4A616D")
    }
    ax_text(ax_wrap(paste0("Counted as tested: a documented COVID visit and qualifying PCR within ",horizon,
      " days of symptom onset, on the first eligible episode in a higher-positivity week. Panel A also reflects which participants had a recorded eligible episode. Panel B compares testing conditional on that episode."),7.72,10.35),.4,8.40,10.35)
    ax_text(ax_wrap("Diamonds fall outside the approximate 95% limits. These are screening signals, particularly with sparse counts or zero events. Exact Fisher results for the eligible participants are reported in V21.",7.72,10.15),.4,9.00,10.15,col="#4A616D")
    for(args in text_queue)do.call(put_text,args)
  }
  ax_png(a,"V20","site_recruitment_testing_ratio",draw,9.70,rows=nrow(m))
}

phmpt_export_audit_assets <- function(root,out,data_dir,audit_tables) {
  a<-ax_begin(root,out,"V18_V27_episodes")
  fp<-ax_read(file.path(out,"anomaly_examples"),"figure_parameters.csv")
  setting<-function(key,default){z<-fp$value[fp$parameter==key];if(length(z)==1)as.character(z) else default}
  cut<-as.numeric(setting("high_positivity_threshold","10"));horizon<-as.numeric(setting("onset_to_pcr_days","7"))
  window<-setting("window","UNKNOWN")
  high_label<-paste0("High positivity >=",cut,"% | Low <",cut,"%")
  gettab<-function(name)if(name%in%names(audit_tables))audit_tables[[name]] else NULL
  specs<-list(
    V18=c("Documented testing after symptomatic episodes","episode_testing_overview"),
    V19=c("Adjusted comparisons during higher background positivity","adjusted_testing_models"),
    V21=c("Sites with nominal testing imbalances during higher positivity","site_testing_screen"),
    V26=c("Sites with opposing changes across positivity periods","site_positivity_interactions"))
  for(id in names(specs)) {
    if(id=="V21")v17_plot(a,data_dir,out,window,cut,horizon)
    s<-specs[[id]];d<-gettab(s[1])
    if(is.null(d)){ax_unavailable(a,id,s[2],s[1],if(id=="V19")"No converged, input-hash-matched adjusted models. Run analyse_site_testing_imbalance.R on this export." else "No matching comparison in this run.");next}
    if(id=="V26") {
      d<-do.call(rbind,lapply(seq_len(nrow(d)),function(i)data.frame(Site=d[i,1],Arm=c("mRNA-1273","Placebo"),
        `Low positivity`=c(d[i,2],d[i,4]),`High positivity`=c(d[i,3],d[i,5]),`Interaction p`=c(d[i,6],""),check.names=FALSE)))
    }
    note<-switch(id,V18="Endpoint: documented COVID visit plus qualifying PCR in the configured linkage window. Repeated episodes can contribute. Differences are vaccine minus placebo; missing background remains in overall totals.",
      V19="Site-clustered uncertainty. Odds ratios compare testing in vaccine versus placebo episodes. Onset recording source is a sensitivity adjustment. Early post-dose symptoms remain in the main comparison.",
      V21="First eligible high-positivity episode per participant. Fisher tests are two-sided; Bonferroni p uses the complete site family. Nominal signals select records for review.",
      V26="Low/high are pooled positivity categories, not consecutive weeks. Every arm-period cell has at least ten episodes. Interaction p values account for repeated participants and are nominal.")
    ax_table(a,id,s[2],s[1],d,paste0("COVID forms plus matched AEs | ",window," | ",high_label),note,
      widths=switch(id,V18=c(.30,.25,.25,.20),V19=c(.45,.16,.23,.16),V21=c(.30,.22,.22,.13,.13),V26=c(.15,.20,.25,.25,.15)))
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
  z<-ax_read(ex,"first_high_episode_site_screen.csv");z$fisher_p_bonferroni<-p.adjust(ifelse(is.finite(z$fisher_p),z$fisher_p,1),"bonferroni");s<-site("coverage","US372");z<-z[z$SITEID==s,]
  if(nrow(z)==1){chart("V22","site_testing_example",paste(s,"Documented testing after symptoms"),paste0("First eligible high-positivity episode per participant | ",high_label),
    paste0("Visit-plus-PCR endpoint. Two-sided Fisher p ",ax_p(z$fisher_p),"; Bonferroni p ",ax_p(z$fisher_p_bonferroni)," across the complete site family."),
    function()bars(c(z$v_y,z$p_y),c(z$v_n,z$p_n),c("mRNA-1273","Placebo")));ax_data(a,"V22",z)
  }else ax_unavailable(a,"V22","site_testing_example","Site testing example","Selected site not available.")
  z<-ax_read(ex,"pcr_episode_review.csv");s<-site("pcr","US325")
  if(nrow(z)) {
    zz<-ax_read(ex,"selected_site_episodes.csv");zz<-zz[zz$SITEID==s & zz$high==1 & zz$onset_source=="COVID_FORM",,drop=FALSE]
    armcount<-function(ar){u<-zz[zz$arm==ar,,drop=FALSE];paste0(sum(u$tested),"/",nrow(u))}
    chart("V23","delayed_pcr_example",paste(s,"PCR timing after COVID-form symptoms"),paste0("High-positivity form-onset episodes meeting the ",horizon,"-day endpoint: vaccine ",armcount("mRNA-1273"),", placebo ",armcount("Placebo"),". Vaccine-arm exceptions shown."),
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
    },height=max(6.6,3.5+.7*nrow(z)));ax_data(a,"V23",z)
  }else ax_unavailable(a,"V23","delayed_pcr_example","PCR timing review","No selected episodes.")
  z<-ax_read(ex,"recording_source_counts.csv");s<-site("source","US352")
  d<-data.frame(Arm=z$arm,`Onset source`=z$source,`Share of episodes`=paste0(z$episodes," (",sprintf("%.1f%%",z$share_pct),")"),
    `Tested within source`=ax_rate(z$tested,z$episodes),check.names=FALSE)
  # A compact table gives both the composition and the within-source denominator.
  ax_table(a,"V24","source_composition_example",paste(s,"Recording source and documented testing"),d,
    paste0("All eligible episodes | Source at first symptom date | ",high_label),
    "AE = AE records only; COVID_FORM = COVID forms only; AE; COVID_FORM = both on the first symptom date. The share column describes source composition; the final column describes visit-plus-PCR testing within that source.",widths=c(.19,.24,.27,.30))
  z<-ax_read(ex,"postdose_timing_counts.csv");s<-site("timing","US372")
  if(nrow(z)) {
    chart("V25","postdose_timing_example",paste(s,"Testing by post-dose symptom timing"),paste0("All eligible episodes | Repeated episodes retained | ",high_label),
      "Each pair shows mRNA-1273 (blue) and placebo (brown). Post-dose days 0-6 remain in the main comparison. The partition localizes the difference without identifying its cause.",function(){
      labs<-rep(c("Vaccine","Placebo"),nrow(z));k<-as.vector(rbind(z$v_y,z$p_y));n<-as.vector(rbind(z$v_n,z$p_n));x<-bars(k,n,labs)
      for(i in seq_len(nrow(z)))mtext(z$group[i],1,at=mean(x[(2*i-1):(2*i)]),line=2.2,cex=.82)
    });ax_data(a,"V25",z)
  }else ax_unavailable(a,"V25","postdose_timing_example","Testing by post-dose timing","No selected episodes.")
  z<-ax_read(ex,"site_high_low_interactions.csv");s<-site("waves","US302");z<-z[z$subset=="ALL_EPISODES"&z$SITEID==s,]
  if(nrow(z)==1) {
    chart("V27","positivity_testing_example",paste(s,"Testing in lower and higher positivity"),paste0("All eligible episodes | Pooled categories | ",high_label),
      paste0("Vaccine change minus placebo change: ",sprintf("%+.1f",z$difference_in_changes_pp)," percentage points. Participant-clustered interaction p ",ax_p(z$interaction_p)," (nominal). Lines connect pooled categories, not successive calendar weeks."),function(){
      plot(c(.65,2.35),c(0,106),type="n",axes=FALSE,xlab="",ylab="");abline(h=seq(0,100,20),col="#E1E7EB")
      axis(2,at=seq(0,100,20),labels=paste0(seq(0,100,20),"%"),tick=FALSE,cex.axis=.85);axis(1,at=1:2,labels=c("Lower positivity","Higher positivity"),tick=FALSE,cex.axis=.85)
      for(i in 1:2){p<-if(i==1)"v" else "p";k<-unlist(z[paste0(p,c("_low_y","_high_y"))],use.names=FALSE);n<-unlist(z[paste0(p,c("_low_n","_high_n"))],use.names=FALSE)
        y<-100*k/n;lines(1:2,y,col=col[i],lwd=2);points(1:2,y,pch=16,col=col[i]);other<-if(i==1)100*c(z$p_low_y/z$p_low_n,z$p_high_y/z$p_high_n) else 100*c(z$v_low_y/z$v_low_n,z$v_high_y/z$v_high_n);text(1:2,y+ifelse(y>=other,8,-8),ax_rate(k,n),col=col[i],cex=.87,font=2)}
      legend("topright",c("mRNA-1273","Placebo"),col=col,lwd=2,bty="n",cex=.86)
    });ax_data(a,"V27",z)
  }else ax_unavailable(a,"V27","positivity_testing_example","Testing across positivity periods","Selected site unavailable.")
  ax_finish(a)
}

if(isTRUE(getOption('phmpt.autorun',TRUE))) {
  args<-commandArgs(trailingOnly=TRUE)
  if(length(args)>=1L)DATA_DIR<-args[1]
  if(length(args)>=2L)OUTPUT_DIR<-args[2]
  PHMPT_SITE_AUDIT_RESULT<-analyse_site_testing_imbalance(DATA_DIR,OUTPUT_DIR,
    PHMPT_REPORT_ONLY,PHMPT_REPORT_START,PHMPT_REFERENCE_START)
}
