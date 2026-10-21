# Moderna P301 - exploratory site testing audit | v1.2.0 | 2026-10-21
# Stand-alone base R >=4.3: select all and run in RStudio from the project root.
# Generates report_insert.html/txt, numbered prose, tables, references and five
# illustrations as PDF/PNG (and SVG with Cairo). Does not modify the draft.
# Primary high-positivity rule: county weekly NAAT positivity >=10%.
# Command line: Rscript analyse_site_testing_imbalance.R <export_folder> <output_folder>
# To load functions without running: options(phmpt.autorun=FALSE).

DATA_DIR <- getOption('phmpt.audit_data_dir','')   # blank: discover latest BLINDED heatmap export
OUTPUT_DIR <- getOption('phmpt.audit_output_dir','') # blank: <export>/positivity_audit
PHMPT_REPORT_START <- getOption('phmpt.report_start_paragraph',25L)
PHMPT_REFERENCE_START <- getOption('phmpt.report_reference_start',32L)
PHMPT_REPORT_ONLY <- getOption('phmpt.report_only',FALSE) # reuse existing, hash-matched model CSVs

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
  DATA_DIR<-phmpt_find_export(data_dir)
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
  repo_url='https://github.com/OpenVaet/moderna_trial_analysis/blob/main') {
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
  body<-character();plain<-character();number<-as.integer(start_paragraph);tables<-list();captions<-list();fig_index<-0L
  addh<-function(s){body<<-c(body,paste0('<h2>',esc(s),'</h2>'));plain<<-c(plain,s,'')}
  addp<-function(s,refs=integer()) {
    marker<-if(length(refs))paste0(' <sup>',paste0('<a href="#ref-',refs,'">',refs,'</a>',collapse=', '),'</sup>') else ''
    body<<-c(body,paste0('<p class="numbered"><span class="number">',number,'.</span> ',esc(s),marker,'</p>'))
    plain<<-c(plain,paste0(number,'. ',s,if(length(refs))paste0(' [',paste(refs,collapse=', '),']') else ''),'')
    number<<-number+1L
  }
  addtable<-function(title,headers,rows,note='') {
    rows<-as.data.frame(rows,stringsAsFactors=FALSE);names(rows)<-headers
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
  message('Copy-ready report section: ',normalizePath(file.path(output_dir,'report_insert.html'),winslash='/'))
  invisible(list(html=file.path(output_dir,'report_insert.html'),paragraph_start=start_paragraph,paragraph_end=number-1L,adjusted_models=!is.null(ma)&&cut==10))
}


# ---- RUN: selecting/running the whole script executes the analysis ----
if(isTRUE(getOption('phmpt.autorun',TRUE))) {
  args<-commandArgs(trailingOnly=TRUE)
  if(length(args)>=1L)DATA_DIR<-args[1]
  if(length(args)>=2L)OUTPUT_DIR<-args[2]
  PHMPT_SITE_AUDIT_RESULT<-analyse_site_testing_imbalance(DATA_DIR,OUTPUT_DIR,
    PHMPT_REPORT_ONLY,PHMPT_REPORT_START,PHMPT_REFERENCE_START)
}
