# Standalone: set this bundle folder as working directory and run the entire script.
# options(phmpt.followup_autorun=FALSE) loads functions only.
local({r<-getOption('phmpt.followup_root',getwd());source(file.path(r,'audit_followup_helpers.R'),local=FALSE)})
analyse_testing_robustness <- function(root=audit_root()) {
 out<-file.path(root,'results','testing');assets<-file.path(root,'audit_assets');dir.create(assets,recursive=TRUE,showWarnings=FALSE);set.seed(4301)
 e<-rd(root,'episodes');q<-rd(root,'required_evaluation_episodes');f<-rd(root,'subject_followup');ck(!anyDuplicated(e$episode_id),'Duplicate episode IDs')
 main<-e[yes(e$incident)&e$classification=='MANDATORY'&yes(e$complete3),];seven<-main[yes(main$complete7)&yes(main$orig_eligible),]
 sets<-list(All_required_evaluations=q,All_mandatory_swab_episodes=main,Complete_seven_day_opportunity=seven,No_linked_positive_PCR_D0_7=seven[text0(seven$pcr_result_d0_7)!='POSITIVE',],PP_D2plus14_complete7=seven[yes(seven$pp14),],PP_D2plus14_no_linked_positive=seven[yes(seven$pp14)&text0(seven$pcr_result_d0_7)!='POSITIVE',])
 # Same complete-7 denominator in both rows isolates the effect of positivity exclusion.
 ans<-list();for(nm in names(sets))for(m in if(nm=='All_required_evaluations')'form_by_d3'else c('swab_d0_3','swab_d0_7','any_sample_d0_7'))ans[[length(ans)+1]]<-cbind(scope=nm,endpoint=m,cluster_compare(sets[[nm]],m))
 a<-do.call(rbind,ans);wr(a,out,'complete_and_selected_comparisons')
 outcomes<-do.call(rbind,lapply(names(sets),function(nm){z<-sets[[nm]];as.data.frame(table(scope=rep(nm,nrow(z)),arm=z$arm,result=ifelse(yes(z$any_sample_d0_7),ifelse(text0(z$pcr_result_d0_7)=='','INDETERMINATE',text0(z$pcr_result_d0_7)),'UNTESTED')))}));wr(outcomes,out,'pcr_outcome_counts')
 fields<-c('site','calendar_month','trigger_family','severity','elapsed_onset_band','trigger_duration_band','recorded_group_count_band','dose_phase','days_since_dose_band')
 q$calendar_month<-format(date0(q$trigger_date),'%Y-%m');q$stratum<-do.call(key,q[fields]);common<-intersect(q$stratum[q$arm=='mRNA-1273'],q$stratum[q$arm=='Placebo']);q$supported<-q$stratum%in%common
 wr(q,out,'matched_and_unmatched_episode_ledger')
 rows<-list();for(g in c('ALL','MATCHED','UNMATCHED'))for(arm in c('mRNA-1273','Placebo')){z<-q[q$arm==arm&(g=='ALL'|q$supported==(g=='MATCHED')),];rows[[length(rows)+1]]<-data.frame(group=g,arm=arm,episodes=nrow(z),subjects=length(unique(z$subject)),sites=length(unique(z$site)),forms3=sum(yes(z$form_by_d3)),forms3_pct=100*mean(yes(z$form_by_d3)),swabs3=sum(yes(z$swab_d0_3)),swabs3_pct=100*mean(yes(z$swab_d0_3)),sample7=sum(yes(z$any_sample_d0_7)),complete7=sum(yes(z$complete7)),positive7=sum(text0(z$pcr_result_d0_7)=='POSITIVE'),untested7=sum(!yes(z$any_sample_d0_7)),ae_trigger=sum(text0(z$trigger_source)=='AE'),ae_origin=sum(text0(z$onset_source)=='AE'))}
 ms<-do.call(rbind,rows);wr(ms,out,'matched_subset_summary')
 profiles<-list();for(col in c('site','onset_source','trigger_source','trigger_family','severity','trigger_duration_band','days_since_dose_band','dose_phase')){z<-as.data.frame(table(arm=q$arm,supported=q$supported,level=text0(q[[col]])));z$variable<-col;profiles[[col]]<-z};wr(do.call(rbind,profiles),out,'matched_unmatched_profiles')
 st<-aggregate(cbind(n=rep(1,nrow(q)),forms=as.integer(yes(q$form_by_d3)),swabs=as.integer(yes(q$swab_d0_3)),positive=as.integer(text0(q$pcr_result_d0_7)=='POSITIVE')),q[c('site','arm','supported')],sum);wr(st,out,'matched_unmatched_by_site')
 c1<-q[q$supported&q$arm=='mRNA-1273',];c0<-q[q$supported&q$arm=='Placebo',]
 agg<-function(z)aggregate(cbind(n=rep(1,nrow(z)),yes=as.integer(yes(z$form_by_d3))),z[c('stratum','site')],sum)
 z<-merge(agg(c1),agg(c0),by=c('stratum','site'));z$w<-z$n.x*z$n.y/(z$n.x+z$n.y);wr(data.frame(vaccine_supported=nrow(c1),placebo_supported=nrow(c0),vaccine_support_pct=100*nrow(c1)/sum(q$arm=='mRNA-1273'),placebo_support_pct=100*nrow(c0)/sum(q$arm=='Placebo'),sites=length(unique(z$site)),strata=nrow(z),vaccine_standardized_pct=100*weighted.mean(z$yes.x/z$n.x,z$w),placebo_standardized_pct=100*weighted.mean(z$yes.y/z$n.y,z$w)),out,'matched_reproduction')
 match_result<-function(d,outcome,stage){
 d$calendar_month<-format(date0(d$trigger_date),'%Y-%m');d$stratum<-do.call(key,d[fields]);co<-intersect(d$stratum[d$arm=='mRNA-1273'],d$stratum[d$arm=='Placebo']);v<-d[d$stratum%in%co&d$arm=='mRNA-1273',];p<-d[d$stratum%in%co&d$arm=='Placebo',]
 estimate<-lo<-hi<-vp<-pp<-NA_real_;ns<-0L
 if(nrow(v)&&nrow(p)){
  ag<-function(x)aggregate(cbind(n=rep(1,nrow(x)),x=as.integer(yes(x[[outcome]]))),x[c('stratum','site')],sum)
  z<-merge(ag(v),ag(p),by=c('stratum','site'));z$w<-z$n.x*z$n.y/(z$n.x+z$n.y);vp<-100*weighted.mean(z$x.x/z$n.x,z$w);pp<-100*weighted.mean(z$x.y/z$n.y,z$w);estimate<-vp-pp
  g<-aggregate(cbind(num=z$w*(z$x.x/z$n.x-z$x.y/z$n.y),den=z$w),list(site=z$site),sum);ns<-nrow(g)
  if(ns>=3){w<-rmultinom(2000,ns,rep(1/ns,ns));ci<-quantile(100*as.vector(crossprod(g$num,w))/as.vector(crossprod(g$den,w)),c(.025,.975),names=FALSE);lo<-ci[1];hi<-ci[2]}
 }
 data.frame(stage=stage,vaccine_eligible=sum(d$arm=='mRNA-1273'),placebo_eligible=sum(d$arm=='Placebo'),vaccine_matched=nrow(v),placebo_matched=nrow(p),sites=ns,vaccine_standardized_pct=vp,placebo_standardized_pct=pp,gap_pp=estimate,lower95=lo,upper95=hi)
 }
 # Boundary scenarios truncate or exclude only; they never invent additional follow-up.
 bs<-boundaries(f);j<-match(e$subject,f$subject);ck(!anyNA(j),'Missing follow-up for an episode');jr<-match(q$subject,f$subject)
 sens<-list();lost<-list();matched_bounds<-list();for(nm in names(bs)){
 end<-bs[[nm]][j];valid3<-!is.na(end)&date0(e$trigger_date)+3<=end;valid7<-!is.na(end)&date0(e$trigger_date)+7<=end
 m<-e[yes(e$incident)&e$classification=='MANDATORY'&yes(e$complete3)&valid3,];s<-e[yes(e$incident)&e$classification=='MANDATORY'&yes(e$complete7)&yes(e$orig_eligible)&valid7,]
 vq<-!is.na(bs[[nm]][jr])&date0(q$trigger_date)+3<=bs[[nm]][jr];r<-q[vq,]
 matched_bounds[[length(matched_bounds)+1]]<-cbind(boundary=nm,match_result(r,'form_by_d3','Checklist_required_evaluations'))
 matched_bounds[[length(matched_bounds)+1]]<-cbind(boundary=nm,match_result(m[yes(m$form_by_d3),],'form_then_swab_by_d3','Swab_given_checklist'))
 for(spec in list(list(r,'form_by_d3','required_evaluation'),list(m,'swab_d0_3','mandatory'),list(s,'swab_d0_7','mandatory_complete7')))sens[[length(sens)+1]]<-cbind(boundary=nm,denominator=spec[[3]],endpoint=spec[[2]],cluster_compare(spec[[1]],spec[[2]]))
 for(arm in c('mRNA-1273','Placebo')){ff<-f$arm==arm;lost[[length(lost)+1]]<-data.frame(boundary=nm,arm=arm,subjects=sum(ff),missing_UBLDDT=sum(ff&is.na(date0(f$UBLDDT))),no_recorded_boundary=sum(ff&!yes(f$has_blind_boundary)),excluded=sum(ff&is.na(bs[[nm]])),earlier_end=sum(ff&!is.na(bs[[nm]])&bs[[nm]]<date0(f$follow_end)))}
 }
 wr(do.call(rbind,matched_bounds),out,'matched_boundary_sensitivity');bsout<-do.call(rbind,sens);wr(bsout,out,'boundary_sensitivity');wr(do.call(rbind,lost),out,'boundary_population')
 core<-a[a$scope%in%c('All_required_evaluations','All_mandatory_swab_episodes','Complete_seven_day_opportunity','No_linked_positive_PCR_D0_7') & ((a$scope=='All_required_evaluations'&a$endpoint=='form_by_d3')|(a$scope=='All_mandatory_swab_episodes'&a$endpoint=='swab_d0_3')|(a$scope%in%c('Complete_seven_day_opportunity','No_linked_positive_PCR_D0_7')&a$endpoint=='swab_d0_7')),]
 labels<-c('Checklist by day 3: all required evaluations','Swab by day 3: all mandatory swab episodes','Swab by day 7: all episodes with full follow-up','Swab by day 7: excluding linked PCR-positive episodes')
 v<-data.frame(Comparison=labels,'mRNA-1273'=rate(core$vaccine_x,core$vaccine_n),Placebo=rate(core$placebo_x,core$placebo_n),check.names=FALSE);wr(v,assets,'V47_complete_episode_comparison');print_table_png(v,file.path(assets,'V47_complete_episode_comparison.png'),'V47 All qualifying episodes and positivity-restricted sensitivity','Day 3 is a calendar-date proxy; rows 3 and 4 share the same follow-up rule.')

 v<-data.frame(Population=c('All eligible episodes: checklist by day 3','Episodes in common support (% of all eligible)','Common support: raw checklist coverage','Common support: standardized checklist coverage','Outside common support: checklist by day 3'),
 'mRNA-1273'=c(rate(ms$forms3[1],ms$episodes[1]),rate(nrow(c1),ms$episodes[1]),rate(ms$forms3[3],ms$episodes[3]),sprintf('%.1f%%',100*weighted.mean(z$yes.x/z$n.x,z$w)),rate(ms$forms3[5],ms$episodes[5])),
 'Placebo'=c(rate(ms$forms3[2],ms$episodes[2]),rate(nrow(c0),ms$episodes[2]),rate(ms$forms3[4],ms$episodes[4]),sprintf('%.1f%%',100*weighted.mean(z$yes.y/z$n.y,z$w)),rate(ms$forms3[6],ms$episodes[6])),check.names=FALSE)
 wr(v,assets,'V48_common_support');print_table_png(v,file.path(assets,'V48_common_support.png'),'V48 Checklist coverage inside and outside the matched subset','Standardization gives both arms the same distribution across 365 shared clinical/site/month strata.')

 html_report(file.path(out,'testing_results.html'),'Testing coverage, selection and common support',c('<p>Complete eligible episodes estimate documented investigation or swabbing among recorded illnesses warranting it. Excluding detected positives conditions on an event downstream of testing. That sensitivity describes the remaining recorded episodes; it cannot establish equal ascertainment. No efficacy estimate is calculated.</p>',html_table(a),'<h2>Common support</h2>',html_table(ms),'<p>Exact matching uses site, month and seven recorded clinical/timing descriptors. Missing grades remain an explicit category. Unmatched episodes have no counterpart under these rules; this is not evidence of clinical misconduct. Source and PCR result do not enter the matching variables.</p>','<h2>Blinding-boundary sensitivity</h2>',html_table(bsout),'<p>The 17 December cap is a deliberately conservative calendar scenario, not an imputed personal unblinding date. All analyses remain conditional on recorded complaints and the original trigger reconstruction.</p>'))
 invisible(a)
}
if(isTRUE(getOption('phmpt.followup_autorun',TRUE)))analyse_testing_robustness()
