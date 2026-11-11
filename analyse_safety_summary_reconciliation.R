# AE -> S3 ADAE -> released CSR Tables 7-12/7-13. Base R, standalone bundle.
local({r<-getOption('phmpt.followup_root',getwd());source(file.path(r,'audit_followup_helpers.R'),local=FALSE)})
analyse_safety_summary_reconciliation <- function(root=audit_root()) {
 out<-file.path(root,'results','safety');assets<-audit_asset_path(root);dir.create(assets,recursive=TRUE,showWarnings=FALSE)
 ae<-rd(root,'S3_AE');ae27<-rd(root,'S27_AE');a<-rd(root,'S3_ADAE');sl<-rd(root,'S3_ADSL');f<-rd(root,'subject_followup');ar<-rd(root,'new_flag_case_ledger');changes<-rd(root,'seriousness_changes');set.seed(4302)
 for(nm in c('ae','ae27','a')){d<-get(nm);d$key<-key(d$SUBJID,d$AESPID);ck(!anyDuplicated(d$key),paste('Nonunique subject/form in',nm));assign(nm,d)}
 j<-match(ae$key,a$key);j27<-match(ae$key,ae27$key);ck(!anyNA(j27),'S3 AE keys not all retained in S27; inspect before proceeding')
 global<-data.frame(subject=ae$SUBJID,form=ae$AESPID,AESEQ_S3=ae$AESEQ,AESEQ_S27=ae27$AESEQ[j27],AE_term=ae$AETERM,AE_serious=ae$AESER,in_S3_ADAE=!is.na(j),ADAE_serious=a$AESER[j],ADAE_TRTEMFL=a$TRTEMFL[j],ADAE_APERIOD=a$APERIOD[j],ADAE_D28FL=a$D28FL[j],ADAE_AR_term=a$ARTERM[j]);flags<-rd(root,'S3_flags');global$S3_REMOVEFL<-key(ae$SUBJID,ae$AESEQ)%in%key(flags$subject,flags$seq);global$in_S3_ADSL<-ae$SUBJID%in%sl$SUBJID
 global$bridge_status<-ifelse(global$in_S3_ADAE,'RETAINED_IN_ADAE',ifelse(global$S3_REMOVEFL,'S3_AR_REMOVE_FLAG',ifelse(!nzchar(text0(ae$AETERM)),'EMPTY_AE_TERM',ifelse(!global$in_S3_ADSL,'NO_S3_ADSL_ROW','UNEXPLAINED_ABSENCE'))))
 wr(global,out,'global_S3_AE_ADAE_bridge');wr(as.data.frame(table(global$bridge_status)),out,'source_exclusion_categories')
 common<-!is.na(j);ck(all(text0(ae$AESER[common])==text0(a$AESER[j[common]])),'AE/ADAE seriousness mismatch; inspect bridge')
 a$partA<-text0(a$SAFFL)=='Y'&text0(a$TRTEMFL)=='Y'&text0(a$APERIOD)=='1'&nzchar(text0(a$AETERM));a$days28<-a$partA&text0(a$D28FL)=='Y'
 # Exact row selectors from released table shell/program, checked against CSR pp290-291.
 masks<-list(All=rep(TRUE,nrow(a)),Serious=text0(a$AESER)=='Y',Fatal=text0(a$AEOUT)=='FATAL',Medically_attended=text0(a$AEMAFL)=='Y',Vaccine_discontinued=grepl('WITHDRAWN',text0(a$AEACN),fixed=TRUE),Study_discontinued=text0(a$AEDISFL)=='Y',Severe=text0(a$AESEV)=='SEVERE')
 expected<-rd(root,'published_CSR_counts');counts<-list();members<-list()
 for(win in c('D28','Part_A'))for(rel in c('Any','Related'))for(cat in names(masks))for(arm in c('mRNA-1273','Placebo')){
 keep<-(if(win=='D28')a$days28 else a$partA)&masks[[cat]]&text0(a$TRT01A)==arm
 if(rel=='Related')keep<-keep&text0(a$ARELN)=='2'
 n<-length(unique(a$SUBJID[keep]));den<-sum(text0(sl$SAFFL)=='Y'&text0(sl$TRT01A)==arm);idx<-which(expected$window==win&expected$relationship==rel&expected$category==cat&expected$arm==arm);ck(length(idx)==1,'Missing/duplicate published table reference')
 counts[[length(counts)+1]]<-data.frame(window=win,relationship=rel,category=cat,arm=arm,subjects=n,denominator=den,published=expected$count[idx],published_N=expected$denominator[idx],match=n==expected$count[idx]&den==expected$denominator[idx])
 if(any(keep))members[[length(members)+1]]<-data.frame(window=win,relationship=rel,category=cat,arm=arm,subject=a$SUBJID[keep],form=a$AESPID[keep],AESEQ=a$AESEQ[keep])
 }
 tab<-do.call(rbind,counts);wr(tab,out,'CSR_summary_reproduction');wr(do.call(rbind,members),out,'table_record_contributions')
 # New removal flags were added after S3. They must not be applied retrospectively.
 ar$key<-key(ar$subject,ar$AE_spid);idx<-match(ar$key,a$key)
 ar$in_S3_ADAE<-!is.na(idx);ar$S3_partA_selected<-a$partA[idx];ar$S3_D28_selected<-a$days28[idx];ar$S3_analysis_serious<-a$AESER[idx]
 ar$only_D28_event_for_subject<-vapply(seq_len(nrow(ar)),function(i){other<-a$days28&a$SUBJID==ar$subject[i]&!a$key%in%ar$key;isTRUE(ar$S3_D28_selected[i])&&!any(other)},logical(1))
 ar$S27_final_analysis_status<-'NOT_VERIFIED_NO_S27_ADAE_OR_FINAL_TABLE_LINEAGE'
 wr(ar,out,'new_removal_flags_ADAE_bridge')
 arsum<-do.call(rbind,lapply(c('mRNA-1273','Placebo'),function(arm){z<-ar[ar$arm==arm,];data.frame(arm=arm,flagged=nrow(z),ongoing_day7=sum(yes(z$CE_ongoing_day7)),in_S3_ADAE=sum(z$in_S3_ADAE),D28_selected=sum(yes(z$S3_D28_selected)),only_D28_event_subjects=sum(z$only_D28_event_for_subject))}));wr(arsum,out,'new_flag_summary')
 # Original row indices are validated against subject, onset and term before use.
 cidx<-as.integer(changes$bla_row);ck(all(ae$SUBJID[cidx]==changes$subject)&all(gsub('[[:space:]]+',' ',trimws(toupper(ae$AETERM[cidx])))==gsub('[[:space:]]+',' ',trimws(toupper(changes$bla_entity))))&all(substr(ae$AESTDTC[cidx],1,10)==substr(changes$bla_date,1,10)),'Seriousness candidate row ordering changed')
 changes$S3_form<-ae$AESPID[cidx];changes$S27_serious<-ae27$AESER[match(ae$key[cidx],ae27$key)];k<-match(ae$key[cidx],a$key)
 changes$in_S3_ADAE<-!is.na(k);changes$S3_ADAE_serious<-a$AESER[k];changes$S3_partA_selected<-a$partA[k];changes$S3_D28_selected<-a$days28[k]
 changes$other_PartA_SAE<-vapply(seq_len(nrow(changes)),function(i)any(a$partA&a$SUBJID==changes$subject[i]&text0(a$AESER)=='Y'&a$key!=ae$key[cidx[i]]),logical(1))
 changes$cross_release_identity_status<-ifelse(changes$match_confidence=='REFERENCE_LINKED','FORM_REFERENCE_LINKED_NOT_INDEPENDENTLY_ADJUDICATED','CANDIDATE_PAIR_REQUIRES_IDENTITY_REVIEW')
 wr(changes,out,'seriousness_changes_ADAE_bridge')
 missing<-data.frame(subject=c('US3172069','US3172205','US3432392'));missing$arm<-sl$TRT01P[match(missing$subject,sl$SUBJID)];missing$in_S3_ADSL<-missing$subject%in%sl$SUBJID
 missing$S3_AE_rows<-vapply(missing$subject,function(s)sum(ae$SUBJID==s),integer(1));missing$S27_AE_rows<-vapply(missing$subject,function(s)sum(ae27$SUBJID==s),integer(1));missing$S3_ADAE_rows<-vapply(missing$subject,function(s)sum(a$SUBJID==s),integer(1));missing$unsolicited_table_contribution<-ifelse(missing$S3_ADAE_rows==0,'NO_EVENT_CONTRIBUTION_THROUGH_S3_ADAE','REVIEW');wr(missing,out,'participants_without_licensing_AE')
 # Sensitivities retain actual treatment for safety; crossover events remain Part B.
 bs<-boundaries(f);fi<-match(a$SUBJID,f$subject);si<-match(sl$SUBJID,f$subject);onset<-date0(a$ASTDT);sens<-list()
 for(bn in names(bs))for(arm in c('mRNA-1273','Placebo'))for(cat in c('All','Serious','Fatal')){
 eligible<-text0(sl$SAFFL)=='Y'&text0(sl$TRT01A)==arm&!is.na(bs[[bn]][si]);ids<-sl$SUBJID[eligible]
 keep<-a$partA&a$SUBJID%in%ids&masks[[cat]]&!is.na(onset)&!is.na(bs[[bn]][fi])&onset<=bs[[bn]][fi]
 sens[[length(sens)+1]]<-data.frame(boundary=bn,arm=arm,category=cat,subjects=length(unique(a$SUBJID[keep])),denominator=length(ids),events=sum(keep),missing_onset_in_source_subset=sum(a$partA&a$SUBJID%in%ids&masks[[cat]]&is.na(onset)),subjects_if_missing_onsets_included=length(unique(a$SUBJID[keep|(a$partA&a$SUBJID%in%ids&masks[[cat]]&is.na(onset))])))
 };wr(do.call(rbind,sens),out,'safety_boundary_sensitivity')
 # Optional later ADAE must be explicitly supplied; no concatenation or guessed linkage.
 later<-getOption('phmpt.S27_ADAE',NULL)
 status<-data.frame(check=c('S3_AE_to_ADAE','S3_CSR_summary_tables','Paired_EUA_ADAE_comparison','S27_ADAE_and_final_summaries'),status=c('COMPLETE_FOR_RELEASED_S3_ROWS',if(all(tab$match))'56_NUMERATOR_AND_DENOMINATOR_CELLS_MATCH'else'COUNT_MISMATCH_REVIEW','NOT_COMPLETED_NO_PAIRED_EUA_ADAE_IN_BUNDLE','PENDING'))
 if(!is.null(later)){
  ck(file.exists(later),'Explicit S27 ADAE not found');l<-if(grepl('[.]csv([.]gz)?$',later,ignore.case=TRUE))read.csv(later,stringsAsFactors=FALSE)else if(requireNamespace('haven',quietly=TRUE))as.data.frame(if(grepl('[.]xpt$',later,ignore.case=TRUE))haven::read_xpt(later)else haven::read_sas(later))else stop('Supply CSV or install haven for the explicit later ADAE.')
  ck(all(c('SUBJID','AESPID')%in%names(l)),'Later ADAE needs SUBJID/AESPID');lk<-key(l$SUBJID,l$AESPID);ck(!anyDuplicated(lk),'Later ADAE has ambiguous form keys')
  z<-ar[c('subject','arm','AE_spid','CE_ongoing_day7')];z$in_supplied_later_ADAE<-ar$key%in%lk;wr(z,out,'optional_S27_ADAE_bridge');status$status[4]<-'SUPPLIED_LATER_ADAE_LINKED_FINAL_TABLE_RULES_STILL_REQUIRE_VALIDATION'
 }
 eua_ae_path<-getOption('phmpt.EUA_AE',NULL);eua_adae_path<-getOption('phmpt.EUA_ADAE',NULL)
 if(!is.null(eua_ae_path)||!is.null(eua_adae_path)){
  ck(!is.null(eua_ae_path)&&!is.null(eua_adae_path),'Supply both phmpt.EUA_AE and phmpt.EUA_ADAE for the paired bridge')
  load_table<-function(p){ck(file.exists(p),paste('Missing',p));if(grepl('[.]csv([.]gz)?$',p,ignore.case=TRUE))return(read.csv(p,stringsAsFactors=FALSE));if(!requireNamespace('haven',quietly=TRUE))stop('Export to CSV or install haven for XPT/SAS7BDAT inputs');if(grepl('[.]xpt$',p,ignore.case=TRUE))as.data.frame(haven::read_xpt(p))else as.data.frame(haven::read_sas(p))}
  ea<-load_table(eua_ae_path);ed<-load_table(eua_adae_path)
  for(nm in c('ea','ed')){z<-get(nm);ck(all(c('SUBJID','AESPID')%in%names(z)),paste(nm,'requires SUBJID/AESPID'));z$key<-key(z$SUBJID,z$AESPID);ck(!anyDuplicated(z$key),paste('Ambiguous keys',nm));assign(nm,z)}
  ej<-match(ea$key,ed$key);sj<-match(ea$key,a$key);z<-data.frame(subject=ea$SUBJID,form=ea$AESPID,in_EUA_ADAE=!is.na(ej),same_form_in_S3_ADAE=!is.na(sj));for(field in intersect(c('AETERM','AESER','TRTEMFL','D28FL','APERIOD'),intersect(names(ed),names(a)))){z[[paste0('EUA_',field)]]<-ed[[field]][ej];z[[paste0('S3_',field)]]<-a[[field]][sj]};wr(z,out,'optional_paired_EUA_S3_ADAE_bridge');status$status[3]<-'PAIRED_FORM_BRIDGE_EXPORTED_EVENT_IDENTITY_AND_EUA_TABLE_RULES_REQUIRE_REVIEW'
 }
 wr(status,out,'completion_status');wr(data.frame(rows_AE=nrow(ae),rows_ADAE=nrow(a),source_rows_absent=sum(!global$in_S3_ADAE),ADAE_keys_without_AE=sum(!a$key%in%ae$key),seriousness_matches=sum(common),seriousness_mismatches=sum(text0(ae$AESER[common])!=text0(a$AESER[j[common]]))),out,'global_bridge_counts')
 v<-data.frame(Trace=c('Newly flagged reaction records','Marked ongoing at Day 7','Present and selected in S3 ADAE for 28-day table','Participants whose only selected event is one of these records'),'mRNA-1273'=c(arsum$flagged[1],arsum$ongoing_day7[1],arsum$D28_selected[1],arsum$only_D28_event_subjects[1]),Placebo=c(arsum$flagged[2],arsum$ongoing_day7[2],arsum$D28_selected[2],arsum$only_D28_event_subjects[2]),check.names=FALSE)
 wr(v,assets,'V41_AR_flag_analysis_trace');print_table_png(v,file.path(assets,'V41_AR_flag_analysis_trace.png'),'V41 Reaction records retained in the S3 safety analysis','These counts establish S3 treatment; later S27 table treatment remains unresolved.')
 html_report(file.path(out,'safety_results.html'),'Source events, analysis records and safety summaries',c('<p>The comparison uses subject plus sponsor form ID within each submission, not AESEQ across submissions. Exact matches to released S3 table totals establish aggregate reproducibility and record eligibility under the table program. They do not independently validate every source event or establish later S27 treatment.</p>',html_table(status),'<h2>New reaction flags</h2>',html_table(arsum),'<h2>Participants without licensing event rows</h2>',html_table(missing),'<h2>CSR Tables 7-12 and 7-13</h2>',html_table(tab),'<p>The explicit-unblinding-only subset excludes deaths or withdrawals before formal unblinding and is a selection diagnostic, not a preferred safety analysis. Missing onset dates are exported separately with counts if they were included. These censored-onset sensitivity counts do not replace the reproduced CSR population or table windows.</p>','<p>Fatal-event counts use event onset in Part A and final AE outcome; they are not counts of deaths occurring during Part A. CSR p291 explicitly explains this distinction.</p>'))
 print_asset_manifest(assets,c('V41'),'V41_safety_followup')
 invisible(status)
}
if(isTRUE(getOption('phmpt.followup_autorun',TRUE)))analyse_safety_summary_reconciliation()
