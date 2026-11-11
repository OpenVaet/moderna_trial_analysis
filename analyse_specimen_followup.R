# Distinguish submission, test-record identity, accession and specimen type.
# Never resolve a mixed-result date by overwriting its individual laboratory rows.
local({r<-getOption('phmpt.followup_root',getwd());source(file.path(r,'audit_followup_helpers.R'),local=FALSE)})
analyse_specimen_followup <- function(root=audit_root()) {
 out<-file.path(root,'results','specimens');dir.create(out,recursive=TRUE,showWarnings=FALSE)
 refs<-rd(root,'pcr_comparison');raw<-rd(root,'raw_PCR_specimen_fields');eua<-rd(root,'raw_EUA_PCR_specimen_fields');p<-rd(root,'pcr_records');n<-rd(root,'serology_records');d<-rd(root,'subjects_BLA_PRIMARY');f<-rd(root,'subject_followup')
 normalize<-function(z){z<-toupper(trimws(text0(z)));ifelse(z%in%c('POSITIVE','DETECTED','Y','YES'),'POSITIVE',ifelse(z%in%c('NEGATIVE','NOT DETECTED','N','NO'),'NEGATIVE','INDETERMINATE'))}
 fix<-function(x,release){x$release<-release;x$date<-substr(text0(x$MBDTC),1,10);x$binary_result<-normalize(ifelse(nzchar(text0(x$MBSTRESC)),text0(x$MBSTRESC),text0(x$MBORRES)));x$date_key<-key(x$subject,x$date);x$accession<-text0(x$MBREFID);x$unique_record_key<-key(release,x$subject,x$sequence);x}
 raw<-fix(raw,'S3');eua<-fix(eua,'EUA');for(col in setdiff(names(raw),names(eua)))eua[[col]]<-NA_character_;all<-rbind(raw,eua[names(raw)])
 ck(!anyDuplicated(all$unique_record_key),'Raw PCR record keys are not unique within submission')
 refs$release<-ifelse(refs$seq=='S6','EUA','S3');refs$date_key<-key(refs$sid,refs$test_date);refs$comparison_key<-key(refs$release,refs$date_key)
 ix<-split(seq_len(nrow(all)),key(all$release,all$date_key));ledger<-list();detail<-list()
 for(i in seq_len(nrow(refs))){r<-refs[i,];z<-all[ix[[r$comparison_key]],,drop=FALSE];if(is.null(z))z<-all[FALSE,];pos<-any(z$binary_result=='POSITIVE');neg<-any(z$binary_result=='NEGATIVE');known<-z[nzchar(z$accession),]
  both_accession<-FALSE;if(nrow(known))both_accession<-any(vapply(split(known$binary_result,key(known$subject,known$accession)),function(x)all(c('POSITIVE','NEGATIVE')%in%x),logical(1)))
  state<-if(both_accession)'SAME_ACCESSION_DISCORDANCE'else if(pos&&neg){if(all(nzchar(z$accession[z$binary_result%in%c('POSITIVE','NEGATIVE')])) )'DISTINCT_ACCESSIONS_MIXED_DATE'else'MIXED_DATE_WITH_UNIDENTIFIED_SPECIMEN'}else if(pos)'POSITIVE_ONLY'else if(neg)'NEGATIVE_ONLY'else'NO_BINARY_RESULT'
  ledger[[i]]<-cbind(r,raw_records=nrow(z),positive_records=sum(z$binary_result=='POSITIVE'),negative_records=sum(z$binary_result=='NEGATIVE'),distinct_known_accessions=length(unique(known$accession)),records_without_accession=sum(!nzchar(z$accession)),opposite_results_in_distinct_known_accessions=any(known$binary_result=='POSITIVE')&&any(known$binary_result=='NEGATIVE')&&!both_accession,normalized_date_status=state,CRF_result_represented=r$CRF_result%in%z$binary_result,all_binary_agree=!(pos&&neg)&r$CRF_result%in%z$binary_result)
  if(!yes(r$agree)||!r$CRF_result%in%z$binary_result){if(nrow(z)){z$CRF_document<-r$doc;z$CRF_page<-r$positive_field_pg;z$CRF_result<-r$CRF_result;detail[[length(detail)+1]]<-z}}
 }
 res<-do.call(rbind,ledger);wr(res,out,'CRF_PCR_submission_specific_reconciliation');if(length(detail))wr(do.call(rbind,detail),out,'flagged_dates_individual_test_records')
 wr(res[!yes(res$agree),],out,'previously_flagged_comparisons')
 wr(as.data.frame(table(previous=refs$same_date_MB,reconciled=res$normalized_date_status,CRF_represented=res$CRF_result_represented)),out,'PCR_resolution_counts')
 # Mixed within a date is not a conflict within a specimen or source-record key.
 all$binary<-all$binary_result%in%c('POSITIVE','NEGATIVE');spec<-all[nzchar(all$accession)&all$binary,];spkey<-key(spec$release,spec$subject,spec$accession)
 bad<-names(which(vapply(split(spec$binary_result,spkey),function(x)length(unique(x))>1,logical(1))));wr(spec[spkey%in%bad,],out,'same_accession_conflicts')
 wr(all[all$subject=='US3712171'&all$date%in%c('2020-09-25','2020-10-01'),],out,'US3712171_cross_snapshot_test_rows')
 wr(p[yes(p$result_issue)|yes(p$date_issue)|!yes(p$raw_present),],out,'PCR_record_level_review_queue')
 wr(n[yes(n$result_issue)|yes(n$date_issue)|!yes(n$raw_present),],out,'serology_record_level_review_queue')
 # Preserve the published run's result-blind Day57 selection. Sensitivities may censor
 # the selected draw; they do not substitute a different, more favorable sample.
 ck(!anyDuplicated(d$subject),'Duplicate Day57 participant');j<-match(d$subject,f$subject);ck(!anyNA(j),'Day57 follow-up linkage missing');bs<-boundaries(f);sumrows<-list();follow<-list()
 si<-split(seq_len(nrow(n)),n$subject);pi<-split(seq_len(nrow(p)),p$subject)
 for(bn in names(bs)){
  end<-bs[[bn]][j];draw_ok<-!is.na(end)&!is.na(date0(d$draw_date))&date0(d$draw_date)<=end;retained<-!is.na(end)
  for(arm in c('mRNA-1273','Placebo')){
   armset<-d$arm==arm&retained;binary<-armset&draw_ok&yes(d$binary);anti<-binary&yes(d$anti_n_positive);without<-anti&text0(d$pcr_status)%in%c('NEGATIVE_ONLY','NO_COMPLETED_PCR');cal<-binary&text0(d$pcr_status)=='POSITIVE'&!is.na(d$pcr_lag)&d$pcr_lag>=14
   sumrows[[length(sumrows)+1]]<-data.frame(boundary=bn,arm=arm,original_due=sum(d$arm==arm),retained_due=sum(armset),binary_draw=sum(binary),antiN_positive=sum(anti),positive_without_recorded_PCR=sum(without),PCR_positive_lag14_denominator=sum(cal),PCR_positive_lag14_antiN_positive=sum(cal&yes(d$anti_n_positive)),draws_censored=sum(d$arm==arm&yes(d$binary)&!draw_ok))
  }
 }
 for(i in seq_len(nrow(d))){
  nidx<-si[[d$subject[i]]];pidx<-pi[[d$subject[i]]];nn<-n[if(is.null(nidx))integer()else nidx,];pp<-p[if(is.null(pidx))integer()else pidx,];dd<-date0(d$draw_date[i]);be<-date0(d$scope_end[i]);base<-date0(d$baseline_date[i])
  # Known PCR infection with no eligible later serology stays visible, never a negative response.
  goodp<-yes(pp$record_usable)&yes(pp$method_ok)&!yes(pp$derived)&!yes(pp$boundary_visit)&pp$result=='POSITIVE'&!is.na(date0(pp$date))&date0(pp$date)>base&date0(pp$date)<=be
  first<-if(any(goodp))min(date0(pp$date[goodp]))else as.Date(NA)
  goodn<-yes(nn$record_usable)&!yes(nn$derived)&!yes(nn$boundary_visit)&nn$result%in%c('POSITIVE','NEGATIVE')&!is.na(date0(nn$date))&date0(nn$date)<=be
  after<-goodn&!is.na(first)&date0(nn$date)>=first+14;prior<-goodn&!is.na(dd)&date0(nn$date)<dd
  follow[[i]]<-data.frame(subject=d$subject[i],arm=d$arm[i],day57_draw=d$draw_date[i],draw_status=d$draw_status[i],PCR_before_draw=d$pcr_status[i],first_followup_PCR_positive=as.character(first),later_serology_ge14_records=sum(after),later_serology_ge14_positive=sum(after&nn$result=='POSITIVE'),last_negative_before_draw=if(any(prior&nn$result=='NEGATIVE'))as.character(max(date0(nn$date[prior&nn$result=='NEGATIVE'])))else'',day57_serology_record_ids=d$serology_record_ids[i],later_serology_record_ids=paste(nn$record_key[after],collapse=';'),day57_source_refs=d$serology_source_refs[i])
 }
 fs<-do.call(rbind,follow);sr<-do.call(rbind,sumrows);wr(sr,out,'day57_boundary_sensitivity');wr(fs,out,'serology_followup_ledger')
 fs$PCR_positive<-nzchar(text0(fs$first_followup_PCR_positive));wr(as.data.frame(table(arm=fs$arm,PCR_positive=fs$PCR_positive,later_sample14=fs$later_serology_ge14_records>0)),out,'PCR_followup_availability')
 html_report(file.path(out,'specimen_results.html'),'Specimen identity and serology follow-up',c('<p>Each CRF is compared with the corresponding EUA or S3 laboratory snapshot. MBSEQ identifies a row within that submission; MBREFID supplies a laboratory reference when present. A different result on the same calendar date is not, by itself, a changed result for one specimen. Blank accession fields prevent complete specimen identity adjudication.</p>',html_table(res[!yes(res$agree),c('sid','seq','test_date','CRF_result','normalized_date_status','CRF_result_represented')]),'<h2>Day57: fixed draw, alternative boundaries</h2>',html_table(sr),'<p>These are follow-up and detection summaries, not corrected infection counts or vaccine efficacy. The calibration remains conditional on PCR-detected infection and adequate later sampling. Serology accession identities are not established by this record-level follow-up ledger.</p>'))
 invisible(res)
}
if(isTRUE(getOption('phmpt.followup_autorun',TRUE)))analyse_specimen_followup()
