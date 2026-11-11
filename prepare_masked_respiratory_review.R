# Produces review packets; does NOT pretend automated classification is clinical adjudication.
# Two reviewers work independently. Only the custodian receives the restricted key.
local({r<-getOption('phmpt.followup_root',getwd());source(file.path(r,'audit_followup_helpers.R'),local=FALSE)})
prepare_masked_respiratory_review <- function(root=audit_root(),mode=getOption('phmpt.review_mode','prepare')) {
 out<-file.path(root,'results','masked_review');rev<-file.path(out,'reviewers_only');res<-file.path(out,'custodian_only');dir.create(rev,recursive=TRUE,showWarnings=FALSE);dir.create(res,recursive=TRUE,showWarnings=FALSE)
 k<-rd(root,'reviewer_key');s<-rd(root,'independent_review_sheet');ex<-rd(root,'clinical_excerpts');ep<-rd(root,'episodes');ck(!anyDuplicated(k$review_id)&&setequal(k$review_id,s$review_id),'Review packet keys disagree')
 ck(all(ex$review_id%in%k$review_id),'Unknown review ID in clinical excerpts');k$tested<-yes(k$any_sample_d0_7)
 wr(k,res,'linkage_key');wr(ex,rev,'coded_clinical_excerpts')
 # Do not distribute the input folder or custodian directory to reviewers.
 s$reviewer_id<-'';s$original_source_available<-'';s$masked_source_file<-'';s$source_locator<-'';s$source_redaction_checked<-'';s$temperature_threshold_established<-'';s$duration_48h_established<-'';s$respiratory_symptom_established<-'';s$protocol_version_confirmed<-'';s$clinical_justification_for_no_swab<-'';s$review_completed_utc<-''
 if(mode=='prepare'){
  for(nm in c('reviewer_A','reviewer_B'))if(!file.exists(file.path(rev,paste0(nm,'.csv'))))wr(s,rev,nm)
  j<-match(k$episode_id,ep$episode_id);m<-k[c('review_id','episode_id','subject','site')];m$trigger_date<-ep$trigger_date[j];m$needed_source<-'Original complaint narrative; contact/report time; temperature; duration; severity; applicable site protocol; dated clinical rationale'
  m$masked_source_file<-'';m$redaction_checked<-'';wr(m,res,'source_retrieval_manifest')
  writeLines(c('Give reviewers only reviewers_only/ after source masking is checked. Do not distribute inputs/ or custodian_only/.','Existing excerpts are coded/limited descriptions, not a substitute for complete source narratives. Obtain relevant source pages and contact logs; remove arm, PCR/serology results, site, subject identifiers and outcome-revealing commentary. Preserve relative timing and required clinical context.','Each reviewer independently uses YES / NO / INDETERMINATE for evaluation_warranted and swab_warranted. Record a reason, source locator, protocol version confirmation, source availability and completion time. Do not treat missing source material as evidence that a swab was unnecessary.','Original-source completion and coded-record completion are reported separately. Resolve disagreements before unmasking. No automated rating is supplied.'),file.path(rev,'READ_FIRST.txt'))
  wr(data.frame(status='PENDING_INDEPENDENT_REVIEW',episodes=nrow(k),untested=sum(!k$tested),tested_controls=sum(k$tested),human_reviews_completed=0),out,'review_status')
 } else if(mode=='validate') {
  readreview<-function(nm){p<-file.path(rev,paste0(nm,'.csv'));ck(file.exists(p),paste('Missing',p));d<-read.csv(p,stringsAsFactors=FALSE,na.strings='');ck(!anyDuplicated(d$review_id)&&setequal(d$review_id,k$review_id),paste('Invalid ID inventory',nm));d<-d[match(k$review_id,d$review_id),];for(c in c('evaluation_warranted','swab_warranted'))ck(all(toupper(text0(d[[c]]))%in%c('YES','NO','INDETERMINATE')),paste('Incomplete/invalid ratings',nm,c));ck(all(nzchar(text0(d$reviewer_id)))&&all(nzchar(text0(d$review_completed_utc))),paste('Reviewer identity/time missing',nm));ck(all(nzchar(text0(d$reason_or_records_needed))),paste('Missing clinical reason',nm));d}
  a<-readreview('reviewer_A');b<-readreview('reviewer_B');ck(all(a$reviewer_id!=b$reviewer_id),'Reviewers must be independent people')
  pair<-data.frame(review_id=k$review_id,A_evaluation=a$evaluation_warranted,B_evaluation=b$evaluation_warranted,A_swab=a$swab_warranted,B_swab=b$swab_warranted,A_source=a$original_source_available,B_source=b$original_source_available)
  pair$agreement<-toupper(pair$A_evaluation)==toupper(pair$B_evaluation)&toupper(pair$A_swab)==toupper(pair$B_swab);wr(pair,rev,'blinded_agreement')
  consensus<-file.path(rev,'consensus.csv')
  if(!file.exists(consensus)){z<-data.frame(review_id=k$review_id,evaluation_warranted=ifelse(pair$agreement,a$evaluation_warranted,''),swab_warranted=ifelse(pair$agreement,a$swab_warranted,''),consensus_reviewer='',reason='',locked_utc='');wr(z,rev,'consensus');stop('Blinded consensus sheet created. Complete and lock it before unmasking.')}
  z<-read.csv(consensus,stringsAsFactors=FALSE,na.strings='');ck(!anyDuplicated(z$review_id)&&setequal(z$review_id,k$review_id),'Consensus ID mismatch');z<-z[match(k$review_id,z$review_id),]
  for(c in c('evaluation_warranted','swab_warranted'))ck(all(toupper(text0(z[[c]]))%in%c('YES','NO','INDETERMINATE')),'Consensus is incomplete')
  ck(all(nzchar(text0(z$locked_utc)))&&all(nzchar(text0(z$consensus_reviewer)))&&all(nzchar(text0(z$reason))),'Consensus must be signed, reasoned and locked before unmasking')
  files<-file.path(rev,c('reviewer_A.csv','reviewer_B.csv','consensus.csv'));wr(data.frame(file=basename(files),MD5=unname(tools::md5sum(files)),locked_at=format(Sys.time(),tz='UTC',usetz=TRUE)),res,'locked_review_fingerprints')
  final<-cbind(k,z[setdiff(names(z),'review_id')]);final$source_review_complete<-yes(a$original_source_available)&yes(b$original_source_available)&yes(a$source_redaction_checked)&yes(b$source_redaction_checked)&nzchar(text0(a$source_locator))&nzchar(text0(b$source_locator))
  wr(final,res,'unmasked_adjudication');wr(as.data.frame(table(arm=final$arm,tested=final$tested,source_complete=final$source_review_complete,swab_warranted=final$swab_warranted)),out,'adjudicated_summary')
  wr(data.frame(status=if(all(final$source_review_complete))'COMPLETED_SOURCE_REVIEW'else'RATINGS_COMPLETE_SOURCE_REVIEW_PARTIAL',episodes=nrow(k),source_complete=sum(final$source_review_complete),human_reviews_completed=nrow(a)+nrow(b)),out,'review_status')
 } else stop('review_mode must be prepare or validate')
 invisible(out)
}
if(isTRUE(getOption('phmpt.followup_autorun',TRUE)))prepare_masked_respiratory_review()
