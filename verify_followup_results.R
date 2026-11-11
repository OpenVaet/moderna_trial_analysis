# Independent result invariants; run after run_audit_followup.R.
verify_followup_results <- function(root=getOption('phmpt.followup_root',getwd())) {
 read<-function(folder,file)read.csv(file.path(root,'results',folder,paste0(file,'.csv')),stringsAsFactors=FALSE)
 checks<-list();check<-function(name,x){checks[[length(checks)+1]]<<-data.frame(check=name,passed=isTRUE(x));if(!isTRUE(x))stop('Verification failed: ',name)}
 z<-read('safety','CSR_summary_reproduction');check('56 CSR numerical cells and denominators reproduced',nrow(z)==56&&all(z$match))
 z<-read('safety','global_bridge_counts');check('Complete S3 bridge dimensions',z$rows_AE==39215&&z$rows_ADAE==37756&&z$source_rows_absent==1459&&z$ADAE_keys_without_AE==0&&z$seriousness_mismatches==0)
 z<-read('safety','new_removal_flags_ADAE_bridge');check('All 62 retained and selected, including 51 ongoing',nrow(z)==62&&all(z$in_S3_ADAE)&all(z$S3_D28_selected)&sum(z$CE_ongoing_day7)==51)
 z<-read('safety','participants_without_licensing_AE');check('Three targeted participants have no S3/S27 AE or S3 ADAE',nrow(z)==3&&all(z$S3_AE_rows+z$S27_AE_rows+z$S3_ADAE_rows==0))
 z<-read('safety','seriousness_changes_ADAE_bridge');check('All licensing-side seriousness records reach ADAE',nrow(z)==31&&all(z$in_S3_ADAE)&&all(z$S3_partA_selected)&&sum(z$S3_ADAE_serious=='Y')==15)
 z<-read('testing','matched_reproduction');check('Original exact-stratum membership and point estimates',z$vaccine_supported==504&&z$placebo_supported==539&&z$strata==365&&abs(z$vaccine_standardized_pct-55.380532)<1e-5&&abs(z$placebo_standardized_pct-56.234972)<1e-5)
 z<-read('testing','matched_subset_summary');z<-z[z$group=='UNMATCHED',];check('Unmatched population retained',identical(z$episodes,c(2391L,3055L))&&identical(z$forms3,c(1134L,1862L))&&all(z$sites==99))
 z<-read('testing','complete_and_selected_comparisons');a<-z[z$scope=='Complete_seven_day_opportunity'&z$endpoint=='swab_d0_7',];b<-z[z$scope=='No_linked_positive_PCR_D0_7'&z$endpoint=='swab_d0_7',]
 check('Complete and positivity-excluded denominator identities',a$vaccine_n==2818&&a$placebo_n==3526&&a$vaccine_n-b$vaccine_n==53&&a$placebo_n-b$placebo_n==689&&b$vaccine_x==1115&&b$placebo_x==1220)
 z<-read('specimens','CRF_PCR_submission_specific_reconciliation');check('640 CRF comparisons, 630 uniform agreements',nrow(z)==640&&sum(z$all_binary_agree)==630)
 z<-read('specimens','previously_flagged_comparisons');m<-z[z$same_date_MB=='MIXED',];check('Nine mixed comparisons represent seven dates and all retain matching result',nrow(m)==9&&length(unique(paste(m$sid,m$test_date)))==7&&all(m$CRF_result_represented))
 z<-read('specimens','day57_boundary_sensitivity');z<-z[z$boundary=='recorded',];check('Original Day57 sample and calibration preserved',identical(z$original_due,c(14745L,14743L))&&identical(z$positive_without_recorded_PCR,c(24L,47L))&&identical(z$PCR_positive_lag14_denominator,c(21L,79L))&&identical(z$PCR_positive_lag14_antiN_positive,c(12L,61L)))
 z<-read('masked_review','review_status');check('No clinical review falsely claimed',z$status=='PENDING_INDEPENDENT_REVIEW'&&z$episodes==1140&&z$untested==570&&z$tested_controls==570&&z$human_reviews_completed==0)
 z<-read.csv(file.path(root,'audit_assets','V02_table_data.csv'));check('V02 six participant/episode rows preserved',nrow(z)==6&&identical(z$vaccine_tested,c(1634L,1624L,1557L,1465L,1143L,863L)))
 ans<-do.call(rbind,checks);write.csv(ans,file.path(root,'results','verification.csv'),row.names=FALSE);message(nrow(ans),' verification checks passed');invisible(ans)
}
if(isTRUE(getOption('phmpt.followup_verify_autorun',TRUE)))verify_followup_results()
