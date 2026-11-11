# Standalone launcher: set the extracted folder as the working directory in RStudio.
# Or options(phmpt.followup_root='D:/Projects/moderna_trial_analysis/audit_followup').
run_audit_followup <- function(root=getOption('phmpt.followup_root',getwd())) {
 root<-normalizePath(root,winslash='/',mustWork=TRUE);old<-options(phmpt.followup_root=root,phmpt.followup_autorun=FALSE,phmpt.autorun=FALSE);on.exit(options(old),add=TRUE)
 scripts<-c('analyse_safety_summary_reconciliation.R','analyse_testing_robustness.R','analyse_specimen_followup.R','prepare_masked_respiratory_review.R')
 for(s in scripts)if(!file.exists(file.path(root,s)))stop('Missing script: ',s)
 source(file.path(root,'audit_followup_helpers.R'));audit_root(root)
 manifest<-file.path(root,'inputs','INPUT_MD5.csv');if(file.exists(manifest)){m<-read.csv(manifest,stringsAsFactors=FALSE);actual<-unname(tools::md5sum(file.path(root,'inputs',m$file)));if(anyNA(actual)||any(actual!=m$MD5))stop('Pinned inputs differ from INPUT_MD5.csv. Re-extract or deliberately rebuild and document the input manifest.')}
 for(s in scripts)source(file.path(root,s),local=FALSE)
 message('[1/5] S3 safety tables and source-event traceability');analyse_safety_summary_reconciliation(root)
 message('[2/5] Complete episode coverage, matching and boundary sensitivity');analyse_testing_robustness(root)
 message('[3/5] PCR specimen records and serology follow-up');analyse_specimen_followup(root)
 message('[4/5] Masked respiratory-review workflow (no clinical adjudication implied)');prepare_masked_respiratory_review(root)
 message('[5/5] V02 interpretation and print layout');source(file.path(root,'analyse_intro_figures.R'),local=FALSE);render_intro_figures(root=root,selected='V02',out=file.path(root,'audit_assets'),inputs=file.path(root,'table_inputs_intro'))
 writeLines(capture.output(sessionInfo()),file.path(root,'results','sessionInfo.txt'));message('Complete. See results/, audit_assets/ and audit_edits.html. Independent source review and unavailable later analysis files remain explicitly pending.')
 invisible(root)
}
if(isTRUE(getOption('phmpt.followup_launcher_autorun',TRUE)))run_audit_followup()
