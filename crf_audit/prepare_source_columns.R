# Optional raw JMP route. Uses the bundled project reader for the named S3 submissions.
# Run from project root (copy paths quoted in CRF_EVIDENCE.md). Source only, no autorun.
crf_source_columns <- function(root=getwd(), adsl_file=NULL, suppdv_file=NULL, out) {
 opts<-options(phmpt.autorun=FALSE);on.exit(options(opts),add=TRUE)
 env<-new.env(parent=environment());sys.source(file.path(root,"crf_audit","read_source_jmp.R"),envir=env)
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 if (!is.null(adsl_file)) {
   cols<-c("SUBJID","TRT01P","DOSE1FL","DOSE2FL","PPROTFL","SAFFL","FASFL","PPEREAS","UBLDDT","PDVDT")
   a<-env$crf_read_jmp(adsl_file,cols);for (n in c("UBLDDT","PDVDT")) if (is.numeric(a[[n]])) a[[n]]<-as.character(as.Date(a[[n]],origin="1960-01-01"))
   write.csv(a,file.path(out,"ADSL_columns.csv"),row.names=FALSE,na="")
 }
 if (!is.null(suppdv_file)) {
   ids<-env$crf_read_jmp(suppdv_file,"USUBJID"); take<-which(endsWith(ids$USUBJID,"US393-2018")); if(!length(take))stop("Case not found in SUPPDV")
   x<-env$crf_read_jmp(suppdv_file,c("USUBJID","IDVAR","IDVARVAL","QNAM","QVAL"),rows=take)
   want<-sub(".*-US","US",x$USUBJID);want<-gsub("-","",want,fixed=TRUE)
   x<-x[want=="US3932018" & x$IDVAR=="DVSEQ" & x$QNAM %in% c("DECSRTNL","DECSRTN1","DECSRTN2","DECSRTN3"),]
   write.csv(x,file.path(out,"case_SUPPDV.csv"),row.names=FALSE,na="")
 }
}
