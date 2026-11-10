# Standalone in RStudio: set the project root as working directory, then select/run.
# Default: offline ledger replay. Set options(phmpt.crf.mode="source") to rescan
# original local PDFs and link against the pinned XPTs. See CRF_EVIDENCE.md.
analyse_crf_evidence <- function(root=getOption("phmpt.root", getwd()),
                          mode=getOption("phmpt.crf.mode","replay"),
                          data_root=getOption("phmpt.crf.data_root",root),
                          out=getOption("phmpt.crf.out",NULL)) {
 root <- normalizePath(root,mustWork=TRUE); data_root <- normalizePath(data_root,mustWork=TRUE)
 if (is.null(out)) dir.create(file.path(root,"analysis"),recursive=TRUE,showWarnings=FALSE)
 if (is.null(out)) out <- tempfile(paste0("crf_",mode,"_"),tmpdir=file.path(root,"analysis"))
 ev <- file.path(root,"evidence"); tool <- file.path(root,"crf_audit")
 if (!file.exists(file.path(tool,"replay_crf.R"))) stop("Use the moderna_trial_analysis project as root.")
 if (dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE))) stop("Output directory must be new or empty.")
 dir.create(out,recursive=TRUE,showWarnings=FALSE);if(!dir.exists(out))stop("Cannot create output: ",out)
 lib <- new.env(parent=environment())
 if (mode=="replay") {
   sys.source(file.path(tool,"replay_crf.R"),envir=lib);sys.source(file.path(tool,"replay_forms.R"),envir=lib)
   got <- lib$crf_summary(root,out); ref <- read.csv(file.path(ev,"crf_source","reference_counts.csv"))
   stopifnot(identical(got$measure,ref$measure),all(got$value==ref$value))
   lib$replay_566_review(root,out)
   readcase <- function(f) read.csv(file.path(ev,f),colClasses="character",na.strings=character(),check.names=FALSE)
   q<-readcase("form_dispositions_566.csv"); s<-readcase("serious_12.csv"); stopifnot(setequal(s$key,q$key[q$EUA_serious=="Y"]),length(unique(s$key))==12L)
   arm<-readcase("crf_source/Moderna_171_dose2_visit_comparison.csv"); flag<-arm[tolower(arm$discrepant_final_D2)=="true",]
   stopifnot(nrow(arm)==171L,!anyDuplicated(arm$subject),nrow(flag)==50L,sum(flag$arm=="Placebo")==48L,all(flag$PPROTFL=="N"),all(flag$SAFFL=="Y"),sum(tolower(flag$US331_explicit_workaround)=="true")==3L)
 } else if (mode=="source") {
   py <- getOption("phmpt.crf.python", unname(Sys.which(if (.Platform$OS.type=="windows") "python" else "python3")))
   if (!nzchar(py)) stop("Install Python 3.10+ with pandas and PyMuPDF or set options(phmpt.crf.python=...).")
   run <- function(name,args) {if (system2(py,shQuote(c(file.path(tool,name),args)))!=0L) stop("Failed stage: ",name)}
   only <- isTRUE(getOption("phmpt.crf.allow_partial",FALSE)); sc <- file.path(out,"scan"); jo <- file.path(out,"joins")
   loc <- function(pattern) {
     a<-list.files(data_root,pattern=pattern,recursive=TRUE,full.names=TRUE,ignore.case=TRUE)
     if (length(a)!=1L) stop("Expected exactly one ",pattern," under ",data_root,"; for duplicate snapshots use the explicit Python CLI.")
     normalizePath(a)
   }
   files <- c("--ex",loc("_125752_S3_M5_mrna-1273-p301_S_ex[.]xpt$"),
    "--ae-s5",loc("_27073_S5_M5_mrna-1273-p301_S_ae[.]xpt$"),
    "--ae-s3",loc("_125752_S3_M5_mrna-1273-p301_S_ae[.]xpt$"),
    "--ae-s27",loc("_125752_S27_M5_mrna-1273-p301_S_ae[.]xpt$"))
   run("extract_fields.py",c("--root",data_root,"--out",sc,"--jobs",as.character(getOption("phmpt.crf.jobs",2L)),if (only) "--allow-partial"))
   run("reconcile_fields.py",c("--scan",sc,"--out",jo,files))
   run("verify_extract.py",c("--scan",sc,"--joins",jo))
   run("dose2_visit_history.py",c("--root",data_root,"--joins",jo,"--ex",files[2],"--out",file.path(out,"dose2_workflow")))
 } else stop("mode must be replay or source")
 message("Results: ",normalizePath(out));invisible(out)
}
if(isTRUE(getOption("phmpt.autorun",TRUE))) analyse_crf_evidence()
