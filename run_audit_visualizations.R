# Moderna COVE audit V01–V49 | standalone launcher | two introductory graphics
# Extract the script bundle into your project root, open the project in RStudio,
# select this entire file and Run, or: Rscript run_audit_visualizations.R
# Include four unnumbered supplemental charts/maps: options(phmpt.figures_include_supplementary=TRUE).
# Preview the sequence only: options(phmpt.figures_plan_only=TRUE); source(this file).
# Each generator runs in a fresh R process. Individual generators also remain standalone.
# V03–V05 accepts analyse_moderna_demographics.R or its earlier name analyse_demographics.R.
figure_plan <- function(include_supplementary=FALSE) {
  numbered<-data.frame(script=paste0(c("analyse_intro_figures","analyse_moderna_demographics","analyse_dose_timeline","analyse_deviations","analyse_audit_tables",
    "analyse_ae_covid_symptoms","analyse_weekly_symptoms_ae","analyse_background_covid",
    "analyse_site_symptom_recording","analyse_pcr_testing_by_lab","analyse_audit_tables","analyse_testing_robustness","analyse_symptom_testing_heatmaps",
    "analyse_site_testing_imbalance","analyse_anti_n_seroconversion","analyse_day57_anti_n_sensitivity",
    "analyse_eua_bla_changes","analyse_safety_summary_reconciliation","analyse_efficacy_cases","analyse_death_reporting","analyse_subject_numbering"),".R"),
    ids=c("V01–V02","V03–V05","V06","V07, V10, V11","V08–V09","V12","V13–V14","background prerequisite","V15","V16–V18","V19","V20–V21",
    "V22–V31: episode data and descriptive plots","V22–V31: final models and plots","V32","V33–V35",
    "V36–V40","V41","V42–V44","V45","V46–V49"),stringsAsFactors=FALSE)
  if(!include_supplementary)return(numbered)
  rbind(numbered,data.frame(script=c("analyse_anti_n_testing.R","analyse_us_weekly_covid_tests.R",
    "analyse_background_covid_map_sites_only.R","analyse_background_covid_map.R"),
    ids=c("Unnumbered: anti-N testing","Unnumbered: national testing","Unnumbered: site counties map","Unnumbered: all US counties map"),stringsAsFactors=FALSE))
}
figure_csv <- function(x,p) write.csv(x,p,row.names=FALSE,na="",fileEncoding="UTF-8")
figure_resolve_plan <- function(plan,scripts_dir) {
  modern<-"analyse_moderna_demographics.R";legacy<-"analyse_demographics.R"
  present<-function(name)file_test("-f",file.path(scripts_dir,name))
  if(!present(modern)&&present(legacy))plan$script[plan$script==modern]<-legacy
  missing<-plan$script[!vapply(plan$script,present,logical(1))]
  if(length(missing)){
    missing[missing==modern]<-paste0(modern," (or ",legacy,")")
    stop("Missing generator(s) in ",scripts_dir,": ",paste(missing,collapse=", "),
      ". Place the missing file(s) there, or set options(phmpt.scripts_dir=...) to the directory containing your generators.",call.=FALSE)
  }
  plan
}
# Older V03–V05 releases hard-code audit_assets and prod_data. For the legacy
# filename, make a run-local copy with only those settings adapted. Leave its
# analysis code and the user's original file alone. Refuse older scripts that
# cannot accept the staged output directory.
figure_legacy_demographics <- function(path,root) {
  s<-readLines(path,warn=FALSE,encoding="UTF-8")
  k<-grep("^\\s*PHMPT_AUDIT_ASSET_DIR\\s*<-",s,perl=TRUE)
  if(length(k)!=1L||!any(grepl("V03_V05_demographics_manifest.csv",s,fixed=TRUE))||
     !any(grepl("^\\s*bla1_export_audit_assets\\s*<-\\s*function",s,perl=TRUE))||
     !any(grepl("audit_asset_dir\\s*=\\s*PHMPT_AUDIT_ASSET_DIR\\b",s,perl=TRUE)))
    stop("Found ",path,", but it does not expose the V03–V05 audit output settings. Place the supplied analyse_moderna_demographics.R in the generator directory and rerun.",call.=FALSE)
  s[k]<-'PHMPT_AUDIT_ASSET_DIR <- getOption("phmpt.audit_asset_dir","audit_assets")'
  folder<-if(dir.exists(file.path(root,"prod_data_moderna")))"prod_data_moderna" else "prod_data"
  # Only the two literal release-folder arguments used in the older generator.
  s<-gsub('file[.]path\\((p|root),[ \\t]*([\x22\x27])prod_data(_moderna)?\\2',
    paste0('file.path(\\1,"',folder,'"'),s,perl=TRUE)
  s
}
figure_manifest <- function(dir) {
  files<-list.files(dir,pattern="^V[0-9]{2}.*_manifest[.]csv$",full.names=TRUE)
  if(!length(files))return(data.frame())
  do.call(rbind,lapply(files,function(f) {
    d<-read.csv(f,fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE,check.names=FALSE)
    if(!all(c("visual","file")%in%names(d))||anyNA(d$file)||any(d$file!=basename(d$file)))stop("Invalid figure manifest: ",f)
    data.frame(visual=d$visual,file=d$file,status=if("status"%in%names(d))d$status else "AVAILABLE",
      md5=if("md5"%in%names(d))d$md5 else NA_character_,stringsAsFactors=FALSE)
  }))
}
# Old IDs are used solely to recognize previous generated assets when publishing.
# Same suffix, older figure ID. Recognize both the immediately preceding
# 44-figure edition and the earlier 41-figure edition. New figures have no predecessor.
figure_previous_name <- function(n) {
  previous <- c(setNames(1:19,sprintf("V%02d",1:19)),V20=48L,V21=47L,
                setNames(20:38,sprintf("V%02d",22:40)),V41=49L,
                setNames(39:46,sprintf("V%02d",42:49)))
  convert <- function(x,edition=0L) {
    hits <- gregexpr("V[0-9]{2}(?![0-9])",x,perl=TRUE)
    tags <- regmatches(x,hits)[[1]]
    if(!length(tags)||any(!tags%in%names(previous)))return(character())
    ids <- unname(previous[tags])
    if(edition>0L) {
      if(any(ids>46L))return(character())
      ids <- ids-2L
      if(any(ids<1L))return(character())
      if(edition==2L) {
        if(any(ids%in%c(6L,7L,17L)))return(character())
        ids <- ifelse(ids>=18L,ids-3L,ifelse(ids>=8L,ids-2L,ids))
      }
    }
    regmatches(x,hits)<-list(sprintf("V%02d",ids));x
  }
  unique(unlist(lapply(n,function(x)c(convert(x),convert(x,1L),convert(x,2L))),use.names=FALSE))
}

figure_publish <- function(stage,dest,run) {
  backup<-file.path(run,"previous_audit_assets")
  if(dir.exists(dest)) {
    # Same filesystem as project; preserve unrelated files and manual assets.
    new_files<-list.files(stage,all.files=TRUE,no..=TRUE)
    old_files<-list.files(dest,all.files=TRUE,no..=TRUE)
    old_manifests<-old_files[old_files%in%unique(c(new_files,figure_previous_name(new_files))) & grepl("_manifest[.]csv$",old_files)]
    old_png<-unlist(lapply(old_manifests,function(f){
      x<-read.csv(file.path(dest,f),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE)
      if("file"%in%names(x))x$file[x$file==basename(x$file)] else character()
    }),use.names=FALSE)
    keep<-setdiff(old_files,unique(c(old_png,new_files,figure_previous_name(new_files))))
    if(length(keep)&&!all(file.copy(file.path(dest,keep),stage,recursive=TRUE,copy.date=TRUE)))
      stop("Could not preserve unrelated existing assets. Previous directory remains active.")
    if(!file.rename(dest,backup))stop("Could not archive the previous audit_assets directory (close files open in other applications).")
  }
  if(!file.rename(stage,dest)){
    if(dir.exists(backup)&&!file.rename(backup,dest))stop("Publish failed. Previous assets are in ",backup)
    stop("Publish failed; the previous audit_assets directory was restored.")
  }
}
run_audit_visualizations <- function(root=getOption("phmpt.root",getwd()),
  scripts_dir=getOption("phmpt.scripts_dir",root), plan_only=isTRUE(getOption("phmpt.figures_plan_only",FALSE)),
  include_supplementary=isTRUE(getOption("phmpt.figures_include_supplementary",FALSE)),
  reuse_background=isTRUE(getOption("phmpt.figures_reuse_background",TRUE)),
  allow_flagged=isTRUE(getOption("phmpt.figures_allow_flagged",FALSE))) {
  root<-normalizePath(root,winslash="/",mustWork=TRUE);scripts_dir<-normalizePath(scripts_dir,winslash="/",mustWork=TRUE)
  plan<-figure_resolve_plan(figure_plan(include_supplementary),scripts_dir);paths<-file.path(scripts_dir,plan$script)
  print(plan,row.names=FALSE)
  for(p in paths)invisible(parse(p,encoding="UTF-8"))
  legacy<-match("analyse_demographics.R",plan$script)
  if(!is.na(legacy)){
    message("V03–V05: using ",paths[legacy]," (legacy filename).")
    adapted<-figure_legacy_demographics(paths[legacy],root)
    invisible(parse(text=adapted,encoding="UTF-8"))
  }
  if(plan_only)return(invisible(plan))
  rscript<-getOption("phmpt.figure_rscript",file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe" else "Rscript"))
  if(!file.exists(rscript))stop("Rscript executable not found: ",rscript)
  base<-file.path(root,"analysis","figure_runs");dir.create(base,recursive=TRUE,showWarnings=FALSE)
  lock<-file.path(base,"launcher.lock")
  if(!dir.create(lock,showWarnings=FALSE))stop("Another run has the launcher lock: ",lock,". If a prior run was interrupted, confirm it has stopped before removing that directory.")
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  run<-tempfile(paste0(format(Sys.time(),"%Y%m%d_%H%M%S"),"_"),tmpdir=base)
  if(!dir.create(run))stop("Cannot create run directory: ",run)
  logs<-file.path(run,"logs");dir.create(logs)
  stage<-file.path(run,"new_audit_assets");dir.create(stage)
  bg<-file.path(root,"analysis","background_covid","site_weekly.csv")
  weekly<-file.path(run,"weekly_symptoms");heat<-file.path(run,"heatmap");site<-file.path(run,"site_models")
  bg_ok<-function() {
    if(!file.exists(bg))return(FALSE)
    d<-tryCatch(read.csv(bg,fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE),error=function(e)NULL)
    !is.null(d)&&all(c("SITEID","week_start","county_fips","county_name","state_abbr","state_name","county_cases",
      "county_positivity_pct","state_cases","state_positivity_pct")%in%names(d))&&length(unique(d$SITEID))==99L&&
      any(is.finite(suppressWarnings(as.numeric(d$county_positivity_pct))))&&!anyDuplicated(d[c("SITEID","week_start")])
  }
  # No command-line arguments are passed through to the child generator (the
  # site-audit script historically interprets any such arguments as input paths).
  opts<-list(phmpt.root=root,phmpt.autorun=TRUE,phmpt.audit_asset_dir=stage,
    phmpt.followup_root=getOption("phmpt.followup_root",scripts_dir),phmpt.followup_autorun=TRUE,
    phmpt.weekly_output=weekly,phmpt.weekly_dir=weekly,phmpt.heatmap_output=heat,
    phmpt.audit_data_dir=heat,phmpt.audit_output_dir=site,phmpt.report_only=FALSE,
    phmpt.background_file=bg,
    phmpt.intro_input_dir=getOption("phmpt.intro_input_dir",file.path(scripts_dir,"table_inputs_intro")),
    phmpt.intro_coverage_file=getOption("phmpt.intro_coverage_file"),
    phmpt.intro_lab_file=getOption("phmpt.intro_lab_file"),
    phmpt.intro_trigger_file=getOption("phmpt.intro_trigger_file"),
    phmpt.deviation_tables_dir=getOption("phmpt.deviation_tables_dir",file.path(scripts_dir,"table_inputs_v022","deviations")),
    phmpt.trigger_table_dir=getOption("phmpt.trigger_table_dir",file.path(scripts_dir,"table_inputs_v022","protocol")),
    phmpt.v17_loess_span=getOption("phmpt.v17_loess_span",0.40),
    phmpt.v17_p0_randomized=getOption("phmpt.v17_p0_randomized"),
    phmpt.v17_p0_eligible=getOption("phmpt.v17_p0_eligible"),
    phmpt.v17_adjustment=match.arg(getOption("phmpt.v17_adjustment","funnel"),c("funnel","nominal","bonferroni")))
  state<-data.frame(step=seq_len(nrow(plan)),plan,status="PENDING",started="",seconds=NA_real_,log="",message="",stringsAsFactors=FALSE)
  save_state<-function()figure_csv(state,file.path(run,"launcher_status.csv"))
  save_state();message("Run directory: ",run)
  exec_paths<-paths
  if(!is.na(legacy)){
    exec_paths[legacy]<-file.path(logs,"analyse_demographics.R")
    writeLines(enc2utf8(adapted),exec_paths[legacy],useBytes=TRUE)
  }
  figure_csv(data.frame(script=plan$script,source=paths,md5=unname(tools::md5sum(paths)),
    executed_file=exec_paths,executed_md5=unname(tools::md5sum(exec_paths))),file.path(run,"generator_checksums.csv"))
  for(i in seq_len(nrow(plan))) {
    s<-plan$script[i];message(sprintf("[%02d/%02d] %s  |  %s",i,nrow(plan),s,plan$ids[i]))
    if(s=="analyse_background_covid.R"&&reuse_background&&bg_ok()){
      state$status[i]<-"REUSED";state$message[i]<-paste("Validated existing site_weekly.csv, MD5",unname(tools::md5sum(bg)));save_state();next
    }
    logfile<-file.path(logs,sprintf("%02d_%s.log",i,tools::file_path_sans_ext(s)))
    state$status[i]<-"RUNNING";state$started[i]<-format(Sys.time(),"%FT%T%z");state$log[i]<-logfile;save_state()
    config<-file.path(logs,sprintf("%02d_config.rds",i));driver<-file.path(logs,sprintf("%02d_run.R",i))
    this_opts<-opts
    if(s=="analyse_audit_tables.R")this_opts$phmpt.audit_table_ids<-if(plan$ids[i]=="V19")"V19" else c("V08","V09")
    saveRDS(list(root=root,script=exec_paths[i],opts=this_opts,libs=.libPaths()),config)
    writeLines(c(paste0("cfg<-readRDS(",encodeString(config,quote='"'),")"),
      ".libPaths(cfg$libs); setwd(cfg$root); options(cfg$opts)",
      "source(cfg$script,local=.GlobalEnv,echo=FALSE,encoding='UTF-8')"),driver)
    started<-Sys.time();ec<-tryCatch(system2(rscript,c("--vanilla",shQuote(driver)),stdout=logfile,stderr=logfile,wait=TRUE),
      error=function(e){message(conditionMessage(e));127L})
    state$seconds[i]<-round(as.numeric(difftime(Sys.time(),started,units="secs")),1)
    if(ec!=0L){
      state$status[i]<-"FAILED";state$message[i]<-paste("Rscript returned",ec);save_state()
      tailmsg<-if(file.exists(logfile))tail(readLines(logfile,warn=FALSE),8L) else character()
      stop(s," failed. No new figure set was published.\n",paste(tailmsg,collapse="\n"),"\nFull log: ",logfile,call.=FALSE)
    }
    state$status[i]<-"DONE";save_state()
  }
  mf<-figure_manifest(stage);expected<-sprintf("V%02d",1:49)
  if(!nrow(mf)||!setequal(expected,unique(mf$visual))||anyDuplicated(mf$file))
    stop("Final figure inventory incomplete or duplicated: ",paste(setdiff(expected,mf$visual),collapse=", "),". Inspect ",run)
  files<-file.path(stage,mf$file)
  if(any(!file.exists(files))||any(file.info(files)$size<40))stop("A figure is missing or empty. Existing collection was not replaced: ",run)
  actual<-unname(tools::md5sum(files));if(any(!is.na(mf$md5)&nzchar(mf$md5)&mf$md5!=actual))stop("A figure checksum changed before publication.")
  mf$md5<-actual;mf<-mf[order(mf$visual,mf$file),];figure_csv(mf,file.path(run,"figure_inventory.csv"))
  flagged<-grepl("UNAVAILABLE|INCOMPLETE|PROVISIONAL|PARTIAL|DIFFER|REFERENCE_ONLY|NO_MATCH",mf$status)
  if(any(flagged)&&!allow_flagged)stop("Flagged output(s): ",paste(unique(mf$visual[flagged]),collapse=", "),
    ". See ",file.path(run,"figure_inventory.csv"),". Previous audit_assets remains active. Set options(phmpt.figures_allow_flagged=TRUE) only to publish explicitly flagged figures.")
  dest<-file.path(root,"audit_assets");figure_publish(stage,dest,run)
  message("Published V01–V49 in ",dest,". Status, checksums and logs: ",run)
  invisible(list(directory=dest,run=run,inventory=mf))
}
if(isTRUE(getOption("phmpt.launcher_autorun",TRUE)))run_audit_visualizations()
