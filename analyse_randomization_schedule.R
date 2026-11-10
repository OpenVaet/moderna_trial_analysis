# Standalone RStudio: extract the whole ZIP and set the working directory to its
# extracted folder. Select/run this script. No R packages required for replay.
# To re-extract primary inputs, see README.txt. Input AS USED, planned treatment and
# actual treatment are different variables. No missing screening ID is assigned an arm.
analyse_randomization_schedule <- function(root=getOption("phmpt.randomization.root",getwd()),
  pdf=getOption("phmpt.randomization.pdf",""), eua=getOption("phmpt.randomization.eua",""),
  bla=getOption("phmpt.randomization.bla",""), gaps=getOption("phmpt.randomization.gaps",""),
  out=getOption("phmpt.randomization.out",file.path(root,"results"))) {
 ref<-file.path(root,"reference"); if (!dir.exists(out)) dir.create(out,recursive=TRUE)
 rd<-function(path)read.csv(path,check.names=FALSE,colClasses="character",na.strings=character(),fileEncoding="UTF-8")
 wr<-function(x,name){write.csv(x,file.path(out,paste0(name,".csv")),row.names=FALSE,na="")}
 if (nzchar(pdf)) {
  exe<-Sys.which("pdftotext");if(!nzchar(exe))stop("Install Poppler pdftotext or leave pdf unset to replay the extracted schedule.")
  txt<-tempfile(fileext=".txt");on.exit(unlink(txt),add=TRUE)
  if(system2(exe,shQuote(c("-layout",normalizePath(pdf),txt)))!=0)stop("PDF extraction failed")
  t<-readLines(txt,warn=FALSE,encoding="UTF-8");pat<-"^\\s*([0-9]+)\\s+([0-9]+)\\s+([0-9]+)\\s+([0-9]{6})\\s+([234])\\s+.+?\\s+([12])\\s+(mRNA1273|Placebo)\\s*$"
  use<-grepl(pat,t,perl=TRUE);pageno<-1L+cumsum(lengths(regmatches(t,gregexpr("\f",t,fixed=TRUE))))
  parts<-regmatches(t[use],regexec(pat,t[use],perl=TRUE)); if(!length(parts))stop("No schedule rows")
  s<-as.data.frame(do.call(rbind,lapply(parts,function(q)q[-1L])),stringsAsFactors=FALSE)
  names(s)<-c("record","block","size","RANDNUM","sched_stratum","code","sched_arm");s$pdf_page<-pageno[use]
  s$sched_arm[s$sched_arm=="mRNA1273"]<-"mRNA-1273"
 } else s<-rd(file.path(ref,"P301_actual_schedule.csv.gz"))
 for(k in c("record","block","size","RANDNUM","sched_stratum","code","pdf_page"))s[[k]]<-as.integer(s[[k]])
 s<-s[order(s$record),];stopifnot(nrow(s)==126000L,identical(s$record,1:126000),identical(s$RANDNUM,100001:226000),all(s$size==6L),all(table(s$block,s$sched_arm)==3L),all(s$block==(s$record-1L)%/%6L+1L),all(s$pdf_page==2L+(s$record-1L)%/%26L),all(s$code==ifelse(s$sched_arm=="mRNA-1273",1L,2L)),all(s$sched_stratum==2L+(s$record-1L)%/%42000L))
 s$slot<-(s$record-1L)%%6L+1L
 g<-rd(if(nzchar(gaps))gaps else file.path(ref,"screening_gaps.csv")); gu<-g[g$evidence=="Untraced",];splitdate<-min(gu$inferred_date[gu$inferred_date!=""])
 if(!length(splitdate)||is.na(splitdate))stop("No dated untraced positions")
 fs<-c("SUBJID","SITEID","RANDFL","RANDNUM","RANDDT","TRT01P","TRT01A","STRATAR","STRATARN","STRATAV","STRATAVN","AGE")
 loadad<-function(path,rel){
  if(!nzchar(path))path<-file.path(ref,paste0(rel,"_ADSL.csv.gz"))
  if(grepl("[.]jmp$",path,ignore.case=TRUE)){env<-new.env(parent=environment());sys.source(file.path(root,"read_source_jmp.R"),envir=env);x<-env$crf_read_jmp(path,fs)}
  else x<-rd(path)
  if(!all(fs%in%names(x)))stop(rel,": missing ",paste(setdiff(fs,names(x)),collapse=", "))
  x<-x[,fs]; x$release<-rel;stopifnot(!anyDuplicated(x$SUBJID),all(x$RANDNUM[x$RANDFL=="Y"]!=""))
  x<-x[x$RANDNUM!="",];stopifnot(all(grepl("^[0-9]{6}$",x$RANDNUM)),!anyDuplicated(x$RANDNUM))
  x$RANDNUM<-as.integer(x$RANDNUM);j<-match(x$RANDNUM,s$RANDNUM)
  x<-cbind(x,s[j,c("block","slot","sched_stratum","sched_arm","pdf_page")]);x$on_schedule<-!is.na(j)
  # File-specific coding: schedule 2/3/4 = ADSL 1/2/3. Check against text too.
  labels<-c(">=18 and <65 years and not at risk",">=18 and <65 years and at risk",">=65 years")
  idx<-match(tolower(trimws(x$STRATAR)),labels)
  x$planned_match<-x$on_schedule & x$TRT01P==x$sched_arm
  x$stratum_match<-x$on_schedule & !is.na(idx) & idx+1L==x$sched_stratum & idx==as.integer(x$STRATARN)
  x$verified_stratum_diff<-ifelse(x$STRATAVN!=""&x$STRATARN!="",x$STRATAVN!=x$STRATARN,NA)
  x$actual_diff<-ifelse(x$TRT01A=="",NA,x$TRT01A!=x$TRT01P)
  x$period<-ifelse(x$RANDDT=="","No date",ifelse(x$RANDDT<splitdate,"Before first gap proxy","On/after first gap proxy"));x
 }
 x<-rbind(loadad(eua,"EUA"),loadad(bla,"BLA"));keep<-x$RANDFL=="Y";y<-x[keep,];wr(x,"subject_schedule_crosswalk");wr(y[y$verified_stratum_diff%in%TRUE,],"verification_stratum_changes");wr(y[y$actual_diff%in%TRUE,],"actual_vs_planned_opposite");wr(x[x$RANDFL!="Y",],"numbered_rows_outside_randomized_set")
 counts<-function(q)data.frame(release=q$release[1],period=q$period[1],people=nrow(q),no_schedule=sum(!q$on_schedule),planned_mismatches=sum(!q$planned_match),as_used_stratum_mismatches=sum(!q$stratum_match),verified_stratum_changes=sum(q$verified_stratum_diff,na.rm=TRUE),actual_opposite=sum(q$actual_diff,na.rm=TRUE),blank_actual=sum(is.na(q$actual_diff)))
 res<-do.call(rbind,lapply(split(y,list(y$release,y$period),drop=TRUE),counts));wr(res,"by_randomization_date")
 # Site-specific timing is also a proxy, not a verified allocation timestamp.
 first_site<-tapply(gu$inferred_date[gu$inferred_date!=""],gu$SITEID[gu$inferred_date!=""],min)
 sitecut<-unname(first_site[y$SITEID]);yp<-y
 yp$period<-ifelse(is.na(sitecut),"No untraced screening gap at site",ifelse(yp$RANDDT<sitecut,"Before site gap proxy","On/after site gap proxy"))
 bysite<-do.call(rbind,lapply(split(yp,interaction(yp$release,yp$SITEID,yp$period,drop=TRUE)),function(q)cbind(SITEID=q$SITEID[1],counts(q))))
 wr(bysite,"site_specific_gap_date_checks")
 yd<-y;yd$period<-yd$RANDDT
 daily<-do.call(rbind,lapply(split(yd,interaction(yd$release,yd$period,drop=TRUE)),counts));names(daily)[2]<-"RANDDT";wr(daily,"daily_checks")
 holes<-blocks<-inv<-list()
 for (rel in c("EUA","BLA")) for (scope in c("randomized","all_numbered")) {
  w<-x[x$release==rel & (scope=="all_numbered"|x$RANDFL=="Y") & x$on_schedule,];w<-w[order(w$RANDNUM),];b<-split(w,w$block)
  blockstats<-do.call(rbind,lapply(b,function(a)data.frame(SITEID=paste(unique(a$SITEID),collapse=";"),stratum=a$sched_stratum[1],block=a$block[1],sites=length(unique(a$SITEID)),subjects=nrow(a),last_slot=max(a$slot))))
  blockstats$terminal<-with(blockstats,ave(block,SITEID,stratum,FUN=max)==block)
  blocks[[length(blocks)+1]]<-cbind(release=rel,scope=scope,blockstats)
  for(a in b) {
   o<-setdiff(seq_len(max(a$slot)),a$slot); if(length(o)){
    nn<-a$RANDNUM[1]-a$slot[1]+o; holes[[length(holes)+1]]<-data.frame(release=rel,scope=scope,SITEID=a$SITEID[1],block=a$block[1],slot=o,RANDNUM=nn,scheduled_arm=s$sched_arm[match(nn,s$RANDNUM)])
   }
   ix<-which(diff(as.numeric(as.Date(a$RANDDT)))<0)+1L
   if(length(ix))inv[[length(inv)+1]]<-data.frame(release=rel,scope=scope,site=a$SITEID[ix],block=a$block[ix],previous_id=a$SUBJID[ix-1],previous_number=a$RANDNUM[ix-1],previous_date=a$RANDDT[ix-1],next_id=a$SUBJID[ix],next_number=a$RANDNUM[ix],next_date=a$RANDDT[ix])
  }
 }
 b<-do.call(rbind,blocks);wr(b,"block_occupancy");hh<-do.call(rbind,holes);if(is.null(hh))hh<-data.frame(RANDNUM=integer())
 for (rel in c("EUA","BLA")){q<-x[x$release==rel,];j<-match(hh$RANDNUM,q$RANDNUM);hh[[paste0(rel,"_subject")]]<-q$SUBJID[j];hh[[paste0(rel,"_RANDFL")]]<-q$RANDFL[j]}
 wr(hh,"unrepresented_slots_before_later_used_slot");wr(if(length(inv))do.call(rbind,inv) else data.frame(block=integer()),"date_order_leads")
 bb<-x[x$release=="BLA",];runs<-gu[!duplicated(gu$run_id),c("run_id","SITEID","run_length","previous_subject_id","next_subject_id","previous_screening_date","next_screening_date","inferred_date")]
 for(side in c("previous","next")){j<-match(runs[[paste0(side,"_subject_id")]],bb$SUBJID);for(f in c("RANDNUM","RANDDT","block","slot","sched_stratum","sched_arm"))runs[[paste0(side,"_",f)]]<-bb[[f]][j]}
 runs$consecutive_in_same_block<-with(runs,previous_block==next_block & previous_slot+1L==next_slot)
 wr(runs,"untraced_runs_and_numbered_neighbors")
 eu<-x[x$release=="EUA",];jg<-match(eu$SUBJID,bb$SUBJID);ac<-data.frame(SUBJID=eu$SUBJID,EUA_number=eu$RANDNUM,BLA_number=bb$RANDNUM[jg],EUA_random=eu$RANDFL,BLA_random=bb$RANDFL[jg])
 wr(ac[is.na(ac$BLA_number)|ac$EUA_number!=ac$BLA_number|ac$EUA_random!=ac$BLA_random,],"shared_roster_changes");wr(bb[!bb$SUBJID%in%eu$SUBJID,],"newly_numbered_in_BLA")
 print(res,row.names=FALSE);message("First untraced screening-date proxy: ",splitdate,". Outputs: ",normalizePath(out));invisible(res)
}
if(isTRUE(getOption("phmpt.autorun",TRUE)))analyse_randomization_schedule()
