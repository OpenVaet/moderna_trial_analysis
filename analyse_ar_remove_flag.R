# Moderna P301: compare the AR Remove Flag in S3 and S27. Base R; no downloads.
# In RStudio set the project root as your working directory and run the whole script.
# Prepare ar_flag_inputs/ as described in docs/audit-reproduction.md (or use the
# companion extracts). Alternative: options(phmpt.ar_inputs="...",phmpt.ar_output="...").
analyse_ar_remove_flag<-function(input=getOption("phmpt.ar_inputs","ar_flag_inputs"),
                                  out=getOption("phmpt.ar_output","ar_flag_results")){
 read<-function(f)read.csv(file.path(input,paste0(f,".csv")),colClasses="character",check.names=FALSE,na.strings=NULL,encoding="UTF-8")
 if(!file.exists(file.path(input,"flag_parent_crosswalk.csv")))stop("Missing AR replay inputs: see docs/audit-reproduction.md or set phmpt.ar_inputs.")
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 write<-function(x,n)write.csv(x,file.path(out,paste0(n,".csv")),row.names=FALSE,na="",fileEncoding="UTF-8")
 key<-function(a,b)paste(a,as.numeric(b),sep=":"); isdate<-function(s)grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}",s)
 day<-function(s)as.Date(ifelse(isdate(s),substr(s,1,10),NA),format="%Y-%m-%d")
 norm<-function(s)toupper(trimws(gsub("[[:space:]]+"," ",s)))
 arms<-read("arms");stopifnot(!anyDuplicated(arms$subject));lk<-read("flag_parent_crosswalk")
 lk$oldkey<-key(lk$subject,lk$old_seq);lk$newkey<-key(lk$subject,lk$new_seq);stopifnot(!anyDuplicated(lk$oldkey),!anyDuplicated(lk$newkey))
 all<-do.call(rbind,lapply(c("S3","S27"),function(v){
  x<-read(paste0(v,"_flags"));x$source<-v;x$arm<-arms$arm[match(x$subject,arms$subject)]
  stopifnot(!anyNA(x$arm),all(x$arm%in%c("mRNA-1273","Placebo")),all(x$QLABEL=="AR Remove Flag"&x$QORIG=="Derived"),!anyDuplicated(key(x$subject,x$seq)));x}))
 matchsets<-lapply(c("S3","S27"),function(v){x<-all[all$source==v,];key(x$subject,x$seq)})
 stopifnot(all(matchsets[[1]]%in%lk$oldkey),all(matchsets[[2]]%in%lk$newkey))
 lk$s3<-lk$oldkey%in%matchsets[[1]];lk$s27<-lk$newkey%in%matchsets[[2]]
 counts<-function(x,v){do.call(rbind,lapply(c("mRNA-1273","Placebo"),function(arm){z<-x[x$arm==arm,,drop=FALSE];data.frame(snapshot=v,arm=arm,records=nrow(z),participants=length(unique(z$subject)))}))}
 adds<-lk[!lk$s3&lk$s27,];loss<-lk[lk$s3&!lk$s27,];stopifnot(nrow(adds)>0)
 tab<-do.call(rbind,c(lapply(c("S3","S27"),function(v)counts(all[all$source==v,],v)),list(counts(adds,"ADDED"),counts(loss,"LOST"))));write(tab,"arm_counts");print(tab,row.names=FALSE)
 data<-lapply(c("S3_focus_subject_AE","S27_focus_subject_AE","S3_focus_SUPPAE","S27_focus_SUPPAE","S27_CE","S27_SUPPCE","S27_FACE","S3_EX","S27_lookup"),read)
 names(data)<-c("old","ae","so","sn","ce","sc","face","ex","lookup")
 bind<-function(x,col1,col2)split(seq_len(nrow(x)),key(x[[col1]],x[[col2]]))
 ai<-bind(data$ae,"SUBJID","AESEQ");bi<-bind(data$old,"SUBJID","AESEQ");stopifnot(max(lengths(ai))==1,max(lengths(bi))==1)
 lkkey<-key(data$lookup$subject,data$lookup$ae_record_number);stopifnot(!anyDuplicated(lkkey))
 yes_categories<-c("DOES NOT INTERFERE WITH ACTIVITY","REPEATED USE OF OVER-THE-COUNTER PAIN RELIEVER > 24 HOURS OR INTERFERES WITH SOME ACTIVITY","ANY USE OF PRESCRIPTION PAIN RELIEVER OR PREVENTS DAILY ACTIVITY")
 status<-function(x){if(!nrow(x))return(logical());(x$FATESTCD=="OCCUR"&x$FASTRESC=="Y")|(x$FATESTCD=="SEV"&(x$FASTRESC%in%yes_categories|grepl("^GRADE [1-4]",x$FASTRESC)))}
 qvals<-function(x,id,seq,q){unique(x$QVAL[x$USUBJID==id&as.numeric(x$IDVARVAL)==as.numeric(seq)&x$QNAM==q])}
 stamp<-function(s){if(!length(s)||!grepl("T[0-9]{2}:[0-9]{2}",s))return(as.POSIXct(NA,tz="UTC"));as.POSIXct(sub("T"," ",substr(s,1,19),fixed=TRUE),tz="UTC")}
 ledger<-do.call(rbind,lapply(seq_len(nrow(adds)),function(i){
  k<-adds[i,];a<-data$ae[ai[[k$newkey]],];b<-data$old[bi[[k$oldkey]],];stopifnot(nrow(a)==1,nrow(b)==1)
  id<-a$USUBJID
  dates<-data$ex[data$ex$USUBJID==id&!is.na(day(data$ex$EXSTDTC)) & day(data$ex$EXSTDTC)<=day(a$AESTDTC),]
  if(!nrow(dates))stop("No prior dose for ",id);dates<-dates[order(dates$EXSTDTC,decreasing=TRUE),];dose<-dates[1,]
  ref<-switch(dose$VISIT,"Visit 1 Day 1"="DOSE 1","Visit 2 Day 29"="DOSE 2",stop("Check dose reference for ",id))
  num<-sub("^AE-","",a$AESPID);ss<-data$lookup[match(key(k$subject,num),lkkey),]
  stopifnot(!is.na(ss$subject),norm(ss$reported_term)==norm(a$AETERM),norm(ss$pt)==norm(a$AEDECOD))
  ce<-data$ce[data$ce$USUBJID==id&data$ce$CELNKGRP==a$AELNKGRP&data$ce$CETPTREF==ref&data$ce$CECAT=="REACTOGENICITY"&norm(data$ce$CETERM)==norm(ss$solicited_term),]
  stopifnot(nrow(ce)<=1);cev<-if(nrow(ce))qvals(data$sc,id,ce$CESEQ,"CEEVAL")else character();refcheck<-if(nrow(ce))qvals(data$sc,id,ce$CESEQ,"AESPIDX")else character()
  fa<-data$face[data$face$USUBJID==id&data$face$FACAT=="REACTOGENICITY"&data$face$FATPTREF==ref&data$face$FALNKGRP==a$AELNKGRP&norm(data$face$FAOBJ)==norm(ss$solicited_term),];pos<-fa[status(fa),]
  at_home<-pos$FAEVAL=="STUDY SUBJECT"&!grepl("30 MINUTES",pos$FATPT,fixed=TRUE);pd<-day(pos$FADTC[at_home]);ig<-day(pos$FADTC[pos$FAEVAL=="INVESTIGATOR"]);on<-day(a$AESTDTC);en<-day(a$AEENDTC);base<-day(dose$EXSTDTC)
  other<-data$ae[data$ae$USUBJID==id&data$ae$AELNKGRP==a$AELNKGRP&data$ae$AESEQ!=a$AESEQ,,drop=FALSE]
  ob<-qvals(data$so,id,b$AESEQ,"AESOFL");nb<-qvals(data$sn,id,a$AESEQ,"AESOFL")
  val<-function(z)if(length(z))paste(z,collapse=";")else ""
  data.frame(subject=k$subject,arm=k$arm,S3_seq=b$AESEQ,S27_seq=a$AESEQ,AE_spid=a$AESPID,verbatim=a$AETERM,PT=a$AEDECOD,serious=a$AESER,S3_solicited_CRFS=val(ob),S27_solicited_CRFS=val(nb),S3_link_group=b$AELNKGRP,S27_link_group=a$AELNKGRP,category=a$AECAT,subcategory=a$AESCAT,
   lookup_symptom=ss$solicited_term,lookup_flag=ss$solicited_flag,dose_reference=ref,dose_given=dose$EXTRT,dose_timestamp=dose$EXSTDTC,ae_onset=a$AESTDTC,ae_end=a$AEENDTC,
   onset_days=as.integer(on-base),end_days=as.integer(en-base),measured_end_hours=round(as.numeric(difftime(stamp(a$AEENDTC),stamp(dose$EXSTDTC),units="hours")),3),
   CE_rows=nrow(ce),CE_seq=val(ce$CESEQ),CE_positive=any(ce$CEOCCUR=="Y"),CE_sponsor_id_matches=a$AESPID%in%refcheck,CE_evaluator=val(cev),CE_onset=val(ce$CESTDTC),CE_end=val(ce$CEENDTC),CE_ongoing_day7=any(ce$CEENRTPT=="ONGOING"&ce$CEENTPT=="DAY 7"),other_AE_same_group=nrow(other),
   diary_positive_same_dose=length(pd)>0,diary_positive_on_onset=any(pd==on,na.rm=TRUE),diary_positive_during_AE=any(pd>=on&pd<=en,na.rm=TRUE),diary_positive_dates=val(sort(unique(as.character(pd)))),
   investigator_positive_on_onset=any(ig==on,na.rm=TRUE))
 }))
 write(ledger,"new_flag_case_ledger")
 counts2<-function(cols){aggregate(list(n=rep(1,nrow(ledger))),ledger[cols],sum)}
 write(counts2(c("arm","PT","lookup_symptom")),"terms_by_arm");write(counts2(c("arm","onset_days")),"onset_by_arm");write(counts2(c("arm","end_days","CE_ongoing_day7")),"end_by_arm")
 measures<-lapply(c("CE_positive","CE_sponsor_id_matches","diary_positive_same_dose","diary_positive_on_onset","diary_positive_during_AE","investigator_positive_on_onset","CE_ongoing_day7"),function(v)do.call(rbind,lapply(c("mRNA-1273","Placebo"),function(a){s<-ledger[ledger$arm==a,];data.frame(measure=v,arm=a,n=sum(s[[v]]),total=nrow(s))})))
 write(do.call(rbind,measures),"record_linkage_counts")
 print(with(ledger,table(lookup_symptom,arm)));print(with(ledger,table(onset_days,arm)));invisible(ledger)
}
if(isTRUE(getOption("phmpt.autorun",TRUE)))analyse_ar_remove_flag()
