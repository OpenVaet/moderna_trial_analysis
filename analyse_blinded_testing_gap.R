# P301: locate the full-blinded diagnostic-PCR coverage difference in v0.30.
# In RStudio run this file from the project root; optionally set:
# options(phmpt.gap_lab=".../unpacked_strict_M3_pcr_export",phmpt.gap_adsl=".../BLA_ADSL_selected.csv")
# Output in analysis/blinded_testing_gap/. No packages required.
analyse_blinded_testing_gap<-function(root=getOption("phmpt.root",getwd()),lab=getOption("phmpt.gap_lab",""),
   adsl=getOption("phmpt.gap_adsl",file.path(root,"reference_v1_12","BLA_ADSL_selected.csv")),
   out=getOption("phmpt.gap_out",file.path(root,"analysis","blinded_testing_gap"))){
  rd<-function(f)read.csv(f,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
  if(!nzchar(lab)){
    f<-list.files(file.path(root,"analysis"),"^M3_mapping_version[.]txt$",full.names=TRUE,recursive=TRUE)
    labs<-dirname(f);labs<-labs[file.exists(file.path(labs,"Prod","pcr_record_audit.csv"))]
    if(length(labs)!=1L)stop("Set phmpt.gap_lab to one unpacked, strict-mapping PCR-by-lab export.")
    lab<-labs
  }
  stopifnot(identical(trimws(readLines(file.path(lab,"M3_mapping_version.txt"),n=1L)), "1.1.0-M3"))
  folder<-file.path(lab,"symptom_links","Prod_BLA_SNAPSHOT_STRICT_BLINDED")
  paths<-c(anchors=file.path(folder,"symptom_evidence_dates.csv"),links=file.path(folder,"symptom_date_to_pcr_links.csv"),
    tests=file.path(lab,"Prod","pcr_record_audit.csv"),adsl=adsl)
  a<-rd(paths["anchors"]);l<-rd(paths["links"]);t<-rd(paths["tests"]);d<-rd(adsl)
  tf<-function(x){x<-toupper(as.character(x));stopifnot(!anyNA(x),all(x%in%c("TRUE","FALSE")));x=="TRUE"}
  l<-l[tf(l$diagnostic),];key<-function(s,r)paste(s,r,sep=":");k<-key(t$subject,t$record_id)
  j<-match(key(l$subject,l$record_id),k);ix<-match(l$anchor_id,a$anchor_id)
  stopifnot(!anyDuplicated(a$anchor_id),!anyDuplicated(k),!anyNA(j),!anyNA(ix),all(a$subject[ix]==l$subject),all(tf(t$usable_dated[j])))
  dt<-as.Date(a$evidence_date);ld<-as.Date(l$test_date)
  stopifnot(all(ld>=dt[ix]),all(ld<=dt[ix]+7),all(ld<=as.Date(a$link_end[ix])))
  a$hit<-a$anchor_id%in%l$anchor_id;a$pos<-a$anchor_id%in%l$anchor_id[tf(t$positive[j])&!tf(t$result_conflict[j])]
  id<-match(a$subject,d$SUBJID);stopifnot(!anyNA(id),!anyDuplicated(d$SUBJID))
  d1<-as.numeric(dt-as.Date(d$TR01SDT[id]));d2<-as.numeric(dt-as.Date(d$DOSE2DT[id]))
  e1<-!is.na(d1)&d1%in%0:6;e2<-!is.na(d2)&d2%in%0:6;a$early<-e1|e2
  stopifnot(all(tf(a$early_postdose)==a$early),!any(e1&e2),all(!is.na(d1)&d1>=0))
  a$membership<-d$PPROTFL[id]=="Y" & toupper(d$SCOV2BL[id])=="NEGATIVE";a$membership[is.na(a$membership)]<-FALSE
  a$late<-!is.na(d2)&d2>=14
  eff<-as.Date(d$EFFCODT[id]);a$eff<-a$membership&a$late&!is.na(eff)&dt<=eff
  stopifnot(all(ld[a$eff[ix]]<=eff[ix[a$eff[ix]]]))
  a$period<-ifelse(e2,"dose2_0_6",ifelse(e1,"dose1_0_6",ifelse(a$late,"from_dose2_day14",
    ifelse(is.na(d2),"no_dose2_other_days","other_days_before_dose2_day14"))))
  stopifnot(all(a$site==substr(a$subject,1L,5L)),all(a$arm%in%c("mRNA-1273","Placebo")),
    !anyDuplicated(unique(a[c("subject","arm")])$subject))
  sites<-sort(unique(a$site));stopifnot(length(sites)==99L)
  arms<-c("mRNA-1273","Placebo");set.seed(30114);B<-4000L;w<-rmultinom(B,length(sites),rep(1/length(sites),length(sites)))
  full<-a[!a$subject%in%a$subject[a$pos],]
  person<-function(x){
    if(!nrow(x))return(data.frame(subject=character(),arm=character(),hit=logical()))
    aggregate(x["hit"],x[c("subject","arm")],any)
  }
  summary<-counts<-draws<-list()
  one<-function(label,x){
    stopifnot(!is.null(x$hit),!anyDuplicated(x$subject))
    i<-match(substr(x$subject,1,5),sites);stopifnot(!anyNA(i))
    cts<-do.call(rbind,lapply(arms,function(arm){
      s<-x$arm==arm;data.frame(row=label,site=sites,arm=arm,n=tabulate(i[s],length(sites)),y=tabulate(i[s&x$hit],length(sites)))
    }))
    wr<-function(v)as.numeric(crossprod(v,w))
    z<-lapply(arms,function(arm){f<-cts[cts$arm==arm,];dn<-wr(f$n);ifelse(dn>0,100*wr(f$y)/dn,NA_real_)})
    y<-vapply(arms,function(arm)sum(cts$y[cts$arm==arm]),integer(1));n<-vapply(arms,function(arm)sum(cts$n[cts$arm==arm]),integer(1))
    g<-z[[1]]-z[[2]];ci<-if(any(is.finite(g)))quantile(g[is.finite(g)],c(.025,.975),names=FALSE) else c(NA,NA)
    pct<-ifelse(n>0,100*y/n,NA_real_)
    summary[[label]]<<-data.frame(row=label,v_y=y[1],v_n=n[1],v_pct=pct[1],p_y=y[2],p_n=n[2],p_pct=pct[2],
      gap=unname(pct[1]-pct[2]),lo=ci[1],hi=ci[2],valid_draws=sum(is.finite(g)),row.names=NULL)
    counts[[label]]<<-cts;draws[[label]]<<-g
  }
  row<-function(label,x)one(label,person(x))
  row("full",full)
  row("early",full[full$early,]);row("outside_early",full[!full$early,])
  for(period in c("dose1_0_6","dose2_0_6","other_days_before_dose2_day14","from_dose2_day14","no_dose2_other_days")){
    x<-full[full$period==period,];row(period,x)
    row(paste0(period,"_pp_baseline_neg"),x[x$membership,]);row(paste0(period,"_outside_pp_baseline_neg"),x[!x$membership,])
  }
  row("pp_baseline_neg_full",full[full$membership,]);row("outside_pp_baseline_neg_full",full[!full$membership,])
  row("pp_baseline_neg_outside_early",full[full$membership&!full$early,])
  row("outside_pp_baseline_neg_outside_early",full[!full$membership&!full$early,])
  row("eff_fixed",full[full$eff,])
  localneg<-function(x)x[!x$subject%in%x$subject[x$pos],]
  row("eff_window_specific",localneg(a[a$eff,]));row("outside_early_window_specific",localneg(a[!a$early,]))
  # Within mutually exclusive timing profiles, keep the ORIGINAL endpoint (tested anywhere in blinded follow-up).
  pr<-person(full);earlyids<-unique(full$subject[full$early]);outids<-unique(full$subject[!full$early])
  pr$timing<-ifelse(!pr$subject%in%earlyids,"outside_only",ifelse(!pr$subject%in%outids,"early_only","both"))
  grps<-c("early_only","outside_only","both")
  for(group in grps)one(group,pr[pr$timing==group,])
  ss<-do.call(rbind,summary); cnt<-do.call(rbind,counts)
  get<-function(ids)ss[match(ids,ss$row),]
  h<-get(c("full","eff_window_specific"));stopifnot(
    identical(h$v_y,c(1557L,1143L)),identical(h$p_y,c(1538L,1176L)),identical(h$v_n,c(3398L,1546L)),identical(h$p_n,c(3061L,1600L)),
    max(abs(h$lo-c(-6.97920305198342,-2.87132982665515)))<1e-8,max(abs(h$hi-c(-1.86002556051843,3.84571405774983)))<1e-8)
  step<-c("full","outside_early","from_dose2_day14","eff_fixed","eff_window_specific")
  change<-do.call(rbind,lapply(2:length(step),function(i){
    from<-step[i-1];to<-step[i];ci<-quantile(draws[[to]]-draws[[from]],c(.025,.975),names=FALSE)
    data.frame(from,to,change_pp=get(to)$gap-get(from)$gap,lo=ci[1],hi=ci[2])}))
  # Symmetric exact decomposition: sum over timing profiles of average risk * difference in shares,
  # plus average share * within-profile difference. Descriptive only; not causal or a split into exclusive infections.
  g<-get(grps);f<-get("full");wv<-g$v_n/f$v_n;wp<-g$p_n/f$p_n
  dec<-data.frame(timing=g$row,v_share=wv,p_share=wp,mix_pp=(wv-wp)*(g$v_pct+g$p_pct)/2,
    within_pp=(wv+wp)/2*(g$v_pct-g$p_pct))
  stopifnot(abs(sum(dec$mix_pp+dec$within_pp)-f$gap)<1e-10,
    all(colSums(g[c("v_n","p_n","v_y","p_y")])==unlist(f[c("v_n","p_n","v_y","p_y")])))
  dir.create(out,showWarnings=FALSE,recursive=TRUE)
  save<-function(x,f)write.csv(x,file.path(out,f),row.names=FALSE,na="")
  save(ss,"coverage.csv");save(cnt,"counts_by_site.csv");save(change,"nested_changes.csv");save(dec,"timing_composition.csv")
  save(data.frame(input=names(paths),file=basename(paths),md5=unname(tools::md5sum(paths))),"inputs.csv")
  message("Results: ",normalizePath(out));invisible(ss)
}
if(isTRUE(getOption("phmpt.autorun",TRUE)))analyse_blinded_testing_gap()
