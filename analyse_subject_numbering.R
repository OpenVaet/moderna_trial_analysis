# Moderna P301: site-numbering audit and V46-V49 | version 1.1.0 | 2026-11-06
# Standalone, base R only. Run the whole script in RStudio from your project root.
# Prerequisite: P301_subject_screening_index.csv from analyse_moderna_demographics.R.
# Optional overrides before running:
# options(phmpt.root="D:/Projects/moderna_trial_analysis")
# options(phmpt.numbering_index="path/to/P301_subject_screening_index.csv")
# options(phmpt.audit_asset_dir="audit_assets", phmpt.audit_dpi=300)
# options(phmpt.autorun=FALSE) loads functions without running them.
# Fixed reference rosters are embedded below: EUA ADSL 25 Nov 2020 and production S3 DM.
# No number is treated as a proven participant unless a checked roster contains it.

ni_csv <- function(d,f)write.csv(d,f,row.names=FALSE,na="",fileEncoding="UTF-8")
ni_read <- function(f)read.csv(f,stringsAsFactors=FALSE,check.names=FALSE,fileEncoding="UTF-8-BOM",na.strings="")
ni_path <- function(root,p)if(grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)",p))p else file.path(root,p)
ni_num <- function(n)format(n,big.mark=",",scientific=FALSE,trim=TRUE)
ni_date <- function(x){
  x<-as.character(x);ok<-!is.na(x)&grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$",x)
  d<-as.Date(rep(NA_character_,length(x)));d[ok]<-suppressWarnings(as.Date(x[ok],format="%Y-%m-%d"))
  d[!is.na(d)&format(d,"%Y-%m-%d")!=x]<-NA;d
}
ni_colors<-c("EUA: not randomized"="#EAA245","BLA DM only"="#9C6BC3","Untraced"="#CB4C51","Mixed evidence"="#637E93")
ni_class<-function(eua,dm)ifelse(eua=="N","EUA: not randomized",ifelse(eua=="Y","EUA: randomized",ifelse(dm,"BLA DM only","Untraced")))
ni_colors<-c(ni_colors,"EUA: randomized"="#D88526")
ni_pick <- function(root,assets,explicit){
  if(nzchar(explicit)){
    f<-ni_path(root,explicit);if(!file.exists(f))stop("Missing explicit index: ",f);return(normalizePath(f,winslash="/"))
  }
  for(d in unique(c(assets,file.path(root,"audit_assets"),root))){
    f<-file.path(d,"P301_subject_screening_index.csv");if(file.exists(f))return(normalizePath(f,winslash="/"))
  }
  candidates<-list.files(file.path(root,"analysis"),pattern="^P301_subject_screening_index[.]csv$",recursive=TRUE,full.names=TRUE)
  candidates<-candidates[!grepl("previous_audit_assets",candidates,fixed=TRUE)]
  if(length(candidates)){
    hashes<-unname(tools::md5sum(candidates))
    if(length(unique(hashes))>1L)stop("Multiple different screening-index exports found. Set options(phmpt.numbering_index=...) to the intended complete export.")
    return(normalizePath(candidates[1],winslash="/"))
  }
  f<-file.path(root,"reference_subject_numbering","P301_subject_screening_index.csv")
  if(file.exists(f)){message("Using supplied frozen screening-index reference: ",f);return(normalizePath(f,winslash="/"))}
  stop("Run analyse_moderna_demographics.R first, or set options(phmpt.numbering_index=...).")
}

ni_reconstruct <- function(z,ref,start=2001L){
  required<-c("SITEID","subject_number","SUBJID","screening_date","population")
  if(!all(required%in%names(z)))stop("Screening index lacks: ",paste(setdiff(required,names(z)),collapse=", "))
  if(!nrow(z)||anyNA(z[c("SITEID","subject_number","SUBJID")])||anyDuplicated(z$SUBJID))stop("Empty, missing, or duplicate observed IDs.")
  z$subject_number<-as.integer(z$subject_number)
  if(anyNA(z$subject_number)||any(z$subject_number<start)||any(!grepl("^US[0-9]{3}$",z$SITEID))||any(z$SUBJID!=paste0(z$SITEID,z$subject_number)))stop("Unexpected site/subject numbering; review before imputing gaps.")
  z$date<-ni_date(z$screening_date);z<-z[order(z$SITEID,z$subject_number),];rownames(z)<-NULL
  all<-runs<-list();rno<-0L
  for(s in unique(z$SITEID)){
    a<-z[z$SITEID==s,];missing<-setdiff(seq.int(start,max(a$subject_number)),a$subject_number)
    if(!length(missing))next
    run<-cumsum(c(TRUE,diff(missing)!=1L))
    for(k in unique(run)){
      n<-missing[run==k];left<-a[a$subject_number<min(n),,drop=FALSE];right<-a[a$subject_number>max(n),,drop=FALSE]
      li<-if(nrow(left))tail(left$subject_number,1) else NA_integer_;ri<-if(nrow(right))right$subject_number[1] else NA_integer_
      ld<-if(nrow(left))tail(left$date,1) else as.Date(NA_character_);rd<-if(nrow(right))right$date[1] else as.Date(NA_character_)
      ds<-c(ld,rd);date<-if(any(!is.na(ds)))min(ds,na.rm=TRUE) else as.Date(NA_character_)
      late<-if(any(!is.na(ds)))max(ds,na.rm=TRUE) else as.Date(NA_character_)
      kind<-if(is.na(li))"LEADING_ASSUMPTION" else "INTERNAL"
      basis<-if(sum(!is.na(ds))==2L)"EARLIER_OF_TWO_NEIGHBORS" else if(sum(!is.na(ds))==1L)"ONE_DATED_NEIGHBOR" else "NO_DATED_NEIGHBOR"
      ids<-paste0(s,n);j<-match(ids,ref$SUBJID);e<-ref$EUA[j];e[is.na(e)|e=="-"]<-"absent";dm<-!is.na(j)&ref$DM[j]
      status<-ni_class(e,dm);rno<-rno+1L;run_id<-sprintf("G%03d",rno)
      d<-data.frame(run_id=run_id,SITEID=s,candidate_id=ids,subject_number=n,kind=kind,run_length=length(n),
        previous_observed=li,next_observed=ri,previous_subject_id=if(is.na(li))""else paste0(s,li),next_subject_id=if(is.na(ri))""else paste0(s,ri),
        previous_screening_date=ld,next_screening_date=rd,inferred_date=date,later_neighbor_date=late,date_basis=basis,
        bracket_days=if(sum(!is.na(ds))==2L)as.integer(late-date)else NA_integer_,
        reversed_neighbor_dates=!is.na(ld)&!is.na(rd)&ld>rd,EUA_RANDFL=e,BLA_DM=dm,evidence=status,stringsAsFactors=FALSE)
      d$week_start<-d$inferred_date-(as.integer(format(d$inferred_date,"%u"))-1L)
      all[[length(all)+1L]]<-d
      run_evidence<-if(length(unique(status))==1L)status[1] else "Mixed evidence"
      runs[[length(runs)+1L]]<-data.frame(run_id=run_id,SITEID=s,kind=kind,first=min(n),last=max(n),length=length(n),evidence=run_evidence,
        inferred_date=date,previous_screening_date=ld,next_screening_date=rd,stringsAsFactors=FALSE)
    }
  }
  if(!length(all))stop("No gaps in the selected index; no missing-identifier figures are required.")
  gaps<-do.call(rbind,all);runs<-do.call(rbind,runs)
  stopifnot(!anyDuplicated(gaps$candidate_id),sum(runs$length)==nrow(gaps),!any(gaps$candidate_id%in%z$SUBJID))
  list(observed=z,gaps=gaps,runs=runs)
}

ni_theme<-function(){par(family="sans",col="#192D3A",col.axis="#3A4A54",col.lab="#192D3A",fg="#667783",cex=1,las=1,bty="n",tcl=-.2)}
ni_draw <- function(out,id,slug,height,draw,dpi){
  stem<-paste0(id,"_",slug);png<-file.path(out,paste0(stem,".png"));pdf<-file.path(out,paste0(stem,".pdf"))
  grDevices::png(png,width=8.5,height=height,units="in",res=dpi,pointsize=12,bg="white",type=if(capabilities("cairo"))"cairo"else getOption("bitmapType"))
  tryCatch(draw(),finally=grDevices::dev.off())
  if(!file.exists(png)||file.info(png)$size<100L)stop("PNG device failed: ",png)
  if(capabilities("cairo"))grDevices::cairo_pdf(pdf,width=8.5,height=height,family="sans",pointsize=12,bg="white")else grDevices::pdf(pdf,width=8.5,height=height,family="Helvetica",pointsize=12,bg="white",useDingbats=FALSE)
  tryCatch(draw(),finally=grDevices::dev.off())
  if(!file.exists(pdf)||file.info(pdf)$size<100L)stop("PDF device failed: ",pdf)
  data.frame(visual=id,file=basename(png),page=1L,pages=1L,status="AVAILABLE",width_px=round(8.5*dpi),height_px=round(height*dpi),dpi=dpi,md5=unname(tools::md5sum(png)),stringsAsFactors=FALSE)
}
ni_wrap<-function(x,n=100)paste(strwrap(x,width=n),collapse="\n")
ni_header<-function(id,title,subtitle){mtext(paste(id,title),side=3,line=4.4,adj=0,cex=1.33,font=2,col="#192D3A");mtext(subtitle,side=3,line=2.95,adj=0,cex=.91,col="#526773")}

analyse_subject_numbering <- function(root=getOption("phmpt.root",getwd()),index_file=getOption("phmpt.numbering_index",""),
  audit_asset_dir=getOption("phmpt.audit_asset_dir","audit_assets"),dpi=getOption("phmpt.audit_dpi",300L)){
  root<-normalizePath(root,winslash="/",mustWork=TRUE);out<-ni_path(root,audit_asset_dir);dir.create(out,recursive=TRUE,showWarnings=FALSE)
  input<-ni_pick(root,out,index_file);z<-ni_read(input)
  if(!all(c("screening_status","screening_source_values")%in%names(z))||!any(z$screening_status=="Recorded"))stop("A complete screening-SV index is required, not an ADSL-only roster.")
  a<-ni_reconstruct(z,ni_reference());z<-a$observed;g<-a$gaps;r<-a$runs
  lev<-names(ni_colors)[names(ni_colors)%in%unique(g$evidence)]
  counts<-table(factor(g$evidence,levels=lev));int<-g$kind=="INTERNAL"
  weeks<-seq(min(g$week_start,na.rm=TRUE),max(g$week_start,na.rm=TRUE),by="week")
  weekly<-expand.grid(week_start=as.character(weeks),evidence=lev,stringsAsFactors=FALSE)
  weekly$positions<-vapply(seq_len(nrow(weekly)),function(i)sum(as.character(g$week_start)==weekly$week_start[i]&g$evidence==weekly$evidence[i],na.rm=TRUE),integer(1))
  weekly$leading_positions<-vapply(seq_len(nrow(weekly)),function(i)sum(as.character(g$week_start)==weekly$week_start[i]&g$evidence==weekly$evidence[i]&!int,na.rm=TRUE),integer(1))
  weekly$one_neighbor_positions<-vapply(seq_len(nrow(weekly)),function(i)sum(as.character(g$week_start)==weekly$week_start[i]&g$evidence==weekly$evidence[i]&g$date_basis=="ONE_DATED_NEIGHBOR",na.rm=TRUE),integer(1))
  weekly_later<-aggregate(rep(1L,nrow(g)),list(week_start=g$later_neighbor_date-(as.integer(format(g$later_neighbor_date,"%u"))-1L),evidence=g$evidence),sum);names(weekly_later)[3]<-"positions"
  ni_csv(g,file.path(out,"V47_numbering_gaps.csv"));ni_csv(r,file.path(out,"V49_gap_runs.csv"));ni_csv(weekly,file.path(out,"V47_weekly.csv"));ni_csv(weekly_later,file.path(out,"V47_later_neighbor_sensitivity.csv"))
  dist<-do.call(rbind,lapply(c("All missing positions","Internal gaps only"),function(scope){
    d<-if(scope=="Internal gaps only")r[r$kind=="INTERNAL",]else r
    x<-aggregate(rep(1L,nrow(d)),list(run_length=d$length,evidence=d$evidence),sum);names(x)[3]<-"runs";x$positions<-x$run_length*x$runs;x$scope<-scope;x
  }))
  ni_csv(dist,file.path(out,"V49_run_distribution.csv"))
  sites<-do.call(rbind,lapply(unique(z$SITEID),function(s){
    d<-g[g$SITEID==s,];rs<-r[r$SITEID==s,];zz<-z[z$SITEID==s,]
    data.frame(SITEID=s,observed=nrow(zz),first_observed=min(zz$subject_number),last_observed=max(zz$subject_number),
      missing_positions=nrow(d),internal=sum(d$kind=="INTERNAL"),assumed_leading=sum(d$kind=="LEADING_ASSUMPTION"),
      EUA_not_randomized=sum(d$EUA_RANDFL=="N"),EUA_randomized=sum(d$EUA_RANDFL=="Y"),DM_only=sum(d$EUA_RANDFL=="absent"&d$BLA_DM),
      untraced=sum(d$evidence=="Untraced"),largest_run=if(nrow(rs))max(rs$length)else 0L,stringsAsFactors=FALSE)
  }))
  ni_csv(sites,file.path(out,"V47_sites.csv"))
  summary<-data.frame(metric=c("observed_ids","observed_sites","randomized_in_index","missing_positions","internal_positions","leading_positions","gap_sites","EUA_not_randomized","EUA_randomized","DM_only","untraced_positions","untraced_sites","two_dated_neighbors","one_dated_neighbor","no_dated_neighbor","reversed_neighbor_dates","same_neighbor_date","internal_runs","all_runs"),
    value=c(nrow(z),length(unique(z$SITEID)),sum(z$population=="Randomized (V05)"),nrow(g),sum(int),sum(!int),length(unique(g$SITEID)),sum(g$EUA_RANDFL=="N"),sum(g$EUA_RANDFL=="Y"),sum(g$EUA_RANDFL=="absent"&g$BLA_DM),sum(g$evidence=="Untraced"),length(unique(g$SITEID[g$evidence=="Untraced"])),sum(g$date_basis=="EARLIER_OF_TWO_NEIGHBORS"),sum(g$date_basis=="ONE_DATED_NEIGHBOR"),sum(g$date_basis=="NO_DATED_NEIGHBOR"),sum(g$reversed_neighbor_dates),sum(g$bracket_days==0,na.rm=TRUE),sum(r$kind=="INTERNAL"),nrow(r)))
  ni_csv(summary,file.path(out,"V47_summary.csv"))
  prov<-data.frame(role=c("Current screening index","Fixed reference EUA ADSL","Fixed reference BLA DM"),
    file=c(normalizePath(input,winslash="/"),"FDA-CBER-2022-1614-6328950-6329713_27073_S5_M5_mrna-1273-p301_A_adsl.jmp","FDA-CBER-2022-1614-3825109-3825910_125752_S3_M5_mrna-1273-p301_S_dm.xpt"),
    md5=c(unname(tools::md5sum(input)),"4a6830e7875e6a853a3123ca46a28e0a","e22e3bfcc741c460652658eb5e761f1b"),
    use=c("ADSL plus screening-SV union; complete original screening dates only","Full embedded membership and RANDFL, 30,571 IDs","Full embedded membership, 32,102 IDs"),stringsAsFactors=FALSE)
  ni_csv(prov,file.path(out,"V47_input_provenance.csv"))
  source_manifest<-file.path(dirname(input),"P301_subject_index_sources.csv")
  if(file.exists(source_manifest))file.copy(source_manifest,file.path(out,"V47_index_source_files.csv"),overwrite=TRUE)
  mf<-list()
  # V46 is calculated from the same reconstructed gaps used by V47-V49.
  eua_n<-sum(g$EUA_RANDFL!="absent");eua_nonrandom<-sum(g$EUA_RANDFL=="N");eua_dm<-sum(g$EUA_RANDFL!="absent"&g$BLA_DM)
  dm_ids<-g$candidate_id[g$EUA_RANDFL=="absent"&g$BLA_DM]
  # Additional DM fields were read directly for these two documented records.
  # Unknown future DM-only candidates are not assigned the same field status.
  known_dm_empty<-c("US3912225","US3982056")
  dm_note<-if(length(dm_ids)&&all(dm_ids%in%known_dm_empty))
    paste0(paste(dm_ids,collapse=" and "),". Planned and actual arms, treatment dates and consent dates are blank.") else
    paste0(if(length(dm_ids))paste(dm_ids,collapse=", ")else"No DM-only identifiers",if(length(dm_ids))". Review the source DM fields for treatment and consent details."else".")
  eua_note<-paste0(if(eua_n>0&&eua_nonrandom==eua_n)"All have RANDFL=N: explicitly nonrandomized."else
    paste0(eua_nonrandom," have RANDFL=N; ",sum(g$EUA_RANDFL=="Y")," have RANDFL=Y.")," ",if(eua_dm==2)"Two"else eua_dm," also appear in BLA DM.")
  untraced<-g$evidence=="Untraced"
  trace<-data.frame(category=c("Found in EUA ADSL","Found only in BLA DM","Found in neither comparison roster"),
    positions=c(eua_n,length(dm_ids),sum(untraced)),
    record_evidence=c(eua_note,dm_note,paste0(sum(untraced&g$kind=="INTERNAL")," internal positions and ",sum(untraced&g$kind=="LEADING_ASSUMPTION")," assumed leading positions; allocation and disposition remain unresolved.")),stringsAsFactors=FALSE)
  ni_csv(trace,file.path(out,"V46_traceability_table.csv"))
  draw38<-function(){
    grid::grid.newpage()
    text_at<-function(txt,x,y,size=11.7,bold=FALSE,col="#192D3A",just=c("left","top"))
      grid::grid.text(txt,x=grid::unit(x,"in"),y=grid::unit(5.15-y,"in"),just=just,
        gp=grid::gpar(fontfamily="sans",fontsize=size,fontface=if(bold)"bold"else"plain",col=col,lineheight=1.2))
    wrap<-function(txt,width,size=11.7,bold=FALSE){
      words<-strsplit(txt,"[[:space:]]+")[[1]];lines<-character();line<-""
      fits<-function(t)grid::convertWidth(grid::grobWidth(grid::textGrob(t,gp=grid::gpar(fontfamily="sans",fontsize=size,fontface=if(bold)"bold"else"plain"))),"in",valueOnly=TRUE)<=width
      for(word in words){candidate<-if(nzchar(line))paste(line,word)else word
        if(nzchar(line)&&!fits(candidate)){lines<-c(lines,line);line<-word}else line<-candidate}
      paste(c(lines,line),collapse="\n")
    }
    text_at("V46  Traceability of missing identifier positions",.35,.28,17,TRUE)
    text_at(paste0(nrow(g)," positions outside the reconstructed screening index"),.35,.72,11.5,col="#526773")
    left<-c(.35,2.48,3.23);width<-c(2.13,.75,4.92);top<-1.15;hh<-.58;rh<-c(1.0,1.03,1.03)
    rect<-function(y,h,fill){grid::grid.rect(x=grid::unit(.35,"in"),y=grid::unit(5.15-y,"in"),width=grid::unit(7.8,"in"),height=grid::unit(h,"in"),just=c("left","top"),gp=grid::gpar(fill=fill,col=NA))}
    rect(top,hh,"#EAF0F3")
    heads<-c("Traceability category","Positions","What the records establish")
    for(j in 1:3)text_at(wrap(heads[j],width[j]-.20,if(j==2)10.4 else 11.5,TRUE),left[j]+.1,top+.14,if(j==2)10.4 else 11.5,TRUE)
    y<-top+hh;fills<-c("#FFE2BA","#E8DAF5","#F8D8DA");stripes<-c("#D78C2C","#8958B2","#C84E55")
    for(i in 1:3){
      rect(y,rh[i],fills[i]);grid::grid.rect(x=grid::unit(.35,"in"),y=grid::unit(5.15-y,"in"),width=grid::unit(.04,"in"),height=grid::unit(rh[i],"in"),just=c("left","top"),gp=grid::gpar(fill=stripes[i],col=NA))
      text_at(wrap(trace$category[i],width[1]-.22,11.8,TRUE),left[1]+.12,y+.19,11.8,TRUE)
      text_at(as.character(trace$positions[i]),left[2]+width[2]/2,y+rh[i]/2,16,TRUE,just="centre")
      text_at(wrap(trace$record_evidence[i],width[3]-.24,11.7),left[3]+.12,y+.16,11.7)
      y<-y+rh[i]
      grid::grid.lines(x=grid::unit(c(.35,8.15),"in"),y=grid::unit(rep(5.15-y,2),"in"),gp=grid::gpar(col="white",lwd=2))
    }
    text_at("Counts describe missing index positions, not a count of missing randomized participants.",.35,4.99,9.4,col="#526773")
  }
  mf[[1]]<-ni_draw(out,"V46","identifier_traceability_table",5.15,draw38,dpi)
  draw39<-function(){
    ni_theme();par(mar=c(7,4.5,8,1.2));m<-matrix(weekly$positions,nrow=length(weeks),ncol=length(lev));m<-t(m)
    totals<-colSums(m);ylim<-max(5,max(totals)*1.18);x<-barplot(m,col=ni_colors[lev],border="white",space=.25,ylim=c(0,ylim),axes=FALSE,names.arg=rep("",length(weeks)),ylab="Identifier positions")
    axis(2);abline(h=0,col="#82909A");axis(1,at=x,labels=format(weeks,"%d %b"),cex.axis=.83,tick=FALSE,gap.axis=-1)
    text(x,totals+.025*ylim,totals,cex=.9,font=2,col="#243B4C")
    for(j in seq_len(ncol(m))){lo<-0;for(i in seq_len(nrow(m))){if(m[i,j]>=3)text(x[j],lo+m[i,j]/2,m[i,j],cex=.83,col="white",font=2);lo<-lo+m[i,j]}}
    ni_header("V47","Timing of missing identifier positions",paste0(nrow(g)," positions | Earlier bordering screening date | Monday-start weeks, 2020"))
    legend("top",inset=c(0,-.15),legend=paste0(lev," (",as.integer(counts),")"),fill=ni_colors[lev],horiz=TRUE,bty="n",cex=.84,xpd=NA)
    mtext("Week beginning",side=1,line=2.1,cex=.9)
    mtext(ni_wrap(paste0(sum(g$date_basis=="EARLIER_OF_TWO_NEIGHBORS")," positions have two dated neighbours; ",sum(g$date_basis=="ONE_DATED_NEIGHBOR")," use one available date, including ",sum(!int)," assumed leading positions.\nTiming is inferred from nearby IDs; these are not documented recruitment dates."),105),side=1,line=4.5,adj=0,cex=.82)
  }
  mf[[2]]<-ni_draw(out,"V47","missing_identifier_timeline",5.9,draw39,dpi)
  # Pick the site with the longest internal run untraced in both references.
  pool<-r[r$kind=="INTERNAL"&r$evidence=="Untraced",];if(!nrow(pool))pool<-r
  chosen<-pool[order(-pool$length,pool$SITEID,pool$first),][1,];site<-chosen$SITEID;ss<-sites[sites$SITEID==site,];zz<-z[z$SITEID==site,];gg<-g[g$SITEID==site,]
  # Full site sequence up to 120 slots; longer sites show a bounded window around the selected run.
  low<-min(zz$subject_number);high<-max(zz$subject_number)
  if(high-low+1L>120L){low<-max(2001L,chosen$first-10L);high<-min(max(zz$subject_number),chosen$last+10L)}
  slots<-seq.int(low,high);jj<-match(paste0(site,slots),gg$candidate_id);state<-ifelse(is.na(jj),"Observed",gg$evidence[jj]);slot<-data.frame(SITEID=site,subject_number=slots,subject_id=paste0(site,slots),evidence=state,stringsAsFactors=FALSE)
  slot$screening_date<-as.character(zz$date[match(slot$subject_id,zz$SUBJID)]);ni_csv(slot,file.path(out,"V48_site_example.csv"))
  draw40<-function(){
    grid::grid.newpage();tx<-function(t,x,y,size=12,bold=FALSE,col="#192D3A",just="left")grid::grid.text(t,x=x,y=y,just=just,gp=grid::gpar(fontfamily="sans",fontsize=size,fontface=if(bold)"bold"else"plain",col=col))
    tx(paste0("V48  ",site,": a block of ",chosen$length," untraced identifiers"),.04,.96,18,TRUE)
    tx(paste0(ss$observed," observed IDs | ",ss$EUA_not_randomized," EUA-confirmed nonrandomized | ",ss$untraced," untraced"),.04,.91,11.5)
    ncol<-10L;nrows<-ceiling(length(slots)/ncol);rh<-min(.052,.57/nrows);top<-.83;cw<-.088
    colors<-c("Observed"="#E6ECEF",ni_colors)
    for(i in seq_along(slots)){
      x<-.06+((i-1L)%%ncol)*cw+cw/2;y<-top-floor((i-1L)/ncol)*rh-rh/2
      grid::grid.roundrect(x,y,width=cw-.006,height=rh-.006,r=grid::unit(.04,"in"),gp=grid::gpar(fill=colors[state[i]],col="white"))
      tx(as.character(slots[i]),x,y,11,TRUE,if(state[i]=="Observed")"#294352"else"white",just="centre")
    }
    tx(paste0("Site ",site,"; suffixes ",low,"-",high,". Each cell represents one number."),.04,.235,10.5)
    leg<-c("Observed in index","EUA: not randomized","Untraced in both references");lc<-c("#E6ECEF",ni_colors[c("EUA: not randomized","Untraced")]);lx<-c(.045,.32,.65)
    for(i in seq_along(leg)){grid::grid.rect(lx[i]+.01,.19,width=.022,height=.02,gp=grid::gpar(fill=lc[i],col=NA));tx(leg[i],lx[i]+.03,.19,9.5)}
    edge<-gg[gg$run_id==chosen$run_id,][1,]
    fmt<-function(d)if(is.na(d))"date unavailable"else format(d,"%d %b %Y")
    tx(paste0(edge$previous_subject_id," screened ",fmt(edge$previous_screening_date),"  |  ",edge$next_subject_id," screened ",fmt(edge$next_screening_date)),.04,.13,10.5)
    tx(paste0("Missing block: ",site,chosen$first,"-",site,chosen$last,". V47 proxy date: ",fmt(edge$inferred_date),"."),.04,.085,10.5,TRUE)
    tx("Numbering order is not necessarily screening-date order.",.04,.04,10,col="#526773")
  }
  mf[[3]]<-ni_draw(out,"V48","site_identifier_example",7.0,draw40,dpi)
  draw41<-function(){
    ni_theme();par(mar=c(5.4,6.6,8,1.4));d<-dist[dist$scope=="All missing positions",];lens<-sort(unique(d$run_length));ls<-names(ni_colors)[names(ni_colors)%in%d$evidence]
    m<-matrix(0L,nrow=length(ls),ncol=length(lens),dimnames=list(ls,lens));for(i in seq_len(nrow(d)))m[d$evidence[i],as.character(d$run_length[i])]<-d$runs[i]
    totals<-colSums(m);lim<-max(totals)*1.16;x<-barplot(m,horiz=TRUE,col=ni_colors[ls],border="white",space=.7,xlim=c(0,lim),axes=FALSE,names.arg=paste(lens,"missing"),cex.names=.95)
    axis(1);text(totals+lim*.02,x,paste0(totals," run",ifelse(totals==1,"","s")),adj=0,cex=.95,font=2)
    for(j in seq_len(ncol(m))){lo<-0;for(i in seq_len(nrow(m))){if(m[i,j]>=8)text(lo+m[i,j]/2,x[j],m[i,j],cex=.95,font=2,col="white");lo<-lo+m[i,j]}}
    ni_header("V49","Distribution of consecutive numbering gaps",paste0(nrow(r)," runs contain ",nrow(g)," missing positions | One run = one uninterrupted block"))
    legend("top",inset=c(0,-.15),legend=ls,fill=ni_colors[ls],horiz=TRUE,bty="n",cex=.84,xpd=NA)
    mtext("Number of runs",side=1,line=2.2,cex=.95)
    lead<-r[r$kind=="LEADING_ASSUMPTION",];di<-table(r$length[r$kind=="INTERNAL"])
    phrase<-paste(paste0(as.integer(di)," run",ifelse(di==1,"","s")," of ",names(di)),collapse=", ")
    mtext(paste0("Includes ",nrow(lead)," assumed leading run",if(nrow(lead)!=1)"s"else"",if(nrow(lead)==1)paste0(" at ",lead$SITEID," (",lead$first,"-",lead$last,")")else"", ".\nInternal-only distribution: ",phrase,"."),side=1,line=4.2,adj=0,cex=.82)
  }
  mf[[4]]<-ni_draw(out,"V49","consecutive_gap_distribution",5.7,draw41,dpi)
  manifest<-do.call(rbind,mf);ni_csv(manifest,file.path(out,"V46_V49_manifest.csv"))
  notes<-c("Moderna P301 numbering audit, V46-V49",paste("Input:",input),
    "V46: evidence classification table. Counts are derived from the current reconstructed gaps. DM-only blank field details are verified for US3912225 and US3982056; new IDs require source review.",
    "Start at 2001 within each site; stop at that site's last observed ID. No trailing gaps are inferred.",
    "The observed index is the union of BLA ADSL and screening-related SV IDs. Reference status is looked up in complete fixed EUA ADSL and BLA DM rosters.",
    "V47: every gap gets the earlier complete Screening date of its two immediately bordering observed IDs. If only one exists it is used and flagged; no wider search or linear interpolation is performed. Missing dates remain missing. Leading positions have only an upper neighbour.",
    "The proxy date describes numerical neighbourhood, not proven recruitment, consent, or randomization. Both original dates, interval widths, reversals and later-neighbour sensitivity counts are exported.",
    "V48: the site with the longest wholly untraced internal run is selected deterministically. All slots are shown for sites up to 120 slots; otherwise a window includes the run and ten surrounding slots where available.",
    "V49: count maximal uninterrupted runs, then show their lengths. All-position and internal-only distributions are both exported. Run counts and position counts are distinct.",
    "Orange is found in EUA ADSL; purple only in BLA DM; red absent from both reference rosters. 'Untraced' never implies absence from every PHMPT file, that a person existed, or that a randomized participant disappeared.",
    "Figures contain no source URLs or footnote markers. Source metadata are in the CSV and accompanying HTML insert.")
  writeLines(notes,file.path(out,"V46_V49_methods.txt"),useBytes=TRUE)
  obsolete<-c("V46_missing_identifier_timeline.png","V46_missing_identifier_timeline.pdf",
    "V47_site_identifier_example.png","V47_site_identifier_example.pdf","V48_consecutive_gap_distribution.png","V48_consecutive_gap_distribution.pdf",
    "V46_V48_manifest.csv","V46_V48_methods.txt","V46_numbering_gaps.csv","V46_weekly.csv","V46_later_neighbor_sensitivity.csv",
    "V46_sites.csv","V46_summary.csv","V46_input_provenance.csv","V46_index_source_files.csv","V47_site_example.csv","V48_gap_runs.csv","V48_run_distribution.csv")
  obsolete<-file.path(out,obsolete);obsolete<-obsolete[file.exists(obsolete)]
  if(length(obsolete)){
    archive<-file.path(out,paste0("previous_numbering_V46_V48_",format(Sys.time(),"%Y%m%d_%H%M%S")))
    dir.create(archive,recursive=TRUE,showWarnings=FALSE)
    if(!all(file.rename(obsolete,file.path(archive,basename(obsolete)))))warning("Some obsolete numbering outputs could not be archived; close them and rerun.")
  }
  message("Wrote V46-V49 PNG and PDF, data CSVs and manifest in ",out)
  message(nrow(g)," positions: ",paste(paste(names(counts),as.integer(counts),sep="="),collapse="; "))
  invisible(list(gaps=g,runs=r,sites=sites,summary=summary,manifest=manifest))
}

ni_reference <- function(){
  compressed_hex <- paste0(c(
    "78da6d9db98e363d929d7d5dcb5f407227cd167a040cd02a43ad31c61ac811204096a4fb873ef2708987994ea0eae0bcb932c960acfffcd7fffe2ffffaf7bffee7fffa3ffff7fffdf5bfff",
    "c71ff12ffff6b7fff8dbdffff98ffff86f7ffbfdfb7ff9c75ffff91f7ffb8fbfffd7fff46fff0ccff3977f1ef747f8f2d7bfffe536e6abc48fc5da1f911278294b585e2a7f448de0d5d405",
    "7975f01a786d8adfbf9e83fdb93e17b3e5b958247e2cf6e79a1d8fe7c66f1b30dfefd7571ccff76be9e2c7627f8ed7ff32bcd07fdbc58fc57c1738477053589eebbfcd0ebcec252c2f872e",
    "0a7955c23c97909b84e59529ecf18a93583cb7def9938865891f8bf577e41a78fe91b03cdf45e0398297b0bc10bac8e41509f0faf88b3c6fea22f3b7b948d8dff667f59407bc3285e5957e",
    "7d25901725c01b63b782579b84e5f5b1ebf8ac5c7f568ecfcaf929c0f35d04f2a204784902bcdc4521af7f0b1ed7e742bf9680fb75a19f2364f2fab715781ff191b0bc3804af2f667eabee",
    "e35b75eb5b4dfc6d4a1296d7e71797f9acfa37e332c695ebdf8ccb3c5e4e12e0f5e3f15b707d1c746179a51faf24f2b204784502bc2a015eeb73139f7d8d1296d7e74957f9fcfadce42a9f",
    "5f9f9bfae46679addf47e3186afd1cad92d7fa248877e997f8778b7909f04217181bbecf75def1786e88425eed0263d2f7efc3738cfb3ec63dc7b8ef63fccf240b5e98c2f2fa9ce3a3072f",
    "0609cbfb336cff9c02cfd9f775cb275e5feed797f1fc7c9f77bbb0bc3eef7ace39be4c015ebfbecae75cbd84e5d57ecd359157b81eb98ff5c8adf5c891e7a6b03c377881bc28015eea02ef",
    "3cf4f71638af853eaf053efbd09f7de0b30ffdd987c4eb4b53585eead7c73922f43922708e087d8e087c97a1bfcb2ec0eb6b683af7e1f73ae8884df163b16bddf21feb96dfeb56012f5609",
    "cb8b43977ac04b7391b3bcfe0c9e4c5e76129697bdd637cb2b53585e9fd79ec267501f09cbab63418be0b5246179adcf890faecf8d49e6c1f53989485e92583acdc0c6f10a7fdbe7ba07d7",
    "ecdc2361797d8c7701de14e00509f0a20478fdfa5c26afcfc5ae92d7b8eefb8f75dfaf753f78f0ba3ed485e5694d2ee455895f8b35adc99617e7c26c8f17bb483c6f0a129697a2165cf0fa",
    "fd26de6f6a5c93fd5e93f9fcb4fef23ef2581bf95cca5cd02caf3609cb6b43794f96e7fb78f146671d589fc73dde475f3a86b0bc302778cbcb4ec2f2f21091bcbe46e5425ee51ae03fd600",
    "bfd68092c92b12e055ae15fe63adf01f6b85d7f66608f0fa391ade876f55c2f2ba7e10f84dbfd619bfd619332786359f9a7d40f8d80784b50f8801bcbe2e74617931697f6979637f69e6f6",
    "8155cd759657a6b0bca1df371eaf1509cb6b557398e1b93591fd582c68d202af700e0b6b0e73093c97252c4ff34b21af72ce096bcef10f786b5f6179da5734f0b40f78c00bd71e3b7ce8ed",
    "61e9ed91c74b8f84e52547fd3e6cfd9ef73be6974cdeadb787a5b7673e973197643e97dc68030873dee8c2f29aa31e1b3ef4d8b0f558479e9700af0b8ff1dcd5da212c6fe8b6fc3ea4b3f2",
    "fbf061e8ac182fbeefadbab0bc387456f2ba3ee4f99ca5db168c973937f1bc2568e201af49585e1d06123ebf313f378c3fdfc68412c94b12e065ce4de1636e0a7b6ec2fd86be068407f71b",
    "fa376875e0f0a103870f1d387ce8c061e9c001f7114292b0bcd0af2514f2aa0478fd3e22c641579d87b0bca13f9bef286e9d35132b123f16eb73a219bb1debf7fb98b13bb03ee1b94a5e93",
    "b0bc61f76909bcfedebab0bc56645330bc6953c0f54d9da68177db19e2d26912ae4f3a48c2f5b9be5773663d8f525386006f2a2ce00509f0a236d596d702f7d871efb123794902bcac7db7",
    "e1ed7df78fc5bc2614cb7353589e1b73099e9fefcfcff3f9f93837bce02509f032f7c471ef89f1de7caa523a2caf4e6179dd96d19504c39b1f279e9f8423cf0509cbd37799c91bdf7925af",
    "6933ba78697d33e6b71d7345e2c762d7fe2d7decdfd2debff11c5da77e7200d6d732bb2f4b1ffbb2b47495f680d79c84e5354ffd2519fde5d76295fbad64f65b3f16cbd423d2871e91961e",
    "112a7863fff147acfd5bdadf2faeefa547a4a547986f307d7c83697e83de81a7b5f68f30e79ddf4726af488057b9ae26b3ae5ade58570b9e9f2f53fc582c53f74e1fba77dab61b8c83ad44",
    "9bf3eafba8e77de435768ddd22cb2430c48fc59c0c01e0458ebffc31fef21c7f768ecd1f736c5e73acb1990f936894b0bcfeacac5d347fd845f3b28b369eb73509c3f372a0e078fdb3c7b7",
    "5fb72feb5cdf70634509cbeb7ad863ec8955db16ec53eaf65744f2926c41e015197fc06b1296d7bf8fc77c1f73ba1902bc2ae38fe58df9c5ecd90716b517026f0af0fab5185fd174910c61",
    "79653a41c04b12e065faf8eaf6f19157337d27f5c3775297efc4e3b97497c6103f161bfb8f04de985f12aee5e54ba8db978077fed20fead60f0a794d63dc621ae33cc7980f1e5cdfcb4651",
    "b78d02efb76f1f86b03ced21f0de5eb68c6a6ce1e05dfb8abaf6151ccf3e2409f032f715f5635f51f7bee2dc479bdfa0f50b35bdb2217e2c96f9dedac77b6bf3bdf5eb333c5d9ff90607d6",
    "afc5ec4d3bd6f7a6ded85407e625c01b9ba4485e929ddff264e7e7f164cf89e45d6b45db6b4501af5609cbebdf87d9bbb8eda33febe0c0a4b7eff1222cca196b797dedeec2f27c809d6660",
    "b2d3f078dd1669ecdec268f7165631ff0d6ccc7589f7913ce635619cd7848d8d44b33cf74cf16331ee0d0676ad5b6e6e0bac0d6f6cf565c34b96e77396a3c0f0bacb6008c30bfd78f1d87c",
    "9d33be841f8b79bd14cb7349c2f2ba2dc8ecc19cf5b3039bc2fe76bccbb3e6094b7a01168b1e363761b4b9094bd06d0736de51c3fd6e3b3f30dac8844d4319785102bc84fda0fbf03bbbed",
    "773e6b9e9b43c0daa98535f87addf1ebf25a2efd45d8dcbf59def0759cbdaeb0a2c9dcf0bc26785c8b1f1f36dfd1f67302cb12f6b7be60de18d8982302de6f0853fc586cd86803795102bc",
    "043b83fbf0b5b9ed6bab7847a1afdda1e11d85be77890ef716fb188f0ef716bb8e1e39c6a39fc2f2fa1c168f7e20ac4880f7e7d9478efb18a7b0bcfe2d44330f75ac8ff168e6a181fd7956",
    "31630cc5ae1fc48c311efbde25e6405e9400efcf338d95e7edf6fb2e2cafefa7a379cec75759880da5f7cc4d63ba7938e7f83de764f20ae762ff3117fb351767f2f2dceb5adeb0979831e9",
    "975dd98c492f17d870728197e0b772c76f759e7df87806613f83f3bd053d92b1365ade581bcdbc1bd6bc5bc8bb628084cd4db5e5f57dca73f629c3ada8350ac77bcd7561cd75e6db0f3b8e",
    "2783e7a9dbba6074a45f8b6519c9172fbef7422ebef742c2127581f8f6d9080b3240585e0b129637ec6167af26d72a7d18c29af6839637f683c7d63cb06e6b76674f2ccc49fc5accc35ee7",
    "acbdcef296bdce9e63d8a4cf7ed56d9b91999bd2c7fa9b967e95f9db1183767c13c2b8071bd8b5077369efc10279117b3061193623f7613372c766e4c9bbde5b9aefcdd8965c7afbc6dcb1",
    "2de13e6e7fbc3087f83597def16bc232d7fdf4b1eea7eddbf6e40509f0a202a92cafef434d0c95b0aa05dcf2c6bd35bcb7192f552d4febb99977d3f61fe13e667c13eec37b271bade5c946",
    "7bc6735218d3d8c0585e72d8cf08f38a5bb23ce9bb3ceff06d738ccbc75c70bfdbdf6379f2f7f078d5c13f2d8cfe696101f1b62ebde36d855539780c6f06d8e2bcd3171dc99b02bc0cdf8e",
    "c298e8db1196e490b1bc3efe02c75f77d7d8f82661fd7895d752a7b0bfad5902bc22015eb769714c8656252caf8fc9e8f1ec63df07445e5f8c41e2c76251e2d7625380d7f5358ebf989c84",
    "e5a5a10fe1f9c53eaf45ce6bb10cc586d7dcc74be47899fa108fd775ccd878bca163b648de14e0fdb98ff4e0b9a4be27e9e2d762b18b4c5e91f8b1d89ff7913cc655eaeb6af21857a9cf75",
    "e9f812dcb6e59eb8cd81adf9f9c76241f381e1693e88e029b6d17cfbf9e3dbcffbdbe7f18a87ff57d874ee5adeface2daf0edbf59913cbd6cd02b128f163b18418b4812906cd8157e9eb70",
    "65f9161ff0e69e13e77de95c65eb5c38deb43156f07c43ccc4c0ae9809610ef1d4c2a2047889eb56f958b7ca5ab7ccf756761c05ef3779ee7fcbc7feb7bce3b406b67c9fd3e7208cb1ce03",
    "bb629ddd0c73b6b6cde15eb96c9bae6cdb269eb3f6b07cceddac69ed7a03bbec7ac28a06bae58df8dd88e72c9f26ef77c66ee139fbe153aabcbe7b3d2a3bff03637c18e88259f7cbcaf588",
    "385e8853585eac887b1dd815f72acc239ed595773cabb0aa6411cbebf716329e4be87b4e932732b02b4f44d8c813f1e40509f0fa391ade5be8736c68786f7d19fba3e6e07e639f2322bfdf",
    "d8bfdff864f2a600afefed1d8f2723028fe73c6d0ae5c3a650b64d21901725c04bb233805725c0eb6b72e07d84292c2f54adbf9677afbf65afbf781f2f1b4559368ac4e3a524015e3f5ee1",
    "7d74fd2f16de47d7ff92b1eb1dbf1ab1dbae5ca72dd21dbf81db7e97e3371056e0af701ffe0ab7fc15c1ecf7e72739c48fc59a82652c6fe866668e1858e237d3b67f3503538ce139477bfb",
    "d00676efdfdadcbf99d8b28eddb165c21c6dd2edc326dd964dda8cbf11de50257e2cd6107b34b02bf64858a31edbde71fac222bf8ff6f17db4f97da488e7d7d59e217e2cf6e75a72c3f172",
    "df73e686e3e5aec315b33f6afd5f27f163b13fe3a5243ce7d2f59cfae0b7f571123f16f35d14f2aa047853585e7f9795efb2f677593d9e41ed6b5ef57806b5bfb776740bbf7d3b67ec0eec",
    "b2cb0b1b7a5305af3409cbab8cd352d816e3b484cdc5d5f07ca8123f166b18a7fec377e2b7efc4e33e421f2fc1e33e421f2f81cfa09bb287b0bc61cf3e3af0c0ae183961236682bca19f56",
    "f2c6bee28c8d8ecdb5e73c2bf78e61f61fb98c7ee7329ef137b02b2ec3cf548df1ed1bde8e2bfdb198c7f73bd269aeef5758836ddd2fdbb089431ee13a571cb25fb914c65fe1fddb5f21ac",
    "21666c6057cc98b0eb9dfb8f77eed73b3fbacfc0d2149637d6370fdebd760bf35a072d6fad8353b71596b0371dd8b5371536f6a63c6f9dfb55cbebe32af1fa924425d6247e0dd6c7780a78",
    "47a9eb0729e01da5fece53c6fb4d3dcf3555625d67cd2786a063b9afbf5d986790fbfbc809ef3cf7b52c67feb6eb8ed99c77dbdbcfbc265763e5f7f111e3efc33b76465892002fc34f3cb0",
    "950bb5ee23bc73a106b672a17e2de6652bb5bc612b3dfb1961639319c0eb7e97e7c4310a4b4a5607af4ad8eb5b712d9637f2cf1dae79e729fd586ce6298137057851fb55cb5bfb55cbd37e",
    "95bc2b5e6f60c3ce9a716f775cbd3f71f578a6d3a68a673a7355cf798f1fa2010b0e3ee681add8801f8b05c43709637c93b0a67160797926c6595e9e2fdef24a95b0bce187383af548c5be",
    "62537cdcba6d218f7b71bfe3953379f7738ef339dbb5227eac1571e520243ce76e06b2f9ab03cb8f0c9c86370c38e1f802fd8ea30d38af7cc767ff31b01424ec6f5384ef58187dc7c2b26a",
    "1d589e6a1d44f0b45fade4cd5a079657e6426f79354858def0a51efbd050619fa875ffc762597b49cb73418bcdc28e5f88d88ac1f8b1d874ae819725c02ba8bf218cf5378435ea1be943df",
    "48734e347ed3e172493307f0c7627faeb91c3bbacf3b76ab5aaca7fb0cf163b01157c5dffa11e87af2397d5ef6a1e30b1c98f2c0ab3d9e77978e993f74cc3c754c63c7f41f764cff61c7f4",
    "e5ed3b1646dfb1b006bfa95e0ffda6c2927474cb1b3abad1918ac23dad9f7d60979f5d58447d1561acafe2cbda9b9ab9ae7ccc7565cf75e45db1bfc2187b3eb0cb8f288cb964fed845c193",
    "7f3fe1bc73fce17e67be69248f3962c20a62c57d79c78afbb262c53dceeb8792ea71bf23ced2e437faf2ce6f14e63977968fb9b3acb9d3ec3feadb77accf99be63bf6d2889bc2b0fdcd7ad",
    "fb54f0963fd9f2863de2f87a0776f97a8545e4368fb0da957bf863b18afc00bfed25278f79602120765058d4b76079e35b38396cc2aa12252c2f67d47511c6ba2ebebdebba08f3d4ebda87",
    "5ed7965e77fc4703bb72dc8565f82184d12f2eac52d76b1fba5e5bba9e19a7ed5dff455894f3c0f24246bd166185df5bfbf8dedafade32eea3ab244358dec8637678bf3be7e2c76251be6d",
    "cb1b1b598ffbf56bf36879572ea3b082fca3815df947c22a6a72f8f6aec9e157ac6ee0fbb873008579f984c19b8e61f098dbe73fec6b7edbd70a9edf6d47173694181eafeba2a1f07865ec",
    "f771bfd14d617871183402c67decdf600c18f7b1efc162e0f1c21496d7f706668f1d4e1cb20376e514861d9b7c749a706a8a35f0243278577c71d8f1c5676d0c1fb5c7c2ae3d9679de2b16",
    "475883cf702c6397cf5058828f4fcb1d7d7cc2127c6803bb7c68c2582f4818f36284315e343cef7c626109be3661d3e1065ec57e6160d77e213ccb77ccebbb63448435f8c686f9eef28d85",
    "5d07edd80064e663deacb00cfbd0c02edf9830af2009cb1bb1b095e71df1020ee78d2318c0e1bc8a493d6b453871dc0d98622acfb5b835ee037fab78f9085e4858b7c2a947e6c11b7ad8d1",
    "2306367ce095bc1a608f1036ed11e035adfb96d766809ee5b519e8055e9500af617d0beebdbe09f3a8af228cf5558425d8320676d932c247ddb2e03ebe41b7be413397b877deb130d63018",
    "d855c320b8ed7bf7e4b16e4a707b5dc57b9b7e76f2467cecd173c2ae4796337905b5bdc2a9ed857130621bfd83eb9bf5b9705eaf6461dcef1d9716dc5ebb13798c470f3b1e3d601cdcb5b8",
    "c2a9c545de8a73b33cc5b544f2686b16367506f08a0478157a44706f3d425845bc80cac1cc822f96d7f71a666f2a8c7b536173830a5e419d1861d45fc2474db1e04cae91e55d7b5d618173",
    "b6fb98b3dd8aebe3fc32e3f5385e6a444e6b70ef9cd6e0763d191e4ff51e787d2d62bf153eea960567f261c163be78702b5ffc49e4b10684b082ba9c03bbea7286531b2d9297b0cf135690",
    "832a8c39a8c2987f3e303941f05c86c1d4f84e84d177222cc9716d79c3ef12797db12077696057eed2c0f24c60b2bc3c05781e763861b4c30d6cc59258deaa390a9e97d20b5e92008f76b8",
    "e0564c87c3f39bb6391c6fc66ae0783d6d64046758de58fc03c6410c0ef193fa372856c3f214abc1e3a54762dabe82dbb1923cde15bf114e8e099e5f4f3919c2f246de49c5fc12fb3718b9",
    "eef7d0491fb9eef730c9212cafcd0047c3536ce383e3a5a749585eff8e12bfa3d4bf23e31b1b989fc2f296e3cbf286cf2bf0bc61040de0f9a5fe3e12df47eaef23155e5ff781a7caebebfe",
    "bcd4787d6d0acbeb7358a6fe92bbfe92a9bfe43e476497c8cb12e0753f1dd7dfdcd7dfccf537fb21c8ebe334735dcd7d5dcd7c7eb93fbfcce797fbf3cb66fdf5ef5883815db106c292127a",
    "2c2fcf82169637ec34272f6660575ecc30078ebc18f34d7fd4e60bbb365fc635df39da61d79b33fb2dffb1dff26bbf55c9bbfc56c20ad7bc8f1a6a61d7503b3ebe81a9eed1035e9b13b7e1",
    "0dc374f0389e0c88668cfb772c8eb0a80042f09af64c96971c7c2c03938f85c7539c60216fd69806af697eb6bcfe9c039f733763585fccc02e5f8cb0a8d878f0666c3c784531ef96376a3f",
    "715cf5700eac017ecff78518e3f58435ad0196b71602cb9348e46589c50bef1ced10de39da61e7a125f22e7f80b08878b3815df16621bc6b3684f0aed910b60fdcec1bb7cfbaf15a2e7bb1",
    "b02863c8e26dbfae396ffc386f9ce775669d89bba667b53c3792fdcd3a13df7945c2e8a31ad8e5a30a71ede932aef9b5df8a6bbf75ec75035b79c29657e85f08f1ed5f1016b437b0bcb537",
    "b0bc7ce9cff1437f1ed82c72667853b72d96a7b84887f7b1ed98bf16cbd463e3d2633d9eb3828d3cdec71da725ccd1c6133f6c3c71d978cc7a14774e8dc3f5c539b9595e7a242c2f39f8a7",
    "c3877f3a6cffb49913b70fd7e8021fb973e1e4ce35cbbbeb10089b45dc2d4fc11901bccba6bf4c16d6a62faca006acb0590316bca6c438cb8b3379d6f2946f90c9635da6b0eb32159ef7b6",
    "4fa6659f7438ef5d1731ec5a4d01cff9652b48cb5670621b437ad71d0c2b0fcde47e858fdcafb073bf12aeef8e371336d6463c97be2c0e6179636dac3cde1563282c20967d98daaf587661",
    "95b6f564e2d72c2f25eadee943f74e4bf776e4751d3379f2fc54982dafef2553c4f853fc2ec741eae320711c74b5db676347df394e46a7c9523987f8b15894f30fbc5917dbf2d22c406179",
    "f29547f2126a3b848fda5461d7a63a397b0353ddfd04ded013cd1c91773f1107de98b33dcefbb2d5e75dfb13c77b7d0b79c7c1f3786321083bb631e45db7a792d724ecf1464dc0c4f35e79",
    "a4c23c6a948de9faaa51262cc0af36b0cbaf266cda3c16af7cd8d6cb876dbd2cdbba79a61f311861c760181b72597132c6865c969f29f178f7182a730c9918c3504c5d35cb53fd834c1e7b",
    "9108a3cf7a6097cf3aecb88c46debd1eed188c07bcbbde7a38397189bc2c011eeb0d87f2ced517c61ae703bb6a9c0b636f13615e02bcc075ab7cac5b3b4f2ef1bcabcea7e5ad3d2278ec13",
    "126cdc087811b5cb8525ea8ee543772cc67f0e1ee3fa84b1164d28bbef48006fcd57e0b1ce58d8f9790de3feae9b37d2e64623a380fbb873ecc2478e5df8c8b10b3bc72ee19b9eb5b3701f",
    "bd9416fc8de5c3df5876ae2a9edfb68f83d750333d9477cd7461957a44f9d023cace67e779af5a5ca1ac3a3b0fdec74bdf284bdff0f82ea7fe8cef5276e080fbd8c2f2465c01af4fb9e1bc",
    "bee1fb8c9c0fee1cbb7072ecf01df554ba212c6fe5d38177d95ecbb2bdba48dee55f2dcbbfcaf960e6dd615cc96ecb71fada63db7cba5f8bb1e6cfc0e2b4e5dadf0e5b6ec6d89d757b782d",
    "39d1275c763d1ef26a920176f176dc97f161d4d53bc4f830eaca8bf6e45d31a461d7ae347ba12af5083ad28919e3f16e9bc2ae176974a4faae771356cd473b8fd78f79bcee3d4e228f756c",
    "c2a9e588fbd839d5bf16bbf6f6f5636f5fb7bf16cffe8ebb16d6503b6660ab768ce5554f5daf7ee87a75eb7a782e2f1bc0aee518717df24b465c5f2f5708ddac1addecd762ac5928ec9a77",
    "eb9a778d2fffe4373ec0ae1843618cef54ba00f3de8531bf31d8fcc65f8b39da54db874db57df8119bb183d8e30d7fb7f16bb4e95bec02bc84b8cdd04cdce6afc5664ef5e4c51ddf74e25a",
    "843589a9030f6c1547b4bc95db6279da9834f0ae986d61cc13898f89ad30c71b3673534323ba77ffc6815df1e3d12dbdf8eccba27befcbe2ae09789e73fc881b892b6ec4d4ed89cef49f01",
    "46bd4ed88ca3b0bc150f6c792b1e18bc19140c5e855e17dd5baf8bab279ca9c11e6defb35f8b79d48491998a356196e9cac623e95fc62345ffce9f8efe9d3f1dfd5a038e1e163ffc3371fb",
    "67ceba10fd7b5d88bba7d459a30676c50a096b88858d1ffd9ea2dfb1b03cde152b14ad0fc8f06ebb59dcbda28ecd3cfa77de53dcfd998e7e15fdbb4e79f46b8de2f5dd359d8435d4d81f21",
    "28578d7d6109fde90676f5a713961187278c7178c22af6f6035b7b7bcb1b7bfb84e777d77589a75f51228fb666610535fb0776d5ec1796b48f5fbc60749a1f8bb18fd8c0ae3e627186b3e1",
    "3b0f1fdff9aa6967f24946599b954ff26331c600c5f08e011256557fd7f246fddde3b78a336cc6eaedc3f576e9edc268cb88c1c46f8097904735b0ab2ea7302f251fbc193c6c799a73701f",
    "77ec87b082188c815d3118c21883218c3118c232eadbc6b0e28b3ddeefcccdade43509cb0b8c251ed8c8cdadb8b71ede3084fd6dd70f52c335a711bbd070cda9cf4df9e84831bef57161d4",
    "c7635cfab823eff2dd099b3628cb930d2a813792a9ce7e26c6a50b9cfd8cb0821ca218df3944c2984314e3db5635b0ab3760dc7990c767ade59db625619e6b40fc5803e2bbe7df285b74d9",
    "a08479e498c4f8ce3189bb36a4c7bdddbe931897efa4e0fa6edb88b00c9d3fc6b7ce1f57bd48cffbb8fbb1c4dd8fe5f82f8531675e1863b6637cdb818525d4298fc71f95815d392671d75e",
    "3cfbb75156fdb271c7b4eb17e378a356b1e3f19c6acb25cbbb6bf9089b9bf4c5dbb98267cf14ad7d7ce9accb3e6efc8863fabafc88c21c7cf9b1bc7b1dc5f2ee7514771e9f59cfcbaef71a",
    "c17389ba72f9d095b70df9c433c4b2e3ae0378a387ef83e3ddbe6d610db5716379d7c68d3b77eee456c56d5330fb94d3e3e1ccb175ed2bccd8f8d8dbc7ddcfc1ec21ea3bce5c58401d9181",
    "5d754484157ebff5e3fbaddbf789f3eefbb5bcd11333e2bcb3af74012f56d41a8a1fb58662ddf98315bc618fcd787ebd2da58d0388f51d07208c3d950776d95e87faa166f38fe5bde6976a",
    "fa2c83c7386e61113ed758df3e576101bd7084cd5e38e055f44a967982bd9285cd465096a7d8e907bc3c03a82d6fc54e8397a8b3d60f9db52e9d95ef63c63af3fa465c86d145ebbb5e69ac",
    "bb8e217957ec96b019c80cdeb43d80c7be55a3bcdbd5b72ad65573d4cc6b75d7c049e079d6798a27d73200bb724c8425d4958ced5d5732b6bd8f6f9637d7ee48aca0b7ac4c701e3104c22a",
    "fa2fc4f6eebf3030f50de2f146df57a39fb6e513693c47aba815364c6657ad3061013949e9d9f9ec8fc5ee1eb9c23ce692655eb4791369db5a8ebf6c60a3f0d9f197096bf04b26f7f64b0e",
    "6cf8254f4cefdaae0f6179953597d3cebdc9b896dbf724ac611e4aee3d0fa5ddeffdcc2f031b73c9d14fd3476e465ab919e65b48dbe6716c0fe9c3f69096edc1e42f24ffce5f10c6beb469",
    "f7a03df692e4df75df8425e44626bff220f9db11db9d13ce91fb0de6632b15d6c5d145f56f92b0bcae8be6f3bd75acf475a11c3b88b05edfec21afdb188b23afbff3726c9bc2b204787fc6",
    "5039f6a6746af92460573d8a14def5285278eb48c29a363096b7fac0589e8ccb01bccbf72e8cfb9914defb99b4ebfb9cd8e41476be7825c6fcf874e227f9db2b6631851db3886b99f3a4b7",
    "bcd9db2492c73e7b296cdf7b03ef8a6d4ce16df34dbb3690e76f2fbb7cda3d7c8f6f5118f3df8455d86807b66cb4f6bc31a01e6d0aef7ab4c23c6cb429bc6db4296c1b2dded1ad870dac3c",
    "884113c61834614d9b38cbab5358ded0d12bdfd115ef93763da346ded5c753d86c5267785bd7fbb198437efcc0aefc78611ebdc05278f7024bbb9f70c07774db0e852519c02c6f2c9011cf",
    "79f6c6c2fdaae7554ae465e87503cb0175755358757573228f7a9db0823cb4b4ed34c73e39b0f18e8e8f79980dd558e05c5f5c7db4237957bf316199e3397e8ce7b8c673c6f5296ea4e1fa",
    "ee3e12c9f67d058fbd97878beeea0b2f8cf546d247ada6b46b351d5bf3c0461ec1b1350bab883518d88a3558df795cb10666de4d6fbb9930dacd8445ecbb0776edbb53dafbee049ed69e02",
    "9e671fb1815d7d1453da7d14237897dd4c5846bf4561ecb798768fd7c4ebbbea0b082bf0e7a58f1a4c29ed758bc7d31ac5e38d98b1c2fbb8d7a3f4eeb798d2ae4bc7e32ddd11bc4be74adb",
    "dfc3777ed904d3897fe67d2cbb14b0865e60c36da7b5b1109b7aac3987d6418f31245b81c7189abd433279acdf2eacc18692d2b229041e6f05415bde322c589e3635bcb7ab9e564a26371c",
    "3cfa3a560939ac97e963bd4c6bbd2c78f6778e534abba644258f7d7406a63e3abc3ead8d3c9e6af8f1b98cfd07c790ab0d79ea03bbf2d48505ae97e963bd4c7b4f52c8a39d35edde2b0ee7",
    "bded2029bded2069f7ffe53c74e7b80b638e7b4aef7cf694de3d19845d7ba6f4ae7793768f61ce39778f4d61456bade569adc5b8925da5f0faee7526ad1c728e839df36d798ad5a8e4cd06",
    "9596d71ee4598fd27f579e754a3b463890c77e2c69f763f1b88fdb2727cc6b91036f16c2b6bc308b095b9efc6f9ebc80fa3ec258df27a5e5a74bbcbe348bd7029bc927f6b743645e4b9eb5",
    "812ca6da40781f33c78ee71d39080d632db404fb7d4aef5afcc3f576c5090af36aca6779aea2d6d0c0962fd0f2fc6cd4075e402e9e30e6e209632edec0ae5c3c6133210fbc20015e549c20",
    "785502bc59f7c8f246dda388b1bb7b0882e71163984ebf9d4c5e9100af4a80d7242c6fe495538f885d8f88d423622eaa696c79a3a631bffd58bc0afc5b5e5f174c9de364fbf280c73e85c2",
    "665c24785380d7afb9f1787d5d888dc76b11bd7ad2e9cb836bbe7360d2c98171e4f55c7317c88b12e00d91c92b4a56b7bc3e9e93e7f1d4fb27935724c0aba8132d6cd689b6bcab4eb430d6",
    "894ef95d8b6b60572d2e6181fa5a7ec7ec086b5cdff2bb5e4bda3eb913fb36b0abaeb3b08c1abfa9ec1abfbb3e75da7918c79f32b02b877ca847ca212f9677e72aa49dab707cef035b3976",
    "bf169b090696371cc3c67658e65ed2f4621acbf1bdee97bdee37f0d2ccbdb1bce4e0d7500b6dfa358445c4b0089b312cbf162bb49fda3878f01a6db465fb2ba2e5cd9cd142deac73025e43",
    "9f9a64fbd458dee85353305e76cc89e58dbcf2c6f3aede6796377a9f193b7f7dc76da66ae2367f0ce61eda18eb878db1ee1ae2e45d793bc23c7a1c0e95fdea7198764cf4f169a6e3d3cce0",
    "45d6725448f3b5d7a86baf6174f4faa1a3d70f1dbd7ee8e875ebe8e75d6edf9323369a701d3fa2b08a5ed3035320a9036f04061a3da7bdfdd3c2e89f1e589d85d88165ee7fdbdcff9abc1d",
    "95dd66dece2ac56de3b8537bc771a7b6ec70c746969fb70f6860970f485893b0bcd5d3dbf2465dd8e3031ad8e503125691139c9f774eb030d62e1f98fce2fbfbc86ed973ce33cda73e4804",
    "76d95584b1266a5e71a0260657e9578cc11596e057cb7ed9452b7e7bc727e6159f68f62479fb35ce7a39b0ab8664debe0973bfe13d867278db507278db50f247ed84fc513b219fda09e0ed",
    "9ee83f168b88a3c8e11d47212ca2b67fb6b5fd2d6f048c27f22edf6c0e269fcef2ae9a2b399afef6ff6e318fd8d0bce3eb8e6e9be3caeb3dbaedc0d283fa1603bbea5b08631dbe31a496ad",
    "65ea11799a5946cc84e12966e2d4f41496d0c74e18fbd8092ba8d336b0f220d640ade4186b206c3655373c25c5069c43f6d8b32ee46d7b3df55e07b6e6921f83f90731e539bd63caf3ae4d",
    "7ffc0603bb6208725a3a9c79cee99d43242c620dcde9bd86e60ffd347fe8a779e9a7c6473fa6ebcb472fcc636e5fd3ba9ddbf52f7b50e4f23147948f39a2ac39e2f81695d2cbba47c23274",
    "cc5cde751cf2d6314fde673e799f1ebc1567f46bb1885e1063c9720e762e6101399ef923c7339f1ccf02de959324ac29df6ff1b67e65be699b87068c794579d72e3f36f35cdf36f37cf4b0",
    "405e82dd3bd76df74ef8edf0dd79fc76ea619518eb8be6938796c00b053ea5815d3ea5bc75b8131b953ff2d0f2c94373e08d9cb31373272ca05e95b082ded579eb7a8dd777f92a85d15799",
    "ebdbf69a4f0c5a22ef9a13eb9a13cd377df4c4079873fccedb8ea7cee4519fcccde852bf166b88a71ed8154f9ddbfea62b79f41f0d4cfe235ef3a5b3e6adb3465ef3e53fcabb86fdd9fbe5",
    "f6cec7c9ed5dcb421863218425f424cca72721cf9b2ae7e2f6311737d31f02bc20011e7bc0e78f5af779d7ba7ff09c77fd26c39b310e786ff2d7ba4c5e411ee9c094478af33a5f11439adb",
    "3b8654d8f4f780379d3ee065c44ce4f68e991056a973b50f9dabedbef5bcbe2bff4d58c03e2fb7f73e4f58e0fcd23ee697f6eef5abb00ad63f1016d003207ff400c8ab078089b557f902c6",
    "da0bcb88d510c6588ddc76ac06cf7bcf6b6dce6b56876b1f3a5c7beb70e5e46e466249e2c762193eeb815d3e6b6109fbb762f76f53ef2ccfbbaf87308ff83f618cff2ba7febd07efeaa551",
    "9edd4be3016fcd57f65a345f15f2aa843dded2c72def9aaf8439cd4d9677d5152fee1d2f50dc3b5ea038d30f1cbcacfbb098ee83bfbdea0509e3bea2b87d1fbcbe6b8e1596e0a31f582e88",
    "2313c6383299811ed8878a7bdb875689773b9f0eec9a4f8badb70e1e6da0631bf43c98ff0676cd7fc5bde7bfe2def35f71eff9af7cd44c2fa7667a01eff2950b9bb9fabf0653306dc36faf",
    "1a2ec5bd6bb814877ee50663ef9ff2515bbd7cd4562fbbb67ac173debe727b0ef9caf17e770e11781eb165c2666c1978d4c3d42278ea61bf169b02bc991f005ec03c59769df72793477dad",
    "38e32bb7bcd5df05bc28ff397805316dc5999836cbbb62da8431a64d5844de6c71266f16c7630c78712606dcf242404f5b615535b12c6fd4c48abcbe95fb005e442df9626bc95b5e627f0d",
    "61b3de8de5e5073172e5d48dafe4cd7af09627db03be99eddf078f35e2cba9118f31b97d09f82d6bea14f7aea9534e2d795e4b7588eb2ba7be7cc479eb0c2e006f16e001af215f63602b5f",
    "c3f29a43fd8872ead0f35dae5a11e02509f0326afe8c34daabe64fd935e71dbe855d076dadc927f7ff7c5bfe5dcb4cd80cf4b2bca14b9975d5bf6d87c2226c1e031bb6ff13ab51668993d1",
    "60d1f086d1a48b5f8b15d4e92d1f757acbaad36b6aab0ff7f2555bbd9cbaba87b7edb6c7de543eeae0965307b78077f5f4293bdebbf01c574efa6ab162fd0bc298ff5666aa9fed1154c2bb",
    "4750d9f1de668d0fef5c4b615593aae1a9bec089672d3364c9c6f1e8d3651c8f30af469686d75b5adabe9625acbe96c77f5e8e2d3711cb123f162bb0e195f8b6e1959d2b6df4c49d9f7c62",
    "8985314ebfc40ffd6adb818f3fa0c4770fb7b26384cdfc1cdfb17465f7643df10c035bb1e2bf168bf0ab95f8f6ab096ba833b15202ad5da59c7c623cfb5d2bc7f29c432dc7818d6236c71e",
    "31b0c85e6ac2185b266c3a9a2d6fd510b2bc15f40d1e63d04a5c6b9ed187e29eef79bc1a38177fd4342ea7a631ce7bf7f028bbceb1c3fbdd7d382c6ff9a7c16b1296b79acd5bde559fa1ec",
    "deb701efedaeff524efde248deecd7015e46dc5789bbbf069ecbee870b5ed2fc0c5e416f8eb27be4165e9f6a40609cde3de684cddaea86a7daea1ccf3b660c3caffebae045ce4df1636e8a",
    "efdab8c2ea6bbeba6bd50973e88b278c7df1063662d5f8cdf410ad217e2d561493657f3b6ce109f3d0dd675ed808ae226f097bbc92d0677e60579f79611e35e84a7cd7a0133682a6cefb38",
    "fea304cc65891f8b15f8b18555f899caf633993936bd6d91c22a6aa296f4ae892aaca23f7b49effeecc20ad7dff4b1fe26936f051e6b9d0a2be84131b0d5abccf21c6b74978f1add25bde3",
    "870676c50f9564fc06e031f6a3ec7adc66dd4fef1c5961ac5b51d23b9f539847efe592debd974b5abd978d9d2bbdfb8309735c7b3ee29acb8e6be6fbbd6bef08637fa672e286711faa237c",
    "6c912b74c8c6c7eadf8c5edd65fb204f7ccec0ae7a2d25af5c94c8df5eb5d1caae596df67979f59e3f3ef581e5877a4efed073f2b623f1bc43f92c3caf941df22ae30a94b278ed0df2fe66",
    "0a79ac5151f2b2df1f7fdec0567ed4afc5d87b65a524d8bc846516865e9c978e14f1dbd718ca6b0cf1f9dd3dd24a5e7b62334fe677dc88b0c4bd647ee7f90b639ebfb0a8a2d786b77b90ad",
    "e7b27dccc616548cafd2f22afb7008637c89b0c6f7f6515ba4ecda22665e2befde3ac202fa320a635f466189f3df876fbb6cdf36eff7ceab149669ab2f6fdfa2b0cbb654def992e5f4f42e",
    "e45509f01a723806163ced3ee55df3a2d8dac7e031ef5dd86c366f798adbc433f59af0f80cee3971c7729a7d777df79529f55d9b4a58a4cfa12eff42e4f1ae1edcc28a8c0d9657e67edaf2",
    "b49f8ee05dfef8b2fcf1d69e5d77fdbf44ec5aa7eb5aa7cd38a8ef5a7f65d7a30d782ed3ee8d6b79d9b8eb878dbb2e1b77e1f1aefa2fa5be6b7497faf68d09a36f6cb9c5ac6f4cff36e465",
    "0d6ce52c2f1db3ae6fcb3cd3ed7b37f697b67d4a67ed1e98a38fe5f8b63d7891b15ba57df858966fdbac3df5c445366057bc63b5fe32f002fcd8f579fbb18579c48cd5e71d3326cc21a66d",
    "98f3af9836611ef5b40676d5d3129610bf21f700e337567a9db5a38fd0b0abaea4b0849ea2d5bdf5f6eade7a7b754b6f3fefbcbaf73b17c6da31f5a33f71ddfd89038fb7e2232c2f720c55",
    "f71e4375d7943db106d5bd630d8445e84875fbda0acf5b1ce286ab7bc70d0b63dc70fde89f5c3ffa27d78ffec975f54f36f9e7d5bdf3cf8565c46d0eec8adb141611e3208cb5b285715f31",
    "b06b5f212c69005bdee51babdb3776f615d5bdf715c266133df066bf72f02a7c63d5bdebeb55b7ebebe17eef75bf6e1f55c0fb983e1bbc0fd9c312deef6dd31216e08ba9eedd7ba01e5f0c",
    "8f77f5ffadbbff6fc173bef3358445e8a2d5195df4d762097e83eadefd538479e4f90feccaf31756d12fab9e386e07ecfe56fdfe56f9db7b3ef56b3e3d3e8ceadf39dac298bb39d213afdc",
    "cdea57ee26cfbbfbc85adeb05f559c77f610c479434dc81d515ae4cc1db1bc913b1271bcdeded5f678adab0e503473d3a99fe3815d35fb85cd02a0e015d43aa8aba68e37f3f86c9535c48f",
    "c51cea7ed48fba1ff5d4fd78c05bf95bbf16f3c869a81f753feaaefb51c85b9b79cbbbf28f6a7cd70caeb6ee2a78538097104f3830e7d023a81e3f04cf7be561086b5c8fd25ba7aebb6e85",
    "9987d2bb9e51ddf528ccbc914c7c9de565fa8fea87fdaa6efb95c3f1eefd9b30d617a81ff5056adab594cff3cbeffc2861cc8f1ad8c88ff2e45d714b035bfa86e585a97480e765c8b33c6d",
    "2e2a7885b914da3e30974258d47ec1f0b45f3836e4baed1b46c7cceffe152b64186b72fe5893f35a9303ae6fee35707d6e35a3b43cc5bf3ce0a5b9305bdeea39647939c3d6523f6c2d75db",
    "5a0acfbb62582c6f2dcce0b16692b0887a37ab34928d2f16c63e442bccc5d6c2d4f6ff411ccac0ae9a0dc22efd397fe8cf79fb4d71bfd3569088b1b683b0c2f9347fcca7a7bf1a8f77c5d5",
    "0b2b524e2c4ff12578f6732ee679af1e4bc292821016afbceb400eecd6e5771e6ee46faf78656195df5bf9f8decafadecc5c5c4c8db21f8b79eac0e543072e3b7789bcfb3b2aeb3b0a993c",
    "eed98555d82c6b79d77412c6dc91fad1cfab96777c582d880f7b36b6368096777f47d33c66e3b44698cb15cf2a2c211fb196773ea2b080589c6afb5bfd1aac3ed431cb878e59968e9971de",
    "db2f292ca8a785e58d5892533b559883bf6c84875c7d2d858de20685bc2a015e53d102cb1b410d05cf59beb6c6e375fd2f9edc61614dfded0d6ff7b737bc595f2082e766a180c5abeffad4",
    "c2667d6af02af50d6bc3b33cc521f378cb4f075ea50da59adc2fcbf301fe0f61d3ff015ea64da19a35def264474ae02d7ba2e58d3927f1fa9407c4eb4b8c1ba9d5e43e8057a87355a37359",
    "9eec0c3cdecadf022f6b5363790ab8e2f1962f66f1b6cdedb9b04b8f6d3bf63c93c7dae5c22a626c0676c5d8d49d3f6374aeb6f4ab13975bdb3bdebb36d357d5f2ae5caddacc3b028f3159",
    "b5bd63b20676c564d5b6ec3995f776db644e7e77248f76e5f6bc9fbd303efbb6e3fecff736b0cb66de9e773f87815dfe41611936b2f6bc6d64edc4e427f232d6dfb663f233cf7be5190a63",
    "5e5b7b8c1de9d7625321b53cada1782e6e359c06af2087bc3d26b765f14e8c7f21368b2afc184c7bb5089e4bb06336f78e3f681fb594dba9a5ecc1bbf6ddcdbdf7ddcdbdf7dd6df7c60abc",
    "beab7e49733bbf8cd777ed0785b1af56fbe8abd53ef2179a7bf73b6e6e7d9767ae6beebdbf6cd6a60a5e94002fc1ceda3e721adac969e079d79c081e7bff0cec8ab36cee1d67d93eecb1cd",
    "bdbf7d61d7b7efb64f89c7bb6a9c0b0be863d2dcee63c2f1d70a6aa137f7ee51d5dcaefb1bc8e33e5e18e35084d1bf25ac20774e1873e79adbb973b8dfdb2e202c40f76e1ff6e776eccf99",
    "bc821c3b61ac0fdc3e7ab8b5ddc32df27eaf3e0dede489f0395f7be2e6b6de1ec99b060df032623edb474fb8767ac2f13eae9cb8b66de685e7bdfc79c232f6d8ede48938f0965f083cee0d",
    "8431665b1863b6db6cb162f7d8c2b8c76e279f04cf65d64fc435dff1a7c26662087801fb6e61ec87d6dcbb1f9a30f6436beedd0fad9dbc0ef2b49fc61c3b6361cf37e8dff5288405e852cd",
    "bf75296101351684b1c6c2c096ffd2f2ae7cbf767af9f1faae18ffe68d6e015e466d75610dfaf8c02e1f9f30877c6761cc771e98749502deca59b63ce5ce45f0ae386e61193d6984cdc634",
    "e0d156dafcc75ae6d75a76f679cdbfe3a08445e448b4695e1c02bc8a78a9815dfd7b8431df59cb09f3f3daae6de3300ee65c8c7774f7446ffedd13bded3e853e9017514baa7953f316bc39",
    "6783575107b77953ebc0f2aebaf1c21c7547ffa13b7a933ff86bb18cbcc0815d7981c21a6abf0fecaafd3eb0ab074ff3ef180c6111b115c2185b218c75175a78e73e0863ee83b022fbb8e5",
    "b987fa6e30b50ec00bc8271696b0671fd89a372c6fcc1be67d84b7ad4f58a27efad17fa1edfe0b27966960574d276119b18d03d3f7cbfbcd0e35abdbe9a1c0e7571853d4c23ba6a89d1c1d",
    "deef95ebd1c2caf5e07bdbfad5afc5e80f10c6ba0103bbf2665bd8b6c84cdeb43bda73c8eec8e35df9fb2dacfcfd4cde95d3dac2f6dbf3de967d12bca26fc6f2caecad63cf3bfcb50d6363",
    "f641c5b39f3a038ef7d20fc2d20f8e6fb17df816dbf62d72acedb839f032fa5e0fecea7b2d8cfefd164cdf6bfc364a809790dbd282c96d99b65c6105fdad5a30f99c385e431de6667b1ed8",
    "6b5ef99ce0317675608a5dc5d898f992bcdfe164ad99bc2261cf5b2bfa632b9495fdb1850dff02afaf3535d331bc1ef23ac48fc522f27684cdbc1dfcb6482c5e7cf79f69f1dd7f4658a1ad",
    "257ed85ae287ad25be631b8515d43a68f15deb4018f31185d1f6d5e2ae25c0e35df58cdae9c519c9bb6c3c71d9221b7957dfcd76fa6ee2faee78ef16f7de349177cd9371cd937cceb37e09",
    "eec3ad26dee0cd0631e045c4dcb5f461334acb667462988535f8fc0776f9fc5b32ba3c78975d2fbd63610776f9b107269d9ad777afa1dbe76ff6c4e91d4f232c70aca58fb196de75805ada",
    "f52878cdb2c9646205fe73614d0ba9e52dfb8be5c9fec2e3d5196c6e79eda1be9bde7ebf76e2142a798d3a70da3a30dedbb45164f20afa64aac53cfb640af3a88529aca21ecfc0e2c3b19b",
    "f6d88de43147a77de4e8b49da39378bc558bd0f2d43406e365aefbe4adfa98e0310e4018e300da472f88964cbc1e78975e9cdebd930676db28767f88c6fbb8eac7b6d3e381e7bdfaae0b6b",
    "a883d13e7a41b493cfe4c90bf09bb6b4fca69c0f760c81e5b948bb457ec700b5fc3197e43d9704f07ca48f20bf7d682d7fecd9f3db87d6f2c7ba958ddf1e3cc6620fecaac3d7b259cbc00b",
    "9c73f2bb2e8eb0c47dc0478e533b394e0ebc7b6f9fb79d9ac7bbe27c5bde71be19bc7b6fbfe3872a7957bd5c61157dafd5baf691b0bc650c068ff5d15bde3d6e305e5eebefceb772186b33",
    "6e18ef77f6c23ed757de715aadbce3b406a67da803efb673959d979fc0bbfa1ab5b2fdab0f7857fcbd30d69113c6d86e6111b5955a79f73f6a65f73f0ae0c9c693c02b197164adbce3c884",
    "35fa12ca872f61e760359e574a56046f4c5a0f9ef3ddf34f187bfe0963cfbf815d3dffda8ea7f13cde55435c58a63da77cd8734e5f694f1e636cda8eb18978e7b30f2a7f7baf6f65c71df2",
    "5aae9a4ec21cf7b5c5d49b03c6bea5ad2cfb10c7c6b4fbf0b9948c3cb4f69187d6561e9aa9edd0cabbb68330d67650298f6bef5cf6ba856f55bd89f8adcefd349ed59d5fdb4eaedb39ef8e",
    "2f31768bba7cfec66e5157cd38b3a67ce479b5139741deb2e55a5eaaf42dd60fdfe28ead68e4dddf51dddf11782f1f5a5ddf91b167d7779dc581a9b6358f77c5c2b6bafc47bc8fad73fd5a",
    "2c51bfb2f95bf678a31612dfd19d7bd8eace3dc47dc887c1fb509c20efe365cba8cb9661ae6fc5a12c7fdeff07f1e303a5"
  ),collapse="")
  at<-seq.int(1L,nchar(compressed_hex),by=2L)
  raw<-as.raw(strtoi(substring(compressed_hex,at,at+1L),16L))
  t<-read.csv(text=rawToChar(memDecompress(raw,type="gzip")),stringsAsFactors=FALSE)
  v<-Map(seq.int,t$first,t$last);i<-rep(seq_len(nrow(t)),lengths(v))
  ref<-data.frame(SUBJID=paste0(t$SITEID[i],unlist(v,use.names=FALSE)),EUA=t$EUA_ADSL_RANDFL[i],DM=t$BLA_DM[i]==1L,stringsAsFactors=FALSE)
  stopifnot(!anyDuplicated(ref$SUBJID),sum(ref$EUA!="-")==30571,sum(ref$DM)==32102)
  ref
}

if(isTRUE(getOption("phmpt.autorun",TRUE)))analyse_subject_numbering()
