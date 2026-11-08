# Moderna P301, report v0.22. Select all and run from the project root in RStudio.
# Default: the named aggregate inputs distributed with this edition. For another run:
# options(phmpt.deviation_tables_dir=".../recomputed", phmpt.trigger_table_dir=".../protocol_trigger_audit")
render_audit_tables <- function(root=getOption("phmpt.root",getwd()),
  deviations=getOption("phmpt.deviation_tables_dir",file.path(root,"table_inputs_v022","deviations")),
  trigger=getOption("phmpt.trigger_table_dir",file.path(root,"table_inputs_v022","protocol")),
  out=getOption("phmpt.audit_asset_dir",file.path(root,"audit_assets"))) {
  selected<-getOption("phmpt.audit_table_ids",c("V06","V07","V17"));stopifnot(length(selected)>0,all(selected%in%c("V06","V07","V17")))
  root<-normalizePath(root,mustWork=TRUE)
  if(!grepl("^(/|[A-Za-z]:|\\\\)",out))out<-file.path(root,out)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  csv<-function(d,n)write.csv(d,file.path(out,n),row.names=FALSE,na="",fileEncoding="UTF-8")
  inputs<-c(file.path(deviations,c("inventory.csv","opportunities.csv","late_reporting_context.csv")),file.path(trigger,"arm_summary.csv"))
  load<-function(i,need){f<-inputs[i];if(!file.exists(f))stop("Missing table input: ",f)
    x<-read.csv(f,check.names=FALSE,fileEncoding="UTF-8-BOM");stopifnot(all(need%in%names(x)));x}
  inv<-load(1,c("code","arm","people"));opp<-load(2,c("positive_index","code","followup","days","arm","people_with_opportunity","affected_in_fixed_window"))
  late<-load(3,c("window","arm","category","Freq"));t<-load(4,c("scope","arm","n","form","swab","form_then_swab"))
  arms<-c("mRNA-1273","Placebo")
  pair<-function(x,col){stopifnot(nrow(x)==2L,setequal(x$arm,arms));as.numeric(x[[col]][match(arms,x$arm)])}
  counts<-function(x){stopifnot(!anyNA(x),all(is.finite(x)&x>=0&x==floor(x)));x}
  frac<-function(y,n){counts(c(y,n));stopifnot(all(y<=n));ifelse(n>0,sprintf("%s/%s (%.1f%%)",formatC(y,format="d",big.mark=","),formatC(n,format="d",big.mark=","),100*y/n),"0/0 (not estimable)")}
  codes<-c("P08.12","P07.08","P08.14");labels<-c("P08.12  Missed follow-up visit","P07.08  Missed saliva collection","P08.14  Missed illness monitoring")
  o<-subset(opp,positive_index=="first_initial_illness_positive_PCR"&followup=="complete_minimum_followup")
  v6<-do.call(rbind,lapply(seq_along(codes),function(i){
    a<-o[o$code==codes[i],];stopifnot(all(a$days==c(28,28,14)[i]));y<-counts(pair(a,"affected_in_fixed_window"));n<-counts(pair(a,"people_with_opportunity"))
    b<-counts(pair(inv[inv$code==codes[i],],"people"));stopifnot(all(y<=n),all(y<=b))
    data.frame(deviation=labels[i],days=c(28,28,14)[i],global_vaccine=b[1],global_placebo=b[2],vaccine_y=y[1],vaccine_n=n[1],placebo_y=y[2],placebo_n=n[2])}))
  catkeys<-c("COVID_AE_WITH_HOSPITAL_EVIDENCE","COVID_AE_WITHOUT_OWN_HOSPITAL_EVIDENCE","OTHER_SERIOUS_AE","SUSPECTED_OR_NONSPECIFIC_CORONAVIRUS_ONLY","NO_COVID_OR_SERIOUS_CANDIDATE")
  latelab<-c("Confirmed or explicitly COVID-coded AE, with hospitalization documented on that event","Confirmed or explicitly COVID-coded AE, without hospitalization documented on that event","Other serious AE","Suspected or nonspecific coronavirus AE only","No qualifying COVID or serious AE identified")
  z<-late[late$window==28,];stopifnot(setequal(z$category,catkeys))
  zz<-t(vapply(catkeys,function(k)counts(pair(z[z$category==k,],"Freq")),numeric(2)))
  total<-counts(pair(inv[inv$code=="I12.01",],"people"));stopifnot(all(colSums(zz)==total))
  v7<-data.frame(category=c(latelab,"Total with I12.01"),rbind(zz,total),row.names=NULL,check.names=FALSE);names(v7)[2:3]<-c("mRNA-1273","Placebo")
  q<-t[t$scope=="Required_investigation",];s<-t[t$scope=="All",]
  rows<-c("Checklist documented by day 3 among episodes warranting investigation","Swab collected by day 3 among all swab-indicated episodes","Swab collected on or after the first checklist date and by day 3, among swab-indicated episodes with a checklist by day 3")
  y<-rbind(pair(q,"form"),pair(s,"swab"),pair(s,"form_then_swab"));n<-rbind(pair(q,"n"),pair(s,"n"),pair(s,"form"))
  counts(c(y,n));stopifnot(all(y<=n),all(y[3,]<=y[2,]))
  v17<-data.frame(decision=rows,vaccine_y=y[,1],vaccine_n=n[,1],placebo_y=y[,2],placebo_n=n[,2])
  printtab<-function(id,slug,title,sub,head,body,widths,footer,data) {
    if(!id%in%selected)return(NULL)
    # Grid coordinates in inches; word wrapping measured in the active device font.
    dpi<-getOption("phmpt.audit_dpi",300L);w<-8.5;ht<-switch(id,V06=4.75,V07=5.55,V17=5.15)
    f<-paste0(id,"_",slug,".png"); grDevices::png(file.path(out,f),width=w,height=ht,units="in",res=dpi,bg="white")
    on.exit(grDevices::dev.off())
    grid::grid.newpage(); left<-.28; usable<-w-2*left
    gx<-function(x)grid::unit(x,"in"); gy<-function(y)grid::unit(ht-y,"in")
    ink<-"#17384b";gp<-function(size,bold=FALSE,col=ink)grid::gpar(fontsize=size,fontface=if(bold)"bold" else "plain",col=col,fontfamily="sans")
    txt<-function(x,y,s,size,bold=FALSE,col=ink){grid::grid.text(s,gx(x),gy(y),just=c("left","top"),gp=gp(size,bold,col))}
    wrap<-function(s,width,size,bold=FALSE){
      measure<-function(v)grid::convertWidth(grid::grobWidth(grid::textGrob(v,gp=gp(size,bold))),"in",valueOnly=TRUE)
      lines<-character();line<-""
      for(word in strsplit(s," ",fixed=TRUE)[[1]]){v<-trimws(paste(line,word));if(nzchar(line)&&measure(v)>width){lines<-c(lines,line);line<-word}else line<-v}
      c(lines,line)}
    txt(left,.24,paste(id,title),17,TRUE);txt(left,.64,sub,10,FALSE,"#385565")
    x<-left+c(0,cumsum(widths*usable)); y0<-1.05
    row<-function(vals,head=FALSE,shade=FALSE,bold=FALSE) {
      size<-if(head)10.5 else 11; l<-Map(function(j)wrap(vals[j],widths[j]*usable-.18,size,head||bold),seq_along(vals))
      hh<-max(lengths(l))*.205+.20
      if(head||shade)grid::grid.rect(gx(left),gy(y0),gx(usable),gx(hh),just=c("left","top"),gp=grid::gpar(fill=if(head)"#dce8ed" else "#f1f5f6",col=NA))
      for(j in seq_along(l))for(k in seq_along(l[[j]]))txt(x[j]+.08,y0+.10+(k-1)*.205,l[[j]][k],size,head||bold)
      y0<<-y0+hh }
    row(head,TRUE)
    for(i in seq_len(nrow(body)))row(as.character(unlist(body[i,],use.names=FALSE)),shade=i%%2==0,bold=id=="V07"&&i==nrow(body))
    lines<-wrap(footer,usable-.12,10);stopifnot(y0+.23+length(lines)*.18<ht)
    for(i in seq_along(lines))txt(left+.04,y0+.20+(i-1)*.18,lines[i],10,FALSE,"#395766")
    csv(data,paste0(id,"_table_data.csv"))
    data.frame(visual=id,file=f,status="AVAILABLE",page=1,pages=1,stringsAsFactors=FALSE)
  }
  m6<-printtab("V06","deviations_after_confirmed_illness","Deviations after confirmed COVID-19","Production | Individual blinded follow-up | Participants with a qualifying initial illness-visit PCR",
    c("Recorded deviation","V05 total: vaccine / placebo","mRNA-1273: affected / eligible","Placebo: affected / eligible"),
    data.frame(v6$deviation,paste(v6$global_vaccine,v6$global_placebo,sep=" / "),frac(v6$vaccine_y,v6$vaccine_n),frac(v6$placebo_y,v6$placebo_n)),
    c(.33,.18,.245,.245),"The V05 counts use all randomized participants. The last two columns require observation through day 28 for visits and saliva collection, or through day 14 for monitoring; day 0 is the index PCR collection date.",v6)
  m7<-printtab("V07","late_reporting_event_context","Events recorded near late safety notifications","Production | Participants with I12.01 | AE onset within 28 days before or after a deviation",
    c("Highest-priority event found","mRNA-1273","Placebo"),v7,c(.66,.17,.17),
    "Each participant appears once, in the first applicable category. Hospitalization must be recorded on the COVID event itself. Proximity in time does not identify which event triggered the notice.",v7)
  m17<-printtab("V17","investigation_and_swabbing","Investigation and swabbing by day 3","Production | Blinded illness episodes | Three further calendar days of observation required",
    c("Decision and eligible episodes","mRNA-1273","Placebo"),data.frame(rows,frac(y[,1],n[,1]),frac(y[,2],n[,2])),c(.49,.255,.255),
    "Day 0 is the relevant protocol trigger. A same-day checklist and swab do not establish their order. Day 3 is a calendar-day measure, not 72 verified hours from first contact.",v17)
  for(m in Filter(Negate(is.null),list(m6,m7,m17))){m$md5<-unname(tools::md5sum(file.path(out,m$file)));csv(m,paste0(m$visual,"_audit_table_manifest.csv"))}
  csv(data.frame(input=basename(inputs),path=normalizePath(inputs,winslash="/"),md5=unname(tools::md5sum(inputs))),"audit_table_inputs.csv")
  message("Wrote ",paste(selected,collapse=", ")," in ",normalizePath(out,winslash="/"))
}
if(isTRUE(getOption("phmpt.autorun",TRUE)))render_audit_tables()
