# Moderna P301, opening figures V01 and V02. Select all and Run from the project
# root in RStudio. Base R only. Both PNGs go into audit_assets/ at 300 dpi.
# Reprint the released audit from the three bundled CSVs, or explicitly set:
# options(phmpt.intro_lab_file=..., phmpt.intro_coverage_file=..., phmpt.intro_trigger_file=...)
# To generate only one: options(phmpt.intro_figures="V01") or "V02".
render_intro_figures <- function(root=getOption("phmpt.root",getwd()),
  selected=getOption("phmpt.intro_figures",c("V01","V02")),
  out=getOption("phmpt.audit_asset_dir","audit_assets"),
  inputs=getOption("phmpt.intro_input_dir","table_inputs_intro"),
  dpi=getOption("phmpt.audit_dpi",300L)) {
  root<-normalizePath(root,mustWork=TRUE);path<-function(x)if(grepl("^(/|[A-Za-z]:|\\\\)",x))x else file.path(root,x)
  out<-path(out);inputs<-path(inputs);dir.create(out,showWarnings=FALSE,recursive=TRUE)
  stopifnot(length(selected)>0L,all(selected%in%c("V01","V02")),length(dpi)==1,is.finite(dpi),dpi>0)
  ink<-"#17384B"; muted<-"#466474"; teal<-"#137268"; blue<-"#235979"
  U<-function(x)grid::unit(x,"inches")
  txt<-function(x,y,label,size=12,bold=FALSE,col=ink,align="left"){
    grid::grid.text(label,U(x),U(y),just=c(align,"centre"),
      gp=grid::gpar(fontfamily="sans",fontsize=size,fontface=if(bold)"bold" else "plain",col=col))
  }
  rect<-function(x,y,w,h,fill,col=NA,lwd=1,lty=1){
    grid::grid.rect(U(x),U(y),U(w),U(h),just=c("left","bottom"),
      gp=grid::gpar(fill=fill,col=col,lwd=lwd,lty=lty))}
  seg<-function(x,y,xx,yy,col=muted,lty=1,arrow=FALSE){
    grid::grid.lines(U(c(x,xx)),U(c(y,yy)),gp=grid::gpar(col=col,lwd=1.35,lty=lty),
       arrow=if(arrow)grid::arrow(length=U(.068),type="closed") else NULL)}
  csv<-function(x,f)utils::write.csv(x,file.path(out,f),row.names=FALSE,na="",fileEncoding="UTF-8")
  draw<-function(id,slug,h,fun){
    name<-paste0(id,"_",slug,".png");dst<-file.path(out,name);tmp<-tempfile(tmpdir=out,fileext=".png")
    on.exit(unlink(tmp),add=TRUE)
    grDevices::png(tmp,width=8.5,height=h,units="in",res=dpi,bg="white")
    tryCatch({grid::grid.newpage();fun()},finally=grDevices::dev.off())
    if(!file.exists(tmp)||file.info(tmp)$size<2000||!file.copy(tmp,dst,overwrite=TRUE))stop("PNG write failed: ",dst)
    csv(data.frame(visual=id,file=name,status="AVAILABLE",md5=unname(tools::md5sum(dst))),
      paste0(id,"_intro_manifest.csv"))
  }
  if("V01"%in%selected)draw("V01","symptom_to_efficacy_case",4.60,function(){
    txt(.25,4.25,"V01  Only a swabbed, PCR-positive illness can become a case",17,TRUE)
    txt(.25,3.89,"How a complaint reported to a trial site could reach the primary efficacy endpoint",11,FALSE,muted)
    node<-function(x,y,w,h,title,body,fill="#f1f5f6",border="#92a7af",dash=1){
      rect(x,y,w,h,fill,border,1.25,dash);cx<-x+w/2;cy<-y+h/2
      txt(cx,cy+.18,title,11.3,TRUE,ink,"centre")
      for(j in seq_along(body))txt(cx,cy-.08-(j-1)*.19,body[j],10.25,FALSE,muted,"centre")}
    X0<-.25; X1<-2.22; X2<-4.32; X3<-6.44; leftend<-1.78;branch<-1.98
    yu<-2.83;yl<-1.31;ym<-(yu+yl)/2; H<-.96
    seg(leftend,ym,branch,ym);seg(branch,yl,branch,yu)
    seg(branch,yu,X1-.025,yu,arrow=TRUE);seg(branch,yl,X1-.025,yl,arrow=TRUE)
    seg(X1+1.77,yu,X2-.025,yu,arrow=TRUE);seg(X2+1.75,yu,X3-.025,yu,arrow=TRUE)
    seg(X1+1.77,yl,X2-.025,yl,lty=2,arrow=TRUE)
    node(X0,ym-H/2,leftend-X0,H,"Symptom reported",c("to the trial site"))
    node(X1,yu-H/2,1.77,H,"COVID checklist",c("illness assessed"))
    node(X2,yu-H/2,1.75,H,"Nasal swab + PCR",c("sampling target:","within 72 hours"))
    node(X3,yu-H/2,1.81,H,"Efficacy case",c("qualifying illness +","positive PCR"),"#e1eef3",blue)
    node(X1,yl-H/2,1.77,H,"Adverse event",c("AE record only"))
    node(X2,yl-H/2,1.75,H,"No swab on record",c("no PCR-confirmed","case documented"),"#fafafa","#8ba1ab",2)
    txt(.25,.52,"A complaint could appear in both records. An AE entry alone does not establish whether a swab was taken.",10,FALSE,muted)
    txt(.25,.25,"For the primary endpoint, qualifying illness had to start at least 14 days after dose 2 in an eligible participant.",10,FALSE,muted)
  })
  if("V02"%in%selected){
    files<-c(lab=getOption("phmpt.intro_lab_file",file.path(inputs,"symptom_summary_global.csv")),
      headline=getOption("phmpt.intro_coverage_file",file.path(inputs,"headline_pcr_coverage.csv")),
      protocol=getOption("phmpt.intro_trigger_file",file.path(inputs,"arm_summary.csv")))
    files[]<-vapply(files,path,character(1))
    read<-function(k,cols){
      if(!file.exists(files[[k]]))stop("Missing ",k," aggregate: ",files[[k]])
      d<-utils::read.csv(files[[k]],check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
      if(length(e<-setdiff(cols,names(d))))stop(k,": missing ",paste(e,collapse=", "))
      d
    }
    arms<-c("mRNA-1273","Placebo")
    lab<-read("lab",c("release","period","scope","stratum","granularity","endpoint","arm","n_eligible","n_tested"))
    lab<-lab[lab$release=="Prod"&lab$period=="BLA_SNAPSHOT"&lab$scope=="STRICT_BLINDED"&lab$granularity=="GLOBAL"&lab$endpoint=="Any PCR",]
    ordered<-function(x){if(nrow(x)!=2L||anyDuplicated(x$arm)||!setequal(x$arm,arms))stop("Expected one record per arm");x[match(arms,x$arm),]}
    headline<-read("headline",c("scope","vaccine_total","vaccine_tested","vaccine_pct","placebo_total","placebo_tested","placebo_pct","gap_pp","low95","high95","vaccine_low95","vaccine_high95","placebo_low95","placebo_high95","bootstrap_draws","n_sites_in_resampling_frame","seed","interval_method"))
    scopes<-c("Participants: blinded | no linked positive",
      "Participants: blinded | no linked positive | excluding post-dose days 0-6",
      "Participants: PP, dose 2 + 14 days | no linked positive")
    if(nrow(headline)!=3||!setequal(headline$scope,scopes))stop("V02 requires all three unique participant scopes. Use the supplied headline_pcr_coverage.csv or rerun the updated analyse_symptom_critique.R.")
    headline<-headline[match(scopes,headline$scope),]
    if(any(!is.finite(as.matrix(headline[c("bootstrap_draws","n_sites_in_resampling_frame","seed")]))) ||
       length(unique(headline$bootstrap_draws))!=1L||length(unique(headline$n_sites_in_resampling_frame))!=1L||
       any(headline$bootstrap_draws<=0 | headline$n_sites_in_resampling_frame<2)||
       !all(grepl("site-cluster",headline$interval_method,fixed=TRUE)))stop("Inconsistent bootstrap metadata.")
    t<-read("protocol",c("scope","arm","form_then_swab","form"))
    t<-ordered(t[t$scope=="All",])
    rows<-list()
    row<-function(label,y,n,unit="Participants",r=NULL){
      if(any(!is.finite(c(y,n)))||any(c(y,n)!=floor(c(y,n)))||any(y<0|y>n|n<=0))stop("Invalid counts: ",label)
      p<-100*y/n;lo<-hi<-c(NA_real_,NA_real_);gap<-gl<-gh<-NA_real_
      if(!is.null(r)){
        lo<-c(r$vaccine_low95,r$placebo_low95);hi<-c(r$vaccine_high95,r$placebo_high95)
        gap<-r$gap_pp;gl<-r$low95;gh<-r$high95
        if(!all(is.finite(c(lo,hi,gap,gl,gh)))||any(lo>hi|lo<0|hi>100)||gl>gh||
           any(abs(p-c(r$vaccine_pct,r$placebo_pct))>1e-7)||abs((p[1]-p[2])-gap)>1e-7)stop("Counts/rates/intervals do not agree: ",label)
      }
      data.frame(label,unit,vaccine_tested=y[1],vaccine_total=n[1],vaccine_pct=p[1],vaccine_low95=lo[1],vaccine_high95=hi[1],
        placebo_tested=y[2],placebo_total=n[2],placebo_pct=p[2],placebo_low95=lo[2],placebo_high95=hi[2],gap_pp=gap,gap_low95=gl,gap_high95=gh)
    }
    for(i in 1:2){
      pop<-c("COVID_FORM_PLUS_AE","COVID_FORM_ONLY")[i]
      s<-ordered(lab[lab$stratum==paste0("SYMPTOMS__",pop,"__BY_CUTOFF__DIAGNOSTIC"),])
      rows[[i]]<-row(c("Symptoms on either record","Symptoms on a COVID checklist")[i],s$n_tested,s$n_eligible)
    }
    for(i in 1:3){r<-headline[i,]; rows[[i+2]]<-row(c("All blinded follow-up; no linked positive PCR",
        "Same full-window no-positive cohort; symptom dates outside days 0-6 after either dose",
        "Per-protocol window; no linked positive PCR")[i],
        c(r$vaccine_tested,r$placebo_tested),c(r$vaccine_total,r$placebo_total),r=r)}
    rows[[6]]<-row("Once a checklist existed: swab by day 3",t$form_then_swab,t$form,"Illness episodes")
    z<-do.call(rbind,rows);csv(z,"V02_table_data.csv")
    csv(data.frame(input=names(files),file=basename(files),md5=unname(tools::md5sum(files))),"V02_inputs.csv")
    fmt<-function(x)formatC(as.integer(x),format="d",big.mark=",");signpc<-function(x)sprintf("%+.1f",x)
    draw("V02","testing_gap_per_protocol_window",7.20,function(){
      txt(.25,6.88,"V02  The testing gap narrowed in the per-protocol window",17,TRUE)
      txt(.25,6.51,"Diagnostic PCR coverage in the blinded licensing data",12,FALSE,muted)
      rect(.25,5.97,8.0,.36,"#dce8ed")
      txt(.37,6.15,"Comparison",11,TRUE);txt(4.25,6.15,"mRNA-1273",11,TRUE);txt(6.28,6.15,"Placebo",11,TRUE)
      labels<-list("Symptoms on either record","Symptoms on a COVID checklist",
        c("All blinded follow-up;","no linked positive PCR"),
        c("Same cohort; set aside symptoms","on days 0-6 after either dose"),
        c("Per-protocol, baseline-negative;","from 14 days after dose 2;","no linked positive PCR"),
        c("Swab by day 3 among swab-indicated","episodes with a checklist by day 3"))
      top<-5.97
      for(i in seq_len(nrow(z))){
        h<-c(.43,.43,.68,.68,.84,.67)[i];mid<-top-h/2
        if(i%in%3:5)rect(.25,top-h,8,h,if(i!=4)"#eef6f5" else "#e3f1ee")
        seg(.25,top-h,8.25,top-h,"#c8d4d7")
        laby<-mid+.104*(length(labels[[i]])-1L)
        for(j in seq_along(labels[[i]]))txt(.37,laby-.208*(j-1),labels[[i]][j],11,i%in%3:5)
        for(a in c("vaccine","placebo")){
          x<-if(a=="vaccine")4.25 else 6.28
          v<-function(s)z[[paste0(a,s)]][i];ci<-is.finite(v("_low95"));y<-if(ci)mid+.11 else mid
          txt(x,y,sprintf("%s/%s (%.1f%%)",fmt(v("_tested")),fmt(v("_total")),v("_pct")),11.4,TRUE)
          if(ci)txt(x,mid-.135,sprintf("95%% CI %.1f%% to %.1f%%",v("_low95"),v("_high95")),10,FALSE,muted)
        }
        top<-top-h
      }
      txt(.30,2.00,"Vaccine minus placebo, percentage points (95% CI)",11.6,TRUE)
      dif<-function(i)sprintf("%s (%s to %s)",signpc(z$gap_pp[i]),signpc(z$gap_low95[i]),signpc(z$gap_high95[i]))
      txt(.30,1.72,paste0("All blinded follow-up: ",dif(3),"     |     Days 0-6 set aside: ",dif(4)),11.1,TRUE,teal)
      txt(.30,1.45,paste0("Per-protocol window: ",dif(5)),11.1,TRUE,teal)
      txt(.30,1.04,"Rows 1-5 count participants; row 6 counts illnesses, from the swab trigger. Untested participants remain",10,FALSE,muted)
      txt(.30,.83,"eligible. Row 4 removes symptom dates, not people with later symptoms; the no-positive rule applies",10,FALSE,muted)
      txt(.30,.62,"across the entire blinded period. Row 5 assesses linked positives within the per-protocol window.",10,FALSE,muted)
      txt(.30,.36,paste0("Intervals resample whole sites, keeping both arms together (",fmt(headline$bootstrap_draws[1])," draws; ",headline$n_sites_in_resampling_frame[1]," sites)."),10,FALSE,muted)
    })
  }
  message("Wrote ",paste(selected,collapse=", ")," in ",normalizePath(out,winslash="/"));invisible(out)
}
if(isTRUE(getOption("phmpt.autorun",TRUE)))render_intro_figures()
