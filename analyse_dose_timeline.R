# Moderna P301 - weekly dose 1 / dose 2 administrations - v1.1.0
# Stand-alone, R >= 4.3: save in your project root, then Select All / Run in RStudio.
# Reads local ADSL: CSV, XPT, SAS7BDAT or the released P301 JMP format.
# CSV/JMP and all exports need only base R. XPT/SAS7BDAT require haven.
# Original blinded study injections (mRNA-1273 OR placebo), not crossover doses.
# ADSL derives TR01SDT from DM.RFXSTDTC and DOSE2DT from dated Day 29 EX records.
# DOSE1FL / DOSE2FL independently identify recorded administration.

# ---- SETTINGS ----
PHMPT_ROOT <- ""                    # Auto-detect; or D:/Projects/moderna_trial_analysis
PHMPT_RELEASE <- "auto"             # "auto" prefers Prod, otherwise EUA; never combines them
PHMPT_ADSL_FILE <- ""               # Optional explicit CSV, XPT, SAS7BDAT or JMP path
PHMPT_POPULATION <- "randomized"    # "randomized" = RANDFL Y; "all_adsl" includes other dosed rows
PHMPT_OUTPUT_DIR <- ""              # Default: analysis/dose_timeline_<timestamp>
PHMPT_CSV_ENCODING <- "auto"        # UTF-8-BOM, then Windows-1252 if decoding fails
PHMPT_NUMERIC_DATE_ORIGIN <- ""     # For unformatted CSV dates only; never guessed
PHMPT_MIN_RANDOMIZED <- 30000L      # Reject subject subsets and incomplete extracts
PHMPT_PNG_DPI <- 300L
PHMPT_AUDIT_ASSET_DIR <- "audit_assets" # Relative to project root; stable V04 PNG
PHMPT_AUDIT_DPI <- 300L                # 8.5-inch width = 2550 px at 300 dpi
PHMPT_SHOW_VALUES <- FALSE          # Optional count labels above individual columns
PHMPT_SHOW_CHART <- interactive()
# Load functions without running: options(phmpt.autorun=FALSE) before source().

dt_text <- function(x) {
  if(inherits(x,"Date")) x<-format(x,"%Y-%m-%d")
  else if(inherits(x,"POSIXt")) x<-format(x,"%Y-%m-%d",tz="UTC")
  x<-as.character(x);x[is.na(x)]<-""
  bad<-!validUTF8(x)
  if(any(bad))x[bad]<-iconv(x[bad],from="windows-1252",to="UTF-8",sub=NA)
  if(anyNA(x))stop("Undecodable input text. Set PHMPT_CSV_ENCODING.")
  enc2utf8(x)
}
dt_csv <- function(x,path) {
  con<-file(path,"wb");on.exit(close(con));writeBin(as.raw(c(239,187,191)),con)
  write.table(x,con,sep=",",row.names=FALSE,quote=TRUE,na="",qmethod="double")
}
dt_path <- function(root,p) if(grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)",p))p else file.path(root,p)
dt_root <- function(root="") {
  if(nzchar(root))return(normalizePath(root,winslash="/",mustWork=TRUE))
  cand<-getwd()
  for(f in sys.frames())if(!is.null(f$ofile))cand<-c(cand,dirname(f$ofile))
  a<-grep("^--file=",commandArgs(),value=TRUE)
  if(length(a))cand<-c(cand,dirname(sub("^--file=","",a[1])))
  if(interactive()&&requireNamespace("rstudioapi",quietly=TRUE)&&rstudioapi::isAvailable()) {
    p<-tryCatch(rstudioapi::getActiveDocumentContext()$path,error=function(e)"")
    if(nzchar(p))cand<-c(cand,dirname(p))
  }
  for(p in unique(cand[dir.exists(cand)]))for(i in seq_len(8L)) {
    if(any(dir.exists(file.path(p,c("prod_data","eua_data")))))return(normalizePath(p,winslash="/",mustWork=TRUE))
    parent<-dirname(p);if(parent==p)break;p<-parent
  }
  stop("Cannot locate eua_data/prod_data. Set PHMPT_ROOT to the project root.")
}
dt_fields <- c("STUDYID","USUBJID","RANDFL","TRT01P","DOSE1FL","DOSE2FL","TR01SDT","DOSE2DT","CUTOFFDT")
dt_date <- function(x,field,origin="") {
  raw<-trimws(dt_text(x));raw[raw %in% c("NA",".")]<-""
  out<-as.Date(rep(NA_character_,length(raw)))
  iso<-grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}($|T| )",raw)
  out[iso]<-suppressWarnings(as.Date(substr(raw[iso],1,10),"%Y-%m-%d"))
  wrong<-iso & (is.na(out)|format(out,"%Y-%m-%d")!=substr(raw,1,10));wrong[is.na(wrong)]<-FALSE
  out[wrong]<-as.Date(NA)
  # SAS-style English dates, without relying on the Windows/French locale.
  sas<-grepl("^[0-9]{1,2}[A-Za-z]{3}[0-9]{4}$",raw)
  if(any(sas)) {
    s<-toupper(raw[sas]);d<-as.integer(sub("[A-Z].*$","",s))
    m<-match(sub("^[0-9]+([A-Z]{3}).*$","\\1",s),toupper(month.abb))
    y<-as.integer(substr(s,nchar(s)-3,nchar(s)))
    v<-sprintf("%04d-%02d-%02d",y,m,d);dd<-suppressWarnings(as.Date(v,"%Y-%m-%d"))
    dd[!is.na(dd)&format(dd,"%Y-%m-%d")!=v]<-as.Date(NA);out[sas]<-dd
  }
  num<-nzchar(raw)&grepl("^-?[0-9]+([.][0-9]+)?$",raw)
  if(any(num)) {
    if(!nzchar(origin))stop(field,": unformatted numeric dates. Set PHMPT_NUMERIC_DATE_ORIGIN from documented metadata, or use the original XPT/JMP.")
    z<-as.numeric(raw[num]);if(any(z!=floor(z)))stop(field,": numeric date has a fractional day.")
    out[num]<-as.Date(z,origin=origin)
  }
  if(any(nzchar(raw)&is.na(out)))stop(field,": partial, invalid or unsupported date(s): ",paste(head(unique(raw[nzchar(raw)&is.na(out)]),5),collapse=", "))
  out
}
dt_read <- function(path,encoding="auto",origin="",min_randomized=30000L) {
  ext<-tolower(tools::file_ext(path))
  if(ext=="jmp")x<-dt_read_jmp(path,dt_fields)
  else if(ext=="csv") {
    x<-NULL;errors<-character()
    for(enc in if(encoding=="auto")c("UTF-8-BOM","windows-1252") else encoding) {
      y<-tryCatch(withCallingHandlers(read.csv(path,colClasses="character",na.strings=c("","NA"),check.names=FALSE,fileEncoding=enc),warning=function(w)stop(conditionMessage(w))),error=function(e)e)
      if(!inherits(y,"error")){x<-y;break};errors<-c(errors,conditionMessage(y))
    }
    if(is.null(x))stop("CSV read failed: ",paste(errors,collapse="; "))
  } else if(ext %in% c("xpt","sas7bdat")) {
    if(!requireNamespace("haven",quietly=TRUE))stop('Install haven for XPT/SAS7BDAT: install.packages("haven")')
    x<-as.data.frame(if(ext=="xpt")haven::read_xpt(path,.name_repair="minimal") else haven::read_sas(path,.name_repair="minimal"))
  } else stop("Unsupported input format: ",ext)
  names(x)<-toupper(sub("^\ufeff","",names(x)))
  if(anyDuplicated(names(x)))stop("Duplicate source column names.")
  missing<-setdiff(dt_fields,names(x));if(length(missing))stop("Missing ADSL fields: ",paste(missing,collapse=", "))
  x<-x[dt_fields];x[]<-lapply(x,dt_text)
  x<-x[toupper(trimws(x$STUDYID))=="MRNA-1273-P301",,drop=FALSE]
  if(!nrow(x))stop("No P301 participants.")
  if(any(!nzchar(trimws(x$USUBJID)))||anyDuplicated(x$USUBJID))stop("Missing or duplicate USUBJID; reconcile the ADSL first.")
  for(n in c("RANDFL","DOSE1FL","DOSE2FL")) {
    x[[n]]<-toupper(trimws(x[[n]]))
    if(any(!x[[n]] %in% c("Y","N")))stop(n,": unexpected or missing flag.")
  }
  if(sum(x$RANDFL=="Y")<min_randomized)stop("ADSL is incomplete or a subject subset: ",sum(x$RANDFL=="Y")," randomized rows.")
  for(n in c("TR01SDT","DOSE2DT","CUTOFFDT"))x[[n]]<-dt_date(x[[n]],n,origin)
  if(anyNA(x$CUTOFFDT)||length(unique(x$CUTOFFDT))!=1L)stop("ADSL must have one non-missing CUTOFFDT.")
  for(pair in list(c("DOSE1FL","TR01SDT"),c("DOSE2FL","DOSE2DT"))) {
    bad<-(x[[pair[1]]]=="Y")!=!is.na(x[[pair[2]]])
    if(any(bad))stop(pair[1]," disagrees with ",pair[2]," for ",sum(bad)," participant(s), including ",paste(head(x$USUBJID[bad],3),collapse=", "))
  }
  if(any(x$DOSE2FL=="Y" & x$DOSE1FL!="Y"))stop("Dose 2 recorded without dose 1.")
  if(any(x$DOSE2DT<x$TR01SDT,na.rm=TRUE))stop("Dose 2 date precedes dose 1.")
  if(any(x$TR01SDT>x$CUTOFFDT|x$DOSE2DT>x$CUTOFFDT,na.rm=TRUE))stop("Dosing dates after CUTOFFDT: investigate the source.")
  x$TRT01P<-trimws(x$TRT01P)
  x$TRT01P[toupper(x$TRT01P)=="PLACEBO"]<-"Placebo"
  x$TRT01P[toupper(x$TRT01P)=="MRNA-1273"]<-"mRNA-1273"
  if(any(x$RANDFL=="Y" & !x$TRT01P %in% c("Placebo","mRNA-1273")))stop("Unexpected randomized treatment arm.")
  x
}
dt_candidates <- function(root,release) {
  paths<-list.files(file.path(root,if(release=="Prod")"prod_data" else "eua_data"),recursive=TRUE,full.names=TRUE,pattern="\\.(csv|xpt|sas7bdat|jmp)$",ignore.case=TRUE)
  paths<-paths[grepl("(^|[_ -])adsl([_. -]|$)",basename(paths),ignore.case=TRUE)]
  paths<-paths[!grepl("p101|p201|dmid",basename(paths),ignore.case=TRUE)]
  paths[order(match(tolower(tools::file_ext(paths)),c("csv","jmp","xpt","sas7bdat")),paths)]
}
dt_signature <- function(x) {
  y<-x[order(x$USUBJID),];rownames(y)<-NULL;p<-tempfile();on.exit(unlink(p))
  saveRDS(y,p,compress=FALSE,version=2);unname(tools::md5sum(p))
}
dt_select <- function(root,release,file,out,encoding,origin,min_randomized) {
  if(!release %in% c("auto","Prod","EUA"))stop("PHMPT_RELEASE must be auto, Prod or EUA.")
  plan<-if(nzchar(file))data.frame(release=if(release=="auto")"Explicit" else release,path=dt_path(root,file)) else {
    rr<-if(release=="auto")c("Prod","EUA") else release
    do.call(rbind,lapply(rr,function(r){p<-dt_candidates(root,r);data.frame(release=rep(r,length(p)),path=p)}))
  }
  if(!nrow(plan))stop("No local ADSL found. Set PHMPT_ADSL_FILE to your complete P301 ADSL CSV, XPT, SAS7BDAT or JMP.")
  audit<-list();valid<-list();hashes<-character()
  for(i in seq_len(nrow(plan))) {
    p<-plan$path[i];rel<-plan$release[i];message("Checking ",rel," ADSL: ",basename(p))
    h<-if(file.exists(p))unname(tools::md5sum(p)) else ""
    key<-paste(rel,h,sep=":");duplicate<-nzchar(h)&&key %in% hashes
    a<-if(duplicate)NULL else tryCatch(dt_read(p,encoding,origin,min_randomized),error=function(e)e)
    ok<-!duplicate&&!inherits(a,"error")
    status<-if(duplicate)"IDENTICAL_COPY" else if(ok)"VALID" else "REJECTED"
    audit[[i]]<-data.frame(release=rel,path=normalizePath(p,winslash="/",mustWork=FALSE),md5=h,
      n_subjects=if(ok)nrow(a) else NA_integer_,n_randomized=if(ok)sum(a$RANDFL=="Y") else NA_integer_,
      cutoff=if(ok)format(a$CUTOFFDT[1]) else "",status=status,detail=if(duplicate)"Identical bytes in this release." else if(ok)"" else conditionMessage(a))
    if(ok)valid[[length(valid)+1L]]<-list(index=i,data=a,signature=dt_signature(a),release=rel)
    hashes<-c(hashes,key)
  }
  audit<-do.call(rbind,audit);log<-file.path(out,"adsl_discovery.csv");dt_csv(audit,log)
  if(!length(valid))stop("No complete usable ADSL. See ",log,". ",paste(unique(audit$detail),collapse=" | "))
  if(!nzchar(file)&&release=="auto") {
    preferred<-if(any(vapply(valid,function(z)z$release=="Prod",logical(1))))"Prod" else "EUA"
    if(preferred=="EUA")warning("No complete usable Prod ADSL found; using EUA. See adsl_discovery.csv.")
    valid<-Filter(function(z)z$release==preferred,valid)
  }
  latest<-max(vapply(valid,function(z)as.numeric(z$data$CUTOFFDT[1]),numeric(1)))
  valid<-Filter(function(z)as.numeric(z$data$CUTOFFDT[1])==latest,valid)
  if(length(unique(vapply(valid,function(z)z$signature,character(1))))!=1L)
    stop("Conflicting ADSL snapshots share the latest cutoff. Set PHMPT_ADSL_FILE explicitly; see ",log)
  selected<-valid[[1]];audit$status[audit$status=="VALID"]<-"NOT_SELECTED"
  audit$status[selected$index]<-"SELECTED";dt_csv(audit,log)
  selected$manifest<-audit[selected$index,,drop=FALSE];selected
}
dt_week <- function(date) date-(as.integer(format(date,"%w"))+6L)%%7L
dt_label <- function(date) paste(format(date,"%d"),month.abb[as.integer(format(date,"%m"))])
dt_n <- function(x) format(x,big.mark=",",scientific=FALSE,trim=TRUE)
dt_chart <- function(weekly,totals,pop_n,population,cutoff,release,show_values=FALSE) {
  old<-par(no.readonly=TRUE);on.exit(par(old))
  par(mar=c(5.5,5.0,4.6,1.0),mgp=c(2.8,.7,0),las=1,family="sans",xaxs="i",yaxs="i",bty="n",cex=.92)
  vals<-rbind(weekly$dose_1,weekly$dose_2);cols<-c("#1E667D","#D27C34")
  ymax<-max(vals);top<-max(pretty(c(0,max(1,ymax)*1.12),n=5));ticks<-pretty(c(0,top),n=5);ticks<-ticks[ticks>=0&ticks<=top]
  positions<-barplot(vals,beside=TRUE,space=c(.12,.65),col=cols,border=NA,axes=FALSE,axisnames=FALSE,ylim=c(0,top),ylab="Participants receiving the dose")
  abline(h=ticks,col="#DFE4E7",lwd=.6)
  # Repaint bars over the guide lines, without opening a second plot.
  barplot(vals,beside=TRUE,space=c(.12,.65),col=cols,border=NA,axes=FALSE,axisnames=FALSE,ylim=c(0,top),add=TRUE)
  axis(2,at=ticks,labels=dt_n(ticks),tick=FALSE,col.axis="#25343D",cex.axis=.85)
  centres<-colMeans(positions);axis(1,at=centres,labels=FALSE,tck=-.012,col="#7E8B93")
  text(centres,par("usr")[3]-.022*top,labels=dt_label(as.Date(weekly$week_start)),srt=55,adj=1,xpd=NA,cex=.76,col="#25343D")
  title(main="Moderna P301: weekly dose 1 and dose 2 administrations",adj=0,line=2.65,cex.main=1.16,col.main="#172F3C")
  endpoints<-c(totals$first_date,totals$last_date)
  years<-range(as.integer(format(as.Date(endpoints[nzchar(endpoints)]),"%Y")))
  period<-if(years[1]==years[2])as.character(years[1]) else paste(years,collapse="-")
  title(sub=paste0("Weeks start on Monday (",period,") | ",release," cutoff ",cutoff," | Original blinded injections: mRNA-1273 and placebo"),line=4.15,cex.sub=.78,col.sub="#435863")
  mtext(paste0(if(population=="randomized")"Randomized population" else "All P301 ADSL participants",": N = ",dt_n(pop_n)," | Zero-dose participants contribute no columns"),side=3,line=1.15,adj=0,cex=.83,col="#435863")
  legend("topright",legend=paste0(c("Dose 1: ","Dose 2: "),dt_n(totals$n_administered)),fill=cols,border=NA,bty="n",cex=.88)
  if(show_values) {
    good<-as.vector(vals)>0
    text(as.vector(positions)[good],as.vector(vals)[good]+.015*top,labels=dt_n(as.vector(vals)[good]),cex=.57,xpd=NA,srt=90,adj=0)
  }
}

# ---- AUDIT PNG EXPORTS: dedicated article-width layouts ----
dt_audit_png <- function(path,draw,width=8.5,height=7,dpi=300L) {
  if(length(dpi)!=1L||!is.finite(dpi)||dpi<150||dpi>600)
    stop("PHMPT_AUDIT_DPI must be between 150 and 600.")
  type<-if(capabilities("cairo"))"cairo" else if(.Platform$OS.type=="windows")"windows" else getOption("bitmapType")
  success<-FALSE
  grDevices::png(path,width=width,height=height,units="in",res=dpi,pointsize=13,bg="white",type=type)
  device<-grDevices::dev.cur()
  on.exit({if(device %in% grDevices::dev.list())grDevices::dev.off(device);if(!success&&file.exists(path))unlink(path)})
  draw();grDevices::dev.off(device)
  bytes<-readBin(path,"raw",n=file.info(path)$size)
  if(length(bytes)<32L || !identical(bytes[1:8],as.raw(c(137,80,78,71,13,10,26,10))) ||
     !identical(tail(bytes,12),as.raw(c(0,0,0,0,73,69,78,68,174,66,96,130))))stop("Incomplete PNG: ",path)
  dims<-vapply(list(bytes[17:20],bytes[21:24]),function(b)sum(as.integer(b)*256^(3:0)),numeric(1))
  if(any(abs(dims-round(c(width,height)*dpi))>1))stop("Unexpected PNG dimensions: ",path)
  success<-TRUE;dims
}
dt_audit_publish <- function(stage,destination,files,owned_pattern) {
  # Rendering/validation finishes before any stable copy is touched. Backups allow
  # rollback if publishing fails partway through (e.g. a file is open in Windows).
  dir.create(destination,recursive=TRUE,showWarnings=FALSE)
  if(!dir.exists(destination))stop("Cannot create audit asset directory: ",destination)
  if(identical(normalizePath(stage,winslash="/"),normalizePath(destination,winslash="/")))return(invisible(destination))
  old<-list.files(destination,pattern=owned_pattern,full.names=FALSE)
  backup<-tempfile("audit_backup_",tmpdir=stage);dir.create(backup)
  if(length(old)&&!all(file.copy(file.path(destination,old),backup,overwrite=TRUE)))
    stop("Cannot back up previous audit assets; stable files were not changed.")
  complete<-FALSE
  on.exit({
    if(!complete) {
      unlink(file.path(destination,setdiff(files,old)))
      if(length(old)&&!all(file.copy(file.path(backup,old),destination,overwrite=TRUE)))
        warning("Could not fully restore earlier audit assets; backup retained at ",backup)
    } else unlink(backup,recursive=TRUE)
  })
  for(name in files) {
    src<-file.path(stage,name);dst<-file.path(destination,name)
    if(!file.copy(src,dst,overwrite=TRUE)||unname(tools::md5sum(src))!=unname(tools::md5sum(dst)))
      stop("Could not publish audit asset: ",dst)
  }
  obsolete<-setdiff(old,files)
  if(length(obsolete)&&unlink(file.path(destination,obsolete))!=0L)
    stop("Could not remove obsolete audit continuation pages.")
  complete<-TRUE;invisible(destination)
}
dt_audit_num <- function(x) format(x,big.mark=",",scientific=FALSE,trim=TRUE)
dt_audit_wrap <- function(text,width,size=11.5,bold=FALSE) {
  # Measure actual glyph widths; avoid clipping long descriptions or long tokens.
  measure<-function(s)grid::convertWidth(grid::grobWidth(grid::textGrob(s,
    gp=grid::gpar(fontfamily="sans",fontsize=size,fontface=if(bold)"bold" else "plain"))),"in",valueOnly=TRUE)
  words<-strsplit(gsub("[[:space:]]+"," ",trimws(text))," ",fixed=TRUE)[[1]]
  lines<-character();line<-""
  for(word in words) {
    candidate<-paste(c(if(nzchar(line))line,word),collapse=" ")
    if(measure(candidate)<=width){line<-candidate;next}
    if(nzchar(line)){lines<-c(lines,line);line<-""}
    while(nchar(word)>1L&&measure(word)>width) {
      n<-nchar(word)-1L
      while(n>1L&&measure(substr(word,1,n))>width)n<-n-1L
      lines<-c(lines,substr(word,1,n));word<-substring(word,n+1L)
    }
    line<-word
  }
  paste(c(lines,line),collapse="\n")
}
dt_audit_text <- function(text,x,y,size=11.5,bold=FALSE,col="#25343D",just=c("left","top")) {
  grid::grid.text(text,x=grid::unit(x,"in"),y=grid::unit(1,"npc")-grid::unit(y,"in"),just=just,
    gp=grid::gpar(fontfamily="sans",fontsize=size,fontface=if(bold)"bold" else "plain",col=col,lineheight=1.14))
}

dt_audit_chart <- function(weekly,totals,pop_n,population,cutoff,release) {
  old<-par(no.readonly=TRUE);on.exit(par(old))
  par(mai=c(2.0,1.02,1.35,.25),mgp=c(2.9,.55,0),family="sans",las=1,bty="n",xaxs="i",yaxs="i")
  weeks<-as.Date(weekly$week_start);x<-as.numeric(weeks)
  vals<-c(weekly$dose_1,weekly$dose_2)
  ymax<-max(pretty(c(0,max(1,vals)*1.15),n=5))
  plot(range(x)+c(-5,5),c(0,ymax),type="n",axes=FALSE,xlab="",ylab="")
  ticks<-pretty(c(0,ymax),n=5);ticks<-ticks[ticks>=0 & ticks<=ymax]
  abline(h=ticks,col="#E1E6E9",lwd=.7)
  rect(x-2.7,0,x-.2,weekly$dose_1,col="#1E667D",border=NA)
  rect(x+.2,0,x+2.7,weekly$dose_2,col="#D27C34",border=NA)
  axis(2,at=ticks,labels=dt_n(ticks),tick=FALSE,cex.axis=.9)
  label<-seq(1,length(x),by=max(1L,ceiling(length(x)/11)))
  axis(1,at=x[label],labels=FALSE,tck=-.015,col="#74838C")
  text(x[label],-.025*ymax,dt_label(weeks[label]),srt=45,adj=1,xpd=NA,cex=.88)
  mtext("Dose administrations",side=2,line=3.35,las=0,cex=.95)
  legend("topright",legend=paste0(c("Dose 1: ","Dose 2: "),dt_n(totals$n_administered)),
    fill=c("#1E667D","#D27C34"),border=NA,bty="n",cex=.95,inset=.015)
  grid::pushViewport(grid::viewport())
  dt_audit_text("V04  Weekly dose 1 and dose 2 administrations",.35,.25,17,TRUE,col="#172F3C")
  dt_audit_text(paste0("Moderna P301 | ",if(population=="randomized")"Randomized population" else "All ADSL participants",
    ": ",dt_n(pop_n)),.35,.69,12.2)
  dt_audit_text(paste0(release," data cutoff: ",cutoff," | mRNA-1273 and placebo combined"),.35,.99,11.5)
  dt_audit_text("Monday-start weeks, labelled by their first day. Original blinded study doses only.",.35,5.64,11.2)
  for(i in seq_len(nrow(totals)))dt_audit_text(paste0("Dose ",totals$dose[i],": ",dt_n(totals$n_administered[i]),
    " administrations | ",totals$first_date[i]," to ",totals$last_date[i]),.35,5.97+.28*(i-1),11.2)
  dt_audit_text("Open-label crossover doses are excluded. Zero-dose participants contribute no columns.",.35,6.6,10.6,col="#52616C")
  grid::popViewport()
}
dt_export_audit <- function(weekly,totals,pop_n,population,cutoff,release,source,root,out,
    asset_dir="audit_assets",dpi=300L) {
  stage<-file.path(out,"audit_assets");dir.create(stage,recursive=TRUE,showWarnings=FALSE)
  name<-"V04_weekly_dose_administrations.png"
  dims<-dt_audit_png(file.path(stage,name),function()dt_audit_chart(weekly,totals,pop_n,population,cutoff,release),height=7,dpi=dpi)
  manifest<-data.frame(visual="V04",file=name,width_px=dims[1],height_px=dims[2],dpi=dpi,
    release=release,cutoff=cutoff,population=population,cohort_n=pop_n,source_adsl_md5=source$md5,
    generated_at=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),run_directory=normalizePath(out,winslash="/"),
    md5=unname(tools::md5sum(file.path(stage,name))))
  dt_csv(manifest,file.path(stage,"V04_dose_manifest.csv"))
  dt_csv(weekly,file.path(stage,"V04_weekly_dose_data.csv"))
  writeLines(c("V04: WEEKLY ORIGINAL BLINDED DOSE ADMINISTRATIONS",
    "8.5 x 7 inches, opaque white PNG; default 300 dpi = 2550 x 2100 pixels.",
    "Use the original PNG for Substack/print, not a screenshot. Counts are identical to weekly_doses.csv.",
    "Bars count dose administrations, grouped by dose; both randomized arms combined. Weeks start Monday.",
    "Original dose dates only: dose 1 TR01SDT/DOSE1FL; dose 2 DOSE2DT/DOSE2FL. Crossover is excluded.",
    "Source file/hash and cohort are recorded in the manifest and selected_adsl.csv in the run directory.",
    "Only this script's named V04 files are refreshed in audit_assets. Other figures are retained."),
    file.path(stage,"V04_dose_README.txt"))
  files<-c(name,"V04_dose_manifest.csv","V04_weekly_dose_data.csv","V04_dose_README.txt")
  destination<-dt_path(root,asset_dir)
  dt_audit_publish(stage,destination,files,"^(V04_weekly_dose_administrations[.]png|V04_dose_(manifest[.]csv|README[.]txt)|V04_weekly_dose_data[.]csv)$")
  message("Audit PNG: ",file.path(destination,name)," (",dpi," dpi)")
  list(directory=destination,manifest=manifest)
}


analyse_dose_timeline <- function(root="",release="auto",adsl_file="",population="randomized",output_dir="",
  encoding="auto",numeric_date_origin="",min_randomized=30000L,png_dpi=300L,show_values=FALSE,show_chart=FALSE,
  audit_asset_dir="audit_assets",audit_dpi=300L) {
  if(!population %in% c("randomized","all_adsl"))stop("Population must be randomized or all_adsl.")
  if(!is.numeric(png_dpi)||length(png_dpi)!=1L||!is.finite(png_dpi)||png_dpi<72||png_dpi>1200)stop("PNG DPI must be between 72 and 1200.")
  root<-dt_root(root)
  out<-if(nzchar(output_dir))dt_path(root,output_dir) else file.path(root,"analysis",paste0("dose_timeline_",format(Sys.time(),"%Y%m%d_%H%M%S")))
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  selected<-dt_select(root,release,adsl_file,out,encoding,numeric_date_origin,min_randomized)
  x<-selected$data;if(population=="randomized")x<-x[x$RANDFL=="Y",,drop=FALSE]
  if(!nrow(x))stop("No participants in the selected population.")
  event<-do.call(rbind,lapply(1:2,function(d) {
    flag<-if(d==1)"DOSE1FL" else "DOSE2FL";field<-if(d==1)"TR01SDT" else "DOSE2DT"
    y<-x[x[[flag]]=="Y",,drop=FALSE]
    data.frame(STUDYID=y$STUDYID,USUBJID=y$USUBJID,RANDFL=y$RANDFL,arm=y$TRT01P,
      dose=rep(d,nrow(y)),date=format(y[[field]]),date_field=rep(field,nrow(y)),flag_field=rep(flag,nrow(y)),stringsAsFactors=FALSE)
  }))
  if(!nrow(event))stop("No administered doses in this population.")
  event$week_start<-format(dt_week(as.Date(event$date)))
  event<-event[order(event$date,event$dose,event$USUBJID),];rownames(event)<-NULL
  weeks<-seq(min(as.Date(event$week_start)),max(as.Date(event$week_start)),by="week")
  make_weekly<-function(e,arm="All arms") {
    data.frame(week_start=format(weeks),week_end=format(weeks+6),arm=arm,
      dose_1=tabulate(match(e$week_start[e$dose==1],format(weeks)),nbins=length(weeks)),
      dose_2=tabulate(match(e$week_start[e$dose==2],format(weeks)),nbins=length(weeks)),stringsAsFactors=FALSE)
  }
  weekly<-make_weekly(event);weekly$total_administrations<-weekly$dose_1+weekly$dose_2
  byarm<-do.call(rbind,lapply(sort(unique(event$arm)),function(a)make_weekly(event[event$arm==a,],a)))
  totals<-do.call(rbind,lapply(1:2,function(d) {
    dates<-as.Date(event$date[event$dose==d]);data.frame(dose=d,n_administered=length(dates),
      first_date=if(length(dates))format(min(dates)) else "",last_date=if(length(dates))format(max(dates)) else "")
  }))
  if(sum(weekly$dose_1)!=totals$n_administered[1]||sum(weekly$dose_2)!=totals$n_administered[2])stop("Weekly totals failed reconciliation.")
  if(anyDuplicated(paste(event$USUBJID,event$dose,sep="|")))stop("Duplicate subject/dose.")
  dt_csv(weekly,file.path(out,"weekly_doses.csv"));dt_csv(byarm,file.path(out,"weekly_doses_by_arm.csv"))
  dt_csv(event,file.path(out,"subject_dose_dates.csv"));dt_csv(totals,file.path(out,"dose_totals_and_dates.csv"))
  dt_csv(selected$manifest,file.path(out,"selected_adsl.csv"))
  dt_csv(x[x$DOSE1FL!="Y",c("STUDYID","USUBJID","RANDFL","TRT01P")],file.path(out,"participants_without_dose1.csv"))
  same<-x$DOSE2FL=="Y" & !is.na(x$DOSE2DT)&!is.na(x$TR01SDT)&x$DOSE2DT==x$TR01SDT
  dt_csv(x[which(same),dt_fields],file.path(out,"same_day_doses_for_review.csv"))
  draw<-function()dt_chart(weekly,totals,nrow(x),population,format(x$CUTOFFDT[1]),selected$release,show_values)
  # A4 landscape vector PDF plus SVG and 300-dpi PNG; no browser/LaTeX required.
  pdfpath<-file.path(out,"weekly_doses.pdf")
  if(capabilities("cairo"))grDevices::cairo_pdf(pdfpath,width=297/25.4,height=210/25.4,family="sans")
  else grDevices::pdf(pdfpath,width=297/25.4,height=210/25.4,paper="special",family="Helvetica")
  tryCatch(draw(),finally=grDevices::dev.off())
  pngpath<-file.path(out,"weekly_doses.png")
  grDevices::png(pngpath,width=297/25.4,height=210/25.4,units="in",res=png_dpi,
    type=if(capabilities("cairo"))"cairo" else getOption("bitmapType"))
  tryCatch(draw(),finally=grDevices::dev.off())
  if(capabilities("cairo")) {
    grDevices::svg(file.path(out,"weekly_doses.svg"),width=297/25.4,height=210/25.4,family="sans")
    tryCatch(draw(),finally=grDevices::dev.off())
  }
  readme<-c("MODERNA P301 - WEEKLY DOSE ADMINISTRATIONS - v1.1.0", "",
    paste("Population:",population,"; ADSL participants:",nrow(x)),
    paste("Release:",selected$release,"; source cutoff:",format(x$CUTOFFDT[1])),
    paste("Selected file:",selected$manifest$path),paste("MD5:",selected$manifest$md5),"",
    "Dose 1 = DOSE1FL Y with TR01SDT; dose 2 = DOSE2FL Y with DOSE2DT. Flag/date inconsistencies stop the run.",
    "These are the original blinded study injections, including placebo. Open-label crossover doses (TR02SDT, OLDOS2DT) and outside vaccinations are not counted.",
    "The default population is the sponsor's RANDFL Y set, matching the BLA randomized population. It is not the safety set. all_adsl includes dated doses in all ADSL rows, including RANDFL N.",
    "Each participant contributes at most one administration per dose. The total across both doses is administrations, not unique participants.",
    "Weeks are Monday-Sunday calendar intervals based on the actual recorded dosing dates. All intervening weeks, including zero-count weeks, are retained.",
    "No dose date is inferred from a planned visit, a 28-day interval, or a participant's last visit. Missing dates are not imputed.",
    "Automatic selection prefers a complete production ADSL and its latest CUTOFFDT. EUA is a fallback and is never pooled with production. Conflicting latest-cutoff inputs require an explicit PHMPT_ADSL_FILE.",
    "CSV dates must be ISO, English SAS-style DDMMMYYYY, or use an explicitly documented numeric origin. Haven dates and JMP date metadata are decoded directly.",
    "CSV/JMP processing and plots use base R. XPT/SAS7BDAT require haven. The embedded JMP reader supports the released P301 fixed-width/Float64 format and stops on unsupported layouts.",
    "AUDIT OUTPUT: audit_assets/V04_weekly_dose_administrations.png, with data CSV, provenance manifest and README; dedicated print/article layout (default 2550 px wide at 300 dpi).",
    "Outputs: weekly_doses.pdf/.png/.svg; weekly_doses.csv; weekly_doses_by_arm.csv; subject_dose_dates.csv; dose_totals_and_dates.csv; source and exception logs.",
    "Protocol/metadata reference: production P301 ADaM define.xml, ADSL TR01SDT, DOSE2DT, DOSE1FL, DOSE2FL and RANDFL derivations.","",
    capture.output(print(totals,row.names=FALSE)))
  writeLines(readme,file.path(out,"README.txt"));capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
  audit_assets<-dt_export_audit(weekly,totals,nrow(x),population,format(x$CUTOFFDT[1]),selected$release,
    selected$manifest,root,out,audit_asset_dir,audit_dpi)
  if(show_chart)draw()
  message("Chart: ",pdfpath);message("Weekly CSV: ",file.path(out,"weekly_doses.csv"));print(totals,row.names=FALSE)
  invisible(list(output_dir=out,weekly=weekly,totals=totals,events=event,selected=selected$manifest,audit_assets=audit_assets))
}

# ---- EMBEDDED P301 JMP READER ----
# Reused from the validated project reader; no external script is sourced.
dt_read_jmp <- function(path, fields) {
  size <- file.info(path)$size
  if (is.na(size) || size < 64 || size > .Machine$integer.max)
    stop("Unsupported JMP file size: ", path)
  b <- readBin(path, "raw", n = size)
  take <- function(buf, at, n) {
    if (length(at) != 1L || length(n) != 1L || !is.finite(at) || !is.finite(n) ||
        at < 1 || n < 0 || at + n - 1 > length(buf)) stop("JMP byte range outside file/payload.")
    if (n == 0) raw() else buf[seq.int(at, length.out = n)]
  }
  uint <- function(buf, at, width) {
    z <- sum(as.integer(take(buf, at, width)) * 256^(seq_len(width) - 1L))
    if (!is.finite(z) || z > 2^53 - 1) stop("JMP integer exceeds exact R representation.")
    z
  }
  txt <- function(z) {
    nonzero <- which(z != as.raw(0))
    if (!length(nonzero)) return("")
    s <- rawToChar(z[seq_len(tail(nonzero, 1L))])
    if (!validUTF8(s)) stop("Invalid UTF-8 in JMP string.")
    Encoding(s) <- "UTF-8"; s
  }
  if (!identical(b[1:8], as.raw(c(255,255,0,0,7,0,0,0))))
    stop("Unsupported JMP signature. Supply an exported CSV/XPT/SAS7BDAT instead.")
  header <- grepRaw(as.raw(c(7,0,8,0,0,0)), b, fixed = TRUE)[1]
  if (is.na(header) || header < 33) stop("JMP dimensions not found.")
  nr <- uint(b, header - 32, 8); nc <- uint(b, header - 24, 4)
  if (nr < 1 || nr > .Machine$integer.max || nc < 1 || nc > 10000) stop("Invalid JMP dimensions.")
  offsets <- NULL
  # Column directory precedes column payloads; validate every pointer and name.
  starts <- grepRaw(as.raw(c(255,255)), b, all = TRUE, fixed = TRUE)
  starts <- starts[starts > 2]
  for (at in starts) {
    p <- at
    while (p <= length(b) && b[p] == as.raw(255)) p <- p + 1L
    p <- p + 10L
    if (p + 15 > length(b)) next
    if (uint(b, p, 4) + uint(b, p + 4, 4) != nc) next
    p <- p + 16 + 6 * nc + 28
    found <- tryCatch({
      while (identical(take(b, p, 2), as.raw(c(253,255))) ||
             identical(take(b, p, 2), as.raw(c(254,255))) ||
             identical(take(b, p, 2), as.raw(c(255,255)))) p <- p + 10 + uint(b, p + 2, 8)
      if (uint(b, p, 4) != nc) stop("Directory column count mismatch.")
      off <- vapply(seq_len(nc), function(i) uint(b, p + 4 + (i - 1) * 8, 8) + 1, numeric(1))
      if (any(diff(off) <= 0) || any(off <= p) || any(off + 8 > length(b))) stop("Invalid JMP column pointers.")
      nm <- vapply(off, function(o) txt(take(b, o + 2, uint(b, o, 2))), character(1))
      if (anyDuplicated(nm) || any(!nzchar(nm))) stop("Invalid JMP column names.")
      list(offsets = off, names = nm)
    }, error = function(e) NULL)
    if (!is.null(found)) { offsets <- found$offsets; column_names <- found$names; break }
  }
  if (is.null(offsets)) stop("JMP column directory not found or unsupported.")
  missing <- setdiff(fields, column_names)
  if (length(missing)) stop("JMP missing required columns: ", paste(missing, collapse = ", "))
  selected <- c(fields, intersect("COUNTRY", setdiff(column_names, fields)))
  values <- setNames(vector("list", length(selected)), selected)
  for (name in selected) {
    i <- match(name, column_names); p <- offsets[i] + 2 + uint(b, offsets[i], 2)
    dt <- as.integer(take(b, p, 6)); p <- p + 6
    end <- if (i < nc) offsets[i + 1] - 1 else length(b)
    if (!dt[1] %in% c(9L,10L)) stop("Unsupported JMP encoding in ", name, "; export this table to CSV.")
    h <- take(b, p, min(2048, end - p + 1))
    marker <- grepRaw(as.raw(c(239,190,254,202)), h, fixed = TRUE)[1]
    if (is.na(marker)) stop("JMP compression header missing in ", name)
    q <- p + marker - 1 + 4
    packed <- uint(b, q, 8); unpacked <- uint(b, q + 8, 8)
    if (unpacked > .Machine$integer.max || q + 16 + packed - 1 > end)
      stop("Invalid JMP compressed block length in ", name)
    payload <- withCallingHandlers(memDecompress(take(b, q + 16, packed), type = "gzip"),
      warning = function(w) stop("JMP decompression warning in ", name, ": ", conditionMessage(w)))
    if (length(payload) != unpacked) stop("JMP decompressed size mismatch in ", name)
    if (dt[1] == 9L && dt[2] %in% c(1L,2L) && all(dt[3:4] == 0L) && dt[5] > 0L) {
      width <- dt[5]
      if (length(payload) != nr * width) stop("Unsupported fixed-width JMP payload in ", name)
      values[[name]] <- vapply(seq_len(nr), function(j) txt(take(payload, (j - 1) * width + 1, width)), character(1))
    } else if (dt[1] == 10L && dt[6] == 8L && dt[4] == dt[5] && dt[4] %in% c(0L,3L,99L,114L)) {
      if (length(payload) != nr * 8) stop("JMP numeric row count mismatch in ", name)
      z <- readBin(payload, "double", n = nr, size = 8, endian = "little")
      if (any(is.infinite(z))) stop("Non-finite JMP number in ", name)
      missing <- is.na(z)
      if (dt[4] == 114L) {
        if (any(abs(z[!is.na(z)] / 86400 - round(z[!is.na(z)] / 86400)) > 1e-8))
          stop("JMP date includes a time component in ", name)
        z <- format(as.Date(z / 86400, origin = "1904-01-01"), "%Y-%m-%d")
      } else z <- as.character(z)
      z[missing | is.na(z)] <- ""; values[[name]] <- z
    } else stop("Unsupported JMP column type in ", name, ": ", paste(dt, collapse = ", "))
  }
  as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
}

# ---- RUN: select and run the entire script in RStudio ----
if(isTRUE(getOption("phmpt.autorun",TRUE))) {
  analyse_dose_timeline(root=PHMPT_ROOT,release=PHMPT_RELEASE,adsl_file=PHMPT_ADSL_FILE,
    population=PHMPT_POPULATION,output_dir=PHMPT_OUTPUT_DIR,encoding=PHMPT_CSV_ENCODING,
    numeric_date_origin=PHMPT_NUMERIC_DATE_ORIGIN,min_randomized=PHMPT_MIN_RANDOMIZED,
    png_dpi=PHMPT_PNG_DPI,show_values=PHMPT_SHOW_VALUES,show_chart=PHMPT_SHOW_CHART,
    audit_asset_dir=PHMPT_AUDIT_ASSET_DIR,audit_dpi=PHMPT_AUDIT_DPI)
}
