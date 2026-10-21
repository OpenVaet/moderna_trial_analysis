# analyse_site_symptom_recording.R | v1.0.0 | 2026-10-21
# Standalone base-R analysis of exports from analyse_weekly_symptoms_ae.R.
# Run the entire file in RStudio from the project root. No package installation.
# Numerator: AE_ONLY + BOTH. Denominator: all distinct subject x symptom-group x week units.
# These are recording-source shares, NOT percentages of subjects or diagnoses.

PHMPT_ROOT <- getOption("phmpt.root", getwd())
PHMPT_WEEKLY_DIR <- getOption("phmpt.weekly_dir", "") # folder or exported ZIP; blank discovers latest complete export
PHMPT_BACKGROUND_FILE <- getOption("phmpt.background_file", "") # site_weekly.csv or site_background_used.csv
PHMPT_REGION_FILE <- getOption("phmpt.region_file", "") # companion site_region_coverage.csv for heatmap exports
PHMPT_OUTPUT_DIR <- getOption("phmpt.site_recording_output", "")
PHMPT_WINDOW <- getOption("phmpt.recording_window", "BLINDED") # chart window; also supports EUA / FULL
PHMPT_SMALL_N <- getOption("phmpt.recording_small_n", 5L)
PHMPT_INDIVIDUAL_PDFS <- isTRUE(getOption("phmpt.individual_site_pdfs", FALSE))
PHMPT_MAKE_PDF <- isTRUE(getOption("phmpt.recording_pdf", TRUE))

sr_csv <- function(x, p) utils::write.csv(x, p, row.names=FALSE, na="", fileEncoding="UTF-8")
sr_read <- function(p) {
  x <- utils::read.csv(p, stringsAsFactors=FALSE, check.names=FALSE, fileEncoding="UTF-8-BOM")
  names(x) <- sub("^\ufeff", "", names(x)); x
}
sr_require <- function(x, cols, where) {
  m <- setdiff(cols, names(x)); if(length(m)) stop(where, ": missing columns: ", paste(m, collapse=", "))
}
sr_path <- function(p, root) if(!nzchar(p)) "" else if(grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)",p)) p else file.path(root,p)
sr_bool <- function(x) {
  y <- toupper(as.character(x)); if(any(!y %in% c("TRUE","FALSE","1","0"))) stop("Invalid boolean values.")
  y %in% c("TRUE","1")
}
sr_date <- function(x, label) {
  d <- as.Date(as.character(x)); if(any(is.na(d) & !is.na(x) & nzchar(as.character(x)))) stop("Invalid dates: ",label)
  d
}
sr_week <- function(d) d - (as.integer(format(d,"%u"))-1L)
sr_p <- function(p) ifelse(is.na(p),"not estimable",ifelse(p<1e-6,"< 0.000001",formatC(p,digits=4,format="f")))

sr_discover <- function(root, explicit, required) {
  p <- sr_path(explicit, root)
  if(nzchar(p)) {
    if(grepl("[.]zip$",p,ignore.case=TRUE)) {
      if(!file.exists(p)) stop("Weekly ZIP not found: ",p)
      listing <- utils::unzip(p,list=TRUE)$Name
      if(any(grepl("(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))",listing))) stop("Unsafe paths in weekly ZIP.")
      dest <- tempfile("weekly_recording_");dir.create(dest)
      utils::unzip(p,files=listing[grepl("[.](csv|txt)$",listing,ignore.case=TRUE)],exdir=dest)
      candidates <- unique(dirname(list.files(dest,pattern="^subject_followup_and_doses[.]csv$",recursive=TRUE,full.names=TRUE)))
    } else candidates <- p
  } else {
    base <- if(dir.exists(file.path(root,"analysis"))) file.path(root,"analysis") else root
    candidates <- unique(dirname(list.files(base,pattern="^subject_followup_and_doses[.]csv$",recursive=TRUE,full.names=TRUE)))
    if(all(file.exists(file.path(root,required)))) candidates <- unique(c(root,candidates))
  }
  candidates <- candidates[vapply(candidates,function(d)all(file.exists(file.path(d,required))),logical(1))]
  if(!length(candidates)) stop("No complete weekly export found. Set PHMPT_WEEKLY_DIR to its folder or ZIP. Required: ",paste(required,collapse=", "))
  stamp <- vapply(candidates,function(d)as.numeric(file.info(file.path(d,"window_totals.csv"))$mtime),numeric(1))
  candidates <- candidates[order(stamp,decreasing=TRUE)]
  message("Weekly export: ",candidates[1L])
  candidates[1L]
}

# Cluster sandwich with HC0 scores and G/(G-1) correction; t reference G-1.
# Site clusters contain all repeated observations from their participants.
sr_model <- function(d, formula, term, label, cluster="SITEID", logistic=FALSE) {
  empty <- data.frame(model=label,term=term,scale=if(logistic)"odds_ratio" else "percentage_points",
    estimate=NA_real_,lower_95=NA_real_,upper_95=NA_real_,p_value=NA_real_,n_units=nrow(d),
    n_subjects=length(unique(d$subject)),n_clusters=length(unique(d[[cluster]])),cluster=cluster,status="not estimable")
  tryCatch({
    if(anyNA(model.frame(formula,d,na.action=na.pass)) || anyNA(d[[cluster]])) stop("Missing model covariates")
    fit <- if(logistic) glm(formula,d,family=binomial(),x=TRUE,y=TRUE,control=glm.control(maxit=100)) else lm(formula,d,x=TRUE,y=TRUE)
    if(logistic && (!fit$converged || fit$boundary)) stop("Logistic convergence/boundary failure")
    b <- coef(fit); keep <- which(!is.na(b)); X <- fit$x[,keep,drop=FALSE]; b <- b[keep]
    if(!term %in% names(b)) stop("Requested coefficient not estimable")
    mu <- fitted(fit); w <- if(logistic) mu*(1-mu) else rep(1,nrow(d))
    bread <- solve(crossprod(X,X*w))
    score <- rowsum(X*as.numeric(fit$y-mu),as.character(d[[cluster]]),reorder=FALSE)
    G <- nrow(score); if(G<3L) stop("Fewer than three clusters")
    V <- bread %*% crossprod(score) %*% bread * G/(G-1)
    se <- sqrt(diag(V)[match(term,names(b))]); z <- unname(b[term]);
    if(!is.finite(se)||se<=0) stop("Zero/invalid robust standard error")
    ci <- z+c(-1,1)*qt(.975,G-1)*se
    empty$estimate <- if(logistic) exp(z) else 100*z
    empty$lower_95 <- if(logistic) exp(ci[1]) else 100*ci[1]
    empty$upper_95 <- if(logistic) exp(ci[2]) else 100*ci[2]
    empty$p_value <- 2*pt(-abs(z/se),G-1);empty$status <- "OK";empty
  },error=function(e){empty$status<-conditionMessage(e);empty})
}

sr_membership <- function(path, roster) {
  d <- sr_read(path);sr_require(d,c("subject","arm","symptom_group","week_start","ae","covid_form","route"),path)
  if(anyDuplicated(d[c("subject","symptom_group","week_start")]))stop("Duplicate subject/group/week units: ",path)
  if(anyNA(d[c("subject","arm","symptom_group","week_start","route")])) stop("Missing membership identifiers")
  d$ae<-sr_bool(d$ae);d$covid_form<-sr_bool(d$covid_form)
  expected<-ifelse(d$ae & d$covid_form,"BOTH",ifelse(d$ae,"AE_ONLY",ifelse(d$covid_form,"COVID_ONLY","INVALID")))
  if(any(d$route!=expected)||any(expected=="INVALID"))stop("Route/flag inconsistency")
  i<-match(d$subject,roster$subject)
  if(anyNA(i)||any(d$arm!=roster$arm[i]))stop("Subject/arm mismatch with selected roster")
  d$SITEID<-roster$SITEID[i];d$week_start<-sr_date(d$week_start,"membership week")
  if(any(d$week_start!=sr_week(d$week_start)))stop("Weeks are not Mondays")
  d$vaccine<-as.integer(d$arm=="mRNA-1273");d$ae_any<-as.integer(d$ae)
  d$ae_only<-as.integer(d$route=="AE_ONLY");d
}
sr_totals <- function(d, window) {
  do.call(rbind,lapply(c("mRNA-1273","Placebo"),function(a){
    z<-d[d$arm==a,];n<-nrow(z)
    data.frame(window=window,arm=a,subjects=length(unique(z$subject)),covid_only=sum(z$route=="COVID_ONLY"),
      ae_only=sum(z$route=="AE_ONLY"),both=sum(z$route=="BOTH"),total_units=n,ae_any=sum(z$ae),
      pct_covid_only=100*sum(z$route=="COVID_ONLY")/n,pct_ae_only=100*sum(z$route=="AE_ONLY")/n,
      pct_both=100*sum(z$route=="BOTH")/n,pct_ae_any=100*sum(z$ae)/n)
  }))
}

sr_background <- function(root, explicit, region_file, sites, weeks) {
  p<-sr_path(explicit,root)
  if(!nzchar(p)) {
    p<-file.path(root,"analysis","background_covid","site_weekly.csv")
    if(!file.exists(p)) {
      base<-if(dir.exists(file.path(root,"analysis")))file.path(root,"analysis") else root
      pp<-list.files(base,pattern="^site_background_used[.]csv$",recursive=TRUE,full.names=TRUE)
      # Generated outputs of THIS script are never treated as source backgrounds.
      pp<-pp[!file.exists(file.path(dirname(pp),"site_recording_parameters.csv"))]
      if(length(pp))p<-pp[order(file.info(pp)$mtime,decreasing=TRUE)][1L]
    }
  }
  if(!file.exists(p))stop("County background missing. Set PHMPT_BACKGROUND_FILE to analysis/background_covid/site_weekly.csv or a heatmap export's site_background_used.csv.")
  message("Background: ",p);b<-sr_read(p);mp<-NULL;rp<-""
  if(all(c("SITEID","week_start","county_positivity_pct") %in% names(b))) {
    sr_require(b,c("county_fips","county_name","state_name"),p)
    b<-unique(b);if(anyDuplicated(b[c("SITEID","week_start")]))stop("Conflicting background site/week duplicates")
    b$week<-sr_date(b$week_start,"background week");b$positivity<-as.numeric(b$county_positivity_pct)
    b$region_key<-paste0("county_",as.character(b$county_fips));b$region_level<-"county"
    b$region_name<-paste(b$county_name,b$state_name,sep=", ")
    mp<-unique(b[c("SITEID","region_key","region_level","region_name")])
  } else {
    sr_require(b,c("SITEID","week","region_key","positivity"),p)
    rp<-sr_path(region_file,root);if(!nzchar(rp))rp<-file.path(dirname(p),"site_region_coverage.csv")
    if(!file.exists(rp))stop("Heatmap background requires companion site_region_coverage.csv")
    mp<-sr_read(rp);sr_require(mp,c("SITEID","region_key","region_name","region_level"),rp)
    mp<-unique(mp[c("SITEID","region_key","region_level","region_name")])
    b$week<-sr_date(b$week,"background week");b$positivity<-as.numeric(b$positivity)
    # The user requested COUNTY positivity. Never relabel state fallback as county.
    k<-match(b$SITEID,mp$SITEID)
    if(anyNA(k)||any(b$region_key!=mp$region_key[k]))stop("Background region-map mismatch")
    b$positivity[mp$region_level[k]!="county"]<-NA_real_
  }
  if(anyDuplicated(mp$SITEID))stop("More than one region per site")
  b<-b[b$SITEID %in% sites & b$week %in% weeks,c("SITEID","region_key","week","positivity")]
  if(anyDuplicated(b[c("SITEID","week")]))stop("Duplicate background site/week values")
  if(any(b$week!=sr_week(b$week)))stop("Background weeks must start Mondays")
  if(any(b$positivity<0|b$positivity>100,na.rm=TRUE))stop("Positivity must be in percent, 0-100")
  geo<-unique(b[c("region_key","week","positivity")]);if(anyDuplicated(geo[c("region_key","week")]))stop("Inconsistent values for a shared county/week")
  grid<-expand.grid(SITEID=sites,week=weeks,stringsAsFactors=FALSE)
  grid<-merge(grid,b,by=c("SITEID","week"),all.x=TRUE,sort=FALSE)
  m<-match(grid$SITEID,mp$SITEID);grid$region_name<-mp$region_name[m];grid$region_level<-mp$region_level[m]
  grid$region_name[is.na(grid$region_name)]<-"County background unavailable"
  list(data=grid,map=mp,file=p,region_file=rp)
}

# Recount the same source records using only the AE onset week. COVID assessments
# retain their actual assessment weeks. This removes repeated interval weeks,
# but is not a reconstruction of illness episodes or a proof of reclassification.
sr_onset <- function(folder, window, original, roster) {
  lp<-file.path(folder,paste0(window,"_source_event_week_links.csv"));ep<-file.path(folder,"symptom_event_evidence.csv")
  if(!all(file.exists(c(lp,ep))))return(NULL)
  l<-sr_read(lp);e<-sr_read(ep)
  sr_require(l,c("event_id","source","subject","arm","symptom_group","week_start"),lp)
  sr_require(e,c("event_id","start_date","subject","source"),ep)
  if(anyDuplicated(e$event_id))stop("Duplicate event IDs in evidence")
  i<-match(l$event_id,e$event_id)
  if(anyNA(i)||any(l$subject!=e$subject[i])||any(l$source!=e$source[i]))stop("Event link/evidence mismatch")
  l$week_start<-sr_date(l$week_start,"event link week")
  onset<-sr_week(sr_date(e$start_date[i],"event onset"))
  l<-l[l$source!="AE" | (!is.na(onset)&l$week_start==onset),]
  l$ae<-as.integer(l$source=="AE");l$covid_form<-as.integer(l$source=="COVID_FORM")
  if(any(l$ae+l$covid_form!=1))stop("Unexpected source in onset recount")
  z<-aggregate(cbind(ae,covid_form)~subject+arm+symptom_group+week_start,l,max)
  z$ae<-z$ae>0;z$covid_form<-z$covid_form>0
  z$route<-ifelse(z$ae&z$covid_form,"BOTH",ifelse(z$ae,"AE_ONLY","COVID_ONLY"))
  k<-match(z$subject,roster$subject);z$SITEID<-roster$SITEID[k]
  z$vaccine<-as.integer(z$arm=="mRNA-1273");z$ae_any<-as.integer(z$ae);z$ae_only<-as.integer(z$route=="AE_ONLY")
  key<-function(q)paste(q$subject,q$symptom_group,q$week_start,sep="|")
  j<-match(key(z),key(original));if(anyNA(j))stop("Onset recount introduced an unknown unit")
  z
}

sr_pdf <- function(path) {
  if(capabilities("cairo"))grDevices::cairo_pdf(path,width=11.69,height=8.27,onefile=TRUE,family="sans")
  else grDevices::pdf(path,width=11.69,height=8.27,onefile=TRUE,family="Helvetica",useDingbats=FALSE)
}
sr_page <- function(site, cells, bg, roster, window, small_n, page, pages) {
  colors<-c("mRNA-1273"="#176E9E","Placebo"="#B65A14")
  z<-cells[cells$SITEID==site,];g<-bg[bg$SITEID==site,];g<-g[order(g$week),]
  weeks<-sort(unique(z$week_start));region<-unique(g$region_name);if(!length(region))region<-"County background unavailable"
  par(mfrow=c(1,1),mar=c(4.7,4.7,7.0,2.2),oma=c(4.0,0,0,0),family="sans",las=1,xaxs="i",yaxs="i")
  plot(weeks,rep(NA_real_,length(weeks)),type="n",xlim=range(weeks)+c(-3,3),ylim=c(0,103),
    xaxt="n",yaxt="n",xlab="",ylab="Percent",bty="n")
  abline(h=seq(0,100,20),col="#E2E7EA",lwd=.7)
  axis(2,at=seq(0,100,20),labels=paste0(seq(0,100,20),"%"),col=NA,col.ticks=NA,cex.axis=.85)
  ticks<-weeks[seq(1,length(weeks),by=3)]
  axis.Date(1,at=ticks,format="%d %b",cex.axis=.85,col="#6B7280")
  lines(g$week,g$positivity,col="#333333",lty=2,lwd=2)
  for(a in names(colors)) {
    q<-z[z$arm==a,];q<-q[order(q$week_start),]
    lines(q$week_start,q$pct_ae_any,col=colors[a],lwd=1.6)
    points(q$week_start,q$pct_ae_any,pch=ifelse(q$total_units<small_n,1,16),col=colors[a],cex=.65)
  }
  mtext(paste(site,region,sep=" - "),3,line=5.4,adj=0,font=2,cex=1.25)
  mtext(paste0(window," window | AE-recorded share of weekly symptom observations"),3,line=3.9,adj=0,cex=.96)
  cohort<-table(factor(roster$arm[roster$SITEID==site],levels=names(colors)))
  mtext(sprintf("Randomized at site: %s mRNA-1273 / %s placebo. All compatible symptoms included, regardless of PCR.",cohort[1],cohort[2]),3,line=2.6,adj=0,cex=.81)
  legend("top",inset=c(0,-.063),legend=c("mRNA-1273: AE share","Placebo: AE share","County test positivity"),
    col=c(colors,"#333333"),lty=c(1,1,2),lwd=c(1.7,1.7,2),pch=c(16,16,NA),horiz=TRUE,bty="n",xpd=TRUE,cex=.84,x.intersp=1.3)
  mtext("Week starting Monday (2020-2021)",1,line=3.2,cex=.86)
  s<-aggregate(cbind(ae_any,total_units)~arm,z,sum)
  summary<-vapply(names(colors),function(a){q<-s[s$arm==a,];sprintf("%s: %s / %s (%.1f%%)",a,q$ae_any,q$total_units,100*q$ae_any/q$total_units)},character(1))
  mtext(paste("Entire window:",paste(summary,collapse="; ")),1,outer=TRUE,line=.5,adj=.055,cex=.81)
  mtext("AE share = (AE-only + both-source units) / all subject x symptom-group x week units. Both-source co-occurrence is not proven reclassification.",1,outer=TRUE,line=1.65,adj=.055,cex=.71)
  mtext(sprintf("Open circles: fewer than %s symptom units; gaps: no observed symptoms. County positivity is an external measure on the same 0-100%% axis.",small_n),1,outer=TRUE,line=2.65,adj=.055,cex=.71)
  mtext(sprintf("Site %s of %s",page,pages),3,line=3.9,adj=1,cex=.7)
}

analyse_site_symptom_recording <- function(root=PHMPT_ROOT,weekly_dir=PHMPT_WEEKLY_DIR,
  background_file=PHMPT_BACKGROUND_FILE,region_file=PHMPT_REGION_FILE,output_dir=PHMPT_OUTPUT_DIR,
  window=PHMPT_WINDOW,small_n=PHMPT_SMALL_N,individual_pdfs=PHMPT_INDIVIDUAL_PDFS,make_pdf=PHMPT_MAKE_PDF) {
  root<-normalizePath(root,mustWork=TRUE);window<-toupper(window)
  if(!window %in% c("BLINDED","EUA","FULL"))stop("Invalid window")
  required<-c("subject_followup_and_doses.csv","window_totals.csv",paste0(window,"_subject_symptom_weeks.csv"))
  folder<-sr_discover(root,weekly_dir,required)
  out<-sr_path(output_dir,root);if(!nzchar(out))out<-file.path(root,"analysis",paste0("site_symptom_recording_",format(Sys.time(),"%Y%m%d_%H%M%S")))
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  roster<-sr_read(file.path(folder,"subject_followup_and_doses.csv"))
  sr_require(roster,c("subject","USUBJID","arm","RANDFL"),"subject roster")
  if(anyDuplicated(roster$subject)||anyNA(roster$subject)||any(roster$RANDFL!="Y"))stop("Invalid randomized roster")
  if(any(!roster$arm %in% c("mRNA-1273","Placebo")))stop("Unrecognized treatment arm")
  # Validate the P301 site convention against the explicit USUBJID site segment.
  explicit_site<-sub("^.*-(US[0-9]{3})-[0-9]+$","\\1",roster$USUBJID)
  prefix<-substr(roster$subject,1,5)
  if(any(!grepl("^US[0-9]{3}$",explicit_site))||any(prefix!=explicit_site))stop("Cannot validate site IDs from subject and USUBJID")
  if("SITEID" %in% names(roster) && any(roster$SITEID!=explicit_site))stop("SITEID/USUBJID mismatch")
  roster$SITEID<-explicit_site;sites<-sort(unique(roster$SITEID))
  windows<-c("BLINDED","EUA","FULL");windows<-windows[file.exists(file.path(folder,paste0(windows,"_subject_symptom_weeks.csv")))]
  all_data<-setNames(lapply(windows,function(w)sr_membership(file.path(folder,paste0(w,"_subject_symptom_weeks.csv")),roster)),windows)
  totals<-do.call(rbind,lapply(windows,function(w)sr_totals(all_data[[w]],w)))
  control<-sr_read(file.path(folder,"window_totals.csv"))
  for(w in windows)for(a in c("mRNA-1273","Placebo")) {
    q<-totals[totals$window==w & totals$arm==a,];r<-control[control$window==w & control$arm==a,]
    if(nrow(r)!=1L || any(as.numeric(q[c("covid_only","ae_only","both","total_units")])!=as.numeric(r[c("covid_only","ae_only","both","symptom_group_weeks")])))stop("Counts do not reconcile to window_totals.csv")
  }
  sr_csv(totals,file.path(out,"recording_source_totals.csv"))
  models<-list();add<-function(x)models[[length(models)+1L]]<<-x
  for(w in intersect(c("BLINDED","EUA"),windows)) {
    d<-all_data[[w]];message("Comparisons: ",w)
    add(sr_model(d,ae_any~vaccine,"vaccine",paste0(w,"_raw_AE_any")))
    add(sr_model(d,ae_only~vaccine,"vaccine",paste0(w,"_raw_AE_only")))
    add(sr_model(d,ae_any~vaccine+factor(SITEID)+factor(week_start)+factor(symptom_group),"vaccine",paste0(w,"_adjusted_AE_any"),logistic=TRUE))
  }
  d<-all_data[[window]]
  # Onset sensitivity is calculated for BLINDED (when available), independent of chart window.
  onset<-if("BLINDED" %in% windows)sr_onset(folder,"BLINDED",all_data$BLINDED,roster) else NULL
  if(!is.null(onset)) {
    sr_csv(sr_totals(onset,"BLINDED_AE_onset_only"),file.path(out,"onset_sensitivity_totals.csv"))
    add(sr_model(onset,ae_any~vaccine,"vaccine","BLINDED_AE_onset_only_raw"))
    add(sr_model(onset,ae_any~vaccine+factor(SITEID)+factor(week_start)+factor(symptom_group),"vaccine","BLINDED_AE_onset_only_adjusted",logistic=TRUE))
  }
  # Subject-weighted sensitivity: each participant with >=1 symptom unit has equal weight.
  if("BLINDED" %in% windows) {
    first<-aggregate(ae_any~subject+arm+vaccine+SITEID,all_data$BLINDED,mean)
    sr_csv(aggregate(ae_any~arm,first,mean),file.path(out,"participant_mean_AE_share.csv"))
    add(sr_model(first,ae_any~vaccine,"vaccine","BLINDED_equal_participant_weight"))
    if("has_blind_boundary" %in% names(roster)) {
      b<-all_data$BLINDED;known<-sr_bool(roster$has_blind_boundary)[match(b$subject,roster$subject)]
      sr_csv(sr_totals(b[known,],"BLINDED_documented_boundary_only"),file.path(out,"documented_boundary_totals.csv"))
      add(sr_model(b[known,],ae_any~vaccine,"vaccine","BLINDED_documented_boundary_only"))
    }
  }
  weekly_control<-file.path(folder,paste0(window,"_weekly_counts.csv"))
  if(file.exists(weekly_control)) {
    wc<-sr_read(weekly_control);sr_require(wc,"week_start",weekly_control)
    ww<-sr_date(wc$week_start,"weekly counts");weeks<-seq(min(ww),max(ww),by="7 days")
    if(any(!d$week_start %in% weeks))stop("Membership outside exported weekly range")
  } else weeks<-seq(min(d$week_start),max(d$week_start),by="7 days")
  bg<-sr_background(root,background_file,region_file,sites,weeks)
  n<-aggregate(cbind(ae_any,ae_only,covid_form)~SITEID+arm+week_start,d,sum)
  den<-aggregate(subject~SITEID+arm+week_start,d,length);names(den)[4]<-"total_units"
  un<-aggregate(subject~SITEID+arm+week_start,d,function(x)length(unique(x)));names(un)[4]<-"unique_subjects"
  cells<-expand.grid(SITEID=sites,arm=c("mRNA-1273","Placebo"),week_start=weeks,stringsAsFactors=FALSE)
  for(x in list(n,den,un))cells<-merge(cells,x,by=c("SITEID","arm","week_start"),all.x=TRUE,sort=FALSE)
  for(k in c("ae_any","ae_only","covid_form","total_units","unique_subjects"))cells[[k]][is.na(cells[[k]])]<-0L
  cells$both<-cells$ae_any-cells$ae_only;cells$covid_only<-cells$total_units-cells$ae_any
  cells$pct_ae_any<-ifelse(cells$total_units>0,100*cells$ae_any/cells$total_units,NA_real_)
  cells$small_denominator<-cells$total_units>0 & cells$total_units<small_n
  bgc<-bg$data;names(bgc)[names(bgc)=="week"]<-"week_start"
  cells<-merge(cells,bgc,by=c("SITEID","week_start"),all.x=TRUE,sort=FALSE)
  cells<-cells[order(cells$SITEID,cells$arm,cells$week_start),];cells$window<-window
  sr_csv(cells,file.path(out,"site_arm_weekly_recording.csv"))
  sr_csv(bg$data,file.path(out,"site_background_used.csv"));sr_csv(bg$map,file.path(out,"site_region_coverage.csv"))
  site_tests<-list();group_tests<-list()
  if("BLINDED" %in% windows) {
    b<-all_data$BLINDED
    for(s in sites) {
      z<-b[b$SITEID==s,];r<-sr_model(z,ae_any~vaccine,"vaccine",s,cluster="subject")
      tt<-sr_totals(z,"BLINDED");r$v_ae_any<-tt$ae_any[1];r$v_units<-tt$total_units[1];r$p_ae_any<-tt$ae_any[2];r$p_units<-tt$total_units[2]
      r$SITEID<-s;site_tests[[length(site_tests)+1L]]<-r
    }
    st<-do.call(rbind,site_tests);st$p_BH<-p.adjust(st$p_value,"BH",n=length(sites));st$p_Holm<-p.adjust(st$p_value,"holm",n=length(sites))
    sr_csv(st,file.path(out,"BLINDED_site_comparisons.csv"))
    for(g in sort(unique(b$symptom_group))) {
      z<-b[b$symptom_group==g,];r<-sr_model(z,ae_any~vaccine,"vaccine",g)
      tt<-sr_totals(z,"BLINDED");r$v_ae_any<-tt$ae_any[1];r$v_units<-tt$total_units[1];r$p_ae_any<-tt$ae_any[2];r$p_units<-tt$total_units[2]
      group_tests[[length(group_tests)+1L]]<-r
    }
    gt<-do.call(rbind,group_tests);gt$p_Holm<-p.adjust(gt$p_value,"holm")
    sr_csv(gt,file.path(out,"BLINDED_symptom_group_comparisons.csv"))
  }
  # Exploratory association: differs by arm per +5 pp county positivity?
  # Run only for the selected BLINDED chart window; covariance clusters by county.
  if(window=="BLINDED") {
    z<-merge(d,bgc,by=c("SITEID","week_start"),all.x=TRUE,sort=FALSE)
    z<-z[is.finite(z$positivity)&!is.na(z$region_key),];z$positivity_5pp<-z$positivity/5
    if(nrow(z))add(sr_model(z,ae_any~vaccine*positivity_5pp+factor(SITEID)+factor(week_start)+factor(symptom_group),
      "vaccine:positivity_5pp","BLINDED_arm_by_county_positivity",cluster="region_key",logistic=TRUE))
  }
  mod<-do.call(rbind,models);mod$p_display<-sr_p(mod$p_value);sr_csv(mod,file.path(out,"arm_comparisons.csv"))
  cohort<-as.data.frame(table(roster$SITEID,roster$arm),stringsAsFactors=FALSE);names(cohort)<-c("SITEID","arm","randomized_subjects")
  sr_csv(cohort,file.path(out,"site_randomized_population.csv"))
  params<-data.frame(parameter=c("script_version","chart_window","unit","AE_numerator","small_n","county_only","onset_sensitivity_available","subjects","sites","missing_blind_boundary","background_missing_site_weeks"),
    value=as.character(c("1.0.0",window,"participant x symptom group x Monday-Sunday week","AE_ONLY + BOTH",small_n,TRUE,!is.null(onset),nrow(roster),length(sites),
      if("has_blind_boundary"%in%names(roster))sum(!sr_bool(roster$has_blind_boundary)) else NA,sum(is.na(bg$data$positivity)))))
  sr_csv(params,file.path(out,"site_recording_parameters.csv"))
  source_files<-c(file.path(folder,required),file.path(folder,paste0(windows,"_subject_symptom_weeks.csv")),file.path(folder,"README.txt"),
    file.path(folder,"BLINDED_source_event_week_links.csv"),file.path(folder,"symptom_event_evidence.csv"),weekly_control,bg$file,bg$region_file)
  source_files<-unique(source_files[nzchar(source_files)&file.exists(source_files)])
  sr_csv(data.frame(file=normalizePath(source_files,winslash="/"),md5=unname(tools::md5sum(source_files))),file.path(out,"input_manifest.csv"))
  if(make_pdf) {
    message("Writing one PDF page for each of ",length(sites)," sites")
    sr_pdf(file.path(out,"site_symptom_recording.pdf"))
    tryCatch(for(i in seq_along(sites))sr_page(sites[i],cells,bg$data,roster,window,small_n,i,length(sites)),finally=dev.off())
    if(individual_pdfs) {
      dir.create(file.path(out,"sites"),showWarnings=FALSE)
      for(i in seq_along(sites)) {
        sr_pdf(file.path(out,"sites",paste0(sites[i],"_symptom_recording.pdf")))
        tryCatch(sr_page(sites[i],cells,bg$data,roster,window,small_n,i,length(sites)),finally=dev.off())
      }
    }
  }
  methods<-c("RECORDING-SOURCE AUDIT - v1.0.0", "",
    "Use: run this file in RStudio from the project root. It discovers a complete analyse_weekly_symptoms_ae.R export under analysis/, or accepts a folder/ZIP in PHMPT_WEEKLY_DIR. The chosen export and source hashes are recorded. No raw XPT reconstruction or network downloads are performed.",
    "Denominator: distinct participant x symptom-group x Monday-Sunday week observed in either source. AE share = (AE_ONLY+BOTH)/all units. Counts are conditional on recorded compatible symptoms, not on all randomized subjects. BOTH means same participant/group/week in both sources; it does not establish the same illness or reclassification. No positivity requirement and no exclusion of post-dose days 0-6.",
    "BLINDED and EUA preserve the input subject-specific follow-up. FULL retains original randomized arm after crossover and is descriptive only. Missing blinding boundaries inherit the input fallback; their roster count is reported. Refer to the source README for exact cutoffs and collection rules.",
    "Main comparisons use all weekly units and site-clustered sandwich standard errors (HC0 scores, G/(G-1) correction, t with G-1 degrees of freedom). The raw linear model estimates a percentage-point difference. The adjusted logistic model includes site, calendar week and symptom-group fixed effects and reports an odds ratio. It is a conditional association, not a causal staff-behavior estimate. Rows are not treated as independent chi-square observations.",
    "Sensitivity: count AE records only in their documented onset week, retaining COVID-form assessment weeks; excludes continuing AE interval weeks. Another sensitivity gives every participant with any symptom unit equal weight. A further sensitivity restricts to participants with a documented blinding boundary. These are alternative summaries, not clinical-episode reconstruction.",
    "Site comparisons use participant-clustered raw differences and both BH and Holm adjustment over every randomized site. Symptom-group comparisons use site-clustered differences with Holm adjustment over all groups. Small or invariant cells may be non-estimable; their status is exported. Site tests do not adjust symptom mix or calendar time and are exploratory.",
    "The exploratory arm x county-positivity logistic interaction uses site/week/symptom-group adjustment and county-clustered uncertainty. It reports the multiplicative change in the arm odds ratio per 5 percentage-point increase in positivity. This is an ecological association; plots and interaction do not establish motive, unblinding or a change in PCR testing.",
    "County background: either the prior analysis/background_covid/site_weekly.csv or heatmap site_background_used.csv plus site_region_coverage.csv. State-only series are not passed off as county data; county gaps stay missing. County positivity is the previous asset's published weekly series, not recalculated from trial tests. Site county is not necessarily participant residence.",
    "Charts: one site per landscape PDF page; shared 0-100% scale for the two AE shares and county positivity; no smoothing or interpolation. Open circles have fewer than PHMPT_SMALL_N symptom units. Weeks without symptoms are missing shares, not zero. Numerical denominators and unique participants are in site_arm_weekly_recording.csv. Set PHMPT_INDIVIDUAL_PDFS=TRUE for an additional individual file per site.",
    "Interpretation: vaccination can change symptom composition and infection occurrence; source forms and interval-expanded AEs use different collection/date rules. Solicited diaries are absent from the input. Even an adjusted arm association cannot by itself show discretionary relabeling, missed COVID infections, faulty blinding or efficacy bias.",
    "Method reference: https://sandwich.r-forge.r-project.org/reference/vcovCL.html",
    "Multiplicity reference: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html")
  writeLines(methods,file.path(out,"README.txt"),useBytes=TRUE)
  capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
  print(totals[,c("window","arm","total_units","ae_any","pct_ae_any")],row.names=FALSE)
  print(mod,row.names=FALSE)
  message("Outputs: ",normalizePath(out,winslash="/"))
  invisible(list(output_dir=out,totals=totals,models=mod,cells=cells))
}

# ---- RUN: starts when the whole script is selected and executed ----
# To load functions only: options(phmpt.autorun=FALSE) before source().
if(isTRUE(getOption("phmpt.autorun",TRUE)))analyse_site_symptom_recording()
