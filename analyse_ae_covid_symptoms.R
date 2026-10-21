# Moderna P301 - observed AE terms matching COVID symptoms - v1.0.0
# R 4.3+. Stand-alone: save in the project root, then Select All / Run in RStudio.
# CSV input and PDF/HTML output use base R + grid. XPT/SAS7BDAT need haven.
# This is a transparent analyst dictionary, NOT an official MedDRA SMQ and NOT
# a diagnosis of COVID-19. It maps symptom compatibility, not visit eligibility.

# ---- SETTINGS ----
PHMPT_ROOT <- ""                  # Auto-detect, or D:/Projects/moderna_trial_analysis
PHMPT_RELEASES <- c("EUA", "Prod")
PHMPT_AE_FILES <- list(EUA=character(), Prod=character()) # Optional explicit paths
PHMPT_OUTPUT_DIR <- ""            # Default: analysis/ae_covid_symptoms_<timestamp>
PHMPT_ASSET_DIR <- ""             # Default: assets/
PHMPT_CSV_ENCODING <- "auto"      # UTF-8-BOM, then Windows-1252 if decoding fails
PHMPT_RULES_FILE <- ""            # Optional reviewed replacement for exported rules CSV
PHMPT_OPEN_REPORT <- interactive()
# To load functions only: options(phmpt.autorun=FALSE) before source().

am_version <- "1.0.0"
am_text <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""
  bad <- !validUTF8(x)
  if(any(bad)) x[bad] <- iconv(x[bad],from="windows-1252",to="UTF-8",sub=NA)
  if(anyNA(x)) stop("Undecodable source text; specify PHMPT_CSV_ENCODING.")
  enc2utf8(x)
}
am_norm <- function(x) toupper(trimws(gsub("[[:space:]\\p{Z}]+"," ",am_text(x),perl=TRUE)))
am_col <- function(x,n) if(n %in% names(x)) am_text(x[[n]]) else rep("",nrow(x))
am_path <- function(root,p) if(grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)",p)) p else file.path(root,p)
am_csv <- function(x,p) {
  dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE)
  con <- file(p,"wb"); on.exit(close(con)); writeBin(as.raw(c(239,187,191)),con)
  write.table(x,con,sep=",",row.names=FALSE,col.names=TRUE,quote=TRUE,na="",qmethod="double")
}
am_root <- function(root="") {
  if(nzchar(root)) return(normalizePath(root,winslash="/",mustWork=TRUE))
  cand <- getwd()
  for(f in sys.frames()) if(!is.null(f$ofile)) cand <- c(cand,dirname(f$ofile))
  a <- grep("^--file=",commandArgs(),value=TRUE)
  if(length(a)) cand <- c(cand,dirname(sub("^--file=","",a[1])))
  if(interactive() && requireNamespace("rstudioapi",quietly=TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path,error=function(e)"")
    if(nzchar(p)) cand <- c(cand,dirname(p))
  }
  for(p in unique(cand[dir.exists(cand)])) for(i in seq_len(8)) {
    if(any(dir.exists(file.path(p,c("eua_data","prod_data")))))
      return(normalizePath(p,winslash="/",mustWork=TRUE))
    parent <- dirname(p); if(parent==p) break; p <- parent
  }
  stop("Cannot find eua_data/prod_data. Set PHMPT_ROOT to the project root.")
}
am_symptoms <- function() {
  data.frame(symptom_code=c("FEVER","CHILLS","COUGH","SHBREATH","DIBREATH","FATIGUE",
    "MYALGIA","BODYACHE","HEADACHE","NLSMELL","NLTASTE","STHROAT","NASALCON","RUNNOSE","NAUSEA","VOMIT","DIARRHEA"),
    symptom_label=c("Fever","Chills","Cough","Shortness of breath","Difficulty breathing","Fatigue",
    "Muscle aches (myalgia)","Body aches","Headache","New loss of smell","New loss of taste","Sore throat",
    "Nasal congestion","Runny nose (rhinorrhea)","Nausea","Vomiting","Diarrhea"),
    symptom_group=c("FEVER","CHILLS","COUGH","DYSPNEA","DYSPNEA","FATIGUE","MYALGIA","MYALGIA",
    "HEADACHE","SMELL_TASTE","SMELL_TASTE","SORE_THROAT","NASAL","NASAL","NAUSEA_VOMITING","NAUSEA_VOMITING","DIARRHEA"),
    min_duration_hours=c(0,0,0,0,0,rep(48,12)),stringsAsFactors=FALSE)
}
am_rules <- function() {
  z <- list()
  add <- function(code,terms,status="DIRECT",note="Explicit symptom or lexical synonym; cause is not inferred.") {
    term <- strsplit(terms,"|",fixed=TRUE)[[1]]
    z[[length(z)+1L]] <<- data.frame(normalized_term=am_norm(term),symptom_code=code,
      match_class=status,include_default=status=="DIRECT",rationale=note,stringsAsFactors=FALSE)
  }
  add("FEVER","Pyrexia|Fever|Post procedural fever|Postoperative fever|Intermittent fever|Fever of unknown origin|Febrile reaction|Fever chills")
  add("FEVER","Body temperature increased|High temperature|Low grade fever|Feverish|Slight temperature","RELATED_REVIEW","Raised temperature or feverishness may be below 38 C; verify measured temperature.")
  add("CHILLS","Chills|Chill|Rigors|Fever chills")
  add("COUGH","Cough|Productive cough|Dry cough|Allergic cough|Cough nonproductive|Coughing|Cough aggravated|Cough increased|Chronic cough|Nocturnal cough|Persistent cough|Drug-induced cough|Smoker's cough")
  add("SHBREATH","Dyspnoea|Dyspnea|Dyspnoea exertional|Dyspnea exertional|Exertional dyspnea|Dyspnoea at rest|Shortness of breath|Short of breath|Increased shortness of breath|Dyspnea exacerbated")
  add("DIBREATH","Difficulty breathing|Breathing difficult|Difficult breathing|Labored breathing|Laboured breathing")
  add("DIBREATH","Respiratory distress|Orthopnoea|Orthopnea","RELATED_REVIEW","Respiratory finding requiring review; do not infer the checklist symptom automatically.")
  add("FATIGUE","Fatigue|Tiredness|Fatigue aggravated|Fatigue extreme")
  add("FATIGUE","Asthenia|Weakness|Weakness generalized|Generalised weakness|Generalized weakness|Muscle fatigue|Mental fatigue|Malaise","RELATED_REVIEW","Weakness, malaise or localized/mental fatigue is broader than the checklist symptom.")
  add("MYALGIA","Myalgia|Muscle ache|Muscle aches|Muscle pain|Muscle soreness|Muscular pain|Muscular pains|Pain muscle|Myalgia aggravated|Generalised muscle aches|Generalized muscle aches|Localised muscle pain|Localized muscle pain|Myalgia of lower extremities|Myalgia upper extremities|Polymyalgia|Tenderness muscle")
  add("BODYACHE","Body ache|Body aches|General body pain|General body aches|Generalized aching|Generalised aching|All over body aches|(Whole body) body aches")
  add("BODYACHE","Generalised pain|Generalized pain","RELATED_REVIEW","General pain is not necessarily the protocol's body-ache symptom; inspect verbatim description.")
  add("HEADACHE","Headache|Head pain|Pain head|Tension headache|Tension headaches|Headache tension|Sinus headache|Headache sinus|Migraine headache|Cluster headache|Cluster headaches|Exertional headache|New daily persistent headache|Post-traumatic headache|Headache post-traumatic|Procedural headache|Headache postoperative|Post procedural headache|Primary cough headache|Primary headache associated with sexual activity|Coital headache|Drug-induced headache|Forehead headache|Frequent headaches|Frontal headache|Headache (excl migraine) aggravated|Headache aggravated|Headache dull|Headache recurrent|Headache temporal|Headache unilateral|Intermittent headache|Occipital headache|Temporal headache")
  add("HEADACHE","Migraine|Migraine aggravated|Complicated migraine|Migraine with aura|Migraine without aura|Migraine aura|Ophthalmic migraine|Vestibular migraine","RELATED_REVIEW","A migraine diagnosis, particularly aura, does not establish headache on a particular day.")
  add("NLSMELL","Anosmia|Loss of smell|New loss of smell|Loss of taste and smell|Loss of smell and taste")
  add("NLSMELL","Hyposmia|Diminished sense of smell|Olfactory disorder|Parosmia|Altered smell sensation|Smell alteration|Smell change|Phantosmia","RELATED_REVIEW","Reduced or distorted smell is not synonymous with new loss; review onset and symptom detail.")
  add("NLTASTE","Ageusia|Loss of taste|New loss of taste|Loss of taste and smell|Loss of smell and taste")
  add("NLTASTE","Hypogeusia|Taste diminished|Dysgeusia|Taste disorder|Impaired taste|Taste alteration|Taste altered|Taste changed|Taste bitter|Taste metallic|Taste salty","RELATED_REVIEW","Reduced or altered taste is not synonymous with new loss; review onset and symptom detail.")
  add("STHROAT","Oropharyngeal pain|Sore throat|Pain throat|Throat pain|Throat sore|Chronic sore throat|Streptococcal sore throat")
  add("STHROAT","Throat irritation|Itchy throat|Pharynx irritated sensation of|Pharyngeal discomfort|Oropharyngeal discomfort|Throat discomfort|Dry throat|Throat tightness|Pharyngeal paraesthesia|Tingling throat","RELATED_REVIEW","Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.")
  add("NASALCON","Nasal congestion|Congestion nasal|Nasal stuffiness|Nose congestion|Chronic nasal congestion|Nasal sinus congestion")
  add("NASALCON","Sinus congestion|Upper respiratory tract congestion|Respiratory tract congestion","RELATED_REVIEW","Congestion site or equivalence to the nasal checklist item needs review.")
  add("RUNNOSE","Rhinorrhoea|Rhinorrhea|Runny nose|Nasal discharge|Nasal mucus increased")
  add("NAUSEA","Nausea|Nausea aggravated|Procedural nausea|Post procedural nausea|Postoperative nausea|Nausea and vomiting|Nausea/vomiting")
  add("VOMIT","Vomiting|Vomited|Emesis|Hyperemesis|Procedural vomiting|Nausea and vomiting|Nausea/vomiting")
  add("DIARRHEA","Diarrhoea|Diarrhea|Acute diarrhea|Diarrhea aggravated|Loose bowel|Loose bowels|Loose stools|Stools watery|Watery diarrhea|Bacterial diarrhoea|Bacterial diarrhea|Viral diarrhoea|Viral diarrhea|Diarrhoea infectious|Traveller's diarrhea")
  out <- unique(do.call(rbind,z)); rownames(out) <- NULL; out
}
am_validate_rules <- function(x) {
  need <- c("normalized_term","symptom_code","match_class","include_default","rationale")
  if(!all(need %in% names(x))) stop("Rules require: ",paste(need,collapse=", "))
  x <- x[need]; x$normalized_term <- am_norm(x$normalized_term)
  x$include_default <- am_norm(x$include_default)
  if(any(!x$include_default %in% c("TRUE","FALSE"))) stop("Rules include_default must be TRUE/FALSE.")
  x$include_default <- x$include_default=="TRUE"
  if(any(!x$symptom_code %in% am_symptoms()$symptom_code) || any(!nzchar(x$normalized_term))) stop("Invalid rule symptom or blank term.")
  if(any(!x$match_class %in% c("DIRECT","RELATED_REVIEW"))) stop("Invalid rule match_class.")
  if(any(x$include_default & x$match_class!="DIRECT")) stop("Review terms cannot be included by default; review and reclassify explicitly.")
  if(anyDuplicated(paste(x$normalized_term,x$symptom_code,sep="\r"))) stop("Duplicate term/symptom rules.")
  x
}
am_read <- function(p,encoding="auto") {
  ext <- tolower(tools::file_ext(p))
  if(ext=="csv") {
    errors <- character(); x <- NULL
    for(enc in if(encoding=="auto") c("UTF-8-BOM","windows-1252") else encoding) {
      y <- tryCatch(withCallingHandlers(read.csv(p,colClasses="character",na.strings=character(),
        check.names=FALSE,fileEncoding=enc),warning=function(w)stop(conditionMessage(w))),error=function(e)e)
      if(!inherits(y,"error")) {x<-y; break}; errors <- c(errors,conditionMessage(y))
    }
    if(is.null(x)) stop("CSV read failed: ",paste(errors,collapse="; "))
  } else if(ext %in% c("xpt","sas7bdat")) {
    if(!requireNamespace("haven",quietly=TRUE)) stop('XPT/SAS7BDAT requires install.packages("haven").')
    x <- as.data.frame(if(ext=="xpt") haven::read_xpt(p,.name_repair="minimal") else haven::read_sas(p,.name_repair="minimal"))
  } else stop("Unsupported input extension: ",ext)
  names(x) <- toupper(sub("^\ufeff","",names(x)))
  if(anyDuplicated(names(x))) stop("Duplicate input column names.")
  x[] <- lapply(x,am_text); x
}
am_candidates <- function(root,release) {
  p <- file.path(root,if(release=="EUA")"eua_data" else "prod_data")
  x <- list.files(p,recursive=TRUE,full.names=TRUE,pattern="\\.(csv|xpt|sas7bdat)$",ignore.case=TRUE)
  x <- x[grepl("(^|[_ -])ae([._ -]|$)",basename(x),ignore.case=TRUE)]
  x <- x[!grepl("p101|p201|dmid",basename(x),ignore.case=TRUE)]
  x[order(match(tolower(tools::file_ext(x)),c("csv","xpt","sas7bdat")),x)]
}
am_prepare <- function(x,snapshot=FALSE) {
  if(!all(c("AETERM","AEDECOD") %in% names(x))) stop("Not a raw AE table: AETERM and AEDECOD are required.")
  n0 <- nrow(x)
  if("STUDYID" %in% names(x)) x <- x[am_norm(x$STUDYID)=="MRNA-1273-P301",,drop=FALSE]
  else if(!snapshot) stop("Missing STUDYID; refusing to guess the study.")
  if(!nrow(x)) stop("No mRNA-1273-P301 AE records.")
  idfield <- intersect(c("USUBJID","SUBJID",if(snapshot)"SUBJECT"),names(x))
  if(!length(idfield)) stop("Missing subject identifier.")
  fields <- c("STUDYID","USUBJID","AESEQ","AETERM","AEDECOD","AELLT","AEPTCD","AELLTCD","AESTDTC","AEENDTC","AESER","AEREL")
  subject <- x[[idfield[1]]]
  seqid <- am_col(x,"AESEQ")
  if(snapshot && "SOURCE_REF" %in% names(x)) seqid <- x$SOURCE_REF
  missingseq <- !nzchar(seqid); seqid[missingseq] <- paste0("SOURCE_ROW_",which(missingseq))
  x <- as.data.frame(setNames(lapply(fields,function(n)am_col(x,n)),fields),stringsAsFactors=FALSE)
  x$USUBJID <- subject; x$AESEQ <- seqid
  if(snapshot) x$STUDYID <- "mRNA-1273-P301"
  if(any(!nzchar(trimws(subject)))) stop("Blank subject identifier.")
  x <- unique(x) # Exact duplicated rows only; distinct AE records are preserved.
  if(anyDuplicated(paste(x$USUBJID,x$AESEQ,sep="\r"))) stop("Conflicting rows share USUBJID/AESEQ; reconcile the source first.")
  attr(x,"dropped_other_study_rows") <- n0-nrow(x); x
}
am_signature <- function(x) {
  y <- x[do.call(order,x[c("USUBJID","AESEQ")]),,drop=FALSE]; rownames(y)<-NULL
  attributes(y) <- attributes(y)[c("names","row.names","class")]
  p <- tempfile(); on.exit(unlink(p)); saveRDS(y,p,compress=FALSE,version=2)
  unname(tools::md5sum(p))
}
am_qualify <- function(code) ifelse(code=="FEVER","Verify measured temperature >=38 C; AE label alone is insufficient.",
  ifelse(code %in% c("NLSMELL","NLTASTE"),"Verify new onset and duration >=48 h.",
  ifelse(code %in% c("CHILLS","COUGH","SHBREATH","DIBREATH"),"Any-duration symptom; review onset and post-vaccination context.",
    "Verify duration >=48 h and post-vaccination context.")))
am_screen <- function(x) grepl(paste(c("\\b(FEVER|PYREXIA|CHILLS?|RIGORS?|COUGH|COUGHING|DYSPNOEA|DYSPNEA|BREATHING|BREATH|FATIGUE|TIREDNESS|ASTHENIA|MYALGIA|HEADACHES?|MIGRAINE|SMELL|TASTE|ANOSMIA|AGEUSIA|HYPOSMIA|HYPOGEUSIA|DYSGEUSIA|PAROSMIA|THROAT|RHINORRHOEA|RHINORRHEA|CONGESTION|NAUSEA|VOMITING|EMESIS|DIARRHOEA|DIARRHEA)\\b","BODY ACHE","MUSCLE ACHE"),collapse="|"),x,perl=TRUE)
am_diagnosis <- function(x) grepl("COVID|CORONAVIRUS|PNEUMONIA|INFLUENZA|RESPIRATORY.*INFECTION|BRONCHITIS|RHINITIS|SINUSITIS|RESPIRATORY FAILURE|ARDS",x,perl=TRUE)
am_escape <- function(x) {
  x<-gsub("&","&amp;",x,fixed=TRUE); x<-gsub("<","&lt;",x,fixed=TRUE); gsub(">","&gt;",x,fixed=TRUE)
}
am_report <- function(tab,pdf,html,subtitle,notes) {
  heading <- c("COVID symptom","Direct preferred terms (AEDECOD)","Other direct terms (AELLT)","Related terms: review only")
  # A4 landscape, vector PDF; no Chrome, LaTeX or external report tools required.
  if(capabilities("cairo")) grDevices::cairo_pdf(pdf,width=297/25.4,height=210/25.4,family="sans",onefile=TRUE)
  else grDevices::pdf(pdf,width=297/25.4,height=210/25.4,paper="special",family="Helvetica",onefile=TRUE)
  on.exit(grDevices::dev.off(),add=TRUE)
  W<-297/25.4*72; H<-210/25.4*72; margin<-28; cw<-c(111,235,260,W-2*margin-606)
  xpos<-margin+c(0,cumsum(cw)[-4]); fontsize<-8.6; lineheight<-11.3
  wrap <- function(s,width,fs=fontsize,bold=FALSE) {
    if(!nzchar(s)) return("-")
    words<-strsplit(s," +")[[1]]; lines<-character(); current<-""
    for(w in words) {
      candidate<-paste(c(if(nzchar(current))current,w),collapse=" ")
      size<-grid::convertWidth(grid::grobWidth(grid::textGrob(candidate,
        gp=grid::gpar(fontsize=fs,fontface=if(bold)"bold" else "plain"))),"points",valueOnly=TRUE)
      if(size>width && nzchar(current)){lines<-c(lines,current);current<-w}else current<-candidate
    }
    c(lines,current)
  }
  textat <- function(label,x,y,size=fontsize,bold=FALSE,color="#182632")
    grid::grid.text(label,x=grid::unit(x,"points"),y=grid::unit(y,"points"),just=c("left","top"),
      gp=grid::gpar(fontsize=size,fontface=if(bold)"bold" else "plain",col=color,lineheight=1.1))
  page<-0L
  newpage <- function() {
    page<<-page+1L; grid::grid.newpage(); grid::pushViewport(grid::viewport(gp=grid::gpar(fontsize=fontsize)))
    textat("MODERNA P301  |  AE TO COVID SYMPTOM MAPPING",margin,H-23,16,TRUE)
    textat(subtitle,margin,H-46,8.5)
    textat("Observed terms only. Direct = symptom compatibility; related = excluded from default matches.",margin,H-61,8.5)
    y<-H-80
    grid::grid.rect(x=grid::unit(margin,"points"),y=grid::unit(y,"points"),width=grid::unit(sum(cw),"points"),height=grid::unit(31,"points"),just=c("left","top"),gp=grid::gpar(fill="#18394A",col=NA))
    for(j in 1:4) textat(paste(wrap(heading[j],cw[j]-12,bold=TRUE),collapse="\n"),xpos[j]+6,y-6,8.6,TRUE,"white")
    textat(paste0("v",am_version,"  |  Full field-specific mapping and verbatim terms: ae_covid_symptom_mapping.csv"),margin,29,7.5)
    textat(paste("Page",page),W-70,29,7.5)
    y-31
  }
  y<-newpage()
  for(i in seq_len(nrow(tab))) {
    cells<-lapply(1:4,function(j)wrap(as.character(tab[i,j]),cw[j]-12,bold=j==1))
    height<-max(lengths(cells))*lineheight+8
    if(height>H-165) stop("A table row exceeds one PDF page. Shorten custom rule labels or use the HTML report.")
    if(y-height<55)y<-newpage()
    grid::grid.rect(x=grid::unit(margin,"points"),y=grid::unit(y,"points"),width=grid::unit(sum(cw),"points"),height=grid::unit(height,"points"),just=c("left","top"),gp=grid::gpar(fill=if(i%%2)"#F0F4F6" else "white",col=NA))
    for(j in 1:4) for(k in seq_along(cells[[j]])) textat(cells[[j]][k],xpos[j]+6,y-4-(k-1)*lineheight,bold=j==1)
    grid::grid.lines(x=grid::unit(c(margin,W-margin),"points"),y=grid::unit(rep(y-height,2),"points"),gp=grid::gpar(col="#CCD5DA",lwd=.5))
    y<-y-height
  }
  # Methods page makes the PDF usable independently of the CSV or this script.
  notespage <- function() {
    page<<-page+1L; grid::grid.newpage(); grid::pushViewport(grid::viewport(gp=grid::gpar(fontsize=10)))
    textat("INTERPRETING AND REUSING THE MAPPING",margin,H-28,16,TRUE)
    textat(paste("Page",page),W-70,29,7.5)
    H-62
  }
  y<-notespage()
  for(s in notes) {
    lines<-wrap(s,W-2*margin,10)
    if(y-length(lines)*14<55)y<-notespage()
    for(line in lines){textat(line,margin,y,10);y<-y-14}
    y<-y-9
  }
  rows<-vapply(seq_len(nrow(tab)),function(i)paste0("<tr>",paste0("<td>",am_escape(as.character(tab[i,])),"</td>",collapse=""),"</tr>"),character(1))
  doc<-paste0('<!doctype html><html lang="en"><meta charset="utf-8"><title>AE to COVID symptom mapping</title><style>',
    '@page{size:A4 landscape;margin:12mm}body{font-family:Arial,sans-serif;color:#182632;margin:25px}h1{font-size:22px}p,li{font-size:12px;line-height:1.5}table{width:100%;border-collapse:collapse;font-size:11px;table-layout:fixed}th{background:#18394a;color:white;text-align:left}td,th{padding:8px;border-bottom:1px solid #cbd4da;vertical-align:top;overflow-wrap:anywhere}tr:nth-child(odd){background:#f0f4f6}thead{display:table-header-group}tr{break-inside:avoid}th:first-child{width:14%}th:nth-child(2){width:30%}th:nth-child(3){width:33%}th:nth-child(4){width:23%}.notes{break-before:page}button{padding:8px}@media print{button{display:none}body{margin:0}}',
    '</style><h1>Moderna P301: AE to COVID symptom mapping</h1><p>',am_escape(subtitle),'</p><p>Observed coded terms only. The CSV also retains exact verbatim matches.</p><button onclick="window.print()">Print / Save as PDF</button><table><thead><tr>',
    paste0('<th>',heading,'</th>',collapse=''),'</tr></thead><tbody>',paste(rows,collapse=''),'</tbody></table><div class="notes"><h2>Methods and interpretation</h2>',
    paste0('<p>',am_escape(notes),'</p>',collapse=''),'</div></html>')
  writeLines(enc2utf8(doc),html,useBytes=TRUE)
}

analyse_ae_covid_symptoms <- function(root="",releases=c("EUA","Prod"),ae_files=list(),output_dir="",asset_dir="",
  encoding="auto",rules_file="",open_report=FALSE,snapshot=FALSE) {
  root<-am_root(root)
  if(any(!releases %in% c("EUA","Prod")))stop("Release must be EUA or Prod.")
  out<-if(nzchar(output_dir))am_path(root,output_dir) else file.path(root,"analysis",paste0("ae_covid_symptoms_",format(Sys.time(),"%Y%m%d_%H%M%S")))
  assetdir<-if(nzchar(asset_dir))am_path(root,asset_dir) else file.path(root,"assets")
  dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(assetdir,recursive=TRUE,showWarnings=FALSE)
  rules<-if(nzchar(rules_file))am_read(am_path(root,rules_file),encoding) else am_rules()
  names(rules)<-tolower(names(rules));rules<-am_validate_rules(rules)
  sy<-am_symptoms(); inventory<-list(); observations<-list(); evidence<-list(); seen<-character(); source_no<-0L
  for(rel in releases) {
    paths<-ae_files[[rel]]
    if(is.null(paths)||!length(paths))paths<-am_candidates(root,rel) else paths<-vapply(paths,function(p)am_path(root,p),character(1))
    if(!length(paths)) {
      inventory[[length(inventory)+1L]]<-data.frame(source_id="",release=rel,path="",md5="",n_records=0L,status="NO_FILES",detail="No candidate AE files found.")
      next
    }
    for(p in unique(paths)) {
      source_no<-source_no+1L; sid<-sprintf("AE%04d",source_no); message("[",rel,"] ",basename(p))
      ans<-tryCatch(am_prepare(am_read(p,encoding),snapshot=snapshot),error=function(e)e)
      good<-!inherits(ans,"error");status<-if(good)"USED" else "REJECTED"; detail<-if(good)"" else conditionMessage(ans)
      if(good) {
        sig<-paste(rel,am_signature(ans),sep=":")
        if(sig %in% seen){status<-"EQUIVALENT_COPY";detail<-"Same projected AE records as an already used source in this release."}else seen<-c(seen,sig)
      }
      inventory[[length(inventory)+1L]]<-data.frame(source_id=sid,release=rel,path=normalizePath(p,winslash="/",mustWork=FALSE),md5=if(file.exists(p))unname(tools::md5sum(p)) else "",n_records=if(good)nrow(ans) else NA_integer_,status=status,detail=detail)
      if(status!="USED")next
      x<-ans
      for(field in c("AEDECOD","AELLT","AETERM")) {
        term<-x[[field]]; norm<-am_norm(term); use<-nzchar(norm)
        if(!any(use))next
        ix<-data.frame(row_id=which(use),match_field=field,ae_term=term[use],normalized_term=norm[use],stringsAsFactors=FALSE)
        if(!nrow(ix))next
        counts<-aggregate(rep(1L,nrow(ix)),ix[c("match_field","ae_term","normalized_term")],sum);names(counts)[4]<-"n_records_in_source"
        counts$source_id<-sid;counts$release<-rel; observations[[length(observations)+1L]]<-counts
        hit<-merge(ix,rules,by="normalized_term",sort=FALSE)
        if(!nrow(hit))next
        hit$source_id<-sid;hit$release<-rel
        for(n in c("USUBJID","AESEQ","AETERM","AEDECOD","AELLT","AEPTCD","AELLTCD","AESTDTC","AEENDTC","AESER","AEREL"))hit[[n]]<-x[[n]][hit$row_id]
        evidence[[length(evidence)+1L]]<-hit
      }
    }
  }
  inv<-do.call(rbind,inventory);am_csv(inv,file.path(out,"ae_input_manifest.csv"))
  if(!length(observations))stop("No usable AE input. See ",file.path(out,"ae_input_manifest.csv"))
  obs<-unique(do.call(rbind,observations));rownames(obs)<-NULL
  # Block silent partial coverage: rejected candidates must be resolved or explicitly excluded.
  if(any(inv$status=="REJECTED"))stop("Some AE candidates could not be read. See ae_input_manifest.csv. Install haven or set PHMPT_AE_FILES to the intended readable inputs, then rerun.")
  terms<-unique(obs[c("match_field","ae_term","normalized_term")])
  mapping<-merge(terms,rules,by="normalized_term",sort=FALSE)
  mapping<-merge(mapping,sy,by="symptom_code",sort=FALSE)
  mapping$qualification_check<-am_qualify(mapping$symptom_code)
  mapping$dictionary_version<-rep(am_version,nrow(mapping))
  mapping<-mapping[order(match(mapping$symptom_code,sy$symptom_code),mapping$match_class,mapping$match_field,mapping$ae_term),]
  mapping<-mapping[c("dictionary_version","match_field","ae_term","normalized_term","symptom_code","symptom_label","symptom_group","match_class","include_default","min_duration_hours","qualification_check","rationale")]
  rownames(mapping)<-NULL
  am_csv(mapping,file.path(out,"ae_covid_symptom_mapping.csv"))
  # Do not silently overwrite an edited asset: keep an immutable copy from each prior run.
  target<-file.path(assetdir,"ae_covid_symptom_mapping.csv")
  if(file.exists(target)) {
    backup<-file.path(assetdir,paste0("ae_covid_symptom_mapping_previous_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid(),".csv"))
    if(!file.copy(target,backup,overwrite=FALSE))stop("Could not back up the existing asset.")
  }
  am_csv(mapping,target);am_csv(rules,file.path(out,"ae_covid_symptom_rules.csv"))
  am_csv(obs,file.path(out,"ae_observed_terms_by_source.csv"))
  if(length(evidence)) {
    ev<-do.call(rbind,evidence);ev<-merge(ev,sy,by="symptom_code",sort=FALSE)
    # Keep each matching field for traceability. This is evidence, NOT event counts.
    am_csv(ev,file.path(out,"ae_covid_match_evidence.csv"))
    # Specific verbatim/LLT mappings take precedence over a more generic PT mapping
    # for the same symptom parameter. Multiple symptom parameters can share one group.
    ev$priority<-match(ev$match_field,c("AETERM","AELLT","AEDECOD"))
    ev<-ev[order(ev$source_id,ev$USUBJID,ev$AESEQ,ev$symptom_code,ev$priority),]
    key<-paste(ev$source_id,ev$USUBJID,ev$AESEQ,ev$symptom_code,sep="\r")
    rec<-ev[!duplicated(key),];rec$qualification_check<-am_qualify(rec$symptom_code)
    am_csv(rec,file.path(out,"ae_covid_record_matches.csv"))
  }
  # No broad regex is promoted to an accepted symptom mapping.
  unmatched<-terms[!terms$normalized_term %in% rules$normalized_term,,drop=FALSE]
  unmatched$lexical_candidate<-am_screen(unmatched$normalized_term)
  unmatched$diagnosis_or_severity_candidate<-am_diagnosis(unmatched$normalized_term)
  unmatched$review_reason<-ifelse(unmatched$diagnosis_or_severity_candidate,"Diagnosis/severity wording: inspect source; no symptom inferred.",ifelse(unmatched$lexical_candidate,"Unmapped lexical candidate; review negation, context and exact wording.","No dictionary match."))
  am_csv(unmatched,file.path(out,"ae_unmapped_terms.csv"))
  am_csv(unmatched[unmatched$lexical_candidate|unmatched$diagnosis_or_severity_candidate,,drop=FALSE],file.path(out,"ae_terms_for_review.csv"))
  jointerms<-function(v)if(length(v))paste(sort(unique(v)),collapse="; ") else "-"
  tab<-do.call(rbind,lapply(seq_len(nrow(sy)),function(i) {
    m<-mapping[mapping$symptom_code==sy$symptom_code[i],];direct<-m[m$include_default,]
    pt<-direct$ae_term[direct$match_field=="AEDECOD"]
    llt<-direct$ae_term[direct$match_field=="AELLT" & !am_norm(direct$ae_term) %in% am_norm(pt)]
    review<-m$ae_term[!m$include_default & m$match_field %in% c("AEDECOD","AELLT")]
    data.frame(symptom=paste0(sy$symptom_label[i]," [",sy$symptom_code[i],"]"),direct_pt=jointerms(pt),direct_llt=jointerms(llt),related_review=jointerms(review),stringsAsFactors=FALSE)
  }))
  used<-inv[inv$status=="USED",];subtitle<-paste0("Source scope: ",paste(unique(used$release),collapse=" + ")," | ",nrow(used)," distinct AE table(s) | ",if(snapshot)"provided production AE export" else "local PHMPT inputs"," | Generated ",format(Sys.Date()))
  notes<-c(
    "The table shows observed coded terms, not a list of proven COVID-19 cases. DIRECT means an explicit symptom or lexical synonym in the selected AE field. RELATED_REVIEW is a broader or ambiguous association and is excluded from default matches. Both the rules and judgments are analyst-authored; this is not an official MedDRA SMQ or sponsor dictionary.",
    "Input: raw P301 AE tables, discovered separately under eua_data and prod_data, or explicit paths. All AE dates and investigator causal attributions are retained. No treatment arm, positivity, efficacy membership, unblinding or post-dose filter is applied. This is a vocabulary inventory, not an incidence or treatment-effect analysis.",
    "Preferred terms (AEDECOD), lower-level terms (AELLT) and verbatim reports (AETERM) are searched by exact equality after case and whitespace normalization. Punctuation is preserved. The printable table shows coded terms; the CSV also contains observed exact verbatim matches. Unmatched terms with suggestive words go to a separate review queue and are never accepted by regex alone.",
    "FEVER: verify measured temperature >=38 C; pyrexia or low-grade fever alone cannot establish this threshold. NLSMELL/NLTASTE: verify new onset. Other qualifying symptoms generally require >=48 h, except chills, cough and breathing symptoms. Calendar dates alone cannot establish 48 elapsed hours. Investigator discretion applied during the first seven days after vaccination, with respiratory symptoms still requiring evaluation.",
    "Generic diagnoses (COVID-19, influenza, pneumonia, bronchitis, rhinitis) are not expanded into unreported symptoms. Pneumonia, respiratory failure and ICU admission belong in a separate severity/endpoint review. Chronic or causally qualified symptoms are retained for burden analysis; their timing and relevance to an illness visit require separate review.",
    "For reuse, join the CSV by match_field plus normalized_term. include_default=TRUE selects direct compatibility only, not protocol eligibility. Use source_id + USUBJID + AESEQ + symptom_code to identify record matches; inspect the evidence file when fields differ. AETERM, then AELLT, then AEDECOD determines the selected evidence for the same parameter. Check review-only evidence before widening the dictionary.",
    "The 17 parameters collapse to 12 symptom groups: breathing, muscle/body aches, smell/taste, nasal symptoms and nausea/vomiting each share a group. Count a participant/day/group once after reconciling AE and COVID-form evidence. One dyspnoea record must not count as two distinct systemic/respiratory symptoms. The 12 groups are an analysis convention, not a replacement endpoint definition.",
    "Multiple releases and production submissions may repeat or revise the same events. Equivalent projected tables are skipped within each release; differing snapshots are inventoried separately. Never add their event counts. Choose one reconciled snapshot per analysis window. Missing releases and source paths/hashes are documented in ae_input_manifest.csv.",
    "References: P301 protocol amendment 4 (30 Sep 2020), sections 8.1.1-8.1.3; production ADaM define.xml, ADSYMP parameter definitions; SDTM AE variables. MedDRA hierarchy and Data Retrieval and Presentation: Points to Consider explain why coded diagnoses and explicit symptoms must be distinguished. https://www.meddra.org/how-to-use/basics/hierarchy",
    paste0("Dictionary version ",am_version,". Actual input tables and file hashes are recorded in the accompanying manifest. R session information and this run's rules are saved with the report.")
  )
  am_csv(tab,file.path(out,"ae_covid_print_table.csv"))
  am_report(tab,file.path(out,"ae_covid_symptom_table.pdf"),file.path(out,"ae_covid_symptom_table.html"),subtitle,notes)
  writeLines(notes,file.path(out,"README.txt"));capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
  message("Report: ",file.path(out,"ae_covid_symptom_table.pdf"));message("Reusable asset: ",target)
  if(open_report)utils::browseURL(file.path(out,"ae_covid_symptom_table.html"))
  invisible(list(output_dir=out,asset=target,mapping=mapping,manifest=inv,table=tab))
}

# ---- RUN ----
if(isTRUE(getOption("phmpt.autorun",TRUE))) {
  analyse_ae_covid_symptoms(root=PHMPT_ROOT,releases=PHMPT_RELEASES,ae_files=PHMPT_AE_FILES,
    output_dir=PHMPT_OUTPUT_DIR,asset_dir=PHMPT_ASSET_DIR,encoding=PHMPT_CSV_ENCODING,
    rules_file=PHMPT_RULES_FILE,open_report=PHMPT_OPEN_REPORT)
}
