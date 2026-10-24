# Moderna P301: EUA death reporting and CRF review | v1.1.0 | 2026-10-24
# Standalone R >=4.3. Run the entire file from the RStudio project root.
# Reads local PHMPT releases only; does not download trial documents.
# XPT/SAS7BDAT: haven. Optional PDF text: pdftools (or pdftotext). ZIP: zip.
# Install once if needed: install.packages(c("haven", "pdftools", "zip"))
# Load without running: options(phmpt.autorun=FALSE) before source().
# The 24-hour clock starts at investigator awareness, not AE onset or death.
# A fatal AE may have been reported promptly before a later death update.

PHMPT_ROOT <- getOption("phmpt.root", "")
PHMPT_OUTPUT_DIR <- getOption("phmpt.deaths_output", "")
PHMPT_SUBMISSIONS <- c(EUA="27073_S5", BLA="125752_S3")
PHMPT_FILES <- list(EUA=list(), BLA=list()) # e.g. PHMPT_FILES$EUA$ADSL <- "eua_data/...adsl.jmp"
PHMPT_EUA_CUTOFF <- as.Date("2020-11-25")
PHMPT_EUA_UPDATE <- as.Date("2020-12-03")
PHMPT_REVIEW_END <- as.Date("2020-12-18")  # authorization date; not the safety-data cutoff
PHMPT_BLA_CUTOFF <- as.Date("2021-03-26")
PHMPT_REVIEW_FILE <- getOption("phmpt.deaths_review", "assets/death_reporting_review.csv")
PHMPT_CSV_DATE_ORIGIN <- ""  # set only for documented numeric CSV dates, e.g. "1960-01-01"
PHMPT_CSV_ENCODING <- "auto"
PHMPT_PDF_ROOTS <- c("eua_data", "prod_data")
PHMPT_PDF_TEXT_SCAN <- TRUE # scans copied subject documents, not the entire PDF collection
PHMPT_COPY_PDFS <- TRUE
PHMPT_MAKE_ZIP <- TRUE
PHMPT_INCLUDE_TEXT_ONLY <- TRUE # separate review queue; never counted as confirmed death evidence
PHMPT_PDF_LINKS_FILE <- "" # optional CSV with subject,path for PDFs whose names lack an ID
PHMPT_EXTRA_SUBJECTS <- character() # optional additional P301 subject IDs to export
PHMPT_SHOW_REPORT <- interactive()
# V36: print-ready Substack table, plus CSV inputs and provenance.
PHMPT_AUDIT_ASSET_DIR <- getOption("phmpt.audit_asset_dir", "audit_assets")
PHMPT_AUDIT_DPI <- getOption("phmpt.audit_dpi", 300L)

DR_VERSION <- "1.1.0"
dr_domains <- c("ADSL","DM","DS","AE","ADAE","DV","SUPPDM","SUPPDS","SUPPAE")
dr_text <- function(x) {
  missing <- is.na(x)
  if(inherits(x,"Date"))x<-format(x,"%Y-%m-%d")
  else if(inherits(x,"POSIXt"))x<-format(x,"%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  else if(is.numeric(x))x<-format(x,scientific=FALSE,trim=TRUE,digits=15)
  x<-as.character(x);x[is.na(x)|missing]<-""
  bad<-!validUTF8(x);if(any(bad))x[bad]<-iconv(x[bad],from="windows-1252",to="UTF-8",sub="byte")
  trimws(x)
}
dr_norm <- function(x)toupper(gsub("[[:space:]]+"," ",dr_text(x)))
dr_col <- function(x,n)if(n%in%names(x))dr_text(x[[n]]) else rep("",nrow(x))
dr_bind <- function(xs) {
  xs<-Filter(function(x)is.data.frame(x)&&ncol(x)>0,xs);if(!length(xs))return(data.frame())
  nm<-unique(unlist(lapply(xs,names)));do.call(rbind,lapply(xs,function(x){for(n in setdiff(nm,names(x)))x[[n]]<-rep(NA,nrow(x));x[nm]}))
}
dr_csv <- function(x,path) {
  if(!ncol(x))x<-data.frame(no_records=character())
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);con<-file(path,"wb");on.exit(close(con))
  writeBin(as.raw(c(239,187,191)),con);write.table(x,con,sep=",",row.names=FALSE,na="",quote=TRUE,qmethod="double")
}
dr_path <- function(root,p)if(grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)",p))p else file.path(root,p)
dr_root <- function(root) {
  if(nzchar(root))return(normalizePath(root,winslash="/",mustWork=TRUE))
  p<-normalizePath(getwd(),winslash="/");repeat {
    if(any(dir.exists(file.path(p,c("eua_data","prod_data")))))return(p)
    up<-dirname(p);if(up==p)break;p<-up
  };stop("Open the RStudio project, or set PHMPT_ROOT to the folder containing eua_data and prod_data.")
}
dr_sid <- function(x) {
  x<-toupper(dr_text(x));m<-regexpr("(?<![A-Z0-9])US[ _-]*[0-9]{3}[ _-]*[0-9]{4}(?![0-9])",x,perl=TRUE)
  ans<-rep("",length(x));ok<-m>0;ans[ok]<-gsub("[^A-Z0-9]","",regmatches(x,m));ans
}
dr_ids <- function(x) {
  z<-regmatches(toupper(x),gregexpr("(?<![A-Z0-9])US[ _-]*[0-9]{3}[ _-]*[0-9]{4}(?![0-9])",toupper(x),perl=TRUE))[[1]]
  unique(gsub("[^A-Z0-9]","",z))
}
dr_arm <- function(x) {
  z<-gsub("[^A-Z0-9]","",dr_norm(x));ifelse(grepl("PLACEBO",z),"Placebo",ifelse(grepl("MRNA1273",z),"mRNA-1273","UNASSIGNED"))
}
dr_projection <- function(domain) switch(domain,
  ADSL=c("STUDYID","USUBJID","SUBJID","SITEID","RANDFL","TRT01P","TRT01A","ARM","SAFFL","FASFL","DTHFL","DTHDT","DTHDTC","DTHDTF","DTHCAUS","EOSDT","EOSSTT","CUTOFFDT","TR01SDT","DOSE2DT","UBLDDT"),
  NULL)

dr_read <- function(path,domain) {
  ext<-tolower(tools::file_ext(path));wanted<-dr_projection(domain)
  if(ext=="jmp") {
    if(domain!="ADSL")stop("JMP supported for ADSL only; use XPT or a CSV export for ",domain)
    x<-ws_read_jmp(path,c("STUDYID","USUBJID"),setdiff(wanted,c("STUDYID","USUBJID")))
  }else if(ext%in%c("xpt","sas7bdat")) {
    if(!requireNamespace("haven",quietly=TRUE))stop("Install haven to read XPT/SAS7BDAT.")
    x<-if(ext=="xpt")haven::read_xpt(path) else haven::read_sas(path)
  }else if(ext=="csv") {
    enc<-if(PHMPT_CSV_ENCODING=="auto")c("UTF-8-BOM","windows-1252") else PHMPT_CSV_ENCODING
    x<-NULL;for(e in enc){q<-tryCatch(withCallingHandlers(read.csv(path,check.names=FALSE,colClasses="character",fileEncoding=e,na.strings=character()),warning=function(w)stop(conditionMessage(w))),error=function(e)NULL);if(!is.null(q)){x<-q;break}}
    if(is.null(x))stop("CSV could not be decoded; set PHMPT_CSV_ENCODING or use XPT.")
  }else stop("Unsupported file type.")
  x<-as.data.frame(x,stringsAsFactors=FALSE);names(x)<-toupper(names(x));if(anyDuplicated(names(x)))stop("Duplicate column names.")
  if(!all(c("STUDYID","USUBJID")%in%names(x)))stop("Missing STUDYID/USUBJID.")
  meta<-data.frame(variable=names(x),label=vapply(x,function(z){a<-attr(z,"label");if(is.null(a))"" else as.character(a)},""),type=vapply(x,function(z)paste(class(z),collapse="/"),""))
  dates<-intersect(names(x),c("DTHDT","EOSDT","CUTOFFDT","TR01SDT","DOSE2DT","UBLDDT","ASTDT","AENDT","ADT"))
  for(n in dates)if(is.numeric(x[[n]])&&!inherits(x[[n]],c("Date","POSIXt")))x[[n]]<-as.Date(x[[n]],origin="1960-01-01")
  x[]<-lapply(x,dr_text)
  for(n in intersect(names(x),c(dates,"DTHDTC","DSSTDTC")))x[[n]][dr_norm(x[[n]])%in%c("NA","N/A","NULL",".")]<-""
  x$.row<-seq_len(nrow(x));x<-x[dr_norm(x$STUDYID)=="MRNA-1273-P301",,drop=FALSE]
  if(!nrow(x))stop("No P301 records.")
  if(ext=="csv")for(n in dates){ii<-grepl("^-?[0-9]+([.]0+)?$",x[[n]]);if(any(ii)){if(!nzchar(PHMPT_CSV_DATE_ORIGIN))stop("Numeric CSV dates in ",n,"; set documented PHMPT_CSV_DATE_ORIGIN.");x[[n]][ii]<-format(as.Date(as.numeric(x[[n]][ii]),origin=PHMPT_CSV_DATE_ORIGIN),"%Y-%m-%d")}}
  x$.subject<-dr_sid(x$USUBJID);if(any(!grepl("^US[0-9]{7}$",x$.subject)))stop("Missing or unrecognized P301 subject ID.")
  if(domain%in%c("ADSL","DM")&&anyDuplicated(x$.subject))stop("Duplicate subject in ",domain)
  required<-switch(domain,AE="AESEQ",ADAE="AESEQ",DS="DSDECOD",DV="DVTERM",SUPPDM=c("QNAM","QVAL"),SUPPDS=c("QNAM","QVAL"),SUPPAE=c("QNAM","QVAL"),character())
  if(!all(required%in%names(x)))stop("Missing expected ",domain," fields: ",paste(setdiff(required,names(x)),collapse=", "))
  list(data=x,metadata=meta,path=normalizePath(path,winslash="/"),md5=unname(tools::md5sum(path)))
}
dr_key <- function(x) {
  nm<-sort(setdiff(names(x),c(".row")));x<-x[nm]
  apply(x,1,function(v)paste0(nchar(v,type="bytes"),":",v,collapse="|"))
}
dr_equivalent <- function(a,b) {
  nm<-sort(setdiff(names(a),".row"));if(!identical(nm,sort(setdiff(names(b),".row")))||nrow(a)!=nrow(b))return(FALSE)
  identical(sort(dr_key(a)),sort(dr_key(b)))
}
dr_pick <- function(root,rel,domain,out,inventory) {
  explicit<-PHMPT_FILES[[rel]][[domain]];if(is.null(explicit))explicit<-""
  paths<-if(nzchar(explicit))dr_path(root,explicit) else inventory[[rel]]
  if(!nzchar(explicit))paths<-paths[grepl(paste0("(^|[_ -])",tolower(domain),"([_. -]|$)"),tolower(basename(paths)))&grepl(paste0(PHMPT_SUBMISSIONS[[rel]],"_"),basename(paths),ignore.case=TRUE)&!grepl("p101|p201|dmid",basename(paths),ignore.case=TRUE)]
  paths<-paths[order(match(tolower(tools::file_ext(paths)),c("xpt","sas7bdat","csv","jmp")),paths)]
  log<-list();selected<-NULL;conflict<-FALSE;seen<-character();successful_stems<-character()
  add<-function(p,status,detail="")log[[length(log)+1L]]<<-data.frame(release=rel,domain=domain,path=p,status=status,detail=detail)
  for(p in paths) {
    ext<-tolower(tools::file_ext(p));stem<-tools::file_path_sans_ext(p)
    if(ext=="csv"&&stem%in%successful_stems){add(p,"CSV_MIRROR_NOT_USED","Readable same-stem original selected; CSV not merged.");next}
    hash<-if(file.exists(p))unname(tools::md5sum(p)) else ""
    if(nzchar(hash)&&hash%in%seen){add(p,"BYTE_IDENTICAL_COPY");next}
    message("Checking ",rel," ",domain,": ",basename(p))
    a<-tryCatch(dr_read(p,domain),error=function(e)e)
    if(inherits(a,"error")){add(p,"UNREADABLE",conditionMessage(a));next}
    if(domain=="ADSL"&&"CUTOFFDT"%in%names(a$data)) {
      dates<-unique(a$data$CUTOFFDT[nzchar(a$data$CUTOFFDT)]);want<-as.character(if(rel=="EUA")PHMPT_EUA_CUTOFF else PHMPT_BLA_CUTOFF)
      if(length(dates)&&!all(substr(dates,1,10)==want)){add(p,"WRONG_CUTOFF",paste(dates,collapse="; "));next}
    }
    seen<-c(seen,a$md5);if(ext%in%c("xpt","sas7bdat"))successful_stems<-c(successful_stems,stem)
    if(is.null(selected)){selected<-a;add(p,"SELECTED")}
    else if(dr_equivalent(selected$data,a$data))add(p,"EQUIVALENT_PROJECTED_VALUES")
    else {conflict<-TRUE;add(p,"CONFLICTING_SNAPSHOT","Set PHMPT_FILES explicitly. This domain is not used.")}
  }
  if(!length(paths))add("","MISSING","No file in the configured submission.")
  lg<-dr_bind(log);dr_csv(lg,file.path(out,paste0("discovery_",rel,"_",domain,".csv")))
  if(conflict)selected<-NULL
  if(!is.null(selected)){selected$release<-rel;selected$domain<-domain}
  list(input=selected,log=lg,status=if(conflict)"CONFLICT" else if(is.null(selected))"UNAVAILABLE" else "AVAILABLE")
}
# Death date bounds use only explicit death dates and DEATH disposition dates.
# AE onset/end dates remain context; neither is silently substituted for death.
dr_day_bounds <- function(x) {
  x<-dr_text(x);lo<-hi<-rep(NA_real_,length(x));precision<-rep("INVALID_OR_MISSING",length(x))
  for(i in seq_along(x)) {
    s<-substr(x[i],1,10)
    if(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$",s)){
      d<-tryCatch(as.Date(s),error=function(e)as.Date(NA));if(!is.na(d)&&format(d,"%Y-%m-%d")==s){lo[i]<-hi[i]<-as.numeric(d);precision[i]<-"DAY"}
    }else if(grepl("^[0-9]{4}-[0-9]{2}$",x[i])){
      d<-tryCatch(as.Date(paste0(x[i],"-01")),error=function(e)as.Date(NA));if(!is.na(d)){lo[i]<-as.numeric(d);hi[i]<-as.numeric(seq(d,by="month",length.out=2)[2]-1);precision[i]<-"MONTH"}
    }else if(grepl("^[0-9]{4}$",x[i])){lo[i]<-as.numeric(as.Date(paste0(x[i],"-01-01")));hi[i]<-as.numeric(as.Date(paste0(x[i],"-12-31")));precision[i]<-"YEAR"}
  };data.frame(lower=lo,upper=hi,precision=precision)
}
dr_signals <- function(input) {
  x<-input$data;dom<-input$domain;rel<-input$release;rows<-list()
  add<-function(ii,reason,strength="STRUCTURED",date_field="",date_value=NULL,imputed=NULL) {
    if(!length(ii))return();if(is.null(date_value))date_value<-if(nzchar(date_field))dr_col(x,date_field)[ii] else rep("",length(ii))
    if(is.null(imputed))imputed<-rep(FALSE,length(ii));bd<-dr_day_bounds(date_value);bd$lower[imputed]<-bd$upper[imputed]<-NA_real_;bd$precision[imputed]<-"IMPUTED_NOT_USED"
    rows[[length(rows)+1L]]<<-data.frame(release=rel,domain=dom,subject=x$.subject[ii],source_row=x$.row[ii],signal=reason,strength=strength,
      date_field=date_field,date_raw=date_value,date_lower=bd$lower,date_upper=bd$upper,date_precision=bd$precision,
      event_id=if(dom%in%c("AE","ADAE"))paste0("AESEQ=",dr_col(x,"AESEQ")[ii],"; AESPID=",dr_col(x,"AESPID")[ii]) else "",
      term=if(dom%in%c("AE","ADAE"))paste(dr_col(x,"AETERM")[ii],dr_col(x,"AEDECOD")[ii],sep=" | ") else if(dom=="DS")dr_col(x,"DSDECOD")[ii] else "",
      ae_start=dr_col(x,"AESTDTC")[ii],ae_end=dr_col(x,"AEENDTC")[ii],source_file=input$path,source_md5=input$md5,stringsAsFactors=FALSE)
  }
  if(dom%in%c("ADSL","DM")) {
    add(which(dr_norm(dr_col(x,"DTHFL"))%in%c("Y","YES")),"DTHFL=Y")
    for(f in intersect(c("DTHDTC","DTHDT"),names(x)))add(which(nzchar(dr_col(x,f))),paste0(f," recorded"),date_field=f,
      imputed=if(f=="DTHDT")nzchar(dr_col(x,"DTHDTF"))[nzchar(dr_col(x,f))] else NULL)
    add(which(dr_norm(dr_col(x,"EOSSTT"))%in%c("DEATH","DEAD","DIED")),"EOSSTT indicates death")
  }
  if(dom=="DS") {ii<-which(dr_norm(dr_col(x,"DSDECOD"))%in%c("DEATH","DEAD","DIED")|dr_norm(dr_col(x,"DSTERM"))%in%c("DEATH","DEAD","DIED"));add(ii,"DEATH disposition",date_field="DSSTDTC")}
  if(dom%in%c("AE","ADAE")) {
    yy<-dr_norm(dr_col(x,"AESDTH"))%in%c("Y","YES")|dr_norm(dr_col(x,"AEOUT"))%in%c("FATAL","DEATH","DEAD","DIED")
    add(which(yy),"AESDTH=Y and/or fatal AE outcome")
    for(f in intersect(c("DTHDTC","DTHDT"),names(x)))add(which(yy&nzchar(dr_col(x,f))),paste0("Fatal AE with ",f),date_field=f)
    kw<-grepl("DEATH|DIED|DECEAS|FATAL|SUICIDE",paste(dr_norm(dr_col(x,"AETERM")),dr_norm(dr_col(x,"AEDECOD"))))
    add(which(kw&!yy),"Death-related term only: participant death not established","TEXT_ONLY")
  }
  if(dom%in%c("SUPPDM","SUPPAE","SUPPDS")) {
    q<-dr_norm(dr_col(x,"QNAM"));v<-dr_col(x,"QVAL")
    add(which(q%in%c("DTHFL","AESDTH")&dr_norm(v)%in%c("Y","YES")),"Supplementary explicit death flag")
    for(f in c("DTHDTC","DTHDT")){ii<-which(q==f&nzchar(v));add(ii,paste("Supplementary",f),date_field=paste0("QVAL:",f),date_value=v[ii])}
  }
  dr_bind(rows)
}
dr_inventory_fields <- function(input) {
  x<-input$data;m<-input$metadata
  m$release<-input$release;m$domain<-input$domain;m$path<-input$path
  m$n_nonblank<-vapply(m$variable,function(n)sum(nzchar(dr_col(x,n))),integer(1))
  m$reporting_date_candidate<-grepl("AWARE|NOTIF|RECEI|REPORT|CREAT|ENTER|ENTRY|AUDIT",paste(m$variable,dr_norm(m$label)))
  if(all(c("QNAM","QVAL")%in%names(x))) {
    q<-unique(data.frame(variable=dr_col(x,"QNAM"),label=dr_col(x,"QLABEL")))
    q<-q[nzchar(q$variable),];q$type<-"SUPPLEMENTARY_QUALIFIER";q$release<-input$release;q$domain<-input$domain;q$path<-input$path
    q$n_nonblank<-vapply(q$variable,function(n)sum(dr_col(x,"QNAM")==n&nzchar(dr_col(x,"QVAL"))),integer(1))
    q$reporting_date_candidate<-grepl("AWARE|NOTIF|RECEI|REPORT|CREAT|ENTER|ENTRY|AUDIT",paste(q$variable,dr_norm(q$label)));m<-dr_bind(list(m,q))
  };m
}

# The uncompressed JMP column layouts below are adapted from JMPReader.jl:
# https://github.com/jaakkor2/JMPReader.jl/blob/main/src/column.jl
# MIT License. Copyright (C) 2024 Jaakko Ruohio.
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
ws_read_jmp <- function(path, fields, optional=character()) {
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
  selected <- unique(c(fields, intersect(optional, column_names)))
  values <- setNames(vector("list", length(selected)), selected)
  for (name in selected) {
    i <- match(name, column_names); p <- offsets[i] + 2 + uint(b, offsets[i], 2)
    dt <- as.integer(take(b, p, 6)); p <- p + 6
    end <- if (i < nc) offsets[i + 1] - 1 else length(b)
    compressed <- dt[1] %in% c(9L,10L)
    if (!dt[1] %in% c(1L,2L,9L,10L))
      stop("Unsupported JMP encoding in ", name, " (descriptor ", paste(dt,collapse=","), ").")
    if (compressed) {
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
    } else payload <- take(b,p,end-p+1)
    fixed_char <- dt[1] %in% c(2L,9L) && dt[2] %in% c(1L,2L) &&
      dt[4] == 0L && dt[3] %in% 0:7 && dt[5] > 0L
    if (fixed_char) {
      width <- dt[5]
      if ((compressed && length(payload) != nr * width) || length(payload) < nr * width)
        stop("Unsupported fixed-width JMP payload in ", name)
      # Uncompressed values occupy the tail of the bounded column block.
      payload <- take(payload,length(payload)-nr*width+1,nr*width)
      values[[name]] <- vapply(seq_len(nr), function(j) txt(take(payload, (j - 1) * width + 1, width)), character(1))
    } else if (dt[1] == 2L && dt[2] %in% c(1L,2L) && all(dt[3:5] == 0L)) {
      # EUA ADSL uses uncompressed variable-length UTF-8 strings. Read the
      # length vector from its metadata, then the string bytes at column end.
      q <- 7L; n1 <- uint(payload,q,8); q <- q+8+n1+2
      n2 <- uint(payload,q,4); q <- q+4+n2+8
      wb <- uint(payload,q,1); maxwidth <- uint(payload,q+1,4); q <- q+5
      if (!wb %in% c(1L,2L,4L)) stop("Unsupported JMP string-length width in ",name,": ",wb)
      lenbytes <- take(payload,q,nr*wb)
      lens <- if(wb==1L)as.integer(lenbytes) else readBin(lenbytes,"integer",n=nr,size=wb,endian="little",signed=TRUE)
      if (anyNA(lens) || any(lens<0) || any(lens>maxwidth)) stop("Invalid JMP string lengths in ",name)
      total <- sum(as.double(lens)); data_at <- length(payload)-total+1
      if (data_at < q+nr*wb) stop("JMP string data overlaps metadata in ",name)
      starts <- data_at+c(0,head(cumsum(as.double(lens)),-1L))
      values[[name]] <- vapply(seq_len(nr),function(j)txt(take(payload,starts[j],lens[j])),character(1))
    } else if (dt[1] %in% c(1L,10L) && dt[6] == 8L && dt[4] == dt[5] && dt[4] %in% c(0L,3L,99L,114L)) {
      if ((compressed && length(payload) != nr * 8) || length(payload) < nr * 8)
        stop("JMP numeric row count mismatch in ", name)
      payload <- take(payload,length(payload)-nr*8+1,nr*8)
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






# ---- Reconcile participant-level evidence without inventing reporting dates ----
dr_iso <- function(x)ifelse(is.na(x),"",as.character(as.Date(x,origin="1970-01-01")))
dr_join <- function(x)paste(sort(unique(x[!is.na(x)&nzchar(x)])),collapse="; ")
dr_subjects <- function(evidence,inputs,coverage) {
  ids<-sort(unique(c(evidence$subject,dr_sid(PHMPT_EXTRA_SUBJECTS))));ids<-ids[nzchar(ids)]
  records<-dr_bind(lapply(inputs,function(a){x<-a$data;if(!a$domain%in%c("ADSL","DM"))return(NULL);data.frame(subject=x$.subject,release=a$release,domain=a$domain,arm=if(a$domain=="ADSL")dr_col(x,"TRT01P") else dr_col(x,"ARM"),rand=dr_col(x,"RANDFL"),dthfl=dr_col(x,"DTHFL"))}))
  core<-coverage[coverage$release=="EUA"&coverage$domain%in%c("ADSL","DM","AE"),]
  complete<-nrow(core)==3L&&all(core$status=="AVAILABLE")&&all(nzchar(core$death_evidence_fields))
  any_eua<-any(coverage$release=="EUA"&coverage$status=="AVAILABLE"&coverage$domain!="DV"&nzchar(coverage$death_evidence_fields))
  dr_bind(lapply(ids,function(id) {
    z<-evidence[evidence$subject==id,,drop=FALSE];s<-z[z$strength=="STRUCTURED",,drop=FALSE]
    eua<-any(s$release=="EUA");bla<-any(s$release=="BLA");dated<-s[!is.na(s$date_lower)&!is.na(s$date_upper),,drop=FALSE]
    conflict<-FALSE;lo<-hi<-NA_real_
    if(nrow(dated)){lo<-max(dated$date_lower);hi<-min(dated$date_upper);conflict<-lo>hi;if(conflict){lo<-min(dated$date_lower);hi<-max(dated$date_upper)}}
    r<-if(nrow(records))records[records$subject==id,,drop=FALSE] else data.frame()
    arms<-if(nrow(r))unique(dr_arm(r$arm)) else character();arms<-setdiff(arms,"UNASSIGNED")
    arm<-if(length(arms)==1L)arms else if(length(arms)>1L)"ARM_CONFLICT" else "UNASSIGNED"
    before<-function(d)!conflict&&!is.na(hi)&&hi<=as.numeric(d)
    missing_eua<-bla&&!eua&&before(PHMPT_EUA_CUTOFF)
    priority<-if(conflict)"DEATH_DATE_CONFLICT" else if(missing_eua&&complete)"PRE_CUTOFF_DEATH_NOT_IDENTIFIED_AT_EUA" else if(missing_eua)"PRE_CUTOFF_EUA_RECONCILIATION_INCOMPLETE_COVERAGE" else if(nrow(s)&&is.na(hi))"DEATH_DATE_UNRESOLVED" else if(!nrow(s)&&nrow(z))"TEXT_ONLY_REVIEW" else "RECONCILIATION"
    flags<-if(nrow(r))dr_norm(r$dthfl) else character()
    data.frame(subject=id,site=substr(id,1,5),arm=arm,randomized_flags=if(nrow(r))dr_join(r$rand) else "",
      structured_death_evidence=nrow(s)>0,EUA_death_evidence=eua,BLA_death_evidence=bla,
      EUA_core_coverage_complete=complete,EUA_any_death_domain_available=any_eua,
      EUA_presence=if(eua)"PRESENT" else if(!any_eua)"EUA_DEATH_DOMAINS_UNAVAILABLE" else if(complete)"NOT_IDENTIFIED_IN_SELECTED_TABLES" else "NOT_IDENTIFIED_IN_AVAILABLE_TABLES_INCOMPLETE_COVERAGE",
      death_date_lower=dr_iso(lo),death_date_upper=dr_iso(hi),death_date_conflict=conflict,
      EUA_recorded_death_dates=dr_join(s$date_raw[s$release=="EUA"]),BLA_recorded_death_dates=dr_join(s$date_raw[s$release=="BLA"]),
      death_by_EUA_cutoff=before(PHMPT_EUA_CUTOFF),death_by_Dec03=before(PHMPT_EUA_UPDATE),death_by_authorization=before(PHMPT_REVIEW_END),
      pre_cutoff_missing_EUA_evidence=if(complete)missing_eua else NA,death_flag_N_present=any(flags%in%c("N","NO"))&&nrow(s)>0,
      evidence_domains=dr_join(paste(s$release,s$domain,sep=":")),priority=priority,
      review_focus=(nrow(s)>0||PHMPT_INCLUDE_TEXT_ONLY||id%in%dr_sid(PHMPT_EXTRA_SUBJECTS))&&(conflict||is.na(lo)||lo<=as.numeric(PHMPT_REVIEW_END)||id%in%dr_sid(PHMPT_EXTRA_SUBJECTS)),
      stringsAsFactors=FALSE)
  }))
}
dr_review_template <- function(subjects) {
  ids<-subjects$subject[subjects$review_focus&subjects$structured_death_evidence]
  if(!length(ids))ids<-character()
  x<-data.frame(subject=rep(ids,each=2),notification_type=rep(c("INITIAL_SAE","DEATH_UPDATE"),length(ids)),stringsAsFactors=FALSE)
  for(n in c("event_key","investigator_awareness","sponsor_notification","awareness_timezone","notification_timezone","verified_awareness","verified_first_notification","verified_event_link","applicable_24h","source_document","source_page","reviewer","review_notes"))x[[n]]<-rep("",nrow(x))
  x$event_key[x$notification_type=="DEATH_UPDATE"]<-"DEATH"
  x
}
# Strict ISO timestamps: timezone-bearing ISO timestamps or a separately specified
# IANA timezone. Date-only records give bounds, not an invented noon timestamp.
dr_clock <- function(x,tz="") {
  fail<-function(reason)list(lo=NA_real_,hi=NA_real_,precision=reason)
  x<-trimws(x);tz<-trimws(tz);if(!nzchar(x))return(fail("MISSING"))
  hasoffset<-grepl("(Z|[+-][0-9]{2}:?[0-9]{2})$",x)
  if(!hasoffset&&(!nzchar(tz)||!tz%in%c("UTC","GMT",OlsonNames())))return(fail("MISSING_OR_INVALID_TIMEZONE"))
  if(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$",x)) {
    d<-tryCatch(as.Date(x),error=function(e)as.Date(NA));if(is.na(d)||as.character(d)!=x)return(fail("INVALID_DATE"))
    a<-as.POSIXct(paste(x,"00:00:00"),tz=tz);b<-as.POSIXct(paste(d+1,"00:00:00"),tz=tz)
    return(list(lo=as.numeric(a),hi=as.numeric(b)-1,precision="DAY"))
  }
  if(!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}[T ][0-9]{2}:[0-9]{2}:[0-9]{2}(Z|[+-][0-9]{2}:?[0-9]{2})?$",x))return(fail("USE_ISO_DATE_OR_TIMESTAMP_WITH_SECONDS"))
  s<-gsub("T"," ",x,fixed=TRUE)
  if(hasoffset&&!grepl("Z$",x)){off<-sub("^.*([+-][0-9]{2}:?[0-9]{2})$","\\1",x);digits<-gsub("[^0-9]","",off);if(as.integer(substr(digits,1,2))>14||as.integer(substr(digits,3,4))>59||(as.integer(substr(digits,1,2))==14&&as.integer(substr(digits,3,4))!=0))return(fail("INVALID_UTC_OFFSET"))}
  if(hasoffset){s<-sub("Z$","+0000",s);s<-sub("([+-][0-9]{2}):([0-9]{2})$","\\1\\2",s);a<-as.POSIXct(strptime(s,"%Y-%m-%d %H:%M:%S%z",tz="UTC"))}
  else a<-as.POSIXct(strptime(s,"%Y-%m-%d %H:%M:%S",tz=tz))
  if(is.na(a))return(fail("INVALID_TIMESTAMP"))
  # Reject invalid calendar dates even if the platform normalizes them.
  d<-tryCatch(as.Date(substr(x,1,10)),error=function(e)as.Date(NA));if(is.na(d)||as.character(d)!=substr(x,1,10))return(fail("INVALID_DATE"))
  if(!hasoffset) {
    candidates<-as.numeric(a)+c(-7200,-3600,0,3600,7200)
    candidates<-candidates[format(as.POSIXct(candidates,origin="1970-01-01",tz=tz),"%Y-%m-%d %H:%M:%S",tz=tz)==s]
    if(!length(candidates))return(fail("NONEXISTENT_LOCAL_TIME"))
    return(list(lo=min(candidates),hi=max(candidates),precision=if(length(candidates)>1)"AMBIGUOUS_LOCAL_TIME" else "SECOND"))
  }
  list(lo=as.numeric(a),hi=as.numeric(a),precision="SECOND")
}
dr_assess <- function(review,subjects) {
  template<-dr_review_template(subjects)
  for(n in setdiff(names(template),names(review)))review[[n]]<-rep("",nrow(review))
  review[]<-lapply(review,dr_text);review$subject<-dr_sid(review$subject)
  review$notification_type<-dr_norm(review$notification_type)
  key<-paste(review$subject,review$notification_type,review$event_key,sep="|")
  dup<-duplicated(key)|duplicated(key,fromLast=TRUE)
  review$arm<-subjects$arm[match(review$subject,subjects$subject)];review$arm[is.na(review$arm)]<-"UNMATCHED_REVIEW_SUBJECT"
  review$review_focus<-subjects$review_focus[match(review$subject,subjects$subject)]
  review$delay_lower_hours<-review$delay_upper_hours<-rep(NA_real_,nrow(review))
  review$assessment<-rep("NOT_ASSESSABLE",nrow(review));review$assessment_reason<-rep("",nrow(review))
  review$awareness_precision<-review$notification_precision<-rep("",nrow(review))
  yes<-function(x)dr_norm(x)%in%c("Y","YES","TRUE")
  for(i in seq_len(nrow(review))) {
    r<-review[i,];reason<-"";status<-"NOT_ASSESSABLE"
    a<-dr_clock(r$investigator_awareness,r$awareness_timezone);b<-dr_clock(r$sponsor_notification,r$notification_timezone)
    review$awareness_precision[i]<-a$precision;review$notification_precision[i]<-b$precision
    if(!nzchar(r$subject)||!r$subject%in%subjects$subject)reason<-"Subject absent from the evidence/review queue; add PHMPT_EXTRA_SUBJECTS if appropriate."
    else if(dup[i]){status<-"CONFLICTING_REVIEW_ROWS";reason<-"Duplicate subject/type/event key; resolve rather than selecting a preferred row."}
    else if(!r$notification_type%in%c("INITIAL_SAE","DEATH_UPDATE"))reason<-"Use INITIAL_SAE or DEATH_UPDATE."
    else if(!nzchar(r$event_key)||!all(yes(c(r$verified_awareness,r$verified_first_notification,r$verified_event_link,r$applicable_24h)))||!nzchar(r$source_document)||!nzchar(r$reviewer))reason<-"Verified awareness, first notification, event link, 24h applicability, event key, source and reviewer are required."
    else if(anyNA(c(a$lo,b$lo)))reason<-paste("Unresolved time:",a$precision,b$precision)
    else {
      lo<-(b$lo-a$hi)/3600;hi<-(b$hi-a$lo)/3600
      if(hi<0){status<-"INCONSISTENT_ORDER";reason<-"Notification precedes awareness; verify the event and timestamp meanings."}
      else {lo<-max(0,lo);review$delay_lower_hours[i]<-lo;review$delay_upper_hours[i]<-hi
        status<-if(lo>24)"EXCEEDS_24H" else if(hi<=24)"WITHIN_24H_BOUND" else "INDETERMINATE_PRECISION"
        reason<-if(a$precision=="DAY"||b$precision=="DAY")"Date-only bounds; lower bound assumes notification followed awareness." else if(a$precision=="AMBIGUOUS_LOCAL_TIME"||b$precision=="AMBIGUOUS_LOCAL_TIME")"Bounds include both offsets for an ambiguous local clock time." else "Verified awareness-to-first-notification timestamps."
      }
    };review$assessment[i]<-status;review$assessment_reason[i]<-reason
  };review
}
dr_h <- function(x){x<-dr_text(x);for(z in list(c("&","&amp;"),c("<","&lt;"),c(">","&gt;"),c('"',"&quot;"),c("'","&#39;")))x<-gsub(z[1],z[2],x,fixed=TRUE);x}
dr_table <- function(x) {
  if(!nrow(x))return("<p class='muted'>No rows available.</p>")
  paste0("<div class='scroll'><table><thead><tr>",paste0("<th>",dr_h(names(x)),"</th>",collapse=""),"</tr></thead><tbody>",paste(apply(x,1,function(z)paste0("<tr>",paste0("<td>",dr_h(z),"</td>",collapse=""),"</tr>")),collapse="\n"),"</tbody></table></div>")
}
# ---- Full document export; page hits are navigation aids, not verified events ----
dr_crfs <- function(root,subjects,export) {
  ids<-subjects$subject[subjects$review_focus]
  files<-unique(unlist(lapply(PHMPT_PDF_ROOTS,function(p)list.files(dr_path(root,p),pattern="[.]pdf$",recursive=TRUE,full.names=TRUE,ignore.case=TRUE))))
  files<-files[!grepl("p101|p201|dmid",basename(files),ignore.case=TRUE)]
  paths_ids<-lapply(basename(files),dr_ids)
  anonymous<-data.frame(path=files[lengths(paths_ids)==0&grepl("CRF|case.report.form",basename(files),ignore.case=TRUE)])
  dr_csv(anonymous,file.path(export,"unassigned_crf_files.csv"))
  links<-dr_bind(lapply(seq_along(files),function(i){s<-intersect(paths_ids[[i]],ids);if(!length(s))return(NULL);data.frame(subject=s,path=files[i],match_method="SUBJECT_ID_IN_FILENAME")}))
  if(nzchar(PHMPT_PDF_LINKS_FILE)) {
    p<-dr_path(root,PHMPT_PDF_LINKS_FILE);x<-read.csv(p,colClasses="character",fileEncoding="UTF-8-BOM",check.names=FALSE,na.strings=character())
    if(!all(c("subject","path")%in%names(x)))stop("PHMPT_PDF_LINKS_FILE needs subject,path columns.")
    x$subject<-dr_sid(x$subject);x$path<-vapply(x$path,function(p)dr_path(root,p),"");x$match_method<-"EXPLICIT_REVIEWER_LINK"
    links<-dr_bind(list(links,x[x$subject%in%ids,c("subject","path","match_method"),drop=FALSE]))
  }
  prototype<-data.frame(subject=character(),document_type=character(),match_method=character(),original_filename=character(),source_path=character(),md5=character(),bytes=numeric(),export_path=character(),copy_status=character(),pages=integer(),text_status=character())
  hits_proto<-data.frame(subject=character(),export_path=character(),pdf_page=integer(),keywords=character(),snippet=character())
  if(!nrow(links))return(list(manifest=prototype,hits=hits_proto,anonymous=nrow(anonymous)))
  links<-unique(links);rows<-hits<-list();cache<-new.env(parent=emptyenv())
  dir.create(file.path(export,"documents"),showWarnings=FALSE)
  for(i in seq_len(nrow(links))) {
    p<-links$path[i];id<-links$subject[i];b<-basename(p);message("CRF/document ",i,"/",nrow(links),": ",b)
    type<-if(grepl("CRF|case.report.form",b,ignore.case=TRUE))"CRF" else "OTHER_SUBJECT_DOCUMENT"
    hash<-if(file.exists(p))unname(tools::md5sum(p)) else "";dest<-"";status<-"SOURCE_MISSING";tx<-character();pages<-NA_integer_;text_status<-"NOT_SCANNED"
    if(nzchar(hash)) {
      dest<-paste0("documents/",if(type=="CRF")"CRF_" else "DOC_",hash,".pdf")
      target<-file.path(export,dest)
      if(PHMPT_COPY_PDFS) {
        ok<-if(file.exists(target))TRUE else file.copy(p,target,overwrite=FALSE)
        status<-if(ok&&identical(unname(tools::md5sum(target)),hash))"COPIED_VERIFIED" else "COPY_OR_HASH_FAILED"
      }else {status<-"COPY_DISABLED";dest<-""}
      if(exists(hash,envir=cache,inherits=FALSE)) {q<-get(hash,envir=cache);tx<-q$text;pages<-q$pages;text_status<-q$status}
      else {
        if(PHMPT_PDF_TEXT_SCAN&&requireNamespace("pdftools",quietly=TRUE)) {
          tx<-tryCatch(pdftools::pdf_text(p),error=function(e){text_status<<-paste("PDF_TEXT_ERROR",conditionMessage(e));character()})
          if(length(tx)){pages<-length(tx);text_status<-if(any(nzchar(trimws(tx))))"TEXT_EXTRACTED" else "NO_TEXT_OCR_REQUIRED"}
        }else if(PHMPT_PDF_TEXT_SCAN&&nzchar(Sys.which("pdftotext"))) {
          tf<-tempfile(fileext=".txt");ef<-tempfile(fileext=".txt")
          rc<-tryCatch(system2(Sys.which("pdftotext"),c("-layout","-enc","UTF-8",shQuote(p),shQuote(tf)),stdout=FALSE,stderr=ef),error=function(e)1L)
          if(identical(rc,0L)&&file.exists(tf)) {
            tx<-strsplit(paste(readLines(tf,warn=FALSE,encoding="UTF-8"),collapse="\n"),"\f",fixed=TRUE)[[1]]
            if(length(tx)&&!nzchar(trimws(tail(tx,1))))tx<-head(tx,-1)
            pages<-length(tx);text_status<-if(any(nzchar(trimws(tx))))"TEXT_EXTRACTED_PDFTOTEXT" else "NO_TEXT_OCR_REQUIRED"
          }else text_status<-"PDFTOTEXT_ERROR"
          unlink(c(tf,ef))
        }else if(PHMPT_PDF_TEXT_SCAN)text_status<-"PDFTOOLS_OR_PDFTOTEXT_UNAVAILABLE"
        assign(hash,list(text=tx,pages=pages,status=text_status),envir=cache)
      }
      if(length(tx))for(j in seq_along(tx)) {
        pat<-c(DEATH="\\bdeath\\b|\\bdied\\b|deceased",FATAL="\\bfatal\\b",AWARENESS="aware|awareness",NOTIFICATION="notif|reported to|report date|date reported|received|receipt|audit trail",SAE="serious adverse|\\bSAE\\b")
        found<-names(pat)[vapply(pat,function(k)grepl(k,tx[j],ignore.case=TRUE,perl=TRUE),logical(1))]
        if(length(found)){
          s<-gsub("[[:space:]]+"," ",tx[j]);at<-regexpr(paste(pat,collapse="|"),s,ignore.case=TRUE,perl=TRUE)[1];start<-max(1,at-100)
          hits[[length(hits)+1L]]<-data.frame(subject=id,export_path=dest,pdf_page=j,keywords=paste(found,collapse="; "),snippet=substr(s,start,start+440))
        }
      }
    }
    rows[[length(rows)+1L]]<-data.frame(subject=id,document_type=type,match_method=links$match_method[i],original_filename=b,source_path=p,md5=hash,bytes=if(file.exists(p))file.info(p)$size else NA_real_,export_path=dest,copy_status=status,pages=pages,text_status=text_status)
  }
  list(manifest=dr_bind(c(list(prototype),rows)),hits=unique(dr_bind(c(list(hits_proto),hits))),anonymous=nrow(anonymous))
}
dr_fatal_changes <- function(inputs,ids) {
  get<-function(rel){a<-Filter(function(a)a$release==rel&&a$domain=="AE",inputs);if(!length(a))return(NULL);x<-a[[1]]$data;x[x$.subject%in%ids,,drop=FALSE]}
  a<-get("EUA");b<-get("BLA");if(is.null(a)||is.null(b))return(data.frame())
  key<-function(x)paste(x$.subject,dr_col(x,"AESEQ"),sep="|")
  ka<-key(a);kb<-key(b);keys<-union(ka,kb);dup<-union(ka[duplicated(ka)],kb[duplicated(kb)])
  rows<-lapply(keys,function(k){
    x<-a[ka==k,,drop=FALSE];y<-b[kb==k,,drop=FALSE]
    fatal<-function(z)any(dr_norm(dr_col(z,"AESDTH"))%in%c("Y","YES")|dr_norm(dr_col(z,"AEOUT"))%in%c("FATAL","DEATH","DIED"))
    if(!fatal(x)&&!fatal(y))return(NULL)
    v<-function(z,n)dr_join(dr_col(z,n));same<-nrow(x)==1L&&nrow(y)==1L&&identical(v(x,"AETERM"),v(y,"AETERM"))&&identical(v(x,"AESTDTC"),v(y,"AESTDTC"))
    data.frame(subject=if(nrow(x))x$.subject[1] else y$.subject[1],AESEQ=if(nrow(x))dr_col(x,"AESEQ")[1] else dr_col(y,"AESEQ")[1],EUA_present=nrow(x)>0,BLA_present=nrow(y)>0,
      pairing_status=if(k%in%dup)"NONUNIQUE_KEY" else if(!nrow(x)||!nrow(y))"UNMATCHED_AESEQ" else if(same)"AESEQ_TERM_ONSET_AGREE" else "AESEQ_ONLY_VERIFY_EVENT_IDENTITY",
      EUA_term=v(x,"AETERM"),BLA_term=v(y,"AETERM"),EUA_onset=v(x,"AESTDTC"),BLA_onset=v(y,"AESTDTC"),
      EUA_AESER=v(x,"AESER"),BLA_AESER=v(y,"AESER"),EUA_AESDTH=v(x,"AESDTH"),BLA_AESDTH=v(y,"AESDTH"),EUA_AEOUT=v(x,"AEOUT"),BLA_AEOUT=v(y,"AEOUT"),
      review_reason=if(!nrow(x)||!nrow(y))"Presence differs; verify follow-up and record lineage" else if(!same)"Verify event identity before comparing flags" else if(v(x,"AESER")!=v(y,"AESER")||v(x,"AESDTH")!=v(y,"AESDTH")||v(x,"AEOUT")!=v(y,"AEOUT"))"Seriousness/fatality/outcome differs; verify dated update" else "No change in selected flags")
  });dr_bind(rows)
}
dr_counts <- function(subjects) {
  if(!nrow(subjects))return(data.frame())
  dr_bind(lapply(c("death_by_EUA_cutoff","death_by_Dec03","death_by_authorization"),function(w)dr_bind(lapply(sort(unique(subjects$arm)),function(arm){z<-subjects[subjects$arm==arm&subjects$structured_death_evidence,,drop=FALSE];d<-z[z[[w]],,drop=FALSE];data.frame(window=switch(w,death_by_EUA_cutoff=as.character(PHMPT_EUA_CUTOFF),death_by_Dec03=as.character(PHMPT_EUA_UPDATE),death_by_authorization=as.character(PHMPT_REVIEW_END)),arm=arm,dated_deaths_reconstructed=nrow(d),with_EUA_death_evidence=if(all(subjects$EUA_any_death_domain_available))sum(d$EUA_death_evidence) else NA_integer_,without_EUA_death_evidence=if(all(subjects$EUA_core_coverage_complete))sum(!d$EUA_death_evidence) else NA_integer_,undated_or_conflicting_deaths=sum(!nzchar(z$death_date_upper)|z$death_date_conflict))}))))
}
dr_methods <- function() c(
  "MODERNA P301: DEATH REPORTING REVIEW",paste("Script version",DR_VERSION),
  "Run from the RStudio project root. Reads local release files; no trial data are downloaded.",
  "Defaults: EUA safety cutoff 2020-11-25; FDA updated deaths through 2020-12-03; authorization 2020-12-18; BLA snapshot 2021-03-26.",
  "The FDA review described 13 deaths (6 vaccine, 7 placebo) as reported through December 3. This is contextual, not an expected count for the November 25 snapshot.",
  "Protocol Amendment 4 (30 September 2020), sections 8.3.6-8.3.10, PDF pages 78-80: SAE notification within 24 hours of awareness; updated SAE information within 24 hours of availability.",
  "Initial SAE notification and notification of the subsequent death/fatal outcome are assessed separately. AE onset and end dates are not used as death or reporting dates.",
  "Death evidence: DTHFL, explicit death dates, death disposition, fatal AE outcome/AESDTH and corresponding supplementary qualifiers. Text-only death/suicide mentions are review leads, not confirmed participant deaths.",
  "All subjects with evidence are retained regardless of causality, randomization or treatment status. Arm uses original planned TRT01P/ARM; discordant labels are flagged. No mortality rate or vaccine efficacy is estimated.",
  "Partial death dates remain intervals. Conflicting dates are not resolved by choosing one source. Imputed DTHDT values are not used to establish timing. AESEQ matches with different term/onset require event-identity review.",
  "Retrospective death dates from BLA identify reconciliation questions, not proof that investigators or FDA knew a death by that date. A difference in outcome may reflect legitimate follow-up. Absence from one released table is not proof of nonreporting.",
  "Review focus includes deaths possibly occurring by authorization, undated/conflicting deaths, text-only review leads and explicitly requested subjects. Later deaths remain in the full subject listing but their documents are not exported by default.",
  "Review workflow: fill assets/death_reporting_review.csv and rerun. The script never overwrites that persistent worksheet. New subjects/types found on subsequent runs are offered in new_review_rows_to_add.csv and included as unverified rows in that run; resolve any duplicate event keys.",
  "For INITIAL_SAE, awareness means awareness that the event met seriousness criteria. For DEATH_UPDATE, it means the death/fatal-outcome information became available to the investigator/site. Sponsor notification is the first valid notification by any permitted route, not necessarily EDC entry.",
  "Required worksheet confirmations (Y): verified_awareness, verified_first_notification, verified_event_link, applicable_24h. Supply event_key, source_document, source_page where available, reviewer. Retain each distinct notification as its own event_key.",
  "Use ISO dates YYYY-MM-DD or timestamps YYYY-MM-DDTHH:MM:SSZ / +/-HH:MM. Local timestamps/date-only entries require an IANA timezone (e.g. America/New_York). Do not infer a timezone solely from trial-site location if the system used another clock.",
  "A date-only entry spans that local day. Delay lower/upper bounds combine the notification and awareness intervals. A lower bound above 24h is flagged EXCEEDS_24H; an upper bound <=24h is WITHIN_24H_BOUND; crossing 24h is INDETERMINATE_PRECISION. Date-only lower bounds are conditioned on notification following awareness.",
  "An unverified date, unclear timestamp meaning, missing timezone, conflicting review row or unresolved event link is not counted as timely. The script does not assess sponsor-to-FDA regulatory deadlines or adjudicate protocol compliance without the source records.",
  "Timestamp candidates in raw tables and SAE reporting deviations are exported as leads. A deviation associated with a deceased participant may relate to another event; it does not automatically establish a death-reporting delay.",
  "CRF matching uses exact P301 subject IDs in filenames. PHMPT_PDF_LINKS_FILE can map subject,path for other documents. All matched versions are retained; byte-identical PDFs share one unchanged copy. Original filenames, source paths and MD5 hashes are recorded. PDFs with no filename subject ID are listed for manual matching, not silently ignored.",
  "PDF page hits are search aids and may be blank form prompts. pdftools (or an installed pdftotext executable) text extraction is optional; image-only PDFs require separate OCR. Missing PDF/text does not imply absent reporting. Complete PDFs are copied without alteration, including potentially redacted/placeholder variants.",
  "Discovery uses configured submissions only; ambiguous snapshots are not concatenated. Same-stem readable XPT/SAS takes priority over CSV; exact copies are deduplicated. Conflicts are logged and the affected domain is withheld until explicitly selected with PHMPT_FILES.",
  "JMP ADSL reading is limited to the declared projection; CSV/XPT preserve all fields. Tables without a usable subject identifier or unexpected schema are logged as unreadable. Missing/ambiguous inputs lead to a partial report, not a claim of no deaths.",
  "V36: audit_assets/V36_death_totals_by_reporting_scope.png, 8.5 inches wide at 300 dpi by default. Published benchmarks and reconstructed BLA ADSL counts are distinct. CSV inputs, per-subject inclusion and source references accompany the PNG; the visual contains no reference markers or source URLs.",
  "Outputs: report.html; review_export/index.html, records/, documents/, death_evidence.csv, deaths_by_subject.csv, fatal_AE_comparison.csv, reporting_date_candidates.csv, review worksheet, timeliness results and CRF manifests. Optional ZIP contains the whole review_export folder.",
  "FDA EUA review: https://www.fda.gov/media/144434/download (pages 42-43). Sponsor briefing: https://www.fda.gov/media/144452/download (section 7.2.4.3).",
  "Released protocol collection: https://phmpt.org/wp-content/uploads/2025/04/125752_S3_M5_5351_mrna1273p301-16-1-1.pdf . Local Amendment 4 filename: 19745_S42_M5_protocol-mrna-1273-p301-amend-4.pdf.",
  "Requested records if unresolved: original SAE notifications, subsequent fatal-outcome updates, awareness/source records, EDC audit trail, follow-up correspondence and sponsor/FDA transmission receipts. These have different meanings and should be reconciled, not substituted for one another."
)
dr_report <- function(subjects,coverage,counts,assessment,crf,export,status,visual_html="") {
  focus<-subjects[subjects$review_focus,,drop=FALSE]
  asummary<-if(nrow(assessment))aggregate(list(notification_records=rep(1L,nrow(assessment))),assessment[c("arm","notification_type","assessment")],sum) else data.frame()
  links<-list.files(export,pattern="[.]csv$",full.names=FALSE)
  nav<-paste0("<ul>",paste0("<li><a href='",dr_h(links),"'>",dr_h(links),"</a></li>",collapse=""),"<li><a href='README.txt'>Methods and worksheet instructions</a></li></ul>")
  m<-crf$manifest;h<-crf$hits
  subjectblocks<-vapply(focus$subject,function(id){
    mm<-m[m$subject==id,,drop=FALSE];hh<-h[h$subject==id,,drop=FALSE];z<-focus[focus$subject==id,,drop=FALSE]
    doclinks<-if(nrow(mm))paste(vapply(seq_len(nrow(mm)),function(j){r<-mm[j,];lab<-paste(r$document_type,r$original_filename,"|",r$copy_status,"|",r$text_status);if(nzchar(r$export_path)&&r$copy_status=="COPIED_VERIFIED")paste0("<li><a href='",dr_h(r$export_path),"'>",dr_h(lab),"</a></li>") else paste0("<li>",dr_h(lab),"</li>")},""),collapse="") else "<li>No subject-linked PDF found. Check the unassigned CRF inventory and local release coverage.</li>"
    pagehtml<-if(nrow(hh))paste0("<details><summary>",nrow(hh)," keyword page hits (navigation only)</summary><ul>",paste(vapply(seq_len(nrow(hh)),function(j){r<-hh[j,];lab<-paste("PDF page",r$pdf_page,"—",r$keywords);p<-if(nzchar(r$export_path))paste0("<a href='",dr_h(r$export_path),"#page=",r$pdf_page,"'>",dr_h(lab),"</a>") else dr_h(lab);paste0("<li>",p,"<br><small>",dr_h(r$snippet),"</small></li>")},""),collapse=""),"</ul></details>") else ""
    paste0("<section class='subject' id='",id,"'><h3>",id," — ",dr_h(z$arm),"</h3><p>",dr_h(z$priority)," | Death date bounds: ",dr_h(z$death_date_lower)," to ",dr_h(z$death_date_upper),"</p><ul>",doclinks,"</ul>",pagehtml,"</section>")
  },"")
  summarycols<-intersect(c("subject","arm","EUA_death_evidence","BLA_death_evidence","death_date_lower","death_date_upper","priority","death_flag_N_present"),names(focus))
  html<-c("<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Moderna P301 — Death reporting review</title><style>body{max-width:1200px;margin:32px auto;padding:0 24px;font:16px/1.5 system-ui,sans-serif;color:#192432}h1,h2,h3{line-height:1.2}h2{margin-top:2em}table{border-collapse:collapse;font-size:13px;width:100%}th,td{border:1px solid #ccd3db;padding:7px;text-align:left;vertical-align:top;overflow-wrap:anywhere}th{background:#eef2f6}.scroll{overflow:auto}.notice{padding:14px 18px;border-left:4px solid #325a7b;background:#f0f5fa}.muted,small{color:#526171}.subject{border-top:1px solid #cbd5df;margin-top:26px;padding-top:8px}a{color:#145b8f}@media print{body{max-width:none;margin:0;padding:0;font-size:10pt}.scroll{overflow:visible}table{font-size:8pt}thead{display:table-header-group}tr{break-inside:avoid}a{color:inherit;text-decoration:none}.downloads,.references,details{display:none}h2,h3{break-after:avoid}}</style></head><body>",
    paste0("<h1>Moderna P301: death reporting at EUA</h1><p>Version ",DR_VERSION," · Generated ",format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z")," · ",dr_h(status),"</p>"),
    paste0("<div class='notice'><strong>Scope.</strong> EUA safety cutoff: ",PHMPT_EUA_CUTOFF,"; updated death count: ",PHMPT_EUA_UPDATE,"; authorization: ",PHMPT_REVIEW_END,". Later BLA records (",PHMPT_BLA_CUTOFF,") provide retrospective reconciliation. A death date, awareness date, first sponsor notification and database entry are distinct. Only verified awareness-to-notification intervals support a 24-hour assessment.</div>"),
    visual_html,
    "<h2>Input coverage</h2>",dr_table(coverage),
    "<p>Unavailable/conflicting domains are not treated as negative evidence. The discovery logs and source manifest identify the selected files and any alternatives withheld.</p>",
    "<h2>Deaths identified in the available records</h2>",dr_table(counts),
    "<p>Blank EUA comparison counts indicate insufficient input coverage, not zero deaths. Counts are unique participants with structured death evidence and a date wholly on or before the stated boundary. Undated/conflicting deaths are shown separately and repeated across windows; they cannot be allocated to a boundary. Evidence in an EUA file does not itself establish when that file was created or transmitted. No arm-specific mortality rate is calculated here.</p>",
    "<h2>Verified notification intervals</h2>",dr_table(asummary),
    if(!nrow(assessment)||!any(assessment$assessment%in%c("EXCEEDS_24H","WITHIN_24H_BOUND")))"<p><strong>No definitive 24-hour classifications are available from verified review entries in this run.</strong> Complete the persistent review worksheet using the original SAE/CRF records and rerun.</p>" else "<p>Counts refer to notification records, not necessarily distinct participants. INITIAL_SAE and DEATH_UPDATE have separate awareness dates. Unassessable entries are not counted as timely.</p>",
    dr_table(assessment[,intersect(c("subject","arm","notification_type","event_key","delay_lower_hours","delay_upper_hours","assessment","assessment_reason"),names(assessment)),drop=FALSE]),
    "<h2>Priority subject review</h2>",dr_table(focus[,summarycols,drop=FALSE]),
    "<p>A pre-cutoff death absent from the selected EUA tables requires reconciliation with the original SAE reports and submission history. A later fatal-outcome update can legitimately differ from an earlier nonfatal snapshot. Death-flag N alongside later death evidence also requires temporal context.</p>",
    paste0("<h2>CRF and subject documents</h2><p>",nrow(m)," subject/source document links; ",length(unique(m$md5[m$copy_status=="COPIED_VERIFIED"]))," unchanged PDFs copied and verified. ",crf$anonymous," CRF-named PDFs lack a subject ID in their filename and require manual matching. Missing text or documents is shown in the export manifest.</p>"),subjectblocks,
    "<div class='downloads'><h2>Review files</h2>",nav,"</div>",
    "<div class='references'><h2>Reference and interpretation</h2><p>Protocol Amendment 4, sections 8.3.6–8.3.10, sets the awareness-based 24-hour criterion. The <a href='https://www.fda.gov/media/144434/download'>FDA EUA review</a> reports 13 deaths (6 vaccine, 7 placebo) through December 3, a later boundary than the November 25 safety cutoff. This contextual total is not imposed on the raw-table reconstruction. See README for the released protocol reference and all decision rules.</p></div></body></html>")
  writeLines(enc2utf8(html),file.path(export,"index.html"),useBytes=TRUE)
  invisible(asummary)
}
analyse_death_reporting <- function(root=PHMPT_ROOT,output_dir=PHMPT_OUTPUT_DIR) {
  root<-dr_root(root)
  out<-if(nzchar(output_dir))dr_path(root,output_dir) else file.path(root,"analysis",paste0("death_reporting_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid()))
  if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop("Output folder is not empty; select a new PHMPT_OUTPUT_DIR to preserve provenance.")
  dir.create(out,recursive=TRUE,showWarnings=FALSE);export<-file.path(out,"review_export");dir.create(export)
  inventory<-list(EUA=list.files(file.path(root,"eua_data"),pattern="[.](csv|xpt|sas7bdat|jmp)$",recursive=TRUE,full.names=TRUE,ignore.case=TRUE),BLA=list.files(file.path(root,"prod_data"),pattern="[.](csv|xpt|sas7bdat|jmp)$",recursive=TRUE,full.names=TRUE,ignore.case=TRUE))
  inputs<-logs<-coverage<-list()
  for(rel in c("EUA","BLA"))for(dom in dr_domains){p<-dr_pick(root,rel,dom,out,inventory);logs[[length(logs)+1L]]<-p$log
    coverage[[length(coverage)+1L]]<-data.frame(release=rel,domain=dom,status=p$status,death_evidence_fields=if(is.null(p$input))"" else dr_join(intersect(names(p$input$data),switch(dom,ADSL=c("DTHFL","DTHDT","DTHDTC","EOSSTT"),DM=c("DTHFL","DTHDTC"),AE=c("AESDTH","AEOUT"),ADAE=c("AESDTH","AEOUT"),DS=c("DSDECOD","DSTERM"),c("QNAM","QVAL")))),n_rows=if(is.null(p$input))NA_integer_ else nrow(p$input$data),file=if(is.null(p$input))"" else p$input$path)
    if(!is.null(p$input))inputs[[paste(rel,dom,sep="_")]]<-p$input
  }
  coverage<-dr_bind(coverage);dr_csv(coverage,file.path(export,"input_coverage.csv"));dr_csv(dr_bind(logs),file.path(export,"input_discovery.csv"))
  if(!length(inputs))stop("No usable P301 inputs. Discovery logs saved in ",out,". Set PHMPT_FILES or resolve reader errors.")
  manifest<-dr_bind(lapply(inputs,function(a)data.frame(release=a$release,domain=a$domain,path=a$path,md5=a$md5,n_rows=nrow(a$data))))
  dr_csv(manifest,file.path(export,"source_manifest.csv"))
  evidence<-dr_bind(lapply(inputs,dr_signals))
  if(!nrow(evidence))evidence<-data.frame(release=character(),domain=character(),subject=character(),source_row=integer(),signal=character(),strength=character(),date_field=character(),date_raw=character(),date_lower=numeric(),date_upper=numeric(),date_precision=character(),event_id=character(),term=character(),ae_start=character(),ae_end=character(),source_file=character(),source_md5=character())
  subjects<-dr_subjects(evidence,inputs,coverage)
  if(!nrow(subjects))subjects<-data.frame(subject=character(),arm=character(),review_focus=logical(),structured_death_evidence=logical())
  # Stable row ordering makes manual reconciliation easier.
  if(nrow(subjects))subjects<-subjects[order(!subjects$review_focus,subjects$priority,subjects$subject),,drop=FALSE]
  dr_csv(subjects,file.path(export,"deaths_by_subject.csv"));dr_csv(evidence,file.path(export,"death_evidence.csv"))
  if(nrow(subjects))dr_csv(subjects[subjects$death_date_conflict,,drop=FALSE],file.path(export,"death_date_conflicts.csv"))
  ids<-subjects$subject[subjects$review_focus];fields<-dr_bind(lapply(inputs,dr_inventory_fields));dr_csv(fields,file.path(export,"field_inventory.csv"))
  candidates<-deviations<-list()
  for(a in inputs){x<-a$data[a$data$.subject%in%ids,,drop=FALSE];dr_csv(x,file.path(export,"records",paste0(a$release,"_",a$domain,".csv")))
    ff<-fields[fields$release==a$release&fields$domain==a$domain&fields$reporting_date_candidate,,drop=FALSE]
    for(n in intersect(ff$variable,names(x))){ii<-which(nzchar(dr_col(x,n)));if(length(ii))candidates[[length(candidates)+1L]]<-data.frame(release=a$release,domain=a$domain,subject=x$.subject[ii],source_row=x$.row[ii],variable=n,value=dr_col(x,n)[ii],source_path=a$path)}
    if(all(c("QNAM","QVAL")%in%names(x))){ii<-which(dr_col(x,"QNAM")%in%ff$variable&nzchar(dr_col(x,"QVAL")));if(length(ii))candidates[[length(candidates)+1L]]<-data.frame(release=a$release,domain=a$domain,subject=x$.subject[ii],source_row=x$.row[ii],variable=dr_col(x,"QNAM")[ii],value=dr_col(x,"QVAL")[ii],source_path=a$path)}
    if(a$domain=="DV") {txt<-apply(x,1,paste,collapse=" ");ii<-grepl("SAE|SERIOUS ADVERSE|AESI",txt,ignore.case=TRUE)&grepl("24|REPORT|NOTIF",txt,ignore.case=TRUE);q<-x[ii,,drop=FALSE];q$release<-rep(a$release,nrow(q));deviations[[length(deviations)+1L]]<-q}
  }
  dr_csv(dr_bind(candidates),file.path(export,"reporting_date_candidates.csv"));dr_csv(dr_bind(deviations),file.path(export,"subject_associated_reporting_deviations.csv"))
  changes<-dr_fatal_changes(inputs,ids);dr_csv(changes,file.path(export,"fatal_AE_comparison.csv"))
  template<-dr_review_template(subjects);dr_csv(template,file.path(export,"death_reporting_review_template.csv"))
  review_path<-dr_path(root,PHMPT_REVIEW_FILE)
  if(file.exists(review_path)) {
    review<-read.csv(review_path,colClasses="character",fileEncoding="UTF-8-BOM",check.names=FALSE,na.strings=character())
    if(!all(c("subject","notification_type","event_key")%in%names(review)))stop("Existing review worksheet lacks subject/notification_type/event_key: ",review_path)
    # Preserve existing reviews. New subject/type pairs are offered in a separate file,
    # avoiding changes to a reviewed worksheet or overriding event-specific keys.
    have<-paste(dr_sid(review$subject),dr_norm(review$notification_type));want<-paste(template$subject,template$notification_type)
    additions<-template[!want%in%have,,drop=FALSE];dr_csv(additions,file.path(export,"new_review_rows_to_add.csv"))
    review<-dr_bind(list(review,additions))
  }else {review<-template;dr_csv(template,review_path)}
  dr_csv(review,file.path(export,"review_input_used.csv"));assessment<-dr_assess(review,subjects);dr_csv(assessment,file.path(export,"notification_timeliness.csv"))
  counts<-dr_counts(subjects);dr_csv(counts,file.path(export,"death_counts_by_boundary.csv"))
  crf<-dr_crfs(root,subjects,export);dr_csv(crf$manifest,file.path(export,"crf_manifest.csv"));dr_csv(crf$hits,file.path(export,"crf_page_index.csv"))
  missing<-data.frame(subject=setdiff(ids,crf$manifest$subject[crf$manifest$document_type=="CRF"&crf$manifest$copy_status=="COPIED_VERIFIED"]))
  dr_csv(missing,file.path(export,"subjects_without_exported_crf.csv"))
  status<-if(all(coverage$status[coverage$domain%in%c("ADSL","DM","DS","AE")]=="AVAILABLE")&&all(nzchar(coverage$death_evidence_fields[coverage$domain%in%c("ADSL","DM","DS","AE")])))"CORE_INPUTS_AVAILABLE" else "PARTIAL_INPUT_COVERAGE"
  writeLines(dr_methods(),file.path(export,"README.txt"),useBytes=TRUE)
  params<-list(script_version=DR_VERSION,root=root,submissions=PHMPT_SUBMISSIONS,explicit_files=PHMPT_FILES,EUA_cutoff=PHMPT_EUA_CUTOFF,EUA_update=PHMPT_EUA_UPDATE,review_end=PHMPT_REVIEW_END,BLA_cutoff=PHMPT_BLA_CUTOFF,review_file=review_path,review_md5=if(file.exists(review_path))unname(tools::md5sum(review_path)) else "",pdf_roots=PHMPT_PDF_ROOTS,pdf_links_file=PHMPT_PDF_LINKS_FILE,extra_subjects=PHMPT_EXTRA_SUBJECTS,pdf_text=PHMPT_PDF_TEXT_SCAN,copy_pdfs=PHMPT_COPY_PDFS,status=status,audit_asset_dir=PHMPT_AUDIT_ASSET_DIR,audit_dpi=PHMPT_AUDIT_DPI)
  dput(params,file.path(export,"run_parameters.R"));writeLines(capture.output(sessionInfo()),file.path(export,"sessionInfo.txt"))
  visual<-dr36_export(root,out,export,inputs)
  asummary<-dr_report(subjects,coverage,counts,assessment,crf,export,status,visual$html);dr_csv(asummary,file.path(export,"timeliness_summary.csv"))
  writeLines("<!doctype html><meta charset='utf-8'><meta http-equiv='refresh' content='0; url=review_export/index.html'><a href='review_export/index.html'>Open the death-reporting review</a>",file.path(out,"report.html"))
  if(PHMPT_MAKE_ZIP) {
    if(requireNamespace("zip",quietly=TRUE)) {
      zpath<-file.path(out,"death_reporting_CRF_review.zip");old<-getwd();on.exit(setwd(old),add=TRUE);setwd(out)
      z<-tryCatch({zip::zipr(zpath,files="review_export",include_directories=TRUE);expected<-list.files(export,recursive=TRUE,full.names=FALSE);actual<-zip::zip_list(zpath)$filename;all(paste0("review_export/",expected)%in%actual)},error=function(e){message("ZIP failed: ",conditionMessage(e));FALSE});setwd(old)
      if(!z)warning("ZIP incomplete/unverified; use the complete review_export folder.")
    }else message("Optional package zip not installed; review_export is ready to copy as a folder.")
  }
  message("Report: ",file.path(export,"index.html"),"\nReview worksheet: ",review_path,"\nCRF export: ",export,"\nStatus: ",status)
  if(PHMPT_SHOW_REPORT&&interactive())browseURL(file.path(export,"index.html"))
  invisible(list(output=out,subjects=subjects,evidence=evidence,assessment=assessment,coverage=coverage,crf=crf,visual=visual))
}
# ---- Embedded V36 renderer: no shared R file or browser required ----
dr36_clean <- function(x) {
  x<-as.character(x);x[is.na(x)]<-"Unavailable"
  x<-gsub("https?://[^[:space:]<>]+","",x,perl=TRUE)
  x<-gsub("(?i)doi:[[:space:]]*[^[:space:]]+","",x,perl=TRUE)
  x<-gsub("\\[[0-9]+([,; -]+[0-9]+)*\\]","",x,perl=TRUE)
  trimws(x)
}
dr36_path <- function(root,p)if(grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)",p))p else file.path(root,p)
dr36_root <- function(root="",out=getwd()) {
  if(nzchar(root))return(normalizePath(root,winslash="/",mustWork=TRUE))
  for(start in unique(c(getwd(),out))) {
    p<-normalizePath(start,winslash="/",mustWork=TRUE)
    repeat {
      if(any(dir.exists(file.path(p,c("prod_data","eua_data")))))return(p)
      if(basename(p)=="analysis")return(dirname(p))
      up<-dirname(p);if(up==p)break;p<-up
    }
  }
  stop("Set PHMPT_ROOT or phmpt.root to the project folder for audit_assets exports.")
}
dr36_read <- function(out,file,required=TRUE) {
  p<-file.path(out,file)
  if(!file.exists(p)){if(required)stop("Missing export: ",p);return(data.frame())}
  read.csv(p,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
}
dr36_csv <- function(x,path) {
  con<-file(path,"wb");on.exit(close(con));writeBin(as.raw(c(239,187,191)),con)
  write.table(x,con,sep=",",row.names=FALSE,quote=TRUE,na="",qmethod="double")
}
dr36_num <- function(x)ifelse(is.na(x),"Unavailable",format(x,big.mark=",",scientific=FALSE,trim=TRUE))
dr36_rate <- function(k,n,digits=1)ifelse(is.na(k)|is.na(n),"Unavailable",ifelse(n>0,
  paste0(dr36_num(k),"/",dr36_num(n),sprintf(paste0(" (%.",digits,"f%%)"),100*k/n)),"0/0 (not estimable)"))
dr36_p <- function(p) {
  z<-ifelse(is.finite(p),formatC(p,digits=4,format="f"),"Not estimable")
  for(i in 2:6)z[is.finite(p)&p<10^-i]<-paste0("< ",format(10^-i,scientific=FALSE,trim=TRUE))
  z
}
dr36_text <- function(s,x,y,size=11.5,bold=FALSE,col="#25343D",just=c("left","top"))
  grid::grid.text(dr36_clean(s),x=grid::unit(x,"in"),y=grid::unit(1,"npc")-grid::unit(y,"in"),just=just,
    gp=grid::gpar(fontfamily="sans",fontsize=size,fontface=if(bold)"bold" else "plain",col=col,lineheight=1.12))
dr36_height <- function(s,size=11.5)grid::convertHeight(grid::grobHeight(grid::textGrob(s,
  gp=grid::gpar(fontfamily="sans",fontsize=size,lineheight=1.12))),"in",valueOnly=TRUE)
dr36_wrap <- function(s,width,size=11.5,bold=FALSE) {
  measure<-function(t)grid::convertWidth(grid::grobWidth(grid::textGrob(t,
    gp=grid::gpar(fontfamily="sans",fontsize=size,fontface=if(bold)"bold" else "plain"))),"in",valueOnly=TRUE)
  paragraphs<-strsplit(dr36_clean(s),"\n",fixed=TRUE)[[1]];all<-character()
  for(p in paragraphs) {
    words<-strsplit(gsub("[[:space:]]+"," ",trimws(p))," ",fixed=TRUE)[[1]];lines<-character();line<-""
    for(w in words) {
      z<-paste(c(if(nzchar(line))line,w),collapse=" ")
      if(measure(z)<=width){line<-z;next}
      if(nzchar(line)){lines<-c(lines,line);line<-""}
      while(nchar(w)>1L&&measure(w)>width) {
        n<-nchar(w)-1L;while(n>1L&&measure(substr(w,1,n))>width)n<-n-1L
        lines<-c(lines,substr(w,1,n));w<-substring(w,n+1L)
      };line<-w
    };all<-c(all,lines,line)
  };paste(all,collapse="\n")
}
dr36_begin <- function(root,out,owner,dpi=PHMPT_AUDIT_DPI,asset_dir=PHMPT_AUDIT_ASSET_DIR) {
  if(length(dpi)!=1||!is.finite(dpi)||dpi<150||dpi>600)stop("Audit DPI must be 150-600.")
  if(length(asset_dir)!=1||is.na(asset_dir)||!nzchar(asset_dir))stop("Audit asset directory is empty.")
  a<-new.env(parent=emptyenv());a$root<-dr36_root(root,out);a$out<-normalizePath(out,winslash="/",mustWork=TRUE)
  a$owner<-owner;a$dpi<-dpi;a$dest<-dr36_path(a$root,asset_dir)
  a$stage<-tempfile(paste0("audit_",owner,"_"),tmpdir=out);dir.create(a$stage)
  a$rows<-list();a$files<-character();a
}
dr36_png <- function(a,id,slug,draw,height=7,page=1L,pages=1L,status="AVAILABLE",rows=NA_integer_) {
  file<-paste0(id,"_",slug,if(page>1)sprintf("_%02d",page) else "",".png");path<-file.path(a$stage,file)
  grDevices::png(path,width=8.5,height=height,units="in",res=a$dpi,pointsize=13,bg="white",
    type=if(capabilities("cairo"))"cairo" else if(.Platform$OS.type=="windows")"windows" else getOption("bitmapType","Xlib"))
  device<-grDevices::dev.cur();ok<-FALSE
  on.exit({if(device%in%grDevices::dev.list())grDevices::dev.off(device);if(!ok)unlink(path)})
  draw();grDevices::dev.off(device)
  bytes<-readBin(path,"raw",n=file.info(path)$size)
  if(length(bytes)<32||!identical(bytes[1:8],as.raw(c(137,80,78,71,13,10,26,10)))||
    !identical(tail(bytes,12),as.raw(c(0,0,0,0,73,69,78,68,174,66,96,130))))stop("Incomplete PNG: ",path)
  dims<-vapply(list(bytes[17:20],bytes[21:24]),function(b)sum(as.integer(b)*256^(3:0)),numeric(1))
  if(any(abs(dims-round(c(8.5,height)*a$dpi))>1))stop("Invalid PNG dimensions.")
  a$rows[[length(a$rows)+1L]]<-data.frame(visual=id,file=file,page=page,pages=pages,status=status,rows=rows,width_px=dims[1],height_px=dims[2],dpi=a$dpi)
  a$files<-c(a$files,file);ok<-TRUE;invisible(path)
}
dr36_data <- function(a,id,x,suffix="data") {
  file<-paste0(id,"_",suffix,".csv");dr36_csv(x,file.path(a$stage,file));a$files<-c(a$files,file)
}
dr36_table <- function(a,id,slug,title,data,subtitle="",note="",widths=NULL,status="AVAILABLE",body_size=11.5,highlight=integer()) {
  data<-as.data.frame(data,check.names=FALSE);dr36_data(a,id,data)
  if(!ncol(data))data<-data.frame(Result="No matching rows",check.names=FALSE)
  if(!nrow(data)){data[1,]<-rep("No matching rows",ncol(data));status<-"NO_MATCHING_ROWS"}
  if(is.null(widths))widths<-if(ncol(data)==3)c(.5,.25,.25) else rep(1/ncol(data),ncol(data))
  widths<-widths/sum(widths)*7.8;left<-.35+c(0,head(cumsum(widths),-1))
  tmp<-tempfile(fileext=".pdf");if(capabilities("cairo"))grDevices::cairo_pdf(tmp,width=8.5,height=13,family="sans") else grDevices::pdf(tmp,width=8.5,height=13)
  dev<-grDevices::dev.cur();on.exit({if(dev%in%grDevices::dev.list())grDevices::dev.off(dev);unlink(tmp)})
  wraprow<-function(z,size=body_size,bold=FALSE)vapply(seq_along(z),function(j)dr36_wrap(z[j],widths[j]-.16,size,bold),character(1))
  heads<-wraprow(names(data),11.2,TRUE);head_h<-max(vapply(heads,dr36_height,numeric(1),size=11.2))+.22
  title<-dr36_wrap(paste(id,title),7.8,17,TRUE);sub<-dr36_wrap(subtitle,7.8,11.2);foot<-dr36_wrap(note,7.8,10.6)
  top<-.25+dr36_height(title,17)+.18+if(nzchar(sub))dr36_height(sub,11.2)+.18 else 0
  bottom<-.42+if(nzchar(foot))dr36_height(foot,10.6) else 0
  budget<-12.8-top-head_h-bottom
  cells<-list();heights<-numeric();origin<-integer()
  # Split a very long source row into labelled continuation rows instead of truncating it.
  for(i in seq_len(nrow(data))) {
    line<-lapply(wraprow(as.character(data[i,])),strsplit,split="\n",fixed=TRUE);line<-lapply(line,`[[`,1)
    cap<-max(1L,floor((min(budget,5.5)-.25)/(body_size/72*1.4)))
    nchunks<-max(1L,ceiling(max(lengths(line))/cap))
    for(k in seq_len(nchunks)) {
      row<-vapply(line,function(z)paste(z[seq_along(z)>=(k-1)*cap+1L & seq_along(z)<=k*cap],collapse="\n"),character(1))
      if(k>1&&!nzchar(row[1]))row[1]<-paste0(as.character(data[i,1])," (continued)")
      row[1]<-dr36_wrap(row[1],widths[1]-.16,body_size)
      h<-max(vapply(row,dr36_height,numeric(1),size=body_size))+.24
      if(h>budget)stop("Audit row is too tall even after wrapping: ",id)
      cells[[length(cells)+1L]]<-row;heights<-c(heights,h);origin<-c(origin,i)
    }
  }
  pages<-list();page<-integer();used<-0
  for(i in seq_along(cells)) {
    if(length(page)&&used+heights[i]>budget){pages[[length(pages)+1L]]<-page;page<-integer();used<-0}
    page<-c(page,i);used<-used+heights[i]
  };if(length(page))pages[[length(pages)+1L]]<-page
  grDevices::dev.off(dev)
  for(p in seq_along(pages)) {
    ix<-pages[[p]];height<-max(3.6,top+head_h+sum(heights[ix])+bottom)
    draw<-function() {
      grid::grid.newpage();dr36_text(title,.35,.25,17,TRUE,col="#17384B")
      sy<-.25+dr36_height(title,17)+.18;if(nzchar(sub))dr36_text(sub,.35,sy,11.2)
      y<-top
      for(r in c(0L,ix)) {
        h<-if(r==0)head_h else heights[r]
        grid::grid.rect(x=grid::unit(.35,"in"),y=grid::unit(1,"npc")-grid::unit(y,"in"),
          width=grid::unit(7.8,"in"),height=grid::unit(h,"in"),just=c("left","top"),
          gp=grid::gpar(fill=if(r==0)"#E8EFF3" else if(origin[r]%in%highlight)"#FFF0D9" else if(origin[r]%%2==0)"#F3F6F8" else "white",col=NA))
        values<-if(r==0)heads else cells[[r]]
        for(j in seq_along(values))dr36_text(values[j],left[j]+.08,y+.11,if(r==0)11.2 else body_size,r==0||(r>0&&j==1&&origin[r]%in%highlight),
          col=if(r!=0&&grepl("mRNA|Vaccine",names(data)[j],ignore.case=TRUE))"#1E667D" else if(r!=0&&grepl("Placebo",names(data)[j]))"#A65B22" else "#25343D")
        y<-y+h
      }
      if(nzchar(foot))dr36_text(foot,.35,y+.16,10.6,col="#52616C")
      dr36_text(paste0("Page ",p," of ",length(pages)),8.15,height-.17,9.5,col="#64737D",just=c("right","bottom"))
    }
    dr36_png(a,id,slug,draw,height,p,length(pages),status,length(unique(origin[ix])))
  };invisible(NULL)
}
dr36_unavailable <- function(a,id,slug,title,reason)
  dr36_table(a,id,slug,title,data.frame(Status="Not available in this run",Reason=reason,check.names=FALSE),widths=c(.3,.7),status="UNAVAILABLE")
dr36_finish <- function(a,notes=character()) {
  if(!length(a$rows))stop("No audit assets were rendered.")
  manifest<-do.call(rbind,a$rows);manifest$renderer_version<-"1.0.0";manifest$generated_at<-format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z")
  manifest$run_directory<-a$out;manifest$md5<-unname(tools::md5sum(file.path(a$stage,manifest$file)))
  mf<-paste0(a$owner,"_manifest.csv");dr36_csv(manifest,file.path(a$stage,mf))
  readme<-paste0(a$owner,"_README.txt")
  writeLines(c("AUDIT PNG ASSETS", "8.5 inches wide; 300 dpi by default (2550 pixels). White background. Upload PNG originals to Substack.",
    "Tables wrap and continue as _02, _03, etc.; no rows are silently omitted. Printed figures have no bibliography, footnote markers, DOI or source URLs.",
    "Statistical definitions, endpoint labels, observation windows and denominators remain on figures. Source references stay in the main report and original provenance/method exports.",
    "A status of UNAVAILABLE is an explicit missing-input/model notice, not a zero finding. Data CSVs preserve the displayed inputs.",notes),file.path(a$stage,readme))
  a$files<-unique(c(a$files,mf,readme));destinations<-unique(c(file.path(a$out,"audit_assets"),a$dest))
  # Complete rendering precedes publishing. Only files listed in this owner's
  # previous manifest and this batch may be replaced or removed.
  for(dest in destinations) {
    dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    oldmf<-file.path(dest,mf);old<-if(file.exists(oldmf))dr36_read(dest,mf)$file else character()
    owned_old<-old[grepl("^V[0-9]{2}_[A-Za-z0-9_]+[.]png$",old)]
    existing<-intersect(unique(c(a$files,owned_old)),list.files(dest))
    backup<-tempfile("backup_",tmpdir=a$stage);dir.create(backup)
    if(length(existing)&&!all(file.copy(file.path(dest,existing),backup,overwrite=TRUE)))stop("Cannot back up previous audit assets.")
    success<-FALSE
    tryCatch({
      for(f in a$files) {
        if(!file.copy(file.path(a$stage,f),file.path(dest,f),overwrite=TRUE))stop("Cannot publish ",f)
        if(unname(tools::md5sum(file.path(a$stage,f)))!=unname(tools::md5sum(file.path(dest,f))))stop("Audit copy mismatch: ",f)
      }
      stale<-setdiff(owned_old,a$files)
      if(length(stale)&&unlink(file.path(dest,stale))!=0L)stop("Cannot remove obsolete continuation pages.")
      success<-TRUE
    },finally={if(!success){unlink(file.path(dest,setdiff(a$files,existing)));file.copy(file.path(backup,existing),dest,overwrite=TRUE)}})
    unlink(backup,recursive=TRUE)
  }
  unlink(a$stage,recursive=TRUE)
  message("Audit PNGs: ",a$dest," | ",paste(unique(manifest$visual),collapse=", "))
  invisible(list(directory=a$dest,manifest=manifest))
}


# V36 keeps published benchmarks separate from locally reconstructed counts.
# No manually reviewed CRF finding is silently treated as an automated result.
dr36_references <- function() data.frame(
  scope_id=c("EUA_PAPER_28D","EUA_FDA_DEC03","BLA_BLINDED","BLA_EXTENDED"),
  scope=c("Initial EUA publication\n28-day fatal-AE window",
    "FDA EUA review\nReports through 3 Dec 2020",
    "Blinded-phase publication\nBlinded follow-up",
    "FDA BLA review\nFollow-up through 4 May 2021"),
  vaccine=c(2L,6L,16L,24L),placebo=c(3L,7L,16L,20L),total=c(5L,13L,32L,44L),
  source=c("Baden et al.; safety results and supplementary Table S8, p. 33",
    "FDA EUA final review memorandum, 18 December 2020, printed p. 44",
    "El Sahly et al.; safety results and Tables S19/S26",
    "FDA BLA clinical review, 28 January 2022, pp. 36, 60, 79-81; Table 37"),
  url=c("https://www.nejm.org/doi/full/10.1056/NEJMoa2035389",
    "https://www.fda.gov/media/144673/download",
    "https://www.nejm.org/doi/full/10.1056/NEJMoa2113017",
    "https://www.fda.gov/media/156342/download"),
  definition=c("Fatal AEs in the 28-day post-injection window; not all deaths known at publication.",
    "Death reports through 3 December; not a count based only on retrospective death dates.",
    "Blinded follow-up, with the 26 March 2021 analysis boundary.",
    "Original randomized arms: 24 vaccine, 17 placebo plus 3 after crossover = 20 original placebo. Extended death surveillance through database lock."),
  stringsAsFactors=FALSE)

dr36_reconstruct <- function(inputs) {
  refs<-dr36_references()[3:4,]
  q<-data.frame(scope_id=refs$scope_id,vaccine=NA_integer_,placebo=NA_integer_,total=NA_integer_,
    reference_vaccine=refs$vaccine,reference_placebo=refs$placebo,reference_total=refs$total,
    status="NOT_RECONSTRUCTED",reason="No usable BLA ADSL selected.",
    adsl_rows=NA_integer_,source_path="",source_md5="",stringsAsFactors=FALSE)
  participants<-data.frame(subject=character(),arm=character(),safety_flag=character(),
    death_date_raw=character(),death_date=character(),unblinding_raw=character(),
    death_date_issue=character(),blinding_issue=character(),
    included_blinded=logical(),included_extended=logical(),stringsAsFactors=FALSE)
  result<-function()list(summary=q,subjects=participants)
  a<-inputs[["BLA_ADSL"]];if(is.null(a))return(result())
  x<-a$data;q$adsl_rows<-nrow(x);q$source_path<-a$path;q$source_md5<-a$md5
  reject<-function(reason){q$reason<<-reason;result()}
  if(as.character(PHMPT_BLA_CUTOFF)!="2021-03-26")return(reject("Configured BLA cutoff differs from the published 26 March 2021 analysis."))
  # The selected complete S3 ADSL has 30,423 unique IDs. A review-export subset
  # can contain every early death yet cannot establish the full death count.
  if(nrow(x)!=30423L||anyDuplicated(x$.subject)||any(!nzchar(x$.subject)))
    return(reject("Full S3 ADSL roster required: 30,423 distinct subject IDs. A subset is not a full reconstruction."))
  required<-c("TRT01P","SAFFL")
  if(!all(required%in%names(x))||!any(c("DTHDT","DTHDTC")%in%names(x)))
    return(reject("ADSL needs original planned arm, safety flag and a death-date column."))
  arm<-dr_arm(dr_col(x,"TRT01P"));safe<-dr_norm(dr_col(x,"SAFFL"))
  ns<-c(sum(safe=="Y"&arm=="mRNA-1273"),sum(safe=="Y"&arm=="Placebo"))
  if(sum(ns)!=30346L||any(safe=="Y"&arm=="UNASSIGNED"))
    return(reject("S3 safety roster does not reconcile to 30,346 participants with identifiable original arms."))
  fields<-intersect(c("DTHDT","DTHDTC"),names(x))
  signal<-dr_norm(dr_col(x,"DTHFL"))%in%c("Y","YES")
  for(f in fields)signal<-signal|nzchar(dr_col(x,f))
  z<-x[signal,,drop=FALSE];n<-nrow(z)
  dd<-rep(NA_real_,n);issues<-rep("",n)
  for(i in seq_len(n)) {
    values<-vapply(fields,function(f)dr_col(z,f)[i],character(1));values<-values[nzchar(values)]
    if("DTHDT"%in%names(values)&&nzchar(dr_col(z,"DTHDTF")[i]))values<-values[names(values)!="DTHDT"]
    bd<-dr_day_bounds(values)
    if(!length(values))issues[i]<-"Missing or imputed death date"
    else if(any(bd$precision!="DAY"))issues[i]<-"Partial or invalid death date"
    else if(length(unique(bd$lower))!=1L)issues[i]<-"Conflicting death dates"
    else dd[i]<-bd$lower[1]
    if(dr_norm(dr_col(z,"DTHFL"))[i]%in%c("N","NO"))issues[i]<-"Death date conflicts with negative death flag"
  }
  sz<-dr_norm(dr_col(z,"SAFFL"));az<-dr_arm(dr_col(z,"TRT01P"))
  rawub<-dr_col(z,"UBLDDT");ub<-dr_day_bounds(rawub)
  ub_issue<-ifelse(nzchar(rawub)&ub$precision!="DAY","Partial or invalid unblinding date","")
  if(!"UBLDDT"%in%names(x))ub_issue[]<-"Unblinding column unavailable"
  eligible<-sz=="Y"
  known<-!nzchar(issues)&!is.na(dd)
  inc_extended<-eligible&known&dd<=as.numeric(as.Date("2021-05-04"))
  inc_blinded<-eligible&known&dd<=as.numeric(as.Date("2021-03-26"))&
    !nzchar(ub_issue)&(!nzchar(rawub)|dd<=ub$lower)
  inc_extended[is.na(inc_extended)]<-FALSE;inc_blinded[is.na(inc_blinded)]<-FALSE
  participants<-data.frame(subject=z$.subject,arm=az,safety_flag=sz,
    death_date_raw=vapply(seq_len(n),function(i)paste(vapply(fields,function(f)paste0(f,"=",dr_col(z,f)[i]),character(1)),collapse="; "),character(1)),
    death_date=dr_iso(dd),unblinding_raw=rawub,death_date_issue=issues,blinding_issue=ub_issue,
    included_blinded=inc_blinded,included_extended=inc_extended,stringsAsFactors=FALSE)
  unknown_safety<-any(!sz%in%c("Y","N"))
  for(j in 1:2) {
    bad<-sum(eligible&nzchar(issues))
    if(j==1L)bad<-bad+sum(eligible&known&dd<=as.numeric(as.Date("2021-03-26"))&nzchar(ub_issue),na.rm=TRUE)
    if(bad||unknown_safety|| (j==1L&&!"UBLDDT"%in%names(x))) {
      q$reason[j]<-"Unresolved death dates, safety membership or required unblinding data; counts withheld."
      next
    }
    take<-if(j==1L)inc_blinded else inc_extended
    q$vaccine[j]<-sum(take&az=="mRNA-1273");q$placebo[j]<-sum(take&az=="Placebo");q$total[j]<-sum(take)
    matches<-identical(as.integer(q[j,c("vaccine","placebo","total")]),as.integer(refs[j,c("vaccine","placebo","total")]))
    q$status[j]<-if(matches)"MATCH" else "DIFFERS"
    q$reason[j]<-if(matches)"Recorded-date counts agree in both original arms and total." else "Computed counts differ from the published reference; inspect V36_subjects.csv."
  }
  result()
}

dr36_export <- function(root,out,export,inputs) {
  a<-dr36_begin(root,out,"V36_deaths");refs<-dr36_references();r<-dr36_reconstruct(inputs);q<-r$summary
  display<-data.frame(`Reporting scope`=refs$scope,Origin="Published",
    `mRNA-1273`=as.character(refs$vaccine),Placebo=as.character(refs$placebo),Total=as.character(refs$total),check.names=FALSE)
  labels<-c("BLA ADSL: blinded reconstruction","BLA ADSL: extended reconstruction")
  status_label<-ifelse(q$status=="MATCH","Matches published counts",ifelse(q$status=="DIFFERS","Differs from published counts","Not reconstructed"))
  shown<-function(v)ifelse(is.na(v),"Unavailable",as.character(v))
  display<-rbind(display,data.frame(`Reporting scope`=paste(labels,status_label,sep="\n"),Origin="Computed",
    `mRNA-1273`=shown(q$vaccine),Placebo=shown(q$placebo),Total=shown(q$total),check.names=FALSE))
  status<-if(any(q$status=="DIFFERS"))"RECONCILIATION_DIFFERENCE" else if(all(q$status=="MATCH"))"RECONCILED" else if(any(q$status=="MATCH"))"PARTIAL_RECONSTRUCTION" else "REFERENCE_ONLY"
  dr36_data(a,"V36",refs,"references")
  dr36_data(a,"V36",q,"reconciliation")
  dr36_data(a,"V36",r$subjects,"subjects")
  dr36_table(a,"V36","death_totals_by_reporting_scope","Deaths by reporting scope",display,
    subtitle="Published totals and reconstruction from released BLA ADSL",
    note="Original randomized arms. Extended follow-up includes crossover. Published rows are reference figures; only the ADSL rows are computed. Death dates do not establish when deaths became known.",
    widths=c(.46,.12,.15,.15,.12),status=status,highlight=5:6)
  saved<-dr36_finish(a,c(
    "V36 shows published reference totals separately from computed BLA ADSL counts. No prior manual CRF finding is hard-coded as a new result.",
    "V36_data.csv contains the displayed table; V36_references.csv carries document links and definitions; V36_reconciliation.csv gives counts/status/reasons and source MD5; V36_subjects.csv contains all death candidates in the selected BLA ADSL and inclusion flags.",
    "Completeness gates: S3 ADSL must contain 30,423 unique P301 IDs and 30,346 safety participants with identifiable original arms. Actual-treatment safety-arm totals are not imposed on the original randomized arms. These are roster checks, not validation of every medical observation. Subsets and incompatible snapshots are not reported as complete reconciliations.",
    "Reconstruction uses SAFFL=Y and original TRT01P. Blinded: recorded death on/before 26 March 2021 and on/before recorded UBLDDT, or UBLDDT blank. Extended: recorded death on/before 4 May 2021. A same-day death/unblinding sequence cannot be resolved from date-only data.",
    "Missing, partial, conflicting or explicitly imputed death dates prevent a full reconstruction unless an exact non-imputed DTHDTC is available. A missing UBLDDT column prevents the blinded comparison; a blank value is treated as no recorded unblinding.",
    "This reproduces recorded-date totals, not the actual date or reporting timeliness of each death. In particular, the CRF-identified uncertainty in US3872496 requires medical review even if its final ADSL date is precise and unflagged.",
    "EUA reference counts are not reconstructed from later death dates: awareness and reporting-cutoff information is needed. The existing verified-review worksheet remains the basis for any 24-hour notification assessment.",
    "Full reconstruction is independent of availability of AE/DM/DS in the wider death audit. A partial-domain audit can still reconcile the complete ADSL counts; that does not certify reporting compliance.",
    paste(q$scope_id,q$status,q$reason,sep=" | ")))
  # The self-contained review ZIP includes V36 and its provenance files.
  target<-file.path(export,"audit_assets");dir.create(target,recursive=TRUE,showWarnings=FALSE)
  paths<-list.files(file.path(out,"audit_assets"),pattern="^V36_",full.names=TRUE)
  if(!all(file.copy(paths,target,overwrite=TRUE)))stop("Cannot copy V36 into the review export.")
  for(p in paths)if(unname(tools::md5sum(p))!=unname(tools::md5sum(file.path(target,basename(p)))))stop("V36 review-export copy mismatch.")
  block<-paste0("<h2>V36: death totals by reporting scope</h2><p>Published reference figures and local ADSL reconstruction are labelled separately. Current status: <strong>",status,
    "</strong>. This table does not establish notification timeliness.</p><a href='audit_assets/V36_death_totals_by_reporting_scope.png'><img src='audit_assets/V36_death_totals_by_reporting_scope.png' alt='V36: death totals by reporting scope' style='width:100%;max-width:1000px;height:auto'></a>",
    dr_table(q[,c("scope_id","vaccine","placebo","total","status","reason")]),
    "<p><a href='audit_assets/V36_data.csv'>Displayed table</a> | <a href='audit_assets/V36_reconciliation.csv'>Reconciliation</a> | <a href='audit_assets/V36_subjects.csv'>Subject-level inclusion</a> | <a href='audit_assets/V36_references.csv'>Sources and definitions</a></p>")
  list(directory=saved$directory,manifest=saved$manifest,summary=q,html=block)
}

# ---- RUN: select all and execute, or source this file in RStudio ----
if(isTRUE(getOption("phmpt.autorun",TRUE)))analyse_death_reporting()
