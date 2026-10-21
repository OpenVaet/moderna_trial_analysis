# =============================================================================
# Moderna P301 | Local / central PCR testing by arm and site | v1.1.1
# 2026-10-21 | Standalone R >= 4.3: select all and run in RStudio.
# Reads local PHMPT eua_data / prod_data; never downloads trial data or installs.
# CSV and supported PHMPT JMP: base R. XPT/SAS7BDAT: install.packages('haven').
# JMP supports compressed production ADSL and uncompressed EUA ADSL, including
# variable-length character columns, fixed-width characters and numeric dates.
# No JMP application, Julia installation, or manual CSV conversion is required.
# PRIMARY: unique subjects with COVID-form symptoms, and the union of subjects
# with COVID-form symptoms or matched symptom-like AEs. Linked PCR within 0-7 days
# of positive assessment / AE onset; cohort-level results remain secondary.
# LOCAL_SOURCE = MBGRPID 003 or LOCALFL=Y (COVID Diagnostic Test reporting route).
# CENTRAL_SOURCE = MBGRPID 001 + Viracor name + no contradictory local evidence.
# These are source-stream classifications, NOT proven physical lab locations.
# Blank LOCALFL alone is NEVER central; Viracor in group 003 remains local-source.
# Distinct MBSEQ rows are records, not necessarily distinct physical assays.
# Load functions only: options(phmpt.autorun=FALSE) before source().
# =============================================================================
PHMPT_ROOT <- getOption('phmpt.root', '')
PHMPT_OUTPUT_DIR <- getOption('phmpt.output_dir', '')
PHMPT_RELEASES <- c('EUA','Prod')
PHMPT_FILES <- list(EUA=list(ADSL='',MB='',ADMB='',SUPPMB='',AE='',ADSYMP='',EVENTS=''),
                    Prod=list(ADSL='',MB='',ADMB='',SUPPMB='',AE='',ADSYMP='',EVENTS=''))
# Optional explicit paths above override automatic discovery; paths can be relative
# to the project root or absolute. The complete prior mb_pcr_metadata.csv,
# admb_pcr_metadata.csv and suppmb_pcr_metadata.csv exports are also supported,
# but require a complete ADSL: a roster inferred from tested subjects is NOT used.
PHMPT_EUA_CUTOFF <- '2020-11-25'
PHMPT_BLA_CUTOFF <- '2021-03-26'
PHMPT_SCOPES <- c('STRICT_BLINDED','CALENDAR')
PHMPT_MIN_RANDOMIZED <- 30000L # safeguard against partial ADSL; lower only for tests
PHMPT_CSV_ENCODING <- 'auto' # UTF-8 first; Windows-1252 fallback
PHMPT_NUMERIC_DATE_ORIGIN <- '' # numeric CSV dates require an explicit origin
PHMPT_SHOW_REPORT <- interactive()
PHMPT_CENTRAL_GROUPS <- '001'
PHMPT_LOCAL_GROUPS <- '003'
PHMPT_CENTRAL_NAME_REGEX <- 'VIRACOR|VIRA[ -]?COR'


PHMPT_ANALYSE_SYMPTOMS <- TRUE
PHMPT_MAPPING_FILE <- '' # existing root/assets AE mapping; otherwise embedded agreed mapping
PHMPT_TEST_LINK_DAYS <- 7L
PHMPT_SYMPTOM_COHORTS <- c('COVID_FORM_ONLY','COVID_FORM_PLUS_AE')
PHMPT_SYMPTOM_FOLLOWUP <- c('BY_CUTOFF','COMPLETE_7D') # COMPLETE_7D means the configured full interval
PHMPT_SYMPTOM_TESTS <- c('DIAGNOSTIC','ANY_NONBASELINE')
# AE + ADSYMP are required for the two-population analysis. EVENTS is an optional
# explicit path to a COMPLETE prior symptom_event_evidence.csv export from the
# same release; default is to reconstruct directly from local source tables.
# BY_CUTOFF counts every subject with usable dated evidence; COMPLETE_7D omits
# evidence dates lacking a full follow-up interval, irrespective of PCR outcome.
# Negative or missing COVID checklist assessments are never positive symptoms.

ws_norm <- function(x) toupper(trimws(gsub("[[:space:]\\p{Z}]+"," ",ws_text(x),perl=TRUE)))
ws_col <- function(x,n) if(n %in% names(x))ws_text(x[[n]]) else rep("",nrow(x))
ws_week <- function(day) floor((as.numeric(day)+3)/7)*7-3
ws_iso <- function(x) {z<-format(as.Date(x,origin="1970-01-01"),"%Y-%m-%d");z[is.na(z)]<-"";z}
ws_n <- function(x)format(x,big.mark=",",scientific=FALSE,trim=TRUE)
ws_sid <- function(x) {
  z<-ws_norm(x);z<-sub("^.*(US[0-9]{3}-?[0-9]{4})$","\\1",z)
  z<-gsub("-","",z,fixed=TRUE)
  if(any(!grepl("^US[0-9]{7}$",z)))stop("Invalid P301 subject ID: ",paste(head(unique(z[!grepl("^US[0-9]{7}$",z)])),collapse=", "))
  z
}
ws_dates <- function(x,field,origin="",strict=FALSE) {
  s<-trimws(ws_text(x));s[s %in% c("NA",".")]<-"";out<-rep(NA_real_,length(s))
  iso<-grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}($|T| )",s)
  v<-substr(s[iso],1,10);d<-suppressWarnings(as.Date(v,"%Y-%m-%d"))
  d[!is.na(d)&format(d,"%Y-%m-%d")!=v]<-as.Date(NA);out[iso]<-as.numeric(d)
  sas<-grepl("^[0-9]{1,2}[A-Za-z]{3}[0-9]{4}$",s)
  if(any(sas)) {
    a<-toupper(s[sas]);dd<-as.integer(sub("[A-Z].*$","",a));mm<-match(sub("^[0-9]+([A-Z]{3}).*$","\\1",a),toupper(month.abb))
    yy<-as.integer(substr(a,nchar(a)-3,nchar(a)));v<-sprintf("%04d-%02d-%02d",yy,mm,dd)
    d<-suppressWarnings(as.Date(v,"%Y-%m-%d"));d[!is.na(d)&format(d,"%Y-%m-%d")!=v]<-as.Date(NA);out[sas]<-as.numeric(d)
  }
  num<-grepl("^-?[0-9]+([.][0-9]+)?$",s)&nzchar(s)
  if(any(num)) {
    if(!nzchar(origin))stop(field,": unformatted numeric dates. Set PHMPT_NUMERIC_DATE_ORIGIN from metadata, or use the original XPT/JMP.")
    v<-as.numeric(s[num]);if(any(v!=floor(v)))stop(field,": fractional numeric dates are unsupported.")
    out[num]<-as.numeric(as.Date(v,origin=origin))
  }
  if(strict&&any(nzchar(s)&!is.finite(out)))stop(field,": invalid/partial dates: ",paste(head(unique(s[nzchar(s)&!is.finite(out)]),4),collapse=", "))
  out
}
ws_text <- function(x) {
  if(inherits(x,"Date")) x<-format(x,"%Y-%m-%d")
  else if(inherits(x,"POSIXt")) x<-format(x,"%Y-%m-%d",tz="UTC")
  x<-as.character(x);x[is.na(x)]<-""
  bad<-!validUTF8(x)
  if(any(bad))x[bad]<-iconv(x[bad],from="windows-1252",to="UTF-8",sub=NA)
  if(anyNA(x))stop("Undecodable input text. Set PHMPT_CSV_ENCODING.")
  enc2utf8(x)
}
ws_csv <- function(x,path) {
  con<-file(path,"wb");on.exit(close(con));writeBin(as.raw(c(239,187,191)),con)
  write.table(x,con,sep=",",row.names=FALSE,quote=TRUE,na="",qmethod="double")
}
ws_path <- function(root,p) if(grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)",p))p else file.path(root,p)
ws_root <- function(root="") {
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






pt_specs <- function(domain) switch(domain,
 ADSL=list(required=c('STUDYID','USUBJID','RANDFL','TRT01P','RANDDT','CUTOFFDT'),optional=c('SUBJID','SITEID','TR01SDT','DOSE2DT','AP01EDT','UBLDDT','PDVDT','AP02SDT','TR02SDT','EUAVACDT','EOSDT','DTHDT')),
 MB=list(required=c('STUDYID','USUBJID','MBSEQ','MBTESTCD','MBSTRESC','MBDTC'),optional=c('MBORRES','MBSTAT','MBREASND','MBGRPID','MBREFID','MBNAM','MBSPEC','MBMETHOD','MBTEST','VISIT','VISITNUM')),
 ADMB=list(required=c('STUDYID','USUBJID','MBSEQ','PARAMCD','AVALC'),optional=c('DTYPE','MBDTC','ADT','VISIT','AVISIT','ABLFL','LOCALFL','MBSTAT','MBMETHOD')),
 AE=list(required=c('STUDYID','USUBJID','AESEQ','AETERM','AEDECOD','AESTDTC'),optional=c('AELLT','AEENDTC','AEENRF','AEOUT','AESER','AEREL')),
 ADSYMP=list(required=c('STUDYID','USUBJID','PARAMCD','AVALC','ADT'),optional=c('ASEQ','PARAM','ASTDT','AENDT','SRCDOM','SRCSEQ','SRCVAR','LNKID','SEV','CUTOFFDT')),
 SUPPMB=list(required=c('STUDYID','USUBJID','IDVAR','IDVARVAL','QNAM','QVAL'),optional=c('QLABEL','QORIG','QEVAL')),
 stop('Unsupported domain'))
pt_read <- function(path,domain) {
 spec<-pt_specs(domain);ext<-tolower(tools::file_ext(path));wanted<-c(spec$required,spec$optional)
 if(ext=='jmp')x<-ws_read_jmp(path,spec$required,spec$optional)
 else if(ext=='csv') {
  x<-NULL;errs<-character()
  for(enc in if(PHMPT_CSV_ENCODING=='auto')c('UTF-8-BOM','windows-1252') else PHMPT_CSV_ENCODING) {
   z<-tryCatch(withCallingHandlers({
    h<-read.csv(path,nrows=0,check.names=FALSE,fileEncoding=enc)
    nm<-toupper(sub('^\ufeff','',names(h)));if(anyDuplicated(nm))stop('Duplicate columns')
    read.csv(path,colClasses=ifelse(nm%in%wanted,'character','NULL'),check.names=FALSE,
             fileEncoding=enc,na.strings=c('','NA'),comment.char='')
   },warning=function(w)stop(conditionMessage(w))),error=function(e)e)
   if(!inherits(z,'error')){x<-z;break};errs<-c(errs,conditionMessage(z))
  }
  if(is.null(x))stop('CSV read failed: ',paste(errs,collapse='; '))
 } else if(ext%in%c('xpt','sas7bdat')) {
  if(!requireNamespace('haven',quietly=TRUE))stop('XPT/SAS7BDAT needs haven: install.packages("haven")')
  fun<-if(ext=='xpt')haven::read_xpt else haven::read_sas
  h<-fun(path,n_max=0,.name_repair='minimal');pick<-names(h)[toupper(names(h))%in%wanted]
  x<-as.data.frame(fun(path,col_select=tidyselect::all_of(pick),.name_repair='minimal'))
 } else stop('Unsupported file type: ',ext)
 names(x)<-toupper(sub('^\ufeff','',names(x)))
 if(anyDuplicated(names(x)))stop('Duplicate columns')
 miss<-setdiff(spec$required,names(x));if(length(miss))stop('Missing fields: ',paste(miss,collapse=', '))
 x[]<-lapply(x,ws_text);x$.source_row<-seq_len(nrow(x))
 x<-x[ws_norm(x$STUDYID)=='MRNA-1273-P301',,drop=FALSE]
 if(!nrow(x))stop('No P301 rows')
 for(n in setdiff(spec$optional,names(x)))x[[n]]<-rep('',nrow(x))
 x$subject<-ws_sid(x$USUBJID)
 for(n in intersect(c('MBSEQ','AESEQ','IDVARVAL','MBGRPID'),names(x))) {
  z<-trimws(x[[n]]);ok<-grepl('^[0-9]+([.]0+)?$',z)
  z[ok]<-as.character(as.numeric(z[ok]));if(n=='MBGRPID')z[ok]<-sprintf('%03d',as.integer(z[ok]));x[[n]]<-z
 }
 if(domain=='AE'&&any(!nzchar(x$AESEQ)))stop('Missing AESEQ')
 x
}
pt_submission <- function(p) {
 z<-tolower(basename(p));m<-regexpr('[0-9]+_s[0-9]+_m[0-9]+_mrna-1273-p301',z)
 v<-rep('',length(p));ok<-m>0;v[ok]<-substring(z[ok],m[ok],m[ok]+attr(m,'match.length')[ok]-1L);v
}
pt_candidates <- function(root,release,domain) {
 folder<-file.path(root,if(release=='Prod')'prod_data' else 'eua_data')
 p<-list.files(folder,recursive=TRUE,full.names=TRUE,pattern='\\.(csv|xpt|sas7bdat|jmp)$',ignore.case=TRUE)
 p<-p[grepl(paste0('(^|[_ -])',tolower(domain),'([_. -]|$)'),basename(p),ignore.case=TRUE)]
 p<-p[!grepl('p101|p201|dmid',basename(p),ignore.case=TRUE)]
 if(release=='Prod'&&any(grepl('_s3_m[0-9]+_',basename(p),ignore.case=TRUE)))p<-p[grepl('_s3_m[0-9]+_',basename(p),ignore.case=TRUE)]
 unique(p[order(match(tolower(tools::file_ext(p)),c('csv','jmp','xpt','sas7bdat')),p)])
}
pt_signature <- function(x) {
 x$.source_row<-NULL;x<-unique(x);x<-x[do.call(order,x),,drop=FALSE];rownames(x)<-NULL
 f<-tempfile();on.exit(unlink(f));saveRDS(x,f,compress=FALSE,version=2);unname(tools::md5sum(f))
}
pt_pick <- function(root,release,domain,explicit,out,submission='',optional=FALSE) {
 paths<-if(nzchar(explicit))ws_path(root,explicit) else pt_candidates(root,release,domain)
 if(!length(paths)){if(optional)return(NULL);stop('No ',release,' ',domain,' found; set PHMPT_FILES$',release,'$',domain)}
 logs<-valid<-list();seen<-character();target<-if(release=='Prod')PHMPT_BLA_CUTOFF else PHMPT_EUA_CUTOFF
 for(i in seq_along(paths)) {
  p<-paths[i];h<-if(file.exists(p))unname(tools::md5sum(p)) else ''
  other<-!nzchar(explicit)&&nzchar(submission)&&pt_submission(p)!=submission
  dup<-nzchar(h)&&h%in%seen
  if(!other&&!dup)message('Checking ',release,' ',domain,': ',basename(p))
  x<-if(other||dup)NULL else tryCatch({
   z<-pt_read(p,domain)
   if(domain=='ADSL') {
    if(anyDuplicated(z$subject))stop('Duplicate ADSL subject')
    if(sum(ws_norm(z$RANDFL)=='Y')<PHMPT_MIN_RANDOMIZED)stop('Incomplete randomized ADSL')
    dd<-ws_dates(z$CUTOFFDT,'CUTOFFDT',PHMPT_NUMERIC_DATE_ORIGIN,TRUE)
    if(anyNA(dd)||any(ws_iso(dd)!=target))stop('Expected ',target,' ADSL cutoff')
   };z
  },error=function(e)e)
  good<-!other&&!dup&&!inherits(x,'error')
  logs[[i]]<-data.frame(release=release,domain=domain,path=normalizePath(p,winslash='/',mustWork=FALSE),md5=h,
   status=if(other)'OTHER_SUBMISSION' else if(dup)'IDENTICAL_BYTES' else if(good)'VALID' else 'REJECTED',
   rows=if(good)nrow(x) else NA_integer_,detail=if(!other&&!dup&&!good)conditionMessage(x) else '')
  if(good)valid[[length(valid)+1L]]<-list(x=x,index=i,signature=pt_signature(x))
  if(!other)seen<-c(seen,h)
 }
 log<-do.call(rbind,logs);logfile<-file.path(out,paste0('discovery_',domain,'.csv'));ws_csv(log,logfile)
 if(!length(valid))stop('No usable ',release,' ',domain,'; see ',logfile,': ',paste(unique(log$detail[nzchar(log$detail)]),collapse=' | '))
 if(length(unique(vapply(valid,function(a)a$signature,character(1))))>1)stop('Different ',domain,' snapshots remain. Select PHMPT_FILES$',release,'$',domain,'; tables are never concatenated.')
 stems<-tolower(tools::file_path_sans_ext(log$path))
 unresolved<-log$status=='REJECTED'&!stems%in%stems[log$status%in%c('VALID','IDENTICAL_BYTES')]
 if(any(unresolved)&&!nzchar(explicit))stop('Unresolved ',domain,' candidates; select one explicit reviewed path. See ',logfile)
 chosen<-valid[[1]];log$status[log$status=='VALID']<-'NOT_SELECTED';log$status[chosen$index]<-'SELECTED';ws_csv(log,logfile)
 list(data=chosen$x,manifest=log[chosen$index,,drop=FALSE])
}
pt_dates_write <- function(x,path) {
 for(n in names(x)[grepl('(^lo$|^hi$|^date$|_date$|^start$|^end$|_end$|^RANDDT$|^TR01SDT$|^DOSE2DT$|^CUTOFFDT$|^AP01EDT$|^UBLDDT$|^PDVDT$|^AP02SDT$|^TR02SDT$|^EUAVACDT$|^EOSDT$|^DTHDT$|^followup_start$)',names(x))])if(is.numeric(x[[n]]))x[[n]]<-ws_iso(x[[n]])
 ws_csv(x,path)
}
pt_context <- function(x,out) {
 if(any(!ws_norm(x$RANDFL)%in%c('Y','N')))stop('Unexpected RANDFL')
 x<-x[ws_norm(x$RANDFL)=='Y',,drop=FALSE]
 x$arm<-ifelse(ws_norm(x$TRT01P)=='MRNA-1273','mRNA-1273',ifelse(ws_norm(x$TRT01P)=='PLACEBO','Placebo','UNKNOWN'))
 if(any(x$arm=='UNKNOWN'))stop('Unexpected randomized treatment arm')
 for(n in c('RANDDT','CUTOFFDT','TR01SDT','DOSE2DT','AP01EDT','UBLDDT','PDVDT','AP02SDT','TR02SDT','EUAVACDT','EOSDT','DTHDT'))x[[n]]<-ws_dates(x[[n]],n,PHMPT_NUMERIC_DATE_ORIGIN,TRUE)
 if(anyNA(x$RANDDT)||any(x$RANDDT>x$CUTOFFDT))stop('Invalid randomization dates')
 x$site<-ws_norm(x$SITEID);x$site_inferred<-!nzchar(x$site);x$site[x$site_inferred]<-substr(x$subject[x$site_inferred],1,5)
 if(any(!grepl('^US[0-9]{3}$',x$site)))stop('Unexpected SITEID; review the source')
 if(any(x$site!=substr(x$subject,1,5)))stop('SITEID disagrees with participant ID')
 x$start<-x$RANDDT;x$followup_start<-ifelse(is.finite(x$TR01SDT),x$TR01SDT,x$RANDDT)+1
 x$calendar_end<-x$CUTOFFDT
 for(n in c('EOSDT','DTHDT'))x$calendar_end<-pmin(x$calendar_end,ifelse(is.finite(x[[n]]),x[[n]],Inf))
 x$blinded_end<-pmin(x$calendar_end,ifelse(is.finite(x$AP01EDT),x$AP01EDT,Inf))
 for(n in c('UBLDDT','PDVDT','AP02SDT','TR02SDT','EUAVACDT'))x$blinded_end<-pmin(x$blinded_end,ifelse(is.finite(x[[n]]),x[[n]]-1,Inf))
 # AP01EDT is an inclusive Part A end. Unblinding/decision/crossover days are
 # excluded when only dates are available; boundary-labelled records excluded too.
 x$missing_blind_boundary<-!is.finite(x$AP01EDT)&!is.finite(x$UBLDDT)&!is.finite(x$PDVDT)&!is.finite(x$AP02SDT)&!is.finite(x$TR02SDT)&!is.finite(x$EUAVACDT)
 x$strict_context_usable<-!x$missing_blind_boundary|(is.finite(x$EOSDT)&x$EOSDT<=ws_dates(PHMPT_EUA_CUTOFF,'EUA'))|(is.finite(x$DTHDT)&x$DTHDT<=ws_dates(PHMPT_EUA_CUTOFF,'EUA'))|x$CUTOFFDT<=ws_dates(PHMPT_EUA_CUTOFF,'EUA')
 pt_dates_write(x,file.path(out,'randomized_roster.csv'));x
}

pt_result <- function(x) {
 z<-ws_norm(x);r<-rep('UNKNOWN',length(z))
 r[z%in%c('NEGATIVE','NOT DETECTED','RNA NOT DETECTED','N')]<-'NEGATIVE'
 r[z%in%c('POSITIVE','DETECTED','RNA DETECTED','P','PRESUMPTIVE POSITIVE')]<-'POSITIVE'
 r[z%in%c('INDETERMINATE','INCONCLUSIVE','EQUIVOCAL','INVALID','SEECOMMENT','SEE COMMENT')]<-'NONBINARY'
 r
}
pt_unique <- function(x,keys,domain,out) {
 if(!nrow(x))return(x)
 cmp<-x;cmp$.source_row<-NULL;keep<-!duplicated(cmp)
 ws_csv(x[!keep,,drop=FALSE],file.path(out,paste0(domain,'_exact_duplicates_removed.csv')));x<-x[keep,,drop=FALSE]
 k<-do.call(paste,c(x[keys],sep='|'))
 if(anyDuplicated(k)) {
  ws_csv(x[duplicated(k)|duplicated(k,fromLast=TRUE),,drop=FALSE],file.path(out,paste0(domain,'_conflicting_keys.csv')))
  stop(domain,' has conflicting rows sharing a key. Resolve the snapshot; see audit output.')
 };x
}
pt_records <- function(m,a,s,ctx,out) {
 m<-m[ws_norm(m$MBTESTCD)=='SARSCOV2',,drop=FALSE]
 if(!nrow(m))stop('No SARSCOV2 MB records')
 if(any(!nzchar(m$MBSEQ)))stop('Missing MBSEQ; cannot safely reconcile PCR records')
 m<-pt_unique(m,c('subject','MBSEQ'),'MB',out);key<-paste(m$subject,m$MBSEQ,sep='|')
 if(!is.null(a)) {
  a<-a[ws_norm(a$PARAMCD)=='SARSCOV2',,drop=FALSE]
  discard<-nzchar(trimws(a$DTYPE))|!nzchar(a$MBSEQ)
  ws_csv(a[discard,,drop=FALSE],file.path(out,'ADMB_derived_or_unkeyed_excluded.csv'));a<-a[!discard,,drop=FALSE]
  a<-pt_unique(a,c('subject','MBSEQ'),'ADMB',out);ak<-paste(a$subject,a$MBSEQ,sep='|')
  ws_csv(a[!ak%in%key,,drop=FALSE],file.path(out,'ADMB_without_matching_MB.csv'));ai<-match(key,ak)
 } else ai<-rep(NA_integer_,nrow(m))
 av<-function(n)if(is.null(a))rep('',nrow(m)) else {z<-a[[n]][ai];z[is.na(z)]<-'';z}
 qvalue<-function(n) {
  if(is.null(s))return(rep('',nrow(m)))
  z<-s[ws_norm(s$IDVAR)=='MBSEQ'&ws_norm(s$QNAM)==n,,drop=FALSE]
  if(!nrow(z))return(rep('',nrow(m)))
  # Aggregate supplemental qualifiers without multiplying MB rows. Conflicting
  # values stay visible (e.g. Y | N), and are not resolved by row order.
  kk<-paste(z$subject,z$IDVARVAL,sep='|')
  vals<-vapply(split(ws_text(z$QVAL),kk),function(v)paste(sort(unique(trimws(v[nzchar(trimws(v))]))),collapse=' | '),character(1))
  r<-unname(vals[match(key,names(vals))]);r[is.na(r)]<-'';r
 }
 local_s<-ws_norm(qvalue('LOCALFL'));local_a<-ws_norm(av('LOCALFL'))
 allowed<-c('','Y','YES','N','NO')
 local_yes<-local_s%in%c('Y','YES')|local_a%in%c('Y','YES')
 local_no<-local_s%in%c('N','NO')|local_a%in%c('N','NO')
 flag_conflict<-!local_s%in%allowed|!local_a%in%allowed|(local_yes&local_no)
 grp<-m$MBGRPID;g_local<-grp%in%PHMPT_LOCAL_GROUPS;g_central<-grp%in%PHMPT_CENTRAL_GROUPS
 conflict<-flag_conflict|(g_central&local_yes)|(g_local&local_no)
 local<-g_local|local_yes
 central_name<-grepl(PHMPT_CENTRAL_NAME_REGEX,ws_norm(m$MBNAM),perl=TRUE)
 route<-ifelse(conflict,'CONFLICT',ifelse(local,'LOCAL_SOURCE',ifelse(g_central&central_name,'CENTRAL_SOURCE','UNCLASSIFIED')))
 basis<-ifelse(conflict,'Conflicting group/LOCALFL evidence',ifelse(local,ifelse(g_local&local_yes,'MBGRPID + LOCALFL',ifelse(g_local,'MBGRPID','LOCALFL')),
        ifelse(g_central&central_name,'Central group + expected central-lab name','Insufficient route evidence')))
 d<-ws_dates(m$MBDTC,'MB.MBDTC',PHMPT_NUMERIC_DATE_ORIGIN)
 da<-ws_dates(av('MBDTC'),'ADMB.MBDTC',PHMPT_NUMERIC_DATE_ORIGIN)
 date_conflict<-is.finite(d)&is.finite(da)&d!=da
 # Raw specimen date owns the cutoff. No partial raw dates are imputed from ADT.
 r1<-pt_result(m$MBSTRESC);r2<-pt_result(m$MBORRES);r3<-pt_result(av('AVALC'))
 pos<-r1=='POSITIVE'|r2=='POSITIVE'|r3=='POSITIVE';neg<-r1=='NEGATIVE'|r2=='NEGATIVE'|r3=='NEGATIVE'
 nonbinary<-r1=='NONBINARY'|r2=='NONBINARY'|r3=='NONBINARY'
 nd<-ws_norm(m$MBSTAT)%in%c('NOT DONE','NOT PERFORMED','CANCELLED','CANCELED')|ws_norm(av('MBSTAT'))%in%c('NOT DONE','NOT PERFORMED','CANCELLED','CANCELED')
 method<-ws_norm(m$MBMETHOD);method_ok<-grepl('PCR|POLYMERASE',method)&!grepl('ANTIGEN|ANTIBODY|SEROLOG',method)
 done<-method_ok&!nd&(pos|neg|nonbinary)
 result_conflict<-(pos&neg)|(nd&(pos|neg|nonbinary))
 visits<-paste(m$VISIT,av('VISIT'),av('AVISIT'),sep=' | ');vv<-ws_norm(visits)
 boundary<-grepl('PARTICIPANT DECISION|OPEN[- ]?LABEL|(^|[^A-Z])OL[- ]',vv)
 baseline<-ws_norm(av('ABLFL'))=='Y'|grepl('(^|[|] *)BASELINE($| *[|])|VISIT 1 DAY 1($| *[|])|(^|[|] *)DAY 1($| *[|])',vv)
 p<-data.frame(record_id=key,subject=m$subject,MBSEQ=m$MBSEQ,MBGRPID=grp,MBREFID=m$MBREFID,
  date=d,original_date=m$MBDTC,analysis_collection_date=av('MBDTC'),date_conflict=date_conflict,
  laboratory=m$MBNAM,method=m$MBMETHOD,raw_result=m$MBSTRESC,original_result=m$MBORRES,
  analysis_result=av('AVALC'),status=m$MBSTAT,reason_not_done=m$MBREASND,comment=qvalue('TSTCOM'),
  local_flag_supp=local_s,local_flag_admb=local_a,route=route,route_basis=basis,
  local_route_named_viracor=local&central_name,completed=done,binary=done&!result_conflict&(pos|neg),
  positive=pos,negative=neg,result_conflict=result_conflict,not_done=nd,method_ok=method_ok,
  visits=visits,baseline=baseline,illness_labelled=grepl("ILLNESS",vv),boundary_visit=boundary,raw_source_row=m$.source_row,
  analysis_source_row=if(is.null(a))NA_integer_ else a$.source_row[ai],stringsAsFactors=FALSE)
 p$usable_dated<-p$completed&is.finite(p$date)&!p$date_conflict
 p$in_randomized_roster<-p$subject%in%ctx$subject
 # MBREFID can flag potential repeated specimens, but cannot justify merging an
 # external result and a central result solely because subject/date agree.
 refkey<-paste(p$subject,p$MBREFID,sep='|');p$repeated_reference<-nzchar(p$MBREFID)&(duplicated(refkey)|duplicated(refkey,fromLast=TRUE))
 pt_dates_write(p,file.path(out,'pcr_record_audit.csv'))
 review<-!p$usable_dated|p$route%in%c('CONFLICT','UNCLASSIFIED')|p$result_conflict|p$repeated_reference|!p$in_randomized_roster
 pt_dates_write(p[review,,drop=FALSE],file.path(out,'pcr_records_for_review.csv'))
 inv<-aggregate(rep(1L,nrow(p)),p[c('route','MBGRPID','local_flag_supp','local_flag_admb','laboratory','method','status')],sum)
 names(inv)[ncol(inv)]<-'records';inv<-inv[order(-inv$records),];ws_csv(inv,file.path(out,'laboratory_classification_inventory.csv'))
 p
}

pt_window <- function(ctx,p,cutoff,scope,stratum,release,label) {
 x<-ctx
 x$start<-if(stratum!='ALL_POST_RANDOMIZATION')x$followup_start else x$RANDDT
 x$end<-pmin(if(scope=='STRICT_BLINDED')x$blinded_end else x$calendar_end,cutoff)
 x$eligible<-x$start<=x$end
 if(scope=='STRICT_BLINDED')x$eligible<-x$eligible&(x$strict_context_usable|cutoff<=ws_dates(PHMPT_EUA_CUTOFF,'EUA'))
 flow<-do.call(rbind,lapply(c('mRNA-1273','Placebo'),function(ar){v<-x[x$arm==ar,];data.frame(release=release,period=label,scope=scope,stratum=stratum,arm=ar,
  randomized=nrow(v),eligible=sum(v$eligible),no_observable_window=sum(v$start>v$end),missing_strict_context=sum(v$start<=v$end&!v$eligible))}))
 x<-x[x$eligible,,drop=FALSE];i<-match(p$subject,x$subject)
 take<-!is.na(i)&p$usable_dated&p$date>=x$start[i]&p$date<=x$end[i]
 if(scope=='STRICT_BLINDED')take<-take&!p$boundary_visit
 if(stratum!='ALL_POST_RANDOMIZATION')take<-take&!p$baseline
 if(stratum=='ILLNESS_LABELLED')take<-take&p$illness_labelled
 take[is.na(take)]<-FALSE;q<-p[take,,drop=FALSE]
 n<-nrow(x);x$days_observed<-x$end-x$start+1
 add<-function(route) {
  z<-if(route=='ANY')q else q[q$route==route,,drop=FALSE]
  ids<-match(z$subject,x$subject)
  list(records=tabulate(ids,nbins=n),days=tabulate(match(unique(z[c('subject','date')])$subject,x$subject),nbins=n))
 }
 for(r in c('ANY','LOCAL_SOURCE','CENTRAL_SOURCE','UNCLASSIFIED','CONFLICT')) {
  v<-add(r);nm<-tolower(r);x[[paste0(nm,'_records')]]<-v$records;x[[paste0(nm,'_days')]]<-v$days
  x[[paste0(nm,'_tested')]]<-v$records>0
 }
 x$both_tested<-x$local_source_tested&x$central_source_tested
 x$unresolved_tested<-x$unclassified_tested|x$conflict_tested
 x$category<-ifelse(x$both_tested,'Both source streams',ifelse(x$local_source_tested,'Local-source only',
    ifelse(x$central_source_tested,'Central-source only',ifelse(x$any_tested,'Unclassified/conflicting only','No completed PCR'))))
 # "only" concerns known route evidence: the accompanying unresolved flag must
 # remain visible; local-only does not mean central testing can be excluded.
 keep<-c('subject','arm','site','start','end','days_observed','category','unresolved_tested',
         names(x)[grepl('_(records|days|tested)$',names(x))])
 x<-x[,unique(keep),drop=FALSE];x$release<-release;x$period<-label;x$scope<-scope;x$stratum<-stratum
 list(subjects=x,flow=flow)
}
pt_ci <- function(k,n) {
 if(n==0)return(c(NA_real_,NA_real_));z<-qnorm(.975);p<-k/n;a<-(p+z*z/(2*n))/(1+z*z/n)
 b<-z*sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n);100*c(a-b,a+b)
}
pt_summary <- function(x,meta=NULL,sites=NULL) {
 if(is.null(meta))meta<-x[1,c('release','period','scope','stratum'),drop=FALSE]
 if(is.null(sites))sites<-sort(unique(x$site))
 endpoints<-c('any_tested','local_source_tested','central_source_tested','both_tested','unresolved_tested')
 recnames<-c('any_records','local_source_records','central_source_records',NA,NA)
 daynames<-c('any_days','local_source_days','central_source_days',NA,NA)
 labels<-c('Any PCR','Local-source PCR','Central-source PCR','Both source streams','Unclassified/conflicting PCR')
 bysite<-setNames(lapply(sites,function(st)x[x$site==st,,drop=FALSE]),sites);groups<-c(list(ALL=x),bysite);rows<-list();j<-0L
 for(site in names(groups))for(ar in c('mRNA-1273','Placebo')) {
  z<-groups[[site]];z<-z[z$arm==ar,,drop=FALSE];n<-nrow(z);pd<-sum(z$days_observed)
  for(e in seq_along(endpoints)) {
   k<-sum(z[[endpoints[e]]]);ci<-pt_ci(k,n);nr<-if(is.na(recnames[e]))NA_integer_ else sum(z[[recnames[e]]]);nd<-if(is.na(daynames[e]))NA_integer_ else sum(z[[daynames[e]]])
   j<-j+1L;rows[[j]]<-data.frame(release=meta$release[1],period=meta$period[1],scope=meta$scope[1],stratum=meta$stratum[1],
    granularity=if(site=='ALL')'GLOBAL' else 'SITE',site=site,arm=ar,endpoint=labels[e],n_eligible=n,n_tested=k,
    pct=if(n)100*k/n else NA_real_,ci95_low=ci[1],ci95_high=ci[2],PCR_records=nr,unique_PCR_days=nd,
    person_days=pd,records_per_100_subjects=if(n)100*nr/n else NA_real_,records_per_100_person_days=if(pd)100*nr/pd else NA_real_,stringsAsFactors=FALSE)
  }
 };do.call(rbind,rows)
}
pt_compare <- function(s) {
 keys<-c('release','period','scope','stratum','granularity','site','endpoint');v<-s[s$arm=='mRNA-1273',];p<-s[s$arm=='Placebo',]
 key<-function(z)do.call(paste,c(z[keys],sep='|'));p<-p[match(key(v),key(p)),]
 rows<-lapply(seq_len(nrow(v)),function(i){a<-v$n_tested[i];b<-v$n_eligible[i]-a;c<-p$n_tested[i];d<-p$n_eligible[i]-c
  mat<-matrix(c(a,b,c,d),nrow=2,byrow=TRUE);method<-'Not estimable';pv<-NA_real_;chi<-NA_real_;emin<-NA_real_
  if(all(rowSums(mat)>0)) {
   if(any(colSums(mat)==0)){method<-'No outcome variation';pv<-1}
   else {ex<-outer(rowSums(mat),colSums(mat))/sum(mat);emin<-min(ex)
    if(all(ex>=5)){r<-chisq.test(mat,correct=FALSE);method<-'Pearson chi-square';pv<-r$p.value;chi<-unname(r$statistic)}
    else {method<-'Fisher exact (two-sided)';pv<-fisher.test(mat)$p.value}
   }
  }
  data.frame(v[i,keys,drop=FALSE],v_tested=a,v_total=a+b,p_tested=c,p_total=c+d,
   v_pct=v$pct[i],p_pct=p$pct[i],difference_pp=v$pct[i]-p$pct[i],
   risk_ratio=if(is.finite(p$pct[i])&&p$pct[i]>0)v$pct[i]/p$pct[i] else NA_real_,
   method=method,min_expected=emin,chi_square=chi,p_value=pv,stringsAsFactors=FALSE)
 });r<-do.call(rbind,rows)
 # Families: all endpoints/sites within release x cutoff x scope x stratum;
 # global comparisons are a separate family. No adjustment across these views.
 fam<-interaction(r$release,r$period,r$scope,r$stratum,r$granularity,drop=TRUE)
 r$p_BH<-NA_real_;r$p_Holm<-NA_real_
 for(g in levels(fam)){i<-which(fam==g);r$p_BH[i]<-p.adjust(r$p_value[i],'BH');r$p_Holm[i]<-p.adjust(r$p_value[i],'holm')}
 r$nominal_p_lt_05<-is.finite(r$p_value)&r$p_value<.05
 r$BH_lt_05<-is.finite(r$p_BH)&r$p_BH<.05;r
}
pt_p <- function(p) {
 v<-rep('NA',length(p));ok<-is.finite(p);v[ok]<-sprintf('%.3f',p[ok])
 # Successively tighter thresholds, through five zeros after the decimal point.
 for(t in c(.01,.001,.0001,.00001,.000001))v[ok&p<t]<-paste0('< ',format(t,scientific=FALSE,trim=TRUE))
 v
}

pt_escape <- function(x) {x<-as.character(x);x[is.na(x)]<-'';x<-gsub('&','&amp;',x,fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE);x<-gsub('>','&gt;',x,fixed=TRUE);gsub('"','&quot;',x,fixed=TRUE)}
pt_table <- function(x) {
 if(!nrow(x))return('<p>No rows.</p>')
 paste0('<div class="scroll"><table><thead><tr>',paste0('<th>',pt_escape(names(x)),'</th>',collapse=''),'</tr></thead><tbody>',
  paste(apply(x,1,function(z)paste0('<tr>',paste0('<td>',pt_escape(z),'</td>',collapse=''),'</tr>')),collapse=''),'</tbody></table></div>')
}
pt_display <- function(t,site=FALSE) {
 fmt<-function(k,n,p)ifelse(n>0,sprintf('%s/%s (%.2f%%)',k,n,p),'0/0 (not estimable)')
 x<-data.frame(Endpoint=t$endpoint,`mRNA-1273 tested`=fmt(t$v_tested,t$v_total,t$v_pct),
  `Placebo tested`=fmt(t$p_tested,t$p_total,t$p_pct),`Difference (pp)`=sprintf('%+.2f',t$difference_pp),
  Test=t$method,p=pt_p(t$p_value),`BH-adjusted p`=pt_p(t$p_BH),check.names=FALSE)
 if(site)x<-cbind(Site=t$site,x);x
}
pt_report <- function(s,t,flow,release_status,out) {
 parts<-c('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Moderna P301 — local and central PCR testing</title>',
 '<style>body{font:16px/1.55 Georgia,serif;color:#17212a;background:#fafaf8;margin:0}main{max-width:1200px;margin:25px auto;padding:30px;background:white}h1,h2,h3{font-family:Arial,sans-serif;line-height:1.25}h1{font-size:29px}h2{font-size:23px;margin-top:32px}h3{font-size:18px}p{max-width:100ch}.note{padding:15px 18px;background:#edf3f7;border-left:4px solid #42637e}.warning{background:#fff5df;border-left-color:#bc8829}table{width:100%;border-collapse:collapse;font:12px/1.4 Arial,sans-serif}th{text-align:left;background:#e5eaee}th,td{padding:8px;border-bottom:1px solid #dbe0e4;vertical-align:top}tr:nth-child(even){background:#f5f7f8}.scroll{overflow:auto}details{margin:20px 0}summary{font:700 17px Arial,sans-serif;cursor:pointer}a{color:#185182}button{padding:8px 12px;margin:8px 8px 8px 0}.muted{color:#52606a;font:13px/1.5 Arial,sans-serif}code{overflow-wrap:anywhere}section.view{border-top:1px solid #c9d0d6;margin-top:28px}@media print{body{background:white}main{max-width:none;margin:0;padding:0}button{display:none}.scroll{overflow:visible}table{font-size:8pt}th,td{padding:5pt}thead{display:table-header-group}tr{break-inside:avoid}h2,h3{break-after:avoid}section.view{break-before:page}p{font-size:10pt}details{display:block}summary{display:none}details>*{display:block!important}a{color:inherit}h1{font-size:20pt}}@page{size:A4 landscape;margin:14mm}</style></head><body><main>',
 '<h1>Local and central PCR testing by randomized arm</h1>',
 paste0('<p class="muted">Moderna P301 · script v1.1.1 · Generated ',pt_escape(format(Sys.time(),tz='UTC',usetz=TRUE)),'.</p>'),
 '<p class="note"><strong>Denominator:</strong> all randomized participants with an observable window, by original treatment assignment. The headline rate is the percentage with at least one completed PCR record in that window. It is not the probability of testing after symptoms: that requires linking the COVID forms and AE episodes. No baseline infection-negative restriction or per-protocol exclusion is applied.</p>',
 '<p><strong>Source-route definitions:</strong> local-source means MBGRPID=003 or LOCALFL=Y. This denotes the COVID Diagnostic Test reporting route. Central-source requires the configured central group (default 001), an expected central-laboratory name (Viracor), and no conflicting local flag. Neither a blank local flag nor a brand name alone establishes the route. Viracor-labelled records in group 003 remain local-source. Unclassified and conflicting records are retained as separate audit categories; these labels do not establish the physical place where an assay was performed.</p>',
 '<p><strong>Completed PCR:</strong> an explicit PCR method plus a documented positive, negative, or nonbinary result, including SeeComment/inconclusive records, without a NOT DONE/cancelled status. The raw MB collection date must be complete and agree with the corresponding ADMB collection date when both exist. Positive/negative discrepancies still document testing but are flagged; a contradiction involving NOT DONE is excluded. Raw MB subject/MBSEQ is the record key. ADMB and SUPPMB enrich it and are never added as additional tests.</p>',
 paste0('<p><strong>Windows:</strong> EUA = ',PHMPT_EUA_CUTOFF,'; BLA = ',PHMPT_BLA_CUTOFF,'. The actual EUA release and a retrospective EUA-date restriction of the production release are labelled separately. STRICT_BLINDED also censors at Part A end, unblinding, participant decision, crossover/outside vaccination, death or study end; it excludes decision/open-label-labelled records. Unblinding/decision/crossover days are excluded when times are unavailable. CALENDAR follows the original randomized groups through the date cutoff and can include post-unblinding/crossover observations; it is descriptive.</p>'),
 '<p><strong>Baseline separation:</strong> ALL_POST_RANDOMIZATION starts on randomization day and includes baseline/day-1 tests. FOLLOW_UP_EXCL_BASELINE starts the day after dose 1 (the day after randomization if undosed) and excludes baseline-labelled records. This still includes scheduled Day 29 testing. ILLNESS_LABELLED restricts the same follow-up window to records explicitly labelled as illness visits in MB or ADMB; records without that label are not assumed to be illness-related. Both denominators require at least one eligible follow-up day, regardless of whether symptoms were reported. For a post-EUA strict window, participants lacking all Part A/unblinding/decision/crossover boundary dates are excluded unless death or study end preceded the EUA date. They are separately inventoried, not silently counted as untested. The calendar view retains them when observable. Person-days use these participant-specific windows; record rates are descriptive.</p>',
 '<p><strong>Inference:</strong> two-sided Pearson chi-square when every expected 2×2 cell is at least five; otherwise two-sided Fisher exact. No-variation comparisons have p=1; a missing arm is not estimable. Tests use participants, not repeated assays. BH and Holm adjustments are exported within each release/cutoff/scope/stratum, with global comparisons separated from the family of all site×endpoint comparisons. Different views are not independent replications. Differences in illness burden, follow-up, visit scheduling or recording route can explain testing differences; these comparisons alone do not establish unequal investigator handling.</p>',
 '<button onclick="window.print()">Print report</button><button onclick="document.querySelectorAll(\'details\').forEach(x=>x.open=true)">Expand site tables</button>',
 '<h2>Inputs and coverage</h2>',pt_table(release_status),
 '<p><a href="summary_global.csv">Global summary CSV</a> · <a href="summary_by_site.csv">Site summary CSV</a> · <a href="comparisons_by_site.csv">Site tests CSV</a> · <a href="population_flow.csv">Population flow CSV</a> · <a href="source_manifest.csv">Source manifest CSV</a> · <a href="parameters.csv">Analysis settings CSV</a></p>',
 '<details><summary>Population flow and excluded windows</summary>',pt_table(flow),'</details>')
 if(any(release_status$status=='SKIPPED'))parts<-c(parts,'<p class="note warning">A requested release was unavailable. Its absence is reported above. A retrospective date restriction of production records does not reproduce an unavailable EUA database snapshot.</p>')
 global<-t[t$granularity=='GLOBAL',];views<-unique(global[c('release','period','scope','stratum')])
 for(i in seq_len(nrow(views))) {
  v<-views[i,];take<-t$release==v$release&t$period==v$period&t$scope==v$scope&t$stratum==v$stratum
  g<-t[take&t$granularity=='GLOBAL',];sites<-t[take&t$granularity=='SITE',]
  sites<-sites[order(sites$endpoint,sites$p_value,sites$site,na.last=TRUE),]
  title<-paste(v,collapse=' · ')
  parts<-c(parts,paste0('<section class="view"><h2>',pt_escape(title),'</h2>'),pt_table(pt_display(g)),
   '<p class="muted">Local and central percentages overlap: a participant may have both. Both-source is a participant history, not proof of paired testing of the same specimen. “Unclassified/conflicting” may overlap either known route. Exact p values and Holm-adjusted values are in the CSV.</p>',
   '<details><summary>Site-level comparisons — all sites and all endpoints</summary>',pt_table(pt_display(sites,TRUE)),'</details></section>')
 }
 parts<-c(parts,'<h2>Audit assets</h2><p>Each input release has its own folder with the randomized roster, one-row-per-MBSEQ PCR audit, quality-review rows, laboratory inventory, exact duplicates, discovery decisions and optional ADMB-only records. The subjects/ folder contains one participant-level file per analysis window. No files from different submissions are silently stacked. PCR record counts may include separate reports of the same physical assay; repeated specimen references are flagged rather than merged without evidence.</p>',
 '<p class="muted">Definition references: P301 Clinical Study Data Reviewer’s Guide, section 3.4.18 (LOCALFL: From COVID Diag Test); P301 Analysis Data Reviewer’s Guide, section 5.2.10 (ADMB derives from MB); reviewed PHMPT PCR metadata export v1.9.0. The operational central-stream rule is editable and should be checked against each release’s Define metadata. <a href="https://www.fda.gov/media/144673/download">FDA EUA review memorandum</a>, section 5.5, identifies Viracor as the contracted RT-qPCR laboratory; this does not imply that every Viracor-labelled result is a central-route record.</p>',
 '<script>let previouslyOpen=[];addEventListener("beforeprint",()=>{previouslyOpen=[...document.querySelectorAll("details")].map(x=>x.open);document.querySelectorAll("details").forEach(x=>x.open=true)});addEventListener("afterprint",()=>document.querySelectorAll("details").forEach((x,i)=>x.open=previouslyOpen[i]));</script></main></body></html>')
 writeLines(enc2utf8(parts),file.path(out,'pcr_testing_report.html'),useBytes=TRUE)
}

ws_mapping <- function(root,path,out) {
  if(!nzchar(path)) {
    p<-file.path(root,c("assets/ae_covid_symptom_mapping.csv","ae_covid_symptom_mapping.csv"))
    path<-if(any(file.exists(p)))p[which(file.exists(p))[1]] else ""
  } else path<-ws_path(root,path)
  if(nzchar(path)) {
    if(!file.exists(path))stop("Mapping asset not found: ",path)
    x<-read.csv(path,colClasses="character",check.names=FALSE,fileEncoding="UTF-8-BOM",na.strings=c("","NA"))
    source<-normalizePath(path,winslash="/");hash<-unname(tools::md5sum(path))
  } else {x<-ws_embedded_mapping();source<-"Embedded reviewed mapping from weekly_symptoms_ae_20261021_055854";hash<-"";message(source)}
  need<-c("match_field","normalized_term","symptom_code","symptom_group","match_class","include_default")
  if(!all(need %in% names(x)))stop("Mapping lacks required fields: ",paste(setdiff(need,names(x)),collapse=", "))
  x[]<-lapply(x,ws_text);x$normalized_term<-ws_norm(x$normalized_term)
  x$match_field<-ws_norm(x$match_field);x$include_default<-ws_norm(x$include_default)
  x$symptom_code<-ws_norm(x$symptom_code);x$symptom_group<-ws_norm(x$symptom_group)
  x$match_class<-ws_norm(x$match_class)
  if(any(!nzchar(x$normalized_term)))stop("Mapping has an empty normalized_term.")
  if(any(!x$include_default %in% c("TRUE","FALSE")))stop("Mapping include_default must be TRUE/FALSE.")
  x$include_default<-x$include_default=="TRUE"
  if(any(!x$match_field %in% c("AETERM","AELLT","AEDECOD")))stop("Unknown mapping match_field.")
  if(any(!x$symptom_code %in% names(ws_param_map)))stop("Unknown mapping symptom_code.")
  if(any(x$symptom_group!=unname(ws_param_map[x$symptom_code])))stop("Mapping symptom groups differ from the agreed 12-group classification.")
  if(any(x$include_default & x$match_class!="DIRECT"))stop("Review-only mapping rows cannot be included without explicit reclassification.")
  # Several verbatim spellings can define the SAME normalized matching rule.
  # Keep the complete source rows; deduplicate by decisions, not by spelling,
  # labels, notes or dictionary provenance. Never silently resolve disagreement.
  key<-paste(x$match_field,x$normalized_term,x$symptom_code,sep="\r")
  group<-match(key,unique(key));groups<-split(seq_len(nrow(x)),group)
  x$mapping_source_row<-seq_len(nrow(x))
  x$mapping_rule_id<-sprintf("MAP%04d",group)
  ws_csv(x,file.path(out,"ae_covid_symptom_mapping_input.csv"))
  # Duration/qualification metadata are retained and checked for agreement even
  # though this burden analysis does not impose protocol-trigger eligibility.
  decisions<-intersect(c("symptom_group","match_class","include_default",
                        "min_duration_hours","qualification_check"),names(x))
  conflicts<-vapply(groups,function(ix) {
    z<-x[ix,decisions,drop=FALSE]
    z[]<-lapply(z,function(v)if(is.logical(v))v else ws_norm(v))
    nrow(unique(z))>1L
  },logical(1))
  bad<-unlist(groups[conflicts],use.names=FALSE)
  conflict_file<-file.path(out,"ae_covid_symptom_mapping_conflicts.csv")
  ws_csv(x[bad,,drop=FALSE],conflict_file)
  if(length(bad))stop("Conflicting mapping decisions for ",sum(conflicts),
    " normalized rule(s). Review ",normalizePath(conflict_file,winslash="/",mustWork=FALSE),
    ". No conflicting rule was selected automatically.")
  keep<-!duplicated(key);used<-x[keep,,drop=FALSE]
  used$mapping_source_rows<-vapply(group[keep],function(g)paste(which(group==g),collapse=";"),character(1))
  if("ae_term" %in% names(x))used$mapping_term_variants<-vapply(group[keep],function(g)
    paste(unique(x$ae_term[group==g]),collapse=" | "),character(1))
  collapsed<-x[duplicated(key)|duplicated(key,fromLast=TRUE),,drop=FALSE]
  ws_csv(collapsed,file.path(out,"ae_covid_symptom_mapping_duplicates.csv"))
  ws_csv(used,file.path(out,"ae_covid_symptom_mapping_used.csv"))
  ws_csv(data.frame(path=source,md5=hash,input_rows=nrow(x),effective_rules=nrow(used),
    equivalent_rows_collapsed=nrow(x)-nrow(used)),file.path(out,"mapping_source.csv"))
  if(nrow(x)>nrow(used))message("Mapping: ",nrow(x)," source rows -> ",nrow(used),
    " effective rules; ",nrow(x)-nrow(used)," equivalent duplicate row(s) collapsed. Original rows retained in mapping audit CSVs.")
  used
}
ws_param_map <- c(FEVER="FEVER",CHILLS="CHILLS",COUGH="COUGH",DIBREATH="DYSPNEA",SHBREATH="DYSPNEA",
  FATIGUE="FATIGUE",BODYACHE="MYALGIA",MYALGIA="MYALGIA",HEADACHE="HEADACHE",NLSMELL="SMELL_TASTE",
  NLTASTE="SMELL_TASTE",STHROAT="SORE_THROAT",NASALCON="NASAL",RUNNOSE="NASAL",NAUSEA="NAUSEA_VOMITING",
  VOMIT="NAUSEA_VOMITING",DIARRHEA="DIARRHEA")
ws_events <- function(ae,sy,mapping,origin,out,ae_timing) {
  evidence<-list();unmatched<-list()
  for(field in c("AETERM","AELLT","AEDECOD")) {
    terms<-ws_norm(ae[[field]]);rules<-mapping[mapping$match_field==field,,drop=FALSE]
    ix<-data.frame(row_id=seq_len(nrow(ae)),normalized_term=terms,stringsAsFactors=FALSE)
    hit<-merge(ix,rules,by="normalized_term",sort=FALSE)
    if(nrow(hit)) {hit$priority<-match(field,c("AETERM","AELLT","AEDECOD"));evidence[[length(evidence)+1L]]<-hit}
    bad<-nzchar(terms)&!terms %in% rules$normalized_term
    if(any(bad)) {
      u<-as.data.frame(table(terms[bad]),stringsAsFactors=FALSE);names(u)<-c("normalized_term","records")
      u$match_field<-field;unmatched[[length(unmatched)+1L]]<-u
    }
  }
  if(length(unmatched))ws_csv(do.call(rbind,unmatched),file.path(out,"unmapped_ae_terms.csv"))
  ev<-data.frame(subject=character(),source=character(),source_id=character(),source_domain=character(),source_sequence=character(),symptom_code=character(),symptom_group=character(),
    lo=numeric(),hi=numeric(),date_quality=character(),term=character(),stringsAsFactors=FALSE)
  if(length(evidence)) {
    h<-do.call(rbind,evidence);h<-h[order(h$row_id,h$symptom_code,h$priority),,drop=FALSE]
    h$selected_evidence<-!duplicated(paste(h$row_id,h$symptom_code,sep="|"))
    h$USUBJID<-ae$USUBJID[h$row_id];h$AESEQ<-ae$AESEQ[h$row_id]
    ws_csv(h,file.path(out,"ae_matching_evidence.csv"))
    # Apply evidence priority BEFORE filtering include_default: a specific
    # review-only term must not be overridden by a generic direct preferred term.
    h<-h[h$selected_evidence & h$include_default,,drop=FALSE];j<-h$row_id
    lo<-ws_dates(ae$AESTDTC,"AESTDTC",origin);ends<-ws_dates(ae$AEENDTC,"AEENDTC",origin)
    quality<-ifelse(!is.finite(lo),"UNKNOWN_START",ifelse(is.finite(ends)&ends<lo,"END_BEFORE_START",
      ifelse(is.finite(ends),"KNOWN_INTERVAL","START_ONLY_END_UNKNOWN")))
    hi<-if(ae_timing=="onset")lo else ifelse(is.finite(ends),ends,lo)
    hi[quality=="END_BEFORE_START"]<-NA_real_
    ev<-data.frame(subject=ae$subject[j],source=rep("AE",length(j)),source_id=if(length(j))paste0("AE:",ae$AESEQ[j]) else character(),
      source_domain=rep("AE",length(j)),source_sequence=ae$AESEQ[j],symptom_code=h$symptom_code,
      symptom_group=h$symptom_group,lo=lo[j],hi=hi[j],date_quality=quality[j],term=ae$AETERM[j],stringsAsFactors=FALSE)
  }
  code<-ws_norm(sy$PARAMCD);value<-ws_norm(sy$AVALC);core<-code %in% names(ws_param_map)
  inventory<-as.data.frame(table(parameter=code,result=value,source_domain=sy$SRCDOM),stringsAsFactors=FALSE);inventory<-inventory[inventory$Freq>0,]
  inventory$in_mapping<-inventory$parameter %in% names(ws_param_map)
  ws_csv(inventory,file.path(out,"symptom_parameter_inventory.csv"))
  positive<-core & value %in% c("Y","YES");j<-which(positive)
  date<-ws_dates(sy$ADT,"ADSYMP.ADT",origin)
  # The same symptom/date can have disagreeing assessments. Retain positive
  # evidence, flag the conflict, and never interpret missing/unknown as negative.
  keys<-paste(sy$subject,code,ws_iso(date),sep="|")
  conflict<-keys %in% keys[core & value %in% c("N","NO")]
  if(any(positive&conflict))ws_csv(sy[positive&conflict,,drop=FALSE],file.path(out,"contradictory_symptom_assessments.csv"))
  se<-data.frame(subject=sy$subject[j],source=rep("COVID_FORM",length(j)),source_id=if(length(j))paste0("ADSYMP:",ifelse(nzchar(sy$ASEQ[j]),sy$ASEQ[j],paste0("ROW_",j))) else character(),
    source_domain=sy$SRCDOM[j],source_sequence=sy$SRCSEQ[j],
    symptom_code=code[j],symptom_group=unname(ws_param_map[code[j]]),lo=date[j],hi=date[j],
    date_quality=ifelse(is.finite(date[j]),ifelse(conflict[j],"POSITIVE_WITH_NEGATIVE_SAME_DAY","POSITIVE_ASSESSMENT_DAY"),"UNKNOWN_ASSESSMENT_DATE"),
    term=sy$PARAM[j],stringsAsFactors=FALSE)
  ev<-unique(rbind(ev,se));ev$event_id<-seq_len(nrow(ev));ev
}
ws_embedded_mapping <- function() read.csv(text=paste(c(
  "\"dictionary_version\",\"match_field\",\"ae_term\",\"normalized_term\",\"symptom_code\",\"symptom_label\",\"symptom_group\",\"match_class\",\"include_default\",\"min_duration_hours\",\"qualification_check\",\"rationale\",\"mapping_source_row\",\"mapping_rule_id\",\"mapping_source_rows\",\"mapping_term_variants\"",
  "\"1.0.0\",\"AEDECOD\",\"Post procedural fever\",\"POST PROCEDURAL FEVER\",\"FEVER\",\"Fever\",\"FEVER\",\"DIRECT\",TRUE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",1,\"MAP0001\",\"1\",\"Post procedural fever\"",
  "\"1.0.0\",\"AEDECOD\",\"Pyrexia\",\"PYREXIA\",\"FEVER\",\"Fever\",\"FEVER\",\"DIRECT\",TRUE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",2,\"MAP0002\",\"2\",\"Pyrexia\"",
  "\"1.0.0\",\"AELLT\",\"Febrile reaction\",\"FEBRILE REACTION\",\"FEVER\",\"Fever\",\"FEVER\",\"DIRECT\",TRUE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",3,\"MAP0003\",\"3\",\"Febrile reaction\"",
  "\"1.0.0\",\"AELLT\",\"Fever\",\"FEVER\",\"FEVER\",\"Fever\",\"FEVER\",\"DIRECT\",TRUE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",4,\"MAP0004\",\"4\",\"Fever\"",
  "\"1.0.0\",\"AELLT\",\"Fever chills\",\"FEVER CHILLS\",\"FEVER\",\"Fever\",\"FEVER\",\"DIRECT\",TRUE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",5,\"MAP0005\",\"5\",\"Fever chills\"",
  "\"1.0.0\",\"AELLT\",\"Fever of unknown origin\",\"FEVER OF UNKNOWN ORIGIN\",\"FEVER\",\"Fever\",\"FEVER\",\"DIRECT\",TRUE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",6,\"MAP0006\",\"6\",\"Fever of unknown origin\"",
  "\"1.0.0\",\"AELLT\",\"Intermittent fever\",\"INTERMITTENT FEVER\",\"FEVER\",\"Fever\",\"FEVER\",\"DIRECT\",TRUE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",7,\"MAP0007\",\"7\",\"Intermittent fever\"",
  "\"1.0.0\",\"AELLT\",\"Post procedural fever\",\"POST PROCEDURAL FEVER\",\"FEVER\",\"Fever\",\"FEVER\",\"DIRECT\",TRUE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",8,\"MAP0008\",\"8\",\"Post procedural fever\"",
  "\"1.0.0\",\"AELLT\",\"Postoperative fever\",\"POSTOPERATIVE FEVER\",\"FEVER\",\"Fever\",\"FEVER\",\"DIRECT\",TRUE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",9,\"MAP0009\",\"9\",\"Postoperative fever\"",
  "\"1.0.0\",\"AELLT\",\"Pyrexia\",\"PYREXIA\",\"FEVER\",\"Fever\",\"FEVER\",\"DIRECT\",TRUE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",10,\"MAP0010\",\"10\",\"Pyrexia\"",
  "\"1.0.0\",\"AETERM\",\"FEVER\",\"FEVER\",\"FEVER\",\"Fever\",\"FEVER\",\"DIRECT\",TRUE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",11,\"MAP0011\",\"11\",\"FEVER\"",
  "\"1.0.0\",\"AETERM\",\"FEVER OF UNKNOWN ORIGIN\",\"FEVER OF UNKNOWN ORIGIN\",\"FEVER\",\"Fever\",\"FEVER\",\"DIRECT\",TRUE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",12,\"MAP0012\",\"12\",\"FEVER OF UNKNOWN ORIGIN\"",
  "\"1.0.0\",\"AETERM\",\"INTERMITTENT FEVER\",\"INTERMITTENT FEVER\",\"FEVER\",\"Fever\",\"FEVER\",\"DIRECT\",TRUE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",13,\"MAP0013\",\"13\",\"INTERMITTENT FEVER\"",
  "\"1.0.0\",\"AETERM\",\"PYREXIA\",\"PYREXIA\",\"FEVER\",\"Fever\",\"FEVER\",\"DIRECT\",TRUE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",14,\"MAP0014\",\"14\",\"PYREXIA\"",
  "\"1.0.0\",\"AEDECOD\",\"Body temperature increased\",\"BODY TEMPERATURE INCREASED\",\"FEVER\",\"Fever\",\"FEVER\",\"RELATED_REVIEW\",FALSE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Raised temperature or feverishness may be below 38 C; verify measured temperature.\",15,\"MAP0015\",\"15\",\"Body temperature increased\"",
  "\"1.0.0\",\"AELLT\",\"Body temperature increased\",\"BODY TEMPERATURE INCREASED\",\"FEVER\",\"Fever\",\"FEVER\",\"RELATED_REVIEW\",FALSE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Raised temperature or feverishness may be below 38 C; verify measured temperature.\",16,\"MAP0016\",\"16\",\"Body temperature increased\"",
  "\"1.0.0\",\"AELLT\",\"Feverish\",\"FEVERISH\",\"FEVER\",\"Fever\",\"FEVER\",\"RELATED_REVIEW\",FALSE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Raised temperature or feverishness may be below 38 C; verify measured temperature.\",17,\"MAP0017\",\"17\",\"Feverish\"",
  "\"1.0.0\",\"AELLT\",\"High temperature\",\"HIGH TEMPERATURE\",\"FEVER\",\"Fever\",\"FEVER\",\"RELATED_REVIEW\",FALSE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Raised temperature or feverishness may be below 38 C; verify measured temperature.\",18,\"MAP0018\",\"18\",\"High temperature\"",
  "\"1.0.0\",\"AELLT\",\"Low grade fever\",\"LOW GRADE FEVER\",\"FEVER\",\"Fever\",\"FEVER\",\"RELATED_REVIEW\",FALSE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Raised temperature or feverishness may be below 38 C; verify measured temperature.\",19,\"MAP0019\",\"19\",\"Low grade fever\"",
  "\"1.0.0\",\"AELLT\",\"Slight temperature\",\"SLIGHT TEMPERATURE\",\"FEVER\",\"Fever\",\"FEVER\",\"RELATED_REVIEW\",FALSE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Raised temperature or feverishness may be below 38 C; verify measured temperature.\",20,\"MAP0020\",\"20\",\"Slight temperature\"",
  "\"1.0.0\",\"AETERM\",\"FEVERISH\",\"FEVERISH\",\"FEVER\",\"Fever\",\"FEVER\",\"RELATED_REVIEW\",FALSE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Raised temperature or feverishness may be below 38 C; verify measured temperature.\",21,\"MAP0021\",\"21\",\"FEVERISH\"",
  "\"1.0.0\",\"AETERM\",\"HIGH TEMPERATURE\",\"HIGH TEMPERATURE\",\"FEVER\",\"Fever\",\"FEVER\",\"RELATED_REVIEW\",FALSE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Raised temperature or feverishness may be below 38 C; verify measured temperature.\",22,\"MAP0022\",\"22\",\"HIGH TEMPERATURE\"",
  "\"1.0.0\",\"AETERM\",\"LOW GRADE FEVER\",\"LOW GRADE FEVER\",\"FEVER\",\"Fever\",\"FEVER\",\"RELATED_REVIEW\",FALSE,\"0\",\"Verify measured temperature >=38 C; AE label alone is insufficient.\",\"Raised temperature or feverishness may be below 38 C; verify measured temperature.\",23,\"MAP0023\",\"23\",\"LOW GRADE FEVER\"",
  "\"1.0.0\",\"AEDECOD\",\"Chills\",\"CHILLS\",\"CHILLS\",\"Chills\",\"CHILLS\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",24,\"MAP0024\",\"24\",\"Chills\"",
  "\"1.0.0\",\"AELLT\",\"Chills\",\"CHILLS\",\"CHILLS\",\"Chills\",\"CHILLS\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",25,\"MAP0025\",\"25\",\"Chills\"",
  "\"1.0.0\",\"AELLT\",\"Fever chills\",\"FEVER CHILLS\",\"CHILLS\",\"Chills\",\"CHILLS\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",26,\"MAP0026\",\"26\",\"Fever chills\"",
  "\"1.0.0\",\"AELLT\",\"Rigors\",\"RIGORS\",\"CHILLS\",\"Chills\",\"CHILLS\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",27,\"MAP0027\",\"27\",\"Rigors\"",
  "\"1.0.0\",\"AETERM\",\"CHILL\",\"CHILL\",\"CHILLS\",\"Chills\",\"CHILLS\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",28,\"MAP0028\",\"28\",\"CHILL\"",
  "\"1.0.0\",\"AETERM\",\"CHILLS\",\"CHILLS\",\"CHILLS\",\"Chills\",\"CHILLS\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",29,\"MAP0029\",\"29\",\"CHILLS\"",
  "\"1.0.0\",\"AETERM\",\"RIGORS\",\"RIGORS\",\"CHILLS\",\"Chills\",\"CHILLS\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",30,\"MAP0030\",\"30\",\"RIGORS\"",
  "\"1.0.0\",\"AEDECOD\",\"Allergic cough\",\"ALLERGIC COUGH\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",31,\"MAP0031\",\"31\",\"Allergic cough\"",
  "\"1.0.0\",\"AEDECOD\",\"Cough\",\"COUGH\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",32,\"MAP0032\",\"32\",\"Cough\"",
  "\"1.0.0\",\"AEDECOD\",\"Productive cough\",\"PRODUCTIVE COUGH\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",33,\"MAP0033\",\"33\",\"Productive cough\"",
  "\"1.0.0\",\"AELLT\",\"Allergic cough\",\"ALLERGIC COUGH\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",34,\"MAP0034\",\"34\",\"Allergic cough\"",
  "\"1.0.0\",\"AELLT\",\"Chronic cough\",\"CHRONIC COUGH\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",35,\"MAP0035\",\"35\",\"Chronic cough\"",
  "\"1.0.0\",\"AELLT\",\"Cough\",\"COUGH\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",36,\"MAP0036\",\"36\",\"Cough\"",
  "\"1.0.0\",\"AELLT\",\"Cough aggravated\",\"COUGH AGGRAVATED\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",37,\"MAP0037\",\"37\",\"Cough aggravated\"",
  "\"1.0.0\",\"AELLT\",\"Cough increased\",\"COUGH INCREASED\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",38,\"MAP0038\",\"38\",\"Cough increased\"",
  "\"1.0.0\",\"AELLT\",\"Cough nonproductive\",\"COUGH NONPRODUCTIVE\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",39,\"MAP0039\",\"39\",\"Cough nonproductive\"",
  "\"1.0.0\",\"AELLT\",\"Coughing\",\"COUGHING\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",40,\"MAP0040\",\"40\",\"Coughing\"",
  "\"1.0.0\",\"AELLT\",\"Drug-induced cough\",\"DRUG-INDUCED COUGH\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",41,\"MAP0041\",\"41\",\"Drug-induced cough\"",
  "\"1.0.0\",\"AELLT\",\"Dry cough\",\"DRY COUGH\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",42,\"MAP0042\",\"42\",\"Dry cough\"",
  "\"1.0.0\",\"AELLT\",\"Nocturnal cough\",\"NOCTURNAL COUGH\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",43,\"MAP0043\",\"43\",\"Nocturnal cough\"",
  "\"1.0.0\",\"AELLT\",\"Persistent cough\",\"PERSISTENT COUGH\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",44,\"MAP0044\",\"44\",\"Persistent cough\"",
  "\"1.0.0\",\"AELLT\",\"Productive cough\",\"PRODUCTIVE COUGH\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",45,\"MAP0045\",\"45\",\"Productive cough\"",
  "\"1.0.0\",\"AELLT\",\"Smoker's cough\",\"SMOKER'S COUGH\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",46,\"MAP0046\",\"46\",\"Smoker's cough\"",
  "\"1.0.0\",\"AETERM\",\"CHRONIC COUGH\",\"CHRONIC COUGH\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",47,\"MAP0047\",\"47\",\"CHRONIC COUGH\"",
  "\"1.0.0\",\"AETERM\",\"COUGH\",\"COUGH\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",48,\"MAP0048\",\"48\",\"COUGH\"",
  "\"1.0.0\",\"AETERM\",\"COUGHING\",\"COUGHING\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",49,\"MAP0049\",\"49\",\"COUGHING\"",
  "\"1.0.0\",\"AETERM\",\"DRY COUGH\",\"DRY COUGH\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",50,\"MAP0050\",\"50\",\"DRY COUGH\"",
  "\"1.0.0\",\"AETERM\",\"PRODUCTIVE COUGH\",\"PRODUCTIVE COUGH\",\"COUGH\",\"Cough\",\"COUGH\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",51,\"MAP0051\",\"51\",\"PRODUCTIVE COUGH\"",
  "\"1.0.0\",\"AEDECOD\",\"Dyspnoea\",\"DYSPNOEA\",\"SHBREATH\",\"Shortness of breath\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",52,\"MAP0052\",\"52\",\"Dyspnoea\"",
  "\"1.0.0\",\"AEDECOD\",\"Dyspnoea exertional\",\"DYSPNOEA EXERTIONAL\",\"SHBREATH\",\"Shortness of breath\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",53,\"MAP0053\",\"53\",\"Dyspnoea exertional\"",
  "\"1.0.0\",\"AELLT\",\"Dyspnea\",\"DYSPNEA\",\"SHBREATH\",\"Shortness of breath\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",54,\"MAP0054\",\"54\",\"Dyspnea\"",
  "\"1.0.0\",\"AELLT\",\"Dyspnea exacerbated\",\"DYSPNEA EXACERBATED\",\"SHBREATH\",\"Shortness of breath\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",55,\"MAP0055\",\"55\",\"Dyspnea exacerbated\"",
  "\"1.0.0\",\"AELLT\",\"Dyspnea exertional\",\"DYSPNEA EXERTIONAL\",\"SHBREATH\",\"Shortness of breath\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",56,\"MAP0056\",\"56\",\"Dyspnea exertional\"",
  "\"1.0.0\",\"AELLT\",\"Exertional dyspnea\",\"EXERTIONAL DYSPNEA\",\"SHBREATH\",\"Shortness of breath\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",57,\"MAP0057\",\"57\",\"Exertional dyspnea\"",
  "\"1.0.0\",\"AELLT\",\"Increased shortness of breath\",\"INCREASED SHORTNESS OF BREATH\",\"SHBREATH\",\"Shortness of breath\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",58,\"MAP0058\",\"58\",\"Increased shortness of breath\"",
  "\"1.0.0\",\"AELLT\",\"Short of breath\",\"SHORT OF BREATH\",\"SHBREATH\",\"Shortness of breath\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",59,\"MAP0059\",\"59\",\"Short of breath\"",
  "\"1.0.0\",\"AELLT\",\"Shortness of breath\",\"SHORTNESS OF BREATH\",\"SHBREATH\",\"Shortness of breath\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",60,\"MAP0060\",\"60\",\"Shortness of breath\"",
  "\"1.0.0\",\"AETERM\",\"DYSPNEA\",\"DYSPNEA\",\"SHBREATH\",\"Shortness of breath\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",61,\"MAP0061\",\"61\",\"DYSPNEA\"",
  "\"1.0.0\",\"AETERM\",\"EXERTIONAL DYSPNEA\",\"EXERTIONAL DYSPNEA\",\"SHBREATH\",\"Shortness of breath\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",62,\"MAP0062\",\"62\",\"EXERTIONAL DYSPNEA\"",
  "\"1.0.0\",\"AETERM\",\"INCREASED SHORTNESS OF BREATH\",\"INCREASED SHORTNESS OF BREATH\",\"SHBREATH\",\"Shortness of breath\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",63,\"MAP0063\",\"63\",\"INCREASED SHORTNESS OF BREATH\"",
  "\"1.0.0\",\"AETERM\",\"SHORT OF BREATH\",\"SHORT OF BREATH\",\"SHBREATH\",\"Shortness of breath\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",64,\"MAP0064\",\"64\",\"SHORT OF BREATH\"",
  "\"1.0.0\",\"AETERM\",\"SHORTNESS OF BREATH\",\"SHORTNESS OF BREATH\",\"SHBREATH\",\"Shortness of breath\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",65,\"MAP0065\",\"65\",\"SHORTNESS OF BREATH\"",
  "\"1.0.0\",\"AELLT\",\"Breathing difficult\",\"BREATHING DIFFICULT\",\"DIBREATH\",\"Difficulty breathing\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",66,\"MAP0066\",\"66\",\"Breathing difficult\"",
  "\"1.0.0\",\"AELLT\",\"Difficulty breathing\",\"DIFFICULTY BREATHING\",\"DIBREATH\",\"Difficulty breathing\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",67,\"MAP0067\",\"67\",\"Difficulty breathing\"",
  "\"1.0.0\",\"AELLT\",\"Labored breathing\",\"LABORED BREATHING\",\"DIBREATH\",\"Difficulty breathing\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",68,\"MAP0068\",\"68\",\"Labored breathing\"",
  "\"1.0.0\",\"AETERM\",\"DIFFICULT BREATHING\",\"DIFFICULT BREATHING\",\"DIBREATH\",\"Difficulty breathing\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",69,\"MAP0069\",\"69\",\"DIFFICULT BREATHING\"",
  "\"1.0.0\",\"AETERM\",\"DIFFICULTY BREATHING\",\"DIFFICULTY BREATHING\",\"DIBREATH\",\"Difficulty breathing\",\"DYSPNEA\",\"DIRECT\",TRUE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",70,\"MAP0070\",\"70\",\"DIFFICULTY BREATHING\"",
  "\"1.0.0\",\"AEDECOD\",\"Respiratory distress\",\"RESPIRATORY DISTRESS\",\"DIBREATH\",\"Difficulty breathing\",\"DYSPNEA\",\"RELATED_REVIEW\",FALSE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Respiratory finding requiring review; do not infer the checklist symptom automatically.\",71,\"MAP0071\",\"71\",\"Respiratory distress\"",
  "\"1.0.0\",\"AELLT\",\"Respiratory distress\",\"RESPIRATORY DISTRESS\",\"DIBREATH\",\"Difficulty breathing\",\"DYSPNEA\",\"RELATED_REVIEW\",FALSE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Respiratory finding requiring review; do not infer the checklist symptom automatically.\",72,\"MAP0072\",\"72\",\"Respiratory distress\"",
  "\"1.0.0\",\"AETERM\",\"RESPIRATORY DISTRESS\",\"RESPIRATORY DISTRESS\",\"DIBREATH\",\"Difficulty breathing\",\"DYSPNEA\",\"RELATED_REVIEW\",FALSE,\"0\",\"Any-duration symptom; review onset and post-vaccination context.\",\"Respiratory finding requiring review; do not infer the checklist symptom automatically.\",73,\"MAP0073\",\"73\",\"RESPIRATORY DISTRESS\"",
  "\"1.0.0\",\"AEDECOD\",\"Fatigue\",\"FATIGUE\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",74,\"MAP0074\",\"74\",\"Fatigue\"",
  "\"1.0.0\",\"AELLT\",\"Fatigue\",\"FATIGUE\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",75,\"MAP0075\",\"75\",\"Fatigue\"",
  "\"1.0.0\",\"AELLT\",\"Fatigue aggravated\",\"FATIGUE AGGRAVATED\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",76,\"MAP0076\",\"76\",\"Fatigue aggravated\"",
  "\"1.0.0\",\"AELLT\",\"Fatigue extreme\",\"FATIGUE EXTREME\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",77,\"MAP0077\",\"77\",\"Fatigue extreme\"",
  "\"1.0.0\",\"AELLT\",\"Tiredness\",\"TIREDNESS\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",78,\"MAP0078\",\"78\",\"Tiredness\"",
  "\"1.0.0\",\"AETERM\",\"FATIGUE\",\"FATIGUE\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",79,\"MAP0079\",\"79\",\"FATIGUE\"",
  "\"1.0.0\",\"AETERM\",\"TIREDNESS\",\"TIREDNESS\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",80,\"MAP0080\",\"80\",\"TIREDNESS\"",
  "\"1.0.0\",\"AEDECOD\",\"Asthenia\",\"ASTHENIA\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Weakness, malaise or localized/mental fatigue is broader than the checklist symptom.\",81,\"MAP0081\",\"81\",\"Asthenia\"",
  "\"1.0.0\",\"AEDECOD\",\"Malaise\",\"MALAISE\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Weakness, malaise or localized/mental fatigue is broader than the checklist symptom.\",82,\"MAP0082\",\"82\",\"Malaise\"",
  "\"1.0.0\",\"AEDECOD\",\"Mental fatigue\",\"MENTAL FATIGUE\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Weakness, malaise or localized/mental fatigue is broader than the checklist symptom.\",83,\"MAP0083\",\"83\",\"Mental fatigue\"",
  "\"1.0.0\",\"AEDECOD\",\"Muscle fatigue\",\"MUSCLE FATIGUE\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Weakness, malaise or localized/mental fatigue is broader than the checklist symptom.\",84,\"MAP0084\",\"84\",\"Muscle fatigue\"",
  "\"1.0.0\",\"AELLT\",\"Asthenia\",\"ASTHENIA\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Weakness, malaise or localized/mental fatigue is broader than the checklist symptom.\",85,\"MAP0085\",\"85\",\"Asthenia\"",
  "\"1.0.0\",\"AELLT\",\"Malaise\",\"MALAISE\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Weakness, malaise or localized/mental fatigue is broader than the checklist symptom.\",86,\"MAP0086\",\"86\",\"Malaise\"",
  "\"1.0.0\",\"AELLT\",\"Mental fatigue\",\"MENTAL FATIGUE\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Weakness, malaise or localized/mental fatigue is broader than the checklist symptom.\",87,\"MAP0087\",\"87\",\"Mental fatigue\"",
  "\"1.0.0\",\"AELLT\",\"Muscle fatigue\",\"MUSCLE FATIGUE\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Weakness, malaise or localized/mental fatigue is broader than the checklist symptom.\",88,\"MAP0088\",\"88\",\"Muscle fatigue\"",
  "\"1.0.0\",\"AELLT\",\"Weakness\",\"WEAKNESS\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Weakness, malaise or localized/mental fatigue is broader than the checklist symptom.\",89,\"MAP0089\",\"89\",\"Weakness\"",
  "\"1.0.0\",\"AELLT\",\"Weakness generalized\",\"WEAKNESS GENERALIZED\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Weakness, malaise or localized/mental fatigue is broader than the checklist symptom.\",90,\"MAP0090\",\"90\",\"Weakness generalized\"",
  "\"1.0.0\",\"AETERM\",\"GENERALIZED WEAKNESS\",\"GENERALIZED WEAKNESS\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Weakness, malaise or localized/mental fatigue is broader than the checklist symptom.\",91,\"MAP0091\",\"91\",\"GENERALIZED WEAKNESS\"",
  "\"1.0.0\",\"AETERM\",\"MALAISE\",\"MALAISE\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Weakness, malaise or localized/mental fatigue is broader than the checklist symptom.\",92,\"MAP0092\",\"92\",\"MALAISE\"",
  "\"1.0.0\",\"AETERM\",\"MENTAL FATIGUE\",\"MENTAL FATIGUE\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Weakness, malaise or localized/mental fatigue is broader than the checklist symptom.\",93,\"MAP0093\",\"93\",\"MENTAL FATIGUE\"",
  "\"1.0.0\",\"AETERM\",\"MUSCLE FATIGUE\",\"MUSCLE FATIGUE\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Weakness, malaise or localized/mental fatigue is broader than the checklist symptom.\",94,\"MAP0094\",\"94\",\"MUSCLE FATIGUE\"",
  "\"1.0.0\",\"AETERM\",\"WEAKNESS\",\"WEAKNESS\",\"FATIGUE\",\"Fatigue\",\"FATIGUE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Weakness, malaise or localized/mental fatigue is broader than the checklist symptom.\",95,\"MAP0095\",\"95\",\"WEAKNESS\"",
  "\"1.0.0\",\"AEDECOD\",\"Myalgia\",\"MYALGIA\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",96,\"MAP0096\",\"96\",\"Myalgia\"",
  "\"1.0.0\",\"AELLT\",\"Generalised muscle aches\",\"GENERALISED MUSCLE ACHES\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",97,\"MAP0097\",\"97\",\"Generalised muscle aches\"",
  "\"1.0.0\",\"AELLT\",\"Generalized muscle aches\",\"GENERALIZED MUSCLE ACHES\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",98,\"MAP0098\",\"98\",\"Generalized muscle aches\"",
  "\"1.0.0\",\"AELLT\",\"Localised muscle pain\",\"LOCALISED MUSCLE PAIN\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",99,\"MAP0099\",\"99\",\"Localised muscle pain\"",
  "\"1.0.0\",\"AELLT\",\"Localized muscle pain\",\"LOCALIZED MUSCLE PAIN\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",100,\"MAP0100\",\"100\",\"Localized muscle pain\"",
  "\"1.0.0\",\"AELLT\",\"Muscle ache\",\"MUSCLE ACHE\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",101,\"MAP0101\",\"101\",\"Muscle ache\"",
  "\"1.0.0\",\"AELLT\",\"Muscle pain\",\"MUSCLE PAIN\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",102,\"MAP0102\",\"102\",\"Muscle pain\"",
  "\"1.0.0\",\"AELLT\",\"Muscle soreness\",\"MUSCLE SORENESS\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",103,\"MAP0103\",\"103\",\"Muscle soreness\"",
  "\"1.0.0\",\"AELLT\",\"Muscular pain\",\"MUSCULAR PAIN\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",104,\"MAP0104\",\"104\",\"Muscular pain\"",
  "\"1.0.0\",\"AELLT\",\"Muscular pains\",\"MUSCULAR PAINS\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",105,\"MAP0105\",\"105\",\"Muscular pains\"",
  "\"1.0.0\",\"AELLT\",\"Myalgia\",\"MYALGIA\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",106,\"MAP0106\",\"106\",\"Myalgia\"",
  "\"1.0.0\",\"AELLT\",\"Myalgia aggravated\",\"MYALGIA AGGRAVATED\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",107,\"MAP0107\",\"107\",\"Myalgia aggravated\"",
  "\"1.0.0\",\"AELLT\",\"Myalgia of lower extremities\",\"MYALGIA OF LOWER EXTREMITIES\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",108,\"MAP0108\",\"108\",\"Myalgia of lower extremities\"",
  "\"1.0.0\",\"AELLT\",\"Myalgia upper extremities\",\"MYALGIA UPPER EXTREMITIES\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",109,\"MAP0109\",\"109\",\"Myalgia upper extremities\"",
  "\"1.0.0\",\"AELLT\",\"Pain muscle\",\"PAIN MUSCLE\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",110,\"MAP0110\",\"110\",\"Pain muscle\"",
  "\"1.0.0\",\"AELLT\",\"Polymyalgia\",\"POLYMYALGIA\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",111,\"MAP0111\",\"111\",\"Polymyalgia\"",
  "\"1.0.0\",\"AELLT\",\"Tenderness muscle\",\"TENDERNESS MUSCLE\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",112,\"MAP0112\",\"112\",\"Tenderness muscle\"",
  "\"1.0.0\",\"AETERM\",\"GENERALIZED MUSCLE ACHES\",\"GENERALIZED MUSCLE ACHES\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",113,\"MAP0113\",\"113\",\"GENERALIZED MUSCLE ACHES\"",
  "\"1.0.0\",\"AETERM\",\"MUSCLE  ACHES\",\"MUSCLE ACHES\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",114,\"MAP0114\",\"114;116\",\"MUSCLE  ACHES | MUSCLE ACHES\"",
  "\"1.0.0\",\"AETERM\",\"MUSCLE ACHE\",\"MUSCLE ACHE\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",115,\"MAP0115\",\"115\",\"MUSCLE ACHE\"",
  "\"1.0.0\",\"AETERM\",\"MUSCLE PAIN\",\"MUSCLE PAIN\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",117,\"MAP0116\",\"117\",\"MUSCLE PAIN\"",
  "\"1.0.0\",\"AETERM\",\"MUSCLE SORENESS\",\"MUSCLE SORENESS\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",118,\"MAP0117\",\"118\",\"MUSCLE SORENESS\"",
  "\"1.0.0\",\"AETERM\",\"MYALGIA\",\"MYALGIA\",\"MYALGIA\",\"Muscle aches (myalgia)\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",119,\"MAP0118\",\"119\",\"MYALGIA\"",
  "\"1.0.0\",\"AELLT\",\"General body pain\",\"GENERAL BODY PAIN\",\"BODYACHE\",\"Body aches\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",120,\"MAP0119\",\"120\",\"General body pain\"",
  "\"1.0.0\",\"AELLT\",\"Generalized aching\",\"GENERALIZED ACHING\",\"BODYACHE\",\"Body aches\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",121,\"MAP0120\",\"121\",\"Generalized aching\"",
  "\"1.0.0\",\"AETERM\",\"(WHOLE BODY) BODY ACHES\",\"(WHOLE BODY) BODY ACHES\",\"BODYACHE\",\"Body aches\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",122,\"MAP0121\",\"122\",\"(WHOLE BODY) BODY ACHES\"",
  "\"1.0.0\",\"AETERM\",\"ALL OVER BODY ACHES\",\"ALL OVER BODY ACHES\",\"BODYACHE\",\"Body aches\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",123,\"MAP0122\",\"123\",\"ALL OVER BODY ACHES\"",
  "\"1.0.0\",\"AETERM\",\"BODY ACHE\",\"BODY ACHE\",\"BODYACHE\",\"Body aches\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",124,\"MAP0123\",\"124\",\"BODY ACHE\"",
  "\"1.0.0\",\"AETERM\",\"BODY ACHES\",\"BODY ACHES\",\"BODYACHE\",\"Body aches\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",125,\"MAP0124\",\"125\",\"BODY ACHES\"",
  "\"1.0.0\",\"AETERM\",\"GENERAL BODY ACHES\",\"GENERAL BODY ACHES\",\"BODYACHE\",\"Body aches\",\"MYALGIA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",126,\"MAP0125\",\"126\",\"GENERAL BODY ACHES\"",
  "\"1.0.0\",\"AETERM\",\"GENERALIZED  PAIN\",\"GENERALIZED PAIN\",\"BODYACHE\",\"Body aches\",\"MYALGIA\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"General pain is not necessarily the protocol's body-ache symptom; inspect verbatim description.\",127,\"MAP0126\",\"127;128\",\"GENERALIZED  PAIN | GENERALIZED PAIN\"",
  "\"1.0.0\",\"AEDECOD\",\"Cluster headache\",\"CLUSTER HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",129,\"MAP0127\",\"129\",\"Cluster headache\"",
  "\"1.0.0\",\"AEDECOD\",\"Exertional headache\",\"EXERTIONAL HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",130,\"MAP0128\",\"130\",\"Exertional headache\"",
  "\"1.0.0\",\"AEDECOD\",\"Headache\",\"HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",131,\"MAP0129\",\"131\",\"Headache\"",
  "\"1.0.0\",\"AEDECOD\",\"New daily persistent headache\",\"NEW DAILY PERSISTENT HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",132,\"MAP0130\",\"132\",\"New daily persistent headache\"",
  "\"1.0.0\",\"AEDECOD\",\"Post-traumatic headache\",\"POST-TRAUMATIC HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",133,\"MAP0131\",\"133\",\"Post-traumatic headache\"",
  "\"1.0.0\",\"AEDECOD\",\"Primary cough headache\",\"PRIMARY COUGH HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",134,\"MAP0132\",\"134\",\"Primary cough headache\"",
  "\"1.0.0\",\"AEDECOD\",\"Primary headache associated with sexual activity\",\"PRIMARY HEADACHE ASSOCIATED WITH SEXUAL ACTIVITY\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",135,\"MAP0133\",\"135\",\"Primary headache associated with sexual activity\"",
  "\"1.0.0\",\"AEDECOD\",\"Procedural headache\",\"PROCEDURAL HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",136,\"MAP0134\",\"136\",\"Procedural headache\"",
  "\"1.0.0\",\"AEDECOD\",\"Sinus headache\",\"SINUS HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",137,\"MAP0135\",\"137\",\"Sinus headache\"",
  "\"1.0.0\",\"AEDECOD\",\"Tension headache\",\"TENSION HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",138,\"MAP0136\",\"138\",\"Tension headache\"",
  "\"1.0.0\",\"AELLT\",\"Cluster headaches\",\"CLUSTER HEADACHES\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",139,\"MAP0137\",\"139\",\"Cluster headaches\"",
  "\"1.0.0\",\"AELLT\",\"Coital headache\",\"COITAL HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",140,\"MAP0138\",\"140\",\"Coital headache\"",
  "\"1.0.0\",\"AELLT\",\"Drug-induced headache\",\"DRUG-INDUCED HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",141,\"MAP0139\",\"141\",\"Drug-induced headache\"",
  "\"1.0.0\",\"AELLT\",\"Exertional headache\",\"EXERTIONAL HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",142,\"MAP0140\",\"142\",\"Exertional headache\"",
  "\"1.0.0\",\"AELLT\",\"Forehead headache\",\"FOREHEAD HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",143,\"MAP0141\",\"143\",\"Forehead headache\"",
  "\"1.0.0\",\"AELLT\",\"Frequent headaches\",\"FREQUENT HEADACHES\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",144,\"MAP0142\",\"144\",\"Frequent headaches\"",
  "\"1.0.0\",\"AELLT\",\"Frontal headache\",\"FRONTAL HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",145,\"MAP0143\",\"145\",\"Frontal headache\"",
  "\"1.0.0\",\"AELLT\",\"Head pain\",\"HEAD PAIN\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",146,\"MAP0144\",\"146\",\"Head pain\"",
  "\"1.0.0\",\"AELLT\",\"Headache\",\"HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",147,\"MAP0145\",\"147\",\"Headache\"",
  "\"1.0.0\",\"AELLT\",\"Headache (excl migraine) aggravated\",\"HEADACHE (EXCL MIGRAINE) AGGRAVATED\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",148,\"MAP0146\",\"148\",\"Headache (excl migraine) aggravated\"",
  "\"1.0.0\",\"AELLT\",\"Headache aggravated\",\"HEADACHE AGGRAVATED\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",149,\"MAP0147\",\"149\",\"Headache aggravated\"",
  "\"1.0.0\",\"AELLT\",\"Headache dull\",\"HEADACHE DULL\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",150,\"MAP0148\",\"150\",\"Headache dull\"",
  "\"1.0.0\",\"AELLT\",\"Headache post-traumatic\",\"HEADACHE POST-TRAUMATIC\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",151,\"MAP0149\",\"151\",\"Headache post-traumatic\"",
  "\"1.0.0\",\"AELLT\",\"Headache postoperative\",\"HEADACHE POSTOPERATIVE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",152,\"MAP0150\",\"152\",\"Headache postoperative\"",
  "\"1.0.0\",\"AELLT\",\"Headache recurrent\",\"HEADACHE RECURRENT\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",153,\"MAP0151\",\"153\",\"Headache recurrent\"",
  "\"1.0.0\",\"AELLT\",\"Headache sinus\",\"HEADACHE SINUS\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",154,\"MAP0152\",\"154\",\"Headache sinus\"",
  "\"1.0.0\",\"AELLT\",\"Headache temporal\",\"HEADACHE TEMPORAL\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",155,\"MAP0153\",\"155\",\"Headache temporal\"",
  "\"1.0.0\",\"AELLT\",\"Headache tension\",\"HEADACHE TENSION\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",156,\"MAP0154\",\"156\",\"Headache tension\"",
  "\"1.0.0\",\"AELLT\",\"Headache unilateral\",\"HEADACHE UNILATERAL\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",157,\"MAP0155\",\"157\",\"Headache unilateral\"",
  "\"1.0.0\",\"AELLT\",\"Intermittent headache\",\"INTERMITTENT HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",158,\"MAP0156\",\"158\",\"Intermittent headache\"",
  "\"1.0.0\",\"AELLT\",\"Migraine headache\",\"MIGRAINE HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",159,\"MAP0157\",\"159\",\"Migraine headache\"",
  "\"1.0.0\",\"AELLT\",\"New daily persistent headache\",\"NEW DAILY PERSISTENT HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",160,\"MAP0158\",\"160\",\"New daily persistent headache\"",
  "\"1.0.0\",\"AELLT\",\"Occipital headache\",\"OCCIPITAL HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",161,\"MAP0159\",\"161\",\"Occipital headache\"",
  "\"1.0.0\",\"AELLT\",\"Pain head\",\"PAIN HEAD\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",162,\"MAP0160\",\"162\",\"Pain head\"",
  "\"1.0.0\",\"AELLT\",\"Post-traumatic headache\",\"POST-TRAUMATIC HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",163,\"MAP0161\",\"163\",\"Post-traumatic headache\"",
  "\"1.0.0\",\"AELLT\",\"Post procedural headache\",\"POST PROCEDURAL HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",164,\"MAP0162\",\"164\",\"Post procedural headache\"",
  "\"1.0.0\",\"AELLT\",\"Primary cough headache\",\"PRIMARY COUGH HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",165,\"MAP0163\",\"165\",\"Primary cough headache\"",
  "\"1.0.0\",\"AELLT\",\"Sinus headache\",\"SINUS HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",166,\"MAP0164\",\"166\",\"Sinus headache\"",
  "\"1.0.0\",\"AELLT\",\"Temporal headache\",\"TEMPORAL HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",167,\"MAP0165\",\"167\",\"Temporal headache\"",
  "\"1.0.0\",\"AELLT\",\"Tension headache\",\"TENSION HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",168,\"MAP0166\",\"168\",\"Tension headache\"",
  "\"1.0.0\",\"AELLT\",\"Tension headaches\",\"TENSION HEADACHES\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",169,\"MAP0167\",\"169\",\"Tension headaches\"",
  "\"1.0.0\",\"AETERM\",\"CLUSTER HEADACHES\",\"CLUSTER HEADACHES\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",170,\"MAP0168\",\"170\",\"CLUSTER HEADACHES\"",
  "\"1.0.0\",\"AETERM\",\"FRONTAL HEADACHE\",\"FRONTAL HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",171,\"MAP0169\",\"171\",\"FRONTAL HEADACHE\"",
  "\"1.0.0\",\"AETERM\",\"HEAD PAIN\",\"HEAD PAIN\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",172,\"MAP0170\",\"172\",\"HEAD PAIN\"",
  "\"1.0.0\",\"AETERM\",\"HEADACHE\",\"HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",173,\"MAP0171\",\"173\",\"HEADACHE\"",
  "\"1.0.0\",\"AETERM\",\"INTERMITTENT HEADACHE\",\"INTERMITTENT HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",174,\"MAP0172\",\"174\",\"INTERMITTENT HEADACHE\"",
  "\"1.0.0\",\"AETERM\",\"MIGRAINE HEADACHE\",\"MIGRAINE HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",175,\"MAP0173\",\"175\",\"MIGRAINE HEADACHE\"",
  "\"1.0.0\",\"AETERM\",\"OCCIPITAL HEADACHE\",\"OCCIPITAL HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",176,\"MAP0174\",\"176\",\"OCCIPITAL HEADACHE\"",
  "\"1.0.0\",\"AETERM\",\"SINUS HEADACHE\",\"SINUS HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",177,\"MAP0175\",\"177\",\"SINUS HEADACHE\"",
  "\"1.0.0\",\"AETERM\",\"TENSION HEADACHE\",\"TENSION HEADACHE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",178,\"MAP0176\",\"178\",\"TENSION HEADACHE\"",
  "\"1.0.0\",\"AETERM\",\"TENSION HEADACHES\",\"TENSION HEADACHES\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",179,\"MAP0177\",\"179\",\"TENSION HEADACHES\"",
  "\"1.0.0\",\"AEDECOD\",\"Migraine\",\"MIGRAINE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"A migraine diagnosis, particularly aura, does not establish headache on a particular day.\",180,\"MAP0178\",\"180\",\"Migraine\"",
  "\"1.0.0\",\"AEDECOD\",\"Migraine with aura\",\"MIGRAINE WITH AURA\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"A migraine diagnosis, particularly aura, does not establish headache on a particular day.\",181,\"MAP0179\",\"181\",\"Migraine with aura\"",
  "\"1.0.0\",\"AEDECOD\",\"Migraine without aura\",\"MIGRAINE WITHOUT AURA\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"A migraine diagnosis, particularly aura, does not establish headache on a particular day.\",182,\"MAP0180\",\"182\",\"Migraine without aura\"",
  "\"1.0.0\",\"AEDECOD\",\"Vestibular migraine\",\"VESTIBULAR MIGRAINE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"A migraine diagnosis, particularly aura, does not establish headache on a particular day.\",183,\"MAP0181\",\"183\",\"Vestibular migraine\"",
  "\"1.0.0\",\"AELLT\",\"Complicated migraine\",\"COMPLICATED MIGRAINE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"A migraine diagnosis, particularly aura, does not establish headache on a particular day.\",184,\"MAP0182\",\"184\",\"Complicated migraine\"",
  "\"1.0.0\",\"AELLT\",\"Migraine\",\"MIGRAINE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"A migraine diagnosis, particularly aura, does not establish headache on a particular day.\",185,\"MAP0183\",\"185\",\"Migraine\"",
  "\"1.0.0\",\"AELLT\",\"Migraine aggravated\",\"MIGRAINE AGGRAVATED\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"A migraine diagnosis, particularly aura, does not establish headache on a particular day.\",186,\"MAP0184\",\"186\",\"Migraine aggravated\"",
  "\"1.0.0\",\"AELLT\",\"Migraine aura\",\"MIGRAINE AURA\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"A migraine diagnosis, particularly aura, does not establish headache on a particular day.\",187,\"MAP0185\",\"187\",\"Migraine aura\"",
  "\"1.0.0\",\"AELLT\",\"Migraine with aura\",\"MIGRAINE WITH AURA\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"A migraine diagnosis, particularly aura, does not establish headache on a particular day.\",188,\"MAP0186\",\"188\",\"Migraine with aura\"",
  "\"1.0.0\",\"AELLT\",\"Migraine without aura\",\"MIGRAINE WITHOUT AURA\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"A migraine diagnosis, particularly aura, does not establish headache on a particular day.\",189,\"MAP0187\",\"189\",\"Migraine without aura\"",
  "\"1.0.0\",\"AELLT\",\"Ophthalmic migraine\",\"OPHTHALMIC MIGRAINE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"A migraine diagnosis, particularly aura, does not establish headache on a particular day.\",190,\"MAP0188\",\"190\",\"Ophthalmic migraine\"",
  "\"1.0.0\",\"AELLT\",\"Vestibular migraine\",\"VESTIBULAR MIGRAINE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"A migraine diagnosis, particularly aura, does not establish headache on a particular day.\",191,\"MAP0189\",\"191\",\"Vestibular migraine\"",
  "\"1.0.0\",\"AETERM\",\"COMPLICATED MIGRAINE\",\"COMPLICATED MIGRAINE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"A migraine diagnosis, particularly aura, does not establish headache on a particular day.\",192,\"MAP0190\",\"192\",\"COMPLICATED MIGRAINE\"",
  "\"1.0.0\",\"AETERM\",\"MIGRAINE\",\"MIGRAINE\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"A migraine diagnosis, particularly aura, does not establish headache on a particular day.\",193,\"MAP0191\",\"193\",\"MIGRAINE\"",
  "\"1.0.0\",\"AETERM\",\"MIGRAINE AURA\",\"MIGRAINE AURA\",\"HEADACHE\",\"Headache\",\"HEADACHE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"A migraine diagnosis, particularly aura, does not establish headache on a particular day.\",194,\"MAP0192\",\"194\",\"MIGRAINE AURA\"",
  "\"1.0.0\",\"AEDECOD\",\"Anosmia\",\"ANOSMIA\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"DIRECT\",TRUE,\"48\",\"Verify new onset and duration >=48 h.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",195,\"MAP0193\",\"195\",\"Anosmia\"",
  "\"1.0.0\",\"AELLT\",\"Anosmia\",\"ANOSMIA\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"DIRECT\",TRUE,\"48\",\"Verify new onset and duration >=48 h.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",196,\"MAP0194\",\"196\",\"Anosmia\"",
  "\"1.0.0\",\"AELLT\",\"Loss of smell\",\"LOSS OF SMELL\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"DIRECT\",TRUE,\"48\",\"Verify new onset and duration >=48 h.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",197,\"MAP0195\",\"197\",\"Loss of smell\"",
  "\"1.0.0\",\"AETERM\",\"ANOSMIA\",\"ANOSMIA\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"DIRECT\",TRUE,\"48\",\"Verify new onset and duration >=48 h.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",198,\"MAP0196\",\"198\",\"ANOSMIA\"",
  "\"1.0.0\",\"AETERM\",\"LOSS OF SMELL\",\"LOSS OF SMELL\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"DIRECT\",TRUE,\"48\",\"Verify new onset and duration >=48 h.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",199,\"MAP0197\",\"199\",\"LOSS OF SMELL\"",
  "\"1.0.0\",\"AETERM\",\"LOSS OF TASTE AND SMELL\",\"LOSS OF TASTE AND SMELL\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"DIRECT\",TRUE,\"48\",\"Verify new onset and duration >=48 h.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",200,\"MAP0198\",\"200\",\"LOSS OF TASTE AND SMELL\"",
  "\"1.0.0\",\"AETERM\",\"NEW LOSS OF SMELL\",\"NEW LOSS OF SMELL\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"DIRECT\",TRUE,\"48\",\"Verify new onset and duration >=48 h.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",201,\"MAP0199\",\"201\",\"NEW LOSS OF SMELL\"",
  "\"1.0.0\",\"AEDECOD\",\"Hyposmia\",\"HYPOSMIA\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or distorted smell is not synonymous with new loss; review onset and symptom detail.\",202,\"MAP0200\",\"202\",\"Hyposmia\"",
  "\"1.0.0\",\"AEDECOD\",\"Parosmia\",\"PAROSMIA\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or distorted smell is not synonymous with new loss; review onset and symptom detail.\",203,\"MAP0201\",\"203\",\"Parosmia\"",
  "\"1.0.0\",\"AELLT\",\"Altered smell sensation\",\"ALTERED SMELL SENSATION\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or distorted smell is not synonymous with new loss; review onset and symptom detail.\",204,\"MAP0202\",\"204\",\"Altered smell sensation\"",
  "\"1.0.0\",\"AELLT\",\"Diminished sense of smell\",\"DIMINISHED SENSE OF SMELL\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or distorted smell is not synonymous with new loss; review onset and symptom detail.\",205,\"MAP0203\",\"205\",\"Diminished sense of smell\"",
  "\"1.0.0\",\"AELLT\",\"Hyposmia\",\"HYPOSMIA\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or distorted smell is not synonymous with new loss; review onset and symptom detail.\",206,\"MAP0204\",\"206\",\"Hyposmia\"",
  "\"1.0.0\",\"AELLT\",\"Phantosmia\",\"PHANTOSMIA\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or distorted smell is not synonymous with new loss; review onset and symptom detail.\",207,\"MAP0205\",\"207\",\"Phantosmia\"",
  "\"1.0.0\",\"AELLT\",\"Smell alteration\",\"SMELL ALTERATION\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or distorted smell is not synonymous with new loss; review onset and symptom detail.\",208,\"MAP0206\",\"208\",\"Smell alteration\"",
  "\"1.0.0\",\"AELLT\",\"Smell change\",\"SMELL CHANGE\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or distorted smell is not synonymous with new loss; review onset and symptom detail.\",209,\"MAP0207\",\"209\",\"Smell change\"",
  "\"1.0.0\",\"AETERM\",\"DIMINISHED SENSE OF SMELL\",\"DIMINISHED SENSE OF SMELL\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or distorted smell is not synonymous with new loss; review onset and symptom detail.\",210,\"MAP0208\",\"210\",\"DIMINISHED SENSE OF SMELL\"",
  "\"1.0.0\",\"AETERM\",\"HYPOSMIA\",\"HYPOSMIA\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or distorted smell is not synonymous with new loss; review onset and symptom detail.\",211,\"MAP0209\",\"211\",\"HYPOSMIA\"",
  "\"1.0.0\",\"AETERM\",\"PHANTOSMIA\",\"PHANTOSMIA\",\"NLSMELL\",\"New loss of smell\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or distorted smell is not synonymous with new loss; review onset and symptom detail.\",212,\"MAP0210\",\"212\",\"PHANTOSMIA\"",
  "\"1.0.0\",\"AEDECOD\",\"Ageusia\",\"AGEUSIA\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"DIRECT\",TRUE,\"48\",\"Verify new onset and duration >=48 h.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",213,\"MAP0211\",\"213\",\"Ageusia\"",
  "\"1.0.0\",\"AELLT\",\"Ageusia\",\"AGEUSIA\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"DIRECT\",TRUE,\"48\",\"Verify new onset and duration >=48 h.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",214,\"MAP0212\",\"214\",\"Ageusia\"",
  "\"1.0.0\",\"AELLT\",\"Loss of taste\",\"LOSS OF TASTE\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"DIRECT\",TRUE,\"48\",\"Verify new onset and duration >=48 h.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",215,\"MAP0213\",\"215\",\"Loss of taste\"",
  "\"1.0.0\",\"AETERM\",\"AGEUSIA\",\"AGEUSIA\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"DIRECT\",TRUE,\"48\",\"Verify new onset and duration >=48 h.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",216,\"MAP0214\",\"216\",\"AGEUSIA\"",
  "\"1.0.0\",\"AETERM\",\"LOSS OF TASTE\",\"LOSS OF TASTE\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"DIRECT\",TRUE,\"48\",\"Verify new onset and duration >=48 h.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",217,\"MAP0215\",\"217\",\"LOSS OF TASTE\"",
  "\"1.0.0\",\"AETERM\",\"LOSS OF TASTE AND SMELL\",\"LOSS OF TASTE AND SMELL\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"DIRECT\",TRUE,\"48\",\"Verify new onset and duration >=48 h.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",218,\"MAP0216\",\"218\",\"LOSS OF TASTE AND SMELL\"",
  "\"1.0.0\",\"AETERM\",\"NEW LOSS OF TASTE\",\"NEW LOSS OF TASTE\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"DIRECT\",TRUE,\"48\",\"Verify new onset and duration >=48 h.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",219,\"MAP0217\",\"219\",\"NEW LOSS OF TASTE\"",
  "\"1.0.0\",\"AEDECOD\",\"Dysgeusia\",\"DYSGEUSIA\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or altered taste is not synonymous with new loss; review onset and symptom detail.\",220,\"MAP0218\",\"220\",\"Dysgeusia\"",
  "\"1.0.0\",\"AEDECOD\",\"Hypogeusia\",\"HYPOGEUSIA\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or altered taste is not synonymous with new loss; review onset and symptom detail.\",221,\"MAP0219\",\"221\",\"Hypogeusia\"",
  "\"1.0.0\",\"AEDECOD\",\"Taste disorder\",\"TASTE DISORDER\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or altered taste is not synonymous with new loss; review onset and symptom detail.\",222,\"MAP0220\",\"222\",\"Taste disorder\"",
  "\"1.0.0\",\"AELLT\",\"Dysgeusia\",\"DYSGEUSIA\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or altered taste is not synonymous with new loss; review onset and symptom detail.\",223,\"MAP0221\",\"223\",\"Dysgeusia\"",
  "\"1.0.0\",\"AELLT\",\"Hypogeusia\",\"HYPOGEUSIA\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or altered taste is not synonymous with new loss; review onset and symptom detail.\",224,\"MAP0222\",\"224\",\"Hypogeusia\"",
  "\"1.0.0\",\"AELLT\",\"Impaired taste\",\"IMPAIRED TASTE\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or altered taste is not synonymous with new loss; review onset and symptom detail.\",225,\"MAP0223\",\"225\",\"Impaired taste\"",
  "\"1.0.0\",\"AELLT\",\"Taste alteration\",\"TASTE ALTERATION\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or altered taste is not synonymous with new loss; review onset and symptom detail.\",226,\"MAP0224\",\"226\",\"Taste alteration\"",
  "\"1.0.0\",\"AELLT\",\"Taste altered\",\"TASTE ALTERED\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or altered taste is not synonymous with new loss; review onset and symptom detail.\",227,\"MAP0225\",\"227\",\"Taste altered\"",
  "\"1.0.0\",\"AELLT\",\"Taste bitter\",\"TASTE BITTER\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or altered taste is not synonymous with new loss; review onset and symptom detail.\",228,\"MAP0226\",\"228\",\"Taste bitter\"",
  "\"1.0.0\",\"AELLT\",\"Taste changed\",\"TASTE CHANGED\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or altered taste is not synonymous with new loss; review onset and symptom detail.\",229,\"MAP0227\",\"229\",\"Taste changed\"",
  "\"1.0.0\",\"AELLT\",\"Taste diminished\",\"TASTE DIMINISHED\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or altered taste is not synonymous with new loss; review onset and symptom detail.\",230,\"MAP0228\",\"230\",\"Taste diminished\"",
  "\"1.0.0\",\"AELLT\",\"Taste metallic\",\"TASTE METALLIC\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or altered taste is not synonymous with new loss; review onset and symptom detail.\",231,\"MAP0229\",\"231\",\"Taste metallic\"",
  "\"1.0.0\",\"AELLT\",\"Taste salty\",\"TASTE SALTY\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or altered taste is not synonymous with new loss; review onset and symptom detail.\",232,\"MAP0230\",\"232\",\"Taste salty\"",
  "\"1.0.0\",\"AETERM\",\"DYSGEUSIA\",\"DYSGEUSIA\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or altered taste is not synonymous with new loss; review onset and symptom detail.\",233,\"MAP0231\",\"233\",\"DYSGEUSIA\"",
  "\"1.0.0\",\"AETERM\",\"HYPOGEUSIA\",\"HYPOGEUSIA\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or altered taste is not synonymous with new loss; review onset and symptom detail.\",234,\"MAP0232\",\"234\",\"HYPOGEUSIA\"",
  "\"1.0.0\",\"AETERM\",\"TASTE ALTERATION\",\"TASTE ALTERATION\",\"NLTASTE\",\"New loss of taste\",\"SMELL_TASTE\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify new onset and duration >=48 h.\",\"Reduced or altered taste is not synonymous with new loss; review onset and symptom detail.\",235,\"MAP0233\",\"235\",\"TASTE ALTERATION\"",
  "\"1.0.0\",\"AEDECOD\",\"Oropharyngeal pain\",\"OROPHARYNGEAL PAIN\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",236,\"MAP0234\",\"236\",\"Oropharyngeal pain\"",
  "\"1.0.0\",\"AELLT\",\"Chronic sore throat\",\"CHRONIC SORE THROAT\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",237,\"MAP0235\",\"237\",\"Chronic sore throat\"",
  "\"1.0.0\",\"AELLT\",\"Oropharyngeal pain\",\"OROPHARYNGEAL PAIN\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",238,\"MAP0236\",\"238\",\"Oropharyngeal pain\"",
  "\"1.0.0\",\"AELLT\",\"Pain throat\",\"PAIN THROAT\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",239,\"MAP0237\",\"239\",\"Pain throat\"",
  "\"1.0.0\",\"AELLT\",\"Sore throat\",\"SORE THROAT\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",240,\"MAP0238\",\"240\",\"Sore throat\"",
  "\"1.0.0\",\"AELLT\",\"Streptococcal sore throat\",\"STREPTOCOCCAL SORE THROAT\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",241,\"MAP0239\",\"241\",\"Streptococcal sore throat\"",
  "\"1.0.0\",\"AELLT\",\"Throat pain\",\"THROAT PAIN\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",242,\"MAP0240\",\"242\",\"Throat pain\"",
  "\"1.0.0\",\"AELLT\",\"Throat sore\",\"THROAT SORE\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",243,\"MAP0241\",\"243\",\"Throat sore\"",
  "\"1.0.0\",\"AETERM\",\"SORE  THROAT\",\"SORE THROAT\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",244,\"MAP0242\",\"244;245\",\"SORE  THROAT | SORE THROAT\"",
  "\"1.0.0\",\"AETERM\",\"THROAT PAIN\",\"THROAT PAIN\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",246,\"MAP0243\",\"246\",\"THROAT PAIN\"",
  "\"1.0.0\",\"AEDECOD\",\"Dry throat\",\"DRY THROAT\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",247,\"MAP0244\",\"247\",\"Dry throat\"",
  "\"1.0.0\",\"AEDECOD\",\"Oropharyngeal discomfort\",\"OROPHARYNGEAL DISCOMFORT\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",248,\"MAP0245\",\"248\",\"Oropharyngeal discomfort\"",
  "\"1.0.0\",\"AEDECOD\",\"Pharyngeal paraesthesia\",\"PHARYNGEAL PARAESTHESIA\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",249,\"MAP0246\",\"249\",\"Pharyngeal paraesthesia\"",
  "\"1.0.0\",\"AEDECOD\",\"Throat irritation\",\"THROAT IRRITATION\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",250,\"MAP0247\",\"250\",\"Throat irritation\"",
  "\"1.0.0\",\"AEDECOD\",\"Throat tightness\",\"THROAT TIGHTNESS\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",251,\"MAP0248\",\"251\",\"Throat tightness\"",
  "\"1.0.0\",\"AELLT\",\"Dry throat\",\"DRY THROAT\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",252,\"MAP0249\",\"252\",\"Dry throat\"",
  "\"1.0.0\",\"AELLT\",\"Itchy throat\",\"ITCHY THROAT\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",253,\"MAP0250\",\"253\",\"Itchy throat\"",
  "\"1.0.0\",\"AELLT\",\"Oropharyngeal discomfort\",\"OROPHARYNGEAL DISCOMFORT\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",254,\"MAP0251\",\"254\",\"Oropharyngeal discomfort\"",
  "\"1.0.0\",\"AELLT\",\"Pharynx irritated sensation of\",\"PHARYNX IRRITATED SENSATION OF\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",255,\"MAP0252\",\"255\",\"Pharynx irritated sensation of\"",
  "\"1.0.0\",\"AELLT\",\"Throat discomfort\",\"THROAT DISCOMFORT\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",256,\"MAP0253\",\"256\",\"Throat discomfort\"",
  "\"1.0.0\",\"AELLT\",\"Throat irritation\",\"THROAT IRRITATION\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",257,\"MAP0254\",\"257\",\"Throat irritation\"",
  "\"1.0.0\",\"AELLT\",\"Throat tightness\",\"THROAT TIGHTNESS\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",258,\"MAP0255\",\"258\",\"Throat tightness\"",
  "\"1.0.0\",\"AELLT\",\"Tingling throat\",\"TINGLING THROAT\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",259,\"MAP0256\",\"259\",\"Tingling throat\"",
  "\"1.0.0\",\"AETERM\",\"DRY THROAT\",\"DRY THROAT\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",260,\"MAP0257\",\"260\",\"DRY THROAT\"",
  "\"1.0.0\",\"AETERM\",\"ITCHY THROAT\",\"ITCHY THROAT\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",261,\"MAP0258\",\"261\",\"ITCHY THROAT\"",
  "\"1.0.0\",\"AETERM\",\"THROAT DISCOMFORT\",\"THROAT DISCOMFORT\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",262,\"MAP0259\",\"262\",\"THROAT DISCOMFORT\"",
  "\"1.0.0\",\"AETERM\",\"THROAT IRRITATION\",\"THROAT IRRITATION\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",263,\"MAP0260\",\"263\",\"THROAT IRRITATION\"",
  "\"1.0.0\",\"AETERM\",\"THROAT TIGHTNESS\",\"THROAT TIGHTNESS\",\"STHROAT\",\"Sore throat\",\"SORE_THROAT\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Irritation, dryness, discomfort or tightness does not necessarily mean sore throat.\",264,\"MAP0261\",\"264\",\"THROAT TIGHTNESS\"",
  "\"1.0.0\",\"AEDECOD\",\"Nasal congestion\",\"NASAL CONGESTION\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",265,\"MAP0262\",\"265\",\"Nasal congestion\"",
  "\"1.0.0\",\"AELLT\",\"Chronic nasal congestion\",\"CHRONIC NASAL CONGESTION\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",266,\"MAP0263\",\"266\",\"Chronic nasal congestion\"",
  "\"1.0.0\",\"AELLT\",\"Congestion nasal\",\"CONGESTION NASAL\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",267,\"MAP0264\",\"267\",\"Congestion nasal\"",
  "\"1.0.0\",\"AELLT\",\"Nasal congestion\",\"NASAL CONGESTION\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",268,\"MAP0265\",\"268\",\"Nasal congestion\"",
  "\"1.0.0\",\"AELLT\",\"Nasal sinus congestion\",\"NASAL SINUS CONGESTION\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",269,\"MAP0266\",\"269\",\"Nasal sinus congestion\"",
  "\"1.0.0\",\"AELLT\",\"Nasal stuffiness\",\"NASAL STUFFINESS\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",270,\"MAP0267\",\"270\",\"Nasal stuffiness\"",
  "\"1.0.0\",\"AELLT\",\"Nose congestion\",\"NOSE CONGESTION\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",271,\"MAP0268\",\"271\",\"Nose congestion\"",
  "\"1.0.0\",\"AETERM\",\"CONGESTION NASAL\",\"CONGESTION NASAL\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",272,\"MAP0269\",\"272\",\"CONGESTION NASAL\"",
  "\"1.0.0\",\"AETERM\",\"NASAL  CONGESTION\",\"NASAL CONGESTION\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",273,\"MAP0270\",\"273;274\",\"NASAL  CONGESTION | NASAL CONGESTION\"",
  "\"1.0.0\",\"AETERM\",\"NOSE CONGESTION\",\"NOSE CONGESTION\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",275,\"MAP0271\",\"275\",\"NOSE CONGESTION\"",
  "\"1.0.0\",\"AEDECOD\",\"Respiratory tract congestion\",\"RESPIRATORY TRACT CONGESTION\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Congestion site or equivalence to the nasal checklist item needs review.\",276,\"MAP0272\",\"276\",\"Respiratory tract congestion\"",
  "\"1.0.0\",\"AEDECOD\",\"Sinus congestion\",\"SINUS CONGESTION\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Congestion site or equivalence to the nasal checklist item needs review.\",277,\"MAP0273\",\"277\",\"Sinus congestion\"",
  "\"1.0.0\",\"AEDECOD\",\"Upper respiratory tract congestion\",\"UPPER RESPIRATORY TRACT CONGESTION\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Congestion site or equivalence to the nasal checklist item needs review.\",278,\"MAP0274\",\"278\",\"Upper respiratory tract congestion\"",
  "\"1.0.0\",\"AELLT\",\"Respiratory tract congestion\",\"RESPIRATORY TRACT CONGESTION\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Congestion site or equivalence to the nasal checklist item needs review.\",279,\"MAP0275\",\"279\",\"Respiratory tract congestion\"",
  "\"1.0.0\",\"AELLT\",\"Sinus congestion\",\"SINUS CONGESTION\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Congestion site or equivalence to the nasal checklist item needs review.\",280,\"MAP0276\",\"280\",\"Sinus congestion\"",
  "\"1.0.0\",\"AELLT\",\"Upper respiratory tract congestion\",\"UPPER RESPIRATORY TRACT CONGESTION\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Congestion site or equivalence to the nasal checklist item needs review.\",281,\"MAP0277\",\"281\",\"Upper respiratory tract congestion\"",
  "\"1.0.0\",\"AETERM\",\"RESPIRATORY TRACT CONGESTION\",\"RESPIRATORY TRACT CONGESTION\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Congestion site or equivalence to the nasal checklist item needs review.\",282,\"MAP0278\",\"282\",\"RESPIRATORY TRACT CONGESTION\"",
  "\"1.0.0\",\"AETERM\",\"SINUS CONGESTION\",\"SINUS CONGESTION\",\"NASALCON\",\"Nasal congestion\",\"NASAL\",\"RELATED_REVIEW\",FALSE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Congestion site or equivalence to the nasal checklist item needs review.\",283,\"MAP0279\",\"283\",\"SINUS CONGESTION\"",
  "\"1.0.0\",\"AEDECOD\",\"Rhinorrhoea\",\"RHINORRHOEA\",\"RUNNOSE\",\"Runny nose (rhinorrhea)\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",284,\"MAP0280\",\"284\",\"Rhinorrhoea\"",
  "\"1.0.0\",\"AELLT\",\"Nasal discharge\",\"NASAL DISCHARGE\",\"RUNNOSE\",\"Runny nose (rhinorrhea)\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",285,\"MAP0281\",\"285\",\"Nasal discharge\"",
  "\"1.0.0\",\"AELLT\",\"Nasal mucus increased\",\"NASAL MUCUS INCREASED\",\"RUNNOSE\",\"Runny nose (rhinorrhea)\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",286,\"MAP0282\",\"286\",\"Nasal mucus increased\"",
  "\"1.0.0\",\"AELLT\",\"Rhinorrhea\",\"RHINORRHEA\",\"RUNNOSE\",\"Runny nose (rhinorrhea)\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",287,\"MAP0283\",\"287\",\"Rhinorrhea\"",
  "\"1.0.0\",\"AELLT\",\"Rhinorrhoea\",\"RHINORRHOEA\",\"RUNNOSE\",\"Runny nose (rhinorrhea)\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",288,\"MAP0284\",\"288\",\"Rhinorrhoea\"",
  "\"1.0.0\",\"AELLT\",\"Runny nose\",\"RUNNY NOSE\",\"RUNNOSE\",\"Runny nose (rhinorrhea)\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",289,\"MAP0285\",\"289\",\"Runny nose\"",
  "\"1.0.0\",\"AETERM\",\"NASAL DISCHARGE\",\"NASAL DISCHARGE\",\"RUNNOSE\",\"Runny nose (rhinorrhea)\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",290,\"MAP0286\",\"290\",\"NASAL DISCHARGE\"",
  "\"1.0.0\",\"AETERM\",\"RHINORRHEA\",\"RHINORRHEA\",\"RUNNOSE\",\"Runny nose (rhinorrhea)\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",291,\"MAP0287\",\"291\",\"RHINORRHEA\"",
  "\"1.0.0\",\"AETERM\",\"RUNNY  NOSE\",\"RUNNY NOSE\",\"RUNNOSE\",\"Runny nose (rhinorrhea)\",\"NASAL\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",292,\"MAP0288\",\"292;293\",\"RUNNY  NOSE | RUNNY NOSE\"",
  "\"1.0.0\",\"AEDECOD\",\"Nausea\",\"NAUSEA\",\"NAUSEA\",\"Nausea\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",294,\"MAP0289\",\"294\",\"Nausea\"",
  "\"1.0.0\",\"AEDECOD\",\"Procedural nausea\",\"PROCEDURAL NAUSEA\",\"NAUSEA\",\"Nausea\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",295,\"MAP0290\",\"295\",\"Procedural nausea\"",
  "\"1.0.0\",\"AELLT\",\"Nausea\",\"NAUSEA\",\"NAUSEA\",\"Nausea\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",296,\"MAP0291\",\"296\",\"Nausea\"",
  "\"1.0.0\",\"AELLT\",\"Nausea aggravated\",\"NAUSEA AGGRAVATED\",\"NAUSEA\",\"Nausea\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",297,\"MAP0292\",\"297\",\"Nausea aggravated\"",
  "\"1.0.0\",\"AELLT\",\"Post procedural nausea\",\"POST PROCEDURAL NAUSEA\",\"NAUSEA\",\"Nausea\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",298,\"MAP0293\",\"298\",\"Post procedural nausea\"",
  "\"1.0.0\",\"AELLT\",\"Postoperative nausea\",\"POSTOPERATIVE NAUSEA\",\"NAUSEA\",\"Nausea\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",299,\"MAP0294\",\"299\",\"Postoperative nausea\"",
  "\"1.0.0\",\"AELLT\",\"Procedural nausea\",\"PROCEDURAL NAUSEA\",\"NAUSEA\",\"Nausea\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",300,\"MAP0295\",\"300\",\"Procedural nausea\"",
  "\"1.0.0\",\"AETERM\",\"NAUSEA\",\"NAUSEA\",\"NAUSEA\",\"Nausea\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",301,\"MAP0296\",\"301\",\"NAUSEA\"",
  "\"1.0.0\",\"AETERM\",\"NAUSEA/VOMITING\",\"NAUSEA/VOMITING\",\"NAUSEA\",\"Nausea\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",302,\"MAP0297\",\"302\",\"NAUSEA/VOMITING\"",
  "\"1.0.0\",\"AETERM\",\"POSTOPERATIVE NAUSEA\",\"POSTOPERATIVE NAUSEA\",\"NAUSEA\",\"Nausea\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",303,\"MAP0298\",\"303\",\"POSTOPERATIVE NAUSEA\"",
  "\"1.0.0\",\"AEDECOD\",\"Procedural vomiting\",\"PROCEDURAL VOMITING\",\"VOMIT\",\"Vomiting\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",304,\"MAP0299\",\"304\",\"Procedural vomiting\"",
  "\"1.0.0\",\"AEDECOD\",\"Vomiting\",\"VOMITING\",\"VOMIT\",\"Vomiting\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",305,\"MAP0300\",\"305\",\"Vomiting\"",
  "\"1.0.0\",\"AELLT\",\"Emesis\",\"EMESIS\",\"VOMIT\",\"Vomiting\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",306,\"MAP0301\",\"306\",\"Emesis\"",
  "\"1.0.0\",\"AELLT\",\"Hyperemesis\",\"HYPEREMESIS\",\"VOMIT\",\"Vomiting\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",307,\"MAP0302\",\"307\",\"Hyperemesis\"",
  "\"1.0.0\",\"AELLT\",\"Procedural vomiting\",\"PROCEDURAL VOMITING\",\"VOMIT\",\"Vomiting\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",308,\"MAP0303\",\"308\",\"Procedural vomiting\"",
  "\"1.0.0\",\"AELLT\",\"Vomited\",\"VOMITED\",\"VOMIT\",\"Vomiting\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",309,\"MAP0304\",\"309\",\"Vomited\"",
  "\"1.0.0\",\"AELLT\",\"Vomiting\",\"VOMITING\",\"VOMIT\",\"Vomiting\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",310,\"MAP0305\",\"310\",\"Vomiting\"",
  "\"1.0.0\",\"AETERM\",\"EMESIS\",\"EMESIS\",\"VOMIT\",\"Vomiting\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",311,\"MAP0306\",\"311\",\"EMESIS\"",
  "\"1.0.0\",\"AETERM\",\"HYPEREMESIS\",\"HYPEREMESIS\",\"VOMIT\",\"Vomiting\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",312,\"MAP0307\",\"312\",\"HYPEREMESIS\"",
  "\"1.0.0\",\"AETERM\",\"NAUSEA/VOMITING\",\"NAUSEA/VOMITING\",\"VOMIT\",\"Vomiting\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",313,\"MAP0308\",\"313\",\"NAUSEA/VOMITING\"",
  "\"1.0.0\",\"AETERM\",\"VOMITED\",\"VOMITED\",\"VOMIT\",\"Vomiting\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",314,\"MAP0309\",\"314\",\"VOMITED\"",
  "\"1.0.0\",\"AETERM\",\"VOMITING\",\"VOMITING\",\"VOMIT\",\"Vomiting\",\"NAUSEA_VOMITING\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",315,\"MAP0310\",\"315\",\"VOMITING\"",
  "\"1.0.0\",\"AEDECOD\",\"Bacterial diarrhoea\",\"BACTERIAL DIARRHOEA\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",316,\"MAP0311\",\"316\",\"Bacterial diarrhoea\"",
  "\"1.0.0\",\"AEDECOD\",\"Diarrhoea\",\"DIARRHOEA\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",317,\"MAP0312\",\"317\",\"Diarrhoea\"",
  "\"1.0.0\",\"AEDECOD\",\"Diarrhoea infectious\",\"DIARRHOEA INFECTIOUS\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",318,\"MAP0313\",\"318\",\"Diarrhoea infectious\"",
  "\"1.0.0\",\"AEDECOD\",\"Viral diarrhoea\",\"VIRAL DIARRHOEA\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",319,\"MAP0314\",\"319\",\"Viral diarrhoea\"",
  "\"1.0.0\",\"AELLT\",\"Acute diarrhea\",\"ACUTE DIARRHEA\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",320,\"MAP0315\",\"320\",\"Acute diarrhea\"",
  "\"1.0.0\",\"AELLT\",\"Bacterial diarrhea\",\"BACTERIAL DIARRHEA\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",321,\"MAP0316\",\"321\",\"Bacterial diarrhea\"",
  "\"1.0.0\",\"AELLT\",\"Diarrhea\",\"DIARRHEA\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",322,\"MAP0317\",\"322\",\"Diarrhea\"",
  "\"1.0.0\",\"AELLT\",\"Diarrhea aggravated\",\"DIARRHEA AGGRAVATED\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",323,\"MAP0318\",\"323\",\"Diarrhea aggravated\"",
  "\"1.0.0\",\"AELLT\",\"Diarrhoea\",\"DIARRHOEA\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",324,\"MAP0319\",\"324\",\"Diarrhoea\"",
  "\"1.0.0\",\"AELLT\",\"Loose bowel\",\"LOOSE BOWEL\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",325,\"MAP0320\",\"325\",\"Loose bowel\"",
  "\"1.0.0\",\"AELLT\",\"Loose bowels\",\"LOOSE BOWELS\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",326,\"MAP0321\",\"326\",\"Loose bowels\"",
  "\"1.0.0\",\"AELLT\",\"Loose stools\",\"LOOSE STOOLS\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",327,\"MAP0322\",\"327\",\"Loose stools\"",
  "\"1.0.0\",\"AELLT\",\"Stools watery\",\"STOOLS WATERY\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",328,\"MAP0323\",\"328\",\"Stools watery\"",
  "\"1.0.0\",\"AELLT\",\"Traveller's diarrhea\",\"TRAVELLER'S DIARRHEA\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",329,\"MAP0324\",\"329\",\"Traveller's diarrhea\"",
  "\"1.0.0\",\"AELLT\",\"Viral diarrhea\",\"VIRAL DIARRHEA\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",330,\"MAP0325\",\"330\",\"Viral diarrhea\"",
  "\"1.0.0\",\"AELLT\",\"Watery diarrhea\",\"WATERY DIARRHEA\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",331,\"MAP0326\",\"331\",\"Watery diarrhea\"",
  "\"1.0.0\",\"AETERM\",\"ACUTE DIARRHEA\",\"ACUTE DIARRHEA\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",332,\"MAP0327\",\"332\",\"ACUTE DIARRHEA\"",
  "\"1.0.0\",\"AETERM\",\"DIARRHEA\",\"DIARRHEA\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",333,\"MAP0328\",\"333\",\"DIARRHEA\"",
  "\"1.0.0\",\"AETERM\",\"LOOSE BOWELS\",\"LOOSE BOWELS\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",334,\"MAP0329\",\"334\",\"LOOSE BOWELS\"",
  "\"1.0.0\",\"AETERM\",\"LOOSE STOOLS\",\"LOOSE STOOLS\",\"DIARRHEA\",\"Diarrhea\",\"DIARRHEA\",\"DIRECT\",TRUE,\"48\",\"Verify duration >=48 h and post-vaccination context.\",\"Explicit symptom or lexical synonym; cause is not inferred.\",335,\"MAP0330\",\"335\",\"LOOSE STOOLS\""
),collapse="\n"),colClasses="character",check.names=FALSE)

# ---- Symptom-conditioned denominators and temporally linked PCR ----
pt_symptom_inputs <- function(root,release,cfg,submission,ctx,out) {
 dest<-file.path(out,'symptoms');dir.create(dest,showWarnings=FALSE)
 if(!is.null(cfg$EVENTS)&&nzchar(cfg$EVENTS)) {
  path<-ws_path(root,cfg$EVENTS)
  e<-read.csv(path,colClasses='character',check.names=FALSE,fileEncoding='UTF-8-BOM',na.strings=c('','NA'))
  e[]<-lapply(e,ws_text)
  req<-c('subject','source','source_id','symptom_code','symptom_group','lo','hi','date_quality','event_id')
  if(!all(req%in%names(e)))stop('EVENTS must be the complete symptom_event_evidence.csv export; missing: ',paste(setdiff(req,names(e)),collapse=', '))
  e$subject<-ws_sid(e$subject);e$source<-ws_norm(e$source)
  if(any(!e$source%in%c('AE','COVID_FORM')))stop('Unexpected symptom export source')
  if(any(!e$symptom_code%in%names(ws_param_map)))stop('Unexpected symptom code in export')
  if(any(e$symptom_group!=unname(ws_param_map[e$symptom_code])))stop('Inconsistent symptom groups in export')
  # Previous weekly exports carry ISO start_date/end_date alongside numeric lo/hi.
  # Heatmap exports use ISO lo/hi. Bare numeric values without a documented origin
  # are rejected rather than guessed.
  e$lo<-ws_dates(if('start_date'%in%names(e))e$start_date else e$lo,'event start',PHMPT_NUMERIC_DATE_ORIGIN)
  e$hi<-ws_dates(if('end_date'%in%names(e))e$end_date else e$hi,'event end',PHMPT_NUMERIC_DATE_ORIGIN)
  if(anyDuplicated(e$event_id))stop('Duplicate event IDs in EVENTS export')
  manifest<-data.frame(release=release,domain='SYMPTOM_EVENTS_EXPORT',path=normalizePath(path,winslash='/'),md5=unname(tools::md5sum(path)),status='SELECTED',rows=nrow(e),detail='Explicit complete prior event export; upstream matching decisions retained, not re-derived')
 } else {
  ae<-pt_pick(root,release,'AE',cfg$AE,dest,submission)
  sy<-pt_pick(root,release,'ADSYMP',cfg$ADSYMP,dest,submission)
  ae$data<-pt_unique(ae$data,c('subject','AESEQ'),'AE',dest)
  mapping<-ws_mapping(root,PHMPT_MAPPING_FILE,dest)
  e<-ws_events(ae$data,sy$data,mapping,PHMPT_NUMERIC_DATE_ORIGIN,dest,'interval')
  manifest<-rbind(ae$manifest,sy$manifest)
 }
 e$evidence_date<-e$lo
 e$anchor_basis<-ifelse(e$source=='AE','AE onset (AESTDTC)','Positive assessment date (ADSYMP ADT)')
 e$usable_anchor<-is.finite(e$lo)&is.finite(e$hi)&e$hi>=e$lo
 e$in_randomized_roster<-e$subject%in%ctx$subject
 pt_dates_write(e,file.path(dest,'symptom_event_evidence_used.csv'))
 pt_dates_write(e[!e$usable_anchor|!e$in_randomized_roster,,drop=FALSE],file.path(dest,'symptom_evidence_for_review.csv'))
 list(events=e,manifest=manifest)
}
pt_symptom_window <- function(ctx,p,e,cutoff,scope,release,period,out) {
 dest<-file.path(out,'symptom_links',paste(release,period,scope,sep='_'));dir.create(dest,recursive=TRUE,showWarnings=FALSE)
 c<-ctx;c$start<-c$RANDDT;c$end<-pmin(cutoff,if(scope=='STRICT_BLINDED')c$blinded_end else c$calendar_end)
 context_ok<-c$start<=c$end
 if(scope=='STRICT_BLINDED')context_ok<-context_ok&(c$strict_context_usable|cutoff<=ws_dates(PHMPT_EUA_CUTOFF,'EUA'))
 c<-c[context_ok,,drop=FALSE]
 j<-match(e$subject,c$subject)
 ok<-!is.na(j)&e$usable_anchor&e$evidence_date>=c$start[j]&e$evidence_date<=c$end[j];ok[is.na(ok)]<-FALSE
 ex<-e[!ok,,drop=FALSE]
 ex$reason<-ifelse(!ex$usable_anchor,'UNKNOWN_OR_INVALID_DATE',ifelse(!ex$subject%in%c$subject,'NO_ELIGIBLE_RANDOMIZED_CONTEXT','OUTSIDE_WINDOW'))
 pt_dates_write(ex,file.path(dest,'excluded_symptom_evidence.csv'))
 v<-e[ok,,drop=FALSE]
 a<-unique(v[c('subject','evidence_date')]);a<-a[order(a$subject,a$evidence_date),,drop=FALSE]
 a$anchor_id<-paste(a$subject,ws_iso(a$evidence_date),sep='|')
 vk<-paste(v$subject,ws_iso(v$evidence_date),sep='|');v$anchor_id<-vk
 a$has_covid_form<-a$anchor_id%in%vk[v$source=='COVID_FORM'];a$has_ae<-a$anchor_id%in%vk[v$source=='AE']
 j<-match(a$subject,c$subject);a$arm<-c$arm[j];a$site<-c$site[j];a$window_end<-c$end[j]
 a$link_end<-pmin(a$evidence_date+PHMPT_TEST_LINK_DAYS,a$window_end)
 a$complete_followup<-a$evidence_date+PHMPT_TEST_LINK_DAYS<=a$window_end
 # Days 0-6 after either dose are retained. Do not impose a PCR-positive condition.
 a$early_postdose<-a$evidence_date>=c$TR01SDT[j]&a$evidence_date<=c$TR01SDT[j]+6 |
  a$evidence_date>=c$DOSE2DT[j]&a$evidence_date<=c$DOSE2DT[j]+6
 a$early_postdose[is.na(a$early_postdose)]<-FALSE
 q<-p[p$usable_dated&!p$baseline,,drop=FALSE]
 if(scope=='STRICT_BLINDED')q<-q[!q$boundary_visit,,drop=FALSE]
 q$diagnostic<-q$illness_labelled|q$MBGRPID%in%PHMPT_LOCAL_GROUPS|
  q$local_flag_supp%in%c('Y','YES')|q$local_flag_admb%in%c('Y','YES')
 ix<-split(seq_len(nrow(q)),q$subject);links<-list();n<-0L
 for(i in seq_len(nrow(a))) {
  k<-ix[[a$subject[i]]];if(!length(k))next
  k<-k[q$date[k]>=a$evidence_date[i]&q$date[k]<=a$link_end[i]]
  if(!length(k))next
  n<-n+1L;links[[n]]<-data.frame(anchor_id=a$anchor_id[i],subject=a$subject[i],evidence_date=a$evidence_date[i],
   record_id=q$record_id[k],test_date=q$date[k],delay_days=q$date[k]-a$evidence_date[i],route=q$route[k],
   diagnostic=q$diagnostic[k],raw_source_row=q$raw_source_row[k],stringsAsFactors=FALSE)
 }
 h<-if(length(links))do.call(rbind,links) else data.frame(anchor_id=character(),subject=character(),evidence_date=numeric(),record_id=character(),test_date=numeric(),delay_days=numeric(),route=character(),diagnostic=logical(),raw_source_row=integer())
 if(nrow(h)){nt<-table(h$record_id);h$anchors_sharing_record<-as.integer(nt[h$record_id])}else h$anchors_sharing_record<-integer()
 pt_dates_write(a,file.path(dest,'symptom_evidence_dates.csv'))
 pt_dates_write(v,file.path(dest,'symptom_record_to_date_links.csv'))
 pt_dates_write(h,file.path(dest,'symptom_date_to_pcr_links.csv'))
 list(anchors=a,links=h,context=c)
}
pt_conditional_subjects <- function(w,source,opportunity,test_selection,release,period,scope,out) {
 a<-w$anchors;c<-w$context;h<-w$links
 base<-if(source=='COVID_FORM_ONLY')a[a$has_covid_form,,drop=FALSE] else a
 selected<-if(opportunity=='COMPLETE_7D')base[base$complete_followup,,drop=FALSE] else base
 h<-h[h$anchor_id%in%selected$anchor_id,,drop=FALSE]
 if(test_selection=='DIAGNOSTIC')h<-h[h$diagnostic,,drop=FALSE]
 subjects<-unique(selected$subject);x<-c[match(subjects,c$subject),,drop=FALSE]
 n<-nrow(x);x$days_observed<-x$end-x$start+1
 for(r in c('ANY','LOCAL_SOURCE','CENTRAL_SOURCE','UNCLASSIFIED','CONFLICT')) {
  z<-if(r=='ANY')h else h[h$route==r,,drop=FALSE]
  # A PCR can lie near multiple assessments. Count the record once per person,
  # never once per linked symptom/AE, and count each person once in inference.
  rec<-unique(z[c('subject','record_id')]);days<-unique(z[c('subject','test_date')]);nm<-tolower(r)
  x[[paste0(nm,'_records')]]<-tabulate(match(rec$subject,x$subject),nbins=n)
  x[[paste0(nm,'_days')]]<-tabulate(match(days$subject,x$subject),nbins=n)
  x[[paste0(nm,'_tested')]]<-x[[paste0(nm,'_records')]]>0
 }
 x$both_tested<-x$local_source_tested&x$central_source_tested
 x$unresolved_tested<-x$unclassified_tested|x$conflict_tested
 x$n_evidence_dates<-tabulate(match(selected$subject,x$subject),nbins=n)
 x$n_complete_dates<-tabulate(match(selected$subject[selected$complete_followup],x$subject),nbins=n)
 x$has_covid_form_evidence<-x$subject%in%selected$subject[selected$has_covid_form]
 x$has_ae_evidence<-x$subject%in%selected$subject[selected$has_ae]
 x$source_category<-ifelse(x$has_covid_form_evidence&x$has_ae_evidence,'Both sources',ifelse(x$has_covid_form_evidence,'COVID forms only','AE only'))
 stratum<-paste('SYMPTOMS',source,opportunity,test_selection,sep='__')
 x$release<-rep(release,n);x$period<-rep(period,n);x$scope<-rep(scope,n);x$stratum<-rep(stratum,n)
 keep<-unique(c('subject','arm','site','start','end','days_observed','source_category','n_evidence_dates','n_complete_dates',
  'has_covid_form_evidence','has_ae_evidence',names(x)[grepl('_(records|days|tested)$',names(x))],'release','period','scope','stratum'))
 x<-x[keep]
 pt_dates_write(x,file.path(out,'symptom_subjects',paste(release,period,scope,paste0(stratum,'.csv'),sep='_')))
 meta<-data.frame(release=release,period=period,scope=scope,stratum=stratum)
 s<-pt_summary(x,meta=meta,sites=sort(unique(w$context$site)))
 s$person_days<-s$records_per_100_person_days<-NA_real_ # no symptom person-time rate is estimated
 s$denominator_source<-source;s$followup_requirement<-opportunity;s$pcr_selection<-test_selection
 f<-do.call(rbind,lapply(c('mRNA-1273','Placebo'),function(ar){z<-base[base$arm==ar,,drop=FALSE];u<-selected[selected$arm==ar,,drop=FALSE]
  data.frame(release=release,period=period,scope=scope,denominator_source=source,followup_requirement=opportunity,pcr_selection=test_selection,arm=ar,
   subjects_with_dated_evidence=length(unique(z$subject)),subjects_with_complete_followup=length(unique(z$subject[z$complete_followup])),
   denominator_subjects=length(unique(u$subject)),evidence_dates=nrow(z),incomplete_dates=sum(!z$complete_followup),selected_dates=nrow(u),
   subjects_excluded_by_followup=length(setdiff(unique(z$subject),unique(u$subject))))}))
 list(summary=s,flow=f)
}
pt_symptom_report <- function(s,t,f,out) {
 header<-c('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>PCR testing among subjects reporting symptoms</title>',
 '<style>body{font:16px/1.55 Georgia,serif;color:#17212a;margin:0;background:#f7f7f4}main{max-width:1200px;padding:30px;margin:auto;background:white}h1,h2,h3,summary{font-family:Arial,sans-serif;line-height:1.25}h1{font-size:29px}h2{font-size:22px}.note{padding:16px;background:#eef4f7;border-left:4px solid #42637e}table{border-collapse:collapse;width:100%;font:12px/1.4 Arial,sans-serif}th{background:#e6ebef;text-align:left}th,td{padding:8px;border-bottom:1px solid #d8dfe5}tr:nth-child(even){background:#f5f7f8}.scroll{overflow:auto}summary{cursor:pointer;font-weight:bold}details{margin:20px 0}section{border-top:1px solid #ccd5de;margin-top:26px;padding-top:8px}a{color:#185182}button{padding:8px;margin:8px}.muted{font:13px/1.5 Arial,sans-serif;color:#54616c}@media print{body{background:white}main{padding:0}button{display:none}section{break-before:page}thead{display:table-header-group}tr{break-inside:avoid}.scroll{overflow:visible}table{font-size:8pt}h2{font-size:15pt}p{font-size:10pt}}@page{size:A4 landscape;margin:14mm}</style><body><main>',
 '<h1>PCR testing among subjects reporting symptoms</h1>',
 paste0('<p class="muted">Moderna P301 · v1.1.1 · EUA cutoff ',PHMPT_EUA_CUTOFF,' · BLA cutoff ',PHMPT_BLA_CUTOFF,'.</p>'),
 '<p class="note"><strong>Two nested populations:</strong> COVID_FORM_ONLY includes every eligible subject with a positive COVID-form symptom assessment, whether or not that subject also has an AE. COVID_FORM_PLUS_AE includes subjects with a positive assessment or an AE matching the agreed symptom dictionary. A subject appearing in both sources is counted once in the combined denominator.</p>',
 paste0('<p><strong>Temporally linked testing:</strong> a subject is counted as tested if at least one qualifying completed PCR falls between a recorded symptom-evidence date and ',PHMPT_TEST_LINK_DAYS,' days later, inclusive, within that subject’s follow-up and the selected cutoff. COVID-form dates are positive ADSYMP assessment dates, not necessarily illness onset. AE dates are AESTDTC onset dates; ongoing AEs do not create a new anchor every day. Repeated assessments may link to the same PCR; participant and record totals are deduplicated. Evidence dates are not reconstructed clinical episodes.</p>'),
 '<p><strong>Two PCR definitions:</strong> DIAGNOSTIC requires an illness label or local COVID Diagnostic Test coding; ANY_NONBASELINE allows any completed non-baseline PCR in the interval, including scheduled follow-up tests. Baseline PCR is excluded from both. This is a temporal association, not proof that symptoms triggered the swab or that an independently documented COVID visit occurred. The earlier heatmap visit-plus-PCR endpoint is different.</p>',
 paste0('<p><strong>Follow-up:</strong> BY_CUTOFF includes all subjects with dated symptom evidence within the usable observation window and reports testing observed by that cutoff. Recent symptom evidence can have less than ',PHMPT_TEST_LINK_DAYS,' days of observation. COMPLETE_7D requires at least one evidence date with a full ',PHMPT_TEST_LINK_DAYS,'-day window; incomplete dates are excluded regardless of whether a test occurred. Both the total reporting population and the eligible denominator are exported. A subject with several dates is counted once per analysis; “tested” means at least one qualifying date was followed by PCR, not that every illness was tested.</p>'),
 '<p><strong>Symptom scope:</strong> the prior exact, field-specific AE matching dictionary is reused. Direct included matches count regardless of investigator causal attribution or PCR result. Days 0–6 after either dose are retained. This includes symptom-compatible AEs, not every AE or solicited reactogenicity diary. Fever/duration thresholds and the complete protocol trigger are not established by these records alone. Missing/partial dates and contradictory records are inventoried.</p>',
 '<p><strong>Local and central:</strong> source-route definitions and uncertainty handling are unchanged. Local-source includes MBGRPID 003 or LOCALFL=Y; central-source requires configured group 001 and the expected laboratory name without contradictory local evidence. These are reporting routes, not confirmed physical laboratory locations. Both-source testing can involve different dates or illnesses in the same subject.</p>',
 '<p><strong>Inference:</strong> chi-square or two-sided Fisher exact compares one tested/not-tested outcome per subject within each population. BH and Holm adjustments cover site×endpoint comparisons within each view. The two populations overlap and are not independent cohorts. Conditional arm differences do not by themselves establish investigator bias; symptom patterns, reporting opportunities and follow-up can differ.</p>',
 '<p><a href="symptom_summary_global.csv">Global symptom-conditioned rates</a> · <a href="symptom_summary_by_site.csv">Site rates</a> · <a href="symptom_comparisons_by_site.csv">Site tests</a> · <a href="symptom_population_flow.csv">Denominators and incomplete follow-up</a> · <a href="pcr_cohort_testing_report.html">Original whole-cohort analysis</a></p>',
 '<button onclick="window.print()">Print displayed tables</button><button onclick="document.querySelectorAll(\'details\').forEach(x=>x.open=true)">Expand all tables</button>',
 '<p class="muted">The default tables below use diagnostic-labelled/local-coded PCR. All non-baseline PCR sensitivity results are expandable. Expand site tables before printing to include them.</p>')
 g<-t[t$granularity=='GLOBAL',];views<-unique(g[c('release','period','scope','stratum')]);parts<-header
 for(i in seq_len(nrow(views))) {
  v<-views[i,];take<-t$release==v$release&t$period==v$period&t$scope==v$scope&t$stratum==v$stratum
  a<-t[take&t$granularity=='GLOBAL',];b<-t[take&t$granularity=='SITE',];b<-b[order(b$endpoint,b$p_value,b$site),]
  title<-paste(v$release,v$period,v$scope,gsub('SYMPTOMS__','',v$stratum,fixed=TRUE),sep=' · ')
  sens<-grepl('ANY_NONBASELINE$',v$stratum)
  if(sens)parts<-c(parts,paste0('<details><summary>Sensitivity: ',pt_escape(title),'</summary>')) else parts<-c(parts,paste0('<section><h2>',pt_escape(title),'</h2>'))
  parts<-c(parts,pt_table(pt_display(a)),'<details><summary>Site-level results</summary>',pt_table(pt_display(b,TRUE)),'</details>',if(sens)'</details>' else '</section>')
 }
 parts<-c(parts,'<details><summary>Population flow for both symptom definitions</summary>',pt_table(f),'</details>',
 '<p class="muted">Audit trail: symptom_subjects/ contains deduplicated subjects and testing flags. symptom_links/ contains source-record → evidence-date → PCR links, including repeated-PCR linkage counts. Each release’s symptoms/ folder contains the mapping, original evidence and quality review. These files permit review without treating repeated symptoms or repeated assays as independent subjects.</p></main></body></html>')
 writeLines(enc2utf8(parts),file.path(out,'pcr_testing_report.html'),useBytes=TRUE)
}

analyse_pcr_testing_by_lab <- function(root=PHMPT_ROOT,output_dir=PHMPT_OUTPUT_DIR,
                                      files=PHMPT_FILES,show_report=PHMPT_SHOW_REPORT) {
 root<-ws_root(root)
 if(any(!PHMPT_RELEASES%in%c('EUA','Prod'))||anyDuplicated(PHMPT_RELEASES))stop('Invalid PHMPT_RELEASES')
 if(any(!PHMPT_SCOPES%in%c('STRICT_BLINDED','CALENDAR'))||anyDuplicated(PHMPT_SCOPES))stop('Invalid PHMPT_SCOPES')
 cuts<-c(EUA=ws_dates(PHMPT_EUA_CUTOFF,'EUA cutoff',strict=TRUE),BLA=ws_dates(PHMPT_BLA_CUTOFF,'BLA cutoff',strict=TRUE))
 if(any(!is.finite(cuts))||cuts[1]>cuts[2])stop('Invalid cutoffs')
 out<-if(nzchar(output_dir))ws_path(root,output_dir) else file.path(root,'analysis',paste0('pcr_testing_by_lab_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_',Sys.getpid()))
 if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop('Output directory must be empty: ',out)
 dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out,winslash='/',mustWork=TRUE)
 dir.create(file.path(out,'subjects'));writeLines('RUNNING',file.path(out,'RUN_STATUS.txt'));complete<-FALSE
 on.exit(if(!complete)writeLines('FAILED - incomplete run. See R console and discovery/audit files.',file.path(out,'RUN_STATUS.txt')),add=TRUE)
 sums<-flows<-manifest<-status<-sym_sums<-sym_flows<-list()
 if(isTRUE(PHMPT_ANALYSE_SYMPTOMS)){
  if(length(PHMPT_TEST_LINK_DAYS)!=1L||!is.finite(PHMPT_TEST_LINK_DAYS)||PHMPT_TEST_LINK_DAYS<0||PHMPT_TEST_LINK_DAYS!=floor(PHMPT_TEST_LINK_DAYS))stop('Invalid PCR linking horizon')
  if(!identical(PHMPT_SYMPTOM_COHORTS,c('COVID_FORM_ONLY','COVID_FORM_PLUS_AE')))stop('Both specified symptom populations must be retained')
  if(any(!PHMPT_SYMPTOM_FOLLOWUP%in%c('BY_CUTOFF','COMPLETE_7D'))||any(!PHMPT_SYMPTOM_TESTS%in%c('DIAGNOSTIC','ANY_NONBASELINE')))stop('Invalid symptom settings')
  dir.create(file.path(out,'symptom_subjects'));dir.create(file.path(out,'symptom_links'))
 }
 for(release in PHMPT_RELEASES) {
  cfg<-files[[release]];if(is.null(cfg))cfg<-list();for(n in c('ADSL','MB','ADMB','SUPPMB','AE','ADSYMP','EVENTS'))if(is.null(cfg[[n]]))cfg[[n]]<-''
  available<-nzchar(cfg$ADSL)||length(pt_candidates(root,release,'ADSL'))>0
  if(!available) {status[[release]]<-data.frame(release=release,status='SKIPPED',detail='No ADSL found; no synthetic release substituted');next}
  message('Processing ',release)
  dest<-file.path(out,release);dir.create(dest)
  ad<-pt_pick(root,release,'ADSL',cfg$ADSL,dest);submission<-pt_submission(ad$manifest$path)
  mb<-pt_pick(root,release,'MB',cfg$MB,dest,submission)
  am<-pt_pick(root,release,'ADMB',cfg$ADMB,dest,submission,optional=TRUE)
  sm<-pt_pick(root,release,'SUPPMB',cfg$SUPPMB,dest,submission,optional=TRUE)
  ctx<-pt_context(ad$data,dest);p<-pt_records(mb$data,if(is.null(am))NULL else am$data,if(is.null(sm))NULL else sm$data,ctx,dest)
  for(z in list(ad,mb,am,sm))if(!is.null(z))manifest[[length(manifest)+1L]]<-z$manifest
  status[[release]]<-data.frame(release=release,status='PROCESSED',detail=paste0(nrow(ctx),' randomized; ',nrow(p),' SARSCOV2 MB records; ADMB ',if(is.null(am))'absent' else 'present','; SUPPMB ',if(is.null(sm))'absent' else 'present'))
  if(isTRUE(PHMPT_ANALYSE_SYMPTOMS)){
   ev<-pt_symptom_inputs(root,release,cfg,submission,ctx,dest);manifest[[length(manifest)+1L]]<-ev$manifest
  }
  periods<-if(release=='Prod')c('EUA','BLA') else 'EUA'
  for(period in periods)for(scope in PHMPT_SCOPES)for(stratum in c('FOLLOW_UP_EXCL_BASELINE','ILLNESS_LABELLED','ALL_POST_RANDOMIZATION')) {
   label<-if(release=='Prod'&&period=='EUA')'EUA_DATE_RETROSPECTIVE' else paste0(period,'_SNAPSHOT')
   message('  ',release,' / ',label,' / ',scope,' / ',stratum)
   z<-pt_window(ctx,p,cuts[period],scope,stratum,release,label)
   flows[[length(flows)+1L]]<-z$flow
   if(!nrow(z$subjects))stop('No eligible subjects in ',release,' / ',label,' / ',scope,' / ',stratum)
   pt_dates_write(z$subjects,file.path(out,'subjects',paste(release,label,scope,paste0(stratum,'.csv'),sep='_')))
   sums[[length(sums)+1L]]<-pt_summary(z$subjects)
  }
  if(isTRUE(PHMPT_ANALYSE_SYMPTOMS))for(period in periods)for(scope in PHMPT_SCOPES) {
   label<-if(release=='Prod'&&period=='EUA')'EUA_DATE_RETROSPECTIVE' else paste0(period,'_SNAPSHOT')
   message('  Symptom-conditioned: ',release,' / ',label,' / ',scope)
   w<-pt_symptom_window(ctx,p,ev$events,cuts[period],scope,release,label,out)
   for(pop in PHMPT_SYMPTOM_COHORTS)for(opp in PHMPT_SYMPTOM_FOLLOWUP)for(sel in PHMPT_SYMPTOM_TESTS) {
    r<-pt_conditional_subjects(w,pop,opp,sel,release,label,scope,out)
    sym_sums[[length(sym_sums)+1L]]<-r$summary;sym_flows[[length(sym_flows)+1L]]<-r$flow
   }
  }
 }
 status<-do.call(rbind,status);ws_csv(status,file.path(out,'release_status.csv'))
 if(!length(sums))stop('No complete release was available. Supply ADSL + MB and, when available, ADMB + SUPPMB.')
 s<-do.call(rbind,sums);f<-do.call(rbind,flows);m<-do.call(rbind,manifest);t<-pt_compare(s)
 stopifnot(all(s$n_tested<=s$n_eligible),all(s$n_tested>=0),!anyDuplicated(s[c('release','period','scope','stratum','granularity','site','arm','endpoint')]))
 # Check that site denominators and participant counts sum to the global results.
 k<-c('release','period','scope','stratum','arm','endpoint');a<-aggregate(s[s$granularity=='SITE',c('n_eligible','n_tested')],s[s$granularity=='SITE',k],sum)
 g<-s[s$granularity=='GLOBAL',];key<-function(z)do.call(paste,c(z[k],sep='|'));a<-a[match(key(g),key(a)),]
 stopifnot(identical(as.integer(a$n_eligible),as.integer(g$n_eligible)),identical(as.integer(a$n_tested),as.integer(g$n_tested)))
 ws_csv(s[s$granularity=='GLOBAL',],file.path(out,'summary_global.csv'));ws_csv(s[s$granularity=='SITE',],file.path(out,'summary_by_site.csv'))
 ws_csv(t[t$granularity=='GLOBAL',],file.path(out,'comparisons_global.csv'));ws_csv(t[t$granularity=='SITE',],file.path(out,'comparisons_by_site.csv'))
 ws_csv(t[t$granularity=='SITE'&t$nominal_p_lt_05,],file.path(out,'sites_nominal_p_lt_05.csv'))
 ws_csv(f,file.path(out,'population_flow.csv'));ws_csv(m,file.path(out,'source_manifest.csv'))
 params<-list(version='1.1.1',EUA_cutoff=PHMPT_EUA_CUTOFF,BLA_cutoff=PHMPT_BLA_CUTOFF,
  central_groups=paste(PHMPT_CENTRAL_GROUPS,collapse=';'),local_groups=paste(PHMPT_LOCAL_GROUPS,collapse=';'),central_name_regex=PHMPT_CENTRAL_NAME_REGEX,
  central_rule='group AND expected lab name AND no contradictory local evidence',local_rule='local group OR LOCALFL Y, except conflict',
  symptom_analysis=PHMPT_ANALYSE_SYMPTOMS,PCR_link_days=PHMPT_TEST_LINK_DAYS,symptom_cohorts=paste(PHMPT_SYMPTOM_COHORTS,collapse=';'),
  scopes=paste(PHMPT_SCOPES,collapse=';'),method='randomized-subject testing proportions; repeated record rates descriptive',numeric_date_origin=PHMPT_NUMERIC_DATE_ORIGIN)
 ws_csv(data.frame(parameter=names(params),value=as.character(unlist(params))),file.path(out,'parameters.csv'))
 capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'));pt_report(s,t,f,status,out)
 if(isTRUE(PHMPT_ANALYSE_SYMPTOMS)) {
  ss<-do.call(rbind,sym_sums);sf<-do.call(rbind,sym_flows);st<-pt_compare(ss)
  # Both definitions share the same date/eligibility rules: the union must contain
  # every COVID-form participant and every COVID-form-linked tested participant.
  keys<-c('release','period','scope','granularity','site','arm','endpoint','followup_requirement','pcr_selection')
  key<-function(z)do.call(paste,c(z[keys],sep='|'))
  one<-ss[ss$denominator_source=='COVID_FORM_ONLY',];both<-ss[ss$denominator_source=='COVID_FORM_PLUS_AE',];both<-both[match(key(one),key(both)),]
  stopifnot(!anyNA(both$n_eligible),all(both$n_eligible>=one$n_eligible),all(both$n_tested>=one$n_tested),all(ss$n_tested<=ss$n_eligible))
  kk<-c('release','period','scope','stratum','arm','endpoint');k<-function(z)do.call(paste,c(z[kk],sep='|'))
  sumsite<-aggregate(ss[ss$granularity=='SITE',c('n_eligible','n_tested')],ss[ss$granularity=='SITE',kk],sum)
  glob<-ss[ss$granularity=='GLOBAL',];sumsite<-sumsite[match(k(glob),k(sumsite)),]
  stopifnot(all(sumsite$n_eligible==glob$n_eligible),all(sumsite$n_tested==glob$n_tested))
  ws_csv(ss[ss$granularity=='GLOBAL',],file.path(out,'symptom_summary_global.csv'));ws_csv(ss[ss$granularity=='SITE',],file.path(out,'symptom_summary_by_site.csv'))
  ws_csv(st[st$granularity=='GLOBAL',],file.path(out,'symptom_comparisons_global.csv'));ws_csv(st[st$granularity=='SITE',],file.path(out,'symptom_comparisons_by_site.csv'))
  ws_csv(sf,file.path(out,'symptom_population_flow.csv'))
  if(!file.rename(file.path(out,'pcr_testing_report.html'),file.path(out,'pcr_cohort_testing_report.html')))stop('Cannot retain cohort report')
  pt_symptom_report(ss,st,sf,out)
 }
 writeLines(c('Site sums match global sums.','No participant numerator exceeds its denominator.','No duplicate summary keys.'),file.path(out,'verification_checks.txt'))
 writeLines('COMPLETE',file.path(out,'RUN_STATUS.txt'));complete<-TRUE
 message('Complete: ',file.path(out,'pcr_testing_report.html'))
 if(isTRUE(show_report))try(utils::browseURL(file.path(out,'pcr_testing_report.html')),silent=TRUE)
 invisible(list(output_dir=out,summary=s,comparisons=t,flow=f,symptom_summary=if(isTRUE(PHMPT_ANALYSE_SYMPTOMS))ss else NULL))
}
# ---- RUN: executes when the entire standalone script is run ----
if(isTRUE(getOption('phmpt.autorun',TRUE)))analyse_pcr_testing_by_lab()
