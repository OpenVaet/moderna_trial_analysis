# Moderna P301 | AE audit: EUA (27073_S5), BLA (125752_S3), amendment (125752_S27).
# Standalone RStudio: set working directory to your PHMPT project, then select all and Run.
# Outputs: analysis/ae_s27_<timestamp>_<id>/report.html plus review CSVs.
# Optional: options(phmpt.root="D:/Projects/moderna_trial_analysis")
# For a named input override, e.g. options(phmpt.s27_files=c(S27_AE="prod_data_moderna/xpt_data/..._S_ae.xpt"))
# Use submission-identified original basenames. S3 and S27 are compared, never pooled.
# XPT: haven (preferred) or foreign. SAS7BDAT requires haven; CSV mirrors supported.
# Load without running: options(phmpt.autorun=FALSE). No network calls or package installation.
PHMPT_EUA_CUTOFF <- "2020-11-25"; PHMPT_BLA_CUTOFF <- "2021-03-26"
PHMPT_CSV_DATE_ORIGIN <- ""; PHMPT_CANDIDATE_DATE_TOLERANCE <- 7L
PHMPT_MAX_CANDIDATES_PER_RECORD <- Inf
sc_spec<-list(AE=list(seq="AESEQ",date="AESTDTC",entity="AETERM",refs=c("AESPID","AEREFID"),clinical=c("AESTDTC","AETERM"),result=character(),serious="AESER"))
sc_text <- function(x) {
 if(inherits(x,'Date'))return(ifelse(is.na(x),'',format(x,'%Y-%m-%d')))
 if(inherits(x,'POSIXt'))return(ifelse(is.na(x),'',format(x,'%Y-%m-%dT%H:%M:%S',tz='UTC')))
 if(is.numeric(x)){
  z<-ifelse(is.na(x),'',sprintf('%.17g',as.numeric(x)))
  if(requireNamespace('haven',quietly=TRUE)){i<-haven::is_tagged_na(x);if(any(i))z[i]<-paste0('<SAS_MISSING_',haven::na_tag(x[i]),'>')}
 }else{z<-as.character(x);z[is.na(z)]<-''}
 bad<-!validUTF8(z);if(any(bad))z[bad]<-iconv(z[bad],from='windows-1252',to='UTF-8',sub=NA)
 if(anyNA(z))stop('Undecodable source text; review the source encoding')
 enc2utf8(trimws(gsub('\r\n','\n',z,fixed=TRUE)))
}
sc_norm <- function(x)toupper(trimws(gsub('[[:space:]]+',' ',sc_text(x))))
sc_col <- function(x,n)if(n%in%names(x))x[[n]] else rep('',nrow(x))
sc_first <- function(x,fields){z<-rep('',nrow(x));for(n in fields){v<-sc_col(x,n);i<-!nzchar(z);z[i]<-v[i]};z}
sc_csv <- function(x,path){con<-file(path,'wb');on.exit(close(con));writeBin(as.raw(c(239,187,191)),con);write.table(x,con,sep=',',row.names=FALSE,quote=TRUE,na='',qmethod='double')}
sc_path <- function(root,path)if(grepl('^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)',path))path else file.path(root,path)
sc_sid <- function(x){z<-sc_norm(x);ok<-grepl('US[0-9]{3}-?[0-9]{4}$',z);z[ok]<-gsub('-','',sub('^.*(US[0-9]{3}-?[0-9]{4})$','\\1',z[ok]));z}
sc_day <- function(x){z<-substr(x,1,10);ok<-grepl('^[0-9]{4}-[0-9]{2}-[0-9]{2}$',z);d<-rep(as.Date(NA),length(z));d[ok]<-suppressWarnings(as.Date(z[ok],format='%Y-%m-%d'));bad<-!is.na(d)&format(d,'%Y-%m-%d')!=z;d[bad]<-as.Date(NA);d}
sc_window <- function(d){ifelse(is.na(d),'UNKNOWN_OR_PARTIAL_DATE',ifelse(d<=as.Date(PHMPT_EUA_CUTOFF),'ON_OR_BEFORE_EUA',ifelse(d<=as.Date(PHMPT_BLA_CUTOFF),'AFTER_EUA_THROUGH_BLA','AFTER_BLA_CUTOFF')))}
sc_key <- function(x,fields,normalize=FALSE,require_all=TRUE){
 if(!length(fields))return(rep('',nrow(x)))
 vs<-lapply(fields,function(n){v<-sc_col(x,n);if(normalize)v<-sc_norm(v);v})
 complete<-if(require_all)Reduce('&',lapply(vs,nzchar)) else Reduce('|',lapply(vs,nzchar))
 z<-do.call(paste,c(lapply(vs,function(v)paste0(nchar(v,type='bytes'),':',v)),sep='|'));z[!complete]<-'';z
}
sc_class <- function(v){
 x<-sc_norm(v);z<-rep('OTHER',length(x));z[!nzchar(x)]<-'MISSING'
 z[x%in%c('POSITIVE','P','DETECTED','RNA DETECTED','SARS-COV-2 DETECTED','PRESUMPTIVE POSITIVE','REACTIVE','YES','Y')]<-'POSITIVE_OR_YES'
 z[x%in%c('NEGATIVE','N','NOT DETECTED','RNA NOT DETECTED','SARS-COV-2 NOT DETECTED','NONREACTIVE','NON-REACTIVE','NO')]<-'NEGATIVE_OR_NO'
 z[x%in%c('NOT DONE','NOT PERFORMED','CANCELLED','CANCELED')]<-'NOT_DONE'
 z[x%in%c('INDETERMINATE','INCONCLUSIVE','EQUIVOCAL','INVALID','SEECOMMENT','SEE COMMENT')]<-'NONBINARY';z
}
sc_read <- function(path,domain){
 ext<-tolower(tools::file_ext(path))
 if(ext%in%c('xpt','sas7bdat')){
  if(requireNamespace('haven',quietly=TRUE))x<-if(ext=='xpt')haven::read_xpt(path) else haven::read_sas(path)
  else if(ext=='xpt'&&requireNamespace('foreign',quietly=TRUE))x<-foreign::read.xport(path)
  else stop('Install haven for SAS files (foreign also supports XPT), or provide a CSV mirror.')
 }else if(ext=='csv'){
  x<-tryCatch(read.csv(path,colClasses='character',check.names=FALSE,fileEncoding='UTF-8-BOM',na.strings=character()),error=function(e)read.csv(path,colClasses='character',check.names=FALSE,fileEncoding='windows-1252',na.strings=character()))
 }else if(ext=='jmp'&&domain=='ADSL'){
  x<-ws_read_jmp(path,fields=c('STUDYID','USUBJID','RANDFL','TRT01P','CUTOFFDT'),optional=c('SITEID'))
 }else stop('Unsupported comparison file: ',path)
 x<-as.data.frame(x,stringsAsFactors=FALSE);names(x)<-toupper(names(x));if(anyDuplicated(names(x)))stop('Duplicate column names: ',path)
 if(domain=='ADSL'){needed<-c('STUDYID','USUBJID','RANDFL','TRT01P','CUTOFFDT');if(!all(needed%in%names(x)))stop('Incomplete ADSL arm/cutoff fields');if(!'SITEID'%in%names(x))x$SITEID<-rep('',nrow(x));x<-x[c(needed,'SITEID')]}
 metadata<-data.frame(variable=names(x),type=vapply(x,function(z)paste(class(z),collapse='/'),''),label=vapply(x,function(z){v<-attr(z,'label');if(is.null(v))'' else as.character(v)},''),stringsAsFactors=FALSE)
 datefields<-intersect(names(x),c('ADT','ASTDT','AENDT','RANDDT','CUTOFFDT','TR01SDT','AP01EDT','DOSE2DT','UBLDDT'))
 for(n in datefields)if(is.numeric(x[[n]])&&!inherits(x[[n]],c('Date','POSIXt')))x[[n]]<-as.Date(as.numeric(x[[n]]),origin='1960-01-01')
 x[]<-lapply(x,sc_text);original_rows<-nrow(x)
 if(!all(c('USUBJID','STUDYID')%in%names(x)))stop('Missing STUDYID / USUBJID: ',path)
 x$.row<-seq_len(nrow(x));x<-x[sc_norm(x$STUDYID)=='MRNA-1273-P301',,drop=FALSE]
 if(!nrow(x))stop('No P301 records: ',path)
 if(ext=='csv')for(n in datefields){i<-grepl('^-?[0-9]+([.]0+)?$',x[[n]]);if(any(i)){if(!nzchar(PHMPT_CSV_DATE_ORIGIN))stop('Numeric CSV dates in ',n,' require PHMPT_CSV_DATE_ORIGIN');x[[n]][i]<-format(as.Date(as.numeric(x[[n]][i]),origin=PHMPT_CSV_DATE_ORIGIN),'%Y-%m-%d')}}
 x$.sid<-sc_sid(x$USUBJID);if(any(!grepl('^US[0-9]{7}$',x$.sid)))stop('Unexpected P301 subject IDs in ',path)
 if(domain!='ADSL'){
  sp<-sc_spec[[domain]];if(!any(sp$date%in%names(x))||!any(sp$entity%in%names(x)))stop('Required date/entity columns absent for ',domain)
  x$.date_raw<-sc_first(x,sp$date);x$.day<-as.character(sc_day(x$.date_raw));x$.day[is.na(x$.day)]<-''
  x$.window<-sc_window(sc_day(x$.date_raw));x$.entity<-sc_norm(sc_first(x,sp$entity));x$.seq<-sc_first(x,sp$seq)
 }
 list(data=x,metadata=metadata,original_rows=original_rows,path=normalizePath(path,winslash='/'),md5=unname(tools::md5sum(path)))
}
sc_equivalent <- function(a,b){
 na<-sort(setdiff(names(a),grep('^\\.',names(a),value=TRUE)));nb<-sort(setdiff(names(b),grep('^\\.',names(b),value=TRUE)))
 identical(na,nb)&&nrow(a)==nrow(b)&&identical(sort(sc_key(a,na,require_all=FALSE)),sort(sc_key(b,nb,require_all=FALSE)))
}
sc_pairs <- function(a,b,domain){
 sp<-sc_spec[[domain]];useda<-rep(FALSE,nrow(a));usedb<-rep(FALSE,nrow(b));pairs<-list()
 add<-function(ia,ib,method,confidence){if(!length(ia))return();stopifnot(!any(useda[ia]),!any(usedb[ib]),!anyDuplicated(ia),!anyDuplicated(ib));useda[ia]<<-TRUE;usedb[ib]<<-TRUE;pairs[[length(pairs)+1L]]<<-data.frame(eua_index=ia,bla_index=ib,match_method=method,match_confidence=confidence,stringsAsFactors=FALSE)}
 unique_pair<-function(ka,kb,method,confidence){
  # Uniqueness is checked on the FULL snapshots, not merely on leftover rows.
  ia<-which(!useda&nzchar(ka)&!duplicated(ka)&!duplicated(ka,fromLast=TRUE));ib<-which(!usedb&nzchar(kb)&!duplicated(kb)&!duplicated(kb,fromLast=TRUE))
  j<-match(ka[ia],kb[ib]);ok<-!is.na(j);add(ia[ok],ib[j[ok]],method,confidence)
 }
 for(ref in intersect(sp$refs,intersect(names(a),names(b)))){
  fields<-c('.sid',ref);if(domain%in%c('MB','LB','ADMB','ADSYMP'))fields<-c(fields,'.entity')
  unique_pair(sc_key(a,fields,TRUE),sc_key(b,fields,TRUE),paste0('REFERENCE:',ref),'REFERENCE_LINKED')
 }
 # Exact clinical-date keys deliberately exclude outcomes and severity, so a
 # result reversal can be detected. They remain provisional record identities.
 fields<-c('.sid','.day','.entity')
 for(n in intersect(setdiff(sp$clinical,c(sp$date,sp$entity)),intersect(names(a),names(b)))){
  if(any(nzchar(a[[n]]))&&any(nzchar(b[[n]])))fields<-c(fields,n)
 }
 clinical_key<-function(x){k<-sc_key(x,fields,TRUE,require_all=FALSE);k[!nzchar(sc_key(x,c('.sid','.day','.entity'),TRUE))]<-'';k}
 unique_pair(clinical_key(a),clinical_key(b),'UNIQUE_CLINICAL_DATE_KEY','PROVISIONAL_CLINICAL_LINK')
 # Equal shared content retains multiplicity, including duplicate identical rows.
 # It can establish record-content retention without proving physical event identity.
 volatile<-unique(c('STUDYID','USUBJID','SUBJID','SITEID',sp$seq,'ASEQ','AESEQ','MBSEQ','LBSEQ','SRCDOM','SRCSEQ','SRCVAR','CUTOFFDT'))
 shared<-setdiff(intersect(names(a),names(b)),c(volatile,grep('^\\.',names(a),value=TRUE)))
 if(length(shared)){
  ka<-sc_key(a,c('.sid',shared),require_all=FALSE);kb<-sc_key(b,c('.sid',shared),require_all=FALSE)
  ia<-which(!useda);ib<-which(!usedb)
  # Occurrence rank makes identical-content multiset matching deterministic.
  rank_key<-function(k)paste0(k,'#OCCURRENCE#',ave(seq_along(k),k,FUN=seq_along))
  if(length(ia)&&length(ib)){j<-match(rank_key(ka[ia]),rank_key(kb[ib]));ok<-!is.na(j)
   add(ia[ok],ib[j[ok]],'EXACT_SHARED_CONTENT_MULTISET','CONTENT_RETENTION_ONLY')}
 }
 ans<-if(length(pairs))do.call(rbind,pairs) else data.frame(eua_index=integer(),bla_index=integer(),match_method=character(),match_confidence=character())
 stopifnot(!anyDuplicated(ans$eua_index),!anyDuplicated(ans$bla_index))
 list(pairs=ans,unmatched_eua=which(!useda),unmatched_bla=which(!usedb))
}
sc_review_candidates <- function(a,b,match,domain){
 ia<-match$unmatched_eua;sp<-sc_spec[[domain]];rows<-list();counts<-integer(length(ia));trunc<-logical(length(ia))
 if(!length(ia))return(list(rows=data.frame(),counts=counts,truncated=trunc))
 # Sequence-only links are review candidates, never automatic identity matches.
 indexes<-list();keys_a<-list()
 addkey<-function(name,fields){ka<-sc_key(a,fields,TRUE);kb<-sc_key(b,fields,TRUE);indexes[[name]]<<-split(which(nzchar(kb)),kb[nzchar(kb)]);keys_a[[name]]<<-ka}
 if(any(nzchar(a$.seq))&&any(nzchar(b$.seq)))addkey('SAME_SEQUENCE_REVIEW',c('.sid','.seq'))
 addkey('SAME_DAY_ENTITY_REVIEW',c('.sid','.day','.entity'))
 addkey('SAME_DAY_OTHER_ENTITY_REVIEW',c('.sid','.day'))
 addkey('NEAR_DATE_SAME_ENTITY_REVIEW',c('.sid','.entity'))
 for(ref in intersect(sp$refs,intersect(names(a),names(b))))addkey(paste0('SHARED_REFERENCE_AMBIGUOUS:',ref),c('.sid',ref))
 matched_b<-rep(FALSE,nrow(b));matched_b[match$pairs$bla_index]<-TRUE
 result_a<-sc_first(a,sp$result);result_b<-sc_first(b,sp$result);class_a<-sc_class(result_a);class_b<-sc_class(result_b)
 day_a<-as.numeric(sc_day(a$.date_raw));day_b<-as.numeric(sc_day(b$.date_raw))
 for(k in seq_along(ia)){
  i<-ia[k];hits<-lapply(names(indexes),function(n){key<-keys_a[[n]][i];if(!nzchar(key))return(integer());v<-indexes[[n]][[key]];if(n=='NEAR_DATE_SAME_ENTITY_REVIEW'){keep<-is.finite(day_a[i])&is.finite(day_b[v])&abs(day_b[v]-day_a[i])<=PHMPT_CANDIDATE_DATE_TOLERANCE;v<-v[keep]};v});names(hits)<-names(indexes)
  j<-sort(unique(unlist(hits,use.names=FALSE)));counts[k]<-length(j);trunc[k]<-length(j)>PHMPT_MAX_CANDIDATES_PER_RECORD
  for(v in j)rows[[length(rows)+1L]]<-data.frame(eua_index=i,bla_index=v,subject=a$.sid[i],eua_row=a$.row[i],bla_row=b$.row[v],reason=paste(names(hits)[vapply(hits,function(z)v%in%z,logical(1))],collapse=';'),bla_already_matched=matched_b[v],
    eua_date=a$.date_raw[i],bla_date=b$.date_raw[v],eua_entity=a$.entity[i],bla_entity=b$.entity[v],eua_result=result_a[i],bla_result=result_b[v],eua_result_class=class_a[i],bla_result_class=class_b[v],stringsAsFactors=FALSE)
 }
 list(rows=if(length(rows))do.call(rbind,rows) else data.frame(),counts=counts,truncated=trunc)
}
ws_text <- sc_text
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







reference_s3 <- read.csv(text='subject,seq,cohort,arm,candidates,s3_row,s3_subject_rows\nUS3132037,1,UNLINKED_SERIOUS,mRNA-1273,2,,11\nUS3132037,2,UNLINKED_SERIOUS,mRNA-1273,2,,11\nUS3172069,1,UNLINKED_SERIOUS,Placebo,0,,0\nUS3172205,1,UNLINKED_SERIOUS,Placebo,0,,0\nUS3222430,2,UNLINKED_SERIOUS,Placebo,1,,3\nUS3252166,3,UNLINKED_SERIOUS,mRNA-1273,0,,2\nUS3272035,1,UNLINKED_SERIOUS,Placebo,1,,1\nUS3352424,2,UNLINKED_SERIOUS,Placebo,2,,4\nUS3432094,3,UNLINKED_SERIOUS,mRNA-1273,2,,8\nUS3432260,3,UNLINKED_SERIOUS,Placebo,1,,6\nUS3432392,1,UNLINKED_SERIOUS,mRNA-1273,0,,0\nUS3552191,3,UNLINKED_SERIOUS,mRNA-1273,2,,9\nUS3552191,8,UNLINKED_SERIOUS,mRNA-1273,1,,9\nUS3572167,1,UNLINKED_SERIOUS,Placebo,2,,4\nUS3752129,4,UNLINKED_SERIOUS,mRNA-1273,1,,3\nUS3752153,5,UNLINKED_SERIOUS,mRNA-1273,1,,5\nUS3762161,4,UNLINKED_SERIOUS,mRNA-1273,6,,13\nUS3862053,3,UNLINKED_SERIOUS,Placebo,1,,5\nUS3872062,15,UNLINKED_SERIOUS,Placebo,6,,21\nUS3912024,1,UNLINKED_SERIOUS,mRNA-1273,1,,19\nUS3912028,2,UNLINKED_SERIOUS,Placebo,1,,3\nUS3942151,1,UNLINKED_SERIOUS,Placebo,2,,3\nUS3952033,5,UNLINKED_SERIOUS,Placebo,2,,6\nUS3042210,1,N_TO_Y,Placebo,,1561,2\nUS3092245,6,Y_TO_N,Placebo,,3033,4\nUS3102052,2,N_TO_Y,Placebo,,3170,9\nUS3182227,1,Y_TO_N,Placebo,,6395,7\nUS3182227,3,Y_TO_N,Placebo,,6396,7\nUS3192486,2,N_TO_Y,mRNA-1273,,7543,3\nUS3192486,3,N_TO_Y,mRNA-1273,,7544,3\nUS3202191,2,Y_TO_N,mRNA-1273,,8718,9\nUS3252166,2,N_TO_Y,mRNA-1273,,10759,2\nUS3252553,1,Y_TO_N,Placebo,,10961,9\nUS3342343,1,N_TO_Y,Placebo,,14833,5\nUS3432094,4,N_TO_Y,mRNA-1273,,17958,8\nUS3532158,2,Y_TO_N,Placebo,,22614,4\nUS3532231,4,Y_TO_N,mRNA-1273,,22719,8\nUS3542199,1,N_TO_Y,Placebo,,23094,1\nUS3552191,6,N_TO_Y,mRNA-1273,,23599,9\nUS3672232,1,Y_TO_N,Placebo,,28568,2\nUS3672243,1,Y_TO_N,Placebo,,28578,2\nUS3712013,1,N_TO_Y,Placebo,,29904,1\nUS3712114,1,N_TO_Y,Placebo,,30020,1\nUS3722066,1,Y_TO_N,Placebo,,30225,8\nUS3722406,2,Y_TO_N,mRNA-1273,,30417,5\nUS3742077,1,N_TO_Y,mRNA-1273,,30984,2\nUS3862053,2,N_TO_Y,Placebo,,35149,5\nUS3872062,2,Y_TO_N,Placebo,,35670,21\nUS3872324,1,Y_TO_N,Placebo,,36349,3\nUS3872324,2,Y_TO_N,Placebo,,36348,3\nUS3872324,3,Y_TO_N,Placebo,,36350,3\nUS3932238,1,Y_TO_N,Placebo,,37944,2\nUS3942001,1,N_TO_Y,Placebo,,38014,1\nUS3942151,3,N_TO_Y,Placebo,,38218,3',colClasses="character",na.strings="",stringsAsFactors=FALSE)
ids<-c(EUA="27073_S5",S3="125752_S3",S27="125752_S27")
colv<-function(z,k)sc_col(z,k)
flag27<-function(x){z<-sc_norm(x);ifelse(z%in%c("Y","YES"),"Y",ifelse(z%in%c("N","NO"),"N",ifelse(!nzchar(z),"MISSING","OTHER")))}
read_aux<-function(p,dom){
  ext<-tolower(tools::file_ext(p))
  fields<-if(dom=="DM")c("STUDYID","USUBJID")else c("STUDYID","USUBJID","RDOMAIN","IDVAR","IDVARVAL","QNAM","QVAL")
  x<-if(ext=="jmp")ws_read_jmp(p,fields,if(dom=="SUPPAE")c("QLABEL","QORIG","QEVAL")else character()) else if(ext=="csv")read.csv(p,colClasses="character",fileEncoding="UTF-8-BOM",check.names=FALSE,na.strings=character()) else if(requireNamespace("haven",quietly=TRUE)){if(ext=="xpt")haven::read_xpt(p) else haven::read_sas(p)}else if(ext=="xpt"&&requireNamespace("foreign",quietly=TRUE))foreign::read.xport(p) else stop("Install haven to read ",p)
  x<-as.data.frame(x);names(x)<-toupper(names(x));stopifnot(all(fields%in%names(x)),!anyDuplicated(names(x)))
  if(dom=="SUPPAE")fields<-union(fields,intersect(c("QLABEL","QORIG","QEVAL"),names(x)))
  x<-x[,fields,drop=FALSE];x[]<-lapply(x,sc_text);x<-x[sc_norm(x$STUDYID)=="MRNA-1273-P301",,drop=FALSE]
  if(!nrow(x))stop("No P301 records: ",p)
  x$.sid<-sc_sid(x$USUBJID);stopifnot(all(grepl("^US[0-9]{7}$",x$.sid)))
  if(dom=="DM"&&anyDuplicated(x$.sid))stop("Duplicate subject in DM.")
  x
}
select_tables<-function(root,over,out) {
  known<-as.vector(outer(names(ids),c("AE","ADSL","DM","SUPPAE"),paste,sep="_"))
  if(length(over)&& (is.null(names(over))||anyDuplicated(names(over))||!all(names(over)%in%known)))stop("Use named phmpt.s27_files, e.g. S27_AE.")
  base<-function(r){b<-if(r=="EUA")"eua_data" else "prod_data";q<-file.path(root,c(paste0(b,"_moderna"),b));q[which(dir.exists(q))[1]]}
  listing<-lapply(c("EUA","S3"),function(r){b<-base(r);if(is.na(b))character() else list.files(b,recursive=TRUE,full.names=TRUE,pattern="[.](xpt|sas7bdat|csv|jmp)$",ignore.case=TRUE)});names(listing)<-c("EUA","Prod")
  log<-data.frame(snapshot=character(),domain=character(),status=character(),path=character(),md5=character(),rows=integer(),message=character())
  add<-function(r,d,st,p="",md="",n=NA_integer_,msg=""){log[nrow(log)+1L,]<<-list(r,d,st,p,md,n,msg);sc_csv(log,file.path(out,"discovery.csv"))}
  get1<-function(r,d){
    key<-paste(r,d,sep="_");override<-key%in%names(over);paths<-if(override)sc_path(root,unname(over[key]))else listing[[if(r=="EUA")"EUA" else "Prod"]]
    paths<-paths[grepl(paste0("(^|[^0-9A-Za-z])",ids[r],"([^0-9A-Za-z]|$)"),basename(paths),ignore.case=TRUE)&grepl(paste0("_",tolower(d),"[.](xpt|sas7bdat|csv|jmp)$"),basename(paths),ignore.case=TRUE)&!grepl("p101|p201|dmid",basename(paths),ignore.case=TRUE)]
    if(d=="AE")paths<-paths[!grepl("[.]jmp$",paths,ignore.case=TRUE)]
    paths<-paths[!grepl("mrna-1273-p[0-9]+",basename(paths),ignore.case=TRUE)|grepl("mrna-1273-p301([^0-9]|$)",basename(paths),ignore.case=TRUE)]
    paths<-unique(paths)
    if(!length(paths)){add(r,d,"MISSING");if(d=="AE"||override)stop("Missing verified ",key," (expected ",ids[r],"). See discovery.csv.");return(NULL)}
    orderfmt<-if(d=="ADSL")c("xpt","sas7bdat","csv","jmp") else c("xpt","sas7bdat","csv","jmp")
    if(!requireNamespace("haven",quietly=TRUE)&&!requireNamespace("foreign",quietly=TRUE)&&any(grepl("[.]csv$",paths,ignore.case=TRUE)))orderfmt<-c("csv","jmp","xpt","sas7bdat")
    ranking<-match(tolower(tools::file_ext(paths)),orderfmt)
    if(!override&&any(ranking!=min(ranking)))for(p in paths[ranking!=min(ranking)])add(r,d,"LOWER_PRIORITY_MIRROR",p)
    paths<-sort(paths[ranking==min(ranking)]);sel<-NULL;hashes<-character()
    for(p in paths){
      if(!file.exists(p))stop("Missing input: ",p)
      md<-unname(tools::md5sum(p));if(md%in%hashes){add(r,d,"IDENTICAL_COPY",p,md);next}
      a<-tryCatch(if(d%in%c("AE","ADSL"))sc_read(p,d)$data else read_aux(p,d),error=function(e){add(r,d,"READ_ERROR",p,md,msg=conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
      if(d=="AE"){need<-c("DOMAIN","AESEQ","AESER","AETERM","AESTDTC","AEDECOD")
        if(!all(need%in%names(a))||any(sc_norm(a$DOMAIN)!="AE")||any(!nzchar(a$AESEQ))||anyDuplicated(a[,c(".sid","AESEQ")]))stop("Invalid subject/AESEQ key or AE schema in ",p)}
      if(d=="ADSL"&&anyDuplicated(a$.sid))stop("Nonunique ADSL roster.")
      if(is.null(sel)){sel<-a;add(r,d,"SELECTED",normalizePath(p,winslash="/"),md,nrow(a))}
      else if(sc_equivalent(sel,a))add(r,d,"EQUIVALENT_VALUES",p,md,nrow(a))
      else {add(r,d,"CONFLICT",p,md,nrow(a));stop("Conflicting ",key," snapshots. Set options(phmpt.s27_files=c(",key,"='...')).")}
      hashes<-c(hashes,md)
    };sel
  }
  data<-list()
  for(r in names(ids))for(domain in c("AE","ADSL","DM","SUPPAE"))data[[paste(r,domain,sep="_")]]<-get1(r,domain)
  list(data=data,manifest=log)
}
# Pairwise S27 comparisons. No cross-submission link uses AESEQ alone.
link27<-function(a,b){
  # Start with identical shared source content except sequence. Repeated copies
  # prove multiset retention only; the physical identity of copies is not assigned.
  cols<-setdiff(intersect(names(a),names(b)),c("AESEQ",grep("^[.]",union(names(a),names(b)),value=TRUE)))
  ka<-sc_key(a,c(".sid",cols),require_all=FALSE);kb<-sc_key(b,c(".sid",cols),require_all=FALSE)
  occ<-function(k)paste(k,ave(seq_along(k),k,FUN=seq_along),sep="#")
  j<-match(occ(ka),occ(kb));ii<-which(!is.na(j));jj<-j[ii]
  dup<-duplicated(ka)|duplicated(ka,fromLast=TRUE)
  dupb<-duplicated(kb)|duplicated(kb,fromLast=TRUE)
  p<-data.frame(i=ii,j=jj,rule=ifelse(dup[ii]|dupb[jj],"IDENTICAL_COPY_GROUP","IDENTICAL_SHARED_CONTENT"),stringsAsFactors=FALSE)
  used<-seq_len(nrow(a))%in%ii;seen<-seq_len(nrow(b))%in%jj
  props<-list()
  for(cols0 in list(c(".sid",".day",".entity"),c(".sid",".day","AEDECOD","AESPID"),c(".sid",".day","AEDECOD","AEREFID"),c(".sid",".entity","AEDECOD","AESPID"),c(".sid",".entity","AEDECOD","AEREFID"))){
    if(!all(cols0%in%cols)&&!all(cols0%in%intersect(names(a),names(b))))next
    x<-sc_key(a,cols0,TRUE);y<-sc_key(b,cols0,TRUE)
    ia<-which(nzchar(x)&!duplicated(x)&!duplicated(x,fromLast=TRUE));ib<-which(nzchar(y)&!duplicated(y)&!duplicated(y,fromLast=TRUE))
    v<-match(x[ia],y[ib]);ok<-!is.na(v);if(!".day"%in%cols0){ k<-which(ok);ld<-a$.day[ia[k]];rd<-b$.day[ib[v[k]]];ok[k]<-!nzchar(ld)|!nzchar(rd)|ld==rd };if(any(ok))props[[length(props)+1L]]<-data.frame(i=ia[ok],j=ib[v[ok]],rule=paste(cols0[-1],collapse="+"))
  }
  if(length(props)){
    g<-do.call(rbind,props);g<-g[!duplicated(g[c("i","j")]),];ok<-!duplicated(g$i)&!duplicated(g$i,fromLast=TRUE)&!duplicated(g$j)&!duplicated(g$j,fromLast=TRUE)&!used[g$i]&!seen[g$j]
    p<-rbind(p,g[ok,])
  }
  stopifnot(!anyDuplicated(p$i),!anyDuplicated(p$j),all(a$.sid[p$i]==b$.sid[p$j]))
  p
}
# Explicit clinical and identifier reasons for every candidate among the same person's rows.
# Do not truncate, and do not exclude rows already allocated to another source event.
candidates<-function(a,i,b,ids_b){
  v<-which(b$.sid==a$.sid[i]);if(!length(v))return(data.frame(j=integer(),evidence=character(),assigned_from=integer()))
  eq<-function(col){x<-sc_norm(colv(a,col)[i]);y<-sc_norm(colv(b,col)[v]);nzchar(x)&x==y}
  delta<-abs(as.numeric(as.Date(b$.day[v],format="%Y-%m-%d")-as.Date(a$.day[i],format="%Y-%m-%d")));same<-!is.na(delta)&delta==0;near<-!is.na(delta)&delta<=7
  # A partial onset may suggest review; it never qualifies as a complete matching date.
  raw<-a$.date_raw[i];tr<-b$.date_raw[v];partial<-grepl("^\\d{4}(-\\d{2})?$",raw)|grepl("^\\d{4}(-\\d{2})?$",tr)
  partial<-partial&nzchar(raw)&nzchar(tr)&(startsWith(tr,raw)|vapply(tr,function(z)startsWith(raw,z),logical(1)))
  ev<-cbind(SAME_AESPID=eq("AESPID"),SAME_AEREFID=eq("AEREFID"),SAME_AESEQ_ONLY_WITHIN_SUBJECT=eq("AESEQ"),SAME_DAY=same,VERBATIM_NEAR_DAY=eq("AETERM")&near,PT_NEAR_DAY=eq("AEDECOD")&near,TERM_PARTIAL_DATE=(eq("AETERM")|eq("AEDECOD"))&partial)
  r<-apply(ev,1,function(v)paste(names(v)[v],collapse=";"));ok<-nzchar(r)
  data.frame(j=v[ok],evidence=r[ok],assigned_from=ids_b[v[ok]])
}
emp<-function(x,cols){if(length(x))do.call(rbind,x) else as.data.frame(setNames(rep(list(character()),length(cols)),cols))}
daytest<-function(s) !is.na(sc_day(s)) & sc_day(s)<=as.Date(PHMPT_EUA_CUTOFF)
link_frame<-function(p,a,b,arms){
  z<-data.frame(subject=a$.sid[p$i],arm=unname(arms[a$.sid[p$i]]),old_row=a$.row[p$i],new_row=b$.row[p$j],old_seq=a$.seq[p$i],new_seq=b$.seq[p$j],rule=p$rule,
    old_onset=a$.date_raw[p$i],new_onset=b$.date_raw[p$j],old_term=a$AETERM[p$i],new_term=b$AETERM[p$j],old_pt=a$AEDECOD[p$i],new_pt=b$AEDECOD[p$j],old_ser=flag27(a$AESER[p$i]),new_ser=flag27(b$AESER[p$j]))
  z$seq_changed<-z$old_seq!=z$new_seq;z$flag_changed<-z$old_ser!=z$new_ser;z
}
supp_index<-function(x,ae,rel,out){
 if(is.null(x))return(NULL)
 good<-sc_norm(x$RDOMAIN)=="AE"&sc_norm(x$IDVAR)=="AESEQ"&nzchar(x$IDVARVAL)&nzchar(x$QNAM)
 ix<-match(paste(x$.sid,x$IDVARVAL),paste(ae$.sid,ae$.seq));x$parent_row<-ae$.row[ix]
 bad<-!good|is.na(ix);if(any(bad))sc_csv(x[bad,],file.path(out,paste0(rel,"_suppae_unlinked.csv")))
 x$parent_index<-ix;x[!bad,intersect(c(".sid","IDVARVAL","parent_row","parent_index","QNAM","QVAL","QLABEL","QORIG","QEVAL"),names(x)),drop=FALSE]
}
compare_supp<-function(s,t,map,out){
 if(is.null(s)||is.null(t))return(NULL)
 # Collapse repeated qualifier values as a multiset; distinguish an empty value from no row.
 prep<-function(x,indices){i<-indices[x$parent_index];k<-!is.na(i);z<-data.frame(link=i[k],QNAM=x$QNAM[k],val=x$QVAL[k])
  z<-z[order(z$link,z$QNAM,z$val),];key<-paste(z$link,z$QNAM,sep="|");duplicates<-duplicated(key)|duplicated(key,fromLast=TRUE)
  if(!any(duplicates)){names(z)[3]<-"value";return(z)}
  vals<-vapply(split(z$val[duplicates],key[duplicates]),function(y)paste0(nchar(y),":",y,collapse=" | "),character(1))
  z$val[duplicates]<-paste0("[repeated] ",vals[key[duplicates]]);z<-z[!duplicated(key),];names(z)[3]<-"value";z}
 x<-prep(s,match(seq_len(max(c(map$i,s$parent_index,0L))),map$i));y<-prep(t,match(seq_len(max(c(map$j,t$parent_index,0L))),map$j))
 v<-merge(x,y,by=c("link","QNAM"),all=TRUE,suffixes=c("_s3","_s27"))
 v$change<-ifelse(is.na(v$value_s3),"QUALIFIER_ADDED",ifelse(is.na(v$value_s27),"QUALIFIER_ABSENT",ifelse(v$value_s3==v$value_s27,"UNCHANGED","VALUE_CHANGED")))
 v<-v[v$change!="UNCHANGED",,drop=FALSE]
 v$s3_index<-map$i[v$link];v$s27_index<-map$j[v$link];v
}
write_html27<-function(out,baseline,ledger,fields,suppsum,pairsum,zerotab,coverage,transitions,remove_flags,naive){
 esc<-function(x){x<-as.character(x);x[is.na(x)]<-"";for(a in list(c("&","&amp;"),c("<","&lt;"),c(">","&gt;"),c('"',"&quot;")))x<-gsub(a[1],a[2],x,fixed=TRUE);x}
 table<-function(x){if(!nrow(x))return("<p>No rows.</p>");paste0('<div class="scroll"><table><thead><tr>',paste0('<th>',esc(names(x)),'</th>',collapse=""),'</tr></thead><tbody>',paste(apply(x,1,function(r)paste0('<tr>',paste0('<td>',esc(r),'</td>',collapse=""),'</tr>')),collapse=""),'</tbody></table></div>')}
 l<-ledger[ledger$cohort=="UNLINKED_SERIOUS",];v<-ledger[ledger$cohort!="UNLINKED_SERIOUS",]
 arm_summary<-function(x,key){if(!nrow(x))return(data.frame());z<-as.data.frame.matrix(table0<-base::table(x[[key]],x$arm));data.frame(outcome=rownames(z),z,row.names=NULL,check.names=FALSE)}
 slim<-c("subject","arm","eua_term","eua_onset","s3_candidates","s3_ae_rows","s27_ae_rows","direct_rule","weak_s27_candidates","assessment")
 swap<-c("subject","arm","eua_term","s3_term","eua_ser","s3_ser","s27_ser_via_s3","s27_ser_direct","trajectory")
 h<-c('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>P301 | S27 adverse-event audit</title><style>body{font:16px/1.48 system-ui;color:#17342b;max-width:1230px;margin:2em auto;padding:0 1.2em}h2{margin-top:1.8em}.note{background:#edf3ee;padding:1em}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:.84em}th,td{border-bottom:1px solid #bdcec6;padding:.35em;text-align:left;vertical-align:top}th{background:#e2ece6}summary{cursor:pointer}code{color:#254b37}@media print{.scroll{overflow:visible}table{font-size:8pt}tr{break-inside:avoid}}</style><h1>Adverse events: EUA, licensing S3 and amendment S27</h1>',
 sprintf('<p class="note">The checked S27 AE table contains %s records in %s participants. The EUA–S3 historical screen identifies %d serious EUA records without an assigned S3 row and %d linked records whose binary seriousness flag changed. Counts describe records unless labelled participants.</p>',coverage$ae_rows[3],coverage$ae_subjects[3],nrow(l),nrow(v)),
 sprintf('<p><b>S3 to S27:</b> %s of %s S3 rows linked to S27; %s changed AESEQ. The linked events showed %d Y→N and %d N→Y seriousness changes. In contrast, pairing on participant and sequence alone would produce %d different seriousness values among %d shared keys (including missing values). Sequence alone is not an event identifier across submissions.</p>',pairsum$linked[1],coverage$ae_rows[2],pairsum$different_AESEQ[1],pairsum$Y_to_N[1],pairsum$N_to_Y[1],naive$discordant_seriousness,naive$shared_sequence_keys),
 sprintf('<p><b>Historical cases:</b> %d/%d originally unlinked serious records acquired a direct clinical/content link to S27. For %d, an earlier S3 candidate has an assigned S27 counterpart; this does not confirm its identity with the EUA event. For %d cases, the participant has no AE rows in S27; for %d, other AE rows are present but no direct candidate was identified. Of the %d historical binary seriousness transitions, %d linked S3 rows retain their S3 flag in S27; %d also have an independent EUA→S27 clinical/content link to the same destination.</p>',sum(!is.na(l$s27_direct_row)),nrow(l),sum(l$assessment=="INDIRECT_FROM_S3_CANDIDATE_ONLY"),sum(l$assessment=="NO_S27_AE_ROWS_FOR_SUBJECT"),sum(l$assessment=="NO_DIRECT_CANDIDATE_OTHER_AE_PRESENT"),nrow(v),sum(v$trajectory=="S3_FLAG_RETAINED_ON_LINKED_ROW"),sum(!is.na(v$s27_direct_row)&!is.na(v$s27_via_s3_row)&v$s27_direct_row==v$s27_via_s3_row)),
 if(!is.null(remove_flags)&&nrow(remove_flags))paste0('<p><b>Supplementary AR Remove Flag:</b> ',nrow(remove_flags),' linked records had REMOVEFL=Y newly added (',esc(paste(names(table0<-base::table(remove_flags$arm)),as.integer(table0),collapse="; ")),'). ',sum(remove_flags$new_ser=="N"),' were non-serious in S27. They remain present in the S27 AE table. These are newly added flags, not the total already flagged. See AR_remove_flag_totals.csv. What effect the derived flag had on an analysis requires the analysis rules and outputs.</p>'),
 '<h2>1. Reproduction of the published audit screen</h2>' ,table(baseline),'<p>The unmatched screen uses known AE onset on or before 25 November 2020. The Y→N and N→Y screens include all linked EUA records, including unknown onset. “Previous” is the reference result, not a value imposed on the new run. See baseline_identity_check.csv for record-level checks.</p>',
 '<h2>2. What S27 contains for the historical unresolved records</h2>',table(arm_summary(l,"assessment")),
 '<p>Direct clinical links require independent agreement beyond sequence alone. Candidates are leads, not confirmed events. The full case ledger shows direct and indirect paths separately. A subject still appearing in demographics is not proof of continuing participation.</p>',table(l[slim]),
 '<h3>Originally without an S3 candidate</h3>',table(zerotab),
 '<h2>3. Previously reclassified seriousness</h2>',table(arm_summary(v,"trajectory")),
 '<p>“S3 flag retained” traces the assigned S3 row to S27. It does not validate the original EUA–S3 match. Source date/term conflicts, competing S27 destinations and missing flags remain explicit.</p>',table(v[swap]),
 '<h2>4. S3 versus S27 across all events</h2>',table(coverage),table(pairsum),table(transitions),
 '<p>First compare identical shared source content excluding AESEQ. Remaining one-to-one links use a unique participant/day/verbatim term or a corroborated reference and clinical descriptor. A partial date never qualifies as a complete date. Competing proposals are not forced. An identical duplicate group preserves content multiplicity only. AESEQ joins qualifiers within a snapshot, not events across snapshots.</p>',table(fields[fields$n_changed>0,]),
 '<h3>Supplementary AE changes on linked parent events</h3>',if(is.null(suppsum))'<p>Not assessed: one or both supplementary tables unavailable.</p>' else table(suppsum),
 if(!is.null(suppsum)&&any(suppsum$qualifier=="AETRTEM"))'<p>For AETRTEM (Treatment Emergent Flag), an absent qualifier means that the S27 SUPPAE table supplies no value on that linked AE; it is not a change to AESER. For individual changes, including the hospital start/end qualifiers, see <code>s3_s27_suppae_changes.csv</code>.</p>',
 '<p>Review files: <code>focus_case_ledger.csv</code>, <code>candidate_paths.csv</code>, <code>focus_subject_events.csv</code> and <code>manual_review_template.csv</code>. Each run retains source paths, checksums and scope in <code>discovery.csv</code>. The comparison establishes what is present in these released AE tables. It does not establish that a record was deleted from the clinical database, was correctly medically reclassified, or was excluded from an FDA analysis.</p></html>')
 writeLines(enc2utf8(h),file.path(out,"report.html"),useBytes=TRUE)
}
analyse_s27_adverse_events<-function(root=getOption("phmpt.root",getwd()),out=getOption("phmpt.s27_output",""),files=getOption("phmpt.s27_files",character())){
 root<-normalizePath(root,mustWork=TRUE,winslash="/")
 if(!nzchar(out)){base<-file.path(root,"analysis");dir.create(base,recursive=TRUE,showWarnings=FALSE);out<-tempfile(paste0("ae_s27_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_"),tmpdir=base)}else out<-sc_path(root,out)
 if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop("Use an empty output folder: ",out)
 dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out,winslash="/")
 cat("S27 review output: ",out,"\n",sep="");status<-file.path(out,"RUN_STATUS.txt");writeLines("RUNNING",status);done<-FALSE;on.exit(if(!done)writeLines("FAILED: incomplete, not a negative result. See discovery.csv and R error.",status),add=TRUE)
 wr<-function(d,f)sc_csv(d,file.path(out,f))
 z<-select_tables(root,files,out);d<-z$data;a<-d$EUA_AE;b<-d$S3_AE;t<-d$S27_AE
 sets<-list(EUA=a,S3=b,S27=t);sub<-sort(unique(unlist(lapply(sets,function(x)x$.sid))))
 roster<-data.frame(subject=sub,stringsAsFactors=FALSE)
 for(n in names(sets)){
  v<-d[[paste0(n,"_ADSL")]];dm<-d[[paste0(n,"_DM")]]
  roster[[paste0(n,"_AE_rows")]]<-tabulate(match(sets[[n]]$.sid,sub),length(sub))
  roster[[paste0(n,"_ADSL_present")]]<-if(is.null(v))NA else sub%in%v$.sid
  roster[[paste0(n,"_DM_present")]]<-if(is.null(dm))NA else sub%in%dm$.sid
  roster[[paste0(n,"_arm")]]<-if(is.null(v))"" else v$TRT01P[match(sub,v$.sid)]
 }
 roster$arm<-apply(roster[c("EUA_arm","S3_arm","S27_arm")],1,function(x){u<-unique(sc_norm(x));u<-u[nzchar(u)];if(!length(u))"UNMAPPED" else if(length(u)>1)"ARM_CONFLICT" else if(u=="MRNA-1273")"mRNA-1273" else if(u=="PLACEBO")"Placebo" else u})
 if(all(roster$arm=="UNMAPPED"))stop("No usable randomized-arm roster. Provide EUA_ADSL or S3_ADSL in phmpt.s27_files.")
 arms<-setNames(roster$arm,roster$subject);wr(roster,"subject_presence.csv")
 coverage<-do.call(rbind,lapply(names(sets),function(n){x<-sets[[n]];data.frame(snapshot=n,submission=unname(ids[n]),ae_rows=nrow(x),ae_subjects=length(unique(x$.sid)),serious=sum(flag27(x$AESER)=="Y"),not_serious=sum(flag27(x$AESER)=="N"),unknown_seriousness=sum(!flag27(x$AESER)%in%c("Y","N")),suppae_available=!is.null(d[[paste0(n,"_SUPPAE")]]))}))
 wr(coverage,"snapshot_coverage.csv")
 message("Replaying original EUA/S3 algorithm; then linking each source independently to S27.")
 m<-sc_pairs(a,b,"AE");q<-sc_review_candidates(a,b,m,"AE");p<-m$pairs
 u<-m$unmatched_eua;ii<-u[daytest(a$.date_raw[u])&flag27(a$AESER[u])=="Y"]
 srcflag<-flag27(a$AESER[p$eua_index]);dstflag<-flag27(b$AESER[p$bla_index])
 moves<-p$eua_index[srcflag%in%c("Y","N")&dstflag%in%c("Y","N")&srcflag!=dstflag]
 old<-c(ii,moves);bidx<-p$bla_index[match(old,p$eua_index)]
 s3key<-paste(a$.sid[old],a$.seq[old]);refkey<-paste(reference_s3$subject,reference_s3$seq)
 co<-c(rep("UNLINKED_SERIOUS",length(ii)),if(length(moves))paste(flag27(a$AESER[moves]),"TO",flag27(b$AESER[bidx[seq_along(moves)+length(ii)]]),sep="_") else character())
 oldq<-q$counts[match(old,u)]
 ver<-data.frame(subject=a$.sid[old],seq=a$.seq[old],cohort=co,arm=unname(arms[a$.sid[old]]),candidates=oldq,s3_subject_rows=roster$S3_AE_rows[match(a$.sid[old],sub)])
 ref<-reference_s3[,setdiff(names(reference_s3),"s3_row")];check<-merge(ver,ref,by=c("subject","seq"),all=TRUE,suffixes=c("_run","_reference"))
 equal<-function(x,y)ifelse(is.na(x)|is.na(y),is.na(x)&is.na(y),as.character(x)==as.character(y))
 check$agrees<-!is.na(check$cohort_run)&!is.na(check$cohort_reference)
 for(f in c("cohort","arm","candidates","s3_subject_rows"))check$agrees<-check$agrees&equal(check[[paste0(f,"_run")]],check[[paste0(f,"_reference")]])
 wr(check,"baseline_identity_check.csv")
 metrics<-function(v){u<-v$cohort=="UNLINKED_SERIOUS";c(sum(u),sum(u&v$candidates>0,na.rm=TRUE),sum(u&v$candidates==0,na.rm=TRUE),length(unique(v$subject[u&v$s3_subject_rows==0])),sum(v$cohort=="Y_TO_N"),sum(v$cohort=="N_TO_Y"))}
 baseline<-do.call(rbind,lapply(c("Both arms","mRNA-1273","Placebo",setdiff(unique(ver$arm),c("mRNA-1273","Placebo"))),function(ar){
  x<-if(ar=="Both arms")ver else ver[ver$arm==ar,];y<-if(ar=="Both arms")ref else ref[ref$arm==ar,]
  data.frame(arm=ar,metric=c("Serious EUA rows not linked to S3","With S3 candidates","Without S3 candidates","Distinct people with zero S3 AE rows","Y to N, all linked EUA dates","N to Y, all linked EUA dates"),previous=metrics(y),current=metrics(x))}))
 wr(baseline,"baseline_counts.csv")
 bt<-link27(b,t);at<-link27(a,t);bc<-link_frame(bt,b,t,arms);ac<-link_frame(at,a,t,arms)
 wr(bc,"s3_s27_links.csv");wr(ac,"eua_s27_links.csv")
 map<-function(p,n,m){x<-rep(NA_integer_,n);y<-rep(NA_integer_,m);x[p$i]<-p$j;y[p$j]<-p$i;list(to=x,from=y)}
 bm<-map(bt,nrow(b),nrow(t));am<-map(at,nrow(a),nrow(t))
 pairsum<-do.call(rbind,lapply(c("S3 / S27","EUA / S27"),function(n){x<-if(n=="S3 / S27")bc else ac;src<-if(n=="S3 / S27")b else a;data.frame(comparison=n,linked=nrow(x),source_unassigned=nrow(src)-nrow(x),s27_unassigned=nrow(t)-nrow(x),different_AESEQ=sum(x$seq_changed),Y_to_N=sum(x$old_ser=="Y"&x$new_ser=="N"),N_to_Y=sum(x$old_ser=="N"&x$new_ser=="Y"))}))
 wr(pairsum,"pair_coverage.csv")
 jk<-match(paste(b$.sid,b$.seq),paste(t$.sid,t$.seq));use<-which(!is.na(jk))
 naive<-data.frame(shared_sequence_keys=length(use),discordant_seriousness=sum(flag27(b$AESER[use])!=flag27(t$AESER[jk[use]])))
 wr(naive,"sequence_only_negative_control.csv")
 trans<-as.data.frame(table(old_ser=bc$old_ser,s27_ser=bc$new_ser),stringsAsFactors=FALSE);wr(trans,"s3_s27_seriousness_matrix.csv")
 shared<-setdiff(intersect(names(b),names(t)),grep("^[.]",names(b),value=TRUE));delta<-list()
 fields<-do.call(rbind,lapply(shared,function(f){oldv<-b[[f]][bt$i];newv<-t[[f]][bt$j];k<-which(oldv!=newv);if(length(k))delta[[length(delta)+1L]]<<-data.frame(subject=bc$subject[k],s3_seq=bc$old_seq[k],s27_seq=bc$new_seq[k],rule=bc$rule[k],field=f,s3_value=oldv[k],s27_value=newv[k]);data.frame(field=f,n_changed=length(k),people=length(unique(bc$subject[k]))) }))
 wr(fields,"s3_s27_field_counts.csv");wr(emp(delta,c("subject","s3_seq","s27_seq","rule","field","s3_value","s27_value")),"s3_s27_field_changes.csv")
 unmatched<-function(x,k,r){v<-x[k,c(".sid",".row",".seq","AETERM","AEDECOD","AESTDTC","AESER"),drop=FALSE];v$snapshot<-rep(r,nrow(v));v}
 wr(rbind(unmatched(b,is.na(bm$to),"S3"),unmatched(t,is.na(bm$from),"S27")),"s3_s27_unassigned.csv")
 # Full participant-level context for the focal rows, including events already paired elsewhere.
 focus<-unique(a$.sid[old]);evcols<-c(".row",".sid","AESEQ","AESPID","AEREFID","AETERM","AEDECOD","AESTDTC","AEENDTC","AESER","AESEV","AESHOSP","AESDTH","AESLIFE","AESDISAB","AESMIE","AEOUT","AESCAT")
 evbind<-function(x,n){v<-as.data.frame(lapply(evcols,function(f)colv(x,f)),stringsAsFactors=FALSE);names(v)<-evcols;v<-v[v$.sid%in%focus,,drop=FALSE];v$snapshot<-rep(n,nrow(v));v}
 wr(do.call(rbind,Map(evbind,sets,names(sets))),"focus_subject_events.csv")
 ss<-lapply(names(sets),function(n)supp_index(d[[paste0(n,"_SUPPAE")]],sets[[n]],n,out));names(ss)<-names(sets)
 for(n in names(ss))if(!is.null(ss[[n]]))wr(ss[[n]][ss[[n]]$.sid%in%focus,],paste0("focus_",n,"_suppae.csv"))
 for(n in names(ss))if(!is.null(ss[[n]]))wr(unique(ss[[n]][intersect(c("QNAM","QLABEL","QORIG"),names(ss[[n]]))]),paste0(n,"_suppae_dictionary.csv"))
 sud<-compare_supp(ss$S3,ss$S27,bt,out);sus<-NULL;remove_flags<-NULL
 if(!is.null(sud)){sud$subject<-b$.sid[sud$s3_index];sud$s3_seq<-b$.seq[sud$s3_index];sud$s27_seq<-t$.seq[sud$s27_index];wr(sud,"s3_s27_suppae_changes.csv");sus<-if(nrow(sud))aggregate(list(parent_qualifiers=rep(1L,nrow(sud))),list(qualifier=sud$QNAM,change=sud$change),sum)else data.frame(qualifier=character(),change=character(),parent_qualifiers=integer());wr(sus,"s3_s27_suppae_summary.csv")}
 if(!is.null(sud)){
   keep<-which(sud$QNAM=="REMOVEFL"&sud$change=="QUALIFIER_ADDED"&sud$value_s27=="Y")
   remove_flags<-bc[sud$link[keep],,drop=FALSE];wr(remove_flags,"s27_added_AR_remove_flag.csv")
   fl<-do.call(rbind,lapply(c("S3","S27"),function(n){v<-ss[[n]];v<-v[v$QNAM=="REMOVEFL"&v$QVAL=="Y",];data.frame(snapshot=n,flagged_events=length(unique(v$parent_index)),participants=length(unique(v$.sid)))}));wr(fl,"AR_remove_flag_totals.csv")
   wr(sud[sud$subject%in%focus,,drop=FALSE],"focus_suppae_changes.csv")
 }
 cmap<-list();ledger<-list(); 
 for(k in seq_along(old)){i<-old[k];bi<-bidx[k];sid<-a$.sid[i];j<-am$to[i];via<-bm$to[bi];isU<-co[k]=="UNLINKED_SERIOUS"
  dct<-candidates(a,i,t,am$from);kk<-if(isU)q$rows[q$rows$eua_index==i,,drop=FALSE] else NULL
  paths<-if(isU)kk$bla_index else bi;pj<-bm$to[paths];possible<-unique(pj[!is.na(pj)])
  addr<-function(route,bi,j,reason,other=NA_integer_){
    n<-max(length(bi),length(j));if(n==0)return()
    cmap[[length(cmap)+1L]]<<-data.frame(subject=sid,eua_row=a$.row[i],eua_seq=a$.seq[i],route=route,s3_row=b$.row[bi],s3_seq=b$.seq[bi],s3_term=b$AETERM[bi],s3_onset=b$.date_raw[bi],s3_ser=ifelse(is.na(bi),NA,flag27(b$AESER[bi])),s27_row=t$.row[j],s27_seq=t$.seq[j],s27_term=t$AETERM[j],s27_onset=t$.date_raw[j],s27_ser=ifelse(is.na(j),NA,flag27(t$AESER[j])),reason=reason,eua_row_also_or_already_assigned=a$.row[other])}
  if(nrow(dct))addr("DIRECT_REVIEW_CANDIDATE",rep(NA_integer_,nrow(dct)),dct$j,dct$evidence,dct$assigned_from)
  if(!is.na(j))addr("DIRECT_ASSIGNED",NA_integer_,j,ac$rule[match(i,at$i)],i)
  if(length(paths))addr(if(isU)"VIA_S3_REVIEW_CANDIDATE" else "VIA_S3_ASSIGNED",paths,pj,if(isU)kk$reason else p$match_method[match(i,p$eua_index)])
  # If an S3 parent has no one-to-one S27 assignment, retain all its weak alternatives.
  if(length(paths))for(bi2 in paths[is.na(pj)]){wc<-candidates(b,bi2,t,bm$from);if(nrow(wc))addr("VIA_S3_REVIEW_ONLY_TO_S27",rep(bi2,nrow(wc)),wc$j,wc$evidence)}
  rn<-match(sid,roster$subject);tn<-roster$S27_AE_rows[rn]
  competing<-!is.na(j)&&length(possible)>0&&any(possible!=j)
  assess<-if(!tn)"NO_S27_AE_ROWS_FOR_SUBJECT" else if(!is.na(j)&&competing)"DIRECT_LINK_WITH_COMPETING_PATH" else if(!is.na(j))"DIRECT_CLINICAL_OR_CONTENT_LINK" else if(length(possible)&&isU)"INDIRECT_FROM_S3_CANDIDATE_ONLY" else if(!is.na(via))"TRACE_VIA_ASSIGNED_S3_ROW" else if(nrow(dct))"WEAK_DIRECT_CANDIDATE_ONLY" else "NO_DIRECT_CANDIDATE_OTHER_AE_PRESENT"
  es<-flag27(a$AESER[i]);bs<-if(is.na(bi))NA_character_ else flag27(b$AESER[bi]);ts<-if(is.na(via))NA_character_ else flag27(t$AESER[via])
  trajectory<-if(isU)"NOT_A_BASELINE_TRANSITION" else if(is.na(via))"S3_ROW_NOT_ASSIGNED_IN_S27" else if(!is.na(j)&&j!=via)"COMPETING_DIRECT_AND_CHAIN_LINK" else if(ts==bs)"S3_FLAG_RETAINED_ON_LINKED_ROW" else if(ts==es)"RETURN_TO_EUA_FLAG_ON_LINKED_ROW" else "OTHER_OR_MISSING_S27_FLAG"
  ledger[[k]]<-data.frame(subject=sid,arm=unname(arms[sid]),cohort=co[k],eua_row=a$.row[i],eua_seq=a$.seq[i],eua_onset=a$.date_raw[i],eua_window=a$.window[i],eua_term=a$AETERM[i],eua_ser=es,
    s3_row=b$.row[bi],s3_seq=b$.seq[bi],s3_term=b$AETERM[bi],s3_onset=b$.date_raw[bi],s3_ser=bs,s3_assign_rule=p$match_method[match(i,p$eua_index)],s3_candidates=oldq[k],s3_ae_rows=roster$S3_AE_rows[rn],s3_ADSL=roster$S3_ADSL_present[rn],s3_DM=roster$S3_DM_present[rn],s27_ae_rows=tn,s27_ADSL=roster$S27_ADSL_present[rn],s27_DM=roster$S27_DM_present[rn],
    s27_direct_row=t$.row[j],s27_direct_seq=t$.seq[j],direct_rule=ac$rule[match(i,at$i)],s27_direct_term=t$AETERM[j],s27_direct_onset=t$.date_raw[j],s27_ser_direct=if(is.na(j))NA_character_ else flag27(t$AESER[j]),
    weak_s27_candidates=nrow(dct),s27_via_s3_row=t$.row[via],s27_via_s3_seq=t$.seq[via],s27_ser_via_s3=ts,s27_rows_via_s3_candidates=if(isU)paste(t$.row[possible],collapse=";")else "",competing_destinations=competing,assessment=assess,trajectory=trajectory)
 }
 L<-emp(ledger,c("subject","arm","cohort","eua_term","eua_onset","s3_candidates","s3_ae_rows","s27_ae_rows","direct_rule","weak_s27_candidates","assessment","s3_term","eua_ser","s3_ser","s27_ser_via_s3","s27_ser_direct","trajectory"))
 if(nrow(L)){
   L$s3_verbatim_agrees<-ifelse(is.na(L$s3_row),NA,sc_norm(L$eua_term)==sc_norm(L$s3_term))
   L$s3_onset_agrees<-ifelse(is.na(L$s3_row)|!nzchar(L$eua_onset)|!nzchar(L$s3_onset),NA,L$eua_onset==L$s3_onset)
   L$direct_chain_agreement<-ifelse(is.na(L$s27_direct_row)|is.na(L$s27_via_s3_row),NA,L$s27_direct_row==L$s27_via_s3_row)
 }
 wr(L,"focus_case_ledger.csv");wr(emp(cmap,c("subject","eua_row","route","reason")),"candidate_paths.csv")
 review<-L;for(f in c("reviewer","review_date","reviewed_S3_row","reviewed_S27_row","decision","supporting_source_and_page","rationale"))review[[f]]<-rep("",nrow(L));wr(review,"manual_review_template.csv")
 zero<-L[L$cohort=="UNLINKED_SERIOUS"&!is.na(L$s3_candidates)&L$s3_candidates==0,c("subject","arm","eua_term","s3_ae_rows","s27_ae_rows","assessment")]
 wr(zero,"four_original_no_candidate_cases.csv")
 write_html27(out,baseline,L,fields,sus,pairsum,zero,coverage,trans,remove_flags,naive)
 wr(data.frame(phase=c("EUA","S3","S27"),submission=ids,scope=c("25 Nov 2020 safety onset boundary for unresolved SAEs","licensing S3 snapshot","later amendment; do not apply a Dec 2021 event-onset cutoff")),"snapshot_scope.csv")
 capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"));done<-TRUE
 ok<-nrow(check)==nrow(reference_s3)&&all(check$agrees)
 writeLines(if(ok)"COMPLETE: reference EUA/S3 case identities reproduced" else "COMPLETE WITH BASELINE DIFFERENCES: inspect baseline_identity_check.csv",status)
 message("Report: ",file.path(out,"report.html"));invisible(out)
}
if(isTRUE(getOption("phmpt.autorun",TRUE)))analyse_s27_adverse_events()
