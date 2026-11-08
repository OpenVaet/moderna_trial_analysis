# P301 audit: US3952070 and the seven randomized-flag removals.
# In RStudio set the working directory to the project and run the entire file.
# Requires original EUA S5 and BLA S3 DM, EX, ADSL, SUPPDM; S3 DV & SUPPDV.
# Optional: EUA DV. Use options(phmpt.root=..., phmpt.enrolment_inputs=characterVector)
# to search specific input directories; phmpt.enrolment_out sets the output directory.
# No downloads. Reads XPT with haven or foreign; includes a reader for the released JMP layouts.
# Matching is a screen for candidates, never a declaration that two IDs are one person.
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
    s <- rawToChar(if(any(z==as.raw(0))) z[seq_len(which(z==as.raw(0))[1]-1L)] else z[seq_len(tail(nonzero, 1L))])
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
audit_enrolment <- function(root=getOption('phmpt.root',getwd()),
 out=getOption('phmpt.enrolment_out',file.path(root,'analysis','enrolment_US3952070'))) {
 id <- 'US3952070'; dir.create(out,showWarnings=FALSE,recursive=TRUE)
 save_csv <- function(x,n) write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='')
 default <- c('prod_data','eua_data','prod_data_moderna','eua_data_moderna')
 dirs <- getOption('phmpt.enrolment_inputs',file.path(root,default))
 dirs <- dirs[dir.exists(dirs)]; if(!length(dirs)) dirs <- root
 files <- unique(unlist(lapply(dirs,list.files,recursive=TRUE,full.names=TRUE)))
 log <- list()
 get <- function(s,dom,fields,optional=FALSE) {
   subm <- if(s=='EUA') '27073_S5' else '125752_S3'
   key <- if(dom=='adsl') 'A(_D)?' else 'S'
   p <- files[grepl(paste0(subm,'_M5_mrna-1273-p301_',key,'_',dom,'\\.(xpt|jmp|csv)$'), basename(files),ignore.case=TRUE)]
   if(!length(p)){if(optional) return(NULL);stop('Missing ',s,' ',dom,'; set phmpt.enrolment_inputs to your extracted data directories.')}
   ext <- tolower(tools::file_ext(p)); p <- p[ext==c('xpt','jmp','csv')[min(match(ext,c('xpt','jmp','csv')))]]
   hh <- tools::md5sum(p); if(length(unique(hh))!=1L) stop('Conflicting copies of ',s,' ',dom,': narrow phmpt.enrolment_inputs.')
   f <- p[1]; message(s,' ',dom,': ',basename(f))
   if(grepl('jmp$',f,ignore.case=TRUE)) z<-ws_read_jmp(f,fields) else if(grepl('xpt$',f,ignore.case=TRUE)) {
     if(requireNamespace('foreign',quietly=TRUE)) z<-foreign::read.xport(f)
     else if(requireNamespace('haven',quietly=TRUE)) z<-as.data.frame(haven::read_xpt(f,encoding='WINDOWS-1252'))
     else stop('Install haven or foreign to read XPT.')
   } else {
     h<-readLines(f,n=1L,warn=FALSE); sep<-if(length(strsplit(h,';',fixed=TRUE)[[1]])>length(strsplit(h,',',fixed=TRUE)[[1]])) ';' else ','
     z<-read.table(f,header=TRUE,sep=sep,quote='"',colClasses='character',comment.char='',check.names=FALSE,fileEncoding='UTF-8-BOM')
   }
   names(z)<-toupper(names(z)); if(length(setdiff(fields,names(z)))) stop('Missing columns in ',basename(f),': ',paste(setdiff(fields,names(z)),collapse=', '))
   z<-z[,fields,drop=FALSE];z[]<-lapply(z,function(x){v<-as.character(x); bad<-!is.na(v)&!validUTF8(v); if(any(bad)) v[bad]<-iconv(v[bad],from='WINDOWS-1252',to='UTF-8',sub='byte');v<-trimws(v);v[!nzchar(v)]<-NA_character_;v})
   log[[paste(s,dom)]] <<- data.frame(snapshot=s,domain=dom,file=basename(f),md5=unname(hh[1]),rows=nrow(z))
   z
 }
 sid<-function(x) {if(any(!grepl('^mRNA-1273-P301-US[0-9]{3}-[0-9]+$',x))) stop('Unrecognized USUBJID');sub('^mRNA-1273-P301-(US[0-9]{3})-','\\1',x)}
 vals<-function(x)paste(sort(unique(na.omit(x))),collapse=';')
 eq<-function(a,b)!is.na(a)&!is.na(b)&a==b
 norm<-function(x)toupper(gsub('[^[:alnum:]]','',x))
 forms<-c('Visit 1 Day 1','Visit 2 Day 29')
 roster<-pres<-expos<-checks<-dc<-list() # entries are keyed by submission
 for(s in c('EUA','S3')){
   dm<-get(s,'dm',c('USUBJID','SUBJID','SITEID','AGE','SEX','RACE','ETHNIC','RFICDTC'))
   a<-get(s,'adsl',c('SUBJID','SITEID','AGE','SEX','RACE','HEIGHTBL','WEIGHTBL','TRT01P','RANDNUM','RANDFL','FASFL','SAFFL','MITTFL','PPROTFL','FASEREAS','SAFEREAS'))
   ex<-get(s,'ex',c('USUBJID','VISIT','EXTRT','EXLNKID','EXSTDTC','EXDOSE','EXDOSU'))
   sup<-get(s,'suppdm',c('USUBJID','QNAM','QVAL')); sup$SUBJID<-sid(sup$USUBJID)
   stopifnot(!anyDuplicated(dm$SUBJID),!anyDuplicated(a$SUBJID),all(dm$SUBJID==sid(dm$USUBJID)))
   ex$SUBJID<-sid(ex$USUBJID);expos[[s]]<-ex
   d<-merge(dm,a,all=TRUE,by='SUBJID',suffixes=c('','.ADSL'))
   extra<-setdiff(unique(ex$SUBJID),d$SUBJID)
   if(length(extra)) d<-merge(d,data.frame(SUBJID=extra),all=TRUE,by='SUBJID')
   for(k in c('SITEID','AGE','SEX','RACE')){
     f<-paste0(k,'.ADSL'); absent<-!d$SUBJID%in%dm$SUBJID;d[absent,k]<-d[absent,f]
   }
   d$AGE[!(grepl('^[0-9]+([.][0-9]+)?$',d$AGE))]<-NA_character_
   for(k in c('HEIGHTBL','WEIGHTBL')) {j<-!is.na(d[[k]]) & !grepl('^[0-9]+([.][0-9]+)?$',d[[k]]);d[[k]][j]<-NA_character_}
   d$in_DM<-d$SUBJID%in%dm$SUBJID;d$in_ADSL<-d$SUBJID%in%a$SUBJID;d$in_EX<-d$SUBJID%in%ex$SUBJID
   ti<-ex[ex$SUBJID==id & ex$VISIT%in%forms & !is.na(ex$EXTRT),]
   stopifnot(nrow(ti)==2L,length(unique(ti$VISIT))==2L)
   for(i in 1:2){
     e<-ex[ex$VISIT%in%forms[i] & !is.na(ex$EXTRT) & !is.na(ex$EXSTDTC),]
     e$date<-substr(e$EXSTDTC,1L,10L)
     dt0<-as.Date(substr(ti$EXSTDTC[ti$VISIT==forms[i]],1,10))
     e$delta<-as.numeric(as.Date(e$date)-dt0)
     sp<-split(seq_len(nrow(e)),e$SUBJID)
     d[[paste0('dose',i)]]<-unname(vapply(sp,function(j)vals(e$date[j]),character(1))[d$SUBJID])
     d[[paste0('offset',i)]]<-unname(vapply(sp,function(j)e$delta[j[which.min(abs(e$delta[j]))]],numeric(1))[d$SUBJID])
   }
   t<-d[d$SUBJID==id,];stopifnot(nrow(t)==1L)
   d$profile<-eq(d$AGE,t$AGE)&eq(d$SEX,t$SEX)&eq(d$RACE,t$RACE)
   d$site<-eq(d$SITEID,t$SITEID);d$other<-d$SUBJID!=id
   u<-d[d$other,];p<-u$profile; st<-u$site
   both<-function(w) !is.na(u$offset1)&!is.na(u$offset2)&abs(u$offset1)<=w&abs(u$offset2)<=w
   cn<-c(other_DM=sum(u$in_DM),other_union=nrow(u),other_site=sum(st),profile_unknown=sum(!complete.cases(u[,c('AGE','SEX','RACE')])),profile=sum(p),site_profile=sum(st&p),both_dates_any_profile=sum(both(0)),site_both_dates_any_profile=sum(st&both(0)),profile_exact=sum(p&both(0)),profile_within7=sum(p&both(7)),profile_within14=sum(p&both(14)),site_profile_within7=sum(st&p&both(7)),site_profile_within14=sum(st&p&both(14)),site_age_within1_sex_race=sum(st & abs(as.numeric(u$AGE)-as.numeric(t$AGE))<=1 & eq(u$SEX,t$SEX)&eq(u$RACE,t$RACE),na.rm=TRUE),profile_height3cm_weight3kg=sum(p & abs(as.numeric(u$HEIGHTBL)-as.numeric(t$HEIGHTBL))<=3 & abs(as.numeric(u$WEIGHTBL)-as.numeric(t$WEIGHTBL))<=3,na.rm=TRUE),other_same_randnum=sum(eq(u$RANDNUM,t$RANDNUM)),other_EX_ID_reusing_target_EXLNKID=length(unique(ex$SUBJID[ex$SUBJID!=id & ex$EXLNKID%in%na.omit(ti$EXLNKID)])))
   checks[[s]]<-data.frame(snapshot=s,screen=names(cn),other_IDs=unname(cn))
   prev<-sup[sup$QNAM=='SUBJID1' & !is.na(sup$QVAL),]
   rev<-prev[norm(prev$QVAL)%in%c(id,paste0('MRNA1273P301',id)) | (norm(prev$QVAL)=='2070' & substr(prev$SUBJID,1,5)=='US395'),]
   pres[[s]]<-data.frame(snapshot=s,target_PREVSCR=vals(sup$QVAL[sup$SUBJID==id & sup$QNAM=='PREVSCR']),target_previous_ID=vals(prev$QVAL[prev$SUBJID==id]),other_ID_naming_target=vals(rev$SUBJID[rev$SUBJID!=id]),total_nonblank_previous_ID_records=nrow(prev))
   roster[[s]]<-d
 }
 save_csv(do.call(rbind,checks),'linkage_counts'); save_csv(do.call(rbind,pres),'previous_screening_links')
 save_csv(do.call(rbind,lapply(names(roster),function(s){d<-roster[[s]];d<-d[d$profile,]; cbind(snapshot=s,d[,c('SUBJID','SITEID','AGE','SEX','RACE','RFICDTC','dose1','dose2','offset1','offset2','TRT01P','HEIGHTBL','WEIGHTBL','SAFFL')])})),'all_profile_candidates')
 save_csv(do.call(rbind,lapply(names(roster),function(s)cbind(snapshot=s,roster[[s]][roster[[s]]$SUBJID==id,c('SUBJID','RANDNUM','RANDFL','SAFFL','FASFL','MITTFL','PPROTFL','FASEREAS','SAFEREAS')]))),'index_flags')
 save_csv(do.call(rbind,lapply(names(expos),function(s)cbind(snapshot=s,expos[[s]][expos[[s]]$SUBJID==id,]))),'index_administrations')
 d<-roster$S3; a<-roster$EUA; ids<-union(d$SUBJID,a$SUBJID)
 save_csv(data.frame(EUA_ids=nrow(a),S3_ids=nrow(d),distinct_across_snapshots=length(ids),EUA_only=length(setdiff(a$SUBJID,d$SUBJID)),S3_only=length(setdiff(d$SUBJID,a$SUBJID))),'roster_scope')
 # Guide section 3.3: also reports US3902125, which is not one of the five Y->N removals.
 pairs<-data.frame(excluded=c('US3522027','US3352166','US3172220','US3162112','US3132036','US3902125'),
 counterpart=c('US3522029','US3172082','US3352442','US3732072','US3562079','US3902081'))
 l<-d[match(pairs$excluded,d$SUBJID),];r<-d[match(pairs$counterpart,d$SUBJID),]
 pairs$EUA_flag<-a$RANDFL[match(pairs$excluded,a$SUBJID)];pairs$S3_flag<-l$RANDFL
 for(k in c('SITEID','AGE','SEX','RACE','dose1','dose2','TRT01P')){
    pairs[[paste0(k,'_excluded')]]<-l[[k]];pairs[[paste0(k,'_counterpart')]]<-r[[k]]
 }
 pairs$profile_agrees<-eq(l$AGE,r$AGE)&eq(l$SEX,r$SEX)&eq(l$RACE,r$RACE)
 pairs$both_dates_agree<-eq(l$dose1,r$dose1)&eq(l$dose2,r$dose2)
 save_csv(pairs,'guide_pairs')
 ds<-get('S3','dv',c('USUBJID','SUBJID','DVSEQ','DVSTDTC','DVTERM','DVSPID'))
 q<-get('S3','suppdv',c('USUBJID','IDVARVAL','QNAM','QVAL'))
 key<-c('USUBJID','DVSEQ');q$DVSEQ<-as.character(as.integer(q$IDVARVAL)); ds$DVSEQ<-as.character(as.integer(ds$DVSEQ))
 stopifnot(!anyDuplicated(ds[key]))
 for(k in c('SEVDES','DECSRTNL','DECSRTN1','DECSRTN2','DECSRTN3','DVSIG','DVTERM1')){
   f<-q[q$QNAM==k,c(key,'QVAL')];stopifnot(!anyDuplicated(f[key]));names(f)[3]<-k
   ds<-merge(ds,f,by=key,all.x=TRUE,sort=FALSE)
 }
 ds$code<-sub(':.*','',ds$DVTERM)
 ds$severity_class<-ifelse(is.na(ds$SEVDES),'missing',ds$SEVDES)
 sc<-c('P07.04','P10.01'); z<-ds[ds$code%in%sc,]
 save_csv(as.data.frame(table(code=z$code,severity=z$severity_class),responseName='records'),'same_code_decisions')
 intersect_ids<-Reduce(intersect,lapply(sc,function(k)z$SUBJID[z$code==k]))
 save_csv(as.data.frame(table(RANDFL=factor(d$RANDFL[match(intersect_ids,d$SUBJID)],levels=c('Y','N')),useNA='ifany'),responseName='IDs'),'both_codes_randomization')
 fchange<-merge(a[,c('SUBJID','RANDFL')],d[,c('SUBJID','RANDFL')],by='SUBJID',suffixes=c('_EUA','_S3'))
 changed<-fchange$SUBJID[which(fchange$RANDFL_EUA=='Y' & fchange$RANDFL_S3=='N')]
 save_csv(ds[ds$SUBJID%in%changed & grepl('^1 *-',ds$severity_class),],'seven_removals_code1_records')
 save_csv(ds[ds$SUBJID==id,],'index_S3_deviations')
 save_csv(ds[ds$SUBJID %in% d$SUBJID[d$profile] & grepl('P04[.]0[23]|two sites|duplicate|enrolled twice|randomized in error',ds$DVTERM,ignore.case=TRUE),],'profile_explicit_duplicate_deviations')
 eu<-get('EUA','dv',c('SUBJID','DVSEQ','DVSTDTC','DVTERM'),TRUE)
 if(!is.null(eu)) save_csv(eu[eu$SUBJID==id,],'index_EUA_deviations')
 save_csv(do.call(rbind,log),'source_manifest')
 message('Completed: ',normalizePath(out))
 invisible(out)
}
if(isTRUE(getOption('phmpt.autorun',TRUE))) audit_enrolment()
