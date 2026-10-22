# =============================================================================
# Moderna P301 - efficacy case reconstruction and interim-analysis timeline
# Version 1.0.0 | 2026-10-22 | Standalone R >= 4.3 (select all / Run in RStudio)
# Run from the project containing eua_data/ and prod_data/, or set PHMPT_ROOT.
# Reads local XPT (preferred), SAS7BDAT or CSV. Never downloads trial datasets.
# For XPT/SAS7BDAT only: install.packages('haven') once if not already installed.
# Base R produces CSV, HTML, PDF and PNG. No ggplot2 / pandoc / browser required.
# Load functions only: options(phmpt.autorun=FALSE) before source().
#
# PURPOSE: recover the sponsor's adjudicated primary-endpoint case membership,
# retain all those cases regardless of PCR-link success, and expose the evidence.
# Case-event dates are NOT result receipt / adjudication / DSMB readiness dates.
# The later EUA/BLA snapshots do NOT recreate the original interim database.
# No alternative vaccine efficacy or political-motive test is calculated.
# =============================================================================
PHMPT_ROOT <- getOption('phmpt.root', '')
PHMPT_OUTPUT_DIR <- getOption('phmpt.output_dir', '')
PHMPT_RELEASES <- c('EUA', 'BLA')
PHMPT_SUBMISSIONS <- c(EUA='27073_S5', BLA='125752_S3')
PHMPT_FILES <- list(EUA=list(TTE='', ADEFF='', MB='', ADMB='', CE='', SUPPCE=''),
                    BLA=list(TTE='', ADEFF='', MB='', ADMB='', CE='', SUPPCE=''))
PHMPT_CUTOFFS <- c(EUA='2020-11-21', BLA='2021-03-26') # efficacy, not safety cutoff
PHMPT_EXPECTED <- data.frame(release=c('EUA','EUA','BLA','BLA'),
  arm=rep(c('mRNA-1273','Placebo'),2), population=c(14134L,14073L,14287L,14164L),
  cases=c(11L,185L,55L,744L), stringsAsFactors=FALSE)
# Metadata-backed exact endpoint titles. Codes are checked against their titles.
PHMPT_PARAMETERS <- c(EUA='Time to COVID-19 14 Days after Dose 2 (Days)',
 BLA='Time to COVID-19 (Prior Infection Censoring) 14 Days after Dose 2 (Days)')
PHMPT_ALLOW_CSV <- TRUE # fallback only when no original XPT/SAS exists
PHMPT_CSV_ENCODING <- 'UTF-8' # change to 'windows-1252' for an old Windows export
PHMPT_CSV_DATE_ORIGIN <- '' # numeric CSV dates require an explicit origin, e.g. '1960-01-01'
PHMPT_CSV_DATE_FORMAT <- '' # ambiguous dates require e.g. '%d/%m/%Y'; never guessed
PHMPT_PCR_CODES <- c('SARSCOV2')
PHMPT_REVIEW_WINDOW_DAYS <- 14L # candidate search only; NOT an endpoint eligibility rule
PHMPT_AVAILABILITY_FILE <- '' # optional reviewed CSV; template is written on every run
PHMPT_SHOW_REPORT <- interactive()
PHMPT_PNG_DPI <- 150L

# ---- Small utilities ---------------------------------------------------------
ec_text <- function(x) {
 if(inherits(x,'Date'))return(ifelse(is.na(x),'',format(x,'%Y-%m-%d')))
 if(inherits(x,'POSIXt'))return(ifelse(is.na(x),'',format(x,'%Y-%m-%dT%H:%M:%S',tz='UTC')))
 z<-as.character(x);z[is.na(z)]<-''
 bad<-!validUTF8(z);if(any(bad))z[bad]<-iconv(z[bad],from='windows-1252',to='UTF-8',sub=NA)
 if(anyNA(z))stop('Unrecoverable text encoding. Set PHMPT_CSV_ENCODING explicitly.')
 trimws(enc2utf8(z))
}
ec_norm <- function(x)toupper(gsub('[[:space:]]+',' ',ec_text(x)))
ec_col <- function(x,n)if(n%in%names(x))ec_text(x[[n]]) else rep('',nrow(x))
ec_first <- function(x,fields){v<-rep('',nrow(x));for(n in fields){z<-ec_col(x,n);i<-!nzchar(v);v[i]<-z[i]};v}
ec_unique <- function(x)paste(sort(unique(x[!is.na(x)&nzchar(x)])),collapse=' | ')
ec_sid <- function(x){x<-ec_norm(x);i<-grepl('US[0-9]{3}-?[0-9]{4}$',x);x[i]<-gsub('-','',sub('^.*(US[0-9]{3}-?[0-9]{4})$','\\1',x[i]));x}
ec_key <- function(...)do.call(paste,c(list(...),sep='|'))
ec_arm <- function(x){z<-ec_norm(x);out<-rep('UNKNOWN',length(z));out[z=='PLACEBO']<-'Placebo';out[z%in%c('MRNA-1273','MRNA-1273 100 MCG','MRNA-1273 100 UG')]<-'mRNA-1273';out}
ec_date <- function(x){
 x<-ec_text(x);d<-rep(as.Date(NA),length(x));iso<-grepl('^[0-9]{4}-[0-9]{2}-[0-9]{2}($|[ T])',x)
 d[iso]<-suppressWarnings(as.Date(substr(x[iso],1,10),format='%Y-%m-%d'))
 bad<-iso&!is.na(d)&format(d,'%Y-%m-%d')!=substr(x,1,10);d[bad]<-as.Date(NA)
 sas<-grepl('^[0-9]{1,2}[- ]?[A-Za-z]{3}[- ]?[0-9]{4}$',x)
 if(any(sas)){s<-toupper(gsub('[- ]','',x[sas]));day<-sub('([0-9]{1,2}).*','\\1',s);mon<-match(sub('^[0-9]{1,2}([A-Z]{3}).*','\\1',s),toupper(month.abb));yr<-sub('.*([0-9]{4})$','\\1',s);d[sas]<-suppressWarnings(as.Date(sprintf('%s-%02d-%02d',yr,mon,as.integer(day)),format='%Y-%m-%d'))}
 num<-grepl('^-?[0-9]+([.]0+)?$',x)&nzchar(x)
 if(any(num)&&nzchar(PHMPT_CSV_DATE_ORIGIN))d[num]<-as.Date(as.numeric(x[num]),origin=PHMPT_CSV_DATE_ORIGIN)
 left<-is.na(d)&nzchar(x)
 if(any(left)&&nzchar(PHMPT_CSV_DATE_FORMAT))d[left]<-suppressWarnings(as.Date(x[left],format=PHMPT_CSV_DATE_FORMAT))
 d
}
ec_csv <- function(x,path){con<-file(path,'wb');on.exit(close(con));writeBin(as.raw(c(239,187,191)),con);write.table(x,con,sep=',',quote=TRUE,row.names=FALSE,na='',qmethod='double')}
ec_bind <- function(xs){xs<-Filter(function(x)!is.null(x)&&nrow(x)>0,xs);if(!length(xs))return(data.frame());ns<-unique(unlist(lapply(xs,names)));xs<-lapply(xs,function(x){for(n in setdiff(ns,names(x)))x[[n]]<-NA;x[ns]});do.call(rbind,xs)}
ec_path <- function(root,p)if(grepl('^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)',p))p else file.path(root,p)
ec_root <- function(root=''){
 if(nzchar(root))return(normalizePath(root,winslash='/',mustWork=TRUE))
 cand<-getwd();for(f in sys.frames())if(!is.null(f$ofile))cand<-c(cand,dirname(f$ofile))
 args<-grep('^--file=',commandArgs(),value=TRUE);if(length(args))cand<-c(cand,dirname(sub('^--file=','',args[1])))
 if(interactive()&&requireNamespace('rstudioapi',quietly=TRUE)&&rstudioapi::isAvailable())cand<-c(cand,dirname(rstudioapi::getActiveDocumentContext()$path))
 for(p in unique(cand[dir.exists(cand)]))for(i in 1:8){if(any(dir.exists(file.path(p,c('eua_data','prod_data')))))return(normalizePath(p,winslash='/'));q<-dirname(p);if(q==p)break;p<-q}
 stop('Project not found. Set PHMPT_ROOT to the folder containing eua_data and/or prod_data.')
}
ec_need <- function(x,vars,what){miss<-setdiff(vars,names(x));if(length(miss))stop(what,' missing: ',paste(miss,collapse=', '))}
ec_keep <- list(
 TTE=c('STUDYID','USUBJID','SUBJID','SITEID','TRT01P','TRT01PN','PPROTFL','SCOV2BL','RTPCRBL','ELECSBL','DOSE2DT','AP01EDT','EFFCODT','CUTOFFDT','PARCAT1','PARAM','PARAMCD','ADT','STARTDT','CNSR','EVNTDESC','APERIOD','CVD53DT','ANL01FL','SRCDOM','SRCSEQ','SRCVAR'),
 ADEFF=c('STUDYID','USUBJID','SUBJID','PARCAT1','PARAM','PARAMCD','ADT','PCRDT','AVALC','ANL01FL','APERIOD','PCRD29FL','EVIDINF','ASEQ','SRCDOM','SRCSEQ','SRCVAR','MBSEQ'),
 MB=c('STUDYID','USUBJID','SUBJID','MBSEQ','MBREFID','MBSPID','MBGRPID','MBTESTCD','MBTEST','MBDTC','MBSTRESC','MBORRES','MBSTAT','MBREASND','MBNAM','MBMETHOD','MBSPEC','VISIT','VISITNUM'),
 ADMB=c('STUDYID','USUBJID','SUBJID','MBSEQ','ASEQ','PARAMCD','PARAM','MBDTC','ADT','AVALC','MBSTAT','MBMETHOD','MBSPEC','AVISIT','VISIT','LOCALFL','DTYPE','SRCDOM','SRCSEQ','SRCVAR'))
# CE/SUPPCE are deliberately read with every column, to retain adjudication detail.
ec_read <- function(path,domain,subjects=NULL){
 ext<-tolower(tools::file_ext(path));fields<-ec_keep[[domain]]
 if(ext%in%c('xpt','sas7bdat')){
  if(!requireNamespace('haven',quietly=TRUE))stop('Install haven once: install.packages("haven")')
  wanted<-unique(c(fields,tolower(fields)))
  if(ext=='xpt')x<-if(is.null(fields))haven::read_xpt(path) else haven::read_xpt(path,col_select=tidyselect::any_of(wanted))
  else x<-if(is.null(fields))haven::read_sas(path) else haven::read_sas(path,col_select=tidyselect::any_of(wanted))
 }else if(ext=='csv'){
  first<-readLines(path,n=1,warn=FALSE,encoding=PHMPT_CSV_ENCODING)
  sep<-if(length(strsplit(first,';',fixed=TRUE)[[1]])>length(strsplit(first,',',fixed=TRUE)[[1]]))';' else ','
  x<-withCallingHandlers(read.table(path,header=TRUE,sep=sep,quote='"',comment.char='',check.names=FALSE,
    colClasses='character',na.strings=character(),fileEncoding=if(PHMPT_CSV_ENCODING=='UTF-8')'UTF-8-BOM' else PHMPT_CSV_ENCODING),
    warning=function(w){if(grepl('invalid|multibyte|EOF within|embedded nul',conditionMessage(w),ignore.case=TRUE))
      stop('CSV read warning: ',conditionMessage(w),'. Check PHMPT_CSV_ENCODING / source integrity; partial input is not accepted.')})
 }else stop('Unsupported file type: ',path)
 x<-as.data.frame(x,stringsAsFactors=FALSE);names(x)<-toupper(sub('^\ufeff','',names(x)))
 if(anyDuplicated(names(x)))stop('Duplicate column names: ',basename(path))
 ec_need(x,c('STUDYID','USUBJID'),domain);x$.row<-seq_len(nrow(x));total<-nrow(x)
 x<-x[ec_norm(x$STUDYID)=='MRNA-1273-P301',,drop=FALSE]
 x$.sid<-ec_sid(x$USUBJID)
 if(any(!grepl('^US[0-9]{7}$',x$.sid)))stop('Unexpected P301 subject identifiers in ',basename(path))
 if(!is.null(subjects))x<-x[x$.sid%in%subjects,,drop=FALSE]
 if(!is.null(fields))x<-x[intersect(c(fields,'.row','.sid'),names(x))]
 meta<-data.frame(variable=names(x),label=vapply(x,function(v){z<-attr(v,'label');if(is.null(z))'' else ec_text(z)},''),r_class=vapply(x,function(v)paste(class(v),collapse='/'),''))
 for(n in intersect(names(x),c('ADT','STARTDT','PCRDT','DOSE2DT','AP01EDT','EFFCODT','CUTOFFDT','CVD53DT'))){
  if(ext!='csv'&&is.numeric(x[[n]])&&!inherits(x[[n]],c('Date','POSIXt')))x[[n]]<-as.Date(as.numeric(x[[n]]),origin='1960-01-01')
 }
 for(n in setdiff(names(x),'.row'))x[[n]]<-ec_text(x[[n]])
 list(data=x,meta=meta,total_rows=total,path=normalizePath(path,winslash='/'),md5=unname(tools::md5sum(path)))
}
ec_pick <- function(root,release,domain,out,subjects=NULL,required=TRUE){
 override<-PHMPT_FILES[[release]][[domain]];if(is.null(override))override<-''
 folder<-file.path(root,if(release=='EUA')'eua_data' else 'prod_data')
 p<-if(nzchar(override))ec_path(root,override) else list.files(folder,recursive=TRUE,full.names=TRUE,pattern='\\.(xpt|sas7bdat|csv)$',ignore.case=TRUE)
 if(!nzchar(override)){
  domain_pat<-if(domain=='TTE')if(release=='EUA')'(adtt|adttea)' else 'adttea' else tolower(domain)
  p<-p[grepl(paste0('_',domain_pat,'\\.(xpt|sas7bdat|csv)$'),basename(p),ignore.case=TRUE)&grepl('mrna-1273-p301_',basename(p),ignore.case=TRUE)]
 }
 log<-data.frame(path=p,status=rep('CANDIDATE',length(p)),md5=rep('',length(p)),stringsAsFactors=FALSE)
 write_log<-function()ec_csv(log,file.path(out,paste0('discovery_',domain,'.csv')))
 if(!length(p)){write_log();if(required)stop('No ',release,' ',domain,' found. Set PHMPT_FILES$',release,'$',domain);return(NULL)}
 use<-rep(TRUE,length(p))
 if(!nzchar(override)){
  use<-grepl(paste0(PHMPT_SUBMISSIONS[[release]],'_'),basename(p),ignore.case=TRUE);log$status[!use]<-'OTHER_SUBMISSION'
  if(!PHMPT_ALLOW_CSV){log$status[use&tolower(tools::file_ext(p))=='csv']<-'CSV_DISABLED';use<-use&tolower(tools::file_ext(p))!='csv'}
  for(ext in c('xpt','sas7bdat'))if(any(use&tolower(tools::file_ext(p))==ext)){z<-use&tolower(tools::file_ext(p))!=ext;log$status[z]<-'PREFER_ORIGINAL';use[z]<-FALSE;break}
 }
 if(!any(use)){write_log();if(required)stop('No ',release,' ',domain,' in selected submission. See discovery log.');return(NULL)}
 ix<-which(use);log$md5[ix]<-unname(tools::md5sum(p[ix]));ids<-unique(log$md5[ix])
 if(anyNA(ids)||any(!nzchar(ids))){write_log();stop('Unreadable input path for ',domain)}
 if(length(ids)>1){log$status[ix]<-'MULTIPLE_DISTINCT_FILES';write_log();stop('Distinct ',release,' ',domain,' snapshots remain. Set PHMPT_FILES$',release,'$',domain,'. Files are never concatenated.')}
 sel<-ix[1];log$status[ix]<-'BYTE_IDENTICAL_COPY';log$status[sel]<-'SELECTED';write_log()
 message('Reading ',release,' ',domain,': ',basename(p[sel]))
 z<-ec_read(p[sel],domain,subjects);ec_csv(z$meta,file.path(out,paste0('metadata_',domain,'.csv')))
 z$manifest<-data.frame(release=release,domain=domain,path=z$path,md5=z$md5,rows_in_file=z$total_rows,rows_retained=nrow(z$data),stringsAsFactors=FALSE);z
}

# ---- Select primary adjudicated population and cases --------------------------
ec_cases <- function(t,release,out){
 ec_need(t,c('PARAM','PARAMCD','PARCAT1','PPROTFL','TRT01P','ADT','CNSR','DOSE2DT'),'TTE')
 inv<-aggregate(rep(1L,nrow(t)),t[c('PARAMCD','PARAM','PARCAT1')],sum);names(inv)[4]<-'rows'
 ec_csv(inv,file.path(out,'endpoint_inventory.csv'))
 keep<-ec_norm(t$PARAM)==ec_norm(PHMPT_PARAMETERS[[release]]) & ec_norm(t$PARCAT1)=='DERIVED PER ADJUDICATION'
 if(release=='BLA'&&'APERIOD'%in%names(t))keep<-keep&ec_norm(t$APERIOD)=='1'
 q<-t[keep,,drop=FALSE]
 if(!nrow(q))stop('Expected adjudicated endpoint title not found. Review endpoint_inventory.csv / PHMPT_PARAMETERS.')
 if(length(unique(q$PARAMCD))!=1)stop('Several codes share the selected endpoint title; review metadata.')
 dup<-duplicated(q$.sid)|duplicated(q$.sid,fromLast=TRUE)
 ec_csv(q[dup,,drop=FALSE],file.path(out,'duplicate_endpoint_subjects.csv'))
 if(any(dup))stop('Selected TTE endpoint has duplicate subject records. Resolve snapshot/period; no first-row selection.')
 q$arm<-ec_arm(q$TRT01P);q$release<-release
 if(any(q$arm=='UNKNOWN'))stop('Unrecognized planned arm label; inspect TTE.TRT01P.')
 if(any(!ec_norm(q$PPROTFL)%in%c('Y','N')))stop('Missing/invalid per-protocol flag in selected TTE endpoint.')
 if(any(!q$CNSR%in%c('0','1')))stop('Unexpected censoring codes; inspect CNSR/EVNTDESC before selecting cases.')
 q$event_date<-ec_date(q$ADT);q$dose2_date<-ec_date(q$DOSE2DT)
 p<-q[ec_norm(q$PPROTFL)=='Y',,drop=FALSE];cases<-p[p$CNSR=='0',,drop=FALSE]
 cases$case_id<-paste(release,cases$.sid,sep=':');cases$site<-sub('^(US[0-9]{3}).*','\\1',cases$.sid)
 cases$days_after_dose2<-as.integer(cases$event_date-cases$dose2_date)
 cases$case_date_valid<-!is.na(cases$event_date)
 cases$within_cutoff<-!is.na(cases$event_date)&cases$event_date<=as.Date(PHMPT_CUTOFFS[[release]])
 cases$after_dose2_14d<-!is.na(cases$days_after_dose2)&cases$days_after_dose2>=14
 cases$baseline_negative<-if('SCOV2BL'%in%names(cases))ec_norm(cases$SCOV2BL)=='NEGATIVE' else NA
 ec_csv(p,file.path(out,'primary_population.csv'));ec_csv(cases,file.path(out,'cases_before_pcr_linkage.csv'))
 bad<-!cases$case_date_valid|!cases$within_cutoff|!cases$after_dose2_14d|(!is.na(cases$baseline_negative)&!cases$baseline_negative)
 ec_csv(cases[bad,,drop=FALSE],file.path(out,'case_membership_checks.csv'))
 if(any(!cases$case_date_valid))stop('A counted case has an unreadable event date. Review CSV date settings; timeline cannot be computed.')
 checks<-PHMPT_EXPECTED[PHMPT_EXPECTED$release==release,,drop=FALSE]
 checks$observed_population<-vapply(checks$arm,function(a)sum(p$arm==a),integer(1));checks$observed_cases<-vapply(checks$arm,function(a)sum(cases$arm==a),integer(1))
 checks$matches_published<-checks$population==checks$observed_population&checks$cases==checks$observed_cases
 checks$case_rule_flags<-vapply(checks$arm,function(a)sum(bad&cases$arm==a),integer(1))
 ec_csv(checks,file.path(out,'published_count_checks.csv'))
 # Preserve source CVD53DT wherever populated, including dedicated 53-case parameters.
 cd<-if('CVD53DT'%in%names(t))unique(t[c('PARAMCD','PARAM','CVD53DT')]) else data.frame(PARAMCD=character(),PARAM=character(),CVD53DT=character())
 cd<-cd[nzchar(cd$CVD53DT),,drop=FALSE];ec_csv(cd,file.path(out,'source_CVD53DT.csv'))
 list(cases=cases,population=p,checks=checks,source53=cd,valid=all(checks$matches_published)&!any(bad))
}

# ---- Reconcile source PCR evidence without collapsing uncertain identities -----
ec_result <- function(x){
 z<-ec_norm(x);v<-rep('OTHER',length(z));v[!nzchar(z)]<-'MISSING'
 v[z%in%c('POSITIVE','P','DETECTED','RNA DETECTED','SARS-COV-2 DETECTED','PRESUMPTIVE POSITIVE','Y')]<-'POSITIVE'
 v[z%in%c('NEGATIVE','N','NOT DETECTED','RNA NOT DETECTED','SARS-COV-2 NOT DETECTED')]<-'NEGATIVE'
 v[z%in%c('NOT DONE','NOT PERFORMED','CANCELLED','CANCELED')]<-'NOT_DONE';v
}
ec_tests <- function(m,a,out){
 ec_need(m,c('MBSEQ','MBTESTCD','MBDTC','MBSTRESC'),'MB')
 m<-m[ec_norm(m$MBTESTCD)%in%PHMPT_PCR_CODES,,drop=FALSE]
 if(is.null(a))a<-m[FALSE,c('STUDYID','USUBJID','.sid','.row'),drop=FALSE]
 else {ec_need(a,c('PARAMCD','MBSEQ','AVALC'),'ADMB');a<-a[ec_norm(a$PARAMCD)%in%PHMPT_PCR_CODES,,drop=FALSE]}
 ignored<-nzchar(ec_col(a,'DTYPE'))|!nzchar(ec_col(a,'MBSEQ'))
 ec_csv(a[ignored,,drop=FALSE],file.path(out,'ADMB_derived_or_unkeyed.csv'));a<-a[!ignored,,drop=FALSE]
 mk<-ec_key(m$.sid,ec_col(m,'MBSEQ'));ak<-ec_key(a$.sid,ec_col(a,'MBSEQ'))
 # Keep duplicate-key rows; never force an arbitrary many-to-many match.
 md<-duplicated(mk)|duplicated(mk,fromLast=TRUE);ad<-duplicated(ak)|duplicated(ak,fromLast=TRUE)
 ec_csv(m[md,,drop=FALSE],file.path(out,'MB_duplicate_keys.csv'));ec_csv(a[ad,,drop=FALSE],file.path(out,'ADMB_duplicate_keys.csv'))
 ai<-match(mk,ak);ai[md|mk%in%ak[ad]|!nzchar(ec_col(m,'MBSEQ'))]<-NA_integer_
 av<-function(n){v<-ec_col(a,n)[ai];v[is.na(v)]<-'';v}
 r<-data.frame(test_id=paste0('MB:',m$.row),subject=m$.sid,USUBJID=m$USUBJID,MBSEQ=ec_col(m,'MBSEQ'),
   MBREFID=ec_col(m,'MBREFID'),MBGRPID=ec_col(m,'MBGRPID'),test_code=ec_col(m,'MBTESTCD'),test_name=ec_col(m,'MBTEST'),
   source_collection=ec_col(m,'MBDTC'),analysis_collection=av('MBDTC'),analysis_date=av('ADT'),
   raw_standard_result=ec_col(m,'MBSTRESC'),raw_original_result=ec_col(m,'MBORRES'),analysis_result=av('AVALC'),
   raw_status=ec_col(m,'MBSTAT'),analysis_status=av('MBSTAT'),method=ec_col(m,'MBMETHOD'),analysis_method=av('MBMETHOD'),
   laboratory=ec_col(m,'MBNAM'),specimen=ec_col(m,'MBSPEC'),visit=ec_col(m,'VISIT'),analysis_visit=av('AVISIT'),local_flag=av('LOCALFL'),
   MB_source_row=m$.row,ADMB_source_row=ifelse(is.na(ai),NA,a$.row[ai]),source_type=rep('MB',nrow(m)),
   key_ambiguous=md|mk%in%ak[ad],stringsAsFactors=FALSE)
 orphan<-!ak%in%mk;ec_csv(a[orphan,,drop=FALSE],file.path(out,'ADMB_without_source_MB.csv'))
 if(any(orphan)){
  z<-a[orphan,,drop=FALSE];extra<-data.frame(test_id=paste0('ADMB_ONLY:',z$.row),subject=z$.sid,USUBJID=z$USUBJID,
   MBSEQ=ec_col(z,'MBSEQ'),MBREFID='',MBGRPID='',test_code=ec_col(z,'PARAMCD'),test_name=ec_col(z,'PARAM'),
   source_collection='',analysis_collection=ec_col(z,'MBDTC'),analysis_date=ec_col(z,'ADT'),raw_standard_result='',raw_original_result='',
   analysis_result=ec_col(z,'AVALC'),raw_status='',analysis_status=ec_col(z,'MBSTAT'),method='',analysis_method=ec_col(z,'MBMETHOD'),
   laboratory='',specimen=ec_col(z,'MBSPEC'),visit='',analysis_visit=ec_col(z,'AVISIT'),local_flag=ec_col(z,'LOCALFL'),
   MB_source_row=NA_integer_,ADMB_source_row=z$.row,source_type='ADMB_ONLY',key_ambiguous=ad[orphan],stringsAsFactors=FALSE)
  r<-rbind(r,extra)
 }
 r$collection_date<-ec_date(r$source_collection);da<-ec_date(r$analysis_collection)
 r$date_source<-ifelse(r$source_type=='MB','MB.MBDTC','ADMB.MBDTC')
 r$collection_date[r$source_type=='ADMB_ONLY']<-da[r$source_type=='ADMB_ONLY']
 # No ADT fallback when a raw collection date is absent/partial.
 r$date_conflict<-!is.na(r$collection_date)&!is.na(da)&r$collection_date!=da
 c1<-ec_result(r$raw_standard_result);c2<-ec_result(r$raw_original_result);c3<-ec_result(r$analysis_result)
 r$positive_any<-c1=='POSITIVE'|c2=='POSITIVE'|c3=='POSITIVE';r$negative_any<-c1=='NEGATIVE'|c2=='NEGATIVE'|c3=='NEGATIVE'
 r$not_done<-ec_norm(r$raw_status)%in%c('NOT DONE','NOT PERFORMED','CANCELLED','CANCELED')|ec_norm(r$analysis_status)%in%c('NOT DONE','NOT PERFORMED','CANCELLED','CANCELED')|c1=='NOT_DONE'|c2=='NOT_DONE'|c3=='NOT_DONE'
 r$result_conflict<-(r$positive_any&r$negative_any)|(r$not_done&(r$positive_any|r$negative_any))
 meth<-paste(ec_norm(r$method),ec_norm(r$analysis_method))
 r$method_documented_pcr<-grepl('PCR|POLYMERASE',meth)&!grepl('ANTIGEN|ANTIBODY|SEROLOG',meth)
 r$raw_positive<-c1=='POSITIVE'|c2=='POSITIVE'
 r$usable_positive<-r$positive_any&r$method_documented_pcr&!r$result_conflict&!r$not_done&!r$key_ambiguous&!r$date_conflict&!is.na(r$collection_date)
 r$raw_positive_supported<-r$usable_positive&r$raw_positive&r$source_type=='MB'
 ec_csv(r,file.path(out,'case_subject_pcr_records.csv'));r
}
ec_link <- function(cases,e,tests,out){
 ec_need(e,c('PARAMCD','PARCAT1','ADT','PCRDT','AVALC'),'ADEFF')
 e$.date<-ec_date(e$ADT);e$.pcr<-ec_date(e$PCRDT)
 ec_csv(e,file.path(out,'case_subject_ADEFF_records.csv'))
 ee<-split(seq_len(nrow(e)),e$.sid);tt<-split(seq_len(nrow(tests)),tests$subject)
 evidence<-vector('list',nrow(cases));edges<-vector('list',nrow(cases));summary<-vector('list',nrow(cases))
 for(i in seq_len(nrow(cases))){
  c<-cases[i,,drop=FALSE];ei<-ee[[c$.sid]];if(is.null(ei))ei<-integer();ti<-tt[[c$.sid]];if(is.null(ti))ti<-integer()
  z<-e[ei,,drop=FALSE];p<-tests[ti,,drop=FALSE]
  yes<-ec_norm(z$AVALC)=='Y';same<-!is.na(z$.date)&z$.date==c$event_date
  adj<-yes&ec_norm(z$PARCAT1)=='DERIVED PER ADJUDICATION'&ec_norm(z$PARAMCD)=='COVIDA'
  sap<-yes&ec_norm(z$PARCAT1)=='DERIVED PER SAP'&ec_norm(z$PARAMCD)=='COVID'
  direct<-unique(z$.pcr[adj&same&!is.na(z$.pcr)]);support<-unique(z$.pcr[sap&same&!is.na(z$.pcr)])
  # An ADEFF source sequence is used only if it explicitly names MB and MBSEQ.
  explicit_seqs<-character()
  ix<-which(adj&same)
  if(length(ix)&&all(c('SRCDOM','SRCSEQ','SRCVAR')%in%names(z)))explicit_seqs<-z$SRCSEQ[ix[ec_norm(z$SRCDOM[ix])=='MB'&ec_norm(z$SRCVAR[ix])=='MBSEQ']]
  seqmatch<-nzchar(p$MBSEQ)&p$MBSEQ%in%explicit_seqs
  dm<-!is.na(p$collection_date)&p$collection_date%in%direct;sm<-!is.na(p$collection_date)&p$collection_date%in%support
  delta<-as.integer(p$collection_date-c$event_date);near<-!is.na(delta)&abs(delta)<=PHMPT_REVIEW_WINDOW_DAYS
  # Include negative/conflicting records on candidate dates to expose disagreements.
  take<-seqmatch|dm|sm|near
  if(any(take)){
   x<-p[take,,drop=FALSE];x$case_id<-c$case_id;x$case_event_date<-c$event_date;x$days_from_case<-delta[take]
   x$link_basis<-ifelse(seqmatch[take],'EXPLICIT_ADJUDICATION_MB_SEQUENCE',ifelse(dm[take],'ADJUDICATED_PCRDT_DATE',ifelse(sm[take],'SAME_EVENT_SAP_PCRDT_DATE','NEARBY_DATE_REVIEW_ONLY')))
   x$qualifying_specimen_verified<-seqmatch[take]&x$raw_positive_supported
   x$link_caveat<-ifelse(x$qualifying_specimen_verified,'Explicit source link plus compatible positive result; verify derivation metadata','Date association is not proof that this specimen qualified the adjudicated case')
   edges[[i]]<-x
  }
  if(nrow(z)){z$case_id<-c$case_id;z$exact_event_date<-same;z$role<-ifelse(adj,'ADJUDICATION',ifelse(sap,'SAP_PRIMARY_DEFINITION','OTHER_PARAMETER'));evidence[[i]]<-z}
  status<-if(any(seqmatch&p$raw_positive_supported))'EXPLICIT_SOURCE_LINK_POSITIVE' else if(any(dm&p$usable_positive))'POSITIVE_ON_ADJUDICATED_PCRDT' else if(any(sm&p$usable_positive))'POSITIVE_ON_SAP_PCRDT' else if(any(near&p$usable_positive))'NEARBY_POSITIVE_REQUIRES_REVIEW' else if(any(p$usable_positive))'POSITIVE_OTHER_DATE_ONLY' else if(any(p$positive_any))'POSITIVE_MENTION_WITH_UNRESOLVED_QUALIFIERS' else 'NO_POSITIVE_IN_RELEASED_MB_ADMB'
  bestdates<-if(any(seqmatch&p$usable_positive))p$collection_date[seqmatch&p$usable_positive] else if(any(dm&p$usable_positive))p$collection_date[dm&p$usable_positive] else if(any(sm&p$usable_positive))p$collection_date[sm&p$usable_positive] else as.Date(character())
  unique_best<-unique(bestdates[!is.na(bestdates)])
  summary[[i]]<-data.frame(case_id=c$case_id,adjudication_rows_same_event=sum(adj&same),sap_rows_same_event=sum(sap&same),
   adjudicated_PCRDT=ec_unique(as.character(direct)),sap_PCRDT_same_event=ec_unique(as.character(support)),
   positive_pcr_records=sum(p$usable_positive),raw_positive_pcr_records=sum(p$raw_positive_supported),
   nearby_positive_records=sum(near&p$usable_positive),candidate_records=sum(take),
   date_linked_positive_records=sum((seqmatch|dm|sm)&p$usable_positive),
   supporting_specimen_date=if(length(unique_best)==1)as.character(unique_best) else '',
   PCR_support_status=status,qualifying_specimen_verified=any(seqmatch&p$raw_positive_supported),stringsAsFactors=FALSE)
 }
 sm<-ec_bind(summary);j<-match(cases$case_id,sm$case_id)
 if(nrow(cases)){cases<-cbind(cases,sm[j,setdiff(names(sm),'case_id'),drop=FALSE])}
 ec_csv(ec_bind(edges),file.path(out,'case_pcr_candidate_links.csv'));ec_csv(ec_bind(evidence),file.path(out,'case_ADEFF_links.csv'))
 ec_csv(cases,file.path(out,'efficacy_cases.csv'));cases
}

# ---- Timelines: explicitly separate event dates from information availability --
ec_sources <- function()data.frame(
 id=c('FDA_EUA','FDA_BLA','BLA_ADRG','EMA','MODERNA_11NOV','MODERNA_16NOV','MODERNA_30NOV','ELECTION'),
 title=c('FDA EUA review memorandum, 18 December 2020, statistical plan / analyses',
 'FDA Summary Basis for Regulatory Action, 30 January 2022, pp.19-20',
 'P301 BLA Analysis Data Reviewer Guide, pp.23-25 and 37',
 'EMA public assessment report, pp.89 and 107-108',
 'Moderna case-accrual announcement, 11 November 2020',
 'Moderna first interim results, 16 November 2020 (SEC exhibit)',
 'Moderna final primary results, 30 November 2020 (SEC exhibit)',
 'FEC official 2020 presidential general-election results'),
 url=c('https://www.fda.gov/media/144673/download','https://www.fda.gov/media/155931/download',
 'https://phmpt.org/wp-content/uploads/2025/05/125752_S3_M5_mrna-1273-p301_A_D_adrg.pdf',
 'https://www.ema.europa.eu/en/documents/assessment-report/spikevax-previously-covid-19-vaccine-moderna-epar-public-assessment-report_en.pdf',
 'https://s29.q4cdn.com/435878511/files/doc_news/2020/11/11/moderna-has-completed-case-accrual-first-planned-interim.pdf',
 'https://www.sec.gov/Archives/edgar/data/1682852/000119312520293735/d53512dex991.htm',
 'https://www.sec.gov/Archives/edgar/data/1682852/000119312520304907/d76031dex991.htm',
 'https://www.fec.gov/documents/2840/2020presgeresults.pdf'),stringsAsFactors=FALSE)
ec_milestones <- function()data.frame(date=as.Date(c('2020-11-03','2020-11-07','2020-11-11','2020-11-16','2020-11-21','2020-11-30','2021-03-26')),
 label=c('US presidential election','Actual first interim efficacy cutoff','Case accrual publicly reported complete; preparation for DSMB ongoing',
 'First interim results announced: 5 vaccine / 90 placebo; 94.5% reported VE','Final EUA efficacy cutoff',
 'Final primary results announced: 11 vaccine / 185 placebo; 94.1% reported VE','BLA blinded-phase efficacy cutoff'),
 source_id=c('ELECTION','FDA_EUA','MODERNA_11NOV','MODERNA_16NOV','FDA_EUA','MODERNA_30NOV','FDA_BLA'),stringsAsFactors=FALSE)
ec_week <- function(d)as.Date(d)-((as.POSIXlt(as.Date(d))$wday+6L)%%7L)
ec_weekly <- function(cases,release){
 weeks<-seq(ec_week(as.Date('2020-07-27')),ec_week(as.Date(PHMPT_CUTOFFS[[release]])),by='7 days')
 out<-expand.grid(week_start=as.character(weeks),arm=c('mRNA-1273','Placebo'),stringsAsFactors=FALSE)
 out$week_start<-as.Date(out$week_start);out$week_end<-out$week_start+6;out$release<-release
 k<-ec_key(ec_week(cases$event_date),cases$arm);tab<-table(k);out$cases<-as.integer(tab[ec_key(out$week_start,out$arm)]);out$cases[is.na(out$cases)]<-0L
 out$partial_last_week<-out$week_end>as.Date(PHMPT_CUTOFFS[[release]]);out
}
ec_thresholds <- function(cases,release){
 out<-list();bases<-list(case_event_date=cases$event_date)
 # A date-linked specimen date is not a result-availability date. No near-date
 # candidate is selected merely to make a complete PCR-time curve.
 if('supporting_specimen_date'%in%names(cases))bases$date_linked_specimen_date<-ec_date(cases$supporting_specimen_date)
 if('analysis_ready_date'%in%names(cases))bases$documented_analysis_ready_date<-ec_date(cases$analysis_ready_date)
 for(b in names(bases)){
  dates<-bases[[b]];valid<-dates[!is.na(dates)];v<-sort(valid)
  for(n in c(53L,106L,151L)){
   d<-if(length(v)>=n)v[n] else as.Date(NA)
   complete<-length(valid)==nrow(cases)
   out[[length(out)+1L]]<-data.frame(release=release,basis=b,threshold=n,
    planned_analysis=if(n==53)'First interim' else if(n==106)'Second interim (not conducted separately)' else 'Final primary',
    threshold_date=as.character(d),dates_available=length(valid),total_cases=nrow(cases),date_coverage_complete=complete,
    n_before=if(is.na(d))NA_integer_ else sum(valid<d),n_on_date=if(is.na(d))NA_integer_ else sum(valid==d),
    n_by_date=if(is.na(d))NA_integer_ else sum(valid<=d),
    days_before_election=if(is.na(d))NA_integer_ else as.integer(as.Date('2020-11-03')-d),
    days_to_interim_announcement=if(is.na(d))NA_integer_ else as.integer(as.Date('2020-11-16')-d),
    interpretation=if(b=='case_event_date')'Retrospective event-date crossing in this snapshot; not an operational trigger date' else if(b=='date_linked_specimen_date')'Date-linked specimen subset; not test receipt, adjudication or readiness; incomplete coverage can shift the crossing' else if(complete)'All counted cases have user-documented readiness dates; still refers to this later snapshot membership' else 'Known-ready subset only; incomplete availability records cannot identify the actual trigger date',stringsAsFactors=FALSE)
  }
 }
 ec_bind(out)
}
ec_availability <- function(cases,root,out){
 template<-data.frame(release=cases$release,subject=cases$.sid,case_event_date=as.character(cases$event_date),
   pcr_result_received_date=rep('',nrow(cases)),adjudication_completed_date=rep('',nrow(cases)),analysis_ready_date=rep('',nrow(cases)),evidence_reference=rep('',nrow(cases)),stringsAsFactors=FALSE)
 ec_csv(template,file.path(out,'case_availability_template.csv'))
 cases$analysis_ready_date<-rep('',nrow(cases));cases$availability_evidence<-rep('',nrow(cases))
 if(!nzchar(PHMPT_AVAILABILITY_FILE))return(cases)
 a<-read.csv(ec_path(root,PHMPT_AVAILABILITY_FILE),colClasses='character',check.names=FALSE,fileEncoding='UTF-8-BOM',na.strings=character())
 ec_need(a,c('release','subject','case_event_date','analysis_ready_date','evidence_reference'),'availability CSV')
 a$subject<-ec_sid(a$subject);key<-ec_key(ec_norm(a$release),a$subject,a$case_event_date)
 if(anyDuplicated(key))stop('Duplicate keys in availability file; retain one reviewed record per release/subject/event date.')
 j<-match(ec_key(cases$release,cases$.sid,as.character(cases$event_date)),key);v<-a$analysis_ready_date[j];v[is.na(v)]<-'';ref<-a$evidence_reference[j];ref[is.na(ref)]<-''
 d<-ec_date(v);bad<-nzchar(v)&(is.na(d)|!nzchar(ref)|d<cases$event_date)
 if(any(bad))stop('Availability dates need a full valid date, evidence reference, and cannot precede the case event. Correct the reviewed availability file.')
 cases$analysis_ready_date<-ifelse(is.na(d),'',as.character(d));cases$availability_evidence<-ref
 ec_csv(a,file.path(out,'availability_input_used.csv'));cases
}

# ---- Print charts (one plot per PDF page / file) -------------------------------
ec_markers <- function(thresholds,include_press=TRUE){
 t<-thresholds[thresholds$basis=='case_event_date'&nzchar(thresholds$threshold_date)&!is.na(thresholds$threshold_date),,drop=FALSE]
 m<-data.frame(date=as.Date(t$threshold_date),code=paste0('T',t$threshold),label=paste0(t$threshold,' cases, retrospective'),color='#505050',lty=3L,stringsAsFactors=FALSE)
 m<-rbind(m,data.frame(date=as.Date('2020-11-03'),code='E',label='US election',color='#BB1E2D',lty=2L))
 if(include_press)m<-rbind(m,data.frame(date=as.Date(c('2020-11-16','2020-11-30')),code=c('P1','P2'),label=c('Interim announcement','Final primary announcement'),color='#74479B',lty=4L))
 m[order(m$date,m$code),,drop=FALSE]
}
ec_plot <- function(cases,weekly,thresholds,release,type='weekly',focus=FALSE,valid=TRUE){
 old<-par(no.readonly=TRUE);on.exit(par(old));par(mar=c(8.0,4.8,4.5,1.3),mgp=c(2.8,.8,0),tcl=-.25,las=1,family='sans',cex=.95)
 lo<-as.Date(if(focus)'2020-10-01' else '2020-07-27');hi<-as.Date(if(release=='EUA'||focus)'2020-12-03' else PHMPT_CUTOFFS[[release]])
 marks<-ec_markers(thresholds);marks<-marks[marks$date>=lo&marks$date<=hi,,drop=FALSE]
 colors<-c('mRNA-1273'='#2166A5','Placebo'='#D5842E');wk<-weekly[weekly$week_start+6>=lo&weekly$week_start<=hi,,drop=FALSE]
 main<-paste0(release,' - ',if(type=='cumulative')'cumulative primary efficacy cases' else 'weekly primary efficacy cases',if(focus)' | election-period detail' else '',if(!valid)' [PROVISIONAL: validation flags]' else '')
 if(type=='weekly'){
  ymax<-max(c(1,wk$cases))*1.25
  plot(NA,xlim=as.numeric(c(lo,hi)),ylim=c(0,ymax),xaxt='n',xlab='',ylab='Counted cases',main=main,bty='l')
  abline(h=pretty(c(0,ymax)),col='#E6E6E6',lwd=.7)
  for(a in names(colors)){
   z<-wk[wk$arm==a,,drop=FALSE];x<-as.numeric(z$week_start);left<-if(a=='mRNA-1273')x+.35 else x+3.6;right<-left+2.8
   rect(left,0,right,z$cases,col=colors[a],border=NA)
   nz<-z$cases>0;text((left[nz]+right[nz])/2,z$cases[nz]+ymax*.022,labels=z$cases[nz],cex=if(focus).85 else .7,col=colors[a])
  }
 }else{
  ymax<-max(1,nrow(cases))*1.13
  plot(NA,xlim=as.numeric(c(lo,hi)),ylim=c(0,ymax),xaxt='n',xlab='',ylab='Cumulative counted cases',main=main,bty='l')
  abline(h=pretty(c(0,ymax)),col='#E6E6E6',lwd=.7)
  for(a in c(names(colors),'Total')){
   d<-if(a=='Total')cases$event_date else cases$event_date[cases$arm==a];u<-sort(unique(d));ct<-vapply(u,function(v)sum(d<=v),integer(1))
   lines(c(as.numeric(lo),as.numeric(u),as.numeric(hi)),c(0,ct,length(d)),type='s',col=if(a=='Total')'#252525' else colors[a],lwd=if(a=='Total')2.3 else 1.8)
  }
  for(n in c(53,106,151))if(n<ymax)abline(h=n,col='#AAAAAA',lty=3)
 }
 ticks<-seq(ec_week(lo),hi,by=if(as.integer(hi-lo)>150)'28 days' else if(focus)'7 days' else '14 days');ticks<-ticks[ticks>=lo]
 axis(1,at=as.numeric(ticks),labels=format(ticks,'%d %b'),cex.axis=.84)
 if(nrow(marks))for(i in seq_len(nrow(marks))){abline(v=as.numeric(marks$date[i]),col=marks$color[i],lty=marks$lty[i],lwd=1.5);text(as.numeric(marks$date[i]),ymax*(.98-.055*((i-1)%%3)),marks$code[i],col=marks$color[i],cex=.83,font=2,pos=if(i%%2)2 else 4,offset=.25)}
 mtext('Adjudicated per-protocol cases; retrospective event dates, not dates available to the DSMB',3,line=.55,cex=.86)
 labs<-c('mRNA-1273','Placebo');cols<-colors
 if(type=='cumulative'){labs<-c(labs,'Total');cols<-c(cols,'#252525')}
 if(type=='weekly')legend('topleft',legend=labs,fill=cols,bty='n',cex=.9,inset=.01)
 else legend('topleft',legend=labs,lty=1,col=cols,lwd=2,bty='n',cex=.9,inset=.01)
 line1<-paste(paste0(marks$code,' = ',format(marks$date,'%d %b %Y')),collapse='; ')
 captions<-c(if(type=='weekly')'Weeks run Monday-Sunday; incomplete cutoff week is retained. Numbers label nonzero bars.' else 'One first qualifying event per counted participant; ties share a calendar date.',
  line1,'T53 / T106 = planned interim thresholds; T151 = planned final threshold. E = election; P1 / P2 = announcements.',
  'Threshold crossings use the later released snapshot. The second interim was not separately conducted.',
  'A pre-election event-date crossing does not establish pre-election analysis readiness or a politically motivated delay.')
 lines<-unlist(lapply(captions,function(s)strwrap(s,width=145)))
 for(i in seq_along(lines))mtext(lines[i],side=1,line=2.2+(i-1)*.8,adj=0,cex=.72,col='#444444')
 invisible(NULL)
}
ec_save_plot <- function(cases,weekly,th,release,out,name,type='weekly',focus=FALSE,valid=TRUE){
 draw<-function()ec_plot(cases,weekly,th,release,type,focus,valid)
 if(capabilities('cairo'))grDevices::cairo_pdf(file.path(out,paste0(name,'.pdf')),width=11.7,height=8.3,family='sans',onefile=TRUE)
 else grDevices::pdf(file.path(out,paste0(name,'.pdf')),width=11.7,height=8.3,useDingbats=FALSE,onefile=TRUE)
 tryCatch(draw(),finally=grDevices::dev.off())
 grDevices::png(file.path(out,paste0(name,'.png')),width=11.7,height=8.3,units='in',res=PHMPT_PNG_DPI)
 tryCatch(draw(),finally=grDevices::dev.off())
}

# ---- HTML results report ------------------------------------------------------
ec_html <- function(x){x<-ec_text(x);x<-gsub('&','&amp;',x,fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE);x<-gsub('>','&gt;',x,fixed=TRUE);gsub('"','&quot;',x,fixed=TRUE)}
ec_table <- function(x,limit=40L){
 if(is.null(x)||!nrow(x))return('<p>No records.</p>')
 total<-nrow(x);x<-head(x,limit)
 header<-paste0('<tr>',paste0('<th>',ec_html(names(x)),'</th>',collapse=''),'</tr>')
 rows<-vapply(seq_len(nrow(x)),function(i){values<-vapply(x,function(v)ec_text(v[i]),'');paste0('<tr>',paste0('<td>',ec_html(values),'</td>',collapse=''),'</tr>')},'')
 paste0('<div class="table-wrap"><table><thead>',header,'</thead><tbody>',paste(rows,collapse=''),'</tbody></table></div>',if(total>limit)paste0('<p>Showing ',limit,' of ',total,' rows. Full detail is in the CSV exports.</p>') else '')
}
ec_report <- function(results,status,out,root){
 s<-ec_sources();ec_csv(s,file.path(out,'sources.csv'));milestones<-ec_milestones();ec_csv(milestones,file.path(out,'documented_milestones.csv'))
 chunks<-c('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Moderna efficacy cases and analysis timeline</title>',
 '<style>body{max-width:1120px;margin:36px auto;padding:0 26px;font:16px/1.55 Arial,sans-serif;color:#20252b}h1{font-size:29px}h2{margin-top:36px;font-size:23px}h3{font-size:19px}.note{border-left:4px solid #2166a5;padding:12px 18px;background:#f1f6fa}.warn{border-left:4px solid #b53432;padding:12px 18px;background:#fff4ee}.table-wrap{overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px;margin:12px 0}th,td{text-align:left;padding:8px;border-bottom:1px solid #d7dde2;vertical-align:top}th{background:#eaf0f5}img{width:100%;height:auto}code{background:#f0f2f3;padding:2px 4px}a{color:#165c91}small{color:#555}@media print{body{max-width:none;margin:0;padding:0;font-size:11pt}h2,h3{break-after:avoid}img{break-inside:avoid}thead{display:table-header-group}.table-wrap{overflow:visible}a{color:inherit}table{font-size:9pt}}@page{size:A4;margin:16mm}</style><body>',
 '<h1>Moderna P301: efficacy cases and interim-analysis timeline</h1>',
 paste0('<p>Generated ',ec_html(format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z')),'; script v1.0.0.</p>'),
 '<div class="note"><strong>Question assessed:</strong> when did the eventually counted cases occur relative to the planned analysis thresholds, the US election and Moderna\'s announcements? Event dates support a retrospective chronology. They do not identify when results arrived, cases were adjudicated, a database was ready, or management learned an efficacy result.</div>',
 '<h2>Processing status</h2>',ec_table(status),
 '<h2>Documented plan and announcements</h2>',
 '<p>Two interim analyses were planned at 53 and 106 cases, followed by a final analysis at 151. The actual first interim analysis included 95 cases, and the final primary analysis included 196. The second interim was not conducted separately [FDA_EUA; EMA].</p>',
 '<p>EMA questioned whether the pandemic acceleration and FDA safety-follow-up explanation adequately justified the changed timing, while considering the impact on its assessment of vaccination benefit low [EMA, pp.89 and 107-108]. This is a documented methodological concern, not a finding of election-related intent.</p>',
 ec_table(milestones),
 '<p>On 11 November Moderna said it had completed accrual and was preparing the data for the independent DSMB, while remaining blinded [MODERNA_11NOV]. On 16 November it announced the interim result, followed by the final primary result on 30 November. An EUA safety-follow-up requirement and readiness to announce an interim efficacy result are separate questions.</p>')
 for(rel in names(results)){
  r<-results[[rel]];c<-r$cases;t<-r$thresholds;et<-t[t$basis=='case_event_date',,drop=FALSE]
  chunks<-c(chunks,paste0('<h2>',rel,': primary efficacy reconstruction</h2>'),ec_table(r$checks),
   if(!r$valid)'<div class="warn">Published counts or membership checks do not reconcile. The following figures are provisional and must not be described as a validated reconstruction.</div>' else '<p>The arm-specific primary population and case totals match the published benchmarks. This validates the selected aggregate case set; it does not independently verify every case definition or specimen linkage.</p>')
  npre<-sum(c$event_date<as.Date('2020-11-03'));nelec<-sum(c$event_date==as.Date('2020-11-03'));n7<-sum(c$event_date<=as.Date('2020-11-07'))
  first<-et[et$threshold==53,,drop=FALSE];d<-if(nrow(first))first$threshold_date[1] else NA_character_
  para<-paste0('In the ',rel,' snapshot, ',npre,' counted cases have an event date before the election and ',nelec,' fall on election day. ',
   if(!is.na(d)&&nzchar(d))paste0('The 53-case event-date threshold is crossed on ',d,', ',as.integer(as.Date('2020-11-16')-as.Date(d)),' calendar days before the interim announcement. ') else 'The 53-case threshold is not reached. ',
   'These are retrospective case-event intervals, not measured operational delays. The later snapshot contains ',n7,' selected cases dated through 7 November; this is not assumed to reproduce the original 95-case interim database.')
  chunks<-c(chunks,paste0('<p>',ec_html(para),'</p>'),ec_table(et[c('threshold','planned_analysis','threshold_date','n_before','n_on_date','n_by_date','days_before_election','days_to_interim_announcement')]),
   paste0('<p><a href="',rel,'/efficacy_cases.csv">Case-level CSV</a> | <a href="',rel,'/thresholds.csv">All threshold calculations</a> | <a href="',rel,'/weekly_cases.csv">Weekly counts</a></p>'),
   paste0('<img src="',rel,'/weekly_cases.png" alt="Weekly counted cases by treatment arm"><p><a href="',rel,'/weekly_cases.pdf">Print PDF</a></p>'))
  if(rel=='EUA')chunks<-c(chunks,paste0('<img src="',rel,'/weekly_cases_election_focus.png" alt="Election period detail"><p><a href="',rel,'/weekly_cases_election_focus.pdf">Print detail PDF</a></p><img src="',rel,'/cumulative_cases.png" alt="Cumulative counted cases"><p><a href="',rel,'/cumulative_cases.pdf">Print cumulative PDF</a></p>'))
  if(rel=='BLA')chunks<-c(chunks,'<p>The BLA plot reflects the updated case definitions, per-protocol population and follow-up through March 2021. Its historical crossings are a sensitivity comparison, not evidence of what was known in October-November 2020.</p>')
  chunks<-c(chunks,'<h3>PCR evidence and linkage</h3>')
  summary<-if(nrow(c))as.data.frame(table(arm=c$arm,PCR_support_status=c$PCR_support_status),stringsAsFactors=FALSE) else data.frame();if(nrow(summary))summary<-summary[summary$Freq>0,,drop=FALSE]
  chunks<-c(chunks,ec_table(summary),paste0('<p>',sum(c$adjudication_rows_same_event==0),' cases have no same-date adjudicated ADEFF row; ',sum(!nzchar(c$adjudicated_PCRDT)),' have no populated positive-PCR date on the matched adjudicated row. ',sum(c$qualifying_specimen_verified),' have an explicit MB source-sequence link with a compatible raw positive result.</p>'),
   '<p>ADJUDICATED_PCRDT and SAP_PCRDT links identify candidate specimen dates. A matching date can contain several tests. NEARBY_POSITIVE means a candidate within the configurable retrieval window, not a verified qualifying test. Blank PCR dates, unmatched records and ambiguous or conflicting results remain visible; none remove a case from the official numerator.</p>',
   paste0('<p><a href="',rel,'/case_pcr_candidate_links.csv">Case-to-test candidate links</a> | <a href="',rel,'/case_subject_pcr_records.csv">All released SARS-CoV-2 records for case subjects</a> | <a href="',rel,'/case_ADEFF_links.csv">ADEFF supporting records</a></p>'))
  cov<-sum(nzchar(c$analysis_ready_date));chunks<-c(chunks,paste0('<h3>Information availability</h3><p>Evidence-backed analysis-ready dates supplied by the analyst: ',cov,'/',nrow(c),'. ',if(cov<nrow(c))'Incomplete coverage prevents identifying the actual operational trigger date from these records.' else 'Availability dates cover the selected later snapshot; confirm that its membership matches the contemporaneous review set.','</p>'),
   '<p>CVD53DT is exported separately as a sponsor-derived 53rd-case date. Its presence does not turn an event date into a database-lock or adjudication-completion date. No CE assessment date is automatically treated as an adjudication-completion timestamp.</p>')
 }
 if(all(c('EUA','BLA')%in%names(results)))chunks<-c(chunks,'<h2>Changes in case membership between submissions</h2><p><a href="eua_bla_case_comparison.csv">Subject-level EUA/BLA case comparison</a> reports retained cases, changed event dates and cases absent from the other selected primary set. Differences require review of population, adjudication and censoring changes; they do not by themselves establish deleted or altered source tests.</p>')
 chunks<-c(chunks,'<h2>What is needed to assess an avoidable announcement delay?</h2><p>The testable questions are the time from specimen collection to a reported result, from result receipt to case adjudication, from adjudication to the prespecified analysis, and from the DSMB recommendation to public disclosure. The necessary evidence includes timestamped laboratory transmissions, case-submission and adjudication logs, database freeze/lock records, DSMB schedules and minutes, and dated communications about analysis timing. The event-date charts identify intervals to investigate; they cannot establish political motivation.</p>',
 '<h2>Methods and audit trail</h2><p>Primary cases are selected from adjudicated ADTT (EUA) or ADTTEA (BLA), using the exact endpoint title, per-protocol flag and uncensored-event code. BLA is restricted to period 1 where APERIOD exists. ADSL is not required because these subject-level variables are already carried into the time-to-event datasets. Planned treatment TRT01P identifies the randomized arm. No AE or symptom-only record is added to the official case numerator.</p><p>ADMB joins MB only on subject and MBSEQ within the same submission; duplicate keys are not collapsed. A derived ADMB row with DTYPE populated is exported separately. All nearby positives, negatives and conflicts are retained. Partial dates are not imputed, PCR result dates are not substituted for collection dates, and clinical timestamps are not converted to information-availability timestamps.</p>',
 '<p>Thresholds count unique primary-endpoint subjects across both arms, use exact calendar dates, and report same-day ties. Weekly bins begin Monday. The primary timing analysis uses EUA; BLA is kept separate. No inferential test of an election effect is produced because one election and an unknown operational timeline do not support such a test.</p>',
 '<p><a href="input_manifest.csv">Input paths / checksums</a> | <a href="parameters.csv">Run settings</a> | <a href="document_inventory.csv">Candidate local documentation</a> | <a href="sources.csv">Sources</a></p><h2>Sources</h2>')
 for(i in seq_len(nrow(s)))chunks<-c(chunks,paste0('<p>[',s$id[i],'] <a href="',ec_html(s$url[i]),'">',ec_html(s$title[i]),'</a></p>'))
 chunks<-c(chunks,'</body></html>');writeLines(chunks,file.path(out,'efficacy_cases_report.html'),useBytes=TRUE)
}

# ---- Main standalone execution ------------------------------------------------
analyse_efficacy_cases <- function(root=PHMPT_ROOT,output_dir=PHMPT_OUTPUT_DIR){
 root<-ec_root(root)
 if(!length(PHMPT_RELEASES)||any(!PHMPT_RELEASES%in%c('EUA','BLA')))stop('PHMPT_RELEASES must contain EUA and/or BLA.')
 if(!is.numeric(PHMPT_REVIEW_WINDOW_DAYS)||length(PHMPT_REVIEW_WINDOW_DAYS)!=1||is.na(PHMPT_REVIEW_WINDOW_DAYS)||PHMPT_REVIEW_WINDOW_DAYS<0)stop('Invalid candidate retrieval window.')
 base<-if(nzchar(output_dir))ec_path(root,output_dir) else file.path(root,'analysis')
 dir.create(base,recursive=TRUE,showWarnings=FALSE)
 out<-tempfile(paste0('efficacy_cases_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_'),tmpdir=base);dir.create(out)
 pars<-data.frame(setting=c('version','root','EUA_submission','BLA_submission','EUA_efficacy_cutoff','BLA_efficacy_cutoff','PCR_review_window_days','CSV_encoding','CSV_date_origin','CSV_date_format','availability_file'),
  value=c('1.0.0',root,PHMPT_SUBMISSIONS,PHMPT_CUTOFFS,PHMPT_REVIEW_WINDOW_DAYS,PHMPT_CSV_ENCODING,PHMPT_CSV_DATE_ORIGIN,PHMPT_CSV_DATE_FORMAT,PHMPT_AVAILABILITY_FILE))
 ec_csv(pars,file.path(out,'parameters.csv'));writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
 manifests<-list();results<-list();status<-list()
 for(rel in PHMPT_RELEASES){
  dest<-file.path(out,rel);dir.create(dest);message('Processing ',rel)
  ans<-tryCatch({
   t<-ec_pick(root,rel,'TTE',dest);manifests[[length(manifests)+1L]]<-t$manifest
   selected<-ec_cases(t$data,rel,dest);rm(t);gc(verbose=FALSE)
   subs<-selected$cases$.sid
   inputs<-list()
   for(dom in c('ADEFF','MB','ADMB','CE','SUPPCE')){
    z<-ec_pick(root,rel,dom,dest,subjects=subs,required=dom%in%c('ADEFF','MB'))
    if(!is.null(z)){inputs[[dom]]<-z$data;manifests[[length(manifests)+1L]]<-z$manifest}
   }
   if(!is.null(inputs$CE)){
    ce<-inputs$CE;adj<-if('CECAT'%in%names(ce))grepl('ADJUDIC',ec_norm(ce$CECAT)) else rep(TRUE,nrow(ce))
    ec_csv(ce[adj,,drop=FALSE],file.path(dest,'case_subject_CE_adjudication.csv'))
   }
   if(!is.null(inputs$SUPPCE))ec_csv(inputs$SUPPCE,file.path(dest,'case_subject_SUPPCE.csv'))
   tests<-ec_tests(inputs$MB,inputs$ADMB,dest);cases<-ec_link(selected$cases,inputs$ADEFF,tests,dest)
   cases<-ec_availability(cases,root,dest);ec_csv(cases,file.path(dest,'efficacy_cases.csv'))
   weekly<-ec_weekly(cases,rel);th<-ec_thresholds(cases,rel)
   ec_csv(weekly,file.path(dest,'weekly_cases.csv'));ec_csv(th,file.path(dest,'thresholds.csv'))
   daily_dates<-seq(as.Date('2020-07-27'),as.Date(PHMPT_CUTOFFS[[rel]]),by='day')
   daily<-data.frame(date=daily_dates,mRNA_1273=vapply(daily_dates,function(d)sum(cases$arm=='mRNA-1273'&cases$event_date<=d),integer(1)),
    Placebo=vapply(daily_dates,function(d)sum(cases$arm=='Placebo'&cases$event_date<=d),integer(1)))
   daily$total<-daily$mRNA_1273+daily$Placebo;ec_csv(daily,file.path(dest,'cumulative_cases.csv'))
   ec_save_plot(cases,weekly,th,rel,dest,'weekly_cases',valid=selected$valid)
   if(rel=='EUA'){
    ec_save_plot(cases,weekly,th,rel,dest,'weekly_cases_election_focus',focus=TRUE,valid=selected$valid)
    ec_save_plot(cases,weekly,th,rel,dest,'cumulative_cases',type='cumulative',valid=selected$valid)
   }
   list(cases=cases,weekly=weekly,thresholds=th,checks=selected$checks,source53=selected$source53,valid=selected$valid)
  },error=function(e)e)
  if(inherits(ans,'error')){status[[rel]]<-data.frame(release=rel,status='FAILED',detail=conditionMessage(ans));message(rel,': ',conditionMessage(ans))}
  else{results[[rel]]<-ans;status[[rel]]<-data.frame(release=rel,status=if(ans$valid)'COUNTS_RECONCILED' else 'PROVISIONAL_COUNT_OR_RULE_MISMATCH',detail=paste(nrow(ans$cases),'cases; PCR linkage is audited separately'))}
  ec_csv(ec_bind(manifests),file.path(out,'input_manifest.csv'));ec_csv(ec_bind(status),file.path(out,'release_status.csv'))
 }
 if(all(c('EUA','BLA')%in%names(results))){
  a<-results$EUA$cases;b<-results$BLA$cases;sid<-union(a$.sid,b$.sid);i<-match(sid,a$.sid);j<-match(sid,b$.sid)
  cmp<-data.frame(subject=sid,EUA_case=!is.na(i),BLA_case=!is.na(j),EUA_arm=a$arm[i],BLA_arm=b$arm[j],EUA_event_date=a$event_date[i],BLA_event_date=b$event_date[j])
  cmp$comparison<-ifelse(!cmp$EUA_case,'BLA_PRIMARY_CASE_ONLY',ifelse(!cmp$BLA_case,'EUA_PRIMARY_CASE_NOT_IN_BLA_PRIMARY_SET',ifelse(cmp$EUA_event_date==cmp$BLA_event_date,'RETAINED_SAME_EVENT_DATE','RETAINED_EVENT_DATE_CHANGED')))
  ec_csv(cmp,file.path(out,'eua_bla_case_comparison.csv'))
 }
 docs<-unique(unlist(lapply(file.path(root,c('eua_data','prod_data')),function(p)list.files(p,recursive=TRUE,full.names=TRUE,pattern='\\.(pdf|xml|txt|sas)$',ignore.case=TRUE))))
 docs<-docs[grepl('adrg|csdrg|define|sap|statistical|protocol|adjudic|dsmb|a40adeff|a50adtt|t14020201030101',basename(docs),ignore.case=TRUE)]
 ec_csv(data.frame(path=docs,filename=basename(docs)),file.path(out,'document_inventory.csv'))
 ec_report(results,ec_bind(status),out,root)
 message('Report: ',normalizePath(file.path(out,'efficacy_cases_report.html'),winslash='/'))
 message('CSV assets and one-chart-per-page PDFs: ',normalizePath(out,winslash='/'))
 if(PHMPT_SHOW_REPORT)utils::browseURL(normalizePath(file.path(out,'efficacy_cases_report.html'),winslash='/'))
 if(any(vapply(status,function(x)x$status=='FAILED',logical(1))))warning('One or more releases failed. The HTML and discovery logs identify the missing or ambiguous inputs.',call.=FALSE)
 invisible(list(output_dir=out,results=results,status=ec_bind(status)))
}

# ---- RUN: executes when the whole standalone script is run ---------------------
if(isTRUE(getOption('phmpt.autorun',TRUE)))analyse_efficacy_cases()
