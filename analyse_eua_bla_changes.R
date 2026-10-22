# =============================================================================
# Moderna P301 | EUA / BLA-targeted production snapshot reconciliation | v1.0.0
# 2026-10-22 | Standalone R >= 4.3: select all and run in RStudio.
# Reads original local XPTs; no trial files are downloaded or modified.
# XPT/SAS7BDAT requires haven. No packages are installed automatically.
# This is a comparison of released snapshots, NOT an electronic audit trail.
# Changed / unmatched records do not establish deletion, misconduct or cause.
# Load only: options(phmpt.autorun=FALSE) before source().
# =============================================================================
PHMPT_ROOT <- getOption('phmpt.root','')
PHMPT_OUTPUT_DIR <- getOption('phmpt.output_dir','')
PHMPT_EUA_CUTOFF <- '2020-11-25'
PHMPT_BLA_CUTOFF <- '2021-03-26'
PHMPT_SHOW_REPORT <- interactive()
# Explicit overrides can be absolute or project-relative. BLA scans prod_data.
PHMPT_FILES <- list(EUA=list(ADSL='',AE='',ADSYMP='',MB='',LB='',ADMB='',ADAE=''),
                    BLA=list(ADSL='',AE='',ADSYMP='',MB='',LB='',ADMB='',ADAE=''))
# These identify the submissions used in the preceding P301 analyses, not URLs.
# If a selected submission is absent, discovery logs alternatives and stops for
# that domain; files from distinct submissions are never silently concatenated.
PHMPT_SUBMISSIONS <- c(EUA='27073_S5',BLA='125752_S3')
PHMPT_DOMAINS <- c('AE','ADSYMP','MB','LB','ADMB','ADAE')
PHMPT_REQUIRED_DOMAINS <- c('AE','ADSYMP','MB')
PHMPT_ALLOW_CSV_FALLBACK <- FALSE # original XPTs preferred; applies to compared domains
PHMPT_CSV_DATE_ORIGIN <- ''       # numeric CSV dates need explicit provenance
PHMPT_CANDIDATE_DATE_TOLERANCE <- 7L # review-only candidates; never automatic matches
PHMPT_MAX_CANDIDATES_PER_RECORD <- 25L # every count retained; detail list may be capped
PHMPT_HTML_EXAMPLES <- 30L

sc_spec <- list(
 AE=list(seq='AESEQ',date='AESTDTC',entity='AETERM',refs=c('AESPID','AEREFID'),
         clinical=c('AESTDTC','AETERM'),result=character(),serious='AESER',kind='Source AE records'),
 ADSYMP=list(seq='ASEQ',date=c('ADT','ASTDT'),entity='PARAMCD',refs='LNKID',
             clinical=c('ADT','PARAMCD'),result='AVALC',serious=character(),kind='Derived symptom assessments (ADSYMP / EUA ADSY alias)'),
 MB=list(seq='MBSEQ',date='MBDTC',entity='MBTESTCD',refs=c('MBREFID','MBSPID'),
         clinical=c('MBDTC','MBTESTCD','MBSPEC','MBMETHOD'),result=c('MBSTRESC','MBORRES'),serious=character(),kind='Source microbiology: PCR and other tests'),
 LB=list(seq='LBSEQ',date='LBDTC',entity='LBTESTCD',refs=c('LBREFID','LBSPID'),
         clinical=c('LBDTC','LBTESTCD','LBSPEC','LBMETHOD'),result=c('LBSTRESC','LBORRES'),serious=character(),kind='Source laboratory results'),
 ADMB=list(seq=c('MBSEQ','ASEQ'),date=c('ADT','MBDTC'),entity='PARAMCD',refs=character(),
           clinical=c('ADT','PARAMCD','DTYPE','AVISIT'),result='AVALC',serious=character(),kind='Derived microbiology analysis records'),
 ADAE=list(seq='AESEQ',date=c('ASTDT','AESTDTC'),entity=c('AETERM','AEDECOD'),refs=c('AESPID','AEREFID'),
           clinical=c('ASTDT','AETERM'),result=character(),serious='AESER',kind='Derived AE analysis records')
)
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
sc_root <- function(root=''){
 if(nzchar(root))return(normalizePath(root,winslash='/',mustWork=TRUE))
 cand<-getwd();for(f in sys.frames())if(!is.null(f$ofile))cand<-c(cand,dirname(f$ofile))
 a<-grep('^--file=',commandArgs(),value=TRUE);if(length(a))cand<-c(cand,dirname(sub('^--file=','',a[1])))
 if(interactive()&&requireNamespace('rstudioapi',quietly=TRUE)&&rstudioapi::isAvailable())cand<-c(cand,dirname(rstudioapi::getActiveDocumentContext()$path))
 for(p in unique(cand[dir.exists(cand)]))for(i in 1:8){if(all(dir.exists(file.path(p,c('eua_data','prod_data')))))return(normalizePath(p,winslash='/'));q<-dirname(p);if(q==p)break;p<-q}
 stop('Cannot locate eua_data and prod_data. Set PHMPT_ROOT.')
}
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
  if(!requireNamespace('haven',quietly=TRUE))stop('Package haven is required to read XPT/SAS7BDAT. Install it once with install.packages("haven").')
  x<-if(ext=='xpt')haven::read_xpt(path) else haven::read_sas(path)
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
sc_candidates <- function(root,release,domain){
 folder<-file.path(root,if(release=='EUA')'eua_data' else 'prod_data')
 ext<-if(domain=='ADSL')'xpt|sas7bdat|csv|jmp' else if(isTRUE(PHMPT_ALLOW_CSV_FALLBACK))'xpt|sas7bdat|csv' else 'xpt|sas7bdat'
 p<-list.files(folder,recursive=TRUE,full.names=TRUE,pattern=paste0('\\.(',ext,')$'),ignore.case=TRUE)
 aliases<-if(domain=='ADSYMP')'adsymp|adsy' else tolower(domain)
 p<-p[grepl(paste0('(^|[_ -])(',aliases,')([_. -]|$)'),basename(p),ignore.case=TRUE)&!grepl('p101|p201|dmid',basename(p),ignore.case=TRUE)]
 p[order(match(tolower(tools::file_ext(p)),c('xpt','sas7bdat','csv','jmp')),p)]
}
sc_equivalent <- function(a,b){
 na<-sort(setdiff(names(a),grep('^\\.',names(a),value=TRUE)));nb<-sort(setdiff(names(b),grep('^\\.',names(b),value=TRUE)))
 identical(na,nb)&&nrow(a)==nrow(b)&&identical(sort(sc_key(a,na,require_all=FALSE)),sort(sc_key(b,nb,require_all=FALSE)))
}
sc_pick <- function(root,release,domain,out){
 explicit<-PHMPT_FILES[[release]][[domain]];if(is.null(explicit))explicit<-''
 paths<-if(nzchar(explicit))sc_path(root,explicit) else sc_candidates(root,release,domain)
 log<-data.frame(release=rep(release,length(paths)),domain=rep(domain,length(paths)),path=paths,status=rep('CANDIDATE',length(paths)),detail=rep('',length(paths)),stringsAsFactors=FALSE)
 if(!length(paths)){sc_csv(log,file.path(out,paste0('discovery_',release,'_',domain,'.csv')));return(NULL)}
 if(!all(file.exists(paths)))stop('Explicit file does not exist: ',paste(paths[!file.exists(paths)],collapse=', '))
 choose<-rep(TRUE,length(paths));sub<-PHMPT_SUBMISSIONS[[release]]
 if(!nzchar(explicit)&&!is.null(sub)&&nzchar(sub)){
  choose<-grepl(paste0(sub,'_'),basename(paths),ignore.case=TRUE)
  log$status[!choose]<-'OTHER_SUBMISSION';log$detail[!choose]<-'Not combined with selected submission'
  if(!any(choose)){sc_csv(log,file.path(out,paste0('discovery_',release,'_',domain,'.csv')));stop('No ',release,' ',domain,' candidate in ',sub,'. Review discovery log and PHMPT_SUBMISSIONS / PHMPT_FILES.')}
 }
 # A CSV mirror does not compete with its available original XPT.
 if(!nzchar(explicit)&&domain!='ADSL'&&any(choose&tolower(tools::file_ext(paths))=='xpt')){i<-choose&tolower(tools::file_ext(paths))!='xpt';choose[i]<-FALSE;log$status[i]<-'PREFER_ORIGINAL_XPT'}
 selected<-NULL;errors<-character()
 for(i in which(choose)){
  message('Reading ',release,' ',domain,': ',basename(paths[i]))
  a<-tryCatch(sc_read(paths[i],domain),error=function(e)e)
  if(inherits(a,'error')){log$status[i]<-'UNREADABLE';log$detail[i]<-conditionMessage(a);errors<-c(errors,paths[i]);next}
  if(is.null(selected)){selected<-a;log$status[i]<-'SELECTED'}
  else if(identical(selected$md5,a$md5)||sc_equivalent(selected$data,a$data)){log$status[i]<-'EQUIVALENT_VALUES';log$detail[i]<-'Same cell values and duplicate multiplicity; row order and import labels ignored'}
  else {log$status[i]<-'CONFLICTING_SNAPSHOT';sc_csv(log,file.path(out,paste0('discovery_',release,'_',domain,'.csv')));stop('Multiple non-equivalent ',release,' ',domain,' tables. Select PHMPT_FILES$',release,'$',domain,' explicitly; no tables were concatenated.')}
 }
 sc_csv(log,file.path(out,paste0('discovery_',release,'_',domain,'.csv')))
 if(is.null(selected))stop('No readable ',release,' ',domain,'; see discovery log')
 if(length(errors)&&!nzchar(explicit))stop('Unreadable competing candidates remain for ',release,' ',domain,'. Select a reviewed explicit path.')
 selected$release<-release;selected$domain<-domain;selected
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
  for(v in head(j,PHMPT_MAX_CANDIDATES_PER_RECORD))rows[[length(rows)+1L]]<-data.frame(eua_index=i,bla_index=v,subject=a$.sid[i],eua_row=a$.row[i],bla_row=b$.row[v],reason=paste(names(hits)[vapply(hits,function(z)v%in%z,logical(1))],collapse=';'),bla_already_matched=matched_b[v],
    eua_date=a$.date_raw[i],bla_date=b$.date_raw[v],eua_entity=a$.entity[i],bla_entity=b$.entity[v],eua_result=result_a[i],bla_result=result_b[v],eua_result_class=class_a[i],bla_result_class=class_b[v],stringsAsFactors=FALSE)
 }
 list(rows=if(length(rows))do.call(rbind,rows) else data.frame(),counts=counts,truncated=trunc)
}
sc_rosters <- function(a,b,out){
 one<-function(x,rel){if(is.null(x))return(data.frame(subject=character(),arm=character(),randomized=character(),cutoff=character()));z<-x$data
  r<-data.frame(subject=z$.sid,arm=sc_first(z,c('TRT01P','ARM')),randomized=sc_col(z,'RANDFL'),cutoff=sc_col(z,'CUTOFFDT'),stringsAsFactors=FALSE)
  if(anyDuplicated(r$subject))stop('Duplicate ADSL subject in ',rel);r}
 x<-one(a,'EUA');y<-one(b,'BLA');ids<-union(x$subject,y$subject);i<-match(ids,x$subject);j<-match(ids,y$subject)
 z<-data.frame(subject=ids,eua_present=!is.na(i),bla_present=!is.na(j),eua_arm=x$arm[i],bla_arm=y$arm[j],eua_randomized=x$randomized[i],bla_randomized=y$randomized[j],eua_cutoff=x$cutoff[i],bla_cutoff=y$cutoff[j],stringsAsFactors=FALSE)
 z$arm<-ifelse(!is.na(z$eua_arm)&nzchar(z$eua_arm),z$eua_arm,ifelse(!is.na(z$bla_arm)&nzchar(z$bla_arm),z$bla_arm,'UNMAPPED'))
 z$arm_source<-ifelse(!is.na(z$eua_arm)&nzchar(z$eua_arm),'EUA','BLA_ONLY_OR_UNMAPPED')
 z$arm_conflict<-!is.na(z$eua_arm)&!is.na(z$bla_arm)&sc_norm(z$eua_arm)!=sc_norm(z$bla_arm)
 z$arm[z$arm_conflict]<-'ARM_CONFLICT';sc_csv(z,file.path(out,'subject_roster_comparison.csv'));z
}
sc_schema <- function(a,b,domain){
 nm<-union(a$metadata$variable,b$metadata$variable);i<-match(nm,a$metadata$variable);j<-match(nm,b$metadata$variable)
 data.frame(domain=domain,variable=nm,in_eua=!is.na(i),in_bla=!is.na(j),eua_type=a$metadata$type[i],bla_type=b$metadata$type[j],eua_label=a$metadata$label[i],bla_label=b$metadata$label[j],stringsAsFactors=FALSE)
}
sc_compare <- function(old,new,domain,roster,out){
 a<-old$data;b<-new$data;sp<-sc_spec[[domain]];folder<-file.path(out,domain);dir.create(folder,showWarnings=FALSE)
 schema<-sc_schema(old,new,domain);sc_csv(schema,file.path(folder,'schema_comparison.csv'))
 m<-sc_pairs(a,b,domain);q<-sc_review_candidates(a,b,m,domain);pairs<-m$pairs
 arm<-function(s){z<-roster$arm[match(s,roster$subject)];z[is.na(z)]<-'UNMAPPED';z}
 common<-setdiff(intersect(names(a),names(b)),grep('^\\.',names(a),value=TRUE))
 administrative<-c('CUTOFFDT','DCO','DCODT','DCOCUTOFF')
 ids<-unique(c('STUDYID','USUBJID','SUBJID','SITEID',sp$seq,'ASEQ','AESEQ','MBSEQ','LBSEQ','SRCSEQ','SRCVAR','SRCDOM'))
 changes<-list();nchanged<-nclinical<-integer(nrow(pairs));priority<-list()
 if(nrow(pairs))for(field in common){
  ov<-a[[field]][pairs$eua_index];nv<-b[[field]][pairs$bla_index];idx<-which(ov!=nv);if(!length(idx))next
  nchanged[idx]<-nchanged[idx]+1L;if(!field%in%c(ids,administrative))nclinical[idx]<-nclinical[idx]+1L
  z<-pairs[idx,];changes[[length(changes)+1L]]<-data.frame(domain=domain,subject=a$.sid[z$eua_index],arm=arm(a$.sid[z$eua_index]),eua_row=a$.row[z$eua_index],bla_row=b$.row[z$bla_index],match_method=z$match_method,match_confidence=z$match_confidence,eua_window=a$.window[z$eua_index],eua_date=a$.date_raw[z$eua_index],bla_date=b$.date_raw[z$bla_index],eua_entity=a$.entity[z$eua_index],bla_entity=b$.entity[z$bla_index],field=field,change_type=if(field%in%ids)'IDENTIFIER_OR_LINEAGE' else if(field%in%administrative)'SNAPSHOT_METADATA' else if(field%in%sp$result)'RESULT' else if(field%in%sp$serious)'SERIOUSNESS' else 'OTHER_VALUE',eua_value=ov[idx],bla_value=nv[idx],value_transition=ifelse(!nzchar(ov[idx]),'MISSING_TO_VALUE',ifelse(!nzchar(nv[idx]),'VALUE_TO_MISSING','VALUE_CHANGED')),stringsAsFactors=FALSE)
  if(field%in%sp$result){
   oc<-sc_class(ov[idx]);nc<-sc_class(nv[idx]);keep<-oc!=nc&(oc%in%c('POSITIVE_OR_YES','NEGATIVE_OR_NO')|nc%in%c('POSITIVE_OR_YES','NEGATIVE_OR_NO'))
   if(any(keep)){r<-changes[[length(changes)]][keep,];r$flag<-paste0(oc[keep],'_TO_',nc[keep]);priority[[length(priority)+1L]]<-r}
  }
  if(field%in%sp$serious){keep<-sc_norm(ov[idx])%in%c('Y','YES')&!sc_norm(nv[idx])%in%c('Y','YES');if(any(keep)){r<-changes[[length(changes)]][keep,];r$flag<-'SERIOUS_YES_TO_OTHER';priority[[length(priority)+1L]]<-r}}
 }
 template<-data.frame(domain=character(),subject=character(),arm=character(),eua_row=integer(),bla_row=integer(),match_method=character(),match_confidence=character(),eua_window=character(),eua_date=character(),bla_date=character(),eua_entity=character(),bla_entity=character(),field=character(),change_type=character(),eua_value=character(),bla_value=character(),value_transition=character())
 delta<-if(length(changes))do.call(rbind,changes) else template
 flags<-if(length(priority))do.call(rbind,priority) else transform(template,flag=character())
 if(nrow(pairs)){
  cross<-a$.window[pairs$eua_index]=='ON_OR_BEFORE_EUA'&b$.window[pairs$bla_index]%in%c('AFTER_EUA_THROUGH_BLA','AFTER_BLA_CUTOFF')
  if(any(cross)){z<-pairs[cross,];r<-data.frame(domain=domain,subject=a$.sid[z$eua_index],arm=arm(a$.sid[z$eua_index]),eua_row=a$.row[z$eua_index],bla_row=b$.row[z$bla_index],match_method=z$match_method,match_confidence=z$match_confidence,eua_window=a$.window[z$eua_index],eua_date=a$.date_raw[z$eua_index],bla_date=b$.date_raw[z$bla_index],eua_entity=a$.entity[z$eua_index],bla_entity=b$.entity[z$bla_index],field=paste(sp$date,collapse='/'),change_type='DATE_WINDOW',eua_value=a$.date_raw[z$eua_index],bla_value=b$.date_raw[z$bla_index],value_transition='VALUE_CHANGED',flag='DATE_MOVED_BEYOND_EUA_CUTOFF');flags<-rbind(flags,r)}
 }
 sc_csv(delta,file.path(folder,'field_changes.csv'))
 if(nrow(delta)){z<-aggregate(delta$eua_row,delta[c('eua_window','field','change_type','value_transition')],function(v)length(unique(v)));names(z)[ncol(z)]<-'records';sc_csv(z,file.path(folder,'field_change_summary.csv'))}else sc_csv(data.frame(eua_window=character(),field=character(),change_type=character(),value_transition=character(),records=integer()),file.path(folder,'field_change_summary.csv'))
 sc_csv(flags,file.path(folder,'priority_changes.csv'))
 crosswalk<-data.frame(domain=rep(domain,nrow(pairs)),subject=a$.sid[pairs$eua_index],arm=arm(a$.sid[pairs$eua_index]),eua_row=a$.row[pairs$eua_index],bla_row=b$.row[pairs$bla_index],eua_sequence=a$.seq[pairs$eua_index],bla_sequence=b$.seq[pairs$bla_index],eua_date=a$.date_raw[pairs$eua_index],bla_date=b$.date_raw[pairs$bla_index],eua_window=a$.window[pairs$eua_index],bla_window=b$.window[pairs$bla_index],match_method=pairs$match_method,match_confidence=pairs$match_confidence,changed_fields=nchanged,changed_content_fields=nclinical,stringsAsFactors=FALSE)
 sc_csv(crosswalk,file.path(folder,'record_crosswalk.csv'))
 # Export source values for every unresolved row, including ALL SAE / AE terms.
 ua<-a[m$unmatched_eua,,drop=FALSE];ua$.arm<-arm(ua$.sid);ua$.candidate_count<-q$counts;ua$.candidate_list_truncated<-q$truncated
 ua$.review_status<-ifelse(q$counts>0,'UNRESOLVED_WITH_POSSIBLE_COUNTERPART','NO_COUNTERPART_IDENTIFIED')
 ua$.original_positive<-Reduce('|',lapply(sp$result,function(n)sc_class(sc_col(ua,n))=='POSITIVE_OR_YES'),init=rep(FALSE,nrow(ua)))
 ua$.original_serious<-sc_norm(sc_first(ua,sp$serious))%in%c('Y','YES')
 j<-match(ua$.sid,roster$subject);ua$.subject_in_bla_adsl<-roster$bla_present[j]
 ub<-b[m$unmatched_bla,,drop=FALSE];ub$.arm<-arm(ub$.sid);ub$.review_status<-rep('UNMATCHED_PRODUCTION_RECORD',nrow(ub))
 sc_csv(ua,file.path(folder,'unmatched_eua_records.csv'));sc_csv(ub,file.path(folder,'unmatched_bla_records.csv'))
 if(!ncol(q$rows))q$rows<-data.frame(eua_index=integer(),bla_index=integer(),subject=character(),eua_row=integer(),bla_row=integer(),reason=character(),bla_already_matched=logical(),eua_date=character(),bla_date=character(),eua_entity=character(),bla_entity=character(),eua_result=character(),bla_result=character(),eua_result_class=character(),bla_result_class=character())
 sc_csv(q$rows,file.path(folder,'unresolved_match_candidates.csv'))
 # One row per record: no repeated symptom/test rows are treated as independent subjects.
 subjects<-union(a$.sid,b$.sid);sc_csv(data.frame(subject=subjects,arm=arm(subjects),eua_records=tabulate(match(a$.sid,subjects),nbins=length(subjects)),bla_records=tabulate(match(b$.sid,subjects),nbins=length(subjects))),file.path(folder,'subject_record_counts.csv'))
 changed_subjects<-unique(crosswalk$subject[crosswalk$changed_content_fields>0&crosswalk$eua_window=='ON_OR_BEFORE_EUA'])
 info<-data.frame(domain=domain,eua_records=nrow(a),bla_records=nrow(b),matched_records=nrow(pairs),matched_reference=sum(pairs$match_confidence=='REFERENCE_LINKED'),matched_provisional=sum(pairs$match_confidence=='PROVISIONAL_CLINICAL_LINK'),matched_content_only=sum(pairs$match_confidence=='CONTENT_RETENTION_ONLY'),matched_records_with_content_changes=sum(nclinical>0),eua_pre_cutoff_records=sum(a$.window=='ON_OR_BEFORE_EUA'),eua_pre_cutoff_changed_records=sum(nclinical>0&a$.window[pairs$eua_index]=='ON_OR_BEFORE_EUA'),eua_pre_cutoff_changed_subjects=length(changed_subjects),unmatched_eua=nrow(ua),unmatched_eua_pre_cutoff=sum(ua$.window=='ON_OR_BEFORE_EUA'),unmatched_eua_with_candidates=sum(q$counts>0),unmatched_bla=nrow(ub),unmatched_bla_pre_eua=sum(ub$.window=='ON_OR_BEFORE_EUA'),unmatched_bla_after_eua=sum(ub$.window=='AFTER_EUA_THROUGH_BLA'),unmatched_bla_after_bla=sum(ub$.window=='AFTER_BLA_CUTOFF'),unmatched_bla_unknown_date=sum(ub$.window=='UNKNOWN_OR_PARTIAL_DATE'))
 stopifnot(nrow(pairs)+nrow(ua)==nrow(a),nrow(pairs)+nrow(ub)==nrow(b))
 byarm<-list();arms<-sort(unique(c(arm(a$.sid),arm(b$.sid))))
 for(ar in arms){oldids<-which(arm(a$.sid)==ar);pre<-oldids[a$.window[oldids]=='ON_OR_BEFORE_EUA'];z<-crosswalk[crosswalk$arm==ar&crosswalk$eua_window=='ON_OR_BEFORE_EUA',,drop=FALSE];u<-ua[ua$.arm==ar&ua$.window=='ON_OR_BEFORE_EUA',,drop=FALSE];f<-flags[flags$arm==ar&flags$eua_window=='ON_OR_BEFORE_EUA',,drop=FALSE]
  byarm[[length(byarm)+1L]]<-data.frame(domain=domain,arm=ar,eua_pre_cutoff_records=length(pre),eua_pre_cutoff_subjects=length(unique(a$.sid[pre])),changed_records=sum(z$changed_content_fields>0),changed_subjects=length(unique(z$subject[z$changed_content_fields>0])),unmatched_records=nrow(u),unmatched_subjects=length(unique(u$.sid)),unmatched_positive_records=sum(u$.original_positive),unmatched_serious_records=sum(u$.original_serious),positive_to_negative_records=length(unique(f$eua_row[f$flag=='POSITIVE_OR_YES_TO_NEGATIVE_OR_NO'])),positive_to_negative_subjects=length(unique(f$subject[f$flag=='POSITIVE_OR_YES_TO_NEGATIVE_OR_NO'])),stringsAsFactors=FALSE)
 }
 list(summary=info,by_arm=do.call(rbind,byarm),flags=flags,schema=schema,crosswalk=crosswalk,unmatched=ua)
}
sc_escape <- function(x){x<-as.character(x);x[is.na(x)]<-'';for(p in list(c('&','&amp;'),c('<','&lt;'),c('>','&gt;'),c('"','&quot;')))x<-gsub(p[1],p[2],x,fixed=TRUE);x}
sc_table <- function(x){if(!nrow(x))return('<p>No records in this category.</p>');paste0('<div class="scroll"><table><thead><tr>',paste0('<th>',sc_escape(names(x)),'</th>',collapse=''),'</tr></thead><tbody>',paste(apply(x,1,function(r)paste0('<tr>',paste0('<td>',sc_escape(r),'</td>',collapse=''),'</tr>')),collapse=''),'</tbody></table></div>')}
sc_report <- function(results,coverage,manifest,roster,out){
 summary<-if(length(results))do.call(rbind,lapply(results,`[[`,'summary')) else data.frame()
 byarm<-if(length(results))do.call(rbind,lapply(results,`[[`,'by_arm')) else data.frame()
 flags<-if(length(results))do.call(rbind,lapply(results,`[[`,'flags')) else data.frame()
 if(ncol(summary))sc_csv(summary,file.path(out,'comparison_summary.csv'))
 if(ncol(byarm))sc_csv(byarm,file.path(out,'changes_by_arm.csv'))
 if(ncol(flags))sc_csv(flags,file.path(out,'priority_changes_all_domains.csv'))
 css<-'body{font:15px/1.55 Arial,sans-serif;color:#192c3b;background:#eff3f6;margin:0}main{max-width:1250px;margin:24px auto;padding:34px 42px;background:white}h1{font-size:29px}h2{font-size:21px;margin-top:34px}h3{font-size:17px}table{border-collapse:collapse;width:100%;font-size:12px}th,td{padding:8px;border-bottom:1px solid #d7e0e7;text-align:left;vertical-align:top;overflow-wrap:anywhere}th{background:#e8eff5}tr:nth-child(even){background:#f6f8fa}.scroll{overflow:auto;margin:18px 0}.note{padding:16px;background:#edf5fb;border-left:4px solid #326887}.warning{background:#fff2dc;border-color:#b58024}.meta{color:#5a6e7d;font-size:13px}a{color:#195b83}details{margin:18px 0}summary{cursor:pointer;font-weight:bold}code{overflow-wrap:anywhere}nav a{margin-right:16px}button{padding:8px 12px}@media print{body{background:white;font-size:10pt}main{margin:0;padding:0}nav,button{display:none}.scroll{overflow:visible}table{font-size:8pt}thead{display:table-header-group}tr{break-inside:avoid}h2{break-after:avoid;font-size:15pt}h1{font-size:20pt}}@page{size:A4 landscape;margin:13mm}'
 parts<-c('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>P301 EUA / BLA snapshot reconciliation</title>',paste0('<style>',css,'</style></head><body><main>'),'<h1>EUA / BLA snapshot reconciliation</h1>',paste0('<p class="meta">Moderna P301 · v1.0.0 · EUA reference ',PHMPT_EUA_CUTOFF,' · BLA target ',PHMPT_BLA_CUTOFF,' · Generated ',format(Sys.time(),tz='UTC',usetz=TRUE),'</p>'),
 '<p class="note">This report identifies differences between selected released tables. A changed value is a snapshot discrepancy; an unmatched record has no confidently identified counterpart under these rules. Neither finding alone establishes deletion, an undocumented correction or misconduct. No electronic data-capture audit trail is supplied by this comparison.</p>',
 '<nav><a href="#coverage">Coverage</a><a href="#overview">Reconciliation</a><a href="#arms">By arm</a><a href="#priority">Result / seriousness changes</a><a href="#domains">Record review</a><a href="#methods">Methods</a></nav><button onclick="window.print()">Print report</button>',
 '<h2 id="coverage">1. Files and snapshot coverage</h2>',sc_table(coverage),'<p>Only the selected submissions are compared. Added follow-up is distinguished from records newly appearing with dates on or before the EUA reference date. Unknown/partial dates remain separate. The BLA label refers to the selected production submission; its exact analysis cutoff must be checked in the manifest and ADSL coverage below.</p>')
 if(!nrow(roster))parts<-c(parts,'<p class="note warning">ADSL arm mapping is unavailable. Record reconciliation remains possible, but treatment-arm summaries are unmapped. Review the ADSL discovery logs.</p>')
 if(nrow(roster)){
  rr<-data.frame(Metric=c('EUA ADSL subjects','BLA ADSL subjects','EUA subjects absent from BLA ADSL','Conflicting arm assignments','Latest EUA ADSL cutoff','Latest BLA ADSL cutoff'),Value=c(sum(roster$eua_present),sum(roster$bla_present),sum(roster$eua_present&!roster$bla_present),sum(roster$arm_conflict),if(any(nzchar(na.omit(roster$eua_cutoff))))max(roster$eua_cutoff,na.rm=TRUE) else 'Unavailable',if(any(nzchar(na.omit(roster$bla_cutoff))))max(roster$bla_cutoff,na.rm=TRUE) else 'Unavailable'))
  parts<-c(parts,sc_table(rr))
 }
 if(any(coverage$status!='COMPARED'))parts<-c(parts,'<p class="note warning">Coverage is incomplete for the domains marked unavailable or failed. Absence of a source table is not interpreted as absence of its events.</p>')
 parts<-c(parts,'<h2 id="overview">2. Record reconciliation</h2>',sc_table(if(nrow(summary))data.frame(Domain=summary$domain,`EUA records`=summary$eua_records,`Production records`=summary$bla_records,Matched=summary$matched_records,`Changed pre-EUA records`=summary$eua_pre_cutoff_changed_records,`Changed pre-EUA subjects`=summary$eua_pre_cutoff_changed_subjects,`Unmatched pre-EUA records`=summary$unmatched_eua_pre_cutoff,check.names=FALSE) else data.frame()),'<p class="meta">Changed-content counts exclude identifier/lineage fields and snapshot cutoff metadata. They still include coding updates and previously empty fields being populated; use each domain’s field-change summary to distinguish these from result reversals.</p><h3>Timing of unmatched production records</h3>',sc_table(if(nrow(summary))data.frame(Domain=summary$domain,`Dated on/before EUA`=summary$unmatched_bla_pre_eua,`After EUA through BLA`=summary$unmatched_bla_after_eua,`After BLA target`=summary$unmatched_bla_after_bla,`Unknown/partial date`=summary$unmatched_bla_unknown_date,check.names=FALSE) else data.frame()),'<p><a href="comparison_summary.csv">Summary CSV</a> · <a href="input_manifest.csv">Input files and MD5 fingerprints</a> · <a href="subject_roster_comparison.csv">Subject/arm roster comparison</a></p>',
 '<h2 id="arms">3. Changes affecting records dated on or before EUA, by arm</h2>',sc_table(byarm),'<p class="meta">The denominator includes all P301 source records, including non-randomized subjects if present. Arm assignment uses EUA ADSL where available; BLA-only assignments and conflicts are identified in the roster. Both record and unique-subject counts are provided. Derived domains overlap the source domains and must not be added together.</p>',
 '<h2 id="priority">4. Result reversals, seriousness changes and dates crossing the EUA cutoff</h2>')
 if(nrow(flags)){
  keep<-flags$eua_window=='ON_OR_BEFORE_EUA';f<-flags[keep,,drop=FALSE]
  if(nrow(f)){
   keys<-c('domain','eua_entity','arm','match_confidence','flag');count<-aggregate(f$eua_row,f[keys],function(x)length(unique(x)));names(count)[ncol(count)]<-'Records'
   parts<-c(parts,sc_table(count),'<p>Reference-linked and provisional clinical links are shown separately. Raw and standardized result fields may both flag one record; the summary counts unique EUA rows within each category.</p>',sc_table(head(f,PHMPT_HTML_EXAMPLES)),paste0('<p class="meta">Showing at most ',PHMPT_HTML_EXAMPLES,' flagged field changes. All rows, including later or undated records, are in <a href="priority_changes_all_domains.csv">priority_changes_all_domains.csv</a>.</p>'))
  }else parts<-c(parts,'<p>No priority changes among successfully matched records dated on or before EUA. Check unmatched and ambiguous candidates before drawing a completeness conclusion.</p>')
 }else parts<-c(parts,'<p>No priority changes among matched records. Unmatched and ambiguous records still require review.</p>')
 parts<-c(parts,'<h2 id="domains">5. Domain-level review files</h2>')
 for(d in names(results)){
  r<-results[[d]];u<-r$unmatched;u<-u[u$.window=='ON_OR_BEFORE_EUA'&(u$.original_positive|u$.original_serious),,drop=FALSE]
  viewcols<-intersect(c('.sid','.arm','.row','.date_raw','.entity','.original_positive','.original_serious','.candidate_count','.review_status'),names(u))
  parts<-c(parts,paste0('<h3>',d,' — ',sc_escape(sc_spec[[d]]$kind),'</h3><p><a href="',d,'/record_crosswalk.csv">Matched records</a> · <a href="',d,'/field_changes.csv">Every changed shared field</a> · <a href="',d,'/field_change_summary.csv">Changes by field</a> · <a href="',d,'/unmatched_eua_records.csv">Unmatched EUA records</a> · <a href="',d,'/unmatched_bla_records.csv">Unmatched production records</a> · <a href="',d,'/unresolved_match_candidates.csv">Ambiguous / sequence-only candidates</a> · <a href="',d,'/schema_comparison.csv">Schema and label differences</a> · <a href="',d,'/subject_record_counts.csv">Participant record counts</a></p>'),
   '<details><summary>Unmatched pre-EUA positive / serious records requiring review</summary>',sc_table(head(u[viewcols],PHMPT_HTML_EXAMPLES)),'</details>')
 }
 parts<-c(parts,'<h2 id="methods">6. Comparison rules and limitations</h2><ul>',
 '<li>Matching starts with unique subject-specific reference IDs (and test identity for findings). Then it uses a unique subject/date/event or test key, without the result value. These clinical matches are provisional. Finally, identical shared content is reconciled as a multiset, preserving duplicate multiplicity.</li>',
 '<li>Uniqueness is checked in the complete snapshots. Sequence numbers alone never establish identity. Unmatched rows with the same sequence, collection/assessment day, nearby date with the same event/test, or ambiguous reference appear in the candidate table; a candidate can already be matched elsewhere, which helps identify splitting, merging or reassignment.</li>',
 '<li>Every shared XPT column is compared after removal of import labels/attributes and outer whitespace. Numeric values retain precision; tagged SAS missing values remain distinct. Record order is immaterial. New or removed columns and their labels/types are reported separately, not turned into changes to every row.</li>',
 '<li>All production dates remain available during matching, so a date moved beyond EUA can be detected. Record timing uses AE onset, microbiology/lab collection date, or symptom assessment date. This is not a data-entry timestamp; a pre-EUA date newly present in BLA cannot establish when the record was entered.</li>',
 '<li>AE compares all adverse events, not only COVID-compatible AEs. ADSYMP/ADSY, ADMB and ADAE are derived tables: changed derivation or source selection can change their rows. Comparison with source domains and the reviewer guides is needed before calling an event removed.</li>',
 '<li>Positive/negative flags use explicit qualitative codes only. Numeric laboratory values are compared as values, without inventing positivity thresholds. Positive/negative changes in symptoms mean yes/no assessment changes, not infection status. Multiple findings can concern one specimen or illness.</li>',
 '<li>No baseline-negative, randomized-only, per-protocol or blinded-window exclusion is imposed on the source-record reconciliation. Records missing from a later release are not silently removed from denominators. No efficacy calculation or causal attribution is made.</li></ul>',
 '<p>Review priority: trace reference-linked result reversals and serious records first, then provisional matches and unmatched historical records; compare raw/standardized values, dates, test identity, specimen/reference IDs and source-versus-derived records. Confirm explanations against source documents or an audit trail.</p></main></body></html>')
 writeLines(enc2utf8(parts),file.path(out,'snapshot_changes_report.html'),useBytes=TRUE)
}
analyse_eua_bla_changes <- function(root=PHMPT_ROOT,output_dir=PHMPT_OUTPUT_DIR,show_report=PHMPT_SHOW_REPORT){
 root<-sc_root(root);if(any(!PHMPT_DOMAINS%in%names(sc_spec)))stop('Unsupported PHMPT_DOMAINS')
 if(!is.finite(PHMPT_CANDIDATE_DATE_TOLERANCE)||PHMPT_CANDIDATE_DATE_TOLERANCE<0)stop('Invalid candidate date tolerance')
 if(!is.finite(PHMPT_MAX_CANDIDATES_PER_RECORD)||PHMPT_MAX_CANDIDATES_PER_RECORD<1)stop('Invalid candidate detail limit')
 if(is.na(as.Date(PHMPT_EUA_CUTOFF))||is.na(as.Date(PHMPT_BLA_CUTOFF))||as.Date(PHMPT_EUA_CUTOFF)>as.Date(PHMPT_BLA_CUTOFF))stop('Invalid cutoffs')
 out<-if(nzchar(output_dir))sc_path(root,output_dir) else file.path(root,'analysis',paste0('snapshot_changes_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_',Sys.getpid()))
 if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop('Output directory must be empty: ',out)
 dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out,winslash='/');writeLines('RUNNING',file.path(out,'RUN_STATUS.txt'))
 complete<-FALSE;on.exit(if(!complete)writeLines('FAILED - see coverage and console',file.path(out,'RUN_STATUS.txt')),add=TRUE)
 manifest<-list();coverage<-list();results<-list();register<-function(z){if(is.null(z))return();manifest[[length(manifest)+1L]]<<-data.frame(release=z$release,domain=z$domain,path=z$path,md5=z$md5,file_bytes=file.info(z$path)$size,original_rows=z$original_rows,p301_rows=nrow(z$data),stringsAsFactors=FALSE)}
 ra<-tryCatch(sc_pick(root,'EUA','ADSL',out),error=function(e){message('EUA arm roster unavailable: ',conditionMessage(e));NULL});rb<-tryCatch(sc_pick(root,'BLA','ADSL',out),error=function(e){message('BLA arm roster unavailable: ',conditionMessage(e));NULL});register(ra);register(rb)
 roster<-sc_rosters(ra,rb,out)
 for(domain in PHMPT_DOMAINS){
  message('Comparing domain ',domain)
  answer<-tryCatch({a<-sc_pick(root,'EUA',domain,out);b<-sc_pick(root,'BLA',domain,out);register(a);register(b)
   if(is.null(a)||is.null(b))list(status='UNAVAILABLE',detail=paste('Missing',paste(c(if(is.null(a))'EUA',if(is.null(b))'BLA'),collapse=' / '),'table'))
   else {z<-sc_compare(a,b,domain,roster,out);list(status='COMPARED',detail=paste(nrow(a$data),'EUA /',nrow(b$data),'production rows'),result=z)}
  },error=function(e)list(status='FAILED',detail=conditionMessage(e)))
  coverage[[domain]]<-data.frame(domain=domain,required=domain%in%PHMPT_REQUIRED_DOMAINS,status=answer$status,detail=answer$detail,stringsAsFactors=FALSE)
  if(answer$status=='COMPARED')results[[domain]]<-answer$result else message(domain,': ',answer$detail)
 }
 cv<-do.call(rbind,coverage);mf<-if(length(manifest))do.call(rbind,manifest) else data.frame(release=character(),domain=character(),path=character(),md5=character(),file_bytes=numeric(),original_rows=integer(),p301_rows=integer())
 sc_csv(cv,file.path(out,'coverage.csv'));sc_csv(mf,file.path(out,'input_manifest.csv'))
 sc_csv(data.frame(parameter=c('version','EUA_cutoff','BLA_target_cutoff','EUA_submission','BLA_submission','allow_csv','candidate_detail_limit','candidate_date_tolerance'),value=c('1.0.0',PHMPT_EUA_CUTOFF,PHMPT_BLA_CUTOFF,PHMPT_SUBMISSIONS['EUA'],PHMPT_SUBMISSIONS['BLA'],PHMPT_ALLOW_CSV_FALLBACK,PHMPT_MAX_CANDIDATES_PER_RECORD,PHMPT_CANDIDATE_DATE_TOLERANCE)),file.path(out,'parameters.csv'))
 capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'));sc_report(results,cv,mf,roster,out)
 failed<-any(cv$status=='FAILED'|cv$required&cv$status!='COMPARED');status<-if(failed)'PARTIAL - required domain missing or failed' else if(any(cv$status!='COMPARED'))'COMPLETE_WITH_OPTIONAL_GAPS' else 'COMPLETE'
 writeLines(status,file.path(out,'RUN_STATUS.txt'));complete<-TRUE
 message('Report: ',file.path(out,'snapshot_changes_report.html'))
 if(isTRUE(show_report))try(utils::browseURL(file.path(out,'snapshot_changes_report.html')),silent=TRUE)
 if(failed)warning('Partial comparison: inspect coverage.csv and discovery logs. Missing/failed domains are not negative findings.')
 invisible(list(output_dir=out,coverage=cv,results=results))
}
# JMP ADSL reader is embedded below; no dependency on another analysis script.

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







# ---- RUN: executes when the complete standalone script is run ----
if(isTRUE(getOption("phmpt.autorun",TRUE)))analyse_eua_bla_changes()
