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
 sc_audit_auxiliary(root,out,ra,rb)
 sc_export_snapshot_assets(root,out)
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
# Audit PNG export extension v1.0.0 | 2026-10-23
# Standalone; stable Vxx filenames in <project>/audit_assets plus this run's copy.
# Long tables continue as _02.png, _03.png, etc. References remain in HTML/provenance.
PHMPT_AUDIT_ASSET_DIR <- getOption("phmpt.audit_asset_dir", "audit_assets")
PHMPT_AUDIT_DPI <- getOption("phmpt.audit_dpi", 300L)

# ---- Standalone audit-asset renderer (V08-V35), version 1.0.0 ----
# Embedded in each script: no shared R file or browser is required.
ax_clean <- function(x) {
  x<-as.character(x);x[is.na(x)]<-"Unavailable"
  x<-gsub("https?://[^[:space:]<>]+","",x,perl=TRUE)
  x<-gsub("(?i)doi:[[:space:]]*[^[:space:]]+","",x,perl=TRUE)
  x<-gsub("\\[[0-9]+([,; -]+[0-9]+)*\\]","",x,perl=TRUE)
  trimws(x)
}
ax_path <- function(root,p)if(grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)",p))p else file.path(root,p)
ax_root <- function(root="",out=getwd()) {
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
ax_read <- function(out,file,required=TRUE) {
  p<-file.path(out,file)
  if(!file.exists(p)){if(required)stop("Missing export: ",p);return(data.frame())}
  read.csv(p,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
}
ax_csv <- function(x,path) {
  con<-file(path,"wb");on.exit(close(con));writeBin(as.raw(c(239,187,191)),con)
  write.table(x,con,sep=",",row.names=FALSE,quote=TRUE,na="",qmethod="double")
}
ax_num <- function(x)ifelse(is.na(x),"Unavailable",format(x,big.mark=",",scientific=FALSE,trim=TRUE))
ax_rate <- function(k,n,digits=1)ifelse(is.na(k)|is.na(n),"Unavailable",ifelse(n>0,
  paste0(ax_num(k),"/",ax_num(n),sprintf(paste0(" (%.",digits,"f%%)"),100*k/n)),"0/0 (not estimable)"))
ax_p <- function(p) {
  z<-ifelse(is.finite(p),formatC(p,digits=4,format="f"),"Not estimable")
  for(i in 2:6)z[is.finite(p)&p<10^-i]<-paste0("< ",format(10^-i,scientific=FALSE,trim=TRUE))
  z
}
ax_text <- function(s,x,y,size=11.5,bold=FALSE,col="#25343D",just=c("left","top"))
  grid::grid.text(ax_clean(s),x=grid::unit(x,"in"),y=grid::unit(1,"npc")-grid::unit(y,"in"),just=just,
    gp=grid::gpar(fontfamily="sans",fontsize=size,fontface=if(bold)"bold" else "plain",col=col,lineheight=1.12))
ax_height <- function(s,size=11.5)grid::convertHeight(grid::grobHeight(grid::textGrob(s,
  gp=grid::gpar(fontfamily="sans",fontsize=size,lineheight=1.12))),"in",valueOnly=TRUE)
ax_wrap <- function(s,width,size=11.5,bold=FALSE) {
  measure<-function(t)grid::convertWidth(grid::grobWidth(grid::textGrob(t,
    gp=grid::gpar(fontfamily="sans",fontsize=size,fontface=if(bold)"bold" else "plain"))),"in",valueOnly=TRUE)
  paragraphs<-strsplit(ax_clean(s),"\n",fixed=TRUE)[[1]];all<-character()
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
ax_begin <- function(root,out,owner,dpi=PHMPT_AUDIT_DPI,asset_dir=PHMPT_AUDIT_ASSET_DIR) {
  if(length(dpi)!=1||!is.finite(dpi)||dpi<150||dpi>600)stop("Audit DPI must be 150-600.")
  if(length(asset_dir)!=1||is.na(asset_dir)||!nzchar(asset_dir))stop("Audit asset directory is empty.")
  a<-new.env(parent=emptyenv());a$root<-ax_root(root,out);a$out<-normalizePath(out,winslash="/",mustWork=TRUE)
  a$owner<-owner;a$dpi<-dpi;a$dest<-ax_path(a$root,asset_dir)
  a$stage<-tempfile(paste0("audit_",owner,"_"),tmpdir=out);dir.create(a$stage)
  a$rows<-list();a$files<-character();a
}
ax_png <- function(a,id,slug,draw,height=7,page=1L,pages=1L,status="AVAILABLE",rows=NA_integer_) {
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
ax_data <- function(a,id,x,suffix="data") {
  file<-paste0(id,"_",suffix,".csv");ax_csv(x,file.path(a$stage,file));a$files<-c(a$files,file)
}
ax_table <- function(a,id,slug,title,data,subtitle="",note="",widths=NULL,status="AVAILABLE",body_size=11.5,highlight=integer()) {
  data<-as.data.frame(data,check.names=FALSE);ax_data(a,id,data)
  if(!ncol(data))data<-data.frame(Result="No matching rows",check.names=FALSE)
  if(!nrow(data)){data[1,]<-rep("No matching rows",ncol(data));status<-"NO_MATCHING_ROWS"}
  if(is.null(widths))widths<-if(ncol(data)==3)c(.5,.25,.25) else rep(1/ncol(data),ncol(data))
  widths<-widths/sum(widths)*7.8;left<-.35+c(0,head(cumsum(widths),-1))
  tmp<-tempfile(fileext=".pdf");if(capabilities("cairo"))grDevices::cairo_pdf(tmp,width=8.5,height=13,family="sans") else grDevices::pdf(tmp,width=8.5,height=13)
  dev<-grDevices::dev.cur();on.exit({if(dev%in%grDevices::dev.list())grDevices::dev.off(dev);unlink(tmp)})
  wraprow<-function(z,size=body_size,bold=FALSE)vapply(seq_along(z),function(j)ax_wrap(z[j],widths[j]-.16,size,bold),character(1))
  heads<-wraprow(names(data),11.2,TRUE);head_h<-max(vapply(heads,ax_height,numeric(1),size=11.2))+.22
  title<-ax_wrap(paste(id,title),7.8,17,TRUE);sub<-ax_wrap(subtitle,7.8,11.2);foot<-ax_wrap(note,7.8,10.6)
  top<-.25+ax_height(title,17)+.18+if(nzchar(sub))ax_height(sub,11.2)+.18 else 0
  bottom<-.42+if(nzchar(foot))ax_height(foot,10.6) else 0
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
      row[1]<-ax_wrap(row[1],widths[1]-.16,body_size)
      h<-max(vapply(row,ax_height,numeric(1),size=body_size))+.24
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
      grid::grid.newpage();ax_text(title,.35,.25,17,TRUE,col="#17384B")
      sy<-.25+ax_height(title,17)+.18;if(nzchar(sub))ax_text(sub,.35,sy,11.2)
      y<-top
      for(r in c(0L,ix)) {
        h<-if(r==0)head_h else heights[r]
        grid::grid.rect(x=grid::unit(.35,"in"),y=grid::unit(1,"npc")-grid::unit(y,"in"),
          width=grid::unit(7.8,"in"),height=grid::unit(h,"in"),just=c("left","top"),
          gp=grid::gpar(fill=if(r==0)"#E8EFF3" else if(origin[r]%in%highlight)"#FFF0D9" else if(origin[r]%%2==0)"#F3F6F8" else "white",col=NA))
        values<-if(r==0)heads else cells[[r]]
        for(j in seq_along(values))ax_text(values[j],left[j]+.08,y+.11,if(r==0)11.2 else body_size,r==0||(r>0&&j==1&&origin[r]%in%highlight),
          col=if(r!=0&&grepl("mRNA|Vaccine",names(data)[j],ignore.case=TRUE))"#1E667D" else if(r!=0&&grepl("Placebo",names(data)[j]))"#A65B22" else "#25343D")
        y<-y+h
      }
      if(nzchar(foot))ax_text(foot,.35,y+.16,10.6,col="#52616C")
      ax_text(paste0("Page ",p," of ",length(pages)),8.15,height-.17,9.5,col="#64737D",just=c("right","bottom"))
    }
    ax_png(a,id,slug,draw,height,p,length(pages),status,length(unique(origin[ix])))
  };invisible(NULL)
}
ax_unavailable <- function(a,id,slug,title,reason)
  ax_table(a,id,slug,title,data.frame(Status="Not available in this run",Reason=reason,check.names=FALSE),widths=c(.3,.7),status="UNAVAILABLE")
ax_finish <- function(a,notes=character()) {
  if(!length(a$rows))stop("No audit assets were rendered.")
  manifest<-do.call(rbind,a$rows);manifest$renderer_version<-"1.0.0";manifest$generated_at<-format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z")
  manifest$run_directory<-a$out;manifest$md5<-unname(tools::md5sum(file.path(a$stage,manifest$file)))
  mf<-paste0(a$owner,"_manifest.csv");ax_csv(manifest,file.path(a$stage,mf))
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
    oldmf<-file.path(dest,mf);old<-if(file.exists(oldmf))ax_read(dest,mf)$file else character()
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

sc_audit_aux_read <- function(path,domain="ADSL") {
  fields<-if(domain=="ADSL")c("STUDYID","USUBJID","RANDFL","TRT01P","CUTOFFDT") else c("STUDYID","USUBJID")
  optional<-c("FASFL","SAFFL","MITTFL","PPROTFL","DOSE1FL","DOSE2FL","TR01SDT","DOSE2DT","RANDDT","FASEREAS","SAFEREAS","SITEID","ARM")
  ext<-tolower(tools::file_ext(path))
  x<-if(ext=="jmp")ws_read_jmp(path,fields,optional) else if(ext=="csv")
    tryCatch(read.csv(path,check.names=FALSE,colClasses="character",fileEncoding="UTF-8-BOM"),error=function(e)read.csv(path,check.names=FALSE,colClasses="character",fileEncoding="windows-1252")) else {
      if(!requireNamespace("haven",quietly=TRUE))stop("Package haven required for ",path)
      if(ext=="xpt")haven::read_xpt(path) else if(ext=="sas7bdat")haven::read_sas(path) else stop("Unsupported auxiliary format.")
    }
  x<-as.data.frame(x);names(x)<-toupper(names(x));if(!all(fields%in%names(x)))stop("Missing auxiliary fields: ",path)
  x<-x[sc_norm(x$STUDYID)=="MRNA-1273-P301",intersect(c(fields,optional),names(x)),drop=FALSE]
  if(!nrow(x))stop("No P301 rows in auxiliary input.")
  for(n in intersect(c("TR01SDT","DOSE2DT","RANDDT","CUTOFFDT"),names(x))) {
    if(inherits(x[[n]],"POSIXt"))x[[n]]<-format(x[[n]],"%Y-%m-%d")
    else if(inherits(x[[n]],"Date"))x[[n]]<-format(x[[n]],"%Y-%m-%d")
    else if(is.numeric(x[[n]]))x[[n]]<-format(as.Date(x[[n]],origin="1960-01-01"),"%Y-%m-%d")
    else {ii<-!is.na(x[[n]])&grepl("^-?[0-9]+([.]0+)?$",x[[n]]);if(any(ii)){if(!nzchar(PHMPT_CSV_DATE_ORIGIN))stop("Set PHMPT_CSV_DATE_ORIGIN for numeric auxiliary CSV dates.");x[[n]][ii]<-format(as.Date(as.numeric(x[[n]][ii]),origin=PHMPT_CSV_DATE_ORIGIN),"%Y-%m-%d")}}
  }
  x$subject<-sc_sid(x$USUBJID);if(anyDuplicated(x$subject))stop("Duplicate auxiliary subject IDs.")
  x
}
sc_audit_auxiliary <- function(root,out,ra,rb) {
  unlink(file.path(out,c("audit_EUA_ADSL_fields.csv","audit_BLA_ADSL_fields.csv","audit_BLA_DM_subjects.csv","audit_auxiliary_inputs.csv")))
  log<-list();add<-function(domain,path,status,detail=""){log[[length(log)+1L]]<<-data.frame(domain=domain,path=path,status=status,detail=detail,md5=if(file.exists(path))unname(tools::md5sum(path)) else "")}
  for(rel in c("EUA","BLA")) {
    item<-if(rel=="EUA")ra else rb
    if(is.null(item)){add(paste(rel,"ADSL"),"","UNAVAILABLE");next}
    x<-tryCatch(sc_audit_aux_read(item$path),error=function(e)e)
    if(inherits(x,"error")){add(paste(rel,"ADSL"),item$path,"UNAVAILABLE",conditionMessage(x));next}
    ax_csv(x,file.path(out,paste0("audit_",rel,"_ADSL_fields.csv")));add(paste(rel,"ADSL"),item$path,"SELECTED")
  }
  explicit<-PHMPT_FILES$BLA$DM;if(is.null(explicit))explicit<-""
  p<-if(nzchar(explicit))sc_path(root,explicit) else {
    all<-list.files(file.path(root,"prod_data"),recursive=TRUE,full.names=TRUE,pattern="[.](csv|xpt|sas7bdat)$",ignore.case=TRUE)
    all[grepl("(^|[_ -])dm([_. -]|$)",basename(all),ignore.case=TRUE)&grepl(paste0(PHMPT_SUBMISSIONS[['BLA']],"_"),basename(all),ignore.case=TRUE)]
  }
  selected<-NULL;conflict<-FALSE
  for(path in p) {
    x<-tryCatch(sc_audit_aux_read(path,"DM"),error=function(e)e)
    if(inherits(x,"error")){add("BLA DM",path,"UNREADABLE",conditionMessage(x));next}
    ids<-sort(x$subject)
    if(is.null(selected)){selected<-ids;add("BLA DM",path,"SELECTED")}
    else if(identical(ids,selected))add("BLA DM",path,"EQUIVALENT_SUBJECT_IDS")
    else {conflict<-TRUE;add("BLA DM",path,"CONFLICTING_SUBJECT_IDS")}
  }
  if(!length(p))add("BLA DM","","UNAVAILABLE","No matching production submission DM")
  if(!is.null(selected)&&!conflict)ax_csv(data.frame(subject=selected),file.path(out,"audit_BLA_DM_subjects.csv"))
  if(conflict)warning("Conflicting BLA DM rosters: DM presence unavailable in audit figures; set PHMPT_FILES$BLA$DM.")
  ax_csv(do.call(rbind,log),file.path(out,"audit_auxiliary_inputs.csv"))
}
sc_export_snapshot_assets <- function(root,out) {
  a<-ax_begin(root,out,"V28_V32_snapshots");r<-ax_read(out,"subject_roster_comparison.csv")
  dm<-ax_read(out,"audit_BLA_DM_subjects.csv",FALSE);old<-ax_read(out,"audit_EUA_ADSL_fields.csv",FALSE);new<-ax_read(out,"audit_BLA_ADSL_fields.csv",FALSE)
  dmcheck<-file.exists(file.path(out,"audit_BLA_DM_subjects.csv"));dmcount<-function(ids)if(dmcheck)as.character(sum(ids%in%dm$subject)) else "Not checked: DM unavailable"
  y<-function(x)!is.na(x)&toupper(x)=="Y";removed<-r$eua_present&!r$bla_present;added<-!r$eua_present&r$bla_present
  d<-data.frame(Measure=c("ADSL subject IDs","Randomized flag = Y","Shared subject IDs","IDs present only in this snapshot","EUA-randomized IDs absent from BLA ADSL","EUA-only ADSL IDs still present in BLA DM","EUA-only ADSL IDs absent from checked BLA DM"),
    EUA=c(ax_num(sum(r$eua_present)),ax_num(sum(y(r$eua_randomized))),ax_num(sum(r$eua_present&r$bla_present)),ax_num(sum(removed)),"Not applicable",ax_num(sum(removed)),"Not applicable"),
    BLA=c(ax_num(sum(r$bla_present)),ax_num(sum(y(r$bla_randomized))),ax_num(sum(r$eua_present&r$bla_present)),ax_num(sum(added)),ax_num(sum(y(r$eua_randomized)&!r$bla_present)),dmcount(r$subject[removed]),if(dmcheck)ax_num(sum(!r$subject[removed]%in%dm$subject)) else "Not checked: DM unavailable"),check.names=FALSE)
  roster_ok<-nrow(r)>0&&any(r$eua_present)&&any(r$bla_present)
  if(roster_ok)ax_table(a,"V28","subject_roster_reconciliation","EUA and BLA subject-roster reconciliation",d,
    "P301 | Selected EUA and production submission snapshots",
    "ADSL membership, randomized flags and DM membership are different checks. An ID missing from ADSL is not necessarily absent from the submission. DM is checked directly only when a compatible DM input is readable.")
  else ax_unavailable(a,"V28","subject_roster_reconciliation","EUA and BLA subject-roster reconciliation","Both readable ADSL rosters are required; unavailable releases are not zero participants.")
  z<-r[roster_ok&((r$eua_present&r$bla_present&y(r$eua_randomized)!=y(r$bla_randomized))|(!r$eua_present&r$bla_present&y(r$bla_randomized))),]
  if(nrow(z)) {
    value<-function(data,ids,col){if(!col%in%names(data))return(rep("Unavailable",length(ids)));v<-as.character(data[[col]][match(ids,data$subject)]);v[is.na(v)]<-"Unavailable";v[v==""]<-"Blank";v}
    # Context is transcribed from the reviewer guide, never inferred from an AE/flag.
    duplicates<-c("US3132036","US3162112","US3172220","US3352166","US3522027")
    context<-ifelse(z$subject%in%duplicates,"Reviewer-guide context: duplicate / invalid randomization","No reviewer-guide classification supplied")
    d<-data.frame(`Participant / arm`=paste(z$subject,z$arm,sep="\n"),
      `RANDFL EUA -> BLA`=paste(ifelse(z$eua_present,z$eua_randomized,"Absent"),z$bla_randomized,sep=" -> "),
      `BLA dose and population fields`=paste0("Dose 1: ",value(new,z$subject,"TR01SDT"),"\nDose 2: ",value(new,z$subject,"DOSE2DT"),
        "\nFAS / safety: ",value(new,z$subject,"FASFL")," / ",value(new,z$subject,"SAFFL"),"\nmITT / PP: ",value(new,z$subject,"MITTFL")," / ",value(new,z$subject,"PPROTFL")),
      `Recorded reasons / separate context`=paste0("FAS: ",value(new,z$subject,"FASEREAS"),"\nSafety: ",value(new,z$subject,"SAFEREAS"),"\n",context),check.names=FALSE)
    ax_table(a,"V29","randomization_flag_changes","Changes in randomized-population flags",d,
      "Changed RANDFL and newly included randomized IDs | Dose dates and analysis flags from BLA ADSL",
      "Blank means the field was present but empty; unavailable means it could not be read. Reviewer-guide context is a separate annotation, not a dataset-derived reason. Dose dates alone do not resolve analysis-population eligibility.",widths=c(.19,.14,.32,.35),highlight=which(z$subject=="US3952070"))
    ax_data(a,"V29",data.frame(subject=z$subject,context=context,context_source="ADRG reviewer-guide discussion of duplicate/invalid randomization; see report references"),"guide_context")
  }else ax_unavailable(a,"V29","randomization_flag_changes","Randomized-population flags","No changed shared-subject flags, or required roster inputs unavailable.")
  f<-ax_read(out,"AE/field_changes.csv",FALSE);u<-ax_read(out,"AE/unmatched_eua_records.csv",FALSE)
  if(file.exists(file.path(out,"AE/field_changes.csv"))&&file.exists(file.path(out,"AE/unmatched_eua_records.csv"))) {
    f<-f[f$field=="AESER",];f<-f[(f$eua_value=="Y"&f$bla_value=="N")|(f$eua_value=="N"&f$bla_value=="Y"),]
    down<-f[f$eua_value=="Y",];up<-f[f$eua_value=="N",];us<-u[u$.window=="ON_OR_BEFORE_EUA"&u$.original_serious,]
    cnt<-function(x,ar,sid="subject",arm="arm")paste0(sum(x[[arm]]==ar)," records / ",length(unique(x[[sid]][x[[arm]]==ar]))," participants")
    rows<-list("Serious -> non-serious (all linked EUA records)"=down,"Non-serious -> serious (all linked EUA records)"=up)
    d<-do.call(rbind,lapply(names(rows),function(n)data.frame(Change=n,`mRNA-1273`=cnt(rows[[n]],"mRNA-1273"),Placebo=cnt(rows[[n]],"Placebo"),check.names=FALSE)))
    d<-rbind(d,data.frame(Change=c("Pre-EUA serious records with no linked counterpart","Subset with no possible review candidate"),
      `mRNA-1273`=c(cnt(us,"mRNA-1273",".sid",".arm"),cnt(us[us$.candidate_count==0,],"mRNA-1273",".sid",".arm")),
      Placebo=c(cnt(us,"Placebo",".sid",".arm"),cnt(us[us$.candidate_count==0,],"Placebo",".sid",".arm")),check.names=FALSE))
    ax_table(a,"V30","ae_seriousness_changes","AE seriousness changes and unmatched records",d,
      "Record linkage across selected snapshots | AESER indicates seriousness, not severity",
      "Linked seriousness changes include undated records; unmatched rows are restricted to onset on/before the EUA cutoff. Record counts and unique participant counts are shown separately. Possible counterparts require review and are not treated as confirmed matches.")
    ax_data(a,"V30",f,"linked_seriousness_changes")
    counts<-ax_read(out,"AE/subject_record_counts.csv",FALSE)
    v<-us[us$.candidate_count==0&us$.sid%in%counts$subject[counts$bla_records==0],]
    d<-if(!nrow(v))data.frame(`Participant / arm`=character(),`EUA serious event / onset`=character(),`BLA records`=character(),check.names=FALSE) else data.frame(`Participant / arm`=paste(v$.sid,v$.arm,sep="\n"),`EUA serious event / onset`=paste(v$AETERM,v$AESTDTC,sep="\n"),
      `BLA records`=paste0("AE: 0\nADSL: ",if(any(r$bla_present))ifelse(v$.sid%in%r$subject[r$bla_present],"present","absent") else rep("not checked",nrow(v)),"\nDM: ",if(dmcheck)ifelse(v$.sid%in%dm$subject,"present","absent") else "not checked"),check.names=FALSE)
    ax_table(a,"V31","unmatched_serious_events","Serious-event participants with no BLA AE rows",d,
      "Pre-EUA serious event | No linked counterpart or possible review candidate",
      "Zero AE rows refers to the selected production AE snapshot. ADSL and DM membership are independently checked. This identifies records requiring reconciliation; it does not establish deletion from the original clinical database.",widths=c(.25,.49,.26))
  }else {ax_unavailable(a,"V30","ae_seriousness_changes","AE seriousness changes","AE comparison unavailable.");ax_unavailable(a,"V31","unmatched_serious_events","Unmatched serious events","AE comparison unavailable.")}
  f<-ax_read(out,"ADSYMP/field_changes.csv",FALSE)
  if(file.exists(file.path(out,"ADSYMP/field_changes.csv"))) {
    f<-f[f$field=="AVALC",];fv<-sc_class(f$eua_value);bv<-sc_class(f$bla_value)
    rows<-list();for(same in c(FALSE,TRUE))for(direction in c("YES -> NO","NO -> YES")) {
      sel<-if(direction=="YES -> NO")fv=="POSITIVE_OR_YES"&bv=="NEGATIVE_OR_NO" else fv=="NEGATIVE_OR_NO"&bv=="POSITIVE_OR_YES"
      if(same)sel<-sel&!is.na(f$eua_date)&nzchar(f$eua_date)&f$eua_date==f$bla_date
      z<-f[which(sel),];rows[[length(rows)+1L]]<-data.frame(`Linked-record scope / transition`=paste(if(same)"Same date" else "All linked dates",direction,sep="\n"),
        `mRNA-1273`=ax_num(sum(z$arm=="mRNA-1273")),Placebo=ax_num(sum(z$arm=="Placebo")),Total=ax_num(nrow(z)),check.names=FALSE)
    }
    ax_table(a,"V32","symptom_response_reversals","Changes in recorded symptom responses",do.call(rbind,rows),
      "ADSYMP AVALC | Counts of linked records, not unique participants",
      "Same-date rows are a subset of all linked dates. Reference and provisional clinical links retain their original classification in the CSV crosswalk. Revised dates and response fields require source-level reconciliation.",widths=c(.43,.19,.19,.19))
  }else ax_unavailable(a,"V32","symptom_response_reversals","Symptom response changes","ADSYMP comparison unavailable.")
  ax_finish(a)
}

if(isTRUE(getOption("phmpt.autorun",TRUE)))analyse_eua_bla_changes()
