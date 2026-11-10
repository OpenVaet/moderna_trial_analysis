# Moderna P301: standalone audit of AE form numbers and PCR row-number gaps.
# Run in the project root in RStudio. Reads original *_S_ae and *_S_mb files only.
# Defaults to XPT (haven or foreign), otherwise CSV. Never pools submissions.
# Optional exact overrides, e.g. options(phmpt.seq.files=c(S5_AE="...", S3_MB="...")).
# Optional extra domains PR/CM/MH/HO from CSV or XPT if available. JMP-only sources
# must first be exported using the project's conversion tools. No data download.
analyse_record_number_gaps <- function(root=getOption("phmpt.root",getwd()),
  overrides=getOption("phmpt.seq.files",character()), dest=getOption("phmpt.seq.out",
  file.path(root,"analysis",paste0("record_number_gaps_",format(Sys.time(),"%Y%m%d_%H%M%S"))))) {
  if(!dir.create(dest,recursive=TRUE,showWarnings=FALSE)&&!dir.exists(dest)) stop("Output directory unavailable")
  paths<-list.files(root,pattern="\\.(csv|xpt)$",recursive=TRUE,full.names=TRUE,ignore.case=TRUE)
  paths<-paths[grepl("mrna-1273-p301_s_",basename(paths),ignore.case=TRUE)]
  statuses<-list(); prof<-list(); gaps<-list(); retained<-list(); analyses<-list()
  write_out<-function(x,n)write.csv(x,file.path(dest,n),row.names=FALSE,na="")
  pick<-function(sn,dm){
    key<-paste(sn,dm,sep="_"); if(key%in%names(overrides)){f<-unname(overrides[key]);if(!file.exists(f))stop("Override not found: ",f);return(f)}
    pat<-paste0("_",sn,"_M5_mrna-1273-p301_S_",dm,"\\.(csv|xpt)$")
    f<-paths[grepl(pat,basename(paths),ignore.case=TRUE)]
    if(!length(f))return(NULL)
    typ<-if(any(grepl("\\.xpt$",f,ignore.case=TRUE))&&
       (requireNamespace("haven",quietly=TRUE)||requireNamespace("foreign",quietly=TRUE)||!any(grepl("\\.csv$",f,ignore.case=TRUE)))) "xpt" else "csv"
    f<-sort(f[grepl(paste0("\\.",typ,"$"),f,ignore.case=TRUE)])
    if(length(f)>1 && length(unique(unname(tools::md5sum(f))))!=1L)stop("Conflicting inputs: choose the exact ",key," using phmpt.seq.files")
    f[1L]
  }
  read1<-function(f,want){
    if(grepl("\\.xpt$",f,ignore.case=TRUE)){
      if(requireNamespace("haven",quietly=TRUE)) x<-as.data.frame(haven::read_xpt(f))
      else if(requireNamespace("foreign",quietly=TRUE)) x<-foreign::read.xport(f)
      else stop("Install 'haven' to read XPT, or specify your CSV export.")
    } else {
      n<-names(read.csv(f,nrows=0L,check.names=FALSE,fileEncoding="UTF-8-BOM"))
      x<-read.csv(f,colClasses=ifelse(toupper(n)%in%want,"character","NULL"),
                  check.names=FALSE,fileEncoding="UTF-8-BOM",na.strings="",comment.char="")
    }
    names(x)<-toupper(names(x));x<-x[intersect(want,names(x))]
    if("STUDYID"%in%names(x)) x<-x[!is.na(x$STUDYID)&trimws(x$STUDYID)=="mRNA-1273-P301",,drop=FALSE]
    if(!"USUBJID"%in%names(x))stop("No USUBJID in ",f)
    bad<-is.na(x$USUBJID)|!grepl("^mRNA-1273-P301-US[0-9]{3}-[0-9]{4}$",x$USUBJID,ignore.case=TRUE)
    if(any(bad))stop("Missing or unexpected participant key in ",f)
    x
  }
  summary1<-function(sn,dm,field,scope,sub,num,keep=FALSE){
    num<-suppressWarnings(as.numeric(num));ok<-is.finite(num)&num>=1&num==floor(num);sub<-as.character(sub)
    v<-unique(data.frame(USUBJID=sub[ok],number=num[ok]));v<-v[order(v$USUBJID,v$number),,drop=FALSE]
    if(!nrow(v))return(invisible(NULL))
    b<-c(TRUE,v$USUBJID[-1L]!=head(v$USUBJID,-1L));last<-c(b[-1L],TRUE)
    span<-v$number[last]-v$number[b]+1; sz<-diff(c(which(b),nrow(v)+1L));h<-span-sz
    key<-paste(sn,dm,scope);prof[[key]]<<-data.frame(snapshot=sn,domain=dm,field=field,scope=scope,rows=length(num),invalid_or_blank=sum(!ok),duplicate_keys=sum(ok)-nrow(v),subjects=sum(b),subjects_with_two_or_more=sum(sz>=2),contiguous_from_one=sum(v$number[b]==1&h==0),subjects_with_internal_gaps=sum(h>0),internal_positions=sum(h),subjects_start_above_one=sum(v$number[b]>1),assumed_leading_positions=sum(v$number[b]-1),max_number=max(v$number))
    if(keep){
      before<-c(0,head(v$number,-1L));before[b]<-0; idx<-which(v$number-before>1)
      if(length(idx))gaps[[key]]<<-do.call(rbind,lapply(idx,function(j)data.frame(snapshot=sn,domain=dm,USUBJID=v$USUBJID[j],number=seq.int(before[j]+1,v$number[j]-1),position=if(b[j])"leading" else "internal",previous=if(b[j])NA_integer_ else before[j],following=v$number[j])))
      retained[[paste(sn,dm,sep="_")]]<<-v
    };invisible(NULL)
  }
  fields<-list(AE=c("AETERM","AESTDTC","AESER","AEGRPID"),MB=c("MBTESTCD","MBTEST","MBMETHOD","MBREFID","MBDTC"),HO="HOSPID",PR="PRCAT",MH="MHCAT",CM="CMCAT")
  for(sn in c("S5","S3","S27"))for(dm in if(sn=="S27")"AE" else names(fields)){
    f<-pick(sn,dm);key<-paste(sn,dm,sep="_");req<-dm%in%c("AE","MB");seqk<-paste0(dm,"SEQ");sp<-paste0(dm,"SPID")
    if(is.null(f)){if(req)stop("Missing ",key,"; set phmpt.seq.files");statuses[[key]]<-data.frame(key,path="",status="optional input not supplied");next}
    x<-read1(f,unique(c("STUDYID","USUBJID",seqk,sp,fields[[dm]])))
    need<-switch(dm,AE=c(seqk,sp,"AETERM","AESTDTC","AESER"),MB=c(seqk,"MBTESTCD","MBTEST","MBMETHOD"),HO=c(seqk,sp),c(seqk,sp,fields[[dm]]))
    miss<-setdiff(need,names(x));if(length(miss)){if(req)stop("Missing fields in ",key,": ",paste(miss,collapse=",")); statuses[[key]]<-data.frame(key,path=f,status=paste("fields unavailable",paste(miss,collapse=",")));next}
    statuses[[key]]<-data.frame(key,path=f,status="read");message("Checking ",key," (",nrow(x)," rows)")
    summary1(sn,dm,seqk,"entire domain",x$USUBJID,x[[seqk]])
    if(dm=="AE"){
      tok<-trimws(as.character(x$AESPID));n<-ifelse(grepl("^AE-[0-9]+$",tok),sub("^AE-","",tok),NA_character_)
      if(anyNA(n)||anyDuplicated(data.frame(x$USUBJID,n)))stop("Unexpected non-unique AE form key; inspect before linkage")
      x$number<-as.integer(n);summary1(sn,dm,sp,"source: all numbered AE forms",x$USUBJID,n,TRUE);analyses[[key]]<-x
    }else if(dm=="MB"){
      method<-toupper(trimws(x$MBMETHOD));m<-!is.na(x$MBTESTCD)&toupper(trimws(x$MBTESTCD))=="SARSCOV2"&!is.na(method)&grepl("PCR|POLYMERASE",method)&!grepl("ANTIGEN|ANTIBODY|SEROLOG",method)
      ss<-as.character(x$USUBJID);n<-as.numeric(x$MBSEQ);if(anyDuplicated(data.frame(ss,n)))stop("Duplicate MB keys")
      summary1(sn,dm,seqk,"SARS-CoV-2 RT-PCR rows (not necessarily completed)",ss[m],n[m])
      lo<-tapply(n[m],ss[m],min);hi<-tapply(n[m],ss[m],max);ix<-!m&!is.na(lo[ss])&n>lo[ss]&n<hi[ss];vv<-as.data.frame(table(x$MBTESTCD[ix],x$MBTEST[ix],useNA="ifany"));names(vv)<-c("MBTESTCD","MBTEST","positions");write_out(vv[vv$positions>0,],paste0(sn,"_PCR_intervening_record_types.csv"))
    }else if(dm%in%c("CM","MH","PR")){
      catval<-switch(dm,MH="GENERAL",PR="CONCOMITANT PROCEDURES",CM="all categories");if(dm!="CM")x<-x[which(x[[paste0(dm,"CAT")]]==catval),,drop=FALSE]
      summary1(sn,dm,sp,paste0("source: ",catval),x$USUBJID,x[[sp]],TRUE)
    }else if(dm=="HO"){analyses[[key]]<-x}
  }
  kfun<-function(a)paste(a$USUBJID,a$number,sep="|")
  g<-do.call(rbind,gaps)
  if(NROW(g)>0L){g$in_EUA_numbered_forms<-NA;for(d in unique(g$domain)){baseline<-retained[[paste0("S5_",d)]];i<-g$domain==d;if(!is.null(baseline))g$in_EUA_numbered_forms[i]<-kfun(g[i,])%in%kfun(baseline)};write_out(g,"source_form_gaps.csv")}
  a<-analyses$S5_AE;b<-analyses$S3_AE;c<-analyses$S27_AE
  l<-a[!kfun(a)%in%kfun(b),];lo<-tapply(b$number,b$USUBJID,min);hi<-tapply(b$number,b$USUBJID,max)
  l$location<-ifelse(is.na(lo[l$USUBJID]),"no numbered form for subject",ifelse(l$number<lo[l$USUBJID],"leading",ifelse(l$number>hi[l$USUBJID],"trailing","internal")))
  write_out(l,"EUA_AE_form_numbers_absent_in_S3.csv")
  j<-match(kfun(b),kfun(c));write_out(data.frame(S3_rows=nrow(b),S27_rows=nrow(c),shared=sum(!is.na(j)),AESEQ_changed=sum(b$AESEQ!=c$AESEQ[j],na.rm=TRUE)),"AE_renumbering.csv")
  for(sn in c("S5","S3")){
    t<-analyses[[paste0(sn,"_HO")]];if(is.null(t))next
    h<-ifelse(is.na(t$HOSPID),"",toupper(t$HOSPID));hits<-regmatches(h,gregexpr("\\bAE[- ]0*[0-9]+\\b",h,perl=TRUE));lens<-lengths(hits);u<-data.frame(USUBJID=rep(t$USUBJID,lens),HOSEQ=rep(t$HOSEQ,lens),number=as.integer(sub("^AE[- ]","",unlist(hits))),stringsAsFactors=FALSE)
    u$present<-kfun(u)%in%kfun(analyses[[paste0(sn,"_AE")]]);write_out(u,paste0(sn,"_HO_AE_references.csv"))
  }
  write_out(do.call(rbind,prof),"sequence_profiles.csv");write_out(do.call(rbind,statuses),"input_files.csv");message("Wrote ",normalizePath(dest));invisible(dest)
}
if(isTRUE(getOption("phmpt.autorun",TRUE)))analyse_record_number_gaps()
