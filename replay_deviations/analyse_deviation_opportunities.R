# Moderna COVE, Production S3 / BLA. Base R; no packages or network required.
# Unzip the evidence bundle. Open this script in RStudio, set the working directory
# to the extracted folder and Run All. Alternatively set options(phmpt.m1_dir=...).
# Reads the bundled extracts; does not modify the project's original source files.
analyse_deviation_opportunities <- function(root=getOption("phmpt.m1_dir",getwd())) {
  rd <- function(s) read.csv(file.path(root,"inputs",paste0(s,".csv")),colClasses="character",
                              na.strings=c("","NA"),check.names=FALSE, fileEncoding="UTF-8")
  out <- file.path(root,"recomputed");dir.create(out,showWarnings=FALSE)
  wr <- function(x,n) write.csv(x,file.path(out,paste0(n,".csv")),row.names=FALSE,na="")
  pop <- rd("population"); arm <- c("mRNA-1273","Placebo")
  for(n in c("RANDDT","CUTOFFDT","UBLDDT","AP01EDT","AP02SDT","EOSDT","DTHDT"))pop[[n]]<-as.Date(pop[[n]])
  pop$end<-with(pop,pmin(CUTOFFDT,UBLDDT-1,AP01EDT,AP02SDT-1,EOSDT,DTHDT,na.rm=TRUE))
  stopifnot(!anyDuplicated(pop$subject),all(pop$end>=pop$RANDDT),all(pop$arm%in%arm))
  id <- function(x) match(x,pop$subject)
  observable <- function(id_,d) {j<-id(id_);!is.na(j)&!is.na(d)&d>=pop$RANDDT[j]&d<=pop$end[j]}
  miss<-is.na(pop$UBLDDT); no3<-miss&is.na(pop$AP01EDT)&is.na(pop$AP02SDT)
  caps<-do.call(rbind,lapply(arm,function(a){z<-pop$arm==a;data.frame(arm=a,N=sum(z),
    missing_UBLDDT=sum(z&miss),alternative_period_date=sum(z&miss&!no3),
    all_three_missing=sum(z&no3),admin_if_all_three_missing=sum(z&no3&pop$end==pop$CUTOFFDT),
    missing_UBLDDT_end_at_cutoff_including_ties=sum(z&miss&pop$end==pop$CUTOFFDT))}));wr(caps,"caps")
  dv<-rd("deviations");dv$date<-as.Date(dv$date);dv<-dv[observable(dv$subject,dv$date),]
  dv$arm<-pop$arm[id(dv$subject)]; stopifnot(!anyDuplicated(dv[,c("subject","DVSEQ")]))
  inventory<-do.call(rbind,lapply(split(dv,dv$code),function(g) do.call(rbind,lapply(arm,function(a){
    z<-g[g$arm==a,];j<-unique(id(z$subject));data.frame(code=g$code[1],arm=a,records=nrow(z),
      people=length(j),missing_UBLDDT=sum(miss[j]),all_three_missing=sum(no3[j]))
  }))));wr(inventory,"inventory")
  # Index: earliest usable nonbaseline positive PCR at the initial illness visit;
  # local-coded diagnostics can also qualify 0-3 days after an independent Day 1 visit.
  pos<-rd("usable_nonbaseline_positive_PCR");pos$date<-as.Date(pos$date)
  pos<-pos[observable(pos$subject,pos$date),]
  vis<-rd("independent_illness_visits");vis$visit_date<-as.Date(vis$visit_date)
  ini<-function(x) grepl("ILLNESS VISIT DAY 1(?![0-9])",toupper(x),perl=TRUE)
  loc<-as.logical(pos$local_coded);di<-as.logical(pos$usable_diagnostic)
  near<-function(v){dd<-split(as.integer(v$visit_date),v$subject)
    vapply(seq_len(nrow(pos)),function(i) any((as.integer(pos$date[i])-dd[[pos$subject[i]]]) %in% 0:3),logical(1))}
  initial<-di & (grepl("^ILLNESS VISIT DAY 1(?![0-9])",toupper(pos$raw_visit),perl=TRUE) | (loc & near(vis[ini(vis$labels),])))
  broad<-di & (as.logical(pos$original_illness_label) | (loc & near(vis)))
  first<-function(z) aggregate(date~subject,z,min)
  defs<-list(first_initial_illness_positive_PCR=first(pos[initial,]),
             first_illness_linked_positive_PCR=first(pos[broad,]),
             first_nonbaseline_positive_PCR=first(pos))
  horizons<-c("P08.12"=28L,"P07.08"=28L,"P08.14"=14L)
  tab<-list()
  for(def in names(defs)) for(code in names(horizons)) for(full in c(FALSE,TRUE)) {
    h<-horizons[[code]];v<-defs[[def]]
    v<-v[v$date + (if(full) h else 0) <= pop$end[id(v$subject)],,drop=FALSE]
    g<-dv[dv$code==code,];lag<-as.numeric(g$date-v$date[match(g$subject,v$subject)])
    ok<-!is.na(lag)&lag>=0&lag<=h
    for(a in arm){den<-sum(pop$arm[id(v$subject)]==a);n<-length(unique(g$subject[ok&g$arm==a]))
      tab[[length(tab)+1]]<-data.frame(positive_index=def,code=code,
        followup=if(full) "complete_minimum_followup" else "any_followup",
        days=h,arm=a,people_with_opportunity=den,affected_in_fixed_window=n,percent=100*n/den)}
  };templ<-do.call(rbind,tab);wr(templ,"opportunities")
  # The deviation record has no verified AE identifier. A temporal match is only a candidate.
  ae<-rd("adverse_events");ae$start<-as.Date(ae$start);ae<-ae[observable(ae$subject,ae$start),]
  ae$conf<-ae$AECVDIAG%in%"Yes";pt<-toupper(ae$AEDECOD)
  ae$cov<-ae$conf | pt %in% c("COVID-19","COVID-19 PNEUMONIA","ASYMPTOMATIC COVID-19","SARS-COV-2 TEST POSITIVE")
  ae$susp<-pt%in%c("SUSPECTED COVID-19","CORONAVIRUS INFECTION");ae$sae<-ae$AESER%in%"Y"
  ae$hosp<-ae$AESHOSP%in%"Y" | !is.na(ae$AEHOSPST) | !is.na(ae$AEHOSPEN)
  late<-dv[dv$code=="I12.01",]
  levels_<-c("COVID_AE_WITH_HOSPITAL_EVIDENCE","COVID_AE_WITHOUT_OWN_HOSPITAL_EVIDENCE","OTHER_SERIOUS_AE",
             "SUSPECTED_OR_NONSPECIFIC_CORONAVIRUS_ONLY","NO_COVID_OR_SERIOUS_CANDIDATE")
  e<-merge(late[,c("subject","arm","DVSEQ","date")],ae,by="subject",all.x=TRUE)
  e$distance<-abs(as.numeric(e$start-e$date))
  ed<-list();sc<-list()
  for(w in c(0,7,14,28)){
    t<-e[!is.na(e$distance)&e$distance<=w,,drop=FALSE]
    people<-unique(late[,c("subject","arm")]); b<-rep(5L,nrow(people))
    priority<-ifelse(t$cov&t$hosp,1L,ifelse(t$cov,2L,ifelse(t$sae,3L,ifelse(t$susp,4L,5L))))
    if(nrow(t)) {k<-tapply(priority,t$subject,min);ii<-match(names(k),people$subject);b[ii]<-k}
    by<-as.data.frame(table(arm=factor(people$arm,levels=arm),category=factor(levels_[b],levels=levels_)))
    by$window<-w;ed[[length(ed)+1]]<-by
    for(a in arm){sa<-ae$subject[ae$conf & pop$arm[id(ae$subject)]==a];cn<-unique(t$subject[t$arm==a & t$conf])
      sc[[length(sc)+1]]<-data.frame(arm=a,days=w,people_with_confirmed_symptomatic_AE=length(unique(sa)),
        people_with_nearby_deviation=length(cn))}
  }
  wr(do.call(rbind,ed),"late_reporting_context");wr(do.call(rbind,sc),"confirmed_AE_opportunities")
  message("Recomputed results: ",normalizePath(out));invisible(templ)
}
if(isTRUE(getOption("phmpt.autorun",TRUE))) analyse_deviation_opportunities()
