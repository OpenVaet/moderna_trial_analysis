crf_summary <- function(root,dest) {
 read <- function(n) read.csv(file.path(root,"evidence",paste0(n,".csv")),colClasses="character",check.names=FALSE, na.strings="",fileEncoding="UTF-8")
 uniq <- function(x) length(unique(na.omit(x))); yes <- function(x) toupper(x)=="TRUE"
 c <- read("corpus_manifest"); f<-c[c$class=="full",]; a<-read("arm_checkbox"); q<-read("dose2_flags"); ae<-read("AE_term_onset_checks"); p<-read("pcr_comparison"); d<-read("contact_lags"); h<-read("history_leads")
 du <- subset(a,!is.na(dose2)); dd<-du[!duplicated(du$sid),]
 conflicts<-subset(dd,dose2 %in% c("Placebo","mRNA-1273") & DOSE2FL=="N" & EX2=="Not Administered")
 assigned<-subset(a,!is.na(assigned)); lg<-subset(d,yes(eligible)); posneg<-c("POSITIVE","NEGATIVE")
 stopifnot(setequal(q$sid,conflicts$sid),!anyDuplicated(q$sid),all(yes(ae$present)),all(yes(ae$term_equal_normalized)),all(yes(ae$start_equal)))
 n<-c(files_total=nrow(c),files_read=nrow(f),subjects=uniq(f$subject),sites=uniq(f$site),pages=sum(as.numeric(f$pages)),placeholders=sum(c$class=="placeholder"),
   assignment_people=uniq(assigned$sid),assignment_discordant=sum(assigned$assigned!=assigned$rand),
   dose2_people=nrow(dd),dose2_placebo=sum(dd$rand=="Placebo"),dose2_mRNA=sum(dd$rand=="mRNA-1273"),
   no_dose2_placebo=sum(dd$rand=="Placebo" & dd$DOSE2FL=="N"),no_dose2_mRNA=sum(dd$rand=="mRNA-1273" & dd$DOSE2FL=="N"),
   dose2_conflicts=nrow(q),dose2_placebo_conflicts=sum(q$rand=="Placebo"),dose2_mRNA_conflicts=sum(q$rand=="mRNA-1273"),dose2_sites=uniq(substr(q$sid,1,5)),dose2_visit_said_no=sum(q$exposure_crf_response=="No (N)"),
   dose1_conflict_people=uniq(subset(a,!is.na(dose1) & dose1!=EX1)$sid),contact_fields=nrow(lg),late30_fields=sum(as.numeric(lg$entry_lag)>=30),late30_people=uniq(lg$sid[as.numeric(lg$entry_lag)>=30]),
   EUA_AE_forms=sum(ae$snapshot=="S5"),BLA_AE_forms_S3=sum(ae$snapshot=="S3"),BLA_AE_forms_S27=sum(ae$snapshot=="S27"),
   PCR_comparisons=nrow(p),PCR_agreements=sum(p$CRF_result==p$same_date_MB),PCR_opposite=sum(p$same_date_MB %in% posneg & p$CRF_result!=p$same_date_MB),PCR_mixed=sum(p$same_date_MB=="MIXED"),PCR_no_exact_date=sum(p$same_date_MB=="NO SAME-DATE BINARY"),
   assignment_histories_both_products=uniq(subset(h,field=="assigned")$sid),PCR_positive_histories_yes_no=uniq(subset(h,field=="positive")$sid),seriousness_histories_yes_no=uniq(subset(h,field=="sae")$sid))
 result<-data.frame(measure=names(n),value=as.numeric(n)); write.csv(result,file=file.path(dest,"publication_counts.csv"),row.names=FALSE);print(result,row.names=FALSE);invisible(result)
}
