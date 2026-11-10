replay_566_review <- function(root,dest) {
 f<-function(n)read.csv(file.path(root,"evidence",n),check.names=FALSE,colClasses="character",na.strings=character(),fileEncoding="UTF-8")
 d<-f("form_dispositions_566.csv");x<-f("ae_form_review/comparison_AE_rows_287_people.csv")
 k<-function(g)paste(g$SUBJID,g$AESPID,sep=":")
 a<-subset(x,snapshot=="S5");b<-subset(x,snapshot=="S3");c<-subset(x,snapshot=="S27")
 stopifnot(!anyDuplicated(k(a)),!anyDuplicated(k(b)),!anyDuplicated(k(c)),setequal(k(b),k(c)))
 num<-function(y)as.integer(sub("^AE-","",y))
 lo<-tapply(num(b$AESPID),b$SUBJID,min);hi<-tapply(num(b$AESPID),b$SUBJID,max)
 test<-!(k(a)%in%k(b)) & num(a$AESPID)>lo[a$SUBJID] & num(a$AESPID)<hi[a$SUBJID]
 g<-a[which(test),];stopifnot(nrow(g)==566L,length(unique(g$SUBJID))==287L,setequal(k(g),d$key),!anyDuplicated(d$key))
 j<-match(d$key,k(g));stopifnot(identical(d$EUA_serious,g$AESER[j]),sum(d$EUA_serious=="Y")==12,all(nzchar(d$classification)),all(nzchar(d$reason)))
 res<-as.data.frame.matrix(xtabs(~classification+randomized_arm,data=d));res$forms<-rowSums(res)
 res$serious<-as.integer(tapply(d$EUA_serious=="Y",d$classification,sum)[rownames(res)])
 old<-f("ae_form_review/summary.csv");ix<-match(rownames(res),old$classification)
 stopifnot(all(res$forms==as.integer(old$forms[ix])),all(res$serious==as.integer(old$serious[ix])))
 print(res);print(table(released_subject_CRF=d$source_coverage!=""))
 write.csv(data.frame(classification=rownames(res),res,check.names=FALSE),file.path(dest,"replayed_disposition_counts.csv"),row.names=FALSE)
 invisible(res)
}
