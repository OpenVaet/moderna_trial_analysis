# Run this whole file in RStudio with the extracted bundle as your working directory.
# Optional: options(phmpt.cm_root = "D:/.../extracted_bundle"). No packages needed.
# Reconstructs the denominator from S5 numbered keys, then checks and summarizes the
# attached review ledgers. It does not try to automate interpretation of safety narratives.
replay_CMSPID_review <- function(root = getOption("phmpt.cm_root", getwd())) {
  rd <- function(name) read.csv(file.path(root,"evidence",name), colClasses="character", na.strings=character(), check.names=FALSE)
  k <- rd("S5_numbered_keys.csv"); x <- rd("all_1170_positions.csv")
  id <- function(s,n) paste(s,as.integer(n),sep=":")
  stopifnot(!anyDuplicated(id(k$USUBJID,k$CMSPID)), !anyDuplicated(id(x$USUBJID,x$CMSPID)))
  pos <- lapply(split(as.integer(k$CMSPID),k$USUBJID),function(v) setdiff(seq.int(min(v),max(v)),v))
  pos <- pos[lengths(pos)>0]
  ks <- id(rep(names(pos),lengths(pos)),unlist(pos,use.names=FALSE))
  stopifnot(setequal(ks,id(x$USUBJID,x$CMSPID)))
  early <- rd("107_S80_original_row_fields.csv"); cv <- rd("CRF_comparison_by_export.csv")
  key <- function(s,n) paste(s,as.integer(n),sep=":")
  hit <- key(x$subject,x$CMSPID) %in% key(cv$subject,cv$gap)
  exact <- key(x$subject,x$CMSPID) %in% key(cv$subject[cv$exact_form!="no numbered form"],cv$gap[cv$exact_form!="no numbered form"])
  named <- key(x$subject,x$CMSPID) %in% key(cv$subject[cv$exact_form=="named"],cv$gap[cv$exact_form=="named"])
  idx <- match(id(x$USUBJID,x$CMSPID),id(early$USUBJID,early$CMSPID)); prior <- !is.na(idx); term <- ifelse(prior,early$CMTRT[idx],"")
  stopifnot(all(hit==(x$available_CRF!="")), all(exact==(tolower(x$exact_form_seen)=="true")))
  # Matched key and non-blank drug name in the earlier submission.
  byclass <- ifelse(named, "Named CRF", ifelse(exact, "Blank CRF only", ifelse(prior,ifelse(term!="","Named S80","Blank S80"), "No direct record")))
  print(data.frame(S5_numbered_subjects=length(unique(k$USUBJID)),gap_subjects=length(pos),gaps=length(ks),CRF_covered=sum(hit)))
  print(table(byclass)); print(table(byclass[hit])); print(table(x$TRT01P))
  od <- file.path(root,"replay"); dir.create(od,showWarnings=FALSE,recursive=TRUE)
  write.csv(data.frame(category=names(table(byclass)), positions=as.integer(table(byclass))),file.path(od,"counts.csv"),row.names=FALSE)
  write.csv(x[prior & term!="" & x$S5_level=="0" & x$S3_level=="0",],file.path(od,"no_later_product_match.csv"),row.names=FALSE)
  invisible(x)
}
if (isTRUE(getOption("phmpt.autorun",TRUE))) replay_CMSPID_review()
