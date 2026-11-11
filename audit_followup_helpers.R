# Shared base-R utilities for Moderna audit follow-up 1.0.0.
audit_root <- function(root=getOption('phmpt.followup_root',getwd())) {
 if(!file.exists(file.path(root,'inputs','episodes.csv.gz')))stop('Set options(phmpt.followup_root=...) to the extracted bundle folder, or set it as the RStudio working directory.')
 normalizePath(root,winslash='/',mustWork=TRUE)
}
rd <- function(root,name) {p<-file.path(root,'inputs',paste0(name,'.csv.gz'));if(!file.exists(p))stop('Missing pinned input: ',p);read.csv(gzfile(p),check.names=FALSE,stringsAsFactors=FALSE,na.strings='',fileEncoding='UTF-8')}
yes <- function(x) !is.na(x)&toupper(as.character(x))%in%c('TRUE','Y','YES','1')
text0 <- function(x) {x<-as.character(x);x[is.na(x)]<-'';x}
date0 <- function(x) as.Date(substr(text0(x),1,10),format='%Y-%m-%d')
key <- function(...) do.call(paste,c(lapply(list(...),text0),sep='|'))
ck <- function(ok,msg)if(!isTRUE(ok))stop(msg,call.=FALSE)
wr <- function(x,out,name) {dir.create(out,recursive=TRUE,showWarnings=FALSE);write.csv(x,file.path(out,paste0(name,'.csv')),row.names=FALSE,na='')}
pct <- function(x,n) ifelse(n>0,100*x/n,NA_real_)
rate <- function(x,n) sprintf('%s/%s (%.1f%%)',format(x,big.mark=',',trim=TRUE),format(n,big.mark=',',trim=TRUE),pct(x,n))
escape <- function(x) {x<-text0(x);for(p in list(c('&','&amp;'),c('<','&lt;'),c('>','&gt;'),c('"','&quot;')))x<-gsub(p[1],p[2],x,fixed=TRUE);x}
html_table <- function(x) paste0('<table><thead><tr>',paste0('<th>',escape(names(x)),'</th>',collapse=''),'</tr></thead><tbody>',paste(apply(x,1,function(r)paste0('<tr>',paste0('<td>',escape(r),'</td>',collapse=''),'</tr>')),collapse=''),' </tbody></table>')
html_report <- function(path,title,body)writeLines(paste0('<!doctype html><html lang="en"><meta charset="utf-8"><title>',escape(title),'</title><style>body{font:17px/1.55 Arial,sans-serif;max-width:1100px;margin:40px auto;padding:0 24px;color:#182c3b}h1{font-size:28px}h2{font-size:22px;margin-top:32px}table{border-collapse:collapse;width:100%;font-size:14px}th,td{border-bottom:1px solid #ccd6de;padding:9px;text-align:left}th{background:#edf2f6}code{font-size:.9em} .note{background:#fff2d9;padding:12px}img{width:100%}</style><body><h1>',title,'</h1>',paste(body,collapse='\n'),'</body></html>'),path,useBytes=TRUE)
# Site-cluster bootstrap retains all episodes of a participant within a selected site.
cluster_compare <- function(z,outcome,B=2000L) {
 out<-lapply(c('mRNA-1273','Placebo'),function(a){q<-z[z$arm==a,,drop=FALSE];c(n=nrow(q),x=sum(yes(q[[outcome]])))})
 n1<-out[[1]]['n'];n0<-out[[2]]['n'];x1<-out[[1]]['x'];x0<-out[[2]]['x'];ci<-c(NA_real_,NA_real_)
 if(n1>0&&n0>0){sites<-sort(unique(z$site));g<-sapply(sites,function(s){q<-z[z$site==s,,drop=FALSE];c(sum(q$arm=='mRNA-1273'&yes(q[[outcome]])),sum(q$arm=='mRNA-1273'),sum(q$arm=='Placebo'&yes(q[[outcome]])),sum(q$arm=='Placebo'))});if(length(sites)>2){w<-rmultinom(B,length(sites),rep(1/length(sites),length(sites)));g<-g%*%w;boot<-100*(g[1,]/g[2,]-g[3,]/g[4,]);ci<-quantile(boot,c(.025,.975),na.rm=TRUE,names=FALSE)}}
 data.frame(vaccine_x=unname(x1),vaccine_n=unname(n1),vaccine_pct=pct(x1,n1),placebo_x=unname(x0),placebo_n=unname(n0),placebo_pct=pct(x0,n0),gap_pp=pct(x1,n1)-pct(x0,n0),lower95=ci[1],upper95=ci[2],row.names=NULL)
}
boundaries <- function(f) {
 ck(!anyDuplicated(f$subject),'Duplicate follow-up subject');n<-nrow(f);b<-date0(f$follow_end);ub<-date0(f$UBLDDT);any<-yes(f$has_blind_boundary)
 uncertain<-is.na(ub);unresolved<-!any & is.na(date0(f$EOSDT)) & is.na(date0(f$DTHDT));cut<-as.Date(getOption('phmpt.uncertain_boundary_cap','2020-12-17'))
 list(recorded=b,explicit_unblinding_only=replace(b,uncertain,as.Date(NA)),recorded_boundary_only=replace(b,!any,as.Date(NA)),uncertain_capped_17Dec2020=as.Date(ifelse(uncertain,pmin(as.numeric(b),as.numeric(cut)),as.numeric(b)),origin='1970-01-01'),unresolved_end_capped_17Dec2020=as.Date(ifelse(unresolved,pmin(as.numeric(b),as.numeric(cut)),as.numeric(b)),origin='1970-01-01'),all_capped_17Dec2020=pmin(b,cut))
}
# Fixed raster output, no footnotes or references embedded in the figure.
print_table_png <- function(x,file,title,subtitle='',width=2400,row_height=115) {
 h<-max(650,340+(nrow(x)+1)*row_height);png(file,width=width,height=h,res=230,type=if(capabilities('cairo'))'cairo'else getOption('bitmapType'))
 on.exit(dev.off());par(mar=c(1,1,1,1),xpd=NA,family='sans');plot.new();plot.window(c(0,1),c(0,1));text(.01,.96,title,adj=c(0,1),cex=1.3,font=2);text(.01,.88,subtitle,adj=c(0,1),cex=.85,col='#435364')
 top<-.77;bottom<-.04;rh<-(top-bottom)/(nrow(x)+1);widths<-if(ncol(x)==3)c(.52,.24,.24)else rep(1/ncol(x),ncol(x));xs<-cumsum(c(0,widths));cols<-names(x)
 for(i in 0:nrow(x)){y<-top-i*rh;rect(0,y-rh,1,y,col=if(i==0)'#e7eef4'else if(i%%2)'#f5f7f9'else'white',border=NA);vals<-if(i==0)cols else as.character(x[i,]);for(j in seq_along(vals))text(xs[j]+.012,y-rh/2,paste(strwrap(vals[j],width=if(j==1&&ncol(x)==3)53 else 29),collapse='\n'),adj=c(0,.5),cex=.82,font=if(i==0)2 else 1);segments(0,y-rh,1,y-rh,col='#d4dce4')}
}
