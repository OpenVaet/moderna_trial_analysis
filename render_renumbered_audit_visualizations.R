# Re-render V20-V49 from the supplied, already calculated audit tables.
# Base R only. Set this extracted folder as the working directory and Run all.
# This preserves the published data; the analyse_*.R files remain the full generators.
# Optional: options(phmpt.figure_root='D:/Projects/moderna_trial_analysis')

render_renumbered_audit_visualizations <- function(root=getOption('phmpt.figure_root',getwd())) {
 root<-normalizePath(root,winslash='/',mustWork=TRUE)
 old<-options(phmpt.autorun=FALSE,phmpt.followup_autorun=FALSE,
   phmpt.figure_autorun=FALSE,phmpt.root=root,phmpt.followup_root=root,
   phmpt.audit_asset_dir=file.path(root,'audit_assets'))
 on.exit(options(old),add=TRUE)
 assets<-file.path(root,'audit_assets');dir.create(assets,recursive=TRUE,showWarnings=FALSE)
 load_script<-function(name){e<-new.env(parent=globalenv());sys.source(file.path(root,name),envir=e);e}
 read_asset<-function(stem)read.csv(file.path(assets,paste0(stem,'.csv')),
   check.names=FALSE,stringsAsFactors=FALSE,fileEncoding='UTF-8-BOM')
 calls<-function(x,head){
   if(!is.call(x)&&!is.expression(x)&&!is.pairlist(x))return(list())
   hit<-if(is.call(x)&&identical(x[[1]],as.name(head)))list(x)else list()
   children<-as.list(x);if(is.call(x))children<-children[-1]
   c(hit,unlist(lapply(children,calls,head=head),recursive=FALSE))
 }
 assignment<-function(fun,name){
   x<-Filter(function(c)length(c)==3L&&identical(c[[2]],as.name(name)),calls(body(fun),'<-'))
   if(length(x)!=1L)stop('Expected one native drawing definition: ',name)
   x[[1]]
 }
 literal_call<-function(fun,head,id){
   pos<-if(head=='chart')2L else 3L
   x<-Filter(function(c)length(c)>=pos&&is.character(c[[pos]])&&identical(c[[pos]],id),calls(body(fun),head))
   if(length(x)!=1L)stop('Expected one native ',head,' call for ',id)
   x[[1]]
 }
 work<-tempfile('figure_render_');dir.create(work);on.exit(unlink(work,recursive=TRUE),add=TRUE)
 site<-load_script('analyse_site_testing_imbalance.R')
 ctx<-new.env(parent=site);ctx$window<-'BLINDED';ctx$cut<-10;ctx$horizon<-7
 ctx$high_label<-'High positivity >=10% | Low <10%'
 ctx$a<-site$ax_begin(root,work,'V22_V31_episodes')
 # Reuse the native table metadata and drawing functions; no raster relabelling.
 eval(assignment(site$phmpt_export_audit_assets,'specs'),ctx)
 table_loop<-Filter(function(c)identical(c[[3]],as.name('id')),
                    calls(body(site$phmpt_export_audit_assets),'ax_table'))[[1]]
 for(id in names(ctx$specs)){
   ctx$id<-id;ctx$s<-ctx$specs[[id]];ctx$d<-read_asset(paste0(id,'_data'))
   eval(assignment(site$phmpt_export_audit_assets,'note'),ctx)
   eval(table_loop,ctx)
 }
 ctx$s<-'US352';ctx$d<-read_asset('V28_data')
 table<-literal_call(site$phmpt_export_audit_assets,'ax_table','V28');table[[6]]<-quote(d)
 eval(table,ctx)
 # Native two-panel funnel, using its stored point, limit and LOESS coordinates.
 funnel<-new.env(parent=site);funnel$a<-ctx$a;funnel$window<-'BLINDED';funnel$cut<-10;funnel$horizon<-7
 tab<-list(sites=read_asset('V24_funnel_sites'),summary=read_asset('V24_funnel_summary'),
           curve=read_asset('V24_funnel_limits'),loess=read_asset('V24_funnel_loess'))
 funnel$panels<-lapply(c('randomized','eligible'),function(p)lapply(tab,function(d)d[d$population==p,,drop=FALSE]))
 funnel$m<-read_asset('V24_data')
 eval(assignment(site$v17_plot,'draw'),funnel)
 eval(literal_call(site$v17_plot,'ax_png','V24'),funnel)
 # Four native site illustrations.
 eval(assignment(site$phmpt_export_audit_assets,'chart'),ctx)
 eval(assignment(site$phmpt_export_audit_assets,'bars'),ctx)
 ctx$col<-c('#1E667D','#A65B22')
 for(spec in list(c('V26','US372'),c('V27','US325'),c('V29','US372'),c('V31','US302'))){
   ctx$z<-read_asset(paste0(spec[1],'_data'));ctx$s<-spec[2]
   ctx$armcount<-function(ar)if(ar=='mRNA-1273')'15/19'else'26/26'
   eval(literal_call(site$phmpt_export_audit_assets,'chart',spec[1]),ctx)
 }
 site$ax_finish(ctx$a)

 # Tables retain native titles, subtitles, notes, column widths and highlights.
 table_group<-function(script,fun,ids,owner,extra=list(),head='ax_table'){
   e<-load_script(script);env<-new.env(parent=e);env$a<-e$ax_begin(root,work,owner)
   env$ax_table<-e$ax_table
   for(n in names(extra))assign(n,extra[[n]],env)
   for(id in ids){
     env$d<-read_asset(paste0(id,'_data'))
     c<-literal_call(e[[fun]],head,id);c[[1]]<-as.name('ax_table');c[[6]]<-quote(d)
     if('highlight'%in%names(c))c$highlight<-which(grepl('US3952070',env$d[[1]],fixed=TRUE))
     eval(c,env)
   }
   e$ax_finish(env$a)
 }
 table_group('analyse_anti_n_seroconversion.R','sc_export_seroconversion_assets','V32','V32_seroconversion',
   list(b=list(contrast=list(
     fisher_p=fisher.test(matrix(c(21,31,605,43),nrow=2,byrow=TRUE))$p.value,
     difference_pp=100*(21/52-605/648)))))
 table_group('analyse_day57_anti_n_sensitivity.R','d57_export_audit_assets',sprintf('V%02d',33:35),'V33_V35_day57')
 table_group('analyse_eua_bla_changes.R','sc_export_snapshot_assets',sprintf('V%02d',36:40),'V36_V40_snapshots')

 eff<-load_script('analyse_efficacy_cases.R');e<-new.env(parent=eff)
 e$a<-eff$ax_begin(root,work,'V42_V44_efficacy');e$valid<-TRUE
 for(id in c('V42','V44')){
   e$d<-read_asset(paste0(id,'_data'));c<-literal_call(eff$ec_export_audit_assets,'ax_table',id)
   c[[6]]<-quote(d);eval(c,e)
 }
 e$w<-read_asset('V43_data');e$w$week_start<-as.Date(e$w$week_start)
 e$weeks<-sort(unique(e$w$week_start));e$arms<-c('mRNA-1273','Placebo')
 e$vals<-sapply(e$arms,function(ar){q<-e$w[e$w$arm==ar,];q$cases[match(e$weeks,q$week_start)]})
 e$cols<-c('#1E667D','#A65B22');e$markers<-read_asset('V43_calendar_markers');e$markers$date<-as.Date(e$markers$date)
 e$height<-8.5;e$z<-list(valid=TRUE)
 eval(assignment(eff$ec_export_audit_assets,'draw'),e)
 eval(literal_call(eff$ec_export_audit_assets,'ax_png','V43'),e);eff$ax_finish(e$a)
 death<-load_script('analyse_death_reporting.R');de<-new.env(parent=death)
 de$a<-death$dr36_begin(root,work,'V45_deaths');de$visual<-read_asset('V45_data');de$status<-'RECONCILED'
 dc<-literal_call(death$dr36_export,'dr36_table','V45');dc[[6]]<-quote(visual)
 eval(dc,de);death$dr36_finish(de$a)

 quick<-load_script('render_V20_V21_V41.R');quick$render_V20_V21_V41(root)
 numbering<-load_script('analyse_subject_numbering.R')
 numbering$analyse_subject_numbering(root,
   index_file=file.path(root,'figure_inputs','P301_subject_screening_index.csv'),
   audit_asset_dir=assets,dpi=300L)
 message('Regenerated V20-V49 from the supplied audit tables and screening index.')
 invisible(assets)
}
if(isTRUE(getOption('phmpt.figure_autorun',TRUE)))render_renumbered_audit_visualizations()
