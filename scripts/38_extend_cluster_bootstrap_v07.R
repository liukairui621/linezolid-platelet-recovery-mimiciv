# Same fitted point weights; extend the original frozen sequence to 2000 draws.
suppressPackageStartupMessages({library(brglm2);library(data.table);library(jsonlite)})
R<-'/root/projects/linezolid_platelet_recovery'
p<-fromJSON(file.path(R,'config/analysis_plan_v0.7.json'))$phase7
m<-readRDS(file.path(R,'cache/phase61_overlap.rds'));raw<-m$raw;x<-m$data
old<-list(data=x);f<-formula(m$fit)
# Reuse the exact frozen preparation and fit function bodies without executing the old pipeline.
expr<-parse(file.path(R,'scripts/31b_bias_reduced_overlap_v061.R'))
for(e in expr)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])&&as.character(e[[2]])%in%c('prep','fitone'))eval(e)
stopifnot(exists('prep'),exists('fitone'),nrow(raw)==1569,length(m$bootstrap)==500)
cl<-split(seq_len(nrow(raw)),raw$subject_id)
pat<-vapply(cl,function(ii)paste(sort(unique(x$A[ii])),collapse=','),'');strata<-split(names(cl),pat)
set.seed(p$bootstrap$seed);B<-p$bootstrap$total_B
boot<-vector('list',B);diag<-vector('list',B)
for(i in seq_len(B)){
 namesdraw<-unlist(lapply(strata,function(z)sample(z,length(z),replace=TRUE)),use.names=FALSE)
 ix<-unlist(cl[namesdraw],use.names=FALSE)
 if(i<=500){
  stopifnot(identical(ix,m$bootstrap[[i]]$ix));boot[[i]]<-m$bootstrap[[i]]
  diag[[i]]<-data.table(replicate=i,valid=TRUE,warnings=0L,reused=TRUE)
 }else{
  z<-prep(raw[ix,,drop=FALSE]);ob<-tryCatch(fitone(z),error=function(e)list(valid=FALSE,warnings=conditionMessage(e)))
  if(ob$valid)boot[[i]]<-list(ix=ix,weights=ob$w)
  diag[[i]]<-data.table(replicate=i,valid=ob$valid,warnings=length(ob$warnings),reused=FALSE)
 }
 if(i%%250==0){cat(i,'/',B,'draws\n');flush.console()}
}
dg<-rbindlist(diag);stopifnot(mean(dg$valid)>=.98)
saveRDS(list(bootstrap=boot,diagnostics=dg),file.path(R,'cache/phase7_bootstrap.rds'))
report<-list(completed_utc=format(Sys.time(),tz='UTC',usetz=TRUE),requested=B,valid=sum(dg$valid),failed=sum(!dg$valid),draws_with_warnings=sum(dg$warnings>0),
 original500_index_sequences_identical=TRUE,original500_weights_reused=TRUE,additional_refits=B-500,seed=p$bootstrap$seed,
 point_weights_unchanged=TRUE,median_fills_recomputed_in_new_resamples=TRUE)
write_json(report,file.path(R,'reports/BOOTSTRAP_EXTENSION_v0.7.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(R,'reports/R_sessionInfo_v0.7.txt'))
cat(toJSON(report,auto_unbox=TRUE),'\n')
