suppressPackageStartupMessages({library(brglm2);library(data.table);library(jsonlite)})
R<-'/root/projects/linezolid_platelet_recovery'
# Execute only setup/definitions before the bootstrap loop; no endpoint reads.
ee<-parse(file.path(R,'scripts/61_comparability_models_v09.R'))
for(e in ee){if(is.call(e)&&identical(e[[1]],as.name('for'))&&identical(e[[2]],as.name('pop')))break;eval(e)}
rows<-list();checks<-list()
for(pop in cfg$populations){
 ix<-which(ct[[pop]]);r<-raw[ix,,drop=FALSE];r$icu_unit<-factor(r$icu_unit);x<-prep(r);ob<-fitone(x);stopifnot(ob$valid)
 mm<-model.matrix(~icu_unit-1,x);a<-x$A==1
 extra<-cbind(mm,ob$mm[,!grepl('^icu_unit|Intercept',colnames(ob$mm)),drop=FALSE])
 rows[[pop]]<-rbindlist(lapply(colnames(extra),function(k){v<-extra[,k];sd0<-sqrt((var(v[a])+var(v[!a]))/2);data.table(population=pop,term=k,mean_LZD_before=mean(v[a]),mean_VAN_before=mean(v[!a]),mean_LZD_after=weighted.mean(v[a],ob$w[a]),mean_VAN_after=weighted.mean(v[!a],ob$w[!a]),SMD_before=if(sd0>0)(mean(v[a])-mean(v[!a]))/sd0 else NA_real_,SMD_after=if(sd0>0)(weighted.mean(v[a],ob$w[a])-weighted.mean(v[!a],ob$w[!a]))/sd0 else NA_real_)}))
 checks[[pop]]<-list(n_LZD=sum(a),n_VAN=sum(!a),max_all_level_SMD=max(abs(rows[[pop]]$SMD_after),na.rm=TRUE),weight_ranges=list(LZD=range(ob$w[a]),VAN=range(ob$w[!a])),ESS=list(LZD=ess(ob$w[a]),VAN=ess(ob$w[!a])),nonintercept_df=ncol(ob$mm)-1,point_fallback=ob$numerical_fallback,point_warning_count=length(ob$warnings),omitted_columns=ob$omitted)
 cat(toJSON(checks[[pop]],auto_unbox=TRUE),'\n')
}
fwrite(rbindlist(rows),file.path(R,'reports/ALL_LEVEL_BALANCE_v0.9.csv'))
write_json(list(utc=format(Sys.time(),tz='UTC',usetz=TRUE),no_outcomes_loaded=TRUE,all_icu_levels_including_reference=TRUE,populations=checks),file.path(R,'reports/COMPARABILITY_POINT_AUDIT_v0.9.json'),pretty=TRUE,auto_unbox=TRUE,digits=12)
