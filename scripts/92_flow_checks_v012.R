suppressPackageStartupMessages({library(data.table);library(jsonlite)})
R<-'/root/projects/linezolid_platelet_recovery';O<-file.path(R,'reports/submission_v0.12');d<-fread(file.path(R,'cache/phase10_baseline.csv'))
# Explicit sequential flow; each filter is cumulative. Keep patient records on server.
stages<-list(known_baseline_low=rep(TRUE,nrow(d)))
stages$physical_location<-stages[[1]]&d$physical_verified
stages$no_known_vre<-stages$physical_location&!d$known_vre7d
stages$routine_culture<-stages$no_known_vre&d$routine_bacterial_culture72h
stages$shared_units<-d$routine_bacterial_culture72h_shared
stages$primary_no_recorded_lmwh<-stages$shared_units&!d$lmwh_available7d
stopifnot(all(stages$shared_units<=stages$routine_culture),all(stages$primary_no_recorded_lmwh<=stages$shared_units))
out<-list();prev<-NULL
for(k in names(stages)){
 ix<-stages[[k]]
 for(dr in c('linezolid','vancomycin')){j<-d$drug==dr;n<-sum(ix&j);removed<-if(is.null(prev))NA_integer_ else sum(prev&j)-n;out[[length(out)+1]]<-data.table(stage=k,drug=dr,n=as.character(n),removed=if(!is.na(removed)&&removed>0&&removed<10)'<10' else as.character(removed))}
 prev<-ix
}
fwrite(rbindlist(out),file.path(O,'FLOW_v0.12.csv'))
ms<-fromJSON(file.path(R,'reports/BASELINE_MODEL_READINESS_v0.10.json'))$models
gate<-rbindlist(lapply(ms,function(z)data.table(population=z$population,n_LZD=z$n_LZD,n_VAN=z$n_VAN,max_SMD=z$max_all_level_SMD,pass=z$readiness_gate)))
stopifnot(sum(gate$pass)==3);fwrite(gate,file.path(O,'READINESS_SELECTION_v0.12.csv'))
# Weight dispersion is descriptive; no endpoint recalculation or changes.
m<-readRDS(file.path(R,'cache/phase11_routine_culture_no_recorded_lmwh_shared.rds'));a<-m$raw$drug=='linezolid'
write_json(list(primary_admissions=nrow(m$raw),primary_patients=length(unique(m$raw$subject_id)),three_of_five_readiness_passed=TRUE,weight_summary_LZD=as.list(summary(m$weights[a])),weight_summary_VAN=as.list(summary(m$weights[!a])),no_new_outcome_models=TRUE),file.path(O,'FLOW_CHECKS_v0.12.json'),pretty=TRUE,auto_unbox=TRUE)
cat('Flow and selection audit complete; three of five candidate populations pass.\n')
