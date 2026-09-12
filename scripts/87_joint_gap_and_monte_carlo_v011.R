# Joint contrast of previously specified endpoints; interval simulation-error diagnostics only.
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
R<-'/root/projects/linezolid_platelet_recovery';p<-fromJSON(file.path(R,'config/analysis_plan_v0.11.json'));res<-list();mc<-list()
set.seed(20260916)
for(pop in c(p$populations$primary,p$populations$sensitivity)){
 a<-readRDS(file.path(R,'cache',paste0('phase11_associations_',pop,'.rds')))
 pp<-a$point['hospital_single__recovery',]-a$point['hospital__recovery',]
 bb<-a$bootstrap[,'hospital_single__recovery',]-a$bootstrap[,'hospital__recovery',]
 q<-quantile(bb[,3],c(.025,.975),na.rm=TRUE,type=7)
 res[[pop]]<-data.table(population=pop,metric='single_minus_confirmed_gap',LZD=pp[1],VAN=pp[2],RD=pp[3],lower=unname(q[1]),upper=unname(q[2]))
 for(metric in c('hospital__recovery','hospital_single__recovery','physical__recovery','physical_timed__recovery','hospital_timed__recovery')){
  v<-a$bootstrap[,metric,3];v<-v[is.finite(v)];qq<-replicate(1000,quantile(sample(v,length(v),replace=TRUE),c(.025,.975),type=7));se<-apply(qq,1,sd)
  mc[[length(mc)+1]]<-data.table(population=pop,metric=metric,lower_endpoint_MCSE=se[1],upper_endpoint_MCSE=se[2],simulations=1000,seed=20260916)
 }
}
fwrite(rbindlist(res),file.path(R,'reports/CONFIRMATION_GAP_v0.11.csv'))
fwrite(rbindlist(mc),file.path(R,'reports/INTERVAL_MONTE_CARLO_v0.11.csv'))
write_json(list(utc=format(Sys.time(),tz='UTC',usetz=TRUE),purpose='Joint algebraic contrast and numerical error of already fixed estimates; not a causal decomposition or additional independent outcome.',simulation_error='Resampling the retained bootstrap statistic vector estimates endpoint Monte Carlo variability, not clinical sampling uncertainty. No additional patient bootstrap draws or model changes.',disclosure='Numerical audit added after first v0.11 estimates; does not change population, outcomes, weights or confidence-interval rule.'),file.path(R,'reports/INTERVAL_DIAGNOSTICS_v0.11.json'),pretty=TRUE,auto_unbox=TRUE)
print(rbindlist(res));print(rbindlist(mc))
