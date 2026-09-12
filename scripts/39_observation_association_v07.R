suppressPackageStartupMessages({library(data.table);library(jsonlite);library(survival)})
R<-'/root/projects/linezolid_platelet_recovery'
m<-readRDS(file.path(R,'cache/phase61_overlap.rds'))
draws<-readRDS(file.path(R,'cache/phase7_bootstrap.rds'))$bootstrap
ep<-as.data.frame(fread(file.path(R,'cache/phase7_endpoints.csv')))
d<-ep[match(m$data$hadm_id,ep$hadm_id),];w<-m$weights
stopifnot(nrow(d)==1569,!anyNA(d$hadm_id),identical(as.character(d$hadm_id),as.character(m$data$hadm_id)),identical(d$drug,m$data$drug))
d$A<-m$data$A
prefixes<-c('hospital','hospital_single','icu','icu_single','hospital_timed','hospital_dated','icu_timed','icu_dated')
met<-list();kind<-character()
add<-function(name,value,type){met[[name]]<<-value;kind[name]<<-type}
for(pr in prefixes){
 s<-d[[paste0(pr,'_status')]]
 for(k in 0:3)add(paste0(pr,'__',c('no_event','recovery','death','live_exit')[k+1]),as.numeric(s==k),'probability')
}
add('hospital__CIF_area_days',(14-d$hospital_day)*(d$hospital_status==1),'days')
for(k in c('confirmed','never_threshold','unconfirmed_observed_low','unconfirmed_insufficient_clock','unconfirmed_no_late_measurement'))add(paste0('reason__',k),as.numeric(d$reason==k),'probability')
add('reason__all_unconfirmed',as.numeric(d$hospital_single_status==1 & d$hospital_status!=1),'probability')
for(k in c('hospital_tests','hospital_days','hospital_tests_per_day','specimen_groups','canonical_rows','icu_tests','icu_days','icu_tests_per_day','hospital_icu_days','hospital_outside_icu_days','hospital_icu_tests','hospital_outside_icu_tests'))add(paste0('process__',k),d[[k]],'process')
add('conditional__confirmation_among_threshold',ifelse(d$hospital_single_status==1,as.numeric(d$hospital_status==1),NA_real_),'conditional_probability')
add('conditional__first_high_day',d$first_high_day,'conditional_days')
mat<-do.call(cbind,met)
stat<-function(ix,w){
 z<-mat[ix,,drop=FALSE];a<-d$A[ix]
 arms<-lapply(c(1,0),function(v){ii<-a==v;y<-z[ii,,drop=FALSE];ww<-w[ii];present<-!is.na(y);y[!present]<-0
  colSums(y*ww)/colSums(present*ww)})
 # Distinguish the weighted mean of individual rates from the pooled rate.
 for(pr in c('hospital','icu'))for(j in 1:2){
  ii<-a==c(1,0)[j];v<-sum(w[ii]*d[[paste0(pr,'_tests')]][ix[ii]])/sum(w[ii]*d[[paste0(pr,'_days')]][ix[ii]])
  arms[[j]][paste0('process__',pr,'_pooled_tests_per_day')]<-v
 }
 cbind(LZD=arms[[1]],VAN=arms[[2]],RD=arms[[1]]-arms[[2]])
}
point<-stat(seq_len(nrow(d)),w);plain<-stat(seq_len(nrow(d)),rep(1,nrow(d)))
kind<-c(kind,process__hospital_pooled_tests_per_day='process',process__icu_pooled_tests_per_day='process')
boot<-array(NA_real_,c(length(draws),nrow(point),3),dimnames=list(NULL,rownames(point),colnames(point)))
for(b in seq_along(draws)){
 ob<-draws[[b]];if(is.null(ob))next
 boot[b,,]<-stat(ob$ix,ob$weights)
}
old<-readRDS(file.path(R,'cache/phase61_clinical_bootstrap.rds'))$boot
check<-max(abs(boot[1:500,'hospital__recovery','RD']-old[,'RD']),abs(boot[1:500,'hospital__CIF_area_days','RD']-old[,'area_difference']))
stopifnot(check<1e-12)
# Check nesting and arithmetic at every draw, not merely point estimates.
err<-max(abs(boot[,'hospital_single__recovery','RD']-boot[,'hospital__recovery','RD']-boot[,'reason__all_unconfirmed','RD']),na.rm=TRUE)
stopifnot(err<1e-12)
ajerr<-0
for(pr in prefixes)for(a in 0:1){
 ii<-d$A==a;tt<-d[[paste0(pr,'_day')]][ii];ss<-d[[paste0(pr,'_status')]][ii]
 fit<-survfit(Surv(tt,factor(ss,levels=0:3,labels=c('censor','recovery','death','live_exit')))~1,weights=w[ii],se.fit=FALSE,time0=TRUE,timefix=FALSE)
 for(k in 1:3){j<-match(c('recovery','death','live_exit')[k],fit$states)
  expected<-vapply(fit$time,function(t)sum(w[ii][ss==k & tt<=t])/sum(w[ii]),0.0)
  ajerr<-max(ajerr,max(abs(fit$pstate[,j]-expected)))
 }
}
stopifnot(ajerr<1e-10)
qfun<-function(z)unname(quantile(z,c(.025,.975),na.rm=TRUE,type=7))
tab<-rbindlist(lapply(seq_len(nrow(point)),function(j){
 qq<-apply(boot[,j,],2,qfun);qq500<-qfun(boot[1:500,j,'RD'])
 data.table(metric=rownames(point)[j],kind=unname(kind[rownames(point)[j]]),LZD=point[j,1],VAN=point[j,2],RD=point[j,3],
  LZD_lower=qq[1,1],LZD_upper=qq[2,1],VAN_lower=qq[1,2],VAN_upper=qq[2,2],lower=qq[1,3],upper=qq[2,3],
  lower_500=qq500[1],upper_500=qq500[2],unweighted_LZD=plain[j,1],unweighted_VAN=plain[j,2])
}))
valid<-apply(boot,1,function(x)all(is.finite(x)));stopifnot(mean(valid)>=.98)
# A second Monte Carlo diagnostic resamples bootstrap statistic vectors, never clinical rows.
set.seed(20260912);mc<-array(NA_real_,c(1000,nrow(point),2))
vals<-boot[valid,,'RD'];B<-nrow(vals)
for(i in 1:1000){ix<-sample.int(B,B,replace=TRUE);mc[i,,]<-t(apply(vals[ix,,drop=FALSE],2,qfun))}
tab[,MCSE_lower:=apply(mc[,,1],2,sd)];tab[,MCSE_upper:=apply(mc[,,2],2,sd)]
fwrite(tab,file.path(R,'reports/OBSERVATION_ASSOCIATION_v0.7.csv'))
counts<-rbindlist(lapply(prefixes,function(pr)data.table(drug=d$drug,status=d[[paste0(pr,'_status')]])[,.(n=.N),by=.(drug,status)][,endpoint:=pr]))
counts[,state:=c('no_event','recovery','death','live_exit')[status+1]];counts[,status:=NULL]
counts[,n:=ifelse(n>0&n<10,'<10',as.character(n))]
fwrite(counts,file.path(R,'reports/OBSERVATION_EVENT_COUNTS_v0.7.csv'))
saveRDS(list(point=point,boot=boot,valid=valid),file.path(R,'cache/phase7_association.rds'))
report<-list(completed_utc=format(Sys.time(),tz='UTC',usetz=TRUE),n=1569,B_requested=length(draws),B_valid=sum(valid),
 primary_first500_bootstrap_max_error=check,algebraic_identity_all_draws_max_error=err,AJ_empirical_max_error=ajerr,
 MCSE_method='SD of interval endpoint from1000 resamplings of2000 bootstrap-statistic vectors; seed20260912; not a sampling standard error.',
 interpretation='All added endpoints exploratory. ICU exits change estimand. No measurement-bias percentage or true-effect bounds identified.',estimates=tab)
write_json(report,file.path(R,'reports/OBSERVATION_ASSOCIATION_v0.7.json'),pretty=TRUE,auto_unbox=TRUE,digits=10)
print(tab[metric%in%c('hospital__recovery','hospital_single__recovery','icu__recovery','icu_single__recovery','hospital_timed__recovery','hospital_dated__recovery','icu_timed__recovery','icu_dated__recovery','hospital__CIF_area_days','reason__unconfirmed_observed_low','reason__unconfirmed_insufficient_clock','reason__unconfirmed_no_late_measurement'),.(metric,LZD,VAN,RD,lower,upper,MCSE_lower,MCSE_upper)])
