# Reweight the already defined v0.6.1/v0.7 interpretation tools; no new outcome definitions.
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
R<-'/root/projects/linezolid_platelet_recovery';pops<-fromJSON(file.path(R,'config/analysis_plan_v0.9.json'))$clinical$populations
old<-as.data.frame(fread(file.path(R,'cache/phase61_endpoints.csv')));ep<-as.data.frame(fread(file.path(R,'cache/phase7_endpoints.csv')))
curves<-list();bins<-list();risks<-list();tabs<-list();identity_error<-0
for(pop in pops){
 m<-readRDS(file.path(R,'cache',paste0('phase9_',pop,'.rds')));d<-ep[match(m$data$hadm_id,ep$hadm_id),];od<-old[match(m$data$hadm_id,old$hadm_id),];a<-m$data$A;w<-m$weights
 stopifnot(identical(as.character(d$hadm_id),as.character(m$data$hadm_id)),identical(as.character(od$hadm_id),as.character(d$hadm_id)),all(d$hospital_status==od$status),all(abs(d$hospital_day-od$event_day)<1e-10))
 for(arm in 0:1){ix<-a==arm;g<-seq(0,14,by=.25);q<-t(vapply(g,function(t)vapply(1:3,function(k)sum(w[ix][d$hospital_status[ix]==k & d$hospital_day[ix]<=t])/sum(w[ix]),0.0),numeric(3)))
  curves[[paste(pop,arm)]]<-data.table(population=pop,drug=if(arm)'linezolid' else 'vancomycin',day=g,recovery=q[,1],death=q[,2],discharge=q[,3],no_first_event=1-rowSums(q))
  rr<-data.table(population=pop,drug=if(arm)'linezolid' else 'vancomycin',day=c(0,3,7,14),n_at_risk=vapply(c(0,3,7,14),function(t)sum(d$hospital_day[ix]>=t),0L));rr[,n_at_risk:=ifelse(n_at_risk>0&n_at_risk<10,'<10',as.character(n_at_risk))];risks[[paste(pop,arm)]]<-rr
 }
 ld<-as.data.table(od[od$live_discharge,]);ld[,weight:=w[od$live_discharge]];ld[,discharge_window:=ifelse(discharge_day<=14,'Within14d','After14d')];ld[,bin:=ifelse(is.na(last_platelets),'Missing / ambiguous',ifelse(last_platelets<50,'Below 50',ifelse(last_platelets<100,'50 to 99','100 or above')))];lb<-ld[,.(n=.N,weight_sum=sum(weight)),by=.(drug,discharge_window,bin)];lb[,fraction:=weight_sum/sum(weight_sum),by=.(drug,discharge_window)];lb[,n:=ifelse(n>0&n<10,'<10',as.character(n))];lb[,population:=pop];bins[[pop]]<-lb
 mat<-cbind(confirmed=as.numeric(d$hospital_status==1),single=as.numeric(d$hospital_single_status==1),CIF_area_days=(14-d$hospital_day)*(d$hospital_status==1),never_threshold=as.numeric(d$reason=='never_threshold'),unconfirmed_observed_low=as.numeric(d$reason=='unconfirmed_observed_low'),unconfirmed_insufficient_clock=as.numeric(d$reason=='unconfirmed_insufficient_clock'),unconfirmed_no_late_measurement=as.numeric(d$reason=='unconfirmed_no_late_measurement'),hospital_tests=d$hospital_tests,hospital_observation_days=d$hospital_days)
 stat<-function(ix,w){z<-mat[ix,,drop=FALSE];aa<-a[ix];v1<-colSums(z[aa==1,,drop=FALSE]*w[aa==1])/sum(w[aa==1]);v0<-colSums(z[aa==0,,drop=FALSE]*w[aa==0])/sum(w[aa==0]);v1<-c(v1,pooled_tests_per_day=unname(v1['hospital_tests']/v1['hospital_observation_days']));v0<-c(v0,pooled_tests_per_day=unname(v0['hospital_tests']/v0['hospital_observation_days']));cbind(LZD=v1,VAN=v0,RD=v1-v0)}
 point<-stat(seq_len(nrow(d)),w);bs<-array(NA_real_,c(length(m$bootstrap),nrow(point),3))
 for(i in seq_along(m$bootstrap)){b<-m$bootstrap[[i]];if(!is.null(b))bs[i,,]<-stat(b$ix,b$weights)}
 # Confirmation-gap decomposition must be an identity at every refitted draw.
 jj<-match(c('unconfirmed_observed_low','unconfirmed_insufficient_clock','unconfirmed_no_late_measurement'),rownames(point));gap<-bs[,match('single',rownames(point)),3]-bs[,match('confirmed',rownames(point)),3];identity_error<-max(identity_error,max(abs(gap-rowSums(bs[,jj,3])),na.rm=TRUE))
 tabs[[pop]]<-rbindlist(lapply(seq_len(nrow(point)),function(j){qq<-quantile(bs[,j,3],c(.025,.975),na.rm=TRUE);data.table(population=pop,metric=rownames(point)[j],LZD=point[j,1],VAN=point[j,2],RD=point[j,3],lower=unname(qq[1]),upper=unname(qq[2]),diagnostic_gate=m$report$diagnostic_gate)}))
}
stopifnot(identity_error<1e-10)
fwrite(rbindlist(curves),file.path(R,'reports/FIRST_EVENT_CURVES_v0.9.csv'));fwrite(rbindlist(risks),file.path(R,'reports/RISK_TABLE_v0.9.csv'));fwrite(rbindlist(bins),file.path(R,'reports/LAST_DISCHARGE_COUNT_BINS_v0.9.csv'));fwrite(rbindlist(tabs),file.path(R,'reports/OBSERVATION_CONTEXT_v0.9.csv'))
write_json(list(utc=format(Sys.time(),tz='UTC',usetz=TRUE),confirmation_gap_identity_max_error=identity_error,grid_step_days=.25,endpoint_definitions_reused=TRUE,conditional_discharge_summaries_descriptive=TRUE,CIF_area_is_not_recovered_alive_time=TRUE),file.path(R,'reports/OBSERVATION_CONTEXT_CHECKS_v0.9.json'),pretty=TRUE,auto_unbox=TRUE,digits=12)
print(rbindlist(tabs))
