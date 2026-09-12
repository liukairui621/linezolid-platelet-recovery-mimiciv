# Independent direct weighted summaries and deterministic draw checks, server only.
suppressPackageStartupMessages({library(brglm2);library(data.table);library(jsonlite)})
R<-'/root/projects/linezolid_platelet_recovery'
for(e in parse(file.path(R,'scripts/76_baseline_model_readiness_v010.R'))){if(is.call(e)&&identical(e[[1]],as.name('for'))&&identical(e[[2]],as.name('pop')))break;eval(e)}
plan<-fromJSON(file.path(R,'config/analysis_plan_v0.11.json'));ep<-as.data.frame(fread(file.path(R,'cache/phase11_endpoints.csv')));results<-list()
for(pop in c(plan$populations$primary,plan$populations$sensitivity)){
 m<-readRDS(file.path(R,'cache',paste0('phase11_',pop,'.rds')));old<-readRDS(file.path(R,'cache',paste0('phase10_readiness_',pop,'.rds')));aout<-readRDS(file.path(R,'cache',paste0('phase11_associations_',pop,'.rds')))
 raw<-m$raw;x<-prep(raw);point<-fitone(x);stopifnot(identical(raw,old$raw),max(abs(point$w-m$weights))<1e-10,max(abs(m$weights-old$weights))==0)
 d<-ep[match(raw$hadm_id,ep$hadm_id),];aa<-x$A
 vv<-cbind(hospital__recovery=as.numeric(d$hospital_status==1),hospital_single__recovery=as.numeric(d$hospital_single_status==1),hospital__death=as.numeric(d$hospital_status==2),hospital__live_exit=as.numeric(d$hospital_status==3),physical__recovery=as.numeric(d$physical_status==1),CIF_area_days=(14-d$hospital_day)*(d$hospital_status==1))
 direct<-function(ix,w){
  z<-vv[ix,,drop=FALSE];a<-aa[ix];v1<-vapply(seq_len(ncol(z)),function(k)weighted.mean(z[a==1,k],w[a==1]),0.0);v0<-vapply(seq_len(ncol(z)),function(k)weighted.mean(z[a==0,k],w[a==0]),0.0)
  p1<-sum(w[a==1]*d$hospital_tests[ix][a==1])/sum(w[a==1]*d$hospital_days[ix][a==1]);p0<-sum(w[a==0]*d$hospital_tests[ix][a==0])/sum(w[a==0]*d$hospital_days[ix][a==0])
  ans<-cbind(LZD=c(v1,p1),VAN=c(v0,p0),RD=c(v1-v0,p1-p0));rownames(ans)<-c(colnames(vv),'hospital_pooled_tests_per_day');ans
 }
 check<-direct(seq_len(nrow(raw)),m$weights);point_err<-max(abs(check-aout$point[rownames(check),]));stopifnot(point_err<1e-10)
 cl<-split(seq_len(nrow(raw)),raw$subject_id);pat<-vapply(cl,function(ii)paste(sort(unique(aa[ii])),collapse=','),'');strata<-split(names(cl),pat);set.seed(plan$bootstrap$seed)
 checked<-0;err<-0;refiterr<-0
 for(i in 1:2000){
  nd<-unlist(lapply(strata,function(z)sample(z,length(z),replace=TRUE)),use.names=FALSE);ix<-unlist(cl[nd],use.names=FALSE);bb<-m$bootstrap[[i]]
  if(is.null(bb)){stopifnot(!m$diagnostics$valid[i]);next}
  stopifnot(identical(ix,bb$ix),all(is.finite(bb$weights)),all(bb$weights>0&bb$weights<1));checked<-checked+1
  if(i<=200)stopifnot(identical(bb,old$bootstrap[[i]]))
  vals<-direct(ix,bb$weights);err<-max(err,max(abs(vals-aout$bootstrap[i,rownames(vals),])))
  if(i%in%c(1,200,500,1000,2000)){z<-fitone(prep(raw[ix,,drop=FALSE]));stopifnot(z$valid);refiterr<-max(refiterr,max(abs(z$w-bb$weights)))}
 }
 stopifnot(err<1e-10,refiterr<1e-10,checked==m$report$valid_draws)
 results[[pop]]<-list(membership_unchanged=TRUE,first200_full_draw_objects_identical=TRUE,valid_draw_indices_verified=checked,point_weight_max_error=max(abs(point$w-m$weights)),selected_refits=c(1,200,500,1000,2000),selected_refit_weight_max_error=refiterr,direct_point_statistic_max_error=point_err,direct_all_bootstrap_statistic_max_error=err,independent_statistics_verified=rownames(check))
}
write_json(list(utc=format(Sys.time(),tz='UTC',usetz=TRUE),passed=TRUE,patient_rows_exported=FALSE,populations=results),file.path(R,'reports/INFERENCE_REPRODUCTION_v0.11.json'),pretty=TRUE,auto_unbox=TRUE,digits=12)
print(results)
