suppressPackageStartupMessages({library(data.table);library(jsonlite);library(survival)})
R<-'/root/projects/linezolid_platelet_recovery';cfg<-fromJSON(file.path(R,'config/analysis_plan_v0.8.json'))$clinical
ep<-as.data.frame(fread(file.path(R,'cache/phase7_endpoints.csv')));res<-list();counts<-list();ajerr<-0
for(pop in cfg$populations){
 m<-readRDS(file.path(R,'cache',paste0('phase8_',pop,'.rds')));d<-ep[match(m$data$hadm_id,ep$hadm_id),]
 stopifnot(!anyNA(d$hadm_id),identical(d$drug,m$data$drug));a<-m$data$A
 prefixes<-c('hospital','hospital_single','icu','icu_single');if(pop!='availability_corrected_all')prefixes<-c(prefixes,'physical','physical_single')
 mat<-do.call(cbind,lapply(prefixes,function(pr){s<-d[[paste0(pr,'_status')]];stopifnot(!anyNA(s));cbind(recovery=as.numeric(s==1),death=as.numeric(s==2),live_exit=as.numeric(s==3))}))
 colnames(mat)<-unlist(lapply(prefixes,function(pr)paste(pr,c('recovery','death','live_exit'),sep='__')))
 stat<-function(ix,w){z<-mat[ix,,drop=FALSE];aa<-a[ix];v1<-colSums(z[aa==1,,drop=FALSE]*w[aa==1])/sum(w[aa==1]);v0<-colSums(z[aa==0,,drop=FALSE]*w[aa==0])/sum(w[aa==0]);cbind(LZD=v1,VAN=v0,RD=v1-v0)}
 point<-stat(seq_len(nrow(d)),m$weights);boot<-array(NA_real_,c(length(m$bootstrap),nrow(point),3))
 for(i in seq_along(m$bootstrap)){b<-m$bootstrap[[i]];if(!is.null(b))boot[i,,]<-stat(b$ix,b$weights)}
 z<-rbindlist(lapply(seq_len(nrow(point)),function(j){q<-quantile(boot[,j,3],c(.025,.975),na.rm=TRUE);data.table(population=pop,metric=rownames(point)[j],LZD=point[j,1],VAN=point[j,2],RD=point[j,3],lower=unname(q[1]),upper=unname(q[2]),diagnostic_gate=m$report$diagnostic_gate,n_LZD=sum(a==1),n_VAN=sum(a==0),ESS_LZD=m$report$ESS_LZD,ESS_VAN=m$report$ESS_VAN)}))
 res[[pop]]<-z
 for(pr in prefixes){
  cnt<-data.table(drug=d$drug,status=as.character(d[[paste0(pr,'_status')]]))[,.(n=.N),by=.(drug,status)];cnt[,`:=`(population=pop,endpoint=pr,n=ifelse(n>0&n<10,'<10',as.character(n)))];counts[[paste(pop,pr)]]<-cnt
  for(arm in 0:1){ii<-a==arm;tt<-d[[paste0(pr,'_day')]][ii];ss<-d[[paste0(pr,'_status')]][ii];w<-m$weights[ii]
   sf<-survfit(Surv(tt,factor(ss,levels=0:3,labels=c('censor','recovery','death','live_exit')))~1,weights=w,se.fit=FALSE,time0=TRUE,timefix=FALSE)
   for(k in 1:3){jj<-match(c('recovery','death','live_exit')[k],sf$states);ex<-vapply(sf$time,function(t)sum(w[ss==k&tt<=t])/sum(w),0.0);ajerr<-max(ajerr,max(abs(sf$pstate[,jj]-ex)))}
  }
 }
 saveRDS(list(point=point,boot=boot),file.path(R,'cache',paste0('phase8_association_',pop,'.rds')))
}
stopifnot(ajerr<1e-10);z<-rbindlist(res);fwrite(z,file.path(R,'reports/CLINICAL_ASSOCIATIONS_v0.8.csv'));fwrite(rbindlist(counts),file.path(R,'reports/CLINICAL_EVENT_COUNTS_v0.8.csv'))
write_json(list(completed_utc=format(Sys.time(),tz='UTC',usetz=TRUE),AJ_empirical_max_error=ajerr,bootstrap_B=cfg$bootstrap$B,estimates=z,interpretation='All post-primary exploratory analyses; no multiplicity-adjusted confirmatory inference. Different population definitions change overlap estimands. Failed diagnostic gates remain diagnostic only.'),file.path(R,'reports/CLINICAL_ASSOCIATIONS_v0.8.json'),pretty=TRUE,auto_unbox=TRUE,digits=10)
print(z[metric%in%c('hospital__recovery','hospital_single__recovery','physical__recovery'),.(population,metric,LZD,VAN,RD,lower,upper,diagnostic_gate)])
