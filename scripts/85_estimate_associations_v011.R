# Joint inference from fixed memberships, point weights and patient-cluster refits.
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(survival)})
R<-'/root/projects/linezolid_platelet_recovery';plan<-fromJSON(file.path(R,'config/analysis_plan_v0.11.json'))
ep<-as.data.frame(fread(file.path(R,'cache/phase11_endpoints.csv')))
prefixes<-c('hospital','hospital_single','hospital72','hospital_timed','hospital_dated','physical','physical_single','physical_timed','physical_dated','icu','icu_single','icu_timed','icu_dated')
res<-list();curveout<-list();risks<-list();binout<-list();lagout<-list();rawcounts<-list();ajerr<-0;gaperr<-0
small<-function(x)ifelse(x>0&x<10,'<10',as.character(x))
for(pop in c(plan$populations$primary,plan$populations$sensitivity)){
 m<-readRDS(file.path(R,'cache',paste0('phase11_',pop,'.rds')));stopifnot(length(m$bootstrap)==2000)
 d<-ep[match(m$raw$hadm_id,ep$hadm_id),];stopifnot(!anyNA(d$hadm_id),identical(as.character(d$hadm_id),as.character(m$raw$hadm_id)),identical(d$drug,m$raw$drug))
 a<-m$data$A;w<-m$weights;n<-nrow(d)
 mat<-do.call(cbind,lapply(prefixes,function(pr){s<-d[[paste0(pr,'_status')]];stopifnot(!anyNA(s));cbind(recovery=as.numeric(s==1),death=as.numeric(s==2),live_exit=as.numeric(s==3))}))
 colnames(mat)<-unlist(lapply(prefixes,function(pr)paste(pr,c('recovery','death','live_exit'),sep='__')))
 mat<-cbind(mat,CIF_area_days=(14-d$hospital_day)*(d$hospital_status==1),never_threshold=as.numeric(d$reason=='never_threshold'),unconfirmed_observed_low=as.numeric(d$reason=='unconfirmed_observed_low'),unconfirmed_insufficient_clock=as.numeric(d$reason=='unconfirmed_insufficient_clock'),unconfirmed_no_late_measurement=as.numeric(d$reason=='unconfirmed_no_late_measurement'),hospital_tests=d$hospital_tests,hospital_days=d$hospital_days,hospital_mean_person_tests_per_day=d$hospital_tests_per_day,physical_tests=d$physical_tests,physical_days=d$physical_days,physical_mean_person_tests_per_day=d$physical_tests_per_day)
 basecols<-colnames(mat);grid<-seq(0,14,by=.25)
 cm<-do.call(cbind,lapply(1:3,function(k)sapply(grid,function(t)as.numeric(d$hospital_status==k&d$hospital_day<=t))))
 cn<-unlist(lapply(c('recovery','death','discharge'),function(k)paste0('curve_',k,'_',seq_along(grid))))
 colnames(cm)<-cn;mat<-cbind(mat,cm)
 # Conditional discharge descriptors: ratios, not causal contrasts or latent recovery.
 strata<-list(All=rep(TRUE,n),ConfirmedBeforeDischarge=d$confirmed_before_discharge,NotConfirmedBeforeDischarge=!d$confirmed_before_discharge)
 bins<-list(Below50=!is.na(d$last_platelets)&d$last_platelets<50,From50to99=!is.na(d$last_platelets)&d$last_platelets>=50&d$last_platelets<100,AtLeast100=!is.na(d$last_platelets)&d$last_platelets>=100,MissingOrAmbiguous=is.na(d$last_platelets))
 windows<-list(Within14d=d$discharge_day<=14,After14d=d$discharge_day>14);cond<-list();binmeta<-list();lagmeta<-list()
 for(ww in names(windows))for(st in names(strata)){
  mask<-d$live_discharge&windows[[ww]]&strata[[st]];key<-paste(ww,st,sep='__');den<-paste0('den_',key);cond[[den]]<-as.numeric(mask)
  for(bn in names(bins)){nm<-paste0('bin_',key,'__',bn);cond[[nm]]<-as.numeric(mask&bins[[bn]]);binmeta[[nm]]<-list(window=ww,stratum=st,bin=bn,den=den,mask=mask&bins[[bn]],denmask=mask)}
  lm<-mask&is.finite(d$last_draw_to_discharge_hours);ld<-paste0('lagden_',key);ln<-paste0('lag_',key);cond[[ld]]<-as.numeric(lm);cond[[ln]]<-ifelse(lm,d$last_draw_to_discharge_hours,0);lagmeta[[ln]]<-list(window=ww,stratum=st,den=ld)
 }
 cmat<-do.call(cbind,cond);colnames(cmat)<-names(cond);mat<-cbind(mat,cmat);storage.mode(mat)<-'double'
 stat<-function(ix,ww){
  sw<-numeric(n);rs<-rowsum(ww,ix,reorder=FALSE);sw[as.integer(rownames(rs))]<-rs[,1]
  v1<-drop(crossprod(sw*(a==1),mat))/sum(sw[a==1]);v0<-drop(crossprod(sw*(a==0),mat))/sum(sw[a==0])
  cbind(LZD=v1,VAN=v0,RD=v1-v0)
 }
 pointall<-stat(seq_len(n),w);bs<-array(NA_real_,c(2000,ncol(mat),3),dimnames=list(NULL,colnames(mat),c('LZD','VAN','RD')))
 for(i in seq_len(2000)){bb<-m$bootstrap[[i]];if(!is.null(bb))bs[i,,]<-stat(bb$ix,bb$weights)}
 # Append arm-specific nonlinear summaries using the same draws.
 basepoint<-pointall[basecols,,drop=FALSE];basebs<-bs[,basecols,,drop=FALSE]
 for(pre in c('hospital','physical')){
  nm<-paste0(pre,'_pooled_tests_per_day');nu<-paste0(pre,'_tests');de<-paste0(pre,'_days')
  pp<-pointall[nu,1:2]/pointall[de,1:2];xx<-bs[,nu,1:2]/bs[,de,1:2]
  basepoint<-rbind(basepoint,setNames(c(pp,pp[1]-pp[2]),c('LZD','VAN','RD')));rownames(basepoint)[nrow(basepoint)]<-nm
  arr<-array(c(xx[,1],xx[,2],xx[,1]-xx[,2]),c(2000,1,3),dimnames=list(NULL,nm,c('LZD','VAN','RD')))
  # Append within each arm block, preserving metric and bootstrap alignment.
  oldcols<-rownames(basepoint)[-nrow(basepoint)];tmp<-array(NA_real_,c(2000,nrow(basepoint),3),dimnames=list(NULL,rownames(basepoint),c('LZD','VAN','RD')))
  for(k in 1:3){tmp[,oldcols,k]<-basebs[,oldcols,k];tmp[,nm,k]<-arr[,1,k]}
  basebs<-tmp
 }
 qs<-function(v)if(all(!is.finite(v)))c(NA_real_,NA_real_)else unname(quantile(v,c(.025,.975),na.rm=TRUE,type=7))
 res[[pop]]<-rbindlist(lapply(seq_len(nrow(basepoint)),function(j){q<-qs(basebs[,j,3]);q1<-qs(basebs[,j,1]);q0<-qs(basebs[,j,2]);data.table(population=pop,metric=rownames(basepoint)[j],LZD=basepoint[j,1],VAN=basepoint[j,2],RD=basepoint[j,3],lower=q[1],upper=q[2],LZD_lower=q1[1],LZD_upper=q1[2],VAN_lower=q0[1],VAN_upper=q0[2],valid_draws=sum(is.finite(basebs[,j,3])),diagnostic_gate=m$report$readiness_gate)}))
 jj<-c('unconfirmed_observed_low','unconfirmed_insufficient_clock','unconfirmed_no_late_measurement');gaperr<-max(gaperr,max(abs(basebs[,'hospital_single__recovery',3]-basebs[,'hospital__recovery',3]-rowSums(basebs[,jj,3])),na.rm=TRUE))
 for(arm in 0:1){k<-if(arm)1 else 2;drug<-if(arm)'linezolid' else 'vancomycin';ix<-a==arm
  cc<-sapply(c('recovery','death','discharge'),function(state)pointall[paste0('curve_',state,'_',seq_along(grid)),k]);cc<-cbind(cc,no_first_event=1-rowSums(cc))
  for(j in seq_along(grid))for(st in colnames(cc)){
   v<-if(st=='no_first_event')1-rowSums(bs[,paste0('curve_',c('recovery','death','discharge'),'_',j),k])else bs[,paste0('curve_',st,'_',j),k]
   q<-qs(v);curveout[[length(curveout)+1]]<-data.table(population=pop,drug=drug,day=grid[j],state=st,estimate=cc[j,st],lower=q[1],upper=q[2])
  }
  rr<-data.table(population=pop,drug=drug,day=c(0,3,7,14));rr[,n_at_risk:=small(vapply(day,function(t)sum(d$hospital_day[ix]>=t),0L))];rr[,weighted_at_risk_fraction:=vapply(day,function(t)sum(w[ix][d$hospital_day[ix]>=t])/sum(w[ix]),0.0)];risks[[length(risks)+1]]<-rr
  for(nm in names(binmeta)){
   z<-binmeta[[nm]];v<-bs[,nm,k]/bs[,z$den,k];q<-qs(v)
   binout[[length(binout)+1]]<-data.table(population=pop,drug=drug,discharge_window=z$window,prior_confirmation=z$stratum,bin=z$bin,n=small(sum(z$mask&ix)),denominator=small(sum(z$denmask&ix)),fraction=pointall[nm,k]/pointall[z$den,k],lower=q[1],upper=q[2],valid_draws=sum(is.finite(v)))
  }
  for(nm in names(lagmeta)){z<-lagmeta[[nm]];v<-bs[,nm,k]/bs[,z$den,k];q<-qs(v);lagout[[length(lagout)+1]]<-data.table(population=pop,drug=drug,discharge_window=z$window,prior_confirmation=z$stratum,mean_last_draw_to_discharge_hours=pointall[nm,k]/pointall[z$den,k],lower=q[1],upper=q[2])}
 }
 for(pr in prefixes){
  for(arm in 0:1){ix<-a==arm;tt<-d[[paste0(pr,'_day')]][ix];ss<-d[[paste0(pr,'_status')]][ix];ww<-w[ix]
   sf<-survfit(Surv(tt,factor(ss,levels=0:3,labels=c('censor','recovery','death','live_exit')))~1,weights=ww,se.fit=FALSE,time0=TRUE,timefix=FALSE)
   for(k in 1:3){jj<-match(c('recovery','death','live_exit')[k],sf$states);ex<-vapply(sf$time,function(t)sum(ww[ss==k&tt<=t])/sum(ww),0.0);ajerr<-max(ajerr,max(abs(sf$pstate[,jj]-ex)))}
   for(k in 0:3)rawcounts[[length(rawcounts)+1]]<-data.table(population=pop,drug=if(arm)'linezolid' else 'vancomycin',endpoint=pr,state=c('no_first_event','recovery','death','live_exit')[k+1],n=small(sum(ss==k)))
  }
 }
 saveRDS(list(point=basepoint,bootstrap=basebs),file.path(R,'cache',paste0('phase11_associations_',pop,'.rds')))
 cat('Finished',pop,'\n');flush.console()
}
stopifnot(ajerr<1e-10,gaperr<1e-10)
fwrite(rbindlist(res),file.path(R,'reports/CLINICAL_ASSOCIATIONS_v0.11.csv'))
fwrite(rbindlist(curveout),file.path(R,'reports/FIRST_EVENT_CURVES_v0.11.csv'))
fwrite(rbindlist(risks),file.path(R,'reports/RISK_TABLE_v0.11.csv'))
fwrite(rbindlist(binout),file.path(R,'reports/LAST_DISCHARGE_COUNT_BINS_v0.11.csv'))
fwrite(rbindlist(lagout),file.path(R,'reports/LAST_DISCHARGE_LAGS_v0.11.csv'))
fwrite(rbindlist(rawcounts),file.path(R,'reports/FIRST_EVENT_COUNTS_v0.11.csv'))
write_json(list(utc=format(Sys.time(),tz='UTC',usetz=TRUE),AJ_max_error=ajerr,confirmation_gap_identity_max_error=gaperr,bootstrap_B=2000,grid_step_days=.25,pointwise_intervals=TRUE,conditional_discharge_descriptive=TRUE,last_confirmation_stratum='Any confirmed recovery before live discharge, including after14d for late discharges',CIF_area_is_not_alive_recovered_time=TRUE),file.path(R,'reports/ASSOCIATION_CHECKS_v0.11.json'),pretty=TRUE,auto_unbox=TRUE,digits=12)
print(rbindlist(res)[metric%in%c('hospital__recovery','hospital_single__recovery','hospital__death','hospital__live_exit','physical__recovery','hospital_pooled_tests_per_day')])
