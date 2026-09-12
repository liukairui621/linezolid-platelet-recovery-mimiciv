# First exploratory clinical association. No outcome-based model selection.
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(survival);library(ggplot2)})
R<-'/root/projects/linezolid_platelet_recovery';V<-'v0.6.1';start<-format(Sys.time(),tz='UTC',usetz=TRUE)
m<-readRDS(file.path(R,'cache/phase61_overlap.rds'));ep<-as.data.frame(fread(file.path(R,'cache/phase61_endpoints.csv')))
d<-ep[match(m$data$hadm_id,ep$hadm_id),];stopifnot(!anyNA(d$hadm_id),identical(as.character(d$hadm_id),as.character(m$data$hadm_id)))
stopifnot(identical(d$drug,m$data$drug),nrow(d)==1569,all(d$event_day[d$status==0]==14),m$report$outcome_free_reporting_gate)
d$A<-as.integer(d$drug=='linezolid');w<-m$weights;grid<-seq(0,14,by=.25)
emp<-function(t,s,w,g){vapply(g,function(tt)vapply(1:3,function(k)sum(w[s==k & t<=tt])/sum(w),0.0),numeric(3))}
stats<-function(z,w,suffix=''){
 t<-z[[paste0('event_day',suffix)]];s<-z[[paste0('status',suffix)]]
 v<-lapply(c(1,0),function(a){ix<-z$A==a;c(recovery=weighted.mean(s[ix]==1,w[ix]),
  area=weighted.mean((14-t[ix])*(s[ix]==1),w[ix]),death=weighted.mean(s[ix]==2,w[ix]),discharge=weighted.mean(s[ix]==3,w[ix]))})
 c(recovery_LZD=v[[1]][1],recovery_VAN=v[[2]][1],RD=v[[1]][1]-v[[2]][1],
   area_LZD=v[[1]][2],area_VAN=v[[2]][2],area_difference=v[[1]][2]-v[[2]][2],
   death_LZD=v[[1]][3],death_VAN=v[[2]][3],discharge_LZD=v[[1]][4],discharge_VAN=v[[2]][4])
}
# Independent survival package AJ verification, including administrative14d ties.
aj_checks<-list()
for(a in 0:1)for(suf in c('','_72h')){
 ix<-d$A==a;z<-d[ix,];tt<-z[[paste0('event_day',suf)]];ss<-factor(z[[paste0('status',suf)]],levels=0:3,labels=c('censor','recovery','death','discharge'))
 fit<-survfit(Surv(tt,ss)~1,weights=w[ix],se.fit=FALSE,time0=TRUE,timefix=FALSE)
 p<-fit$pstate[,match(c('recovery','death','discharge'),fit$states),drop=FALSE]
 q<-t(emp(tt,as.integer(ss)-1,w[ix],fit$time));err<-max(abs(p-q));stopifnot(err<1e-10)
 aj_checks[[paste(a,suf)]]<-err
}
point<-stats(d,w);plain<-stats(d,rep(1,nrow(d)));gap<-stats(d,w,'_72h')
names(point)<-names(plain)<-names(gap)<-sub('\\..*','',names(point))
boot<-matrix(NA_real_,length(m$bootstrap),length(point),dimnames=list(NULL,names(point)));boot72<-boot
curveboot<-array(NA_real_,c(length(m$bootstrap),2,3,length(grid)))
for(b in seq_along(m$bootstrap)){
 ob<-m$bootstrap[[b]];if(is.null(ob))next
 z<-d[ob$ix,,drop=FALSE];boot[b,]<-unname(stats(z,ob$weights));boot72[b,]<-unname(stats(z,ob$weights,'_72h'))
 for(a in 0:1){ii<-z$A==a;curveboot[b,a+1,,]<-emp(z$event_day[ii],z$status[ii],ob$weights[ii],grid)}
}
ci<-function(mat) t(apply(mat,2,quantile,probs=c(.025,.975),na.rm=TRUE))
intervals<-ci(boot);ints72<-ci(boot72)
tab<-data.table(estimand=names(point),estimate=unname(point),lower=intervals[,1],upper=intervals[,2],unweighted=unname(plain),
 estimate_gap72=unname(gap),lower_gap72=ints72[,1],upper_gap72=ints72[,2])
# Same corrected cohort, same formula, MLE only as the prespecified point sensitivity.
ml<-glm(formula(m$fit),data=m$data,family=binomial(),control=glm.control(maxit=100,epsilon=1e-10))
wml<-ifelse(d$A==1,1-fitted(ml),fitted(ml));mlpoint<-stats(d,wml);tab[,unpenalized_point:=unname(mlpoint)]
fwrite(tab,file.path(R,paste0('reports/CLINICAL_ASSOCIATION_',V,'.csv')))
curves<-list();risk<-list();bounds<-list()
for(a in 0:1){
 ix<-d$A==a;g<-sort(unique(c(0,d$event_day[ix],14)));p<-t(emp(d$event_day[ix],d$status[ix],w[ix],g))
 curves[[a+1]]<-data.table(drug=if(a)'linezolid' else 'vancomycin',day=g,recovery=p[,1],death=p[,2],discharge=p[,3],no_first_event=1-rowSums(p))
 risk[[a+1]]<-data.table(drug=if(a)'linezolid' else 'vancomycin',day=c(0,3,7,14),
  n_at_risk=vapply(c(0,3,7,14),function(t)sum(d$event_day[ix]>=t),0L),
  weight_at_risk=vapply(c(0,3,7,14),function(t)sum(w[ix][d$event_day[ix]>=t]),0.0))
 bounds[[a+1]]<-rbindlist(lapply(1:3,function(k){qq<-apply(curveboot[,a+1,k,],2,quantile,probs=c(.025,.975),na.rm=TRUE)
  data.table(drug=if(a)'linezolid' else 'vancomycin',state=c('recovery','death','discharge')[k],day=grid,lower=qq[1,],upper=qq[2,])}))
}
cc<-rbindlist(curves);bb<-rbindlist(bounds);rr<-rbindlist(risk)
fwrite(cc,file.path(R,paste0('reports/FIRST_EVENT_CURVES_',V,'.csv')));fwrite(bb,file.path(R,paste0('reports/CIF_POINTWISE_INTERVALS_',V,'.csv')))
rr[,n_at_risk:=ifelse(n_at_risk>0&n_at_risk<10,'<10',as.character(n_at_risk))];fwrite(rr,file.path(R,paste0('reports/RISK_TABLE_',V,'.csv')))
ld<-as.data.table(d[d$live_discharge,]);ld[,weight:=w[d$live_discharge]]
ld[,discharge_window:=ifelse(discharge_day<=14,'Discharge within 14 days','Discharge after 14 days')]
ld[,bin:=factor(ifelse(is.na(last_platelets),'Missing / ambiguous',ifelse(last_platelets<50,'Below 50',ifelse(last_platelets<100,'50 to 99','100 or above'))),levels=c('Below 50','50 to 99','100 or above','Missing / ambiguous'))]
lastbins<-ld[,.(n=.N,weight_sum=sum(weight)),by=.(drug,discharge_window,bin)]
lastbins[,fraction:=weight_sum/sum(weight_sum),by=.(drug,discharge_window)]
lastsummary<-ld[,.(n=.N,known=sum(!is.na(last_platelets)),median=median(last_platelets,na.rm=TRUE),q25=quantile(last_platelets,.25,na.rm=TRUE),q75=quantile(last_platelets,.75,na.rm=TRUE),
 median_draw_to_discharge_hours=median(last_draw_to_discharge_hours,na.rm=TRUE),ambiguous=sum(last_ambiguous),no_post_t0_sample=sum(!last_post_t0)),by=.(drug,discharge_window,confirmed_before_discharge)]
lastsummary[n<10,c('median','q25','q75','median_draw_to_discharge_hours'):=NA_real_]
for(k in c('n','known','ambiguous','no_post_t0_sample'))lastsummary[,(k):=ifelse(get(k)>0&get(k)<10,'<10',as.character(get(k)))]
lastbins[,n:=ifelse(n>0&n<10,'<10',as.character(n))]
fwrite(lastbins,file.path(R,paste0('reports/LAST_DISCHARGE_COUNT_BINS_',V,'.csv')));fwrite(lastsummary,file.path(R,paste0('reports/LAST_DISCHARGE_COUNT_SUMMARY_',V,'.csv')))
fdir<-file.path(R,'figures',V);dir.create(fdir,recursive=TRUE,showWarnings=FALSE)
theme_set(theme_classic(base_size=12,base_family='sans')+theme(legend.position='bottom',plot.title.position='plot'))
colors<-c('Confirmed recovery'='#3B8D83','Death before recovery'='#B85A5A','Live discharge before recovery'='#D4B86B','No first event yet'='#DEE4E8')
long<-melt(cc,id.vars=c('drug','day'),variable.name='state',value.name='probability')
long[,state:=factor(state,levels=c('recovery','death','discharge','no_first_event'),labels=names(colors))]
# Repeat each x-boundary with its left and right value to preserve exact step geometry.
step<-long[,{
 z<- .SD[order(day)];n<-nrow(z);if(n==1)z else data.table(day=c(z$day[1],rep(z$day[-1],each=2)),probability=c(z$probability[1],as.vector(rbind(z$probability[-n],z$probability[-1]))))
},by=.(drug,state)]
p1<-ggplot(step,aes(day,probability,fill=state))+geom_area(position='stack',stat='identity')+facet_wrap(~drug)+
 scale_fill_manual(values=colors,drop=FALSE)+scale_y_continuous(labels=function(z)paste0(round(z*100),'%'),limits=c(0,1))+
 labs(x='Days from first recorded dose',y='Weighted first-event probability',fill=NULL,title='Recorded recovery, competing exits and no first event',subtitle='First-event categories; recovery remains absorbing after subsequent death or discharge.')
p2<-ggplot(cc,aes(day,recovery,color=drug))+geom_ribbon(data=bb[state=='recovery'],aes(x=day,ymin=lower,ymax=upper,fill=drug),inherit.aes=FALSE,alpha=.16,color=NA)+
 geom_step(linewidth=.8)+scale_color_manual(values=c(linezolid='#B85A5A',vancomycin='#287A9F'))+scale_fill_manual(values=c(linezolid='#B85A5A',vancomycin='#287A9F'))+
 scale_y_continuous(labels=function(z)paste0(round(z*100),'%'),limits=c(0,max(bb[state=='recovery']$upper)*1.05))+
 labs(x='Days from first recorded dose',y='Cumulative confirmed recovery',color=NULL,fill=NULL,title='Exploratory association after overlap weighting',subtitle='Shading: pointwise 95% patient-cluster bootstrap intervals; 500 refits.')
p3<-ggplot(lastbins,aes(drug,fraction,fill=bin))+geom_col(width=.65)+facet_wrap(~discharge_window)+
 scale_fill_manual(values=c('#B85A5A','#D4B86B','#3B8D83','#DEE4E8'),drop=FALSE)+scale_y_continuous(labels=function(z)paste0(round(z*100),'%'),limits=c(0,1))+
 labs(x=NULL,y='Weighted fraction among live discharges',fill='Last count (10^9/L)',title='Last recorded platelet count before live discharge',subtitle='Descriptive post-treatment subsets; a single high count is not confirmed recovery.')
for(nm in c('Fig1A_first_event_stack','Fig1B_recovery_CIF','Fig1C_last_discharge_count')){
 pp<-list(Fig1A_first_event_stack=p1,Fig1B_recovery_CIF=p2,Fig1C_last_discharge_count=p3)[[nm]]
 ggsave(file.path(fdir,paste0(nm,'.pdf')),pp,width=9.2,height=6.2);ggsave(file.path(fdir,paste0(nm,'.png')),pp,width=9.2,height=6.2,dpi=180)
}
report<-list(started_utc=start,completed_utc=format(Sys.time(),tz='UTC',usetz=TRUE),plan_version='0.6.1',
 stage='First exploratory clinical association; no causal efficacy/safety conclusion',n_LZD=sum(d$A==1),n_VAN=sum(d$A==0),
 bootstrap_valid=sum(complete.cases(boot)),bootstrap_requested=nrow(boot),AJ_vs_empirical_CIF_max_error=max(unlist(aj_checks)),
 estimates=tab,interpretation=c('RD is LZD minus VAN; lower recovery can reflect biological, indication, discharge and monitoring mechanisms.',
 'CIF area difference is not days alive with maintained recovery.',
 'Weights use recorded baseline covariates and median-fill indicators. Missing information, common indication and source coverage remain unresolved.',
 'Last-count plots condition on live discharge and are descriptive, not a causal subgroup contrast.',
 'Transfusion-restricted and ICU-restricted sensitivity analyses are pending.'))
write_json(report,file.path(R,paste0('reports/CLINICAL_ASSOCIATION_',V,'.json')),pretty=TRUE,auto_unbox=TRUE,digits=10)
saveRDS(list(boot=boot,boot72=boot72,curveboot=curveboot,grid=grid),file.path(R,'cache/phase61_clinical_bootstrap.rds'))
capture.output(sessionInfo(),file=file.path(R,'reports/R_sessionInfo_v0.6.1_clinical.txt'))
cat(toJSON(report,pretty=TRUE,auto_unbox=TRUE,digits=7),'\n')
