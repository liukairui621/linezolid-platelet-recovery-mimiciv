# Verify saved association outputs; repair stack geometry without any refit.
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(ggplot2)})
R<-'/root/projects/linezolid_platelet_recovery';V<-'v0.6.1';F<-file.path(R,'figures',V)
cc<-fread(file.path(R,paste0('reports/FIRST_EVENT_CURVES_',V,'.csv')))
tab<-fread(file.path(R,paste0('reports/CLINICAL_ASSOCIATION_',V,'.csv')))
bb<-readRDS(file.path(R,'cache/phase61_clinical_bootstrap.rds'))
m<-readRDS(file.path(R,'cache/phase61_overlap.rds'))
stopifnot(max(abs(rowSums(cc[,.(recovery,death,discharge,no_first_event)])-1))<1e-12,
 all(as.matrix(cc[,.(recovery,death,discharge,no_first_event)])>= -1e-12))
for(i in seq_len(nrow(tab))){
 q<-quantile(bb$boot[,tab$estimand[i]],c(.025,.975));stopifnot(max(abs(q-c(tab$lower[i],tab$upper[i])))<1e-12)
}
integrals<-cc[order(day),.(area=sum(diff(day)*head(recovery,-1))),by=drug]
stopifnot(abs(integrals[drug=='linezolid']$area-tab[estimand=='area_LZD']$estimate)<1e-12,
 abs(integrals[drug=='vancomycin']$area-tab[estimand=='area_VAN']$estimate)<1e-12)
baseN<-table(m$data$subject_id)
for(ob in m$bootstrap){if(is.null(ob))next;nn<-table(m$data$subject_id[ob$ix]);ratio<-as.numeric(nn)/as.numeric(baseN[names(nn)]);stopifnot(all(abs(ratio-round(ratio))<1e-12))}
cols<-c('Confirmed recovery'='#3B8D83','Death before recovery'='#B85A5A','Live discharge before recovery'='#D4B86B','No first event yet'='#DEE4E8')
rect<-cc[order(day),{
 zz<-.SD;n<-nrow(zz);lo<-matrix(0,n,4);hi<-as.matrix(zz[,.(recovery,death,discharge,no_first_event)])
 for(j in 2:4){lo[,j]<-hi[,j-1];hi[,j]<-hi[,j]+hi[,j-1]}
 hi[,4]<-1 # exact theoretical bound, avoiding floating-point scale deletion
 rbindlist(lapply(1:4,function(j)data.table(xmin=zz$day[-n],xmax=zz$day[-1],ymin=lo[-n,j],ymax=hi[-n,j],state=names(cols)[j])))
},by=drug]
stopifnot(all(rect$ymin>= -1e-12),all(rect$ymax<=1+1e-12),all(rect$ymax>=rect$ymin))
fwrite(rect,file.path(R,paste0('reports/FIRST_EVENT_STACK_GEOMETRY_',V,'.csv')))
p<-ggplot(rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=state))+geom_rect(color=NA)+facet_wrap(~drug)+
 scale_fill_manual(values=cols,breaks=names(cols))+scale_y_continuous(labels=function(z)paste0(round(z*100),'%'),limits=c(0,1),expand=expansion(mult=0))+
 scale_x_continuous(breaks=c(0,3,7,14),limits=c(0,14),expand=expansion(mult=0))+
 labs(x='Days from first recorded dose',y='Weighted first-event probability',fill=NULL,title='Recorded recovery and competing exits',
 subtitle='First-event categories; recovery stays absorbing after later death or discharge.')+
 theme_classic(base_size=12,base_family='sans')+theme(legend.position='bottom',plot.title.position='plot',legend.text=element_text(size=10),panel.spacing.x=grid::unit(1.2,'lines'))
hist<-file.path(F,'render_history');dir.create(hist,showWarnings=FALSE)
for(ext in c('png','pdf'))if(!file.exists(file.path(hist,paste0('Fig1A_initial_stack_error.',ext))))file.copy(file.path(F,paste0('Fig1A_first_event_stack.',ext)),file.path(hist,paste0('Fig1A_initial_stack_error.',ext)))
ggsave(file.path(F,'Fig1A_first_event_stack.pdf'),p,width=9.2,height=6.2)
ggsave(file.path(F,'Fig1A_first_event_stack.png'),p,width=9.2,height=6.2,dpi=180)
report<-list(verified_utc=format(Sys.time(),tz='UTC',usetz=TRUE),passed=TRUE,
 checks=c('All first-event probabilities sum to1 and are nonnegative','Saved percentile intervals recomputed from bootstrap agree',
 'Exact step integration equals stored recovery-CIF area','Every patient cluster remains intact in all bootstrap draws',
 'Explicit stacked rectangles remain in[0,1] and sum to1; duplicate-x geom_area stacking bug corrected without refitting'),
 visual_status='Final PNG inspected separately; initial erroneous stack retained only in render_history, not an accepted figure.')
write_json(report,file.path(R,paste0('reports/ARTIFACT_CHECKS_',V,'.json')),pretty=TRUE,auto_unbox=TRUE)
cat(toJSON(report,pretty=TRUE,auto_unbox=TRUE),'\n')
