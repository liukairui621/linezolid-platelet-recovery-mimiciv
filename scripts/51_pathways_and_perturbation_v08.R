# Frozen eight-program test; existing outcome/DE results informed the follow-up plan.
suppressPackageStartupMessages({library(edgeR);library(limma);library(data.table);library(jsonlite)})
R<-'/root/projects/linezolid_platelet_recovery';O<-file.path(R,'public_data/pathways_v0.8');dir.create(O,recursive=TRUE,showWarnings=FALSE)
cfg<-fromJSON(file.path(R,'config/pathway_sets_v0.8.json'),simplifyVector=FALSE)$sets
setids<-vapply(cfg,function(z)z$id,'');sets<-setNames(lapply(cfg,function(z)unlist(z$genes)),setids)
setnames<-setNames(vapply(cfg,function(z)z$name,''),setids)
checks<-list();allres<-list();coverage<-list();gene_sets<-list()
testsets<-function(v,design,coef,symbol,dataset,branch){
 stopifnot(length(symbol)==nrow(v),ncol(design)<nrow(design),qr(design)$rank==ncol(design))
 avg<-rowMeans(v$E);ok<-!is.na(symbol)&nzchar(symbol)
 oo<-order(!ok,-avg,rownames(v));oo<-oo[ok[oo]];oo<-oo[!duplicated(symbol[oo])]
 vv<-v[oo,];sym<-symbol[oo];ix<-lapply(sets,function(s)which(sym%in%s));eligible<-lengths(ix)>=10
 cov<-data.table(dataset=dataset,branch=branch,set_id=setids,set_name=unname(setnames),deposited_N=lengths(sets),expressed_N=lengths(ix),background_N=length(sym),testable=eligible)
 coverage[[paste(dataset,branch,sep='__')]]<<-cov
 if(!any(eligible))stop('No testable frozen sets')
 cam<-camera(vv,index=ix[eligible],design=design,contrast=coef,inter.gene.cor=NA,allow.neg.cor=FALSE,sort=FALSE)
 fixed<-camera(vv,index=ix[eligible],design=design,contrast=coef,inter.gene.cor=0.01,allow.neg.cor=FALSE,sort=FALSE)
 fr<-fry(vv,index=ix[eligible],design=design,contrast=coef,sort=FALSE)
 z<-copy(cov);z[,`:=`(camera_direction=NA_character_,camera_P=1,camera_correlation=NA_real_,fixed01_P=1,fry_direction=NA_character_,fry_P=1,fry_mixed_P=1,mean_log2FC=NA_real_,mean_moderated_t=NA_real_)]
 ff<-eBayes(lmFit(vv,design),robust=TRUE,trend=isTRUE(attr(v,'trend')))
 for(k in which(eligible)){
  nm<-setids[k];z[k,`:=`(camera_direction=cam[nm,'Direction'],camera_P=cam[nm,'PValue'],camera_correlation=cam[nm,'Correlation'],fixed01_P=fixed[nm,'PValue'],fry_direction=fr[nm,'Direction'],fry_P=fr[nm,'PValue'],fry_mixed_P=fr[nm,'PValue.Mixed'],mean_log2FC=mean(ff$coefficients[ix[[k]],coef]),mean_moderated_t=mean(ff$t[ix[[k]],coef]))]
 }
 z[,`:=`(camera_FDR=p.adjust(camera_P,'BH'),fixed01_FDR=p.adjust(fixed01_P,'BH'),fry_FDR=p.adjust(fry_P,'BH'),fry_mixed_FDR=p.adjust(fry_mixed_P,'BH'))]
 allres[[paste(dataset,branch,sep='__')]]<<-z
 if(!grepl('^LOO_',branch)){
  gene_sets[[paste(dataset,branch,sep='__')]]<<-rbindlist(lapply(which(eligible),function(k){jj<-ix[[k]];data.table(dataset=dataset,branch=branch,set_id=setids[k],gene_id=rownames(vv)[jj],symbol=sym[jj],logFC=ff$coefficients[jj,coef],moderated_t=ff$t[jj,coef])}))
 }
 z
}
# Original disease contrast, reconstructed without modifying archived files.
P<-file.path(R,'public_data/GSE252275_processed');dt<-fread(file.path(P,'total_exon_counts.tsv.gz'));ids<-dt[[1]]
ss<-as.data.frame(fread(file.path(P,'samples_and_qc.tsv')));ss<-ss[ss$group%in%c('sepsis_preserve','sepsis_decrease'),]
ss$group<-factor(ss$group,levels=c('sepsis_preserve','sepsis_decrease'));cts<-as.matrix(dt[,ss$sample,with=FALSE]);rownames(cts)<-ids
ann<-fread(file.path(P,'gene_annotations.tsv'));symbol<-ann$symbol[match(ids,ann$gene_id_as_supplied)]
globins<-c('HBA1','HBA2','HBB','HBD','HBE1','HBG1','HBG2','HBM','HBQ1','HBZ')
ss$globin_fraction<-colSums(cts[symbol%in%globins,,drop=FALSE])/colSums(cts)
ss$globin_logit<-qlogis(pmax(1e-6,pmin(1-1e-6,ss$globin_fraction)));flag<-ss$globin_fraction>.2
keep<-rowSums(t(t(cts)/colSums(cts))*1e6>=1)>=9;stopifnot(sum(keep)==11681)
dbfit<-function(ix,adjust=FALSE){
 d<-ss[ix,,drop=FALSE];d$globin_centered<-d$globin_logit-mean(d$globin_logit)
 X<-model.matrix(if(adjust)~group+globin_centered else ~group,d)
 y<-calcNormFactors(DGEList(cts[,ix,drop=FALSE]),method='TMM');y<-y[keep,,keep.lib.sizes=TRUE]
 v<-voom(y,X,plot=FALSE,normalize.method='none');list(v=v,X=X)
}
specs<-list(all20_group=list(ix=1:20,adjust=FALSE),all20_group_globin=list(ix=1:20,adjust=TRUE),exclude3_group=list(ix=which(!flag),adjust=FALSE),exclude3_group_globin=list(ix=which(!flag),adjust=TRUE))
for(nm in names(specs)){
 z<-do.call(dbfit,specs[[nm]]);tt<-topTable(eBayes(lmFit(z$v,z$X),robust=TRUE),coef=2,number=Inf,sort.by='none')
 old<-fread(file.path(R,'public_data/GSE252275_DE_v0.5',paste0(nm,'.tsv.gz')))
 stopifnot(identical(rownames(tt),old$gene_id));err<-max(abs(tt$logFC-old$logFC),abs(tt$t-old$t));stopifnot(err<1e-8)
 checks[[nm]]<-list(DE_max_error=err,genes=nrow(tt));testsets(z$v,z$X,2,symbol[keep],'GSE252275',nm)
 cat('Disease',nm,'complete\n');flush.console()
}
for(i in 1:20){z<-dbfit(setdiff(1:20,i));testsets(z$v,z$X,2,symbol[keep],'GSE252275',paste0('LOO_',ss$sample[i]));if(i%%5==0){cat('LOO',i,'/20\n');flush.console()}}
# Approved HGNC exact Ensembl mapping; ambiguous Ensembl mappings are excluded.
h<-fread(file.path(R,'public_data/reference_sources_v0.8/hgnc_complete_set.txt'),select=c('symbol','status','ensembl_gene_id'))
h<-unique(h[status=='Approved' & !is.na(ensembl_gene_id) & nzchar(ensembl_gene_id),.(symbol,ensembl_gene_id)])
amb<-h[,.(n=uniqueN(symbol)),by=ensembl_gene_id][n>1]$ensembl_gene_id;h<-h[!ensembl_gene_id%in%amb]
map_symbol<-function(ids)h$symbol[match(sub('\\.[0-9]+$','',ids),h$ensembl_gene_id)]
drugreports<-list()
for(cell in c('fibroblast','macrophage')){
 path<-file.path(R,'public_data',if(cell=='fibroblast')'GSE310202_fibroblast_raw_count_data.txt.gz' else 'GSE310202_macrophage_raw_count_data.tabular.gz')
 dd<-read.delim(gzfile(path),row.names=1,check.names=FALSE);cc<-as.matrix(dd);storage.mode(cc)<-'double';stopifnot(ncol(cc)==6,all(is.finite(cc)),all(cc>=0),!anyDuplicated(rownames(cc)))
 sm<-data.frame(sample=colnames(cc),LZD=as.integer(grepl('^TGFL',colnames(cc))),block=sub('^TGFL?','',colnames(cc)))
 sm$block<-factor(sm$block);stopifnot(all(table(sm$block,sm$LZD)==1))
 sy<-map_symbol(rownames(cc));k<-rowSums(t(t(cc)/colSums(cc))*1e6>=1)>=3
 y<-calcNormFactors(DGEList(cc),method='TMM');y<-y[k,,keep.lib.sizes=TRUE]
 for(branch in c('blocked_voom','unblocked_voom','blocked_trend')){
  X<-model.matrix(if(branch=='unblocked_voom')~LZD else ~block+LZD,sm);coef<-match('LZD',colnames(X))
  if(branch=='blocked_trend'){v<-new('EList',list(E=cpm(y,log=TRUE,prior.count=.5),design=X));attr(v,'trend')<-TRUE}
  else v<-voom(y,X,plot=FALSE,normalize.method='none')
  fit<-eBayes(lmFit(v,X),robust=TRUE,trend=branch=='blocked_trend');tt<-topTable(fit,coef=coef,number=Inf,sort.by='none')
  fwrite(data.table(gene_id=rownames(tt),symbol=sy[k],tt),file.path(O,paste0('GSE310202_',cell,'_',branch,'_DE.tsv.gz')),sep='\t',compress='gzip')
  drugreports[[paste(cell,branch)]]<-list(n_samples=6,blocks=3,residual_df=nrow(X)-ncol(X),genes=sum(k),FDR05=sum(tt$adj.P.Val<.05),fractional_values=sum(cc!=floor(cc)),mapped=sum(!is.na(sy[k])))
  testsets(v,X,coef,sy[k],paste0('GSE310202_',cell),branch)
 }
 fwrite(sm,file.path(O,paste0(cell,'_sample_design.tsv')),sep='\t');cat('Drug',cell,'complete\n');flush.console()
}
# Independent sepsis platelet context, not thrombocytopenia/recovery replication.
dd<-fread(file.path(R,'public_data/reference_sources_v0.8/GSE210797_filtered_counts.txt.gz'))
samples<-c(paste('Patient',1:6),paste('Donor',1:6));stopifnot(all(samples%in%names(dd)))
cc<-as.matrix(dd[,samples,with=FALSE]);rownames(cc)<-dd$id;stopifnot(all(cc>=0),all(cc==floor(cc)),!anyDuplicated(rownames(cc)))
sm<-data.frame(sample=samples,sepsis=as.integer(grepl('Patient',samples)))
sy<-dd$gene_name;k<-rowSums(t(t(cc)/colSums(cc))*1e6>=1)>=6
sm$globin_fraction<-colSums(cc[sy%in%globins,,drop=FALSE])/colSums(cc)
sm$globin_centered<-qlogis(pmax(1e-6,pmin(1-1e-6,sm$globin_fraction)));sm$globin_centered<-sm$globin_centered-mean(sm$globin_centered)
y<-calcNormFactors(DGEList(cc),method='TMM');y<-y[k,,keep.lib.sizes=TRUE]
for(branch in c('sepsis_group','sepsis_group_globin')){
 X<-model.matrix(if(branch=='sepsis_group')~sepsis else ~sepsis+globin_centered,sm);v<-voom(y,X,plot=FALSE,normalize.method='none')
 tt<-topTable(eBayes(lmFit(v,X),robust=TRUE),coef=2,number=Inf,sort.by='none')
 fwrite(data.table(gene_id=rownames(tt),symbol=sy[k],tt),file.path(O,paste0('GSE210797_',branch,'_DE.tsv.gz')),sep='\t',compress='gzip')
 testsets(v,X,2,sy[k],'GSE210797',branch)
}
fwrite(sm,file.path(O,'GSE210797_sample_design.tsv'),sep='\t')
res<-rbindlist(allres);fwrite(res,file.path(R,'reports/PATHWAY_TESTS_v0.8.tsv'),sep='\t')
fwrite(rbindlist(coverage),file.path(R,'reports/PATHWAY_COVERAGE_v0.8.tsv'),sep='\t')
fwrite(rbindlist(gene_sets),file.path(O,'pathway_member_statistics.tsv.gz'),sep='\t',compress='gzip')
loo<-res[dataset=='GSE252275' & grepl('^LOO_',branch)]
ref<-res[dataset=='GSE252275' & branch=='all20_group']
ls<-loo[,.(LOO_count=.N,LOO_camera_FDR05=sum(camera_FDR<.05),LOO_fixed01_FDR05=sum(fixed01_FDR<.05),LOO_direction_Up=sum(camera_direction=='Up'),LOO_min_FDR=min(camera_FDR),LOO_max_FDR=max(camera_FDR)),by=set_id]
fwrite(ls,file.path(R,'reports/PATHWAY_LOO_v0.8.tsv'),sep='\t')
report<-list(completed_utc=format(Sys.time(),tz='UTC',usetz=TRUE),disease_reproduction=checks,drug=drugreports,
 method='Competitive camera estimated nonnegative intergene correlation; fixed0.01 sensitivity; self-contained fry. All8 sets in each BH family; untestable P=1.',
 independent=list(dataset='GSE210797',n=12,contrast='sepsis6 minus healthy6',supplied_genes=nrow(cc),tested_genes=sum(k),upstream_filtered=TRUE,thrombocytopenia_validation=FALSE),
 limits=c('Frozen after prior gene-level diagnostics. Exploratory, not independently preregistered.', 'mRNA direction is not translation flux or mitochondrial activity.', 'Drug cells are skin fibroblasts and iPSC macrophages under TGF stimulation, not platelets or megakaryocytes.', 'Drug suffix pairing verified against filenames/GEO replicate labels; no PCSS analyzed.', 'Self-contained and competitive tests answer different questions; their significance is not interchangeable.', 'Independent sepsis-versus-healthy study does not validate a thrombocytopenia or drug-specific interaction.'))
write_json(report,file.path(R,'reports/OMICS_BRIDGE_v0.8.json'),pretty=TRUE,auto_unbox=TRUE,digits=10)
capture.output(sessionInfo(),file=file.path(R,'reports/R_sessionInfo_v0.8_omics.txt'))
print(res[!grepl('LOO_',branch),.(dataset,branch,set_id,expressed_N,camera_direction,camera_FDR,fixed01_FDR,fry_FDR,mean_log2FC)])
