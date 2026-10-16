# B9A conditional limited historical compartment validation.
# No formal effect can be evaluated without the new-source S6 receipt and freeze.
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else '/mnt/h/SCI2/MRQI'
.libPaths(c(file.path(root,'0_admin/environment/rp241_b9_microarray/R_library'),.libPaths()))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(jsonlite))
date <- '2026-10-14'
paired_fit <- function(y, m) {
    stopifnot(all(is.finite(y)),!anyDuplicated(rownames(y)),identical(colnames(y),m$analysis_sample_id))
    m$patient <- factor(m$donor_key)
    m$compartment <- factor(m$role,levels=c('HOST_BAT','TUMOR'))
    design <- model.matrix(~patient+compartment,m)
    stopifnot(qr(design)$rank==ncol(design),nrow(design)>ncol(design))
    fit <- eBayes(lmFit(y,design),trend=FALSE,robust=TRUE)
    co <- which(colnames(design)=='compartmentTUMOR');stopifnot(length(co)==1)
    se <- fit$stdev.unscaled[,co]*sqrt(fit$s2.post)
    dt <- rep_len(fit$df.total,nrow(y))
    b <- fit$coefficients[,co];tt <- fit$t[,co]
    table <- data.frame(HGNC_ID=rownames(y),beta_TUMOR_minus_HOST=b,moderated_SE=se,
             moderated_t=tt,df_total=dt,source_raw_residual_df=fit$df.residual,
             P_two_sided=fit$p.value[,co],CI95_lower=b-qt(.975,dt)*se,CI95_upper=b+qt(.975,dt)*se,
             source_BH_full12405=p.adjust(fit$p.value[,co],method='BH'))
    list(table=table,fit=fit,design=design,metadata=m)
}
validate81 <- function(table, directions) {
    stopifnot(nrow(directions)==81,!anyDuplicated(directions$HGNC_ID),all(directions$expected_GBM_direction %in% c('UP','DOWN')))
    index <- match(directions$HGNC_ID,table$HGNC_ID)
    out <- directions
    for(col in setdiff(colnames(table),'HGNC_ID'))out[[col]]<-table[[col]][index]
    expected <- ifelse(out$expected_GBM_direction=='UP',1,-1)
    out$P_one_sided <- pt(expected*out$moderated_t,df=out$df_total,lower.tail=FALSE)
    out$BH_full81 <- p.adjust(out$P_one_sided,method='BH',n=81)
    out$BY_full81 <- p.adjust(out$P_one_sided,method='BY',n=81)
    out$direction_correct <- expected*out$beta_TUMOR_minus_HOST>0
    out$status <- ifelse(is.na(out$P_one_sided),'NOT_TESTABLE_NO_UNIQUE_SOURCE_IDENTITY',
                 ifelse(!out$direction_correct,'POINT_DIRECTION_OPPOSED',
                 ifelse(out$BH_full81<.05,'LIMITED_COMPARTMENT_DIRECTION_SUPPORTED_BH','DIRECTION_CONCORDANT_NOT_FDR_SUPPORTED')))
    out$family_denominator <- 81
    out$population_GBM_confirmation <- FALSE
    out$causal_or_Tier1 <- FALSE
    out
}
aggregate_host <- function(gx, m, variant='MEAN_ALL_REGIONS') {
    agg<-list();md<-list()
    for(donor in sort(unique(m$donor_key))) for(role in c('HOST_BAT','TUMOR')) {
        ix<-which(m$donor_key==donor & m$role==role)
        if(role=='HOST_BAT' && variant=='FIRST_NON_R') ix<-ix[!grepl('R$',m$source_sample_suffix[ix])]
        if(role=='HOST_BAT' && variant=='ALTERNATE_R' && any(grepl('R$',m$source_sample_suffix[ix]))) ix<-ix[grepl('R$',m$source_sample_suffix[ix])]
        stopifnot(length(ix)>=1)
        key<-paste(donor,role,sep='__')
        agg[[key]]<-rowMeans(gx[,ix,drop=FALSE])
        md[[key]]<-data.frame(analysis_sample_id=key,donor_key=donor,role=role,source_gsms=paste(m$gsm[ix],collapse=';'),source_arrays=length(ix))
    }
    mm<-do.call(rbind,md);rownames(mm)<-NULL
    list(y=do.call(cbind,agg),m=mm)
}

if('--synthetic' %in% args) {
    # Real source expression and candidate gene rows are not opened in this branch.
    set.seed(20261014)
    m<-data.frame(analysis_sample_id=paste0('s',1:10),donor_key=rep(paste0('p',1:5),each=2),role=rep(c('HOST_BAT','TUMOR'),5))
    y<-matrix(rnorm(1205*10),1205,10,dimnames=list(paste0('TEST_GENE',1:1205),m$analysis_sample_id))
    known<-seq(-1,1,length.out=1205)
    y[,m$role=='TUMOR']<-y[,m$role=='TUMOR']+known
    result<-paired_fit(y,m)
    # Independent direct donor difference formula, not a duplicate lmFit.
    diff<-y[,m$role=='TUMOR']-y[,m$role=='HOST_BAT']
    stopifnot(max(abs(result$table$beta_TUMOR_minus_HOST-rowMeans(diff)))<1e-12)
    stopifnot(max(abs(result$table$P_two_sided-2*pt(abs(result$table$moderated_t),result$table$df_total,lower.tail=FALSE)))<1e-12)
    dir<-data.frame(HGNC_ID=c(paste0('TEST_GENE',1:60),paste0('MISSING',1:21)),expected_GBM_direction=rep(c('UP','DOWN'),length.out=81))
    out<-validate81(result$table,dir)
    stopifnot(nrow(out)==81,sum(is.na(out$P_one_sided))==21,all(out$family_denominator==81),
      all(out$P_one_sided[out$direction_correct%in%FALSE]>.5))
    # Independent complete-family BH calculation with missing mapped to1.
    p<-out$P_one_sided;p[is.na(p)]<-1;ord<-order(p);q<-pmin(1,rev(cummin(rev(p[ord]*81/(1:81)))))
    qb<-numeric(81);qb[ord]<-q
    stopifnot(max(abs(qb[!is.na(out$P_one_sided)]-out$BH_full81[!is.na(out$P_one_sided)]))<1e-12)
    fake<-data.frame(gsm=paste0('a',1:12),donor_key=c(rep('p1',2),rep('p2',2),rep('p3',3),rep('p4',3),rep('p5',2)),
        role=c('HOST_BAT','TUMOR','HOST_BAT','TUMOR','HOST_BAT','HOST_BAT','TUMOR','HOST_BAT','HOST_BAT','TUMOR','HOST_BAT','TUMOR'),
        source_sample_suffix=c('1','1','2','2','3','3R','3','4','4R','4','5','5'))
    fx<-matrix(rep(1:12,each=1205),1205,12,dimnames=list(rownames(y),fake$gsm))
    a<-aggregate_host(fx,fake);a1<-aggregate_host(fx,fake,'FIRST_NON_R');a2<-aggregate_host(fx,fake,'ALTERNATE_R')
    stopifnot(a$y[1,'p3__HOST_BAT']==5.5,a1$y[1,'p3__HOST_BAT']==5,a2$y[1,'p3__HOST_BAT']==6,ncol(a$y)==10)
    path<-file.path(root,'0_admin/validation',paste0('RP241_AD_GBM_B9A_Synthetic_Statistical_QA_',date,'.json'))
    stopifnot(!file.exists(path))
    write_json(list(status='PASS',real_expression_opened=FALSE,seed=20261014,source_genes=1205,paired_patients=5,
        candidates=81,tested_identity=60,missing_retained=21,checks=c('Independent paired mean beta','Independent t tail','Full81 BH including missing','Opposite sign one-sided P above0.5','Mean/nonR/R host nested aggregation')),
        path,pretty=TRUE,auto_unbox=TRUE)
    cat('Synthetic B9A statistical QA PASS; no formal source contrast.\n');quit(status=0)
}

admission_path<-file.path(root,'0_admin/decision',paste0('RP241_AD_GBM_B9A_GSE13276_Source_Admission_',date,'.json'))
if(!file.exists(admission_path)) {
    cat('HOLD_B9A_NEW_SOURCE_S6_REVIEW_REQUIRED__NO_CANDIDATE_CONTRAST_OPENED\n');quit(status=3)
}
admission<-fromJSON(admission_path)
stopifnot(identical(admission$source,'GSE13276'),identical(admission$status,'GO_LIMITED_HISTORICAL_PAIRED_COMPARTMENT_VALIDATION'),
    identical(admission$evidence_kind,'USER_REPORTED_COAUTHOR_SOURCE_REVIEW'),nzchar(admission$reviewer_name),nzchar(admission$review_date),
    isTRUE(admission$accepted_5_patient_paired_scope),isTRUE(admission$accepted_excluded_cavernoma_controls),
    isTRUE(admission$accepted_GEO_R3_R4_nesting),!isTRUE(admission$population_normal_claim))
# The PowerShell launcher verifies the frozen bytes before invoking R.
qc<-file.path(root,'3_results/rp241_ad_gbm_b9_blind_qc');src<-file.path(root,'3_results/rp241_ad_gbm_b9_source_audit')
gx<-readRDS(file.path(qc,paste0('GSE13276_RMA12_All12405_UniqueHGNC_Expression_',date,'.rds')))
m<-read.delim(file.path(src,paste0('GSE13276_Primary12_Metadata_',date,'.tsv')),stringsAsFactors=FALSE,check.names=FALSE)
d<-read.delim(file.path(src,paste0('All81_GSE13276_PreEffect_Identity_Missing_Ledger_',date,'.tsv')),stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(gx)==12405,ncol(gx)==12,nrow(d)==81,identical(colnames(gx),m$gsm),sum(d$unique_probe_n>0)==60)
outdir<-file.path(root,'3_results/rp241_ad_gbm_b9a_limited_validation');stopifnot(!dir.exists(outdir));dir.create(outdir)
base<-aggregate_host(gx,m)
models<-list(M0_ALL5_HOST_MEAN=base,S1_FIRST_NON_R=aggregate_host(gx,m,'FIRST_NON_R'),S2_ALTERNATE_R=aggregate_host(gx,m,'ALTERNATE_R'))
for(patient in sort(unique(base$m$donor_key))) {
    ok<-base$m$donor_key!=patient
    models[[paste0('LOPO_',patient)]]<-list(y=base$y[,ok,drop=FALSE],m=base$m[ok,,drop=FALSE])
}
stopifnot(length(models)==8)
summaries<-list()
for(name in names(models)) {
    fit<-paired_fit(models[[name]]$y,models[[name]]$m);r<-validate81(fit$table,d)
    write.table(fit$table,file.path(outdir,paste0(name,'_All12405_Source_Diagnostics.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
    write.table(r,file.path(outdir,paste0(name,'_All81_Frozen_Direction_Results.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
    write.table(fit$design,file.path(outdir,paste0(name,'_Design.tsv')),sep='\t',col.names=NA,quote=FALSE)
    saveRDS(fit,file.path(outdir,paste0(name,'_limma_Fit.rds')))
    summaries[[name]]<-data.frame(model=name,source_family=12405,frozen_candidate_family=81,identity_tested=sum(!is.na(r$P_one_sided)),missing=sum(is.na(r$P_one_sided)),
        supported_BH=sum(r$status=='LIMITED_COMPARTMENT_DIRECTION_SUPPORTED_BH'),supported_BY=sum(r$BY_full81<.05 & r$direction_correct,na.rm=TRUE),
        point_direction_opposed=sum(r$status=='POINT_DIRECTION_OPPOSED'),paired_patients=length(unique(fit$metadata$donor_key)),residual_df=nrow(fit$design)-ncol(fit$design),population_confirmation_complete=FALSE)
}
s<-do.call(rbind,summaries);rownames(s)<-NULL
write.table(s,file.path(outdir,paste0('All8_Models_Execution_Summary_',date,'.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
write_json(list(phase='B9A_LIMITED_HISTORICAL_PAIRED_COMPARTMENT',date=date,primary_supported=s$supported_BH[1],
    population_GBM_confirmation_complete=FALSE,modern_IDHwt_confirmation_complete=FALSE,new_Tier1_causal_targets=0,
    statistical_completion_not_source_claim_upgrade=TRUE,review_receipt=admission_path),file.path(outdir,paste0('Execution_Receipt_',date,'.json')),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(outdir,paste0('R_sessionInfo_',date,'.txt')))
cat('All8 B9A models completed; independent adjudication required before claims.\n')
