args <- commandArgs(trailingOnly=TRUE)
argval <- function(key) { i <- match(key,args); if(is.na(i)||i==length(args))stop(key); args[i+1] }
root <- normalizePath(argval('--root'),mustWork=TRUE)
stamp <- argval('--stamp')
qa_only <- '--qa-only' %in% args
Sys.setenv(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
suppressPackageStartupMessages(library(DESeq2))
stopifnot(as.character(getRversion())=='4.1.2',as.character(packageVersion('DESeq2'))=='1.34.0')
suppressPackageStartupMessages(library(jsonlite))
dir.create(file.path(root,'3_results/rp241_ad_b2e_exploratory'),recursive=TRUE,showWarnings=FALSE)
out <- file.path(root,'3_results/rp241_ad_b2e_exploratory')
write_json_atomic <- function(x,path) {
  writeLines(jsonlite::toJSON(x,pretty=TRUE,auto_unbox=TRUE,na='null',digits=16),path,useBytes=TRUE)
}
write_table <- function(x,path) write.table(x,path,sep='\t',quote=FALSE,row.names=FALSE,na='NA')

# Synthetic QA checks orientation and family accounting before human effects.
if(qa_only) {
  set.seed(20261013)
  n <- 40L; g <- 200L
  diagnosis <- factor(rep(c('Control','AD'),each=n/2),levels=c('Control','AD'))
  meta <- data.frame(age_centered_per_10y=rep(seq(-1,1,length.out=20),2),sex_M=rep(c(0,1),20),diagnosis=diagnosis)
  mu <- matrix(rep(seq(70,250,length.out=g),n),nrow=g)
  mu[1:5,diagnosis=='AD'] <- mu[1:5,diagnosis=='AD']*4
  mu[6:10,diagnosis=='AD'] <- mu[6:10,diagnosis=='AD']/4
  k <- matrix(rnbinom(g*n,mu=as.vector(mu),size=20),g,n)
  storage.mode(k) <- 'integer'
  rownames(k)<-paste0('QA_',seq_len(g));colnames(k)<-paste0('QA_sample_',seq_len(n));rownames(meta)<-colnames(k)
  d <- DESeqDataSetFromMatrix(k,meta,~age_centered_per_10y+sex_M+diagnosis)
  d <- DESeq(d,test='Wald',fitType='parametric',betaPrior=FALSE,minReplicatesForReplace=Inf,parallel=FALSE,quiet=TRUE)
  r <- results(d,contrast=c('diagnosis','AD','Control'),independentFiltering=FALSE)
  reverse <- results(d,contrast=c('diagnosis','Control','AD'),independentFiltering=FALSE)
  ok <- is.finite(r$log2FoldChange)&is.finite(reverse$log2FoldChange)
  stopifnot(qr(model.matrix(~age_centered_per_10y+sex_M+diagnosis,meta))$rank==4,
            all(r$log2FoldChange[1:5]>0),all(r$log2FoldChange[6:10]<0),
            max(abs(r$log2FoldChange[ok]+reverse$log2FoldChange[ok]))<1e-8,
            identical(counts(d),k))
  bh <- rep(NA_real_,g); p <- is.finite(r$pvalue);bh[p]<-p.adjust(r$pvalue[p],method='BH',n=g)
  stopifnot(all(bh[p]>=r$pvalue[p]),all(bh[p]<=1))
  receipt <- list(status='PASS',R_version=as.character(getRversion()),DESeq2_version=as.character(packageVersion('DESeq2')),
                  genes=g,samples=n,design_rank=4,positive_fixture_genes=5,negative_fixture_genes=5,
                  contrast_swap_sign_verified=TRUE,count_replacement=FALSE,full_BH_family=g,
                  actual_dispersion_fit=attr(dispersionFunction(d),'fitType'),human_data_effects_tested=FALSE,
                  completed_at=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'))
  write_json_atomic(receipt,file.path(root,paste0('0_admin/validation/RP241_AD_B2E_Synthetic_Statistical_QA_',stamp,'.json')))
  cat('SYNTHETIC_COUNT_AWARE_QA_PASS\n');quit(status=0)
}

qa_path<-file.path(root,paste0('0_admin/validation/RP241_AD_B2E_Synthetic_Statistical_QA_',stamp,'.json'))
stopifnot(file.exists(qa_path),jsonlite::fromJSON(qa_path)$status=='PASS')
matrix_file <- file.path(root,'3_results/rp241_ad_b1_source_use_gate/GSE125583_Entrez_Raw_Integer_Counts.tsv.gz')
raw<-read.delim(gzfile(matrix_file),check.names=FALSE,colClasses=c('character',rep('integer',289)))
stopifnot(nrow(raw)==30727,ncol(raw)==290,!anyDuplicated(raw$Entrez_ID),all(as.matrix(raw[,-1])>=0))
all_counts<-as.matrix(raw[,-1]);rownames(all_counts)<-raw$Entrez_ID;rm(raw)
universe_file<-file.path(root,paste0('3_results/rp241_ad_b2_annotation_model/GSE125583_PreEffect_Primary_Measurable_Universe_',stamp,'.tsv'))
universe<-read.delim(universe_file,colClasses='character',check.names=FALSE)
gene_ids<-universe$original_Entrez_ID
stopifnot(length(gene_ids)==19108,!anyDuplicated(gene_ids),all(gene_ids %in% rownames(all_counts)))
counts19<-all_counts[gene_ids,,drop=FALSE];rm(all_counts);gc()
models<-c('M0_CLINICAL_STUDY_LABEL_PRIMARY','S1_LOW_TANGLE_CONTROL_VS_HIGH_TANGLE_AD','S2_COMMON_OBSERVED_AGE_RANGE','S3_REMOVE_EXPLICIT_GSE95587_OVERLAPPING_DONORS')
expected <- list(c(218,70),c(193,13),c(201,70),c(157,42))
summaries<-list()
for(j in seq_along(models)) {
 name<-models[j]
 mf<-file.path(root,paste0('3_results/rp241_ad_b1_source_use_gate/GSE125583_',name,'_Design_Matrix_2026-10-13.tsv'))
 frozen<-read.delim(mf,check.names=FALSE)
 m<-data.frame(age_centered_per_10y=frozen$age_centered_per_10y,sex_M=frozen$sex_M,
               diagnosis=factor(ifelse(frozen$diagnosis_AD==1,'AD','Control'),levels=c('Control','AD')))
 rownames(m)<-frozen$gsm;k<-counts19[,frozen$gsm,drop=FALSE]
 stopifnot(sum(m$diagnosis=='AD')==expected[[j]][1],sum(m$diagnosis=='Control')==expected[[j]][2],
           qr(model.matrix(~age_centered_per_10y+sex_M+diagnosis,m))$rank==4,nrow(k)==19108,
           identical(rownames(k),gene_ids),!anyNA(k))
 cat(format(Sys.time()),'START',name,ncol(k),'libraries\n');flush.console()
 d<-DESeqDataSetFromMatrix(k,m,~age_centered_per_10y+sex_M+diagnosis)
 d<-DESeq(d,test='Wald',fitType='parametric',betaPrior=FALSE,minReplicatesForReplace=Inf,parallel=FALSE,quiet=FALSE)
 actual_fit<-attr(dispersionFunction(d),'fitType')
 if(!identical(actual_fit,'parametric')) {
   write_json_atomic(list(status='HOLD_DISPERSION_FIT_TECHNICAL_AUDIT_REQUIRED',model=name,actual_fit=actual_fit),file.path(out,paste0(name,'_Technical_Hold.json')))
   stop('Frozen parametric fit changed; no effect output published')
 }
 stopifnot(identical(counts(d),k))
 r<-results(d,contrast=c('diagnosis','AD','Control'),alpha=0.05,independentFiltering=FALSE)
 stopifnot(identical(rownames(r),gene_ids))
 p_ok<-is.finite(r$pvalue);bh<-rep(NA_real_,19108);bh[p_ok]<-p.adjust(r$pvalue[p_ok],method='BH',n=19108)
 maxcook<-mcols(d)$maxCooks
 if(is.null(maxcook))maxcook<-rep(NA_real_,19108)
 reason<-ifelse(p_ok,'TESTED_FINITE_P',ifelse(r$baseMean==0,'ALL_ZERO_IN_THIS_FROZEN_SENSITIVITY','P_MISSING_COOKS_OR_MODEL_DIAGNOSTIC'))
 rr<-data.frame(original_Entrez_ID=gene_ids,baseMean=r$baseMean,log2FoldChange_AD_minus_Control=r$log2FoldChange,
                  lfcSE=r$lfcSE,Wald_statistic=r$stat,pvalue=r$pvalue,BH_FDR_full19108=bh,
                  CI95_lower_log2=r$log2FoldChange-qnorm(.975)*r$lfcSE,
                  CI95_upper_log2=r$log2FoldChange+qnorm(.975)*r$lfcSE,
                  maxCooks=maxcook,p_value_status=reason,
                  effect_direction=ifelse(is.na(r$log2FoldChange),'NOT_ESTIMABLE',ifelse(r$log2FoldChange>0,'AD_UP',ifelse(r$log2FoldChange<0,'AD_DOWN','ZERO'))),
                  HGNC_ID=universe$HGNC_ID,HGNC_symbol_frozen=universe$HGNC_symbol_frozen,
                  Ensembl_ID_frozen=universe$HGNC_Ensembl_gene_ID,GENCODE27_ID=universe$GENCODE27_versioned_ID,
                  strict_cross_disease_identity_eligible=universe$strict_cross_disease_identity_eligible,
                  evidence_scope='SINGLE_SOURCE_AD_TISSUE_EXPLORATORY__NO_AD_GBM_ANTAGONISTIC_CLASSIFICATION')
 write_table(rr,file.path(out,paste0(name,'_All19108_Results_',stamp,'.tsv')))
 sf<-data.frame(gsm=colnames(k),size_factor=sizeFactors(d),donor=frozen$donor,diagnosis=as.character(m$diagnosis))
 write_table(sf,file.path(out,paste0(name,'_Size_Factors_',stamp,'.tsv')))
 # Retain the fitted object for diagnostics without placing it in a lightweight ZIP.
 saveRDS(d,file.path(out,paste0(name,'_DESeq2_Model_',stamp,'.rds')),compress=TRUE)
 sig<-is.finite(bh)&bh<.05
 summaries[[j]]<-data.frame(model=name,n_AD=expected[[j]][1],n_control=expected[[j]][2],design_rank=4,
                    genes_input=19108,genes_finite_P=sum(p_ok),P_missing=sum(!p_ok),BH_family=19108,
                    FDR_lt_005=sum(sig),AD_UP_FDR_lt_005=sum(sig&r$log2FoldChange>0,na.rm=TRUE),
                    AD_DOWN_FDR_lt_005=sum(sig&r$log2FoldChange<0,na.rm=TRUE),actual_dispersion_fit=actual_fit,
                    count_replacement=FALSE,automatic_sample_exclusion=FALSE,evidence_scope='EXPLORATORY_ONLY')
 write_table(do.call(rbind,summaries),file.path(out,paste0('GSE125583_Exploratory_Models_Execution_Summary_',stamp,'.tsv')))
 cat(format(Sys.time()),'COMPLETE',name,'\n');flush.console();rm(d,r,k);gc()
}
write_json_atomic(list(status='FOUR_FROZEN_COUNT_AWARE_MODELS_COMPLETE__SINGLE_SOURCE_EXPLORATORY_ONLY',
                models=4,genes_per_model=19108,DESeq2_version=as.character(packageVersion('DESeq2')),
                independent_replication_admitted=FALSE,AD_GBM_opposite_gene_classification_unlocked=FALSE,
                completed_at=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z')),
                file.path(out,paste0('Execution_Receipt_',stamp,'.json')))
writeLines(capture.output(sessionInfo()),file.path(out,paste0('R_sessionInfo_',stamp,'.txt')))
