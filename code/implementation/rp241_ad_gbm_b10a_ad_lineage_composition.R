# B10A executable contract. Formal mode is blocked without scope-specific S6.
suppressPackageStartupMessages({library(DESeq2);library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else '/mnt/h/SCI2/MRQI'
synthetic <- '--synthetic' %in% args
date <- '2026-10-14'
stopifnot(as.character(packageVersion('DESeq2'))=='1.34.0')
write_tsv <- function(x,p) write.table(x,p,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
fixed_fdr <- function(p,n,method='BH') {stopifnot(length(p)==n);p[!is.finite(p)]<-1;p.adjust(p,method=method,n=n)}
one_sided <- function(stat,sign) pnorm(-sign*stat)
fit_counts <- function(counts,meta,form) {
 stopifnot(all(counts>=0),all(counts==floor(counts)),identical(colnames(counts),rownames(meta)))
 X<-model.matrix(form,meta);stopifnot(qr(X)$rank==ncol(X),nrow(X)-ncol(X)>=3)
 keep<-rowSums(counts>=10)>=3
 warnings<-character()
 fit<-withCallingHandlers(DESeq(DESeqDataSetFromMatrix(counts[keep,,drop=FALSE],meta,form),
    sfType='poscounts',fitType='parametric',betaPrior=FALSE,minReplicatesForReplace=Inf,quiet=TRUE),
    warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 rs<-as.data.frame(results(fit,contrast=c('Diagnosis','AD','Control'),independentFiltering=FALSE,cooksCutoff=FALSE))
 stopifnot('betaConv' %in% colnames(mcols(fit)))
 rs$beta_converged<-mcols(fit)$betaConv
 cook<-assays(fit)[['cooks']];rs$max_raw_Cook<-if(is.null(cook)) NA_real_ else apply(cook,1,max,na.rm=TRUE)
 rs$Ensembl_stable_ID<-rownames(rs)
 list(fit=fit,result=rs,keep=keep,warning=warnings,rank=qr(X)$rank,df=nrow(X)-ncol(X),
      fit_type=attr(dispersionFunction(fit),'fitType'))
}
linear_hc3 <- function(y,meta,form) {
 X<-model.matrix(form,meta);stopifnot(qr(X)$rank==ncol(X),nrow(X)>ncol(X))
 inv<-solve(crossprod(X));beta<-inv%*%crossprod(X,y);res<-y-as.numeric(X%*%beta)
 h<-rowSums((X%*%inv)*X);V<-inv%*%crossprod(X,X*as.numeric((res/(1-h))^2))%*%inv
 j<-which(colnames(X)=='DiagnosisAD');se<-as.numeric(sqrt(V[j,j]));df<-nrow(X)-ncol(X);b<-as.numeric(beta[j,1])
 c(beta=b,SE=se,CI_lower=b-qt(.975,df)*se,CI_upper=b+qt(.975,df)*se,P_two_sided=2*pt(-abs(b/se),df),df=df)
}
if(synthetic) {
 set.seed(20261014);n<-18;g<-1500
 meta<-data.frame(Batch=factor(rep(1:3,6)),Age=seq(75,92),Sex=factor(rep(c('F','M'),9)),Diagnosis=factor(c(rep('Control',7),rep('AD',11)),levels=c('Control','AD')))
 rownames(meta)<-paste0('syn',seq_len(n));means<-exp(runif(g,3,6));mu<-outer(means,rep(1,n));mu[1,meta$Diagnosis=='AD']<-mu[1,meta$Diagnosis=='AD']*8;mu[2,meta$Diagnosis=='AD']<-mu[2,meta$Diagnosis=='AD']/8
 meta$Age<-(meta$Age-mean(meta$Age))/10
 cc<-matrix(rnbinom(g*n,mu=as.numeric(mu),size=10),g,n,dimnames=list(paste0('gene',seq_len(g)),rownames(meta)))
 a<-fit_counts(cc,meta,~Batch+Age+Sex+Diagnosis);r<-a$result
 stopifnot(r['gene1','log2FoldChange']>0,r['gene2','log2FoldChange']<0)
 p<-one_sided(c(r['gene1','stat'],r['gene2','stat']),c(1,-1))
 stopifnot(all(is.finite(p)),all(p<.05))
 pv<-c(.002,.01,rep(1,256));o<-order(pv);q<-pmin(1,rev(cummin(rev(pv[o]*258/seq_len(258)))));back<-numeric(258);back[o]<-q
 stopifnot(max(abs(back-fixed_fdr(pv,258)))<1e-12)
 hc<-linear_hc3(as.numeric(cc[1,]),meta,~Batch+Age+Sex+Diagnosis)
 stopifnot(hc['beta']>0,all(is.finite(hc)))
 out<-file.path(root,'0_admin/validation',paste0('RP241_AD_GBM_B10A_CountAware_Synthetic_QA_',date,'.json'))
 write_json(list(status='PASS',synthetic_only=TRUE,real_gene_contrasts=0,DESeq2=as.character(packageVersion('DESeq2')),
  controlled_positive_and_negative_signs=TRUE,manual_fixed258_BH_match=TRUE,HC3_finite=TRUE,fit_type=a$fit_type,warnings=a$warning),out,pretty=TRUE,auto_unbox=TRUE)
 cat('PASS synthetic count-aware/HC3/fixed-family QA; no real source effects.\n');quit(status=0)
}

# Read attestation and frozen SHA list BEFORE actual expression data.
approval<-file.path(root,'0_admin/audit',paste0('RP241_AD_GBM_B10_S6_User_Reported_Review_',date,'.json'))
if(!file.exists(approval))stop('HOLD_B10A_S6_SOURCE_REVIEW_REQUIRED__NO_FORMAL_EFFECTS')
review<-fromJSON(approval)
stopifnot(review$evidence_type=='USER_REPORTED_COAUTHOR_SOURCE_REVIEW',nzchar(review$reviewer),nzchar(review$review_date),
 'B10_GSE174367_AD_EXPLORATORY_SIX_LINEAGES' %in% review$accepted_scopes)
contract<-fromJSON(file.path(root,'0_admin/protocol/RP241',paste0('RP241_AD_GBM_B10_Execution_Input_Code_Freeze_',date,'.json')),simplifyVector=FALSE)
for(r in contract$inputs) {
 p<-file.path(root,r$relative_path)
 got<-substr(system2('sha256sum',shQuote(p),stdout=TRUE),1,64)
 if(length(got)!=1||got!=r$sha256)stop('FROZEN_INPUT_HASH_MISMATCH ',r$relative_path)
}
out<-file.path(root,'3_results/rp241_ad_gbm_b10a_ad_lineage_composition')
if(dir.exists(out)&&length(list.files(out)))stop('OUTPUT_EXISTS__PRESERVE_PRIOR_RUN_NO_OVERWRITE')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
base<-file.path(root,'3_results/rp241_ad_gbm_b10_source_model_freeze')
cc<-as.matrix(read.delim(gzfile(file.path(base,paste0('GSE174367_All58676_StableID_DonorLineage_IntegerCounts_',date,'.tsv.gz'))),row.names=1,check.names=FALSE))
md<-read.delim(file.path(base,paste0('GSE174367_All126_DonorLineage_Pseudobulk_Metadata_',date,'.tsv')),check.names=FALSE)
cand<-read.delim(file.path(base,paste0('All81_Immutable_Identity_Direction_Input_',date,'.tsv')),check.names=FALSE)
cand$primary43<-as.character(cand$primary43) %in% c('TRUE','True','true','1')
md$eligible_50_cells<-as.character(md$eligible_50_cells) %in% c('TRUE','True','true','1')
stopifnot(nrow(cand)==81,sum(cand$primary43)==43,ncol(cc)==126,nrow(cc)==58676)
elig<-read.delim(file.path(base,paste0('GSE174367_OutcomeBlind_Model_Eligibility_',date,'.tsv')))
elig$eligible_for_declared_primary_or_sensitivity<-as.character(elig$eligible_for_declared_primary_or_sensitivity) %in% c('TRUE','True','true','1')
lineages<-c('ASC','EX','INH','MG','ODC','OPC');allrows<-list();progrows<-list();warns<-list();dq<-list();influence<-list();score_influence<-list()
forms<-list(M0=~Batch+Age+Sex+Diagnosis,S1_RIN=~Batch+Age+Sex+RIN+Diagnosis,S2_PMI_CC=~Batch+Age+Sex+PMI+Diagnosis)
for(model in names(forms))for(lin in lineages) {
 b<-md[md$author_lineage==lin & md$eligible_50_cells,,drop=FALSE]
 b<-b[complete.cases(b[,all.vars(forms[[model]]),drop=FALSE]),,drop=FALSE]
 b$Batch<-factor(b$Batch);b$Sex<-factor(b$Sex);b$Diagnosis<-factor(b$Diagnosis,levels=c('Control','AD'));rownames(b)<-b$pseudobulk_id
 b$Age_raw<-b$Age;b$Age<-(b$Age-mean(b$Age))/10
 for(k in c('RIN','PMI')) {b[[paste0(k,'_raw')]]<-b[[k]];b[[k]]<-b[[k]]-mean(b[[k]],na.rm=TRUE)}
 write_tsv(b,file.path(out,paste0(model,'_',lin,'_Raw_And_Centered_Model_Metadata.tsv')))
 er<-elig[elig$model==model & elig$author_lineage==lin,,drop=FALSE]
 rr<-cand;rr$lineage<-lin;rr$model<-model;rr$n_AD<-sum(b$Diagnosis=='AD');rr$n_Control<-sum(b$Diagnosis=='Control')
 rr$log2FoldChange<-rr$lfcSE<-rr$stat<-rr$pvalue<-rr$max_raw_Cook<-NA_real_;rr$beta_converged<-FALSE
 rr$test_status<-'NOT_RUN_SMALL_N_SENSITIVITY';rr$normalization<-'DESeq2_poscounts'
 if(isTRUE(er$eligible_for_declared_primary_or_sensitivity)) {
  fit<-fit_counts(cc[,rownames(b),drop=FALSE],b,forms[[model]]);rs<-fit$result
  saveRDS(fit$fit,file.path(out,paste0(model,'_',lin,'_DESeq2_fit.rds')))
  write_tsv(rs,file.path(out,paste0(model,'_',lin,'_AllSourceGenes_Diagnostics.tsv')))
  ix<-match(rr$AD_external_Ensembl_stable_ID,rs$Ensembl_stable_ID)
  for(k in c('log2FoldChange','lfcSE','stat','pvalue','max_raw_Cook','beta_converged'))rr[[k]]<-rs[[k]][ix]
  rr$test_status<-ifelse(is.na(ix),'LOW_COUNT_OR_SOURCE_NOT_MEASURABLE',ifelse(!rr$beta_converged|!is.finite(rr$pvalue),'MODEL_P_NOT_FINITE_OR_NOT_CONVERGED','TESTABLE'))
  rr$design_rank<-fit$rank;rr$residual_df<-fit$df;rr$dispersion_fit_type<-fit$fit_type
  warns[[paste(model,lin,sep='_')]]<-fit$warning
  if(model=='M0') {
   # Fast conditional influence: fixed main gene filter, size factors and
   # dispersions. Does NOT purport to re-estimate the full-source dispersion.
   rid<-which(rownames(fit$fit) %in% cand$AD_external_Ensembl_stable_ID)
   for(drop in rownames(b)) {
    bk<-b[rownames(b)!=drop,,drop=FALSE];X<-model.matrix(forms$M0,bk)
    ir<-cand;ir$lineage<-lin;ir$removed_donor<-b[drop,'donor_id'];ir$n_AD<-sum(bk$Diagnosis=='AD');ir$n_Control<-sum(bk$Diagnosis=='Control')
    ir$method<-'CONDITIONAL_FIXED_MAIN_FILTER_SIZEFACTORS_DISPERSIONS';ir$beta<-ir$SE<-ir$P2<-NA_real_
    ir$fit_status<-'RANK_NOT_IDENTIFIABLE';ir$below5_after_drop<-ir$n_AD<5|ir$n_Control<5
    if(qr(X)$rank==ncol(X)&&nrow(X)-ncol(X)>=3) {
     di<-fit$fit[rid,rownames(bk)];design(di)<-forms$M0
     di<-nbinomWaldTest(di,betaPrior=FALSE,quiet=TRUE)
     ri<-as.data.frame(results(di,contrast=c('Diagnosis','AD','Control'),independentFiltering=FALSE,cooksCutoff=FALSE))
     ix2<-match(ir$AD_external_Ensembl_stable_ID,rownames(ri));ir$beta<-ri$log2FoldChange[ix2];ir$SE<-ri$lfcSE[ix2];ir$P2<-ri$pvalue[ix2]
     ir$fit_status<-ifelse(is.na(ix2),'MAIN_GENE_NOT_MEASURABLE',ifelse(is.finite(ir$beta)&is.finite(ir$P2),'CONDITIONAL_INFLUENCE_FIT','NONFINITE'))
    }
    ir$direction_correct<-ifelse(grepl('^R1',ir$direction),1,-1)*ir$beta>0;influence[[length(influence)+1]]<-ir
   }
  }
  # Equal-weight new43 score. Never use DE P, coefficients, or diagnosis-specific weights.
  norm<-log2(counts(fit$fit,normalized=TRUE)+1);primary<-cand[cand$primary43,,drop=FALSE]
  ii<-match(primary$AD_external_Ensembl_stable_ID,rownames(norm));z<-matrix(NA_real_,nrow(primary),ncol(norm))
  for(i in seq_len(nrow(primary)))if(!is.na(ii[i])) {
   v<-norm[ii[i],];s<-sd(v);if(is.finite(s)&&s>0)z[i,]<-(v-mean(v))/s
  }
  r1<-grepl('^R1',primary$direction);nc1<-sum(is.finite(z[r1,1]));nc2<-sum(is.finite(z[!r1,1]))
  prog<-list(R1=if(nc1>=23)colMeans(z[r1,,drop=FALSE],na.rm=TRUE) else rep(NA_real_,ncol(z)),
             R2=if(nc2>=12)colMeans(z[!r1,,drop=FALSE],na.rm=TRUE) else rep(NA_real_,ncol(z)))
  prog$AMS43<-prog$R1-prog$R2
  for(component in names(prog)) {
   vals<-prog[[component]];pr<-data.frame(model=model,lineage=lin,component=component,covered_R1=nc1,denominator_R1=28,covered_R2=nc2,denominator_R2=15,status='INSUFFICIENT_COMPONENT_COVERAGE',beta=NA_real_,SE=NA_real_,CI_lower=NA_real_,CI_upper=NA_real_,P_two_sided=NA_real_)
   if(all(is.finite(vals))) {a<-linear_hc3(vals,b,forms[[model]]);pr$status<-'COHORT_WITHIN_LINEAGE_STANDARDIZED_EXPLORATORY';for(k in names(a)[names(a)!='df'])pr[[k]]<-a[k]}
   if(model=='M0'&&all(is.finite(vals)))for(i in seq_len(nrow(b))) {
    bk<-b[-i,,drop=FALSE];X<-model.matrix(forms$M0,bk)
    if(qr(X)$rank==ncol(X)) {aa<-linear_hc3(vals[-i],bk,forms$M0);score_influence[[length(score_influence)+1]]<-data.frame(lineage=lin,component=component,removed_donor=b$donor_id[i],t(aa),normalization='FIXED_FULL_COHORT_SCORE_SCALING_CONDITIONAL_INFLUENCE')}
   }
   progrows[[length(progrows)+1]]<-pr
   write_tsv(data.frame(donor=b$donor_id,model=model,lineage=lin,component=component,score=vals),file.path(out,paste0(model,'_',lin,'_',component,'_DonorScores.tsv')))
  }
  members<-data.frame(primary,score_usable=is.finite(z[,1]));write_tsv(members,file.path(out,paste0(model,'_',lin,'_All43_Program_Members_Coverage.tsv')))
 }else {
  rr$design_rank<-er$rank;rr$residual_df<-er$residual_df;rr$dispersion_fit_type<-NA_character_
  for(component in c('R1','R2','AMS43'))progrows[[length(progrows)+1]]<-data.frame(model=model,lineage=lin,component=component,covered_R1=NA_integer_,denominator_R1=28,covered_R2=NA_integer_,denominator_R2=15,status='NOT_RUN_SMALL_N_SENSITIVITY',beta=NA_real_,SE=NA_real_,CI_lower=NA_real_,CI_upper=NA_real_,P_two_sided=NA_real_)
 }
 rr$CI95_lower<-rr$log2FoldChange-1.96*rr$lfcSE;rr$CI95_upper<-rr$log2FoldChange+1.96*rr$lfcSE
 signs<-ifelse(grepl('^R1',rr$direction),1,-1);rr$P_one_sided<-one_sided(rr$stat,signs);rr$direction_correct<-signs*rr$log2FoldChange>0
 rr$P_one_sided[rr$test_status!='TESTABLE']<-NA_real_;allrows[[length(allrows)+1]]<-rr
 cat(model,lin,'complete\n');gc()
}
res<-do.call(rbind,allrows);res$BH_primary258<-res$BY_primary258<-NA_real_;res$BH_secondary486<-NA_real_
for(model in names(forms)) {
 ix<-which(res$model==model & res$primary43);stopifnot(length(ix)==258)
 res$BH_primary258[ix]<-fixed_fdr(res$P_one_sided[ix],258);res$BY_primary258[ix]<-fixed_fdr(res$P_one_sided[ix],258,'BY')
 ix<-which(res$model==model);stopifnot(length(ix)==486);res$BH_secondary486[ix]<-fixed_fdr(res$P_one_sided[ix],486)
}
write_tsv(res,file.path(out,'All81_AllSixLineages_AllThreeModels_Complete_Ledger.tsv'))
write_tsv(res[res$primary43,,drop=FALSE],file.path(out,'All43_AllSixLineages_AllThreeModels_Complete_Ledger.tsv'))
pr<-do.call(rbind,progrows);pr$BH_fixed18<-NA_real_;for(model in names(forms)){ix<-which(pr$model==model);stopifnot(length(ix)==18);pr$BH_fixed18[ix]<-fixed_fdr(pr$P_two_sided[ix],18)}
write_tsv(pr,file.path(out,'R1_R2_AMS43_Complete_Program_Family.tsv'))
write_tsv(do.call(rbind,influence),file.path(out,'All81_M0_Conditional_LeaveOneDonorOut_Complete_Ledger.tsv'))
write_tsv(do.call(rbind,score_influence),file.path(out,'Program_M0_Conditional_LeaveOneDonorOut.tsv'))
composition<-list();donors<-unique(md$donor_id)
for(lin in c(lineages,'PER.END')) {
 b<-md[md$author_lineage==lin,,drop=FALSE];stopifnot(nrow(b)==18)
 b$Batch<-factor(b$Batch);b$Sex<-factor(b$Sex);b$Diagnosis<-factor(b$Diagnosis,levels=c('Control','AD'))
 b$Age<-(b$Age-mean(b$Age))/10
 total<-tapply(md$cells,md$donor_id,sum)[b$donor_id];y<-log((b$cells+.5)/(total-b$cells+.5))
 a<-linear_hc3(y,b,forms$M0);composition[[lin]]<-data.frame(lineage=lin,t(a),inference='DONOR_EQUAL_WEIGHT_RELATIVE_CAPTURE_LOGIT_NOT_ABSOLUTE_TISSUE')
}
co<-do.call(rbind,composition);co$BH_fixed7<-fixed_fdr(co$P_two_sided,7);write_tsv(co,file.path(out,'SevenLineages_Donor_RelativeCapture_Composition.tsv'))
write_json(list(status='B10A_EXECUTED__EXPLORATORY_REANALYSIS_NOT_NEW_INDEPENDENT_CONFIRMATION',new43_source='GSE174367',gene_family_primary=258,gene_family_secondary=486,program_family=18,composition_family=7,
   LODO_conditional_executed=TRUE,LODO_full_dispersion_refit_executed=FALSE,GBM_tumor_state_executed=FALSE,population_GBM_confirmation=FALSE,cell_intrinsic_antagonism_established=FALSE,Tier1=0,fit_warnings=warns,
   next_required='INDEPENDENT_NUMERICAL_VALIDATION_BEFORE_BIOLOGICAL_CLOSURE'),file.path(out,'Execution_Receipt.json'),pretty=TRUE,auto_unbox=TRUE)
