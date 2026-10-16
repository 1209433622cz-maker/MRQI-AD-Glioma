"""Fixed-family, source-qualified reciprocal expression discovery; no external AD effects."""
from pathlib import Path
from datetime import datetime
import hashlib,importlib.util,json,sys
import numpy as np,pandas as pd
root=Path(sys.argv[1] if len(sys.argv)>1 else r'H:\SCI2\MRQI');date='2026-10-14'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()
freeze=json.loads((root/f'0_admin/protocol/RP241/RP241_AD_GBM_B7_Directional_Execution_Freeze_{date}.json').read_text())
assert all(sha(root/p)==v for p,v in freeze['sha256'].items())
gbmdir=root/'3_results/rp241_gbm_b1_discovery'
receipt=json.loads((gbmdir/'Execution_Receipt.json').read_text())
assert receipt['status']=='FULL17_GBM_SOURCE_MODELS_COMPLETE__EXPLORATORY_ONLY' and receipt['family']==21442 and receipt['models']==17
summary=pd.read_csv(gbmdir/'All17_Models_Execution_Summary.tsv',sep='\t')
assert len(summary)==17 and summary['family'].eq(21442).all() and summary.actual_fit.eq('parametric').all()
spec=importlib.util.spec_from_file_location('fixed_iut',root/'2_code/rp241/rp241_reciprocal_direction_preflight.py')
iut=importlib.util.module_from_spec(spec);spec.loader.exec_module(iut)
dest=root/'3_results/rp241_ad_gbm_b7_reciprocal_discovery';assert not dest.exists();dest.mkdir()
inputs=[]
def read(p,**kwargs):
    inputs.append({'relative_path':p.relative_to(root).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)})
    return pd.read_csv(p,sep='\t',**kwargs)
g=read(root/f'3_results/rp241_ad_b6_fullqc_joint_universe/JointG_Frozen_PreEffect_Unique_Genes_{date}.tsv',dtype=str)
family=read(root/f'3_results/rp241_ad_b6_fullqc_joint_universe/All2G_Directional_Hypotheses_PreEffect_Freeze_{date}.tsv',dtype=str)
ad=read(root/'3_results/rp241_ad_b2e_exploratory/M0_CLINICAL_STUDY_LABEL_PRIMARY_All19108_Results_2026-10-13.tsv',dtype={'original_Entrez_ID':str})
gbm=read(gbmdir/'M0_PAIRED15_FULL30_AllPrimaryGene_Results.tsv',dtype={'gene_id':str})
assert len(g)==16068 and len(family)==32136 and len(ad)==19108 and len(gbm)==21442
assert ad.original_Entrez_ID.is_unique and gbm.gene_id.is_unique
a=ad.set_index('original_Entrez_ID').loc[g.original_Entrez_ID].reset_index()
b=gbm.set_index('gene_id').loc[g.GENCODE27_gene_ID].reset_index()
assert a.HGNC_ID.tolist()==g.HGNC_ID_AD.tolist() and a.GENCODE27_ID.tolist()==g.GENCODE27_gene_ID.tolist()
wide=g[['original_Entrez_ID','GENCODE27_gene_ID','HGNC_ID_AD','HGNC_symbol_frozen']].copy()
for key in ['log2FoldChange_AD_minus_Control','lfcSE','Wald_statistic','pvalue','CI95_lower_log2','CI95_upper_log2','p_value_status']:
    wide['AD_'+key]=a[key].to_numpy()
for key in ['log2FoldChange_Tumor_minus_Normal','lfcSE','Wald_statistic','pvalue','CI95_lower_log2','CI95_upper_log2','p_value_status']:
    wide['TUMOR_HOST_'+key]=b[key].to_numpy()
rows=[]
for j,r in g.iterrows():
    ab=float(a.iloc[j].log2FoldChange_AD_minus_Control);ap=float(a.iloc[j].pvalue)
    tb=float(b.iloc[j].log2FoldChange_Tumor_minus_Normal);tp=float(b.iloc[j].pvalue)
    p1,p2=iut.reciprocal(ab,ap,tb,tp)
    for direction,p in [('R1_AD_UP_HOST_TUMOR_DOWN',p1),('R2_AD_DOWN_HOST_TUMOR_UP',p2)]:
        rows.append({'HGNC_ID':r.HGNC_ID_AD,'HGNC_symbol':r.HGNC_symbol_frozen,'AD_original_Entrez_ID':r.original_Entrez_ID,
         'GBM_original_gene_id':r.GENCODE27_gene_ID,'direction':direction,'AD_beta':ab,'AD_SE':a.iloc[j].lfcSE,'AD_P2':ap,
         'TUMOR_HOST_beta':tb,'TUMOR_HOST_SE':b.iloc[j].lfcSE,'TUMOR_HOST_P2':tp,
         'IUT_P':p,'testability':'TESTED' if np.isfinite(p) else 'NOT_TESTABLE_SOURCE_P_OR_EFFECT_MISSING',
         'full_directional_family':32136,'evidence_scope':'TWO_SOURCE_DIFFERENT_ESTIMANDS_EXPLORATORY_RECIPROCITY__NOT_CONFIRMED_CAUSAL_ANTAGONISM'})
alltests=pd.DataFrame(rows)
assert alltests[['HGNC_ID','direction']].to_records(index=False).tolist()==family[['HGNC_ID','direction']].to_records(index=False).tolist()
alltests['IUT_BH_full2G']=iut.adjust_full_family(alltests.IUT_P.to_numpy())
alltests['IUT_BY_full2G']=iut.adjust_full_family(alltests.IUT_P.to_numpy(),'BY')
alltests['selected_for_reserved_validation']=alltests.IUT_BH_full2G.lt(.05)
primary_beta=b.log2FoldChange_Tumor_minus_Normal.to_numpy()
lo=[];lop=[]
for r in summary.itertuples():
    if r.model=='M0_PAIRED15_FULL30':continue
    xx=read(gbmdir/(r.model+'_AllPrimaryGene_Results.tsv'),dtype={'gene_id':str})
    assert len(xx)==21442 and xx.gene_id.tolist()==gbm.gene_id.tolist()
    x=xx.set_index('gene_id').loc[g.GENCODE27_gene_ID].reset_index()
    if r.model=='S1_GRADE4_PAIRED11':
        wide['GRADE4_beta']=x.log2FoldChange_Tumor_minus_Normal.to_numpy();wide['GRADE4_P2']=x.pvalue.to_numpy()
    else:
        lo.append(x.log2FoldChange_Tumor_minus_Normal.to_numpy());lop.append(x.pvalue.to_numpy())
lo=np.array(lo);lop=np.array(lop);assert lo.shape==(15,16068)
wide['LOPO_sign_matches_primary_count']=((np.sign(lo)==np.sign(primary_beta))&np.isfinite(lo)).sum(axis=0)
wide['LOPO_finite_beta_count']=np.isfinite(lo).sum(axis=0);wide['LOPO_finite_P_count']=np.isfinite(lop).sum(axis=0)
wide['LOPO_beta_min']=np.nanmin(lo,axis=0);wide['LOPO_beta_max']=np.nanmax(lo,axis=0)
for model in ['S1_LOW_TANGLE_CONTROL_VS_HIGH_TANGLE_AD','S2_COMMON_OBSERVED_AGE_RANGE','S3_REMOVE_EXPLICIT_GSE95587_OVERLAPPING_DONORS']:
    x=read(root/('3_results/rp241_ad_b2e_exploratory/'+model+'_All19108_Results_2026-10-13.tsv'),dtype={'original_Entrez_ID':str})
    assert len(x)==19108;x=x.set_index('original_Entrez_ID').loc[g.original_Entrez_ID].reset_index()
    wide[model+'_beta']=x.log2FoldChange_AD_minus_Control.to_numpy();wide[model+'_P2']=x.pvalue.to_numpy()
alltests.to_csv(dest/f'All32136_Directional_IUT_BH_BY_Exploratory_{date}.tsv',sep='\t',index=False,na_rep='NA')
wide.to_csv(dest/f'All16068_Primary_Effects_Predeclared_Source_Sensitivity_{date}.tsv',sep='\t',index=False,na_rep='NA')
selected=alltests[alltests.selected_for_reserved_validation].copy()
assert not selected.HGNC_ID.duplicated().any()
selected.to_csv(dest/f'All_Primary_BH_Selected_Directions_Reserved_Validation_Candidates_{date}.tsv',sep='\t',index=False,na_rep='NA')
# Coverage is recorded only after selection, without recomputing G or opening validation effects.
coverage=read(root/'3_results/rp241_ad_b5_source_admission/GSE159699_AllGene_ExactSymbol_Identity_Coverage_2026-10-14.tsv',dtype=str).fillna('')
valid=coverage.identity_status.eq('EXACT_CURRENT_SYMBOL_UNIQUE_ENSEMBL')&coverage.ensembl_source_collision.eq('False')
validmap=coverage[valid].set_index('ensembl_exact').source_gene
assert validmap.index.is_unique
lookup=g.set_index('HGNC_ID_AD').HGNC_Ensembl_gene_ID
selected['Ensembl_stable_ID']=selected.HGNC_ID.map(lookup)
selected['GSE159699_exact_source_row']=selected.Ensembl_stable_ID.map(validmap)
selected['reserved_coverage_status']=np.where(selected.GSE159699_exact_source_row.notna(),'REFERENCE_IDENTITY_MAPPED__COUNTS_MEASURABILITY_AND_EFFECT_NOT_EXECUTED','NOT_TESTABLE_EXACT_UNIQUE_IDENTITY_ABSENT_OR_CONFLICT')
selected['reserved_effect_opened']=False;selected['full_selected_direction_validation_family']=len(selected)
selected.to_csv(dest/f'All_Selected_Directions_Reserved_AD_Identity_Coverage_Ledger_{date}.tsv',sep='\t',index=False,na_rep='NA')
pd.DataFrame(inputs).to_csv(root/f'0_admin/manifest/RP241_AD_GBM_B7_Real_Input_Receipt_{date}.tsv',sep='\t',index=False)
result={'executed_at':datetime.now().astimezone().isoformat(),'status':'FIXED2G_RECIPROCAL_EXPLORATORY_DISCOVERY_COMPLETE__INDEPENDENT_VALIDATION_NOT_EXECUTED',
 'joint_G':16068,'directional_family':32136,'finite_directional_tests':int(alltests.IUT_P.notna().sum()),'missing_hypotheses_retained':int(alltests.IUT_P.isna().sum()),
 'BH_R1':int((alltests.direction.eq('R1_AD_UP_HOST_TUMOR_DOWN')&alltests.IUT_BH_full2G.lt(.05)).sum()),
 'BH_R2':int((alltests.direction.eq('R2_AD_DOWN_HOST_TUMOR_UP')&alltests.IUT_BH_full2G.lt(.05)).sum()),
 'BY_R1':int((alltests.direction.eq('R1_AD_UP_HOST_TUMOR_DOWN')&alltests.IUT_BY_full2G.lt(.05)).sum()),
 'BY_R2':int((alltests.direction.eq('R2_AD_DOWN_HOST_TUMOR_UP')&alltests.IUT_BY_full2G.lt(.05)).sum()),
 'selected_validation_directions':len(selected),'reserved_exact_identity_mapped':int(selected.GSE159699_exact_source_row.notna().sum()),
 'reserved_identity_unmapped_or_conflict':int(selected.GSE159699_exact_source_row.isna().sum()),
 'AD_P2_numeric_zero':int(a.pvalue.eq(0).sum()),'TUMOR_HOST_P2_numeric_zero':int(b.pvalue.eq(0).sum()),
 'GBM_models':17,'AD_sensitivity_models':3,'independent_AD_validation_effect_tests':0,
 'healthy_population_comparison_claimed':False,'causal_antagonism_claimed':False,'confirmed_targets':0}
(dest/'Execution_Receipt.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
candidate_path=dest/f'All_Primary_BH_Selected_Directions_Reserved_Validation_Candidates_{date}.tsv'
(root/f'0_admin/protocol/RP241/RP241_AD_GBM_B7_Reserved_Validation_Candidates_Freeze_{date}.json').write_text(json.dumps({
 'frozen_at':datetime.now().astimezone().isoformat(),'discovery_rule':'PRIMARY_IUT_BH_FULL32136_LT0.05',
 'selected_direction_family':len(selected),'independent_source':'GSE159699','independent_effects_opened':False,
 'validation_coverage_cannot_reduce_candidate_family':True,'human_second_author_source_review_claimed':False,
 'sha256':{candidate_path.relative_to(root).as_posix():sha(candidate_path),
 (dest/f'All32136_Directional_IUT_BH_BY_Exploratory_{date}.tsv').relative_to(root).as_posix():sha(dest/f'All32136_Directional_IUT_BH_BY_Exploratory_{date}.tsv')
 }},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2))
