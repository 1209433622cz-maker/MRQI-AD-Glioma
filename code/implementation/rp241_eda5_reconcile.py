"""EDA5 public-source evidence reconciliation. No disease refit or new FDR."""
from pathlib import Path
import hashlib,json,platform,sys
import numpy as np
import pandas as pd

R=Path('H:/SCI2/MRQI');D='2026-10-16'
O=R/f'3_results/rp241_eda5_evidence_reconciliation_{D}'
Q=R/f'0_admin/validation/rp241_eda5_{D}'
I={
 'identity':'3_results/rp241_ad_gbm_b10_source_model_freeze/All81_Immutable_Identity_Direction_Input_2026-10-14.tsv',
 'AD159':'3_results/rp241_ad_gbm_b8_evidence_integration/All508_Candidate_Complete_AD_GBM_Evidence_Ledger_2026-10-14.tsv',
 'discovery_sensitivity':'3_results/rp241_ad_gbm_b7_reciprocal_discovery/All16068_Primary_Effects_Predeclared_Source_Sensitivity_2026-10-14.tsv',
 'AD125_model_config':'0_admin/protocol/RP241/RP241_AD_B2E_Frozen_Execution_Config_2026-10-13.json',
 'GBM132':'3_results/rp241_ad_gbm_b9a_limited_validation/M0_ALL5_HOST_MEAN_All81_Frozen_Direction_Results.tsv',
 'GBM132_context':'3_results/rp241_ad_gbm_b9a_evidence_integration/All81_Complete_AD_GBM_Limited_Compartment_Evidence_Ledger_2026-10-14.tsv',
 'GBM132_probe_context':'3_results/rp241_ad_gbm_b9a_evidence_integration/All81_Probe_Suffix_Specificity_Annotations_No_Prune_2026-10-14.tsv',
 'BA17':'3_results/rp241_eda4_gse203206_2026-10-16/All81_All7Models_Expression_Validation_2026-10-16.tsv',
 'BA17_LODO':'3_results/rp241_eda4_gse203206_2026-10-16/All81_M0_LeaveOneDonor_Descriptive_Closure_2026-10-16.tsv',
 'ADcell':'3_results/rp241_ad_gbm_b10a_ad_lineage_composition/All81_AllSixLineages_AllThreeModels_Complete_Ledger.tsv',
 'ADprogram':'3_results/rp241_ad_gbm_b10a_ad_lineage_composition/R1_R2_AMS43_Complete_Program_Family.tsv',
 'GBMstate':'3_results/rp241_ad_gbm_b10b_evidence_integration/All81_Complete_CrossSource_Context_Ledger_2026-10-14.tsv',
 'GBMstate_donor':'3_results/rp241_eda1_public_evidence_2026-10-16/GBM_All81_CrossSource_Mean_And_Donor_Context_2026-10-16.tsv',
 'GBM_CNV':'3_results/rp241_ad_gbm_b10c_cnv_description/All81_TwoMask_vs_AuthorPrimary_vs_Source211_Context_2026-10-14.tsv',
 'GBM_regional':'3_results/rp241_eda3_regional_new_AD_QC_2026-10-16/All81_Source_ThreeMethod_Fixed_And_Matched_Patient_Means_2026-10-16.tsv',
 'GBM_region_influence':'3_results/rp241_eda3_regional_new_AD_QC_2026-10-16/All17_All81_Region_Influence_Summary_2026-10-16.tsv',
 'GBMprotein':'3_results/rp241_ad_gbm_b10d_p2_protein_direction/GBM_PRIMARY_All81_NoPrune.tsv',
 'GBMprotein_unshared':'3_results/rp241_ad_gbm_b10d_p2_protein_direction/GBM_UNSHARED_All81_NoPrune.tsv',
 'spatial':'3_results/rp241_ad_gbm_b10d_p4_patient_spatial_description/All81_All3Pair_EqualPatient_Descriptions_NoP_2026-10-15.tsv',
 'spatial_MP':'3_results/rp241_ad_gbm_b10d_p3_spatial_source/All81_Author_MP_Gene_Overlap_Circularity_Audit_2026-10-15.tsv',
 'RP241':'9.Version/AD_GBM_研究计划书_RP2.4.1_完整版_2026-10-09.docx',
 'scope':f'0_admin/protocol/RP241/RP241_EDA5_CrossSource_Reconciliation_Scope_{D}.md',
 'controlled_boundary':'0_admin/protocol/RP241/RP241_AD_GBM_P5A_Controlled_Data_AI_Boundary_Addendum_2026-10-15.md',
 'protein_application_status':f'0_admin/audit/RP241_DeepConsensus_Application_Status_User_Confirmation_{D}.json',
 'BA17_S6':f'0_admin/audit/RP241_EDA4_GSE203206_S6_Human_Review_{D}.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(df,name):
 O.mkdir(parents=True,exist_ok=True)
 (df if isinstance(df,pd.DataFrame) else pd.DataFrame(df)).to_csv(O/f'{name}_{D}.tsv',sep='\t',index=False,na_rep='NA',float_format='%.12g')
def write(p,v):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def flag(v):
 if pd.isna(v):return None
 if isinstance(v,(bool,np.bool_)):return bool(v)
 if str(v).strip().lower() in ['true','1','1.0']:return True
 if str(v).strip().lower() in ['false','0','0.0']:return False
 raise ValueError(f'Invalid boolean {v!r}')
def load(k):return pd.read_csv(R/I[k],sep='\t')
def state(beta,q,expected,p=None):
 if not np.isfinite(beta):return 'MISSING_OR_NOT_ESTIMABLE'
 direction=np.sign(beta)==expected
 if np.isfinite(q) and q<.05:return 'FDR_DIRECTION_SUPPORT' if direction else 'FDR_OPPOSITE_EFFECT'
 if beta==0:return 'ZERO_POINT_NO_FDR_SUPPORT'
 return 'DIRECTION_ONLY_NO_FDR' if direction else 'POINT_OPPOSITE_NO_FDR'

def main():
 O.mkdir(parents=True,exist_ok=True);Q.mkdir(parents=True,exist_ok=True)
 assert sha(R/I['identity'])=='93dff5698e60a1b03d2a232ab2963ffe5cd8fd60699bda479d11ed7f1eec24e0'
 assert sha(R/I['RP241'])=='f752314510725261c71e6ac72557b9c94b5f81d5428c8a644d403d660fe5203b'
 ad125_config=json.loads((R/I['AD125_model_config']).read_text(encoding='utf-8'))
 assert ad125_config['formula']=='~ age_centered_per_10y + sex_M + diagnosis'
 receipt=[dict(key=k,path=p,bytes=(R/p).stat().st_size,sha256=sha(R/p)) for k,p in I.items()]
 save(receipt,'Input_Source_Provenance_Receipt')
 b=load('identity');b.primary43=b.primary43.map(flag);assert len(b)==81 and b.primary43.sum()==43 and b.HGNC_ID.is_unique
 ids=list(b.HGNC_ID);expected=dict(zip(b.HGNC_ID,np.where(b.direction.str.startswith('R1'),1,-1)))
 a=load('AD159');g=load('GBM132');gc=load('GBM132_context');probe=load('GBM132_probe_context').set_index('HGNC_ID');ds=load('discovery_sensitivity').rename(columns={'HGNC_ID_AD':'HGNC_ID','HGNC_symbol_frozen':'HGNC_symbol'})
 ba=load('BA17');bl=load('BA17_LODO');cell=load('ADcell');prog=load('ADprogram');ctx=load('GBMstate');donor=load('GBMstate_donor')
 cnv=load('GBM_CNV');reg=load('GBM_regional');ri=load('GBM_region_influence');prot=load('GBMprotein');pu=load('GBMprotein_unshared');sp=load('spatial');mp=load('spatial_MP')
 # Enforce one-to-one identity without fuzzy names, substring or outcome selection.
 for name,df,multi in [('AD159',a,False),('GBM132',g,False),('GBM132_context',gc,False),('discovery',ds,False),('BA17',ba,True),('BA17_LODO',bl,False),('ADcell',cell,True),('GBMstate',ctx,False),('GBM_CNV',cnv,True),('protein',prot,False),('protein_unshared',pu,False),('spatial',sp,True),('spatial_MP',mp,False)]:
  assert set(ids)<=set(df.HGNC_ID),(name,'missing frozen gene')
  if not multi:assert df.HGNC_ID.is_unique,(name,'duplicate key')
  if 'HGNC_symbol' in df:
   check=df[df.HGNC_ID.isin(ids)].merge(b[['HGNC_ID','HGNC_symbol']],on='HGNC_ID',suffixes=('_source','_frozen'))
   assert (check.HGNC_symbol_source==check.HGNC_symbol_frozen).all(),name
 a=a.set_index('HGNC_ID');g=g.set_index('HGNC_ID');gc=gc.set_index('HGNC_ID');ds=ds.set_index('HGNC_ID');ctx=ctx.set_index('HGNC_ID');donor=donor.set_index('HGNC_ID');prot=prot.set_index('HGNC_ID');pu=pu.set_index('HGNC_ID');mp=mp.set_index('HGNC_ID');bl=bl.set_index('HGNC_ID')
 ba=ba.set_index(['HGNC_ID','model']);cell=cell.set_index(['HGNC_ID','model','lineage']);sp=sp.set_index(['HGNC_ID','zone_a','zone_b'])
 assert set(a[a.BH_full508.lt(.05)].index)==set(ids)
 assert set(g[g.BH_full81.lt(.05)].index)==set(b[b.primary43].HGNC_ID)
 assert (a.loc[ids,'Ensembl_stable_ID'].to_numpy()==b.AD_external_Ensembl_stable_ID.to_numpy()).all()
 masters=[];long=[];sens=[];joinlog=[]
 def ev(hid,source,model,role,beta,se,p,q,qfamily,denom,expected_sign,units,clinical,lo=np.nan,hi=np.nan,p_sides='two-sided',n=np.nan,source_status='source_finite_result_or_original_NA_retained'):
  row=dict(HGNC_ID=hid,gene=b.set_index('HGNC_ID').loc[hid,'HGNC_symbol'],primary43=bool(b.set_index('HGNC_ID').loc[hid,'primary43']),source=source,model=model,role=role,beta=beta,SE=se,CI95_low=lo,CI95_high=hi,P=p,P_semantics=p_sides,q=q,q_family=qfamily,q_denominator=denom,expected_sign=expected_sign,status=state(beta,q,expected_sign),original_source_status=source_status,units=units,clinical_scope=clinical,n=n)
  long.append(row);return row['status']
 for r in b.itertuples():
  hid=r.HGNC_ID;adsg=expected[hid];gbsg=-adsg;v=a.loc[hid];z=g.loc[hid];d=ds.loc[hid];zz=gc.loc[hid]
  row=dict(HGNC_ID=hid,gene=r.HGNC_symbol,Ensembl_stable_ID=r.AD_external_Ensembl_stable_ID,direction=r.direction,primary43=bool(r.primary43),selected_by_AD159=True,selected_by_GBM132=bool(r.primary43))
  assert np.sign(v.discovery_AD_beta)==adsg and np.sign(v.discovery_TUMOR_HOST_beta)==gbsg
  row.update(AD125_beta=v.discovery_AD_beta,AD125_SE=v.discovery_AD_SE,AD125_P2=v.discovery_AD_P2,AD125_CI_low=d.AD_CI95_lower_log2,AD125_CI_high=d.AD_CI95_upper_log2,GBM165_beta=v.discovery_TUMOR_HOST_beta,GBM165_SE=v.discovery_TUMOR_HOST_SE,GBM165_P2=v.discovery_TUMOR_HOST_P2,IUT_q32136=v.discovery_IUT_BH_full2G,IUT_BY32136=v.discovery_IUT_BY_full2G,discovery_role='SELECTED_DISCOVERY_NOT_EXTERNAL_VALIDATION',GBM165_grade4_beta=d.GRADE4_beta,GBM165_LOPO_sign_stable_n=d.LOPO_sign_matches_primary_count,AD125_no95587_beta=d.S3_REMOVE_EXPLICIT_GSE95587_OVERLAPPING_DONORS_beta)
  ev(hid,'GSE125583','M0','selected_discovery',v.discovery_AD_beta,v.discovery_AD_SE,v.discovery_AD_P2,np.nan,'not_a_source_q__see_IUT_32136',32136,adsg,'DESeq2_log2FC','FuG;218AD/70Control',d.AD_CI95_lower_log2,d.AD_CI95_upper_log2,n=288)
  ev(hid,'GSE165595','M0','selected_discovery',v.discovery_TUMOR_HOST_beta,v.discovery_TUMOR_HOST_SE,v.discovery_TUMOR_HOST_P2,np.nan,'not_a_source_q__see_IUT_32136',32136,gbsg,'DESeq2_log2FC','15 historical tumour-host pairs;11 grade4+4 grade3',d.TUMOR_HOST_CI95_lower_log2,d.TUMOR_HOST_CI95_upper_log2,n=15)
  row.update(AD159_beta=v.log2FoldChange_AD_minus_Old,AD159_SE=v.lfcSE,AD159_P1=v.one_sided_P,AD159_q508=v.BH_full508,AD159_BY508=v.BY_full508,AD159_LODO_same_direction_n=v.LODO_correct_direction_count,AD159_LODO_q508_support_n=v.LODO_BH_supported_count,AD159_common_age_q508=v.common_age_BH_full508,AD159_raw_Cook_warning=flag(v.raw_Cook_above_default_cutoff))
  row['AD159_state']=ev(hid,'GSE159699','M0','external_AD_then_81_selection',v.log2FoldChange_AD_minus_Old,v.lfcSE,v.one_sided_P,v.BH_full508,'one_sided_BH508',508,adsg,'DESeq2_log2FC','BA20/21;12AD/10Old',v.CI95_lower_log2,v.CI95_upper_log2,p_sides='one-sided frozen direction',n=22)
  row.update(GBM132_beta=z.beta_TUMOR_minus_HOST,GBM132_P1=z.P_one_sided,GBM132_q81=z.BH_full81,GBM132_BY81=z.BY_full81,GBM132_unique_probe_n=z.unique_probe_n,GBM132_original_status=z.status,GBM132_identity_status=z.pre_effect_identity_status,GBM132_LOPO_same_direction_n=zz.LOPO_correct_direction_n if np.isfinite(z.beta_TUMOR_minus_HOST) else np.nan,GBM132_LOPO_q81_support_n=zz.LOPO_BH_supported_n if np.isfinite(z.beta_TUMOR_minus_HOST) else np.nan)
  row.update(GBM132_exact_probe_ids=probe.loc[hid,'exact_GPL96_probes'],GBM132_x_at_flag=flag(probe.loc[hid,'contains_x_at_probe']),GBM132_s_at_flag=flag(probe.loc[hid,'contains_s_at_probe']),GBM132_gene_ID_is_not_sequence_specificity_proof=True)
  row['GBM132_state']=ev(hid,'GSE13276','M0','limited_core_adjacent_WM_then_43_selection',z.beta_TUMOR_minus_HOST,z.moderated_SE,z.P_one_sided,z.BH_full81,'one_sided_BH81',81,gbsg,'RMA_log2_intensity_difference','5 adult historical GBM;core-adjacent WM',z.CI95_lower,z.CI95_upper,p_sides='one-sided frozen direction',n=5,source_status=z.status)
  for model in ['M0_primary','M1_mapping','M2_PMI','M3_age_common_support','M4_EOAD','M5_LOAD','M6_relative_scale']:
   t=ba.loc[(hid,model)];qval=t.q_BH_43 if r.primary43 else t.q_BH_81;qn=43 if r.primary43 else 81
   ss=ev(hid,'GSE203206',model,'bounded_BA17_primary' if model=='M0_primary' else 'supportive_only_no_primary_rescue',t.beta,t.SE_HC3,t.P_HC3,qval,f'two_sided_HC3_BH{qn}',qn,adsg,'log2_provided_plus1' if model!='M6_relative_scale' else 'log2_relative_provided_plus1',f'BA17;{int(t.AD)}AD/{int(t.Control)}Control',t.CI95_low,t.CI95_high,n=int(t.n))
   sens.append(dict(HGNC_ID=hid,gene=r.HGNC_symbol,primary43=bool(r.primary43),model=model,beta=t.beta,P_HC3=t.P_HC3,q43=t.q_BH_43,q81=t.q_BH_81,status=ss,interpretation_role='primary' if model=='M0_primary' else 'sensitivity_no_rescue'))
  t=ba.loc[(hid,'M0_primary')];m6=ba.loc[(hid,'M6_relative_scale')];l=bl.loc[hid]
  qcol='q_BH_43' if r.primary43 else 'q_BH_81';qt=float(t[qcol]);qm=float(m6[qcol])
  row.update(BA17_beta=t.beta,BA17_HC3_SE=t.SE_HC3,BA17_CI_low=t.CI95_low,BA17_CI_high=t.CI95_high,BA17_P2=t.P_HC3,BA17_q43=t.q_BH_43,BA17_q81=t.q_BH_81,BA17_state=state(t.beta,qt,adsg),BA17_q_scope='primary43' if r.primary43 else 'background81',BA17_M6_beta=m6.beta,BA17_M6_q43=m6.q_BH_43,BA17_M6_q81=m6.q_BH_81,BA17_scale_FDR_status_differs=bool(np.isfinite(qt) and np.isfinite(qm) and ((qt<.05)!=(qm<.05))),BA17_scale_point_sign_differs=bool(np.isfinite(t.beta) and np.isfinite(m6.beta) and np.sign(t.beta)!=np.sign(m6.beta)),BA17_LODO_sign_flip_n=l.LODO_sign_flips_from_M0,BA17_LODO_beta_min=l.LODO_beta_min,BA17_LODO_beta_max=l.LODO_beta_max,BA17_missing_reason='EXACT_ENSEMBL_ID_ABSENT' if not flag(t.measured) else 'MEASURED')
  for model in ['M0','S1_RIN','S2_PMI_CC']:
   hits=[];by=[];meas=0;opposite=[];unmeasured=[]
   for lineage in ['ASC','EX','INH','MG','ODC','OPC']:
    c=cell.loc[(hid,model,lineage)];cv=c.BH_primary258 if r.primary43 else c.BH_secondary486;cd=258 if r.primary43 else 486
    st=ev(hid,'GSE174367',model+'_'+lineage,'AD_lineage_reanalysis',c.log2FoldChange,c.lfcSE,c.P_one_sided,cv,f'one_sided_BH{cd}',cd,adsg,'DESeq2_log2FC',f'{lineage};{int(c.n_AD)}AD/{int(c.n_Control)}Control',c.CI95_lower,c.CI95_upper,p_sides='one-sided frozen direction',n=int(c.n_AD+c.n_Control),source_status=c.test_status)
    if np.isfinite(c.log2FoldChange):meas+=1
    if st=='FDR_DIRECTION_SUPPORT':hits.append(lineage)
    if np.isfinite(c.BY_primary258) and c.BY_primary258<.05:by.append(lineage)
    if st=='POINT_OPPOSITE_NO_FDR':opposite.append(lineage)
    if not np.isfinite(c.log2FoldChange):unmeasured.append(lineage+':'+c.test_status)
   row[f'ADcell_{model}_testable_lineages']=meas;row[f'ADcell_{model}_BH_support_lineages']=';'.join(hits)
   row[f'ADcell_{model}_BY_primary_support_lineages']=';'.join(by);row[f'ADcell_{model}_point_opposite_lineages']=';'.join(opposite)
   row[f'ADcell_{model}_untestable_reasons']=';'.join(unmeasured)
  c=ctx.loc[hid];dd=donor.loc[hid]
  row.update(GBM211_peak=c.GSE211376_dominant,GBM182_peak=c.GSE182109_dominant,GBM211_common_patients=c.GSE211376_common_N,GBM182_common_patients=c.GSE182109_common_N,GBM_state_cross_source=c.cross_source_status,GBM_state_is_disease_replication=False,GBM211_donor_agreement_n=dd.source211_strict_donor_agree_n,GBM182_donor_agreement_n=dd.source182_strict_donor_agree_n,GBM_both_all_donors_agree=flag(dd.same_mean_and_all_donors_both_sources),GBM_both_all_LODO_same_peak=flag(dd.same_mean_and_all_LODO_both_sources))
  cv=cnv[cnv.HGNC_ID==hid];row['GBM_CNV_mask_changes_point_label_n']=sum(flag(v) is True for v in cv.author_primary_point_label_changed)
  rr=ri[ri.HGNC_ID==hid];row['GBM_regional_patient_peak_changes_n']=int(rr.donor_peak_changed_n.sum());row['GBM_regional_group_peak_changes_n']=int(rr.group_peak_changed_n.sum())
  rt=reg[(reg.HGNC_ID==hid)&reg.scope.eq('FIXED_ORIGINAL_PATIENTS')&reg.source.eq('GSE182109')].set_index('method')
  assert set(rt.index)=={'A_all_regions_pooled_count','B_common_regions_pooled_count','C_common_regions_equal_logCPM'}
  pa=rt.loc['A_all_regions_pooled_count','peak'];pb=rt.loc['B_common_regions_pooled_count','peak'];pc=rt.loc['C_common_regions_equal_logCPM','peak']
  row['GBM182_region_A_peak']=pa;row['GBM182_region_B_peak']=pb;row['GBM182_region_C_peak']=pc
  row['GBM182_AB_label_change']=False if pd.isna(pa) or pd.isna(pb) else pa!=pb
  row['GBM182_BC_label_change']=False if pd.isna(pb) or pd.isna(pc) else pb!=pc
  for label,p in [('main',prot.loc[hid]),('unshared',pu.loc[hid])]:
   qv=p.OLS_Q_primary43 if r.primary43 else p.OLS_Q_all81_exploratory
   actualn=int(prot.loc[b[b.primary43].HGNC_ID,'OLS_P'].notna().sum()) if r.primary43 else int(prot.OLS_P.notna().sum())
   row[f'GBMprotein_{label}_beta']=p.OLS_beta;row[f'GBMprotein_{label}_OLS_q']=qv;row[f'GBMprotein_{label}_HC3_primary_q']=p.HC3_Q_primary43;row[f'GBMprotein_{label}_status']=state(p.OLS_beta,qv,gbsg)
   row[f'GBMprotein_{label}_original_status']=p.status
   for fld in ['n','cases','controls','rank','df','max_leverage']:row[f'GBMprotein_{label}_{fld}']=p[fld]
   row[f'GBMprotein_{label}_OLS_q_actual_n']=int(pu.loc[b[b.primary43].HGNC_ID,'OLS_P'].notna().sum()) if (label=='unshared' and r.primary43) else (int(pu.OLS_P.notna().sum()) if label=='unshared' else actualn)
   ev(hid,'CPTAC_PDC000204',label,'public_GBM_protein_limited_source_confounded',p.OLS_beta,p.OLS_SE,p.OLS_P,qv,'two_sided_OLS_BH_actual_estimable_primary_subset' if r.primary43 else 'two_sided_OLS_BH_actual_estimable_background',row[f'GBMprotein_{label}_OLS_q_actual_n'],gbsg,'PDC_log_ratio','11 sequence-reported-negative historical GBM /8 shared-plex GTEx; gene complete-case',p.OLS_CI_low,p.OLS_CI_high,n=p.n,source_status=p.status)
  p=prot.loc[hid];row['GBMprotein_HC3_state']=state(p.OLS_beta,p.HC3_Q_primary43,gbsg) if r.primary43 else 'NO_PRIMARY43_Q_FOR_BACKGROUND'
  overlap=mp.loc[hid];row['spatial_author_MP_overlap_n']=int(overlap.membership_n);row['spatial_author_label_defining_MP_memberships']=overlap.author_label_defining_MP_exact_memberships
  for xa,xb,key in [('st_mal','st_norm','stmal_stnorm'),('dis_mal','dis_norm','dismal_disnorm'),('st_mal','dis_mal','stmal_dismal')]:
   s=sp.loc[(hid,xa,xb)];row[f'spatial_{key}_N']=int(s.eligible_patients);row[f'spatial_{key}_mean']=s.equal_patient_mean;row[f'spatial_{key}_min']=s.range_min;row[f'spatial_{key}_max']=s.range_max;row[f'spatial_{key}_positive_patients']=s.positive_point_n;row[f'spatial_{key}_negative_patients']=s.negative_point_n
  row.update(spatial_inferential_replication=False,new_81_target_genetics='NOT_EXECUTED__LEGACY_PROGRAM_NOT_TRANSFERRED',real_functional_validation='NOT_EXECUTED',causal_Tier1=False,population_normal_GBM_validation_established=False,same_cell_AD_GBM_antagonism_established=False)
  masters.append(row)
 master=pd.DataFrame(masters);save(master,'All81_Complete_CrossSource_Evidence_Ledger');save(master[master.primary43],'Original43_Complete_CrossSource_Evidence_Ledger')
 save(pd.DataFrame(long),'All81_Source_Model_Quantitative_Evidence_Long');save(pd.DataFrame(sens),'All81_BA17_SevenModel_Status_NoReselection');save(prog,'AD174367_Original43_Program_M0_and_Sensitivity_Unchanged')
 compact=['HGNC_ID','gene','direction','primary43','AD159_q508','GBM132_state','GBM132_q81','BA17_state','BA17_q43','BA17_q81','BA17_scale_FDR_status_differs','BA17_LODO_sign_flip_n','ADcell_M0_BH_support_lineages','ADcell_S1_RIN_BH_support_lineages','GBM_state_cross_source','GBM_regional_patient_peak_changes_n','GBMprotein_main_status','GBMprotein_HC3_state','spatial_author_MP_overlap_n']
 save(master[compact],'All81_Readable_Status_Review');save(master.loc[master.primary43,compact],'Original43_Readable_Status_Review')
 counts=[]
 for subset,mm in [('all81',master),('original43',master[master.primary43])]:
  for col in ['GBM132_state','BA17_state','GBM_state_cross_source','GBMprotein_main_status','GBMprotein_HC3_state']:
   for stateval,n in mm[col].value_counts(dropna=False).items():counts.append(dict(subset=subset,dimension=col,status=stateval,n=int(n),fixed_subset_n=len(mm),is_independent_replication_rate=False))
 save(counts,'Status_Counts_With_Selection_And_Denominator_Notes')
 summary=dict(date=D,status='EDA5_COMPLETE__SELECTION_HISTORY_AND_HETEROGENEOUS_SOURCE_EVIDENCE_RECONCILED__NO_UNIVERSAL_OR_CAUSAL_CLAIM',identities=81,original43=43,original43_R1=int(master[master.primary43].direction.str.startswith('R1').sum()),AD159_support_selected81=int(master.AD159_state.eq('FDR_DIRECTION_SUPPORT').sum()),GBM132_support_selected43=int(master.GBM132_state.eq('FDR_DIRECTION_SUPPORT').sum()),BA17_primary43_point_concordant=int(master.loc[master.primary43,'BA17_state'].eq('DIRECTION_ONLY_NO_FDR').sum()),BA17_primary43_point_opposite=int(master.loc[master.primary43,'BA17_state'].eq('POINT_OPPOSITE_NO_FDR').sum()),BA17_primary43_missing=int(master.loc[master.primary43,'BA17_state'].eq('MISSING_OR_NOT_ESTIMABLE').sum()),BA17_primary43_FDR_support=0,BA17_all81_point_opposite=int(master.BA17_state.eq('POINT_OPPOSITE_NO_FDR').sum()),BA17_scale_FDR_changed_genes=master.loc[master.BA17_scale_FDR_status_differs,'gene'].tolist(),ADcell_M0_primary43_genes_with_BH258_context=int(master.loc[master.primary43,'ADcell_M0_BH_support_lineages'].ne('').sum()),ADcell_S1_RIN_primary43_genes_with_BH258_context=int(master.loc[master.primary43,'ADcell_S1_RIN_BH_support_lineages'].ne('').sum()),GBM_primary43_same_state_point_peaks=int(master.loc[master.primary43,'GBM_state_cross_source'].eq('SAME_DOMINANT_DESCRIPTIVE_CONTEXT').sum()),GBM_primary43_protein_OLS_FDR_same=int(master.loc[master.primary43,'GBMprotein_main_status'].eq('FDR_DIRECTION_SUPPORT').sum()),GBM_primary43_protein_OLS_FDR_opposite_genes=master.loc[master.primary43&master.GBMprotein_main_status.eq('FDR_OPPOSITE_EFFECT'),'gene'].tolist(),GBM_primary43_protein_HC3_FDR_same=int(master.loc[master.primary43,'GBMprotein_HC3_state'].eq('FDR_DIRECTION_SUPPORT').sum()),GBM_primary43_region_peak_sensitive_genes=master.loc[master.primary43&master.GBM_regional_patient_peak_changes_n.gt(0),'gene'].tolist(),new_disease_models_fitted=0,new_multiple_testing_adjustments=0,controlled_files_or_restricted_derivatives_read=False,causal_Tier1=0,legacy_genetics_relabelled_as_current81_null=False,source_inputs=receipt)
 write(O/f'EDA5_Reconciliation_Summary_{D}.json',summary)
 assert all(sha(R/rec['path'])==rec['sha256'] for rec in receipt)
 write(Q/f'Builder_Integrity_Receipt_{D}.json',dict(status='PASS',all81_and_original43_preserved=True,source_gene_symbols_and_keys_exact=True,AD159_selected_set_equals_original81=True,GBM132_selected_set_equals_original43=True,all_input_hashes_unchanged=True,code_sha256=sha(Path(__file__)),environment=dict(python=sys.version,executable=sys.executable,pandas=pd.__version__,numpy=np.__version__,platform=platform.platform())))
 print(json.dumps({k:v for k,v in summary.items() if k!='source_inputs'},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
