"""Integrate complete context ledgers, do not manufacture causal support."""
import argparse,json
from pathlib import Path
import pandas as pd,numpy as np
DATE='2026-10-14'
def run(root):
 for name in ['Independent_Raw_Replay_Validation','Independent_Descriptive_Validation']:
  assert json.loads((root/f'0_admin/validation/RP241_AD_GBM_B10B_{name}_{DATE}.json').read_text())['status']=='PASS'
 out=root/'3_results/rp241_ad_gbm_b10b_evidence_integration'
 c=pd.read_csv(out/f'All81_Complete_CrossSource_Context_Ledger_{DATE}.tsv',sep='\t')
 prior=pd.read_csv(root/f'3_results/rp241_ad_gbm_b10a_evidence_integration/All81_Complete_CellContext_Evidence_Ledger_{DATE}.tsv',sep='\t')
 assert len(c)==81 and len(prior)==81
 merged=prior.merge(c,on=['HGNC_ID','HGNC_symbol','direction'],how='left',validate='one_to_one')
 assert len(merged)==81 and np.array_equal(merged.B9A_primary43,merged.primary43)
 merged['C_GBM_evidence_semantics']='TUMOR_ONLY_DONOR_STATE_DESCRIPTIVE_CONTEXT__NOT_DISEASE_REPLICATION'
 merged['C_identity_uncertainty']='AUTHOR_LABELS_PRIMARY__COPYKAT_NOT_GROUND_TRUTH__SENSITIVITY_PENDING'
 merged['composition_adjusted_bulk_inverse_expression_established']=False
 merged['AD_programme_substate_composition_adjustment_executed']=False
 merged['G_current43_projection_executed']=False
 merged['P_protein_or_spatial_executed']=False;merged['F_real_function_executed']=False;merged['Tier1_causal_target']=False
 merged.to_csv(out/f'All81_AD_GBM_Complete_Multidimension_Context_Ledger_{DATE}.tsv',sep='\t',index=False,na_rep='NA')
 merged[merged.primary43].to_csv(out/f'All43_AD_GBM_Complete_Multidimension_Context_Ledger_NoPrune_{DATE}.tsv',sep='\t',index=False,na_rep='NA')
 prof=pd.read_csv(out/f'All81_ThreeCompartment_CommonDonor_Profile_{DATE}.tsv',sep='\t')
 margins=[]
 for r in prof.itertuples():
  vals=[getattr(r,s+'_mean_logCPM') for s in ['malignant','myeloid','oligodendrocyte']]
  if all(np.isfinite(vals)):
   v=sorted(vals,reverse=True);gap=v[0]-v[1]
  else:gap=np.nan
  margins.append({'source':r.source,'HGNC_ID':r.HGNC_ID,'HGNC_symbol':r.HGNC_symbol,'primary43':r.primary43,
   'dominant_compartment':r.dominant_compartment,'point_top_minus_second_logCPM':gap,'common_N':r.common_three_donors50,
   'margin_is_significance_or_gene_selection_rule':False})
 pd.DataFrame(margins).to_csv(out/f'All81_Point_Profile_Margins_Descriptive_NoThreshold_{DATE}.tsv',sep='\t',index=False,na_rep='NA')
 # Keep all CNV categories, including literal NA, in the annotation audit.
 qc=pd.read_csv(root/f'3_results/rp241_ad_gbm_b10b_state_localization/GSE182109_Author_Assignment_CopyKat_Audit_{DATE}.tsv',sep='\t',keep_default_na=False)
 q=qc.groupby(['role','Assignment','CopyKatPrediction'],as_index=False).cells.sum()
 q['within_author_assignment_fraction']=q.cells/q.groupby(['role','Assignment']).cells.transform('sum')
 q['prediction_is_histologic_truth']=False
 q.to_csv(out/f'GSE182109_AllRoles_AuthorAssignment_CopyKat_Fractions_{DATE}.tsv',sep='\t',index=False)
 prim=c[c.primary43];s=prim[prim.cross_source_status=='SAME_DOMINANT_DESCRIPTIVE_CONTEXT']
 strat=s.groupby(['direction','GSE211376_dominant']).size().reset_index(name='genes')
 strat.to_csv(out/f'Primary43_Point_Context_Concordance_Summary_{DATE}.tsv',sep='\t',index=False)
 dec={'date':DATE,'status':'COMPLETE_B10B_DONOR_STATE_LOCALIZATION__GO_B10C_COMPOSITION_AND_CELL_IDENTITY_SENSITIVITY',
  'S6_reused':'USER_REPORTED_QI_TENG_2026-10-14_LIMITED_GBM_SCOPES__NO_RECONFIRMATION_NEEDED',
  'candidate_primary43_kept':43,'full81_kept':81,'main_GBM_sources':{'GSE211376':11,'GSE182109':14},
  'profile_common_three_donors':{'GSE211376':10,'GSE182109':7},
  'primary43_same_point_dominant_among_three':39,'primary43_different_point_dominant':2,'primary43_identity_not_comparable':2,
  'R1_same_oligodendrocyte_point_context':23,'R2_same_malignant_point_context':13,
  'discordant_primary_genes':['APEX1','RNF130'],'missing_source211_primary_genes':['FLCN','OR7A5'],
  'full81_cross_status':c.cross_source_status.value_counts().to_dict(),
  'all_new_state_metrics_descriptive':True,'new_P_values':0,'new_GBM_AMS43_scores':0,
  'normal_GBM_controls_present':False,'modern_IDHwt_population_confirmation':False,'C_cell_intrinsic_AD_GBM_antagonism_established':False,
  'G_current43_projection_executed':False,'P_protein_or_spatial_executed':False,'F_real_function_executed':False,'Tier1_causal_targets':0,
  'bulk_composition_proven_to_explain_inverse_expression':False,'old310_and_genetic_nulls_unchanged':True,
  'next_goal':'B10C_FIXED_AD_NEURONAL_PROGRAMME_COMPOSITION_CONDITIONAL_MODEL_PLUS_GBM_ANNOTATION_CNV_SENSITIVITY__KEEP43_AND81'}
 (root/f'0_admin/decision/RP241_AD_GBM_B10B_Formal_Decision_{DATE}.json').write_text(json.dumps(dec,ensure_ascii=False,indent=2),encoding='utf-8')
 print(json.dumps(dec,ensure_ascii=False,indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--root',default=r'H:\SCI2\MRQI');a=ap.parse_args();run(Path(a.root))
