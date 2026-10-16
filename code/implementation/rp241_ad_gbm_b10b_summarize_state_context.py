"""Frozen descriptive donor bootstrap / paired state / cross-source context.
No P values, no programme or gene selection; complete81 and complete43.
"""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
DATE='2026-10-14';SEED=20261014;B=2000
COMP=['malignant','myeloid','oligodendrocyte']
PAIRS=[('malignant','myeloid'),('malignant','oligodendrocyte'),('myeloid','oligodendrocyte')]
def ci_mean(x,idx):
 x=np.asarray(x,float)
 if len(x)<2:return np.nan,np.nan
 b=x[idx].mean(axis=1);return tuple(np.quantile(b,[.025,.975]))
def run(root):
 base=root/'3_results/rp241_ad_gbm_b10b_state_localization'
 out=root/'3_results/rp241_ad_gbm_b10b_evidence_integration';out.mkdir(parents=True,exist_ok=True)
 d=pd.read_csv(base/f'All81_TwoSource_DonorState_Complete_Localization_{DATE}.tsv',sep='\t')
 cand=pd.read_csv(root/f'3_results/rp241_ad_gbm_b10_source_model_freeze/All81_Immutable_Identity_Direction_Input_{DATE}.tsv',sep='\t')
 assert len(cand)==81 and cand.primary43.sum()==43
 rng=np.random.default_rng(SEED);summ=[];pairs=[];dom=[];commonqc=[];lodo=[]
 for (source,role,level,state),xx in d.groupby(['source','role','level','state'],sort=True):
  donors=sorted(xx.loc[xx.eligible_50,'donor_id'].unique());n=len(donors)
  idx=rng.integers(0,n,(B,n)) if n>=2 else None
  for c in cand.itertuples():
   z=xx[(xx.HGNC_ID==c.HGNC_ID)&xx.donor_id.isin(donors)].sort_values('donor_id')
   assert len(z)==n
   identity=bool(xx[xx.HGNC_ID==c.HGNC_ID].identity_present.all());ok=n>0 and identity
   v=z.log2_CPM_plus1.to_numpy(float) if ok else np.array([])
   lo,hi=ci_mean(v,idx) if ok else (np.nan,np.nan)
   summ.append({'source':source,'role':role,'level':level,'state':state,'HGNC_ID':c.HGNC_ID,'HGNC_symbol':c.HGNC_symbol,'primary43':bool(c.primary43),
    'eligible_donors50':n,'eligible_cells':int(z.cells.sum()),'identity_present':identity,
    'mean_log2_CPM_plus1':float(v.mean()) if ok else np.nan,'median_CPM':float(z.CPM.median()) if ok else np.nan,
    'median_detection_fraction':float(z.detection_fraction.median()) if ok else np.nan,'bootstrap_mean_logCPM_CI95_lower':lo,'bootstrap_mean_logCPM_CI95_upper':hi,
    'status':'IDENTITY_NOT_TESTABLE' if not identity else 'NO_ELIGIBLE_DONORS' if not n else 'LOW_N_DESCRIPTIVE_ONLY' if n<5 else 'DONOR_LOCALIZATION_DESCRIPTIVE',
    'CI_family_adjusted':False,'P_value_computed':False,'disease_effect':False})
 # All roles remain visible in summary. Paired descriptions main14/11 only,
 # plus predefined newly/recurrent source182 splits; no effect-driven subsets.
 for source,src in d.groupby('source',sort=True):
  main=src[(src.role=='MAIN_GBM')&(src.level=='broad_compartment')]
  subgroups=[('MAIN_GBM_ALL',main)]
  if source=='GSE182109':
   subgroups += [('MAIN_NEW_GBM',main[main.clinical_context=='newly_diagnosed_GBM']),('MAIN_RECURRENT_GBM',main[main.clinical_context=='recurrent_GBM'])]
  for subgroup,m in subgroups:
   for a,b in PAIRS:
    aa=set(m.loc[(m.state==a)&m.eligible_50,'donor_id']);bb=set(m.loc[(m.state==b)&m.eligible_50,'donor_id']);ds=sorted(aa&bb);n=len(ds)
    idx=rng.integers(0,n,(B,n)) if n>=2 else None
    for c in cand.itertuples():
     ca=m[(m.state==a)&(m.HGNC_ID==c.HGNC_ID)&m.donor_id.isin(ds)].set_index('donor_id').reindex(ds)
     cb=m[(m.state==b)&(m.HGNC_ID==c.HGNC_ID)&m.donor_id.isin(ds)].set_index('donor_id').reindex(ds)
     identity=bool(m[m.HGNC_ID==c.HGNC_ID].identity_present.all());v=(ca.log2_CPM_plus1-cb.log2_CPM_plus1).to_numpy(float) if n and identity else np.array([])
     lo,hi=ci_mean(v,idx) if len(v) else (np.nan,np.nan)
     lm=np.array([(v.sum()-v[i])/(len(v)-1) for i in range(len(v))]) if len(v)>=2 else np.array([])
     pairs.append({'source':source,'subgroup':subgroup,'pair':a+'_minus_'+b,'HGNC_ID':c.HGNC_ID,'HGNC_symbol':c.HGNC_symbol,'primary43':bool(c.primary43),
      'paired_donors50':n,'identity_present':identity,'mean_paired_logCPM_delta':float(v.mean()) if len(v) else np.nan,
      'median_paired_logCPM_delta':float(np.median(v)) if len(v) else np.nan,'bootstrap_CI95_lower':lo,'bootstrap_CI95_upper':hi,
      'LODO_mean_delta_min':float(lm.min()) if len(lm) else np.nan,'LODO_mean_delta_max':float(lm.max()) if len(lm) else np.nan,
      'status':'IDENTITY_NOT_TESTABLE' if not identity else 'NO_PAIRED_DONORS' if n==0 else 'LOW_N_DESCRIPTIVE_ONLY' if n<5 else 'PAIRED_TUMOR_STATE_DESCRIPTION',
      'CI_family_adjusted':False,'P_value_computed':False,'normal_control_comparison':False})
     for i,donor in enumerate(ds):
      lodo.append({'source':source,'subgroup':subgroup,'pair':a+'_minus_'+b,'HGNC_ID':c.HGNC_ID,'HGNC_symbol':c.HGNC_symbol,'removed_donor':donor,
       'remaining_N':n-1,'mean_delta_after_remove':float(lm[i]) if len(lm) else np.nan,'paired_value':float(v[i]) if len(v) else np.nan})
  ds=set(main.donor_id.unique())
  for state in COMP:ds &= set(main.loc[(main.state==state)&main.eligible_50,'donor_id'])
  ds=sorted(ds);commonqc.append({'source':source,'common_three_donors50':len(ds),'donors':';'.join(ds),'paired_compartments':';'.join(COMP),'normal_controls':0})
  for c in cand.itertuples():
   z=main[(main.HGNC_ID==c.HGNC_ID)&main.state.isin(COMP)&main.donor_id.isin(ds)]
   identity=bool(main[main.HGNC_ID==c.HGNC_ID].identity_present.all());row={'source':source,'HGNC_ID':c.HGNC_ID,'HGNC_symbol':c.HGNC_symbol,'primary43':bool(c.primary43),'common_three_donors50':len(ds),'identity_present':identity}
   if identity and ds:
    pv=z.pivot(index='donor_id',columns='state',values='log2_CPM_plus1').reindex(index=ds,columns=COMP);assert np.isfinite(pv.values).all()
    mean=pv.mean(axis=0);mx=mean.max()
    dominant='UNDETECTED_ALL_THREE' if mx==0 else ';'.join(mean.index[np.abs(mean-mx)<=1e-12])
    for state in COMP:
     row[state+'_mean_logCPM']=float(mean[state]);row[state+'_donor_top_fraction']=float((np.abs(pv[state]-pv.max(axis=1))<=1e-12).mean())
    row['dominant_compartment']=dominant;row['status']='DESCRIPTIVE_COMMON_DONOR_PROFILE' if len(ds)>=5 else 'LOW_N_COMMON_DONOR_PROFILE'
   else:
    for state in COMP:row[state+'_mean_logCPM']=np.nan;row[state+'_donor_top_fraction']=np.nan
    row['dominant_compartment']='NOT_TESTABLE';row['status']='IDENTITY_OR_COMMON_DONOR_NOT_TESTABLE'
   dom.append(row)
 df=pd.DataFrame(dom);allrows=[];ranks=[]
 sources=['GSE211376','GSE182109']
 for c in cand.itertuples():
  zs=[df[(df.source==s)&(df.HGNC_ID==c.HGNC_ID)].iloc[0] for s in sources]
  labels=[z.dominant_compartment for z in zs]
  valid=all(z.identity_present and z.common_three_donors50>=5 and z.dominant_compartment not in ['NOT_TESTABLE','UNDETECTED_ALL_THREE'] and ';' not in z.dominant_compartment for z in zs)
  status='SAME_DOMINANT_DESCRIPTIVE_CONTEXT' if valid and labels[0]==labels[1] else 'DIFFERENT_DOMINANT_DESCRIPTIVE_CONTEXT' if valid else 'NOT_COMPARABLE_OR_TIED_OR_UNDETECTED'
  allrows.append({'HGNC_ID':c.HGNC_ID,'HGNC_symbol':c.HGNC_symbol,'direction':c.direction,'primary43':bool(c.primary43),
   'GSE211376_dominant':labels[0],'GSE182109_dominant':labels[1],'GSE211376_common_N':int(zs[0].common_three_donors50),
   'GSE182109_common_N':int(zs[1].common_three_donors50),'cross_source_status':status,'C_disease_replication':False,'C_cell_intrinsic_AD_GBM_antagonism':False,'Tier1':False})
 for universe in ['primary43','full81']:
  ids=set(cand.loc[cand.primary43,'HGNC_ID'] if universe=='primary43' else cand.HGNC_ID)
  for state in COMP:
   x=df[(df.source==sources[0])&df.HGNC_ID.isin(ids)].set_index('HGNC_ID');y=df[(df.source==sources[1])&df.HGNC_ID.isin(ids)].set_index('HGNC_ID')
   common=sorted(set(x.index)&set(y.index));xx=x.loc[common,state+'_mean_logCPM'];yy=y.loc[common,state+'_mean_logCPM'];mask=np.isfinite(xx)&np.isfinite(yy)
   rho=float(spearmanr(xx[mask],yy[mask]).statistic) if mask.sum()>=3 else np.nan
   ranks.append({'universe':universe,'broad_compartment':state,'gene_pairs':int(mask.sum()),'Spearman_rho_descriptive_only':rho,'P_value_computed':False,'genes_are_patients':False})
 tables=[(summ,'All81_AllRoles_State_EqualDonor_Localization_Summary'),(pairs,'All81_Predeclared_PairedState_Descriptions'),(lodo,'All81_PairedState_LeaveOneDonor_Descriptions'),
  (dom,'All81_ThreeCompartment_CommonDonor_Profile'),(commonqc,'ThreeCompartment_CommonDonor_Eligibility'),(allrows,'All81_Complete_CrossSource_Context_Ledger'),(ranks,'FixedUniverse_CrossSource_GeneRank_Description')]
 for rows,name in tables:pd.DataFrame(rows).to_csv(out/f'{name}_{DATE}.tsv',sep='\t',index=False,na_rep='NA')
 ef=pd.DataFrame(allrows);ef[ef.primary43].to_csv(out/f'All43_Complete_CrossSource_Context_Ledger_NoPrune_{DATE}.tsv',sep='\t',index=False)
 assert len(ef)==81 and ef.primary43.sum()==43
 assert len(pd.DataFrame(pairs).query("subgroup=='MAIN_GBM_ALL'"))==486
 rec={'status':'COMPLETE_DESCRIPTIVE_LOCALIZATION__NO_DISEASE_OR_P_INFERENCE','date':DATE,'bootstrap_seed':SEED,'bootstrap_B':B,
  'state_summary_rows':len(summ),'paired_description_rows':len(pairs),'main_two_source_pair81_slots':486,'lodo_description_rows':len(lodo),
  'common_donor_counts':commonqc,'primary43_cross_status':ef[ef.primary43].cross_source_status.value_counts().to_dict(),
  'full81_cross_status':ef.cross_source_status.value_counts().to_dict(),'P_values_generated':0,'new_GBM_program_scores':0,'Tier1':0}
 (out/f'Descriptive_Integration_Receipt_{DATE}.json').write_text(json.dumps(rec,indent=2),encoding='utf-8')
 print(json.dumps(rec,indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--root',default=r'H:\SCI2\MRQI');a=ap.parse_args();run(Path(a.root))
