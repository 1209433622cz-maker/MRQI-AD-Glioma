"""Fixed M0 programme scores, all author neuronal ILRs, HC3 + complete LODO."""
from pathlib import Path
import json,argparse
import numpy as np,pandas as pd
from rp241_ad_gbm_b10c_common import ROOT,DATE,check_freeze,ilr,design,hc3,fixed_bh
def run(root):
 check_freeze(root);out=root/'3_results/rp241_ad_gbm_b10c_ad_composition'
 if out.exists() and any(out.iterdir()):raise RuntimeError('Existing outputs preserved')
 out.mkdir(parents=True)
 src=root/'3_results/rp241_ad_gbm_b10a_ad_lineage_composition';metaonly=root/'3_results/rp241_ad_gbm_b10b_evidence_integration'
 old=pd.read_csv(src/'R1_R2_AMS43_Complete_Program_Family.tsv',sep='\t');res=[];loo=[];scores=[];cov=[]
 for lin,n in [('EX',16),('INH',17)]:
  md=pd.read_csv(src/f'M0_{lin}_Raw_And_Centered_Model_Metadata.tsv',sep='\t').sort_values('donor_id').set_index('donor_id');assert len(md)==n
  counts=pd.read_csv(metaonly/f'B10C_MetadataOnly_{lin}_AuthorSubstate_Counts_{DATE}.tsv',sep='\t').set_index('donor_id').loc[md.index]
  assert np.array_equal(counts.sum(axis=1).values,md.cells.values);Z=ilr(counts.values)
  for model,rin,comp in [('M0_reproduced',False,False),('Mcomp_primary',False,True),('M0_RIN_fixedScore',True,False),('Mcomp_RIN_fixedScore',True,True)]:
   take=np.isfinite(md.RIN.to_numpy(float)) if rin else np.ones(len(md),bool);m=md.loc[take];z=Z[take]
   X,names=design(m,z,rin,comp);j=names.index('Diagnosis_AD')
   if model=='Mcomp_primary':assert np.linalg.matrix_rank(X)==(10 if lin=='EX' else 9)
   for i,donor in enumerate(m.index):cov.append({'lineage':lin,'model':model,'donor_id':donor,**dict(zip(names,X[i]))})
   for part in ['R1','R2','AMS43']:
    s=pd.read_csv(src/f'M0_{lin}_{part}_DonorScores.tsv',sep='\t');assert s.donor.is_unique and set(s.donor)==set(md.index)
    y=s.set_index('donor').loc[m.index,'score'].to_numpy(float);assert np.isfinite(y).all()
    fit=hc3(X,y,j);fit.update(lineage=lin,component=part,model=model,AD_N=int((m.Diagnosis=='AD').sum()),Control_N=int((m.Diagnosis=='Control').sum()),
      score_origin='UNCHANGED_B10A_M0_DONOR_SCORES__NO_NEW_SCALING',family='PRIMARY_CONDITIONAL6' if model=='Mcomp_primary' else 'RIN_CONDITIONAL6_SENSITIVITY' if comp else 'BASELINE_REPLAY_NOT_NEW_FAMILY')
    res.append(fit)
    if model=='M0_reproduced':
     ref=old[(old.model=='M0')&(old.lineage==lin)&(old.component==part)].iloc[0]
     for k in ['beta','SE','CI_lower','CI_upper','P_two_sided']:assert abs(fit[k]-ref[k])<1e-10,(lin,part,k)
    for i,donor in enumerate(m.index):scores.append({'lineage':lin,'model':model,'component':part,'donor_id':donor,'fixed_score':y[i]})
    if model in ['M0_reproduced','Mcomp_primary']:
     for i,donor in enumerate(m.index):
      keep=np.arange(len(m))!=i;a=hc3(X[keep],y[keep],j)
      a.update(lineage=lin,component=part,model=model,removed_donor=donor,
       AD_N=int((m.iloc[np.where(keep)[0]].Diagnosis=='AD').sum()),Control_N=int((m.iloc[np.where(keep)[0]].Diagnosis=='Control').sum()),
       score_scaling_fixed=True,ILR_definition_fixed=True)
      a['below5_after_drop']=a['AD_N']<5 or a['Control_N']<5;loo.append(a)
 df=pd.DataFrame(res);df['BH_conditional_fixed6']=np.nan
 for model in ['Mcomp_primary','Mcomp_RIN_fixedScore']:
  ix=df.model==model;assert ix.sum()==6;df.loc[ix,'BH_conditional_fixed6']=fixed_bh(df.loc[ix,'P_two_sided'])
 df.to_csv(out/f'All6_Programmes_All4Models_Complete_Conditional_Ledger_{DATE}.tsv',sep='\t',index=False,na_rep='NA')
 pd.DataFrame(loo).to_csv(out/f'All6_M0_Mcomp_Complete_LeaveOneDonor_{DATE}.tsv',sep='\t',index=False,na_rep='NA')
 pd.DataFrame(cov).to_csv(out/f'AllModels_Fixed_Metadata_ILR_Designs_{DATE}.tsv',sep='\t',index=False,na_rep='NA')
 pd.DataFrame(scores).to_csv(out/f'AllModels_Fixed_Score_Identity_Replay_{DATE}.tsv',sep='\t',index=False)
 pairs=[]
 for lin in ['EX','INH']:
  for part in ['R1','R2','AMS43']:
   a=df[(df.lineage==lin)&(df.component==part)&(df.model=='M0_reproduced')].iloc[0];b=df[(df.lineage==lin)&(df.component==part)&(df.model=='Mcomp_primary')].iloc[0]
   pairs.append({'lineage':lin,'component':part,'M0_beta':a.beta,'M0_SE':a.SE,'Mcomp_beta':b.beta,'Mcomp_SE':b.SE,
    'point_beta_difference_not_formal_difference_test':b.get('beta',np.nan)-a.beta,
    'same_score_same_donor':True,'difference_P_computed':False,'conditional_is_causal_direct_effect':False})
 pd.DataFrame(pairs).to_csv(out/f'All6_SameScore_SameDonor_PointComparison_{DATE}.tsv',sep='\t',index=False,na_rep='NA')
 rec={'date':DATE,'status':'AD_FIXED_PROGRAMME_COMPOSITION_EXECUTED__INDEPENDENT_VERIFICATION_PENDING','main_family':6,
  'total_model_slots':24,'LODO_slots':len(loo),'all_M0_numerical_replay_match':True,
  'new_programme_scaling_or_member_selection':False,'new_single_gene_DE_executed':False,'causal_direct_effect_established':False}
 (out/f'Execution_Receipt_{DATE}.json').write_text(json.dumps(rec,indent=2),encoding='utf-8');print(df[['lineage','component','model','beta','SE','BH_conditional_fixed6']].to_string(index=False))
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--root',default=str(ROOT));run(Path(a.parse_args().root))
