"""Complete fixed CNV sensitivity descriptions, no disease P or new scoring."""
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd
from rp241_ad_gbm_b10c_common import ROOT,DATE,COMP,MASKS,check_freeze
PAIRS=[('malignant','myeloid'),('malignant','oligodendrocyte'),('myeloid','oligodendrocyte')]
B=2000;SEED=20261014
def CI(v,idx):return tuple(np.quantile(v[idx].mean(axis=1),[.025,.975])) if len(v)>=2 else (np.nan,np.nan)
def run(root):
 check_freeze(root)
 extra=json.loads((root/f'0_admin/protocol/RP241/RP241_AD_GBM_B10C_Descriptive_Code_Registration_{DATE}.json').read_text())
 from rp241_ad_gbm_b10c_common import sha
 for x in extra['inputs']:assert sha(root/x['relative_path'])==x['sha256']
 out=root/'3_results/rp241_ad_gbm_b10c_cnv_description'
 if out.exists() and any(out.iterdir()):raise RuntimeError('Preserve existing outputs')
 out.mkdir(parents=True)
 c=pd.read_csv(root/f'3_results/rp241_ad_gbm_b10_source_model_freeze/All81_Immutable_Identity_Direction_Input_{DATE}.tsv',sep='\t');rng=np.random.default_rng(SEED)
 summaries=[];pairs=[];lodo=[];profiles=[];common=[]
 for mask in MASKS:
  d=pd.read_csv(root/f'3_results/rp241_ad_gbm_b10c_cnv_count/{mask}/GSE182109_All81_DonorState_Complete_Localization_{DATE}.tsv',sep='\t')
  assert len(d)==8748 and not d.duplicated(['donor_id','level','state','HGNC_ID']).any()
  for (role,level,state),z in d.groupby(['role','level','state'],sort=True):
   ds=sorted(z.loc[z.eligible_50,'donor_id'].unique());n=len(ds);idx=rng.integers(0,n,(B,n)) if n>=2 else None
   for x in c.itertuples():
    v=z[(z.HGNC_ID==x.HGNC_ID)&z.donor_id.isin(ds)].sort_values('donor_id');identity=bool(z[z.HGNC_ID==x.HGNC_ID].identity_present.all())
    yy=v.log2_CPM_plus1.values if n and identity else np.array([]);lo,hi=CI(yy,idx)
    summaries.append({'mask':mask,'role':role,'level':level,'state':state,'HGNC_ID':x.HGNC_ID,'HGNC_symbol':x.HGNC_symbol,'primary43':bool(x.primary43),
     'donors50':n,'cells_eligible':int(v.cells.sum()),'identity_present':identity,'mean_logCPM':float(yy.mean()) if len(yy) else np.nan,
     'median_CPM':float(v.CPM.median()) if len(yy) else np.nan,'median_detection':float(v.detection_fraction.median()) if len(yy) else np.nan,
     'CI95_lower':lo,'CI95_upper':hi,'status':'IDENTITY_NA' if not identity else 'NO_DONORS' if n==0 else 'SMALL_N_DESCRIPTION' if n<5 else 'DESCRIPTIVE_ONLY',
     'P_value_computed':False,'disease_normal_contrast':False,'CI_family_adjusted':False})
  main=d[(d.role=='MAIN_GBM')&(d.level=='broad_compartment')]
  for label,sub in [('MAIN_ALL',main),('MAIN_NEW',main[main.clinical_context=='newly_diagnosed_GBM']),('MAIN_RECURRENT',main[main.clinical_context=='recurrent_GBM'])]:
   for a,b in PAIRS:
    ds=sorted(set(sub.loc[(sub.state==a)&sub.eligible_50,'donor_id'])&set(sub.loc[(sub.state==b)&sub.eligible_50,'donor_id']));n=len(ds);idx=rng.integers(0,n,(B,n)) if n>=2 else None
    for x in c.itertuples():
     aa=sub[(sub.HGNC_ID==x.HGNC_ID)&(sub.state==a)].set_index('donor_id').reindex(ds)
     bb=sub[(sub.HGNC_ID==x.HGNC_ID)&(sub.state==b)].set_index('donor_id').reindex(ds)
     identity=bool(sub[sub.HGNC_ID==x.HGNC_ID].identity_present.all());v=(aa.log2_CPM_plus1-bb.log2_CPM_plus1).values if n and identity else np.array([])
     lo,hi=CI(v,idx);lv=np.array([(v.sum()-u)/(n-1) for u in v]) if n>=2 and len(v) else np.array([])
     pairs.append({'mask':mask,'subgroup':label,'pair':a+'_minus_'+b,'HGNC_ID':x.HGNC_ID,'HGNC_symbol':x.HGNC_symbol,'primary43':bool(x.primary43),'paired_N':n,
      'identity_present':identity,'mean_delta':float(v.mean()) if len(v) else np.nan,'median_delta':float(np.median(v)) if len(v) else np.nan,
      'CI95_lower':lo,'CI95_upper':hi,'LODO_min':float(lv.min()) if len(lv) else np.nan,'LODO_max':float(lv.max()) if len(lv) else np.nan,
      'status':'IDENTITY_NA' if not identity else 'NO_PAIRED_DONORS' if n==0 else 'SMALL_N_DESCRIPTION' if n<5 else 'DESCRIPTIVE_ONLY','CI_family_adjusted':False,'P_computed':False})
     for i,donor in enumerate(ds):lodo.append({'mask':mask,'subgroup':label,'pair':a+'_minus_'+b,'HGNC_ID':x.HGNC_ID,'removed_donor':donor,'remaining_N':n-1,'mean_delta_after_drop':float(lv[i]) if len(lv) else np.nan,'value_removed':float(v[i]) if len(v) else np.nan})
   ds=set(sub.donor_id.unique())
   for s in COMP:ds &= set(sub.loc[(sub.state==s)&sub.eligible_50,'donor_id'])
   ds=sorted(ds);common.append({'mask':mask,'subgroup':label,'common_three50_N':len(ds),'donors':';'.join(ds)})
   for x in c.itertuples():
    identity=bool(sub[sub.HGNC_ID==x.HGNC_ID].identity_present.all());row={'mask':mask,'subgroup':label,'HGNC_ID':x.HGNC_ID,'HGNC_symbol':x.HGNC_symbol,'primary43':bool(x.primary43),'common_N':len(ds),'identity_present':identity}
    if identity and ds:
     p=sub[(sub.HGNC_ID==x.HGNC_ID)&sub.donor_id.isin(ds)].pivot(index='donor_id',columns='state',values='log2_CPM_plus1').loc[ds,COMP];assert np.isfinite(p.values).all()
     m=p.mean(axis=0);top=m.max();labeltop='UNDETECTED' if top==0 else ';'.join(m.index[np.abs(m-top)<=1e-12])
     row['dominant']=labeltop;row['status']='SMALL_N_PROFILE' if len(ds)<5 else 'DESCRIPTIVE_PROFILE'
     for s in COMP:row[s+'_mean_logCPM']=float(m[s]);row[s+'_donor_top_fraction']=float((abs(p[s]-p.max(axis=1))<=1e-12).mean())
    else:
     row['dominant']='NOT_TESTABLE';row['status']='IDENTITY_OR_COMMON_N_NA'
     for s in COMP:row[s+'_mean_logCPM']=np.nan;row[s+'_donor_top_fraction']=np.nan
    profiles.append(row)
 tables=[(summaries,'All81_AllRoles_TwoMask_State_Descriptions'),(pairs,'All81_AllSubgroups_TwoMask_Paired_Descriptions'),(lodo,'All81_AllSubgroups_TwoMask_LODO'),(profiles,'All81_AllSubgroups_TwoMask_CommonProfiles'),(common,'TwoMask_CommonDonor_Identity')]
 for rows,name in tables:pd.DataFrame(rows).to_csv(out/f'{name}_{DATE}.tsv',sep='\t',index=False,na_rep='NA')
 old=pd.read_csv(root/f'3_results/rp241_ad_gbm_b10b_evidence_integration/All81_ThreeCompartment_CommonDonor_Profile_{DATE}.tsv',sep='\t')
 old211=old[old.source=='GSE211376'].set_index('HGNC_ID');old182=old[old.source=='GSE182109'].set_index('HGNC_ID');pr=pd.DataFrame(profiles);cross=[]
 for mask in MASKS:
  pp=pr[(pr['mask']==mask)&(pr.subgroup=='MAIN_ALL')].set_index('HGNC_ID')
  for x in c.itertuples():
   n=pp.loc[x.HGNC_ID];a=old211.loc[x.HGNC_ID];b=old182.loc[x.HGNC_ID]
   comparable=bool(n.identity_present and a.identity_present and n.common_N>=5 and a.common_three_donors50>=5 and n.dominant not in ['UNDETECTED','NOT_TESTABLE'] and ';' not in n.dominant)
   cross.append({'mask':mask,'HGNC_ID':x.HGNC_ID,'HGNC_symbol':x.HGNC_symbol,'direction':x.direction,'primary43':bool(x.primary43),
    'unchanged_source211_dominant':a.dominant_compartment,'source182_author_primary_dominant':b.dominant_compartment,
    'source182_mask_dominant':n.dominant,'source182_mask_common_N':int(n.common_N),
    'cross_status':'SAME_POINT_CONTEXT' if comparable and n.dominant==a.dominant_compartment else 'DIFFERENT_POINT_CONTEXT' if comparable else 'NOT_COMPARABLE',
    'author_primary_point_label_changed':bool(n.dominant!=b.dominant_compartment),'author_primary_not_replaced':True,'C_causal_antagonism':False})
 cf=pd.DataFrame(cross);cf.to_csv(out/f'All81_TwoMask_vs_AuthorPrimary_vs_Source211_Context_{DATE}.tsv',sep='\t',index=False)
 cf[cf.primary43].to_csv(out/f'All43_TwoMask_Complete_Context_NoPrune_{DATE}.tsv',sep='\t',index=False)
 rec={'date':DATE,'status':'TWO_MASK_DESCRIPTIONS_EXECUTED__INDEPENDENT_VERIFICATION_PENDING','state_rows':len(summaries),'pair_rows':len(pairs),'LODO_rows':len(lodo),'profile_rows':len(profiles),
  'P_values_or_GBM_programmes':0,'complete43_and81_retained':True,'bootstrap_B':B,'bootstrap_seed':SEED,'primary43_by_mask':{m:cf[(cf['mask']==m)&cf.primary43].cross_status.value_counts().to_dict() for m in MASKS}}
 (out/f'Execution_Receipt_{DATE}.json').write_text(json.dumps(rec,indent=2),encoding='utf-8');print(json.dumps(rec,indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--root',default=str(ROOT));run(Path(ap.parse_args().root))
