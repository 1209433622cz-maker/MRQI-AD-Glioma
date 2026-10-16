"""Public count-derived region sensitivity; fixed patients/genes; no tests."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
R=Path('H:/SCI2/MRQI');D='2026-10-16';OLD='2026-10-14';TOL=1e-12
O=R/f'3_results/rp241_eda3_regional_new_AD_QC_{D}';Q=R/f'0_admin/validation/rp241_eda3_{D}'
for p in (O,Q):p.mkdir(parents=True,exist_ok=True)
COMP=['malignant','myeloid','oligodendrocyte'];SOURCES=['GSE211376','GSE182109']
METHODS=['A_all_regions_pooled_count','B_common_regions_pooled_count','C_common_regions_equal_logCPM']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def save(x,name):
 z=x if isinstance(x,pd.DataFrame) else pd.DataFrame(x);z.to_csv(O/f'{name}_{D}.tsv',sep='\t',index=False,na_rep='NA',float_format='%.12g');return z
def yes(x):return str(x).lower() in ['true','1']
def pk(v):
 if not np.isfinite(v).all():return 'NOT_ESTIMABLE'
 if max(v)<=TOL:return 'UNDETECTED'
 ix=np.flatnonzero(abs(v-max(v))<=TOL)
 return COMP[ix[0]] if len(ix)==1 else 'TIED'
def logs(count,lib):return np.log2(1+count/lib[:,None]*1e6)
paths={
 'base':R/f'3_results/rp241_ad_gbm_b10_source_model_freeze/All81_Immutable_Identity_Direction_Input_{OLD}.tsv',
 'donor':R/f'3_results/rp241_ad_gbm_b10b_state_localization/All81_TwoSource_DonorState_Complete_Localization_{OLD}.tsv',
 'common':R/f'3_results/rp241_ad_gbm_b10b_evidence_integration/ThreeCompartment_CommonDonor_Eligibility_{OLD}.tsv',
 'EDA2_profile':R/f'3_results/rp241_eda2_donor_depth_clinical_{D}/All162_Gene_Donor_Depth_Influence_Ledger_{D}.tsv',
 'scope':R/f'0_admin/protocol/RP241/RP241_EDA3_Regional_and_GSE203206_Technical_Exploration_Scope_{D}.md'}
for s in SOURCES:
 for name in ['All81_RegionState_Complete_Localization','Region_Library_Conservation']:paths[s+'_'+name]=R/f'3_results/rp241_ad_gbm_b10b_state_localization/{s}_{name}_{OLD}.tsv'
freeze=dict(date=D,post_EDA2_exploration=True,inputs=[dict(key=k,path=p.relative_to(R).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for k,p in paths.items()],code_sha256=sha(Path(__file__)))
write(Q/f'Regional_Input_Code_Freeze_{D}.json',freeze)
b=pd.read_csv(paths['base'],sep='\t');b.primary43=b.primary43.map(yes);assert len(b)==81 and b.primary43.sum()==43 and b.HGNC_ID.is_unique
gids=b.HGNC_ID.tolist();base=pd.read_csv(paths['donor'],sep='\t')
for c in ['primary43','identity_present','eligible_50']:base[c]=base[c].map(yes)
base=base[(base.role=='MAIN_GBM')&(base.level=='broad_compartment')]
common={r.source:sorted(r.donors.split(';')) for r in pd.read_csv(paths['common'],sep='\t').itertuples()}
coverage=[];structural=[];methodrows=[];region_gene=[];LOR=[];patient_summary=[];group_rows=[];summaries=[];conservation=0
for s in SOURCES:
 sub=pd.read_csv(paths[s+'_All81_RegionState_Complete_Localization'],sep='\t')
 sub=sub[(sub.role=='MAIN_GBM')&(sub.level=='broad_compartment')].copy();sub.identity_present=sub.identity_present.map(yes)
 assert not sub.duplicated(['unit','state','HGNC_ID']).any()
 reg=pd.read_csv(paths[s+'_Region_Library_Conservation'],sep='\t');reg=reg[reg.donor_id.isin(base.loc[base.source==s,'donor_id'])]
 assert reg.unit.is_unique
 tech=sub[['unit','donor_id','state','cells','source_library_UMI']].drop_duplicates()
 assert not tech.duplicated(['unit','state']).any()
 check=tech.groupby('unit')[['cells','source_library_UMI']].sum().reindex(reg.unit)
 assert np.array_equal(check.cells,reg.joined_cells) and np.array_equal(check.source_library_UMI,reg.source_library_UMI);conservation+=len(reg)
 ident=base[base.source==s].groupby('HGNC_ID').identity_present.agg(['min','max']).reindex(gids)
 assert (ident['min']==ident['max']).all();present=ident['min'].to_numpy(bool)
 # Every unit is retained. A compartment absent from author-assigned cells is structural zero, not missing gene identity.
 data={}
 for row in reg.itertuples():
  counts=np.zeros((3,81));counts[:,~present]=np.nan;cells=np.zeros(3,dtype=int);lib=np.zeros(3,dtype=np.int64)
  states=sub[sub.unit==row.unit]
  info=dict(source=s,donor_id=row.donor_id,unit=row.unit,clinical_context=states.clinical_context.iloc[0],fixed_common_patient=row.donor_id in common[s])
  for j,st in enumerate(COMP):
   z=states[states.state==st]
   if len(z):
    z=z.set_index('HGNC_ID').reindex(gids);assert len(z)==81 and z.cells.nunique()==1 and z.source_library_UMI.nunique()==1
    cells[j]=int(z.cells.iloc[0]);lib[j]=int(z.source_library_UMI.iloc[0]);counts[j]=z.raw_UMI.to_numpy(float)
    assert np.array_equal(z.identity_present.to_numpy(bool),present)
   else:structural.append(dict(source=s,donor_id=row.donor_id,unit=row.unit,state=st,cells=0,source_library_UMI=0,source_identity_missing_not_zero_filled=True))
   info[st+'_cells']=int(cells[j]);info[st+'_UMI']=int(lib[j]);info[st+'_eligible50']=bool(cells[j]>=50)
  qual=bool((cells>=50).all());info['common_region_qualified']=qual;info['ineligible_compartments']=';'.join(COMP[j] for j in range(3) if cells[j]<50)
  coverage.append(info);data[row.unit]=dict(donor=row.donor_id,count=counts,cells=cells,lib=lib,info=info,qual=qual)
 methodmap={};qualified_by_donor={};patient_units={}
 for don in common[s]:
  units=sorted(k for k,v in data.items() if v['donor']==don);patient_units[don]=units
  good=[k for k in units if data[k]['qual']];qualified_by_donor[don]=good
  count=sum((data[k]['count'] for k in units));lib=sum((data[k]['lib'] for k in units));cells=sum((data[k]['cells'] for k in units))
  assert (cells>=50).all();orig=logs(count,lib)
  old=base[(base.source==s)&(base.donor_id==don)&base.state.isin(COMP)]
  for j,st in enumerate(COMP):
   z=old[old.state==st].set_index('HGNC_ID').reindex(gids)
   assert np.allclose(orig[j],z.log2_CPM_plus1,equal_nan=True,rtol=0,atol=1e-10)
   assert np.allclose(count[j],z.raw_UMI,equal_nan=True,rtol=0,atol=0) and int(z.source_library_UMI.iloc[0])==int(lib[j])
  if good:
   pooled=logs(sum((data[k]['count'] for k in good)),sum((data[k]['lib'] for k in good)))
   equal=np.mean([logs(data[k]['count'],data[k]['lib']) for k in good],axis=0)
  else:pooled=equal=np.full((3,81),np.nan)
  methodmap[don]=dict(zip(METHODS,[orig,pooled,equal]))
  for meth,arr in methodmap[don].items():
   for gi,g in enumerate(b.itertuples()):
    v=arr[:,gi];row=dict(source=s,donor_id=don,clinical_context=data[units[0]]['info']['clinical_context'],HGNC_ID=g.HGNC_ID,gene=g.HGNC_symbol,primary43=g.primary43,method=meth,
     total_original_regions=len(units),common_qualified_regions=len(good),qualified_units=';'.join(good),identity_present=bool(present[gi]),
     technical_estimable=bool(meth==METHODS[0] or good),peak=pk(v),top_margin=float(np.sort(v)[-1]-np.sort(v)[-2]) if np.isfinite(v).all() else np.nan)
    row['status']='IDENTITY_MISSING' if not present[gi] else 'NO_COMMON_QUALIFIED_REGION' if meth!=METHODS[0] and not good else 'DESCRIPTIVE_ESTIMATE'
    for j,st in enumerate(COMP):row[st+'_logCPM']=v[j]
    methodrows.append(row)
  for unit in units:
   rg=data[unit];vals=logs(rg['count'],rg['lib']) if rg['qual'] else np.full((3,81),np.nan)
   for gi,g in enumerate(b.itertuples()):
    row=dict(source=s,donor_id=don,unit=unit,HGNC_ID=g.HGNC_ID,gene=g.HGNC_symbol,primary43=g.primary43,common_region_qualified=rg['qual'],identity_present=bool(present[gi]),peak=pk(vals[:,gi]))
    for j,st in enumerate(COMP):row[st+'_logCPM']=vals[j,gi]
    region_gene.append(row)
 all_ids=common[s];matched=[d for d in all_ids if qualified_by_donor[d]]
 group_orig=np.mean([methodmap[don][METHODS[0]] for don in all_ids],axis=0)
 for scope,scope_ids in [('FIXED_ORIGINAL_PATIENTS',all_ids),('MATCHED_ESTIMABLE_PATIENTS',matched)]:
  for meth in METHODS:
   available=[don for don in scope_ids if meth==METHODS[0] or qualified_by_donor[don]]
   complete=len(available)==len(scope_ids) and bool(scope_ids)
   arr=np.mean([methodmap[don][meth] for don in scope_ids],axis=0) if complete else np.full((3,81),np.nan)
   for gi,g in enumerate(b.itertuples()):
    row=dict(source=s,scope=scope,method=meth,HGNC_ID=g.HGNC_ID,gene=g.HGNC_symbol,primary43=g.primary43,original_fixed_patients=len(all_ids),
     patients_in_scope=len(scope_ids),technical_available_patients=len(available),complete_scope=complete,patient_ids=';'.join(scope_ids),identity_present=bool(present[gi]),
     peak=pk(arr[:,gi]),status='IDENTITY_MISSING' if not present[gi] else 'FIXED_SCOPE_INCOMPLETE_NO_MEAN' if not complete else 'DESCRIPTIVE_MEAN')
    for j,st in enumerate(COMP):row[st+'_mean_logCPM']=arr[j,gi]
    group_rows.append(row)
 for don in all_ids:
  units=patient_units[don];arr=methodmap[don][METHODS[0]];group_old_labels=[pk(group_orig[:,gi]) for gi in range(81)]
  deletion_by_gene={gi:[] for gi in range(81)}
  for unit in units:
   remaining=[k for k in units if k!=unit]
   cells=sum((data[k]['cells'] for k in remaining),start=np.zeros(3,dtype=int));lib=sum((data[k]['lib'] for k in remaining),start=np.zeros(3,dtype=np.int64))
   eligible=bool(remaining and (cells>=50).all())
   if eligible:
    new=logs(sum((data[k]['count'] for k in remaining)),lib);group_new=(group_orig*len(all_ids)-arr+new)/len(all_ids)
   else:new=group_new=np.full((3,81),np.nan)
   for gi,g in enumerate(b.itertuples()):
    oldpk,npk,gpk=pk(arr[:,gi]),pk(new[:,gi]),pk(group_new[:,gi]);v=new[:,gi]
    rr=dict(source=s,donor_id=don,removed_unit=unit,HGNC_ID=g.HGNC_ID,gene=g.HGNC_symbol,primary43=g.primary43,identity_present=bool(present[gi]),
      original_regions=len(units),remaining_regions=len(remaining),coverage_eligible_after_removal=eligible,
      status='IDENTITY_MISSING' if not present[gi] else 'NO_REGION_REMAINS' if not remaining else 'LOW_COMPARTMENT_COVERAGE' if not eligible else 'DESCRIPTIVE_DELETION',
      original_patient_peak=oldpk,deletion_patient_peak=npk,original_group_peak=group_old_labels[gi],deletion_group_peak=gpk,
      donor_unique_peak_changed=int(oldpk!=npk) if oldpk in COMP and npk in COMP else np.nan,
      group_unique_peak_changed=int(group_old_labels[gi]!=gpk) if group_old_labels[gi] in COMP and gpk in COMP else np.nan,
      max_abs_patient_logCPM_change=float(max(abs(v-arr[:,gi]))) if np.isfinite(v).all() else np.nan)
    for j,st in enumerate(COMP):rr[st+'_remaining_cells']=int(cells[j]);rr[st+'_remaining_UMI']=int(lib[j]);rr[st+'_patient_logCPM']=new[j,gi];rr[st+'_group_logCPM']=group_new[j,gi]
    LOR.append(rr);deletion_by_gene[gi].append(rr)
  for gi,g in enumerate(b.itertuples()):
   z=deletion_by_gene[gi];good=[r for r in z if r['coverage_eligible_after_removal'] and r['identity_present']];dp=[r['donor_unique_peak_changed'] for r in good if np.isfinite(r['donor_unique_peak_changed'])];gp=[r['group_unique_peak_changed'] for r in good if np.isfinite(r['group_unique_peak_changed'])]
   patient_summary.append(dict(source=s,donor_id=don,HGNC_ID=g.HGNC_ID,gene=g.HGNC_symbol,primary43=g.primary43,original_patient_peak=pk(arr[:,gi]),
    regions=len(units),common_qualified_regions=len(qualified_by_donor[don]),attempted_deletions=len(z),coverage_eligible_deletions=sum(r['coverage_eligible_after_removal'] for r in z),identity_present=bool(present[gi]),
    unique_donor_peak_comparisons=len(dp),donor_peak_changed_n=sum(dp) if dp else np.nan,unique_group_peak_comparisons=len(gp),group_peak_changed_n=sum(gp) if gp else np.nan,
    max_abs_logCPM_change=max(r['max_abs_patient_logCPM_change'] for r in good) if good else np.nan,
    status='IDENTITY_MISSING' if not present[gi] else 'NO_ESTIMABLE_DELETION' if not good else 'DESCRIPTIVE_INFLUENCE'))
 summaries.append(dict(source=s,fixed_patients=len(all_ids),matched_estimable_patients=len(matched),missing_common_region_patients=';'.join(d for d in all_ids if d not in matched),
  common_patient_regions=sum(len(patient_units[d]) for d in all_ids),qualified_common_regions=sum(len(qualified_by_donor[d]) for d in all_ids)))
save(coverage,'All_MAIN_Region_ThreeCompartment_Coverage')
save(pd.DataFrame(structural,columns=['source','donor_id','unit','state','cells','source_library_UMI','source_identity_missing_not_zero_filled']),'Absent_Region_Compartment_Structural_Zero_Audit')
save(methodrows,'All17_All81_Patient_ThreeMethod_Descriptions');save(region_gene,'All_CommonRegion_All81_Peak_Descriptions')
group_rows=save(group_rows,'All81_Source_ThreeMethod_Fixed_And_Matched_Patient_Means')
LOR=save(LOR,'All81_LeaveOneRegionOut_Patient_And_Group_Descriptions');save(patient_summary,'All17_All81_Region_Influence_Summary');save(summaries,'Source_Region_Method_Eligibility')
ref=pd.read_csv(paths['EDA2_profile'],sep='\t').set_index(['source','HGNC_ID'])
original=group_rows[(group_rows.scope=='FIXED_ORIGINAL_PATIENTS')&(group_rows.method==METHODS[0])]
for row in original.itertuples():
 expected=ref.loc[(row.source,row.HGNC_ID),'mean_peak'];assert row.peak==expected or (row.peak=='NOT_ESTIMABLE' and expected=='IDENTITY_MISSING')
counts=[]
for (s,scope),z in group_rows.groupby(['source','scope']):
 for universe,genes in [('all81',gids),('primary43',b.loc[b.primary43,'HGNC_ID'].tolist())]:
  zz=z[z.HGNC_ID.isin(genes)];wide=zz.pivot(index='HGNC_ID',columns='method',values='peak')
  for first,last in [(METHODS[0],METHODS[1]),(METHODS[1],METHODS[2]),(METHODS[0],METHODS[2])]:
   valid=wide[first].isin(COMP)&wide[last].isin(COMP)
   counts.append(dict(source=s,scope=scope,universe=universe,first_method=first,last_method=last,panel_genes=len(genes),comparable_unique_peaks=int(valid.sum()),same_peak=int((wide.loc[valid,first]==wide.loc[valid,last]).sum()),different_peak=int((wide.loc[valid,first]!=wide.loc[valid,last]).sum())))
save(counts,'Source_Method_Peak_Comparison_Summary')
influence=[]
for s in SOURCES:
 for universe,zz in [('all81',LOR[LOR.source==s]),('primary43',LOR[(LOR.source==s)&LOR.primary43])]:
  for key in ['donor_unique_peak_changed','group_unique_peak_changed']:
   valid=zz[np.isfinite(zz[key])]
   influence.append(dict(source=s,universe=universe,comparison=key,scenario_gene_records=len(zz),unique_peak_comparisons=len(valid),changed_records=int(valid[key].sum()),
    genes_changed=int(valid.loc[valid[key]==1,'HGNC_ID'].nunique()),patients_contributing_changes=int(valid.loc[valid[key]==1,'donor_id'].nunique()),independent_patient_N_increased=False))
save(influence,'LeaveOneRegionOut_Complete_Influence_Summary')
summary=dict(date=D,status='REGION_SENSITIVITY_COMPLETE__ORIGINAL_PATIENTS_AND_GENES_RETAINED',all81=81,primary43=43,region_units_total=len(coverage),
 region_gene_records=len(region_gene),patient_gene_method_records=len(methodrows),patient_deletion_gene_records=len(LOR),fixed_and_matched_group_records=len(group_rows),
 source_unit_conservation_checks=conservation,structurally_absent_region_compartments=len(structural),eligibility=summaries,
 new_P_FDR=False,new_disease_or_recurrence_effect=False,old_evidence_overwritten=False,controlled_data_read=False)
write(O/f'EDA3_Region_Analysis_Summary_{D}.json',summary)
assert all(sha(p)==next(x['sha256'] for x in freeze['inputs'] if x['key']==k) for k,p in paths.items())
print(json.dumps(summary,ensure_ascii=False,indent=2))
