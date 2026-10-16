"""Metadata-only source join. No candidate expression contrast is fitted."""
from pathlib import Path
import json, hashlib, gzip
import numpy as np
import pandas as pd
import h5py

R=Path('H:/SCI2/MRQI');D='2026-10-15'
O=R/'3_results/rp241_ad_gbm_b10d_p3_spatial_source';O.mkdir(parents=True,exist_ok=True)
B=R/'1_data/spatial/GSE237183'
q=pd.read_csv(R/f'3_results/rp241_ad_gbm_b10d_p1_preparation/GSE237183_All13_Counts_Barcode_Position_QC_{D}.tsv',sep='\t')
codes=['MGH258','ZH881inf','ZH881T1','ZH916bulk','ZH916inf','ZH916T1','ZH1007inf','ZH1007nec','ZH1019inf','ZH1019T1','ZH8811Abulk','ZH8811Bbulk','ZH8812bulk']
mapping={f'GSM{7596587+i}':x for i,x in enumerate(codes)}
assert set(q.GSM)==set(mapping)
meta=pd.read_csv(B/'author_inputs_metadata/Inputs/general/visium_metadata.csv',keep_default_na=False)
assert not meta.duplicated(['sample','spot_id']).any()
final=pd.read_csv(O/'Author_final_zones_All26.tsv',sep='\t',keep_default_na=False).rename(columns={'label':'final_zone'})
org=pd.read_csv(O/'Author_spatial_zonesv3_All26.tsv',sep='\t',keep_default_na=False).rename(columns={'label':'organisation_zone'})
mp=pd.read_csv(O/'Author_MP_Assignments_All13.tsv',sep='\t',keep_default_na=False).rename(columns={'mp':'MP_RDS'})
author=meta.merge(final,on=['sample','spot_id'],how='outer',validate='1:1',indicator='metadata_final_join')
author=author.merge(org,on=['sample','spot_id'],how='outer',validate='1:1')
author=author.merge(mp,on=['sample','spot_id'],how='outer',validate='1:1')
readiness=[];allspots=[]
for _,r in q.iterrows():
 sample=mapping[r.GSM];p=B/'processed'/r.H5_name
 with h5py.File(p) as f:bc=[x.decode() for x in f['matrix/barcodes'][:]]
 assert len(set(bc))==len(bc)==r.filtered_spots
 posfile=next((B/'processed').glob(r.GSM+'*tissue_positions_list.csv.gz'))
 positions=pd.read_csv(posfile,header=None,names=['spot_id','in_tissue','array_row','array_col','pixel_row','pixel_col'])
 raw=pd.DataFrame({'sample':sample,'spot_id':bc,'GSM':r.GSM,'patient':r.patient,'specimen_region':r.author_region})
 raw=raw.merge(positions,on='spot_id',how='left',validate='1:1')
 aa=author.loc[author['sample'].eq(sample)]
 joined=raw.merge(aa,on=['sample','spot_id'],how='left',validate='1:1',indicator='local_author_join')
 joined['author_metadata_present']=joined.metadata_final_join.eq('both')
 joined['author_final_zone_present']=joined.final_zone.notna() & joined.final_zone.ne('')
 joined['author_MP_RDS_present']=joined.MP_RDS.notna() & joined.MP_RDS.ne('')
 have=joined.author_metadata_present
 coords=(joined.loc[have,'centroid_x'].eq(joined.loc[have,'array_col']) & joined.loc[have,'centroid_y'].eq(-joined.loc[have,'array_row']))
 zmatch=joined.loc[have,'org2'].eq(joined.loc[have,'final_zone'])
 mm=joined.loc[have & joined.author_MP_RDS_present,'mp'].eq(joined.loc[have & joined.author_MP_RDS_present,'MP_RDS'])
 readiness.append(dict(GSM=r.GSM,sample=sample,patient=r.patient,specimen_region=r.author_region,
   local_filtered_spots=len(bc),author_metadata_spots=int(have.sum()),author_final_zone_spots=int(joined.author_final_zone_present.sum()),
   author_MP_RDS_spots=int(joined.author_MP_RDS_present.sum()),local_unannotated_spots=int((~have).sum()),
   author_only_barcodes=int((~aa.spot_id.isin(bc)).sum()),coordinate_transform='centroid_x=array_col;centroid_y=-array_row',
   coordinate_disagreements=int((~coords).sum()),CSV_org2_vs_final_RDS_disagreements=int((~zmatch).sum()),
   CSV_mp_vs_MP_RDS_disagreements=int((~mm).sum()),age='NOT_RESOLVED',individual_treatment='NOT_RESOLVED',
   old_GBM_patient_identity_bridge='NOT_RESOLVED',new_candidate_spatial_effects=False))
 allspots.append(joined)
rr=pd.DataFrame(readiness);ss=pd.concat(allspots,ignore_index=True)
assert len(ss)==25488 and len(rr)==13 and rr.patient.nunique()==5
rr.to_csv(O/f'All13_Exact_Author_Annotation_Join_Readiness_{D}.tsv',sep='\t',index=False)
ss.to_csv(O/f'All25488_Local_Spots_Author_Join_NoDrop_{D}.tsv',sep='\t',index=False,na_rep='')
presence=(ss.groupby(['patient','sample','specimen_region','final_zone'],dropna=False).size().rename('spots').reset_index())
presence.to_csv(O/f'All5_Patient_Specimen_AuthorZone_Presence_{D}.tsv',sep='\t',index=False,na_rep='UNANNOTATED')
genes=pd.read_csv(O/'Author_Label_Defining_MP_Genes.tsv',sep='\t')
ident=pd.read_csv(R/'3_results/rp241_ad_gbm_b10_source_model_freeze/All81_Immutable_Identity_Direction_Input_2026-10-14.tsv',sep='\t')
overlap=[]
for _,g in ident.iterrows():
 matches=genes.loc[genes.gene.eq(g.HGNC_symbol)]
 overlap.append(dict(HGNC_ID=g.HGNC_ID,HGNC_symbol=g.HGNC_symbol,primary43=g.primary43,direction=g.direction,
   author_label_defining_MP_exact_memberships=';'.join(matches.metaprogram),membership_n=len(matches),
   label_derived_from_same_RNA_platform=True,independent_spatial_gene_validation_established=False))
ov=pd.DataFrame(overlap);ov.to_csv(O/f'All81_Author_MP_Gene_Overlap_Circularity_Audit_{D}.tsv',sep='\t',index=False)
sm=dict(all_local_spots=25488,author_metadata_all26_rows=len(meta),local_samples=13,local_patients=5,
 annotated_local_spots=int(rr.author_metadata_spots.sum()),unannotated_local_spots=int(rr.local_unannotated_spots.sum()),
 author_only_barcodes=int(rr.author_only_barcodes.sum()),coordinate_disagreements=int(rr.coordinate_disagreements.sum()),
 org2_vs_final_RDS_disagreements=int(rr.CSV_org2_vs_final_RDS_disagreements.sum()),
 mp_CSV_vs_RDS_disagreements=int(rr.CSV_mp_vs_MP_RDS_disagreements.sum()),
 main43_MP_overlap_genes=ov.loc[ov.primary43 & ov.membership_n.gt(0),'HGNC_symbol'].tolist(),
 full81_MP_overlap_genes=ov.loc[ov.membership_n.gt(0),'HGNC_symbol'].tolist(),
 author_spot_annotation_resolved=True,age_treatment_closed=False,patient_overlap_closed=False,
 spatial_S6_received=False,new_candidate_spatial_effects=False)
(O/f'Spatial_Source_Closure_Summary_{D}.json').write_text(json.dumps(sm,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(sm,ensure_ascii=False,indent=2))
