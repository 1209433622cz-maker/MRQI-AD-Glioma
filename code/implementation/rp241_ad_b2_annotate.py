"""Freeze all-row numerical gene identities and outcome-blind DE universe."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import argparse
import csv
import ctypes
import gzip
import hashlib
import json
import re
from collections import Counter,defaultdict
from pathlib import Path

import numpy as np

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()

def read_tsv(p):
 with p.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def table(p,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,keys,delimiter='\t');w.writeheader();w.writerows(rows)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--stamp',required=True);a=ap.parse_args()
 root=a.root.resolve();stamp=a.stamp
 if stamp=='2026-10-13' and (root/'6_release/phase_packages/MRQI_RP241_AD_B2_B2E_Annotation_CountAware_Exploratory_delta_2026-10-14.zip').exists():
  raise RuntimeError('Published B2 identity/model preparation is frozen. New annotation requires a new stage/version.')
 if hasattr(ctypes,'windll'):ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
 out=root/'3_results/rp241_ad_b2_annotation_model';out.mkdir(parents=True,exist_ok=True)
 hgnc_path=root/'1_data/reference/gene_identity/hgnc_complete_set_2026-09-04.tsv'
 ncbi_path=root/f'1_data/reference/gene_identity/ncbi_gene_{stamp}/Homo_sapiens.gene_info.gz'
 hgnc=read_tsv(hgnc_path);by_entrez=defaultdict(list)
 for r in hgnc:
  if r['status']=='Approved' and r['entrez_id']:by_entrez[r['entrez_id']].append(r)
 with gzip.open(ncbi_path,'rt',encoding='utf-8') as f:
  all_ncbi_rows=list(csv.DictReader(f,delimiter='\t'))
 source_taxa=dict(Counter(r['#tax_id'] for r in all_ncbi_rows))
 ncbi_rows=[r for r in all_ncbi_rows if r['#tax_id']=='9606']
 assert len({r['GeneID'] for r in ncbi_rows})==len(ncbi_rows)
 ncbi={r['GeneID']:r for r in ncbi_rows}
 # Identity index is independent of current GBM quantitative outcomes.
 gtf_path=root/'1_data/reference/gencode_v27_GRCh38/gencode.v27.annotation.gtf.gz'
 gencode=defaultdict(list);gtf_records=[]
 with gzip.open(gtf_path,'rt',encoding='utf-8') as f:
  for line in f:
   if '\tgene\t' not in line:continue
   cols=line.rstrip('\n').split('\t');attrs=dict(re.findall(r'(\w+) "([^"]+)"',cols[8]))
   versioned=attrs['gene_id'];base=re.sub(r'\.\d+(?=_PAR_Y$|$)','',versioned)
   row={'Ensembl_base_ID':base,'GENCODE27_gene_ID':versioned,'gene_name_source':attrs.get('gene_name',''),
        'gene_type_source':attrs.get('gene_type',''),'chromosome':cols[0],'start':int(cols[3]),'end':int(cols[4]),
        'strand':cols[6],'PAR_Y_alternative':versioned.endswith('_PAR_Y')}
   gencode[base].append(row);gtf_records.append(row)
 assert len(gtf_records)==58288
 table(out/f'GENCODE27_Exact_Gene_Identity_Index_{stamp}.tsv',gtf_records)
 matrix_path=root/'3_results/rp241_ad_b1_source_use_gate/GSE125583_Entrez_Raw_Integer_Counts.tsv.gz'
 with gzip.open(matrix_path,'rt',encoding='utf-8') as f:
  rd=csv.reader(f,delimiter='\t');gsms=next(rd)[1:];data=list(rd)
 genes=[r[0] for r in data];counts=np.array([r[1:] for r in data],dtype=np.int64);del data
 assert counts.shape==(30727,289) and len(set(genes))==30727
 members=read_tsv(root/'3_results/rp241_ad_b1_source_use_gate/GSE125583_M0_CLINICAL_STUDY_LABEL_PRIMARY_Design_Matrix_2026-10-13.tsv')
 lookup={gsm:i for i,gsm in enumerate(gsms)};indices=[lookup[r['gsm']] for r in members]
 assert len(indices)==288
 detected=(counts[:,indices]>=10).sum(axis=1)
 measured=detected>=70
 # A missing history reference never changes absent/current IDs into a retired assertion.
 history_path=root/f'1_data/reference/gene_identity/ncbi_gene_{stamp}/Homo_sapiens.gene_history.tsv.gz'
 history=defaultdict(list)
 if history_path.exists():
  with gzip.open(history_path,'rt',encoding='utf-8') as f:
   for r in csv.DictReader(f,delimiter='\t'):history[r['Discontinued_GeneID']].append(r)
 annotations=[]
 for i,eid in enumerate(genes):
  hh=by_entrez.get(eid,[]);nn=ncbi.get(eid);current=eid in ncbi
  h=hh[0] if len(hh)==1 else None
  current_hgnc=sorted(set(re.findall(r'HGNC:HGNC:(\d+)',nn['dbXrefs']))) if nn else []
  current_hgnc=['HGNC:'+x for x in current_hgnc]
  current_ens=sorted(set(re.findall(r'Ensembl:(ENSG\d+)',nn['dbXrefs']))) if nn else []
  ens=h['ensembl_gene_id'] if h else ''
  ens_valid=bool(re.fullmatch(r'ENSG\d+',ens))
  gg=gencode.get(ens,[]) if ens_valid else []
  identity_issues=[]
  if not current:identity_issues.append('NOT_IN_CURRENT_NCBI_HUMAN_INFO')
  if len(hh)==0:identity_issues.append('NO_EXACT_HGNC_ENTREZ')
  elif len(hh)>1:identity_issues.append('MULTIPLE_HGNC_ENTREZ_HITS')
  if h and (not current_hgnc or current_hgnc!=[h['hgnc_id']]):identity_issues.append('NCBI_HGNC_MISSING_OR_CONFLICT')
  if not ens_valid:identity_issues.append('NO_SINGLE_HGNC_ENSEMBL')
  if current_ens and current_ens!=[ens]:identity_issues.append('NCBI_ENSEMBL_MULTIPLE_OR_CONFLICT')
  if len(gg)!=1:identity_issues.append('GENCODE27_NOT_UNIQUE_OR_NOT_PRESENT')
  raw_history=history.get(eid,[])
  if current:history_status='CURRENT_EXACT_NCBI_GENE_ID'
  elif not history_path.exists():history_status='HISTORY_NOT_FROZEN__RETIREMENT_NOT_ASSERTED'
  elif len(raw_history)==1:history_status='DISCONTINUED_RECORD__NO_COUNTS_ALIAS_OR_MERGE_APPLIED'
  elif len(raw_history)>1:history_status='MULTIPLE_HISTORY_RECORDS__NO_ALIAS_APPLIED'
  else:history_status='NO_CURRENT_OR_HISTORY_HIT'
  r={'original_Entrez_ID':eid,'raw_row_index_1based':i+1,'NCBI_current_human_record':current,
     'NCBI_symbol':nn['Symbol'] if nn else '', 'NCBI_type':nn['type_of_gene'] if nn else '',
     'HGNC_exact_hits':len(hh),'HGNC_ID':h['hgnc_id'] if h else '', 'HGNC_symbol_frozen':h['symbol'] if h else '',
     'HGNC_status':h['status'] if h else '', 'HGNC_locus_type':h['locus_type'] if h else '',
     'HGNC_Ensembl_gene_ID':ens,'NCBI_HGNC_IDs':'|'.join(current_hgnc),'NCBI_Ensembl_IDs':'|'.join(current_ens),
     'GENCODE27_records_for_base_ID':len(gg),'GENCODE27_versioned_ID':gg[0]['GENCODE27_gene_ID'] if len(gg)==1 else '',
     'GENCODE27_gene_type':gg[0]['gene_type_source'] if len(gg)==1 else '',
     'symbol_name_difference_preserved':bool(h and nn and h['symbol']!=nn['Symbol']),
     'history_status':history_status,'history_replacements':'|'.join(x['GeneID'] for x in raw_history),
     'numeric_alias_applied':False,'counts_merged':False,
     'primary_samples_count_ge10':int(detected[i]),'DE_measurable_primary_count_rule':bool(measured[i]),
     'raw_total_read_count':int(counts[i,:].sum()),'identity_issues':'|'.join(identity_issues),
     'strict_cross_disease_identity_eligible':not identity_issues,'HGNC_or_Ensembl_collision':False}
  annotations.append(r)
 hgnc_targets=defaultdict(list);ens_targets=defaultdict(list)
 for r in annotations:
  if r['HGNC_ID']:hgnc_targets[r['HGNC_ID']].append(r['original_Entrez_ID'])
  if r['HGNC_Ensembl_gene_ID']:ens_targets[r['HGNC_Ensembl_gene_ID']].append(r['original_Entrez_ID'])
 for r in annotations:
  collision=bool((r['HGNC_ID'] and len(hgnc_targets[r['HGNC_ID']])>1) or
                 (r['HGNC_Ensembl_gene_ID'] and len(ens_targets[r['HGNC_Ensembl_gene_ID']])>1))
  if collision:
   r['HGNC_or_Ensembl_collision']=True;r['strict_cross_disease_identity_eligible']=False
   r['identity_issues']+='|MULTIPLE_RAW_ENTREZ_ROWS_SHARE_HGNC_OR_ENSEMBL'
 table(out/f'GSE125583_All_Raw_Row_Identity_Ledger_{stamp}.tsv',annotations)
 table(out/f'GSE125583_Identity_Unresolved_or_Conflict_Rows_{stamp}.tsv',[r for r in annotations if not r['strict_cross_disease_identity_eligible']])
 table(out/f'GSE125583_PreEffect_Primary_Measurable_Universe_{stamp}.tsv',[r for r in annotations if r['DE_measurable_primary_count_rule']])
 table(out/f'GSE125583_Strict_GENCODE27_Bridge_Measurable_Rows_{stamp}.tsv',[r for r in annotations if r['DE_measurable_primary_count_rule'] and r['strict_cross_disease_identity_eligible']])
 summary={'date':stamp,'raw_ID_rows':len(annotations),'NCBI_source_taxa':source_taxa,'NCBI_taxonomy_filter':'#tax_id==9606 exactly','NCBI_human_reference_rows':len(ncbi_rows),'HGNC_approved_reference_rows':len(hgnc),
          'GENCODE27_gene_records':len(gtf_records),'GENCODE27_PAR_Y_records':sum(r['PAR_Y_alternative'] for r in gtf_records),
          'GENCODE27_basic_ID_conflicts':sum(len(v)>1 for v in gencode.values()),
          'original_ID_in_current_NCBI':sum(r['NCBI_current_human_record'] for r in annotations),
          'original_ID_absent_current_NCBI':sum(not r['NCBI_current_human_record'] for r in annotations),
          'HGNC_exact_one_hit':sum(r['HGNC_exact_hits']==1 for r in annotations),
          'strict_cross_disease_rows_all':sum(r['strict_cross_disease_identity_eligible'] for r in annotations),
          'primary_measurable_universe_rows':sum(r['DE_measurable_primary_count_rule'] for r in annotations),
          'strict_bridge_primary_measurable_rows':sum(r['DE_measurable_primary_count_rule'] and r['strict_cross_disease_identity_eligible'] for r in annotations),
          'original_ID_collision_rows':sum(r['HGNC_or_Ensembl_collision'] for r in annotations),
          'history_source_frozen':history_path.exists(),'history_aliases_applied':0,'counts_rows_merged_or_dropped':0,
          'filter':'count>=10 in>=70 of288 complete-covariate primary samples',
          'gene_effects_or_P_values_used_in_filter':False,'primary_DE_family_raw_Entrez_rows':int(measured.sum()),
          'DE_run':False,'reference_SHA256':{'HGNC':sha(hgnc_path),'NCBI_gene_info':sha(ncbi_path),'GENCODE27_GTF':sha(gtf_path),'B1_counts':sha(matrix_path)}}
 (out/f'GSE125583_Annotation_Measurable_Universe_Receipt_{stamp}.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
 print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
