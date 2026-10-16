"""Author-FASTA global peptide specificity + outcome-blind AD aggregation/QC.

No gene disease contrasts, ANOVA outputs, or candidate selection by P.
"""
from pathlib import Path
import csv,json,re,bisect,collections,hashlib,openpyxl,numpy as np,datetime
R=Path(r'H:\SCI2\MRQI');O=R/'3_results/rp241_ad_gbm_b10d_p1_preparation';O.mkdir(parents=True,exist_ok=True)
def read(p):return list(csv.DictReader(p.open(encoding='utf-8'),delimiter='\t'))
def write(n,rows):
 rows=list(rows)
 with (O/n).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
m=read(R/'3_results/rp241_ad_gbm_b10d_p0_source_audit/All81_Fixed_HGNC_UniProt_Identity_2026-10-15.tsv');candidate={x['HGNC_ID'] for x in m}
h=read(R/'1_data/reference/gene_identity/hgnc_complete_set_2026-09-04.tsv');approved={x['symbol']:x['hgnc_id'] for x in h if x['status']=='Approved'};uni=collections.defaultdict(set)
for x in h:
 if x['status']!='Approved':continue
 for u in x['uniprot_ids'].split('|'):
  if u:uni[u].add(x['hgnc_id'])
fa=R/'1_data/proteomics/AD/PXD020296/IPI_HUMAN_042015_maptmtbrmod_abeta_apoe.fasta'
assert fa.stat().st_size==45363996 and hashlib.sha1(fa.read_bytes()).hexdigest()=='1f1bd1333c3ffbc68495ff007668106a7bc4246a'
entries=[];header=None;seq=[]
def add(header,seq):
 if header is None:return
 parts=header.split('|');acc=parts[1] if len(parts)>=3 else header.split()[0]
 gn=re.search(r'\bGN=(\S+)',header);gn=gn[1] if gn else ''
 stable=re.sub(r'-\d+$','',acc);by_uni=uni.get(stable,set());by_gn={approved[gn]} if gn in approved else set()
 conflict=bool(by_uni and by_gn and by_uni!=by_gn);gs=by_uni or by_gn
 gene=next(iter(gs)) if len(gs)==1 and not conflict else None
 entries.append({'accession':acc,'GN':gn,'gene':gene,'conflict':conflict,'sequence':''.join(seq).upper().replace('I','J').replace('L','J')})
for line in fa.open(encoding='utf-8'):
 if line.startswith('>'):add(header,seq);header=line[1:].strip();seq=[]
 else:seq.append(line.strip())
add(header,seq);acc_gene={e['accession']:e['gene'] for e in entries}
starts=[];parts=[];offset=0
for e in entries:starts.append(offset);parts.append(e['sequence']);offset+=len(e['sequence'])+1
database='|'.join(parts)
channels=read(R/'3_results/rp241_ad_gbm_b10d_p0_source_audit/AD_PXD020296_All33_Channel_Donor_Join_2026-10-15.tsv')
w=openpyxl.load_workbook(R/'1_data/proteomics/AD/PXD020296/dataanalysis.xlsx',read_only=True,data_only=True)
rows=[];all_values=[];all_GIS=[];qc_raw=collections.Counter();potential=collections.defaultdict(set)
for row in w['total peptide anova input '].iter_rows(min_row=5,values_only=True):
 if not row[0]:continue
 values=np.asarray([float(v) if v is not None else np.nan for v in row[41:74]])
 if np.isfinite(values).all():all_GIS.append([values[b*11]-values[b*11+10] for b in range(3)])
 accs={a.strip() for col in [row[13],row[18]] for a in str(col or '').split(';') if a.strip()};genes={acc_gene.get(a) for a in accs};targets=genes&candidate
 if not targets:continue
 sequence=str(row[3]).upper().replace('I','J').replace('L','J');mods=str(row[5] or '');q=float(row[8]) if row[8] is not None else np.nan
 rows.append({'peptide_id':row[0],'seq':sequence,'mods':mods,'potential_genes':targets,'source_accessions':';'.join(sorted(accs)),'q':q,'values':values})
 for g in targets:potential[g].add(sequence)
w.close();gis=np.asarray(all_GIS);mean=gis.mean(axis=0);sd=gis.std(axis=0,ddof=1)
specificity={};qs=sorted({x['seq'] for x in rows})
print('FASTA',len(entries),'candidate_rows',len(rows),'unique_queries',len(qs),flush=True)
for i,q in enumerate(qs):
 genes=set();unknown=False;hit_proteins=[];pos=database.find(q)
 while pos>=0:
  j=bisect.bisect_right(starts,pos)-1;gene=entries[j]['gene'];hit_proteins.append(entries[j]['accession'])
  if gene is None:unknown=True
  else:genes.add(gene)
  pos=database.find(q,pos+1)
 unique=next(iter(genes)) if len(genes)==1 and not unknown else None
 specificity[q]={'gene':unique,'all_genes':genes,'unknown':unknown,'hits':sorted(set(hit_proteins))}
 if (i+1)%100==0:print('FASTA sequence queries',i+1,'/',len(qs),flush=True)
ledger=[];retained=collections.defaultdict(lambda:collections.defaultdict(list));eligible_rows=0
for x in rows:
 s=specificity[x['seq']];v=x['values'];finite=np.isfinite(v).all();d=np.asarray([v[b*11]-v[b*11+10] for b in range(3)]);qc=bool(finite and np.all(abs(d-mean)<=2*sd));nonphos=not bool(re.search('phospho',x['mods'],re.I));fdr=bool(np.isfinite(x['q']) and x['q']<=.01)
 good=s['gene'] in candidate and qc and nonphos and fdr
 if good:retained[s['gene']][x['seq']].append(v);eligible_rows+=1
 ledger.append({'peptide_id':x['peptide_id'],'I_L_equivalent_sequence':x['seq'],'modifications':x['mods'],'source_potential_HGNC':';'.join(sorted(x['potential_genes'])),'full_FASTA_matched_accessions':';'.join(s['hits']),'full_FASTA_HGNC':';'.join(sorted(s['all_genes'])),'unknown_competing_entries':s['unknown'],'unique_HGNC':s['gene'] or '', 'nonphosphorylated':nonphos,'finite_all33':finite,'source_qvalue_le_001':fdr,'all3_GIS_within2SD':qc,'retained_for_gene_aggregation':good})
write('AD_All_Candidate_Peptide_FullFASTA_Specificity_QC_2026-10-15.tsv',ledger)
matrix=[];coverage=[]
for x in m:
 g=x['HGNC_ID'];seqs=retained[g];n=len(seqs);vals=np.median(np.asarray([np.median(np.asarray(vs),axis=0) for vs in seqs.values()]),axis=0) if n else np.full(33,np.nan)
 status='MEASURABLE_MULTI_PEPTIDE' if n>=2 else 'SINGLE_PEPTIDE_DESCRIPTIVE_ONLY' if n==1 else 'AMBIGUOUS_OR_TECHNICAL_FILTERED' if potential[g] else 'NO_ASSOCIATED_INPUT_PEPTIDE'
 coverage.append(dict(x,potential_sequence_n=len(potential[g]),retained_gene_specific_sequence_n=n,status=status,main_inference_eligible=n>=2))
 matrix.append({'HGNC_ID':g,'HGNC_symbol':x['HGNC_symbol'],'primary43':x['primary43'],'status':status,**{f'B{c["Batch"]}_{c["Channel"]}':float(vals[i]) if n>=2 else '' for i,c in enumerate(channels)}})
write('AD_All81_Prepared_Protein_Log2_GIS_Matrix_2026-10-15.tsv',matrix);write('AD_All81_Measurability_NoPrune_2026-10-15.tsv',coverage)
bio=[x for x in channels if x['Diagnosis'] in ['AD','control']];n=len(bio);design=np.asarray([[1,float(c['Diagnosis']=='AD'),float(c['Age']),float(c['Sex']=='female'),float(c['Batch']=='2'),float(c['Batch']=='3')] for c in bio]);pmides=np.column_stack([design,[float(c['PMI']) for c in bio]])
summary={'recorded_at_local':datetime.datetime.now().astimezone().isoformat(),'FASTA_entries':len(entries),'FASTA_gene_unresolved_entries':sum(e['gene'] is None for e in entries),'FASTA_annotation_conflicts':sum(e['conflict'] for e in entries),'candidate_associated_rows':len(rows),'global_sequence_queries':len(qs),'unique_nonphospho_QC_rows':eligible_rows,'GIS_all_source_finite_rows':len(gis),'GIS_batch_difference_mean':mean.tolist(),'GIS_batch_difference_SD_ddof1':sd.tolist(),'all81_measurable_multi_peptide':sum(x['main_inference_eligible'] for x in coverage),'primary43_measurable_multi_peptide':sum(x['main_inference_eligible'] and x['primary43']=='True' for x in coverage),'AD_primary_model_N':n,'AD_primary_model_columns':6,'AD_primary_model_rank':int(np.linalg.matrix_rank(design)),'AD_primary_model_residual_df':n-int(np.linalg.matrix_rank(design)),'AD_PMI_sensitivity_rank':int(np.linalg.matrix_rank(pmides)),'AD_PMI_sensitivity_residual_df':n-int(np.linalg.matrix_rank(pmides)),'normalization':'SOURCE_ALREADY_LOG2_NORMALIZED_ABUNDANCE_OVER_ARITHMETIC_GIS__NO_REPEAT_LOG','new_AD_contrast_fitted':False,'source_S6_received':False}
(O/'AD_Preparation_Specificity_Normalization_Model_QC_2026-10-15.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))
