"""Read-only standard-library validator; content verification, not model replay."""
from pathlib import Path
import csv,hashlib,re,sys
r=Path(__file__).resolve().parents[1]
def rows(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def bad(k):return 'flip_units' in k or k.endswith('largest_change_deleted_donor') or bool(re.fullmatch(r'EDA53_(Pearson|bias|RMSE|centered_RMSE|sd_ratio)_(min|max|worst_departure)_donors',k))
def main():
 manifest=rows(r/'PUBLIC_MANIFEST.tsv');names=[x['path'] for x in manifest]
 assert len(names)==len(set(names));assert all(not Path(n).is_absolute() and '..' not in Path(n).parts for n in names)
 assert set(names)=={p.relative_to(r).as_posix() for p in r.rglob('*') if p.is_file() and p.name!='PUBLIC_MANIFEST.tsv' and '.git' not in p.relative_to(r).parts and '__pycache__' not in p.relative_to(r).parts},'exact file-set mismatch'
 for x in manifest:
  p=r/x['path'];assert p.stat().st_size==int(x['bytes']) and digest(p)==x['sha256'],x['path']
 genes={}
 for stem,n in [('Complete81',81),('Original43',43)]:
  a=rows(r/f'results/{stem}_GroupLevel_AllEvidence.tsv');assert len(a)==n and len(a[0])==1207
  assert not any(bad(k) for k in a[0]),'linkable identity field'
  assert len({x['HGNC_ID'] for x in a})==n;genes[stem]={x['HGNC_ID']:x for x in a}
 for k,x in genes['Original43'].items():assert genes['Complete81'][k]['direction']==x['direction']
 redact=rows(r/'docs/IDENTITY_COLUMN_REDACTIONS.tsv');assert len(redact)==74
 for stem in genes:assert sum(x['table']==stem for x in redact)==37
 f=rows(r/'docs/All57_Figure_Disclosure_Status.tsv');assert len(f)==57
 selected=[x for x in f if x['public_image_path']];assert len(selected)==43
 assert len([x for x in f if not x['public_image_path']])==14
 assert {p.relative_to(r).as_posix() for p in (r/'figures').rglob('*.png')}=={x['public_image_path'] for x in selected}
 for x in selected:assert x['single_image_review']=='True' and x['original_changed']=='False' and digest(r/x['public_image_path'])==x['sha256']
 t=rows(r/'docs/All211_Actual_Payload_Index_And_Local_Disclosure.tsv');assert len(t)==211
 selected=[x for x in t if x['public_payload_path']];assert len(selected)==37
 literal=[x for x in selected if x['transformation']=='LITERAL_BYTES'];assert len(literal)==35
 assert {p.relative_to(r).as_posix() for p in (r/'supplements').rglob('*') if p.is_file()}=={x['public_payload_path'] for x in literal}
 for x in selected:
  p=r/x['public_payload_path'];assert p.stat().st_size==int(x['public_bytes']) and digest(p)==x['public_sha256']
  if x in literal:assert x['public_sha256']==x['source_sha256']
 assert len(rows(r/'reproduction/All17_Route_Scope.tsv'))==17
 assert len(rows(r/'docs/All211_Supplement_Index_NOT_ALL_PAYLOAD.tsv'))==211
 assert len(rows(r/'docs/All57_Figure_Index.tsv'))==57
 assert len(rows(r/'diagnostics/All12_GBM_R2_N_EAF_PrimarySource_Adjudication.tsv'))==12
 assert len(rows(r/'docs/MR_Paragraph_To_Source_And_Disclosure_Crosswalk.tsv'))==10
 for x in rows(r/'PUBLIC_ASSET_PROVENANCE.tsv'):
  p=r/x['public_path'];assert p.stat().st_size==int(x['public_bytes']) and digest(p)==x['public_sha256']
 print('PASS: exact byte/SHA/file set; full81/43 1207cols;37identityfields excluded;43images;211scope/35literal+2derivatives;17route scope. NO model replay or new P/FDR.')
if __name__=='__main__':main()
