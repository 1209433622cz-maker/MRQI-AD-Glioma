"""Read-only EDA91 document/attachment checks. No scientific models are run."""
from pathlib import Path
import csv,hashlib,json,re,zipfile,xml.etree.ElementTree as ET
R=Path(__file__).resolve().parents[1]
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def paragraphs(j):
 out=j.get('abstract_zh',[])+j.get('abstract_en',[])
 def walk(ss):
  for s in ss:out.extend(s.get('paragraphs',[]));walk(s.get('subsections',[]))
 walk(j['sections']);return out
def main():
 v=json.loads((R/'PUBLICATION_VERSION.json').read_text('utf-8'));assert v['version'].startswith('RP241_EDA91_') and not v['new_scientific_analysis'] and not v['causal_upgrade']
 j=json.loads((R/'docs/CURRENT_MANUSCRIPT_STRUCTURED.json').read_text('utf-8'));p=paragraphs(j);assert len(p)==104 and len({x['id'] for x in p})==104
 md=(R/'docs/CURRENT_MANUSCRIPT_WORKING_DRAFT.md').read_text('utf-8');doc=R/v['public_manuscript_word'];assert sha(doc)==v['manuscript_word_sha256']
 with zipfile.ZipFile(doc) as z:
  assert z.testzip() is None;root=ET.fromstring(z.read('word/document.xml'))
 ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
 texts=[''.join(n.itertext()) for n in root.findall('./w:body/w:p',ns)]
 for a in p:assert md.count(a['text'])==1 and sum(t.endswith(a['text']) for t in texts)==1,a['id']
 assert len(j['references'])==29
 seen=[]
 for a in p:
  for token in re.findall(r'\[([\d,]+)\]',a['text']):
   for x in map(int,token.split(',')):
    assert 1<=x<=29
    if x not in seen:seen.append(x)
 assert seen==list(range(1,30))
 f=rows(R/'docs/EDA91_ACTUAL_FIGURE_SET.tsv');ix={x['figure_id']:x for x in rows(R/'docs/All57_Figure_Index.tsv')};d={x['figure_id']:x for x in rows(R/'docs/All57_Figure_Disclosure_Status.tsv')}
 assert len(f)==55 and {x['figure_id'] for x in f}==set(ix)-{'S1','S5'}
 assert sum(x['figure_id'].startswith('M') for x in f)==12
 assert sum(bool(x['public_payload_path']) for x in f)==41 and sum(not x['public_payload_path'] for x in f)==14
 for x in f:
  assert x['body_paragraphs'] and all(a in {p['id'] for p in p} for a in x['body_paragraphs'].split(';'))
  if x['public_payload_path']:assert sha(R/x['public_payload_path'])==x['sha256']==d[x['figure_id']]['sha256']
  else:assert x['release_scope']=='LOCAL_REVIEW_ONLY_NOT_AUTO_SUBMITTED_OR_PUBLIC'
 assert len(rows(R/'docs/EDA91_HISTORICAL_BACKGROUND_FIGURES.tsv'))==2
 fix={'S36':'E48M03;E48R01','S37':'E48M03;E48R02','S38':'E48M04;E48R03','S39':'E48M03;E48R02','S33':'E48M01;E48R04'}
 assert all(ix[k]['manuscript_paragraph_ids']==t for k,t in fix.items())
 t=rows(R/'docs/EDA91_ACTUAL_TABLE_SET.tsv');original=rows(R/'docs/All211_Actual_Payload_Index_And_Local_Disclosure.tsv')
 assert len(t)==211 and {x['table_id'] for x in t}=={x['table_id'] for x in original}
 assert sum(bool(x['public_payload_path']) for x in t)==37 and sum(x['public_transform']=='LITERAL_BYTES' for x in t)==35
 for x in t:
  if x['public_payload_path']:assert sha(R/x['public_payload_path'])==x['public_sha256']
 bindings=rows(R/'docs/EDA91_UNCHANGED_SCIENTIFIC_PAYLOADS.tsv');assert len(bindings)==94
 for x in bindings:assert (R/x['path']).stat().st_size==int(x['bytes']) and sha(R/x['path'])==x['sha256']
 legends=(R/'docs/EDA91_CURRENT_SUBMISSION_LEGENDS.md').read_text('utf-8');assert len(re.findall(r'^## ',legends,re.M))==55 and not re.search(r'^## 补图S(?:1|5) ',legends,re.M)
 print('PASS EDA91:104 Word/MD/JSON paragraphs;29 references;55 actual figures=41public+14local;211 table roles=35literal+2derivative+174local;94 scientific payloads unchanged. NO journal submission, model replay or causal upgrade.')
if __name__=='__main__':main()
