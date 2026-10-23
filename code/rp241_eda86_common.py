from pathlib import Path
import csv, hashlib, json, re, io
from datetime import datetime, timezone
import os
R=Path(os.environ.get('MRQI_ROOT', str(Path.cwd()))); D='2026-10-24'
O=R/f'3_results/rp241_eda86_core_gap_public_archive_{D}'
V=R/f'0_admin/validation/rp241_eda86_core_gap_public_archive_{D}'
A=R/f'0_admin/source_audit/rp241_eda86_public_metadata_{D}'
P=R/'0_admin/protocol/RP241'
S=R/f'6_release/public_archive_candidates/RP241_EDA86_{D}'
for p in (O,V,A): p.mkdir(parents=True,exist_ok=True)
def utc(): return datetime.now(timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''): h.update(b)
 return h.hexdigest()
def rel(p):return Path(p).relative_to(R).as_posix()
def rec(p):p=Path(p);return dict(path=rel(p),bytes=p.stat().st_size,sha256=sha(p))
def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def put(p,s):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 if p.exists():
  assert p.read_text(encoding='utf-8-sig')==s,('Preserve differing artifact',str(p));return
 p.write_text(s,encoding='utf-8')
def js(p,x):put(p,json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def tab(p,x,fields=None):
 x=list(x);fields=fields or list(x[0]);b=io.StringIO(newline='');w=csv.DictWriter(b,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(x);put(p,b.getvalue())

