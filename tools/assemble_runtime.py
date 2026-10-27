"""Verify both release manifests and assemble an offline runtime in a new folder."""
from pathlib import Path,PurePosixPath
import argparse,json,hashlib,shutil
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(base):
 m=json.loads((base/'RELEASE_MANIFEST.json').read_text(encoding='utf-8'))
 for r in m['files']:
  rel=PurePosixPath(r['path'])
  if rel.is_absolute() or '..' in rel.parts:raise ValueError('Unsafe release path')
  p=base/r['path']
  if not p.is_file() or p.stat().st_size!=r['bytes'] or digest(p)!=r['sha256']:raise ValueError('RELEASE_HASH_MISMATCH: '+r['path'])
 return len(m['files'])
def main():
 p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 code=Path(__file__).resolve().parents[1];data=a.data.resolve();out=a.out.resolve()
 if out.exists() or any(out==r or r in out.parents for r in [code,data]):raise ValueError('Use a new output outside both input records')
 counts={'software':verify(code),'data':verify(data)}
 out.mkdir(parents=True,exist_ok=False)
 for src in [code/'runtime',data/'runtime']:
  for old in sorted(src.rglob('*')):
   if not old.is_file():continue
   new=out/old.relative_to(src)
   if new.exists():
    if digest(new)!=digest(old):raise ValueError('Conflicting component files')
   else:new.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(old,new)
 for comp in ['baseline','crossprediction','figures']:
  m=json.loads((out/comp/'MANIFEST.json').read_text(encoding='utf-8'))
  for r in m['files']:
   f=out/comp/r['path']
   if digest(f)!=r['sha256']:raise ValueError('ASSEMBLED_MANIFEST_MISMATCH: '+r['path'])
 result={'status':'PASS','source_records_verified':counts,'new_fits':0,'new_predictions':0}
 (out/'ASSEMBLY_RECEIPT.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))
if __name__=='__main__':main()
