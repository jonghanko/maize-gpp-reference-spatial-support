"""Verify current publication figure sources and render to a fresh external folder."""
from pathlib import Path,PurePosixPath
import argparse,json,hashlib,subprocess,sys
p=argparse.ArgumentParser();p.add_argument('--runtime',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
module=a.runtime.resolve()/'publication_figures12';out=a.out.resolve()
if out.exists() or out==module or module in out.parents:raise ValueError('Use a fresh output outside publication figure inputs')
for r in json.loads((module/'SHA256_MANIFEST.json').read_text(encoding='utf-8'))['files']:
 rel=PurePosixPath(r['path'])
 if rel.is_absolute() or '..' in rel.parts:raise ValueError('Unsafe figure source path')
 f=module/r['path']
 if not f.is_file() or f.stat().st_size!=r['bytes'] or hashlib.sha256(f.read_bytes()).hexdigest()!=r['sha256']:raise ValueError('FIGURE_SOURCE_HASH_MISMATCH: '+r['path'])
subprocess.run([sys.executable,str(module/'run_figures.py'),'--out',str(out)],check=True)
