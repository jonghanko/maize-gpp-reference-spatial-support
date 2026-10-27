"""Render frozen numerical figure inputs only; no model fitting or prediction."""
from pathlib import Path
import argparse,subprocess,sys
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--groups',nargs='+',choices=['main','contrasts','supp'],default=['main','contrasts','supp'])
p.add_argument('--out',type=Path,default=ROOT/'output')
p.add_argument('--verify-contrasts',action='store_true',help='Requires pypdf and pdfplumber')
a=p.parse_args();out=a.out.resolve()
commands={'main':[ROOT/'code/main12/plot_main12.py','--source',ROOT/'data/main12','--out',out/'main'],
'contrasts':[ROOT/'code/contrasts12/plot_contrasts12.py','--source',ROOT/'data/contrasts12','--out',out/'contrasts'],
'supp':[ROOT/'code/supp12/build_supp12.py','--source',ROOT/'data/supp12','--out',out/'supp']}
for name in a.groups:
    subprocess.run([sys.executable,*map(str,commands[name])],check=True)
if a.verify_contrasts:
    subprocess.run([sys.executable,str(ROOT/'code/contrasts12/verify_contrasts12.py'),'--source',str(ROOT/'data/contrasts12'),'--out',str(out/'contrasts')],check=True)
print('PASS: figure rendering only; no model calls')
