"""Reproduce all11 article figures from source-bound geometry and aggregate data."""
from pathlib import Path
from types import SimpleNamespace
import argparse,hashlib,importlib.util,json,shutil,sys
ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name):
 spec=importlib.util.spec_from_file_location(name,ROOT/'code'/(name+'.py'));mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();out=a.out.resolve()
 assert not out.exists() and not out.is_relative_to(ROOT),'Output must be new and outside immutable figure module'
 manifest=json.loads((ROOT/'FIGURE_INPUT_MANIFEST.json').read_text(encoding='utf-8'))
 for r in manifest['files']:assert sha(ROOT/r['path'])==r['sha256']
 before={r['path']:r['sha256'] for r in manifest['files']};out.mkdir(parents=True)
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 def reset():plt.rcdefaults();plt.rcParams['pdf.fonttype']=42
 reset();m=load('plot_design_geometry');sys.argv=['plot_design_geometry','--data',str(ROOT/'data/FIGURE_DATA.json'),'--out',str(out/'design')];m.main()
 reset();m=load('plot_conditional');m.main(SimpleNamespace(results=ROOT/'data',out=out/'conditional'))
 reset();m=load('plot_aggregate');m.OUT=out/'aggregate';m.OUT.mkdir();m.AUDIT=out/'aggregate_inputs';m.AUDIT.mkdir();m.FIGS=m.OUT
 for src,dst in [('TableS6_aggregation_estimands.csv','aggregation_metric_contrasts.csv'),('TableS7_site_season_contrasts.csv','LOYO_site_season_contrasts_all_families.csv'),('TableS3_reference_input_interactions.csv','own_reference_interactions.csv')]:shutil.copy2(ROOT/'data/tables'/src,m.AUDIT/dst)
 m.build_plots()
 reset();m=load('plot_omissions');m.ROOT=ROOT/'data';m.SOURCE=ROOT/'data/OMISSION_EXACT.json';m.CSV=ROOT/'data/tables/TableS13_fixed_OOF_summary_omissions.csv';m.HERE=out/'omission';m.HERE.mkdir();m.main()
 reset();load('plot_hindcast').main(ROOT/'data/tables',out/'hindcast')
 reset();aliases=out/'component_inputs';aliases.mkdir();shutil.copy2(ROOT/'data/tables/TableS22_component_summary_estimands.csv',aliases/'TableS23_component_summary_estimands.csv');load('plot_components').main(aliases,out/'components')
 reset();aliases=out/'support_inputs';aliases.mkdir()
 for src,dst in [('TableS23_monthly_weight_concentration.csv','TableS24_monthly_weight_concentration.csv'),('TableS25_site_season_input_summary.csv','TableS26_site_season_input_summary.csv')]:shutil.copy2(ROOT/'data/tables'/src,aliases/dst)
 load('plot_input_support').main(aliases,out/'support')
 final=out/'figures';final.mkdir()
 for folder in ['design','conditional','aggregate','omission','hindcast','components','support']:
  for f in (out/folder).glob('Figure*'):
   if f.suffix in ['.png','.pdf','.svg']:shutil.copy2(f,final/f.name)
 assert len(list(final.glob('*.png')))==len(list(final.glob('*.pdf')))==len(list(final.glob('*.svg')))==11
 assert all(sha(ROOT/rel)==h for rel,h in before.items())
 (out/'FIGURE_REPLAY_RECEIPT.json').write_text(json.dumps({'status':'pass','logical_figures':11,'formats':['png','pdf','svg'],'source_files_unchanged':len(before),'outputs':[{'path':p.relative_to(out).as_posix(),'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(final.iterdir())],'scope':'Original plotting bodies with public inputs and preserved priorities/styles/geometry; not a pixel-identity guarantee on arbitrary platforms. No model fit, prediction or score recomputation.','new_fits':0,'new_predictions':0},indent=2)+'\n',encoding='utf-8')
 print(json.dumps({'status':'pass','figures':11,'formats':3,'new_fits':0,'new_predictions':0}))
if __name__=='__main__':main()
