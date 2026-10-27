"""Render 20 already-computed summary sensitivities; never fit or rescore predictions."""
from pathlib import Path
from fractions import Fraction
from decimal import Decimal,localcontext
import csv,hashlib,json,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
SOURCE=ROOT/'FIXED_OOF_OMISSION_SENSITIVITY.json'
CSV=ROOT/'summary_sensitivity.csv'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def exact(o):return Fraction(int(o['numerator']),int(o['denominator']))
def decimal17(q):
 with localcontext() as c:
  c.prec=40
  return format(Decimal(q.numerator)/Decimal(q.denominator),'.17g')
def main():
 names=['FigureS1_fixed_OOF_summary_omissions.png','FigureS1_fixed_OOF_summary_omissions.svg','FigureS1_plot_values.csv','FIGURE_VERIFICATION.json']
 assert not any((HERE/n).exists() for n in names),'Figure outputs already exist; preserve the first render.'
 before={p.name:sha(p) for p in [SOURCE,CSV]}
 data=json.loads(SOURCE.read_text())['all_summaries']
 ix={(r['model'],r['reference'],r['scenario']):r for r in data if r['contrast']=='C_minus_B'}
 order=[('ridge','VUT_REF'),('ridge','CUT_REF'),('random_forest','VUT_REF'),('random_forest','CUT_REF')]
 scenarios=['full_four_site_seasons','omit_US-Ne1-2012','omit_US-Ne1-2013','omit_US-Ne3-2013','omit_US-Ro1-2015']
 rownames=['Ridge · VUT','Ridge · CUT','RF · VUT','RF · CUT']
 colnames=['Full four\nsite-seasons','Omit\nNe1 2012','Omit\nNe1 2013','Omit\nNe3 2013','Omit Ro1 2015\nMead evaluation only']
 values=[];records=[]
 for model,ref in order:
  row=[]
  for scenario in scenarios:
   r=ix[(model,ref,scenario)];q=exact(r['relative_difference_percent']);row.append(float(q))
   records.append({'model':model,'reference':ref,'contrast':'C_minus_B','scenario':scenario,
     'relative_difference_percent':decimal17(q),'exact_numerator':q.numerator,'exact_denominator':q.denominator,
     'n_retained_site_seasons':r['n_retained_site_seasons'],'n_retained_regions':r['n_retained_regions']})
  values.append(row)
 with CSV.open(encoding='utf-8-sig',newline='') as f:scsv=list(csv.DictReader(f))
 ci={(r['model'],r['reference'],r['scenario']):r for r in scsv if r['contrast']=='C_minus_B'}
 assert all(abs(float(r['relative_difference_percent'])-float(ci[(r['model'],r['reference'],r['scenario'])]['relative_difference_percent']))<1e-12 for r in records)
 assert len(records)==20 and all(len(row)==5 for row in values)
 assert all(ix[(m,r,scenarios[-1])]['n_retained_regions']==1 for m,r in order)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
 fig=plt.figure(figsize=(12.6,5.3))
 ax=fig.add_axes([.14,.27,.71,.48]);cax=fig.add_axes([.88,.29,.021,.44])
 # One symmetric colour scale across all rows and columns; zero is fixed.
 extent=math.ceil(max(abs(x) for row in values for x in row)*10)/10
 cmap=LinearSegmentedColormap.from_list('C_B_signed',['#2166ac','#fafafa','#d6601f'])
 image=ax.imshow(values,cmap=cmap,norm=TwoSlopeNorm(vmin=-extent,vcenter=0,vmax=extent),aspect='auto')
 ax.set_xticks(range(5),colnames,fontsize=10);ax.set_yticks(range(4),rownames,fontsize=11)
 ax.tick_params(axis='both',length=0,pad=9)
 ax.set_xticks([-.5,.5,1.5,2.5,3.5,4.5],minor=True);ax.set_yticks([-.5,.5,1.5,2.5,3.5],minor=True)
 ax.grid(which='minor',color='white',linewidth=2);ax.tick_params(which='minor',length=0)
 for i,row in enumerate(values):
  for j,val in enumerate(row):
   ax.text(j,i,f'{val:+.3f}%',ha='center',va='center',fontsize=11,
     color='white' if abs(val)>.65*extent else '#101820')
 for spine in ax.spines.values():spine.set_visible(False)
 colorbar=fig.colorbar(image,cax=cax);colorbar.set_label('Relative RMSE difference (%)\nC minus B',labelpad=11,fontsize=10)
 fig.suptitle('C − B sensitivity to site-season omissions from the error summary',fontsize=15,y=.95)
 fig.text(.5,.84,'Fixed LOYO out-of-fold outputs · blue: lower C error · orange: higher C error',ha='center',fontsize=10)
 fig.text(.5,.095,'Evaluation-summary omissions only: no retraining, new validation or uncertainty intervals.\nOmitting Ro1 leaves only Mead in evaluation; its labels remain in the original training folds.',ha='center',va='center',fontsize=10)
 fig.savefig(HERE/names[0],dpi=220,facecolor='white');fig.savefig(HERE/names[1],facecolor='white');fig.savefig(HERE/'FigureS1_fixed_OOF_summary_omissions.pdf',facecolor='white');plt.close(fig)
 with (HERE/names[2]).open('x',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
 after={p.name:sha(p) for p in [SOURCE,CSV]};assert before==after
 verification={'status':'pass','rows':4,'columns':5,'cells':20,'source_sha256_before':before,'source_sha256_after':after,
  'plot_code_sha256':sha(Path(__file__)),'all_exact_JSON_and_csv_values_agree':True,'shared_zero_centred_colour_extent':[-extent,extent],
  'negative_colour':'blue; favours C','positive_colour':'orange; favours B','new_model_fits':0,'new_validation':False,
  'rendering_only_float_conversion':True,'statistical_intervals':False,
  'outputs':{n:{'sha256':sha(HERE/n),'bytes':(HERE/n).stat().st_size} for n in names[:3]}}
 with (HERE/names[3]).open('x',encoding='utf-8') as f:json.dump(verification,f,indent=2);f.write('\n')
 print(json.dumps({'status':'pass','cells':20,'new_fits':0,'extent':extent}))

if __name__=='__main__':main()
