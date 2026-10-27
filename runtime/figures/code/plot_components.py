"""Plot exact fixed-output components; reads aggregate CSV only."""
from pathlib import Path
import argparse,json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator,FormatStrFormatter
import pandas as pd
import numpy as np
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main(run,out):
 run=Path(run);out=Path(out);out.mkdir(exist_ok=False)
 path=run/'TableS23_component_summary_estimands.csv';df=pd.read_csv(path)
 d=df.loc[df.scope.eq('equal_site_season')&df.fit_reference.eq(df.score_reference)]
 panels=[('principal_642','C_core_weighted__minus__B_uniform_pixels','(a) Principal cohort: C − B'),('principal_642','B_uniform_pixels__minus__A_shared_scene','(b) Principal cohort: B − A shared'),('T1_437','C_core_weighted__minus__B_uniform_pixels','(c) Earlier-year Mead: C − B')]
 order=[(m,r) for m in ['ridge','random_forest','hist_gradient_boosting'] for r in ['VUT_REF','CUT_REF']]
 vals=[];frames=[]
 for cohort,contrast,title in panels:
  q=d.loc[d.cohort.eq(cohort)&d.contrast.eq(contrast)].set_index(['model','fit_reference']).loc[order]
  frames.append(q);vals.extend(q[['mean_bias_rmse_component','mean_centered_rmse_component','mean_delta_rmse']].to_numpy().ravel())
 lo=min(vals)*1.15;hi=max(.02,max(vals)*1.6)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.5,'axes.titlesize':11.5,'svg.fonttype':'none','savefig.facecolor':'white'})
 fig,axes=plt.subplots(1,3,figsize=(15.3,5.2),sharex=True)
 labels=[{'ridge':'Ridge','random_forest':'RF','hist_gradient_boosting':'HGB'}[m]+' · '+r.replace('_REF','') for m,r in order]
 for ax,q,(_,_,title) in zip(axes,frames,panels):
  y=np.arange(len(q));ax.barh(y-.14,q.mean_bias_rmse_component,height=.23,color='#D55E00',zorder=3)
  ax.barh(y+.14,q.mean_centered_rmse_component,height=.23,color='#0072B2',zorder=3)
  ax.scatter(q.mean_delta_rmse,y,s=30,marker='D',color='#202020',zorder=4)
  ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlim(lo,hi);ax.axvline(0,color='#777777',lw=.8)
  ax.set_title(title,pad=14);ax.grid(axis='x',alpha=.2,zorder=0);ax.spines[['top','right']].set_visible(False)
  ax.xaxis.set_major_locator(MaxNLocator(nbins=4));ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
  ax.set_xlabel('Contribution to mean RMSE difference\n(g C m$^{-2}$ d$^{-1}$)',labelpad=9)
 legend=[Patch(facecolor='#D55E00',label='Squared-mean-error component'),Patch(facecolor='#0072B2',label='Centered-error component'),Line2D([0],[0],color='#202020',marker='D',linestyle='none',label='Total RMSE difference')]
 fig.legend(handles=legend,loc='lower center',ncol=3,bbox_to_anchor=(.5,.015),frameon=False)
 fig.subplots_adjust(left=.075,right=.99,top=.9,bottom=.25,wspace=.43)
 for ext in ['png','svg','pdf']:fig.savefig(out/('FigureS3_fixed_error_components.'+ext),dpi=220)
 plt.close(fig)
 receipt={'status':'pass','aggregate_input_sha256':sha(path),'plot_code_sha256':sha(__file__),'panels':3,'displayed_rows':18,'selection':'Ridge/RF prespecified primary pair and HGB existing secondary; both own-reference paths,all3declared contrasts/cohorts','all_five_families_and_four_reference_paths_in_tables':True,'new_model_fits':0,'files':{p.name:sha(p) for p in out.iterdir() if p.is_file()}}
 (out/'PLOT_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
 print(json.dumps({'status':'pass','panels':3,'rows':18}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--run',required=True);p.add_argument('--out',required=True);a=p.parse_args();main(a.run,a.out)
