from pathlib import Path
import argparse,json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main(run,out):
 run=Path(run);out=Path(out);out.mkdir(exist_ok=False)
 m=pd.read_csv(run/'TableS24_monthly_weight_concentration.csv');s=pd.read_csv(run/'TableS26_site_season_input_summary.csv')
 order=[('US-Ne1',2012),('US-Ne1',2013),('US-Ne3',2013),('US-Ro1',2015),('US-Ne1',2010),('US-Ne1',2011),('US-Ne3',2011)]
 labels=[('P' if i<4 else 'T1')+' · '+site.replace('US-','')+' '+str(year) for i,(site,year) in enumerate(order)]
 matrix=np.full((7,12),np.nan)
 for i,(site,year) in enumerate(order):
  q=m.loc[m.site_id.eq(site)&m.year.eq(year)&m.used_by_retained_case_anchor.eq(True)]
  for _,r in q.iterrows():matrix[i,int(r.month)-1]=r.concentration_cell_count
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
 fig,(a,b)=plt.subplots(1,2,figsize=(14.5,5.6),gridspec_kw={'width_ratios':[1.45,1]})
 cmap=plt.get_cmap('cividis').copy();cmap.set_bad('#eeeeee')
 im=a.imshow(np.ma.masked_invalid(matrix),vmin=1,vmax=3,cmap=cmap,aspect='auto')
 for i in range(7):
  for j in range(12):
   if np.isfinite(matrix[i,j]):a.text(j,i,f'{matrix[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if matrix[i,j]<1.7 else '#111111')
 a.set_yticks(range(7),labels);a.set_xticks(range(12),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],rotation=45,ha='right')
 a.set_title('(a) Weight concentration in used monthly maps',pad=13);a.set_xlabel('Actual composite-start month of a retained-date anchor',labelpad=12)
 a.axhline(3.5,color='white',lw=2);fig.colorbar(im,ax=a,fraction=.035,pad=.025,label='Concentration cell count, $N_c$')
 q=s.set_index(['site_id','year']).loc[order];y=np.arange(7)
 b.barh(y,q.mean_absolute_NDVI_difference,color=['#0072B2']*4+['#D55E00']*3,height=.45)
 b.scatter(q.max_absolute_NDVI_difference,y,marker='o',s=45,facecolors='none',edgecolors='#222222',zorder=4)
 b.set_yticks(y,labels);b.invert_yaxis();b.axhline(3.5,color='#aaaaaa',lw=.8);b.grid(axis='x',alpha=.2);b.set_axisbelow(True)
 b.spines[['top','right']].set_visible(False);b.set_xlim(-.002,max(q.max_absolute_NDVI_difference)*1.13)
 b.set_xlabel('Absolute NDVI difference between C and B',labelpad=12);b.set_title('(b) Realized input differences on retained dates',pad=13)
 b.legend(handles=[Line2D([0],[0],color='#555555',lw=7,label='Mean absolute difference'),Line2D([0],[0],marker='o',color='#222222',markerfacecolor='none',linestyle='none',label='Maximum absolute difference')],loc='upper center',bbox_to_anchor=(.5,-.30),frameon=False,fontsize=8.5)
 fig.subplots_adjust(left=.10,right=.985,top=.88,bottom=.32,wspace=.48)
 fig.text(.10,.045,'P: principal 642 dates; T1: separate 437-date hindcast. Blank, light-grey months are unused by retained-date anchors.',fontsize=9)
 for ext in ['png','svg']:fig.savefig(out/('FigureS4_conditional_input_support.'+ext),dpi=220)
 plt.close(fig)
 (out/'PLOT_RECEIPT.json').write_text(json.dumps({'status':'pass','panels':2,'site_seasons':7,'used_map_cells_shown':int(np.isfinite(matrix).sum()),'input_hashes':{p.name:sha(p) for p in [run/'TableS24_monthly_weight_concentration.csv',run/'TableS26_site_season_input_summary.csv']},'code_sha256':sha(__file__),'files':{p.name:sha(p) for p in out.iterdir() if p.is_file()},'new_predictions':0},indent=2)+'\n',encoding='utf-8')
 print(json.dumps({'status':'pass','panels':2,'site_seasons':7}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--run',required=True);p.add_argument('--out',required=True);a=p.parse_args();main(a.run,a.out)
