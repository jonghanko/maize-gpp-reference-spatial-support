from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
def main(data,out):
 data=Path(data);out=Path(out);out.mkdir(exist_ok=False)
 macro=pd.read_csv(data/"TableS17a_T1_all_reference_macro.csv");season=pd.read_csv(data/"TableS18_T1_site_season.csv")
 own=macro.loc[macro.fit_reference.eq(macro.score_reference)].copy();ownseason=season.loc[season.fit_reference.eq(season.score_reference)].copy()
 assert len(own)==20 and len(ownseason)==60
 rows=[]
 seasons=[('US-Ne1',2010),('US-Ne1',2011),('US-Ne3',2011)]
 for site,year in seasons:
  row=[]
  for family,ref in [('random_forest','VUT_REF'),('random_forest','CUT_REF'),('ridge','VUT_REF'),('ridge','CUT_REF')]:
   z=ownseason.loc[ownseason.site_id.eq(site)&ownseason.harvest_year.eq(year)&ownseason.model.eq(family)&ownseason.fit_reference.eq(ref)]
   b=float(z.loc[z.arm.eq('B_uniform_pixels'),'rmse'].iloc[0]);c=float(z.loc[z.arm.eq('C_core_weighted'),'rmse'].iloc[0]);row.append(100*(c-b)/b)
  rows.append(row)
 values=np.array(rows)
 fig=plt.figure(figsize=(9.4,7.0),layout='constrained');grid=fig.add_gridspec(2,2,height_ratios=[1,1.05])
 for col,(family,label) in enumerate([('random_forest','Random forest (primary)'),('ridge','Ridge (secondary)')]):
  ax=fig.add_subplot(grid[0,col]);x=np.arange(2)
  for offset,arm,color,title in [(-.18,'B_uniform_pixels','#687d8e','Uniform cells (B)'),(.18,'C_core_weighted','#b55b35','Core weights (C)')]:
   y=[float(own.loc[own.model.eq(family)&own.fit_reference.eq(ref)&own.arm.eq(arm),'macro_site_season_rmse'].iloc[0]) for ref in ['VUT_REF','CUT_REF']]
   bars=ax.bar(x+offset,y,width=.34,color=color,label=title)
   ax.bar_label(bars,fmt='%.3f',padding=3,fontsize=8)
  ax.set_xticks(x,['VUT fit/score','CUT fit/score']);ax.set_title(label,fontsize=11);ax.set_ylabel('Mean site-season RMSE (g C m⁻² d⁻¹)',fontsize=9);ax.set_ylim(0,ax.get_ylim()[1]*1.40);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
  if col==0:ax.legend(fontsize=8,loc='upper left')
 ax=fig.add_subplot(grid[1,:]);limit=max(float(abs(values).max()),0.05)
 im=ax.imshow(values,cmap='RdBu_r',vmin=-limit,vmax=limit,aspect='auto')
 ax.set_xticks(range(4),['RF · VUT','RF · CUT','Ridge · VUT','Ridge · CUT']);ax.set_yticks(range(3),['Ne1 2010 (n=128)','Ne1 2011 (n=148)','Ne3 2011 (n=161)'])
 for i in range(3):
  for j in range(4):ax.text(j,i,f'{values[i,j]:+.2f}%',ha='center',va='center',color='white' if abs(values[i,j])>limit*.62 else '#141414',fontsize=10)
 ax.set_title('Site-season C − B: negative values indicate lower RMSE',fontsize=10)
 fig.colorbar(im,ax=ax,label='RMSE difference (%)',shrink=.86)
 fig.suptitle('Earlier-year Mead hindcast: train 2012/2013 → test 2010/2011\n480 training days; 437 test days; one existing region',fontsize=12)
 fig.supxlabel('Dependent processed GPP references · fixed monthly 80% DAY cores · retrospective, not a regional/forecast test',fontsize=8)
 for suffix in ['png','svg','pdf']:fig.savefig(out/('FigureS2_T1_hindcast.'+suffix),dpi=220)
 plt.close(fig)
