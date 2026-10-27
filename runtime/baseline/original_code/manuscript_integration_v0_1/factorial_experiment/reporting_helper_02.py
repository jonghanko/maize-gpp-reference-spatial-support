"""Postfit aggregate reporting only; never mutate frozen plan/code/run files."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ARMS=['A_existing','A_shared_scene','B_uniform_pixels','C_core_weighted']
VARS=['VUT_REF','CUT_REF'];MODELS=['train_mean','seasonal_harmonic','ridge','random_forest','hist_gradient_boosting']
CONTRASTS=[('C_minus_B','C_core_weighted','B_uniform_pixels'),('C_minus_A_shared','C_core_weighted','A_shared_scene'),('A_existing_minus_A_shared','A_existing','A_shared_scene')]
NAMES={'train_mean':'Mean','seasonal_harmonic':'Seasonal','ridge':'Ridge','random_forest':'RF','hist_gradient_boosting':'HGB'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 def clean(x):
  if isinstance(x,dict):return {k:clean(z) for k,z in x.items()}
  if isinstance(x,list):return [clean(z) for z in x]
  if isinstance(x,float) and not np.isfinite(x):return None
  if isinstance(x,np.generic):return clean(x.item())
  return x
 with Path(p).open('x') as f:json.dump(clean(v),f,indent=2,allow_nan=False)
def main(root):
 freeze=json.loads((root/'plan01/freeze.json').read_text());manifest=json.loads((root/'runs/run_manifest.json').read_text())['files']
 def check():
  for group in ['inputs','code','plan_files']:assert all(sha(p)==h for p,h in freeze[group].items())
  assert all(sha(root/'runs'/p)==d['sha256'] for p,d in manifest.items())
 check();out=root/'reports_02';out.mkdir(exist_ok=False);s=root/'runs/summary'
 c=pd.read_csv(s/'input_contrasts.csv');i=pd.read_csv(s/'reference_input_interactions.csv');m=pd.read_csv(s/'cross_scores_macro_summary.csv');ss=pd.read_csv(s/'cross_scores_site_season.csv');ranks=pd.read_csv(s/'family_rankings.csv')
 own=c.loc[c.fit_reference.eq(c.score_reference)]
 mainc=own.loc[own.axis.eq('year_loyo')&own.metric.eq('macro_site_season_rmse')]
 mainc.to_csv(out/'all_family_LOYO_input_effects.csv',index=False,float_format='%.17g')
 c.loc[c.primary].to_csv(out/'primary_input_effects.csv',index=False,float_format='%.17g')
 maini=i.loc[i.axis.eq('year_loyo')&i.metric.eq('macro_site_season_rmse')&i.interaction_type.eq('CUT_own_minus_VUT_own')]
 maini.to_csv(out/'all_family_LOYO_reference_interactions.csv',index=False,float_format='%.17g')
 ownm=m.loc[m.fit_reference.eq(m.score_reference)]
 ownm.to_csv(out/'all_own_reference_macro_scores.csv',index=False,float_format='%.17g')
 partial=own.loc[own.axis.eq('cluster_loso')&own.metric.eq('macro_site_season_rmse')]
 partial.to_csv(out/'partial_Ro1_input_effects.csv',index=False,float_format='%.17g')
 seasons=[]
 q=ss.loc[ss.axis.eq('year_loyo')&ss.fit_reference.eq(ss.score_reference)]
 for key,g in q.groupby(['model','fit_reference','site_id','spatial_cluster','harvest_year'],sort=True):
  for contrast,first,second in CONTRASTS:
   a=float(g.loc[g.arm.eq(first),'rmse'].iloc[0]);b=float(g.loc[g.arm.eq(second),'rmse'].iloc[0])
   seasons.append({'model':key[0],'variant':key[1],'site_id':key[2],'region':key[3],'year':int(key[4]),'n':int(g.n.iloc[0]),'contrast':contrast,'first_RMSE':a,'second_RMSE':b,'absolute_RMSE_difference':a-b,'relative_percent':100*(a-b)/b})
 seasonal=pd.DataFrame(seasons);seasonal.to_csv(out/'all_LOYO_site_season_input_effects.csv',index=False,float_format='%.17g')
 rankdiff=[]
 for (axis,arm,metric),g in ranks.loc[ranks.fit_reference.eq(ranks.score_reference)].groupby(['axis','arm','metric'],sort=True):
  p=g.pivot(index='model',columns='fit_reference',values='rank')
  rankdiff.append({'axis':axis,'arm':arm,'metric':metric,'all_family_ranks_same_between_own_references':bool(p.VUT_REF.eq(p.CUT_REF).all()),'changed_models':p.index[~p.VUT_REF.eq(p.CUT_REF)].tolist()})
 save(out/'RESULT_SUMMARY.json',{'cohort':freeze['summary'],'execution':json.loads((root/'runs/execution.json').read_text()),'guard_reason_counts':json.loads((root/'audit_metadata/GUARD_REASON_WORDING_ADDENDUM_01.json').read_text()),'primary_input_effects':c.loc[c.primary].to_dict('records'),'primary_interactions':i.loc[i.primary_descriptive].to_dict('records'),'all_five_family_LOYO_effects':mainc.to_dict('records'),'partial_Ro1_effects':partial.to_dict('records'),'rank_reference_comparisons':rankdiff,'scientific_primary_approved':False,'ET_fit':False,'all_frozen_files_unchanged':True})
 plt.rcParams.update({'font.size':10,'axes.titlesize':11,'axes.labelsize':10})
 colors=['#2369a0','#db8733'];models=['ridge','random_forest'];ct=['C_minus_B','C_minus_A_shared'];labs=['C − B','C − A_shared']
 # Both primary models on one absolute-error-difference scale.
 fig,axs=plt.subplots(1,2,figsize=(10.5,4.9));p=c.loc[c.primary];limit=max(abs(p.absolute_RMSE_difference))*1.35
 for ax,model in zip(axs,models):
  for k,var in enumerate(VARS):
   z=p.loc[p.model.eq(model)&p.fit_reference.eq(var)].set_index('contrast').reindex(ct)
   bars=ax.bar(np.arange(2)+(k-.5)*.34,z.absolute_RMSE_difference,.34,label=var[:3],color=colors[k])
   ax.bar_label(bars,labels=[f'{v:+.4f}' for v in z.absolute_RMSE_difference],padding=3,fontsize=9)
  ax.set_xticks([0,1],labs);ax.set_ylim(-limit,limit);ax.axhline(0,color='.3',lw=.8);ax.set_title(NAMES[model]);ax.set_ylabel('RMSE difference (g C m$^{-2}$ d$^{-1}$)');ax.legend(frameon=False)
 fig.suptitle('Complete year LOYO: own-reference input contrasts, same 642 cases',fontsize=13)
 fig.text(.5,.03,'Negative values indicate lower C error. Four site-seasons, two regions; descriptive estimates without uncertainty bars.',ha='center',fontsize=9)
 fig.tight_layout(rect=(0,.08,1,.93));fig.savefig(out/'primary_LOYO_input_effects.png',dpi=220);fig.savefig(out/'primary_LOYO_input_effects.svg');plt.close(fig)
 # All five families retained in the dependent own-reference interaction.
 fig,axs=plt.subplots(1,2,figsize=(11.7,4.9));limit=max(abs(maini.interaction_RMSE_units))*1.4
 limit=max(limit,.001)
 for ax,contrast,label in zip(axs,ct,labs):
  z=maini.loc[maini.contrast.eq(contrast)].set_index('model').reindex(MODELS)
  bars=ax.bar(np.arange(5),z.interaction_RMSE_units,color=['#bfc4c8','#9da6ad','#647b8d','#2369a0','#80bdd8'])
  ax.bar_label(bars,labels=[f'{v:+.4f}' for v in z.interaction_RMSE_units],padding=3,fontsize=8)
  ax.set_xticks(np.arange(5),[NAMES[k] for k in MODELS]);ax.set_ylim(-limit,limit);ax.axhline(0,color='.3',lw=.8);ax.set_title(label);ax.set_ylabel('CUT contrast − VUT contrast\n(g C m$^{-2}$ d$^{-1}$)')
 fig.suptitle('Reference × input interaction: own-fit/own-score LOYO contrasts',fontsize=13)
 fig.text(.5,.03,'Both training and scoring path change here; full fixed-fit/fixed-score matrices are retained. No causal or population inference.',ha='center',fontsize=9)
 fig.tight_layout(rect=(0,.08,1,.93));fig.savefig(out/'reference_input_interactions.png',dpi=220);fig.savefig(out/'reference_input_interactions.svg');plt.close(fig)
 # Every season is shown, including the unchanged-X Ro1 test; no selected counterexample.
 order=[('US-Ne1',2012),('US-Ne1',2013),('US-Ne3',2013),('US-Ro1',2015)]
 seasonlabs=['Ne1 2012\nn=156','Ne1 2013\nn=167','Ne3 2013\nn=157','Ro1 2015\nn=162']
 p=seasonal.loc[seasonal.model.isin(models)&seasonal.contrast.isin(ct)];limit=max(abs(p.absolute_RMSE_difference))*1.35
 fig,axs=plt.subplots(2,2,figsize=(12.4,7.9));fig.subplots_adjust(left=.09,right=.98,bottom=.16,top=.87,wspace=.28,hspace=.5)
 for row,contrast in enumerate(ct):
  for col,model in enumerate(models):
   ax=axs[row,col];ax.axvspan(2.55,3.45,color='#ececec',zorder=0)
   for k,var in enumerate(VARS):
    z=p.loc[p.model.eq(model)&p.contrast.eq(contrast)&p.variant.eq(var)].set_index(['site_id','year']).reindex(order)
    bars=ax.bar(np.arange(4)+(k-.5)*.32,z.absolute_RMSE_difference,.32,label=var[:3],color=colors[k])
    ax.bar_label(bars,labels=[f'{v:+.3f}' for v in z.absolute_RMSE_difference],padding=2,fontsize=7)
   ax.set_ylim(-limit,limit);ax.axhline(0,color='.3',lw=.8);ax.set_xticks(np.arange(4),seasonlabs);ax.set_title(NAMES[model]+' · '+labs[row]);ax.set_ylabel('RMSE difference (g C m$^{-2}$ d$^{-1}$)')
   if row==0:ax.legend(frameon=False,fontsize=8,ncol=2)
 fig.suptitle('All four LOYO site-seasons: heterogeneous input effects',fontsize=13,y=.97)
 fig.text(.5,.045,'Shaded Ro1 has identical test inputs in all four arms: changes arise through Mead training representation.\nRo1 LOYO2015, partial region LOSO and joint Ro1/2015 have identical train/test membership.',ha='center',fontsize=9)
 fig.savefig(out/'all_LOYO_site_season_effects.png',dpi=220);fig.savefig(out/'all_LOYO_site_season_effects.svg');plt.close(fig)
 check();print(json.dumps({'status':'aggregate_reporting_complete','frozen_files_unchanged':True}))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args().root)
