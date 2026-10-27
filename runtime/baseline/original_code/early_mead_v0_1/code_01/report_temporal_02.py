"""Report all fixed outputs; no fitting, selection, new metrics or source edits."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path);args=p.parse_args();root=args.root
 run=root/'run01';plan=root/'plan01';out=args.output or root/'reports_01';assert not out.exists()
 post=json.loads((root/'independent_postfit/INDEPENDENT_POSTFIT_VERIFICATION.json').read_text())
 assert post['status']=='pass'
 protocol=json.loads((plan/'protocol.json').read_text());execution=json.loads((run/'execution.json').read_text())
 assert execution['protocol_sha256']==sha(plan/'protocol.json')==protocol.get('_protocol_sha256',execution['protocol_sha256'])
 assert post['protocol_sha256']==execution['protocol_sha256'] and Path(post['run_root']).resolve()==run.resolve()
 macro=pd.read_csv(run/'scores_macro.csv');season=pd.read_csv(run/'scores_site_season.csv');contrast=pd.read_csv(run/'input_contrasts.csv');pooled=pd.read_csv(run/'scores_pooled.csv')
 own=macro.loc[macro.fit_reference.eq(macro.score_reference)].copy()
 ownseason=season.loc[season.fit_reference.eq(season.score_reference)].copy()
 assert len(own)==20 and len(ownseason)==60 and len(contrast)==20
 out.mkdir();pd.DataFrame(protocol['test_site_season_counts']).to_csv(out/'TableS17_T1_case_counts.csv',index=False)
 for name,data in [('TableS18_T1_all_reference_macro',macro),('TableS18b_T1_pooled',pooled),('TableS19_T1_site_season',season),('TableS20_T1_C_minus_B',contrast)]:data.to_csv(out/(name+'.csv'),index=False,float_format='%.17g')
 owncontrast=contrast.loc[contrast.fit_reference.eq(contrast.score_reference)].copy()
 controls=owncontrast.loc[owncontrast.model.isin(['train_mean','seasonal_harmonic'])]
 assert np.allclose(controls.absolute_difference,0,rtol=0,atol=0)
 pairs={}
 for family in ['random_forest','ridge']:
  for ref in ['VUT_REF','CUT_REF']:
   row=owncontrast.loc[owncontrast.model.eq(family)&owncontrast.fit_reference.eq(ref)].iloc[0]
   pairs[family,ref]=row.to_dict()
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
 for suffix in ['png','svg']:fig.savefig(out/('FigureS2_T1_hindcast.'+suffix),dpi=220)
 plt.close(fig)
 RF=[pairs['random_forest',ref] for ref in ['VUT_REF','CUT_REF']];R=[pairs['ridge',ref] for ref in ['VUT_REF','CUT_REF']]
 signs={}
 for family in ['random_forest','ridge']:
  signs[family]={}
  for ref in ['VUT_REF','CUT_REF']:
   z=ownseason.loc[ownseason.model.eq(family)&ownseason.fit_reference.eq(ref)].pivot(index=['site_id','harvest_year'],columns='arm',values='rmse')
   signs[family][ref]=int((z.C_core_weighted<z.B_uniform_pixels).sum())
 result={'status':'reported_after_independent_verification','protocol_sha256':sha(plan/'protocol.json'),
  'n_train':480,'n_test':437,'test_counts':protocol['test_site_season_counts'],'model_tasks':20,'learned_estimators':16,'training_mean_controls':4,
  'own_prediction_records':execution['own_prediction_records'],'all_two_scoring_reference_records':execution['cross_prediction_records'],
  'RF_own_reference_contrasts':RF,'Ridge_own_reference_contrasts':R,'favorable_site_season_counts':signs,
  'figure_heatmap_values_percent':values.tolist(),'new_independent_region':False,'forecasting':False,
  'source_dependence':'VUT/CUT are dependent processed NT references','additional_model_fits':0}
 (out/'RESULT_SUMMARY.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 lines=['# Conditional earlier-year Mead GPP extension','',
  'This separately frozen retrospective experiment trains on the 480 existing Mead cases from Ne1 2012/2013 and Ne3 2013, then evaluates Ne1 2010/2011 and Ne3 2011. It uses the same five fixed model recipes and both dependent processed NT references, with uniform-cell B and conditional-core C optical inputs. No point-arm contrast, tuning, reference selection or new-region/forecast claim is made.','',
  'Source preparation: 92 Terra MOD09A1 items and six assets per item were acquired using exact source metadata, native geometry and the inherited QA. Shared point-valid/all-positive-core-valid scene anchors are followed by reflectance-before-NDVI and internal same-year interpolation only across brackets ≤24 days. Full May–October optical coverage is 522/552 dates. The separate fixed publication-window weather set contains 476 dates; it is not the result of sequential filtering of those 522 optical dates. Within the 476-date window, daily source/native QC gives 450 paired reference dates. Their intersection with complete B/C optical and seven weather predictors contains 437 (128/148/161). The thirteen additional exclusions are missing optical coverage: Ne1 17–26 October 2011 and Ne3 17–19 October 2011, selected before performance scores. No target values are imputed or clipped.','',
  'Twenty model tasks (sixteen learned estimators and four training-mean controls) produce 8,740 own-reference records and 17,480 records across both scoring references, including the own-reference directions. These are repeated predictions on 437 physical days, not independent field replications. All source/plan hashes remain unchanged; serialized-model predictions and metrics were independently replayed without fitting additional estimators.','',
  '| Family | Fit/scoring reference | B mean season RMSE | C mean season RMSE | C − B (%) | Seasons with lower C error |',
  '|---|---|---:|---:|---:|---:|']
 for family,label in [('random_forest','RF (primary)'),('ridge','Ridge (secondary)')]:
  for ref in ['VUT_REF','CUT_REF']:
   x=pairs[family,ref];lines.append(f"| {label} | {ref} | {x['B_RMSE']:.6f} | {x['C_RMSE']:.6f} | {x['relative_percent']:+.6f} | {signs[family][ref]}/3 |")
 lines += ['', 'Negative C − B values indicate lower RMSE for C. All five families, all site-seasons, both own/cross-reference scoring paths and pooled summaries are retained in Tables S18–S20, regardless of direction. The primary estimand equally averages the three test site-season RMSEs; these are dependent observations in one already represented region.', '',
  'The prespecified secondary HGB has lower C-arm mean site-season RMSE than RF under both own-reference paths: 3.104956 versus 3.112224 g C m⁻² d⁻¹ under VUT and 3.053804 versus 3.091772 under CUT. Its C − B differences are −2.976237% and −2.958710%, respectively. RF therefore remains the prespecified primary comparator but is not the best-performing family in this extension.', '',
  '![Figure S2](FigureS2_T1_hindcast.png)','',
  'The test years precede training years and both sites were seen in training. This adds conditional earlier-year evidence, not spatial transfer, chronological forecasting or sealed confirmation. The published normalized monthly 80% DAY cores are not full daily footprints. Ne1 physical mast registration and historical human exposure remain unresolved; fixed point-valid anchor selection makes B/C support conditional even though no A-point contrast is evaluated. Annual maize and operation dates are publication-supported, not original field logs. The DD/HR CUT quality-fraction mismatch remains a resolution-specific diagnostic; it was not used to invent an exclusion threshold. Original uStar run logs were not recovered and malformed micro-prefix text in HR metadata was retained as a source caveat; targets retain the verified DD g C unit and unchanged source values. Neither processed reference is independent GPP truth.','',
  'This does not settle Ro1 field dates, German source clock/NT/uStar/LE history, the four-crop GPP/ET study or incremental agent value. Existing principal 642-case results and every prior fit remain preserved. The source of any model-direction change must not be attributed causally to crop, management or footprints from this one-region comparison.']
 (out/'TEMPORAL_EXTENSION_EN_v0_1.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 print(json.dumps({'status':'reported','RF_percent':[x['relative_percent'] for x in RF],'Ridge_percent':[x['relative_percent'] for x in R]}))

if __name__=='__main__':main()
