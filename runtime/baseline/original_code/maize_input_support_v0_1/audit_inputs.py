"""Target-free audit of existing map weights, optical inputs and interpolation."""
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
import json,hashlib,argparse,platform
import numpy as np
import pandas as pd
FEATURES=['NDVI_modis','EA_kPa_agera5','WS_ms_agera5','P_mm_agera5','TA_max_C_agera5','TA_min_C_agera5','SW_IN_MJ_m2_d_agera5','VPD_kPa_agera5']
IDS=['site_id','date','record_id','harvest_year']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def normalize(d):
 d=d.copy();d['date']=pd.to_datetime(d.date).dt.strftime('%Y-%m-%d');return d
def weight_metrics(weights):
 p=np.asarray(weights,dtype=float);assert len(p) and np.isfinite(p).all() and (p>0).all()
 total=float(p.sum());assert abs(total-1)<1e-7;p=p/total;n=len(p)
 nc=float(1/(p@p));tv=float(np.abs(p-1/n).sum()/2)
 assert 1-1e-12<=nc<=n+1e-12 and -1e-12<=tv<=1
 return dict(n_cells=n,source_conditional_weight_sum=total,concentration_cell_count=nc,max_conditional_weight=float(p.max()),distance_from_uniform=tv)
def boundary_checks():
 a=weight_metrics([1]);assert a['n_cells']==a['concentration_cell_count']==1 and a['distance_from_uniform']==0
 a=weight_metrics([.5,.5]);assert a['concentration_cell_count']==2 and a['distance_from_uniform']==0
 a=weight_metrics([.75,.25]);assert np.isclose(a['concentration_cell_count'],1.6) and a['distance_from_uniform']==.25
 return 3
def run(protocol_path,expected,out):
 assert sha(protocol_path)==expected;plan=read(protocol_path);assert sha(__file__)==plan['code_sha256']
 out=Path(out);out.mkdir(exist_ok=False);tests=boundary_checks();sources=plan['files']
 for s in sources.values():assert sha(s['path'])==s['sha256']
 for s in plan['manifest_pins']+plan['code_pins']:assert sha(s['path'])==s['sha256']
 for s in plan['map_pins']:assert sha(s['path'])==s['sha256']
 def file(key):return sources[key]['path']
 monthly=[];month_index={}
 for cohort,key in [('principal_642','primary_monthly_weights'),('T1_test437','T1_monthly_weights_in_protocol')]:
  for site in read(file(key))['sites']:
   for m in site['maps']:
    assert sha(m['path'])==m['sha256']
    pixels=m['pixels'];assert len({(p['row'],p['col']) for p in pixels})==len(pixels)
    row=dict(cohort=cohort,site_id=site['site_id'],year=int(site['year']),month=int(m['month']),map_path=m['path'],map_sha256=m['sha256'])
    row.update(weight_metrics([p['conditional_core_weight'] for p in pixels]));monthly.append(row);month_index[(row['site_id'],row['year'],row['month'])]=row
 assert len(monthly)==84 and len(month_index)==84
 # Read predictors/keys only from snapshots that also contain target columns.
 frames={}
 for cohort,prefix in [('principal_642','primary642'),('T1_test437','T1')]:
  for arm in ['B','C']:
   for ref in ['VUT','CUT']:
    key=prefix+'_'+arm+'_'+ref;d=normalize(pd.read_parquet(file(key),columns=IDS+FEATURES))
    assert len(d)==(642 if cohort=='principal_642' else 917)
    assert not d.duplicated(['site_id','date']).any() and np.isfinite(d[FEATURES].to_numpy()).all()
    frames[(cohort,arm,ref)]=d
  base=frames[(cohort,'B','VUT')]
  for arm in ['B','C']:
   for ref in ['VUT','CUT']:
    d=frames[(cohort,arm,ref)];pd.testing.assert_frame_equal(d[IDS],base[IDS],check_dtype=False)
    assert np.array_equal(d[FEATURES[1:]].to_numpy(),base[FEATURES[1:]].to_numpy())
   assert np.array_equal(frames[(cohort,arm,'VUT')][FEATURES].to_numpy(),frames[(cohort,arm,'CUT')][FEATURES].to_numpy())
 fold=read(file('T1_frozen_fold'));assert fold['n_train']==480 and fold['n_test']==437
 train_indices=fold['train_indices'];test_indices=fold['test_indices']
 assert not set(train_indices)&set(test_indices) and set(train_indices)|set(test_indices)==set(range(917))
 for arm in ['B','C']:
  train=frames[('T1_test437',arm,'VUT')].iloc[train_indices].sort_values(['site_id','date']).reset_index(drop=True)
  old=frames[('principal_642',arm,'VUT')];old=old.loc[old.site_id.isin(['US-Ne1','US-Ne3'])].sort_values(['site_id','date']).reset_index(drop=True)
  assert len(train)==len(old)==480
  pd.testing.assert_frame_equal(train[['site_id','date','harvest_year']+FEATURES],old[['site_id','date','harvest_year']+FEATURES],check_exact=True,check_dtype=False)
 scenes={
  'principal_642':pd.read_parquet(file('primary_original_accepted_scene_audit')),
  'T1_test437':pd.DataFrame(read(file('T1_accepted_scene_audit'))['scenes'])}
 for cohort,a in scenes.items():
  a['anchor']=pd.to_datetime(a.anchor).dt.strftime('%Y-%m-%d')
  assert not a.duplicated(['site_id','year','anchor']).any()
  for r in a.to_dict('records'):
   assert int(r['positive_core_pixels'])==month_index[(r['site_id'],int(r['year']),pd.Timestamp(r['anchor']).month)]['n_cells']
 daily=[];anchor_use=Counter();map_use=Counter();largest_replay_difference=0.
 for cohort in ['principal_642','T1_test437']:
  b=frames[(cohort,'B','VUT')];c=frames[(cohort,'C','VUT')]
  if cohort=='T1_test437':b=b.iloc[test_indices];c=c.iloc[test_indices]
  b=b.set_index(['site_id','date']).sort_index();c=c.set_index(['site_id','date']).sort_index();assert b.index.equals(c.index)
  anchor=scenes[cohort].loc[scenes[cohort].shared_anchor_usable.eq(True)].copy()
  assert np.isfinite(anchor[['B_uniform_pixels','C_core_weighted']].to_numpy()).all()
  for (site,date),r in b.iterrows():
   year=int(r['harvest_year']);d=pd.Timestamp(date);a=anchor.loc[anchor.site_id.eq(site)&anchor.year.eq(year)].sort_values('anchor').reset_index(drop=True)
   dates=pd.DatetimeIndex(a.anchor);j=dates.searchsorted(d,side='right')-1
   assert j>=0,(site,date,'no left anchor')
   left=a.iloc[j];ld=dates[j]
   if ld==d:right=left;rd=ld;fraction=0.
   else:
    assert j+1<len(a),(site,date,'no right anchor')
    right=a.iloc[j+1];rd=dates[j+1];assert (rd-ld).days<=24
    fraction=(d-ld).days/(rd-ld).days
   assert ld.year==rd.year==year and ld<=d<=rd
   bv=float(r['NDVI_modis']);cv=float(c.loc[(site,date),'NDVI_modis'])
   for arm,actual in [('B_uniform_pixels',bv),('C_core_weighted',cv)]:
    replay=(1-fraction)*float(left[arm])+fraction*float(right[arm]);difference=abs(replay-actual)
    assert difference<1e-12,(cohort,site,date,arm,actual,replay)
    largest_replay_difference=max(largest_replay_difference,difference)
   role='T1_test' if cohort=='T1_test437' else ('T1_training_reuse' if site!='US-Ro1' else 'not_used_in_T1')
   daily.append(dict(cohort=cohort,site_id=site,year=year,date=date,T1_role=role,B_NDVI=bv,C_NDVI=cv,C_minus_B_NDVI=cv-bv,absolute_NDVI_difference=abs(cv-bv),left_anchor=ld.strftime('%Y-%m-%d'),right_anchor=rd.strftime('%Y-%m-%d'),left_item_id=str(left['item_id']),right_item_id=str(right['item_id']),left_map_month=ld.month,right_map_month=rd.month,bracket_days=(rd-ld).days,right_anchor_fraction=fraction,direct_anchor=(ld==rd),right_anchor_after_case=(rd>d),cross_month_bracket=(ld.month!=rd.month)))
   for s in {ld.strftime('%Y-%m-%d'),rd.strftime('%Y-%m-%d')}:anchor_use[(cohort,site,year,s)]+=1
   for m in {ld.month,rd.month}:map_use[(site,year,m)]+=1
 d=pd.DataFrame(daily);assert len(d)==1079 and not d.duplicated(['site_id','year','date']).any()
 expected_cases=plan['case_counts']
 assert {f'{s}:{int(y)}':int(n) for (s,y),n in d.groupby(['site_id','year']).size().items()}==expected_cases
 assert (d.T1_role.eq('T1_training_reuse')).sum()==480
 # Check the original saved daily curves independently of model snapshots.
 for cohort,key in [('principal_642','primary_original_daily_optics'),('T1_test437','T1_daily_full_year_optics')]:
  optical=normalize(pd.read_parquet(file(key))).set_index(['site_id','date'])
  selected=d.loc[d.cohort.eq(cohort)].set_index(['site_id','date'])
  for outfield,srcfield in [('B_NDVI','B_uniform_pixels'),('C_NDVI','C_core_weighted')]:
   assert np.max(np.abs(selected[outfield].to_numpy()-optical.loc[selected.index,srcfield].to_numpy()))<1e-12
 for r in monthly:
  r['retained_dates_using_map']=map_use[(r['site_id'],r['year'],r['month'])]
  r['used_by_retained_case_anchor']=bool(r['retained_dates_using_map'])
 usage=[]
 for cohort,a in scenes.items():
  for r in a.to_dict('records'):
   key=(cohort,r['site_id'],int(r['year']),r['anchor']);month=pd.Timestamp(r['anchor']).month
   usage.append(dict(cohort=cohort,site_id=r['site_id'],year=int(r['year']),anchor=r['anchor'],item_id=r['item_id'],map_month=month,shared_anchor_usable=bool(r['shared_anchor_usable']),retained_dates_using_anchor=anchor_use[key],B_anchor_NDVI=r['B_uniform_pixels'],C_anchor_NDVI=r['C_core_weighted']))
 summaries=[]
 for (cohort,site,year,role),g in d.groupby(['cohort','site_id','year','T1_role'],sort=True):
  active=[r for r in monthly if r['site_id']==site and r['year']==year and r['used_by_retained_case_anchor']]
  summaries.append(dict(cohort=cohort,site_id=site,year=int(year),T1_role=role,n_days=len(g),mean_C_minus_B_NDVI=float(g.C_minus_B_NDVI.mean()),mean_absolute_NDVI_difference=float(g.absolute_NDVI_difference.mean()),median_absolute_NDVI_difference=float(g.absolute_NDVI_difference.median()),max_absolute_NDVI_difference=float(g.absolute_NDVI_difference.max()),n_exact_equal_BC=int((g.C_minus_B_NDVI==0).sum()),n_difference_above_tolerance=int((g.absolute_NDVI_difference>1e-12).sum()),active_map_months=len(active),minimum_active_cell_count=min(r['n_cells'] for r in active),maximum_active_cell_count=max(r['n_cells'] for r in active),minimum_active_concentration=min(r['concentration_cell_count'] for r in active),maximum_active_concentration=max(r['concentration_cell_count'] for r in active),direct_anchor_days=int(g.direct_anchor.sum()),right_anchor_after_case_days=int(g.right_anchor_after_case.sum()),cross_month_bracket_days=int(g.cross_month_bracket.sum()),maximum_bracket_days=int(g.bracket_days.max())))
 for name,frame in [('TableS24_monthly_weight_concentration.csv',pd.DataFrame(monthly)),('TableS25_retained_daily_input_lineage.csv',d),('TableS26_site_season_input_summary.csv',pd.DataFrame(summaries)),('TableS27_scene_anchor_usage.csv',pd.DataFrame(usage))]:
  frame.to_csv(out/name,index=False,float_format='%.17g')
 assert all(sha(s['path'])==s['sha256'] for s in sources.values())
 save(out/'EXECUTION_RECEIPT.json',{'status':'pass','created_utc':datetime.now(timezone.utc).isoformat(),'protocol_sha256':expected,'code_sha256':sha(__file__),'runtime':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__},'boundary_checks':tests,'monthly_map_rows':len(monthly),'retained_daily_rows':len(d),'site_season_summary_rows':len(summaries),'scene_audit_rows':len(usage),'active_monthly_maps':sum(r['used_by_retained_case_anchor'] for r in monthly),'maximum_NDVI_replay_difference':largest_replay_difference,'T1_training_reuse_days':480,'target_columns_read':[],'new_predictions':0,'new_model_fits':0,'new_FFP_runs':0,'new_eligibility_approvals':0,'sources_unchanged':True,'Cloud_calls':0,'scope':'Post-review target-free construction audit; concentration is a property of conditional raster-cell weights, not field/crop purity, independent sample size or full daily flux support.'})
 print(json.dumps({'status':'pass','maps':len(monthly),'daily':len(d),'summaries':len(summaries),'scene_rows':len(usage),'max_replay_difference':largest_replay_difference}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--protocol',required=True);p.add_argument('--expected-sha',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(a.protocol,a.expected_sha,a.out)
