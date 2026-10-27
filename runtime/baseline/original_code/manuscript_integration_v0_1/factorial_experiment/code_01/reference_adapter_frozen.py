"""Paired NT VUT/CUT diagnostic; source evidence and exact cases precede fits.

The numerical engine is a byte-identical copy of the fixed rice baseline.
This adapter owns selection, all supported/unsupported strata, and cross-scoring.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
import hashlib,json,time
from pathlib import Path
import numpy as np
import pandas as pd
import engine

VARIANTS=['VUT_REF','CUT_REF']
ENGINE_SHA='c6dab8c357fd0432dbf8b06fa263d04e19b4527ad41984982862ad69409b2e9d'
BASE_PROTOCOL_SHA='2ea92e66fe137de3579d1b8ff14294d93e8972fa79439d604d829afa0a322d75'
CLUSTERS={'US-Ne1':'Mead_NE','US-Ne3':'Mead_NE','US-Ro1':'Rosemount_MN','US-IAC':'Ames_IA','US-UiC':'EnergyFarm_IL'}
TARGET='gpp_gC_m2_d'

def now():return datetime.now(timezone.utc).isoformat()
def write_json(path,value):
 with Path(path).open('x',encoding='utf-8') as f:json.dump(value,f,indent=2,allow_nan=False,default=lambda z:z.item() if isinstance(z,np.generic) else str(z));f.write('\n')
def strict_flag(d,c,v=True):
 if c not in d or not pd.api.types.is_bool_dtype(d[c].dtype) or d[c].isna().any():raise ValueError('Invalid boolean gate: '+c)
 return d[c].eq(v)
def normalize(d):
 d=d.copy();d['date']=pd.to_datetime(d.date).dt.strftime('%Y-%m-%d')
 if d.duplicated(['site_id','date']).any() or d.record_id.duplicated().any():raise ValueError('Duplicate identity')
 return d.sort_values(['site_id','date']).reset_index(drop=True)
def same(a,b,cols):
 pd.testing.assert_frame_equal(a[cols].reset_index(drop=True),b[cols].reset_index(drop=True),check_exact=True,check_dtype=True)

def construct_pair(sources,guards):
 """No ET eligibility condition, score-dependent selector or spatial relocation."""
 source={k:normalize(sources[k]) for k in VARIANTS};guard={k:normalize(guards[k]) for k in VARIANTS}
 key=['site_id','date'];same(source['VUT_REF'],source['CUT_REF'],key)
 for var in VARIANTS:
  d,g=source[var],guard[var]
  same(d,g,key)
  if not d.record_id.equals(g.record_id):raise ValueError('Guard source record identity mismatch')
  if not d.gpp_reference_variant.eq(var).all() or not g.variant.eq(var).all():raise ValueError('Reference identity mismatch')
  if not d.gpp_reference_method.eq('NT').all():raise ValueError('NT reference changed')
  if not strict_flag(d,'training_release_approved',False).all() or not strict_flag(d,'optical_representativeness_verified',False).all():raise ValueError('Original scientific approval changed')
  if not strict_flag(g,'strict_fullflux_approved',False).all():raise ValueError('Guard promoted strict approval')
  if not d.spatial_split_group.eq(d.site_id.map(CLUSTERS)).all():raise ValueError('Dependency group changed')
  dates=pd.to_datetime(d.date)
  if not dates.dt.year.eq(d.calendar_year).all():raise ValueError('Source year mismatch')
  if not g.calendar_year.equals(d.calendar_year) or not g.spatial_split_group.equals(d.spatial_split_group):raise ValueError('Guard split metadata changed')
  a=dates<pd.to_datetime(g.sowing_date_guard,errors='coerce');b=dates>pd.to_datetime(g.harvest_date_guard,errors='coerce')
  outside=a|b
  if not a.equals(strict_flag(g,'known_pre_sowing')) or not b.equals(strict_flag(g,'known_post_harvest')):raise ValueError('Documented bounds do not reproduce flags')
  if not outside.equals(strict_flag(g,'known_outside_documented_bounds')) or not (~outside).equals(strict_flag(g,'approved_for_new_calendar_proxy')):raise ValueError('Guard acceptance changed')
  # A malformed non-null date cannot quietly be treated as an unknown bound.
  for c in ['sowing_date_guard','harvest_date_guard']:
   if (g[c].notna()&pd.to_datetime(g[c],errors='coerce').isna()).any():raise ValueError('Malformed source date')
  native=strict_flag(d,'crop_year_maize_confirmed')&strict_flag(d,'calendar_may_oct_screen')&strict_flag(d,'input_alignment_complete')&strict_flag(d,'qc_gpp_native_audited_candidate')
  if not native.equals(strict_flag(d,'analysis_calendar_gpp')):raise ValueError('Original registered mask changed')
  if not d.gpp_gC_m2_d.equals(d.gpp_nt_gC_m2_d):raise ValueError('NT target identity changed')
 # Equality includes missing entries outside the eligible cohort.
 same(source['VUT_REF'],source['CUT_REF'],key+engine.FEATURES+['et_mm_d','et_qc','qc_et_native_audited_candidate','analysis_calendar_et','crop_year_maize_confirmed','calendar_may_oct_screen','input_alignment_complete'])
 same(guard['VUT_REF'],guard['CUT_REF'],key+['sowing_date_guard','harvest_date_guard','calendar_evidence_status','calendar_sowing_conflict','known_outside_documented_bounds','approved_for_new_calendar_proxy'])
 eligible=strict_flag(source['VUT_REF'],'analysis_calendar_gpp')&strict_flag(source['CUT_REF'],'analysis_calendar_gpp')&strict_flag(guard['VUT_REF'],'approved_for_new_calendar_proxy')&strict_flag(guard['CUT_REF'],'approved_for_new_calendar_proxy')
 result={}
 guard_cols=['sowing_date_guard','harvest_date_guard','calendar_evidence_status','calendar_sowing_conflict','known_outside_documented_bounds','approved_for_new_calendar_proxy','strict_fullflux_approved']
 for var,d in source.items():
  x=d.loc[eligible].copy().reset_index(drop=True);g=guard[var].loc[eligible].reset_index(drop=True)
  x['source_record_id']=x.record_id
  x['observation_key']=x.site_id+'|'+x.date
  x['record_id']='paired_reference|'+x.observation_key
  x['spatial_cluster']=x.spatial_split_group;x['harvest_year']=x.calendar_year.astype(int)
  x['season_split_group']=x.site_id+'|'+x.calendar_year.astype(str)
  x['system_id']='us_maize_processed_flux'
  x['benchmark_permission']='paired_processing_path_calendar_proxy_diagnostic'
  x['scientific_primary_approved']=False;x['actual_crop_window_applied']=False
  x['documented_bound_quarantine_applied']=True;x['all_actual_crop_windows_verified']=False
  for c in guard_cols:x[c]=g[c]
  if not len(x) or not np.isfinite(x[[TARGET,*engine.FEATURES]].to_numpy(float)).all():raise ValueError('Empty/incomplete paired cohort')
  if not pd.to_datetime(x.date).dt.month.between(5,10).all():raise ValueError('Calendar screen changed')
  result[var]=x
 same(result['VUT_REF'],result['CUT_REF'],['record_id',*key,*engine.FEATURES])
 return result

def load_inputs(config):
 if config['expected_source_rows']!=12422 or config['min_train_rows']!=300 or config['min_prior_years']!=3:raise ValueError('Scope/support threshold changed')
 if config['target']!='GPP' or config['variants']!=VARIANTS:raise ValueError('Target/reference roster changed')
 paths={}
 for kind in ['sources','guards']:
  for var,info in config[kind].items():
   p=Path(info['path']);paths[str(p)]=info['sha256']
 for p,h in config.get('frozen_evidence',{}).items():paths[p]=h
 for p,h in paths.items():
  if engine.sha(p)!=h:raise ValueError('Source/evidence hash changed: '+p)
 sources={v:pd.read_parquet(config['sources'][v]['path']) for v in VARIANTS}
 guards={v:pd.read_parquet(config['guards'][v]['path']) for v in VARIANTS}
 if any(len(d)!=12422 for d in [*sources.values(),*guards.values()]):raise ValueError('Unexpected source/guard row count')
 result=construct_pair(sources,guards)
 if any(len(d)!=config['expected_common_rows'] for d in result.values()):raise ValueError('Expected common case count changed')
 return result,paths

def plan_splits(d,min_train=300,min_prior=3):
 """Register every stratum; unsupported strata never generate a fit task."""
 clusters=sorted(d.spatial_cluster.unique());years=sorted(map(int,d.harvest_year.unique()))
 specs=[('cluster_loso',c,None) for c in clusters]+[('year_loyo',None,y) for y in years]
 specs += [('cluster_year_disjoint',c,int(y)) for c,y in sorted(set(zip(d.spatial_cluster,d.harvest_year)))]
 specs += [('past_years',None,y) for y in years]
 all_folds=[];member=[]
 for axis,c,y in specs:
  if axis=='cluster_loso':test=d.spatial_cluster.eq(c);train=~test
  elif axis=='year_loyo':test=d.harvest_year.eq(y);train=~test
  elif axis=='cluster_year_disjoint':test=d.spatial_cluster.eq(c)&d.harvest_year.eq(y);train=~d.spatial_cluster.eq(c)&~d.harvest_year.eq(y)
  else:test=d.harvest_year.eq(y);train=d.harvest_year.lt(y)
  tr,te=d.loc[train],d.loc[test];nprior=tr.harvest_year.nunique() if axis=='past_years' else None
  reasons=[]
  if not len(te):reasons.append('no_test_cases')
  if len(tr)<min_train:reasons.append('training_rows_below_300')
  if axis=='past_years' and nprior<min_prior:reasons.append('fewer_than_3_prior_years')
  fold_id=f'f{sum(f["axis"]==axis for f in all_folds):03d}'
  f={'axis':axis,'fold_id':fold_id,'held_cluster':c,'held_year':y,'train_indices':np.flatnonzero(train).tolist(),'test_indices':np.flatnonzero(test).tolist(),'n_train':len(tr),'n_test':len(te),'n_excluded':int((~(train|test)).sum()),'train_clusters':sorted(tr.spatial_cluster.unique()),'test_clusters':sorted(te.spatial_cluster.unique()),'train_years':list(map(int,sorted(tr.harvest_year.unique()))),'test_years':list(map(int,sorted(te.harvest_year.unique()))),'train_date_max':str(tr.date.max()) if len(tr) else None,'test_date_min':str(te.date.min()),'test_novel_sites':sorted(set(te.site_id)-set(tr.site_id)),'test_seen_sites':sorted(set(te.site_id)&set(tr.site_id)),'train_record_ids_sha256':engine.ids_hash(tr.record_id),'test_record_ids_sha256':engine.ids_hash(te.record_id),'combined_domain_shift':False,'prior_year_count':int(nprior) if nprior is not None else None,'supported':not reasons,'unsupported_reasons':reasons}
  if (train&test).any():raise ValueError('Train/test overlap')
  if axis in ['cluster_loso','cluster_year_disjoint'] and set(tr.spatial_cluster)&set(te.spatial_cluster):raise ValueError('Cluster leakage')
  if axis in ['year_loyo','cluster_year_disjoint','past_years'] and set(tr.harvest_year)&set(te.harvest_year):raise ValueError('Year leakage')
  all_folds.append(f);member.append(pd.DataFrame({'axis':axis,'fold_id':fold_id,'record_id':d.record_id,'role':np.where(test,'test',np.where(train,'train','excluded'))}))
 return all_folds,pd.concat(member,ignore_index=True)

def summarize_cohort(pair,folds):
 d=pair['VUT_REF'];delta=pair['CUT_REF'][TARGET]-d[TARGET]
 return {'cohort_rows':len(d),'sites':int(d.site_id.nunique()),'clusters':int(d.spatial_cluster.nunique()),'site_seasons':len(d[['site_id','harvest_year']].drop_duplicates()),'scope':'calendar_proxy_minus_documented_outside_bounds','strict_fullflux_approved':False,'scientific_primary_approved':False,'target_specific_no_ET_QC_condition':True,'ET_equal_all_source_rows':True,'calendar_conflict_selected_cases':int(d.calendar_sowing_conflict.sum()),'guarded_known_outside_selected':int(d.known_outside_documented_bounds.sum()),'reference_delta':{'n_exact_equal':int((delta==0).sum()),'n_different':int((delta!=0).sum()),'mean_cut_minus_vut':float(delta.mean()),'rms_reference_difference':float(np.sqrt(np.mean(delta**2))),'max_abs_difference':float(delta.abs().max())},'fold_inventory':{axis:{'planned':sum(f['axis']==axis for f in folds),'supported':sum(f['axis']==axis and f['supported'] for f in folds),'test_rows':sum(f['n_test'] for f in folds if f['axis']==axis and f['supported']),'min_train_rows':min([f['n_train'] for f in folds if f['axis']==axis and f['supported']],default=None)} for axis in engine.AXES},'expected_fits':2*5*sum(f['supported'] for f in folds)}

def protocol(config):
 if engine.sha(Path(__file__).parent/'base_model_protocol.json')!=BASE_PROTOCOL_SHA:raise ValueError('Fixed model protocol changed')
 p=json.loads((Path(__file__).parent/'base_model_protocol.json').read_text())
 p.update(protocol_version='3.0.0',study='Paired US NT VUT/CUT processing-path sensitivity on common fixed point inputs',target='GPP',variants=VARIANTS,reference_question='composite_processing_path_not_isolated_ustar',expected_rows=config['expected_common_rows'],source_config=config,primary_comparisons=[{'axis':'cluster_loso','fit_reference':v,'score_reference':v,'contrast':'random_forest versus ridge','metric':m} for v in VARIANTS for m in ['macro_site_season_rmse','macro_cluster_rmse']],cross_scoring='Complete 2 train-reference by 2 scoring-reference matrix; dependent diagnostics only',family_ranking='All five fixed families, each axis and macro metric, ascending RMSE; exact ties receive minimum rank; no selection',flag_meanings={'actual_crop_window_applied':'Legacy alias: False means no fully verified actual-season window was applied to every case; individual known bounds were quarantined','documented_bound_quarantine_applied':True,'all_actual_crop_windows_verified':False},interpretation=['Exact same physical cases, X and folds for both fits; no ET QC conditioning','Original scientific approval flags false; calendar proxy plus conservative documented-bound quarantine','Ne1/Ne3 one Mead cluster; Ne3 2013 sowing conflict retained','No crop-pure/full dynamic footprint, exact full season, independent truth, forecast, physiology or agent superiority claim','No actual DD AUXNEE member metadata in acquired archives; VUT/CUT is composite processing-path sensitivity','All five models and four axes retained; no outcome-based reference/model/split/site/QC selection'])
 for c in ['expected_reference_files','expected_fold_counts','expected_test_rows']:p.pop(c,None)
 return p

def prepare(args,preview=False):
 config=json.loads(args.config.read_text());pair,inputs=load_inputs(config)
 folds,membership=plan_splits(pair['VUT_REF']);summary=summarize_cohort(pair,folds)
 if preview:
  if args.output.exists():raise FileExistsError(args.output)
  args.output.mkdir(parents=True);write_json(args.output/'candidate_summary.json',summary)
  pair['VUT_REF'].groupby(['site_id','spatial_cluster','harvest_year']).size().rename('n').reset_index().to_csv(args.output/'cohort_summary.csv',index=False)
  pd.DataFrame([{k:v for k,v in f.items() if k not in ['train_indices','test_indices']} for f in folds]).to_csv(args.output/'fold_support.csv',index=False)
  print(json.dumps(summary));return
 if config.get('root_authorized_prefit_freeze') is not True:raise ValueError('Root prefit review authorization pending')
 receipt=config.get('root_independent_selection_receipt')
 if not isinstance(receipt,dict) or engine.sha(receipt['path'])!=receipt['sha256']:raise ValueError('Root independent selection receipt absent/changed')
 inputs[receipt['path']]=receipt['sha256']
 if args.output.exists():raise FileExistsError(args.output)
 if engine.sha(Path(__file__).parent/'engine.py')!=ENGINE_SHA:raise ValueError('Fixed numerical engine changed')
 args.output.mkdir(parents=True);p=protocol(config);write_json(args.output/'protocol.json',p)
 for var,d in pair.items():d.to_parquet(args.output/(var+'_snapshot.parquet'),index=False)
 membership.to_parquet(args.output/'all_membership.parquet',index=False);write_json(args.output/'all_folds.json',folds)
 write_json(args.output/'supported_folds.json',[f for f in folds if f['supported']]);write_json(args.output/'cohort_summary.json',summary)
 pair['VUT_REF'].groupby(['site_id','spatial_cluster','harvest_year']).size().rename('n').reset_index().to_csv(args.output/'cohort_summary.csv',index=False)
 pd.DataFrame([{k:v for k,v in f.items() if k not in ['train_indices','test_indices']} for f in folds]).to_csv(args.output/'fold_support.csv',index=False)
 code={str(p.resolve()):engine.sha(p) for p in sorted(Path(__file__).parent.iterdir()) if p.suffix in ['.py','.json']}
 files={str(p.resolve()):engine.sha(p) for p in args.output.iterdir() if p.is_file()}
 write_json(args.output/'freeze.json',{'created_utc':now(),'status':'frozen_before_any_fit_after_root_selection_review','inputs':inputs,'code':code,'plan_files':files,'protocol_sha256':engine.sha(args.output/'protocol.json'),'summary':summary,'historical_score_exposure':True,'genuinely_unseen_claim':False})
 print(json.dumps(summary))

def require_frozen(plan):
 f=json.loads((plan/'freeze.json').read_text())
 for name in ['inputs','code','plan_files']:
  for path,digest in f[name].items():
   if engine.sha(path)!=digest:raise ValueError('Frozen '+name+' changed: '+path)
 if f['status']!='frozen_before_any_fit_after_root_selection_review':raise ValueError('Missing prefit freeze')
 return f

def cross_score(preds,pair):
 result={};long=[];key=['axis','fold_id','model','target','seed','record_id']
 for fit_ref,values in preds.items():
  if values.duplicated(key).any():raise ValueError('Duplicate prediction scope')
  for score_ref in VARIANTS:
   ref=pair[score_ref][['record_id',TARGET]].rename(columns={TARGET:'scoring_truth'})
   q=values.merge(ref,on='record_id',validate='many_to_one');q['y_true']=q.pop('scoring_truth')
   tables=engine.score_tables(q)
   for name,t in tables.items():
    t=t.copy();t['fit_reference']=fit_ref;t['score_reference']=score_ref
    result.setdefault(name,[]).append(t)
   q['fit_reference']=fit_ref;q['score_reference']=score_ref;long.append(q)
 return {k:pd.concat(v,ignore_index=True) for k,v in result.items()},pd.concat(long,ignore_index=True)

def contrast_table(macro):
 rows=[]
 for axis in engine.AXES:
  for fit_ref in VARIANTS:
   for score_ref in VARIANTS:
    part=macro.loc[macro.axis.eq(axis)&macro.fit_reference.eq(fit_ref)&macro.score_reference.eq(score_ref)]
    if part.empty:continue
    for metric in ['macro_site_season_rmse','macro_cluster_rmse']:
     a=float(part.loc[part.model.eq('random_forest'),metric].iloc[0]);b=float(part.loc[part.model.eq('ridge'),metric].iloc[0])
     rows.append({'axis':axis,'fit_reference':fit_ref,'score_reference':score_ref,'contrast':'RF_minus_Ridge','metric':metric,'random_forest':a,'ridge':b,'absolute_difference':a-b,'relative_percent':100*(a-b)/b,'primary':axis=='cluster_loso' and fit_ref==score_ref})
 return pd.DataFrame(rows)

def execute(args):
 if args.workers!=4:raise ValueError('Fixed four-worker execution required')
 frozen=require_frozen(args.plan)
 if args.output.exists():raise FileExistsError(args.output)
 p=json.loads((args.plan/'protocol.json').read_text());p['_protocol_sha256']=frozen['protocol_sha256']
 folds=json.loads((args.plan/'supported_folds.json').read_text());args.output.mkdir(parents=True)
 preds={};pair={};allfits=[];start=time.perf_counter()
 for var in VARIANTS:
  root=args.output/var;root.mkdir();snapshot=args.plan/(var+'_snapshot.parquet');pair[var]=pd.read_parquet(snapshot)
  tasks=[(f,m,'GPP',p,str(root)) for f in folds for m in engine.MODELS]
  saved=[];fits=[]
  with ProcessPoolExecutor(max_workers=args.workers,initializer=engine.worker_init,initargs=(snapshot,)) as pool:
   futures=[pool.submit(engine.fit_task,t) for t in tasks]
   for future in as_completed(futures):
    values,fit=future.result();saved.append(values);fits.append(fit)
  values=pd.concat(saved,ignore_index=True).sort_values(['axis','fold_id','model','record_id']).reset_index(drop=True)
  values.to_parquet(root/'predictions.parquet',index=False);values.to_csv(root/'predictions.csv.gz',index=False,float_format='%.17g')
  for name,t in engine.score_tables(values).items():t.to_csv(root/('scores_'+name+'.csv'),index=False,float_format='%.17g')
  write_json(root/'fit_log.json',sorted(fits,key=lambda z:(z['axis'],z['fold_id'],z['model'])))
  preds[var]=values;allfits+=fits
 tables,cross=cross_score(preds,pair)
 summary=args.output/'summary';summary.mkdir()
 for name,t in tables.items():t.to_csv(summary/('cross_scores_'+name+'.csv'),index=False,float_format='%.17g')
 cross.to_parquet(args.output/'cross_predictions.parquet',index=False)
 contrast_table(tables['macro_summary']).to_csv(summary/'prespecified_contrasts.csv',index=False,float_format='%.17g')
 ranks=[]
 for metric in ['macro_site_season_rmse','macro_cluster_rmse']:
  z=tables['macro_summary'][['axis','fit_reference','score_reference','model',metric]].copy()
  z['metric']=metric;z['rank']=z.groupby(['axis','fit_reference','score_reference'])[metric].rank(method='min',ascending=True)
  z=z.rename(columns={metric:'rmse'});ranks.append(z)
 pd.concat(ranks,ignore_index=True).to_csv(summary/'model_rankings.csv',index=False,float_format='%.17g')
 key=['axis','fold_id','model','target','seed','record_id']
 a=preds['VUT_REF'];b=preds['CUT_REF'];q=a.merge(b,on=key,suffixes=('_vut','_cut'),validate='one_to_one')
 q['prediction_delta_cut_minus_vut']=q.y_pred_cut-q.y_pred_vut;q['reference_delta_cut_minus_vut']=q.y_true_cut-q.y_true_vut
 q.to_parquet(args.output/'paired_prediction_deltas.parquet',index=False)
 delta=[]
 for group,g in q.groupby(['axis','model'],sort=True):
  delta.append({'axis':group[0],'model':group[1],'n':len(g),'mean_prediction_delta':float(g.prediction_delta_cut_minus_vut.mean()),'rms_prediction_delta':float(np.sqrt(np.mean(g.prediction_delta_cut_minus_vut**2))),'mean_reference_delta':float(g.reference_delta_cut_minus_vut.mean()),'rms_reference_delta':float(np.sqrt(np.mean(g.reference_delta_cut_minus_vut**2)))})
 pd.DataFrame(delta).to_csv(summary/'reference_and_prediction_delta.csv',index=False,float_format='%.17g')
 require_frozen(args.plan)
 if len(allfits)!=frozen['summary']['expected_fits']:raise ValueError('Fit count mismatch')
 write_json(args.output/'execution.json',{'created_utc':now(),'status':'complete','fit_count':len(allfits),'prediction_records_per_reference':{k:len(v) for k,v in preds.items()},'cross_prediction_records':len(cross),'elapsed_seconds':time.perf_counter()-start,'workers':args.workers,'estimator_threads':1,'gpu_used':False,'no_ET_fit':True,'frozen_inputs_rechecked_unchanged':True})
 manifest={str(p.relative_to(args.output)):{'bytes':p.stat().st_size,'sha256':engine.sha(p)} for p in sorted(args.output.rglob('*')) if p.is_file()}
 write_json(args.output/'run_manifest.json',{'created_utc':now(),'files':manifest})
 print(json.dumps({'status':'complete','fit_count':len(allfits),'elapsed_seconds':time.perf_counter()-start}))

def cli():
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
 for name in ['preview','prepare']:
  q=sub.add_parser(name);q.add_argument('--config',type=Path,required=True);q.add_argument('--output',type=Path,required=True)
 q=sub.add_parser('execute');q.add_argument('--plan',type=Path,required=True);q.add_argument('--output',type=Path,required=True);q.add_argument('--workers',type=int,default=4)
 args=p.parse_args()
 if args.command=='execute':execute(args)
 else:prepare(args,preview=args.command=='preview')

if __name__=='__main__':cli()
