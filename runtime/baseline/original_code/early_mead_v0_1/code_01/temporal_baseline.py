"""Prepare and execute a source-frozen, same-region retrospective temporal diagnostic.

Preparation reads approved inputs/processed references but never predictions.
Execution needs the exact freeze and an independent prefit verification receipt.
The prior 642-case release and earlier preparation contracts remain unchanged.
"""
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse, hashlib, json, platform, sys, time
import numpy as np
import pandas as pd
import sklearn

ARMS=['B_uniform_pixels','C_core_weighted']
REFS=['VUT_REF','CUT_REF']
FEATURES=['NDVI_modis','EA_kPa_agera5','WS_ms_agera5','P_mm_agera5','TA_max_C_agera5','TA_min_C_agera5','SW_IN_MJ_m2_d_agera5','VPD_kPa_agera5']
TRAIN_SEASONS={('US-Ne1',2012),('US-Ne1',2013),('US-Ne3',2013)}
TEST_BOUNDS={('US-Ne1',2010):('2010-04-19','2010-09-21'),('US-Ne1',2011):('2011-05-18','2011-10-26'),('US-Ne3',2011):('2011-05-03','2011-10-19')}
ENGINE_SHA='c6dab8c357fd0432dbf8b06fa263d04e19b4527ad41984982862ad69409b2e9d'
BASE_PROTOCOL_SHA='a558444a518d5201af429b858fd4533ce1e6ad138056ccc10a81b0872af1accb'

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def now(): return datetime.now(timezone.utc).isoformat()
def read(p): return json.loads(Path(p).read_text(encoding='utf-8'))
def save(p,v):
 with Path(p).open('x',encoding='utf-8') as f: json.dump(v,f,indent=2,ensure_ascii=False,allow_nan=False);f.write('\n')
def ids(values): return hashlib.sha256('\n'.join(map(str,values)).encode()).hexdigest()
def runtime():return {'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'sklearn':sklearn.__version__}
def load_engine(path):
 import importlib.util
 path=Path(path); assert sha(path)==ENGINE_SHA
 spec=importlib.util.spec_from_file_location('engine',path); mod=importlib.util.module_from_spec(spec)
 sys.modules['engine']=mod;spec.loader.exec_module(mod);return mod
def norm(frame):
 d=frame.copy();d['date']=pd.to_datetime(d.date).dt.strftime('%Y-%m-%d')
 assert not d[['site_id','date']].duplicated().any()
 return d.sort_values(['site_id','date']).reset_index(drop=True)
def physical(d): return list(zip(d.site_id,d.date))
def check_train_test(train,test):
 assert len(train)==480
 assert set(zip(train.site_id,pd.to_datetime(train.date).dt.year))==TRAIN_SEASONS
 assert set(zip(test.site_id,pd.to_datetime(test.date).dt.year))==set(TEST_BOUNDS)
 assert not set(physical(train))&set(physical(test)), 'Physical train/test overlap'
 assert not set(train.record_id)&set(test.record_id), 'Record train/test overlap'
 assert pd.to_datetime(test.date).max()<pd.to_datetime(train.date).min(), 'Expected retrospective hindcast ordering'
 assert train.spatial_cluster.eq('Mead_NE').all() and test.spatial_cluster.eq('Mead_NE').all()
 assert np.isfinite(train[FEATURES+['gpp_gC_m2_d']].to_numpy(float)).all()
 assert np.isfinite(test[FEATURES+['gpp_gC_m2_d']].to_numpy(float)).all()

def prepare(args):
 assert not args.output.exists(), 'Use a new plan directory'
 c=read(args.config); sources={}
 for item in c['source_files']:
  p=Path(item['path']);assert sha(p)==item['sha256'],p;sources[str(p)]=item['sha256']
 sources[str(args.config.resolve())]=sha(args.config)
 required=[c[k] for k in ['base_protocol','engine','optical_daily','weather7','paired_target_payload','native_qc_receipt','optical_qa_receipt']]
 required += [c['training_snapshots'][a][r] for a in ARMS for r in REFS]
 assert {str(Path(p).resolve()) for p in required}<={str(Path(p).resolve()) for p in sources}, 'Incomplete source freeze'
 assert c['receipt_checks'], 'Exact final QA receipt checks required'
 for item in c['receipt_checks']:
  assert item['path'] in [c['native_qc_receipt'],c['optical_qa_receipt']]
  receipt=read(item['path'])
  for key,value in item['required_equal'].items():assert receipt[key]==value,(item['path'],key)
 assert {x['path'] for x in c['receipt_checks']}=={c['native_qc_receipt'],c['optical_qa_receipt']}
 assert sha(c['base_protocol'])==BASE_PROTOCOL_SHA
 engine=load_engine(c['engine']); base=read(c['base_protocol'])
 optical=norm(pd.read_csv(c['optical_daily']))
 assert len(optical)==552
 assert set(zip(optical.site_id,optical.year))==set(TEST_BOUNDS)
 assert np.array_equal(optical.B_uniform_pixels.notna(),optical.C_core_weighted.notna())
 assert np.array_equal(optical.complete_BC_optical.astype(str).str.lower().eq('true'),optical[ARMS].notna().all(axis=1))
 weather=norm(pd.read_csv(c['weather7'])); assert len(weather)==476
 targets=norm(pd.read_parquet(c['paired_target_payload']))
 ycols=c['target_columns']; assert set(ycols)==set(REFS)
 assert len(targets)<=450 and np.isfinite(targets[list(ycols.values())].to_numpy(float)).all()
 assert set(physical(targets))<=set(physical(weather))
 common=weather[['site_id','date',*FEATURES[1:]]].merge(optical[['site_id','date',*ARMS]],on=['site_id','date'],validate='one_to_one')
 assert len(common)==476
 common=common.merge(targets[['site_id','date',*ycols.values()]],on=['site_id','date'],validate='one_to_one')
 eligibility=common[['site_id','date']].copy()
 eligibility['complete_8X_both_arms']=np.isfinite(common[[*FEATURES[1:],*ARMS]].to_numpy(float)).all(axis=1)
 selected=common.loc[eligibility.complete_8X_both_arms].copy().reset_index(drop=True)
 assert len(selected)>0
 assert set(zip(selected.site_id,pd.to_datetime(selected.date).dt.year))==set(TEST_BOUNDS), 'All three preselected seasons required'
 for site,year in TEST_BOUNDS:
  d=selected.loc[selected.site_id.eq(site)&pd.to_datetime(selected.date).dt.year.eq(year)]
  lo,hi=TEST_BOUNDS[(site,year)]
  assert len(d)>0 and d.date.between(lo,hi).all() and pd.to_datetime(d.date).dt.month.between(5,10).all()
 selected['record_id']=['T1:'+s+':'+d for s,d in physical(selected)]
 selected['observation_key']=selected.record_id
 selected['system_id']=selected.site_id
 selected['spatial_cluster']='Mead_NE'
 selected['harvest_year']=pd.to_datetime(selected.date).dt.year.astype(int)
 selected['spatial_split_group']='Mead_NE'
 selected['season_split_group']=selected.site_id+'__'+selected.harvest_year.astype(str)
 frames={};train_frames={};test_frames={}
 for arm in ARMS:
  for ref in REFS:
   path=c['training_snapshots'][arm][ref]
   tr=norm(pd.read_parquet(path));assert len(tr)==480
   te=selected[engine.KEYS].copy()
   te['gpp_gC_m2_d']=selected[ycols[ref]].to_numpy(float)
   for feature in FEATURES:te[feature]=selected[arm if feature=='NDVI_modis' else feature].to_numpy(float)
   tr=tr[[*engine.KEYS,'gpp_gC_m2_d',*FEATURES]].copy()
   check_train_test(tr,te)
   train_frames[arm,ref]=tr;test_frames[arm,ref]=te
   frames[arm,ref]=pd.concat([tr,te],ignore_index=True)
 baseline=train_frames[ARMS[0],REFS[0]]
 for arm in ARMS:
  for ref in REFS:
   d=train_frames[arm,ref]
   pd.testing.assert_frame_equal(baseline[engine.KEYS],d[engine.KEYS],check_exact=True)
   np.testing.assert_array_equal(baseline[FEATURES[1:]].to_numpy(),d[FEATURES[1:]].to_numpy())
   np.testing.assert_array_equal(train_frames[arm,REFS[0]][FEATURES].to_numpy(),d[FEATURES].to_numpy())
   np.testing.assert_array_equal(test_frames[arm,REFS[0]][FEATURES].to_numpy(),test_frames[arm,ref][FEATURES].to_numpy())
 for ref in REFS:
  for collection in [train_frames,test_frames]:
   np.testing.assert_array_equal(collection[ARMS[0],ref].gpp_gC_m2_d,collection[ARMS[1],ref].gpp_gC_m2_d)
 assert not args.output.exists();args.output.mkdir(parents=True)
 eligibility.to_csv(args.output/'input_eligibility.csv',index=False)
 selected.to_parquet(args.output/'paired_test_source_snapshot.parquet',index=False)
 for (arm,ref),d in frames.items():d.to_parquet(args.output/(arm+'__'+ref+'_snapshot.parquet'),index=False)
 n=len(selected);fold={
  'axis':'temporal_hindcast','fold_id':'f000','n_train':480,'n_test':n,
  'train_indices':list(range(480)),'test_indices':list(range(480,480+n)),
  'train_record_ids_sha256':ids(baseline.record_id),'test_record_ids_sha256':ids(selected.record_id),
  'train_years':[2012,2013],'test_years':[2010,2011],
  'train_clusters':['Mead_NE'],'test_clusters':['Mead_NE'],
  'independent_new_region':False,'chronological_forecasting':False}
 save(args.output/'fold.json',fold)
 protocol={
  'created_utc':now(),'version':'T1_temporal_hindcast_v0_1','scope':'Processed NT-reference same-Mead retrospective temporal extension',
  'selection_basis':'All three publication-supported early maize seasons with same-year cores; exact finite B/C intersection after fixed source QC; no new score-driven selection',
  'seed':base['seed'],'models':base['models'],'features':FEATURES,'arms':ARMS,'references':REFS,
  'axis':'temporal_hindcast','primary_model':'random_forest','primary_contrast':'C_core_weighted_minus_B_uniform_pixels',
  'primary_metric':'equal_test_site_season_mean_RMSE','primary_paths':'VUT own and CUT own reported separately; both dependent processed products',
  'secondary':'All fixed five families, per-season and pooled metrics, both cross-reference scoring directions; no outcome-based choice',
  'test_site_season_counts':selected.groupby(['site_id','harvest_year']).size().rename('n').reset_index().to_dict('records'),
  'n_train':480,'n_test':n,'expected_fits':20,'expected_own_prediction_records':20*n,'expected_cross_prediction_records':40*n,
  'scientific_scope_limit':[
   'One existing region and correlated site-seasons; no population uncertainty/significance claim',
   'Test years precede training years; hindcast, not forecasting',
   'Unknown historical human exposure; retrospective role, never sealed confirmation',
   'Published normalized monthly 80% DAY cores, not full daily flux footprints',
   'Ne1 mast/point authority unresolved; fixed point QA makes selection conditional; no A point contrast',
   'No ET/cross-crop/agent-advantage or full-season-total claim'],
  'training_weights':'Uniform daily; no imputation/target clipping; train-only transforms',
  'tuning':'None; seed42, no model/threshold/coordinate/variant selection from scores',
  'resource_budget':{'model_fits':20,'workers':4,'estimator_threads':1,'GPU':False,'LLM_or_Cloud_benchmark_calls':0,'refit_after_scores':False},
  'runtime':runtime(),
  'parent_S1_regional_execution_released':False,'source_config_sha256':sha(args.config),
  'native_QC_receipt':c['native_qc_receipt'],'optical_QA_receipt':c['optical_qa_receipt']}
 save(args.output/'protocol.json',protocol)
 files={str(p.resolve()):sha(p) for p in args.output.iterdir() if p.is_file()}
 save(args.output/'freeze.json',{'status':'frozen_before_fit_pending_independent_prefit_verification','created_utc':now(),
  'inputs':sources,'code':{str(Path(__file__).resolve()):sha(__file__),str(Path(c['engine']).resolve()):sha(c['engine'])},
  'plan_files':files,'protocol_sha256':sha(args.output/'protocol.json'),'execution_allowed':False})
 print(json.dumps({'status':'prepared','n_train':480,'n_test':n,'fits_proposed':20,'protocol_sha256':sha(args.output/'protocol.json')}),flush=True)

def frozen(plan):
 f=read(plan/'freeze.json')
 for group in ['inputs','code','plan_files']:
  for path,h in f[group].items():assert sha(path)==h,(group,path)
 return f

def execute(args):
 assert not args.output.exists()
 f=frozen(args.plan);a=read(args.authorization)
 assert a['status']=='authorized_conditional_retrospective_T1'
 assert a['protocol_sha256']==f['protocol_sha256']
 q=read(a['independent_prefit_receipt'])
 assert sha(a['independent_prefit_receipt'])==a['independent_prefit_receipt_sha256'] and q['status']=='pass'
 assert q['protocol_sha256']==f['protocol_sha256']
 c=read(args.config);engine=load_engine(c['engine']);p=read(args.plan/'protocol.json');p['_protocol_sha256']=f['protocol_sha256']
 assert sha(args.config)==p['source_config_sha256']
 assert runtime()==p['runtime'], 'Numerical runtime changed after freeze'
 assert p['expected_fits']==20 and p['resource_budget']['workers']==4
 fold=read(args.plan/'fold.json');args.output.mkdir(parents=True)
 start=time.perf_counter();cross=[];fits=[];own=[]
 for arm in ARMS:
  for ref in REFS:
   snapshot=args.plan/(arm+'__'+ref+'_snapshot.parquet');out=args.output/arm/ref;out.mkdir(parents=True)
   tasks=[(fold,m,'GPP',p,str(out)) for m in engine.MODELS]
   rows=[]
   with ProcessPoolExecutor(max_workers=4,initializer=engine.worker_init,initargs=(snapshot,)) as pool:
    pending=[pool.submit(engine.fit_task,t) for t in tasks]
    for future in as_completed(pending):
     preds,fit=future.result();rows.append(preds);fits.append({**fit,'arm':arm,'fit_reference':ref})
   pred=pd.concat(rows).sort_values(['model','record_id']).reset_index(drop=True)
   pred.to_parquet(out/'own_predictions.parquet',index=False);pred['arm']=arm;pred['fit_reference']=ref;own.append(pred)
   for score_ref in REFS:
    truth=pd.read_parquet(args.plan/(arm+'__'+score_ref+'_snapshot.parquet')).iloc[fold['test_indices']][['record_id','gpp_gC_m2_d']]
    d=pred.drop(columns='y_true').merge(truth,on='record_id',validate='many_to_one').rename(columns={'gpp_gC_m2_d':'y_true'})
    d['score_reference']=score_ref;cross.append(d)
 own=pd.concat(own,ignore_index=True);cross=pd.concat(cross,ignore_index=True)
 assert len(fits)==20 and len(own)==p['expected_own_prediction_records'] and len(cross)==p['expected_cross_prediction_records']
 controls=own.loc[own.model.isin(['train_mean','seasonal_harmonic'])].pivot(index=['fit_reference','model','record_id'],columns='arm',values='y_pred')
 np.testing.assert_array_equal(controls[ARMS[0]],controls[ARMS[1]])
 own.to_parquet(args.output/'own_predictions.parquet',index=False);cross.to_parquet(args.output/'cross_predictions.parquet',index=False)
 groups=['arm','model','fit_reference','score_reference']
 pooled=[];season=[]
 for key,d in cross.groupby(groups,sort=True):
  pooled.append({**dict(zip(groups,key)),**engine.metric_values(d.y_true,d.y_pred)})
  for (site,year),s in d.groupby(['site_id','harvest_year'],sort=True):
   season.append({**dict(zip(groups,key)),'site_id':site,'harvest_year':int(year),**engine.metric_values(s.y_true,s.y_pred)})
 pooled=pd.DataFrame(pooled);season=pd.DataFrame(season)
 macro=season.groupby(groups,as_index=False).agg(macro_site_season_rmse=('rmse','mean'),macro_site_season_mae=('mae','mean'),n_site_seasons=('site_id','size'))
 contrasts=[]
 for (model,fit_ref,score_ref),g in macro.groupby(['model','fit_reference','score_reference']):
  b=float(g.loc[g.arm.eq(ARMS[0]),'macro_site_season_rmse'].iloc[0]);cval=float(g.loc[g.arm.eq(ARMS[1]),'macro_site_season_rmse'].iloc[0])
  contrasts.append({'model':model,'fit_reference':fit_ref,'score_reference':score_ref,'contrast':'C_minus_B','B_RMSE':b,'C_RMSE':cval,
   'absolute_difference':cval-b,'relative_percent':100*(cval-b)/b if b else None,
   'primary':model=='random_forest' and fit_ref==score_ref})
 for name,data in [('scores_pooled',pooled),('scores_site_season',season),('scores_macro',macro),('input_contrasts',pd.DataFrame(contrasts))]:data.to_csv(args.output/(name+'.csv'),index=False,float_format='%.17g')
 save(args.output/'fit_log.json',sorted(fits,key=lambda x:(x['arm'],x['fit_reference'],x['model'])))
 frozen(args.plan)
 save(args.output/'execution.json',{'created_utc':now(),'status':'complete','protocol_sha256':f['protocol_sha256'],
  'authorization_sha256':sha(args.authorization),'fits':20,'own_prediction_records':len(own),'cross_prediction_records':len(cross),
  'n_train':480,'n_test':p['n_test'],'seconds':time.perf_counter()-start,'controls_identical':True,
  'all_frozen_sources_unchanged':True,'new_region_or_forecast_claim':False,'Cloud_calls':0,'GPU':False})
 save(args.output/'run_manifest.json',{'files':{str(p.relative_to(args.output)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(args.output.rglob('*')) if p.is_file()}})
 print(json.dumps({'status':'complete','fits':20,'n_test':p['n_test'],'seconds':time.perf_counter()-start}),flush=True)

if __name__=='__main__':
 parser=argparse.ArgumentParser();sub=parser.add_subparsers(dest='command',required=True)
 q=sub.add_parser('prepare');q.add_argument('--config',type=Path,required=True);q.add_argument('--output',type=Path,required=True)
 q=sub.add_parser('execute');q.add_argument('--config',type=Path,required=True);q.add_argument('--plan',type=Path,required=True);q.add_argument('--output',type=Path,required=True);q.add_argument('--authorization',type=Path,required=True)
 args=parser.parse_args();prepare(args) if args.command=='prepare' else execute(args)
