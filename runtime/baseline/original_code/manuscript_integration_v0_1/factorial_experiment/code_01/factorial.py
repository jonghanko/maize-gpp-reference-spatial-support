"""Fixed same-case optical-input by processed-reference factorial diagnostic.

New source-evidence quarantine, not tuning. Earlier releases are read-only.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
import argparse,json,time
import numpy as np
import pandas as pd
import engine
import reference_adapter_frozen as ref

ARMS=['A_existing','A_shared_scene','B_uniform_pixels','C_core_weighted']
VARS=['VUT_REF','CUT_REF']
CONTRASTS=[('C_minus_B','C_core_weighted','B_uniform_pixels'),('C_minus_A_shared','C_core_weighted','A_shared_scene'),('A_existing_minus_A_shared','A_existing','A_shared_scene')]
METRICS=['macro_site_season_rmse','macro_cluster_rmse']
def now():return datetime.now(timezone.utc).isoformat()
def save(p,v):return ref.write_json(p,v)
def slug(arm,var):return arm+'__'+var

def load(config):
 inputs={}
 for key in ['source_pair_config','old_spatial_freeze','approved_selection_receipt','source_audit_summary','independent_source_receipt']:
  info=config[key];assert engine.sha(info['path'])==info['sha256'],key;inputs[info['path']]=info['sha256']
 pairconfig=json.loads(Path(config['source_pair_config']['path']).read_text())
 pair,sourceinputs=ref.load_inputs(pairconfig);inputs.update(sourceinputs)
 old=json.loads(Path(config['old_spatial_freeze']['path']).read_text())
 for family in ['inputs','code']:
  for path,h in old[family].items():assert engine.sha(path)==h,(family,path);inputs[path]=h
 optical={}
 for arm in ARMS:
  info=old['snapshots'][arm];assert engine.sha(info['path'])==info['sha256'];inputs[info['path']]=info['sha256']
  optical[arm]=ref.normalize(pd.read_parquet(info['path']))
  assert len(optical[arm])==681
 base=optical['A_existing'];keys=['site_id','date']
 for arm,d in optical.items():
  ref.same(base,d,keys+['record_id','gpp_gC_m2_d',*engine.FEATURES[1:]])
  assert np.isfinite(d[engine.FEATURES].to_numpy(float)).all()
 # Intersection is determined by native QC and the documented-bound source guard.
 paired=pair['VUT_REF'];common=base[keys].merge(paired[keys],on=keys,validate='one_to_one').sort_values(keys).reset_index(drop=True)
 accepted=ref.normalize(pd.read_parquet(config['approved_selection_receipt']['path']))
 ref.same(common,accepted,keys)
 assert len(common)==config['expected_rows']==642
 result={}
 for arm in ARMS:
  oldselected=common.merge(optical[arm],on=keys,validate='one_to_one')
  for var in VARS:
   x=common.merge(pair[var],on=keys,validate='one_to_one')
   if var=='VUT_REF':ref.same(x,oldselected,keys+['gpp_gC_m2_d',*engine.FEATURES[1:]])
   ref.same(x,oldselected,keys+engine.FEATURES[1:])
   point=common.merge(base,on=keys,validate='one_to_one')
   np.testing.assert_array_equal(x.NDVI_modis,point.NDVI_modis)
   x['record_id']='factorial|'+x.site_id+'|'+x.date
   x['old_spatial_record_id']=oldselected.record_id
   for name in ARMS:
    view=common.merge(optical[name],on=keys,validate='one_to_one')
    x[name]=view.NDVI_modis.to_numpy()
   x['NDVI_modis']=oldselected.NDVI_modis.to_numpy()
   x['optical_arm']=arm;x['benchmark_permission']='matched_input_reference_factorial_calendar_proxy_diagnostic'
   x['documented_bound_quarantine_applied']=True;x['all_actual_crop_windows_verified']=False
   x['scientific_primary_approved']=False;x['actual_crop_window_applied']=False
   assert len(x)==642 and np.isfinite(x[['gpp_gC_m2_d',*engine.FEATURES]].to_numpy(float)).all()
   result[(arm,var)]=x
 for arm in ARMS:ref.same(result[(arm,'VUT_REF')],result[(arm,'CUT_REF')],['record_id',*keys,*engine.FEATURES])
 for var in VARS:
  for arm in ARMS:ref.same(result[('A_existing',var)],result[(arm,var)],['record_id',*keys,'gpp_gC_m2_d',*engine.FEATURES[1:]])
 return result,inputs

def summarize(frames,folds):
 x=frames[('A_existing','VUT_REF')]
 return {'rows':len(x),'sites':int(x.site_id.nunique()),'regions':int(x.spatial_cluster.nunique()),'site_seasons':len(x[['site_id','harvest_year']].drop_duplicates()),'Ne3_sowing_conflict_cases':int(x.calendar_sowing_conflict.sum()),'Ro1_unknown_bounds_cases':int((x.site_id.eq('US-Ro1')&x.sowing_date_guard.isna()).sum()),'known_outside_selected':int(x.known_outside_documented_bounds.sum()),'scientific_primary_approved':False,'all_actual_crop_windows_verified':False,'source_old681_known_outside_quarantine':39,'fold_inventory':{axis:{'planned':sum(f['axis']==axis for f in folds),'supported':sum(f['axis']==axis and f['supported'] for f in folds),'test_rows_per_arm_reference_family':sum(f['n_test'] for f in folds if f['axis']==axis and f['supported']),'minimum_training':min([f['n_train'] for f in folds if f['axis']==axis and f['supported']],default=None)} for axis in engine.AXES},'expected_fits':4*2*5*sum(f['supported'] for f in folds),'expected_own_prediction_records':4*2*5*sum(f['n_test'] for f in folds if f['supported'])}

def build_protocol(config,summary):
 assert engine.sha(Path(__file__).parent/'base_model_protocol.json')==ref.BASE_PROTOCOL_SHA
 p=json.loads((Path(__file__).parent/'base_model_protocol.json').read_text())
 for key in ['expected_reference_files','expected_fold_counts','expected_test_rows','primary_comparisons']:p.pop(key,None)
 p.update(protocol_version='4.0.0',study='Evidence-updated fully matched optical-input by NT-reference factorial',target='GPP',arms=ARMS,variants=VARS,expected_rows=642,source_config=config,primary_axis='year_loyo',primary_metric='macro_site_season_rmse',primary_models=['ridge','random_forest'],primary_input_contrasts=['C_minus_B','C_minus_A_shared'],primary_reference_interaction='[RMSE_C−RMSE_B]_CUT_own−[RMSE_C−RMSE_B]_VUT_own; gC m−2 d−1; descriptive only',cross_scoring='Every arm/family/supported axis: full2fit_reference×2score_reference on exact same cases',controls=['A_existing−A_shared; expected identical inputs on source audit','Training mean and seasonal curve must be identical across optical arms within a reference'],summary=summary,reason_for_new_experiment='39 later-documented outside-bound dates and identifiability on exact common native-QC/input cases; not model-score tuning',interpretation=['Only complete LOYO is primary; regional OOF partial; all chronological strata unsupported','Ne1/Ne3 one Mead region; Ro1 inputs identical across arms, so arm-dependent Ro1 predictions arise through changed Mead training representation','Ne3 sowing conflict157 and Ro1 unknown bounds162 retained; no exact full season or strict crop/fullflux approval','B/C NDVI use weighted red/NIR first, then NDVI; published retained80%core conditional weights and shared scene/pixel QA unchanged from frozen input release','Reference×input interaction descriptive and dependent; not causal or population inference; neither NT path independent truth','No new acquisition, optical interpolation, pixel relocation, ET fit, tuning or clipping','All five fixed families retained; no outer-score-driven model/ref/site/threshold selection','VUT/CUT composite processing-path contrast; exact DD AUXNEE member metadata absent'])
 return p

def prepare(args):
 if args.output.exists():raise FileExistsError(args.output)
 config=json.loads(args.config.read_text());frames,inputs=load(config)
 frame=frames[('A_existing','VUT_REF')];folds,membership=ref.plan_splits(frame,300,3)
 summary=summarize(frames,folds)
 assert {k:v['supported'] for k,v in summary['fold_inventory'].items()}=={'cluster_loso':1,'year_loyo':3,'cluster_year_disjoint':1,'past_years':0}
 assert engine.sha(Path(__file__).parent/'engine.py')==ref.ENGINE_SHA
 args.output.mkdir(parents=True);p=build_protocol(config,summary);save(args.output/'protocol.json',p)
 for (arm,var),x in frames.items():x.to_parquet(args.output/(slug(arm,var)+'_snapshot.parquet'),index=False)
 save(args.output/'all_folds.json',folds);save(args.output/'supported_folds.json',[f for f in folds if f['supported']])
 membership.to_parquet(args.output/'all_membership.parquet',index=False);save(args.output/'cohort_summary.json',summary)
 frame.groupby(['site_id','spatial_cluster','harvest_year']).size().rename('n').reset_index().to_csv(args.output/'site_season_counts.csv',index=False)
 pd.DataFrame([{k:v for k,v in f.items() if k not in ['train_indices','test_indices']} for f in folds]).to_csv(args.output/'fold_support.csv',index=False)
 code={str(p.resolve()):engine.sha(p) for p in sorted(Path(__file__).parent.iterdir()) if p.suffix in ['.py','.json']}
 files={str(p.resolve()):engine.sha(p) for p in args.output.iterdir() if p.is_file()}
 save(args.output/'freeze.json',{'created_utc':now(),'status':'frozen_before_any_fit_pending_root_final_review','inputs':inputs,'code':code,'plan_files':files,'protocol_sha256':engine.sha(args.output/'protocol.json'),'summary':summary,'historical_score_exposure':True,'genuinely_unseen_claim':False})
 print(json.dumps(summary))

def frozen(plan):
 f=json.loads((plan/'freeze.json').read_text())
 for group in ['inputs','code','plan_files']:
  for path,h in f[group].items():assert engine.sha(path)==h,(group,path)
 return f

def input_contrasts(macro):
 rows=[]
 for (axis,model,fit,score),g in macro.groupby(['axis','model','fit_reference','score_reference'],sort=True):
  for name,first,second in CONTRASTS:
   for metric in METRICS:
    a=float(g.loc[g.arm.eq(first),metric].iloc[0]);b=float(g.loc[g.arm.eq(second),metric].iloc[0])
    rows.append({'axis':axis,'model':model,'fit_reference':fit,'score_reference':score,'contrast':name,'metric':metric,'first_arm':first,'second_arm':second,'first_RMSE':a,'second_RMSE':b,'absolute_RMSE_difference':a-b,'relative_percent':100*(a-b)/b,'primary':axis=='year_loyo' and metric==METRICS[0] and model in ['ridge','random_forest'] and fit==score and name in ['C_minus_B','C_minus_A_shared']})
 return pd.DataFrame(rows)

def interactions(contrasts):
 rows=[]
 for (axis,model,name,metric),g in contrasts.groupby(['axis','model','contrast','metric'],sort=True):
  get=lambda fit,score:float(g.loc[g.fit_reference.eq(fit)&g.score_reference.eq(score),'absolute_RMSE_difference'].iloc[0])
  primary=axis=='year_loyo' and model in ['ridge','random_forest'] and name=='C_minus_B' and metric==METRICS[0]
  rows.append({'axis':axis,'model':model,'contrast':name,'metric':metric,'interaction_type':'CUT_own_minus_VUT_own','fixed_reference':None,'interaction_RMSE_units':get('CUT_REF','CUT_REF')-get('VUT_REF','VUT_REF'),'primary_descriptive':primary})
  for var in VARS:
   rows.append({'axis':axis,'model':model,'contrast':name,'metric':metric,'interaction_type':'training_path_effect_fixed_scoring_reference','fixed_reference':var,'interaction_RMSE_units':get('CUT_REF',var)-get('VUT_REF',var),'primary_descriptive':False})
   rows.append({'axis':axis,'model':model,'contrast':name,'metric':metric,'interaction_type':'scoring_label_effect_fixed_training_reference','fixed_reference':var,'interaction_RMSE_units':get(var,'CUT_REF')-get(var,'VUT_REF'),'primary_descriptive':False})
 return pd.DataFrame(rows)

def execute(args):
 if args.workers!=4:raise ValueError('Fixed4workers')
 f=frozen(args.plan);assert not args.output.exists();args.output.mkdir(parents=True)
 p=json.loads((args.plan/'protocol.json').read_text());p['_protocol_sha256']=f['protocol_sha256']
 folds=json.loads((args.plan/'supported_folds.json').read_text());start=time.perf_counter();fitcount=0;allcross=[];allmetrics={};predcounts={}
 for arm in ARMS:
  pairs={};preds={}
  for var in VARS:
   out=args.output/arm/var;out.mkdir(parents=True);snapshot=args.plan/(slug(arm,var)+'_snapshot.parquet');pairs[var]=pd.read_parquet(snapshot)
   tasks=[(fold,m,'GPP',p,str(out)) for fold in folds for m in engine.MODELS];rows=[];fits=[]
   with ProcessPoolExecutor(max_workers=4,initializer=engine.worker_init,initargs=(snapshot,)) as pool:
    futures=[pool.submit(engine.fit_task,t) for t in tasks]
    for future in as_completed(futures):
     q,fit=future.result();rows.append(q);fits.append(fit)
   q=pd.concat(rows,ignore_index=True).sort_values(['axis','fold_id','model','record_id']).reset_index(drop=True)
   q.to_parquet(out/'predictions.parquet',index=False);q.to_csv(out/'predictions.csv.gz',index=False,float_format='%.17g')
   for name,t in engine.score_tables(q).items():t.to_csv(out/('scores_'+name+'.csv'),index=False,float_format='%.17g')
   save(out/'fit_log.json',sorted(fits,key=lambda z:(z['axis'],z['fold_id'],z['model'])))
   preds[var]=q;predcounts[slug(arm,var)]=len(q);fitcount+=len(fits)
  metrics,cross=ref.cross_score(preds,pairs);cross['arm']=arm;allcross.append(cross)
  for name,t in metrics.items():t['arm']=arm;allmetrics.setdefault(name,[]).append(t)
 # Invariant calendar/mean controls have no optical input and must match exactly.
 long=pd.concat(allcross,ignore_index=True)
 for var in VARS:
  controls=long.loc[long.fit_reference.eq(var)&long.score_reference.eq(var)&long.model.isin(['train_mean','seasonal_harmonic'])]
  pivot=controls.pivot(index=['axis','fold_id','model','record_id'],columns='arm',values='y_pred')
  for arm in ARMS:np.testing.assert_array_equal(pivot[arm],pivot.A_existing)
 long.to_parquet(args.output/'cross_predictions.parquet',index=False)
 tables={name:pd.concat(t,ignore_index=True) for name,t in allmetrics.items()};summary=args.output/'summary';summary.mkdir()
 for name,t in tables.items():t.to_csv(summary/('cross_scores_'+name+'.csv'),index=False,float_format='%.17g')
 contrast=input_contrasts(tables['macro_summary']);contrast.to_csv(summary/'input_contrasts.csv',index=False,float_format='%.17g')
 interactions(contrast).to_csv(summary/'reference_input_interactions.csv',index=False,float_format='%.17g')
 ranks=[]
 for metric in METRICS:
  z=tables['macro_summary'][['axis','arm','model','fit_reference','score_reference',metric]].copy();z['metric']=metric
  z['rank']=z.groupby(['axis','arm','fit_reference','score_reference'])[metric].rank(method='min');z=z.rename(columns={metric:'RMSE'});ranks.append(z)
 pd.concat(ranks,ignore_index=True).to_csv(summary/'family_rankings.csv',index=False,float_format='%.17g')
 frozen(args.plan);assert fitcount==f['summary']['expected_fits']==200
 save(args.output/'execution.json',{'created_utc':now(),'status':'complete','fits':fitcount,'own_predictions_by_arm_reference':predcounts,'own_prediction_records':sum(predcounts.values()),'cross_prediction_records':len(long),'elapsed_seconds':time.perf_counter()-start,'workers':4,'estimator_threads':1,'unsupported_fits':0,'ET_fits':0,'controls_identical_across_arms':True,'all_frozen_files_unchanged':True})
 save(args.output/'run_manifest.json',{'created_utc':now(),'files':{str(p.relative_to(args.output)):{'bytes':p.stat().st_size,'sha256':engine.sha(p)} for p in sorted(args.output.rglob('*')) if p.is_file()}})
 print(json.dumps({'status':'complete','fits':fitcount,'seconds':time.perf_counter()-start}))

if __name__=='__main__':
 p=argparse.ArgumentParser();s=p.add_subparsers(dest='command',required=True)
 q=s.add_parser('prepare');q.add_argument('--config',type=Path,required=True);q.add_argument('--output',type=Path,required=True)
 q=s.add_parser('execute');q.add_argument('--plan',type=Path,required=True);q.add_argument('--output',type=Path,required=True);q.add_argument('--workers',type=int,default=4)
 args=p.parse_args();prepare(args) if args.command=='prepare' else execute(args)
