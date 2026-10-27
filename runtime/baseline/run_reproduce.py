"""Offline reproduction from supplied processed inputs and frozen outputs.

Default replay requires numpy/pandas, not PyArrow or scikit-learn. Optional refit
uses the immutable archived engine, original snapshots/folds and fixed recipes.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,hashlib,os,sys,importlib.util,time
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
from portable_io import read_frame,NATIVE_MAP_CHECKS

ROOT=Path(__file__).resolve().parent
GROUP=['axis','arm','model','target','seed','fit_reference','score_reference']
SORT_KEYS=['cohort','scope','axis','arm','model','fit_reference','score_reference','target','seed','fold_id','site_id','harvest_year','spatial_cluster','unit','contrast','first_arm','second_arm','year','month','date','anchor','item_id']
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def load(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def save(path,data):Path(path).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def install_source_tree_guard(roots,analysis_root,allowed):
    from public_guard import install
    install(allowed)
def validate():
    man=load(ROOT/'MANIFEST.json');count=0
    for f in man['files']:
        p=ROOT/f['path']
        if not p.is_file() or p.stat().st_size!=f['bytes'] or sha(p)!=f['sha256']:
            raise ValueError('MANIFEST_MISMATCH: '+f['path'])
        count+=1
    return count

def metrics(y,p):
    y=np.asarray(y,dtype=float);p=np.asarray(p,dtype=float)
    if not len(y) or y.shape!=p.shape or not np.isfinite(y).all() or not np.isfinite(p).all():raise ValueError('Incomplete paired predictions')
    e=p-y;sse=float(np.dot(e,e));sst=float(np.sum((y-y.mean())**2))
    skill=1-sse/sst if len(y)>=3 and np.ptp(y)>0 and sst>0 else None
    r=float(np.corrcoef(y,p)[0,1]) if len(y)>=3 and np.ptp(y)>0 and np.ptp(p)>0 else None
    return dict(n=len(y),rmse=float(np.sqrt(np.mean(e**2))),mae=float(np.mean(abs(e))),bias_pred_minus_obs=float(np.mean(e)),nse=skill,r2_reference_variance=skill,pearson_r=r,squared_pearson_r=r*r if r is not None else None)

def score_tables(frame):
    out={}
    for name,extra in [('pooled',[]),('fold',['fold_id']),('cluster',['spatial_cluster']),('site_season',['site_id','spatial_cluster','harvest_year','season_split_group'])]:
        groups=GROUP+extra;rows=[]
        for keys,g in frame.groupby(groups,sort=True):rows.append({**dict(zip(groups,keys)),**metrics(g.y_true,g.y_pred)})
        out[name]=pd.DataFrame(rows)
    rows=[]
    for keys,g in out['site_season'].groupby(GROUP,sort=True):
        info=dict(zip(GROUP,keys));mask=np.ones(len(out['cluster']),dtype=bool)
        for k,v in info.items():mask&=out['cluster'][k].eq(v).to_numpy()
        c=out['cluster'].loc[mask];valid=g.nse.dropna()
        rows.append({**info,'n_site_seasons':len(g),'n_clusters':len(c),'macro_site_season_rmse':float(g.rmse.mean()),'median_site_season_rmse':float(g.rmse.median()),'macro_site_season_mae':float(g.mae.mean()),'macro_site_season_bias':float(g.bias_pred_minus_obs.mean()),'macro_site_season_nse':float(valid.mean()) if len(valid) else None,'n_finite_season_nse':len(valid),'macro_cluster_rmse':float(c.rmse.mean())})
    out['macro_summary']=pd.DataFrame(rows)
    return out

def compare_frames(actual,expected,label):
    missing=set(expected.columns)-set(actual.columns)
    if missing:raise ValueError(label+' missing columns '+repr(sorted(missing)))
    actual=actual[expected.columns].copy();expected=expected.copy()
    if len(actual)!=len(expected):raise ValueError(label+' row count mismatch '+str((len(actual),len(expected))))
    keys=[k for k in SORT_KEYS if k in expected.columns]
    if keys:
        for f in [actual,expected]:
            for k in keys:
                if k in ['seed','harvest_year','year','month']:f[k]=pd.to_numeric(f[k],errors='raise')
        actual=actual.sort_values(keys,kind='stable',na_position='last').reset_index(drop=True);expected=expected.sort_values(keys,kind='stable',na_position='last').reset_index(drop=True)
    else:actual=actual.reset_index(drop=True);expected=expected.reset_index(drop=True)
    cells=0;maxdiff=0.
    for col in expected.columns:
        a=actual[col];b=expected[col]
        aa=pd.to_numeric(a,errors='coerce');bb=pd.to_numeric(b,errors='coerce')
        numeric=bool((a.isna()|aa.notna()).all() and (b.isna()|bb.notna()).all())
        if numeric:
            av=aa.to_numpy(dtype=float,na_value=np.nan);bv=bb.to_numpy(dtype=float,na_value=np.nan)
            if not np.isclose(av,bv,rtol=1e-10,atol=1e-10,equal_nan=True).all():raise ValueError(label+' numeric mismatch '+col)
            finite=np.isfinite(av)&np.isfinite(bv)
            if finite.any():maxdiff=max(maxdiff,float(np.max(np.abs(av[finite]-bv[finite]))))
        else:
            av=a.fillna('<NA>').astype(str).to_numpy();bv=b.fillna('<NA>').astype(str).to_numpy()
            if not np.array_equal(av,bv):raise ValueError(label+' text mismatch '+col)
        cells+=len(a)
    return {'label':label,'rows':len(expected),'compared_cells':cells,'maximum_numeric_absolute_difference':maxdiff,'status':'pass'}

def replay(layout,out):
    out.mkdir(parents=True,exist_ok=False);checks=[];matrices={}
    for cohort,record in layout['cohorts'].items():
        data=read_frame(ROOT/record['cross_predictions']);key=GROUP+['fold_id','record_id']
        if data.duplicated(key).any():raise ValueError('Duplicate scoped prediction key '+cohort)
        if len(data)!=record['cross_prediction_rows']:raise ValueError('Prediction representation count mismatch')
        tables=score_tables(data);folder=out/'scores'/cohort;folder.mkdir(parents=True)
        for name,frame in tables.items():frame.to_csv(folder/(name+'.csv'),index=False,float_format='%.17g')
        for name,expected in record['expected_scores'].items():checks.append(compare_frames(tables[name],pd.read_csv(ROOT/expected,float_precision='round_trip'),cohort+'/'+name))
        matrices[cohort]=str((folder/'macro_summary.csv').relative_to(out))
    for component in ['decomposition','input_support']:
        p=ROOT/layout['portable_protocols'][component];adapter=ROOT/layout['portable_adapters'][component]
        mod=module(adapter,'portable_'+component)
        if component=='decomposition':mod.run(str(p),sha(p),str(out/component))
        else:mod.run(str(p),sha(p),str(out/component))
        for expected in layout['expected_supplements'][component]:
            name=Path(expected).name;actual=pd.read_csv(out/component/name,float_precision='round_trip');ref=pd.read_csv(ROOT/expected,float_precision='round_trip')
            checks.append(compare_frames(actual,ref,component+'/'+name))
    native=list(NATIVE_MAP_CHECKS.values())
    return {'status':'pass','scores_compared':checks,'full_fit_by_score_reference_matrices':matrices,'Figure2_numeric_source':'scores/principal/macro_summary.csv','native_map_preflight':{'unique_references':len(native),'performed':sum(x['status']=='verified_byte_identity_only' for x in native),'not_performed':sum(x['status'].startswith('not_performed') for x in native),'default':'off; archived selected-cell weights are supplied, native rasters are not'},'figures':'Original plotting code and archived figure references are supplied; this default replay regenerates numeric tables, not all manuscript figures.','new_scientific_fits':0,'new_predictions':0,'raw_acquisition_reproduced':False}

_ENGINE=None;_FRAMES={};_WORKER_GUARD_PROBED=False
def refit_worker_init(layout,output_root):
    global _WORKER_GUARD_PROBED
    os.chdir(ROOT)
    install_source_tree_guard(layout['forbidden_original_source_roots'],layout['original_analysis_root'],[str(ROOT),str(output_root)])
    try:
        with open(layout['original_analysis_root']+'/manuscript_integration_v0_1/factorial_experiment/plan01/protocol.json','rb'):
            raise RuntimeError('Worker original-source guard positive control unexpectedly opened')
    except PermissionError:_WORKER_GUARD_PROBED=True
def worker_guard_report(_):
    import multiprocessing
    if not _WORKER_GUARD_PROBED:raise RuntimeError('Worker guard initializer did not run')
    return {'start_method':multiprocessing.get_start_method(),'guard_initializer_ran':True,'original_source_positive_control_open_blocked':True}
def spawned_guard_test(layout,out):
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor
    out.mkdir(parents=True,exist_ok=False)
    with ProcessPoolExecutor(max_workers=1,mp_context=multiprocessing.get_context('spawn'),initializer=refit_worker_init,initargs=(layout,str(out))) as pool:
        result=list(pool.map(worker_guard_report,[0]))[0]
    if result['start_method']!='spawn':raise RuntimeError('Fresh spawn was not tested')
    return {'status':'pass','fresh_worker':result,'fits':0}
def fit_worker(job):
    global _ENGINE
    if not _WORKER_GUARD_PROBED:raise RuntimeError('Refit worker has no verified source-tree guard')
    if _ENGINE is None:_ENGINE=module(ROOT/job['engine'],'immutable_original_engine')
    path=job['snapshot']
    if path not in _FRAMES:_FRAMES[path]=read_frame(ROOT/path)
    _ENGINE._FRAME=_FRAMES[path]
    protocol=dict(job['protocol']);protocol['_protocol_sha256']=job['original_protocol_sha256']
    values,fit=_ENGINE.fit_task((job['fold'],job['model'],'GPP',protocol,job['output_root']))
    values['arm']=job['arm'];values['fit_reference']=job['reference'];values['score_reference']=job['reference']
    expected=job.get('expected_model_sha256')
    return values,{'cohort':job['cohort'],'arm':job['arm'],'reference':job['reference'],**fit,'expected_source_model_sha256':expected,'model_byte_hash_matches_source':fit['model_sha256']==expected if expected else None,'worker_guard_positive_control_blocked':_WORKER_GUARD_PROBED}

def refit(layout,out,workers):
    from concurrent.futures import ProcessPoolExecutor
    out.mkdir(parents=True,exist_ok=False);jobs=[]
    for cohort,info in layout['cohorts'].items():
        protocol=load(ROOT/info['protocol']);folds=load(ROOT/info['folds']);folds=folds if isinstance(folds,list) else [folds]
        folds=[f for f in folds if f.get('supported',True)]
        for s in info['snapshots']:
            for fold in folds:
                for model in protocol['models']:
                    suffix='models/'+fold['axis']+'/'+fold['fold_id']+'/GPP/'+model+'.joblib'
                    candidates=[x for x in layout['expected_model_hashes'][cohort] if x['key'].endswith(suffix) and s['arm'] in x['key'] and s['reference'] in x['key']]
                    jobs.append(dict(cohort=cohort,engine=layout['original_engine'],snapshot=s['path'],arm=s['arm'],reference=s['reference'],protocol=protocol,original_protocol_sha256=info['original_protocol_sha256'],fold=fold,model=model,output_root=str(out/'models'/cohort/(s['arm']+'__'+s['reference'])),expected_model_sha256=candidates[0]['sha256'] if len(candidates)==1 else None))
    if len(jobs)!=220:raise ValueError('Expected200principal+20T1reproduction tasks')
    outputs=[];fits=[]
    if workers==1:
        refit_worker_init(layout,str(out));results=map(fit_worker,jobs)
    else:
        pool=ProcessPoolExecutor(max_workers=workers,initializer=refit_worker_init,initargs=(layout,str(out)));results=pool.map(fit_worker,jobs)
    for values,fit in results:outputs.append((fit['cohort'],values));fits.append(fit)
    if workers!=1:pool.shutdown()
    comparisons=[]
    for cohort,info in layout['cohorts'].items():
        actual=pd.concat([v for c,v in outputs if c==cohort],ignore_index=True)
        expected=read_frame(ROOT/info['cross_predictions']);expected=expected.loc[expected.fit_reference.eq(expected.score_reference)].copy()
        fields=GROUP+['fold_id','record_id','observation_key','date','site_id','system_id','spatial_cluster','harvest_year','spatial_split_group','season_split_group','y_true','y_pred']
        # CSV dates are normalized only for identity comparison; no timezone conversion.
        actual['date']=pd.to_datetime(actual.date).dt.strftime('%Y-%m-%d');expected['date']=pd.to_datetime(expected.date).dt.strftime('%Y-%m-%d')
        comparisons.append(compare_frames(actual[fields],expected[fields],cohort+'/refit_own_predictions'))
        actual.to_csv(out/(cohort+'_reproduction_predictions.csv.gz'),index=False,float_format='%.17g',compression={'method':'gzip','mtime':0})
    save(out/'FIT_TASKS.json',fits)
    return {'status':'pass','reproducibility_refit_tasks':len(fits),'principal_tasks':sum(f['cohort']=='principal' for f in fits),'T1_tasks':sum(f['cohort']=='T1' for f in fits),'prediction_comparisons':comparisons,'model_byte_hash_matches_source':sum(f['model_byte_hash_matches_source'] is True for f in fits),'model_byte_hash_comparisons_available':sum(f['expected_source_model_sha256'] is not None for f in fits),'tasks_with_verified_worker_source_guard':sum(f['worker_guard_positive_control_blocked'] for f in fits),'workers':workers,'estimator_threads':1,'new_scientific_experiments':0,'new_scientific_eligibility':0,'note':'These are reproduction refits of the original fixed tasks, not new experiments or new independent cases. Numeric comparison tolerance1e-10; model byte hashes reported only when actually compared.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['validate','replay','refit','all','guardtest'],default='replay');p.add_argument('--out',type=Path);p.add_argument('--workers',type=int,default=4);a=p.parse_args()
    os.chdir(ROOT);layout=load(ROOT/'ARTIFACT_MAP.json');install_source_tree_guard(layout['forbidden_original_source_roots'],layout['original_analysis_root'],[str(ROOT),str((a.out or ROOT/'outputs').resolve())])
    count=validate();print(json.dumps({'manifest':'pass','immutable_files':count}),flush=True)
    if a.mode=='validate':return
    try:
        with open(layout['original_analysis_root']+'/manuscript_integration_v0_1/factorial_experiment/plan01/protocol.json','rb') as blocked_probe:
            raise RuntimeError('Original source-tree guard did not block its positive-control probe')
    except PermissionError:guard_probe_blocked=True
    out=(a.out or ROOT/'outputs'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')).resolve()
    if out==ROOT or ROOT in out.parents:raise ValueError('Output must be new and outside immutable public component')
    out.mkdir(parents=True,exist_ok=False)
    receipt={'status':'pass','created_utc':datetime.now(timezone.utc).isoformat(),'package_root':str(ROOT),'mode':a.mode,'manifest_files_verified':count,'original_source_tree_reads':'blocked_by_Python_audit_hook','original_source_tree_positive_control_probe_blocked':guard_probe_blocked,'runtime':{'python':sys.version.split()[0],'numpy':np.__version__,'pandas':pd.__version__},'raw_Flux_MODIS_AgERA5_acquisition_reproduced':False}
    if a.mode in ['replay','all']:receipt['replay']=replay(layout,out/'replay')
    if a.mode in ['refit','all']:receipt['refit']=refit(layout,out/'refit',a.workers)
    if a.mode=='guardtest':receipt['spawned_worker_guard_test']=spawned_guard_test(layout,out/'spawned_guard_test')
    save(out/'REPRODUCTION_RECEIPT.json',receipt);print(json.dumps({'status':'pass','receipt':str(out/'REPRODUCTION_RECEIPT.json')}),flush=True)
if __name__=='__main__':
    try:main()
    except Exception as exc:
        print(json.dumps({'status':'FAIL','error_type':type(exc).__name__,'error':str(exc)}),file=sys.stderr,flush=True);raise
