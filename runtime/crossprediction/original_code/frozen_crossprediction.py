"""Frozen saved-model train-arm by test-input crossprediction; never fits models."""
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ProcessPoolExecutor
import argparse, hashlib, json, os, platform, sys, time
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('MKL_NUM_THREADS','1')
import joblib
import numpy as np
import pandas as pd
import sklearn
from threadpoolctl import threadpool_limits

FEATURES=['NDVI_modis','EA_kPa_agera5','WS_ms_agera5','P_mm_agera5','TA_max_C_agera5','TA_min_C_agera5','SW_IN_MJ_m2_d_agera5','VPD_kPa_agera5']
KEYS=['record_id','observation_key','date','site_id','system_id','spatial_cluster','harvest_year','spatial_split_group','season_split_group']
GROUP=['cohort','axis','model','fit_reference','score_reference','train_arm','test_arm']
PAIR=['cohort','axis','model','fit_reference','score_reference']
TARGET='gpp_gC_m2_d'
_FRAMES={}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ids_hash(x):return hashlib.sha256('\n'.join(map(str,x)).encode()).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def now():return datetime.now(timezone.utc).isoformat()
def assert_same(a,b,label):
    pd.testing.assert_frame_equal(a.reset_index(drop=True),b.reset_index(drop=True),check_exact=True,obj=label)
def calendar(dates):
    d=pd.to_datetime(dates).dt.dayofyear.to_numpy(dtype=float)
    return np.column_stack([f(2*np.pi*k*d/365.25) for k in (1,2) for f in (np.sin,np.cos)])
def read_frame(path):
    if path not in _FRAMES:_FRAMES[path]=pd.read_parquet(path)
    return _FRAMES[path]
def worker_init():
    os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
def worker(job):
    model=joblib.load(job['model_path']);fold=job['fold'];base=read_frame(job['snapshots'][job['train_arm']])
    tr=base.iloc[fold['train_indices']];tb=read_frame(job['snapshots']['B_uniform_pixels']).iloc[fold['test_indices']]
    tc=read_frame(job['snapshots']['C_core_weighted']).iloc[fold['test_indices']]
    for key in ['axis','fold_id','train_record_ids_sha256','test_record_ids_sha256']:
        assert model[key]==fold[key],(job['model_path'],key)
    assert model['model']==job['model'] and model['target']=='GPP' and model['seed']==42
    assert model['protocol_sha256']==job['original_protocol_sha256']
    assert model['n_train']==len(tr) and model['n_test']==len(tb)
    assert ids_hash(tr.record_id)==model['train_record_ids_sha256']
    assert ids_hash(tb.record_id)==model['test_record_ids_sha256']
    assert not set(tr.record_id)&set(tb.record_id)
    assert_same(tb[KEYS+FEATURES[1:]+[TARGET]],tc[KEYS+FEATURES[1:]+[TARGET]],'paired test identity/weather/reference')
    outputs=[]
    for test_arm,frame in [('B_uniform_pixels',tb),('C_core_weighted',tc)]:
        with threadpool_limits(limits=1):
            if job['model']=='train_mean':
                assert model['feature_columns']==[]
                pred=np.full(len(frame),model['training_mean'])
            elif job['model']=='seasonal_harmonic':
                assert model['feature_columns']==['sin_doy_1','cos_doy_1','sin_doy_2','cos_doy_2']
                pred=model['estimator'].predict(calendar(frame.date))
            else:
                assert model['feature_columns']==FEATURES
                pred=model['estimator'].predict(frame[FEATURES])
        assert np.isfinite(pred).all()
        out=frame[KEYS].copy()
        for k in ['cohort','train_arm','fit_reference','model']:out[k]=job[k]
        out['test_arm']=test_arm;out['axis']=fold['axis'];out['fold_id']=fold['fold_id'];out['target']='GPP';out['seed']=42
        out['y_pred']=pred;outputs.append(out)
    ro=tb.site_id.eq('US-Ro1').to_numpy()
    if ro.any():
        assert_same(tb.loc[ro,FEATURES],tc.loc[ro,FEATURES],'Ro1 paired test features')
        assert np.array_equal(outputs[0].y_pred.to_numpy()[ro],outputs[1].y_pred.to_numpy()[ro]),'Ro1 per-fixed-model test-input prediction equality'
    return pd.concat(outputs,ignore_index=True),{'cohort':job['cohort'],'axis':fold['axis'],'fold_id':fold['fold_id'],'model':job['model'],'train_arm':job['train_arm'],'fit_reference':job['fit_reference'],'model_path':job['model_path'],'model_sha256':job['model_sha256'],'prediction_passes':2,'new_fit_tasks':0,'n_train':len(tr),'n_test':len(tb),'ro1_pair_cases':int(ro.sum())}

def metric(y,p):
    y=np.asarray(y,dtype=float);p=np.asarray(p,dtype=float)
    assert len(y) and y.shape==p.shape and np.isfinite(y).all() and np.isfinite(p).all()
    e=p-y;sse=float(np.dot(e,e));sst=float(np.sum((y-y.mean())**2))
    nse=1-sse/sst if len(y)>=3 and np.ptp(y)>0 and sst>0 else None
    r=float(np.corrcoef(y,p)[0,1]) if len(y)>=3 and np.ptp(y)>0 and np.ptp(p)>0 else None
    return dict(n=len(y),rmse=float(np.sqrt(np.mean(e**2))),mae=float(np.mean(np.abs(e))),bias_pred_minus_reference=float(np.mean(e)),nse=nse,pearson_r=r,squared_pearson_r=r*r if r is not None else None)

def numeric_tables(pred):
    native=[]
    for scope,extra in [('per_fold',['fold_id']),('per_site_season',['site_id','spatial_cluster','harvest_year']),('per_spatial_cluster',['spatial_cluster']),('pooled_days',[])]:
        keys=GROUP+extra
        for values,g in pred.groupby(keys,sort=True):
            row=dict(zip(keys,values));row['scope']=scope
            row['unit']=':'.join(str(row[k]) for k in extra) if extra else 'ALL'
            row.update(metric(g.y_true,g.y_pred));native.append(row)
    native=pd.DataFrame(native)
    cells=[]
    for scope,source in [('equal_site_season','per_site_season'),('equal_region_pooled','per_spatial_cluster'),('pooled_days','pooled_days')]:
        for values,g in native.loc[native.scope.eq(source)].groupby(GROUP,sort=True):
            row=dict(zip(GROUP,values));row.update(estimand=scope,n_units=len(g),n_days=int(g.n.sum()),rmse=float(g.rmse.mean()),mae=float(g.mae.mean()),bias_pred_minus_reference=float(g.bias_pred_minus_reference.mean()))
            cells.append(row)
    cells=pd.DataFrame(cells)
    def add_contrasts(frame,units):
        rows=[]
        for values,g in frame.groupby(PAIR+units,sort=True):
            row=dict(zip(PAIR+units,values))
            assert len(g)==4
            get=lambda t,e:float(g.loc[g.train_arm.eq(t)&g.test_arm.eq(e),'rmse'].iloc[0])
            bb=get('B_uniform_pixels','B_uniform_pixels');bc=get('B_uniform_pixels','C_core_weighted');cb=get('C_core_weighted','B_uniform_pixels');cc=get('C_core_weighted','C_core_weighted')
            row.update(R_BB=bb,R_BC=bc,R_CB=cb,R_CC=cc,
                training_given_test_B=cb-bb,training_given_test_C=cc-bc,
                test_given_training_B=bc-bb,test_given_training_C=cc-cb,
                interaction=cc-cb-bc+bb,original_pipeline_contrast=cc-bb)
            row['n_units']=int(g.n_units.iloc[0]) if 'n_units' in g else 1
            row['n_days']=int(g.n_days.iloc[0]) if 'n_days' in g else int(g.n.iloc[0])
            identity=[row['interaction']-(row['training_given_test_C']-row['training_given_test_B']),row['interaction']-(row['test_given_training_C']-row['test_given_training_B']),row['original_pipeline_contrast']-row['training_given_test_B']-row['test_given_training_C'],row['original_pipeline_contrast']-row['test_given_training_B']-row['training_given_test_C']]
            assert max(map(abs,identity))<1e-10
            row['maximum_identity_residual']=max(map(abs,identity))
            for field in ['training_given_test_B','training_given_test_C','test_given_training_B','test_given_training_C','interaction','original_pipeline_contrast']:
                row[field+'_common_BB_pp']=100*row[field]/bb if bb else 0.
            rows.append(row)
        return pd.DataFrame(rows)
    units=add_contrasts(native,['scope','unit'])
    contrasts=add_contrasts(cells,['estimand'])
    # Linear unit averaging must reproduce the explicitly labelled cell estimands.
    for _,r in contrasts.iterrows():
        scope={'equal_site_season':'per_site_season','equal_region_pooled':'per_spatial_cluster','pooled_days':'pooled_days'}[r.estimand]
        subset=units.loc[units.scope.eq(scope)]
        for k in PAIR:subset=subset.loc[subset[k].eq(r[k])]
        for field in ['training_given_test_B','training_given_test_C','test_given_training_B','test_given_training_C','interaction','original_pipeline_contrast']:
            assert abs(float(subset[field].mean())-r[field])<1e-10
    is_primary=lambda d:(d.cohort.eq('principal')&d.axis.eq('year_loyo'))|(d.cohort.eq('T1')&d.axis.eq('temporal_hindcast'))
    compact=cells.loc[is_primary(cells)&cells.estimand.eq('equal_site_season')&cells.fit_reference.eq(cells.score_reference)].copy()
    compact_contrasts=contrasts.loc[is_primary(contrasts)&contrasts.estimand.eq('equal_site_season')&contrasts.fit_reference.eq(contrasts.score_reference)].copy()
    return {'CELL_SCORES_NATIVE_UNITS.csv':native,'CELL_SCORES_ESTIMANDS.csv':cells,'CONDITIONAL_CONTRASTS_NATIVE_UNITS.csv':units,'CONDITIONAL_CONTRASTS_ESTIMANDS.csv':contrasts,'MAIN_CELLS_OWN_PATHS.csv':compact,'MAIN_CONDITIONALS_OWN_PATHS.csv':compact_contrasts}

def run(args):
    started=time.monotonic();assert sha(args.protocol)==args.expected_protocol_sha
    protocol=load(args.protocol);assert protocol['version']=='JAG_BC_crossprediction_v0.1'
    assert sha(args.bindings)==args.expected_bindings_sha
    binding=load(args.bindings);assert len(binding['models'])==120
    assert protocol['source_authority']['original_model_manifests']==binding['model_manifest_hashes']
    for r in binding['sources']+binding['models']:
        p=Path(r['path']);assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],r['path']
    args.out.mkdir(parents=True,exist_ok=False)
    execution={'status':'frozen_before_any_new_cross_input_prediction','frozen_utc':now(),'protocol_sha256':args.expected_protocol_sha,'bindings_sha256':args.expected_bindings_sha,'code_sha256':sha(__file__),'sources_verified':len(binding['sources']),'models_verified':len(binding['models']),'historical_exposure':protocol['exposure']}
    save(args.out/'EXECUTION_FREEZE.json',execution)
    frames={};foldmaps={};jobs=[];pair_checks=[]
    for cohort,record in binding['cohorts'].items():
        folds=load(record['folds']);folds=folds if isinstance(folds,list) else [folds]
        folds=[f for f in folds if f.get('supported',True)];foldmaps[cohort]={(f['axis'],f['fold_id']):f for f in folds}
        assert len(folds)==protocol['cohorts'][cohort]['supported_fold_instances']
        for arm in protocol['arms']:
            for ref in protocol['fit_references']:
                frame=pd.read_parquet(record['snapshots'][arm][ref]);frames[(cohort,arm,ref)]=frame
                assert len(frame)==(642 if cohort=='principal' else 917)
                assert not frame.record_id.duplicated().any() and not frame.observation_key.duplicated().any()
                assert np.isfinite(frame[FEATURES+[TARGET]].to_numpy()).all()
                dates=pd.to_datetime(frame.date);assert dates.dt.tz is None
        for ref in protocol['fit_references']:
            b=frames[(cohort,'B_uniform_pixels',ref)];c=frames[(cohort,'C_core_weighted',ref)]
            assert_same(b[KEYS+FEATURES[1:]+[TARGET]],c[KEYS+FEATURES[1:]+[TARGET]],'cohort arm pair')
        for arm in protocol['arms']:
            assert_same(frames[(cohort,arm,'VUT_REF')][KEYS+FEATURES],frames[(cohort,arm,'CUT_REF')][KEYS+FEATURES],'reference-invariant features')
        for fold in folds:
            f=frames[(cohort,'B_uniform_pixels','VUT_REF')];tr=f.iloc[fold['train_indices']];te=f.iloc[fold['test_indices']]
            assert not set(tr.record_id)&set(te.record_id)
            assert ids_hash(tr.record_id)==fold['train_record_ids_sha256'] and ids_hash(te.record_id)==fold['test_record_ids_sha256']
        primary_axis='year_loyo' if cohort=='principal' else 'temporal_hindcast'
        f=frames[(cohort,'B_uniform_pixels','VUT_REF')]
        ix=[i for fold in folds if fold['axis']==primary_axis for i in fold['test_indices']]
        test=f.iloc[ix];assert not test.record_id.duplicated().any()
        observed={str(s)+':'+str(int(y)):int(n) for (s,y),n in test.groupby(['site_id','harvest_year']).size().items()}
        assert observed==protocol['cohorts'][cohort]['site_seasons']
        pair_checks.append({'cohort':cohort,'primary_cases':len(test),'site_seasons':observed,'snapshot_rows':len(f),'positional_folds_checked':len(folds)})
    for model in binding['models']:
        cohort=model['cohort'];parts=model['key'].split('/');arm,ref,_,axis,fold_id,target,filename=parts
        assert target=='GPP' and filename[:-7] in protocol['models']
        jobs.append({'cohort':cohort,'train_arm':arm,'fit_reference':ref,'model':filename[:-7],'model_path':model['path'],'model_sha256':model['sha256'],'original_protocol_sha256':binding['cohorts'][cohort]['original_protocol_sha256'],'fold':foldmaps[cohort][(axis,fold_id)],'snapshots':{a:binding['cohorts'][cohort]['snapshots'][a][ref] for a in protocol['arms']}})
    outputs=[];logs=[]
    with ProcessPoolExecutor(max_workers=4,initializer=worker_init) as pool:
        for i,(frame,log) in enumerate(pool.map(worker,jobs),1):
            outputs.append(frame);logs.append(log)
            if i%40==0:print(json.dumps({'stored_models_processed':i,'total':120}),flush=True)
    base=pd.concat(outputs,ignore_index=True)
    assert len(base)==protocol['prediction_policy']['expected_unique_prediction_roles_before_score_expansion']
    allpred=[]
    for (cohort,test_arm),g in base.groupby(['cohort','test_arm'],sort=True):
        for score_ref in protocol['score_references']:
            truth=frames[(cohort,test_arm,score_ref)][['record_id',TARGET]].rename(columns={TARGET:'y_true'})
            q=g.merge(truth,on='record_id',how='left',validate='many_to_one');q['score_reference']=score_ref
            assert q.y_true.notna().all();allpred.append(q)
    pred=pd.concat(allpred,ignore_index=True)
    key=protocol['prediction_policy']['row_key'];assert not pred.duplicated(key).any()
    pred['date']=pd.to_datetime(pred.date).dt.strftime('%Y-%m-%d')
    pred=pred.sort_values(key,kind='stable').reset_index(drop=True)
    assert len(pred)==protocol['prediction_policy']['expected_scored_rows']
    diagonal=[]
    for cohort,record in binding['cohorts'].items():
        q=pred.loc[pred.cohort.eq(cohort)];assert len(q)==protocol['cohorts'][cohort]['expected_scored_rows_all_axes']
        actual=q.loc[q.train_arm.eq(q.test_arm)].copy();actual['arm']=actual.train_arm
        old=pd.read_parquet(record['original_predictions']);old=old.loc[old.arm.isin(protocol['arms'])].copy();old['date']=pd.to_datetime(old.date).dt.strftime('%Y-%m-%d')
        compare=['axis','fold_id','arm','model','target','seed','fit_reference','score_reference','record_id']
        a=actual.set_index(compare).sort_index();b=old.set_index(compare).sort_index();assert a.index.equals(b.index)
        for field in ['y_true','y_pred']:assert np.array_equal(a[field].to_numpy(),b[field].to_numpy()),(cohort,field,float(np.max(np.abs(a[field]-b[field]))))
        for field in KEYS[1:]:assert np.array_equal(a[field].to_numpy(),b[field].to_numpy()),(cohort,field)
        diagonal.append({'cohort':cohort,'scored_diagonal_rows':len(a),'maximum_prediction_difference':0.,'exact_prediction_and_reference_equality':True})
    # Scoring-path expansion must not change predictions.
    p=pred.pivot(index=[k for k in key if k!='score_reference'],columns='score_reference',values='y_pred')
    assert np.array_equal(p['VUT_REF'],p['CUT_REF'])
    control=pred.loc[pred.model.isin(['train_mean','seasonal_harmonic'])]
    p=control.pivot(index=[k for k in key if k not in ['train_arm','test_arm']],columns=['train_arm','test_arm'],values='y_pred')
    assert np.all(p.to_numpy()==p.iloc[:,[0]].to_numpy()),'Optical-free four-cell equality'
    ro=pred.loc[pred.site_id.eq('US-Ro1')]
    p=ro.pivot(index=[k for k in key if k!='test_arm'],columns='test_arm',values='y_pred')
    assert np.array_equal(p['B_uniform_pixels'],p['C_core_weighted'])
    pred.to_parquet(args.out/'cross_predictions.parquet',index=False)
    pred.to_csv(args.out/'cross_predictions.csv.gz',index=False,float_format='%.17g',compression={'method':'gzip','mtime':0})
    tables=numeric_tables(pred)
    for name,frame in tables.items():frame.to_csv(args.out/name,index=False,float_format='%.17g')
    save(args.out/'MODEL_REPLAY_LOG.json',logs)
    for r in binding['sources']+binding['models']:assert sha(r['path'])==r['sha256'],r['path']
    assert sha(args.protocol)==args.expected_protocol_sha and sha(args.bindings)==args.expected_bindings_sha and sha(__file__)==execution['code_sha256']
    receipt={'status':'PASS','created_utc':now(),'protocol_sha256':args.expected_protocol_sha,'bindings_sha256':args.expected_bindings_sha,'code_sha256':sha(__file__),'runtime':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'sklearn':sklearn.__version__,'joblib':joblib.__version__},'models_loaded':len(logs),'prediction_passes':sum(l['prediction_passes'] for l in logs),'new_fit_tasks':0,'scored_prediction_records':len(pred),'base_prediction_role_records':len(base),'diagonal_verification':diagonal,'snapshot_pairing':pair_checks,'optical_free_four_cell_prediction_equality':'PASS','fixed_model_Ro1_BC_test_prediction_equality':'PASS','scoring_reference_prediction_invariance':'PASS','source_files_unchanged_after':len(binding['sources'])+len(binding['models']),'table_rows':{k:len(v) for k,v in tables.items()},'maximum_unit_identity_residual':float(tables['CONDITIONAL_CONTRASTS_NATIVE_UNITS.csv'].maximum_identity_residual.max()),'maximum_estimand_identity_residual':float(tables['CONDITIONAL_CONTRASTS_ESTIMANDS.csv'].maximum_identity_residual.max()),'new_raw_acquisition':False,'new_FFP_runs':0,'new_scientific_eligibility':0,'independent_test_claim':False,'seconds':time.monotonic()-started}
    save(args.out/'EXECUTION_RECEIPT.json',receipt)
    files=[{'path':p.name,'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(args.out.iterdir()) if p.is_file()]
    save(args.out/'RESULT_MANIFEST.json',{'status':'PASS','files':files,'manifest_self_excluded':True})
    print(json.dumps(receipt),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--protocol',type=Path,required=True);p.add_argument('--expected-protocol-sha',required=True);p.add_argument('--bindings',type=Path,required=True);p.add_argument('--expected-bindings-sha',required=True);p.add_argument('--out',type=Path,required=True);run(p.parse_args())
