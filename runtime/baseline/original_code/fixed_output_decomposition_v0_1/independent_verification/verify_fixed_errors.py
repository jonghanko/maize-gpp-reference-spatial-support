"""Independent fixed-output accounting verification. No production imports/fits.

Run on CPU with --protocol FILE --run PRODUCTION_OUTPUT --out OWN_NEW_DIRECTORY.
Numerical reductions deliberately use math.fsum and explicit centered residuals.
"""
import argparse
import ast
import hashlib
import itertools
import json
import math
import platform
import time
from collections import Counter
from pathlib import Path
import pandas as pd

EXPECTED_PROTOCOL_SHA = 'b976269c2d1d7fb6043276faa5f980784c561772a9e00adb017f86a3c2650669'
EXPECTED_PRODUCTION_SHA = '60ae41160d8848e34b3edc0c214b0b7982670d312ba9ac6f6e0614c7cdc25fe2'
SCOPES = ['equal_site_season','equal_region_pooled','pooled_days']
METRIC_IDS = ['cohort','axis','scope','unit','arm','model','fit_reference','score_reference']
CONTRAST_IDS = ['cohort','axis','scope','unit','model','fit_reference','score_reference','contrast']
SUMMARY_IDS = ['cohort','axis','scope','model','fit_reference','score_reference','contrast']
ABS_TOL = 1e-10

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))

def dump(path,value):
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def average(values):
    values = list(values)
    return math.fsum(values)/len(values)

def residual_components(values):
    errors = [float(v) for v in values]
    if not errors or not all(math.isfinite(v) for v in errors):
        raise ValueError('Nonfinite/empty paired residuals')
    b = math.fsum(errors)/len(errors)
    centered = math.fsum((v-b)*(v-b) for v in errors)/len(errors)
    mse = math.fsum(v*v for v in errors)/len(errors)
    return dict(n=len(errors),bias=b,bias_squared=b*b,centered_mse=centered,mse=mse,rmse=math.sqrt(mse))

def pair_components(a,b):
    den = a['rmse']+b['rmse']
    db = a['bias_squared']-b['bias_squared']
    dc = a['centered_mse']-b['centered_mse']
    if den == 0:
        if any(a[k] != 0 or b[k] != 0 for k in ['mse','bias_squared','centered_mse']):
            raise ValueError('Zero RMSE with nonzero component')
        rb,rc = 0.0,0.0
    else:
        rb,rc = db/den,dc/den
    return dict(delta_mse=a['mse']-b['mse'],delta_bias_squared=db,delta_centered_mse=dc,
                delta_rmse=a['rmse']-b['rmse'],bias_rmse_component=rb,
                centered_rmse_component=rc,normalizing_denominator=den)

def run(protocol_path,production_run,out_path):
    start=time.monotonic()
    protocol_path,production_run,out_path=map(Path,(protocol_path,production_run,out_path))
    if not out_path.resolve().is_relative_to(Path(__file__).resolve().parent):
        raise ValueError('Verifier may write only below its own folder')
    if out_path.exists():
        raise FileExistsError('Preserve prior verification outputs')
    if sha(protocol_path)!=EXPECTED_PROTOCOL_SHA:
        raise ValueError('Frozen protocol mismatch')
    plan=load(protocol_path)
    production_code=protocol_path.parent/'decompose_frozen_errors.py'
    if sha(production_code)!=EXPECTED_PRODUCTION_SHA or plan['analysis_code_sha256']!=EXPECTED_PRODUCTION_SHA:
        raise ValueError('Production source mismatch')
    filenames=['TableS21_error_components_by_unit.csv','TableS22_paired_component_contrasts.csv',
               'TableS23_component_summary_estimands.csv','EXECUTION_RECEIPT.json']
    prod_pins={f:sha(production_run/f) for f in filenames}
    receipt=load(production_run/'EXECUTION_RECEIPT.json')
    out_path.mkdir(parents=True,exist_ok=False)
    counts=Counter();failures=[];worst=Counter();source_pins={};original_paths=[]
    metrics=[];contrasts=[];summaries=[];support={};original_comparisons=0

    def check(category,name,ok):
        counts[category]+=1
        if not bool(ok): failures.append({'category':category,'name':str(name)})
    def close(category,name,a,b,tol=ABS_TOL):
        a,b=float(a),float(b)
        diff=abs(a-b)
        worst[category]=max(worst[category],diff)
        check(category,name,math.isfinite(a) and math.isfinite(b) and diff<=tol)
    def pin(path,expected=None):
        path=Path(path);actual=sha(path)
        if expected is not None: check('source_hash',str(path),actual==expected)
        source_pins[str(path)]=actual
    def read_pinned_csv(source,manifest,relative):
        path=Path(source['path']).parent/relative
        expected=manifest['files'][relative]['sha256']
        pin(path,expected);original_paths.append(str(path))
        return pd.read_csv(path)

    # Independent boundary cases include the zero denominator and a macro trap.
    z=pair_components(residual_components([0,0]),residual_components([0,0]))
    check('boundary','zero denominator',z['normalizing_denominator']==0 and z['bias_rmse_component']==0 and z['centered_rmse_component']==0)
    p=pair_components(residual_components([1,1]),residual_components([0,0]))
    close('boundary','pure bias',p['bias_rmse_component'],1)
    close('boundary','pure bias centered',p['centered_rmse_component'],0)
    p=pair_components(residual_components([-2,2]),residual_components([-1,1]))
    close('boundary','pure centered',p['centered_rmse_component'],1)
    p=pair_components(residual_components([-1,3]),residual_components([2,2]))
    check('boundary','offsetting components',p['bias_rmse_component']<0<p['centered_rmse_component'] and p['delta_rmse']>0)
    a=pair_components(residual_components([1,1]),residual_components([2,2]))
    b=pair_components(residual_components([-3,3]),residual_components([-4,4]))
    close('boundary','macro uses individual denominators',average([a['bias_rmse_component'],b['bias_rmse_component']]),-.5)
    # Squared macro bias is not mean squared unit bias.
    check('boundary','squared macro bias trap',average([1.0**2,(-1.0)**2])!=average([1.0,-1.0])**2)
    check('protocol','three estimands exact',plan['aggregations']==SCOPES)
    check('protocol','no production fit call',not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr in {'fit','fit_transform','partial_fit'} for n in ast.walk(ast.parse(production_code.read_text()))))

    for source in plan['sources']:
        cohort=source['cohort'];pin(source['path'],source['sha256']);pin(source['manifest_path'],source['manifest_sha256'])
        manifest=load(source['manifest_path'])
        check('source_hash',cohort+' array manifest pin',manifest['files']['cross_predictions.parquet']['sha256']==source['sha256'])
        full=pd.read_parquet(source['path'])
        check('support',cohort+' full representations',len(full)==source['all_source_rows'])
        df=full.loc[full['axis']==source['axis']].copy()
        check('support',cohort+' only GPP',set(df['target'])=={'GPP'})
        check('support',cohort+' families',set(df['model'])==set(plan['models']))
        check('support',cohort+' arms',set(df['arm'])==set(source['arms']))
        check('support',cohort+' ref paths',set(df['fit_reference'])==set(plan['references']) and set(df['score_reference'])==set(plan['references']))
        check('support',cohort+' unique identities',not df.duplicated(['arm','model','fit_reference','score_reference','site_id','date']).any())
        check('support',cohort+' finite pairs',all(math.isfinite(float(x)) for col in ['y_true','y_pred'] for x in df[col]))
        physical=df[['site_id','date','harvest_year','spatial_cluster','fold_id']].drop_duplicates().sort_values(['site_id','date']).reset_index(drop=True)
        check('support',cohort+' one metadata binding per day',not physical.duplicated(['site_id','date']).any())
        site_counts={f'{s}:{int(y)}':int(n) for (s,y),n in physical.groupby(['site_id','harvest_year']).size().items()}
        check('support',cohort+' exact season counts',site_counts==source['site_season_cases'])
        physical_n=sum(site_counts.values())
        expected_groups=list(itertools.product(source['arms'],plan['models'],plan['references'],plan['references']))
        check('support',cohort+' retained full matrix count',len(df)==physical_n*len(expected_groups))
        support[cohort]={'physical_days':physical_n,'site_seasons':site_counts,'prediction_representations':len(df),'region_count':physical.spatial_cluster.nunique()}
        target_by_score={};pred_by_arm_model_fit={};controls={};metrics_this=[]
        for arm,model,fit,score in expected_groups:
            g=df.loc[(df.arm==arm)&(df.model==model)&(df.fit_reference==fit)&(df.score_reference==score)].sort_values(['site_id','date'])
            binding=g[['site_id','date','harvest_year','spatial_cluster','fold_id']].reset_index(drop=True)
            check('binding',(cohort,arm,model,fit,score),binding.equals(physical))
            target=g.y_true.tolist();prediction=g.y_pred.tolist()
            if score not in target_by_score: target_by_score[score]=target
            check('binding',(cohort,'one target per scoring reference',arm,model,fit,score),target==target_by_score[score])
            pkey=(arm,model,fit)
            if pkey not in pred_by_arm_model_fit:pred_by_arm_model_fit[pkey]=prediction
            check('binding',(cohort,'scoring path does not alter prediction',pkey,score),prediction==pred_by_arm_model_fit[pkey])
            if model in ['train_mean','seasonal_harmonic']:
                ck=(model,fit,score)
                if ck not in controls:controls[ck]=prediction
                check('controls',(cohort,arm,ck),prediction==controls[ck])
            for scope in SCOPES:
                if scope=='equal_site_season': units=[(f'{s}:{int(y)}',u) for (s,y),u in g.groupby(['site_id','harvest_year'],sort=True)]
                elif scope=='equal_region_pooled':units=[(str(region),u) for region,u in g.groupby('spatial_cluster',sort=True)]
                else:units=[('ALL',g)]
                for unit,u in units:
                    row=dict(cohort=cohort,axis=source['axis'],scope=scope,unit=unit,arm=arm,model=model,fit_reference=fit,score_reference=score)
                    row.update(residual_components(float(p)-float(y) for p,y in zip(u.y_pred,u.y_true)))
                    check('identity',(cohort,scope,unit,arm,model,fit,score,'nonnegative centered'),row['centered_mse']>=0)
                    close('identity',(cohort,scope,unit,arm,model,fit,score,'MSE sum'),row['mse'],row['bias_squared']+row['centered_mse'])
                    metrics.append(row);metrics_this.append(row)
        if 'A_existing' in source['arms']:
            for model,fit in itertools.product(plan['models'],plan['references']):
                check('controls',(cohort,'point anchor identity',model,fit),pred_by_arm_model_fit[('A_existing',model,fit)]==pred_by_arm_model_fit[('A_shared_scene',model,fit)])

        # Original per-unit metrics are independently pinned by the original run manifest.
        relative_by_scope={'equal_site_season':'summary/cross_scores_site_season.csv','equal_region_pooled':'summary/cross_scores_cluster.csv','pooled_days':'summary/cross_scores_pooled.csv'} if cohort=='principal_642' else {'equal_site_season':'scores_site_season.csv','equal_region_pooled':'scores_pooled.csv','pooled_days':'scores_pooled.csv'}
        for scope,relative in relative_by_scope.items():
            table=read_pinned_csv(source,manifest,relative)
            for row in [r for r in metrics_this if r['scope']==scope]:
                t=table.loc[(table.arm==row['arm'])&(table.model==row['model'])&(table.fit_reference==row['fit_reference'])&(table.score_reference==row['score_reference'])]
                if 'axis' in table:t=t.loc[t.axis==source['axis']]
                if scope=='equal_site_season':
                    site,year=row['unit'].split(':');t=t.loc[(t.site_id==site)&(t.harvest_year==int(year))]
                elif scope=='equal_region_pooled' and 'spatial_cluster' in table:t=t.loc[t.spatial_cluster==row['unit']]
                check('original_metrics',(cohort,scope,row['unit'],row['arm'],row['model'],row['fit_reference'],row['score_reference'],'unique row'),len(t)==1)
                if len(t)!=1:continue
                for ours,original in [('rmse','rmse'),('bias','bias_pred_minus_obs'),('n','n')]:
                    close('original_metrics',(cohort,scope,row['unit'],row['arm'],row['model'],row['fit_reference'],row['score_reference'],ours),row[ours],t.iloc[0][original])
                    original_comparisons+=1
        for ref in source['score_checks']:
            pin(ref['path'],ref['sha256']);table=pd.read_csv(ref['path'])
            for arm,model,fit,score in expected_groups:
                rows=[r for r in metrics_this if r['scope']==ref['scope'] and (r['arm'],r['model'],r['fit_reference'],r['score_reference'])==(arm,model,fit,score)]
                t=table.loc[(table.arm==arm)&(table.model==model)&(table.fit_reference==fit)&(table.score_reference==score)]
                if 'axis' in table:t=t.loc[t.axis==source['axis']]
                check('original_macros',(cohort,ref['scope'],arm,model,fit,score,'unique'),len(t)==1)
                if len(t)==1:
                    field={'equal_region_pooled':'macro_cluster_rmse','equal_site_season':'macro_site_season_rmse','pooled_days':'rmse'}[ref['scope']]
                    close('original_macros',(cohort,ref['scope'],arm,model,fit,score),average(r['rmse'] for r in rows),t.iloc[0][field])

    mlookup={tuple(row[k] for k in METRIC_IDS):row for row in metrics}
    check('output_keys','independent metric keys unique',len(mlookup)==len(metrics))
    for source in plan['sources']:
        for first,second in source['contrasts']:
            for row in [r for r in metrics if r['cohort']==source['cohort'] and r['arm']==first]:
                key=[row[k] for k in METRIC_IDS];key[4]=second;b=mlookup[tuple(key)]
                out={k:row[k] for k in ['cohort','axis','scope','unit','model','fit_reference','score_reference']}
                out.update(first_arm=first,second_arm=second,contrast=first+'__minus__'+second,n=row['n'],first_rmse=row['rmse'],second_rmse=b['rmse'],bias_shift=row['bias']-b['bias'])
                out.update(pair_components(row,b));contrasts.append(out)
                close('identity',tuple(out[k] for k in CONTRAST_IDS)+('contrast MSE',),out['delta_mse'],out['delta_bias_squared']+out['delta_centered_mse'])
                close('identity',tuple(out[k] for k in CONTRAST_IDS)+('contrast RMSE',),out['delta_rmse'],out['bias_rmse_component']+out['centered_rmse_component'])
    groups={}
    for row in contrasts:groups.setdefault(tuple(row[k] for k in SUMMARY_IDS),[]).append(row)
    for key,rows in groups.items():
        out=dict(zip(SUMMARY_IDS,key));out.update(first_arm=rows[0]['first_arm'],second_arm=rows[0]['second_arm'],n_units=len(rows),n_days=sum(r['n'] for r in rows))
        for field in ['first_rmse','second_rmse','delta_rmse','bias_rmse_component','centered_rmse_component','delta_mse','delta_bias_squared','delta_centered_mse']:
            out['mean_'+field]=average(r[field] for r in rows)
        den=out['mean_second_rmse']
        if den<=0:raise ValueError('Actual zero-comparator macro percentage requires explicit missing-value protocol')
        out['delta_rmse_percent']=100*out['mean_delta_rmse']/den
        out['bias_rmse_percentage_points']=100*out['mean_bias_rmse_component']/den
        out['centered_rmse_percentage_points']=100*out['mean_centered_rmse_component']/den
        out['favorable_units']=sum(r['delta_rmse']<0 for r in rows)
        out['unfavorable_units']=sum(r['delta_rmse']>0 for r in rows)
        out['unchanged_units']=sum(r['delta_rmse']==0 for r in rows)
        close('identity',key+('macro RMSE',),out['mean_delta_rmse'],out['mean_bias_rmse_component']+out['mean_centered_rmse_component'])
        close('identity',key+('macro MSE',),out['mean_delta_mse'],out['mean_delta_bias_squared']+out['mean_delta_centered_mse'])
        close('identity',key+('percentage points',),out['delta_rmse_percent'],out['bias_rmse_percentage_points']+out['centered_rmse_percentage_points'])
        summaries.append(out)

    def compare_output(filename,expected,ids):
        table=pd.read_csv(production_run/filename)
        check('output_keys',filename+' no duplicate keys',not table.duplicated(ids).any())
        actual={tuple(row[k] for k in ids):row for row in table.to_dict('records')}
        wanted={tuple(row[k] for k in ids):row for row in expected}
        check('output_keys',filename+' exact complete keys',set(actual)==set(wanted))
        check('output_keys',filename+' exact rows',len(table)==len(expected))
        for key,row in wanted.items():
            if key not in actual:continue
            for field,value in row.items():
                if field in ids:continue
                if isinstance(value,(int,float)):
                    close('production_numeric',key+(field,),actual[key][field],value,tol=1e-9 if 'percent' in field else ABS_TOL)
                else:check('production_labels',key+(field,),actual[key][field]==value)
    compare_output(filenames[0],metrics,METRIC_IDS)
    compare_output(filenames[1],contrasts,CONTRAST_IDS)
    compare_output(filenames[2],summaries,SUMMARY_IDS)
    check('receipt','frozen protocol',receipt['protocol_sha256']==EXPECTED_PROTOCOL_SHA)
    check('receipt','frozen production code',receipt['code_sha256']==EXPECTED_PRODUCTION_SHA)
    check('receipt','counts',receipt['unit_metric_rows']==len(metrics) and receipt['paired_contrast_rows']==len(contrasts) and receipt['summary_rows']==len(summaries))
    check('receipt','source physical support',all(receipt['physical_case_support'][k]['cases']==v['physical_days'] and receipt['physical_case_support'][k]['site_seasons']==v['site_seasons'] for k,v in support.items()))
    for field in ['new_model_fits','new_predictions','new_FFP_runs','new_scientific_eligibility','Cloud_calls']:
        check('receipt',field,receipt[field]==0)
    for path,expected in source_pins.items():check('unchanged_sources',path,sha(path)==expected)
    check('unchanged_sources','protocol',sha(protocol_path)==EXPECTED_PROTOCOL_SHA)
    check('unchanged_sources','production source',sha(production_code)==EXPECTED_PRODUCTION_SHA)
    for filename,pinned in prod_pins.items():check('unchanged_outputs',filename,sha(production_run/filename)==pinned)

    for name,rows in [('INDEPENDENT_METRICS.csv',metrics),('INDEPENDENT_CONTRASTS.csv',contrasts),('INDEPENDENT_SUMMARIES.csv',summaries)]:
        pd.DataFrame(rows).to_csv(out_path/name,index=False,float_format='%.17g')
    result={'status':'PASS' if not failures else 'FAIL','protocol_sha256':EXPECTED_PROTOCOL_SHA,'production_code_sha256':EXPECTED_PRODUCTION_SHA,
            'verifier_code_sha256':sha(__file__),'production_run':str(production_run),'checks':sum(counts.values()),'checks_by_category':dict(counts),
            'failed_count':len(failures),'failures':failures,'maximum_absolute_difference_by_category':dict(worst),'support':support,
            'independent_rows':{'metrics':len(metrics),'contrasts':len(contrasts),'summaries':len(summaries)},
            'original_per_unit_metric_comparisons':original_comparisons,'source_pins':source_pins,'production_output_pins':prod_pins,
            'runtime':{'python':platform.python_version(),'pandas':pd.__version__,'seconds':time.monotonic()-start},
            'numerics':'math.fsum of explicit squared residuals and explicitly centered squared residuals, ddof=0; no production-function imports',
            'limits':'Post-review algebraic accounting; no bias correction, causal attribution, confidence interval, field approval or independent confirmation.',
            'model_fits':0,'new_predictions':0,'new_cases':0,'external_messages_sent':0,'Cloud_calls':0,
            'model_binary_hash_scope':'Original run manifests verified unchanged; individual serialized model binaries not reread because only frozen predictions are used.',
            'preflight_note':'Initial pre-freeze read had incomplete pooled_days original-metric field map; corrected final frozen code confirmed before CPU execution.'}
    dump(out_path/'INDEPENDENT_VERIFICATION.json',result)
    text=f"""# Independent fixed-output decomposition verification

**{result['status']}** — {result['checks']} checks, {len(failures)} failures.

Recomputed {len(metrics)} unit metrics, {len(contrasts)} paired contrasts and {len(summaries)} summaries from the original pinned prediction arrays using `math.fsum` and explicit centered residual sums. All case/arm/model/reference bindings, optical-free controls, original per-unit metrics, original macro metrics, exact unit and macro RMSE contributions, separate aggregation scopes, and source/output hashes were checked. Complete details and source pins are in `INDEPENDENT_VERIFICATION.json`.

The principal cohort retains 642 physical days and the separate hindcast 437; repeated reference/arm rows are not new observations. Equal-site-season RMSE, region-pooled RMSE and all-days RMSE remain distinct estimands. The zero-denominator accounting rule was tested independently. Source arrays/manifests and production outputs remained unchanged throughout verification. No model was fitted and no new prediction or case generated.

This verifies post-review numerical accounting only. It does not resolve field dates, mast registration, displacement evidence, independent regional support or historical exposure; it does not establish bias-corrected predictive skill or physical causes.
"""
    (out_path/'README_EN.md').write_text(text,encoding='utf-8')
    print(json.dumps({'status':result['status'],'checks':result['checks'],'failed':len(failures),'rows':result['independent_rows'],'maximum_difference':dict(worst)}))
    if failures:raise SystemExit(1)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--protocol',required=True);parser.add_argument('--run',required=True);parser.add_argument('--out',required=True)
    args=parser.parse_args();run(args.protocol,args.run,args.out)
