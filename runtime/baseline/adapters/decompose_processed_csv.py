"""Post-review, fixed-prediction error accounting. No estimator fitting."""
import argparse,hashlib,json,platform,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
from portable_io import read_frame

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def components(errors):
 e=np.asarray(errors,dtype=float);assert len(e) and np.isfinite(e).all()
 bias=float(e.mean());mse=float(np.mean(e*e));centered=float(np.mean((e-bias)**2))
 assert abs(mse-bias*bias-centered)<1e-10
 return dict(n=len(e),bias=bias,bias_squared=bias*bias,centered_mse=centered,mse=mse,rmse=float(np.sqrt(mse)))
def contrast(first,second):
 den=first['rmse']+second['rmse'];db=first['bias_squared']-second['bias_squared'];dc=first['centered_mse']-second['centered_mse']
 out=dict(delta_mse=first['mse']-second['mse'],delta_bias_squared=db,delta_centered_mse=dc,delta_rmse=first['rmse']-second['rmse'],bias_rmse_component=db/den if den else 0.,centered_rmse_component=dc/den if den else 0.,normalizing_denominator=den)
 assert abs(out['delta_mse']-db-dc)<1e-10
 assert abs(out['delta_rmse']-out['bias_rmse_component']-out['centered_rmse_component'])<1e-10
 return out
def boundary_checks():
 checks=[]
 def check(name,ok):assert ok,name;checks.append(name)
 q=contrast(components([1,1]),components([0,0]));check('pure_bias',q['bias_rmse_component']==1 and q['centered_rmse_component']==0)
 q=contrast(components([-2,2]),components([-1,1]));check('pure_centered_error',q['bias_rmse_component']==0 and q['centered_rmse_component']==1)
 q=contrast(components([-1,3]),components([2,2]));check('opposing_components',q['bias_rmse_component']<0<q['centered_rmse_component'] and q['delta_rmse']>0)
 q=contrast(components([0,0]),components([0,0]));check('two_zero_rmse',q['normalizing_denominator']==0 and q['bias_rmse_component']==q['centered_rmse_component']==0)
 a=contrast(components([1,1]),components([2,2]));b=contrast(components([-3,3]),components([-4,4]))
 check('unit_denominators_before_averaging',np.isclose((a['bias_rmse_component']+b['bias_rmse_component'])/2,-.5) and np.isclose((a['centered_rmse_component']+b['centered_rmse_component'])/2,-.5))
 return checks

def run(protocol_path,expected_protocol_sha,out):
 start=time.monotonic();assert sha(protocol_path)==expected_protocol_sha
 plan=load(protocol_path);assert sha(__file__)==plan['analysis_code_sha256']
 out=Path(out);out.mkdir(parents=True,exist_ok=False)
 tests=boundary_checks();metrics=[];contrasts=[];source_records=[];source_tables={};physical_counts={}
 expected_models=set(plan['models']);refs=set(plan['references'])
 for source in plan['sources']:
  path=Path(source['path']);assert sha(path)==source['sha256']
  manifest=load(source['manifest_path']);assert sha(source['manifest_path'])==source['manifest_sha256']
  assert manifest['files'][source['portable_manifest_key']]['sha256']==source['sha256']
  df=read_frame(path);assert len(df)==source['all_source_rows']
  df=df.loc[df.axis.eq(source['axis'])].copy()
  assert set(df.arm)==set(source['arms']) and set(df.model)==expected_models
  assert set(df.fit_reference)==refs and set(df.score_reference)==refs and set(df.target)=={'GPP'}
  assert np.isfinite(df[['y_true','y_pred']].to_numpy()).all()
  key=['arm','fit_reference','score_reference','model','site_id','date']
  assert not df.duplicated(key).any()
  actual_units={f'{s}:{int(y)}':int(n) for (s,y),n in df[['site_id','harvest_year','date']].drop_duplicates().groupby(['site_id','harvest_year']).size().items()}
  assert actual_units==source['site_season_cases']
  physical_counts[source['cohort']]={'cases':sum(actual_units.values()),'site_seasons':actual_units,'filtered_prediction_representations':len(df)}
  grouped=['arm','model','fit_reference','score_reference']
  baseline=None
  for keys,g in df.groupby(grouped,sort=True):
   arm,model,fit,score=keys;g=g.sort_values(['site_id','date'])
   binding=g[['site_id','date','harvest_year','spatial_cluster','fold_id']].reset_index(drop=True)
   if baseline is None:baseline=binding
   pd.testing.assert_frame_equal(binding,baseline,check_dtype=False)
   assert len(g)==sum(actual_units.values())
   for scope,columns in [('equal_site_season',['site_id','harvest_year']),('equal_region_pooled',['spatial_cluster']),('pooled_days',[])]:
    units=g.groupby(columns,sort=True) if columns else [('ALL',g)]
    for unit,u in units:
     if not isinstance(unit,tuple):unit=(unit,)
     unit=':'.join(map(str,unit))
     row=dict(cohort=source['cohort'],axis=source['axis'],scope=scope,unit=unit,arm=arm,model=model,fit_reference=fit,score_reference=score)
     row.update(components(u.y_pred.to_numpy()-u.y_true.to_numpy()));metrics.append(row)
  for model in sorted(expected_models):
   for fit in sorted(refs):
    sub=df.loc[df.model.eq(model)&df.fit_reference.eq(fit)]
    paths=[]
    for score in sorted(refs):
     paths.append(sub.loc[sub.score_reference.eq(score)].sort_values(['arm','site_id','date']).y_pred.to_numpy())
    assert np.array_equal(*paths),'Changing scoring labels must not change prediction'
  for first,second in source['contrasts']:
   k=['model','fit_reference','score_reference','site_id','date']
   f=df.loc[df.arm.eq(first)].set_index(k).sort_index();s=df.loc[df.arm.eq(second)].set_index(k).sort_index()
   assert f.index.equals(s.index) and np.array_equal(f.y_true.to_numpy(),s.y_true.to_numpy())
  # Fixed optical-free controls must still have exactly equal predictions.
  for model in ['train_mean','seasonal_harmonic']:
   piv=df.loc[df.model.eq(model)].pivot(index=['fit_reference','score_reference','site_id','date'],columns='arm',values='y_pred')
   assert np.all(piv.to_numpy()==piv.iloc[:,[0]].to_numpy())
  if 'A_existing' in source['arms']:
   k=['model','fit_reference','score_reference','site_id','date'];a=df.loc[df.arm.eq('A_existing')].set_index(k).sort_index();b=df.loc[df.arm.eq('A_shared_scene')].set_index(k).sort_index()
   assert a.index.equals(b.index) and np.array_equal(a.y_pred.to_numpy(),b.y_pred.to_numpy())
  source_records.append({'cohort':source['cohort'],'path':str(path),'sha256_before':source['sha256'],'sha256_after':sha(path),'manifest_sha256':sha(source['manifest_path'])})
  for ref in source['score_checks']:
   assert sha(ref['path'])==ref['sha256'];source_tables[(source['cohort'],ref['scope'])]=pd.read_csv(ref['path'])
 m=pd.DataFrame(metrics)
 for source in plan['sources']:
  for first,second in source['contrasts']:
   ids=['cohort','axis','scope','unit','model','fit_reference','score_reference']
   a=m.loc[m.cohort.eq(source['cohort'])&m.arm.eq(first)].set_index(ids)
   b=m.loc[m.cohort.eq(source['cohort'])&m.arm.eq(second)].set_index(ids)
   assert a.index.equals(b.index)
   for k,firstrow in a.iterrows():
    secondrow=b.loc[k];assert firstrow['n']==secondrow['n']
    row=dict(zip(ids,k));row.update(first_arm=first,second_arm=second,contrast=first+'__minus__'+second,n=int(firstrow['n']),first_rmse=firstrow['rmse'],second_rmse=secondrow['rmse'],bias_shift=firstrow['bias']-secondrow['bias'])
    row.update(contrast(firstrow,secondrow));contrasts.append(row)
 c=pd.DataFrame(contrasts);summaries=[]
 ids=['cohort','axis','scope','model','fit_reference','score_reference','contrast','first_arm','second_arm']
 for k,g in c.groupby(ids,sort=True):
  row=dict(zip(ids,k));row.update(n_units=len(g),n_days=int(g.n.sum()))
  for field in ['first_rmse','second_rmse','delta_rmse','bias_rmse_component','centered_rmse_component','delta_mse','delta_bias_squared','delta_centered_mse']:
   row['mean_'+field]=float(g[field].mean())
  row['delta_rmse_percent']=100*row['mean_delta_rmse']/row['mean_second_rmse']
  row['bias_rmse_percentage_points']=100*row['mean_bias_rmse_component']/row['mean_second_rmse']
  row['centered_rmse_percentage_points']=100*row['mean_centered_rmse_component']/row['mean_second_rmse']
  row['favorable_units']=int((g.delta_rmse<0).sum());row['unfavorable_units']=int((g.delta_rmse>0).sum());row['unchanged_units']=int((g.delta_rmse==0).sum())
  assert abs(row['mean_delta_rmse']-row['mean_bias_rmse_component']-row['mean_centered_rmse_component'])<1e-10
  summaries.append(row)
 s=pd.DataFrame(summaries)
 # Reproduce every available original primary metric, not merely contrast sums.
 original_metric_checks=0
 for (cohort,scope),table in source_tables.items():
  target=m.loc[m.cohort.eq(cohort)&m.scope.eq(scope)]
  for keys,g in target.groupby(['arm','model','fit_reference','score_reference']):
   arm,model,fit,score=keys
   t=table.loc[table.arm.eq(arm)&table.model.eq(model)&table.fit_reference.eq(fit)&table.score_reference.eq(score)]
   if 'axis' in table:t=t.loc[t.axis.eq(str(g.axis.iloc[0]))]
   assert len(t)==1,(cohort,scope,keys)
   field='macro_cluster_rmse' if scope=='equal_region_pooled' else ('rmse' if scope=='pooled_days' else 'macro_site_season_rmse')
   assert np.isclose(g.rmse.mean(),t.iloc[0][field],rtol=0,atol=1e-10)
   original_metric_checks+=1
 for name,frame in [('TableS21_error_components_by_unit.csv',m),('TableS22_paired_component_contrasts.csv',c),('TableS23_component_summary_estimands.csv',s)]:
  frame.to_csv(out/name,index=False,float_format='%.17g')
 save(out/'EXECUTION_RECEIPT.json',{'status':'pass','created_utc':datetime.now(timezone.utc).isoformat(),'protocol_sha256':expected_protocol_sha,'code_sha256':sha(__file__),'runtime':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__},'boundary_checks':tests,'physical_case_support':physical_counts,'unit_metric_rows':len(m),'paired_contrast_rows':len(c),'summary_rows':len(s),'original_metric_checks':original_metric_checks,'maximum_unit_MSE_identity_residual':float(np.max(np.abs(m.mse-m.bias_squared-m.centered_mse))),'maximum_contrast_RMSE_identity_residual':float(np.max(np.abs(c.delta_rmse-c.bias_rmse_component-c.centered_rmse_component))),'source_records':source_records,'new_model_fits':0,'new_predictions':0,'new_FFP_runs':0,'new_scientific_eligibility':0,'Cloud_calls':0,'seconds':time.monotonic()-start,'scope':'Post-review descriptive fixed-output algebra; not bias-corrected predictive skill, causal attribution, confidence interval or independent validation.'})
 print(json.dumps({'status':'pass','metrics':len(m),'contrasts':len(c),'summaries':len(s),'original_metric_checks':original_metric_checks}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--protocol',required=True);p.add_argument('--expected-protocol-sha',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(a.protocol,a.expected_protocol_sha,a.out)
