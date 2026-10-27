"""Independent target-free input/anchor audit. No production imports/model reads."""
import argparse
import bisect
import hashlib
import json
import math
import statistics
import time
from collections import Counter,defaultdict
from datetime import date
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq

PROTOCOL_SHA='839a20c09e8940e084a3dc555c6afbab1c53ce695333e310991d66d995c1f986'
CODE_SHA='5e5222e2ea6199342e45599a4d5aee991fc98eb42b652b2a2b648ba1b4ad8591'
FEATURES=['NDVI_modis','EA_kPa_agera5','WS_ms_agera5','P_mm_agera5','TA_max_C_agera5','TA_min_C_agera5','SW_IN_MJ_m2_d_agera5','VPD_kPa_agera5']
IDS=['site_id','date','record_id','harvest_year']
TABLES={
 'monthly':('TableS24_monthly_weight_concentration.csv',['cohort','site_id','year','month']),
 'daily':('TableS25_retained_daily_input_lineage.csv',['cohort','site_id','year','date']),
 'summary':('TableS26_site_season_input_summary.csv',['cohort','site_id','year','T1_role']),
 'scenes':('TableS27_scene_anchor_usage.csv',['cohort','site_id','year','anchor'])}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def mean(x):x=list(x);return math.fsum(x)/len(x)
def day(x):return pd.Timestamp(x).date()
def date_text(x):return day(x).isoformat()
def native(x):return x.item() if hasattr(x,'item') else x
def weight_values(w):
    w=[float(x) for x in w]
    if not w or not all(math.isfinite(x) and x>0 for x in w):raise ValueError('Invalid original positive weights')
    total=math.fsum(w);p=[x/total for x in w];n=len(p)
    return dict(n_cells=n,source_conditional_weight_sum=total,concentration_cell_count=1/math.fsum(x*x for x in p),max_conditional_weight=max(p),distance_from_uniform=math.fsum(abs(x-1/n) for x in p)/2)
def bracket(anchors,current):
    i=bisect.bisect_left(anchors,current)
    if i<len(anchors) and anchors[i]==current:return i,i,0.0
    if i==0 or i==len(anchors):raise ValueError('Extrapolation would be required')
    left,right=anchors[i-1],anchors[i]
    gap=(right-left).days
    if not 0<gap<=24:raise ValueError('Unsupported interpolation bracket')
    if left.year!=current.year or right.year!=current.year:raise ValueError('Cross-year bracket')
    return i-1,i,(current-left).days/gap

def run(protocol,production_run,out):
    started=time.monotonic();protocol,production_run,out=map(Path,(protocol,production_run,out))
    if not out.resolve().is_relative_to(Path(__file__).resolve().parent):raise ValueError('Write scope must be own verifier folder')
    if out.exists():raise FileExistsError('Preserve existing verification attempt')
    if sha(protocol)!=PROTOCOL_SHA:raise ValueError('Frozen protocol hash mismatch')
    plan=read(protocol);code=protocol.parent/'audit_inputs.py'
    if sha(code)!=CODE_SHA or plan['code_sha256']!=CODE_SHA:raise ValueError('Frozen production code mismatch')
    output_pins={name:sha(production_run/name) for name,_ in TABLES.values()}
    output_pins['EXECUTION_RECEIPT.json']=sha(production_run/'EXECUTION_RECEIPT.json')
    receipt=read(production_run/'EXECUTION_RECEIPT.json')
    out.mkdir(parents=True,exist_ok=False)
    checks=Counter();failed=[];maxdiff=Counter();source_pins={};columns_read={}
    def check(category,name,condition):
        checks[category]+=1
        if not bool(condition):failed.append({'category':category,'name':str(name)})
    def close(category,name,a,b,tol=1e-12):
        a,b=float(a),float(b);difference=abs(a-b);maxdiff[category]=max(maxdiff[category],difference)
        check(category,name,math.isfinite(a) and math.isfinite(b) and difference<=tol)
    def pin(path,expected):
        path=str(path);actual=sha(path);check('source_hash',path,actual==expected);source_pins[path]=expected
    def f(key):return plan['files'][key]['path']
    def projected(key,fields):
        p=f(key);columns_read[p]=fields
        return pd.read_parquet(p,columns=fields)
    for item in list(plan['files'].values())+plan['manifest_pins']+plan['code_pins']+plan['map_pins']:pin(item['path'],item['sha256'])
    for w,nc,tv in [([1],1,0),([.5,.5],2,0),([.75,.25],1.6,.25)]:
        v=weight_values(w);close('boundary',str(w)+' concentration',v['concentration_cell_count'],nc);close('boundary',str(w)+' TV',v['distance_from_uniform'],tv)
    anchors_test=[date(2012,4,25),date(2012,5,19)]
    check('boundary','exact reuse',bracket(anchors_test,anchors_test[0])==(0,0,0.0))
    check('boundary','24day cross-month bracket',bracket(anchors_test,date(2012,5,1))[:2]==(0,1))
    for label,aa,dd in [('25day',[date(2012,4,24),date(2012,5,19)],date(2012,5,1)),('left extrapolation',anchors_test,date(2012,4,24)),('right extrapolation',anchors_test,date(2012,5,20)),('cross year',[date(2011,12,31),date(2012,1,8)],date(2012,1,2))]:
        try:bracket(aa,dd);rejected=False
        except ValueError:rejected=True
        check('boundary',label,rejected)

    months=[];month_lookup={}
    for cohort,key in [('principal_642','primary_monthly_weights'),('T1_test437','T1_monthly_weights_in_protocol')]:
        for site in read(f(key))['sites']:
            check('monthly_binding',(cohort,site['site_id'],site['year'],'all months'),{int(m['month']) for m in site['maps']}==set(range(1,13)))
            for m in site['maps']:
                pin(m['path'],m['sha256'])
                pixels=m['pixels'];cellids=[(p['row'],p['col']) for p in pixels]
                check('monthly_binding',(cohort,site['site_id'],site['year'],m['month'],'unique positive cells'),len(set(cellids))==len(cellids))
                r=dict(cohort=cohort,site_id=site['site_id'],year=int(site['year']),month=int(m['month']),map_path=m['path'],map_sha256=m['sha256'])
                r.update(weight_values(p['conditional_core_weight'] for p in pixels))
                check('monthly_values',(r['site_id'],r['year'],r['month'],'source normalized'),abs(r['source_conditional_weight_sum']-1)<1e-7)
                check('monthly_values',(r['site_id'],r['year'],r['month'],'bounds'),1-1e-12<=r['concentration_cell_count']<=r['n_cells']+1e-12 and 0<=r['distance_from_uniform']<=1)
                months.append(r);month_lookup[(r['site_id'],r['year'],r['month'])]=r
    check('monthly_binding','all84 unique maps',len(months)==len(month_lookup)==84)

    frames={}
    for cohort,prefix,n in [('principal_642','primary642',642),('T1_test437','T1',917)]:
        for arm in ['B','C']:
            for ref in ['VUT','CUT']:
                d=projected(prefix+'_'+arm+'_'+ref,IDS+FEATURES);d['date']=d['date'].map(date_text)
                check('snapshot_support',(cohort,arm,ref,'row count'),len(d)==n)
                check('snapshot_support',(cohort,arm,ref,'unique days'),not d.duplicated(['site_id','date']).any())
                check('snapshot_support',(cohort,arm,ref,'finite predictors'),all(math.isfinite(float(v)) for field in FEATURES for v in d[field]))
                frames[(cohort,arm,ref)]=d
        base=frames[(cohort,'B','VUT')]
        for arm in ['B','C']:
            for ref in ['VUT','CUT']:
                d=frames[(cohort,arm,ref)]
                check('snapshot_binding',(cohort,arm,ref,'ordered IDs'),list(d[IDS].itertuples(index=False,name=None))==list(base[IDS].itertuples(index=False,name=None)))
                check('snapshot_binding',(cohort,arm,ref,'weather exact'),list(d[FEATURES[1:]].itertuples(index=False,name=None))==list(base[FEATURES[1:]].itertuples(index=False,name=None)))
            check('snapshot_binding',(cohort,arm,'reference X exact'),list(frames[(cohort,arm,'VUT')][FEATURES].itertuples(index=False,name=None))==list(frames[(cohort,arm,'CUT')][FEATURES].itertuples(index=False,name=None)))

    fold=read(f('T1_frozen_fold'));train_indices=fold['train_indices'];test_indices=fold['test_indices']
    check('fold','480 training and437 test',fold['n_train']==len(train_indices)==480 and fold['n_test']==len(test_indices)==437)
    check('fold','full disjoint partition',not set(train_indices)&set(test_indices) and set(train_indices)|set(test_indices)==set(range(917)))
    train_keys=set(frames[('T1_test437','B','VUT')].iloc[train_indices][['site_id','date']].itertuples(index=False,name=None))
    for arm in ['B','C']:
        train=frames[('T1_test437',arm,'VUT')].iloc[train_indices].sort_values(['site_id','date'])
        old=frames[('principal_642',arm,'VUT')];old=old[old.site_id.isin(['US-Ne1','US-Ne3'])].sort_values(['site_id','date'])
        fields=['site_id','date','harvest_year']+FEATURES
        check('training_reuse',arm+' exact480 keys/predictors',len(train)==len(old)==480 and list(train[fields].itertuples(index=False,name=None))==list(old[fields].itertuples(index=False,name=None)))
    for cohort,key in [('principal_642','primary642_selected_keys'),('T1_test437','T1_test_source_keys')]:
        # Schema metadata inspected only to choose predeclared identifier columns.
        names=pq.read_schema(f(key)).names
        fields=[x for x in ['site_id','date','harvest_year','record_id'] if x in names]
        q=projected(key,fields);q['date']=q['date'].map(date_text)
        expected=frames[(cohort,'B','VUT')]
        if cohort=='T1_test437':expected=expected.iloc[test_indices]
        check('selected_keys',cohort+' exact frozen selected physical keys',set(q[['site_id','date']].itertuples(index=False,name=None))==set(expected[['site_id','date']].itertuples(index=False,name=None)) and len(q)==len(expected))

    primary_columns=['site_id','year','anchor','item_id','positive_core_pixels','shared_anchor_usable','B_uniform_pixels','C_core_weighted']
    scenes={'principal_642':projected('primary_original_accepted_scene_audit',primary_columns).to_dict('records'),'T1_test437':read(f('T1_accepted_scene_audit'))['scenes']}
    accepted={};all_scene_records=[]
    for cohort,rows in scenes.items():
        scene_keys=[]
        for raw in rows:
            r=dict(raw);r['anchor']=date_text(r['anchor']);year=int(r['year']);anchor=day(r['anchor']);site=r['site_id']
            scene_keys.append((site,year,r['anchor']))
            check('scene_binding',(cohort,site,year,r['anchor'],'year'),anchor.year==year)
            check('scene_binding',(cohort,site,year,r['anchor'],'positive cells'),int(r['positive_core_pixels'])==month_lookup[(site,year,anchor.month)]['n_cells'])
            all_scene_records.append((cohort,r))
            if r['shared_anchor_usable']:
                check('scene_binding',(cohort,site,year,r['anchor'],'finite accepted NDVI'),all(math.isfinite(float(r[x])) for x in ['B_uniform_pixels','C_core_weighted']))
                accepted.setdefault((cohort,site,year),[]).append(r)
        check('scene_binding',cohort+' unique all-year scene keys',len(scene_keys)==len(set(scene_keys)))
    for key in accepted:accepted[key].sort(key=lambda x:x['anchor'])

    daily=[];anchor_users=defaultdict(set);month_users=defaultdict(set);max_replay=0.0
    for cohort in ['principal_642','T1_test437']:
        bb=frames[(cohort,'B','VUT')];cc=frames[(cohort,'C','VUT')]
        if cohort=='T1_test437':bb=bb.iloc[test_indices];cc=cc.iloc[test_indices]
        blookup={(r.site_id,r.date):r for r in bb.itertuples(index=False)}
        clookup={(r.site_id,r.date):r for r in cc.itertuples(index=False)}
        check('daily_binding',cohort+' B/C case equality',set(blookup)==set(clookup))
        for (site,dt),r in sorted(blookup.items()):
            year=int(r.harvest_year);current=day(dt);aa=accepted[(cohort,site,year)];anchor_dates=[day(a['anchor']) for a in aa]
            il,ir,fraction=bracket(anchor_dates,current);left,right=aa[il],aa[ir];ld,rd=anchor_dates[il],anchor_dates[ir]
            bv,cv=float(r.NDVI_modis),float(clookup[(site,dt)].NDVI_modis)
            for field,actual in [('B_uniform_pixels',bv),('C_core_weighted',cv)]:
                # Alternate arithmetic form to production's convex combination.
                replay=float(left[field])+fraction*(float(right[field])-float(left[field]))
                difference=abs(replay-actual);max_replay=max(max_replay,difference)
                check('NDVI_replay',(cohort,site,dt,field),difference<1e-12)
            role='T1_test' if cohort=='T1_test437' else ('T1_training_reuse' if (site,dt) in train_keys else 'not_used_in_T1')
            row=dict(cohort=cohort,site_id=site,year=year,date=dt,T1_role=role,B_NDVI=bv,C_NDVI=cv,C_minus_B_NDVI=cv-bv,absolute_NDVI_difference=abs(cv-bv),left_anchor=ld.isoformat(),right_anchor=rd.isoformat(),left_item_id=str(left['item_id']),right_item_id=str(right['item_id']),left_map_month=ld.month,right_map_month=rd.month,bracket_days=(rd-ld).days,right_anchor_fraction=fraction,direct_anchor=il==ir,right_anchor_after_case=rd>current,cross_month_bracket=ld.month!=rd.month)
            daily.append(row);case=(cohort,site,year,dt)
            for a in {ld,rd}:anchor_users[(cohort,site,year,a.isoformat())].add(case)
            for month in {ld.month,rd.month}:month_users[(site,year,month)].add(case)
    case_counts=Counter(f"{r['site_id']}:{r['year']}" for r in daily)
    check('daily_binding','1079 exact season support',len(daily)==1079 and dict(case_counts)==plan['case_counts'])
    check('daily_binding','no physical duplicate across cohorts',len({(r['site_id'],r['year'],r['date']) for r in daily})==1079)
    check('training_reuse','480 principal rows marked reused',sum(r['T1_role']=='T1_training_reuse' for r in daily)==480)
    for cohort,key in [('principal_642','primary_original_daily_optics'),('T1_test437','T1_daily_full_year_optics')]:
        fields=['site_id','date','B_uniform_pixels','C_core_weighted'];optical=projected(key,fields);optical['date']=optical['date'].map(date_text)
        lookup={(r.site_id,r.date):r for r in optical.itertuples(index=False)}
        check('original_optics',cohort+' unique original curve keys',len(lookup)==len(optical))
        for r in [x for x in daily if x['cohort']==cohort]:
            source=lookup[(r['site_id'],r['date'])]
            close('original_optics',(cohort,r['site_id'],r['date'],'B'),r['B_NDVI'],source.B_uniform_pixels)
            close('original_optics',(cohort,r['site_id'],r['date'],'C'),r['C_NDVI'],source.C_core_weighted)
    for r in months:
        n=len(month_users[(r['site_id'],r['year'],r['month'])]);r.update(retained_dates_using_map=n,used_by_retained_case_anchor=n>0)
    usage=[]
    for cohort,r in all_scene_records:
        usage.append(dict(cohort=cohort,site_id=r['site_id'],year=int(r['year']),anchor=r['anchor'],item_id=r['item_id'],map_month=day(r['anchor']).month,shared_anchor_usable=bool(r['shared_anchor_usable']),retained_dates_using_anchor=len(anchor_users[(cohort,r['site_id'],int(r['year']),r['anchor'])]),B_anchor_NDVI=r['B_uniform_pixels'],C_anchor_NDVI=r['C_core_weighted']))
    groups=defaultdict(list)
    for r in daily:groups[(r['cohort'],r['site_id'],r['year'],r['T1_role'])].append(r)
    summaries=[]
    for (cohort,site,year,role),rows in sorted(groups.items()):
        active=[r for r in months if r['site_id']==site and r['year']==year and r['used_by_retained_case_anchor']]
        absolute=[r['absolute_NDVI_difference'] for r in rows];signed=[r['C_minus_B_NDVI'] for r in rows]
        summaries.append(dict(cohort=cohort,site_id=site,year=year,T1_role=role,n_days=len(rows),mean_C_minus_B_NDVI=mean(signed),mean_absolute_NDVI_difference=mean(absolute),median_absolute_NDVI_difference=statistics.median(absolute),max_absolute_NDVI_difference=max(absolute),n_exact_equal_BC=sum(v==0 for v in signed),n_difference_above_tolerance=sum(v>1e-12 for v in absolute),active_map_months=len(active),minimum_active_cell_count=min(r['n_cells'] for r in active),maximum_active_cell_count=max(r['n_cells'] for r in active),minimum_active_concentration=min(r['concentration_cell_count'] for r in active),maximum_active_concentration=max(r['concentration_cell_count'] for r in active),direct_anchor_days=sum(r['direct_anchor'] for r in rows),right_anchor_after_case_days=sum(r['right_anchor_after_case'] for r in rows),cross_month_bracket_days=sum(r['cross_month_bracket'] for r in rows),maximum_bracket_days=max(r['bracket_days'] for r in rows)))
    check('summary','seven unique season summaries',len(summaries)==7)
    expected_tables={'monthly':months,'daily':daily,'summary':summaries,'scenes':usage}
    for label,expected in expected_tables.items():
        filename,ids=TABLES[label];actual=pd.read_csv(production_run/filename)
        check('output_keys',filename+' unique keys',not actual.duplicated(ids).any())
        alookup={tuple(native(r[k]) for k in ids):r for r in actual.to_dict('records')}
        elookup={tuple(native(r[k]) for k in ids):r for r in expected}
        check('output_keys',filename+' exact complete keys',set(alookup)==set(elookup) and len(actual)==len(expected))
        for key,row in elookup.items():
            if key not in alookup:continue
            for field,value in row.items():
                if field in ids:continue
                value=native(value);av=alookup[key][field]
                if value is None or (isinstance(value,float) and math.isnan(value)):
                    check('production_missing',key+(field,),pd.isna(av))
                elif isinstance(value,bool):check('production_boolean',key+(field,),bool(av)==value)
                elif isinstance(value,(int,float)):close('production_numeric',key+(field,),av,value)
                else:check('production_labels',key+(field,),av==value)
    check('receipt','protocol/code pins',receipt['protocol_sha256']==PROTOCOL_SHA and receipt['code_sha256']==CODE_SHA)
    check('receipt','output counts',receipt['monthly_map_rows']==84 and receipt['retained_daily_rows']==1079 and receipt['site_season_summary_rows']==7 and receipt['scene_audit_rows']==len(usage))
    check('receipt','active maps',receipt['active_monthly_maps']==sum(r['used_by_retained_case_anchor'] for r in months))
    check('receipt','T1 training480',receipt['T1_training_reuse_days']==480)
    check('receipt','no target columns',receipt['target_columns_read']==[])
    for field in ['new_predictions','new_model_fits','new_FFP_runs','new_eligibility_approvals','Cloud_calls']:check('receipt',field,receipt[field]==0)
    check('projection','only identifier/predictor/optical columns decoded',all(set(cols).issubset(set(IDS+FEATURES+primary_columns+['B_uniform_pixels','C_core_weighted'])) for cols in columns_read.values()))
    for path,pinned in source_pins.items():check('unchanged_sources',path,sha(path)==pinned)
    check('unchanged_sources','protocol',sha(protocol)==PROTOCOL_SHA);check('unchanged_sources','production source',sha(code)==CODE_SHA)
    for filename,pinned in output_pins.items():check('unchanged_outputs',filename,sha(production_run/filename)==pinned)
    for label,rows in expected_tables.items():pd.DataFrame(rows).to_csv(out/('INDEPENDENT_'+label.upper()+'.csv'),index=False,float_format='%.17g')
    result={'status':'PASS' if not failed else 'FAIL','checks':sum(checks.values()),'failed_count':len(failed),'failures':failed,'checks_by_category':dict(checks),'maximum_difference_by_category':dict(maxdiff),'maximum_NDVI_replay_difference':max_replay,'protocol_sha256':PROTOCOL_SHA,'production_code_sha256':CODE_SHA,'verifier_code_sha256':sha(__file__),'production_run':str(production_run),'rows':{k:len(v) for k,v in expected_tables.items()},'source_pins':source_pins,'production_output_pins':output_pins,'projected_columns_read':columns_read,'target_columns_read':[],'model_files_read':[],'new_model_fits':0,'new_predictions':0,'new_FFP_runs':0,'seconds':time.monotonic()-started,'methods':'Own math.fsum weight arithmetic, bisect-left exact/nearest accepted-anchor bracketing, alternate linear-interpolation arithmetic; no production imports.','scope':'Target-free construction verification only; concentration does not establish independent sample size, crop purity, full daily footprint or model skill.'}
    (out/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (out/'README_EN.md').write_text(f"# Independent input-support verification\n\n**{result['status']}** — {result['checks']} checks; {len(failed)} failures.\n\nIndependently reconstructed all 84 monthly weight metrics, 1,079 retained case lineages, seven season summaries and {len(usage)} scene-usage records from pinned target-free source projections. Own `math.fsum` arithmetic and `bisect_left` anchor lookup use every original accepted site/year anchor. Exact reuse, cross-month brackets, source map months, <=24-day gaps, retained-case/map/anchor counters, all B/C and reference-X bindings, original optical curves and T1's exact 480 training-case reuse were checked. Source and output hashes remained unchanged.\n\nNo flux target column or serialized model was decoded, and no model was fitted, prediction generated or physical footprint computed. This verifies input construction; it supplies no crop-purity, independent-sample, field-authority or model-skill approval.\n",encoding='utf-8')
    print(json.dumps({'status':result['status'],'checks':result['checks'],'failed':len(failed),'rows':result['rows'],'max_replay':max_replay,'max_difference':dict(maxdiff)}))
    if failed:raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--protocol',required=True);p.add_argument('--run',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(a.protocol,a.run,a.out)
