"""Guarded portable adapter for the byte-identical original inference runner.

Only active path rebasing, worker I/O controls and reference comparisons are
added. The original inference, scoring, grouping and contrast code is retained.
"""
from pathlib import Path
import argparse,hashlib,json,multiprocessing,os,sys
sys.dont_write_bytecode=True
os.environ['PYTHONDONTWRITEBYTECODE']='1'
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'original_code'))
import frozen_crossprediction as original
import numpy as np
import pandas as pd
_ORIGINAL_WORKER=original.worker
_ORIGINAL_WORKER_INIT=original.worker_init
_GUARD_INSTALLED=False
_GUARD_PROBED=False
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def install_guard():
    global _GUARD_INSTALLED
    if _GUARD_INSTALLED:return
    from public_guard import install
    install([ROOT,os.environ['MAIZE_CROSSPREDICTION_OUTPUT_ROOT']])
    _GUARD_INSTALLED=True
def canary():
    source=load(ROOT/'PORTABLE_PATH_MAP.json')['original_source_probe']
    try:
        with open(source,'rb'):raise RuntimeError('Source-read canary unexpectedly succeeded')
    except PermissionError:return True
def guarded_worker_init():
    global _GUARD_PROBED
    os.chdir(ROOT);_ORIGINAL_WORKER_INIT();install_guard();_GUARD_PROBED=canary()
def guarded_worker(job):
    if not _GUARD_PROBED:raise RuntimeError('No verified per-worker I/O guard')
    values,record=_ORIGINAL_WORKER(job)
    record['portable_worker_guard_positive_control']=True;record['portable_worker_pid']=os.getpid()
    return values,record
def guard_report(_):
    if not _GUARD_PROBED:raise RuntimeError('Fresh worker guard was not initialized')
    return {'start_method':multiprocessing.get_start_method(),'old_source_read_blocked':canary(),'fits':0}
def validate():
    m=load(ROOT/'MANIFEST.json')
    for r in m['files']:
        p=ROOT/r['path']
        if not p.is_file() or p.stat().st_size!=r['bytes'] or sha(p)!=r['sha256']:raise ValueError('MANIFEST_MISMATCH: '+r['path'])
    binding=load(ROOT/'RELATIVE_INPUT_BINDINGS.json')
    for r in binding['sources']+binding['models']:
        p=Path(r['path'])
        if p.is_absolute() or '..' in p.parts:raise ValueError('Execution binding is not package-relative: '+str(p))
    return len(m['files'])
def compare_scientific(actual):
    expected=ROOT/'reference_results';checks=[]
    names=['cross_predictions.csv.gz','CELL_SCORES_NATIVE_UNITS.csv','CELL_SCORES_ESTIMANDS.csv','CONDITIONAL_CONTRASTS_NATIVE_UNITS.csv','CONDITIONAL_CONTRASTS_ESTIMANDS.csv','MAIN_CELLS_OWN_PATHS.csv','MAIN_CONDITIONALS_OWN_PATHS.csv']
    for name in names:
        a=pd.read_csv(actual/name,float_precision='round_trip');b=pd.read_csv(expected/name,float_precision='round_trip')
        pd.testing.assert_frame_equal(a,b,check_exact=True)
        checks.append({'file':name,'rows':len(a),'exact_frame_equality':True,'byte_identity':sha(actual/name)==sha(expected/name)})
    jobs=load(actual/'MODEL_REPLAY_LOG.json');assert len(jobs)==120 and all(j['portable_worker_guard_positive_control'] for j in jobs)
    return {'status':'PASS','scientific_outputs':checks,'all120_worker_task_guards_verified':True,'worker_pids':sorted(set(j['portable_worker_pid'] for j in jobs)),'maximum_numeric_difference':0.}
def main():
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['validate','replay','guardtest'],default='replay');p.add_argument('--out',type=Path);p.add_argument('--start-method',choices=['fork','spawn','forkserver'],default='spawn');a=p.parse_args()
    os.chdir(ROOT);count=validate()
    if a.mode=='validate':print(json.dumps({'status':'PASS','immutable_files_verified':count}));return
    if a.out is None:raise ValueError('Use a new explicit --out folder')
    out=a.out.resolve()
    if out==ROOT or ROOT in out.parents:raise ValueError('Output must be new and outside immutable public component')
    os.environ['MAIZE_CROSSPREDICTION_OUTPUT_ROOT']=str(out)
    install_guard();assert canary();multiprocessing.set_start_method(a.start_method,force=True)
    if a.mode=='guardtest':
        from concurrent.futures import ProcessPoolExecutor
        out.mkdir(parents=True,exist_ok=False)
        with ProcessPoolExecutor(max_workers=1,mp_context=multiprocessing.get_context('spawn'),initializer=guarded_worker_init) as pool:
            result=list(pool.map(guard_report,[0]))[0]
        save(out/'PORTABLE_GUARD_RECEIPT.json',{'status':'PASS','fresh_worker':result,'new_fits':0});print(json.dumps(result));return
    original.worker_init=guarded_worker_init;original.worker=guarded_worker
    original.run(argparse.Namespace(protocol=ROOT/'BC_CROSSPREDICTION_PROTOCOL_v0_1.json',expected_protocol_sha=sha(ROOT/'BC_CROSSPREDICTION_PROTOCOL_v0_1.json'),bindings=ROOT/'RELATIVE_INPUT_BINDINGS.json',expected_bindings_sha=sha(ROOT/'RELATIVE_INPUT_BINDINGS.json'),out=out))
    comparison=compare_scientific(out)
    receipt={'status':'PASS','immutable_files_verified_before':count,'package_manifest_sha256':sha(ROOT/'MANIFEST.json'),'adapter_sha256':sha(__file__),'original_runner_sha256':sha(ROOT/'original_code/frozen_crossprediction.py'),'protocol_sha256':sha(ROOT/'BC_CROSSPREDICTION_PROTOCOL_v0_1.json'),'relative_bindings_sha256':sha(ROOT/'RELATIVE_INPUT_BINDINGS.json'),'start_method':a.start_method,'main_source_guard_positive_control':True,'all_worker_controls':comparison,'source_access':'Generic offline data/config guard active; private original source locations are not distributed.','new_fits':0,'new_scientific_experiment':False}
    save(out/'PORTABLE_REPLAY_RECEIPT.json',receipt)
    assert validate()==count
    print(json.dumps({'status':'PASS','portable_scientific_comparison':comparison,'new_fits':0}),flush=True)
if __name__=='__main__':main()
