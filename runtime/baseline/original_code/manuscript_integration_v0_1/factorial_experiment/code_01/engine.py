"""Fixed CPU rice baselines on immutable v0.2 references; grouped OOF only.

No hyperparameter search, imputation, target clipping, prediction feedback or
new blind-test claim. The saved protocol/folds precede every model fit.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import math
import multiprocessing
import os
from pathlib import Path
import platform
import sys
import time
from zoneinfo import ZoneInfo

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.ensemble import HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

FEATURES = ["NDVI_modis", "EA_kPa_agera5", "WS_ms_agera5", "P_mm_agera5",
            "TA_max_C_agera5", "TA_min_C_agera5", "SW_IN_MJ_m2_d_agera5", "VPD_kPa_agera5"]
TARGETS = {"GPP":"gpp_gC_m2_d", "ET":"et_mm_d"}
KEYS = ["record_id", "observation_key", "date", "site_id", "system_id", "spatial_cluster",
        "harvest_year", "spatial_split_group", "season_split_group"]
MODELS = ["train_mean", "seasonal_harmonic", "ridge", "random_forest", "hist_gradient_boosting"]
AXES = ["cluster_loso", "year_loyo", "cluster_year_disjoint", "past_years"]
_FRAME = None


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda:stream.read(2**20), b""): h.update(chunk)
    return h.hexdigest()


def ids_hash(values):
    return hashlib.sha256("\n".join(map(str, values)).encode()).hexdigest()


def json_write(path, value):
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)


def metric_values(y_true, y_pred):
    y, p = np.asarray(y_true, dtype=float), np.asarray(y_pred, dtype=float)
    if y.ndim != 1 or y.shape != p.shape or not len(y) or not np.isfinite(y).all() or not np.isfinite(p).all():
        raise ValueError("Metrics require complete finite paired predictions")
    error = p-y
    sse = float(np.dot(error, error))
    sst = float(np.sum((y-y.mean())**2))
    skill = 1-sse/sst if len(y) >= 3 and np.ptp(y) > 0 and sst > 0 else None
    r = float(np.corrcoef(y, p)[0,1]) if len(y) >= 3 and np.ptp(y) > 0 and np.ptp(p) > 0 else None
    return {"n":len(y), "rmse":float(np.sqrt(np.mean(error**2))),
            "mae":float(np.mean(abs(error))), "bias_pred_minus_obs":float(np.mean(error)),
            "nse":skill, "r2_reference_variance":skill, "pearson_r":r,
            "squared_pearson_r":r*r if r is not None else None}


def calendar_matrix(dates):
    doy = pd.to_datetime(dates).dt.dayofyear.to_numpy(dtype=float)
    return np.column_stack([f(2*np.pi*k*doy/365.25) for k in (1,2) for f in (np.sin,np.cos)])


def load_cohort(reference, protocol):
    reference = Path(reference)
    for name, expected in protocol["expected_reference_files"].items():
        if sha(reference/name) != expected: raise ValueError("Frozen reference mismatch: "+name)
    contract = json.loads((reference/"feature_contract.json").read_text())
    if contract["primary_predictors_allowlist"] != FEATURES: raise ValueError("Predictor contract mismatch")
    panel = pd.read_parquet(reference/"nt_without_declared_holds.parquet").sort_values("record_id").reset_index(drop=True)
    predictors = pd.read_parquet(reference/"common_predictors.parquet")
    if panel.record_id.duplicated().any() or panel.observation_key.duplicated().any() or predictors.record_id.duplicated().any():
        raise ValueError("Duplicate primary identity")
    if len(panel) != protocol["expected_rows"]: raise ValueError("Unexpected cohort size")
    if panel[KEYS+["crop","gpp_reference_method"]].isna().any().any(): raise ValueError("Missing required reference identity")
    if not panel.crop.eq("rice").all() or not panel.gpp_reference_method.eq("NT").all(): raise ValueError("Cohort target method changed")
    for flag in ["qc_paired_dd_good_pass", "qc_calendar_pass", "common_inputs_complete", "panel_nt_without_declared_holds"]:
        if not pd.api.types.is_bool_dtype(panel[flag].dtype) or panel[flag].isna().any() or not panel[flag].eq(True).all():
            raise ValueError("Required source selection failed: "+flag)
    for flag in ["qc_provider_lineage_hold", "qc_naju_operator_hold", "qc_mixed_crop_footprint_hold"]:
        if not pd.api.types.is_bool_dtype(panel[flag].dtype) or panel[flag].isna().any() or not panel[flag].eq(False).all():
            raise ValueError("Scientific hold entered modeling: "+flag)
    if not panel.gpp_gC_m2_d.equals(panel.gpp_nt_gC_m2_d): raise ValueError("NT identity mismatch")
    if not panel.spatial_split_group.equals(panel.spatial_cluster): raise ValueError("Cluster-group mismatch")
    x = panel[KEYS+list(TARGETS.values())].merge(predictors[["record_id", *FEATURES]], on="record_id", validate="one_to_one")
    x = x.sort_values("record_id").reset_index(drop=True)
    pd.testing.assert_frame_equal(x[FEATURES], panel[FEATURES], check_exact=True, check_dtype=True)
    if not np.isfinite(x[FEATURES+list(TARGETS.values())].to_numpy(dtype=float)).all(): raise ValueError("Incomplete modeling cohort")
    date = pd.to_datetime(x.date)
    if not date.dt.year.eq(x.harvest_year).all() or not date.dt.month.between(5,10).all(): raise ValueError("Rice calendar mismatch")
    if x.spatial_cluster.nunique() != 4 or len(x[["spatial_cluster","harvest_year"]].drop_duplicates()) != 24:
        raise ValueError("Spatial/season support changed")
    return x


def design_splits(frame, min_prior_years=3, min_train_rows=300):
    """Return complete memberships; excluded rows cannot silently enter fit."""
    clusters, years = sorted(frame.spatial_cluster.unique()), sorted(frame.harvest_year.unique())
    specs = [("cluster_loso", c, None) for c in clusters]
    specs += [("year_loyo", None, int(y)) for y in years]
    specs += [("cluster_year_disjoint", c, int(y)) for c,y in
              sorted(set(zip(frame.spatial_cluster, frame.harvest_year)))]
    specs += [("past_years", None, int(y)) for i,y in enumerate(years) if i >= min_prior_years]
    folds, membership = [], []
    for axis, cluster, year in specs:
        if axis == "cluster_loso":
            test, train = frame.spatial_cluster.eq(cluster), ~frame.spatial_cluster.eq(cluster)
        elif axis == "year_loyo":
            test, train = frame.harvest_year.eq(year), ~frame.harvest_year.eq(year)
        elif axis == "cluster_year_disjoint":
            test = frame.spatial_cluster.eq(cluster) & frame.harvest_year.eq(year)
            train = ~frame.spatial_cluster.eq(cluster) & ~frame.harvest_year.eq(year)
        else:
            test, train = frame.harvest_year.eq(year), frame.harvest_year.lt(year)
        if not test.any() or int(train.sum()) < min_train_rows: raise ValueError("Insufficient prespecified fold support")
        if (train & test).any(): raise ValueError("Train/test row overlap")
        tr, te = frame.loc[train], frame.loc[test]
        if axis in {"cluster_loso", "cluster_year_disjoint"} and set(tr.spatial_cluster)&set(te.spatial_cluster):
            raise ValueError("Spatial overlap")
        if axis in {"year_loyo", "cluster_year_disjoint", "past_years"} and set(tr.harvest_year)&set(te.harvest_year):
            raise ValueError("Year overlap")
        if axis == "past_years" and not (tr.harvest_year.max() < te.harvest_year.min()): raise ValueError("Future training label")
        fold_id = f"f{sum(f['axis']==axis for f in folds):03d}"
        tri, tei = np.flatnonzero(train), np.flatnonzero(test)
        novel = sorted(set(te.site_id)-set(tr.site_id))
        fold = {"axis":axis, "fold_id":fold_id, "held_cluster":cluster, "held_year":year,
                "train_indices":tri.tolist(), "test_indices":tei.tolist(),
                "n_train":len(tr), "n_test":len(te), "n_excluded":int((~(train|test)).sum()),
                "train_clusters":sorted(tr.spatial_cluster.unique().tolist()), "test_clusters":sorted(te.spatial_cluster.unique().tolist()),
                "train_years":list(map(int,sorted(tr.harvest_year.unique()))), "test_years":list(map(int,sorted(te.harvest_year.unique()))),
                "train_date_max":str(tr.date.max()), "test_date_min":str(te.date.min()), "test_novel_sites":novel,
                "test_seen_sites":sorted(set(te.site_id)&set(tr.site_id)),
                "train_record_ids_sha256":ids_hash(tr.record_id), "test_record_ids_sha256":ids_hash(te.record_id),
                "combined_domain_shift": bool("JP-Mse" in set(te.site_id) and "JP-Mse" not in set(tr.site_id))}
        folds.append(fold)
        membership.append(pd.DataFrame({"axis":axis,"fold_id":fold_id,"record_id":frame.record_id,
                                        "role":np.where(test,"test",np.where(train,"train","excluded"))}))
    return folds, pd.concat(membership, ignore_index=True)


def score_tables(predictions):
    results = {}
    common = ["axis","model","target","seed"]
    for name, extra in {"pooled":[], "fold":["fold_id"], "cluster":["spatial_cluster"],
                        "site_season":["site_id","spatial_cluster","harvest_year","season_split_group"]}.items():
        groups = common+extra
        results[name] = pd.DataFrame([{**dict(zip(groups,key)), **metric_values(part.y_true,part.y_pred)}
                                      for key,part in predictions.groupby(groups, sort=True)])
    macro = []
    for key, part in results["site_season"].groupby(common, sort=True):
        row = {**dict(zip(common,key)), "n_site_seasons":len(part), "n_clusters":int(part.spatial_cluster.nunique()),
               "macro_site_season_rmse":float(part.rmse.mean()), "median_site_season_rmse":float(part.rmse.median()),
               "macro_site_season_mae":float(part.mae.mean()), "macro_site_season_bias":float(part.bias_pred_minus_obs.mean()),
               "macro_site_season_nse":float(part.nse.mean()) if part.nse.notna().any() else None,
               "n_finite_season_nse":int(part.nse.notna().sum())}
        c = results["cluster"]
        mask = np.ones(len(c), dtype=bool)
        for name,value in zip(common,key): mask &= c[name].eq(value)
        row["macro_cluster_rmse"] = float(c.loc[mask,"rmse"].mean())
        macro.append(row)
    results["macro_summary"] = pd.DataFrame(macro)
    return results


def worker_init(snapshot):
    global _FRAME
    _FRAME = pd.read_parquet(snapshot)


def fit_task(task):
    fold, model_name, target, protocol, root = task
    tr, te = _FRAME.iloc[fold["train_indices"]], _FRAME.iloc[fold["test_indices"]]
    y = tr[TARGETS[target]].to_numpy(dtype=float)
    seed = protocol["seed"]
    estimator, train_mean = None, float(y.mean())
    features = []
    start = time.perf_counter()
    with threadpool_limits(limits=1):
        if model_name == "train_mean":
            prediction = np.full(len(te), train_mean)
        elif model_name == "seasonal_harmonic":
            features = ["sin_doy_1","cos_doy_1","sin_doy_2","cos_doy_2"]
            estimator = Ridge(**protocol["models"][model_name])
            estimator.fit(calendar_matrix(tr.date),y)
            prediction = estimator.predict(calendar_matrix(te.date))
        else:
            features = FEATURES.copy()
            hp = dict(protocol["models"][model_name])
            if model_name == "ridge": estimator = make_pipeline(StandardScaler(),Ridge(**hp))
            elif model_name == "random_forest": estimator = RandomForestRegressor(**hp, random_state=seed, n_jobs=1)
            elif model_name == "hist_gradient_boosting": estimator = HistGradientBoostingRegressor(**hp,random_state=seed)
            else: raise ValueError("Unapproved model")
            estimator.fit(tr[FEATURES],y)
            prediction = estimator.predict(te[FEATURES])
    if not np.isfinite(prediction).all(): raise ValueError("Nonfinite model predictions")
    path = Path(root)/"models"/fold["axis"]/fold["fold_id"]/target/(model_name+".joblib")
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists(): raise FileExistsError(path)
    model_record = {k:fold[k] for k in ["axis","fold_id","train_record_ids_sha256","test_record_ids_sha256"]}
    model_record.update(estimator=estimator, training_mean=train_mean, feature_columns=features, model=model_name,
                        target=target, seed=seed, protocol_sha256=protocol["_protocol_sha256"],
                        n_train=len(tr), n_test=len(te), hyperparameters=protocol["models"][model_name])
    joblib.dump(model_record,path,compress=3)
    output = te[KEYS].copy()
    output["axis"],output["fold_id"],output["model"],output["target"],output["seed"] = fold["axis"],fold["fold_id"],model_name,target,seed
    output["y_true"],output["y_pred"] = te[TARGETS[target]].to_numpy(dtype=float),prediction
    fit = {k:fold[k] for k in ["axis","fold_id","n_train","n_test"]}
    fit.update(model=model_name,target=target,seed=seed,seconds=time.perf_counter()-start,
               model_path=str(path),model_sha256=sha(path),features=features)
    return output,fit


def run(args):
    if args.output.exists(): raise FileExistsError("New baseline run directory required")
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    if protocol["features"] != FEATURES or list(protocol["models"]) != MODELS or protocol["axes"] != AXES:
        raise ValueError("Unapproved protocol feature/model/split roster")
    if protocol["models"]["hist_gradient_boosting"].get("early_stopping") is not False: raise ValueError("Random-row early stopping forbidden")
    args.output.mkdir(parents=True)
    protocol_hash = sha(args.protocol)
    frame = load_cohort(args.reference,protocol)
    frame.to_parquet(args.output/"data_snapshot.parquet",index=False)
    folds,membership = design_splits(frame,protocol["min_prior_years"],protocol["min_train_rows"])
    membership.to_parquet(args.output/"split_membership.parquet",index=False)
    json_write(args.output/"folds.json",folds)
    counts = {axis:sum(f["axis"]==axis for f in folds) for axis in AXES}
    if counts != protocol["expected_fold_counts"]: raise ValueError("Unexpected fold structure")
    for axis in AXES:
        test = membership.loc[(membership.axis==axis)&membership.role.eq("test")]
        if test.record_id.duplicated().any(): raise ValueError("Repeated OOF case")
        if axis != "past_years" and set(test.record_id) != set(frame.record_id): raise ValueError("Incomplete OOF coverage")
        if len(test) != protocol["expected_test_rows"][axis]: raise ValueError("Unexpected OOF sample support")
    code_hashes = {str(p.resolve()):sha(p) for p in sorted(Path(__file__).parent.glob("*.py"))}
    input_hashes = {str((args.reference/name).resolve()):sha(args.reference/name) for name in protocol["expected_reference_files"]}
    freeze = {"frozen_utc":datetime.now(timezone.utc).isoformat(),"frozen_kst":datetime.now(ZoneInfo("Asia/Seoul")).isoformat(),
              "protocol_sha256":protocol_hash,"protocol_source":str(args.protocol.resolve()),"inputs":input_hashes,"code":code_hashes,
              "data_snapshot_sha256":sha(args.output/"data_snapshot.parquet"),"split_membership_sha256":sha(args.output/"split_membership.parquet"),
              "folds_sha256":sha(args.output/"folds.json"),"fold_counts":counts,"cohort_rows":len(frame),
              "status":"frozen_before_any_fit","historical_test_exposure":True,"new_blind_test_claim":False}
    json_write(args.output/"protocol.json",protocol)
    json_write(args.output/"protocol_freeze.json",freeze)
    summary = frame.groupby(["site_id","spatial_cluster","harvest_year"],sort=True).agg(n=("record_id","size"),date_min=("date","min"),date_max=("date","max")).reset_index()
    summary.to_csv(args.output/"cohort_summary.csv",index=False)
    protocol["_protocol_sha256"] = protocol_hash
    tasks = [(fold,model,target,protocol,str(args.output)) for fold in folds for model in MODELS for target in TARGETS]
    print(json.dumps({"status":"frozen","cohort_rows":len(frame),"fold_counts":counts,"fits":len(tasks),"protocol_sha256":protocol_hash}),flush=True)
    started = time.perf_counter()
    predictions,fit_log = [],[]
    with ProcessPoolExecutor(max_workers=protocol["workers"],mp_context=multiprocessing.get_context("spawn"),
                             initializer=worker_init,initargs=(str(args.output/"data_snapshot.parquet"),)) as pool:
        futures = [pool.submit(fit_task,task) for task in tasks]
        for i,future in enumerate(as_completed(futures),1):
            pred,log = future.result()
            predictions.append(pred);fit_log.append(log)
            if i % 50 == 0 or i == len(tasks): print(json.dumps({"fits_completed":i,"total":len(tasks),"wall_seconds":round(time.perf_counter()-started,2)}),flush=True)
    pred = pd.concat(predictions,ignore_index=True).sort_values(["axis","fold_id","model","target","record_id"]).reset_index(drop=True)
    if pred.duplicated(["axis","model","target","record_id"]).any(): raise ValueError("Duplicate OOF predictions")
    pred.to_parquet(args.output/"predictions.parquet",index=False)
    pred.to_csv(args.output/"predictions.csv.gz",index=False,float_format="%.17g")
    pd.testing.assert_frame_equal(pred,pd.read_parquet(args.output/"predictions.parquet"),check_exact=True)
    for name,table in score_tables(pred).items():
        table.to_csv(args.output/("macro_summary.csv" if name=="macro_summary" else f"scores_{name}.csv"),index=False,float_format="%.17g")
    json_write(args.output/"fit_log.json",sorted(fit_log,key=lambda r:(r["axis"],r["fold_id"],r["target"],r["model"])))
    for path,expected in input_hashes.items():
        if sha(path) != expected: raise ValueError("Original input changed during training")
    if sha(args.protocol) != protocol_hash: raise ValueError("Protocol changed during training")
    if any(sha(path)!=expected for path,expected in code_hashes.items()): raise ValueError("Code changed during training")
    verification = {"status":"passed","targets_and_features_unchanged":True,"source_inputs_unchanged":True,
                    "split_checks_passed":True,"expected_fit_count":len(tasks),"actual_fit_count":len(fit_log),
                    "prediction_rows":len(pred),"rows":len(frame),"fold_counts":counts,"retrospective":True,
                    "hyperparameter_search":False,"new_blind_test_claim":False,"new_agent_value_claim":False,
                    "model_training_wall_seconds":time.perf_counter()-started}
    json_write(args.output/"verification.json",verification)
    manifest = {"created_utc":datetime.now(timezone.utc).isoformat(),"created_kst":datetime.now(ZoneInfo("Asia/Seoul")).isoformat(),
                "protocol_sha256":protocol_hash,"inputs":input_hashes,"code":code_hashes,
                "environment":{"python":sys.version,"pandas":pd.__version__,"numpy":np.__version__,"sklearn":sklearn.__version__,"host":platform.node()},
                "outputs":{str(p.relative_to(args.output).as_posix()):{"bytes":p.stat().st_size,"sha256":sha(p)}
                           for p in sorted(args.output.rglob("*")) if p.is_file()}}
    json_write(args.output/"run_manifest.json",manifest)
    print(json.dumps(verification),flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference",type=Path,required=True)
    p.add_argument("--protocol",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())
