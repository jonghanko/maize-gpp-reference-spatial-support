"""Independent T1 post-fit verification: saved-model replay only, never fitting.

Prepared separately from the production adapter. No engine metric/calendar/
prepare/fit function is imported or called. Writes are confined to this folder.
Run only after parent confirms fitrun01 completion.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time

sys.dont_write_bytecode = True
import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parent
CPU = Path('source://analysis/early_mead_v0_1')
PROTOCOL_SHA = '233f0bc5ad008bbe87a588180e796f9dedf9d2e934e9d077b710c681b529248d'
ENGINE_SHA = 'c6dab8c357fd0432dbf8b06fa263d04e19b4527ad41984982862ad69409b2e9d'
ARMS = ['B_uniform_pixels', 'C_core_weighted']
REFS = ['VUT_REF', 'CUT_REF']
MODELS = ['train_mean', 'seasonal_harmonic', 'ridge', 'random_forest', 'hist_gradient_boosting']
FEATURES = ['NDVI_modis', 'EA_kPa_agera5', 'WS_ms_agera5', 'P_mm_agera5', 'TA_max_C_agera5',
            'TA_min_C_agera5', 'SW_IN_MJ_m2_d_agera5', 'VPD_kPa_agera5']
KEYS = ['record_id', 'observation_key', 'date', 'site_id', 'system_id', 'spatial_cluster',
        'harvest_year', 'spatial_split_group', 'season_split_group']
BOUNDS = {('US-Ne1', 2010): ('2010-04-19', '2010-09-21'),
          ('US-Ne1', 2011): ('2011-05-18', '2011-10-26'),
          ('US-Ne3', 2011): ('2011-05-03', '2011-10-19')}
SEASONS = {('US-Ne1', 2010): 128, ('US-Ne1', 2011): 148, ('US-Ne3', 2011): 161}
CHECKS = []
FIT_CALL_ATTEMPTS = 0

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for chunk in iter(lambda: f.read(2**20), b''):
            h.update(chunk)
    return h.hexdigest()

def read(p):
    return json.loads(Path(p).read_text(encoding='utf-8'))

def check(name, ok, details=None):
    CHECKS.append({'name': name, 'passed': bool(ok), 'details': details})
    if not ok:
        raise AssertionError(name)

def close(a, b, atol=5e-12, rtol=2e-12):
    return np.allclose(np.asarray(a, dtype=float), np.asarray(b, dtype=float),
                       atol=atol, rtol=rtol, equal_nan=True)

def normalized(d):
    d = d.copy()
    d['date'] = pd.to_datetime(d.date).dt.strftime('%Y-%m-%d')
    assert not d.duplicated(['site_id', 'date']).any()
    return d.sort_values(['site_id', 'date']).reset_index(drop=True)

def idhash(series):
    return hashlib.sha256('\n'.join(str(v) for v in series).encode()).hexdigest()

def physical(d):
    return set(zip(d.site_id, d.date))

def forbid_fit(*args, **kwargs):
    global FIT_CALL_ATTEMPTS
    FIT_CALL_ATTEMPTS += 1
    raise RuntimeError('Post-fit verifier forbids estimator.fit')

def install_fit_guards():
    for cls in [Ridge, RandomForestRegressor, HistGradientBoostingRegressor, StandardScaler, Pipeline]:
        cls.fit = forbid_fit

def calendar(dates):
    # Independent implementation; the production calendar_matrix is never called.
    return np.asarray([[fun(2 * math.pi * k * datetime.fromisoformat(str(d)[:10]).timetuple().tm_yday / 365.25)
        for k in (1, 2) for fun in (math.sin, math.cos)] for d in dates], dtype=float)

def metrics(truth, predicted):
    # Python fsum/explicit centered products, independent of engine.metric_values.
    y, p = [float(x) for x in truth], [float(x) for x in predicted]
    assert len(y) == len(p) and y and all(math.isfinite(x) for x in y+p)
    n = len(y)
    error = [b-a for a, b in zip(y, p)]
    ybar, pbar = math.fsum(y)/n, math.fsum(p)/n
    sse = math.fsum(e*e for e in error)
    sy = math.fsum((a-ybar)**2 for a in y)
    sp = math.fsum((a-pbar)**2 for a in p)
    skill = 1-sse/sy if n >= 3 and max(y) > min(y) and sy > 0 else math.nan
    corr = math.fsum((a-ybar)*(b-pbar) for a, b in zip(y, p))/math.sqrt(sy*sp) \
        if n >= 3 and max(y) > min(y) and max(p) > min(p) and sy > 0 and sp > 0 else math.nan
    return {'n': n, 'rmse': math.sqrt(sse/n), 'mae': math.fsum(abs(e) for e in error)/n,
        'bias_pred_minus_obs': math.fsum(error)/n, 'nse': skill, 'r2_reference_variance': skill,
        'pearson_r': corr, 'squared_pearson_r': corr*corr if math.isfinite(corr) else math.nan}

def table_equal(name, actual, expected, key_columns):
    check(name+'_columns', set(actual.columns) == set(expected.columns))
    check(name+'_unique_groups', not actual.duplicated(key_columns).any())
    a = actual.sort_values(key_columns).reset_index(drop=True)
    e = expected.sort_values(key_columns).reset_index(drop=True)
    check(name+'_group_keys', len(a) == len(e) and a[key_columns].equals(e[key_columns]))
    for col in expected.columns:
        if col in key_columns:
            continue
        if col in ['contrast', 'primary']:
            check(name+'_'+col, a[col].equals(e[col]))
        else:
            check(name+'_'+col, close(a[col], e[col]))

def verify(args):
    started = time.perf_counter()
    output = args.output.resolve()
    check('output_confined_to_independent_postfit', output == ROOT or ROOT in output.parents)
    check('verification_outputs_not_replaced', not any((output/name).exists() for name in
        ['SAVED_MODEL_REPLAY.csv', 'INDEPENDENT_INPUT_CONTRASTS.csv', 'INDEPENDENT_POSTFIT_VERIFICATION.json']))
    plan, run = args.plan.resolve(), args.run.resolve()
    p = read(plan/'protocol.json')
    f = read(plan/'freeze.json')
    check('exact_parent_protocol_SHA', sha(plan/'protocol.json') == f['protocol_sha256'] == PROTOCOL_SHA)
    check('frozen_task_budget_20_with_437_tests', p['n_train'] == 480 and p['n_test'] == 437
          and p['expected_fits'] == 20 and p['expected_own_prediction_records'] == 8740
          and p['expected_cross_prediction_records'] == 17480)
    check('feature_model_arm_reference_rosters', p['features'] == FEATURES and list(p['models']) == MODELS
          and p['arms'] == ARMS and p['references'] == REFS and p['seed'] == 42)
    check('primary_estimand_is_frozen_RF_own_macro', p['primary_model'] == 'random_forest'
          and p['primary_metric'] == 'equal_test_site_season_mean_RMSE'
          and p['primary_contrast'] == 'C_core_weighted_minus_B_uniform_pixels')
    freeze_paths = {str(plan/'freeze.json'): sha(plan/'freeze.json')}
    for group in ['inputs', 'code', 'plan_files']:
        for path, digest in f[group].items():
            check('frozen_'+group+'_'+str(len(freeze_paths)), sha(path) == digest)
            freeze_paths[path] = digest
    config = read(args.config)
    check('source_config_matches_protocol', sha(args.config) == p['source_config_sha256'])
    check('frozen_engine_SHA', sha(config['engine']) == ENGINE_SHA)
    run_manifest_path = run/'run_manifest.json'
    manifest = read(run_manifest_path)['files']
    run_paths = {str(run_manifest_path): sha(run_manifest_path)}
    for rel, r in manifest.items():
        path = run/rel
        check('run_manifest_'+rel, path.is_file() and sha(path) == r['sha256'] and path.stat().st_size == r['bytes'])
        run_paths[str(path)] = r['sha256']
    execution = read(run/'execution.json')
    check('completed_execution_receipt', execution['status'] == 'complete' and execution['protocol_sha256'] == PROTOCOL_SHA
          and execution['fits'] == 20 and execution['n_train'] == 480 and execution['n_test'] == 437
          and execution['own_prediction_records'] == 8740 and execution['cross_prediction_records'] == 17480)
    check('CPU_only_no_Cloud_no_forecast_claim', execution['Cloud_calls'] == 0 and execution['GPU'] is False
          and execution['new_region_or_forecast_claim'] is False and p['parent_S1_regional_execution_released'] is False)
    authorization = read(args.authorization)
    check('authorization_hash_and_exact_scope', sha(args.authorization) == execution['authorization_sha256']
          and authorization['status'] == 'authorized_conditional_retrospective_T1'
          and authorization['protocol_sha256'] == PROTOCOL_SHA)
    prefit = read(authorization['independent_prefit_receipt'])
    check('independent_prefit_receipt_hash_and_status', sha(authorization['independent_prefit_receipt']) ==
          authorization['independent_prefit_receipt_sha256'] and prefit['status'] == 'pass'
          and prefit['protocol_sha256'] == PROTOCOL_SHA)
    freeze_paths[str(args.authorization.resolve())] = sha(args.authorization)
    freeze_paths[authorization['independent_prefit_receipt']] = authorization['independent_prefit_receipt_sha256']
    check('runtime_consistent', p['runtime'] == {'python': platform.python_version(), 'numpy': np.__version__,
          'pandas': pd.__version__, 'sklearn': sklearn.__version__})
    fold = read(plan/'fold.json')
    check('exact_fold_indices', fold['train_indices'] == list(range(480)) and fold['test_indices'] == list(range(480, 917))
          and fold['axis'] == 'temporal_hindcast' and fold['fold_id'] == 'f000')

    # Reconstruct candidate keys and values directly from the three frozen sources.
    weather = normalized(pd.read_csv(config['weather7'], float_precision='round_trip'))
    optical = normalized(pd.read_csv(config['optical_daily'], float_precision='round_trip'))
    target = normalized(pd.read_parquet(config['paired_target_payload']))
    ycols = config['target_columns']
    check('source_cardinalities_476_552_450', len(weather) == 476 and len(optical) == 552 and len(target) == 450)
    expected = weather[['site_id', 'date', *FEATURES[1:]]].merge(
        optical[['site_id', 'date', *ARMS]], on=['site_id', 'date'], validate='one_to_one')
    expected = expected.merge(target[['site_id', 'date', *ycols.values()]], on=['site_id', 'date'], validate='one_to_one')
    mask = np.isfinite(expected[[*FEATURES[1:], *ARMS]].to_numpy(float)).all(axis=1)
    expected = normalized(expected.loc[mask])
    check('independently_reconstructed_437_physical_keys', len(expected) == 437 and len(physical(expected)) == 437)
    check('source_season_counts_128_148_161', expected.assign(year=pd.to_datetime(expected.date).dt.year)
          .groupby(['site_id', 'year']).size().to_dict() == SEASONS)
    for (site, year), (lo, hi) in BOUNDS.items():
        d = expected.loc[expected.site_id.eq(site) & pd.to_datetime(expected.date).dt.year.eq(year)]
        check(f'published_field_bounds_{site}_{year}', d.date.between(lo, hi).all()
              and pd.to_datetime(d.date).dt.month.between(5, 10).all())
    expected['record_id'] = ['T1:'+s+':'+d for s, d in zip(expected.site_id, expected.date)]
    expected['observation_key'] = expected.record_id
    expected['system_id'] = expected.site_id
    expected['spatial_cluster'] = 'Mead_NE'
    expected['harvest_year'] = pd.to_datetime(expected.date).dt.year.astype(int)
    expected['spatial_split_group'] = 'Mead_NE'
    expected['season_split_group'] = expected.site_id+'__'+expected.harvest_year.astype(str)
    saved_selected = normalized(pd.read_parquet(plan/'paired_test_source_snapshot.parquet'))
    check('stored_paired_snapshot_exact_keys', saved_selected[KEYS].equals(expected[KEYS]))
    check('stored_paired_snapshot_source_values', close(saved_selected[[*FEATURES[1:], *ARMS, *ycols.values()]],
          expected[[*FEATURES[1:], *ARMS, *ycols.values()]], atol=1e-14, rtol=1e-14))
    frames = {}
    for arm in ARMS:
        for ref in REFS:
            frame = pd.read_parquet(plan/(arm+'__'+ref+'_snapshot.parquet'))
            frame['date'] = pd.to_datetime(frame.date).dt.strftime('%Y-%m-%d')
            tr, te = frame.iloc[:480].reset_index(drop=True), frame.iloc[480:].reset_index(drop=True)
            source_train = normalized(pd.read_parquet(config['training_snapshots'][arm][ref]))
            check(f'snapshot_{arm}_{ref}_917', len(frame) == 917 and len(tr) == 480 and len(te) == 437)
            check(f'train_source_keys_{arm}_{ref}', tr[KEYS].equals(source_train[KEYS]))
            check(f'train_source_values_{arm}_{ref}', close(tr[[*FEATURES, 'gpp_gC_m2_d']], source_train[[*FEATURES, 'gpp_gC_m2_d']], atol=0, rtol=0))
            check(f'test_keys_source_{arm}_{ref}', te[KEYS].equals(expected[KEYS]))
            check(f'test_truth_source_{arm}_{ref}', close(te.gpp_gC_m2_d, expected[ycols[ref]], atol=0, rtol=0))
            check(f'test_optical_source_{arm}_{ref}', close(te.NDVI_modis, expected[arm], atol=1e-14, rtol=1e-14))
            check(f'test_weather_source_{arm}_{ref}', close(te[FEATURES[1:]], expected[FEATURES[1:]], atol=1e-14, rtol=1e-14))
            check(f'no_train_test_physical_overlap_{arm}_{ref}', not physical(tr)&physical(te)
                  and not set(tr.record_id)&set(te.record_id) and max(te.date) < min(tr.date))
            check(f'saved_membership_hashes_{arm}_{ref}', idhash(tr.record_id) == fold['train_record_ids_sha256']
                  and idhash(te.record_id) == fold['test_record_ids_sha256'])
            check(f'train_roster_Mead_only_{arm}_{ref}', set(zip(tr.site_id, tr.harvest_year)) ==
                  {('US-Ne1', 2012), ('US-Ne1', 2013), ('US-Ne3', 2013)} and frame.spatial_cluster.eq('Mead_NE').all())
            frames[arm, ref] = (tr, te)

    own = pd.read_parquet(run/'own_predictions.parquet')
    cross = pd.read_parquet(run/'cross_predictions.parquet')
    fit_log = read(run/'fit_log.json')
    task_keys = ['arm', 'fit_reference', 'model']
    check('20_unique_saved_model_tasks', len(fit_log) == 20 and len({tuple(r[k] for k in task_keys) for r in fit_log}) == 20)
    check('8740_unique_own_predictions', len(own) == 8740 and not own.duplicated([*task_keys, 'record_id']).any())
    check('17480_unique_both_score_reference_records', len(cross) == 17480 and
          not cross.duplicated([*task_keys, 'score_reference', 'record_id']).any())
    check('prediction_protocol_fields', own.axis.eq('temporal_hindcast').all() and own.fold_id.eq('f000').all()
          and own.target.eq('GPP').all() and own.seed.eq(42).all() and np.isfinite(own[['y_true', 'y_pred']]).all().all())
    install_fit_guards()
    replay_rows = []
    for log in fit_log:
        arm, ref, model = (log[k] for k in task_keys)
        tr, te = frames[arm, ref]
        key = arm+'__'+ref+'__'+model
        path = run/arm/ref/'models/temporal_hindcast/f000/GPP'/(model+'.joblib')
        check(key+'_model_path_hash', str(path) == log['model_path'] and sha(path) == log['model_sha256'])
        record = joblib.load(path)
        check(key+'_saved_record_identity', record['model'] == model and record['target'] == 'GPP'
              and record['seed'] == 42 and record['axis'] == 'temporal_hindcast' and record['fold_id'] == 'f000'
              and record['protocol_sha256'] == PROTOCOL_SHA and record['n_train'] == 480 and record['n_test'] == 437)
        check(key+'_saved_membership_hashes', record['train_record_ids_sha256'] == fold['train_record_ids_sha256']
              and record['test_record_ids_sha256'] == fold['test_record_ids_sha256'])
        check(key+'_hyperparameters', record['hyperparameters'] == p['models'][model])
        check(key+'_training_mean', close(record['training_mean'], math.fsum(tr.gpp_gC_m2_d)/480, atol=1e-13, rtol=1e-14))
        estimator = record['estimator']
        if model == 'train_mean':
            check(key+'_mean_control_has_no_estimator', estimator is None and record['feature_columns'] == [])
            prediction = np.full(437, record['training_mean'])
        else:
            if model == 'seasonal_harmonic':
                check(key+'_calendar_four_feature_contract', isinstance(estimator, Ridge) and estimator.n_features_in_ == 4
                      and record['feature_columns'] == ['sin_doy_1', 'cos_doy_1', 'sin_doy_2', 'cos_doy_2'])
                check(key+'_calendar_params', all(estimator.get_params()[k] == v for k, v in p['models'][model].items()))
                design = calendar(te.date)
            else:
                check(key+'_eight_feature_order', record['feature_columns'] == FEATURES
                      and estimator.n_features_in_ == 8 and list(estimator.feature_names_in_) == FEATURES)
                design = te[FEATURES]
                if model == 'ridge':
                    check(key+'_ridge_pipeline', isinstance(estimator, Pipeline) and list(estimator.named_steps) == ['standardscaler', 'ridge'])
                    scaler, ridge = estimator.named_steps['standardscaler'], estimator.named_steps['ridge']
                    train_x = tr[FEATURES].to_numpy(float)
                    check(key+'_train_only_scaler_mean', close(scaler.mean_, train_x.mean(axis=0)))
                    check(key+'_train_only_scaler_variance', close(scaler.var_, train_x.var(axis=0)))
                    check(key+'_train_only_scaler_scale', close(scaler.scale_, np.where(train_x.var(axis=0) == 0, 1, np.sqrt(train_x.var(axis=0)))))
                    check(key+'_train_only_scaler_sample_count', np.asarray(scaler.n_samples_seen_ == 480).all()
                          and scaler.n_features_in_ == 8 and list(scaler.feature_names_in_) == FEATURES)
                    check(key+'_ridge_params', all(ridge.get_params()[k] == v for k, v in p['models'][model].items()))
                else:
                    cls = RandomForestRegressor if model == 'random_forest' else HistGradientBoostingRegressor
                    check(key+'_estimator_class_params_seed', isinstance(estimator, cls) and estimator.random_state == 42
                          and all(estimator.get_params()[k] == v for k, v in p['models'][model].items()))
                    if model == 'random_forest':
                        tree_seeds = np.random.RandomState(42).randint(np.iinfo(np.int32).max, size=200)
                        check(key+'_RF200_singlethread_seed_sequence', len(estimator.estimators_) == 200 and estimator.n_jobs == 1
                              and np.array_equal([t.random_state for t in estimator.estimators_], tree_seeds))
                        check(key+'_RF_training_sample_totals_and_leaf_guard', all(t.n_features_in_ == 8
                              and t.tree_.weighted_n_node_samples[0] == 480
                              and np.min(t.tree_.n_node_samples[t.tree_.children_left == -1]) >= 5 for t in estimator.estimators_))
                    else:
                        check(key+'_HGB_fixed200_no_earlystop', estimator.n_iter_ == 200 and estimator.do_early_stopping_ is False)
            with threadpool_limits(limits=1):
                prediction = estimator.predict(design)
        pred = own.loc[own.arm.eq(arm) & own.fit_reference.eq(ref) & own.model.eq(model)].sort_values('record_id').reset_index(drop=True)
        test_sorted = te.sort_values('record_id').reset_index(drop=True)
        replay_sorted = pd.DataFrame({'record_id': te.record_id, 'prediction': prediction}).sort_values('record_id').prediction.to_numpy()
        check(key+'_exact_437_prediction_keys_and_identity', len(pred) == 437 and pred[KEYS].equals(test_sorted[KEYS]))
        check(key+'_own_truth_exact', close(pred.y_true, test_sorted.gpp_gC_m2_d, atol=0, rtol=0))
        check(key+'_saved_model_prediction_replay', close(pred.y_pred, replay_sorted, atol=2e-12, rtol=2e-13))
        replay_rows.append({'arm': arm, 'fit_reference': ref, 'model': model, 'records': 437,
                            'maximum_absolute_replay_difference': float(np.max(abs(pred.y_pred.to_numpy()-replay_sorted)))})
    for arm in ARMS:
        for ref in REFS:
            saved = pd.read_parquet(run/arm/ref/'own_predictions.parquet').sort_values(['model', 'record_id']).reset_index(drop=True)
            overall = own.loc[own.arm.eq(arm) & own.fit_reference.eq(ref), saved.columns].sort_values(['model', 'record_id']).reset_index(drop=True)
            check(f'per_taskfile_equals_global_{arm}_{ref}', saved.equals(overall))

    # Cross-reference pairing independently keys each own prediction to source truth.
    source_truth = {r: expected.set_index('record_id')[ycols[r]].to_dict() for r in REFS}
    own_lookup = own.set_index([*task_keys, 'record_id']).sort_index()
    check('all_cross_reference_rosters', set(cross.score_reference) == set(REFS))
    cross_ok, max_cross_error = True, 0.0
    for key, part in cross.groupby([*task_keys, 'score_reference'], sort=True):
        arm, fit_ref, model, score_ref = key
        base = own_lookup.loc[(arm, fit_ref, model)]
        cross_ok &= len(part) == 437 and set(part.record_id) == set(source_truth[score_ref])
        base_meta = base.reset_index().sort_values('record_id').reset_index(drop=True)
        cross_meta = part.sort_values('record_id').reset_index(drop=True)
        cross_ok &= cross_meta[[*KEYS, 'axis', 'fold_id', 'target', 'seed']].equals(
            base_meta[[*KEYS, 'axis', 'fold_id', 'target', 'seed']])
        desired = part.record_id.map(source_truth[score_ref]).to_numpy(float)
        predictions = part.record_id.map(base.y_pred.to_dict()).to_numpy(float)
        cross_ok &= close(part.y_true, desired, atol=0, rtol=0) and close(part.y_pred, predictions, atol=0, rtol=0)
        max_cross_error = max(max_cross_error, float(np.max(abs(part.y_true.to_numpy()-desired))))
    check('all_40n_manual_cross_reference_truth_and_prediction_pairs', cross_ok)
    groups = ['arm', 'model', 'fit_reference', 'score_reference']
    pooled, season = [], []
    for key, d in cross.groupby(groups, sort=True):
        pooled.append({**dict(zip(groups, key)), **metrics(d.y_true, d.y_pred)})
        for (site, year), part in d.groupby(['site_id', 'harvest_year'], sort=True):
            season.append({**dict(zip(groups, key)), 'site_id': site, 'harvest_year': int(year), **metrics(part.y_true, part.y_pred)})
    pooled, season = pd.DataFrame(pooled), pd.DataFrame(season)
    macro = season.groupby(groups, as_index=False).agg(macro_site_season_rmse=('rmse', 'mean'),
        macro_site_season_mae=('mae', 'mean'), n_site_seasons=('site_id', 'size'))
    check('all_40_metric_paths_have_three_test_seasons', len(pooled) == 40 and len(season) == 120
          and len(macro) == 40 and macro.n_site_seasons.eq(3).all())
    contrasts = []
    for (model, fit_ref, score_ref), group in macro.groupby(['model', 'fit_reference', 'score_reference']):
        b = float(group.loc[group.arm.eq(ARMS[0]), 'macro_site_season_rmse'].iloc[0])
        c = float(group.loc[group.arm.eq(ARMS[1]), 'macro_site_season_rmse'].iloc[0])
        contrasts.append({'model': model, 'fit_reference': fit_ref, 'score_reference': score_ref,
            'contrast': 'C_minus_B', 'B_RMSE': b, 'C_RMSE': c, 'absolute_difference': c-b,
            'relative_percent': 100*(c-b)/b if b else np.nan,
            'primary': model == 'random_forest' and fit_ref == score_ref})
    contrasts = pd.DataFrame(contrasts)
    for name, expected_table, table_keys in [('scores_pooled', pooled, groups),
        ('scores_site_season', season, [*groups, 'site_id', 'harvest_year']),
        ('scores_macro', macro, groups), ('input_contrasts', contrasts, ['model', 'fit_reference', 'score_reference'])]:
        table_equal(name, pd.read_csv(run/(name+'.csv'), float_precision='round_trip'), expected_table, table_keys)
    check('exact_two_primary_RF_own_contrasts', int(contrasts.primary.sum()) == 2
          and set(contrasts.loc[contrasts.primary, 'fit_reference']) == set(REFS))
    controls = own.loc[own.model.isin(['train_mean', 'seasonal_harmonic'])].pivot(
        index=['fit_reference', 'model', 'record_id'], columns='arm', values='y_pred')
    check('all_mean_and_calendar_control_predictions_identical_BC', np.array_equal(controls[ARMS[0]], controls[ARMS[1]]))
    check('all_control_contrasts_zero', (contrasts.loc[contrasts.model.isin(['train_mean', 'seasonal_harmonic']), 'absolute_difference'] == 0).all())
    check('no_estimator_fit_call_attempts', FIT_CALL_ATTEMPTS == 0)
    check('all_frozen_sources_and_plans_unchanged_after_replay', all(sha(path) == digest for path, digest in freeze_paths.items()))
    check('all_run_files_unchanged_after_replay', all(sha(path) == digest for path, digest in run_paths.items()))
    output.mkdir(exist_ok=True)
    pd.DataFrame(replay_rows).to_csv(output/'SAVED_MODEL_REPLAY.csv', index=False, float_format='%.17g')
    contrasts.to_csv(output/'INDEPENDENT_INPUT_CONTRASTS.csv', index=False, float_format='%.17g')
    primary = contrasts.loc[contrasts.primary].to_dict('records')
    result = {'created_utc': datetime.now(timezone.utc).isoformat(), 'status': 'pass',
        'protocol_sha256': PROTOCOL_SHA, 'verifier_sha256': sha(__file__), 'plan_root': str(plan), 'run_root': str(run),
        'checks_total': len(CHECKS), 'checks_passed': sum(c['passed'] for c in CHECKS), 'checks': CHECKS,
        'physical_train_keys': 480, 'physical_test_keys': 437, 'saved_model_tasks': 20,
        'mean_control_records': 4, 'fitted_estimator_records': 16, 'saved_prediction_replay_rows': 8740,
        'own_prediction_rows': 8740, 'both_score_reference_rows_including_own': 17480,
        'maximum_replay_difference': max(r['maximum_absolute_replay_difference'] for r in replay_rows),
        'maximum_cross_truth_difference': max_cross_error, 'primary_contrasts': primary,
        'frozen_sources_and_plan_files_checked': len(freeze_paths), 'run_files_checked': len(run_paths),
        'new_estimator_fit_calls': FIT_CALL_ATTEMPTS, 'new_predictive_experiments': 0,
        'new_model_or_threshold_selection': False, 'new_Cloud_or_LLM_calls': 0, 'GPU': False,
        'scope': 'Verification of a same-Mead retrospective processed-reference hindcast; no sealed, regional, forecast, causal or agent advantage claim',
        'elapsed_seconds': time.perf_counter()-started}
    with (output/'INDEPENDENT_POSTFIT_VERIFICATION.json').open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: result[k] for k in ['status', 'checks_total', 'physical_test_keys',
        'saved_model_tasks', 'maximum_replay_difference', 'primary_contrasts']}), flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, default=CPU/'plan01')
    parser.add_argument('--run', type=Path, default=CPU/'run01')
    parser.add_argument('--config', type=Path, default=CPU/'code_01/source_config.json')
    parser.add_argument('--authorization', type=Path, default=CPU/'EXECUTION_AUTHORIZATION.json')
    parser.add_argument('--output', type=Path, default=ROOT)
    verify(parser.parse_args())
