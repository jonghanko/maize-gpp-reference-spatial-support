"""Source-pinned, predictors-only early Mead MODIS preparation on CPU WS3.

The four original pure functions are AST-selected without executing the old
module's target-bearing align/build_weights entry points. New outputs only.
"""
from __future__ import annotations
import argparse
import ast
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import threading
from urllib.parse import urlsplit

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
import planetary_computer
import pyproj
import rasterio
from rasterio.windows import Window
import shapely
from shapely.geometry import box

BASE = Path('source://research')
ANALYSIS = BASE / 'analysis_v3_20261026'
OUT = ANALYSIS / 'early_mead_v0_1/optical_inputs'
OLD = ANALYSIS / 'matched_support_pilot_v0_1/spatial_experiment'
SPATIAL = ANALYSIS / 'support_resolution_v0_1/us_spatial'
PRIOR = ANALYSIS / 'candidate_readiness_v0_1'
TARGETS = [('US-Ne1', 2010), ('US-Ne1', 2011), ('US-Ne3', 2011)]
BANDS = ['sur_refl_b01', 'sur_refl_b02', 'sur_refl_b04', 'sur_refl_b06',
         'sur_refl_b07', 'sur_refl_state_500m']

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(2**20), b''):
            h.update(b)
    return h.hexdigest()

def arrsha(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()

def load(p):
    return json.loads(Path(p).read_text(encoding='utf-8'))

def save(p, value):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('x', encoding='utf-8') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')

def now():
    return datetime.now(timezone.utc).isoformat()

def pure_functions(source):
    names = {'cells', 'qa', 'scale', 'interpolate'}
    tree = ast.parse(Path(source).read_text(encoding='utf-8'))
    functions = [x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name in names]
    assert {x.name for x in functions} == names
    module = ast.Module(body=functions, type_ignores=[])
    ns = {'np': np, 'pd': pd, 'pyproj': pyproj, 'shapely': shapely}
    exec(compile(module, str(source), 'exec'), ns)
    return {n: ns[n] for n in names}, {x.name: hashlib.sha256(
        ast.dump(x, include_attributes=False).encode()).hexdigest() for x in functions}

def source_check(frozen):
    return all(Path(p).is_file() and sha(p) == h for p, h in frozen.items())

def freeze():
    OUT.mkdir(parents=True, exist_ok=True)
    assert not (OUT / 'FROZEN_OPTICAL_PROTOCOL.json').exists()
    old_audit = load(PRIOR / 'optical_preparation/audit02/OPTICAL_CACHE_AUDIT.json')
    expected = dict(old_audit['sources_before'])
    expected.update(old_audit['map_hashes'])
    assert source_check(expected), 'Old source/map changed since readiness audit'
    extra = [PRIOR / 'optical_preparation/TERRA_STAC_ITEM_MANIFEST.json',
             PRIOR / 'optical_preparation/metadata_01/STAC_2010_h10v04_response.json',
             PRIOR / 'optical_preparation/metadata_01/STAC_2011_h10v04_response.json',
             OUT / 'frozen_reference/T1_INPUT_PREPARATION_CONTRACT.json',
             OUT / 'frozen_reference/INPUT_REUSE_GUIDE_EN.md',
             SPATIAL / 'run01/monthly_core_overlap.csv', Path(__file__)]
    expected.update({str(p): sha(p) for p in extra})
    functions, function_hashes = pure_functions(OLD / 'code/inputs.py')
    geometry = {s['site_id']: s for s in load(SPATIAL / 'pixel_geometry_01.json')['sites']}
    members = load(SPATIAL / 'acquisition_manifest_01.json')['selected_members']
    previous_overlap = pd.read_csv(SPATIAL / 'run01/monthly_core_overlap.csv')
    terra = load(PRIOR / 'optical_preparation/TERRA_STAC_ITEM_MANIFEST.json')
    items = terra['selected_items']
    assert len(items) == len({i['item_id'] for i in items}) == 92
    stac_items = {}
    for year in [2010, 2011]:
        response_path = PRIOR / f'optical_preparation/metadata_01/STAC_{year}_h10v04_response.json'
        assert sha(response_path) == terra['response_sha256'][response_path.name]
        for item in load(response_path)['features']:
            stac_items[item['id']] = item
    selected = []
    for site, year in TARGETS:
        geom = geometry[site]
        affine = rasterio.Affine(*geom['transform'])
        site_maps = []
        for month in range(1, 13):
            hits = [m for m in members if Path(m['path']).name == f'{site}_{year}_{month:02d}_DAY_fpt_weight.tif']
            assert len(hits) == 1
            member = hits[0]
            assert sha(member['path']) == member['sha256'] == expected[member['path']]
            with rasterio.open(member['path']) as ds:
                a = ds.read(1, masked=True).filled(0).astype(float)
                assert np.isfinite(a).all() and (a >= 0).all()
                rr, cc = np.nonzero(a > 0)
                assert len(rr) > 0
                weights = a[rr, cc]
                polygons = functions['cells'](ds, rr, cc, geom['native_crs_wkt'])
            total = math.fsum(weights.tolist())
            areas = shapely.area(polygons)
            assert np.isfinite(areas).all() and (areas > 0).all()
            bounds = shapely.total_bounds(polygons)
            c0, r1 = (~affine) * (bounds[0], bounds[1])
            c1, r0 = (~affine) * (bounds[2], bounds[3])
            entries = []
            for row in range(math.floor(r0), math.floor(r1) + 1):
                for col in range(math.floor(c0), math.floor(c1) + 1):
                    assert 0 <= row < 2400 and 0 <= col < 2400
                    x0, y0 = affine * (col, row)
                    x1, y1 = affine * (col + 1, row + 1)
                    fraction = shapely.area(shapely.intersection(polygons, box(x0, y1, x1, y0))) / areas
                    mass = math.fsum((weights * fraction).tolist())
                    if mass > 0:
                        entries.append({'row': row, 'col': col, 'stored_core_mass': mass,
                                        'conditional_core_weight': mass / total})
            assert abs(sum(e['conditional_core_weight'] for e in entries) - 1) < 1e-7
            point = next(e for e in entries if e['row'] == geom['row'] and e['col'] == geom['column'])
            old = previous_overlap.loc[previous_overlap.site_id.eq(site) & previous_overlap.year.eq(year)
                & previous_overlap.month.eq(month) & previous_overlap.period.eq('DAY')]
            assert len(old) == 1
            assert abs(point['conditional_core_weight'] - old.iloc[0].fraction_of_published80pct_core_inside_MODIS_pixel) < 1e-7
            site_maps.append({'month': month, 'path': member['path'], 'sha256': member['sha256'],
                              'stored_weight_sum': total, 'pixels': entries})
        pixels = {(e['row'], e['col']) for m in site_maps for e in m['pixels']}
        r0, c0 = min(r for r, c in pixels), min(c for r, c in pixels)
        r1, c1 = max(r for r, c in pixels), max(c for r, c in pixels)
        selected.append({'site_id': site, 'year': year, 'geometry': geom, 'maps': site_maps,
            'window': {'row_off': r0, 'col_off': c0, 'height': r1-r0+1, 'width': c1-c0+1}})
    unions = {}
    for year in [2010, 2011]:
        sites = [s for s in selected if s['year'] == year]
        r0 = min(s['window']['row_off'] for s in sites)
        c0 = min(s['window']['col_off'] for s in sites)
        r1 = max(s['window']['row_off'] + s['window']['height'] for s in sites)
        c1 = max(s['window']['col_off'] + s['window']['width'] for s in sites)
        unions[str(year)] = {'row_off': r0, 'col_off': c0, 'height': r1-r0, 'width': c1-c0}
    for i in items:
        assert i['item_id'].startswith(f'MOD09A1.A{i["year"]}') and '.h10v04.061.' in i['item_id']
        item = stac_items[i['item_id']]
        assert item['properties']['start_datetime'][:10] == i['composite_start']
        assert rasterio.Affine(*item['properties']['proj:transform']).almost_equals(
            rasterio.Affine(*selected[0]['geometry']['transform']), precision=1e-10)
        assert all(item['assets'][b]['href'] == i['required_unsigned_hrefs'][b] for b in BANDS)
    config = {'created_utc': now(), 'status': 'predictors_only_frozen_before_raster_access',
        'recipe': 'Exact original cells/qa/scale/interpolate functions; no target-bearing entry point executed',
        'staged_reference_origin': 'Workstation candidate-readiness T1 contract and input-reuse guide; byte-identical copies in new CPU stage, old CPU release not amended',
        'pure_function_ast_sha256': function_hashes, 'sources_before': expected,
        'adapter_sha256': sha(__file__), 'sites': selected, 'year_union_windows': unions,
        'items': items, 'raw_window_policy': 'One common union window per Terra item; lossless per-site subwindows',
        'bands': BANDS, 'monthly_maps': 36, 'unique_items': 92, 'site_year_windows': 138,
        'point_anchor_policy': 'Reconstructed point must pass same state/five-band QA and finite NDVI before B/C shared screen',
        'core_policy': 'All strictly positive geometric core cells retained; no after-QA removal or renormalization',
        'numeric_weight_policy': 'Original unconditional pre-QA normalization mass/mass.sum() only',
        'date_policy': 'Composite START month; same-year internal linear interpolation with brackets <=24 days; no extrapolation',
        'field_bounds_applied': False, 'A_existing_exact_comparison': False,
        'target_reads': 0, 'model_fits': 0, 'model_execution_allowed': False,
        'coordinate_changes': 0, 'Ne1_period_valid_registration_authority_resolved': False,
        'full_daily_footprint_claim': False, 'new_region_or_sealed_claim': False}
    assert source_check(expected)
    save(OUT / 'FROZEN_OPTICAL_PROTOCOL.json', config)
    print(json.dumps({'status': config['status'], 'map_count': 36, 'items': 92,
        'union_windows': unions, 'monthly_pixel_counts': [{'site': s['site_id'], 'year': s['year'],
        'counts': [len(m['pixels']) for m in s['maps']]} for s in selected]}), flush=True)

def acquire_item(item, config, stop):
    if stop.is_set():
        return None
    item_id = item['item_id']
    directory = OUT / 'raw_windows' / str(item['year'])
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / (item_id + '.npz')
    meta = directory / (item_id + '.json')
    assert not path.exists() and not meta.exists(), 'Do not silently replace an acquired window'
    win = config['year_union_windows'][str(item['year'])]
    window = Window(win['col_off'], win['row_off'], win['width'], win['height'])
    geom = config['sites'][0]['geometry']
    expected = rasterio.Affine(*geom['transform'])
    arrays, assets = {}, []
    for band in BANDS:
        if stop.is_set():
            return None
        href = item['required_unsigned_hrefs'][band]
        assert urlsplit(href).hostname == 'modiseuwest.blob.core.windows.net'
        # The provider's ordinary public SAS signing method is inherited from acquire_one.
        # A denial stops the campaign; no alternate hostname/identity/route is attempted.
        signed = planetary_computer.sign(href)
        with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR', GDAL_HTTP_MAX_RETRY='0',
                GDAL_HTTP_TIMEOUT='45', CPL_VSIL_CURL_ALLOWED_EXTENSIONS='.tif'), rasterio.open(signed) as ds:
            assert ds.transform.almost_equals(expected, precision=1e-10)
            assert ds.crs == rasterio.crs.CRS.from_wkt(geom['native_crs_wkt'])
            assert (ds.height, ds.width) == (2400, 2400) and ds.count == 1
            a = ds.read(1, window=window)
            assert a.shape == (win['height'], win['width'])
            assert a.dtype.kind in 'iu'
            arrays[band] = a
            assets.append({'band': band, 'unsigned_url': href, 'array_sha256': arrsha(a),
                'dtype': str(a.dtype), 'shape': list(a.shape), 'raster_transform': list(ds.transform)[:6],
                'raster_crs_wkt': ds.crs.to_wkt(), 'raster_shape': [ds.height, ds.width],
                'affine_crs_shape_verified': True, 'window': win, 'read_success': True})
    with path.open('xb') as f:
        np.savez_compressed(f, **arrays)
    result = {'item_id': item_id, 'year': item['year'], 'anchor': item['composite_start'],
        'path': str(path), 'sha256': sha(path), 'bytes': path.stat().st_size, 'window': win,
        'assets': assets, 'frozen_item_manifest_sha256': config['sources_before'][str(PRIOR / 'optical_preparation/TERRA_STAC_ITEM_MANIFEST.json')]}
    save(meta, result)
    return result

def acquire():
    config = load(OUT / 'FROZEN_OPTICAL_PROTOCOL.json')
    assert source_check(config['sources_before'])
    assert sha(__file__) == config['adapter_sha256']
    assert not (OUT / 'ACQUISITION_RECEIPT.json').exists()
    results, failures = [], []
    stop = threading.Event()
    # Probe the exact planned first item before parallelizing any public reads.
    remaining = list(config['items'])
    first = remaining.pop(0)
    try:
        results.append(acquire_item(first, config, stop))
        print(json.dumps({'completed': 1, 'total': 92, 'failures': 0}), flush=True)
    except Exception as exc:
        stop.set()
        failures.append({'item_id': first['item_id'], 'exception_type': type(exc).__name__,
            'denial_403': '403' in str(exc), 'error_detail_redacted': True,
            'error_classification': 'ordinary_provider_access_failed_no_bypass'})
    if not stop.is_set():
        with ThreadPoolExecutor(max_workers=4) as pool:
            futures = {pool.submit(acquire_item, item, config, stop): item for item in remaining}
            for future in as_completed(futures):
                item = futures[future]
                try:
                    value = future.result()
                    if value is not None:
                        results.append(value)
                except Exception as exc:
                    stop.set()
                    failures.append({'item_id': item['item_id'], 'exception_type': type(exc).__name__,
                        'denial_403': '403' in str(exc), 'error_detail_redacted': True,
                        'error_classification': 'ordinary_provider_access_failed_no_bypass'})
                if len(results) % 5 == 0 or failures:
                    print(json.dumps({'completed': len(results), 'total': 92, 'failures': len(failures)}), flush=True)
    results.sort(key=lambda x: (x['year'], x['anchor']))
    assert source_check(config['sources_before'])
    receipt = {'created_utc': now(), 'status': 'passed' if len(results) == 92 and not failures else 'blocked_incomplete',
        'scenes': results, 'failures': failures, 'expected_unique_Terra_items': 92,
        'complete_unique_items': len(results), 'complete_asset_reads': 6 * len(results),
        'raw_window_file_bytes': sum(r['bytes'] for r in results),
        'http_request_count': None, 'upstream_transfer_bytes': None,
        'bytes_note': 'Local compressed native-window bytes only; COG HTTP tile/header traffic is not instrumented',
        'method': 'Frozen unsigned hrefs, ordinary planetary_computer.sign, rasterio native range-window reads',
        'denial_bypass_attempts': 0, 'new_metadata_searches': 0, 'target_reads': 0, 'model_fits': 0,
        'source_hashes_unchanged': True, 'protocol_sha256': sha(OUT / 'FROZEN_OPTICAL_PROTOCOL.json')}
    save(OUT / 'ACQUISITION_RECEIPT.json', receipt)
    print(json.dumps({k: receipt[k] for k in ['status', 'complete_unique_items', 'complete_asset_reads', 'raw_window_file_bytes', 'failures']}), flush=True)
    if receipt['status'] != 'passed':
        raise RuntimeError('Ordinary provider acquisition incomplete; campaign stopped without bypass')

def ndvi(red, nir):
    denominator = float(nir) + float(red)
    return float((nir-red) / denominator) if np.isfinite(denominator) and denominator > 0 else math.nan

def scene_values(data, site, date, config, functions):
    win = config['year_union_windows'][str(site['year'])]
    pixels = site['maps'][pd.Timestamp(date).month - 1]['pixels']
    refl = {b: functions['scale'](data[b]) for b in BANDS[:-1]}
    usable = functions['qa'](data[BANDS[-1]])
    for a in refl.values():
        usable &= np.isfinite(a)
    pr = site['geometry']['row'] - win['row_off']
    pc = site['geometry']['column'] - win['col_off']
    red, nir = refl['sur_refl_b01'], refl['sur_refl_b02']
    point_usable = bool(usable[pr, pc])
    point = ndvi(red[pr, pc], nir[pr, pc]) if point_usable else math.nan
    point_usable = bool(point_usable and np.isfinite(point))
    rr = np.asarray([p['row']-win['row_off'] for p in pixels])
    cc = np.asarray([p['col']-win['col_off'] for p in pixels])
    mass = np.asarray([p['conditional_core_weight'] for p in pixels])
    assert abs(mass.sum()-1) < 1e-7
    # Inherited unconditional floating-point normalization, before any QA selection.
    mass = mass / mass.sum()
    core_usable = bool(usable[rr, cc].all())
    shared = point_usable and core_usable
    b = c = math.nan
    if shared:
        b = ndvi(red[rr, cc].mean(), nir[rr, cc].mean())
        c = ndvi(float(mass @ red[rr, cc]), float(mass @ nir[rr, cc]))
        shared = bool(np.isfinite(b) and np.isfinite(c))
    w = site['window']
    r0, c0 = w['row_off']-win['row_off'], w['col_off']-win['col_off']
    per_site_hashes = {band: arrsha(data[band][r0:r0+w['height'], c0:c0+w['width']]) for band in BANDS}
    return {'site_id': site['site_id'], 'year': site['year'], 'anchor': date,
        'map_month': pd.Timestamp(date).month, 'positive_core_pixels': len(pixels),
        'point_anchor_usable': point_usable, 'positive_core_qa_usable': core_usable,
        'shared_anchor_usable': shared, 'required_pixel_failures': int((~usable[rr, cc]).sum()),
        'point_source_ndvi': float(point) if point_usable else None,
        'B_uniform_pixels': float(b) if shared else None,
        'C_core_weighted': float(c) if shared else None,
        'site_native_window': w, 'site_native_array_sha256': per_site_hashes}

def derive():
    config = load(OUT / 'FROZEN_OPTICAL_PROTOCOL.json')
    acq = load(OUT / 'ACQUISITION_RECEIPT.json')
    assert acq['status'] == 'passed' and len(acq['scenes']) == 92
    assert source_check(config['sources_before']) and sha(__file__) == config['adapter_sha256']
    functions, hashes = pure_functions(OLD / 'code/inputs.py')
    assert hashes == config['pure_function_ast_sha256']
    scenes, daily, audit = [], [], []
    for site in config['sites']:
        comp = []
        for item in [i for i in acq['scenes'] if i['year'] == site['year']]:
            assert sha(item['path']) == item['sha256']
            with np.load(item['path'], allow_pickle=False) as data:
                assert sorted(data.files) == sorted(BANDS)
                assert all(arrsha(data[b]) == a['array_sha256'] for b, a in zip(BANDS, item['assets']))
                row = scene_values(data, site, item['anchor'], config, functions)
            row.update(item_id=item['item_id'], raw_path=item['path'], raw_sha256=item['sha256'])
            scenes.append(row)
            comp.append(row)
        comp = pd.DataFrame(comp).sort_values('anchor')
        dates = pd.date_range(f'{site["year"]}-01-01', f'{site["year"]}-12-31')
        frame = pd.DataFrame({'site_id': site['site_id'], 'year': site['year'], 'date': dates})
        anchors = pd.to_datetime(comp.anchor)
        for arm in ['B_uniform_pixels', 'C_core_weighted']:
            frame[arm] = functions['interpolate'](anchors, pd.to_numeric(comp[arm], errors='coerce'), dates)
        frame['complete_BC_optical'] = np.isfinite(frame[['B_uniform_pixels', 'C_core_weighted']]).all(axis=1)
        assert np.isfinite(frame.B_uniform_pixels).equals(np.isfinite(frame.C_core_weighted))
        keep = frame.date.dt.month.between(5, 10)
        assert int(keep.sum()) == 184
        audit.append({'site_id': site['site_id'], 'year': site['year'], 'anchor_candidates': len(comp),
            'point_usable_anchors': int(comp.point_anchor_usable.sum()),
            'all_positive_core_qa_anchors': int(comp.positive_core_qa_usable.sum()),
            'shared_usable_anchors': int(comp.shared_anchor_usable.sum()),
            'MayOct_calendar_days': 184, 'MayOct_complete_BC_optical_days': int(frame.loc[keep, 'complete_BC_optical'].sum()),
            'full_calendar_days': len(frame), 'full_calendar_complete_BC_optical_days': int(frame.complete_BC_optical.sum())})
        daily.append(frame)
    combined = pd.concat(daily, ignore_index=True)
    assert len(combined) == 1095 and not combined.duplicated(['site_id', 'date']).any()
    combined.to_parquet(OUT / 'daily_BC_optical_full_year.parquet', index=False)
    combined.loc[combined.date.dt.month.between(5, 10)].to_csv(OUT / 'daily_BC_optical_MayOct.csv', index=False, date_format='%Y-%m-%d')
    pd.DataFrame(audit).to_csv(OUT / 'OPTICAL_COVERAGE_SUMMARY.csv', index=False)
    save(OUT / 'SCENE_QA_AUDIT.json', {'status': 'completed', 'scenes': scenes})
    source_unchanged = source_check(config['sources_before'])
    assert source_unchanged
    outputs = {p.name: sha(p) for p in [OUT / 'daily_BC_optical_full_year.parquet', OUT / 'daily_BC_optical_MayOct.csv',
        OUT / 'OPTICAL_COVERAGE_SUMMARY.csv', OUT / 'SCENE_QA_AUDIT.json']}
    save(OUT / 'OPTICAL_CONSTRUCTION_RECEIPT.json', {'created_utc': now(), 'status': 'passed',
        'site_year_window_candidates': len(scenes), 'unique_Terra_items': len(acq['scenes']),
        'native_raster_assets_verified': 552, 'monthly_maps_verified': 36, 'calendar_rows': len(combined),
        'MayOct_rows': 552, 'coverage': audit, 'outputs': outputs,
        'protocol_sha256': sha(OUT / 'FROZEN_OPTICAL_PROTOCOL.json'),
        'acquisition_sha256': sha(OUT / 'ACQUISITION_RECEIPT.json'),
        'source_hashes_before': config['sources_before'], 'source_hashes_unchanged': source_unchanged,
        'python': sys.version, 'platform': platform.platform(),
        'libraries': {'numpy': np.__version__, 'pandas': pd.__version__, 'rasterio': rasterio.__version__,
            'shapely': shapely.__version__, 'pyproj': pyproj.__version__},
        'target_reads': 0, 'weather_joins': 0, 'model_fits': 0, 'new_scores': 0,
        'model_execution_allowed': False, 'coordinate_changes': 0,
        'scope': 'Earlier-year same-region retrospective conditional DAY-core B/C predictors only; field-window/native target QC and fit protocol separate',
        'A_point_comparison_approved': False, 'Ne1_registration_resolved': False,
        'full_daily_flux_footprint_claim': False, 'sealed_or_new_region_claim': False})
    print(json.dumps({'status': 'passed', 'coverage': audit, 'outputs': outputs}), flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=['freeze', 'acquire', 'derive'])
    args = parser.parse_args()
    {'freeze': freeze, 'acquire': acquire, 'derive': derive}[args.phase]()
