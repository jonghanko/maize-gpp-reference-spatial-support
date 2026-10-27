"""Build four prespecified optical regimes on frozen same-year core/QC support."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
from urllib.request import urlopen
import numpy as np
import pandas as pd
import planetary_computer
import pyproj
import rasterio
from rasterio.windows import Window
import shapely
from shapely.geometry import box

FEATURES=['NDVI_modis','EA_kPa_agera5','WS_ms_agera5','P_mm_agera5','TA_max_C_agera5','TA_min_C_agera5','SW_IN_MJ_m2_d_agera5','VPD_kPa_agera5']
ARMS=['A_existing','A_shared_scene','B_uniform_pixels','C_core_weighted']
BANDS=['sur_refl_b01','sur_refl_b02','sur_refl_b04','sur_refl_b06','sur_refl_b07','sur_refl_state_500m']
SITE_YEARS=[('US-Ne1',2012),('US-Ne1',2013),('US-Ne3',2013),('US-Ro1',2015)]

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()

def save(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False)

def roots(root):
    base=root.parents[1]
    return base,base/'support_resolution_v0_1/us_spatial',base/'us_baseline_v0_1/runs/gpp_vut/data_snapshot.parquet'

def cells(ds,rr,cc,target):
    offsets=np.asarray([(0,0),(.25,0),(.5,0),(.75,0),(1,0),(1,.25),(1,.5),(1,.75),(1,1),(.75,1),(.5,1),(.25,1),(0,1),(0,.75),(0,.5),(0,.25)])
    x=ds.transform.c+(cc[:,None]+offsets[None,:,0])*ds.transform.a
    y=ds.transform.f+(rr[:,None]+offsets[None,:,1])*ds.transform.e
    xx,yy=pyproj.Transformer.from_crs(ds.crs,target,always_xy=True).transform(x,y)
    return shapely.polygons(np.stack([xx,yy],axis=-1))

def build_weights(root):
    base,spatial,source=roots(root)
    geometry=json.loads((spatial/'pixel_geometry_01.json').read_text());geometries={s['site_id']:s for s in geometry['sites']}
    maps=json.loads((spatial/'acquisition_manifest_01.json').read_text())
    discovery=json.loads((spatial/'source_discovery_01.json').read_text())
    prior=pd.read_csv(spatial/'run01/monthly_core_overlap.csv')
    selected=[];sources={str(source):sha(source),str(root/'code/INPUT_PROTOCOL.json'):sha(root/'code/INPUT_PROTOCOL.json')}
    for site,year in SITE_YEARS:
        geometry=geometries[site];affine=rasterio.Affine(*geometry['transform'])
        attrs=[a for a in discovery['attrs'] if a['metadata']['site_id']==site and a['metadata']['year']==year]
        assert len(attrs)==1;attr=attrs[0]
        sources[attr['path']]=sha(attr['path']);sources[attr['metadata']['cache_path']]=sha(attr['metadata']['cache_path'])
        site_maps=[]
        for month in range(1,13):
            member=next(m for m in maps['selected_members'] if Path(m['path']).name==f'{site}_{year}_{month:02d}_DAY_fpt_weight.tif')
            sources[member['path']]=sha(member['path'])
            with rasterio.open(member['path']) as ds:
                a=ds.read(1,masked=True).filled(0).astype(float);assert np.isfinite(a).all() and (a>=0).all()
                rr,cc=np.nonzero(a>0);weights=a[rr,cc];polygons=cells(ds,rr,cc,geometry['native_crs_wkt'])
            total=math.fsum(weights.tolist());areas=shapely.area(polygons);bounds=shapely.total_bounds(polygons)
            c0,r1=(~affine)*(bounds[0],bounds[1]);c1,r0=(~affine)*(bounds[2],bounds[3]);entries=[]
            for row in range(math.floor(r0),math.floor(r1)+1):
                for col in range(math.floor(c0),math.floor(c1)+1):
                    assert 0<=row<2400 and 0<=col<2400
                    x0,y0=affine*(col,row);x1,y1=affine*(col+1,row+1)
                    fraction=shapely.area(shapely.intersection(polygons,box(x0,y1,x1,y0)))/areas
                    mass=math.fsum((weights*fraction).tolist())
                    if mass>0:entries.append({'row':row,'col':col,'stored_core_mass':mass,'conditional_core_weight':mass/total})
            assert abs(sum(e['conditional_core_weight'] for e in entries)-1)<1e-7
            point=next(e for e in entries if e['row']==geometry['row'] and e['col']==geometry['column'])
            old=prior.loc[prior.site_id.eq(site)&prior.year.eq(year)&prior.month.eq(month)&prior.period.eq('DAY')].iloc[0]
            assert abs(point['conditional_core_weight']-old.fraction_of_published80pct_core_inside_MODIS_pixel)<1e-7
            site_maps.append({'month':month,'path':member['path'],'sha256':member['sha256'],'stored_weight_sum':total,'pixels':entries})
        allpixels={(e['row'],e['col']) for m in site_maps for e in m['pixels']}
        rmin,cmin=min(r for r,c in allpixels),min(c for r,c in allpixels)
        rmax,cmax=max(r for r,c in allpixels),max(c for r,c in allpixels)
        selected.append({'site_id':site,'year':year,'geometry':geometry,'attrs_path':attr['path'],'attrs_sha256':attr['sha256'],'cache_path':attr['metadata']['cache_path'],
                         'granule_ids':attr['metadata']['granule_ids'],'qa_policy':attr['metadata']['qa_policy'],'maps':site_maps,
                         'window':{'row_off':rmin,'col_off':cmin,'height':rmax-rmin+1,'width':cmax-cmin+1}})
    save(root/'weights.json',{'created_utc':datetime.now(timezone.utc).isoformat(),'sites':selected,'sources':sources,'code_sha256':sha(__file__),
                            'source_policy':'Full positive retained-core MODIS pixel set; no after-QA pixel removal. Equal native pixel area in B; conditional core mass in C.'})
    print(json.dumps({'status':'weights_frozen','site_years':len(selected),'monthly_maps':48,'pixel_counts':[{ 'site':s['site_id'],'year':s['year'],'counts':[len(m['pixels']) for m in s['maps']],'window':s['window']} for s in selected]}),flush=True)

def acquire_one(root,site,item_id):
    directory=root/'raw_windows'/site['site_id']/str(site['year']);directory.mkdir(parents=True,exist_ok=True)
    path=directory/(item_id+'.npz');meta=directory/(item_id+'.json')
    if path.exists() and meta.exists():return json.loads(meta.read_text())
    url='https://planetarycomputer.microsoft.com/api/stac/v1/collections/modis-09A1-061/items/'+item_id
    with urlopen(url,timeout=60) as response:raw=response.read()
    item=json.loads(raw);assert item['id']==item_id
    metadata_path=directory/(item_id+'_STAC.json');metadata_path.write_bytes(raw)
    props=item['properties'];expected=rasterio.Affine(*site['geometry']['transform'])
    assert rasterio.Affine(*props['proj:transform']).almost_equals(expected,precision=1e-10)
    win=site['window'];window=Window(win['col_off'],win['row_off'],win['width'],win['height'])
    arrays={};assets=[]
    for band in BANDS:
        href=item['assets'][band]['href']
        with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',GDAL_HTTP_MAX_RETRY='4',GDAL_HTTP_RETRY_DELAY='1',GDAL_HTTP_TIMEOUT='60'),rasterio.open(planetary_computer.sign(href)) as ds:
            assert ds.transform.almost_equals(expected,precision=1e-10) and ds.crs==rasterio.crs.CRS.from_wkt(site['geometry']['native_crs_wkt'])
            arrays[band]=ds.read(1,window=window)
            assert arrays[band].shape==(win['height'],win['width'])
            assets.append({'band':band,'unsigned_url':href,'array_sha256':hashlib.sha256(arrays[band].tobytes()).hexdigest(),'dtype':str(arrays[band].dtype)})
    with path.open('xb') as f:np.savez_compressed(f,**arrays)
    result={'site_id':site['site_id'],'year':site['year'],'item_id':item_id,'anchor':props['start_datetime'][:10],
            'path':str(path),'sha256':sha(path),'bytes':path.stat().st_size,'STAC_path':str(metadata_path),'STAC_sha256':sha(metadata_path),'assets':assets,'window':win}
    save(meta,result);return result

def acquire(root):
    config=json.loads((root/'weights.json').read_text());tasks=[(site,g) for site in config['sites'] for g in site['granule_ids']]
    results=[];failures=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        jobs={pool.submit(acquire_one,root,s,g):(s['site_id'],g) for s,g in tasks}
        for n,future in enumerate(as_completed(jobs),1):
            try:results.append(future.result())
            except Exception as exc:failures.append({'site':jobs[future][0],'granule':jobs[future][1],'error':str(exc)})
            if n%20==0 or n==len(tasks):print(json.dumps({'completed':n,'total':len(tasks),'failures':len(failures)}),flush=True)
    results.sort(key=lambda r:(r['site_id'],r['year'],r['anchor']))
    save(root/'acquisition_manifest.json',{'status':'passed' if not failures else 'failed','created_utc':datetime.now(timezone.utc).isoformat(),'scenes':results,'failures':failures,'attempted':len(tasks),
                                         'assets_per_scene':BANDS,'raw_window_bytes':sum(r['bytes'] for r in results)})
    if failures:raise ValueError('Source acquisition failures require resolution before alignment')

def qa(word):
    word=np.asarray(word,dtype=np.int64)
    return (word!=65535)&~np.isin(word&3,[1,2])&(((word>>2)&1)==0)&(((word>>8)&3)<1)&(((word>>10)&1)==0)&(((word>>12)&1)==0)&(((word>>15)&1)==0)

def scale(dn):
    raw=np.asarray(dn,dtype=float);out=np.where(np.isfinite(raw)&(raw!=-28672)&(raw>=-100)&(raw<=16000),raw*.0001,np.nan)
    return np.where(out>1,np.nan,np.where(out<0,0,out))

def interpolate(anchors,values,dates):
    valid=np.isfinite(values);anchors=pd.DatetimeIndex(anchors)[valid];values=np.asarray(values)[valid]
    result=np.full(len(dates),np.nan)
    for i,date in enumerate(dates):
        j=anchors.searchsorted(date,side='right')-1
        if j>=0 and date==anchors[j]:result[i]=values[j]
        elif j>=0 and j+1<len(anchors) and (anchors[j+1]-anchors[j]).days<=24:
            ratio=(date-anchors[j]).days/(anchors[j+1]-anchors[j]).days
            result[i]=(1-ratio)*values[j]+ratio*values[j+1]
    return result

def align(root):
    base,spatial,source=roots(root);settings=json.loads((root/'weights.json').read_text());acq=json.loads((root/'acquisition_manifest.json').read_text());assert acq['status']=='passed'
    original=pd.read_parquet(source);initial=original.loc[[ (s,int(y)) in SITE_YEARS for s,y in zip(original.site_id,original.harvest_year)]].copy();assert len(initial)==692
    daily=[];scenes=[];sources={str(source):sha(source),str(root/'weights.json'):sha(root/'weights.json'),str(root/'acquisition_manifest.json'):sha(root/'acquisition_manifest.json')}
    for site in settings['sites']:
        s,y=site['site_id'],site['year'];cache=pd.read_parquet(site['cache_path']);cache['time']=pd.to_datetime(cache.time)
        composite=[]
        for item in [r for r in acq['scenes'] if r['site_id']==s and r['year']==y]:
            assert sha(item['path'])==item['sha256'];sources[item['path']]=item['sha256'];sources[item['STAC_path']]=item['STAC_sha256']
            date=pd.Timestamp(item['anchor']);weights=site['maps'][date.month-1]['pixels'];win=site['window']
            with np.load(item['path'],allow_pickle=False) as data:
                reflectance={band:scale(data[band]) for band in BANDS[:-1]};usable=qa(data[BANDS[-1]])
            for values in reflectance.values():usable &=np.isfinite(values)
            pr=site['geometry']['row']-win['row_off'];pc=site['geometry']['column']-win['col_off'];assert usable[pr,pc]
            red,nir=reflectance['sur_refl_b01'],reflectance['sur_refl_b02'];point=(nir[pr,pc]-red[pr,pc])/(nir[pr,pc]+red[pr,pc])
            old=cache.loc[cache.time.eq(date),'NDVI_modis'];assert len(old)==1
            if np.isfinite(old.iloc[0]):assert abs(point-float(old.iloc[0]))<1e-12
            rr=np.array([w['row']-win['row_off'] for w in weights]);cc=np.array([w['col']-win['col_off'] for w in weights]);mass=np.array([w['conditional_core_weight'] for w in weights])
            assert abs(mass.sum()-1)<1e-7;mass=mass/mass.sum()
            shared=bool(usable[rr,cc].all() and np.isfinite(point))
            b=c=math.nan
            if shared:
                red_b,nir_b=red[rr,cc].mean(),nir[rr,cc].mean();red_c,nir_c=float(mass@red[rr,cc]),float(mass@nir[rr,cc])
                b=(nir_b-red_b)/(nir_b+red_b);c=(nir_c-red_c)/(nir_c+red_c)
                shared=bool(np.isfinite(b) and np.isfinite(c))
            row={'site_id':s,'year':y,'anchor':item['anchor'],'item_id':item['item_id'],'positive_core_pixels':len(weights),'shared_anchor_usable':shared,
                 'point_source_ndvi':float(point) if np.isfinite(point) else None,'A_shared_scene':float(point) if shared else None,
                 'B_uniform_pixels':float(b) if shared else None,'C_core_weighted':float(c) if shared else None,
                 'required_pixel_failures':int((~usable[rr,cc]).sum()),'raw_path':item['path'],'raw_sha256':item['sha256']}
            scenes.append(row);composite.append(row)
        dates=pd.date_range(f'{y}-01-01',f'{y}-12-31');frame=pd.DataFrame({'site_id':s,'date':dates,'year':y})
        comp=pd.DataFrame(composite).sort_values('anchor')
        for arm in ARMS[1:]:frame[arm]=interpolate(pd.to_datetime(comp.anchor),pd.to_numeric(comp[arm],errors='coerce'),dates)
        daily.append(frame)
    daily=pd.concat(daily,ignore_index=True);assert len(daily)==1461 and not daily.duplicated(['site_id','date']).any()
    initial['date']=pd.to_datetime(initial.date)
    common=initial.merge(daily[['site_id','date',*ARMS[1:]]],on=['site_id','date'],how='left',validate='one_to_one')
    common['A_existing']=common.NDVI_modis
    keep=np.isfinite(common[ARMS+FEATURES[1:]+['gpp_gC_m2_d']].to_numpy(float)).all(axis=1)
    rejected=common.loc[~keep,['record_id','site_id','date','harvest_year',*ARMS]]
    common=common.loc[keep].sort_values('record_id').reset_index(drop=True)
    assert not common.empty
    assert common.training_release_approved.eq(False).all() and common.scientific_primary_approved.eq(False).all()
    output=root/'inputs';output.mkdir(exist_ok=False)
    daily.to_parquet(output/'daily_new_optical_curves.parquet',index=False)
    pd.DataFrame(scenes).to_parquet(output/'scene_audit.parquet',index=False)
    common.to_parquet(output/'common_cohort.parquet',index=False)
    rejected.to_parquet(output/'excluded_candidates.parquet',index=False)
    summary=common.groupby(['site_id','spatial_cluster','harvest_year']).size().rename('rows').reset_index();summary.to_csv(output/'cohort_summary.csv',index=False)
    for p,h in settings['sources'].items():assert sha(p)==h
    save(output/'input_receipt.json',{'status':'passed','created_utc':datetime.now(timezone.utc).isoformat(),'initial_cases':len(initial),'common_cases':len(common),'excluded_cases':len(rejected),
                                     'scenes':len(scenes),'shared_usable_scenes':sum(r['shared_anchor_usable'] for r in scenes),'source_point_samples_verified':len(scenes),
                                     'source_hashes':sources,'scope':'Four arms, same target cases/weather; shared scenes for A_shared/B/C; A_existing exact source; no actual-season/full-footprint approval',
                                     'outputs':{p.name:sha(p) for p in sorted(output.iterdir()) if p.is_file()},'input_protocol_sha256':sha(root/'code/INPUT_PROTOCOL.json')})
    print(summary.to_json(orient='records'));print(json.dumps({'initial':len(initial),'common':len(common),'excluded':len(rejected),'scenes':len(scenes),'shared_scenes':sum(r['shared_anchor_usable'] for r in scenes)}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('phase',choices=['weights','acquire','align']);p.add_argument('--root',type=Path,required=True);a=p.parse_args()
    {'weights':build_weights,'acquire':acquire,'align':align}[a.phase](a.root)
