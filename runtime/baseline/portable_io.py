"""Portable processed-table I/O and explicit optional native-map preflight."""
from pathlib import Path
import os,hashlib
import pandas as pd

NATIVE_MAP_CHECKS={}
def read_frame(path,columns=None):
    frame=pd.read_csv(path,usecols=columns,float_precision='round_trip')
    if columns is not None:frame=frame.loc[:,columns]
    if 'date' in frame:frame['date']=pd.to_datetime(frame['date'])
    return frame

def optional_map_check(record):
    name=Path(record['path']).name
    expected=record.get('sha256',record.get('expected_sha256'))
    root=os.environ.get('MAIZE_NATIVE_MAP_ROOT')
    if not root:
        NATIVE_MAP_CHECKS[name]={'status':'not_performed_source_rasters_not_packaged','expected_sha256':expected}
        return False
    candidate=Path(root)/name
    if not candidate.is_file():raise FileNotFoundError('Explicit native-map root lacks '+name)
    actual=hashlib.sha256(candidate.read_bytes()).hexdigest()
    if actual!=expected:raise ValueError('Native-map SHA mismatch '+name)
    NATIVE_MAP_CHECKS[name]={'status':'verified_byte_identity_only','expected_sha256':expected,'actual_sha256':actual}
    return True
