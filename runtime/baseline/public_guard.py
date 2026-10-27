"""Offline analysis I/O guard, independent of private source-server locations."""
from pathlib import Path
import os,sys
DATA_SUFFIXES=('.json','.csv','.csv.gz','.parquet','.joblib','.tif','.tiff','.npy','.npz','.yaml','.yml')
def install(roots):
    allow=[os.path.realpath(p).rstrip(os.sep)+os.sep for p in roots]+[os.path.realpath(sys.prefix).rstrip(os.sep)+os.sep,os.path.realpath(sys.base_prefix).rstrip(os.sep)+os.sep]
    def hook(event,args):
        if event in ('socket.connect','socket.getaddrinfo'):
            raise PermissionError('Network access is disabled during offline analysis')
        if event!='open' or not args or not isinstance(args[0],(str,bytes,os.PathLike)):return
        p=os.path.realpath(os.fsdecode(args[0]));allowed=any(p==a[:-1] or p.startswith(a) for a in allow)
        if not allowed and p.lower().endswith(DATA_SUFFIXES):
            raise PermissionError('Data/config access outside public bundle/output/runtime is blocked')
    sys.addaudithook(hook)
