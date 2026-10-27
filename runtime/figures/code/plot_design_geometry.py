"""Publication illustrations from frozen cohort definitions and native geometry."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,hashlib,math,platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,Rectangle,FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from pyproj import Transformer,CRS,Geod
import rasterio

ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def box(ax,x,y,w,h,text,face='#F3F7FA',size=8.4):
    p=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.006,rounding_size=0.015',facecolor=face,edgecolor='#718899',lw=.8)
    ax.add_patch(p);ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=size,linespacing=1.4);return p
def arrow(ax,x1,y1,x2,y2,color='#4D6576'):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=10,lw=.9,color=color))
def export(fig,out,stem):
    outputs=[]
    for ext in ['png','svg','pdf']:
        p=out/(stem+'.'+ext);fig.savefig(p,dpi=320,facecolor='white',metadata={'Creator':'Matplotlib; exact sources and plot code are archived'})
        outputs.append({'file':p.name,'bytes':p.stat().st_size,'sha256':sha(p)})
    plt.close(fig);return outputs

def design(out):
    fig=plt.figure(figsize=(7.4,7.6));gs=fig.add_gridspec(2,1,height_ratios=[1.07,1],left=.04,right=.98,bottom=.025,top=.98,hspace=.13)
    a=fig.add_subplot(gs[0]);a.set(xlim=(0,1),ylim=(0,1));a.axis('off')
    a.text(.0,.99,'a  Frozen common-case experiments',fontsize=11,weight='bold',va='top')
    box(a,.02,.66,.46,.24,'Principal: 642 daily cases\n4 site-seasons · 3 towers · 2 regions\nPrimary LOYO: equal weight to 4 seasons\nOriginal supported fold roles retained',size=8.6)
    box(a,.53,.66,.45,.24,'Separate Mead hindcast: 437 test dates\nTrain: 480 Mead dates in 2012/2013\nTest: 2010/2011 at already-seen sites\nB/C only; later-year training labels',size=8.6)
    a.plot([.25,.25,.75,.75],[.65,.61,.61,.65],color='#4D6576',lw=.9)
    a.plot([.25,.125],[.61,.61],color='#4D6576',lw=.9);arrow(a,.125,.61,.125,.575)
    box(a,.02,.27,.215,.30,'Optical representations\nPrincipal: four arms\nA: existing point\nA: shared-scene point\nB: uniform; C: core\nT1: B/C only',size=8.2)
    box(a,.275,.27,.215,.30,'Seven AgERA5 inputs\nSame cases, weather,\nfold roles and recipes\nAll five fixed families\nOriginal priorities retained',size=8.2)
    box(a,.53,.27,.19,.30,'Fitting reference\nVUT or CUT\nSeparate pipelines\nfor each input arm',size=8.5)
    box(a,.76,.27,.22,.30,'Scoring reference\nScore each prediction\nagainst VUT and CUT\nShared observations;\ndependent products',size=8.2)
    for x1,x2 in [(.237,.274),(.491,.529),(.722,.758)]:arrow(a,x1,.42,x2,.42)
    a.text(.5,.14,'B/C use identical cells and accepted scenes; cases and weather are fixed.\nEqual-season and equal-region pooled summaries remain distinct.',ha='center',va='center',fontsize=8.5,linespacing=1.45)
    b=fig.add_subplot(gs[1]);b.set(xlim=(0,1),ylim=(0,1));b.axis('off')
    b.text(.0,.99,'b  Post-review cross-prediction from unchanged saved models',fontsize=11,weight='bold',va='top')
    b.text(.455,.84,'TEST INPUT',ha='center',fontsize=8.6,weight='bold')
    b.text(.335,.77,'B: uniform',ha='center',fontsize=8.5);b.text(.58,.77,'C: core weighted',ha='center',fontsize=8.5)
    b.text(.08,.46,'TRAINING INPUT',rotation=90,ha='center',va='center',fontsize=8.6,weight='bold')
    b.text(.15,.60,'B',ha='center',va='center',fontsize=10,weight='bold');b.text(.15,.35,'C',ha='center',va='center',fontsize=10,weight='bold')
    for x,y,label,detail,diag in [(.22,.49,'$f_B(X_B)$  /  $R_{BB}$','Original diagonal',True),(.465,.49,'$f_B(X_C)$  /  $R_{BC}$','New cross-prediction',False),(.22,.24,'$f_C(X_B)$  /  $R_{CB}$','New cross-prediction',False),(.465,.24,'$f_C(X_C)$  /  $R_{CC}$','Original diagonal',True)]:
        box(b,x,y,.23,.215,label+'\n'+detail,face='#DDEBF2' if diag else '#FFF0D8',size=8.8)
    box(b,.755,.24,.225,.465,'100 principal + 20 T1\nsaved B/C models\n\n240 prediction passes\n0 new fits\n\nEach fit/score path\nis reported separately',size=8.1)
    b.text(.48,.115,'Conditional training and test-input contrasts + interaction',ha='center',fontsize=8.8,weight='bold')
    b.text(.5,.025,'Same cases and reference paths · Exposed historical domains\nComputational sensitivity; no physical causal attribution or independent confirmation',ha='center',va='center',fontsize=8.15,linespacing=1.4)
    return export(fig,out,'Figure1_controlled_design')

def native_xy(record):
    g=record['geometry'];tr=Transformer.from_crs('EPSG:4326',CRS.from_wkt(g['native_crs_wkt']),always_xy=True)
    a=g['registered_coordinate_native'];calculated=tr.transform(g['lon'],g['lat'])
    assert np.max(np.abs(np.array(a)-np.array(calculated)))<1e-6
    o=record['coordinate_comparison'];published=tr.transform(o['published_lon'],o['published_lat'])
    return tr,np.array(a),np.array(published)

def draw_map(ax,record,base,norm,cmap,zoom=False):
    g=record['geometry'];tr,origin,published=native_xy(record)
    tf=g['transform'];sx,sy=tf[0],tf[4];x0,y0=tf[2],tf[5]
    r0,c0=g['row'],g['column'];coords=[]
    pixels=record['map']['pixels'];weights={(v['row'],v['col']):v['conditional_core_weight'] for v in pixels}
    for r in range(r0-2,r0+3):
        for c in range(c0-2,c0+3):
            x=x0+sx*c-origin[0];y=y0+sy*(r+1)-origin[1]
            mass=weights.get((r,c))
            ax.add_patch(Rectangle((x,y),sx,-sy,facecolor=cmap(norm(mass)) if mass is not None else '#FFFFFF',edgecolor='#BFC8CE',lw=.65,zorder=0))
            if mass is not None and not zoom:
                ax.text(x+sx*.5,y+(-sy)*.82,f'{mass:.3f}',ha='center',va='center',fontsize=8.1,color='#16354D',bbox=dict(boxstyle='round,pad=.12',facecolor='white',alpha=.85,edgecolor='none'),zorder=5)
    with rasterio.open(base/record['local_map']) as ds:
        arr=ds.read(1);mask=(np.isfinite(arr)&(arr>0)).astype(float)
        pad=np.pad(mask,1);rr,cc=np.meshgrid(np.arange(-1,ds.height+1),np.arange(-1,ds.width+1),indexing='ij')
        lon=ds.transform.c+(cc+.5)*ds.transform.a;lat=ds.transform.f+(rr+.5)*ds.transform.e
        xx,yy=tr.transform(lon,lat);xx-=origin[0];yy-=origin[1]
        ax.contourf(xx,yy,pad,levels=[.5,1.5],colors=['#E5A339'],alpha=.5,zorder=2)
        ax.contour(xx,yy,pad,levels=[.5],colors=['#30333A'],linewidths=.75,zorder=3)
    x=x0+sx*c0-origin[0];y=y0+sy*(r0+1)-origin[1]
    ax.add_patch(Rectangle((x,y),sx,-sy,facecolor='none',edgecolor='#8A4A94',lw=1.25,ls=(0,(4,2)),zorder=4))
    ax.scatter([0],[0],marker='o',s=27 if not zoom else 28,color='#131D2A',edgecolor='white',linewidth=.6,zorder=8)
    px,py=published-origin
    ax.scatter([px],[py],marker='+',s=54 if not zoom else 55,color='#C24B25',linewidth=1.3,zorder=9)
    ax.set_aspect('equal');ax.tick_params(labelsize=7.1,length=2.5,width=.6)
    for sp in ax.spines.values():sp.set_color('#8E9AA3');sp.set_linewidth(.65)
    if not zoom:
        ax.set(xlim=(-510,510),ylim=(-510,510),xticks=[-400,0,400],yticks=[-400,0,400])
        # True geographic north in the local Sinusoidal display.
        north=np.array(tr.transform(g['lon'],g['lat']+.001))-origin;north=north/np.linalg.norm(north)*125
        ax.annotate('',xy=(-405+north[0],330+north[1]),xytext=(-405,330),arrowprops=dict(arrowstyle='-|>',color='#424D57',lw=.7),zorder=10)
        ax.text(-405+north[0],330+north[1]+35,'N',fontsize=7.5,ha='center',color='#424D57')
    return {'registered_native_xy':origin.tolist(),'published_native_xy':published.tolist(),'published_offset_native_xy':(published-origin).tolist(),'cell_width_m':sx,'cell_height_m':-sy,'pixel_weights':pixels}

def spatial(out,data,base):
    fig=plt.figure(figsize=(7.4,5.25));gs=fig.add_gridspec(1,3,left=.075,right=.985,bottom=.40,top=.86,wspace=.17)
    norm=Normalize(0,1);cmap=plt.get_cmap('Blues');meta=[]
    for i,record in enumerate(data['records']):
        ax=fig.add_subplot(gs[0,i]);meta.append(draw_map(ax,record,base,norm,cmap))
        g=record['geometry'];sid=record['site_id'];n=len(record['map']['pixels']);cw=[x['conditional_core_weight'] for x in record['map']['pixels']]
        ax.set_title(f'{chr(97+i)}  {sid} · July {record["year"]}',fontsize=9.5,weight='bold',pad=10)
        if i>0:ax.set_yticklabels([])
        B=' / '.join(f'{1/n:.3f}' for _ in cw);C=' / '.join(f'{w:.3f}' for w in cw)
        ax.text(.5,-.15,f'B: {B}\nC: {C}',transform=ax.transAxes,ha='center',va='top',fontsize=8.1,linespacing=1.5)
        if i==0:
            inset=ax.inset_axes([.06,.07,.42,.42]);draw_map(inset,record,base,norm,cmap,zoom=True)
            inset.set(xlim=(-16,12),ylim=(-15,13),xticks=[],yticks=[]);inset.set_title('Ne1 grid-edge zoom',fontsize=6.8,pad=3)
            inset.text(.5,-.07,'5.57 m geodesic separation',ha='center',va='top',transform=inset.transAxes,fontsize=6.0)
        if i==2:ax.text(.5,.05,'One selected cell\nB = C',transform=ax.transAxes,ha='center',fontsize=8.3,bbox=dict(boxstyle='round,pad=.25',facecolor='white',alpha=.92,edgecolor='#BBC4CD'))
    fig.suptitle('Published monthly core support on the frozen MODIS grid',fontsize=11.1,weight='bold',y=.98)
    fig.text(.5,.928,'Normalized monthly retained 80% cores; actual field boundaries unavailable',ha='center',fontsize=8.0,color='#46545F')
    fig.supxlabel('MODIS Sinusoidal x offset from registered input coordinate (m)',fontsize=8.3,y=.273)
    fig.supylabel('Sinusoidal y offset (m)',fontsize=8.3,x=.01,y=.63)
    cax=fig.add_axes([.30,.196,.40,.020]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal');cb.ax.tick_params(labelsize=7,length=2);cb.set_label('C: conditional retained-core cell weight',fontsize=8,labelpad=2)
    handles=[Rectangle((0,0),1,1,fc='#E5A339',alpha=.5,ec='#30333A',lw=.8,label='Retained monthly core'),Line2D([0],[0],color='#8A4A94',ls='--',lw=1.2,label='A containing cell'),Line2D([0],[0],marker='o',color='none',mfc='#131D2A',mec='white',markersize=5,label='Registered input coordinate'),Line2D([0],[0],marker='+',color='#C24B25',ls='none',markersize=6,label='Published footprint coordinate')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.018),ncol=2,frameon=False,fontsize=7.2,columnspacing=1.5,handlelength=1.7)
    return export(fig,out,'Figure2_spatial_support'),meta

def main():
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();data=json.loads(a.data.read_text(encoding='utf-8'));base=a.data.parent
    a.out.mkdir(parents=True,exist_ok=False)
    for r in data['records']:
        assert sha(base/r['local_map'])==r['map']['sha256']
        assert abs(sum(v['conditional_core_weight'] for v in r['map']['pixels'])-1)<1e-7
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.linewidth':.7,'svg.fonttype':'none','pdf.fonttype':42,'savefig.bbox':None})
    files=design(a.out);generated,geometry=spatial(a.out,data,base);files+=generated
    save(a.out/'PLOT_RECEIPT.json',{'status':'PASS','created_utc':datetime.now(timezone.utc).isoformat(),'code_sha256':sha(__file__),'data_sha256':sha(a.data),'selection_rule':data['selection_rule'],'runtime':{'python':platform.python_version(),'numpy':np.__version__,'matplotlib':matplotlib.__version__,'rasterio':rasterio.__version__},'outputs':files,'geometry':geometry,'native_maps_read':3,'new_fits':0,'new_model_predictions':0,'changed_old_artifacts':0,'scope':'Diagram plus visualization of existing native monthly conditional cores, archived cell mass and source-coordinate ambiguity; no field boundaries or new physical footprint calculation.'})
    print(json.dumps({'status':'PASS','outputs':len(files),'native_maps':3}))
if __name__=='__main__':main()
