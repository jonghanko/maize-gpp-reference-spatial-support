"""12-point exact-width rendering of frozen design roles and native geometry.
No source data, fitted model, prediction, spatial coordinate or weight is modified.
"""
from pathlib import Path
import argparse,hashlib,json,platform,importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,FancyBboxPatch,FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
P=argparse.ArgumentParser();P.add_argument('--source',type=Path,required=True);P.add_argument('--out',type=Path,required=True);a=P.parse_args()
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('frozen_geometry',ROOT/'original_plot_design_geometry.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
data=json.loads((a.source/'FIGURE_DATA.json').read_text(encoding='utf-8'))
for r in data['records']:assert sha(a.source/r['local_map'])==r['map']['sha256']
a.out.mkdir(parents=True,exist_ok=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.titlesize':12,'axes.labelsize':12,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':12,'svg.fonttype':'none','pdf.fonttype':42,'savefig.bbox':None})
outputs=[];design_roles=[]
def export(fig,stem):
    assert abs(fig.get_figwidth()-7.16)<1e-12
    fig.canvas.draw()
    texts=[t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    assert min(t.get_fontsize() for t in texts)>=12
    for ext in ['pdf','svg','png']:
        p=a.out/(stem+'.'+ext);fig.savefig(p,dpi=600,facecolor='white',metadata={'Creator':'Matplotlib; frozen scientific sources; 12 pt exact 7.16 in style'})
        outputs.append({'file':p.name,'sha256':sha(p),'bytes':p.stat().st_size})
    plt.close(fig)
def box(ax,x,y,w,h,text,face='#F3F7FA'):
    patch=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.006,rounding_size=0.015',facecolor=face,edgecolor='#718899',lw=.8)
    ax.add_patch(patch);ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=12,linespacing=1.25)
    design_roles.append(text);return patch
def arrow(ax,x1,y1,x2,y2):ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=12,lw=.9,color='#4D6576'))
fig=plt.figure(figsize=(7.16,6.0));gs=fig.add_gridspec(2,1,left=.04,right=.985,bottom=.04,top=.98,hspace=.08,height_ratios=[1,1.12])
x=fig.add_subplot(gs[0]);x.set(xlim=(0,1),ylim=(0,1));x.axis('off');x.text(0,.99,'a',ha='left',va='top',weight='bold',fontsize=12)
box(x,.03,.70,.44,.20,'Principal LOYO');box(x,.53,.70,.44,.20,'Mead hindcast')
x.plot([.25,.25,.75,.75],[.69,.62,.62,.69],color='#4D6576',lw=.9);x.plot([.25,.125],[.62,.62],color='#4D6576',lw=.9);arrow(x,.125,.62,.125,.53)
box(x,.02,.15,.215,.37,'Optical inputs\nA · A_shared\nB · C')
box(x,.275,.15,.215,.37,'AgERA5\nmeteorology')
box(x,.53,.15,.19,.37,'Fitting\nreference\nVUT / CUT')
box(x,.76,.15,.22,.37,'Scoring\nreference\nVUT / CUT')
for s,e in [(.237,.274),(.491,.529),(.722,.758)]:arrow(x,s,.335,e,.335)
y=fig.add_subplot(gs[1]);y.set(xlim=(0,1),ylim=(0,1));y.axis('off');y.text(0,.99,'b',ha='left',va='top',weight='bold',fontsize=12)
y.text(.49,.92,'Test input',ha='center',va='center',fontsize=12,weight='bold')
y.text(.35,.82,'B: uniform',ha='center',fontsize=12);y.text(.60,.82,'C: core weighted',ha='center',fontsize=12)
y.text(.06,.46,'Training\ninput',ha='center',va='center',fontsize=12,weight='bold')
y.text(.17,.60,'B',ha='center',va='center',fontsize=12,weight='bold');y.text(.17,.30,'C',ha='center',va='center',fontsize=12,weight='bold')
for px,py,text,diagonal in [(.23,.48,'f_B(X_B)\nR_BB',True),(.48,.48,'f_B(X_C)\nR_BC',False),(.23,.18,'f_C(X_B)\nR_CB',False),(.48,.18,'f_C(X_C)\nR_CC',True)]:box(y,px,py,.23,.24,text,'#DDEBF2' if diagonal else '#FFF0D8')
box(y,.77,.18,.21,.54,'Saved B/C\nmodels')
export(fig,'Figure1_controlled_design')
# All scientific map drawing is called from the byte-exact frozen source.
fig=plt.figure(figsize=(7.16,4.8));gs=fig.add_gridspec(1,3,left=.095,right=.985,bottom=.43,top=.94,wspace=.19)
norm=Normalize(0,1);cmap=plt.get_cmap('Blues');geometry=[]
for i,r in enumerate(data['records']):
    ax=fig.add_subplot(gs[0,i]);geometry.append(old.draw_map(ax,r,a.source,norm,cmap))
    ax.tick_params(labelsize=12,length=2.5,width=.6)
    for t in ax.texts:
        t.set_fontsize(12)
        if t.get_text()=='N':t.set_position((-430,480));t.set_va('top')
    ax.set_title(chr(97+i),loc='left',fontsize=12,weight='bold',pad=6)
    if i>0:ax.tick_params(labelleft=False)
    if i==0:
        inset=ax.inset_axes([.05,.06,.42,.42]);old.draw_map(inset,r,a.source,norm,cmap,zoom=True);inset.set(xlim=(-16,12),ylim=(-15,13),xticks=[],yticks=[]);inset.tick_params(labelsize=12)
fig.text(.53,.405,'MODIS Sinusoidal x offset (m)',ha='center',fontsize=12)
fig.text(.095,.966,'Sinusoidal y offset (m)',ha='left',va='center',fontsize=12)
cax=fig.add_axes([.25,.30,.55,.032]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal');cb.solids.set_rasterized(False);cb.ax.tick_params(labelsize=12,length=3);cb.set_label('C: conditional core cell weight',fontsize=12,labelpad=5)
handles=[Rectangle((0,0),1,1,fc='#E5A339',alpha=.5,ec='#30333A',lw=.8,label='Retained monthly core'),Line2D([0],[0],color='#8A4A94',ls='--',lw=1.2,label='A containing cell'),Line2D([0],[0],marker='o',color='none',mfc='#131D2A',mec='white',markersize=6,label='Registered coordinate'),Line2D([0],[0],marker='+',color='#C24B25',ls='none',markersize=7,label='Published coordinate')]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.53,.035),ncol=2,frameon=False,fontsize=12,columnspacing=1.2,handlelength=1.5)
export(fig,'Figure2_spatial_support')
receipt={'status':'PASS','width_in':7.16,'heights_in':{'Figure1':6.0,'Figure2':4.8},'minimum_artist_font_pt':12,'mathtext_used':False,'source_JSON_sha256':sha(a.source/'FIGURE_DATA.json'),'frozen_geometry_script_sha256':sha(ROOT/'original_plot_design_geometry.py'),'style_script_sha256':sha(__file__),'geometry':geometry,'source_selection_rule':data['selection_rule'],'design_essential_labels':design_roles,'outputs':outputs,'new_fits':0,'new_predictions':0,'source_data_modified':False,'runtime':{'python':platform.python_version(),'matplotlib':matplotlib.__version__,'numpy':np.__version__}}
(a.out/'STYLE12_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8');print(json.dumps({'status':'PASS','outputs':len(outputs),'min_font':12,'new_fits':0}))
