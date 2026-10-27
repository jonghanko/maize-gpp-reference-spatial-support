"""Display-only redraw at the exact final insertion width; source CSVs read-only."""
from pathlib import Path
import hashlib,json,platform,itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
import numpy as np
import pandas as pd
import argparse
PACKAGE=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source',type=Path,default=PACKAGE/'data/contrasts12')
parser.add_argument('--out',type=Path,default=PACKAGE/'output/contrasts12')
args=parser.parse_args()
SOURCE=args.source.resolve()
ROOT=args.out.resolve()
ROOT.mkdir(parents=True,exist_ok=True)
W=7.16
FIELDS=['training_given_test_B','training_given_test_C','test_given_training_B','test_given_training_C','interaction','original_pipeline_contrast']
LABELS=['Train C−B | B','Train C−B | C','Test C−B | B','Test C−B | C','Interaction','Pipeline']
REFS=['VUT_REF','CUT_REF']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(name,obj):(ROOT/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.labelsize':12,'axes.titlesize':12,'legend.fontsize':12,'xtick.labelsize':12,'ytick.labelsize':12,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'savefig.bbox':None})
receipts=[]
def bbox_audit(fig,stem):
 fig.canvas.draw();renderer=fig.canvas.get_renderer();width,height=fig.canvas.get_width_height()
 texts=[]
 for obj in fig.findobj(match=Text):
  if not obj.get_visible() or not obj.get_text().strip():continue
  b=obj.get_window_extent(renderer)
  if not b.width or not b.height:continue
  texts.append({'text':obj.get_text(),'font_pt':obj.get_fontsize(),'bbox_px':list(b.extents)})
 outside=[x for x in texts if x['bbox_px'][0]<-0.1 or x['bbox_px'][1]<-0.1 or x['bbox_px'][2]>width+0.1 or x['bbox_px'][3]>height+0.1]
 collisions=[]
 for a,b in itertools.combinations(texts,2):
  x0,y0,x1,y1=a['bbox_px'];u0,v0,u1,v1=b['bbox_px'];dx=min(x1,u1)-max(x0,u0);dy=min(y1,v1)-max(y0,v0)
  if dx>1 and dy>1:collisions.append({'first':a['text'],'second':b['text'],'intersection_px':[dx,dy]})
 axes=[list(ax.get_window_extent(renderer).extents) for ax in fig.axes]
 axis_collisions=[]
 for i,j in itertools.combinations(range(len(axes)),2):
  a,b=axes[i],axes[j]
  if min(a[2],b[2])-max(a[0],b[0])>1 and min(a[3],b[3])-max(a[1],b[1])>1:axis_collisions.append([i,j])
 return {'stem':stem,'font_minimum_pt':min(x['font_pt'] for x in texts),'text_objects':texts,'outside_figure':outside,'text_bbox_collisions':collisions,'axes_bbox_collisions':axis_collisions}
def export(fig,stem,values):
 audit=bbox_audit(fig,stem)
 for ext in ['pdf','svg','png']:fig.savefig(ROOT/(stem+'.'+ext),dpi=600,facecolor='white')
 fig.savefig(ROOT/(stem+'_final_width_100dpi.png'),dpi=100,facecolor='white')
 pd.DataFrame(values).to_csv(ROOT/(stem+'_plot_values.csv'),index=False,float_format='%.17g')
 receipts.append({'stem':stem,'width_in':W,'height_in':float(fig.get_figheight()),'bbox_audit':audit,'outputs':[{'file':stem+'.'+e,'sha256':sha(ROOT/(stem+'.'+e))} for e in ['pdf','svg','png']]})
 plt.close(fig)
def conditional():
 src=SOURCE/'MAIN_CONDITIONALS_OWN_PATHS.csv';before=sha(src);frame=pd.read_csv(src,float_precision='round_trip')
 for stem,panels,h in [
  ('Figure3_BC_conditional_contrasts',[('principal','ridge'),('principal','random_forest'),('T1','ridge'),('T1','random_forest')],6.7),
  ('FigureS7_HGB_conditional_contrasts',[('principal','hist_gradient_boosting'),('T1','hist_gradient_boosting')],4.2)]:
  rows=len(panels)//2;fig,axs=plt.subplots(rows,2,figsize=(W,h),sharex=True,sharey=True);values=[]
  for pi,(ax,(cohort,model)) in enumerate(zip(np.ravel(axs),panels)):
   for ri,ref in enumerate(REFS):
    r=frame.loc[frame.cohort.eq(cohort)&frame.model.eq(model)&frame.fit_reference.eq(ref)].iloc[0]
    x=[float(r[c]) for c in FIELDS];y=np.arange(6)+(ri-.5)*.19
    ax.scatter(x,y,c=['#0072B2','#D55E00'][ri],marker=['o','^'][ri],s=58,edgecolors='white',linewidths=.45,zorder=4)
    for key,value in zip(FIELDS,x):values.append({'panel':chr(97+pi),'cohort':cohort,'model':model,'fit_reference':ref,'score_reference':ref,'contrast':key,'value':value,'n_days':int(r.n_days),'n_units':int(r.n_units)})
   ax.axvline(0,color='#222222',lw=1,zorder=2)
   for boundary in [1.5,3.5,4.5]:ax.axhline(boundary,color='#d6d6d6',lw=.7,zorder=1)
   ax.axhspan(3.5,4.5,color='#f2f2f2',zorder=0)
   ax.set_xlim(-.15,.10);ax.set_ylim(5.5,-.5);ax.set_xticks([-.15,-.10,-.05,0,.05,.10],['−.15','−.10','−.05','0','.05','.10']);ax.grid(axis='x',color='#dedede',lw=.65,zorder=0)
   ax.set_yticks(np.arange(6));ax.set_yticklabels(LABELS)
   ax.tick_params(axis='x',rotation=0,labelbottom=True)
   if pi%2:ax.tick_params(axis='y',labelleft=False)
   ax.text(.015,.97,'('+chr(97+pi)+')',transform=ax.transAxes,ha='left',va='top',fontsize=12,bbox={'facecolor':'white','edgecolor':'none','pad':.3},zorder=8)
  handles=[Line2D([0],[0],marker=['o','^'][i],linestyle='none',color=['#0072B2','#D55E00'][i],markersize=8,label=ref[:3]+' fit / '+ref[:3]+' score') for i,ref in enumerate(REFS)]
  fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.58,.995),ncol=2,frameon=False,handletextpad=.4,columnspacing=1.1)
  fig.supxlabel('Signed RMSE contrast (g C m⁻² d⁻¹)',fontsize=12,y=.02)
  fig.subplots_adjust(left=.23,right=.96,top=.90 if rows==2 else .83,bottom=.12 if rows==2 else .21,hspace=.35,wspace=.22)
  export(fig,stem,values)
 assert sha(src)==before
 return {'file':src.name,'sha256':before}
def figure4():
 src=SOURCE/'TableS6_aggregation_estimands.csv';before=sha(src);agg=pd.read_csv(src)
 primary=['ridge','random_forest'];cts=['C_minus_B','B_minus_A_shared','C_minus_A_shared'];cols=['equal_site_season_absolute_difference','equal_region_pooled_absolute_difference']
 mx=max(abs(agg.loc[agg.model.isin(primary)&agg.contrast.isin(cts),cols[0]]).max(),abs(agg.loc[agg.model.isin(primary)&agg.contrast.isin(cts),cols[1]]).max())*1.45
 fig,axs=plt.subplots(2,2,figsize=(W,7.0),sharey=True);values=[];groups=[]
 for row,column in enumerate(cols):
  for col,model in enumerate(primary):
   pi=row*2+col;ax=axs[row,col]
   for ri,var in enumerate(REFS):
    z=agg.loc[agg.model.eq(model)&agg.variant.eq(var)].set_index('contrast').reindex(cts)
    bars=ax.bar(np.arange(3)+(ri-.5)*.34,z[column],.34,color=['#2166a1','#c87820'][ri],label=var[:3])
    groups.append((ax,ri,list(bars),z[column].to_numpy(dtype=float)))
    for ct,v,b in zip(cts,z[column],bars):values.append({'panel':chr(97+pi),'row':row,'model':model,'variant':var,'contrast':ct,'metric_column':column,'value':float(v),'bar_height':float(b.get_height()),'display':f'{v:+.4f}'})
   ax.set_xticks(range(3),['C − B','B − A\nshared','C − A\nshared']);ax.set_ylim(-mx,mx);ax.set_yticks([-.2,-.1,0,.1,.2]);ax.axhline(0,color='.3',lw=.8)
   ax.text(.02,.97,'('+chr(97+pi)+')',transform=ax.transAxes,ha='left',va='top',fontsize=12)
   ax.grid(axis='y',color='#e4e4e4',lw=.5,zorder=0);ax.set_axisbelow(True)
 fig.legend(handles=[Line2D([0],[0],color=c,lw=8,label=r[:3]) for c,r in zip(['#2166a1','#c87820'],REFS)],loc='upper center',bbox_to_anchor=(.55,.995),ncol=2,frameon=False)
 fig.text(.55,.932,'RMSE difference (g C m⁻² d⁻¹)',ha='center',va='center',fontsize=12)
 fig.subplots_adjust(left=.12,right=.98,top=.88,bottom=.18,wspace=.14,hspace=.75)
 fig.canvas.draw()
 # Exact four-decimal values form two reference-colored readout rows below each
 # axis rather than colliding above adjacent paired bars. The bar geometry stays exact.
 for ax,ri,bars,vals in groups:
  for ci,val in enumerate(vals):
   ax.text(ci,-.34-ri*.13,f'{val:+.4f}',transform=ax.get_xaxis_transform(),ha='center',va='center',fontsize=12,color=['#2166a1','#c87820'][ri],clip_on=False)
 export(fig,'Figure4_primary_LOYO_input_effects',values)
 assert sha(src)==before
 return {'file':src.name,'sha256':before,'ylim':[-float(mx),float(mx)]}
if __name__=='__main__':
 sources=[conditional(),figure4()]
 dump('PLOT_RECEIPT.json',{'status':'GENERATED_REQUIRES_INDEPENDENT_AND_VISUAL_QA','sources':sources,'plots':receipts,'all_target_font_pt':12,'physical_width_in':W,'new_fits':0,'new_predictions':0,'source_csvs_unchanged':True,'runtime':{'python':platform.python_version(),'matplotlib':matplotlib.__version__,'pandas':pd.__version__,'numpy':np.__version__}})
 print(json.dumps({'figures':3,'rows':[48,24,24],'font_minimums':[x['bbox_audit']['font_minimum_pt'] for x in receipts],'bbox_collisions':[len(x['bbox_audit']['text_bbox_collisions']) for x in receipts]}))
