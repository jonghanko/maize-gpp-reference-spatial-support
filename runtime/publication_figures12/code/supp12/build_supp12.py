"""Layout-only supplement plots from immutable public aggregate inputs."""
from pathlib import Path
from fractions import Fraction
from decimal import Decimal, localcontext
import csv, io, json, hashlib, math, itertools, sys, os, argparse
SCRIPT_ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description='Layout-only S1-S6 figures from frozen aggregate inputs; no model fitting or prediction.')
parser.add_argument('--source',type=Path,required=True,help='Immutable frozen figure module, with data/ and FIGURE_INPUT_MANIFEST.json')
parser.add_argument('--out',type=Path,required=True,help='Writable output directory outside the immutable input module')
args=parser.parse_args()
ROOT=args.out.resolve();SOURCE=args.source.resolve()
assert not ROOT.is_relative_to(SOURCE),'Keep outputs outside the immutable source module'
ROOT.mkdir(parents=True,exist_ok=True)
os.environ['MPLCONFIGDIR']=str(ROOT/'_mplconfig')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.ticker import MaxNLocator, FormatStrFormatter

DATA=SOURCE/'data';TABLES=DATA/'tables'
OUT=ROOT/'figures';OUT.mkdir(parents=True,exist_ok=True)
QA=ROOT/'visual_qa';QA.mkdir(exist_ok=True)
W=7.16
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.labelsize':12,'axes.titlesize':12,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':12,'figure.titlesize':12,'figure.labelsize':12,'svg.fonttype':'none','pdf.fonttype':42,'axes.formatter.use_mathtext':False,'savefig.facecolor':'white'})
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
manifest=json.loads((SOURCE/'FIGURE_INPUT_MANIFEST.json').read_text(encoding='utf-8'))
for row in manifest['files']:assert sha(SOURCE/row['path'])==row['sha256'],row['path']
source_before={row['path']:row['sha256'] for row in manifest['files']}
receipts=[];plotvalues={};caption_moves=[]

def caption(name,panels,moved):caption_moves.append({'figure':name,'panel_mapping':panels,'moved_to_caption':moved})
def panel(ax,label):ax.text(-.13,1.025,label,transform=ax.transAxes,ha='left',va='bottom',fontsize=12)
def save(fig,name,values,identity):
 fig.canvas.draw()
 renderer=fig.canvas.get_renderer();texts=[]
 for t in fig.findobj(match=Text):
  if t.get_visible() and t.get_text().strip():
   assert t.get_fontsize()>=12,(name,t.get_text(),t.get_fontsize())
   box=t.get_window_extent(renderer=renderer)
   if box.width>0 and box.height>0:texts.append({'text':t.get_text(),'size_pt':t.get_fontsize(),'bbox_px':[float(box.x0),float(box.y0),float(box.x1),float(box.y1)]})
 collisions=[]
 for a,b in itertools.combinations(texts,2):
  x=max(0,min(a['bbox_px'][2],b['bbox_px'][2])-max(a['bbox_px'][0],b['bbox_px'][0])); y=max(0,min(a['bbox_px'][3],b['bbox_px'][3])-max(a['bbox_px'][1],b['bbox_px'][1]))
  if x>1 and y>1:collisions.append({'first':a['text'],'second':b['text'],'intersection_px':[x,y]})
 extent=fig.get_window_extent(renderer);outside=[]
 for t in texts:
  x0,y0,x1,y1=t['bbox_px']
  if x0<-.5 or y0<-.5 or x1>extent.width+.5 or y1>extent.height+.5:outside.append(t)
 for ext in ['pdf','svg','png']:fig.savefig(OUT/(name+'.'+ext),dpi=600,bbox_inches=None,pad_inches=0)
 fig.savefig(QA/(name+'.png'),dpi=150,bbox_inches=None,pad_inches=0)
 receipt={'figure':name,'canvas_width_in':W,'canvas_height_in':float(fig.get_size_inches()[1]),'minimum_artist_font_pt':min(t['size_pt'] for t in texts),'all_artist_fonts_at_least_12':True,'text_artist_count':len(texts),'text_bounding_boxes':texts,'text_pair_intersections':collisions,'outside_canvas_text':outside,'axis_bounds':[list(a.get_position().bounds) for a in fig.axes],'value_identity':identity,'new_fits':0,'new_predictions':0,'files':{p.name:{'bytes':p.stat().st_size,'sha256':sha(p)} for p in OUT.glob(name+'.*')}}
 receipts.append(receipt);plotvalues[name]=values
 plt.close(fig)

# S1: same exact-rational cells and zero-centered color range.
data=json.loads((DATA/'OMISSION_EXACT.json').read_text(encoding='utf-8'))['all_summaries']
ix={(r['model'],r['reference'],r['scenario']):r for r in data if r['contrast']=='C_minus_B'}
order=[('ridge','VUT_REF'),('ridge','CUT_REF'),('random_forest','VUT_REF'),('random_forest','CUT_REF')]
scenarios=['full_four_site_seasons','omit_US-Ne1-2012','omit_US-Ne1-2013','omit_US-Ne3-2013','omit_US-Ro1-2015']
values=[];exactvalues=[]
for m,r in order:
 row=[]
 for s in scenarios:
  q=ix[m,r,s]['relative_difference_percent'];frac=Fraction(int(q['numerator']),int(q['denominator']));row.append(float(frac));exactvalues.append({'model':m,'reference':r,'scenario':s,'numerator':frac.numerator,'denominator':frac.denominator})
 values.append(row)
check=pd.read_csv(TABLES/'TableS13_fixed_OOF_summary_omissions.csv').set_index(['model','reference','scenario','contrast'])
for (m,r),row in zip(order,values):
 for s,v in zip(scenarios,row):assert abs(v-float(check.loc[(m,r,s,'C_minus_B'),'relative_difference_percent']))<1e-12
extent=math.ceil(np.max(np.abs(values))*10)/10
fig=plt.figure(figsize=(W,3.75));ax=fig.add_axes([.16,.24,.68,.69]);cax=fig.add_axes([.87,.24,.025,.69])
im=ax.imshow(values,cmap=LinearSegmentedColormap.from_list('C_B_signed',['#2166ac','#fafafa','#d6601f']),norm=TwoSlopeNorm(vmin=-extent,vcenter=0,vmax=extent),aspect='auto')
ax.set_xticks(range(5),['Full\n4 seasons','Omit\nNe1 2012','Omit\nNe1 2013','Omit\nNe3 2013','Omit\nRo1 2015']);ax.set_yticks(range(4),['Ridge · VUT','Ridge · CUT','RF · VUT','RF · CUT']);ax.tick_params(length=0,pad=7)
ax.set_xticks(np.arange(-.5,5,.5)[::2],minor=True);ax.set_yticks(np.arange(-.5,4,.5)[::2],minor=True);ax.grid(which='minor',color='white',linewidth=2);ax.tick_params(which='minor',length=0)
for i,row in enumerate(values):
 for j,v in enumerate(row):ax.text(j,i,f'{v:+.3f}%',ha='center',va='center',fontsize=12,color='white' if abs(v)>.65*extent else '#101820')
for spine in ax.spines.values():spine.set_visible(False)
cb=fig.colorbar(im,cax=cax,ticks=[-3,0,3]);cb.set_label('RMSE difference (%)',fontsize=12,labelpad=7)
save(fig,'FigureS1_fixed_OOF_summary_omissions',{'cells':values,'exact_cells':exactvalues,'color_limits':[-extent,extent]},'20 exact-rational JSON cells agree with frozen table within 1e-12; colors and groups unchanged')
caption('S1',{'single_panel':'20-cell C-minus-B omission heatmap; same row and column order'},['C − B sensitivity to site-season omissions from the error summary','Fixed LOYO out-of-fold outputs; blue means lower C error and orange higher C error.','Evaluation-summary omissions only: no retraining, new validation or uncertainty intervals.','Omitting Ro1 leaves only Mead in evaluation; its labels remain in the original training folds.','The full column contains four site-seasons; each omission leaves three.'])

# S2: same top macro scores and bottom per-season percentages.
macro=pd.read_csv(TABLES/'TableS17a_T1_all_reference_macro.csv');seasons=pd.read_csv(TABLES/'TableS18_T1_site_season.csv')
own=macro.loc[macro.fit_reference.eq(macro.score_reference)];ownseason=seasons.loc[seasons.fit_reference.eq(seasons.score_reference)]
vals=[];bars_values={}
for site,year in [('US-Ne1',2010),('US-Ne1',2011),('US-Ne3',2011)]:
 row=[]
 for model,ref in [('random_forest','VUT_REF'),('random_forest','CUT_REF'),('ridge','VUT_REF'),('ridge','CUT_REF')]:
  z=ownseason.loc[ownseason.site_id.eq(site)&ownseason.harvest_year.eq(year)&ownseason.model.eq(model)&ownseason.fit_reference.eq(ref)]
  b=float(z.loc[z.arm.eq('B_uniform_pixels'),'rmse'].iloc[0]);c=float(z.loc[z.arm.eq('C_core_weighted'),'rmse'].iloc[0]);row.append(100*(c-b)/b)
 vals.append(row)
vals=np.array(vals);limit=max(float(abs(vals).max()),.05)
fig=plt.figure(figsize=(W,6.65));axs=[fig.add_axes([.12,.57,.36,.34]),fig.add_axes([.61,.57,.36,.34])]
handles=[]
for col,(model,label) in enumerate([('random_forest','RF'),('ridge','Ridge')]):
 ax=axs[col];top=[]
 for k,(off,arm,color,title) in enumerate([(-.18,'B_uniform_pixels','#687d8e','Uniform (B)'),(.18,'C_core_weighted','#b55b35','Core weighted (C)')]):
  y=[float(own.loc[own.model.eq(model)&own.fit_reference.eq(ref)&own.arm.eq(arm),'macro_site_season_rmse'].iloc[0]) for ref in ['VUT_REF','CUT_REF']]
  top+=y;bars=ax.bar(np.arange(2)+off,y,width=.34,color=color,label=title)
  for bar,v in zip(bars,y):ax.annotate(f'{v:.3f}',(bar.get_x()+bar.get_width()/2,v),xytext=(0,3+14*k),textcoords='offset points',ha='center',va='bottom',fontsize=12)
  if col==0:handles.append(Patch(facecolor=color,label=title))
  bars_values[model+'/'+arm]=y
 ax.set_xticks([0,1],['VUT','CUT']);ax.set_ylim(0,max(top)*1.05*1.4);ax.set_yticks(np.arange(0,int(ax.get_ylim()[1])+1));ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True);panel(ax,'('+chr(97+col)+')')
axs[0].set_ylabel('RMSE (g C m⁻² d⁻¹)',fontsize=12)
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.53,1.005),frameon=False,ncol=2,fontsize=12)
ax=fig.add_axes([.19,.13,.64,.30]);cax=fig.add_axes([.86,.13,.025,.30]);im=ax.imshow(vals,cmap='RdBu_r',vmin=-limit,vmax=limit,aspect='auto')
ax.set_xticks(range(4),['RF · VUT','RF · CUT','Ridge · VUT','Ridge · CUT'],rotation=20,ha='right');ax.set_yticks(range(3),['Ne1 2010','Ne1 2011','Ne3 2011']);panel(ax,'(c)')
for i in range(3):
 for j in range(4):ax.text(j,i,f'{vals[i,j]:+.2f}%',ha='center',va='center',fontsize=12,color='white' if abs(vals[i,j])>limit*.62 else '#141414')
cb=fig.colorbar(im,cax=cax,ticks=[-3,0,3]);cb.set_label('Difference (%)',fontsize=12,labelpad=5)
save(fig,'FigureS2_T1_hindcast',{'macro_bars':bars_values,'relative_cells':vals.tolist(),'color_limits':[-limit,limit]},'Same own-fit/own-score macro and site-season rows as original frozen plot; original per-family y limits retained')
caption('S2',{'(a)':'Random forest, primary','(b)':'Ridge, secondary','(c)':'Per-site-season C-minus-B percentages, same 3×4 order'},['Earlier-year Mead hindcast: train 2012/2013 and test 2010/2011.','480 training days; 437 test days; one existing region.','Panel (c) season counts: Ne1 2010 n=128, Ne1 2011 n=148, Ne3 2011 n=161.','VUT/CUT represent paired fit/score reference paths.','Negative percentages indicate lower RMSE; dependent processed GPP references, fixed monthly 80% DAY cores; retrospective, not a regional or forecast test.'])

# S3: same 18 rows and algebraic components, stacked for 12-point readability.
df=pd.read_csv(TABLES/'TableS22_component_summary_estimands.csv');d=df.loc[df.scope.eq('equal_site_season')&df.fit_reference.eq(df.score_reference)]
panels=[('principal_642','C_core_weighted__minus__B_uniform_pixels'),('principal_642','B_uniform_pixels__minus__A_shared_scene'),('T1_437','C_core_weighted__minus__B_uniform_pixels')]
order=[(m,r) for m in ['ridge','random_forest','hist_gradient_boosting'] for r in ['VUT_REF','CUT_REF']]
frames=[d.loc[d.cohort.eq(cohort)&d.contrast.eq(contrast)].set_index(['model','fit_reference']).loc[order] for cohort,contrast in panels]
allvals=np.concatenate([x[['mean_bias_rmse_component','mean_centered_rmse_component','mean_delta_rmse']].to_numpy().ravel() for x in frames]);lo=float(min(allvals)*1.15);hi=float(max(.02,max(allvals)*1.6))
fig=plt.figure(figsize=(W,7.5));labels=[{'ridge':'Ridge','random_forest':'RF','hist_gradient_boosting':'HGB'}[m]+' · '+ref[:3] for m,ref in order]
for i,(bottom,q) in enumerate(zip([.745,.49,.235],frames)):
 ax=fig.add_axes([.20,bottom,.77,.21]);y=np.arange(6)
 ax.barh(y-.14,q.mean_bias_rmse_component,height=.23,color='#D55E00',zorder=3);ax.barh(y+.14,q.mean_centered_rmse_component,height=.23,color='#0072B2',zorder=3);ax.scatter(q.mean_delta_rmse,y,s=35,marker='D',color='#202020',zorder=4)
 ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlim(lo,hi);ax.axvline(0,color='#777777',lw=.8);ax.grid(axis='x',alpha=.2,zorder=0);ax.spines[['top','right']].set_visible(False);ax.set_xticks([v for v in [-.15,-.10,-.05,0,.05] if lo<=v<=hi]);ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'));panel(ax,'('+chr(97+i)+')')
 if i<2:ax.tick_params(axis='x',labelbottom=False)
 else:ax.set_xlabel('Contribution to mean RMSE difference (g C m⁻² d⁻¹)',fontsize=12,labelpad=8)
leg=[Patch(facecolor='#D55E00',label='Squared-mean-error component'),Patch(facecolor='#0072B2',label='Centered-error component'),Line2D([0],[0],marker='D',color='#202020',linestyle='none',label='Total RMSE difference')]
fig.legend(handles=leg,loc='lower center',bbox_to_anchor=(.55,.005),frameon=False,ncol=1,fontsize=12)
save(fig,'FigureS3_fixed_error_components',{'panels':[q[['mean_bias_rmse_component','mean_centered_rmse_component','mean_delta_rmse']].to_dict('records') for q in frames],'x_limits':[lo,hi]},'Same frozen component matrix, 18 displayed model/reference rows, unchanged shared numeric axis')
caption('S3',{'(a)':'Principal cohort C minus B','(b)':'Principal cohort B minus shared-scene point','(c)':'Earlier-year Mead hindcast C minus B'},['The three original panels are stacked vertically; scientific order and shared limits are unchanged.','Orange bars are squared-mean-error contributions, blue bars centered-error contributions, and diamonds their sum.','Both own-reference paths are retained for Ridge, random forest and secondary histogram gradient boosting.'])

# S4: same 45/84 monthly support cells and 7 realized-input rows.
m=pd.read_csv(TABLES/'TableS23_monthly_weight_concentration.csv');s=pd.read_csv(TABLES/'TableS25_site_season_input_summary.csv')
order=[('US-Ne1',2012),('US-Ne1',2013),('US-Ne3',2013),('US-Ro1',2015),('US-Ne1',2010),('US-Ne1',2011),('US-Ne3',2011)]
labels=[('P' if i<4 else 'T1')+' · '+site.replace('US-','')+' '+str(year) for i,(site,year) in enumerate(order)]
matrix=np.full((7,12),np.nan)
for i,(site,year) in enumerate(order):
 for _,row in m.loc[m.site_id.eq(site)&m.year.eq(year)&m.used_by_retained_case_anchor.eq(True)].iterrows():matrix[i,int(row.month)-1]=row.concentration_cell_count
assert np.isfinite(matrix).sum()==45
fig=plt.figure(figsize=(W,7.5));a=fig.add_axes([.20,.58,.69,.36]);cax=fig.add_axes([.91,.58,.018,.36]);cmap=plt.get_cmap('cividis').copy();cmap.set_bad('#eeeeee');im=a.imshow(np.ma.masked_invalid(matrix),vmin=1,vmax=3,cmap=cmap,aspect='auto')
for i in range(7):
 for j in range(12):
  if np.isfinite(matrix[i,j]):a.text(j,i,f'{matrix[i,j]:.2f}',ha='center',va='center',fontsize=12,color='white' if matrix[i,j]<1.7 else '#111111')
a.set_yticks(range(7),labels);a.set_xticks(range(12),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],rotation=45,ha='right');a.set_xticks(np.arange(-.5,12,.5)[::2],minor=True);a.set_yticks(np.arange(-.5,7,.5)[::2],minor=True);a.grid(which='minor',color='white',linewidth=.6);a.tick_params(which='minor',length=0);a.axhline(3.5,color='white',lw=2);a.set_xlabel('Month',fontsize=12,labelpad=5);panel(a,'(a)');cb=fig.colorbar(im,cax=cax,ticks=[1,2,3]);cb.set_label('Nc',fontsize=12,labelpad=5)
q=s.set_index(['site_id','year']).loc[order];b=fig.add_axes([.22,.18,.71,.29]);y=np.arange(7);b.barh(y,q.mean_absolute_NDVI_difference,color=['#0072B2']*4+['#D55E00']*3,height=.45);b.scatter(q.max_absolute_NDVI_difference,y,marker='o',s=45,facecolors='none',edgecolors='#222222',zorder=4)
b.set_yticks(y,labels);b.invert_yaxis();b.axhline(3.5,color='#aaaaaa',lw=.8);b.grid(axis='x',alpha=.2);b.set_axisbelow(True);b.spines[['top','right']].set_visible(False);b.set_xlim(-.002,float(max(q.max_absolute_NDVI_difference)*1.13));b.set_xticks([v for v in [0,.02,.04,.06,.08,.10,.12] if b.get_xlim()[0]<=v<=b.get_xlim()[1]]);b.set_xlabel('Absolute NDVI difference',fontsize=12,labelpad=7);panel(b,'(b)')
fig.legend(handles=[Line2D([0],[0],color='#555555',lw=7,label='Mean absolute difference'),Line2D([0],[0],marker='o',color='#222222',markerfacecolor='none',linestyle='none',label='Maximum absolute difference')],loc='lower center',bbox_to_anchor=(.58,.008),frameon=False,ncol=1,fontsize=12)
save(fig,'FigureS4_conditional_input_support',{'monthly_grid':[[None if not np.isfinite(v) else float(v) for v in row] for row in matrix],'mean_absolute_difference':q.mean_absolute_NDVI_difference.tolist(),'maximum_absolute_difference':q.max_absolute_NDVI_difference.tolist(),'color_limits':[1,3],'x_limits':list(b.get_xlim())},'Same monthly support and realized-input source rows; 45 finite map cells, fixed color limits and input-difference x limits unchanged')
caption('S4',{'(a)':'Weight concentration Nc in monthly maps used by retained-date anchors','(b)':'Realized input differences between C and B on retained dates'},['Panels are stacked vertically for readability.','Month is the actual composite-start month of a retained-date anchor.','P denotes principal 642 dates; T1 denotes separate 437-date hindcast.','Blank light-grey months are unused by retained-date anchors.','The quantity Nc is conditional concentration cell count, not an independent sample size or field-purity metric.'])

# S5: same 15 own-reference processing-path contrasts and shared y limits.
interactions=pd.read_csv(TABLES/'TableS3_reference_input_interactions.csv');models=['train_mean','seasonal_harmonic','ridge','random_forest','hist_gradient_boosting'];names={'train_mean':'Mean','seasonal_harmonic':'Seasonal','ridge':'Ridge','random_forest':'RF','hist_gradient_boosting':'HGB'}
cts=['C_minus_B','B_minus_A_shared','C_minus_A_shared'];lim=float(max(abs(interactions.CUT_minus_VUT_own_interaction_RMSE_units).max()*1.4,.001));fig=plt.figure(figsize=(W,7.5));vlist=[]
for i,(bottom,contrast) in enumerate(zip([.735,.425,.115],cts)):
 ax=fig.add_axes([.17,bottom,.79,.215]);z=interactions.loc[interactions.contrast.eq(contrast)].set_index('model').reindex(models);vals=z.CUT_minus_VUT_own_interaction_RMSE_units.to_numpy();vlist.append(vals.tolist());bars=ax.bar(range(5),vals,color=['#bbc1c7','#98a2aa','#647b8d','#2166a1','#80b4cc'])
 for bar,v in zip(bars,vals):ax.annotate(f'{v:+.4f}',(bar.get_x()+bar.get_width()/2,v),xytext=(0,3 if v>=0 else -3),textcoords='offset points',ha='center',va='bottom' if v>=0 else 'top',fontsize=12)
 ax.set_xticks(range(5),[names[x] for x in models]);ax.set_ylim(-lim,lim);ax.set_yticks([-.04,-.02,0,.02,.04]);ax.axhline(0,color='.3',lw=.8);panel(ax,'('+chr(97+i)+')')
fig.text(.014,.535,'Reference-path contrast difference (g C m⁻² d⁻¹)',rotation=90,ha='center',va='center',fontsize=12)
save(fig,'FigureS5_reference_input_interactions',{'bar_values':vlist,'shared_y_limits':[-lim,lim],'contrasts':cts,'models':models},'Same 15 source bar coordinates, model order, own-reference path comparison and shared y limits')
caption('S5',{'(a)':'C minus B input contrast','(b)':'B minus shared-scene point input contrast','(c)':'C minus shared-scene point input contrast'},['Dependent reference paths: own-fit/own-score input interactions.','All five families are retained under the equal-site-season estimand.','The B-minus-shared-point contrast was added after review.','Both fitting and scoring reference change together; this does not isolate uStar or establish a causal interaction.'])

# S6: same 48 bars, fixed common y range, 3×2 family/contrast matrix.
seasons=pd.read_csv(TABLES/'TableS7_site_season_contrasts.csv');primary=['ridge','random_forest'];p=seasons.loc[seasons.model.isin(primary)&seasons.contrast.isin(cts)];lim=float(abs(p.absolute_RMSE_difference).max()*1.38)
order=[('US-Ne1',2012),('US-Ne1',2013),('US-Ne3',2013),('US-Ro1',2015)];labels=['Ne1\n2012','Ne1\n2013','Ne3\n2013','Ro1\n2015'];fig=plt.figure(figsize=(W,7.5));bar_data=[]
for row,(bottom,contrast) in enumerate(zip([.79,.48,.17],cts)):
 for col,model in enumerate(primary):
  ax=fig.add_axes([.12+.47*col,bottom,.38,.16]);ax.axvspan(2.55,3.45,color='#ededed',zorder=0)
  for k,ref in enumerate(['VUT_REF','CUT_REF']):
   z=p.loc[p.model.eq(model)&p.contrast.eq(contrast)&p.variant.eq(ref)].set_index(['site_id','harvest_year']).reindex(order);vals=z.absolute_RMSE_difference.to_numpy();bars=ax.bar(np.arange(4)+(k-.5)*.32,vals,.32,color=['#2166a1','#c87820'][k])
   for x,v in enumerate(vals):
    ax.text(x,-.52-.19*k,f'{v:+.3f}',transform=ax.get_xaxis_transform(),ha='center',va='center',fontsize=12,color=['#2166a1','#a65d0a'][k],clip_on=False)
   bar_data.append({'panel':row*2+col,'contrast':contrast,'model':model,'reference':ref,'values':vals.tolist()})
  ax.set_ylim(-lim,lim);ax.set_yticks([-.6,-.3,0,.3,.6]);ax.axhline(0,color='.3',lw=.8);ax.set_xticks(range(4),labels);panel(ax,'('+chr(97+row*2+col)+')')
  if col:ax.tick_params(axis='y',labelleft=False)
fig.text(.015,.535,'RMSE difference (g C m⁻² d⁻¹)',rotation=90,ha='center',va='center',fontsize=12)
fig.legend(handles=[Patch(facecolor='#2166a1',label='VUT'),Patch(facecolor='#c87820',label='CUT')],loc='upper center',bbox_to_anchor=(.55,1.005),frameon=False,ncol=2,fontsize=12)
save(fig,'FigureS6_all_LOYO_site_season_effects',{'bars':bar_data,'shared_y_limits':[-lim,lim],'season_order':order},'Same 48 bar heights in the same six model/contrast panels; original common y limits and Ro1 shading retained')
caption('S6',{'(a)':'Ridge C minus B','(b)':'Random forest C minus B','(c)':'Ridge B minus shared-scene point','(d)':'Random forest B minus shared-scene point','(e)':'Ridge C minus shared-scene point','(f)':'Random forest C minus shared-scene point'},['All four LOYO site-seasons: directions vary by model and comparator.','Numeric readouts beneath each panel show VUT in the upper blue row and CUT in the lower orange row; all original 48 numbers and bar heights are retained.','Site-season sample counts, shared across all panels: Ne1 2012 n=156; Ne1 2013 n=167; Ne3 2013 n=157; Ro1 2015 n=162.','Shaded Ro1 has identical test inputs across arms; its exact 2015 operation dates remain unknown.','The same Ro1 train/test membership recurs in LOYO-2015, partial LOSO and joint exclusion.','No uncertainty intervals are estimated.'])

assert all(sha(SOURCE/path)==v for path,v in source_before.items())
(ROOT/'SOURCE_AND_VALUES.json').write_text(json.dumps({'source_files_verified_before_and_after':source_before,'plot_values':plotvalues,'new_fits':0,'new_predictions':0,'changed_scientific_input_bytes':0},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
(ROOT/'FIGURE_FONT_AND_BBOX_QA.json').write_text(json.dumps({'target_width_in':W,'minimum_required_font_pt':12,'figures':receipts},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
lines=['# Supplementary figure text moved to captions','', 'All descriptive titles and notes listed below have been removed from plotted artwork. Essential model/site/arm identifiers, panel letters, axes, legends and numeric labels remain at least 12 pt at exactly 7.16 inches width. No model fits, predictions, source values, axis limits, color limits or case selections changed.','']
for c in caption_moves:
 lines += ['## Figure '+c['figure'],'','Panel mapping:']
 lines += ['- '+k+': '+v+'.' for k,v in c['panel_mapping'].items()]
 lines += ['','Move or retain these statements in the caption:']
 lines += ['- '+x for x in c['moved_to_caption']]
 lines += ['']
(ROOT/'CAPTION_MOVES_EN.md').write_text('\n'.join(lines),encoding='utf-8')
print(json.dumps({'figures':len(receipts),'font_min':min(x['minimum_artist_font_pt'] for x in receipts),'outside_text_counts':{x['figure']:len(x['outside_canvas_text']) for x in receipts},'bbox_intersections':{x['figure']:len(x['text_pair_intersections']) for x in receipts}},ensure_ascii=False))
