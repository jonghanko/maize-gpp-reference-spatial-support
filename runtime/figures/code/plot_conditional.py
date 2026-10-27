"""Data-derived conditional-contrast plots; source results are read-only."""
from pathlib import Path
import argparse,hashlib,json,platform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

FIELDS=['training_given_test_B','training_given_test_C','test_given_training_B','test_given_training_C','interaction','original_pipeline_contrast']
LABELS=['Training C−B | test B','Training C−B | test C','Test C−B | training B','Test C−B | training C','Interaction (RMSE scale)','Full pipeline CC−BB']
REFS=['VUT_REF','CUT_REF']
COLORS=['#0072B2','#D55E00']
MARKERS=['o','^']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def main(args):
    source=args.results/'MAIN_CONDITIONALS_OWN_PATHS.csv';before=sha(source)
    frame=pd.read_csv(source,float_precision='round_trip')
    assert len(frame)==20 and not frame.duplicated(['cohort','model','fit_reference']).any()
    assert frame.fit_reference.eq(frame.score_reference).all()
    args.out.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','svg.fonttype':'none','pdf.fonttype':42})
    exports=[];figures=[]
    for stem,models,shape in [('Figure3_BC_conditional_contrasts',['ridge','random_forest'],(2,2)),('FigureS7_HGB_conditional_contrasts',['hist_gradient_boosting'],(1,2))]:
        if len(models)==2:
            panels=[('principal','ridge'),('principal','random_forest'),('T1','ridge'),('T1','random_forest')]
            fig,axes=plt.subplots(2,2,figsize=(8.3,7.2),sharex=True,sharey=True)
        else:
            panels=[('principal','hist_gradient_boosting'),('T1','hist_gradient_boosting')]
            fig,axes=plt.subplots(1,2,figsize=(8.3,4.4),sharex=True,sharey=True)
        values=[]
        for pi,(ax,(cohort,model)) in enumerate(zip(np.ravel(axes),panels)):
            for ri,ref in enumerate(REFS):
                r=frame.loc[frame.cohort.eq(cohort)&frame.model.eq(model)&frame.fit_reference.eq(ref)].iloc[0]
                assert r.n_days==(642 if cohort=='principal' else 437)
                x=[float(r[c]) for c in FIELDS];y=np.arange(6)+(ri-.5)*.19
                ax.scatter(x,y,c=COLORS[ri],marker=MARKERS[ri],s=47,edgecolors='white',linewidths=.45,zorder=4)
                for key,value in zip(FIELDS,x):values.append({'cohort':cohort,'model':model,'fit_reference':ref,'score_reference':ref,'contrast':key,'value':value,'units':'g C m^-2 d^-1','estimand':'equal_site_season','n_days':int(r.n_days),'n_units':int(r.n_units)})
            ax.axvline(0,color='#222222',lw=1,zorder=2)
            for boundary in [1.5,3.5,4.5]:ax.axhline(boundary,color='#d6d6d6',lw=.7,zorder=1)
            ax.axhspan(3.5,4.5,color='#f2f2f2',zorder=0)
            ax.set_xlim(-.15,.10);ax.set_ylim(5.5,-.5);ax.set_xticks([-.15,-.10,-.05,0,.05,.10]);ax.grid(axis='x',color='#dedede',lw=.65,zorder=0)
            ax.set_yticks(np.arange(6));ax.set_yticklabels(LABELS)
            if pi%2:ax.tick_params(axis='y',labelleft=False)
            short={'ridge':'Ridge','random_forest':'RF','hist_gradient_boosting':'HGB'}[model]
            label='Principal LOYO' if cohort=='principal' else 'Mead hindcast'
            role=' (secondary)' if cohort=='T1' and model=='ridge' else ''
            ax.set_title('('+chr(97+pi)+') '+label+' · '+short+role,loc='left',pad=23,fontsize=9.3)
            if cohort=='principal':subtitle='642 days · 4 site-seasons · 2 regions'
            else:subtitle='437 days · 3 site-seasons · 1 existing region'
            ax.text(0,1.015,subtitle,transform=ax.transAxes,ha='left',va='bottom',fontsize=7.2,color='#525252')
        fig.suptitle('Conditional representation contrasts from saved models',x=.5,y=.985,fontsize=11,fontweight='bold')
        handles=[Line2D([0],[0],marker=MARKERS[i],linestyle='none',color=COLORS[i],markersize=7,label=ref.replace('_REF','')+' fit / '+ref.replace('_REF','')+' score') for i,ref in enumerate(REFS)]
        fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.62,.941),ncol=2,frameon=False)
        fig.supxlabel('Signed RMSE contrast (g C m$^{-2}$ d$^{-1}$)',fontsize=9,y=.045 if len(models)==2 else .06)
        fig.text(.62,.012,'Negative conditional/pipeline contrasts indicate lower error under C. Interaction is not an error gain.',ha='center',fontsize=6.8,color='#444444')
        fig.subplots_adjust(left=.25,right=.985,top=.81 if len(models)==2 else .77,bottom=.125 if len(models)==2 else .19,hspace=.43,wspace=.15)
        for ext in ['png','svg','pdf']:
            path=args.out/(stem+'.'+ext);fig.savefig(path,dpi=300,metadata={'Title':stem} if ext=='pdf' else None)
            figures.append({'path':path.name,'bytes':path.stat().st_size,'sha256':sha(path)})
        plt.close(fig)
        data=pd.DataFrame(values);data.to_csv(args.out/(stem+'_plot_values.csv'),index=False,float_format='%.17g');exports.append({'path':stem+'_plot_values.csv','rows':len(data),'sha256':sha(args.out/(stem+'_plot_values.csv'))})
    assert sha(source)==before
    caption='''# Figure captions: saved-model B/C crossprediction

Figure 3. Conditional representation contrasts from the original saved models under the equal-site-season RMSE estimand. Rows separate the principal 642-case LOYO experiment and the 437-date earlier-year Mead hindcast; columns show Ridge and RF. Symbols retain the two own-fit/own-score processed reference paths, VUT and CUT. Training contrasts compare C-trained with B-trained models at a fixed test representation; test-input contrasts compare C with B inputs for a fixed trained model. The RMSE-scale interaction is R_CC−R_CB−R_BC+R_BB, and the full-pipeline contrast is the original R_CC−R_BB. First indices denote training and second indices test representation. Negative conditional or pipeline differences indicate lower error for the C condition; the sign of interaction does not itself express an error gain. All cells use the same original cases and fitted models, with no refitting. Ridge and RF are the original principal comparators; Ridge remains a secondary comparator in the separately frozen hindcast. Symbols are descriptive finite-domain comparisons, without population uncertainty bars. The interaction can include the nonlinearity of RMSE and its unit aggregation; it does not identify a unique additive causal decomposition or a physical footprint mechanism.

Figure S7. The same fixed conditional comparisons for the prespecified secondary HGB model, keeping the two populations and both own-reference paths separate. The axes and contrast definitions match Figure 3. All five fixed families and all fitting/scoring paths remain in the complete numeric tables; the displayed three learned families were not selected by their new crossprediction scores.
'''
    (args.out/'FIGURE_CAPTIONS_EN.md').write_text(caption,encoding='utf-8')
    save(args.out/'PLOT_SOURCE_LEDGER.json',{'status':'PASS','source_csv':str(source),'source_csv_sha256':before,'plot_code_sha256':sha(__file__),'runtime':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__},'plots':figures,'numeric_exports':exports,'source_results_mutated':False,'data_selection':'Frozen primary Ridge/RF displays plus prespecified HGB secondary; both own-reference paths; principal LOYO and separate T1; six predeclared RMSE contrasts','x_scale':'Fixed common -0.15 to +0.10 g C m^-2 d^-1, no clipped points','new_fits':0,'new_predictions':0})
    print(json.dumps({'status':'PASS','figures':len(figures),'numeric_rows':[e['rows'] for e in exports],'source_unchanged':True}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
