from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
MODELS = ['train_mean','seasonal_harmonic','ridge','random_forest','hist_gradient_boosting']
PRIMARY = ['ridge','random_forest']
NAMES = {'train_mean':'Mean','seasonal_harmonic':'Seasonal','ridge':'Ridge','random_forest':'RF','hist_gradient_boosting':'HGB'}
VARS = ['VUT_REF','CUT_REF']
CT = ['C_minus_B','B_minus_A_shared','C_minus_A_shared']
LABELS = ['C − B','B − A_shared','C − A_shared']
COLORS = ['#2166a1','#c87820']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return pd.read_csv(p)
def dump(p,v):p.write_text(json.dumps(v,indent=2)+"\n",encoding="utf-8")
def build_plots():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    FIGS.mkdir(parents=True,exist_ok=True)
    plot_sources={str(p):sha(p) for p in [
        AUDIT/'aggregation_metric_contrasts.csv',
        AUDIT/'LOYO_site_season_contrasts_all_families.csv',
        AUDIT/'own_reference_interactions.csv']}
    agg=read(AUDIT/'aggregation_metric_contrasts.csv')
    seasons=read(AUDIT/'LOYO_site_season_contrasts_all_families.csv')
    interactions=read(AUDIT/'own_reference_interactions.csv')
    plt.rcParams.update({'font.size':10,'axes.titlesize':11,'axes.labelsize':10,'svg.fonttype':'none'})
    # Both estimands displayed. Region-pooled is not a pure reweighting of season RMSEs.
    fig,axs=plt.subplots(2,2,figsize=(12.4,8.6),sharey=True)
    mx=max(abs(agg.loc[agg.model.isin(PRIMARY)&agg.contrast.isin(CT),'equal_site_season_absolute_difference']).max(),abs(agg.loc[agg.model.isin(PRIMARY)&agg.contrast.isin(CT),'equal_region_pooled_absolute_difference']).max())*1.45
    for row,(column,title) in enumerate([
        ('equal_site_season_absolute_difference','Mean of four site-season RMSEs'),
        ('equal_region_pooled_absolute_difference','Mean of two region-pooled RMSEs')]):
        for col,model in enumerate(PRIMARY):
            ax=axs[row,col]
            for k,var in enumerate(VARS):
                z=agg.loc[agg.model.eq(model)&agg.variant.eq(var)].set_index('contrast').reindex(CT)
                bars=ax.bar(np.arange(3)+(k-.5)*.34,z[column],.34,label=var[:3],color=COLORS[k])
                ax.bar_label(bars,labels=[f'{v:+.4f}' for v in z[column]],padding=3,fontsize=8)
            ax.set_xticks(range(3),LABELS);ax.set_ylim(-mx,mx);ax.axhline(0,color='.3',lw=.8)
            ax.set_title(NAMES[model]+' · '+title);ax.set_ylabel('RMSE difference (g C m$^{-2}$ d$^{-1}$)')
            if row==0: ax.legend(frameon=False,ncol=2)
    fig.suptitle('Complete year exclusion: input contrasts under two aggregation estimands',fontsize=13)
    fig.text(.5,.025,'Same 642 cases. Negative values favour the first arm. B − A_shared is a post-review descriptive contrast.\nThe alternative changes within-region pooling and regional weighting; C − A_shared remains conditional on point registration.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.08,1,.95))
    for ext in ['png','svg','pdf']:fig.savefig(FIGS/('Figure4_primary_LOYO_input_effects.'+ext),dpi=210)
    plt.close(fig)
    # All five families retained; both reference fit and scoring paths change together.
    fig,axs=plt.subplots(1,3,figsize=(14.3,4.9),sharey=True)
    lim=max(abs(interactions.CUT_minus_VUT_own_interaction_RMSE_units).max()*1.4,.001)
    for ax,contrast,label in zip(axs,CT,LABELS):
        z=interactions.loc[interactions.contrast.eq(contrast)].set_index('model').reindex(MODELS)
        vals=z.CUT_minus_VUT_own_interaction_RMSE_units
        bars=ax.bar(range(5),vals,color=['#bbc1c7','#98a2aa','#647b8d','#2166a1','#80b4cc'])
        ax.bar_label(bars,labels=[f'{v:+.4f}' for v in vals],padding=3,fontsize=8)
        ax.set_xticks(range(5),[NAMES[x] for x in MODELS]);ax.set_ylim(-lim,lim);ax.axhline(0,color='.3',lw=.8)
        ax.set_title(label)
    axs[0].set_ylabel('CUT contrast − VUT contrast\n(RMSE-difference units: g C m$^{-2}$ d$^{-1}$)')
    fig.suptitle('Dependent reference paths: own-fit / own-score input interactions',fontsize=13)
    fig.text(.5,.025,'All five families, equal-site-season estimand. B − A_shared was added after review.\nThis changes fitting and scoring reference together; it does not isolate uStar or establish a causal interaction.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.09,1,.94))
    for ext in ['png','svg','pdf']:fig.savefig(FIGS/('FigureS5_reference_input_interactions.'+ext),dpi=210)
    plt.close(fig)
    order=[('US-Ne1',2012),('US-Ne1',2013),('US-Ne3',2013),('US-Ro1',2015)]
    labels=['Ne1 2012\nn=156','Ne1 2013\nn=167','Ne3 2013\nn=157','Ro1 2015\nn=162']
    p=seasons.loc[seasons.model.isin(PRIMARY)&seasons.contrast.isin(CT)]
    lim=abs(p.absolute_RMSE_difference).max()*1.38
    fig,axs=plt.subplots(3,2,figsize=(12.6,11.2),sharey=True)
    for row,(contrast,label) in enumerate(zip(CT,LABELS)):
        for col,model in enumerate(PRIMARY):
            ax=axs[row,col];ax.axvspan(2.55,3.45,color='#ededed',zorder=0)
            for k,var in enumerate(VARS):
                z=p.loc[p.model.eq(model)&p.contrast.eq(contrast)&p.variant.eq(var)].set_index(['site_id','harvest_year']).reindex(order)
                bars=ax.bar(np.arange(4)+(k-.5)*.32,z.absolute_RMSE_difference,.32,label=var[:3],color=COLORS[k])
                ax.bar_label(bars,labels=[f'{v:+.3f}' for v in z.absolute_RMSE_difference],padding=2,fontsize=7)
            ax.set_ylim(-lim,lim);ax.axhline(0,color='.3',lw=.8);ax.set_xticks(range(4),labels)
            ax.set_title(NAMES[model]+' · '+label);ax.set_ylabel('RMSE difference (g C m$^{-2}$ d$^{-1}$)')
            if row==0:ax.legend(frameon=False,ncol=2,fontsize=9)
    fig.suptitle('All four LOYO site-seasons: directions vary by model and comparator',fontsize=13)
    fig.text(.5,.025,'Shaded Ro1 has identical test inputs across arms; its exact 2015 operation dates remain unknown.\nThe same Ro1 train/test membership recurs in LOYO-2015, partial LOSO and joint exclusion. No uncertainty intervals are estimated.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.065,1,.96))
    for ext in ['png','svg','pdf']:fig.savefig(FIGS/('FigureS6_all_LOYO_site_season_effects.'+ext),dpi=210)
    plt.close(fig)
    assert all(sha(p)==h for p,h in plot_sources.items())
    assert len(list(FIGS.glob('Figure*.png')))==3
    assert len(list(FIGS.glob('Figure*.svg')))==3
    dump(OUT/'PLOT_BUILD_VERIFICATION.json',{
        'status':'pass','source_csvs_unchanged':len(plot_sources),
        'scientific_model_fits':0,'matplotlib_version':matplotlib.__version__,
        'figure3_layout':'3 contrast rows by 2 primary model-family columns; 6 panels',
        'source_csv_hashes':plot_sources,
        'output_hashes':{str(p.relative_to(OUT)):sha(p) for p in FIGS.glob('Figure*.*')}})
