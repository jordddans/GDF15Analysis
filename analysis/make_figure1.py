"""November Figure 1 layout, evaluated from the archived finite bootstrap fits.
Preserves draw(results_dir, figures_dir, workbook=None) for the existing pipeline.
No new seed, fit, or resampling is introduced by plotting.
"""
import json
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import data

def draw(results_dir, figures_dir, workbook=None):
    results_dir, figures_dir = Path(results_dir), Path(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)
    entry=json.loads((results_dir/"results.json").read_text())["models"]["IPAH+HPAH"]
    meta=json.loads((results_dir/"figure1_bootstrap.json").read_text())
    assert (meta["seed"],meta["requested"],meta["finite"],meta["separated"])==(20260907,4000,3991,9)
    coefpath=results_dir/"figure1_bootstrap_coefficients.csv"
    co=pd.read_csv(coefpath,float_precision="round_trip").to_numpy()
    assert co.shape==(3991,2) and np.isfinite(co).all()
    frame=data.load_gdf15(workbook)
    groups=[frame.loc[frame.cohort==g,"gdf15"].to_numpy() for g in data.GROUP1]
    pv,ip,hp=groups; comp=np.r_[ip,hp]; x=np.r_[pv,comp]
    y=np.r_[np.ones(len(pv),bool),np.zeros(len(comp),bool)]
    assert (len(pv),len(ip),len(hp))==(26,98,42)
    low,high=pv.min(),comp.max()
    offset=entry["sampling_log_odds"]
    b=np.array([entry["intercept"],entry["slope_log10_gdf15"]])
    lr=lambda v:np.exp(b[0]+b[1]*np.log10(v)-offset)
    post=lambda v:100*v/(9+v)
    cut=1658.
    sens=(pv>=cut).mean();spec=(comp<cut).mean()
    lp,ln=sens/(1-spec),(1-sens)/spec
    own=post(lr(x));binary=np.where(x>=cut,post(lp),post(ln))
    delta=np.abs(own-binary)
    overlap=(x>=low)&(x<=high)
    counts=[int((delta>10).sum()),int(((delta>10)&overlap).sum()),int(((delta>10)&(x>high)).sum())]
    assert counts==[34,17,17]
    # Exact observed bounds govern classification; labels round to 861 and 1809.5.
    grid=np.unique(np.r_[np.geomspace(70,7000,900),low,high])
    bands=np.percentile(np.exp(co[:,0,None]+co[:,1,None]*np.log10(grid)-offset),[2.5,97.5],axis=0)
    blue,gray,orange="#0072B2","#808080","#D55E00"
    with plt.rc_context({"font.family":"DejaVu Sans","font.size":8,"axes.spines.top":False,
                         "axes.spines.right":False,"axes.edgecolor":"#999999","pdf.fonttype":42}):
        fig,axes=plt.subplots(3,1,figsize=(7.2,9.3),sharex=True,
                             gridspec_kw={"height_ratios":[1,1.38,1.38]})
        fig.subplots_adjust(left=.125,right=.985,top=.96,bottom=.065,hspace=.29)
        a,bb,c=axes
        for ax in axes:
            ax.set_xscale("log");ax.set_xlim(70,7000)
            ax.axvspan(low,high,color="#E7EEF5",zorder=0)
        rng=np.random.default_rng(7)
        for name,v,color,mark,row in [("PVOD",pv,blue,"D",2),("IPAH",ip,gray,"o",1),("HPAH",hp,gray,"o",0)]:
            a.scatter(v,row+rng.uniform(-.16,.16,len(v)),s=18,marker=mark,c=color,edgecolors="white",linewidths=.4,zorder=4)
            a.vlines(np.median(v),row-.25,row+.25,color="black",lw=2,zorder=5)
        a.set_yticks([0,1,2],["HPAH\nn = 42","IPAH\nn = 98","PVOD\nn = 26"])
        a.tick_params(axis="y",length=0,pad=9);a.set_ylim(-.55,2.7)
        a.axvline(cut,color=orange,lw=1.5,ls=(0,(5,3)))
        a.text(cut*.96,2.55,"5 of 26 PVOD below",ha="right",va="top",fontsize=7.5)
        a.text(cut*1.04,2.55,"published cutoff 1,658 pg/ml",color=orange,va="top",fontsize=7.5)
        a.text(cut*1.04,-.48,"4 of 140 comparators at or above",fontsize=7.3)
        a.set_title("A  Plasma GDF15 by diagnosis (bars: group medians)",loc="left",fontsize=10,fontweight="bold",pad=9)
        bb.set_yscale("log");bb.set_ylim(.001,10000)
        bb.fill_between(grid,bands[0],bands[1],color=blue,alpha=.14,lw=0,zorder=1)
        for mask,style in [(grid<=low,(0,(5,3))),((grid>=low)&(grid<=high),"-"),(grid>=high,(0,(5,3)))]:
            bb.plot(grid[mask],lr(grid[mask]),color=blue,lw=2,ls=style,zorder=3)
        bb.plot([70,cut,cut,7000],[ln,ln,lp,lp],color=orange,lw=1.7,ls=(0,(5,3)),zorder=4)
        for v,label in [(1,"LR = 1: no information"),(10,"LR = 10")]:
            bb.axhline(v,color="#CCCCCC",lw=.8)
            bb.text(6800,v*1.12,label,ha="right",color="#666666",fontsize=7.3)
        bb.legend(handles=[
            Line2D([],[],color=blue,lw=2,label="continuous LR (solid within observed overlap)"),
            Line2D([],[],color=blue,lw=2,ls=(0,(5,3)),label="continuous LR, extrapolated"),
            Line2D([],[],color=orange,lw=2,ls=(0,(5,3)),label="LR implied by the 1,658 pg/ml cutoff")
        ],loc="upper left",frameon=False,fontsize=7.1)
        bb.text(880,2500,"observed overlap of\nboth groups\n861–1,809.5 pg/ml",va="top",color="#527391",fontsize=7.4)
        bb.text(80,ln*1.35,f"cutoff-implied LR below 1,658 pg/ml: {ln:.2f}",color=orange,fontsize=7.2)
        bb.text(cut*.95,lp*1.35,f"cutoff-implied LR at or above 1,658 pg/ml: {lp:.1f} →",ha="right",color=orange,fontsize=7.1)
        bb.vlines(comp,.001,.00165,color=gray,lw=.6)
        bb.vlines(pv,.00165,.0027,color=blue,lw=.7)
        bb.set_yticks([.001,.01,.1,1,10,100,1000,10000],["0.001","0.01","0.1","1","10","100","1,000","10,000"])
        bb.set_ylabel("Likelihood ratio for PVOD (log scale)")
        bb.set_title("B  Continuous likelihood ratio and the likelihood ratio implied by the cutoff",loc="left",fontsize=8.8,fontweight="bold",pad=9)
        c.plot(grid,post(lr(grid)),color=blue,lw=1,alpha=.5)
        for mask,color,mark in [(~y,gray,"o"),(y,blue,"D")]:
            c.vlines(x[mask],np.minimum(own[mask],binary[mask]),np.maximum(own[mask],binary[mask]),color=color,lw=.6,alpha=.7)
            c.scatter(x[mask],binary[mask],facecolors="white",edgecolors=orange,s=20,marker=mark,lw=1,zorder=3)
            c.scatter(x[mask],own[mask],c=color,s=20,marker=mark,zorder=4)
        c.axvline(cut,color=orange,ls=(0,(5,3)),lw=1.5)
        c.text(80,96,"Filled: model-estimated probability for the measured value\n"
               "Open: probability implied by the cutoff label\n"
               f"Differ by >10 points in {counts[0]} of {len(x)} patients:\n"
               f"{counts[1]} within the observed overlap, {counts[2]} above the highest comparator",
               va="top",fontsize=7.15,linespacing=1.35)
        pair=np.array([x[x<cut].max(),x[x>=cut].min()])
        c.annotate(f"two patients {np.diff(pair)[0]:.0f} pg/ml apart:\n"
                   f"model estimate {post(lr(pair)).mean():.0f}% for both;\n"
                   f"the labels imply {post(ln):.0f}% and {post(lp):.0f}%",
                   xy=(pair.mean(),post(lr(pair)).mean()),xytext=(2300,43),fontsize=7.5,
                   arrowprops={"arrowstyle":"->","color":"#666666","lw":1})
        c.set_ylim(-4,104);c.set_yticks([0,20,40,60,80,100],["0%","20%","40%","60%","80%","100%"])
        c.set_ylabel("Post-test probability of PVOD (%)\nat 10% pretest probability")
        c.set_xlabel("Plasma GDF15 (pg/ml, log scale)")
        c.set_title("C  Post-test probability for each patient",loc="left",fontsize=10,fontweight="bold",pad=9)
        c.set_xticks([100,200,500,1000,2000,5000],["100","200","500","1,000","2,000","5,000"])
        for ax in axes:ax.tick_params(axis="x",which="minor",labelbottom=False)
        # Keep the pipeline's historical basename, plus submission filenames.
        fig.savefig(figures_dir/"figure1_lr_curve.pdf",metadata={"CreationDate":None})
        fig.savefig(figures_dir/"figure1_lr_curve.png",dpi=180)
        fig.savefig(figures_dir/"GDF15_Figure1.tif",dpi=600,pil_kwargs={"compression":"tiff_lzw"})
        fig.savefig(figures_dir/"GDF15_Figure1.pdf",metadata={"CreationDate":None})
        plt.close(fig)
    report={"seed":meta["seed"],"requested":4000,"finite":len(co),"separated":9,
            "coefficient_file_sha256":hashlib.sha256(coefpath.read_bytes()).hexdigest(),
            "observed_bounds":[float(low),float(high)],"discrepancies_gt10":counts,
            "bootstrap":"Archived ordinary logistic finite fits; no new resampling",
            "lr1_ci":np.percentile(10**((offset-co[:,0])/co[:,1]),[2.5,97.5]).tolist(),
            "lr10_ci":np.percentile(10**((np.log(10)+offset-co[:,0])/co[:,1]),[2.5,97.5]).tolist()}
    (figures_dir/"figure1_validation.json").write_text(json.dumps(report,indent=2)+"\n")

if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument("--results-dir",type=Path,default=Path(__file__).resolve().parents[1]/"results")
    p.add_argument("--figures-dir",type=Path,default=Path(__file__).resolve().parents[1]/"figures")
    p.add_argument("--workbook",type=Path,default=None)
    args=p.parse_args()
    draw(args.results_dir,args.figures_dir,args.workbook)
