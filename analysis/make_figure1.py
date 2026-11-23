"""Production adaptation of Giorgio Lagna's supplied three-panel Figure 1.

The unchanged supplied implementation is archived in qc/figure1_bootstrap/giorgio_original.py.
Public-workbook loading, verified point estimates, validated finite-fit bands and
computed annotations replace the original standalone input/fitting blocks.
"""
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import figure1_bootstrap as fb

def draw(results_dir, figures_dir, workbook=None):
    with mpl.rc_context():
        PVOD_C, PAH_C, CUT_C = "#0072B2", "#8A8A85", "#D55E00"
        INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
        mpl.rcParams.update({"font.family": "sans-serif", "font.size": 7.2,
            "axes.edgecolor": "#c3c2b7", "axes.linewidth": .6,
            "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2,
            "ytick.labelcolor": INK2, "xtick.major.width": .6, "ytick.major.width": .6,
            "axes.labelcolor": INK, "text.color": INK, "figure.facecolor": "white",
            "axes.facecolor": "white", "svg.fonttype": "none", "pdf.fonttype": 42})

        results = json.loads((results_dir / "results.json").read_text())
        entry = results["models"]["IPAH+HPAH"]
        x, y, (PV, IP, HP) = fb.cohort(workbook)
        COMP = np.concatenate([IP, HP])
        CUT = entry["threshold_comparison"]["published_cutoff"]
        PRE = .10
        b = np.array([entry["intercept"], entry["slope_log10_gdf15"]])
        off = entry["sampling_log_odds"]
        LR = lambda g: np.exp(b[0] + b[1] * np.log10(g) - off)
        sens = (PV >= CUT).mean()
        spec = (COMP < CUT).mean()
        LRP, LRN = sens / (1 - spec), (1 - sens) / spec
        post = lambda lr: 100 * (lr * PRE / (1 - PRE)) / (1 + lr * PRE / (1 - PRE))
        PP_POS, PP_NEG = post(LRP), post(LRN)
        XLO, XHI = 70., 7000.
        OBS_LO_PVOD, OBS_HI_COMP = PV.min(), COMP.max()
        coefficients = pd.read_csv(results_dir / "figure1_bootstrap_coefficients.csv", float_precision="round_trip").to_numpy()
        accounting = json.loads((results_dir / "figure1_bootstrap.json").read_text())
        assert accounting["seed"] == 20260907 and accounting["separated"] == 9
        assert coefficients.shape == (3991, 2)
        grid = np.unique(np.r_[np.geomspace(XLO, XHI, 700), OBS_LO_PVOD, OBS_HI_COMP])
        lo, hi = np.percentile(fb.lr(coefficients, grid, off), [2.5, 97.5], axis=0)

        fig = plt.figure(figsize=(7.1, 7.9))
        gs = fig.add_gridspec(3, 1, height_ratios=[.80, 1.25, 1.15], hspace=.42,
                              left=.115, right=.975, top=.965, bottom=.062)
        axA, axB, axC = (fig.add_subplot(g) for g in gs)

        # ---------------------------------------------------------------- A
        groups = [("PVOD", PV, PVOD_C, "D", 2), ("IPAH", IP, PAH_C, "o", 1), ("HPAH", HP, PAH_C, "o", 0)]
        r = np.random.default_rng(7)
        for name, v, c, mk, ypos in groups:
            axA.scatter(v, ypos + r.uniform(-.15, .15, len(v)), s=11, marker=mk, facecolor=c,
                        edgecolor="white", linewidth=.35, alpha=.95, zorder=4)
            axA.plot([np.median(v)] * 2, [ypos - .30, ypos + .30], color=INK, lw=1.4, zorder=6)
        axA.axvline(CUT, color=CUT_C, lw=1.0, ls=(0, (4, 2)), zorder=3)
        axA.text(CUT * 1.04, 2.8, f"published cutoff {CUT:,.0f} pg/mL", fontsize=7, color=CUT_C, va="top")
        axA.text(CUT * .95, 2.8, f"{int((PV < CUT).sum())} of {len(PV)} PVOD below", fontsize=7, ha="right", va="top")
        axA.text(CUT * 1.04, -.35, f"{int((COMP >= CUT).sum())} of {len(COMP)} comparators at or above", fontsize=7, va="center")
        axA.set_ylim(-.55, 2.9)
        axA.set_yticks([2, 1, 0], [f"PVOD\nn = {len(PV)}", f"IPAH\nn = {len(IP)}", f"HPAH\nn = {len(HP)}"])
        axA.tick_params(axis="y", length=0)
        axA.set_title("A", loc="left", fontsize=10, fontweight="bold", pad=4)
        for s in ("left", "right", "top"): axA.spines[s].set_visible(False)

        # ---------------------------------------------------------------- B
        sup = (grid >= OBS_LO_PVOD) & (grid <= OBS_HI_COMP)
        axB.fill_between(grid[sup], lo[sup], hi[sup], color=PVOD_C, alpha=.16, lw=0, zorder=3)
        for tail in (grid < OBS_LO_PVOD, grid > OBS_HI_COMP):
            axB.fill_between(grid[tail], lo[tail], hi[tail], color=PVOD_C, alpha=.07, lw=0, zorder=3)
        axB.plot(grid[sup], LR(grid[sup]), color=PVOD_C, lw=2, zorder=5,
                 label="continuous LR (solid where both groups were observed)")
        for tail in (grid <= OBS_LO_PVOD, grid >= OBS_HI_COMP):
            axB.plot(grid[tail], LR(grid[tail]), color=PVOD_C, lw=1.4, ls=(0, (2, 2)), alpha=1, zorder=5, label="continuous LR, extrapolated" if grid[tail][0] == grid[0] else None)
        step_x = np.array([XLO, CUT, CUT, XHI])
        step_y = np.array([LRN, LRN, LRP, LRP])
        axB.plot(step_x, step_y, color=CUT_C, lw=1.6, ls=(0, (5, 2)), zorder=6,
                 label="LR implied by the 1,658 pg/mL cutoff")
        for L, tx in ((10, "LR = 10"), (1, "LR = 1: no information")):
            axB.axhline(L, color=GRID, lw=.7, zorder=1)
            axB.text(XHI * .97, L * .85, tx, fontsize=6.4, color=INK2, ha="right", va="top")
        axB.set_ylim(.001, 10000)
        for v, c, base in ((COMP, PAH_C, .0012), (PV, PVOD_C, .0022)):
            axB.scatter(v, np.full(v.size, base), marker="|", s=26, color=c, lw=.7, zorder=4)
        axB.text(OBS_LO_PVOD * 1.02, 2700, "observed overlap of\nboth groups\n861–1,809.5 pg/mL", fontsize=6.3, color="#38608a", va="top")
        axB.text(80, LRN * 1.35, f"cutoff-implied LR below 1,658 pg/mL: {LRN:.2f}", fontsize=6.5, color=CUT_C)
        axB.text(CUT * .95, LRP * 1.35, f"cutoff-implied LR at or above 1,658 pg/mL: {LRP:.1f} →", ha="right", fontsize=6.5, color=CUT_C)
        axB.set_yscale("log")
        axB.set_ylabel("Likelihood ratio for PVOD (log scale)", fontsize=7.2)
        axB.set_yticks([.001, .01, .1, 1, 10, 100, 1000, 10000], ["0.001", "0.01", "0.1", "1", "10", "100", "1,000", "10,000"])
        axB.set_title("B", loc="left", fontsize=10, fontweight="bold", pad=4)
        axB.legend(loc="upper left", fontsize=6.3, frameon=False, borderpad=0, handlelength=2.4,
                   bbox_to_anchor=(.005, .995))

        # ---------------------------------------------------------------- C
        own_pp = post(LR(x))
        cut_pp = np.where(x >= CUT, PP_POS, PP_NEG)
        is_pv = y == 1
        for m, c, mk in ((~is_pv, PAH_C, "o"), (is_pv, PVOD_C, "D")):
            axC.vlines(x[m], np.minimum(own_pp[m], cut_pp[m]), np.maximum(own_pp[m], cut_pp[m]),
                       color=c, lw=.55, alpha=.55, zorder=3)
            axC.scatter(x[m], own_pp[m], s=13, marker=mk, facecolor=c, edgecolor="white", lw=.35, zorder=5)
            axC.scatter(x[m], cut_pp[m], s=15, marker=mk, facecolor="white", edgecolor=CUT_C, lw=.8, zorder=4)
        axC.plot(grid, post(LR(grid)), color=PVOD_C, lw=1.2, alpha=.65, zorder=2)
        axC.axvline(CUT, color=CUT_C, lw=1.0, ls=(0, (4, 2)), zorder=1)
        big = np.abs(own_pp - cut_pp) > 10
        near_cutoff = big & (np.abs(x - CUT) <= 200)
        above_comparator = big & (x > OBS_HI_COMP)
        axC.text(78, 63, f"Differ by >10 points in {int(big.sum())} of {len(x)} patients:\n"
                         f"{int(near_cutoff.sum())} near the cutoff, "
                         f"{int(above_comparator.sum())} above the highest comparator",
                 fontsize=6.4, color=INK2, va="top", linespacing=1.5)
        pair = np.array([x[x < CUT].max(), x[x >= CUT].min()])
        pair_pp = post(LR(pair))
        axC.annotate(f"two patients {np.diff(pair)[0]:.0f} pg/mL apart:\n"
                     f"model estimate {pair_pp.mean():.0f}% for both,\n"
                     f"the labels imply {PP_NEG:.0f}% and {PP_POS:.0f}%",
                     xy=(1658, 39), xytext=(2350, 45), fontsize=6.4, color=INK, linespacing=1.45,
                     arrowprops=dict(arrowstyle="->", color=MUTED, lw=.7, connectionstyle="arc3,rad=.2"))
        axC.set_ylim(-4, 104); axC.set_ylabel("Post-test probability of PVOD (%)\nat 10% pretest probability",
                                              fontsize=7.2, linespacing=1.5)
        axC.set_xlabel("Plasma GDF15 (pg/mL, log scale)", fontsize=7.8)
        axC.set_title("C", loc="left", fontsize=10, fontweight="bold", pad=4)

        key = [
            Line2D([], [], marker="D", ls="", mfc=PVOD_C, mec="white", mew=.35,
                   ms=4.4, label="PVOD, model-estimated probability"),
            Line2D([], [], marker="o", ls="", mfc=PAH_C, mec="white", mew=.35,
                   ms=4.4, label="IPAH/HPAH, model-estimated probability"),
            Line2D([], [], marker="D", ls="", mfc="white", mec=CUT_C, mew=.8,
                   ms=4.4, label="PVOD, probability implied by cutoff label"),
            Line2D([], [], marker="o", ls="", mfc="white", mec=CUT_C, mew=.8,
                   ms=4.4, label="IPAH/HPAH, probability implied by cutoff label"),
        ]
        axC.legend(handles=key, title="Key", loc="upper left", fontsize=6.1,
                   title_fontsize=6.4, alignment="left", frameon=False, borderpad=0,
                   handletextpad=.45, labelspacing=.35)

        for ax in (axA, axB, axC):
            ax.axvspan(OBS_LO_PVOD, OBS_HI_COMP, color="#e7eef5", zorder=0)
            ax.set_xscale("log"); ax.set_xlim(XLO, XHI)
            ax.set_xticks([100, 200, 500, 1000, 2000, 5000])
            ax.set_xticklabels(["100", "200", "500", "1,000", "2,000", "5,000"])
            ax.tick_params(axis="x", length=3)
        for ax in (axA, axB): ax.set_xticklabels([])
        for s in ("right", "top"): axB.spines[s].set_visible(False); axC.spines[s].set_visible(False)

        figures_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(figures_dir / "figure1_lr_curve.pdf", dpi=600,
                    metadata={"CreationDate": None})
        fig.savefig(figures_dir / "figure1_lr_curve.tif", dpi=600,
                    pil_kwargs={"compression": "tiff_lzw"})
        # Retain the repository's PNG preview alongside the journal files.
        fig.savefig(figures_dir / "figure1_lr_curve.png", dpi=300)
        plt.close(fig)
