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
        band = pd.read_csv(results_dir / "figure1_band.csv", float_precision="round_trip")
        grid, lo, hi = (band[k].to_numpy() for k in ("gdf15", "lr_ci_low", "lr_ci_high"))
        XLO, XHI = 70., 7000.
        OBS_LO_PVOD, OBS_HI_COMP = PV.min(), COMP.max()

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
            axA.text(XLO * 1.05, ypos, f"{name}\nn = {len(v)}", fontsize=6.8, va="center",
                     color=PVOD_C if name == "PVOD" else INK2,
                     fontweight="bold" if name == "PVOD" else "normal", linespacing=1.3)
        axA.axvline(CUT, color=CUT_C, lw=1.0, ls=(0, (4, 2)), zorder=3)
        axA.text(CUT * 1.05, 3.25, f"published cutoff\n{CUT:,.0f} pg/mL", fontsize=6.4, color=CUT_C,
                 ha="left", va="top", linespacing=1.35, fontweight="bold")
        axA.annotate(f"{int((PV < CUT).sum())} of {len(PV)} PVOD cases fall below the cutoff",
                     xy=(1560, 2.45), fontsize=6.3, color=INK2, ha="right")
        axA.annotate(f"{int((COMP >= CUT).sum())} of {len(COMP)} comparators fall at or above it",
                     xy=(1760, -.42), fontsize=6.3, color=INK2, ha="left", va="center")
        axA.set_ylim(-.65, 3.35); axA.set_yticks([])
        axA.set_title("A   Observed plasma GDF15 by diagnosis. Vertical bars are group medians.",
                      loc="left", fontsize=8, fontweight="bold", pad=4)
        for s in ("left", "right", "top"): axA.spines[s].set_visible(False)

        # ---------------------------------------------------------------- B
        sup = (grid >= OBS_LO_PVOD) & (grid <= OBS_HI_COMP)
        axB.fill_between(grid[sup], lo[sup], hi[sup], color=PVOD_C, alpha=.16, lw=0, zorder=3)
        for tail in (grid < OBS_LO_PVOD, grid > OBS_HI_COMP):
            axB.fill_between(grid[tail], lo[tail], hi[tail], color=PVOD_C, alpha=.07, lw=0, zorder=3)
        axB.plot(grid[sup], LR(grid[sup]), color=PVOD_C, lw=2, zorder=5,
                 label="evidence from the measured value (model)")
        for tail in (grid < OBS_LO_PVOD, grid > OBS_HI_COMP):
            axB.plot(grid[tail], LR(grid[tail]), color=PVOD_C, lw=1.4, ls=(0, (2, 2)), alpha=.55, zorder=5)
        step_x = np.array([XLO, CUT, CUT, XHI])
        step_y = np.array([LRN, LRN, LRP, LRP])
        axB.plot(step_x, step_y, color=CUT_C, lw=1.6, ls=(0, (5, 2)), zorder=6,
                 label="evidence assigned by the 1,658 pg/mL rule")
        for L, tx in ((10, "strong evidence for PVOD"), (1, "no information"), (.1, "strong evidence against")):
            axB.axhline(L, color=GRID, lw=.9, zorder=1)
            axB.text(XHI * .97, L * 1.30, tx, fontsize=6, color=MUTED, ha="right")
        ylo, yhi = 3e-4, 4e3
        axB.set_ylim(ylo, yhi)
        for v, c, base in ((COMP, PAH_C, ylo * 1.25), (PV, PVOD_C, ylo * 2.6)):
            axB.scatter(v, np.full(v.size, base), marker="|", s=26, color=c, lw=.7, zorder=4, alpha=.85)
        axB.text(XLO * 1.10, ylo * 4.6, "observed values", fontsize=6, color=MUTED)
        axB.axvline(OBS_LO_PVOD, color=INK2, lw=.6, ls=(0, (1.5, 2)), zorder=2)
        axB.text(OBS_LO_PVOD * .90, .022, f"{OBS_LO_PVOD:,.0f} pg/mL: lowest observed\nPVOD value. Left of this the\ncurve is extrapolation and is\ndrawn broken.",
                 fontsize=6.1, color=INK2, ha="right", va="top", linespacing=1.4)
        axB.text(XLO * 1.25, LRN * 1.9, f"the rule treats {x.min():,.1f} and {x[x < CUT].max():,.1f} pg/mL\nas the same measurement",
                 fontsize=6.3, color=CUT_C, ha="left", va="bottom", linespacing=1.4)
        axB.annotate(f"the rule assigns every value above the cutoff\nthe same likelihood ratio, {LRP:.1f}",
                     xy=(3100, LRP), xytext=(1050, 300), fontsize=6.3, color=CUT_C, linespacing=1.4,
                     arrowprops=dict(arrowstyle="->", color=CUT_C, lw=.6,
                                     connectionstyle="arc3,rad=-.25"))
        axB.axvline(OBS_HI_COMP, color=INK2, lw=.6, ls=(0, (1.5, 2)), zorder=2)
        axB.text(OBS_HI_COMP * 1.07, .018,
                 f"{OBS_HI_COMP:,.1f} pg/mL: highest comparator\n"
                 "value. Curve is extrapolated to the right.",
                 fontsize=6.1, color=INK2, ha="left", va="top", linespacing=1.4)
        axB.set_yscale("log"); axB.set_ylabel("Likelihood ratio for PVOD", fontsize=7.2)
        axB.set_yticks([1e-3, 1e-2, .1, 1, 10, 100, 1000])
        axB.set_yticklabels(["0.001", "0.01", "0.1", "1", "10", "100", "1000"])
        axB.set_title("B   What a measured value is worth, against what the rule assigns it",
                      loc="left", fontsize=8, fontweight="bold", pad=4)
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
        axC.text(78, 92, f"Filled: probability supported by the patient's own value.\n"
                         f"Open: probability implied by the cutoff label.\n"
                         f"The two differ by more than 10 percentage points in\n"
                         f"{int(big.sum())} of {len(x)} patients.",
                 fontsize=6.4, color=INK2, va="top", linespacing=1.5)
        pair = np.array([x[x < CUT].max(), x[x >= CUT].min()])
        pair_pp = post(LR(pair))
        axC.annotate(f"two patients {np.diff(pair)[0]:.0f} pg/mL apart:\n"
                     f"both values support {pair_pp.mean():.0f}%,\n"
                     f"the rule reports {PP_NEG:.0f}% and {PP_POS:.0f}%",
                     xy=(1658, 39), xytext=(2350, 45), fontsize=6.4, color=INK, linespacing=1.45,
                     arrowprops=dict(arrowstyle="->", color=MUTED, lw=.7, connectionstyle="arc3,rad=.2"))
        axC.set_ylim(-4, 100); axC.set_ylabel("Probability of PVOD after the test (%)\nat 10% pre-test suspicion",
                                              fontsize=7.2, linespacing=1.5)
        axC.set_xlabel("Plasma GDF15 (pg/mL)", fontsize=7.8)
        axC.set_title("C   The consequence for every patient in the cohort",
                      loc="left", fontsize=8, fontweight="bold", pad=4)

        for ax in (axA, axB, axC):
            ax.set_xscale("log"); ax.set_xlim(XLO, XHI)
            ax.set_xticks([100, 200, 500, 1000, 2000, 5000])
            ax.set_xticklabels(["100", "200", "500", "1,000", "2,000", "5,000"])
            ax.tick_params(axis="x", length=3)
        for ax in (axA, axB): ax.set_xticklabels([])
        for s in ("right", "top"): axB.spines[s].set_visible(False); axC.spines[s].set_visible(False)

        leg = [Line2D([], [], marker="D", ls="", mfc=PVOD_C, mec="white", ms=4.4, label="PVOD"),
               Line2D([], [], marker="o", ls="", mfc=PAH_C, mec="white", ms=4.4, label="IPAH or HPAH")]
        axA.legend(handles=leg, loc="lower left", fontsize=6.4, frameon=False, ncol=2,
                   bbox_to_anchor=(.005, -.06), handletextpad=.3, columnspacing=1.4)

        figures_dir.mkdir(parents=True, exist_ok=True)
        for suffix in ("pdf", "png"):
            fig.savefig(figures_dir / f"figure1_lr_curve.{suffix}", dpi=300,
                        metadata={"CreationDate": None} if suffix == "pdf" else None)
        plt.close(fig)
