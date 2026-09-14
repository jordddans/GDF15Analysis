"""
Figure 1, v2, for the GDF15 / PVOD evidence-curve paper.

Every element in this figure is observed data. There are no hypothetical patients.

Usage:      python make_figure1_v2.py [path/to/supporting_data.xlsx]

Input:      the Prabhakar et al. supporting-data workbook. Only sheet GDF15-1 is read,
            which holds raw plasma GDF15 as six group columns. The original and the
            cleaned copy give identical numbers (audited: 0 differences in 4,609 cells),
            so either file works.

Outputs:    Figure1_GDF15_evidence_v2.pdf   (vector, for the journal)
            Figure1_GDF15_evidence_v2.png   (300 dpi)
            and it prints the key numbers to stdout so they can be checked against an
            independent reconstruction.

Requires:   numpy, scipy, matplotlib, openpyxl.
"""
import sys, glob
import numpy as np, openpyxl, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import minimize

PVOD_C, PAH_C, CUT_C = "#0072B2", "#8A8A85", "#D55E00"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
mpl.rcParams.update({"font.family": "sans-serif", "font.size": 7.2,
    "axes.edgecolor": "#c3c2b7", "axes.linewidth": .6,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2,
    "ytick.labelcolor": INK2, "xtick.major.width": .6, "ytick.major.width": .6,
    "axes.labelcolor": INK, "text.color": INK, "figure.facecolor": "white",
    "axes.facecolor": "white", "svg.fonttype": "none", "pdf.fonttype": 42})

# ---------------------------------------------------------------- data
def find_workbook():
    """Use the path given on the command line, else the first .xlsx here that has GDF15-1."""
    if len(sys.argv) > 1:
        return sys.argv[1]
    for f in sorted(glob.glob("*.xlsx")):
        try:
            if "GDF15-1" in openpyxl.load_workbook(f, read_only=True).sheetnames:
                return f
        except Exception:
            continue
    sys.exit("No workbook with a GDF15-1 sheet found here. "
             "Pass the path: python make_figure1_v2.py <workbook.xlsx>")

WB = find_workbook()
print(f"reading {WB}, sheet GDF15-1")
ws = openpyxl.load_workbook(WB, data_only=True)["GDF15-1"]
hdr = [c for c in next(ws.iter_rows(min_row=3, max_row=3, values_only=True)) if c]
d = {h: [] for h in hdr}
for r in ws.iter_rows(min_row=4, values_only=True):
    for i, h in enumerate(hdr):
        if isinstance(r[i], (int, float)):
            d[h].append(float(r[i]))
PV = np.array(d["PVOD"]); IP = np.array(d["IPAH"]); HP = np.array(d["HPAH"])
COMP = np.concatenate([IP, HP])
x = np.concatenate([PV, COMP]); y = np.concatenate([np.ones(len(PV)), np.zeros(len(COMP))])
CUT, PRE = 1658.0, 0.10

def fit(xx, yy):
    lx = np.log10(xx)
    f = lambda b: -np.sum(yy * (b[0] + b[1] * lx) - np.log1p(np.exp(b[0] + b[1] * lx)))
    return minimize(f, [0, 1], method="BFGS").x

b = fit(x, y)
off = np.log((len(PV) / len(x)) / (1 - len(PV) / len(x)))
LR = lambda g, bb=b, o=off: np.exp(bb[0] + bb[1] * np.log10(g) - o)

sens = (PV >= CUT).mean(); spec = (COMP < CUT).mean()
LRP, LRN = sens / (1 - spec), (1 - sens) / spec
post = lambda lr, pre=PRE: 100 * (lr * pre / (1 - pre)) / (1 + lr * pre / (1 - pre))
PP_POS, PP_NEG = post(LRP), post(LRN)

# bootstrap band, 4000 stratified case replicates
rng = np.random.default_rng(20260907)
grid = np.logspace(np.log10(700), np.log10(6500), 400)
B = np.empty((4000, grid.size))
for i in range(4000):
    xs = np.concatenate([rng.choice(PV, len(PV), True), rng.choice(COMP, len(COMP), True)])
    B[i] = LR(grid, fit(xs, y))
lo, hi = np.percentile(B, [2.5, 97.5], axis=0)

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
axA.text(CUT * 1.05, 3.25, "published cutoff\n1,658 pg/mL", fontsize=6.4, color=CUT_C,
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
sup = grid >= OBS_LO_PVOD
axB.fill_between(grid[sup], lo[sup], hi[sup], color=PVOD_C, alpha=.16, lw=0, zorder=3)
axB.fill_between(grid[~sup], lo[~sup], hi[~sup], color=PVOD_C, alpha=.07, lw=0, zorder=3)
axB.plot(grid[sup], LR(grid[sup]), color=PVOD_C, lw=2, zorder=5,
         label="evidence from the measured value (model)")
axB.plot(grid[~sup], LR(grid[~sup]), color=PVOD_C, lw=1.4, ls=(0, (2, 2)), alpha=.55, zorder=5)
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
axB.text(OBS_LO_PVOD * .90, .022, "861 pg/mL: lowest observed\nPVOD value. Left of this the\ncurve is extrapolation and is\ndrawn broken.",
         fontsize=6.1, color=INK2, ha="right", va="top", linespacing=1.4)
axB.text(XLO * 1.25, LRN * 1.9, "the rule treats 91 and 1,657 pg/mL\nas the same measurement",
         fontsize=6.3, color=CUT_C, ha="left", va="bottom", linespacing=1.4)
axB.annotate(f"the rule assigns every value above the cutoff\nthe same likelihood ratio, {LRP:.1f}",
             xy=(3100, LRP), xytext=(1050, 300), fontsize=6.3, color=CUT_C, linespacing=1.4,
             arrowprops=dict(arrowstyle="->", color=CUT_C, lw=.6,
                             connectionstyle="arc3,rad=-.25"))
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
pair = np.array([1657.81, 1658.81])
axC.annotate("two patients 1 pg/mL apart:\nboth values support 22%,\nthe rule reports 2% and 76%",
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

fig.savefig("Figure1_GDF15_evidence_v2.pdf")
fig.savefig("Figure1_GDF15_evidence_v2.png", dpi=300)

print("\n--- numbers to check against an independent build ---")
print(f"n: PVOD {len(PV)}, IPAH {len(IP)}, HPAH {len(HP)}, comparators {len(COMP)}")
print(f"model: logit P(PVOD) = {b[0]:.3f} + {b[1]:.3f} * log10(GDF15)")
print(f"sens {sens:.4f} spec {spec:.4f}  LR+ {LRP:.2f}  LR- {LRN:.3f}")
print(f"post-test at 10%: positive label {PP_POS:.1f}%  negative label {PP_NEG:.1f}%")
for v in (1657.81, 1658.81):
    print(f"  {v}: own LR {LR(v):.3f} -> {post(LR(v)):.1f}% ; rule -> {PP_POS if v>=CUT else PP_NEG:.1f}%")
print("differ >10 pp:", int(big.sum()), "of", len(x), "; >20 pp:", int((np.abs(own_pp-cut_pp)>20).sum()))
print("PVOD below cutoff:", int((PV < CUT).sum()), "comparators >= cutoff:", int((COMP >= CUT).sum()))

# the adjacent pair that Table 1 rests on
order = np.argsort(np.abs(x - CUT))
print("\nfive values closest to the cutoff:")
for k in order[:5]:
    print(f"  {x[k]:10.5f}  {'PVOD' if y[k] else 'comparator':10s}"
          f"  own LR {LR(x[k]):7.3f} -> {post(LR(x[k])):5.1f}%"
          f"   rule -> {PP_POS if x[k] >= CUT else PP_NEG:5.1f}%")
gaps = np.abs(PV[:, None] - COMP[None, :])
i, j = np.unravel_index(gaps.argmin(), gaps.shape)
print(f"minimum case-to-comparator gap anywhere: {gaps.min():.4f} pg/mL "
      f"(PVOD {PV[i]}, comparator {COMP[j]}); cutoff {CUT} falls inside it: "
      f"{min(PV[i], COMP[j]) < CUT <= max(PV[i], COMP[j])}")
