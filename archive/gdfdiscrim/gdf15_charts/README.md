# GDF15 cutoff charts

All 322 patients plotted by their plasma GDF15 level, with two
horizontal lines marking where the diagnosis becomes certain.

## The two lines

- **861 pg/mL** — the lowest PVOD patient. Nobody below this line is
  PVOD. That clears **217 patients**, about 67% of everyone.
- **1,809 pg/mL** — the highest patient who is *not* PVOD (an IPAH case).
  Everyone above this line is PVOD: **20 patients**, which is 77% of all PVOD cases.
- **In between** sit **85 patients** (PVOD 6, IPAH 41, COPD 15, HPAH 10, Other PH 13) where GDF15 alone does not
  settle the diagnosis. That grey band is the real finding.

Both lines are just the highest and lowest actual patients — not thresholds
fitted by a model. Nothing here implies they are precise.

## The five charts

| File | What it shows | Points labelled | Violin-like shape |
|---|---|---|---|
| `Fig1_gdf15_cutoffs` | one group per diagnosis, with medians | no | yes |
| `Fig2_gdf15_cutoffs_symbols` | everyone pooled, symbols + legend | no | yes |
| `Fig3_gdf15_cutoffs_ranked` | everyone in order, lowest to highest | no | no |
| `Fig4_gdf15_cutoffs_scatter` | everyone scattered at random | no | no |
| `FigS1_gdf15_cutoffs_labelled` | every patient named individually | yes | no |

Each comes in `_log` and `_linear`, as PNG, PDF and SVG.

`Fig4_gdf15_cutoffs_scatter_linear_bold` is a second copy of the linear
scatter with the legend, the zone labels and the line captions enlarged,
for showing on a screen or a slide. The points and the lines are identical
to the plain version — only the labelling is louder.

- **log** spreads out the crowded 200–900 region so the grey band is readable.
- **linear** is true spacing, but squashes most patients into the bottom strip.

On every chart the **only real measurement is vertical** — height is the GDF15
level. Horizontal position is layout only, except in `Fig1` (where it is the
diagnosis group) and `Fig3` (where it is rank order).

## Also here

- `CAPTIONS.md` — formal figure captions with group sizes and statistics.
- `gdf15_ladder_points.csv` — the 322 plotted values, so any dot can be
  traced back to its number.

## Regenerating

```
venv/bin/python gdf15_ladder.py            # these screen-sized charts
venv/bin/python gdf15_ladder.py --print    # adds 183 mm / 89 mm journal sizes
```

Source data: Prabhakar et al., *J Clin Invest* 2026;136(2):e199159, Supplementary Source Data, sheet `Fig_1a` (panel Fig. 1a).
