"""Verify the production Figure 1 from public observations and finite MLEs."""
import json
import numpy as np
import pandas as pd
import figure1_bootstrap as fb


def verify(results_dir, check):
    report = json.loads((results_dir / 'figure1_bootstrap.json').read_text())
    main = json.loads((results_dir / 'results.json').read_text())['models']['IPAH+HPAH']
    curve = pd.read_csv(results_dir / 'figure1_band.csv', float_precision='round_trip')
    saved = pd.read_csv(results_dir / 'figure1_bootstrap_coefficients.csv', float_precision='round_trip').to_numpy()
    x, y, groups = fb.cohort()
    repeated, records = fb.bootstrap(x, y)
    check('Figure 1: 4000 requested = 3991 finite + 9 separated',
          report['requested'] == 4000 and report['finite'] == len(saved) == 3991 and report['separated'] == 9)
    check('Figure 1: repeated finite coefficients byte-identical', np.array_equal(saved, repeated))
    rejected = [r for r in records if not r['accepted']]
    check('Figure 1: rejected identities match production accounting',
          [r['replicate'] for r in rejected] == report['separated_replicates_zero_based'])
    check('Figure 1: every excluded sample strictly separated', all(
          x[np.array(r['indices'])][y[np.array(r['indices'])] == 1].min() >
          x[np.array(r['indices'])][y[np.array(r['indices'])] == 0].max() for r in rejected))
    check('Figure 1: all accepted ordinary MLEs pass score/information validation',
          all(r['standardized_score_max_abs'] < 1e-7 and r['information_min_eigenvalue'] > 0
              for r in records if r['accepted']))
    expected = np.percentile(fb.lr(saved, curve.gdf15.to_numpy(), main['sampling_log_odds']), [2.5, 97.5], axis=0)
    check('Figure 1: LR-scale percentile band uses exactly finite fits',
          np.allclose(expected, curve[['lr_ci_low', 'lr_ci_high']].to_numpy().T, rtol=1e-13, atol=0))
    check('Figure 1: point curve uses unchanged production coefficients',
          np.allclose(curve.lr, fb.lr([main['intercept'], main['slope_log10_gdf15']],
              curve.gdf15.to_numpy(), main['sampling_log_odds']), rtol=1e-13, atol=0))
    check('Figure 1: exactly 166 observed group-1 values', len(x) == 166 and list(map(len, groups)) == [26, 98, 42])
    cutoff = main['threshold_comparison']['published_cutoff']
    pair = [x[x < cutoff].max(), x[x >= cutoff].min()]
    observed = json.loads((results_dir / 'centerpiece_verification.json').read_text())
    check('Figure 1: computed adjacent pair is Table 1 observed pair',
          pair == [float(r['raw_concentration']) for r in observed['adjacent_observations']])
    check('Figure 1: computed support boundaries match verified model',
          groups[0].min() == main['support']['pvod_min_gdf15'] and
          np.concatenate(groups[1:]).max() == main['support']['comparator_max_gdf15'])
    check('Figure 1: approved Methods reports finite fits and narrows in tail',
          report['methods_wording'] == fb.METHODS and 'narrows' in fb.METHODS)
