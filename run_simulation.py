"""Run the full simulation study."""
import sys, io, os, time, warnings, contextlib  # noqa: E401
warnings.filterwarnings('ignore')
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
sys.path.insert(0, os.path.dirname(__file__))

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

from ipd_qma import IPDQMA, simulate_location_scale

OUT_DIR = os.path.join(os.path.dirname(__file__), 'output')
os.makedirs(OUT_DIR, exist_ok=True)

N_SIM = 1000

if __name__ != "__main__":
    raise SystemExit("This module should be run directly, not imported.")

scenarios = []

# Type I error (null): VR = 1.0
for dist in ['normal', 'exponential', 'lognormal', 't5']:
    scenarios.append({'name': 'Null_' + dist, 'type': 'null', 'distribution': dist,
                      'variance_ratio': 1.0, 'K': 10, 'n_range': (80, 200)})

# Power by VR
for vr in [1.5, 2.0, 3.0, 5.0]:
    scenarios.append({'name': 'Power_VR' + str(vr), 'type': 'power', 'distribution': 'normal',
                      'variance_ratio': vr, 'K': 10, 'n_range': (80, 200)})

# Power by K
for K in [3, 5, 10, 20]:
    scenarios.append({'name': 'Power_K' + str(K), 'type': 'power', 'distribution': 'normal',
                      'variance_ratio': 2.0, 'K': K, 'n_range': (80, 200)})

# Power by N
for n_max in [50, 100, 200]:
    scenarios.append({'name': 'Power_N' + str(n_max), 'type': 'power', 'distribution': 'normal',
                      'variance_ratio': 2.0, 'K': 10, 'n_range': (max(30, n_max // 2), n_max)})

# Mixed location + scale shift
scenarios.append({'name': 'Mixed_MD0.5_VR2', 'type': 'mixed', 'distribution': 'normal',
                  'variance_ratio': 2.0, 'K': 10, 'n_range': (80, 200), 'mean_shift': 0.5})

all_results = []
t_total = time.time()

for sc_idx, sc in enumerate(scenarios):
    print('[%d/%d] %s (dist=%s, VR=%s, K=%d, N=%s)' %
          (sc_idx + 1, len(scenarios), sc['name'], sc['distribution'],
           sc['variance_ratio'], sc['K'], sc['n_range']))
    t0 = time.time()

    reject_md = reject_q90 = reject_slope = reject_lnvr = 0
    valid_sims = 0
    n_failed = 0

    for sim in range(N_SIM):
        seed = 10000 * sc_idx + sim
        try:
            studies_data, labels, _ = simulate_location_scale(
                K=sc['K'], n_range=sc['n_range'], variance_ratio=sc['variance_ratio'],
                mean_shift=sc.get('mean_shift', 0.0),
                distribution=sc['distribution'], seed=seed)

            model = IPDQMA(quantiles=[0.1, 0.5, 0.9], n_boot=200, conf_level=0.95,
                           use_hksj=True, seed=seed + 1_000_000)
            with contextlib.redirect_stdout(io.StringIO()):
                results = model.fit(studies_data, labels=labels)

            valid_sims += 1
            if results['md_test']['p_value'] < 0.05:
                reject_md += 1
            q90_row = results['profile'][results['profile']['Quantile'] == 0.9]
            q90_p = q90_row['P'].values[0]
            if q90_p < 0.05:
                reject_q90 += 1
            if results['slope_test']['p_value'] < 0.05:
                reject_slope += 1
            if np.isfinite(results['lnvr_test']['p_value']) and results['lnvr_test']['p_value'] < 0.05:
                reject_lnvr += 1
        except Exception as e:
            n_failed += 1
            if n_failed == 1:
                print('  First failure: %s' % e)

    n_failed = N_SIM - valid_sims
    elapsed = time.time() - t0
    if n_failed > 0:
        print('  WARNING: %d/%d simulations failed (%.1f%%)' %
              (n_failed, N_SIM, n_failed / N_SIM * 100))
    n = max(valid_sims, 1)
    rates = {'MD': reject_md / n, 'Q90': reject_q90 / n,
             'Slope': reject_slope / n, 'lnVR': reject_lnvr / n}
    mcse = {k: np.sqrt(v * (1 - v) / n) for k, v in rates.items()}

    row = {
        'Scenario': sc['name'], 'Type': sc['type'], 'Distribution': sc['distribution'],
        'VR': sc['variance_ratio'], 'K': sc['K'],
        'N_range': '%d-%d' % (sc['n_range'][0], sc['n_range'][1]),
        'N_sim': valid_sims,
        'Reject_MD': rates['MD'], 'Reject_Q90': rates['Q90'],
        'Reject_Slope': rates['Slope'], 'Reject_lnVR': rates['lnVR'],
        'MCSE_MD': mcse['MD'], 'MCSE_Q90': mcse['Q90'],
        'MCSE_Slope': mcse['Slope'], 'MCSE_lnVR': mcse['lnVR'],
        'Time_sec': elapsed,
    }
    all_results.append(row)

    print('  %ds | MD=%.3f Q90=%.3f Slope=%.3f lnVR=%.3f' %
          (elapsed, rates['MD'], rates['Q90'], rates['Slope'], rates['lnVR']))

# Save results
sim_df = pd.DataFrame(all_results)
sim_df.to_csv(os.path.join(OUT_DIR, 'table_simulation_results.csv'), index=False)
print('\nTotal time: %.1f min' % ((time.time() - t_total) / 60))

# ---- Generate figures ----

methods = ['Reject_MD', 'Reject_Q90', 'Reject_Slope', 'Reject_lnVR']
labels_m = ['Mean Diff', 'IPD-QMA (Q90)', 'Slope Test', 'lnVR']
colors = ['#3498DB', '#E74C3C', '#2ECC71', '#F39C12']
markers = ['o', 's', '^', 'D']

# Fig 4: Type I error
null_df = sim_df[sim_df['Type'] == 'null']
fig4, ax4 = plt.subplots(figsize=(10, 5))
x = np.arange(len(null_df))
width = 0.2
for i, (m, l, c) in enumerate(zip(methods, labels_m, colors)):
    ax4.bar(x + i * width, null_df[m].values, width, label=l, color=c, alpha=0.8)
ax4.axhline(0.05, color='red', linestyle='--', linewidth=2, label='Nominal alpha=0.05')
ax4.axhspan(0.025, 0.075, color='red', alpha=0.1)
ax4.set_xlabel('Distribution')
ax4.set_ylabel('Rejection Rate')
ax4.set_title('Type I Error Under True Null (VR=1.0, K=10)')
ax4.set_xticks(x + 1.5 * width)
ax4.set_xticklabels(null_df['Distribution'].values)
ax4.legend()
ax4.grid(True, alpha=0.3, axis='y')
fig4.savefig(os.path.join(OUT_DIR, 'fig4_type_i_error.png'), dpi=150, bbox_inches='tight')
plt.close('all')

# Fig 5: Power by VR
vr_df = sim_df[sim_df['Scenario'].str.startswith('Power_VR')]
fig5, ax5 = plt.subplots(figsize=(10, 6))
vrs = vr_df['VR'].values
for m, l, c, mk in zip(methods, labels_m, colors, markers):
    ax5.plot(vrs, vr_df[m].values, '-' + mk, color=c, label=l, linewidth=2, markersize=8)
ax5.axhline(0.80, color='gray', linestyle=':', alpha=0.5, label='80% power')
ax5.axhline(0.05, color='red', linestyle='--', alpha=0.3)
ax5.set_xlabel('Variance Ratio')
ax5.set_ylabel('Power')
ax5.set_title('Power by Variance Ratio (Normal, K=10, N=80-200)')
ax5.legend(loc='lower right')
ax5.set_ylim(-0.02, 1.05)
ax5.grid(True, alpha=0.3)
fig5.savefig(os.path.join(OUT_DIR, 'fig5_power_by_vr.png'), dpi=150, bbox_inches='tight')
plt.close('all')

# Fig 6: Power by K
k_df = sim_df[sim_df['Scenario'].str.startswith('Power_K')]
fig6, ax6 = plt.subplots(figsize=(10, 6))
ks = k_df['K'].values
for m, l, c, mk in zip(methods, labels_m, colors, markers):
    ax6.plot(ks, k_df[m].values, '-' + mk, color=c, label=l, linewidth=2, markersize=8)
ax6.axhline(0.80, color='gray', linestyle=':', alpha=0.5)
ax6.set_xlabel('Number of Studies (K)')
ax6.set_ylabel('Power')
ax6.set_title('Power by Number of Studies (Normal, VR=2.0, N=80-200)')
ax6.legend(loc='lower right')
ax6.set_ylim(-0.02, 1.05)
ax6.grid(True, alpha=0.3)
fig6.savefig(os.path.join(OUT_DIR, 'fig6_power_by_k.png'), dpi=150, bbox_inches='tight')
plt.close('all')

# Fig 7: Power by N
n_df = sim_df[sim_df['Scenario'].str.startswith('Power_N')]
fig7, ax7 = plt.subplots(figsize=(10, 6))
ns = [int(r.split('-')[1]) for r in n_df['N_range'].values]
for m, l, c, mk in zip(methods, labels_m, colors, markers):
    ax7.plot(ns, n_df[m].values, '-' + mk, color=c, label=l, linewidth=2, markersize=8)
ax7.axhline(0.80, color='gray', linestyle=':', alpha=0.5)
ax7.set_xlabel('Max Sample Size per Study')
ax7.set_ylabel('Power')
ax7.set_title('Power by Sample Size (Normal, VR=2.0, K=10)')
ax7.legend(loc='lower right')
ax7.set_ylim(-0.02, 1.05)
ax7.grid(True, alpha=0.3)
fig7.savefig(os.path.join(OUT_DIR, 'fig7_power_by_n.png'), dpi=150, bbox_inches='tight')
plt.close('all')

print('All simulation figures saved to output/')
