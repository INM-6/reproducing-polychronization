import argparse
import numpy as np
import os
import sys
import matplotlib
matplotlib.use('Agg')
import json
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mpl_toolkits.axes_grid.inset_locator
import helper as hf
import plot_helper as phf
import seaborn as sns
import scipy.stats as stat
from matplotlib import mlab
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("-sl", '--spikelist', help="list of spike files", nargs="+")
parser.add_argument("-gl", '--groupstatlist', help="list of group stat files", nargs="+")
parser.add_argument("-cl", '--connectivitylist', help="list of connectivity files", nargs="+")

parser.add_argument('--output', type=str)
args = parser.parse_args()
NFFT=512
noverlap=256
bin_ms=1.

#low_gamma
table_low=None
table_high=None
l=0
h=0
N=len(args.spikelist)
exc_Pxx_tab = np.zeros((noverlap+1,N))
reps=[]
gamma_peak=[]
N_grps=[]

c_high = 'C5'
c_low = 'C6'
for rep,(spk_fl,grp_stat_fl,con_fl) in enumerate(zip(args.spikelist,args.groupstatlist,args.connectivitylist)):

    connectivity = pd.read_json(con_fl)
    connecitivty_e = connectivity.loc[connectivity['pre'] < 800]
    connecitivty_e_e = connectivity.loc[connectivity['post'] < 800]

    connecitivty_e_e['bin_w'] = pd.cut(connecitivty_e_e['weight'], np.arange(0, 10.5, 0.5))

    times, senders = hf.read_spikefile(spk_fl)
    exc_times, exc_sender, inh_times, inh_sender = hf.split_in_ex(times, senders)
    exc_rate, exc_bins = hf.bin_pop_rate(exc_times, exc_sender, bin_ms)
    exc_Pxx, exc_freqs = mlab.psd(exc_rate - np.mean(exc_rate), NFFT=NFFT, Fs=1000. /
                                                                              (exc_bins[1] - exc_bins[0]),
                                  noverlap=noverlap)

    exc_Pxx_tab[:, rep] = exc_Pxx
    idx = np.argmax(exc_Pxx[exc_freqs > 20])
    cut_freqs = exc_freqs[exc_freqs > 20]
    max_freq = cut_freqs[idx]
    if max_freq < 50:
        if table_low is None:
            table_low = pd.pivot_table(connecitivty_e_e, columns='bin_w', index='delay', values='weight', aggfunc=len)
        else:
            table_low = table_low.add(
                pd.pivot_table(connecitivty_e_e, columns='bin_w', index='delay', values='weight', aggfunc=len))
        excolor = c_high
        l += 1
        gamma_peak.append('low')

    else:
        if table_high is None:
            table_high = pd.pivot_table(connecitivty_e_e, columns='bin_w', values='weight', index='delay', aggfunc=len)
        else:
            table_high = table_high.add(
                pd.pivot_table(connecitivty_e_e, values='weight', columns='bin_w', index='delay', aggfunc=len))
        excolor = c_low
        h += 1
        gamma_peak.append('high')


    with open(grp_stat_fl, "r") as f:
        stats = json.load(f)
    reps.append(rep)
    N_grps.append(len(stats['N_fired']))
df=pd.DataFrame(dict(reps=reps,n_groups=N_grps,spektral_peak=gamma_peak))


phf.latexify(fig_height=2.5, columns=2)
fig = plt.figure()

gs0 = gridspec.GridSpec(1, 3)

ax_psd = plt.subplot(gs0[0, 0])

ax_weights = plt.subplot(gs0[0, 1])
ax_groups = plt.subplot(gs0[0, 2])

print(df.columns)
ax_psd.plot(exc_freqs, exc_Pxx_tab[:, df.loc[df['spektral_peak']=='low','reps']], color=c_low, linewidth=1.0)
ax_psd.plot(exc_freqs, exc_Pxx_tab[:, df.loc[df['spektral_peak']=='high','reps']], color=c_high, linewidth=1.0)

ax_psd.set_xlim((0, 100))
ax_psd.set_yscale('log')
ax_psd.set_ylim((np.min(exc_Pxx_tab) * 0.5, np.max(exc_Pxx_tab) * 100))
ax_groups = sns.swarmplot(x="spektral_peak", y="n_groups", edgecolor='k', color='k', data=df, zorder=20,
                          size=2.2, order=['low', 'high'])

ax_groups.set_ylim((0, 8000))
sns.boxplot(x="spektral_peak",
            y="n_groups",
            data=df,
            showcaps=True,
            palette=[c_low, c_high],
            order=['low', 'high'],
            showfliers=False,
            ax=ax_groups,
            linewidth=1.25,
            width=0.6
            )
# sns.boxplot(x="gamma peak", y="Number of Groups", data=df_sub,
#                  showcaps=False,boxprops={'facecolor':'None'},
#                  showfliers=False,whiskerprops={'linewidth':0}, ax=ax_groups)

delays = range(1, 21)
ax_weights.step(delays, table_low.values[:, 0] / l, color='C1', linewidth=2., where='mid', alpha=0.9)
ax_weights.step(delays, table_high.values[:, 0] / h, color='C1', linewidth=2., where='mid', alpha=0.5)

ax_weights.step(delays, table_low.values[:, -1] / l, color=c_low, linewidth=2., where='mid')
ax_weights.step(delays, table_high.values[:, -1] / h, color=c_high, linewidth=2., where='mid')

axin1 = mpl_toolkits.axes_grid.inset_locator.inset_axes(ax_psd,
                                                        width=0.65,  # width = 30% of parent_bbox
                                                        height=0.65,  # height : 1 inch
                                                        borderpad=1.5,
                                                        loc=2
                                                        )

# labels = 'low gamma', 'high gamma',
fracs = [int(l), int(h)]
explode = (0.1, 0)

axin1.pie(fracs, explode=explode, colors=[c_low, c_high], labels=fracs,
          # autopct='%1.1f%%',
          shadow=True, startangle=-95)

ax_psd.set_xlabel('Frequency [Hz]')
ax_psd.set_ylabel('Power [a.u.]')

ax_psd.set_yticks((0.01, 1, 100))

ax_weights.set_xlabel('Delay [ms]')
ax_weights.set_ylabel('Frequency')
# ax_weights.set_yticks(())
ax_weights.set_xticks((1, 10, 20))

ax_groups.set_xlabel('Spectral peak')
ax_groups.set_ylabel(r'Number of groups [$10^3$]')
ax_groups.set_yticks((0, 2500, 5000, 7500))
ax_groups.set_yticklabels((0, 2.5, 5.0, 7.500))
for ax, letter in [(ax_psd, 'A'), (ax_weights, 'B'), (ax_groups, 'C')]:
    ax.annotate(r'\textbf{{{letter}}}'.format(letter=letter), xy=(-0.3, 0.99), xycoords='axes fraction', fontsize=10,
                horizontalalignment='left', verticalalignment='top', annotation_clip=False)

gs0.update(left=0.1, right=0.99, top=0.95, bottom=0.2, hspace=0.2, wspace=0.35)


plt.savefig(args.output)
