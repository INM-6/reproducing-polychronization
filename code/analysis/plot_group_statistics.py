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


# sns.set_palette(sns.light_palette("blue"))
parser = argparse.ArgumentParser()
parser.add_argument("-glo", '--groupstatlist_original', help="list of group stat files", nargs="+")
parser.add_argument("-glp", '--groupstatlist_python', help="list of group stat files", nargs="+")

parser.add_argument('--output', type=str)
args = parser.parse_args()


def return_val(stats):
    if 'Failed' in stats.keys():
        ngr='Failed'
    elif 'Failed' not in stats.keys():

        ngr = len(stats['N_fired'])
    else:
        ngr = np.nan
    return ngr

exp,reps=[],[]
ngr_o=[]
ngr_p=[]
rate_exc=[]
rate_inh=[]
spek_peak=[]
experiments=[i.split('/')[2] for i in args.groupstatlist_original]+[i.split('/')[2] for i in args.groupstatlist_python]

for experiment in np.unique(experiments):
    repetitions = [i.split('/')[3] for i in args.groupstatlist_original  if i.split('/')[2]==experiment] + \
                  [i.split('/')[3] for i in args.groupstatlist_python    if i.split('/')[2]==experiment]
    repetitions=np.sort(np.unique(repetitions))
    for repetition in repetitions:


        print(repetition,experiment)
        file_p='data/NEST_model/{e}/{r}/stats.json'.format(e=experiment,r=repetition)
        if os.path.isfile(file_p):
            with open(file_p, "r") as f_p:
                stats_p = json.load(f_p)
            ngr_p_ = return_val(stats_p)

        else:
            ngr_p_=np.nan


        #print(stats_p)

        file_o='data/NEST_model/{e}/{r}/stats_orig.json'.format(e=experiment,r=repetition )
        #print(file)
        if os.path.isfile(file_o):
            with open(file_o, "r") as f_o:
                stats_o = json.load(f_o)
            ngr_o_ = return_val(stats_o)
        else:
            ngr_o_ = np.nan

        spk_fl = 'data/NEST_model/{e}/{r}/spikes-1001.gdf'.format(e=experiment,r=repetition )
        data = np.loadtxt(spk_fl)
        senders, times = data[:, 0], data[:, 1]

        mean_ex, mean_inh, max_freq = phf.get_rates(times, senders)

        if max_freq<50:
            spek_peak.append('low')
        else:
            spek_peak.append('high')

        ngr_o.append(ngr_o_)
        ngr_p.append(ngr_p_)
        exp.append(experiment.replace('_',' '))
        reps.append(repetition)
        rate_exc.append(mean_ex)
        rate_inh.append(mean_inh)

df=pd.DataFrame({'Number of groups':ngr_o,
                'Number of groups (nest)':ngr_p,
                 'Experiment':exp,
                  'reps':reps,
                 'exc_rate':rate_exc,
                 'inh_rate': rate_inh,
                 'spektral peak':spek_peak
                 })
def iqr(df):
    return df.quantile(.75)-df.quantile(.25)

df_latex=df.replace(value=np.nan,to_replace='Failed').groupby(['Experiment'])['Number of groups', 'Number of groups (nest)','exc_rate','inh_rate','spektral peak'].agg([np.median,iqr,'min','max','count']) #.agg([np.median,iqr])
print(df_latex.to_latex())
print(df_latex)
df_latex_spek=df.replace(value=np.nan,to_replace='Failed').groupby(['Experiment','spektral peak'])['Number of groups', 'Number of groups (nest)','exc_rate','inh_rate','spektral peak'].agg([np.median,iqr,'min','max','count']) #.agg([np.median,iqr])
print(df_latex_spek.to_latex())
print(df_latex_spek)

phf.latexify(fig_height=6., columns=1)
fig = plt.figure()
N = 9

N_bot = 5
M = 4
gs0 = gridspec.GridSpec(N, M)

ax_orig = plt.subplot(gs0[:N_bot, :M - 1])
ax_nest = plt.subplot(gs0[N_bot:, 0:M - 1])

ax_orig_broken = fig.add_subplot(gs0[:N_bot, M - 1])  # , sharey=ax_orig)
ax_nest_broken = fig.add_subplot(gs0[N_bot:, M - 1])  # , sharey=ax_nest)

orig_pal = ['C2', 'C1', 'C0', 'C5', 'C4', 'C4', 'C4', 'C4', 'C4', 'C4', 'C4']
orig_exp_order = ['initial reproduction',
                  'bitwise reproduction',
                  'qualitative model',
                  'poisson stimulus',
                  'stdp window match',
                  'const add value 0p0',
                  'synapse update interval 0p1s',
                  'synapse update interval 10s',
                  'time driven additive 1s',
                  'tau syn update interval 2s',
                  'tau syn update interval 1000s']

orig_names = ['Initial model',
              'Bitwise reproduction',
              'Qualitative model',
              'Poisson stimulus',
              'STDP window match',
              'No additive factor',
              'Buffer length $0.1\;\mathrm{s}$',
              'Buffer length $10\;\mathrm{s}$',
              'No elig. trace',
              'Elig. trace $2\;\mathrm{s}$',
              'Elig. trace $1000\;\mathrm{s}$',
              ]
width = 1.25
ax_orig = sns.boxplot(data=df.replace(value=np.nan,to_replace='Failed'),
                      y='Experiment',
                      x='Number of groups',
                      order=orig_exp_order,
                      palette=orig_pal,
                      fliersize=0,
                      ax=ax_orig,
                      linewidth=width,
                      width=0.6)

ax_orig_broken = sns.boxplot(data=df.replace(value=np.nan,to_replace='Failed'),
                             y='Experiment',
                             x='Number of groups',
                             order=orig_exp_order,
                             palette=orig_pal,
                             fliersize=0,
                             ax=ax_orig_broken,
                             linewidth=width,
                             width=0.6)

ax_orig.set_yticklabels(orig_names)

nest_pal = ['C1', 'C0', 'C3', 'C3', 'C3', 'C3','C4', 'C4']
nest_exp_order = ['bitwise reproduction',
                  'qualitative model',
                  'delay distribution 20',
                  'delay distribution 15',
                  'delay distribution 10',
                  'delay distribution 5',
                  'resolution 0p1 W pspmatched',
                  'qualitative model high res',
                  ]

name_order = ['Bitwise reproduction',
              'Qualitative model',
              r'Delay $\in \left[1,20\right]\;\mathrm{ms}$',
              r'Delay $\in \left[1,15\right]\;\mathrm{ms}$',
              r'Delay $\in \left[1,10\right]\;\mathrm{ms}$',
              r'Delay $\in \left[1,5\right]\;\mathrm{ms}$',
              r'Resolution $0.1\;\mathrm{ms}$',
              r'Improved integration',
              ]
ax_nest = sns.boxplot(data=df.replace(value=np.nan,to_replace='Failed'),
                      y='Experiment',
                      x='Number of groups (nest)',
                      order=nest_exp_order,
                      palette=nest_pal,
                      fliersize=0,
                      ax=ax_nest, linewidth=width,
                      width=0.5)
ax_nest_broken = sns.boxplot(data=df.replace(value=np.nan,to_replace='Failed'),
                             y='Experiment',
                             x='Number of groups (nest)',
                             order=nest_exp_order,
                             palette=nest_pal,
                             fliersize=0,
                             ax=ax_nest_broken,
                             linewidth=width,
                             width=0.6)
print(ax_nest.get_yticks())
ax_nest.set_yticklabels(name_order)

for ax in [ax_orig_broken, ax_nest_broken]:
    # ax.axis('off')
    # if ax !=ax_delay:
    # ax.set_xscale('log')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks(())
    ax.set_ylabel('')

ax_orig_broken.set_xlabel('')
ax_nest_broken.set_xlabel('')

ax_orig_broken.set_xlim((6000, 70000))
ax_orig_broken.set_xticks((10000, 70000))
ax_orig_broken.set_xticklabels(('10k', '70k'))

ax_nest_broken.set_xlim((10000, 40000))
ax_nest_broken.set_xticks((20000, 40000))
ax_nest_broken.set_xticklabels(('20k', '40k'))

ax_orig.set_xlim((-500, 6000))
ax_orig.set_xticks((0, 2500, 5000))
ax_nest.set_xlim((-500, 8500))

for ax in [ax_orig, ax_nest]:
    # ax.axis('off')
    # if ax !=ax_delay:
    # ax.set_xscale('log')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.set_yticks(())

#     if ax !=ax_delay:
#         ax.set_xticks(())
#         ax.set_xlabel('')
#         ax.spines['bottom'].set_visible(False)


# ax_orig.set_ylabel('original group finding algorithm')
# ax_nest.set_ylabel('NEST group finding algorithm')

ax_orig.set_ylabel('')
ax_nest.set_ylabel('')

ax_orig.set_xlabel('Number of groups')

ax_nest.set_xlabel('Number of groups (Python)')

# iqr=df.loc[df['Experiment']=='qualitative model','Number of Groups'].quantile(0.75)-df.loc[df['Experiment']=='qualitative model','Number of Groups'].quantile(0.25)
min_line = df.loc[df['Experiment'] == 'bitwise reproduction', 'Number of groups'].quantile(0.25)  # -1.5*iqr
max_line = df.loc[df['Experiment'] == 'bitwise reproduction', 'Number of groups'].quantile(0.75)  # +1.5*iqr
min_line
ax_orig.axvline(min_line, zorder=0, linestyle='--', color='C1')
ax_orig.axvline(max_line, zorder=0, linestyle='--', color='C1')
min_line = df.loc[df['Experiment'] == 'bitwise reproduction', 'Number of groups (nest)'].quantile(0.25)
max_line = df.loc[df['Experiment'] == 'bitwise reproduction', 'Number of groups (nest)'].quantile(0.75)
min_line
ax_nest.axvline(min_line, zorder=0, linestyle='--', color='C1')
ax_nest.axvline(max_line, zorder=0, linestyle='--', color='C1')

ax_orig.annotate(r'\textbf{A}', xy=(-0.95, 1.05), xycoords='axes fraction',
                 horizontalalignment='left', verticalalignment='top', annotation_clip=False)
ax_nest.annotate(r'\textbf{B}', xy=(-0.95, 1.06), xycoords='axes fraction',
                 horizontalalignment='left', verticalalignment='top', annotation_clip=False)

xy = (0, ax_nest.get_yticks()[-1])

#ax_nest.annotate(xy=xy, xytext=xy, s=r'\textbf{X}', ha='center', va='center')

# xy=(ax_orig_broken.get_xticks()[-1],ax_orig.get_yticks()[-1])
# ax_orig_broken.annotate(xy=xy,xytext=xy,s=r'$\rip$',ha='center',va='center',fontsize=20)

# xy=(ax_orig_broken.get_xticks()[-1],ax_orig.get_yticks()[-1])
# ax_orig_broken.annotate(xy=xy,xytext=xy,s=r'$\rip$',ha='center',va='center',fontsize=20)

xy = (ax_orig_broken.get_xticks()[-1], ax_orig.get_yticks()[-4])
ax_orig_broken.annotate(xy=xy, xytext=xy, s=r'$\rip$', ha='center', va='center', fontsize=20)

# xy=(ax_orig_broken.get_xticks()[-1],ax_orig.get_yticks()[-4])
# ax_orig_broken.annotate(xy=xy,xytext=xy,s=r'$\rip$',ha='center',va='center',fontsize=20)

gs0.update(left=0.4, right=0.95, top=0.97, bottom=0.07, hspace=1.99, wspace=0.35)


plt.savefig(args.output)