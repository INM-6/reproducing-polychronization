import numpy as np
import json
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import pylab as plt
import helper as hf
import plot_helper as phf
import argparse
import matplotlib.patches as mpatches
import mpl_toolkits.axes_grid.inset_locator
import seaborn as sns
import plot_helper as phf
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--spikefile', type=str)
parser.add_argument('-w', '--weightfile', type=str)
parser.add_argument('-fn', '--filename', type=str)


args = parser.parse_args()


phf.latexify(columns=2)
excolor='C0'
incolor='C1'


spikefile = args.spikefile
times, senders = hf.read_spikefile(spikefile)


weights = hf.read_weightfile(args.weightfile)


fig = plt.figure()
gs0 = gridspec.GridSpec(2, 2)
gs0.update(left=0.1, right=0.97, top=0.97, bottom=0.06, hspace=0.15)

gs1 = gridspec.GridSpecFromSubplotSpec(5, 1, subplot_spec=gs0[0, :])

ax01 = plt.subplot(gs1[:4, 0])
ax02 = plt.subplot(gs1[4, 0])

idxes_subsample=senders %4==0
idxes_times=times>np.max(times)-5000
senders=senders[idxes_subsample&idxes_times]
times=times[idxes_subsample&idxes_times]

phf.plot_raster_rate(times, senders, ax01, ax02, incolor=incolor, excolor=excolor,bin_ms=5.,linewidth=.75)

ax1 = plt.subplot(gs0[1, 0])
ax2 = plt.subplot(gs0[1, 1])


phf.plot_weights(weights, ax1, 'k', ylim=[150, 55000], alpha=0.5)
axin1 = mpl_toolkits.axes_grid.inset_locator.inset_axes(ax1,
                                                        width="35   %",  # width = 30% of parent_bbox
                                                        height=.75,  # height : 1 inch
                                                        borderpad=2.0

                                                        )

boxplot_kwargs = dict(positions=range(4),
                      bootstrap=1000,
                      showmeans=True,
                      labels=['Inh', 'Exc']
                      )

idxes_subsample=senders %4==0
idxes_times=times>np.max(times)-5000
senders=senders[idxes_subsample&idxes_times]
times=times[idxes_subsample&idxes_times]

exc_times, exc_sender, inh_times, inh_sender = hf.split_in_ex(times, senders)
inh_rate, inh_bins = hf.bin_pop_rate(inh_times, inh_sender, 1.)
exc_rate, exc_bins = hf.bin_pop_rate(exc_times, exc_sender, 1.)


def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)


bpexc = axin1.boxplot([exc_rate],
                      sym='',
                      widths=0.6,
                      labels=['exc'])
bpinh = axin1.boxplot([inh_rate],
                      sym='',
                      widths=0.6,
                      labels=['inh'])
set_box_color(bpexc, excolor)  # colors are from http://colorbrewer2.org/
set_box_color(bpinh, incolor)
axin1.set_xlim([-1, 3])
axin1.set_yticks([0, 25, 50, 75])
axin1.set_ylabel('rate distribution')


phf.plot_psd(times, senders, ax2, incolor=incolor, excolor=excolor)

# exc = plt.Line2D((0, 1), (0, 0), color='k', linestyle='-')
# inh = plt.Line2D((0, 1), (0, 0), color='b', linestyle='-')
# statistical = plt.Line2D((0, 1), (0, 0), color='k',  linestyle='--')
# bitwise = plt.Line2D((0, 1), (0, 0), color='k', linestyle='-')
#
# # Create legend from custom artist/label lists
# ax2.legend([exc,inh,statistical,bitwise],
#             ['Exc','Inh','Naive', 'Bitwise'], loc=1,
#            prop={'size': 12})
gs0.update(left=0.1, right=0.97, top=0.97, bottom=0.1, hspace=0.25)

for ax, letter in [(ax1, 'B'), (ax2, 'C')]:
    ax.annotate(r'\textbf{{{letter}}}'.format(letter=letter), xy=(-0.08, 0.99), xycoords='axes fraction', fontsize=10,
                horizontalalignment='left', verticalalignment='top', annotation_clip=False)
ax01.annotate(r'\textbf{A}', xy=(-0.03, 0.99), xycoords='axes fraction', fontsize=10,
              horizontalalignment='left', verticalalignment='top', annotation_clip=False)

plt.savefig(args.filename)
