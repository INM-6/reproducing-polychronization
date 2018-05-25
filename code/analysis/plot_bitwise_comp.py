import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import pylab as plt
import helper as hf
import plot_helper as phf
import argparse
import mpl_toolkits.axes_grid.inset_locator



parser = argparse.ArgumentParser()
parser.add_argument('-cs', '--comp_spikefile', type=str)
parser.add_argument('-bs', '--bitwise_spikefile', type=str)

parser.add_argument('-cw', '--comp_weightfile', type=str)
parser.add_argument('-bw', '--bitwise_weightfile', type=str)

parser.add_argument('-fn', '--filename', type=str)


args = parser.parse_args()

comp_spikefile = args.comp_spikefile
comp_times, comp_senders = hf.read_spikefile(comp_spikefile)
bitwise_spikefile = args.bitwise_spikefile
bitwise_times, bitwise_senders = hf.read_spikefile(bitwise_spikefile)

comp_weights = hf.read_weightfile(args.comp_weightfile)
bitwise_weights = hf.read_weightfile(args.bitwise_weightfile)


bin_ms=5.
phf.latexify(columns=2)
excolor='C0'
incolor='C1'


fig = plt.figure()
gs0 = gridspec.GridSpec(2, 2)
gs0.update(left=0.1, right=0.97, top=0.97, bottom=0.1, hspace=0.25)

gs1 = gridspec.GridSpecFromSubplotSpec(7, 1, subplot_spec=gs0[0, :])

ax01 = plt.subplot(gs1[:5, 0])
ax02 = plt.subplot(gs1[5:, 0])
ax1 = plt.subplot(gs0[1, 0])
ax2 = plt.subplot(gs0[1, 1])

#only plot every 10th sender
idxes_subsample=comp_senders %4==0
idxes_times=comp_times>np.max(comp_times)-5000
senders=comp_senders[idxes_subsample&idxes_times]
times=comp_times[idxes_subsample&idxes_times]

phf.plot_raster_rate(times,senders, ax01, ax02, incolor=incolor, excolor=excolor,bin_ms=bin_ms,linewidth=.75)

phf.plot_weights(comp_weights, ax1, 'k', ylim=[150, 55000], alpha=0.5)
phf.plot_weights(bitwise_weights, ax1, 'k', ylim=[150, 55000], alpha=0.3)
axin1 = mpl_toolkits.axes_grid.inset_locator.inset_axes(ax1,
                                                        width="35   %",  # width = 30% of parent_bbox
                                                        height=.75,  # height : 1 inch
                                                        borderpad=2.0

                                                        )

boxplot_kwargs = dict(positions=range(4),
                      bootstrap=1000,
                      showmeans=True,
                      labels=['Inh', 'Exc', 'Inh', 'Exc']
                      )
comp_exc_times, comp_exc_sender, comp_inh_times, comp_inh_sender = hf.split_in_ex(comp_times, comp_senders)
comp_inh_rate, comp_inh_bins = hf.bin_pop_rate(comp_inh_times, comp_inh_sender, bin_ms)
comp_exc_rate, comp_exc_bins = hf.bin_pop_rate(comp_exc_times, comp_exc_sender, bin_ms)
bitwise_exc_times, bitwise_exc_sender, bitwise_inh_times, bitwise_inh_sender = hf.split_in_ex(bitwise_times, bitwise_senders)
bitwise_inh_rate, bitwise_inh_bins = hf.bin_pop_rate(bitwise_inh_times, bitwise_inh_sender, bin_ms)
bitwise_exc_rate, bitwise_exc_bins = hf.bin_pop_rate(bitwise_exc_times, bitwise_exc_sender, bin_ms)





bpexc = axin1.boxplot([bitwise_exc_rate, comp_exc_rate],
                      positions=np.array(range(2)) * 2.0 - 0.4,
                      sym='',
                      widths=0.6,
                      labels=['original', 'NEST'])
bpinh = axin1.boxplot([bitwise_inh_rate, comp_inh_rate],
                      positions=np.array(range(2)) * 2.0 + 0.4,
                      sym='',
                      widths=0.6,
                      labels=['original', 'NEST'])
phf.set_box_color(bpexc, excolor)  # colors are from http://colorbrewer2.org/
phf.set_box_color(bpinh, incolor)
axin1.set_xlim([-1, 3])
axin1.set_yticks([0, 25, 50, 75])
axin1.set_ylabel('rate distribution')

phf.plot_psd(comp_times, comp_senders, ax2, incolor=None, excolor=excolor)
phf.plot_psd(bitwise_times, bitwise_senders, ax2, incolor=None, excolor='C2')

# exc = plt.Line2D((0, 1), (0, 0), color='k', linestyle='-')
# inh = plt.Line2D((0, 1), (0, 0), color='b', linestyle='-')
# comp = plt.Line2D((0, 1), (0, 0), color='k',  linestyle='--')
# bitwise = plt.Line2D((0, 1), (0, 0), color='k', linestyle='-')
#
# # Create legend from custom artist/label lists
# ax2.legend([exc,inh,comp,bitwise],
#             ['Exc','Inh','Naive', 'Bitwise'], loc=1,
#            prop={'size': 12})
gs0.update(left=0.1, right=0.97, top=0.97, bottom=0.1, hspace=0.25)

for ax, letter in [(ax1, 'B'), (ax2, 'C')]:
    ax.annotate(r'\textbf{{{letter}}}'.format(letter=letter), xy=(-0.08, 0.99), xycoords='axes fraction', fontsize=10,
                horizontalalignment='left', verticalalignment='top', annotation_clip=False)
ax01.annotate(r'\textbf{A}', xy=(-0.03, 0.99), xycoords='axes fraction', fontsize=10,
              horizontalalignment='left', verticalalignment='top', annotation_clip=False)

plt.savefig(args.filename)
