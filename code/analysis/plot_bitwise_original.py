import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import pylab as plt
import helper as hf
import plot_helper as phf
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-bs', '--bitwise_spikefile', type=str)
parser.add_argument('-os', '--original_spikefile', type=str)
parser.add_argument('-bmem', '--bitwise_mem_pop_file', type=str)


parser.add_argument('-fn', '--filename', type=str)


args = parser.parse_args()

original_spikefile = args.original_spikefile
original_times, original_senders = hf.read_spikefile(original_spikefile)
bitwise_spikefile = args.bitwise_spikefile
bitwise_times, bitwise_senders = hf.read_spikefile(bitwise_spikefile)

bitwise_mem_pop = np.loadtxt(args.bitwise_mem_pop_file)



phf.latexify(columns=2)
excolor='C0'
incolor='C1'

fig = plt.figure()
gs0 = gridspec.GridSpec(2, 2)
gs0.update(left=0.1, right=0.97, top=0.97, bottom=0.1, hspace=0.25)

gs1 = gridspec.GridSpecFromSubplotSpec(7, 1, subplot_spec=gs0[0, :])

ax01 = plt.subplot(gs1[:5, 0])
ax02 = plt.subplot(gs1[5:, 0])
#only plot every 10th sender
idxes_subsample=bitwise_senders %4==0
idxes_times=bitwise_times>np.max(bitwise_times)-5000
senders=bitwise_senders[idxes_subsample&idxes_times]
times=bitwise_times[idxes_subsample&idxes_times]

phf.plot_raster_rate(times, senders, ax01, ax02, incolor=incolor, excolor=excolor,bin_ms=5.,linewidth=.75)
ax2 = plt.subplot(gs0[1, 1])


ax0, ax1 = phf.mem_spk_plot(bitwise_mem_pop, original_times, original_senders,
                            gs0[1, 0], mem_color='k', spk_exc_color=excolor, spk_inh_color=incolor)

phf.plot_psd(original_times, original_senders, ax2, excolor='C2', incolor=None,linewidth=3)
phf.plot_psd(bitwise_times[bitwise_times>=np.min(original_times)], bitwise_senders[bitwise_times>=np.min(original_times)], ax2, excolor=excolor, incolor=None)

# NEST = plt.Line2D((0, 1), (0, 0), color=excolor, linestyle='-')
# original = plt.Line2D((0, 1), (0, 0), color='C2', linestyle='-')
#
# # Create legend from custom artist/label lists
# ax2.legend([original,NEST],
#             ['NEST','original'], loc=1,
#            prop={'size': 12})

for ax, letter in [(ax0, 'B'), (ax1, 'C'), (ax2, 'D')]:
    ax.annotate(r'\textbf{{{letter}}}'.format(letter=letter), xy=(-0.08, 0.99), xycoords='axes fraction', fontsize=10,
                horizontalalignment='left', verticalalignment='top', annotation_clip=False)
ax01.annotate(r'\textbf{A}', xy=(-0.03, 0.99), xycoords='axes fraction', fontsize=10,
              horizontalalignment='left', verticalalignment='top', annotation_clip=False)


plt.savefig(args.filename)
