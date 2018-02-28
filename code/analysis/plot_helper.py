import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import pylab as plt
import helper as hf
import numpy as np
import matplotlib.mlab as mlab
import json
import os


def mem_spk_plot(data, times, sender, subplotspec, mem_color, spk_inh_color, spk_exc_color):
    id = data[:, 0]
    t = data[:, 1]
    v = data[:, 2]
    t_max = np.max(t)
    unique_sender = np.unique(id)

    gs1 = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=subplotspec,  hspace=0.05)

    in_neurons = unique_sender[unique_sender > 800]
    ex_neurons = unique_sender[unique_sender <= 800]
    if in_neurons.size > 0 and ex_neurons.size > 0:
        in_idx = in_neurons[0]
        ex_idx = ex_neurons[0]

        ax0 = plt.subplot(gs1[0, 0])
        ax1 = plt.subplot(gs1[1, 0])

        ax0.plot(t[id == in_idx], v[id == in_idx], mem_color, linewidth=0.8)
        ax0.plot(times[sender == in_idx - 1], sender[sender == in_idx - 1] * 0,
                 color=spk_inh_color, marker='*', markersize=10, linestyle=' ')
        ax1.plot(t[id == ex_idx], v[id == ex_idx], mem_color, linewidth=0.8)
        ax1.plot(times[sender == ex_idx - 1], sender[sender == ex_idx - 1] * 0,
                 color=spk_exc_color, marker='*', markersize=10, linestyle=' ')
        for ax in [ax0, ax1]:
            ax.set_xlim([t_max - 1000, t_max - 000])
            ax.set_ylim([-100, 100])
            ax.set_yticklabels([-50, 0, 50])
            ax.set_yticks([-50, 0, 50])
            xticks = np.linspace(0, 1000, num=5,
                                 endpoint=True)
        ax1.set_xticks(xticks + t_max - 1000)

        ax1.set_xticklabels(xticks)

        ax0.set_xticks([])
        ax1.set_xlabel('Time [s]')
        ax1.set_ylabel('Membrane Potential [mv]', y=1.1)
    return ax0, ax1


def plot_raster_rate(times, senders, ax01, ax02, incolor='b', excolor='k'):
    exc_times, exc_sender, inh_times, inh_sender = hf.split_in_ex(
        times, senders)

    inh_rate, inh_bins = hf.bin_pop_rate(inh_times, inh_sender, 1.)
    exc_rate, exc_bins = hf.bin_pop_rate(exc_times, exc_sender, 1.)

    ax01.plot(exc_times, exc_sender, excolor,
              marker='.', linestyle='', markersize=1)
    ax01.plot(inh_times, inh_sender, incolor,
              marker='.', linestyle='', markersize=1)

    ax01.set_xlim([np.min(times), np.max(times)])
    ax01.set_ylim([0, 1000])
    ax01.set_xticks([])
    ax01.set_yticks([250, 500, 750])

    # ax01.set_xticklabels(np.linspace(0, 5, num=6,endpoint=True))
    # ax01.set_xlabel('Time [s]')
    ax01.set_ylabel('Neuron Id')

    ax02.plot(inh_bins - np.min(times), inh_rate, incolor)
    ax02.plot(exc_bins - np.min(times), exc_rate, excolor)
    ax02.set_xlabel('Time [s]')
    ax02.set_ylabel(r'$f_{pop}$ [spk/s]')
    ax02.set_ylim([0, 100])

    ax02.set_yticks([0, 50, 100])
    ax02.set_yticklabels([0, 50, 100])

    ax02.set_xlim([np.min(times) - np.min(times),
                   np.max(times) - np.min(times)])
    xticks = np.linspace(np.min(times) - np.min(times), np.max(times) - np.min(times) + 1, num=5,
                         endpoint=True)
    ax02.set_xticks(xticks)
    ax02.set_xticklabels([0, 2.5, 5, 7.5, 10])

    ax01.set_yticks([250, 500, 750])


def plot_weights(weights, ax, c='b', bins=40, normed=False, xlim=[0., 10.], ylim=[150., 50000], scale='log', linestyle='-', alpha=0.5):
    print(np.min(weights), np.max(weights))
    if np.max(weights) > 10:
        xlim = [0, 100]
    ax.hist(weights, bins=np.arange(xlim[0], xlim[1] + xlim[1] * 1. / bins, xlim[1]
                                    * 1. / bins), normed=normed, color=c, linestyle=linestyle, alpha=alpha)
    # ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    ax.set_xlabel('Synaptic weight [mV]')
    ax.set_ylabel('Frequency')
    # ax.set_yscale(scale)


def plot_psd(times, senders, ax, NFFT=512, noverlap=256, xlim=[0., 150.], ylim=[1e-3, 1e2], scale='log', incolor='C0', excolor='C1', linewidth=2):

    exc_times, exc_sender, inh_times, inh_sender = hf.split_in_ex(
        times, senders)
    if incolor is not None:
        inh_rate, inh_bins = hf.bin_pop_rate(inh_times, inh_sender, 1.)
        inh_Pxx, inh_freqs = mlab.psd(inh_rate - np.mean(inh_rate), NFFT=NFFT, Fs=1000. /
                                      (inh_bins[1] - inh_bins[0]), noverlap=noverlap)
        ax.plot(inh_freqs, inh_Pxx, color=incolor, linewidth=linewidth)
    if excolor is not None:
        exc_rate, exc_bins = hf.bin_pop_rate(exc_times, exc_sender, 1.)
        exc_Pxx, exc_freqs = mlab.psd(exc_rate - np.mean(exc_rate), NFFT=NFFT, Fs=1000. /
                                      (exc_bins[1] - exc_bins[0]), noverlap=noverlap)
        ax.plot(exc_freqs, exc_Pxx, color=excolor, linewidth=linewidth)

    # ax3.plot(bins,freqs,Pxx)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Power [a.u.]')
    ax.set_yscale(scale)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_yticks([])


def plot_2D_weights(weights, delays, counts, ax, range=None, cmap='gray_r'):
    if range is None:
        range = [[-0.5, 10.5], [0.5, 20.5]]
    #sns.jointplot(x=weights, y=delays, data=np.array([weights,delays]), kind="kde");
    ax.pcolor(weights, delays, counts, cmap=cmap)
    ax.set_xlim(range[0])
    ax.set_ylim(range[1])
    ax.set_ylabel('Delay [ms]')
    ax.set_xlabel('Weight [mV]')


def plot_group(group, ax, LP=False, numbers=True):
    N_fired = group['N_fired']
    times, senders = hf.get_t_s(group)
    idx_in = senders > 800
    idx_ex = senders < 800
    ax1 = ax.twinx()

    ax.scatter(times[idx_ex], senders[idx_ex], s=80,
               facecolors='none', edgecolors='r')
    ax1.scatter(times[idx_in], (senders[idx_in] - 800)
                * 4, s=80, facecolors='k')

    for link in group['links']:
        if link['post'] < 799:
            post_i = link['post']
        else:
            post_i = (link['post'] - 800) * 4
        pre_times = times[link['pre'] == senders]
        if pre_times.size > 1:
            pre_times = [pre_times[0]]
        post_times = times[link['post'] == senders]
        if post_times.size > 1:
            post_times = [post_times[0]]

        ax.plot([pre_times, post_times],
                [link['pre'], post_i], 'k')
        if numbers:
            ax.text(times[link['pre'] == senders] + 1,
                    link['pre'] - 10,
                    s=str(link['pre']),
                    fontsize=12)

            ax.text(times[link['post'] == senders] + 1,
                    post_i - 10,
                    s=str(link['post']),
                    fontsize=12)

    if LP:

        i = len(group['links']) - 1
        while group['links'][i]['layer'] >= 2:
            i_l = i
            if group['links'][i]['post'] < 799:
                post_i = group['links'][i]['post']
            else:
                post_i = (group['links'][i]['post'] - 800) * 4
            ax.plot([times[group['links'][i]['pre'] == senders],
                     times[group['links'][i]['post'] == senders]], [group['links'][i]['pre'],
                                                                    post_i], 'k', linewidth=4, alpha=0.5)
            if group['links'][i]['layer'] == group['L_max']:
                ax.text(times[group['links'][i]['post'] == senders] - 20,
                        post_i + 20,
                        s='longest path',
                        fontsize=11)
            for j in range(len(group['links'])):
                if group['links'][j]['post'] == group['links'][i]['pre']:
                    i_l = j
            if group['links'][i]['layer'] == 2:
                break
            i = i_l

    ax.set_ylim(-20, 820)
    ax1.set_ylim(-20, 820)

    ax.set_yticks([0, 400, 800])

    ax1.set_yticks([0, 800])
    ax1.set_yticklabels([800, 1000])

    ax.set_xlim(np.min(times) - 10, np.max(times) + 10)
    ax1.set_xlim(np.min(times) - 10, np.max(times) + 10)


def return_NTL(groupfile):
    if os.path.getsize(groupfile) < 422099208 * 1.1:
        groups = hf.read_group_file(groupfile)
    else:
        f = open(groupfile)
        groups = ijson.items(f, "item")
    N_list = []
    L_list = []
    T_list = []
    i = 0
    for i, g in enumerate(groups):
        times, senders = hf.get_t_s(g)

        N_list.append(int(g["N_fired"]))

        T_list.append(max(times))  # time span

        L_list.append(int(g["L_max"]))  # longest path
    return N_list, T_list, L_list
# def return_NTL(group_data):
#     N_list = []
#     L_list = []
#     T_list = []
#     for i, g in enumerate(group_data):
#         times, senders = hf.get_t_s(g)
#
#         N_list.append(int(g["N_fired"]))
#
#         T_list.append(max(times))  # time span
#
#         L_list.append(int(g["L_max"]))  # longest path
#     return N_list,T_list,L_list


def plot_8(group_data, filename, outname):
    fig = plt.figure("plot 8")
    ax0 = fig.add_subplot(221)
    ax1 = fig.add_subplot(222)
    ax2 = fig.add_subplot(223)
    ax3 = fig.add_subplot(224)
    N_list, T_list, L_list = return_NTL(filename)

    ax1.hist(N_list, np.arange(0, 100, 5))
    ax1.set_title('# of neurons, total {}'.format(len(group_data)))
    ax2.hist(T_list,  np.arange(0, 200, 10))
    ax2.set_title('time span[ms]')
    ax3.hist(L_list, np.arange(0, 30, 2))
    ax3.set_title('length of longest path')

    if len(group_data) > 1:
        plot_group(group_data[0], ax0, LP=False, numbers=False)

    plt.savefig(outname)
    return


def plot_combined_groups_statstics(group_files_list, outname):
    fig = plt.figure("plot 8")
    ax0 = fig.add_subplot(221)
    ax1 = fig.add_subplot(222)
    ax2 = fig.add_subplot(223)
    ax3 = fig.add_subplot(224)
    N_combined = []
    L_combined = []
    T_combined = []
    for group_data in group_files_list:
        N_list, T_list, L_list = return_NTL(group_data)
        N_combined += N_list
        T_combined += T_list
        L_combined += L_list

    fig.suptitle('# groups in {} iterations: {}'.format(
        len(group_files_list), len(N_combined)))
    ax1.hist(N_combined, np.arange(0, 100, 5))
    ax1.set_xlabel('# of neurons')
    ax2.hist(T_combined, np.arange(0, 200, 10))
    ax2.set_xlabel('Time span[ms]')
    ax3.hist(L_combined, np.arange(0, 30, 2))
    ax3.set_xlabel('Longest path')
    ax0.set_xticks([])
    boxplot_kwargs = dict(positions=range(3),
                          bootstrap=1000,
                          showmeans=True,
                          labels=['#neurons', 'timespan', 'longest path']
                          )

    def set_box_color(bp, color):
        plt.setp(bp['boxes'], color=color)
        plt.setp(bp['whiskers'], color=color)
        plt.setp(bp['caps'], color=color)
        plt.setp(bp['medians'], color=color)

    bpexc = ax0.boxplot([N_combined, T_combined, L_combined],
                        **boxplot_kwargs)

    set_box_color(bpexc, 'b')  # colors are from http://colorbrewer2.org/

    plt.savefig(outname)
