import numpy as np
import matplotlib.pyplot as plt
import json
import yaml
import os


def conn_dist(data, delay_or_weight):
    delay = np.array([i[delay_or_weight] for i in data])
    pre = np.array([i['pre'] for i in data])
    post = np.array([i['post'] for i in data])

    idx_pre_ex = pre < 800

    idx_post_ex = post < 800
    idx_post_in = post >= 800
    ex_ex = delay[idx_pre_ex & idx_post_ex]
    ex_in = delay[idx_pre_ex & idx_post_in]

    all = delay[idx_pre_ex]

    return all, ex_ex, ex_in


def parse_config(config_file):
    """
    Parses config file for usage in model.py and find_groups_nest
    """
    with open(config_file, 'r') as ymlfile:
        cfg = yaml.load(ymlfile)
    base = os.path.basename(config_file)
    base_and_extless = os.path.splitext(base)[0]
    cfg['simulation-params']['data-prefix'] = base_and_extless
    return cfg



def bin_pop_rate(times, senders, binwidth=1.):
    """
        times: spike times
        senders: associated sender
        binwidth: how many ms to bin in
    """

    t_min = np.min(times)
    t_max = np.max(times)
    N_senders = len(np.unique(senders))
    rate, bin_edges = np.histogram(times, np.arange(t_min, t_max, binwidth))
    return rate * 1000 / (N_senders * binwidth), bin_edges[:-1]


def split_in_ex(times, senders):
    """
    Args:  Splits into two arrays at neuron ID 800
        times: spike times
        senders: associated sender
    """
    exc_sender = senders[senders <= 800]
    exc_times = times[senders <= 800]
    inh_sender = senders[senders > 800]
    inh_times = times[senders > 800]
    return exc_times, exc_sender, inh_times, inh_sender


def calc_specgram(time, rate, NFFT=1024, noverlap=900):
    """
    time bin borders to calculate the frequencies
    rate envelope
    NFFT: FFT window in datapoints
    noverlap: Overlap window in datapoints
    """
    Pxx, freqs, bins, im = plt.specgram(rate, NFFT=NFFT, Fs=1000. / (time[2] - time[1]), noverlap=noverlap)
    return freqs, Pxx, bins


def load_json(fname):
    f = open(fname, 'r')
    data = json.load(f)
    return data


def convert_line(line):
    """
    takes a line from a groupfile and converts it to a proper data format
    """
    split_line = line.split(',')
    N, L = split_line.pop(0).split()  # Number of neurons in group and max Layer in the first line
    fired = []  # save all time/sender pairs of the group, these items are in pairs of two
    while len(split_line[0].split()) < 4:
        s, t = split_line.pop(0).split()
        fired.append({'neuron_id': int(s), 't_fired': int(t)})
    links = []
    while len(split_line[0].split()) == 4:  # save all links, writte nin pairs of 4, pre post delay and layer
        pr, po, d, l = split_line.pop(0).split()
        links.append({'pre': int(pr), 'post': int(po), 'delay': int(d), 'layer': int(l)})

    group = dict(
        N_fired=N,
        L_max=L,
        fired=fired,
        links=links
    )
    return group

def get_t_s(group):
    #get the firing patern form a group
    times = []
    senders = []
    for i in group['fired']:
        times.append(i['t_fired'])
        senders.append(i['neuron_id'])

    return np.array(times).astype(float), np.array(senders).astype(float)


def format_spiketrains(times, senders):
    '''

    Parameters
    ----------
    times   : array of float
            spike times
    senders : array of float
            id corresponding to spike times

    Returns
    spiketrains: dictionary
        dictionary of spiketrains,
        keys are ids
        correspondint items are spike times of that id
    -------

    '''
    spiketrains = dict()
    for id in np.unique(senders):
        spiketrains[id - 1] = times[senders == id].tolist()
    return spiketrains


def read_spikefile(filename):
    """
    read in spike files dependent on the file ending (Izhikevic orders the sender/times differently from nest)
    """
    if '.dat' in filename:
        if 'bitwise' in filename:
            spikes = np.loadtxt(filename)
            times = spikes[:, 1]
            senders = spikes[:, 0]

        else:
            spikes = np.loadtxt(filename)
            times = spikes[:, 0]
            senders = spikes[:, 1]
    else:

        spikes = np.loadtxt(filename)
        times = spikes[:, 1]
        senders = spikes[:, 0]
    return times, senders


def read_weightfile(filename):
    with open(filename, "r") as f:
        all_data = json.load(f)
    weights = [i['weight'] for i in all_data]

    return np.array(weights)


def read_group_file(filename):
    with open(filename, "r") as f:
        if '.json' in filename:
            groups = json.load(f)
        else:
            groups = []
            for line in f.readlines():
                groups.append(convert_line(line))
    return groups


def save_json(fname, json_data):
    with open(fname, "w") as f:
        json.dump(json_data, f)
