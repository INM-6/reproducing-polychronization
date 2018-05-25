import json
import numpy as np
import itertools
from multiprocessing import Pool, TimeoutError
import time
import sys
import nest
import copy
import argparse
import helper

# global parameter used through the whole script
global t_sim, N, Ne, Ni, M, min_delay, max_delay, sim_resolution
global exc_conns, exc_pre, exc_post, exc_weight, exc_delay
global inh_conns, inh_pre, inh_post, inh_weight, inh_delay


parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config', type=str)        # Experiment file defining the network structure and dynamics
parser.add_argument('-o', '--output', type=str)        # output groups file
parser.add_argument('-n', '--num_threads', type=int)   # Number of threads
parser.add_argument('-i', '--connectivity', type=str)  # Connectivity file the group finding is based on
parser.add_argument('-s', '--statistics', type=str)    # group statistics write out in a separate file so they have to be computet every time

args = parser.parse_args()

# load config file
cfg = helper.parse_config(args.config)


in_fn = args.connectivity
max_num_processes = args.num_threads
sim_resolution = cfg["simulation-params"]["resolution"]

Wmax = cfg["network-params"]["plasticity"]["Wmax"]
Winh = cfg["network-params"]["plasticity"]["W_inh"]

out_fn = args.output
out_stat_fn = args.statistics

# load connectivity data
with open(in_fn, "r+") as f:
    final_stdw = json.load(f)

# and select the relevant connections
# only strong exc. and all inh connections
exc_conns = np.array([final_stdw[i] for i in range(len(final_stdw)) if final_stdw[i]['weight'] > Wmax * 0.95])
inh_conns = np.array([final_stdw[i] for i in range(len(final_stdw)) if final_stdw[i]['weight'] == Winh])

# dissamble connections into components
exc_pre = np.array([int(c['pre']) for c in exc_conns])
exc_post = np.array([int(c['post']) for c in exc_conns])
exc_weight = np.array([float(c['weight']) for c in exc_conns])
exc_delay = np.array([float(c['delay']) for c in exc_conns])
min_delay, max_delay = np.min(exc_delay), np.max(exc_delay)

inh_pre = np.array([int(c['pre']) for c in inh_conns])
inh_post = np.array([int(c['post']) for c in inh_conns])
inh_weight = np.array([float(c['weight']) for c in inh_conns])
inh_delay = np.array([float(c['delay']) for c in inh_conns])

t_sim = 1000.0  # simulate only the first second

N = 1000  # total number of neurons
Ne = 800  # number of excitatory neurons
Ni = 200  # number of inhibitory neurons
M = 100  # number of outgoing connections per neuron


def create_network():
    #
    # This function resets and recreates the network. 
    # 

    global t_sim, N, Ne, Ni, M, min_delay, max_delay, sim_resolution
    global exc_conns, exc_pre, exc_post, exc_weight, exc_delay
    global inh_conns, inh_pre, inh_post, inh_weight, inh_delay

    nest.ResetKernel()
    nest.sr("M_ERROR setverbosity")
    nest.SetKernelStatus({'resolution': sim_resolution})

    # build all neurons but only selected connections exc_conns, inh_conns

    # set default values of izhikevich neuron model
    # that excitatory and inhibitory neurons have in common
    nest.SetDefaults('izhikevich', {'b': 0.2,
                                    'c': -65.0,
                                    'V_m': -70.0,
                                    'U_m': -70.0 * 0.2,
                                    'V_th': 30.0,
                                    'consistent_integration': False,
                                    'integration_steps': cfg['simulation-params']['neuron-integration-steps']})

    # create excitatory and inhibitory neurons and set the parameters
    # that excitatory and inhibitory neurons do not have in common
    exc_neurons = nest.Create('izhikevich', Ne, {'a': 0.02, 'd': 8.0})
    inh_neurons = nest.Create('izhikevich', Ni, {'a': 0.1, 'd': 2.0})

    # create list with global ids of all neurons
    all_neurons = exc_neurons + inh_neurons

    # create spike detectors for excitatory and inhibitory neurons
    sd = nest.Create('spike_detector', 1, {'to_file': False, 'label': 'spikes'})

    # create excitatory connections
    nest.Connect(exc_pre, exc_post, {'rule': 'one_to_one'},
                 syn_spec={'model': 'static_synapse', 'weight': exc_weight, 'delay': exc_delay})

    # create inhibitory connections
    nest.Connect(inh_pre, inh_post, {'rule': 'one_to_one'},
                 syn_spec={'model': 'static_synapse', 'weight': inh_weight, 'delay': inh_delay})

    # get all connections
    connections = nest.GetConnections(all_neurons, all_neurons)

    # connect neurons to spike detector
    nest.Connect(all_neurons, sd, 'all_to_all')

    return [exc_neurons, inh_neurons, sd]


# for every excitatory neuron and each possible triplet of excitatory presynaptic neurons
# define the connections that are initially activated and trigger the simulation
def worker(pivot_neuron):
    
    # starting from a pivot neuron (mother neuron)
    import nest

    global t_sim, N, Ne, Ni, M, min_delay, max_delay, sim_resolution
    global exc_conns, exc_pre, exc_post, exc_weight, exc_delay
    global inh_conns, inh_pre, inh_post, inh_weight, inh_delay

    # reset the json data and spike generators
    local_json_data = []
    sgs = None

    # get all incoming excitatory connections to the pivot neuron
    # only strong connections are considered
    inc_exc_conns = exc_conns[np.where(exc_post == pivot_neuron)[0]]

    # convinient function returns triplets of those connections
    # iterate over all triplets (anchor neurons)
    for stim_triplet_id, stim_triplet in enumerate(itertools.combinations(inc_exc_conns, 3)):

        # as creating the network takes some time we usually would create the network only once
        # but as we need to create spike generators for each triplet, which slows down the simulation after a while,
        # we reset and re-create the network after 100 triplets

        if stim_triplet_id % 100 == 0:
            #print("progress", pivot_neuron, stim_triplet_id, float(stim_triplet_id) /
            #      len(list(itertools.combinations(inc_exc_conns, 3))))
            # reset and recreate network and reset spike generators
            exc_neurons, inh_neurons, sd = create_network()
            sgs = None

        # ready data fields for this triplet
        stim_target_gids = []
        stim_times = []
        stim_weights = []
        group = []
        t_fired = []
        group_delays = []

        # calculate the maximal delay between the anchor neurons
        # This is done because in the original group finder only post synaptic neurons
        # after the mother neuron receive the synaptic events from the anchor neurons
        max_delay_triplet = np.max(np.array([c['delay'] for c in stim_triplet]))

        # determine the initially activated neurons 
        for st in sorted(stim_triplet, key=lambda x: x['delay'], reverse=True):
            # iterate over the anchor neurons sorted by delay to the mother neuron

            # by selecting all excitatory connections from the current anchor neuron
            # having a delay larger than the delay to the mother neuron
            stim_conns = exc_conns[np.where(np.all([exc_pre == st['pre'], exc_delay >= st['delay']], axis=0))[0]]

            # disassemble these connections to their components
            stim_pre = np.array([int(c['pre']) for c in stim_conns])
            stim_post = np.array([int(c['post']) for c in stim_conns])
            stim_weight = np.array([float(c['weight']) for c in stim_conns])
            stim_delay = np.array([float(c['delay']) for c in stim_conns])

            # calculate the offet of stimulation relative to the mother neuron
            stim_offset = max_delay_triplet - st['delay']
            stim_target_gids.extend(stim_post)
            stim_times.extend(stim_offset + stim_delay)
            stim_weights.extend(stim_weight)

            # store anchor neurons in this group
            group.append(int(st['pre']))
            t_fired.append(stim_offset + 1)
            group_delays.append(st['delay'])

            # reorder the group (just for cosmetic reasons)
            order = np.argsort(group)
            sorted_group = np.array(group)[order].tolist()
            sorted_t_fired = np.array(t_fired)[order].tolist()
            sorted_stim_delay = np.array(group_delays)[order].tolist()

            # reset neuron parameters
            nest.SetStatus(exc_neurons + inh_neurons, {'b': 0.2,
                                                       'c': -65.0,
                                                       'V_m': -70.0,
                                                       'U_m': -70.0 * 0.2,
                                                       'V_th': 30.0,
                                                       'consistent_integration': False,
                                                       'integration_steps': cfg['simulation-params'][
                                                           'neuron-integration-steps']
                                                       })

            nest.SetStatus(exc_neurons, {'a': 0.02, 'd': 8.0})
            nest.SetStatus(inh_neurons, {'a': 0.1, 'd': 2.0})

        # convert to numpy arrays
        stim_target_gids = np.array(stim_target_gids)
        stim_times = np.array(stim_times)
        stim_weights = np.array(stim_weights)
        stim_delay = np.array(sorted_stim_delay)

        group = sorted_group
        t_fired = sorted_t_fired

        # at this point we created the network and determined the neurons and timepoints when to activate them
        # next, start the simulation and activate the neurons at the calculated times


        ######################
        #    FIND GROUPS     #
        ######################

        nest.ResetNetwork()
        nest.SetKernelStatus({"time": 0.0})

        # every spike_generator is responsible for a specific stimulation time in [sim_resolution, max_delay]
        # create a spike generator for each stim time
        if not sgs == None:
            nest.SetStatus(list(sgs.values()), {"spike_times": []})

        # create spike generators
        stim_times_unique = np.unique(stim_times)
        sgs = dict(zip(stim_times_unique, nest.Create('spike_generator', len(stim_times_unique))))

        # and connect them to the respective targets using the corresponding weights
        for t_stim in sgs.keys():
            nest.SetStatus([sgs[t_stim]], {"spike_times": [t_stim]})
            idxs = np.where(stim_times == t_stim)[0]
            if idxs.size:
                nest.Connect(np.array(len(idxs) * [sgs[t_stim]], np.int64), stim_target_gids[idxs],
                             {'rule': 'one_to_one'}, syn_spec={'model': 'static_synapse',
                                                               "weight": stim_weights[idxs], 
                                                               "delay": len(idxs) * [min_delay]})

        # simulate for sim_time 
        nest.Simulate(t_sim)

        # extract spike times and corresponding global ids from spike detectors
        # and add them to the group
        t_fired.extend(nest.GetStatus(sd, 'events')[0]['times'])
        group.extend(nest.GetStatus(sd, 'events')[0]['senders'])

        # find the links and determine the layers
        N_fired = len(group)
        L_max = 0
        links = np.array([])
        json_group = None

        # drop huge groups for computational reasons
        if N_fired > 1000:
            continue

        # get connectivity of the group
        exc_is_in_group = np.in1d(exc_pre, group) & np.in1d(exc_post, group)
        inh_is_in_group = np.in1d(inh_pre, group) & np.in1d(inh_post, group)

        # disassable these connections to their components
        group_exc_pre = exc_pre[exc_is_in_group]
        group_exc_post = exc_post[exc_is_in_group]
        group_exc_delay = exc_delay[exc_is_in_group]
        group_inh_pre = inh_pre[inh_is_in_group]
        group_inh_post = inh_post[inh_is_in_group]
        group_inh_delay = inh_delay[inh_is_in_group]

        # iterate over all spikes in the group after the three anchor neurons
        for i in range(3, N_fired):
            # i is index of Nth neuron in group
            # group[i] is index of neuron

            # loop through all neurons that fired before this one because they are candidates for linked neurons
            for j in range(i):

                # calculate the delay from neuron j to neuron i
                if group[j] <= Ne:
                    delays = group_exc_delay[
                        np.where(np.all([group_exc_pre == group[j], group_exc_post == group[i]], axis=0))[0]]
                else:
                    delays = group_inh_delay[
                        np.where(np.all([group_inh_pre == group[j], group_inh_post == group[i]], axis=0))[0]]

                # there should only be one delay, except there are multapses in the network which we exclude
                # interate over all found connections between j to j
                for d in delays:
                    # anchor neurons have layer = 1, all others have minimum layer = 2 
                    layer = 2

                    # calculate the layer of a link as the maximal layer of its presnaptic neuron + 1
                    if links.size > 0:
                        # get all links of the presynaptic partner
                        idxs = np.where(np.array(links)[:, 1] == group[j])[0]
                        # get their max layer
                        if idxs.size:
                            layer = int(np.max(np.array(links)[idxs, 3]) + 1)

                            # set the global maximal layer of the current group
                            if layer > L_max:
                                L_max = layer

                        # add link
                        links = np.append(links, [[group[j], group[i], d, layer]], axis=0)
                    else:
                        # if there are no links in this group yet, this link must have layer 2
                        links = np.array([[group[j], group[i], d, layer]])

        
        # viability test as in the original algorithm
        # a group is only then considered if all anchor neurons project to at least one
        # active, excitatory neurons in the group except the mother neuron
        discard = 0
        used = np.zeros(3)
        for i in range(len(used)):
            for j in range(len(links)):
                if links[j][0] == group[i] and links[j][1] <= Ne:
                    used[i] += 1
            if used[i] == 1:
                discard = 1

        # consider only groups with a maximal layer larger than 7
        #
        # as the calculation of layers is slightly different than in the original algorithm
        # we find more groups
        # We believe that in the original code the last spike is sometimes not recorded
        # therefore the threshold of 7 layers is not crossed so the group in the original one is not counted
        if L_max >= 7 and discard == 0:
            # save group in JSON format
            json_group = {}
            json_group["N_fired"] = N_fired
            json_group["L_max"] = L_max

            json_fired = []
            for i in range(N_fired):
                json_fire = {}
                json_fire["neuron_id"] = int(group[i]) - 1
                json_fire["t_fired"] = float(t_fired[i]) - 1
                json_fired.append(json_fire)
            json_group["fired"] = json_fired

            json_links = []
            for j in range(len(links)):
                json_link = {}
                json_link["pre"] = int(links[j][0]) - 1
                json_link["post"] = int(links[j][1]) - 1
                json_link["delay"] = float(links[j][2]) - 1
                json_link["layer"] = int(links[j][3])
                json_links.append(json_link)
            json_group["links"] = json_links

            # print("group found", pivot_neuron, json_group["N_fired"], json_group["L_max"])

        if not json_group == None:
            local_json_data.append(json_group)

    return local_json_data


# generate a field for all groups
json_data = []

# parallelize the algorithm
pool = Pool(processes=max_num_processes)

# This eusures a file is created even when the job fails(which sometimes happens)
# mainly needed to keep snakemake happy because it works based on output / input files

with open(out_fn, "w+") as f:
    json.dump({'Failed':1}, f)
with open(out_stat_fn, "w+") as f:
    json.dump({'Failed':1}, f)


# execute the algorithm for each excitatory neuron as mother neuron
for found_groups in pool.imap_unordered(worker, range(1, Ne + 1)):
    json_data += found_groups

    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    print('!!!!!!!',len(json_data),' groups found so far !!!!!!!!!!')
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')

# calculate statistics about the groups for later use
N_list = []
L_list = []
T_list = []
i = 0

# iterate over all found groups 
for i, g in enumerate(json_data):
    # get times of group
    times, _ = helper.get_t_s(g)

    # add stats for each group
    N_list.append(int(g["N_fired"]))
    T_list.append(max(times))  # time span
    L_list.append(int(g["L_max"]))  # longest path

stats = dict(N_fired=N_list,
             longest_path=L_list,
             time_span=T_list
             )

# save groups file
with open(out_fn, "w+") as f:
    json.dump(json_data, f)

# save statistics file
with open(out_stat_fn, "w+") as fs:
    json.dump(stats, fs)
