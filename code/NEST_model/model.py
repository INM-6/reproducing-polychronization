import nest
import json
import os
import sys
import argparse
import params
import numpy as np

# ugly but not sure how to otherwise handle this
sys.path.insert(0, 'code/analysis/')
import helper


parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config', type=str)
parser.add_argument('-n', '--num_threads', type=int)
parser.add_argument('-r', '--repetition', type=int)

args = parser.parse_args()

cfg = helper.parse_config(args.config)


def write_weights(neuron, fname):
    json_data = []
    for n in neuron:
        conns = nest.GetConnections([n], neuron)

        for c in conns:
            json_conn = {}
            json_conn["pre"] = nest.GetStatus([c], "source")[0]
            json_conn["post"] = nest.GetStatus([c], "target")[0]
            json_conn["delay"] = nest.GetStatus([c], "delay")[0]
            json_conn["weight"] = nest.GetStatus([c], "weight")[0]
            json_data.append(json_conn)

    with open(fname, "w") as f:
        json.dump(json_data, f)


def connect_network(ex_neuron, inh_neuron, cfg):
    conf = cfg["network-params"]["connectivity"]
    if conf["type"] == "reproduce":
        file = helper.load_json(conf["from-file"].format(rep=args.repetition))

        delay = np.array([float(i['delay']) for i in file])
        pre = np.array([i['pre'] for i in file])
        post = np.array([i['post'] for i in file])

        for pre_neuron in np.unique(pre):
            idxes = np.where(pre_neuron == pre)[0]
            if pre_neuron in ex_neuron:

                nest.Connect([pre_neuron], post[idxes].tolist(),
                             conn_spec='all_to_all', syn_spec='EX')

            elif pre_neuron in inh_neuron:
                nest.Connect([pre_neuron], post[idxes].tolist(),
                             conn_spec='all_to_all', syn_spec='II')

        for pr, po, de in zip(pre, post, delay):
            conn = nest.GetConnections(source=[pr], target=[po])
            nest.SetStatus(conn, 'delay', de)

    elif conf["type"] == "generate":
        conn_dict_inh = {'rule': 'fixed_outdegree',
                         'outdegree': params.N_syn, 'multapses': False, 'autapses': False}
        conn_dict_ex = {'rule': 'fixed_outdegree',
                        'outdegree': params.N_syn, 'multapses': False, 'autapses': False}
        nest.Connect(inh_neuron, ex_neuron,
                     conn_spec=conn_dict_inh, syn_spec='II')
        nest.Connect(ex_neuron, neurons, conn_spec=conn_dict_ex, syn_spec='EX')

        if conf["delay-distribution"] == "uniform-non-random":
            delay_range = range(
                int(conf["delay-range"][0]), int(conf["delay-range"][1]))
            delay_list = [
                [j for i in range(int(100 / len(delay_range)))] for j in delay_range]
            delay_list = np.reshape(
                np.array(delay_list).astype(float), (1, 100))[0]
            for n in ex_neuron:
                conns = nest.GetConnections(
                    source=[n], target=np.random.permutation(ex_neuron + inh_neuron).tolist())

                nest.SetStatus(conns, 'delay', delay_list)

        elif conf["delay-distribution"] == "uniform-random":

            for n in ex_neuron:
                delay_list = np.random.uniform(int(conf["delay-range"][0]), int(conf["delay-range"][1] - 1), 100).astype(
                    float)
                delay_list = np.around(
                    delay_list, decimals=int(-np.log10(cfg["simulation-params"]["resolution"])))
                nest.SetStatus(nest.GetConnections(
                    source=[n], target=ex_neuron + inh_neuron), 'delay', delay_list)
                print(delay_list, np.min(delay_list), np.max(delay_list))
    else:
        pass


def set_stimulus(neurons, conf, sim_time):
    def set_stimulus_times(stim_t, stim_id):

        nest.SetStatus([neurons[stim_id[0] - 1]], 'I', 20.)
        # otherwise stimulus must occur at -1ms
        stim_id = stim_id[stim_t > 0]
        stim_t = stim_t[stim_t > 0]
        random_input = nest.Create('spike_generator', len(neurons))
        nest.Connect(random_input, neurons, 'one_to_one')
        nest.SetStatus(nest.GetConnections(random_input), 'weight', 20.)
        for i in np.unique(stim_id):
            idx = stim_id == i
            times = stim_t[idx]
            nest.SetStatus([random_input[int(i - 1)]], {'spike_times': times})
        del stim_id
        del stim_t

    if conf["type"] == "reproduce":
        stimulus = np.loadtxt(conf["from-file"].format(rep=args.repetition))
        stim_id = stimulus[:, 0].astype(int)
        stim_t = stimulus[:, 1] - 1
        del stimulus
        set_stimulus_times(stim_t, stim_id)
        print('reproduce', conf)
        # first neuron gets current manually rest via spike generator

    elif conf["type"] == "generate":

        if conf["distribution"] == "poisson":
            print('generate poisson', conf)

            random_input = nest.Create('poisson_generator')
            nest.SetStatus(random_input, params={'rate': conf["rate"]})
            nest.Connect(random_input, neurons, 'all_to_all',
                         {'weight': conf["weight"]})
        elif conf["distribution"] == "original":
            print('generate original', conf)

            stim_id, stim_t = np.random.choice(1000, int(sim_time)), np.array(
                np.linspace(0, int(sim_time), int(sim_time) + 1))
            stim_t = stim_t[:-1]
            set_stimulus_times(stim_t, stim_id)

        else:
            pass

    else:
        pass


def set_initial_conditions(neurons, conf):
    if conf["type"] == "reproduce":
        initials = np.loadtxt(conf["from-file"].format(rep=args.repetition))
        stim_id = initials[:, 0]
        stim_v = initials[:, 1]
        stim_u = initials[:, 2]
        # first neuron gets current manually rest via spike generator
        nest.SetStatus(neurons, 'V_m', stim_v)
        nest.SetStatus(neurons, 'U_m', stim_u)
    elif conf["type"] == "generate":
        if conf["distribution"] == "uniform":

            vms = np.random.uniform(
                conf["V_m-range"][0], conf["V_m-range"][1], len(neurons))
            ums = np.array([0.2 * v for v in vms])

            nest.SetStatus(neurons, 'V_m', vms)
            nest.SetStatus(neurons, 'U_m', ums)
        else:
            pass
    else:
        pass


nest.ResetKernel()

seed = [np.random.randint(0, 9999999)] * args.num_threads

nest.SetKernelStatus({'resolution': cfg["simulation-params"]["resolution"],
                      'print_time': True,
                      'rng_seeds': seed,
                      'local_num_threads': args.num_threads,
                      'overwrite_files': True,
                      'syn_update_interval': cfg["simulation-params"]["synapse-update-interval"]})

nest.CopyModel('izhikevich', 'ex', params.ex_neuron_model)

nest.CopyModel('izhikevich', 'inh', params.inh_neuron_model)

if cfg["simulation-params"]["resolution"] == 1.0:
    inh_weight = -5.
else:
    inh_weight = -35.

nest.CopyModel("static_synapse", "II", {'weight': inh_weight, 'delay': 1.0})

if cfg["network-params"]["plasticity"]["synapse-model"] == 'stdp_izh_synapse':
    nest.CopyModel(cfg["network-params"]["plasticity"]["synapse-model"], "EX", {
        'weight': cfg["network-params"]["plasticity"]['W_init'],
        'Wmax': cfg["network-params"]["plasticity"]['Wmax'],
        'LTP': cfg["network-params"]["plasticity"]['LTP'],
        'LTD': cfg["network-params"]["plasticity"]['LTD'],
        "tau_syn_update_interval": cfg["network-params"]["plasticity"]["tau_syn_update_interval"],
        "constant_additive_value": cfg["network-params"]["plasticity"]["constant_additive_value"],
        "reset_weight_change_after_update": cfg["network-params"]["plasticity"]["reset_weight_change_after_update"]

    })

elif cfg["network-params"]["plasticity"]["synapse-model"] == 'stdp_synapse':
    mu_plus = None
    mu_minus = None
    if cfg["network-params"]["plasticity"]["stdp-type"] == 'additive':
        mu_plus = .0
        mu_minus = 0.0
    elif cfg["network-params"]["plasticity"]["stdp-type"] == 'multiplicative':
        mu_plus = 1.0
        mu_minus = 1.0
    else:
        pass

    nest.CopyModel(cfg["network-params"]["plasticity"]["synapse-model"], "EX", {
        'weight': cfg["network-params"]["plasticity"]['W_init'],
        'tau_plus': params.tau_plus,
        'lambda': cfg["network-params"]["plasticity"]['lambda'],
        'alpha': cfg["network-params"]["plasticity"]['alpha'],
        'mu_plus': mu_plus,
        'mu_minus': mu_minus,
        'Wmax': cfg["network-params"]["plasticity"]['Wmax']
    })

else:
    nest.CopyModel(cfg["network-params"]["plasticity"]["synapse-model"], "EX", {
        'weight': 6.,
        'consistent_integration': False,
    })

ex_neuron = nest.Create('ex', params.N_ex)
inh_neuron = nest.Create('inh', params.N_inh)
neurons = ex_neuron + inh_neuron
label = os.path.join(cfg["simulation-params"]["data-path"],
                     cfg["simulation-params"]["data-prefix"])
label = os.path.join(label, str(args.repetition))

spikedetector = nest.Create("spike_detector", params={
    'start': cfg["simulation-params"]["sim-time"] - cfg["simulation-params"]["rec_spikes"],
    'withgid': True,
    'withtime': True,
    'to_memory': False,
    'to_file': True,
    'label': os.path.join(label, 'spikes')})

nest.Connect(neurons, spikedetector, 'all_to_all')

mm = nest.Create("multimeter", params={
    'record_from': ['V_m', 'U_m'],
    'withgid': True,
    'withtime': True,
    'to_memory': False,
    'to_file': True,
    'precision': 17,
    'start': cfg["simulation-params"]["sim-time"] - cfg["simulation-params"]["rec_mem"],
    # 'stop': 1000. * 100.,
    'label': os.path.join(label, 'membrane_potential')})
nest.Connect(mm, neurons, 'all_to_all')

set_initial_conditions(neurons, cfg["network-params"]["initial-state"])

connect_network(ex_neuron, inh_neuron, cfg)

set_stimulus(neurons, cfg["network-params"]["stimulus"],
             cfg["simulation-params"]["sim-time"])


nest.Simulate(cfg["simulation-params"]["sim-time"])
write_weights(neurons, os.path.join(label, 'connectivity.json'))
