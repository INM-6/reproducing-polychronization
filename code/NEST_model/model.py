import nest
import json
import os
import sys
import argparse

# ugly but not sure how to otherwise handle this
sys.path.insert(0, 'code/analysis/')
from params import *
import helper

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config', type=str)
parser.add_argument('-n', '--num_threads', type=int)
parser.add_argument('-r', '--repetition', type=int)

args = parser.parse_args()

cfg = helper.parse_config(args.config)


def write_weights(neuron, fname):
    """

    Args:
        neuron: list of neuron Ids
        fname:  filename for output

    Returns:

    """

    #initialize empty data storage
    json_data = []
    #loop through all neurons
    for n in neuron:
        #get all connections of that neuron
        conns = nest.GetConnections([n], neuron)
        # extract relevant parameters of the connection
        for c in conns:
            json_conn = {}
            json_conn["pre"] = nest.GetStatus([c], "source")[0]
            json_conn["post"] = nest.GetStatus([c], "target")[0]
            json_conn["delay"] = nest.GetStatus([c], "delay")[0]
            json_conn["weight"] = nest.GetStatus([c], "weight")[0]
            json_data.append(json_conn)
    #write out in specified output file
    with open(fname, "w") as f:
        json.dump(json_data, f)




def connect_network(ex_neuron, inh_neuron, cfg):
    """

    Args:
        ex_neuron: list o excitatory neuron ids
        inh_neuron:  list of inhibitory neuron ids
        cfg:  config file that defines the network layout, neuron dynamics and plasticity algorithm
            connecitivtes can be:

            "reproduce"

          For comparison we wrote out the exact connections from the original code
          and use them here so we can compare the resulting dynamics
          this step is necessary since we can"t reproduce all aspects of the connectivity drawing algorithm in nest

            "generate"

          The NEST implementation is used, we further distinguish:

          "uniform-non-random"
            where we draw the arbitrary not really random way of the original code
            Here every neuron has M/D connections with delay d where D is max delay and M the amount of post synaptic targets
            Works only when M/D is integer

          "uniform random"
            draws the a uniform distribution between min delay and max delay
            This woeuld be closest to what is usually understood of random disribution of delays
    Returns:

    """
    conf=cfg["network-params"]["connectivity"]
    if conf["type"] == "reproduce":
        ##########################
        #
        #  For comparison we wrote out the exact connections from the original code
        #  and use them here so we can compare the resulting dynamics
        #  this step is necessary since we can"t reproduce all aspects of the connectivity drawing algorithm in nest
        #
        ##########################
        file = helper.load_json(conf["from-file"].format(rep=args.repetition))

        #read out all connections with respective delays
        delay = np.array([float(i['delay']) for i in file])
        pre = np.array([i['pre'] for i in file])
        post = np.array([i['post'] for i in file])

        #loop through all presynaptic neurons and connect them to their targets

        #first connect
        for pre_neuron in np.unique(pre):
            idxes = np.where(pre_neuron == pre)[0]
            if pre_neuron in ex_neuron:

                nest.Connect([pre_neuron], post[idxes].tolist(),
                             conn_spec='all_to_all', syn_spec='EX')

            elif pre_neuron in inh_neuron:
                nest.Connect([pre_neuron], post[idxes].tolist(),
                             conn_spec='all_to_all', syn_spec='II')
        #then set delay
        for pr, po, de in zip(pre, post, delay):
            conn = nest.GetConnections(source=[pr], target=[po])
            nest.SetStatus(conn, 'delay', de)



    elif conf["type"] == "generate":
        ###############################
        #
        #  Use the NEST implementation for drawing connections
        #
        ###############################
        conn_dict_inh = {'rule': 'fixed_outdegree',
                         'outdegree': N_syn, 'multapses': False, 'autapses': False}
        conn_dict_ex = {'rule': 'fixed_outdegree',
                        'outdegree': N_syn, 'multapses': False, 'autapses': False}
        nest.Connect(inh_neuron, ex_neuron,
                     conn_spec=conn_dict_inh, syn_spec='II')
        nest.Connect(ex_neuron, neurons, conn_spec=conn_dict_ex, syn_spec='EX')

        if conf["delay-distribution"] == "uniform-non-random":
            delay_range = range(int(conf["delay-range"][0]), int(conf["delay-range"][1]))
            delay_list = [[j for i in range(int(100 / len(delay_range)))] for j in delay_range]
            delay_list = np.reshape(
                np.array(delay_list).astype(float), (1, 100))[0]
            for n in ex_neuron:
                conns = nest.GetConnections(source=[n], target=np.random.permutation(ex_neuron + inh_neuron).tolist())

                nest.SetStatus(conns, 'delay', delay_list)

        elif conf["delay-distribution"] == "uniform-random":

            for n in ex_neuron:
                delay_list = np.random.uniform(int(conf["delay-range"][0]), int(conf["delay-range"][1]-1), 100).astype(
                    float)
                delay_list=np.around(delay_list,decimals=int(-np.log10(cfg["simulation-params"]["resolution"])))
                nest.SetStatus(nest.GetConnections(
                    source=[n], target=ex_neuron + inh_neuron), 'delay', delay_list)
                print(delay_list,np.min(delay_list),np.max(delay_list))
    else:
        pass


def set_stimulus(neurons, conf, sim_time,Wmax):
    """

    Args:
        neurons: neuron id list
        sim_time: how long the simulation runs in ms
        conf:    stimulus configuration string
        Options are :

        "reproduce"
          stimulus times in original code are written out and read in so we can exaclty reproducce the spiking

        "generate"
          create a stimulus input
            "poisson"
              input statistics according to a poisson input with a certain rate

            "original"
              randomly select one neuron to make it spike every ms of the simulation

    Returns:

    """
    def set_stimulus_times(stim_t, stim_id,Wmax):
        #First ms must be set individually
        # otherwise stimulus spike must have occured at -1ms
        nest.SetStatus([neurons[stim_id[0] - 1]], 'I', 2*Wmax)
        stim_id = stim_id[stim_t > 0]
        stim_t = stim_t[stim_t > 0]

        #create spike generator for every neuron that receives the spike times read out from orginal
        random_input = nest.Create('spike_generator', len(neurons))
        nest.Connect(random_input, neurons, 'one_to_one')
        #20 mv are able to make a neuron spike
        nest.SetStatus(nest.GetConnections(random_input), 'weight', 2*Wmax)

        #set spike times
        for i in np.unique(stim_id):
            idx = stim_id == i
            times = stim_t[idx]
            nest.SetStatus([random_input[int(i - 1)]], {'spike_times': times})
        del stim_id
        del stim_t


    if conf["type"] == "reproduce":
        #load in stimulus/id pairs
        stimulus = np.loadtxt(conf["from-file"].format(rep=args.repetition))
        stim_id = stimulus[:, 0].astype(int)
        stim_t = stimulus[:, 1] - 1#stimulus arrives 1 ms later due to min_delay

        del stimulus

        # first neuron gets current manually rest via spike generator
        set_stimulus_times(stim_t, stim_id,Wmax)
        print('reproduce',conf)

    elif conf["type"] == "generate":

        if conf["distribution"] == "poisson":
            #use poisson statistic input: every neuron gets a individual background input with rate "rate"
            print('generate poisson', conf)

            random_input = nest.Create('poisson_generator')
            nest.SetStatus(random_input, params={'rate': conf["rate"]})
            nest.Connect(random_input, neurons, 'all_to_all', {'weight': conf["weight"]})
        elif conf["distribution"] == "original":
            # Draw inputs according to original code but with numpy and inject them via spike generators
            print('generate original', conf)

            stim_id, stim_t = np.random.choice(N, int(sim_time)), np.array(
                                                                        np.linspace(0, int(sim_time), int(sim_time) + 1))
            stim_t = stim_t[:-1]
            set_stimulus_times(stim_t, stim_id,Wmax)

        else:
            pass

    else:
        pass


def set_initial_conditions(neurons, conf):
    """

    Args:
        neurons: list of neurons
        conf:    defines how the initial distribution is set
            "reproduct"
              uses the initial distributions generated by the priginal code and imports them
            "generate uniform"
               draws accoridng to the same statistics initual distributions
    Returns:

    """
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

            vms = np.random.uniform(conf["V_m-range"][0], conf["V_m-range"][1], len(neurons))
            ums = np.array([0.2 * v for v in vms])

            nest.SetStatus(neurons, 'V_m', vms)
            nest.SetStatus(neurons, 'U_m', ums)
        else:
            pass
    else:
        pass


nest.ResetKernel()

seed = [np.random.randint(0, 9999999)] * args.num_threads

#set resolution random seeds and with how many threads the simulation is run as well as after how many steps the synapse is updated
nest.SetKernelStatus({'resolution': cfg["simulation-params"]["resolution"],
                      'print_time': True,
                      'rng_seeds': seed,
                      'local_num_threads': args.num_threads,
                      'overwrite_files': True,
                      'syn_update_interval': cfg["simulation-params"]["synapse-update-interval"]})

#Neuron model is the original izhikevic neuron as described in the original manuscript
neuron_model = 'izhikevich'

#define parameters for Inh neurons
#Parameters defined in params.py
inh_neuron_model.update({'integration_steps':cfg['simulation-params']['neuron-integration-steps']})
nest.CopyModel(neuron_model, 'inh_Izhi', inh_neuron_model)

#define parameters for excitatory neurons (regular fast spiking)
exc_neuron_model.update({'integration_steps':cfg['simulation-params']['neuron-integration-steps']})
nest.CopyModel(neuron_model, 'ex_Izhi', exc_neuron_model)


#set synapse realted parameters
nest.CopyModel("static_synapse", "II", {'weight': cfg["network-params"]["plasticity"]['W_inh'], 'delay': 1.0})
if cfg["network-params"]["plasticity"]["synapse-model"] == 'stdp_izh_synapse':
    nest.CopyModel(cfg["network-params"]["plasticity"]["synapse-model"], "EX", {
        'weight': cfg["network-params"]["plasticity"]['W_init'],   #defines initial weight distribution (curently only delta distributions)
        'Wmax': cfg["network-params"]["plasticity"]['Wmax'],       # defines cut opff weights for the additive stdp like synapses
        'LTP': cfg["network-params"]["plasticity"]['LTP'],         #LTP max value
        'LTD': cfg["network-params"]["plasticity"]['LTD'],         #LTD max value
        "tau_syn_update_interval": cfg["network-params"]["plasticity"]["tau_syn_update_interval"],   #how much of the synaoptic events are kept in memory
        "constant_additive_value": cfg["network-params"]["plasticity"]["constant_additive_value"],   # constant addivite value that is added after each update interval
        "reset_weight_change_after_update": cfg["network-params"]["plasticity"]["reset_weight_change_after_update"]  # reset weight after update or not (original is not reset)

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
        'tau_plus': tau_plus,
        'lambda': cfg["network-params"]["plasticity"]['lambda'],
        'alpha': cfg["network-params"]["plasticity"]['alpha'],
        'mu_plus': mu_plus,
        'mu_minus': mu_minus,
        'Wmax': cfg["network-params"]["plasticity"]['Wmax']
    })

else:
    nest.CopyModel(cfg["network-params"]["plasticity"]["synapse-model"], "EX", {
        'weight':  cfg["network-params"]["plasticity"]['W_init'],
        'consistent_integration': False,
    })


# create the neurons
ex_neuron = nest.Create('ex_Izhi', N_ex)
inh_neuron = nest.Create('inh_Izhi', N_inh)
neurons = ex_neuron + inh_neuron

#define the labels for the spike/membrane potential file
label = os.path.join(cfg["simulation-params"]["data-path"], cfg["simulation-params"]["data-prefix"])
label = os.path.join(label, str(args.repetition))


#create spike detector for recording
spikedetector = nest.Create("spike_detector", params={
    'start': cfg["simulation-params"]["sim-time"] - cfg["simulation-params"]["rec_spikes"],
    'withgid': True,
    'withtime': True,
    'to_memory': False,
    'to_file': True,
    'label': os.path.join(label, 'spikes')})


#connect spiek detector to the network
nest.Connect(neurons, spikedetector, 'all_to_all')


#create multimeter
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
#connect to the neurons
nest.Connect(mm, neurons, 'all_to_all')




#set all initial conditions
set_initial_conditions(neurons, cfg["network-params"]["initial-state"])

#connect the network accorung to conf
connect_network(ex_neuron, inh_neuron,cfg )

#set stimulus according to conf
set_stimulus(neurons, cfg["network-params"]["stimulus"], cfg["simulation-params"]["sim-time"],cfg["network-params"]["plasticity"]['Wmax'])



nest.Simulate(cfg["simulation-params"]["sim-time"])


write_weights(neurons, os.path.join(label, 'connectivity.json'))
