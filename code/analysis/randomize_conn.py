import json
import numpy as np
import sys
import helper
import copy

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-o', '--out', type=str)        # Experiment file defining the network structure and dynamics
parser.add_argument('-i', '--connectivity', type=str)  # Connectivity file the group finding is based on
parser.add_argument('-r', '--random', type=float)      # ratio of strong synapses 
parser.add_argument('-c', '--config', type=str)      # config file for the network/resolution etc.
parser.add_argument('-e', '--EE', type=int)      # EE synapses only 

args = parser.parse_args()

# load config file
cfg = helper.parse_config(args.config)

Wmax = cfg["network-params"]["plasticity"]["Wmax"]
Winh = cfg["network-params"]["plasticity"]["W_inh"]

# load connectivity data
with open(args.connectivity, "r+") as f:
    conns_ = json.load(f)

conns = copy.deepcopy(conns_)

# and select the relevant connections (EE)
if args.EE == 1:
    exc_conns = np.array([conns[i] for i in range(len(conns)) if conns[i]['weight'] >= 0. and conns[i]['post'] <= 800])
else:
    exc_conns = np.array([conns[i] for i in range(len(conns)) if conns[i]['weight'] >= 0.])

# set all exc_conns to 0 
for c in exc_conns:
    c['weight'] = 0 

# randomly pick synapses with given probability
strong_conns = np.random.choice(exc_conns, int(len(exc_conns) * args.random), replace=False)

# make them strong
for c in strong_conns:
    c['weight'] = Wmax


with open(args.out, "w+") as f:
    json.dump(conns, f)




