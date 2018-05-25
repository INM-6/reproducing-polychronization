import matplotlib
matplotlib.use('Agg')
import numpy as np
import json
import sys
import pylab as plt
import seaborn as sbn
import os
import helper
import argparse

sbn.set_context('paper', font_scale=3.5, rc={"lines.linewidth": 3.5})
sbn.set_style('whitegrid', {"axes.linewidth": 3.5})
plt.rcParams['font.weight'] = 'bold'

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_izh', nargs='+', type=str)
parser.add_argument('-g', '--groups_izh', nargs='+', type=str)
parser.add_argument('-r', '--conn_rand', nargs='+', type=str)
parser.add_argument('-k', '--groups_rand', nargs='+', type=str)
parser.add_argument('-s', '--conn_rand_EE', nargs='+', type=str)
parser.add_argument('-t', '--groups_rand_EE', nargs='+', type=str)
parser.add_argument('-o', '--output', type=str)
parser.add_argument('-c', '--config', type=str)        # Experiment file defining the network structure and dynamics
parser.add_argument('-e', '--EE', type=int)      # EE synapses only 

args = parser.parse_args()


# load config file
cfg = helper.parse_config(args.config)
Wmax = cfg["network-params"]["plasticity"]["Wmax"]

ratios_izh = []
num_groups_izh = []


def analyzeForRatio(conn_fn, group_fn, Wmax, EE):
    # calculate ratio of strong synapses
    # and number of groups

    # load data from experiments
    with open(conn_fn, "r+") as f:
        conns = json.load(f)

    with open(group_fn, "r+") as f:
        groups = json.load(f)


    exc_weights = np.array([c['weight'] for c in conns if c['weight'] >= 0])

    ratio_strong = len(np.where(exc_weights > Wmax * 0.95)[0]) / float(len(exc_weights))

    return([ratio_strong, len(groups)])



# for all experiments
for i, in_fn in enumerate(args.input_izh): 

    try:
        ratio_strong, num_groups = analyzeForRatio(in_fn, args.groups_izh[i], Wmax, args.EE)

        ratios_izh.append(ratio_strong)
        num_groups_izh.append(num_groups)

    except:
        print("no data available")
        pass



ratios_rand = []
num_groups_rand = []


# for all experiments
for i, c_rand_fn in enumerate(args.conn_rand): 

    try:
        ratio_strong, num_groups = analyzeForRatio(c_rand_fn, args.groups_rand[i], Wmax, args.EE)
    
        num_groups_rand.append(num_groups)
        ratios_rand.append(ratio_strong)

    except:
        
        print("no data available random", c_rand_fn, args.groups_rand[i])
        pass



ratios_rand_EE = []
num_groups_rand_EE = []


# for all experiments
for i, c_rand_fn in enumerate(args.conn_rand_EE): 

    try:
        ratio_strong, num_groups = analyzeForRatio(c_rand_fn, args.groups_rand_EE[i], Wmax, args.EE)
    
        num_groups_rand_EE.append(num_groups)
        ratios_rand_EE.append(ratio_strong)

    except:
        
        print("no data available random", c_rand_fn, args.groups_rand[i])
        pass



plt.figure(figsize=[16, 10])
plt.plot(ratios_rand_EE, num_groups_rand_EE)
plt.plot(ratios_rand, num_groups_rand)
plt.plot(ratios_izh, num_groups_izh, 'k*', markersize=30)
plt.xlabel('ratio strong synapses')
plt.ylabel('#groups')
#plt.ylim([0, 10000])

plt.savefig(args.output)


