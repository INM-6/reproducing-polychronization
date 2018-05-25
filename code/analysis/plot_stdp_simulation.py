import numpy as np
import sys
import helper
import json
import argparse
import pylab as plt
import seaborn as sbn

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', nargs='+', type=str)
parser.add_argument('-o', '--output', type=str)

args = parser.parse_args()

colors = ['red', 'green', 'blue', 'black', 'yellow', 'orange', 'white']

plt.figure(figsize=[16,10])

idx = 0

for i in args.input:
    with open(i, 'r+') as f:
        data = json.load(f)
        if "additive_stdp" in data['label'] or 'quali' in data['label'] or 'bitwise' in data['label'] or 'naive' in data['label'] or 'multi' in data['label']:
            #sbn.tsplot(data['w'], color=colors[idx])

            for w in data['w']:
                plt.plot(range(len(w)), w, color=colors[idx], alpha=0.1) #, label=data['label'])
            idx += 1

idx = 0
for i in args.input:
    with open(i, 'r+') as f:
        data = json.load(f)
        if "additive_stdp" in data['label'] or 'quali' in data['label'] or 'bitwise' in data['label'] or 'naive' in data['label'] or 'multi' in data['label']:
            sbn.tsplot(data['w'], color=colors[idx])
            plt.plot(0, 1, label=data['label'], color=colors[idx])

            #plt.fill_between(range(len(w)), np.mean(data['w'], axis=0) - np.var(data['w']), np.mean(data['w'], axis=0) + np.var(data['w']), alpha=0.5, color=colors[idx]) 
            #plt.plot(range(len(w)), np.mean(data['w'], axis=0), linewidth=2., color=colors[idx], label=data['label'], alpha=1.0) #, label=data['label'])
            idx += 1

plt.legend()
plt.xlabel("t (s)")
plt.ylabel("w")
plt.savefig(args.output)


