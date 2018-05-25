import matplotlib
matplotlib.use('Agg')
import numpy as np
import sys
import helper
import json
import argparse
import pylab as plt
import seaborn as sns 
import plot_helper as phf
phf.latexify(columns=1)


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', nargs='+', type=str)
parser.add_argument('-o', '--output', type=str)

args = parser.parse_args()

plt.figure()

for i in args.input:
    with open(i, 'r+') as f:
        data = json.load(f)
        if 'qualitative_model' == data['label'] or  'bitwise' in data['label'] or 'initial' in data['label']: # or 'reso' in data['label']:
            if 'qualitative' in data['label']:
                plt.plot(data['dt'], data['dw'], label=data['label'], linewidth=1.,zorder=10)
            else:
                plt.plot(data['dt'], data['dw'], label=data['label'],linewidth=3.)
plt.subplots_adjust(left=0.25, right=0.99, top=0.95, bottom=0.2, hspace=0.2, wspace=0.25)
plt.axhline(0, color='k', linestyle='--',zorder=0)
plt.axvline(0, color='k', linestyle='--',zorder=0)
#plt.legend(loc=2, prop={'size': 12})
plt.xlabel(r"$\Delta t [ms]$")
plt.ylabel(r"$\frac{\Delta w}{w} [mV]$")
plt.savefig(args.output)


