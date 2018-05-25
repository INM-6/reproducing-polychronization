import os
import numpy as np
import sys
import matplotlib.pyplot as plt
sys.path.append('code/analysis')
import helper as hf
import plot_helper as phf
import pandas as pd
import seaborn as sns
import json
import ijson
import matplotlib.mlab as mlab
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--groupfile', type=str)
parser.add_argument('-o', '--outfolder', type=str)
args = parser.parse_args()




# if os.path.getsize(os.path.join(result_folder,file)) < 422099208 * 1.1:

NTL = phf.return_NTL(args.groupfile)
if NTL:
    N, T, L=NTL[0],NTL[1],NTL[2]
    print(len(N), np.median(N))
    N_groups = len(N)
    N_fired = np.median(N)
    longest_path = np.median(L)
    time_span = np.median(T)
    stats = dict(N_fired=N,
                 longest_path=L,
                 time_span=T
                 )
    with open(args.outfolder, 'w+') as fs:
        json.dump(stats, fs)

else:
    with open(args.outfolder, "w+") as f:
        json.dump({'Failed': 1}, f)

