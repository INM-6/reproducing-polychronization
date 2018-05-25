import matplotlib
matplotlib.use('Agg')
import nest
from nest import voltage_trace as vtr
import numpy as np
import pylab as plt
import seaborn as sns
import argparse
import plot_helper as phf
import matplotlib.gridspec as gridspec

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', nargs='+', type=str)
parser.add_argument('-o', '--output', type=str)

args = parser.parse_args()

neuron_params_low_res = {'consistent_integration': False,
                         'V_th': 30.,
                         'U_m': -0.2 * 65.0,
                         'b': 0.2,
                         'c': -65.0,
                         'a': 0.02,
                         'd': 8.0,
                         'integration_steps': 1}

neuron_params_high_res = {'consistent_integration': False,
                          'V_th': 30.,
                          'U_m': -0.2 * 65.0,
                          'b': 0.2,
                          'c': -65.0,
                          'a': 0.02,
                          'd': 8.0,
                          'integration_steps': 10}


# run low resolution simulation
nest.ResetKernel()
nest.SetKernelStatus({'resolution': 1.0})

neuron_low_res = nest.Create("izhikevich", 1, neuron_params_low_res) 
neuron_high_res = nest.Create("izhikevich", 1, neuron_params_high_res) 

nest.SetStatus(neuron_low_res, {"I_e": 4.})
nest.SetStatus(neuron_high_res, {"I_e": 4.})

mm = nest.Create('multimeter', 1, {'record_from': ['V_m', 'U_m'], 'interval': 1.0})
sd = nest.Create('spike_detector', 1)

nest.Connect(neuron_low_res, sd)
nest.Connect(neuron_high_res, sd)

nest.Connect(mm, neuron_low_res + neuron_high_res)

nest.Simulate(1000.)

ev = nest.GetStatus(sd, 'events')[0]

neuron0_mask = np.where(ev['senders'] == 1)[0]
neuron1_mask = np.where(ev['senders'] == 2)[0]

neuron0_times = ev['times'][neuron0_mask] 
neuron1_times = ev['times'][neuron1_mask]

isi0 = np.diff(neuron0_times)
isi1 = np.diff(neuron1_times)

var0 = np.var(isi0)
var1 = np.var(isi1)

mean0 = np.mean(isi0)
mean1 = np.mean(isi1)

print("isi 0", var0/mean0, 1./mean0)
print("isi 1", var1/mean1, 1./mean1)

ev_mm = nest.GetStatus(mm, 'events')[0]
senders = ev_mm['senders']
times = ev_mm['times']
vms = ev_mm['V_m']
ums = ev_mm['U_m']

mask_low = np.where(senders == neuron_low_res)[0]
mask_high = np.where(senders == neuron_high_res)[0]

Vs_current_low_low = vms[mask_low] 
Vs_current_low_high = vms[mask_high] 

Us_current_low_low = ums[mask_low] 
Us_current_low_high = ums[mask_high] 

times_current_low_low = times[mask_low]
times_current_low_high = times[mask_high]




# run high resolution simulation
nest.ResetKernel()
nest.SetKernelStatus({'resolution': 0.1})

neuron_low_res = nest.Create("izhikevich", 1, neuron_params_low_res) 

nest.SetStatus(neuron_low_res, {"I_e": 4.})

mm = nest.Create('multimeter', 1, {'record_from': ['V_m', 'U_m'], 'interval': 0.1})
sd = nest.Create('spike_detector', 1)

nest.Connect(neuron_low_res, sd)

nest.Connect(mm, neuron_low_res)

nest.Simulate(1000.)

ev = nest.GetStatus(sd, 'events')[0]

neuron0_mask = np.where(ev['senders'] == 1)[0]

neuron0_times = ev['times'][neuron0_mask] 

isi0 = np.diff(neuron0_times)

var0 = np.var(isi0)

mean0 = np.mean(isi0)

print("isi 0", var0/mean0, 1./mean0)

ev_mm = nest.GetStatus(mm, 'events')[0]
senders = ev_mm['senders']
times = ev_mm['times']
vms = ev_mm['V_m']
ums = ev_mm['U_m']

mask_low = np.where(senders == neuron_low_res)[0]

Vs_current_high_low = vms[mask_low] 
Us_current_high_low = ums[mask_low] 

times_current_high_low = times[mask_low]







# run high resolution simulation
nest.ResetKernel()
nest.SetKernelStatus({'resolution': 1.0})

neuron_low_res = nest.Create("izhikevich", 1, neuron_params_low_res) 
neuron_high_res = nest.Create("izhikevich", 1, neuron_params_high_res) 

sg = nest.Create("spike_generator", 1, {"spike_times": [50., 50.]})

nest.Connect(sg, neuron_low_res + neuron_high_res, 'all_to_all', {'weight': 10.})

mm = nest.Create('multimeter', 1, {'record_from': ['V_m', 'U_m'], 'interval': 1.0})
sd = nest.Create('spike_detector', 1)

nest.Connect(neuron_low_res, sd)
nest.Connect(neuron_high_res, sd)

nest.Connect(mm, neuron_low_res)
nest.Connect(mm, neuron_high_res)

nest.Simulate(100.)

ev = nest.GetStatus(sd, 'events')[0]
print(ev)


neuron0_mask = np.where(ev['senders'] == 1)[0]

neuron0_times = ev['times'][neuron0_mask] 

isi0 = np.diff(neuron0_times)

var0 = np.var(isi0)

mean0 = np.mean(isi0)

ev = nest.GetStatus(mm, 'events')[0]

print("isi 0", var0/mean0, 1./mean0)


ev_mm = nest.GetStatus(mm, 'events')[0]
senders = ev_mm['senders']
times = ev_mm['times']
vms = ev_mm['V_m']
ums = ev_mm['U_m']

mask_low = np.where(senders == neuron_low_res)[0]
mask_high = np.where(senders == neuron_high_res)[0]

Vs_spike_low_low = vms[mask_low] 
Vs_spike_low_high = vms[mask_high] 

Us_spike_low_low = ums[mask_low] 
Us_spike_low_high = ums[mask_high] 

times_spike_low_low = times[mask_low]
times_spike_low_high = times[mask_high]





# run high resolution simulation
nest.ResetKernel()
nest.SetKernelStatus({'resolution': 0.1})

neuron_low_res = nest.Create("izhikevich", 1, neuron_params_low_res) 

sg = nest.Create("spike_generator", 1, {"spike_times": [50., 50.]})

nest.Connect(sg, neuron_low_res, 'all_to_all', {'weight': 85.})

mm = nest.Create('multimeter', 1, {'record_from': ['V_m', 'U_m'], 'interval': 0.1})
sd = nest.Create('spike_detector', 1)

nest.Connect(neuron_low_res, sd)

nest.Connect(mm, neuron_low_res)

nest.Simulate(100.)

ev = nest.GetStatus(sd, 'events')[0]
print(ev)

neuron0_mask = np.where(ev['senders'] == 1)[0]

neuron0_times = ev['times'][neuron0_mask] 

isi0 = np.diff(neuron0_times)

var0 = np.var(isi0)

mean0 = np.mean(isi0)

print("isi 0", var0/mean0, 1./mean0)

ev_mm = nest.GetStatus(mm, 'events')[0]
senders = ev_mm['senders']
times = ev_mm['times']
vms = ev_mm['V_m']
ums = ev_mm['U_m']

mask_low = np.where(senders == neuron_low_res)[0]

Vs_spike_high_low = vms[mask_low] 
Us_spike_high_low = ums[mask_low] 

times_spike_high_low = times[mask_low]

phf.latexify(columns=2)


fig = plt.figure()
gs0 = gridspec.GridSpec(2, 2)

ax0 = plt.subplot(gs0[0, 0])
ax1 = plt.subplot(gs0[1, 0])
ax2 = plt.subplot(gs0[0, 1])
ax3 = plt.subplot(gs0[1, 1])


ax0.plot(times_current_low_low, Vs_current_low_low, label="low low")
ax0.plot(times_current_low_high, Vs_current_low_high, label="low high")
ax0.plot(times_current_high_low, Vs_current_high_low, label="high low")
ax0.set_ylabel("V [mV]")
ax0.set_xticks([])

ax1.plot(times_current_low_low, Us_current_low_low, label="low low")
ax1.plot(times_current_low_high, Us_current_low_high, label="low high")
ax1.plot(times_current_high_low, Us_current_high_low, label="high low")
ax1.set_ylabel("U [a.u.]")
ax1.set_xlabel("Time [ms]")

ax2.plot(times_spike_low_low, Vs_spike_low_low, label="low low")
ax2.plot(times_spike_low_high, Vs_spike_low_high, label="low high")
ax2.plot(times_spike_high_low, Vs_spike_high_low, label="high low")
ax2.set_xticks([])

ax3.plot(times_spike_low_low, Us_spike_low_low, label="low low")
ax3.plot(times_spike_low_high, Us_spike_low_high, label="low high")
ax3.plot(times_spike_high_low, Us_spike_high_low, label="high low")
ax3.set_xlabel("Time [ms]")

gs0.update(left=0.08, right=0.99, top=0.95, bottom=0.1, hspace=0.2, wspace=0.25)
for ax, letter in [(ax0, 'A'), (ax1, 'B'), (ax2, 'C'),(ax3,'D')]:
    ax.annotate(r'\textbf{{{letter}}}'.format(letter=letter), xy=(-0.15, 0.99), xycoords='axes fraction', fontsize=10,
                horizontalalignment='left', verticalalignment='top', annotation_clip=False)


plt.savefig(args.output)







