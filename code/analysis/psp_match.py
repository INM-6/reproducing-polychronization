import nest
import numpy as np




def build_model(weight=10.,resolution=1.0):
    nest.ResetKernel()
    print(resolution)
    nest.SetKernelStatus({'resolution': resolution,
                          'print_time': True })
    print(nest.GetStatus([0]))
    neuron_model = 'izhikevich'
    nest.CopyModel(neuron_model, 'inh_Izhi', {'consistent_integration': False,
                                              'U_m': -0.2 * 65.0,
                                              'b': 0.2,
                                              'c': -65.0,
                                              'a': 0.1,
                                              'd': 2.0,
                                              'tau_minus': 20.})
    nest.CopyModel(neuron_model, 'ex_Izhi', {'consistent_integration': False,
                                             'U_m': -0.2 * 65.0,
                                             'b': 0.2,
                                             'c': -65.0,
                                             'a': 0.02,
                                             'd': 8.0,
                                             'tau_minus': 20.})

    nest.CopyModel("static_synapse", "EX", {'weight': weight, 'delay': 1.0})

    ex_neuron = nest.Create('ex_Izhi', 1)
    spike_gen = nest.Create('spike_generator', 1)
    nest.SetStatus(spike_gen, {'spike_times': [250.]})

    mm = nest.Create("multimeter", params={
        'record_from': ['V_m'],
        'withgid': True,
        'withtime': True,
        'to_memory': True,
        'to_file': False,
        'label': 'mem_pot'})
    nest.SetStatus(ex_neuron,'V_m',-70.)
    nest.SetStatus(ex_neuron,'U_m',-70.*.2)

    nest.Connect(mm, ex_neuron, 'all_to_all')
    nest.Connect(spike_gen, ex_neuron, 'all_to_all','EX')

    nest.Simulate(500.)
    return nest.GetStatus(mm)[0]['events']

def synchrony_influence(weight=10.,resolution=0.1,dt=20):
    nest.ResetKernel()
    print(resolution)
    nest.SetKernelStatus({'resolution': resolution,
                          'print_time': True })
    print(nest.GetStatus([0]))
    neuron_model = 'izhikevich'
    nest.CopyModel(neuron_model, 'inh_Izhi', {'consistent_integration': False,
                                              'U_m': -0.2 * 65.0,
                                              'b': 0.2,
                                              'c': -65.0,
                                              'a': 0.1,
                                              'd': 2.0,
                                              'tau_minus': 20.})
    nest.CopyModel(neuron_model, 'ex_Izhi', {'consistent_integration': False,
                                             'U_m': -0.2 * 65.0,
                                             'b': 0.2,
                                             'c': -65.0,
                                             'a': 0.02,
                                             'd': 8.0,
                                             'tau_minus': 20.})

    nest.CopyModel("static_synapse", "EX", {'weight': weight, 'delay': 1.0})

    ex_neuron = nest.Create('ex_Izhi', dt)
    spike_gen = nest.Create('spike_generator', dt)
    for i in range(dt):
        nest.SetStatus([spike_gen[i]], {'spike_times': [250.,250+i*resolution]})

    mm = nest.Create("multimeter", params={
        'record_from': ['V_m'],
        'interval':resolution,
        'withgid': True,
        'withtime': True,
        'to_memory': True,
        'to_file': False,
        'label': 'mem_pot'})

    spk_det = nest.Create("spike_detector", params={
        'to_memory': True,
        'to_file': False,
        'label': 'spk_det'})

    nest.SetStatus(ex_neuron,'V_m',-70.)
    nest.SetStatus(ex_neuron,'U_m',-70.*.2)

    nest.Connect(mm, ex_neuron, 'all_to_all')
    nest.Connect(spike_gen, ex_neuron, 'one_to_one','EX')
    nest.Connect(ex_neuron, spk_det, 'all_to_all','EX')

    nest.Simulate(500.)
    return nest.GetStatus(mm)[0]['events'],nest.GetStatus(spk_det)[0]['events']


#############################################
# First we fit the weight to match the psp for the neurons at the different resolutions
#
#
###############################################
#build_model(weight=10.,resolution=0.1)
res0p1=build_model(weight=85.,resolution=0.1)
res1p0=build_model(weight=10.,resolution=1.0)
res0p1_spike=build_model(weight=300.,resolution=0.1)
res1p0_spike=build_model(weight=20.,resolution=1.0)
res0p1_in=build_model(weight=-35.,resolution=0.1)
res1p0_in=build_model(weight=-5.,resolution=1.0)




import matplotlib.pyplot as plt
plt.plot(res0p1['times'],res0p1['V_m']+70.)
#plt.plot(res0p1['times'],np.cumsum(res0p1['V_m']+70.))

plt.plot(res1p0['times'],res1p0['V_m']+70.)
#plt.plot(res1p0['times'],np.cumsum(res1p0['V_m']+70.))


plt.xlim([240,350])

plt.savefig('wo_spike.png')
plt.close()

plt.plot(res0p1_spike['times'],res0p1_spike['V_m'])
#plt.plot(res0p1_spike['times'],np.cumsum(res0p1_spike['V_m']+70.))

plt.plot(res1p0_spike['times'],res1p0_spike['V_m'])
#plt.plot(res1p0_spike['times'],np.cumsum(res1p0_spike['V_m']+70.))
plt.xlim([240,350])
plt.savefig('w_spike.png')
plt.close()

plt.plot(res0p1_in['times'],res0p1_in['V_m'])
plt.plot(res0p1_in['times'],np.cumsum(res0p1_in['V_m']+70.))

plt.plot(res1p0_in['times'],res1p0_in['V_m'])
plt.plot(res1p0_in['times'],np.cumsum(res1p0_in['V_m']+70.))
plt.xlim([240,350])
plt.savefig('inhib.png')
plt.close()



#############################################
# Now we see the influence on synchrony
#
#
###############################################

#build_model(weight=10.,resolution=0.1)
res0p1_sync,spk=synchrony_influence(weight=85.,resolution=0.1)
print(spk)
Blues = plt.get_cmap('rainbow')
f, axarr = plt.subplots(2)

import matplotlib.pyplot as plt
for sender in np.unique(res0p1_sync['senders']):
    sender_idxes=res0p1_sync['senders']==sender

    dt=(sender-min(res0p1_sync['senders']))*0.1
    c=dt*10/(max(res0p1_sync['senders'])-min(res0p1_sync['senders']))
    print(dt,c,max((res0p1_sync['V_m'][sender_idxes])))
    axarr[0].plot(res0p1_sync['times'][sender_idxes],res0p1_sync['V_m'][sender_idxes],color=Blues(dt))
    axarr[1].plot(dt, max((res0p1_sync['V_m'][sender_idxes])), color=Blues(dt),marker='o')
    if spk['times'][spk['senders']==sender]:
        axarr[0].plot(spk['times'][spk['senders']==sender], max((res0p1_sync['V_m'][sender_idxes]))+10, color=Blues(dt),marker='*',zorder=10)
        axarr[1].plot(dt, max((res0p1_sync['V_m'][sender_idxes])) + 10, color=Blues(dt),
              marker='*', zorder=10)

axarr[0].set_xlim([240,300])

plt.savefig('synchrony_influence_0p1.pdf')
plt.close()


#build_model(weight=10.,resolution=0.1)
res0p1_sync,spk=synchrony_influence(weight=10.,resolution=1.0)
print(spk)
Blues = plt.get_cmap('rainbow')
f, axarr = plt.subplots(2)

for sender in np.unique(res0p1_sync['senders']):
    sender_idxes=res0p1_sync['senders']==sender

    dt=(sender-min(res0p1_sync['senders']))
    c=dt/(max(res0p1_sync['senders'])-min(res0p1_sync['senders']))
    print(dt,c,max((res0p1_sync['V_m'][sender_idxes])))
    axarr[0].plot(res0p1_sync['times'][sender_idxes],res0p1_sync['V_m'][sender_idxes],color=Blues(dt))
    axarr[1].plot(dt, max((res0p1_sync['V_m'][sender_idxes])), color=Blues(dt),marker='o')
    if spk['times'][spk['senders']==sender]:
        axarr[0].plot(spk['times'][spk['senders']==sender], max((res0p1_sync['V_m'][sender_idxes]))+10, color=Blues(dt),marker='*',zorder=10)
        axarr[1].plot(dt, max((res0p1_sync['V_m'][sender_idxes])) + 10, color=Blues(dt),
              marker='*', zorder=10)

axarr[0].set_xlim([240,300])

plt.savefig('synchrony_influence_1p0.pdf')
plt.close()