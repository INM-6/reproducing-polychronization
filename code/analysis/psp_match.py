import nest
import numpy as np


def build_model(weight=10., resolution=1.0):
    nest.ResetKernel()
    print(resolution)
    nest.SetKernelStatus({'resolution': resolution,
                          'print_time': True})
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
    nest.SetStatus(ex_neuron, 'V_m', -70.)
    nest.SetStatus(ex_neuron, 'U_m', -70. * .2)

    nest.Connect(mm, ex_neuron, 'all_to_all')
    nest.Connect(spike_gen, ex_neuron, 'all_to_all', 'EX')

    nest.Simulate(500.)
    return nest.GetStatus(mm)[0]['events']


# build_model(weight=10.,resolution=0.1)
res0p1 = build_model(weight=85., resolution=0.1)
res1p0 = build_model(weight=10., resolution=1.0)
res0p1_spike = build_model(weight=170., resolution=0.1)
res1p0_spike = build_model(weight=20., resolution=1.0)
res0p1_in = build_model(weight=-35., resolution=0.1)
res1p0_in = build_model(weight=-5., resolution=1.0)


import matplotlib.pyplot as plt
plt.plot(res0p1['times'], res0p1['V_m'] + 70.)
plt.plot(res0p1['times'], np.cumsum(res0p1['V_m'] + 70.))

plt.plot(res1p0['times'], res1p0['V_m'] + 70.)
plt.plot(res1p0['times'], np.cumsum(res1p0['V_m'] + 70.))


plt.xlim([240, 350])

plt.savefig('wo_spike.png')
plt.close()

plt.plot(res0p1_spike['times'], res0p1_spike['V_m'] + 70.)
plt.plot(res0p1_spike['times'], np.cumsum(res0p1_spike['V_m'] + 70.))

plt.plot(res1p0_spike['times'], res1p0_spike['V_m'] + 70.)
plt.plot(res1p0_spike['times'], np.cumsum(res1p0_spike['V_m'] + 70.))
plt.xlim([240, 350])
plt.savefig('w_spike.png')
plt.close()

plt.plot(res0p1_in['times'], res0p1_in['V_m'] + 70.)
plt.plot(res0p1_in['times'], np.cumsum(res0p1_in['V_m'] + 70.))

plt.plot(res1p0_in['times'], res1p0_in['V_m'] + 70.)
plt.plot(res1p0_in['times'], np.cumsum(res1p0_in['V_m'] + 70.))
plt.xlim([240, 350])
plt.savefig('inhib.png')
plt.close()
