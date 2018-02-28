
N_inh = 200
N_ex = 800
N_syn = 100
tau_plus = 20
tau_minus = 20


ex_neuron_model = {'consistent_integration': False,
                   'V_th': 30.,
                   'U_m': -0.2 * 65.0,
                   'b': 0.2,
                   'c': -65.0,
                   'a': 0.02,
                   'd': 8.0,
                   'tau_minus': tau_minus}

inh_neuron_model = {'consistent_integration': False,
                    'V_th': 30.,
                    'U_m': -0.2 * 65.0,
                    'b': 0.2,
                    'c': -65.0,
                    'a': 0.1,
                    'd': 2.0,
                    'tau_minus': tau_minus}
