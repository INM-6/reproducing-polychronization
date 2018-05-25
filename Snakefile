import os,glob,sys
import numpy as np

sys.path.insert(0, 'code/NEST_model/') #ugly but not sure how to otherwise handle this

import socket
if "cluster" in socket.gethostname():
    shell.prefix('module load autotools; module load mpi/openmpi/1.10.0;')
    NUM_THREADS=1
else:
    NUM_THREADS=1

#Define folders:
CUR_DIR=os.getcwd()
CODE_DIR='code'
DATA_DIR='data'

nest_prefix='NEST_model'
izhi_prefix='original_model'

NEST_CODE_DIR=os.path.join(CODE_DIR,nest_prefix)
NEST_DATA_DIR=os.path.join(DATA_DIR,nest_prefix)

IZHI_CODE_DIR=os.path.join(CODE_DIR,izhi_prefix)
IZHI_DATA_DIR=os.path.join(DATA_DIR,izhi_prefix)

#compile poly_spnet into a folder because it outputs into ../
IZHI_EXEC_DIR=os.path.join(IZHI_CODE_DIR,'exec')
ANA_DIR=os.path.join(CODE_DIR,'analysis')
NEST_SRC_DIR=os.path.join(CUR_DIR,os.path.join(
            CODE_DIR,'nest/nest-simulator'))

PLOT_FILES = ['dynamic_measures.png']
MAN_DIR='manuscript/8538120cqhctwxyjvvn'
FIG_DIR='figures'
LOG_DIR='logs'
CONFIG_DIR=os.path.join(NEST_CODE_DIR,'experiments')

high_CONFIG_FILES=[file[:-5] for file in os.listdir(CONFIG_DIR) if ('bitwise' in file) or ('qualitative_model.yaml' == file)]
CONFIG_FILES=[file[:-5] for file in os.listdir(CONFIG_DIR)]

CONFIG_FILES_group_finder_orig = [file[:-5] for file in os.listdir(CONFIG_DIR) if not ('delay' in file) and not ('resolution' in file) and not ('high_res' in file)]
CONFIG_FILES_group_finder_nest = [file[:-5] for file in os.listdir(CONFIG_DIR) if ('delay' in file)  or ('qualitative' in file)  or ('resolution' in file) or ('bitwise' in file)]

EXPERIMENTS_FOR_STDP_WINDOW = [file[:-5] for file in os.listdir(CONFIG_DIR)]

NUM_REP=range(10)
high_NUM_REP=range(100)

#RANDOM_RATIOS = np.round(np.linspace(0.1, 0.7, 7), 4)
RANDOM_RATIOS = [0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5]


group_stat_list_orig=  np.unique(expand("{folder}/{experiment}/{rep}/stats_orig.json",
                            folder=NEST_DATA_DIR, experiment=high_CONFIG_FILES, rep=high_NUM_REP) \
                            +expand("{folder}/{experiment}/{rep}/stats_orig.json",
                            folder=NEST_DATA_DIR, experiment=CONFIG_FILES_group_finder_orig, rep=NUM_REP))

group_stat_list_python=group_finder_nest = expand("{folder}/{experiment}/{rep}/stats.json",
                            folder=NEST_DATA_DIR, experiment=CONFIG_FILES_group_finder_nest, rep=NUM_REP)

import numpy as np
include: "Izhikevic.rules"
include: "nest.rules"

rule all:
    input:
        stat_files_orig = group_stat_list_orig,

        stat_files_python = group_stat_list_python,

	    spikes = expand("{folder}/{experiment}/{rep}/spikes-1001.gdf",folder=NEST_DATA_DIR,experiment=CONFIG_FILES,rep=NUM_REP),
	    high_spikes = expand("{folder}/{experiment}/{rep}/spikes-1001.gdf",folder=NEST_DATA_DIR,experiment=high_CONFIG_FILES,rep=high_NUM_REP),



        stdp_plot = "figures/stdp_windows.pdf",
        neuron_dynamics = "figures/neuron_dynamics.pdf",
        group_statistics = "figures/group_stats.pdf",

        rand_bitwise = "figures/bitwise_reproduction/random_groups.pdf",
        rand_resolution = "figures/resolution_0p1_W_pspmatched/random_groups_nest.pdf",

        bimodal_bitwise = "figures/bitwise_reproduction_bimodalgamma.pdf",
        bimodal_qualitative = "figures/qualitative_model_bimodalgamma.pdf",

        bitwise_original_comp='figures/bitwise_original_0.pdf',
        bitwise_initial_comp='figures/bitwise_initial_reproduction_0.pdf',
        bitwise_qualitative_comp='figures/bitwise_qualitative_model_0.pdf',


        plots = expand("{folder}/{experiment}/{rep}/plot_dynamics_{experiment}.pdf",
                            folder=FIG_DIR,experiment=CONFIG_FILES,rep=[0]),


rule clean:
    shell:
        """
        rm -rf {data}/*
        rm -rf {code}/*.dat
        rm -rf {exec}/*
        rm -rf {fig}/*
        rm -rf {logs}/*
        """.format(exec=IZHI_EXEC_DIR,fig=FIG_DIR,data=DATA_DIR,code=IZHI_CODE_DIR,logs=LOG_DIR)

rule compile_find_polychronous_groups:
	output:
	    expand('{folder}/find_polychronous_groups',folder=ANA_DIR)
	input:
	    expand('{folder}/find_polychronous_groups.cpp',folder=ANA_DIR)
	shell:
	    'g++ -o {output} {input} -ljsoncpp'


rule find_groups:
    output:
        "{folder}/{experiment}/{rep}/groups.json"
    input:
        connectivity="{folder}/{experiment}/{rep}/connectivity.json",
        program=rules.compile_find_polychronous_groups.output,
    log: 'logs/find_groups_{experiment}_{rep}.log'
    shell:
        '{input.program} {input.connectivity} {output} &>{log} || true'


rule find_groups_random:
    output:
        "{folder}/{experiment}/random/{r}/groups_random.json"
    input:
        connectivity="{folder}/{experiment}/random/{r}/connectivity_random.json",
        program=rules.compile_find_polychronous_groups.output,
    shell:
        '{input.program} {input.connectivity} {output} || true'



rule find_groups_random_EE:
    output:
        "{folder}/{experiment}/random/{r}/groups_random_EE.json"
    input:
        connectivity="{folder}/{experiment}/random/{r}/connectivity_random_EE.json",
        program=rules.compile_find_polychronous_groups.output,
    shell:
        '{input.program} {input.connectivity} {output} || true'


rule calc_stats:
    output:
        "{folder}/{experiment}/{rep}/stats_orig.json"
    input:
        groups="{folder}/{experiment}/{rep}/groups.json",
    log: 'logs/calculate_stats_{experiment}_{rep}.log'
    shell:
        'python {ANA_DIR}/gather_stats.py -g {{input.groups}} -o {{output}} &>{{log}}'.format(ANA_DIR=ANA_DIR)


rule plot_group_stats:
    #comp stands for the experiemtn we want to compare with, i.e. in our case initial and qualitative
    input:
        original_stats=group_stat_list_orig,
        python_stats=group_stat_list_python
    output:
        'figures/group_stats.pdf',
    priority: 9
    shell:
        'python3 {ANA_DIR}/plot_group_statistics.py -glo {{input.original_stats}} -glp {{input.python_stats}} --output {{output}}'.format(ANA_DIR=ANA_DIR,fig_dir=FIG_DIR)


rule plot_bitwise_comp:
    #comp stands for the experiemtn we want to compare with, i.e. in our case initial and qualitative
    input:
        comp_con=expand('{folder}/{{experiment}}/{{rep}}/connectivity.json',folder=NEST_DATA_DIR),
        bit_con=expand('{folder}/bitwise_reproduction/{{rep}}/connectivity.json',folder=NEST_DATA_DIR),
        comp_spk=expand('{folder}/{{experiment}}/{{rep}}/spikes-1001.gdf',folder=NEST_DATA_DIR),
        bit_spk=expand('{folder}/bitwise_reproduction/{{rep}}/spikes-1001.gdf',folder=NEST_DATA_DIR),
    output:
        'figures/bitwise_{experiment}_{rep}.pdf',
    priority: 9
    shell:
        'python3 {ANA_DIR}/plot_bitwise_comp.py -bs {{input.bit_spk}} -cs {{input.comp_spk}} -bw {{input.bit_con}} -cw {{input.comp_con}} -fn {{output}}'.format(ANA_DIR=ANA_DIR,fig_dir=FIG_DIR)


rule plot_bitwise_original:
    input:
        original_spk=expand('{folder}/bitwise_reproduction/{{rep}}/spikes.dat',folder=IZHI_DATA_DIR),
        nest_mem=expand('{folder}/bitwise_reproduction/{{rep}}/membrane_potential-1002.dat',folder=NEST_DATA_DIR),
        nest_spk=expand('{folder}/bitwise_reproduction/{{rep}}/spikes-1001.gdf',folder=NEST_DATA_DIR),
    output:
        'figures/bitwise_original_{rep}.{ext,(eps|png|pdf|jpg)}',
    priority: 10
    shell:
        'python3 {ANA_DIR}/plot_bitwise_original.py -bs {{input.nest_spk}} -os {{input.original_spk}} -bmem {{input.nest_mem}} -fn {{output}}'.format(ANA_DIR=ANA_DIR,fig_dir=FIG_DIR)

rule plot_bimodal_gamma:
    output:
        outfile=expand('{folder}/{{experiment}}_bimodalgamma.{{ext,(eps|png|pdf|jpg)}}',folder=FIG_DIR),
    input:
        connectivity=expand('{folder}/{{experiment}}/{rep}/connectivity.json',folder=NEST_DATA_DIR,rep=high_NUM_REP),
        spikes=expand('{folder}/{{experiment}}/{rep}/spikes-1001.gdf',folder=NEST_DATA_DIR,rep=high_NUM_REP),
        groups=expand('{folder}/{{experiment}}/{rep}/stats_orig.json',folder=NEST_DATA_DIR,rep=high_NUM_REP),

    priority: 2
    run:
        shell("""
        python3 code/analysis/plot_bimodal_gamma.py \
        -cl {input.connectivity}\
        -sl {input.spikes}\
        -gl {input.groups}\
        --output {output.outfile}
        """)


rule plot_dynamics:
    output:
        file=expand('{folder}/{{experiment}}/{{rep}}/plot_dynamics_{{experiment}}.{{ext,(eps|png|pdf)}}',folder=FIG_DIR),

    input:
        connectivity=expand('{folder}/{{experiment}}/{{rep}}/connectivity.json',folder=NEST_DATA_DIR),
        spikes=expand('{folder}/{{experiment}}/{{rep}}/spikes-1001.gdf',folder=NEST_DATA_DIR),
    priority: 2
    run:
        shell("""
        python3 code/analysis/plot_dynamics.py \
        --spikefile {input.spikes}\
        --weightfile {input.connectivity}\
        --filename {output}
        """)

rule plot_stdp_window:
    input:
        weights=expand('{folder}/{{experiment}}/stdp_window.json'.format(folder=NEST_DATA_DIR), experiment=EXPERIMENTS_FOR_STDP_WINDOW),
        program='{folder}/plot_stdp_window.py'.format(folder=ANA_DIR),
    output:
        plot="{folder}/stdp_windows.pdf".format(folder=FIG_DIR),
    shell:
        """
        python {input.program} -i {input.weights} -o {output.plot}
        """

rule randomize_conn:
    input:
        conf = '{nest_folder}/experiments/{{experiment}}.yaml'.format(nest_folder=NEST_CODE_DIR),
        conns = 'data/NEST_model/{experiment}/0/connectivity.json',
    output:
        fn = 'data/NEST_model/{experiment}/random/{r}/connectivity_random.json',
    shell:
        'python code/analysis/randomize_conn.py -i {input.conns} -c {input.conf} -r {wildcards.r} -o {output.fn} -e 0'
    

rule plot_random_groups:
    input:
        conns=expand('{folder}/{{experiment}}/{rep}/connectivity.json', folder=NEST_DATA_DIR, rep=NUM_REP),
        groups=expand('{folder}/{{experiment}}/{rep}/groups.json', folder=NEST_DATA_DIR, rep=NUM_REP),
        conn_rand=expand('{folder}/{{experiment}}/random/{rand}/connectivity_random.json', folder=NEST_DATA_DIR, rand=RANDOM_RATIOS),
        groups_rand=expand('{folder}/{{experiment}}/random/{rand}/groups_random.json', folder=NEST_DATA_DIR, rand=RANDOM_RATIOS),
        conn_rand_EE=expand('{folder}/{{experiment}}/random/{rand}/connectivity_random_EE.json', folder=NEST_DATA_DIR, rand=RANDOM_RATIOS),
        groups_rand_EE=expand('{folder}/{{experiment}}/random/{rand}/groups_random_EE.json', folder=NEST_DATA_DIR, rand=RANDOM_RATIOS),
        conf='{nest_folder}/experiments/{{experiment}}.yaml'.format(nest_folder=NEST_CODE_DIR),
    output:
        fn = 'figures/{experiment}/random_groups.pdf',
    shell:
        'python code/analysis/plot_random_conn.py -i {input.conns} -g {input.groups} -k {input.groups_rand} -r {input.conn_rand} -s {input.conn_rand_EE} -t {input.groups_rand_EE} -c {input.conf} -o {output.fn} -e 0'


rule plot_random_groups_nest:
    input:
        conns=expand('{folder}/{{experiment}}/{rep}/connectivity.json', folder=NEST_DATA_DIR, rep=NUM_REP),
        groups=expand('{folder}/{{experiment}}/{rep}/groups_nest.json', folder=NEST_DATA_DIR, rep=NUM_REP),
        conn_rand=expand('{folder}/{{experiment}}/random/{rand}/connectivity_random.json', folder=NEST_DATA_DIR, rand=RANDOM_RATIOS),
        groups_rand=expand('{folder}/{{experiment}}/random/{rand}/groups_nest_random.json', folder=NEST_DATA_DIR, rand=RANDOM_RATIOS),
        conn_rand_EE=expand('{folder}/{{experiment}}/random/{rand}/connectivity_random_EE.json', folder=NEST_DATA_DIR, rand=RANDOM_RATIOS),
        groups_rand_EE=expand('{folder}/{{experiment}}/random/{rand}/groups_nest_random_EE.json', folder=NEST_DATA_DIR, rand=RANDOM_RATIOS),
        conf='{nest_folder}/experiments/{{experiment}}.yaml'.format(nest_folder=NEST_CODE_DIR),
    output:
        fn = 'figures/{experiment}/random_groups_nest.pdf',
    shell:
        'python code/analysis/plot_random_conn.py -i {input.conns} -g {input.groups} -k {input.groups_rand} -r {input.conn_rand} -s {input.conn_rand_EE} -t {input.groups_rand_EE} -c {input.conf} -o {output.fn} -e 0'



rule randomize_conn_EE:
    input:
        conf = '{nest_folder}/experiments/{{experiment}}.yaml'.format(nest_folder=NEST_CODE_DIR),
        conns = 'data/NEST_model/{experiment}/0/connectivity.json',
    output:
        fn = 'data/NEST_model/{experiment}/random/{r}/connectivity_random_EE.json',
    shell:
        'python code/analysis/randomize_conn.py -i {input.conns} -c {input.conf} -r {wildcards.r} -o {output.fn} -e 1'
    

