import os,glob,sys
sys.path.insert(0, 'code/NEST_model/') #ugly but not sure how to otherwise handle this

import socket
if "cluster" in socket.gethostname():
    shell.prefix('module load autotools; module load mpi/openmpi/1.10.0;source activate poly-python3;')
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

high_CONFIG_FILES=[file[:-5] for file in os.listdir(CONFIG_DIR) if ('bitwise' in file) or ('qualitative' in file)]
CONFIG_FILES=[file[:-5] for file in os.listdir(CONFIG_DIR) if ('bitwise' not in file) and ('qualitative' not in file) and ('resolution' not in file) and ('synapse_update_interval_10s' not in file)]
CONFIG_FILES_group_finder_nest=[file[:-5] for file in os.listdir(CONFIG_DIR) if ('delay' in file)  or ('qualitative' in file)]
NUM_REP=range(10)
high_NUM_REP=range(10) # use 100 for more exact statistics

include: "original.rules"
include: "nest.rules"

rule all:
    input:
        group_finder_orig=expand("{folder}/{experiment}/{rep}/groups.json",
                            folder=NEST_DATA_DIR,experiment=CONFIG_FILES,rep=NUM_REP),
        group_finder_orig_high=expand("{folder}/{experiment}/{rep}/groups.json",
                            folder=NEST_DATA_DIR,experiment=high_CONFIG_FILES,rep=high_NUM_REP),
        stats_orig_high=expand("{folder}/{experiment}/{rep}/stats_orig.json",
                            folder=NEST_DATA_DIR,experiment=high_CONFIG_FILES,rep=high_NUM_REP),
        stats_orig=expand("{folder}/{experiment}/{rep}/stats_orig.json",
                            folder=NEST_DATA_DIR,experiment=CONFIG_FILES,rep=NUM_REP),
        polytest_data_full_nest=expand("{folder}/{experiment}/{rep}/groups_nest.json",
                            folder=NEST_DATA_DIR,experiment=CONFIG_FILES_group_finder_nest,rep=NUM_REP),
	    spikes=expand("{folder}/{experiment}/{rep}/spikes-1001.gdf",folder=NEST_DATA_DIR,experiment=CONFIG_FILES,rep=NUM_REP),
	    high_spikes=expand("{folder}/{experiment}/{rep}/spikes-1001.gdf",folder=NEST_DATA_DIR,experiment=high_CONFIG_FILES,rep=high_NUM_REP)





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
        '{input.program} {input.connectivity} {output} &>{log}'

rule calc_stats:
    output:
        "{folder}/{experiment}/{rep}/stats_orig.json"
    input:
        groups="{folder}/{experiment}/{rep}/groups.json",
    log: 'logs/calculate_stats_{experiment}_{rep}.log'
    shell:
        'python {ANA_DIR} {{input.groups}} {output} &>{log}'




rule plot_test_statistical_reproduction:
    input:
        stat_con=expand('{folder}/{{experiment}}_reproduction/{{rep}}/connectivity.json',folder=NEST_DATA_DIR),
        bit_con=expand('{folder}/bitwise_reproduction/{{rep}}/connectivity.json',folder=NEST_DATA_DIR),
        stat_spk=expand('{folder}/{{experiment}}_reproduction/{{rep}}/spikes-1001.gdf',folder=NEST_DATA_DIR),
        bit_spk=expand('{folder}/bitwise_reproduction/{{rep}}/spikes-1001.gdf',folder=NEST_DATA_DIR),
    output:
        'figures/bitwise_{experiment}_{rep}.{ext,(eps|png)}',
    priority: 9
    shell:
        'python3 {ANA_DIR}/plot_statistical_reproduction.py -bs {{input.bit_spk}} -ss {{input.stat_spk}} -bw {{input.bit_con}} -sw {{input.stat_con}} -fn {{output}}'.format(ANA_DIR=ANA_DIR,fig_dir=FIG_DIR)

rule plot_test_bitwise_reproduction:
    input:
        original_spk=expand('{folder}/bitwise_reproduction/{{rep}}/spikes.dat',folder=IZHI_DATA_DIR),
        nest_mem=expand('{folder}/bitwise_reproduction/{{rep}}/membrane_potential-1002.dat',folder=NEST_DATA_DIR),
        nest_spk=expand('{folder}/bitwise_reproduction/{{rep}}/spikes-1001.gdf',folder=NEST_DATA_DIR),
    output:
        'figures/bitwise_reproduction_{rep}.{ext,(eps|png)}',
    priority: 10
    shell:
        'python3 {ANA_DIR}/plot_bitwise_reproduction.py -bs {{input.nest_spk}} -os {{input.original_spk}} -bmem {{input.nest_mem}} -fn {{output}}'.format(ANA_DIR=ANA_DIR,fig_dir=FIG_DIR)

rule plot_groups:
    output:
        plot_7=expand('{folder}/{{experiment}}/{{rep}}/plot_7.{{ext,(eps|png)}}',folder=FIG_DIR),
        plot_8=expand('{folder}/{{experiment}}/{{rep}}/plot_8.{{ext,(eps|png)}}',folder=FIG_DIR),

    input:
        groups=expand('{folder}/{{experiment}}/{{rep}}/groups.json',folder=NEST_DATA_DIR),
    priority: 2
    run:
        shell("""
        python3 code/analysis/plot_group_statistics.py \
        --groupfile {input.groups}\
        --outfolder figures/{wildcards.experiment}/{wildcards.rep}\
        """)

rule plot_combined_groups:
    output:
        plot_8=expand('{folder}/{{experiment}}/{{experiment}}_combined_groups.{{ext,(eps|png)}}',folder=FIG_DIR),

    input:
        groups=expand('{folder}/{{experiment}}/{rep}/groups.json',folder=NEST_DATA_DIR,rep=NUM_REP),
    priority: 2
    run:
        shell("""
        python3 code/analysis/plot_combined_group_statistics.py \
        -gl {input.groups}\
        -fn {output.plot_8}\
        """)

rule plot_bimodal_gamma:
    output:
        weight=expand('{folder}/{{experiment}}/{{experiment}}_bimodalgamma_weight_delay.{{ext,(eps|png)}}',folder=FIG_DIR),
        groups=expand('{folder}/{{experiment}}/{{experiment}}_bimodalgamma_groups.{{ext,(eps|png)}}',folder=FIG_DIR),

    input:
        connectivity=expand('{folder}/{{experiment}}/{rep}/connectivity.json',folder=NEST_DATA_DIR,rep=NUM_REP),
        spikes=expand('{folder}/{{experiment}}/{rep}/spikes-1001.gdf',folder=NEST_DATA_DIR,rep=NUM_REP),
        groups=expand('{folder}/{{experiment}}/{rep}/groups.json',folder=NEST_DATA_DIR,rep=NUM_REP),

    priority: 2
    run:
        shell("""
        python3 code/analysis/plot_bimodal_gamma.py \
        -cl {input.connectivity}\
        -sl {input.spikes}\
        -gl {input.groups}\
        --group_plot {output.groups}\
        --gamma_plot {output.weight}\

        """)

rule plot_bimodal_gamma_nest:
    output:
        weight=expand('{folder}/{{experiment}}/{{experiment}}_bimodalgamma_weight_delay.{{ext,(eps|png)}}',folder=FIG_DIR),
        groups=expand('{folder}/{{experiment}}/{{experiment}}_bimodalgamma_groups_nest.{{ext,(eps|png)}}',folder=FIG_DIR),

    input:
        connectivity=expand('{folder}/{{experiment}}/{rep}/connectivity.json',folder=NEST_DATA_DIR,rep=NUM_REP),
        spikes=expand('{folder}/{{experiment}}/{rep}/spikes-1001.gdf',folder=NEST_DATA_DIR,rep=NUM_REP),
        groups=expand('{folder}/{{experiment}}/{rep}/groups_nest.json',folder=NEST_DATA_DIR,rep=NUM_REP),

    priority: 2
    run:
        shell("""
        python3 code/analysis/plot_bimodal_gamma.py \
        -cl {input.connectivity}\
        -sl {input.spikes}\
        -gl {input.groups}\
        --group_plot {output.groups}\
        --gamma_plot {output.weight}\

        """)

rule test_weights_and_delay:
    input:
        nest=expand('{folder}/{{experiment}}/{{rep}}/connectivity.json',folder=NEST_DATA_DIR),
    output:
        weight=expand('{folder}/{{experiment}}/{{rep}}/weight_distribution.{{ext,(eps|png)}}',folder=FIG_DIR),
    priority: 10
    shell:
        'python3 {ANA_DIR}/weight_and_delay_distribution.py -c {{input.nest}} -o {{output.weight}} '.format(ANA_DIR=ANA_DIR)

rule plot_dynamics:
    output:
        file=expand('{folder}/{{experiment}}/{{rep}}/dynamic_measures.{{ext,(eps|png)}}',folder=FIG_DIR),
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
        weights=expand('{folder}/{{experiment}}/stdp_window.json'.format(folder=NEST_DATA_DIR), experiment=CONFIG_FILES),
        program='{folder}/plot_stdp_window.py'.format(folder=ANA_DIR),
    output:
        plot="{folder}/stdp_windows.pdf".format(folder=FIG_DIR),
    shell:
        """
        python {input.program} -i {input.weights} -o {output.plot}
        """

