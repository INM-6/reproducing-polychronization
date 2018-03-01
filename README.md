# Re-Polychronization-Computation-With-Spikes

## Description

This repository is meant for reconstructing the scientific article
"Izhikevich, Eugene M. "Polychronization: computation with spikes." Neural computation 18.2 (2006): 245-282." [1]
using [NEST](http://nest-simulator.org/).

## Usage

To run the original model, the NEST implementation of the model, the analysis, the creation of the figures and the compiling of the manuscript run the following commands:

git clone git@github.com:INM-6/reproducing-polychronization.git

cd reproducing-polychronization

git submodule init

git submodule update --remote

snakemake

## Remarks

In order to reconstruct the network model described in [1], a special implementation of an STDP synapse model is needed.
This model is provided in https://github.com/weidel-p/nest-simulator/tree/feature/stdp_izh_connection and is added 
here as submodule.
When cloning this repository, make sure to run `git submodule init` afterwards.

