#!/usr/bin/env bash
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --time=48:00:00
snakemake --jobs 400\
          --cluster "sbatch -N {cluster.N} -n {cluster.n} --time {cluster.time} --mail-type=FAIL"\
          --cluster-config cluster.json\
          --jobname "{jobid}.{rulename}"\
          --keep-going\
          --rerun-incomplete &>snakelog.log

