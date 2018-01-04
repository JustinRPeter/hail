#!/bin/bash
#PBS -N cluster_prediction
#PBS -P ICACS_1001593
#PBS -l walltime=23:00:00,mem=128GB,vmem=128GB
cd $PBS_O_WORKDIR
module load r/3.2.3-gnu
#Rscript /home/u1066578/jpeter/rfiles/era/cluster_lightning_best_no_clusts.R
Rscript cluster_lightning_best_no_clusts.R

