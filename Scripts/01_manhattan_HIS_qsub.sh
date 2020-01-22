#!/bin/bash
#PBS -N TopMed_Manhattan_HIS
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR


Rscript /home/ashley/TopMed/Scripts/01_man_HIS_run.R

