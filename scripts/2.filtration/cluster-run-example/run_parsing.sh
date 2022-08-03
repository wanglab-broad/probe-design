#!/bin/bash -l
ppath="/stanley/WangLab/Documents/probe/Human_Intron_Probe"
#$ -o /stanley/WangLab/Documents/probe/Human_Intron_Probe/log/qsub_log_o.$JOB_ID.$TASK_ID
#$ -e /stanley/WangLab/Documents/probe/Human_Intron_Probe/log/qsub_log_e.$JOB_ID.$TASK_ID

source "/broad/software/scripts/useuse"
reuse Anaconda3
source activate /stanley/WangLab/envs/probe
now=$(date +"%T")
echo "Current time : $now"
python $ppath/code/picky_parsing.py $ppath

echo "Finished"
now=$(date +"%T")
echo "Current time : $now"
