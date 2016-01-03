#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -o /agbs/cluster/oaskaris/im/mine/script_outputs
#$ -j y
#$ -p 10

matlab -nodisplay nodesktop -nosplash -r "addpath('/agbs/cluster/oaskaris/im/mb/'); RunMain('ServerResults1',$1),exit"

