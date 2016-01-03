#!/bin/bash

for i in {1..6}
do
	qsub -l h_vmem=5G run.sh $i
done

#for i in {14..16}
#do
#	qsub -l h_vmem=25G run.sh $i
#done

#for i in {17..19}
#do
#	qsub -l h_vmem=100G run.sh $i
#done
