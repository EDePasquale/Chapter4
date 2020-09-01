#!/bin/bash


#Dropout
EXP=("VanGalen_CON_MFGENES" "VanGalen_CON_ALLGENES")
VALUE=(0.75 0.50 0.25)
WD="/data/salomonis2/LabFiles/Erica-data/Chapter4_Figures_redo/Figure_3/Sampling_and_Dropout"

#For loop
for a in "${VALUE[@]}" ; do
	for b in "${EXP[@]}" ; do
		JOBNAME="D_${a}_${b}"
		bsub -M32000 -W 24:00 -J $JOBNAME -n 4 -R "span[hosts=1]" -o log/D_${a}_${b}.out <<EOF
module load R/3.6.1

Rscript CL_Dropout.R ${b}_CPTT_strip.txt ${b}_COUNTS_strip.txt $a $WD
EOF
		
	done
done	
		
	