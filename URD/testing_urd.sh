#!/bin/bash

##########################
# Script used to generate URD results in Figure 2 and 5 

# Input: all input except KNN is mandatory 

# Output: Generates URD object and plot files of mains steps shown in URD tutorial at: 
# https://github.com/farrellja/URD/blob/master/Analyses/QuickStart/URD-QuickStart-AxialMesoderm.md

# Refer to UnfoldURDTree.R for unfolding the URD tree plot for Figure 5 D 

##########################


#Fixed arguments
EXPRESSION="VanGalen_CON_MFGENES_CPTT_strip.txt" #log transformed expression matrix file
METADATA="groups.VanGalen_CON_MFGENES_labs.txt" #metadata of the cells in the expression matrix
WD="/data/salomonis2/LabFiles/Erica-data/Chapter4_Figures_redo/Figure_3/URDResults" #working directory to access input files
FILENAME="URD_VanGalen_CHHi_100" #Name of new folder generated in WD with results
ROOTCLUSTER="HSC_c13" #Set root cluster name or number from metadata
ENDCLUSTER="end_states.txt" #Text file with no header and each end cluster within quotes
LOGVALUES="TRUE" #Whether input expression matrix is log transformed (TRUE) or not 

#Optional Arguments 
KNN=100 #KNN for computation of diffusion maps



#For loop

bsub -M32000 -W 24:00 -J $FILENAME -n 4 -R "span[hosts=1]" -o log/urd_VanGalen_CHHi_100_output.out <<EOF
module load R/3.6.1

Rscript CL_RunURD.R $EXPRESSION $METADATA $WD $FILENAME $ROOTCLUSTER $ENDCLUSTER $LOGVALUES $KNN
EOF
		