#!/bin/bash

################################
# Script used to generate results in Figure 2 and 3 

# Outputs Monocle CDS object and the main plot files shown in Monocle tutorials at: 
# http://cole-trapnell-lab.github.io/monocle-release/tutorial_data/Olsson_dataset_analysis_final.html
# http://cole-trapnell-lab.github.io/monocle-release/docs/


################################


#Fixed arguments
DIR="/data/salomonis2/LabFiles/Erica-data/Chapter4_Figures_redo/Figure_3/MonocleResults" #Working directory with access to input files
PATHTORNA="$DIR/VanGalen_CON_MFGENES_COUNTS_strip.txt" #Path and filename to the expression matrix 
PATHTOGROUPS="$DIR/groups.VanGalen_CON_MFGENES_labs.txt" #Path and filename to the metadata of cells in expression matrix
#EXPRESSIONFAMILY=("tobit" "negbinomial.size") 
EXPRESSIONFAMILY="negbinomial.size" #expressionFamily to make CDS object
#FULLMODELFORMULASTR=("Cluster" "Groups") 
FULLMODELFORMULASTR="Groups" #Input to differentialGeneTest #Groups are cellHarmony clusters 
#INITIALMETHOD=("PCA" "destiny_diffusionMaps") 
INITIALMETHOD="PCA" #Input to reduceDimension parameter "initial_method"
#REDUCTIONMETHOD=("DDRTree" "SimplePPT" "L1-graph" "SGL-tree") 
REDUCTIONMETHOD="DDRTree"  #Input to reduceDimension parameter "reduction_method" (RGE technique)
FILENAME="Fig3_100_Monocle" #Name of new folder name generated in WD with results 
LOGVALUES=FALSE #Whether the input expression matrix is log transformed (TRUE) or not 
ROOT= "auto" #Number of state that needs to be set as the root ; auto allows Monocle to find its root


#For loop

bsub -M32000 -W 48:00 -J $FILENAME -n 4 -R "span[hosts=1]" -o log/monocle_output_100.out <<EOF
module load R/3.6.1

Rscript CL_RunMonocle.R $DIR $PATHTORNA $PATHTOGROUPS $EXPRESSIONFAMILY $FULLMODELFORMULASTR $INITIALMETHOD $REDUCTIONMETHOD $FILENAME $LOGVALUES $ROOT
EOF


