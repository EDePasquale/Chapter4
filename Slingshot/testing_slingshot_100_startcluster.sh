#!/bin/bash


#########################
# Script used to generate Slingshot results in Figure 2 and Figure 4 
# Input: for all fixed and optional arguments taking text files except ENDCLUSTER, provide headers to the vectors/matrix

# Output: 
# 1)DR plot 
# 2)Lineage plot (in black and color) 
# 3)Principal curve (in black and color)
# 4)Combined lineage and curve plot 
# 5)Text file containing lineage order information from slingshot object

#########################

#Fixed arguments
LABELS="VG_D_labs.txt" #Labels file providing cluster numbers (no non numeric values)
PLOT="VG_full_coords.txt" #DR plot coordinates 
WD="/data/salomonis2/LabFiles/Erica-data/Chapter4_Figures_redo/Figure_3/SlingshotResults" #working directory to access input files
FILENAME="Fig3_100_Slingshot_startcluster" #Folder name generated within WD with results

#Optional arguments
STARTCLUSTER="14" 
ENDCLUSTER="endclusters.txt" #text file with no header and each end cluster within quotes
LABELSNAME="NULL" # Use a text file with labels corresponding to the colors and labels file. example:"VG_CHHi_SPRING_Labs_K_Names.txt"
COLS="NULL" # Use a text file with color vector corresponding to plot file. example:"VG_CHHi_SPRING_cols.txt"
PCH=16 #Type/shape of each point on the DR plot
CEX=0.5 #Size of each point on the DR plot
PERC=70 ##percent of transparency in DR plot




#For loop

bsub -M32000 -W 24:00 -J $FILENAME -n 4 -R "span[hosts=1]" -o log/slingshot_output_100.out <<EOF
module load R/3.6.1

Rscript CL_RunSlingshot_startcluster.R $LABELS $PLOT $WD $FILENAME $STARTCLUSTER $ENDCLUSTER $LABELSNAME $COLS $PCH $CEX $PERC
EOF
			
		
	