#!/bin/bash

###################################
# Quantifying Assembly Statistics
###################################

scratch=$1
contigs=$2
report=$3
tool_dir=$4
threads=$5

echo $scratch
cd $scratch
date

source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/

# Process scaffold stats with Quast
# commandline output goes to log file
# 2016.08.18 Updated to use quast-3.2, -T changed to -t option
python $tool_dir"/quast-3.2/quast.py" -t $threads $contigs -o quast_output 

mv quast_output/report.txt $report

echo "Quast Completed"
date
exit


