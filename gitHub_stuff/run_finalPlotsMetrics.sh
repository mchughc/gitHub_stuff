#!/bin/bash

# Directory containing R script to run
code_dir="/projects/cidr/Ambrosone/sarahcn/"


# R --vanilla < ${code_dir}/R/run_finalPlotsMetrics.R > ${code_dir}/R/run_finalPlotsMetrics_v1.out

####### version to run plot subset
R --vanilla < ${code_dir}/R/run_finalPlotsMetrics_subset.R > ${code_dir}/plots/subsetPlot/run_finalPlotsMetrics_exons_v1.out

# run_finalPlotsMetrics_customVars_v1.out
# run_finalPlotsMetrics_exons_v1.out

