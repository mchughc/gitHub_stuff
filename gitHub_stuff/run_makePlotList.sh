#!/bin/bash

chr=$1
seg=$2

# Directory containing R script to run
code_dir="/projects/cidr/Ambrosone/sarahcn/R"

R --vanilla --args $chr $seg < ${code_dir}/run_makePlotList.R > \
/projects/cidr/Ambrosone/sarahcn/plots/subsetPlot/variant_lists_logFiles/run_makePlotList_chr${chr}_seg${seg}.out
