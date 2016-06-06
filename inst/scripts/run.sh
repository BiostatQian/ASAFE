#!/bin/bash

jobID=$RANDOM

echo "jobID is ${jobID}"

echo $1 $2

# n_ind, n_markers
Rscript simulate.R ${jobID} $1 $2

