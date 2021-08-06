#!/bin/bash

# Set path
CWD=$(pwd)

# Set host node to the one specified by the user
HOST=b459

# Set the scratch directory
export GAUSS_SCRDIR=/scratch/$USER                                

# SSH into the Blues node and run Gaussian
ssh -n $HOST " module load gaussian/09-e.01                                                       ; 
               mkdir -p $GAUSS_SCRDIR                                                      ; 
               g09 -scrdir=$GAUSS_SCRDIR < $CWD/input.dat > $CWD/output.dat                   "                &

