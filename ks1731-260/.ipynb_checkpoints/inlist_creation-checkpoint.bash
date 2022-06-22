#!/bin/bash

neutron_diffusion_LOGS_directory="neutron_diffusion_LOGS"
dStar_directory="/home/justin/dStar"

num_models=$(python3 inlist_creation.py)

# Run all training models.
for i in $(seq 0 $(($num_models - 1)))
    do
        inlist=$(printf $neutron_diffusion_LOGS_directory"/inlist%04d" $i)
        ./run_dStar -D $dStar_directory -I $inlist
    done