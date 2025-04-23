#!/bin/bash

echo "*****************************************"
echo "*    MPM ANALYSIS BATCH RUNNER (LAB)    *"
echo "*           by Richard Noh              *"
echo "*****************************************"

echo "MPM_ANALYSIS | Initializing parameters"
sim_code_directory="analysis_new/catenoid"
mesh_file_base="meshes/3D/catenoid_msh"
mesh_file_end=".msh"

mp_per_cell_sqrts=("2") # "0.5" "1" "2" "3"
grid_densities=("32") # "4" "6" "8" "12" "16" "24" "32" "48" "64"
grid_length="2"
# consider removing 1/h = 24, 48 from pool.

grid_min_active="1.24" # "0.4" "0.74" "1.24"
grid_max_active="2.26" # "2.8" "2.51" "2.26"

dest="output/analysis/catenoid/convergence/"

JULIA_THREADS="9"

echo "MPM ANALYSIS | IO Files/Directories"
echo "MPM_ANALYSIS |     Input Mesh = ${mesh_file_base}X${mesh_file_end}"
echo "MPM_ANALYSIS |     Destination Directory = ${dest}"

echo "MPM_ANALYSIS | Julia Parameters"
echo "MPM_ANALYSIS |     Num Threads = ${JULIA_THREADS}"

echo "MPM_ANALYSIS | Specified Parameters"
echo "MPM_ANALYSIS |     Grid Active Radius: ${grid_min_active} to ${grid_max_active}"

echo "MPM_ANALYSIS | Parameter Array List"
echo "MPM_ANALYSIS |     SP Per Cell (Sqrt) = ${mp_per_cell_sqrts[@]}"
echo "MPM_ANALYSIS |     Grid Densities     = ${grid_densities[@]}"

function round() {
    printf "%.0f" "$1"
}

function run_simulation() {
    local grid_density=$1
    local mp_per_cell_sqrt=$2
    
    local mp_per_cell=$(echo "$mp_per_cell_sqrt * $mp_per_cell_sqrt / 1" | bc -l)
    local num_mp_side=$(echo "$mp_per_cell_sqrt * $grid_density / 1" | bc -l)
    num_mp_side=$(round "$num_mp_side")

    local mesh_file=${mesh_file_base}${num_mp_side}${mesh_file_end}
    local grid_num_cells=$(($grid_density*$grid_length))

    local output_name="h${grid_density}_mp${mp_per_cell_sqrt}"

    if ! [ -f $mesh_file ]; then
        echo "MPM_ANALYSIS |     Creating missing mesh file for $num_mp_side..."
        python meshers/3D/catenoid.py -m $num_mp_side
    fi

    echo "MPM_ANALYSIS |     Running simulation for grid_density=${grid_density} and mp_per_cell=${mp_per_cell}"
    julia -t $JULIA_THREADS -- ${sim_code_directory}/Main.jl \
        --grid_num_cells $grid_num_cells $grid_num_cells $grid_num_cells \
        --mesh_file $mesh_file \
        --grid_r_min_active $grid_min_active \
        --grid_r_max_active $grid_max_active \
        --output_name $output_name/ \
        --dest $dest
    echo "MPM_ANALYSIS |     Finished this simulation"
}

echo "MPM_ANALYSIS | Beginning batch run..."
sim_count=1
for mp_per_cell_sqrt in "${mp_per_cell_sqrts[@]}"
do
    for grid_density in "${grid_densities[@]}"
    do
        echo "MPM_ANALYSIS | Simulation #${sim_count}:"
        run_simulation $grid_density $mp_per_cell_sqrt
        ((sim_count+=1))
    done
done

echo "MPM_ANALYSIS | Done. :)"
