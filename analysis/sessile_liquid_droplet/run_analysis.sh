#!/bin/bash

echo "*****************************************"
echo "*    MPM ANALYSIS BATCH RUNNER (LAB)    *"
echo "*           by Richard Noh              *"
echo "*****************************************"

echo "MPM_ANALYSIS | Initializing parameters"
sim_code_directory="analysis_new/sessile_liquid_droplet/"
mesh_file_base="meshes/2D/square_msh"
mesh_file_end=".msh"

mp_per_cell_sqrts=("0.5" "1" "2" "3") # "0.5" "1" "2" "3"
grid_densities=("4" "6" "8" "12" "16" "24" "32" "48" "64" "128") # "4" "6" "8" "12" "16" "24" "32" "48" "64" "128"
grid_length="1"

dest="output/analysis/sessile_liquid_2d/convergence/"

JULIA_THREADS="6"

echo "MPM_ANALYSIS | IO Files/Directories"
echo "MPM_ANALYSIS |     Input Mesh = ${mesh_file_base}X${mesh_file_end}"
echo "MPM_ANALYSIS |     Destination Directory = ${dest}"

echo "MPM_ANALYSIS | Julia Parameters"
echo "MPM_ANALYSIS |     Num Threads = ${JULIA_THREADS}"

echo "MPM_ANALYSIS | Parameter Array List"
echo "MPM_ANALYSIS |     MP Per Cell (Sqrt) = ${mp_per_cell_sqrts[@]}"
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
        python meshers/2D/square.py -m $num_mp_side
    fi

    # echo "MPM_ANALYSIS |     Running simulation for grid_density=${grid_density} and mp_per_cell=${mp_per_cell}"
    # julia -t $JULIA_THREADS -- ${sim_code_directory}/Main.jl \
    #     --grid_num_cells $grid_num_cells $grid_num_cells \
    #     --output_name $output_name/ \
    #     --dest $dest \
    #     --mesh_file $mesh_file
    # echo "MPM_ANALYSIS |     Finished this simulation"
}

echo "MPM_ANALYSIS | Beginning batch run..."
sim_count=1

for grid_density in "${grid_densities[@]}"
do
    for mp_per_cell_sqrt in "${mp_per_cell_sqrts[@]}"
    do
        echo "MPM_ANALYSIS | Simulation #${sim_count}:"
        run_simulation $grid_density $mp_per_cell_sqrt
        ((sim_count+=1))
    done
done

echo "MPM_ANALYSIS | Done. :)"
