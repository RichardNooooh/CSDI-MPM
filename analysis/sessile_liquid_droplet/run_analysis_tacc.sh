#!/bin/bash

echo "*****************************************"
echo "*    MPM ANALYSIS BATCH RUNNER (TACC)   *"
echo "*           by Richard Noh              *"
echo "*****************************************"

echo "MPM_ANALYSIS | Initializing parameters"
mesh_file_base="meshes/2D/square_msh"
mesh_file_end=".msh"

mp_per_cell_sqrts=("2" "3")
grid_densities=("64" "128")
JULIA_THREADS=("16" "32")
TACC_RUNTIME=("03:00:00" "06:00:00")
grid_length="1"

ALLOCATION="BCS20003"

dest="output/analysis/sessile_liquid_2d/convergence/"

echo "MPM_ANALYSIS | IO Files/Directories"
echo "MPM_ANALYSIS |     Input Mesh = ${mesh_file_base}X${mesh_file_end}"
echo "MPM_ANALYSIS |     Destination Directory = ${dest}"

echo "MPM_ANALYSIS | TACC Parameters"
echo "MPM_ANALYSIS |     Num Threads for Julia = ${JULIA_THREADS[@]}"
echo "MPM_ANALYSIS |     Requested Runtime     = ${TACC_RUNTIME[@]}"
echo "MPM_ANALYSIS |     Allocation Code       = ${ALLOCATION}"

echo "MPM_ANALYSIS | Parameter Array List"
echo "MPM_ANALYSIS |     MP Per Cell (Sqrt) = ${mp_per_cell_sqrts[@]}"
echo "MPM_ANALYSIS |     Grid Densities     = ${grid_densities[@]}"
echo "MPM_ANALYSIS | Current working directory: $PWD"

function round() {
    printf "%.0f" "$1"
}

function send_batch() {
    local grid_density=$1
    local mp_per_cell_sqrt=$2
    local threads=$3
    local run_time=$4

    local batch_name="batch_ex2_p${mp_per_cell_sqrt}_h${grid_density}.sh"

    local mp_per_cell=$(echo "$mp_per_cell_sqrt * $mp_per_cell_sqrt / 1" | bc -l)
    local num_mp_side=$(echo "$mp_per_cell_sqrt * $grid_density / 1" | bc -l)
    num_mp_side=$(round "$num_mp_side")

    local mesh_file=${mesh_file_base}${num_mp_side}${mesh_file_end}
    local grid_num_cells=$(($grid_density*$grid_length))

    local output_name="p${mp_per_cell_sqrt}_h${grid_density}"

    echo "MPM_ANALYSIS |     Creating batch file for grid_density=${grid_density} and mp_per_cell=${mp_per_cell}"
    echo "MPM_ANALYSIS |         Threads = ${threads}, Runtime = ${runtime}"
    touch $batch_name

    echo "#!/bin/bash" > $batch_name
    echo "#SBATCH -J 2_${output_name}" >> $batch_name
    echo "#SBATCH -o 2_${output_name}.o%j" >> $batch_name
    echo "#SBATCH -e 2_${output_name}.e%j" >> $batch_name
    echo "#SBATCH -p normal" >> $batch_name
    echo "#SBATCH -N 1" >> $batch_name
    echo "#SBATCH -n 1" >> $batch_name
    echo "#SBATCH -t $run_time" >> $batch_name
    echo "#SBATCH --mail-type=all" >> $batch_name
    echo "#SBATCH --mail-user=rjnoh@utexas.edu" >> $batch_name
    echo "" >> $batch_name

    echo "module list" >> $batch_name
    echo "pwd" >> $batch_name
    echo "date" >> $batch_name

    echo "" >> $batch_name
    
    echo "LD_LIBRARY_PATH=\"\" julia -t $threads -- analysis_new/sessile_liquid_droplet/Main.jl \
--grid_num_cells $grid_num_cells $grid_num_cells \
--output_name $output_name/ \
--dest $dest \
--mesh_file $mesh_file" >> $batch_name
    echo "" >> $batch_name
    echo "date" >> $batch_name

    # echo "MPM_ANALYSIS |     Submitting batch file..."
    # sbatch $batch_name
    # rm $batch_name
    # echo "MPM_ANALYSIS |     Removed temporary batch file"
}

echo "MPM_ANALYSIS | Beginning batch run..."
sim_count=1
array_length=${#grid_densities[@]}
for (( i = 0; i < array_length; i++ ))
do
    grid_density=${grid_densities[i]}
    threads=${JULIA_THREADS[i]}
    runtime=${TACC_RUNTIME[i]}
    for mp_per_cell_sqrt in "${mp_per_cell_sqrts[@]}"
    do
        # Create directory for this simulation if it doesn't exist
        echo "MPM_ANALYSIS | Simulation #${sim_count}:"
        send_batch $grid_density $mp_per_cell_sqrt $threads $runtime
        sleep 1
        ((sim_count+=1))
    done
done

echo "MPM_ANALYSIS | Done. :)"
