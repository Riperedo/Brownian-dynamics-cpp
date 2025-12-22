#!/bin/bash

# Ensure we run from the project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT" || exit 1

echo "Starting Production Runs for Hard Spheres..."

# Compile ensuring latest version
make

# Main Loop
for phi in 0.55 0.56
#for phi in 0.1
do
    phi_dir="HardSphere/phi_${phi}"
    echo "Processing Phi = $phi ..."
    
    for run in {6..10}
    #for run in {1..5}
    #for run in 1
    do
        run_dir="${phi_dir}/run_${run}"
        mkdir -p "$run_dir"
        
        echo "  - Run $run/10 -> $run_dir"
        
        # Run Simulation
        # Explicitly set Suffix to empty to get clean filenames "MSD.dat" if logical,
        # or let main generate one. User wanted "saved in each folder".
        # If we just pass output dir, main genertaes "MSD_d3_HS...".
        # User might prefer clean "MSD.dat".
        # We can accept the suffix, or modify main. 
        # But for now, let's keep suffix to avoid overwriting risk if logic fails.
        # Actually user said: "en cada carpeta deberia ir guardado cada archivo dat correspondiente a g(r) y a MSD"
        # Since they are in separate folders, distinct filenames are not strictly necessary, but helpful.
        
        ./bin/brownian_dynamics 3D HARD_SPHERE PHI=$phi OUT_DIR="$run_dir" EQUIL=100000 STEPS=10000000 DT=0.00001 SAVE_FREQ=10 > "$run_dir/log.txt"
    done
done

echo "Production runs complete."
echo "Running averaging script..."
python3 HardSphere/average_results.py
