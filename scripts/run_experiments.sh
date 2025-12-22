#!/bin/bash

# Ensure we run from the project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT" || exit 1

echo "Running from Project Root: $PROJECT_ROOT"

# Compile first
make

echo "Starting 3D Hard Sphere Sweep..."

for phi in 0.1 0.2 0.3 0.4 0.5
do
    echo "Running Phi = $phi ..."
    ./bin/brownian_dynamics 3D HARD_SPHERE PHI=$phi
done

echo "All simulations complete."
echo "All simulations complete."
echo "Results:"
ls -lh scripts/*.dat
