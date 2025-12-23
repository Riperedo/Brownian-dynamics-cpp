# Brownian Dynamics (C++)

This project is a modern, object-oriented C++ implementation of Brownian Dynamics simulations, originally ported from Fortran. It supports 2D and 3D systems and a modular catalog of interaction potentials (Lennard-Jones, Square Well, Hertzian).

## Features
*   **Multi-dimensional**: Native support for 2D and 3D simulations.
*   **Modular Potentials**: Easily extensible architecture for new interactions.
*   **High Performance**: Uses GSL (GNU Scientific Library) for state-of-the-art random number generation (Mersenne Twister).
*   **Analysis**: Automatic calculation of Mean Square Displacement (MSD) and Radial Distribution Function $g(r)$.

## Documentation
*   [**Architecture**](docs/ARCHITECTURE.md): Technical details on class structure and design patterns.
*   [**AI Developer Guide**](docs/AI_DEVELOPER_GUIDE.md): **Start Here for AI Assistance**. A guide explaining the physics-code mapping and tutorials for extending the project.

## Project Structure
*   `src/`: Source code (`.cpp`).
*   `include/`: Header files (`.hpp`).
*   `bin/`: Compiled executables.
*   `scripts/`: Analysis and plotting scripts (`gnuplot`).

## Requirements
*   C++17 compliant compiler (`g++`, `clang++`)
*   GNU Scientific Library (GSL)
*   Make

## Compilation
To build the project:
```bash
make
```
This generates the executable `bin/brownian_dynamics`.

### Default Run (2D, Lennard-Jones)
Run the simulation providing parameters via command line (order doesn't matter):
```bash
./bin/brownian_dynamics [3D|2D] [PotentialType] [KEY=VALUE...]
```

**Parameters**:
*   `N=...`: Number of particles (Default: 121)
*   `PHI=...`: Volume Fraction (Default: 0.5)
*   `T=...`: Temperature (Default: 1.0)
*   `DT=...`: Time Step (Default: 0.00001)
*   `EQUIL=...`: Equilibration Steps (Default: 0)
*   `SAVE_FREQ=...`: Output Frequency (Default: 1000)
*   `CUT=...`: Cutoff Radius (Default: 2.5)
*   `STEPS=...`: Simulation Steps (Default: 1000000)
*   `OUT_DIR=...`: Output Directory (Default: "scripts")
*   `INIT_FILE=...`: Path to initial config file (Format: `Type X Y Z Label`)
*   `CELL_LIST=...`: 1 for $O(N)$ (Default), 0 for $O(N^2)$ (Legacy/Long Range)

**Example**:
```bash
./bin/brownian_dynamics 3D HARD_SPHERE N=500 PHI=0.45 CELL_LIST=1 OUT_DIR=./results
```

## Credits
This project is a C++ port and modernization of the original Fortran implementation:
*   [**Brownian-Dynamics (Fortran)**](https://github.com/Alpixels/Brownian-Dynamics) by Alpixels.

## Visualization
Install `gnuplot` to visualize results:
```bash
sudo apt install gnuplot
cd scripts
gnuplot plot_results.gp
```
This creates `results_summary.png` showing the dynamics and structure of the system. Outputs are `MSD.dat` and `RadialDist.dat`.
