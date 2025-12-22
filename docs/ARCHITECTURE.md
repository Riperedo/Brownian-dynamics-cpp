# System Architecture

## Overview
This project simulates the Brownian motion of particles in a fluid using the **Ermak-McCammon algorithm**. It is designed with modern C++ (C++17) principles, focusing on modularity, performance, and extensibility.

## Core Components

### 1. `Simulation.hpp` (The Engine)
*   **Role**: Manages the main loop, state, and integration.
*   **Design**: Templated class `Simulation<Dim>`.
    *   **Templates**: `Dim` (2 or 3) allows the compiler to generate optimized code for 2D or 3D vector operations at compile-time, avoiding runtime branching.
    *   **State**: Stores `positions` (wrapped in periodic box) and `realPositions` (unwrapped for MSD).
    *   **Integrator**: Implements `integrate()` which updates positions based on forces and random thermal noise (Gaussian).

### 2. `Potentials.hpp` (The Physics)
*   **Role**: Defines how particles interact.
*   **Design**: Strategy Pattern (Runtime Polymorphism).
    *   **Base Class**: `InteractionPotential` (Abstract). Defines `calculateForceMagnitude(r)` and `calculateEnergy(r)`.
    *   **Derived Classes**: `LennardJonesPotential`, `SquareWellPotential`, `HertzianPotential`, `HardSpherePotential`.
*   **Extensibility**: To add a new potential, simply inherit from `InteractionPotential` and implement the force/energy functions. The `Simulation` class handles the rest via a `std::shared_ptr`.

### 3. `Observables.hpp` (The Analysis)
*   **Role**: Measures physical properties during the simulation.
*   **Components**:
    *   `MeanSquareDisplacement`: Tracks particle diffusion over time.
    *   `MeanSquareDisplacement`: Tracks particle diffusion over time.
    *   `RadialDistributionFunction`: Uses `gsl_histogram` to compute the structure of the fluid $g(r)$. Correctly handles normalization for 2D (disks) and 3D (spheres).
    *   `SelfIntermediateScatteringFunction` ($F_s(q,t)$): Computes dynamic relaxation for a specific wave vector $q$. Uses isotropic averaging ($\frac{\sin(qr)}{qr}$ in 3D).

### 4. `Vector.hpp` (The Math)
*   **Role**: Small, header-only struct for N-dimensional vector algebra.
*   **Design**: Templated `Vector<Dim>`. Supports operator overloading (`+`, `-`, `*`) for natural math syntax (e.g., `forces[i] += fVec`).

## Data Flow
1.  **Initialization**: `main.cpp` parses inputs, configures the `Simulation<Dim>` object, and instantiates the chosen `InteractionPotential`.
2.  **Step Loop**:
    *   `calculateForces()`: Computes pair-wise forces (O(N^2)). Finds minimum image distance. Accumulated in `forces` vector.
        *   **Parallelization**: Uses **OpenMP** to distribute the outer loop across threads.
        *   **Thread Safety**: Updates to force vectors use `#pragma omp atomic` to safely handle Newton's 3rd Law updates (`forces[j] -= f`).
        *   **Observables**: Sampling of $g(r)$ is protected by `#pragma omp critical`.
    *   `integrate()`: Updates positions using Langevin dynamics: $\Delta r = \frac{D F}{kT} \Delta t + \xi$. Enforces Periodic Boundary Conditions (PBC).
    *   `sample()`: Periodically records data to `Observables` (`MSD`, `RDF`, `Self-ISF`).
3.  **Equilibration**: If `EQUIL` > 0, the system runs for specified steps *without sampling*. Observables are reset ($t=0$, $r(0)$) after this phase.
4.  **Output**: Writes `.dat` files to user-specified `OUT_DIR` for analysis.

## Directory Structure
*   `src/`: Main entry point.
*   `include/`: Core logic headers.
*   `inputs/`: Configuration files.
*   `scripts/`: Analysis and automation.
*   `bin/`: Compiled artifacts.
