# AI Developer & User Guide: Brownian Dynamics

## Purpose
This document is designed to provide context for AI assistants (like Gemini) or new developers. It bridges the gap between the Physics of Brownian Motion and this specific C++ implementation. Use this context to request features or "courses" on the topic.

---

## 1. Physics Crash Course
**Brownian Dynamics (BD)** simplifies Molecular Dynamics (MD) by simulating the solvent (fluid) implicitly.

### The Langevin Equation
The motion of a particle $i$ is described by:
$$ m \ddot{r}_i = -\gamma \dot{r}_i + F_i(r) + R(t) $$
*   $-\gamma \dot{r}_i$: Friction/Drag from the fluid.
*   $F_i(r)$: Conservative forces (e.g., Lennard-Jones interaction).
*   $R(t)$: Random force from collisions with solvent molecules.

### The Algorithm (Ermak-McCammon)
In the overdamped limit (viscosity dominates inertia), we update positions as:
$$ r(t + \Delta t) = r(t) + \frac{D}{k_B T} F(t) \Delta t + \xi $$
Where $\xi$ is a random Gaussian displacement with variance $\langle \xi^2 \rangle = 2 D \Delta t$.

**In this Code (`Simulation.hpp`):**
*   `calculateForces()` computes $F(t)$.
*   `integrate()` applies the update rule.
*   `gsl_ran_gaussian` generates $\xi$.

---

## 2. Code Mapping
Use this map to navigate the project:

| Physics Concept       | Code Location             | Class/Function                                     |
| :-------------------- | :------------------------ | :------------------------------------------------- |
| **Particle Position** | `include/Simulation.hpp`  | `positions` (wrapped), `realPositions` (unwrapped) |
| **Interaction Force** | `include/Potentials.hpp`  | `InteractionPotential::calculateForceMagnitude`    |
| **Random Noise**      | `include/Simulation.hpp`  | `gsl_ran_gaussian` inside `integrate()`            |
| **Structure g(r)**    | `include/Observables.hpp` | `RadialDistributionFunction`                       |
| **Diffusion (MSD)**   | `include/Observables.hpp` | `MeanSquareDisplacement`                           |
| **Self-ISF Fs(q,t)**  | `include/Observables.hpp` | `SelfIntermediateScatteringFunction`               |
| **Equilibration**     | `include/Simulation.hpp`  | `Simulation::run()` (Phase 1 logic)                |

---

## 3. How to Extend (Tutorials)

### Tutorial A: Adding a New Potential
**Goal**: Add a "Yukawa" (screened Coulomb) potential: $V(r) = A \frac{e^{-\kappa r}}{r}$.

1.  **Open `include/Potentials.hpp`**.
2.  **Create a Class**: Inherit from `InteractionPotential`.
    ```cpp
    class YukawaPotential : public InteractionPotential {
    private:
        double A, kappa;
    public:
        YukawaPotential(double temp, double a_val, double k_val) 
            : A(a_val), kappa(k_val) {}
            
        double calculateForceMagnitude(double r) const override {
            // F = -dV/dr
            // Implementation of derivative here...
        }
        
        double calculateEnergy(double r) const override { ... }
        std::string getName() const override { return "Yukawa"; }
        double getCutoff() const override { return 2.5 / kappa; } // Example
    };
    ```
3.  **Register**: Go to `src/main.cpp`. Add an `else if` block in the potential factory section to instantiate your new class when `p.potentialType == "YUKAWA"`.

### Tutorial B: Adding an Observable
**Goal**: Calculate "Radius of Gyration".

> [!WARNING]
> **Thread Safety**: The force calculation loop is parallelized with OpenMP. If you call `sample()` inside `calculateForces()`, you **MUST** protect shared data structures (like histograms) with `#pragma omp critical` or use atomic operations.

1.  **Open `include/Observables.hpp`**.
2.  **Create Class**: `class RadiusOfGyration`.
3.  **Implement**: Add a `sample(std::vector<Vector<Dim>> positions)` method.
4.  **Integrate**: Add a `std::unique_ptr<RadiusOfGyration>` to `Simulation.hpp` and call `sample()` in the loop.

---

## 4. Prompts for Learning
To learn more from an AI using this codebase, try these prompts:
*   *"Explain how `integrate()` implements the Overdamped Langevin equation step-by-step."*
*   *"Why do we need Minimum Image Convention in `calculateForces()`?"*
*   *"Generate a derived class for a Harmonic Trap potential centered at the origin."*
