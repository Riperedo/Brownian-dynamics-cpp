#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <memory>
#include <omp.h>
#include <string>
#include <vector>

#include "Observables.hpp"
#include "Potentials.hpp"
#include "Vector.hpp"

template <size_t Dim> class Simulation {
private:
  // System State
  std::vector<Vector<Dim>> positions;     // Wrapped positions (in box)
  std::vector<Vector<Dim>> realPositions; // Unwrapped positions (for MSD)
  std::vector<Vector<Dim>> forces;

  // Parameters
  int numParticles;
  double density; // Volume fraction usually input, but we convert to number
                  // density? Input is "Volumen Fraction" (phi). rho* = phi * 4
                  // / pi in 2D?
  double temperature;
  double cutOff;
  double dt; // Time step (H)
  double boxSize;

  // Engine
  gsl_rng *rng;
  std::shared_ptr<InteractionPotential> potential;

  // Observables
  std::unique_ptr<MeanSquareDisplacement<Dim>> msd;

  std::unique_ptr<SelfIntermediateScatteringFunction<Dim>> isf;
  std::unique_ptr<RadialDistributionFunction> rdf;

  // Counters
  long long step;
  long long totalSteps;
  long long equilSteps;
  int saveInterval;

  void applyBoundaryConditions(Vector<Dim> &pos, Vector<Dim> &realPos,
                               const Vector<Dim> &displacement) {
    // Update wrapped position
    for (size_t k = 0; k < Dim; ++k) {
      pos[k] += displacement[k];

      // Periodic wrap
      if (pos[k] > boxSize / 2.0)
        pos[k] -= boxSize;
      if (pos[k] < -boxSize / 2.0)
        pos[k] += boxSize;
    }

    // Update unwrapped position (just add displacement)
    realPos = realPos + displacement;
  }
  std::string outputDir;
  std::string outputSuffix;
  std::vector<double> waveVectors;

public:
  Simulation(int nPart, double phi, double temp, double rc, double timeStep,
             long long steps, std::shared_ptr<InteractionPotential> pot,
             std::string suffix = "", std::string outDir = "scripts",
             std::vector<double> qs = {7.14}, long long eqSteps = 0,
             int saveFreq = 1000)
      : numParticles(nPart), density(0), temperature(temp), cutOff(rc),
        dt(timeStep), boxSize(0), potential(pot), step(0), totalSteps(steps),
        equilSteps(eqSteps), saveInterval(saveFreq), outputDir(outDir),
        outputSuffix(suffix), waveVectors(qs) {

    // Initialize RNG
    gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(NULL)); // Seed with time or fixed

    // Calculate Box Size
    // 2D: phi = rho * pi * sigma^2 / 4  => rho = 4 * phi / (pi * sigma^2)
    // Box L = sqrt(N / rho)
    // Assuming Sigma = 1.0
    double PI = 3.14159265358979323846;
    double rho;
    if (Dim == 2) {
      rho = phi * 4.0 / PI; // From original code: RHOSTAR=PHI*4.0/PI
      boxSize = std::sqrt(numParticles / rho);
    } else {
      // 3D: phi = rho * pi * sigma^3 / 6
      rho = phi * 6.0 / PI;
      boxSize = std::cbrt(numParticles / rho);
    }

    // Init arrays
    positions.resize(numParticles);
    realPositions.resize(numParticles);
    forces.resize(numParticles);

    // Initialize Lattice Grid
    int side = std::ceil(std::pow(numParticles, 1.0 / Dim));
    double del = boxSize / side;
    double offset = -boxSize / 2.0 + del / 2.0; // Center

    // Simple Grid Init (Recursive or nested loops? Nested is hard for template
    // Dim) We'll just place them on a line/plane/cube roughly Better: simple 1D
    // loop mapping index to coordinates
    for (int i = 0; i < numParticles; ++i) {
      // Primitive mapping, robust enough for fluid equilibration
      int iz = (Dim == 3) ? (i / (side * side)) : 0;
      int rem = (Dim == 3) ? (i % (side * side)) : i;
      int iy = rem / side;
      int ix = rem % side;

      if (Dim >= 1) {
        positions[i][0] = offset + ix * del;
      }
      if (Dim >= 2) {
        positions[i][1] = offset + iy * del;
      }
      if (Dim >= 3) {
        positions[i][2] = offset + iz * del;
      }

      realPositions[i] = positions[i];
    }

    // Setup Observables
    msd = std::make_unique<MeanSquareDisplacement<Dim>>(numParticles);
    msd->setInitialPositions(realPositions);

    // Self ISF
    isf = std::make_unique<SelfIntermediateScatteringFunction<Dim>>(
        numParticles, waveVectors);
    isf->setInitialPositions(realPositions);

    rdf = std::make_unique<RadialDistributionFunction>(boxSize, numParticles,
                                                       Dim, 0.01);
  }

  ~Simulation() { gsl_rng_free(rng); }

  void calculateForces() {
    // Reset forces
    for (auto &f : forces)
      for (size_t k = 0; k < Dim; ++k)
        f[k] = 0.0;

// double potentialEnergy = 0.0; // Unused for now

// Pair loop
// OpenMP Parallelization
// Use dynamic scheduling to balance load (inner loop shrinks)
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < numParticles - 1; ++i) {
      for (int j = i + 1; j < numParticles; ++j) {
        Vector<Dim> dr = positions[i] - positions[j];

        // Minimum Image Convention
        for (size_t k = 0; k < Dim; ++k) {
          dr[k] -= std::round(dr[k] / boxSize) * boxSize;
        }

        double r2 = dr.norm2();
        if (r2 < cutOff * cutOff) {
          double r = std::sqrt(r2);
          double fMag = potential->calculateForceMagnitude(r);

          // Force Vector: F = fMag * (dr / r)
          Vector<Dim> fVec = fMag * (1.0 / r) * dr;

          // Update forces[i] - Thread safe (i is unique per thread iteration)
          // Wait, 'forces[i]' is effectively a reduction variable?
          // No, 'i' is the outer loop index. ONLY THIS THREAD processes 'i'.
          // SO forces[i] write is safe?
          // YES, if we only accumulate fVec onto forces[i].
          // BUT we also update forces[j]. 'j' can be accessed by multiple
          // threads (as 'i' in another thread, or 'j' in another). Actually,
          // 'j' varies. Multiple 'i' threads could pick the same 'j'. So
          // forces[j] needs ATOMIC. forces[i] needs protection? No, i is unique
          // to the thread. BUT wait, forces[i] is updated by OTHER threads when
          // THEY are 'i' and THIS 'i' is their 'j'? Yes. If thread 1 handles
          // i=1, j=5. It updates f[5]. Thread 5 handles i=5, j=... (j > 5). BUT
          // thread 2 handles i=2, j=5. It updates f[5]. So 'forces[i]' implies
          // 'forces[some_j]' elsewhere? No. The loop is "for i". Thread owns
          // 'i'. But when 'i' is 'j' in another thread's loop? Yes. Since we do
          // forces[j] -= fVec, any index can be a 'j'. Therefore, ALL force
          // updates must be atomic to be safe.

          for (size_t k = 0; k < Dim; ++k) {
#pragma omp atomic
            forces[i][k] += fVec[k];

#pragma omp atomic
            forces[j][k] -= fVec[k];
          }

#pragma omp critical
          {
            rdf->sample(r);
          }
        }
      }
    }
    rdf->frameComputed();
  }

  void integrate() {
    // Ermak-McCammon (Brownian Dynamics)
    // r(t+dt) = r(t) + F(t)*dt/gamma + random_disp
    // Here assuming gamma/mobility = 1 or normalized?
    // Original code: RXN=RX(I) + FX(I)*H + RDX
    // So mobility is 1.0 implicitly (Mass=1, Gamma=1 -> D = kT/gamma = T)
    // Wait, original params: TSTAR=1.0.
    // Diffusion Coeff D?
    // SIGGMA=SQRT(2.0*H).
    // Random displacement variance is 2*D*dt.
    // If sqrt(2*H) is sigma_gauss, then 2*D*dt = 2*H -> D=1.
    // D = kT / gamma. If T=1 and D=1, then gamma=1.
    // Force term: F * dt / gamma -> F * H. Correct.

    double sigma_gauss = std::sqrt(2.0 * dt);

    for (int i = 0; i < numParticles; ++i) {
      Vector<Dim> displacement;
      for (size_t k = 0; k < Dim; ++k) {
        double noise =
            gsl_ran_gaussian(rng, sigma_gauss); // Stdev = sigma_gauss
        displacement[k] = forces[i][k] * dt + noise;
        // Note: Original code uses "SIGGMA*GASDEV" where GASDEV returns N(0,1).
        // So noise approx N(0, 2*dt). Correct.
      }
      applyBoundaryConditions(positions[i], realPositions[i], displacement);
    }
  }

  void run() {
    std::cout << "Starting Simulation: " << numParticles << " particles, "
              << (Dim == 3 ? "3D" : "2D") << std::endl;

    std::cout << "Starting Simulation: " << numParticles << " particles, "
              << (Dim == 3 ? "3D" : "2D") << std::endl;

    // saveInterval = 1000; // Now set in constructor

    // Equilibration Phase
    if (equilSteps > 0) {
      std::cout << "Equilibrating for " << equilSteps << " steps..."
                << std::endl;
      for (step = 1; step <= equilSteps; ++step) {
        calculateForces();
        integrate();
      }
      std::cout << "Equilibration complete. Resetting clocks." << std::endl;

      // Reset Observables t=0 state
      msd->setInitialPositions(realPositions);
      isf->setInitialPositions(realPositions);
      rdf->reset();
    }

    std::cout << "Running Production for " << totalSteps << " steps..."
              << std::endl;

    for (step = 1; step <= totalSteps; ++step) {
      calculateForces();
      integrate();

      if (step % saveInterval == 0) {
        // std::cout << "Step " << step << " done." << std::endl;
        msd->sample(step * dt, realPositions);
        isf->sample(step * dt, realPositions);
      }
    }

    msd->save(outputDir + "/MSD" +
              (outputSuffix.empty() ? "" : "_" + outputSuffix) + ".dat");
    isf->save(outputDir + "/SelfISF" +
              (outputSuffix.empty() ? "" : "_" + outputSuffix) + ".dat");
    rdf->save(outputDir + "/RadialDist" +
              (outputSuffix.empty() ? "" : "_" + outputSuffix) + ".dat");

    std::cout << "Simulation Complete. Results saved." << std::endl;
  }
};

#endif // SIMULATION_HPP
