#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <fstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <memory>
#include <omp.h>
#include <sstream>
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

  // Cell Lists
  std::vector<int> cellList;    // Next particle index
  std::vector<int> cellHead;    // Head particle index for each cell
  std::array<int, Dim> gridDim; // Cells in each dimension
  Vector<Dim> cellSize;
  int numCellsTotal;

  int getCellIndex(const Vector<Dim> &pos) {
    int index = 0;
    int stride = 1;
    for (size_t k = 0; k < Dim; ++k) {
      // Wrap position to [0, boxSize) for binning
      double p = pos[k] + boxSize / 2.0;
      // Handle small numerical boundary errors
      if (p < 0.0)
        p += boxSize;
      if (p >= boxSize)
        p -= boxSize;

      int idx = static_cast<int>(p / cellSize[k]);
      if (idx >= gridDim[k])
        idx = gridDim[k] - 1;
      if (idx < 0)
        idx = 0;

      index += idx * stride;
      stride *= gridDim[k];
    }
    return index;
  }

  void updateCellList() {
    // 1. Calculate grid dimensions if box changed or first run
    // Assuming isotropic box for now, but general logic is safer
    // Ensure at least 3 cells OR 1 cell (if box < cutoff? No, box > cutoff
    // implied)
    for (size_t k = 0; k < Dim; ++k) {
      gridDim[k] = static_cast<int>(boxSize / cutOff);
      if (gridDim[k] < 1)
        gridDim[k] = 1;
      cellSize[k] = boxSize / gridDim[k];
    }

    numCellsTotal = 1;
    for (size_t k = 0; k < Dim; ++k)
      numCellsTotal *= gridDim[k];

    cellHead.assign(numCellsTotal, -1);
    cellList.resize(numParticles);

    for (int i = 0; i < numParticles; ++i) {
      int c = getCellIndex(positions[i]);
      cellList[i] = cellHead[c];
      cellHead[c] = i;
    }
  }

public:
  Simulation(int nPart, double phi, double temp, double rc, double timeStep,
             long long steps, std::shared_ptr<InteractionPotential> pot,
             std::string suffix = "", std::string outDir = "scripts",
             std::vector<double> qs = {7.14}, long long eqSteps = 0,
             int saveFreq = 1000, std::string initFile = "")
      : numParticles(nPart), density(0), temperature(temp), cutOff(rc),
        dt(timeStep), boxSize(0), potential(pot), step(0), totalSteps(steps),
        equilSteps(eqSteps), saveInterval(saveFreq), outputDir(outDir),
        outputSuffix(suffix), waveVectors(qs) {

    // Initialize RNG
    gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(NULL));

    // Handle Initialization
    if (!initFile.empty()) {
      loadInitialConfiguration(initFile);
    } else {
      // Standard Lattice Init
      positions.resize(numParticles);
      realPositions.resize(numParticles);
      forces.resize(numParticles);
      initializeLattice();
    }

    // Calculate Box Size based on actual numParticles (which might have changed
    // if loaded from file)
    calculateBoxSize(phi);

    // Apply PBC to ensure loaded positions are in range
    for (int i = 0; i < numParticles; ++i) {
      for (size_t k = 0; k < Dim; ++k) {
        // Primitive wrapping to [-L/2, L/2]
        while (positions[i][k] > boxSize / 2.0)
          positions[i][k] -= boxSize;
        while (positions[i][k] < -boxSize / 2.0)
          positions[i][k] += boxSize;
      }
      realPositions[i] = positions[i]; // Reset real positions to wrapped start
    }

    // Forces resize (just in case)
    forces.resize(numParticles);

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
    // Update Cell Lists
    updateCellList();

    // Reset forces
    for (auto &f : forces)
      for (size_t k = 0; k < Dim; ++k)
        f[k] = 0.0;

// Parallelize over particles
// We use dynamic schedule
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < numParticles; ++i) {
      // Identify my cell
      int myCell = getCellIndex(positions[i]);
      Vector<Dim> r_i = positions[i];

      // Recover explicit grid indices kx, ky, kz from myCell
      // This is needed to iterate neighbors correctly handling PBC
      // (Alternatively, pre-compute neighbor list for each cell, but doing it
      // on fly is cheap enough)
      int idx[3] = {0, 0, 0};
      int temp = myCell;
      idx[0] = temp % gridDim[0];
      temp /= gridDim[0];
      if (Dim >= 2) {
        idx[1] = temp % gridDim[1];
        temp /= gridDim[1];
      }
      if (Dim >= 3) {
        idx[2] = temp % gridDim[2];
      }

      // Loop neighbors (3x3x3 block around my cell)
      // D-dimensional loop? Hard with template.
      // We can limit offsets to {-1, 0, 1} per dimension.

      int kx_start = -1, kx_end = 1;
      int ky_start = (Dim >= 2) ? -1 : 0;
      int ky_end = (Dim >= 2) ? 1 : 0;
      int kz_start = (Dim >= 3) ? -1 : 0;
      int kz_end = (Dim >= 3) ? 1 : 0;

      for (int dz = kz_start; dz <= kz_end; ++dz) {
        for (int dy = ky_start; dy <= ky_end; ++dy) {
          for (int dx = kx_start; dx <= kx_end; ++dx) {

            // Apply PBC to neighbor cell index
            int nx = (idx[0] + dx + gridDim[0]) % gridDim[0];
            int ny = (Dim >= 2) ? (idx[1] + dy + gridDim[1]) % gridDim[1] : 0;
            int nz = (Dim >= 3) ? (idx[2] + dz + gridDim[2]) % gridDim[2] : 0;

            // Reconstruct linear neighbor cell index
            int neighborCell = nx;
            if (Dim >= 2)
              neighborCell += ny * gridDim[0];
            if (Dim >= 3)
              neighborCell += nz * gridDim[0] * gridDim[1];

            // Iterate particles in neighbor cell
            int j = cellHead[neighborCell];
            while (j != -1) {
              // Newton's 3rd Law: only if j > i
              // This prevents double counting and self-interaction
              if (j > i) {
                Vector<Dim> dr = r_i - positions[j];

                // Minimum Image Convention
                for (size_t k = 0; k < Dim; ++k) {
                  dr[k] -= std::round(dr[k] / boxSize) * boxSize;
                }

                double r2 = dr.norm2();
                if (r2 < cutOff * cutOff) {
                  double r = std::sqrt(r2);
                  double fMag = potential->calculateForceMagnitude(r);
                  Vector<Dim> fVec = fMag * (1.0 / r) * dr;

                  // Atomic Updates
                  for (size_t k = 0; k < Dim; ++k) {
#pragma omp atomic
                    forces[i][k] += fVec[k];
#pragma omp atomic
                    forces[j][k] -= fVec[k];
                  }

#pragma omp critical
                  rdf->sample(r);
                }
              }
              j = cellList[j]; // Next particle in cell
            }
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
  void calculateBoxSize(double phi) {
    double PI = 3.14159265358979323846;
    double rho;
    if (Dim == 2) {
      rho = phi * 4.0 / PI;
      boxSize = std::sqrt(numParticles / rho);
    } else {
      // 3D
      rho = phi * 6.0 / PI;
      boxSize = std::cbrt(numParticles / rho);
    }
  }

  void initializeLattice() {
    int side = std::ceil(std::pow(numParticles, 1.0 / Dim));
    double del = boxSize / side; // Use the class member boxSize
    double offset = -boxSize / 2.0 + del / 2.0;

    for (int i = 0; i < numParticles; ++i) {
      int iz = (Dim == 3) ? (i / (side * side)) : 0;
      int rem = (Dim == 3) ? (i % (side * side)) : i;
      int iy = rem / side;
      int ix = rem % side;

      if (Dim >= 1)
        positions[i][0] = offset + ix * del;
      if (Dim >= 2)
        positions[i][1] = offset + iy * del;
      if (Dim >= 3)
        positions[i][2] = offset + iz * del;

      realPositions[i] = positions[i];
    }
  }

  void loadInitialConfiguration(const std::string &filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
      std::cerr << "Error: Could not open init file " << filename
                << ". Exiting." << std::endl;
      exit(1);
    }

    positions.clear();
    realPositions.clear();

    std::string line;
    while (std::getline(infile, line)) {
      if (line.empty())
        continue;
      std::stringstream ss(line);
      int type, label;
      double x, y, z;

      // Format: Type X Y Z Label
      if (!(ss >> type >> x >> y >> z >> label)) {
        continue;
      }

      Vector<Dim> p;
      if (Dim >= 1)
        p[0] = x;
      if (Dim >= 2)
        p[1] = y;
      if (Dim >= 3)
        p[2] = z;

      positions.push_back(p);
      realPositions.push_back(p);
    }

    numParticles = positions.size();
    std::cout << "Loaded " << numParticles << " particles from " << filename
              << std::endl;

    if (numParticles == 0) {
      std::cerr << "Error: No particles loaded from " << filename << std::endl;
      exit(1);
    }
  }
};

#endif // SIMULATION_HPP
