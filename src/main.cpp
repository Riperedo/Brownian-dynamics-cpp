#include "Potentials.hpp"
#include "Simulation.hpp"
#include <fstream>
#include <iostream>
#include <memory> // Required for shared_ptr
#include <omp.h>
#include <string>

// Simple parser for input file
struct Parameters {
  int nPart = 121;
  double phi = 0.5;
  double temp = 1.0;
  double cutoff = 5.0;
  long long steps = 1000000;
  int dim = 3; // Default 3D
  std::string potentialType = "LENNARD_JONES";
  // Extra params for other potentials
  double lambda = 1.5; // for Square Well

  std::string suffix = "";
  std::string outputDir = "scripts";
};

Parameters readInput(const std::string &filename) {
  Parameters p;
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    std::cerr << "Warning: Could not open " << filename << ". Using defaults."
              << std::endl;
    return p;
  }

  std::string line;
  // Simple line-based reader.
  // Original input format is just values on separate lines.
  // We can try to follow that or a smarter key-value.
  // Let's assume the ORIGINAL format for backward compatibility,
  // but extended.
  // Original Order: NPART, PHI, TSTAR, CUT, BDSTEP, SSTAT

  // We'll read the first few lines as standard.
  // Extended lines: Dimension, Potential.

  double val;
  if (infile >> val)
    p.nPart = static_cast<int>(val);
  if (infile >> val)
    p.phi = val;
  if (infile >> val)
    p.temp = val;
  if (infile >> val)
    p.cutoff = val;
  if (infile >> val)
    p.steps = static_cast<long long>(val);
  // Skip SSTAT for now or use it. (Input file has 6th line)
  // if (infile >> val) ...

  // Check for our new extensions (optional keys)
  // But since original file is rigid, let's look for separate config or args.
  // For now, let's hardcode checking if more data exists or defaults.
  // Let's allow command line args to override or a "cpp_config.inp"

  infile.close();
  return p;
}

int main(int argc, char *argv[]) {
  // Defaults
  int nPart = 121;
  double phi = 0.5;
  double temp = 1.0;
  double cutoff = 5.0;
  long long steps = 1000000;
  long long equilSteps = 0; // Default 0 for backwards compat
  int saveFreq = 1000;
  int dim = 3;
  std::string potentialType = "HARD_SPHERE";
  std::string outputDir = "scripts";
  std::string initFile = "";
  double lambda = 1.5;
  double qVal = 7.14; // Default Wave Vector magnitude
  double dt = 0.001;  // Default Time Step

  // Process CLI arguments (Simple "Key=Value" or Flags)
  // Usage: ./brownian_dynamics 3D PHI=0.4 N=500 T=1.5
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "3D")
      dim = 3;
    else if (arg == "2D")
      dim = 2;
    else if (arg == "SQUARE_WELL")
      potentialType = "SQUARE_WELL";
    else if (arg == "HERTZIAN")
      potentialType = "HERTZIAN";
    else if (arg == "LENNARD_JONES" || arg == "LJ")
      potentialType = "LENNARD_JONES";
    else if (arg == "HARD_SPHERE" || arg == "HS")
      potentialType = "HARD_SPHERE";
    else if (arg.find("PHI=") == 0)
      phi = std::stod(arg.substr(4));
    else if (arg.find("N=") == 0)
      nPart = std::stoi(arg.substr(2));
    else if (arg.find("T=") == 0)
      temp = std::stod(arg.substr(2));
    else if (arg.find("CUT=") == 0)
      cutoff = std::stod(arg.substr(4));
    else if (arg.find("STEPS=") == 0)
      steps = std::stoll(arg.substr(6));
    else if (arg.find("OUT_DIR=") == 0)
      outputDir = arg.substr(8);
    else if (arg.find("LAMBDA=") == 0)
      lambda = std::stod(arg.substr(7));
    else if (arg.find("Q=") == 0)
      qVal = std::stod(arg.substr(2));
    else if (arg.find("DT=") == 0)
      dt = std::stod(arg.substr(3));
    else if (arg.find("EQUIL=") == 0)
      equilSteps = std::stoll(arg.substr(6));
    else if (arg.find("SAVE_FREQ=") == 0)
      saveFreq = std::stoi(arg.substr(10));
    else if (arg.find("INIT_FILE=") == 0)
      initFile = arg.substr(10);
  }

  // Auto-generate suffix
  std::string suffix = "d" + std::to_string(dim) + "_" + potentialType +
                       "_phi" + std::to_string(phi);

  std::cout << "Configuration (CLI):" << std::endl;
  std::cout << "  Particles: " << nPart << std::endl;
  std::cout << "  Dimension: " << dim << "D" << std::endl;
  std::cout << "  Potential: " << potentialType << std::endl;
  std::cout << "  Volume Fraction (Phi): " << phi << std::endl;
  std::cout << "  Temperature: " << temp << std::endl;
  std::cout << "  Time Step (dt): " << dt << std::endl;
  std::cout << "  Equilibration: " << equilSteps << std::endl;
  std::cout << "  Steps: " << steps << std::endl;
  std::cout << "  Save Frequency: " << saveFreq << std::endl;
  std::cout << "  Wave Vector (q): " << qVal << std::endl;
  std::cout << "  Output Directory: " << outputDir << std::endl;
  std::cout << "  OpenMP Threads: " << omp_get_max_threads() << std::endl;
  if (!initFile.empty())
    std::cout << "  Init File: " << initFile << std::endl;

  // Factory
  std::shared_ptr<InteractionPotential> pot;
  if (potentialType == "SQUARE_WELL") {
    pot = std::make_shared<SquareWellPotential>(temp, 1.0, lambda);
  } else if (potentialType == "HERTZIAN") {
    pot = std::make_shared<HertzianPotential>(temp, 1.0);
  } else if (potentialType == "HARD_SPHERE") {
    pot = std::make_shared<HardSpherePotential>(temp, 1.0);
  } else {
    pot = std::make_shared<LennardJonesPotential>(temp, 1.0, cutoff);
  }

  // Dispatch simulation based on Dimension
  std::vector<double> qs = {
      qVal}; // We pass the single requested Q. Could extend to list later.

  if (dim == 2) {
    Simulation<2> sim(nPart, phi, temp, cutoff, dt, steps, pot, suffix,
                      outputDir, qs, equilSteps, saveFreq, initFile);
    sim.run();
  } else if (dim == 3) {
    Simulation<3> sim(nPart, phi, temp, cutoff, dt, steps, pot, suffix,
                      outputDir, qs, equilSteps, saveFreq, initFile);
    sim.run();
  } else {
    std::cerr << "Error: Only 2D or 3D supported." << std::endl;
    return 1;
  }
  return 0;
}
