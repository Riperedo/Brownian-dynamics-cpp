#ifndef OBSERVABLES_HPP
#define OBSERVABLES_HPP

#include "Vector.hpp"
#include <cmath>
#include <fstream>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>
// #include <memory> // Unused here
#include <string>
#include <vector>

// Forward declaration
template <size_t Dim> struct Particle;

// Mean Square Displacement
template <size_t Dim> class MeanSquareDisplacement {
private:
  std::vector<Vector<Dim>> initialPositions;
  std::vector<double> msdX;
  std::vector<double> msdY;
  std::vector<double> msdZ;
  std::vector<double> msdTotal;
  std::vector<double> timePoints;
  std::vector<int>
      counts; // Number of samples per time point if averaging multiple T0

  // We stick to the simpler logic: single T0 at start?
  // Or multiple T0 origins like original code (CORRELA routine)?
  // Original code: "TAKE A NEW t=0 FOR CORRELATION" every IT0 steps.
  // This is "Multiple Time Origins" to improve statistics. Assumes equilibrium.
  // Implementing full Multiple Time Origin logic is complex but worth it.

  // For simplicity in this first port, let's implement the Simple MSD (from
  // t=0). Original CODE does Multiple Time Origins. We should try to support it
  // or warn. Let's implement Single Origin first for clarity, or a simple
  // version. Wait, the user wants "Viability". I will implement Single Origin
  // (Time = 0 to T) first.

  int numParticles;

public:
  MeanSquareDisplacement(int nPart) : numParticles(nPart) {}

  void setInitialPositions(const std::vector<Vector<Dim>> &positions) {
    initialPositions = positions;
  }

  void sample(double time, const std::vector<Vector<Dim>> &currentPositions) {
    double dx2_sum = 0, dy2_sum = 0, dz2_sum = 0;

    for (int i = 0; i < numParticles; ++i) {
      Vector<Dim> diff = currentPositions[i] - initialPositions[i];
      // Note: For MSD, we usually use UNWRAPPED coordinates.
      // If positions are wrapped to box, MSD will be wrong (bounded).
      // The simulation class needs to provide UNWRAPPED positions or we
      // accumulate displacements. Original code: RX0 is stored. RX(I) is
      // typically wrapped? Let's look at original code INTEGRA: RX(I)=RXN! -
      // ANINT(RXN*IBOXX)*BOXX       !X PERIODIC BOUNDARY CONDITION The boundary
      // condition is COMMENTED OUT in original code INTEGRA? "RX(I)=RXN! -
      // ANINT..." -> The ! makes it a comment. So original code uses UNWRAPPED
      // coordinates? "RX(I)=RXN! ... " -> Wait, if it is commented out,
      // particles fly away? Ah, usually for BD one keeps 'real' position for
      // MSD and 'box' position for Force. I will assume specific implementation
      // details later.

      dx2_sum += diff[0] * diff[0];
      dy2_sum += diff[1] * diff[1];
      if (Dim > 2)
        dz2_sum += diff[2] * diff[2];
    }

    timePoints.push_back(time);
    msdX.push_back(dx2_sum / numParticles);
    msdY.push_back(dy2_sum / numParticles);
    if (Dim > 2)
      msdZ.push_back(dz2_sum / numParticles);
    msdTotal.push_back((dx2_sum + dy2_sum + dz2_sum) / numParticles);
  }

  void save(const std::string &filename) {
    std::ofstream outfile(filename);
    outfile << "# Time MSD_X MSD_Y MSD_Z MSD_Total\n";
    for (size_t i = 0; i < timePoints.size(); ++i) {
      outfile << timePoints[i] << " " << msdX[i] << " " << msdY[i] << " ";
      if (Dim > 2)
        outfile << msdZ[i] << " ";
      outfile << msdTotal[i] << "\n";
    }
    outfile.close();
  }
};

// Radial Distribution Function g(r)
class RadialDistributionFunction {
private:
  gsl_histogram *histogram;
  double binWidth;
  double maxDistance;
  int numBins;
  double boxSize;
  int numParticles;
  int numSamples;
  int dimension;
  double rho; // Number density

public:
  RadialDistributionFunction(double boxL, int nPart, int dim, double dr = 0.01)
      : binWidth(dr), maxDistance(boxL / 2.0),
        numBins(static_cast<int>(maxDistance / dr)), boxSize(boxL),
        numParticles(nPart), numSamples(0), dimension(dim) {

    maxDistance = boxSize / 2.0;
    numBins = static_cast<int>(maxDistance / binWidth);

    histogram = gsl_histogram_alloc(numBins);
    gsl_histogram_set_ranges_uniform(histogram, 0.0, maxDistance);

    // Density = N / V
    double volume =
        (dim == 2) ? (boxSize * boxSize) : (boxSize * boxSize * boxSize);
    rho = numParticles / volume;
  }

  ~RadialDistributionFunction() {
    if (histogram)
      gsl_histogram_free(histogram);
  }

  void sample(double distance) {
    if (distance < maxDistance) {
      gsl_histogram_increment(histogram, distance);
    }
  }

  // Call this once per frame, after looping over all pairs
  void frameComputed() { numSamples++; }

  void reset() {
    gsl_histogram_reset(histogram);
    numSamples = 0;
  }

  void save(const std::string &filename) {
    std::ofstream outfile(filename);
    outfile << "# r g(r)\n";
    double PI = 3.14159265358979323846;

    for (int i = 0; i < numBins; ++i) {
      double count = gsl_histogram_get(histogram, i);
      double r_lower, r_upper;
      gsl_histogram_get_range(histogram, i, &r_lower, &r_upper);
      double r_mid = (r_lower + r_upper) / 2.0;

      // Normalization
      // dN(r) = rho * g(r) * dV(r)
      // g(r) = count / (rho * dV(r) * numSamples * numParticles)
      // Wait, usually pair loop counts each pair twice? Or once?
      // If we loop i=1..N-1, j=i+1..N, we count each pair once.
      // Then count needs factor of 2? Or we normalize by local density relative
      // to 1 particle? Standard: count(r) / (N * numSamples *
      // IdealVolumeElement) IdealVolumeElement: 2D: 2 * pi * r * dr 3D: 4 * pi
      // * r^2 * dr

      double dV;
      if (dimension == 2) {
        dV = 2.0 * PI * r_mid * binWidth;
      } else {
        dV = 4.0 * PI * r_mid * r_mid * binWidth;
      }

      // If we looped unique pairs (i<j), we found 'count' pairs.
      // Total pairs at distance r in full N*N matrix would be 2 * count.
      // Normalization is usually:
      // g(r) = (2 * count) / (N * numSamples * rho * dV)

      double gr = (2.0 * count) / (numParticles * numSamples * rho * dV);

      outfile << r_mid << " " << gr << "\n";
    }
    outfile.close();
  }
};

// Self-Intermediate Scattering Function Fs(q, t)
// Fs(q, t) = < exp(i * q * (r(t) - r(0))) >
// Isotropic average:
// 3D: sin(qr)/qr
// 2D: J0(qr) (Bessel function of first kind order 0)
template <size_t Dim> class SelfIntermediateScatteringFunction {
private:
  std::vector<Vector<Dim>> initialPositions;
  std::vector<double> qValues; // Moduli of q vectors to compute
  // Storage: time -> [val_q1, val_q2, ...]
  std::vector<double> timePoints;
  std::vector<std::vector<double>> isfValues; // isfValues[t][q]
  int numParticles;

public:
  SelfIntermediateScatteringFunction(
      int nPart,
      std::vector<double> qs = {7.14}) // Default q approx 2*pi/sigma (peak)
      : qValues(qs), numParticles(nPart) {}

  void setInitialPositions(const std::vector<Vector<Dim>> &positions) {
    initialPositions = positions;
  }

  void sample(double time, const std::vector<Vector<Dim>> &currentPositions) {
    std::vector<double> currentISF(qValues.size(), 0.0);

    // Loop over particles
    // Note: This loop can be parallelized but accumulation needs care.
    // Given Observables are usually light compared to forces, serial might be
    // fine. But let's protect if called from parallel region (which it is).
    // Actually, caller puts critical section around sample calls.

    for (int i = 0; i < numParticles; ++i) {
      Vector<Dim> diff = currentPositions[i] - initialPositions[i];
      double r = std::sqrt(diff.norm2());

      for (size_t k = 0; k < qValues.size(); ++k) {
        double q = qValues[k];
        double val = 0.0;
        double qr = q * r;

        if (qr < 1e-8) {
          val = 1.0;
        } else if (Dim == 3) {
          val = std::sin(qr) / qr;
        } else {
          val = gsl_sf_bessel_J0(qr); // GSL Bessel J0
        }
        currentISF[k] += val;
      }
    }

    // Average
    for (size_t k = 0; k < qValues.size(); ++k) {
      currentISF[k] /= numParticles;
    }

    timePoints.push_back(time);
    isfValues.push_back(currentISF);
  }

  void save(const std::string &filename) {
    std::ofstream outfile(filename);
    outfile << "# Time";
    for (double q : qValues)
      outfile << " Fs(q=" << q << ")";
    outfile << "\n";

    for (size_t i = 0; i < timePoints.size(); ++i) {
      outfile << timePoints[i];
      for (double val : isfValues[i]) {
        outfile << " " << val;
      }
      outfile << "\n";
    }
    outfile.close();
  }
};

#endif // OBSERVABLES_HPP
