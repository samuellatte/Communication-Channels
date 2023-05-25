#include <iostream>
#include <stdio.h>
#include <math.h>
#include <complex>

using namespace std;

/* Communication parameters */

// Vertical half-wave dipoles as transmit and receive antennas
constexpr double F = 27 * 1e9;
constexpr double H = 2;
constexpr double EIRPmax = 1; // Watt
constexpr double SnrUE = 5; // dB
constexpr double RXnoiseFig = 15; // dB
constexpr double BW = 200 * 1e6; // Hz
constexpr double maxTau = 2 * 1e-6;

/* Ray-tracing parameters */
constexpr double Resistance = 73;
constexpr double Pi = 3.14159265359;
constexpr double W = 2 * Pi * F;
constexpr double C = 299792458;
constexpr double EpsilonR = 4.5;
constexpr double Beta = W/C;
constexpr double Lambda = C/F;

/* Environment parameters */
constexpr double Temp = 298;
constexpr double kBoltz = 1.379 *1e-23;

