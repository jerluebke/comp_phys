#ifndef SIMU_SIMU_H
#define SIMU_SIMU_H

#include "./array.h"

struct Particle
{
  double x; // position
  double v; // velocity
};

struct Simu
{
  // Simulation parameters:
  static int nParticle; // # particles
  static double interval; // time interval to integrate between data output
  static int nInterval; // # of time intervals in total
  static double deltaT; // maximum time step to use

  static double Lx; // System length
  static int Nx; // # grid cells in x-direction
  static double dx; // grid spacing

  static Array1 x; // x-positions of cell centers
  static Array1 n; // density in cells
  static Array1 u; // velocity in cells
  static Array1 T; // temperature in cells
  static Array1 Ex; // electric field

  static Particle* particles;

  /* Init data. */
  static void init();

  /** Start and run simulation. */
  static void run();

  /** Clear allocated data. */
  static void finish();

  /** Compute particle moments on grid. */
  static void moments();

  /* Integrate for one time interval. */
  static void doInterval(int seq);

  /* Write particles to file with given sequence number. */
  static void output(int seq);
};

#endif
