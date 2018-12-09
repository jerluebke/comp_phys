#include <cstdlib>
#include <cmath>

#include "./simu.h"

static double zuf() { return drand48(); }

void Simu::init()
{
  // set parameters:
  nParticle = 100000;

  // iteration:
  interval = .5;
  nInterval = 100;
  deltaT = 0.05;

  Lx = 3.;
  Nx = 150;
  dx = Lx / Nx;

  double vtherm = 1.;

  // initialize particles
  particles = new Particle[nParticle];
  for(int ip=0; ip<nParticle; ++ip){
    Particle* p = particles + ip;
    p -> x = Lx * ip / nParticle;
    //p -> x += 1e-2 * sin( p->x * 2*M_PI / Lx );
    p -> x = Lx * zuf();
    //    p -> v = 0.02*sin(2.*M_PI* p-> x / Lx);
    //p -> v = 0.3 + 0.8 * (ip % 2);
    // Box-Muller
     double phi = 2*M_PI * zuf();
     double z = zuf();
     double rho = vtherm * sqrt( -log( z ) );
     p -> v = 1. + rho*cos(phi);
     //     p -> v = 0.;
  }

  // allocate grid and fields
  x.allocate(Nx);
  n.allocate(Nx);
  u.allocate(Nx);
  T.allocate(Nx);
  Ex.allocate(Nx);

  // set x-positions of cell centers
  for(int i=0; i<Nx; ++i){
    x(i) = (i + 0.5)*dx;
  }
}

void Simu::finish()
{
  // clear grid and fields
  x.free();
  n.free();
  u.free();
  T.free();
  Ex.free();

  delete[] particles; // delete particles
}

