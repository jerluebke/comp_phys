#include <iostream>
#include <cmath>

#include "./simu.h"

static void calcEx()
{
  const int Nx = Simu::Nx;
  double* ex = Simu::Ex.ptr;
  double* n = Simu::n.ptr;
  double dx = Simu::Lx / Nx;
  ex[0] = 0;
  double sumE = 0.;
  for(int i=1; i<Nx; ++i){
    ex[i] = ex[i-1] + dx* ( n[i-1] - 1. );
    sumE += ex[i];
  }

  for(int i=0; i<Nx; ++i){
    ex[i] -= sumE / Nx;
    //ex[i] = 0.;
  }
}

void Simu::moments()
{
  n.fill(0.);
  u.fill(0.);
  T.fill(0.);

  for(int j=0; j<nParticle; ++j){
    double xj = particles[j].x;
    int iCell = lround( xj/dx + 0.5 ) - 1; // cell index of particle center
    if(iCell < 0){
      iCell += Nx;
      xj += Lx;
    }else if(iCell >= Nx){
      iCell -= Nx;
      xj -= Lx;
    }

    // weight for neighbor cell: Distance of partcle center from cell center / dx
    double wNeighb = fabs( xj - x(iCell) ) / dx;
    // right or left neighbor cell with periodicity
    int iNeighb = ( xj > x(iCell) ) ? (iCell+1) % Nx : (iCell+Nx-1) % Nx;

    // add particle velocity moment in neighbor and own cell
    n(iNeighb) += wNeighb;
    n(iCell) += (1.-wNeighb);
    const double v = particles[j].v;
    u(iNeighb) += wNeighb * v;
    u(iCell) += (1.-wNeighb) * v;
    T(iNeighb) += wNeighb * v*v;
    T(iCell) += (1.-wNeighb) * v*v;
  }

  // compute proper moments
  const double factor = ((double)Nx)/nParticle;
  const double klein = 1e-10;
  for(int i=0; i<Nx; ++i){
    n(i) *= factor;
    u(i) *= factor/(n(i) + klein);
    T(i) = factor*T(i)/(n(i) + klein) - u(i)*u(i);
  }
  calcEx();
}

static void drift(Particle* parts, int n, double dt)
{
  // push particles
  for(int j=0; j<n; ++j){
    double newX = parts[j].x + dt*parts[j].v;
    // periodic conditions
    if(newX < 0.){
      newX += Simu::Lx;
    }else if(newX > Simu::Lx){
      newX -= Simu::Lx;
    }
    parts[j].x = newX;
  }
}

static void kick(Particle* parts, int n, double dt)
{
  Simu::moments();

  const double dx = Simu::dx;
  const int Nx = Simu::Nx;
  for(int j=0; j<n; ++j){
    double xj = parts[j].x;

    int iBnd = lround( xj/dx ); // index of left cell boundary closest to particle
    double dist = xj - iBnd*dx; // particle's distance to cell boundary

    // periodic conditions for index
    iBnd = ( iBnd + Nx ) % Nx; // % is modulo-operator in C

    // weight for neighbor boundary: Distance of particle center from cell bndry / dx
    double wNeighb = fabs( dist ) / dx;
    // right or left neighbor bndry with periodicity
    int iNeighb = ( dist > 0. ) ? (iBnd+1) % Nx : (iBnd+Nx-1) % Nx;

    // compute average electric field in particle cloud and advance v
    parts[j].v += dt * ( (1-wNeighb)*Simu::Ex( iBnd ) + wNeighb*Simu::Ex( iNeighb ) );
  }
}

void Simu::doInterval(int i)
{
  std::cout << "Starting interval " << i << " up to t=" << i*interval << std::endl;

  double dt = fmin(deltaT, interval); // time step for regular leapfrogging
  kick(particles, nParticle, -0.5*dt); // sync in: kick v backwards by dt/2
  double toDo = interval;
  while( toDo >= dt ){ // do regular kick-drift sequence with dt
    kick(particles, nParticle, dt);
    drift(particles, nParticle, dt);
    toDo -= dt;
  }
  // Now, there's still a remainder dt > toDo >= 0 for x,
  // and v still lags behind by dt/2
  // sync out of Leapfrog: push v by dt/2.
  // Plus: Do a toDo/2 kick followed by toDo drift, followed by toDo/2 kick
  kick(particles, nParticle, 0.5*dt + 0.5*toDo);
  drift(particles, nParticle, toDo); // last drift
  kick(particles, nParticle, 0.5*toDo); // last k/2
}
