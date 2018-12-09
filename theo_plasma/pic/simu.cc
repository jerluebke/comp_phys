#include <iostream>
#include <cmath>

#include "./simu.h"

// Definition of static (class-) variables
int Simu::nParticle;
double Simu::interval;
int Simu::nInterval;
double Simu::deltaT;

double Simu::Lx; // System length
int Simu::Nx; // System length
double Simu::dx;

// grid arrays
Array1 Simu::x;
Array1 Simu::n;
Array1 Simu::u;
Array1 Simu::T;
Array1 Simu::Ex;

Particle* Simu::particles;

/* Main entry function. */
int main(int argc, char** argv)
{
  std::cout << "Hello, particles!" << std::endl;

  Simu::init();
  Simu::output(0); // initial output (seq=0)
  Simu::run();
  Simu::finish();

  return 0;
}

void Simu::run()
{
  std::cout << "Time interval=" << interval
            << ", deltaT=" << deltaT
            << ", subSteps:" << ceil(interval/deltaT) << std::endl;

  // iteration loop
  for(int iter = 1; iter <= nInterval; ++iter){
    doInterval(iter);
    output(iter);
  }
}

