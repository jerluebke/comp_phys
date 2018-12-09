#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "./simu.h"

void Simu::output(int seq)
{
  // write particle data: Create output file name from seq
  std::ostringstream oss;
  oss << "data/partData." << std::setw(4) << std::setfill('0') << seq;

 // restrict #particles to write out
  const int maxOut = 1000;
  const int nOut = nParticle > maxOut ? maxOut : nParticle;

  FILE* fp = fopen( oss.str().c_str(), "wb"); // open file
  fprintf(fp, "# np = %d time = %g\n*", nOut, seq*interval);
  for(int ip=0; ip<nOut; ++ip){
    int rp = nParticle * drand48();
    fwrite( & particles[rp].x, sizeof(double), 1, fp);
    fwrite( & particles[rp].v, sizeof(double), 1, fp);
  }
  fclose(fp);

  // write field data, compute moments first
  moments();

  std::ostringstream foss;
  foss << "data/fieldData." << std::setw(4) << std::setfill('0') << seq;

  fp = fopen( foss.str().c_str(), "wb");
  fprintf(fp, "# nx = %d time = %g\n*", Nx, seq*interval);

  // write x, n, u, T, Ex to file
  fwrite( x.ptr, sizeof(double), Nx, fp);
  fwrite( n.ptr, sizeof(double), Nx, fp);
  fwrite( u.ptr, sizeof(double), Nx, fp);
  fwrite( T.ptr, sizeof(double), Nx, fp);
  fwrite( Ex.ptr, sizeof(double), Nx, fp);

  fclose(fp);
}
