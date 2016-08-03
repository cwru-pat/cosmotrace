#include "raytrace.h"
#include <iostream>
#include <string>
#include <cstring>
#include <map>  
#include <cmath>
#include <fstream>
#include <sys/stat.h>
#include <iomanip>
#include <limits>
#include <iostream>
#include <zlib.h>
#include <cstdlib>

typedef double real_t;

int main(int argc, char **argv)
{
  // file for writing
  gzFile datafile;
  std::string filename = "output/RayTrace.dat.gz";
  mkdir("output", 0755);
  if(std::ifstream(filename))
  {
    int s=1;
    while(std::ifstream(
        "output/RayTrace." + std::to_string(s) + ".dat.gz"
      ))
      s += 1;
    filename = "output/RayTrace." + std::to_string(s) + ".dat.gz";
    datafile = gzopen(filename.c_str(), "ab");
  }
  else
  {
    datafile = gzopen(filename.c_str(), "ab");
  }

  // Type of spacetime for testing
  std::map<int, std::string> test_type_list = {
    {1, "Standard FRW"},
    {2, "Interpolated FRW"},
    {3, "Sinusoid spacetime"},
    {4, "Interpolated Sinusoid"},
    {5, "Kasner spacetime"},
    {6, "Interpolated Kasner"}
  };
  int test_type = 0;
  if(argc > 1)
  {
    test_type = atoi(argv[1]);

    if(test_type_list.find(test_type) == test_type_list.end() /* not found */)
    {
      std::cout << "Unrecognized spacetime type! Exiting.\n";
      return EXIT_FAILURE;
    }

    std::cout << "\nRunning '" << test_type_list[test_type]
              << "' Simulation.\n";
  }
  else
  {
    std::cout << "\n Enter test simulation type: \n";
    for(auto i=test_type_list.begin(); i!=test_type_list.end(); ++i)
      std::cout << "  " << i->first << ": " << i->second << '\n';
    std::cout << " > ";
    std::cin >> test_type;
  }

  // run time
  real_t t_start = 4.0;
  real_t t_end = 0.0;
  real_t t_ref = 1.0;
  real_t dt = -0.00025; // units tbd
  real_t dx = -10.0*dt;

  // Initial conditions
  cosmo::RaytraceData<real_t> rd = {0};
    // Direction of propagation (will be normalized during sim.)
    rd.V[0] = 1.0/std::sqrt(2.0);
    rd.V[1] = 0.0;
    real_t V2 = rd.V[0]*rd.V[0] + rd.V[1]*rd.V[1];
    if(V2 > 1.0)
      V2 = 1.0;
    rd.V[2] = std::sqrt( 1.0 - V2 );
    // energy in arb. untis
    rd.E = 1.0;
    // Initial beam parameters
    rd.Phi = 1.0;
    rd.ell = 0.0;
    rd.sig_Re = 0.0;
    rd.sig_Im = 0.0;

  // create "ray", initialize with above ICs
  cosmo::RayTrace<real_t, int> * ray;
  ray = new cosmo::RayTrace<real_t, int> (dt, dx, rd);

  // evolve ray
  std::cout << "Running...";
  for(real_t t = t_start; t >= t_end; t += dt)
  {

    // set background spacetime properties according to test_type

    // FRW parameters
    real_t a = std::pow( t, 2.0/3.0 );
    real_t H = 2.0/3.0/t;
    // Sinusoid Spacetime parameters
    real_t L = 0.2;
    real_t eps0 = 0.2;
    // Kasner parameters
    real_t px = 2.0/3.0, py = 2.0/3.0, pz = -1.0/3.0;

    // Ray position
    real_t x = ray->getRayX(1);

    cosmo::RaytracePrimitives<real_t> rp = {0};
    struct cosmo::RaytracePrimitives<real_t> corner_rp[2][2][2];
    switch (test_type)
    {
      case 6:
        cosmo::setKasnerRayCornerPrimitives(px, py, pz, t, corner_rp);
        ray->copyInCornerPrimitives(corner_rp);
        ray->interpolatePrimitives();
        break;
      case 5:
        rp = cosmo::getKasnerRayData(px, py, pz, t);
        ray->setPrimitives(rp);
        break;
      case 4:
        cosmo::setSinusoidRayCornerPrimitives(x, dx, L, eps0, corner_rp);
        ray->copyInCornerPrimitives(corner_rp);
        ray->interpolatePrimitives();
        break;
      case 3:
        rp = cosmo::getSinusoidRayData(x, L, eps0);
        ray->setPrimitives(rp);
        break;
      case 2:
        cosmo::setFRWRayCornerPrimitives(a, H, corner_rp);
        ray->copyInCornerPrimitives(corner_rp);
        ray->interpolatePrimitives();
        break;
      case 1:
      default:
        rp = cosmo::getFRWRayData(a, H);
        ray->setPrimitives(rp);
        break;
    }

    // evolve ray
    ray->setDerivedQuantities();
    ray->evolveRay();


    // print ray information
    rd = ray->getRaytraceData();
    char data_str[35];

    sprintf(data_str, "%.15g\t", (double) t);
    gzwrite(datafile, data_str, strlen(data_str));
    sprintf(data_str, "%.15g\t", (double) a);
    gzwrite(datafile, data_str, strlen(data_str));
    sprintf(data_str, "%.15g\t", (double) (rd.V[iIDX(1)]*rd.V[iIDX(1)] + rd.V[iIDX(2)]*rd.V[iIDX(2)] + rd.V[iIDX(3)]*rd.V[iIDX(3)]) );
    gzwrite(datafile, data_str, strlen(data_str));
    sprintf(data_str, "%.15g\t", (double) rd.E);
    gzwrite(datafile, data_str, strlen(data_str));
    sprintf(data_str, "%.15g\t", (double) rd.ell);
    gzwrite(datafile, data_str, strlen(data_str));
    sprintf(data_str, "%.15g\t", (double) rd.Phi);
    gzwrite(datafile, data_str, strlen(data_str));
    sprintf(data_str, "%.15g\t", (double) rd.sig_Re);
    gzwrite(datafile, data_str, strlen(data_str));
    sprintf(data_str, "%.15g\n", (double) rd.sig_Im);
    gzwrite(datafile, data_str, strlen(data_str));
  }
  std::cout << " done.\n";
  gzclose(datafile);

  return 0;
}
