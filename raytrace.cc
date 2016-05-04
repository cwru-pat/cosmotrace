#include "raytrace.h"
#include <iostream>
#include <string>
#include <map>  
#include <cmath>
#include <fstream>
#include <sys/stat.h>
#include <iomanip>
#include <limits>
#include <iostream>

typedef double real_t;

int main(int argc, char **argv)
{
  // run time
  real_t t_start = 4.0;
  real_t t_end = 0.0;
  real_t t_ref = 1.0;
  real_t dt = -0.005; // units tbd

  // file for writing
  std::ofstream fout;
  mkdir("output", 0755);
  if(std::ifstream("output/RayTrace.dat"))
  {
    int s=1;
    while(std::ifstream(
        "output/RayTrace." + std::to_string(s) + ".dat"
      ))
      s += 1;
    fout.open("output/RayTrace." + std::to_string(s) + ".dat");
  }
  else
  {
    fout.open("output/RayTrace.dat");
  }

  // Type of spacetime for testing
  std::map<int, std::string> test_type_list = {
    {1, "Standard FRW"},
    {2, "Interpolated FRW"},
    {3, "Sinusoid spacetime"},
    {4, "Interpolated Sinusoid"}
  };
  int test_type = 0;
  std::cout << "\n Enter test simulation type: \n";
  for(auto i=test_type_list.begin(); i!=test_type_list.end(); ++i)
    std::cout << "  " << i->first << ": " << i->second << '\n';
  std::cout << " > ";
  std::cin >> test_type;

  // Initial conditions
  cosmo::RaytraceData<real_t> rd = {0};
    // Direction of propagation (will be normalized during sim.)
    real_t a_start = std::pow( t_start, 2.0/3.0 );
    rd.V[0] = 0.2 / a_start;
    rd.V[1] = 0.6 / a_start;
    rd.V[2] = std::sqrt( 1.0 - rd.V[0]*rd.V[0] - rd.V[1]*rd.V[1] ) / a_start;
    // energy in arb. untis
    rd.E = 1.0;
    // Initial 
    rd.Phi = 1.0;
    rd.ell = 0.0;
    rd.sig_Re = 0.0;
    rd.sig_Im = 0.0;

  // create "ray", initialize with above ICs
  cosmo::RayTrace<real_t, int> * ray;
  ray = new cosmo::RayTrace<real_t, int> (dt, rd);

  // evolve ray
  std::cout << "Running...";
  for(real_t t = t_start; t >= t_end; t += dt)
  {

    // set background spacetime properties according to test_type
    real_t L = 1.0;
    real_t dx = 0.0002;
    real_t eps0 = 0.1;
    real_t x = ray->getRayX(1);
    real_t x0 = dx*std::floor(x/dx);
    real_t x1 = dx*std::ceil(x/dx);
    real_t a = std::pow( t, 2.0/3.0 );
    real_t H = 2.0/3.0/t;

    cosmo::RaytracePrimitives<real_t> rp = {0};
    struct cosmo::RaytracePrimitives<real_t> corner_rp[2][2][2];
    switch (test_type)
    {
      case 4:
        cosmo::setSinusoidRayCorners(x0, x1, L, eps0, corner_rp);
        ray->copyInCornerPrimitives(corner_rp);
        ray->setRayX_d(dx);
        ray->interpolatePrimitives();
        break;
      case 3:
        rp = cosmo::getSinusoidRayData(x, L, eps0);
        ray->setPrimitives(rp);
        break;
      case 2:
        cosmo::setFRWRayCorners(a, H, corner_rp);
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
    // std::cout << "Ray is at X = ("
    //           << ray->getRayX(1) << ", "
    //           << ray->getRayX(2) << ", "
    //           << ray->getRayX(3)
    //           << ") with velocity V = ("
    //           << rd.V[0] << ", "
    //           << rd.V[1] << ", "
    //           << rd.V[2]
    //           << ") and E = "
    //           << rd.E
    //           << "\n";
    // std::cout << "{" << ray->getRayX(1) << "," << rd.V[0] << "},";
    
    fout << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
         << rd.ell << '\t'
         << rd.Phi << '\t'
         << rd.E << '\t'
         << rd.V[iIDX(1)]*rd.V[iIDX(1)] + rd.V[iIDX(2)]*rd.V[iIDX(2)] + rd.V[iIDX(3)]*rd.V[iIDX(3)] << '\t'
         << a << '\t'
         << t << '\t'
         << ray->WeylLensingScalarSum_Re() << '\t'
         << ray->WeylLensingScalarSum_Im() << '\n';

  }
  std::cout << " done.\n";
  fout.close();

  return 0;
}
