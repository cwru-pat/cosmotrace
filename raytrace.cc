
#include "raytrace.h"
#include <iostream>
#include <cmath>

typedef double real_t;

int main(int argc, char **argv)
{

  real_t t_start = 1.0;
  real_t t_end = 10.0;
  real_t dt = 0.04; // units tbd


  // Initial conditions
  cosmo::RaytraceData<real_t> rd = {0};
    // Direction of propagation
    rd.V[0] = 0.2;
    rd.V[1] = 0.4;
    rd.V[2] = 0.8944272;
    // energy in arb. untis
    rd.E = 1.0;
    // screen in x-y dirs
    rd.S1[0] = 1.0;
    rd.S2[1] = 1.0;
    // Initial 
    rd.D_A = 1.0;
    rd.Q_Re = 0.0;
    rd.Q_Im = 0.0;

  // create "ray", initialize with above ICs
  cosmo::RayTrace<real_t, int> * ray;
  ray = new cosmo::RayTrace<real_t, int> (dt, rd);

  std::cout << "Running...";
  for(real_t t = t_start; t <= t_end; t += dt)
  {
    // use an FRW spacetime
      real_t a = std::pow( t, 2.0/3.0 );
      real_t H = 2.0/3.0/t;
      // set primitives directly
        cosmo::RaytracePrimitives<real_t> rp = {0};
        rp = cosmo::getFRWRayData(a, H);
        ray->setPrimitives(rp);
      // Or: set primitives nearby for interpolation
        // struct cosmo::RaytracePrimitives<real_t> corner_rp[2][2][2];
        // cosmo::setFRWRayCorners(a, H, corner_rp);
        // ray->copyInCornerPrimitives(corner_rp);
        // ray->interpolatePrimitives();

    // use a static sinusoid spacetime
      real_t L = 1.0;
      real_t dx = 0.2;
      real_t eps0 = 0.1;
      real_t x = ray->getRayX(1);
      // set primitives directly
        // cosmo::RaytracePrimitives<real_t> rp = {0};
        // rp = cosmo::getSinusoidRayData(x, L, eps0);
        // ray->setPrimitives(rp);
      // or: set primitives nearby for interpolation
        // struct cosmo::RaytracePrimitives<real_t> corner_rp[2][2][2];
        // real_t x0 = dx*std::floor(x/dx);
        // real_t x1 = dx*std::ceil(x/dx);
        // real_t x_d = (x - x0) / dx;
        // cosmo::setSinusoidRayCorners(x0, x1, L, eps0, corner_rp);
        // ray->copyInCornerPrimitives(corner_rp);
        // ray->setRayX_d_1(x_d);
        // ray->interpolatePrimitives();


    ray->setDerivedQuantities();
    ray->evolveRay();

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

    std::cout << "{" << ray->getRayX(1) << "," << rd.V[0] << "},";
  }
  std::cout << " done.\n";

  return 0;
}
