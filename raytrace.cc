
#include "raytrace.h"
#include <iostream>
#include <cmath>

typedef double real_t;

int main(int argc, char **argv)
{

  real_t t_start = 1.0;
  real_t t_end = 2.0;
  real_t dt = 0.001; // units tbd


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
    // Metric ray is propagating on
    cosmo::RaytracePrimitives<real_t> rp = {0};

    real_t a = std::pow( t, 2.0/3.0 );
    real_t H = 2.0/3.0/t;

    rp.g[aIDX(1,1)] = a*a;
    rp.g[aIDX(2,2)] = a*a;
    rp.g[aIDX(3,3)] = a*a;

    rp.gi[aIDX(1,1)] = 1.0/a/a;
    rp.gi[aIDX(2,2)] = 1.0/a/a;
    rp.gi[aIDX(3,3)] = 1.0/a/a;

    rp.K[aIDX(1,1)] = -H*a*a;
    rp.K[aIDX(2,2)] = -H*a*a;
    rp.K[aIDX(3,3)] = -H*a*a;

    rp.trK = -3.0*H;

    // set primitives directly
    ray->setPrimitives(rp);


    ray->setDerivedQuantities();
    ray->evolveRay();

    rd = ray->getRaytraceData();

    std::cout << "Ray is at X = ("
              << ray->getRayX(1) << ", "
              << ray->getRayX(2) << ", "
              << ray->getRayX(3)
              << ") with velocity V = ("
              << rd.V[0] << ", "
              << rd.V[1] << ", "
              << rd.V[2]
              << ") and E*a = "
              << rd.E*a
              << "\n";
  }
  std::cout << " done.\n";

  return 0;
}
