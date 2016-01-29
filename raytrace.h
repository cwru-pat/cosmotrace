#ifndef RAYTRACE_H
#define RAYTRACE_H

#include "raytrace_macros.h"
#include "raytrace_data.h"
#include <algorithm>
#include <cmath>

namespace cosmo
{

/*
 * Class to perform raytracing through an arbitrary EdS spacetime, but
 * with a couple caveats:
 *  - Angular diameter distance tracing must be done in synchronous gauge
 *  - observers are assumed to be at rest, eg. U^\mu = {1,0,0,0}
 * To use this class
 * for(t = t_i; t <= t_f; t += dt)
 * {
 *   1) position for ray instance is given by RayTrace::getRayX
 *      (or RayTrace::getRayIDX)
 *   2) Determine metric at position. Then:
 *      call `copyInCornerPrimitives` and `interpolatePrimitives`
 *      (or call `setPrimitives`)
 *   3) Calculate derived components, `setDerivedQuantities`
 *   4) Evolve ray instance from t -> t + dt: `evolveRay`
 * }
 */
template <typename RT, typename IT>
class RayTrace
{
  /* Simulation timestep */
    RT sim_dt;

  /* Evolved variables */
    RaytraceData<RT> rd;

  /* "Primitive" variables used to calculate Riemann tensor components */
    // at point where ray is
    RaytracePrimitives<RT> rp;
    // at adjacent cells from a lattice
    struct RaytracePrimitives<RT> corner_rp[2][2][2];

  /* variables used for calculations at a given step, calculated from primitives */
    // Riemann tensor components
    RT R_1001, R_2002, R_3003, R_1212, R_1313, R_2323, R_1002,
       R_1003, R_1012, R_1013, R_1023, R_2003, R_2012, R_2013,
       R_2023, R_3012, R_3013, R_3023, R_1213, R_1223, R_1323;

    // screen tensors
    // indices run over {01, 02, 03, 12, 13, 23}
    RT Sig1[6], Sig2[6];

  public:
    // initialize ray with some parameters
    RayTrace(RT sim_dt_in, RaytraceData<RT> rd_in)
    {
      sim_dt = sim_dt_in;
      rd = rd_in;
      return;
    }

    ~RayTrace()
    {
      return;
    }

IT X2IDX(RT x, RT sim_dx, IT sim_N)
{
  return ((IT) (x/sim_dx)) % sim_N; // periodic; assumes cubic lattice
}

IT getRayIDX(int dir /* x/y/z (= 1, 2, 3) direction */,
             RT sim_dx, IT sim_N)
{
  // Get index immediately before ray position
  return X2IDX(iIDX(dir), sim_dx, sim_N);
}

void setRayX_d_1(RT sim_x_d_1)
{
  rd.x_d[iIDX(1)] = sim_x_d_1;
}

    RT getRayX(int dir)
    {
      return rd.x[iIDX(dir)];
    }

    RaytraceData<RT> getRaytraceData()
    {
      return rd;
    }

    // Set primitives at corners
    void copyInCornerPrimitives(struct RaytracePrimitives<RT> corner_rp_in[2][2][2])
    {
      std::copy(&corner_rp_in[0][0][0], &corner_rp_in[0][0][0] + 2*2*2, &corner_rp[0][0][0]);
    }

    // Directly set primitives (interpolation or equivalent presumed done externally)
    void setPrimitives(RaytracePrimitives<RT> rp_in)
    {
      rp = rp_in;
    }

    // interpolate primitives
    void interpolatePrimitives()
    {
      // interpolate all quantities

      // rho (doesn't have [6] components)
      rp.rho = LINEAR_INTERPOLATION(rho);

      for(int a=0; a<6; a++)
      {
        rp.g[a] = LINEAR_INTERPOLATION(g[a]);
        rp.gi[a] = LINEAR_INTERPOLATION(gi[a]);
        rp.K[a] = LINEAR_INTERPOLATION(K[a]);
        rp.Ricci[a] = LINEAR_INTERPOLATION(Ricci[a]);

        for(int i=0; i<3; ++i)
        {
          rp.G[i][a] = LINEAR_INTERPOLATION(G[i][a]);
          rp.dK[i][a] = LINEAR_INTERPOLATION(dK[i][a]);
          rp.dg[i][a] = LINEAR_INTERPOLATION(dg[i][a]);
          rp.GL[i][a] = LINEAR_INTERPOLATION(GL[i][a]);
        }

        for(int b=0; b<6; ++b)
        {
          rp.ddg[b][a] = LINEAR_INTERPOLATION(ddg[b][a]);
        }
      }

      int a = 0;
    }

    RT linearInterpolation(
      RT C000, RT C001, RT C010, RT C011, /* values at "C"orners, C_{x,y,z} */
      RT C100, RT C101, RT C110, RT C111  /* "binary" order */
      )
    {
      RT C00 = C000*(1.0 - rd.x_d[0]) + C100*rd.x_d[0];
      RT C01 = C001*(1.0 - rd.x_d[0]) + C101*rd.x_d[0];
      RT C10 = C010*(1.0 - rd.x_d[0]) + C110*rd.x_d[0];
      RT C11 = C011*(1.0 - rd.x_d[0]) + C111*rd.x_d[0];

      RT C0 = C00*(1.0 - rd.x_d[1]) + C10*rd.x_d[1];
      RT C1 = C01*(1.0 - rd.x_d[1]) + C11*rd.x_d[1];

      // return "C"
      return C0*(1.0 - rd.x_d[2]) + C1*rd.x_d[2];
    }

    void setDerivedQuantities()
    {
      setScreenTensors();
      setRiemannComponents();
    }

    void setScreenTensors()
    {
      // k^\mu = (E, E*V^i)
      Sig1[0 /*[01]*/] = -0.5*rd.S1[iIDX(1)]*rd.E;
      Sig1[1 /*[02]*/] = -0.5*rd.S1[iIDX(2)]*rd.E;
      Sig1[2 /*[03]*/] = -0.5*rd.S1[iIDX(3)]*rd.E;
      Sig1[3 /*[12]*/] = 0.5*( rd.S1[iIDX(1)]*rd.E*rd.V[iIDX(2)] - rd.S1[iIDX(2)]*rd.E*rd.V[iIDX(1)] );
      Sig1[4 /*[13]*/] = 0.5*( rd.S1[iIDX(1)]*rd.E*rd.V[iIDX(3)] - rd.S1[iIDX(3)]*rd.E*rd.V[iIDX(1)] );
      Sig1[5 /*[23]*/] = 0.5*( rd.S1[iIDX(2)]*rd.E*rd.V[iIDX(3)] - rd.S1[iIDX(3)]*rd.E*rd.V[iIDX(2)] );

      Sig2[0 /*[01]*/] = -0.5*rd.S2[iIDX(1)]*rd.E;
      Sig2[1 /*[02]*/] = -0.5*rd.S2[iIDX(2)]*rd.E;
      Sig2[2 /*[03]*/] = -0.5*rd.S2[iIDX(3)]*rd.E;
      Sig2[3 /*[12]*/] = 0.5*( rd.S2[iIDX(1)]*rd.E*rd.V[iIDX(2)] - rd.S2[iIDX(2)]*rd.E*rd.V[iIDX(1)] );
      Sig2[4 /*[13]*/] = 0.5*( rd.S2[iIDX(1)]*rd.E*rd.V[iIDX(3)] - rd.S2[iIDX(3)]*rd.E*rd.V[iIDX(1)] );
      Sig2[5 /*[23]*/] = 0.5*( rd.S2[iIDX(2)]*rd.E*rd.V[iIDX(3)] - rd.S2[iIDX(3)]*rd.E*rd.V[iIDX(2)] );
    }

    void setRiemannComponents()
    {
      R_1212 = RIEMANN_ILMJ(1, 2, 1, 2);
      R_1313 = RIEMANN_ILMJ(1, 3, 1, 3);
      R_2323 = RIEMANN_ILMJ(2, 3, 2, 3);
      R_1213 = RIEMANN_ILMJ(1, 2, 1, 3);
      R_1223 = RIEMANN_ILMJ(1, 2, 2, 3);
      R_1323 = RIEMANN_ILMJ(1, 3, 2, 3);

      R_1012 = RIEMANN_I0MJ(1, 1, 2);
      R_1013 = RIEMANN_I0MJ(1, 1, 3);
      R_1023 = RIEMANN_I0MJ(1, 2, 3);
      R_2012 = RIEMANN_I0MJ(2, 1, 2);
      R_2013 = RIEMANN_I0MJ(2, 1, 3);
      R_2023 = RIEMANN_I0MJ(2, 2, 3);
      R_3012 = RIEMANN_I0MJ(3, 1, 2);
      R_3013 = RIEMANN_I0MJ(3, 1, 3);
      R_3023 = RIEMANN_I0MJ(3, 2, 3);

      R_1001 = RIEMANN_I00J(1, 1);
      R_1002 = RIEMANN_I00J(1, 2);
      R_1003 = RIEMANN_I00J(1, 3);
      R_2002 = RIEMANN_I00J(2, 2);
      R_2003 = RIEMANN_I00J(2, 3);
      R_3003 = RIEMANN_I00J(3, 3);
    }

    RT RicciLensingScalarSum()
    {
      // w=0 fluid in synchronous gauge *only*!
      return -4.0*RAY_PI*rp.rho*rd.E*rd.E;
    }

    RT WeylLensingScalarSum_Re()
    {
      // "W_11 - W_22"
      return WeylLensingScalarSum(1,1) - WeylLensingScalarSum(2, 2);
    }

    RT WeylLensingScalarSum_Im()
    {
      // symmetrize to help reduce roundoff
      return -0.5*2.0*( WeylLensingScalarSum(1, 2) + WeylLensingScalarSum(2, 1) );
    }

    // evolve ray; simple eulerian integration for now.
    void evolveRay()
    {
      // General raytracing; non-angular-diameter-distance calcs

      // evolve Energy
      rd.E += sim_dt*rd.E*(
          /* K_{ij} * V^{i} * V^{j} */
          rp.K[aIDX(1,1)]*rd.V[iIDX(1)]*rd.V[iIDX(1)] + rp.K[aIDX(2,2)]*rd.V[iIDX(2)]*rd.V[iIDX(2)] + rp.K[aIDX(3,3)]*rd.V[iIDX(3)]*rd.V[iIDX(3)]
          + 2.0*(rp.K[aIDX(1,2)]*rd.V[iIDX(1)]*rd.V[iIDX(2)] + rp.K[aIDX(1,3)]*rd.V[iIDX(1)]*rd.V[iIDX(3)] + rp.K[aIDX(2,3)]*rd.V[iIDX(2)]*rd.V[iIDX(3)])
        );

      // evolve Velocity vector
      for(int i=1; i<=3; ++i)
      {
        rd.x[iIDX(i)] += sim_dt*rd.V[iIDX(i)];

        rd.V[iIDX(i)] += sim_dt*(
            ( rd.V[iIDX(1)]*rp.K[aIDX(1,1)] + rd.V[iIDX(2)]*rp.K[aIDX(2,1)] + rd.V[iIDX(3)]*rp.K[aIDX(3,1)] )
              *( 2.0*rp.gi[aIDX(i,1)] - rd.V[iIDX(i)]*rd.V[iIDX(1)] )
            + ( rd.V[iIDX(1)]*rp.K[aIDX(1,2)] + rd.V[iIDX(2)]*rp.K[aIDX(2,2)] + rd.V[iIDX(3)]*rp.K[aIDX(3,2)] )
              *( 2.0*rp.gi[aIDX(i,2)] - rd.V[iIDX(i)]*rd.V[iIDX(2)] )
            + ( rd.V[iIDX(1)]*rp.K[aIDX(1,3)] + rd.V[iIDX(2)]*rp.K[aIDX(2,3)] + rd.V[iIDX(3)]*rp.K[aIDX(3,3)] )
              *( 2.0*rp.gi[aIDX(i,3)] - rd.V[iIDX(i)]*rd.V[iIDX(3)] )
            - rp.G[iIDX(i)][aIDX(1,1)]*rd.V[iIDX(1)]*rd.V[iIDX(1)]
              - rp.G[iIDX(i)][aIDX(1,2)]*rd.V[iIDX(1)]*rd.V[iIDX(2)]
              - rp.G[iIDX(i)][aIDX(1,3)]*rd.V[iIDX(1)]*rd.V[iIDX(3)]
            - rp.G[iIDX(i)][aIDX(2,1)]*rd.V[iIDX(2)]*rd.V[iIDX(1)]
              - rp.G[iIDX(i)][aIDX(2,2)]*rd.V[iIDX(2)]*rd.V[iIDX(2)]
              - rp.G[iIDX(i)][aIDX(2,3)]*rd.V[iIDX(2)]*rd.V[iIDX(3)]
            - rp.G[iIDX(i)][aIDX(3,1)]*rd.V[iIDX(3)]*rd.V[iIDX(1)]
              - rp.G[iIDX(i)][aIDX(3,2)]*rd.V[iIDX(3)]*rd.V[iIDX(2)]
              - rp.G[iIDX(i)][aIDX(3,3)]*rd.V[iIDX(3)]*rd.V[iIDX(3)]
          );

        // Evolution of angular diameter distance quantities
        // necessary vars for evolving
        RT R_optical = RicciLensingScalarSum();
        RT W_optical_Re = WeylLensingScalarSum_Re();
        RT W_optical_Im = WeylLensingScalarSum_Im();

        rd.D_A += sim_dt*rd.P;
        rd.P += sim_dt*(
            rd.D_A / rd.E / rd.E * (
              R_optical - (rd.Q_Re*rd.Q_Re + rd.Q_Im*rd.Q_Im)/rd.D_A/rd.D_A/rd.D_A/rd.D_A
            ) + rd.P * (
              // K_ij*Vi*Vj
              rp.K[aIDX(1,1)]*rd.V[iIDX(1)]*rd.V[iIDX(1)] + rp.K[aIDX(2,2)]*rd.V[iIDX(2)]*rd.V[iIDX(2)] + rp.K[aIDX(3,3)]*rd.V[iIDX(3)]*rd.V[iIDX(3)]
              + 2.0*(rp.K[aIDX(1,2)]*rd.V[iIDX(1)]*rd.V[iIDX(2)] + rp.K[aIDX(1,3)]*rd.V[iIDX(1)]*rd.V[iIDX(3)] + rp.K[aIDX(2,3)]*rd.V[iIDX(2)]*rd.V[iIDX(3)])
            )
          );

        rd.Q_Re += -sim_dt/rd.E/2.0*rd.D_A*rd.D_A*W_optical_Re;
        rd.Q_Im += sim_dt/rd.E*rd.D_A*rd.D_A*W_optical_Im;

        // partial parallel transport equations - TODO
        // S_A is spatial, _|_ V
        rd.S1[iIDX(i)] += sim_dt*0.0;
        rd.S2[iIDX(i)] += sim_dt*0.0;
      }
    } // evolveRay

};


template <typename RT>
RaytracePrimitives<RT> getFRWRayData(RT a, RT H)
{
  RaytracePrimitives<RT> rp = {0};

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

  rp.rho = 3.0*H*H/8.0/RAY_PI;

  return rp;
}


template <typename RT>
void setFRWRayCorners(RT a, RT H, struct RaytracePrimitives<RT> corner_rp[2][2][2])
{
  RaytracePrimitives<RT> rp = {0};
  rp = getFRWRayData(a, H);

  corner_rp[0][0][0] = rp; corner_rp[0][0][1] = rp; corner_rp[0][1][0] = rp; corner_rp[0][1][1] = rp;
  corner_rp[1][0][0] = rp; corner_rp[1][0][1] = rp; corner_rp[1][1][0] = rp; corner_rp[1][1][1] = rp;

  return;
}


template <typename RT>
RaytracePrimitives<RT> getSinusoidRayData(RT x, RT L, RT eps0)
{
  RaytracePrimitives<RT> rp = {0};

  // f = sin( 2 pi x / L )

  // g_xx = 1 + eps0*f
  rp.g[aIDX(1,1)] = 1.0 + eps0*std::sin(2.0*RAY_PI*x/L);
  rp.g[aIDX(2,2)] = 1.0;
  rp.g[aIDX(3,3)] = 1.0;

  // g^xx = 1 / g_xx
  rp.gi[aIDX(1,1)] = 1.0 / rp.g[aIDX(1,1)];
  rp.gi[aIDX(2,2)] = 1.0;
  rp.gi[aIDX(3,3)] = 1.0;

  // only one non-zero first derivative
  rp.dg[iIDX(1)][aIDX(1,1)] = 2.0*RAY_PI/L*eps0*std::cos(2.0*RAY_PI*x/L);
  // only non-zero second derivative
  rp.ddg[aIDX(1,1)][aIDX(1,1)] = -2.0*RAY_PI/L*2.0*RAY_PI/L*eps0*std::sin(2.0*RAY_PI*x/L);

  // 3-Christoffel symbols
  rp.G[iIDX(1)][aIDX(1,1)] = RAY_PI/L*eps0*std::cos(2.0*RAY_PI*x/L);
  rp.GL[iIDX(1)][aIDX(1,1)] = RAY_PI/L*eps0*std::cos(2.0*RAY_PI*x/L)/(1.0 + eps0*std::sin(2.0*RAY_PI*x/L));

  // need Ricci and rho for D_A integration
  // but not for geodesic integration

  return rp;
}

template <typename RT>
void setSinusoidRayCorners(RT x0, RT x1, RT L, RT eps0, struct RaytracePrimitives<RT> corner_rp[2][2][2])
{
  RaytracePrimitives<RT> rp = {0};

  rp = getSinusoidRayData(x0, L, eps0);
  corner_rp[0][0][0] = rp;
  corner_rp[0][0][1] = rp;
  corner_rp[0][1][0] = rp;
  corner_rp[0][1][1] = rp;

  rp = getSinusoidRayData(x1, L, eps0);
  corner_rp[1][0][0] = rp;
  corner_rp[1][0][1] = rp;
  corner_rp[1][1][0] = rp;
  corner_rp[1][1][1] = rp;

  return;
}


} // namespace cosmo

#endif
