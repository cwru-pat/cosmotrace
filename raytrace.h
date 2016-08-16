#ifndef RAYTRACE_H
#define RAYTRACE_H

#include "raytrace_macros.h"
#include "raytrace_data.h"
#include <algorithm>
#include <cmath>

#include <iomanip>
#include <limits>
#include <iostream>
#include <random>

namespace cosmo
{

/**
 * @brief      Class to perform raytracing through an arbitrary EdS spacetime,
 * but with a couple caveats:
 *  - Angular diameter distance tracing must be done in synchronous gauge
 *  - observers are assumed to be at rest, eg. U^\mu = {1,0,0,0}
 * To use this class
 * for(t = t_i; t <= t_f; t += dt) (but might want to run "in reverse")
 * {
 *   1) position for ray instance is given by RayTrace::getRayX
 *      (or RayTrace::getRayIDX)
 *   2) Determine metric at position. Then:
 *      call `copyInCornerPrimitives` and `interpolatePrimitives`
 *      (or call `setPrimitives`)
 *   3) Calculate derived components, `setDerivedQuantities`
 *   4) Evolve ray instance from t -> t + dt: `evolveRay`
 * }
 *
 * @tparam     RT    Real type
 * @tparam     IT    Index type
 */
template <typename RT, typename IT>
class RayTrace
{
  RT sim_dt, ///< Simulation timestep (can be negative)
     sim_dx; ///< simulation grid spacing

  RaytraceData<RT> rd = {}; ///< Structure storing evolved variables

  /* "Primitive" variables used to calculate Riemann tensor components */
  RaytracePrimitives<RT> rp = {}; ///< Metric primitives at ray location
  struct RaytracePrimitives<RT> corner_rp[2][2][2]; ///< Metric primitives at adjacent gridpoints

  /* variables used for calculations at a given step, calculated from primitives */
  RT R_1001, R_2002, R_3003, R_1212, R_1313, R_2323, R_1002,
     R_1003, R_1012, R_1013, R_1023, R_2003, R_2012, R_2013,
     R_2023, R_3012, R_3013, R_3023, R_1213, R_1223, R_1323;

  VarSig<RT> varSig; ///< particular instance of screen tensor

  public:
    /**
     * @brief      initialize ray with some parameters
     *
     * @param[in]  sim_dt_in  simulation timestep dt
     * @param[in]  rd_in      RaytraceData containing ICs
     */
    RayTrace(RT sim_dt_in, RT sim_dx_in, RaytraceData<RT> rd_in)
    {
      sim_dt = sim_dt_in;
      sim_dx = sim_dx_in;
      rd = rd_in;
      initializeScreenVectors();
      return;
    }

    ~RayTrace()
    {
      return;
    }

    /**
     * @brief      position to index
     * assumes periodic lattice with period L = sim_N*sim_dx
     *
     * @param[in]  x       { parameter_description }
     * @param[in]  sim_dx  The sim dx
     * @param[in]  sim_N   The sim n
     *
     * @return     { description_of_the_return_value }
     */
    IT X2IDX(RT x, RT sim_dx, IT sim_N)
    {
      if(x < 0)
      {
        return sim_N - 1 - ( ((IT) (-1.0*x/sim_dx)) % sim_N );
      }

      return ((IT) (x/sim_dx)) % sim_N;
    }

    IT getRayIDX(int dir /* x/y/z (= 1, 2, 3) direction */,
                 RT sim_dx, IT sim_N)
    {
      // return index immediately before ray position
      return X2IDX(getRayX(dir), sim_dx, sim_N);
    }

    void setRayX_d()
    {
      rd.x_d[iIDX(1)] = std::fmod(getRayX(1) / sim_dx, 1.0);
      rd.x_d[iIDX(2)] = std::fmod(getRayX(2) / sim_dx, 1.0);
      rd.x_d[iIDX(3)] = std::fmod(getRayX(3) / sim_dx, 1.0);
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
      // fractional coordinates of ray position between corners
      setRayX_d();

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
      varSig.SK1[0 /*[01]*/] = -0.5*rd.S1[iIDX(1)]*rd.E;
      varSig.SK1[1 /*[02]*/] = -0.5*rd.S1[iIDX(2)]*rd.E;
      varSig.SK1[2 /*[03]*/] = -0.5*rd.S1[iIDX(3)]*rd.E;
      varSig.SK1[3 /*[12]*/] = 0.5*( rd.S1[iIDX(1)]*rd.E*rd.V[iIDX(2)] - rd.S1[iIDX(2)]*rd.E*rd.V[iIDX(1)] );
      varSig.SK1[4 /*[13]*/] = 0.5*( rd.S1[iIDX(1)]*rd.E*rd.V[iIDX(3)] - rd.S1[iIDX(3)]*rd.E*rd.V[iIDX(1)] );
      varSig.SK1[5 /*[23]*/] = 0.5*( rd.S1[iIDX(2)]*rd.E*rd.V[iIDX(3)] - rd.S1[iIDX(3)]*rd.E*rd.V[iIDX(2)] );

      varSig.SK2[0 /*[01]*/] = -0.5*rd.S2[iIDX(1)]*rd.E;
      varSig.SK2[1 /*[02]*/] = -0.5*rd.S2[iIDX(2)]*rd.E;
      varSig.SK2[2 /*[03]*/] = -0.5*rd.S2[iIDX(3)]*rd.E;
      varSig.SK2[3 /*[12]*/] = 0.5*( rd.S2[iIDX(1)]*rd.E*rd.V[iIDX(2)] - rd.S2[iIDX(2)]*rd.E*rd.V[iIDX(1)] );
      varSig.SK2[4 /*[13]*/] = 0.5*( rd.S2[iIDX(1)]*rd.E*rd.V[iIDX(3)] - rd.S2[iIDX(3)]*rd.E*rd.V[iIDX(1)] );
      varSig.SK2[5 /*[23]*/] = 0.5*( rd.S2[iIDX(2)]*rd.E*rd.V[iIDX(3)] - rd.S2[iIDX(3)]*rd.E*rd.V[iIDX(2)] );
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

    /**
     * @brief      Real part of Weyl optical scalar (difference of weyl
     * lensing scalar summation terms)
     *
     * @return     \f$\mathcal{W}_{11}^{\Sigma} - \mathcal{W}_{22}^{\Sigma}\f$
     */
    RT WeylLensingScalarSum_Re()
    {
      return 0.5 * ( WEYL_LENSING_SCALAR_SUM(1, 1) - WEYL_LENSING_SCALAR_SUM(2, 2) );
    }

    /**
     * @brief      Imaginary part of Weyl optical scalar, (weyl lensing scalar
     * summation terms symmetrized to help reduce roundoff)
     *
     * @return     \f$\frac{1}{2}(2\mathcal{W}_{12}^{\Sigma} + 2\mathcal{W}_{21}^{\Sigma})\f$
     */
    RT WeylLensingScalarSum_Im()
    {
      return - 0.5 * ( WEYL_LENSING_SCALAR_SUM(1, 2) + WEYL_LENSING_SCALAR_SUM(2, 1) );
    }

    /**
     * @brief      Enforce normalization of photon velocity vector
     */
    void normalizeVelocity()
    {
      RT magV = std::sqrt(
            rp.g[aIDX(1,1)]*rd.V[iIDX(1)]*rd.V[iIDX(1)] + rp.g[aIDX(2,2)]*rd.V[iIDX(2)]*rd.V[iIDX(2)] + rp.g[aIDX(3,3)]*rd.V[iIDX(3)]*rd.V[iIDX(3)]
            + 2.0*(rp.g[aIDX(1,2)]*rd.V[iIDX(1)]*rd.V[iIDX(2)] + rp.g[aIDX(1,3)]*rd.V[iIDX(1)]*rd.V[iIDX(3)] + rp.g[aIDX(2,3)]*rd.V[iIDX(2)]*rd.V[iIDX(3)])
        );

      for(int i=1; i<=3; i++)
      {
        rd.V[iIDX(i)] /= magV;
      }
    }

    /**
     * @brief      try to make a good guess for initial screen vectors
     */
    void initializeScreenVectors()
    {
      // S_A is spatial, _|_ V
      RT s1b[3] = {0.0, 1.0, 0.0}; // basis vector for S1
      RT s2b[3] = {1.0, 0.0, 0.0}; // basis vector for S2

      // special case if in x- or y-direction; or x-y plane
      if(std::abs(rd.V[2]) <= 1e-10 && std::abs(rd.V[1]) <= 1e-10)
      {
        s2b[2] = 1.0;
        s2b[0] = 0.0;
      }
      else if(std::abs(rd.V[2]) <= 1e-10 && std::abs(rd.V[0]) <= 1e-10)
      {
        s1b[2] = 1.0;
        s1b[1] = 0.0;
      }
      else
      {
        // if V is in x-y plane, pick vectors out of that plane
        if(std::abs(rd.V[2]) <= 1e-10) { s1b[2] = 1.0; s2b[2] = 1.0; }
        // if V is in x-z plane, pick vectors out of that plane
        if(std::abs(rd.V[1]) <= 1e-10) { s2b[1] = 1.0; }
        // if V is in y-z plane, pick vectors out of that plane
        if(std::abs(rd.V[0]) <= 1e-10) { s1b[0] = 1.0; }
      }

      // store result in S1, S2
      for(int i=0; i<3; i++)
      {
        rd.S1[i] = s1b[i];
        rd.S2[i] = s2b[i];
      }

      return;
    }

    void orthonormalizeScreenVectors()
    {
      RT new_S1[3], new_S2[3];
      RT s1b[3] = {rd.S1[0], rd.S1[1], rd.S1[2]}; // temporary basis vector for S1
      RT s2b[3] = {rd.S2[0], rd.S2[1], rd.S2[2]}; // temporary basis vector for S2

      // subtract out part of s1 along V & normalize
      for(int i=1; i<=3; i++)
      {
        new_S1[iIDX(i)] = s1b[iIDX(i)] - (
            rp.g[aIDX(1,1)]*s1b[iIDX(1)]*rd.V[iIDX(1)] + rp.g[aIDX(2,2)]*s1b[iIDX(2)]*rd.V[iIDX(2)] + rp.g[aIDX(3,3)]*s1b[iIDX(3)]*rd.V[iIDX(3)]
            + rp.g[aIDX(1,2)]*s1b[iIDX(1)]*rd.V[iIDX(2)] + rp.g[aIDX(1,3)]*s1b[iIDX(1)]*rd.V[iIDX(3)] + rp.g[aIDX(2,3)]*s1b[iIDX(2)]*rd.V[iIDX(3)]
            + rp.g[aIDX(1,2)]*s1b[iIDX(2)]*rd.V[iIDX(1)] + rp.g[aIDX(1,3)]*s1b[iIDX(3)]*rd.V[iIDX(1)] + rp.g[aIDX(2,3)]*s1b[iIDX(3)]*rd.V[iIDX(2)]
          )*rd.V[iIDX(i)];
      }
      RT mags1 = std::sqrt(
            rp.g[aIDX(1,1)]*new_S1[iIDX(1)]*new_S1[iIDX(1)] + rp.g[aIDX(2,2)]*new_S1[iIDX(2)]*new_S1[iIDX(2)] + rp.g[aIDX(3,3)]*new_S1[iIDX(3)]*new_S1[iIDX(3)]
            + 2.0*(rp.g[aIDX(1,2)]*new_S1[iIDX(1)]*new_S1[iIDX(2)] + rp.g[aIDX(1,3)]*new_S1[iIDX(1)]*new_S1[iIDX(3)] + rp.g[aIDX(2,3)]*new_S1[iIDX(2)]*new_S1[iIDX(3)])
        );
      for(int i=1; i<=3; i++)
      {
        new_S1[iIDX(i)] /= mags1;
      }

      // subtract out part of s2 along V & s1
      for(int i=1; i<=3; i++)
      {
        new_S2[iIDX(i)] = s2b[iIDX(i)] - (
            rp.g[aIDX(1,1)]*s2b[iIDX(1)]*rd.V[iIDX(1)] + rp.g[aIDX(2,2)]*s2b[iIDX(2)]*rd.V[iIDX(2)] + rp.g[aIDX(3,3)]*s2b[iIDX(3)]*rd.V[iIDX(3)]
            + 2.0*(rp.g[aIDX(1,2)]*s2b[iIDX(1)]*rd.V[iIDX(2)] + rp.g[aIDX(1,3)]*s2b[iIDX(1)]*rd.V[iIDX(3)] + rp.g[aIDX(2,3)]*s2b[iIDX(2)]*rd.V[iIDX(3)])
          )*rd.V[iIDX(i)] - (
            rp.g[aIDX(1,1)]*s2b[iIDX(1)]*new_S1[iIDX(1)] + rp.g[aIDX(2,2)]*s2b[iIDX(2)]*new_S1[iIDX(2)] + rp.g[aIDX(3,3)]*s2b[iIDX(3)]*new_S1[iIDX(3)]
            + 2.0*(rp.g[aIDX(1,2)]*s2b[iIDX(1)]*new_S1[iIDX(2)] + rp.g[aIDX(1,3)]*s2b[iIDX(1)]*new_S1[iIDX(3)] + rp.g[aIDX(2,3)]*s2b[iIDX(2)]*new_S1[iIDX(3)])
          )*new_S1[iIDX(i)];
      }
      RT mags2 = std::sqrt(
            rp.g[aIDX(1,1)]*new_S2[iIDX(1)]*new_S2[iIDX(1)] + rp.g[aIDX(2,2)]*new_S2[iIDX(2)]*new_S2[iIDX(2)] + rp.g[aIDX(3,3)]*new_S2[iIDX(3)]*new_S2[iIDX(3)]
            + 2.0*(rp.g[aIDX(1,2)]*new_S2[iIDX(1)]*new_S2[iIDX(2)] + rp.g[aIDX(1,3)]*new_S2[iIDX(1)]*new_S2[iIDX(3)] + rp.g[aIDX(2,3)]*new_S2[iIDX(2)]*new_S2[iIDX(3)])
        );
      for(int i=1; i<=3; i++)
      {
        new_S2[iIDX(i)] /= mags2;
      }

      // store result in S1, S2
      for(int i=0; i<3; i++)
      {
        rd.S1[i] = new_S1[i];
        rd.S2[i] = new_S2[i];
      }

      return;
    }

    // evolve ray; simple eulerian integration for now.
    void evolveRay()
    {
      // General raytracing; non-angular-diameter-distance calcs
      RT ev_E, ev_x[3], ev_V[3], ev_Phi, ev_ell, ev_sig_Re, ev_sig_Im, ev_S1[3], ev_S2[3];
      normalizeVelocity();
      orthonormalizeScreenVectors();

      // Energy
      ev_E = rd.E*(
          /* K_{ij} * V^{i} * V^{j} */
          rp.K[aIDX(1,1)]*rd.V[iIDX(1)]*rd.V[iIDX(1)] + rp.K[aIDX(2,2)]*rd.V[iIDX(2)]*rd.V[iIDX(2)] + rp.K[aIDX(3,3)]*rd.V[iIDX(3)]*rd.V[iIDX(3)]
          + 2.0*(rp.K[aIDX(1,2)]*rd.V[iIDX(1)]*rd.V[iIDX(2)] + rp.K[aIDX(1,3)]*rd.V[iIDX(1)]*rd.V[iIDX(3)] + rp.K[aIDX(2,3)]*rd.V[iIDX(2)]*rd.V[iIDX(3)])
        );

      // evolve Velocity vector
      for(int i=1; i<=3; ++i)
      {
        ev_x[iIDX(i)] = rd.V[iIDX(i)];

        ev_V[iIDX(i)] = (
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

        // lowered velocity vectors
        RT V_1 = rp.g[aIDX(1,1)]*rd.V[iIDX(1)] + rp.g[aIDX(1,2)]*rd.V[iIDX(2)] + rp.g[aIDX(1,3)]*rd.V[iIDX(3)];
        RT V_2 = rp.g[aIDX(2,1)]*rd.V[iIDX(1)] + rp.g[aIDX(2,2)]*rd.V[iIDX(2)] + rp.g[aIDX(2,3)]*rd.V[iIDX(3)];
        RT V_3 = rp.g[aIDX(3,1)]*rd.V[iIDX(1)] + rp.g[aIDX(3,2)]*rd.V[iIDX(2)] + rp.g[aIDX(3,3)]*rd.V[iIDX(3)];
        // V_j V^k \Gamma^j_{ki} for a free index i
        RT VjVkGjk1 = (
              V_1*rp.G[iIDX(1)][aIDX(1,1)]*rd.V[iIDX(1)] + V_2*rp.G[iIDX(2)][aIDX(1,1)]*rd.V[iIDX(1)] + V_3*rp.G[iIDX(3)][aIDX(1,1)]*rd.V[iIDX(1)]
              + V_1*rp.G[iIDX(1)][aIDX(2,1)]*rd.V[iIDX(2)] + V_2*rp.G[iIDX(2)][aIDX(2,1)]*rd.V[iIDX(2)] + V_3*rp.G[iIDX(3)][aIDX(2,1)]*rd.V[iIDX(2)]
              + V_1*rp.G[iIDX(1)][aIDX(3,1)]*rd.V[iIDX(3)] + V_2*rp.G[iIDX(2)][aIDX(3,1)]*rd.V[iIDX(3)] + V_3*rp.G[iIDX(3)][aIDX(3,1)]*rd.V[iIDX(3)]
            );
        RT VjVkGjk2 = (
              V_1*rp.G[iIDX(1)][aIDX(1,2)]*rd.V[iIDX(1)] + V_2*rp.G[iIDX(2)][aIDX(1,2)]*rd.V[iIDX(1)] + V_3*rp.G[iIDX(3)][aIDX(1,2)]*rd.V[iIDX(1)]
              + V_1*rp.G[iIDX(1)][aIDX(2,2)]*rd.V[iIDX(2)] + V_2*rp.G[iIDX(2)][aIDX(2,2)]*rd.V[iIDX(2)] + V_3*rp.G[iIDX(3)][aIDX(2,2)]*rd.V[iIDX(2)]
              + V_1*rp.G[iIDX(1)][aIDX(3,2)]*rd.V[iIDX(3)] + V_2*rp.G[iIDX(2)][aIDX(3,2)]*rd.V[iIDX(3)] + V_3*rp.G[iIDX(3)][aIDX(3,2)]*rd.V[iIDX(3)]
            );
        RT VjVkGjk3 = (
              V_1*rp.G[iIDX(1)][aIDX(1,3)]*rd.V[iIDX(1)] + V_2*rp.G[iIDX(2)][aIDX(1,3)]*rd.V[iIDX(1)] + V_3*rp.G[iIDX(3)][aIDX(1,3)]*rd.V[iIDX(1)]
              + V_1*rp.G[iIDX(1)][aIDX(2,3)]*rd.V[iIDX(2)] + V_2*rp.G[iIDX(2)][aIDX(2,3)]*rd.V[iIDX(2)] + V_3*rp.G[iIDX(3)][aIDX(2,3)]*rd.V[iIDX(2)]
              + V_1*rp.G[iIDX(1)][aIDX(3,3)]*rd.V[iIDX(3)] + V_2*rp.G[iIDX(2)][aIDX(3,3)]*rd.V[iIDX(3)] + V_3*rp.G[iIDX(3)][aIDX(3,3)]*rd.V[iIDX(3)]
            );

        // screen vectors (re-use transport term from V evolution)
        ev_S1[iIDX(i)] = rd.V[iIDX(i)]*(
            rd.S1[iIDX(1)]*rp.g[aIDX(1,1)]*ev_V[iIDX(1)] + rd.S1[iIDX(2)]*rp.g[aIDX(2,1)]*ev_V[iIDX(1)] + rd.S1[iIDX(3)]*rp.g[aIDX(3,1)]*ev_V[iIDX(1)]
            + rd.S1[iIDX(1)]*rp.g[aIDX(1,2)]*ev_V[iIDX(2)] + rd.S1[iIDX(2)]*rp.g[aIDX(2,2)]*ev_V[iIDX(2)] + rd.S1[iIDX(3)]*rp.g[aIDX(3,2)]*ev_V[iIDX(2)]
            + rd.S1[iIDX(1)]*rp.g[aIDX(1,3)]*ev_V[iIDX(3)] + rd.S1[iIDX(2)]*rp.g[aIDX(2,3)]*ev_V[iIDX(3)] + rd.S1[iIDX(3)]*rp.g[aIDX(3,3)]*ev_V[iIDX(3)]
            -3.0*(
              rd.S1[iIDX(1)]*rp.K[aIDX(1,1)]*rd.V[iIDX(1)] + rd.S1[iIDX(2)]*rp.K[aIDX(2,1)]*rd.V[iIDX(1)] + rd.S1[iIDX(3)]*rp.K[aIDX(3,1)]*rd.V[iIDX(1)]
              + rd.S1[iIDX(1)]*rp.K[aIDX(1,2)]*rd.V[iIDX(2)] + rd.S1[iIDX(2)]*rp.K[aIDX(2,2)]*rd.V[iIDX(2)] + rd.S1[iIDX(3)]*rp.K[aIDX(3,2)]*rd.V[iIDX(2)]
              + rd.S1[iIDX(1)]*rp.K[aIDX(1,3)]*rd.V[iIDX(3)] + rd.S1[iIDX(2)]*rp.K[aIDX(2,3)]*rd.V[iIDX(3)] + rd.S1[iIDX(3)]*rp.K[aIDX(3,3)]*rd.V[iIDX(3)]
            )
            + rd.S1[iIDX(1)]*VjVkGjk1 + rd.S1[iIDX(2)]*VjVkGjk2 + rd.S1[iIDX(3)]*VjVkGjk3
          ) - (
            rp.G[iIDX(i)][aIDX(1,1)]*rd.S1[iIDX(1)]*rd.V[iIDX(1)] + rp.G[iIDX(i)][aIDX(1,2)]*rd.S1[iIDX(1)]*rd.V[iIDX(2)] + rp.G[iIDX(i)][aIDX(1,3)]*rd.S1[iIDX(1)]*rd.V[iIDX(3)]
            + rp.G[iIDX(i)][aIDX(2,1)]*rd.S1[iIDX(2)]*rd.V[iIDX(1)] + rp.G[iIDX(i)][aIDX(2,2)]*rd.S1[iIDX(3)]*rd.V[iIDX(2)] + rp.G[iIDX(i)][aIDX(2,3)]*rd.S1[iIDX(2)]*rd.V[iIDX(3)]
            + rp.G[iIDX(i)][aIDX(3,1)]*rd.S1[iIDX(3)]*rd.V[iIDX(1)] + rp.G[iIDX(i)][aIDX(3,2)]*rd.S1[iIDX(2)]*rd.V[iIDX(2)] + rp.G[iIDX(i)][aIDX(3,3)]*rd.S1[iIDX(3)]*rd.V[iIDX(3)]
          ) + (
            rp.gi[aIDX(i,1)]*rp.K[aIDX(1,1)]*rd.S1[iIDX(1)] + rp.gi[aIDX(i,1)]*rp.K[aIDX(1,2)]*rd.S1[iIDX(2)] + rp.gi[aIDX(i,1)]*rp.K[aIDX(1,3)]*rd.S1[iIDX(3)]
            + rp.gi[aIDX(i,2)]*rp.K[aIDX(2,1)]*rd.S1[iIDX(1)] + rp.gi[aIDX(i,2)]*rp.K[aIDX(2,2)]*rd.S1[iIDX(2)] + rp.gi[aIDX(i,2)]*rp.K[aIDX(2,3)]*rd.S1[iIDX(3)]
            + rp.gi[aIDX(i,3)]*rp.K[aIDX(3,1)]*rd.S1[iIDX(1)] + rp.gi[aIDX(i,3)]*rp.K[aIDX(3,2)]*rd.S1[iIDX(2)] + rp.gi[aIDX(i,3)]*rp.K[aIDX(3,3)]*rd.S1[iIDX(3)]
          );
        ev_S2[iIDX(i)] = rd.V[iIDX(i)]*(
            rd.S2[iIDX(1)]*rp.g[aIDX(1,1)]*ev_V[iIDX(1)] + rd.S2[iIDX(2)]*rp.g[aIDX(2,1)]*ev_V[iIDX(1)] + rd.S2[iIDX(3)]*rp.g[aIDX(3,1)]*ev_V[iIDX(1)]
            + rd.S2[iIDX(1)]*rp.g[aIDX(1,2)]*ev_V[iIDX(2)] + rd.S2[iIDX(2)]*rp.g[aIDX(2,2)]*ev_V[iIDX(2)] + rd.S2[iIDX(3)]*rp.g[aIDX(3,2)]*ev_V[iIDX(2)]
            + rd.S2[iIDX(1)]*rp.g[aIDX(1,3)]*ev_V[iIDX(3)] + rd.S2[iIDX(2)]*rp.g[aIDX(2,3)]*ev_V[iIDX(3)] + rd.S2[iIDX(3)]*rp.g[aIDX(3,3)]*ev_V[iIDX(3)]
            -3.0*(
              rd.S2[iIDX(1)]*rp.K[aIDX(1,1)]*rd.V[iIDX(1)] + rd.S2[iIDX(2)]*rp.K[aIDX(2,1)]*rd.V[iIDX(1)] + rd.S2[iIDX(3)]*rp.K[aIDX(3,1)]*rd.V[iIDX(1)]
              + rd.S2[iIDX(1)]*rp.K[aIDX(1,2)]*rd.V[iIDX(2)] + rd.S2[iIDX(2)]*rp.K[aIDX(2,2)]*rd.V[iIDX(2)] + rd.S2[iIDX(3)]*rp.K[aIDX(3,2)]*rd.V[iIDX(2)]
              + rd.S2[iIDX(1)]*rp.K[aIDX(1,3)]*rd.V[iIDX(3)] + rd.S2[iIDX(2)]*rp.K[aIDX(2,3)]*rd.V[iIDX(3)] + rd.S2[iIDX(3)]*rp.K[aIDX(3,3)]*rd.V[iIDX(3)]
            )
            + rd.S2[iIDX(1)]*VjVkGjk1 + rd.S2[iIDX(2)]*VjVkGjk2 + rd.S2[iIDX(3)]*VjVkGjk3
          ) - (
            rp.G[iIDX(i)][aIDX(1,1)]*rd.S2[iIDX(1)]*rd.V[iIDX(1)] + rp.G[iIDX(i)][aIDX(1,2)]*rd.S2[iIDX(1)]*rd.V[iIDX(2)] + rp.G[iIDX(i)][aIDX(1,3)]*rd.S2[iIDX(1)]*rd.V[iIDX(3)]
            + rp.G[iIDX(i)][aIDX(2,1)]*rd.S2[iIDX(2)]*rd.V[iIDX(1)] + rp.G[iIDX(i)][aIDX(2,2)]*rd.S2[iIDX(3)]*rd.V[iIDX(2)] + rp.G[iIDX(i)][aIDX(2,3)]*rd.S2[iIDX(2)]*rd.V[iIDX(3)]
            + rp.G[iIDX(i)][aIDX(3,1)]*rd.S2[iIDX(3)]*rd.V[iIDX(1)] + rp.G[iIDX(i)][aIDX(3,2)]*rd.S2[iIDX(2)]*rd.V[iIDX(2)] + rp.G[iIDX(i)][aIDX(3,3)]*rd.S2[iIDX(3)]*rd.V[iIDX(3)]
          ) + (
            rp.gi[aIDX(i,1)]*rp.K[aIDX(1,1)]*rd.S2[iIDX(1)] + rp.gi[aIDX(i,1)]*rp.K[aIDX(1,2)]*rd.S2[iIDX(2)] + rp.gi[aIDX(i,1)]*rp.K[aIDX(1,3)]*rd.S2[iIDX(3)]
            + rp.gi[aIDX(i,2)]*rp.K[aIDX(2,1)]*rd.S2[iIDX(1)] + rp.gi[aIDX(i,2)]*rp.K[aIDX(2,2)]*rd.S2[iIDX(2)] + rp.gi[aIDX(i,2)]*rp.K[aIDX(2,3)]*rd.S2[iIDX(3)]
            + rp.gi[aIDX(i,3)]*rp.K[aIDX(3,1)]*rd.S2[iIDX(1)] + rp.gi[aIDX(i,3)]*rp.K[aIDX(3,2)]*rd.S2[iIDX(2)] + rp.gi[aIDX(i,3)]*rp.K[aIDX(3,3)]*rd.S2[iIDX(3)]
          );
      }

      if(rd.ell == 0)
      {
        ev_Phi = 0;
      }
      else
      {
        RT R_optical = RicciLensingScalarSum();
        ev_Phi = 1.0/rd.E*rd.ell*(
            R_optical - ( rd.sig_Re*rd.sig_Re + rd.sig_Im*rd.sig_Im )/rd.ell/rd.ell/rd.ell/rd.ell
          );
      }

      ev_ell = 1.0/rd.E*rd.Phi;
      
      RT W_optical_Re = WeylLensingScalarSum_Re();
      ev_sig_Re = 1.0/rd.E*rd.ell*rd.ell*W_optical_Re;

      RT W_optical_Im = WeylLensingScalarSum_Im();
      ev_sig_Im = 1.0/rd.E*rd.ell*rd.ell*W_optical_Im;

      // compute evolved values (simple Eulerian integration for now...)
      rd.E += sim_dt*ev_E;
      rd.Phi += sim_dt*ev_Phi;
      rd.ell += sim_dt*ev_ell;
      rd.sig_Re += sim_dt*ev_sig_Re;
      rd.sig_Im += sim_dt*ev_sig_Im;
      for(int i=0; i<3; i++)
      {
        rd.x[i] += sim_dt*ev_x[i];
        rd.V[i] += sim_dt*ev_V[i];
        rd.S1[i] += sim_dt*ev_S1[i];
        rd.S2[i] += sim_dt*ev_S2[i];
      }

      rd.rho = rp.rho;
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
void setFRWRayCornerPrimitives(RT a, RT H, struct RaytracePrimitives<RT> corner_rp[2][2][2])
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
void setSinusoidRayCornerPrimitives(RT x, RT dx, RT L, RT eps0, struct RaytracePrimitives<RT> corner_rp[2][2][2])
{
  RaytracePrimitives<RT> rp = {0};

  RT x0 = dx*std::floor(x/dx);
  RT x1 = dx*std::ceil(x/dx);

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


template <typename RT>
RaytracePrimitives<RT> getKasnerRayData(RT px, RT py, RT pz, RT t)
{
  RaytracePrimitives<RT> rp = {0};

  RT ax = -std::pow( t, px );
  RT ay = -std::pow( t, py );
  RT az = -std::pow( t, pz );
  RT Kx = -ax*px*std::pow( t, px-1.0 );
  RT Ky = -ay*py*std::pow( t, py-1.0 );
  RT Kz = -az*pz*std::pow( t, pz-1.0 );

  rp.g[aIDX(1,1)] = ax*ax;
  rp.g[aIDX(2,2)] = ay*ay;
  rp.g[aIDX(3,3)] = az*az;

  rp.gi[aIDX(1,1)] = 1.0/ax/ax;
  rp.gi[aIDX(2,2)] = 1.0/ay/ay;
  rp.gi[aIDX(3,3)] = 1.0/az/az;

  rp.K[aIDX(1,1)] = Kx;
  rp.K[aIDX(2,2)] = Ky;
  rp.K[aIDX(3,3)] = Kz;

  rp.trK = Kx/ax/ax + Ky/ay/ay + Kz/az/az;

  return rp;
}

template <typename RT>
void setKasnerRayCornerPrimitives(RT px, RT py, RT pz, RT t,
  struct RaytracePrimitives<RT> corner_rp[2][2][2])
{
  // Kasner Parameters
  RaytracePrimitives<RT> rp = {0};

  rp = getKasnerRayData(px, py, pz, t);
  corner_rp[0][0][0] = rp;  corner_rp[0][0][1] = rp;  corner_rp[0][1][0] = rp;  corner_rp[0][1][1] = rp;
  corner_rp[1][0][0] = rp;  corner_rp[1][0][1] = rp;  corner_rp[1][1][0] = rp;  corner_rp[1][1][1] = rp;

  return;
}

} // namespace cosmo

#endif
