#ifndef RAYTRACE_DATA
#define RAYTRACE_DATA

namespace cosmo
{

// To calculate Riemann tensor components, we need all christoffel
//  symbols (raised and lowered), and all first and second (including
// mixed) derivatives of the metric and extrinsic curvature.
// indices follow the convention: g_{ij} = g_{a}, where a=(7-i)*i/2-4+j
// for i<=j, and i,j are in [1,3]. (Eg., g_{23} = g_{4}.)
// metric:
template<typename RT>
struct RaytracePrimitives {
  RT g[6];      // 6 independent metric components
  RT gi[6];     // inverse metric
  RT dg[3][6];  // 3 first derivatives (first index is deriv. index)
  RT ddg[6][6]; // 6 second derivatives (first index is deriv. index)
  // extrinsic curvature:
  RT K[6];      // 6 independent extrinsic curvature components
  RT dK[3][6];  // 3 first derivatives (first index is deriv. index)
  // 3-Christoffel symbols
  // raised first index
  RT G[3][6];
  // lower first index
  RT GL[3][6];
  // Ricci tensor
  RT Ricci[6];
  RT rho;
  RT trK;
};

template<typename RT>
struct RaytraceData {
  // ray position
  RT x[3];
  RT x_d[3]; // normalized ray position for interpolation.
             // given a gridpoint x0 < x < x0 + dx
             // x_d = (x - x0) / dx
  // ray velocity vector
  RT V[3];
  // ray Energy
  RT E;
  // screen vectors
  RT S1[3];
  RT S2[3];
  // beam root area b, and lambda-derivative
  RT b, Omega;
  // sachs scalar sigma
  RT sig_Re;
  RT sig_Im;
};

} // namespace cosmo

#endif
