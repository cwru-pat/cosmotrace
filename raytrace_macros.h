#ifndef RAYTRACE_MACROS
#define RAYTRACE_MACROS

#define RAY_PI 3.14159265358979323846

// map spatial (i,j) to array index
#define aIDX(i,j) ( i <= j  ? (7-i)*i/2 - 4 + j : (7-j)*j/2 - 4 + i )
// map spatial (i) to array index
#define iIDX(i) ( i - 1 )

#define RIEMANN_ILMJ(i, l, m, j) ( \
    0.5*( \
      rp.ddg[aIDX(l,m)][aIDX(i,j)] + rp.ddg[aIDX(i,j)][aIDX(m,l)] \
      - rp.ddg[aIDX(i,m)][aIDX(j,l)] - rp.ddg[aIDX(l,j)][aIDX(i,m)] \
    ) \
    + rp.GL[iIDX(1)][aIDX(i,j)]*rp.G[iIDX(1)][aIDX(l,m)] \
      + rp.GL[iIDX(2)][aIDX(i,j)]*rp.G[iIDX(2)][aIDX(l,m)] \
      + rp.GL[iIDX(3)][aIDX(i,j)]*rp.G[iIDX(3)][aIDX(l,m)] \
    - rp.GL[iIDX(1)][aIDX(i,m)]*rp.G[iIDX(1)][aIDX(j,l)] \
      - rp.GL[iIDX(2)][aIDX(i,m)]*rp.G[iIDX(2)][aIDX(j,l)] \
      - rp.GL[iIDX(3)][aIDX(i,m)]*rp.G[iIDX(3)][aIDX(j,l)] \
    + rp.K[aIDX(i,m)]*rp.K[aIDX(j,l)] - rp.K[aIDX(i,j)]*rp.K[aIDX(m,l)] \
  )

#define RIEMANN_I0MJ(i, m, j) ( \
    rp.dK[iIDX(j)][aIDX(i,m)] - rp.dK[iIDX(m)][aIDX(i,j)] \
    + rp.K[aIDX(j,1)]*rp.G[iIDX(1)][aIDX(i,m)] + rp.K[aIDX(j,1)]*rp.G[iIDX(1)][aIDX(i,m)] + rp.K[aIDX(j,1)]*rp.G[iIDX(1)][aIDX(i,m)] \
    - rp.K[aIDX(m,1)]*rp.G[iIDX(1)][aIDX(i,j)] - rp.K[aIDX(m,1)]*rp.G[iIDX(1)][aIDX(i,j)] - rp.K[aIDX(m,1)]*rp.G[iIDX(1)][aIDX(i,j)] \
  )

#define RIEMANN_I00J(i, j) \
  rp.g[aIDX(1,1)]*rp.K[aIDX(i,1)]*rp.K[aIDX(1,j)] + rp.g[aIDX(1,2)]*rp.K[aIDX(i,1)]*rp.K[aIDX(2,j)] + rp.g[aIDX(1,3)]*rp.K[aIDX(i,1)]*rp.K[aIDX(3,j)] \
  + rp.g[aIDX(2,1)]*rp.K[aIDX(i,2)]*rp.K[aIDX(1,j)] + rp.g[aIDX(2,2)]*rp.K[aIDX(i,2)]*rp.K[aIDX(2,j)] + rp.g[aIDX(2,3)]*rp.K[aIDX(i,2)]*rp.K[aIDX(3,j)] \
  + rp.g[aIDX(3,1)]*rp.K[aIDX(i,3)]*rp.K[aIDX(1,j)] + rp.g[aIDX(3,2)]*rp.K[aIDX(i,3)]*rp.K[aIDX(2,j)] + rp.g[aIDX(3,3)]*rp.K[aIDX(i,3)]*rp.K[aIDX(3,j)] \
  - rp.Ricci[aIDX(i,j)] - rp.trK*rp.K[aIDX(i,j)] + 4.0*RAY_PI*rp.g[aIDX(i,j)]*rp.rho

#define LINEAR_INTERPOLATION(var) \
  linearInterpolation( \
    corner_rp[0][0][0].var, corner_rp[0][0][1].var, corner_rp[0][1][0].var, corner_rp[0][1][1].var, \
    corner_rp[1][0][0].var, corner_rp[1][0][1].var, corner_rp[1][1][0].var, corner_rp[1][1][1].var \
  )

#define WEYL_LENSING_SCALAR_SUM(A, B) ( \
  2.0*( \
    - R_1001 * (varSig.SK##A[0]) * (varSig.SK##B[0]) \
    - R_2002 * (varSig.SK##A[1]) * (varSig.SK##B[1]) \
    - R_3003 * (varSig.SK##A[2]) * (varSig.SK##B[2]) \
    + R_1212 * (varSig.SK##A[3]) * (varSig.SK##B[3]) \
    + R_1313 * (varSig.SK##A[4]) * (varSig.SK##B[4]) \
    + R_2323 * (varSig.SK##A[5]) * (varSig.SK##B[5]) \
    - R_1002*( (varSig.SK##A[0])*(varSig.SK##B[1]) + (varSig.SK##A[1])*(varSig.SK##B[0]) ) \
      - R_1003*( (varSig.SK##A[0])*(varSig.SK##B[2]) + (varSig.SK##A[2])*(varSig.SK##B[0]) ) \
      - R_1012*( (varSig.SK##A[0])*(varSig.SK##B[3]) + (varSig.SK##A[3])*(varSig.SK##B[0]) ) \
    - R_1013*( (varSig.SK##A[0])*(varSig.SK##B[4]) + (varSig.SK##A[4])*(varSig.SK##B[0]) ) \
      - R_1023*( (varSig.SK##A[0])*(varSig.SK##B[5]) + (varSig.SK##A[5])*(varSig.SK##B[0]) ) \
      - R_2003*( (varSig.SK##A[1])*(varSig.SK##B[2]) + (varSig.SK##A[2])*(varSig.SK##B[1]) ) \
    - R_2012*( (varSig.SK##A[1])*(varSig.SK##B[3]) + (varSig.SK##A[3])*(varSig.SK##B[1]) ) \
      - R_2013*( (varSig.SK##A[1])*(varSig.SK##B[4]) + (varSig.SK##A[4])*(varSig.SK##B[1]) ) \
      - R_2023*( (varSig.SK##A[1])*(varSig.SK##B[5]) + (varSig.SK##A[5])*(varSig.SK##B[1]) ) \
    - R_3012*( (varSig.SK##A[2])*(varSig.SK##B[3]) + (varSig.SK##A[3])*(varSig.SK##B[2]) ) \
      - R_3013*( (varSig.SK##A[2])*(varSig.SK##B[4]) + (varSig.SK##A[4])*(varSig.SK##B[2]) ) \
      - R_3023*( (varSig.SK##A[2])*(varSig.SK##B[5]) + (varSig.SK##A[5])*(varSig.SK##B[2]) ) \
    + R_1213*( (varSig.SK##A[3])*(varSig.SK##B[4]) + (varSig.SK##A[4])*(varSig.SK##B[3]) ) \
      + R_1223*( (varSig.SK##A[3])*(varSig.SK##B[5]) + (varSig.SK##A[5])*(varSig.SK##B[3]) ) \
      + R_1323*( (varSig.SK##A[4])*(varSig.SK##B[5]) + (varSig.SK##A[5])*(varSig.SK##B[4]) ) \
  ) )

#endif
