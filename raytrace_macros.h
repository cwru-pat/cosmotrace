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
    - rp.K[aIDX(i,m)]*rp.K[aIDX(j,l)] + rp.K[aIDX(i,j)]*rp.K[aIDX(m,l)] \
  )

#define RIEMANN_I0MJ(i, m, j) ( \
    rp.dK[j][aIDX(i,m)] - rp.dK[m][aIDX(i,j)] \
    + rp.gi[aIDX(1,1)]*( rp.dg[iIDX(m)][aIDX(1,i)] + rp.GL[iIDX(i)][aIDX(m,1)] )*rp.K[aIDX(j,1)] \
      + rp.gi[aIDX(1,2)]*( rp.dg[iIDX(m)][aIDX(1,i)] + rp.GL[iIDX(i)][aIDX(m,1)] )*rp.K[aIDX(j,2)] \
      + rp.gi[aIDX(1,3)]*( rp.dg[iIDX(m)][aIDX(1,i)] + rp.GL[iIDX(i)][aIDX(m,1)] )*rp.K[aIDX(j,3)] \
    + rp.gi[aIDX(2,1)]*( rp.dg[iIDX(m)][aIDX(2,i)] + rp.GL[iIDX(i)][aIDX(m,2)] )*rp.K[aIDX(j,1)] \
      + rp.gi[aIDX(2,2)]*( rp.dg[iIDX(m)][aIDX(2,i)] + rp.GL[iIDX(i)][aIDX(m,2)] )*rp.K[aIDX(j,2)] \
      + rp.gi[aIDX(2,3)]*( rp.dg[iIDX(m)][aIDX(2,i)] + rp.GL[iIDX(i)][aIDX(m,2)] )*rp.K[aIDX(j,3)] \
    + rp.gi[aIDX(3,1)]*( rp.dg[iIDX(m)][aIDX(3,i)] + rp.GL[iIDX(i)][aIDX(m,3)] )*rp.K[aIDX(j,1)] \
      + rp.gi[aIDX(3,2)]*( rp.dg[iIDX(m)][aIDX(3,i)] + rp.GL[iIDX(i)][aIDX(m,3)] )*rp.K[aIDX(j,2)] \
      + rp.gi[aIDX(3,3)]*( rp.dg[iIDX(m)][aIDX(3,i)] + rp.GL[iIDX(i)][aIDX(m,3)] )*rp.K[aIDX(j,3)] \
    - rp.gi[aIDX(1,1)]*( rp.dg[iIDX(j)][aIDX(1,i)] + rp.GL[iIDX(i)][aIDX(j,1)] )*rp.K[aIDX(m,1)] \
      - rp.gi[aIDX(1,2)]*( rp.dg[iIDX(j)][aIDX(1,i)] + rp.GL[iIDX(i)][aIDX(j,1)] )*rp.K[aIDX(m,2)] \
      - rp.gi[aIDX(1,3)]*( rp.dg[iIDX(j)][aIDX(1,i)] + rp.GL[iIDX(i)][aIDX(j,1)] )*rp.K[aIDX(m,3)] \
    - rp.gi[aIDX(2,1)]*( rp.dg[iIDX(j)][aIDX(2,i)] + rp.GL[iIDX(i)][aIDX(j,2)] )*rp.K[aIDX(m,1)] \
      - rp.gi[aIDX(2,2)]*( rp.dg[iIDX(j)][aIDX(2,i)] + rp.GL[iIDX(i)][aIDX(j,2)] )*rp.K[aIDX(m,2)] \
      - rp.gi[aIDX(2,3)]*( rp.dg[iIDX(j)][aIDX(2,i)] + rp.GL[iIDX(i)][aIDX(j,2)] )*rp.K[aIDX(m,3)] \
    - rp.gi[aIDX(3,1)]*( rp.dg[iIDX(j)][aIDX(3,i)] + rp.GL[iIDX(i)][aIDX(j,3)] )*rp.K[aIDX(m,1)] \
      - rp.gi[aIDX(3,2)]*( rp.dg[iIDX(j)][aIDX(3,i)] + rp.GL[iIDX(i)][aIDX(j,3)] )*rp.K[aIDX(m,2)] \
      - rp.gi[aIDX(3,3)]*( rp.dg[iIDX(j)][aIDX(3,i)] + rp.GL[iIDX(i)][aIDX(j,3)] )*rp.K[aIDX(m,3)] \
  )

#define RIEMANN_I00J(i, j) \
  /* todo - set tr(K) */ \
  -rp.Ricci[aIDX(i,j)] - rp.trK*rp.K[aIDX(i,j)] + 4.0*RAY_PI*rp.g[aIDX(i,j)]*rp.rho

#define LINEAR_INTERPOLATION(var) \
  linearInterpolation( \
    corner_rp[0][0][0].var, corner_rp[0][0][1].var, corner_rp[0][1][0].var, corner_rp[0][1][1].var, \
    corner_rp[1][0][0].var, corner_rp[1][0][1].var, corner_rp[1][1][0].var, corner_rp[1][1][1].var \
  )

#define WeylLensingScalarSum(A, B) ( \
  -4.0*( \
    - R_1001 * varSig##A[0 /*01*/] * varSig##B[0 /*01*/] \
    - R_2002 * varSig##A[1 /*02*/] * varSig##B[1 /*02*/] \
    - R_3003 * varSig##A[2 /*03*/] * varSig##B[2 /*03*/] \
    + R_1212 * varSig##A[3 /*12*/] * varSig##B[3 /*12*/] \
    + R_1313 * varSig##A[4 /*13*/] * varSig##B[4 /*13*/] \
    + R_2323 * varSig##A[5 /*23*/] * varSig##B[5 /*23*/] \
    - R_1002*( varSig##A[0 /*01*/]*varSig##B[1 /*02*/] + varSig##A[1 /*02*/]*varSig##B[0 /*01*/] ) \
      - R_1003*( varSig##A[0 /*01*/]*varSig##B[2 /*03*/] + varSig##A[2 /*03*/]*varSig##B[0 /*01*/] ) \
      - R_1012*( varSig##A[0 /*01*/]*varSig##B[3 /*12*/] + varSig##A[3 /*12*/]*varSig##B[0 /*01*/] ) \
    - R_1013*( varSig##A[0 /*01*/]*varSig##B[4 /*13*/] + varSig##A[4 /*13*/]*varSig##B[0 /*01*/] ) \
      - R_1023*( varSig##A[0 /*01*/]*varSig##B[5 /*23*/] + varSig##A[5 /*23*/]*varSig##B[0 /*01*/] ) \
      - R_2003*( varSig##A[1 /*02*/]*varSig##B[2 /*03*/] + varSig##A[2 /*03*/]*varSig##B[1 /*02*/] ) \
    - R_2012*( varSig##A[1 /*02*/]*varSig##B[3 /*12*/] + varSig##A[3 /*12*/]*varSig##B[1 /*02*/] ) \
      - R_2013*( varSig##A[1 /*02*/]*varSig##B[4 /*13*/] + varSig##A[4 /*13*/]*varSig##B[1 /*02*/] ) \
      - R_2023*( varSig##A[1 /*02*/]*varSig##B[5 /*23*/] + varSig##A[5 /*23*/]*varSig##B[1 /*02*/] ) \
    - R_3012*( varSig##A[2 /*03*/]*varSig##B[3 /*12*/] + varSig##A[3 /*12*/]*varSig##B[2 /*03*/] ) \
      - R_3013*( varSig##A[2 /*03*/]*varSig##B[4 /*13*/] + varSig##A[4 /*13*/]*varSig##B[2 /*03*/] ) \
      - R_3023*( varSig##A[2 /*03*/]*varSig##B[5 /*23*/] + varSig##A[5 /*23*/]*varSig##B[2 /*03*/] ) \
    + R_1213*( varSig##A[3 /*12*/]*varSig##B[4 /*13*/] + varSig##A[4 /*13*/]*varSig##B[3 /*12*/] ) \
      + R_1223*( varSig##A[3 /*12*/]*varSig##B[5 /*23*/] + varSig##A[5 /*23*/]*varSig##B[3 /*12*/] ) \
      + R_1323*( varSig##A[4 /*13*/]*varSig##B[5 /*23*/] + varSig##A[5 /*23*/]*varSig##B[4 /*13*/] ) \
  ) )

#endif
