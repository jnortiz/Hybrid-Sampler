#ifndef LIBE_ALGEBRA_H
#define LIBE_ALGEBRA_H

#include <NTL/vec_RR.h>
#include "params.h"

ZZX Cyclo();
ZZX FastMod(const ZZX& f);
void FastMod(CC * const out, CC const * const in, const int N);
void FastMod(vec_RR& f);
ZZ SquaredNorm(const ZZX& f, const unsigned int degree);
void ValidPair(ZZ& PGCD, ZZ& Alpha, ZZ& Beta, ZZX& rho_f, ZZX& rho_g, const ZZX& f, const ZZX& g);
ZZX Reverse(const ZZX& f);
ZZX ReductionCoefficient(const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G, unsigned int & mb);
ZZX FastReductionCoefficient(const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G);
mat_ZZ AnticircularMatrix(const ZZX& f);
mat_ZZ BasisFromPolynomials(const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G);
ZZ_pX Inverse(const ZZX& f);
ZZ_pX Quotient(const ZZX& f, const ZZX& g);
void GS_Norm(const ZZX fx, const ZZX gx, int& flag);
void GenerateBasis(ZZX& f, ZZX& g, ZZX& F, ZZX& G, const ZZ& Norme);
RR_t DotProduct(const RR_t * x1, const RR_t * x2);
void Rotate(RR_t * const dest, RR_t const * const src);
void ClassicMGS(RR_t Bstar[2*N0][2*N0], const RR_t B[2*N0][2*N0]);
void FastMGS(RR_t Bst[2*N0][2*N0], const RR_t B[2*N0][2*N0]);
void FFTMulMod(CC * const c, const vec_RR a, const vec_RR b);
void FFTMul(CC * const c, const vec_RR& a, const vec_RR& b);

void InnerProduct(vec_RR& innerp, const vec_RR& a, const vec_RR& b);
void Inverse(vec_RR& inverse, const vec_RR& a);
void SquareRoot(vec_RR& sqr, const vec_RR& a);
CC Reciprocal(const CC c);

#endif
