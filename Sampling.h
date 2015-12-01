#ifndef LIBE_SAMPLING_H
#define LIBE_SAMPLING_H

#include <NTL/RR.h>
#include <NTL/mat_RR.h>
#include <NTL/matrix.h>

#include "params.h"

/* Sampling algorithms of Thomas Prest */
unsigned int Sample0(unsigned long alea);
unsigned int Sample1(const unsigned int k);
signed int Sample2(const unsigned int k);
signed int Sample3(const RR_t sigma128);
signed int Sample4(RR_t c, RR_t sigma);

/* Gram-Schmidt orthogonalization for block isometric basis of NTRU lattices */
void BlockGSO(mat_RR& BTilde, const mat_RR& B, int n, int precision);
void FasterIsometricGSO(mat_RR& BTilde, const mat_RR& B);

/* Auxiliary functions of Block-GSO and Faster-Isometric-GSO */
RR NormOfBasis(const mat_RR& B);
void rot(mat_RR& out, const vec_RR& b, int n);
RR NormOfBasis(const mat_RR& B);
vec_RR Isometry(const vec_RR& b);

/* Ring-Peikert for sampling from the ring R = Z[x]/<x^n + 1>*/
ZZX RingPeikert(vec_RR& c, const vec_RR& X, RR b, RR eta, RR factor, RR v, long precision, int m, int tailcut);
void OfflineRingPeikert(vec_RR& X, vec_RR& c, RR& b, RR& v, long precision, int m, int n, int tailcut, int q);

/* Continuous Gaussian sampling */
RR Ziggurat(const vec_RR& X, int m, RR sigma, long precision, RR v);
RR ZCreatePartition(vec_RR& X, int m, RR sigma, long n, RR tail);
RR ZRecursion(Vec<RR>& X, int m, RR r, RR sigma, RR& v);
RR NewMarsagliaTailMethod(RR r);

/* Discrete Gaussian sampling from integers */
int KnuthYao(const Vec< Vec<int> >& P, const Vec<int>& begin, int tailcut, RR sigma, RR c);
void BuildProbabilityMatrix(Vec< Vec<int> >& P, Vec<int>& begin, int precision, int tailcut, RR sigma, RR c);
RR Probability(RR x, RR sigma, RR c);
int Select(int a, int b, unsigned bit);

void MulPoly(vec_RR& x, const vec_RR& a, const vec_RR& b);

#endif
