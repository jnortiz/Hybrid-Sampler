#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>

#include "Sampling.h"
#include "Algebra.h"
#include "params.h"

using namespace std;
using namespace NTL;

//==============================================================================
// Takes in input a random value and samples from distribution D_{\sigma_2}^+,   
// Samples an element in Z+ with probability proportionnal to 2^{-x^2}       
//==============================================================================
unsigned int Sample0(unsigned long alea)
{
    if((alea&1UL)==0UL)
    {
        return 0;
    }
    unsigned int i;
    unsigned int k = 1;
    unsigned long mask=0;
    unsigned long aux;
    for(i=1; i<1000;)
    {
        aux = (alea&mask);
        alea = (alea>>k);     
        if(aux)
        {
            return Sample0(alea);
        }
        else
        {
            if((alea&1UL)==0UL)
            {
                return i;
            }
        }
        i++;
        k += 2;
        mask = (mask<<2)|6UL;
    }
    cout << "ERROR" << endl;
    return 999999;
}

//==============================================================================
// Samples from distribution D_{k\sigma_2}^+, ie 
// Samples an element in Z+ with probability proportionnal to 2^{-(x/k)^2} 
//==============================================================================
unsigned int Sample1(const unsigned int k)
{
    unsigned int x, y, z;
    unsigned long alea = rand();

    x = Sample0(alea);
    y = rand()%k;
    z = k*x + y;
    RR_t w = y*( (z<<1) - y );
    RR_t borne =  LDRMX / exp( w*log_2/(k*k) );
    alea = rand();
    if(alea>borne)
    {
        return Sample1(k);
    }
    else
    {
        return z;
    }
    cout << "ERROR" << endl;
    return 999999;
}




//==============================================================================
// Samples from distribution D_{k\sigma_2}, ie                
// Samples an element in Z with probability proportionnal to 2^{-(x/k)^2} 
//==============================================================================
signed int Sample2(const unsigned int k)
{
    signed int signe;
    signed int x;
    unsigned long alea = rand();
    while(1)
    {
        x = Sample1(k);
        if( (x!=0) || ((alea&1)==1) )
        {
            alea >>= 1;
            signe = 1 - 2*(alea&1);
            x *= signe;
            return x;
        }
        alea >>= 1;
    }
}


//==============================================================================
// Samples from distribution D_{sigma}, ie                                       
// Samples an element in Z with probability proportionnal to e^{-x^2/2*(sigma^2)}
//==============================================================================
signed int Sample3(const RR_t sigma128)
{
    signed int x;
    double alea, borne;

    const RR_t sigma = sigma128;
    const unsigned long k = ( (unsigned long) ceil( (RR_t) sigma/sigma_1 ) );
    while(1)

    {
        x = Sample2(k);
        alea = ((RR_t)rand()) / LDRMX;
        borne = exp( -x*x*( 1/(2*sigma*sigma) - 1/(2*k*k*sigma_1*sigma_1) )   );
        assert(borne<=1);
        if(alea<borne)
        {
            return x;
        }
    }
}


//==============================================================================
// Samples from distribution D_{c,sigma}, ie                                              
// Samples an element in Z with probability proportionnal to e^{-(c-x)^2/2*(sigma^2)}    
//==============================================================================
signed int Sample4(RR_t c, RR_t sigma)
{
    RR_t alea, borne;
    signed int x;
    unsigned int coin;

    const signed int intc = ( (signed int) floor(c) );
    const RR_t fracc = c-intc;
    coin = rand();
    const RR_t denom = 1/(2*sigma*sigma);

    while(1)
    {
        x = Sample3(sigma);
        x += (coin&1);
        if(abs(x)>8){cout << x << endl;}
        coin >>= 1;
        borne = exp(-(x-fracc)*(x-fracc)*denom)/ ( exp(-x*x*denom) + exp(-(x-1)*(x-1)*denom) );

        assert(borne<1);
        alea = ( (RR_t)rand() ) / LDRMX;
        if(alea<borne)
        {
            return (x+intc);
        }
    }
}

/*****************************************************************************/

/*
 * For the hybrid Gaussian sampling, the calling procedures sequence is
 * PrecomputationRingKlein(), PrecomputationRingPeikert(), and then RingKlein().
 */

void PrecomputationRingKlein(vec_RR& sigma_squared, mat_RR& innerp, 
        const mat_RR& BTilde, const ZZX& f, const ZZX& g) {

    vec_RR f_RR, g_RR;
    f_RR.SetLength(N0); g_RR.SetLength(N0);
    
    int i;
    for(i = 0; i < N0; i++) {
        f_RR[i] = to_RR(f[i]);
        g_RR[i] = to_RR(g[i]);
    }//end-for
    
    CC f_FFT[N0], g_FFT[N0];    
    FFTStep(f_FFT, f_RR, N0, omega_CC, N0/2);
    FFTStep(g_FFT, g_RR, N0, omega_CC, N0/2);
    
    CC abs_f, abs_g;
    RR D_omega, over_D, q_RR, norm;
    q_RR = to_RR(q0);
    
    abs_f = abs(f_FFT[0]);
    abs_g = abs(g_FFT[0]);
    D_omega = sqrt(abs_f*abs_f + abs_g*abs_g).real();
    over_D = q_RR/D_omega;
    
    if((D_omega - over_D) > 0.0)
        norm = D_omega;
    else
        norm = over_D;
    
    for(i = 1; i < N0; i++) {        
        abs_f = abs(f_FFT[i]);
        abs_g = abs(g_FFT[i]);
        D_omega = sqrt(abs_f*abs_f + abs_g*abs_g).real();
        over_D = q_RR/D_omega;

        if((D_omega - norm) > 0.0)
            norm = D_omega;
        else
            if((over_D - norm) > 0.0)
                norm = over_D;
        
    }//end-for
    
    RR epsilon, eta, overPi, overTwo;
    div(overPi, 1, ComputePi_RR());
    div(overTwo, to_RR(1), to_RR(2));
    div(epsilon, to_RR(1), to_RR(4)); // \epsilon \in (0, 1/2)
    eta = overPi*sqrt(overTwo*(log(2*N0*(1+1/epsilon)))); // See Lemma 3 from (Ducas and Prest, 2015).
    sigma_squared.SetLength(N0);
    sigma_squared[0] = (norm*eta)*(norm*eta);

    innerp.SetDims(N0, N0);
    for(i = 0; i < N0; i++)
        InnerProduct(innerp[i], BTilde[i], BTilde[i]);
    
}//end-PrecomputationRingKlein()

vec_RR RingKlein(const mat_RR& innerp, const mat_RR& B, const mat_RR& BTilde, 
        const mat_RR& b, const vec_RR& X, const vec_RR& sigma_squared, 
        const vec_RR& c, const RR& sigma0, const RR& eta, const RR& v, 
        const long precision) {
        
    cout << "[!] Hybrid Gaussian sampler status: ";
    
    vec_RR center;
    center = c; // 2N0-length
    
    vec_RR V;
    V.SetLength(2*N0);
    
    CC mult[N0];
    vec_RR d_i, inverse, num, sigma_i, H_mul;
    vec_RR R_sample;
    d_i.SetLength(N0);
    sigma_i.SetLength(N0);
    
    int j;
    for(int i = N0-1; i >= 0; i--) {
        
        InnerProduct(num, center, BTilde[i]);
        
        Inverse(inverse, innerp[i]);
        
        FFTMulMod(mult, num, inverse);
        for(j = 0; j < N0; j++)
            d_i[j] = mult[j].real();
        
        FFTMulMod(mult, sigma_squared, inverse);
        for(j = 0; j < N0; j++)
            sigma_i[j] = mult[j].real();
        
        R_sample = RingPeikert(b[i], d_i, X, v, eta, sigma0, precision);
        
        KByHMultiplication(H_mul, R_sample, B[i]);
        
        for(j = 0; j < (2*N0); j++) {
            center[j] = center[j] - H_mul[j];
            V[j] = V[j] + H_mul[j];
        }//end-for
        
    }//end-for
    
    cout << "Pass!" << endl;
    
    return V;
    
}//end-RingKlein()

vec_RR RingPeikert(const vec_RR& b, const vec_RR c, const vec_RR& X, const RR v, 
              const RR eta, const RR sigma0, const long precision) {

//    cout << "[!] Ring_Peikert status: ";
    
    vec_RR cont_poly;
    cont_poly.SetLength(N0);
    int i;
    for(i = 0; i < N0; i++)
        cont_poly[i] = Ziggurat(X, (int)64, sigma0, precision, v);

    /* Step 1: sampling from \bar K */
    CC p_FFT[N0];    
    FFTMulMod(p_FFT, b, cont_poly);

    Vec< Vec<int> > P;
    Vec<int> begin;

    vec_RR center;
    center.SetLength(N0);
    
    for(i = 0; i < N0; i++)
        center[i] = c[i] + p_FFT[i].real();
    
    RR sigma = sqrt(sigma0)*eta;
    vec_RR output;
    output.SetLength(N0);
    
    BuildProbabilityMatrix(P, begin, precision, 13, sigma, to_RR(0));
    
    /* Step 2: sampling from ring R */
    for(i = 0; i < N0; i++) {
        output[i] = to_RR(KnuthYao(P, begin, 13, sigma, to_RR(0))) + center[i];        
    }//end-for    
    
//    cout << "Pass!" << endl;

    return output;
    
}//end-RingPeikert()

void PrecomputationRingPeikert(mat_RR& b, vec_RR& X, RR& v, RR& eta, 
        const mat_RR& innerp, const vec_RR& sigma_squared, 
        const RR sigma0, const long precision, const RR tailcut) {

//    cout << "[!] Precomputation phase of Ring_Peikert sampler status: ";
    
    /**
     * @param BTilde - Gram-Schmidt orthogonalization of short basis B. 2N0 times 2N0 matrix;
     * @param sigma - Standard deviation of Ring_Klein sampler;
     */
    
    /* Calculus of sigma_0 times \eta */
    vec_RR sigmaEta;
    sigmaEta.SetLength(N0);
    RR epsilon, overPi, overTwo;
    div(overPi, 1, ComputePi_RR());
    div(overTwo, to_RR(1), to_RR(2));
    div(epsilon, to_RR(1), to_RR(4)); // \epsilon \in (0, 1/2)
    eta = overPi*sqrt(overTwo*(log(2*N0*(1+1/epsilon)))); // See Lemma 3 from (Ducas and Prest, 2015).
    mul(sigmaEta[0], sigma0, eta);    

    /* Embedding of sigmaEta in K into C^n */
    CC sigmaEta_FFT[N0];    
    FFTStep(sigmaEta_FFT, sigmaEta, N0, omega_CC, N0/2);
    sigmaEta.kill();
    
    /* Embedding of \sigma^2 in K into C^n */
    CC sigma_FFT[N0];
    FFTStep(sigma_FFT, sigma_squared, N0, omega_CC, N0/2);
    
    /* Computation of parameter b for all \tilde b_i in the Gram-Schmidt basis */
    b.SetDims(N0, N0);    
    CC b_FFT[N0], BTilde_FFT[N0];    
    int j;
    
    for(int i = 0; i < N0; i++) {        
        FFTStep(BTilde_FFT, innerp[i], N0, omega_CC, N0/2);        
        for(j = 0; j < N0; j++)
            b_FFT[j] = sqrt(sigmaEta_FFT[j] - sigma_FFT[j]*Reciprocal(BTilde_FFT[j]));
        ReverseFFTStep(BTilde_FFT, b_FFT, N0, omega_1_CC, N0/2); // Reusing BTilde variable
        for(j = 0; j < N0; j++)
            b[i][j] = BTilde_FFT[j].real();        
    }//end-for
    
    v = ZCreatePartition(X, (int)64, sigma0, precision, tailcut);    
    
//    cout << "Pass!" << endl;
    
}//end-PrecomputationRingPeikert()   

void KByHMultiplication(vec_RR& output, const vec_RR& a, const vec_RR& b) {
    /**
     * 
     * @param output - in H = K^m
     * @param a - in K
     * @param b - in H
     */
    
//    cout << "[!] KByHMultiplication status: ";
    
    vec_RR b1, b2;
    b1.SetLength(N0);
    b2.SetLength(N0);
    
    int i;
    for(i = 0; i < N0; i++) {
        b1[i] = b[i];
        b2[i] = b[i+N0];
    }//end-for
    
    CC mult1[N0], mult2[N0];
    FFTMulMod(mult1, a, b1);
    FFTMulMod(mult2, a, b2);
    
    output.SetLength(2*N0);
    for(i = 0; i < N0; i++) {
        output[i] = mult1[i].real();
        output[i+N0] = mult2[i].real();
    }//end-for
    
//    cout << "Pass!" << endl;
    
}//KByHMultiplication

/*****************************************************************************/

// Gram-Schmidt orthogonalization procedure for block isometric basis
void BlockGSO(mat_RR& BTilde, const mat_RR& B, int n, int precision) {
    
    // Precision of floating point operations
    RR::SetPrecision(precision);    
    
    mat_RR outRot, outGSO; // (n x n)-dimension matrices
    vec_RR ortho, mult;
    RR mu, innerpr, norm;
    int i, j, k;
    
    BTilde.SetDims(B.NumRows(), B.NumCols());
    
    k = B.NumRows()/n; // Number of isometric blocks
    
    ortho.SetLength(n);
    mult.SetLength(n);
    
    for(i = 0; i < k; i++) { // For each isometric block i in [1, ..., k]

        clear(ortho);

         /* For each expansion of a vector makes B[i*n] orthogonal 
          * to the previous vectors */
        for(j = 0; j < i*n; j++) {
            if(!IsZero(BTilde[j])) {
                NTL::InnerProduct(innerpr, B[i*n], BTilde[j]);
                NTL::InnerProduct(norm, BTilde[j], BTilde[j]); // Square of norm is the inner product
                div(mu, innerpr, norm);           
                mul(mult, BTilde[j], mu);
                add(ortho, ortho, mult);
            }//end-if
        }//end-for
        
        sub(BTilde[i*n], B[i*n], ortho);            
                
        // Expansion of the vector that is already orthogonal to its predecessors
        // In (Lyubashevsky, and Prest, 2015), it uses B[i*n] instead
        rot(outRot, BTilde[i*n], n);
                
        // Computes its orthogonal basis
        FasterIsometricGSO(outGSO, outRot);
        
        // Copying the orthogonal basis of B[i*n] to the output
        for(j = 0; j < n; j++)
            BTilde[i*n + j] = outGSO[j];
        
    }//end-for    
    
}//end-BlockGSO()

// Gram-Schmidt orthogonalization procedure for isometric basis
void FasterIsometricGSO(mat_RR& BTilde, const mat_RR& B) {
    
    /* This implementation follows the Algorithm 3 
     * of (Lyubashevsky, and Prest, 2015) */
    if(IsZero(B)) {
//        cout << "[!] Warning: the input basis is null." << endl;
        BTilde = B;
        return;
    }//end-if
    
    vec_RR C, D, isometry, mult, V;
    RR CD;
    int n;
    
    n = B.NumCols(); // B is a square matrix
    
    BTilde.SetDims(n, n);
    C.SetLength(n);
    D.SetLength(n);
        
    isometry.SetMaxLength(n);
    mult.SetMaxLength(n);
    V.SetMaxLength(n);
    
    BTilde[0] = B[0];
    V = B[0];
    
    NTL::InnerProduct(C[0], V, Isometry(BTilde[0]));
    NTL::InnerProduct(D[0], BTilde[0], BTilde[0]);
    
    for(int i = 0; i < n-1; i++) {        
        div(CD, C[i], D[i]);        
        isometry = Isometry(BTilde[i]);         
        mul(mult, V, CD);                
        sub(BTilde[i+1], isometry, mult);            
        mul(mult, isometry, CD);        
        sub(V, V, mult);
        NTL::InnerProduct(C[i+1], B[0], Isometry(BTilde[i+1]));
        sub(D[i+1], D[i], CD*C[i]);                
    }//end-for
    
}//end-FasterIsometricGSO()

RR NormOfBasis(const mat_RR& B) {
    
    RR innerp, norm, normB;    
    
    normB = to_RR(0);
    
    for(int i = 0; i < B.NumRows(); i++) {
        NTL::InnerProduct(innerp, B[i], B[i]);
        norm = sqrt(innerp);
        if(norm > normB)
            normB = norm;
    }//end-for
    
    return normB;
    
}//end-NormOfBasis()

void rot(mat_RR& out, const vec_RR& b, int n) {
        
    vec_RR isometry;
    isometry.SetLength(n);
    
    out.SetDims(n, n);    
    
    out[0] = b;    
    isometry = b;    
    
    for(int i = 1; i < n; i++) {
        isometry = Isometry(isometry);
        out[i] = isometry;
    }//end-for
    
}//end-rot()

vec_RR Isometry(const vec_RR& b) {
    
    int n = b.length();    
    vec_RR out;    
    out.SetLength(n);
    
        
    out[0] = -b[n-1];    
    for(int i = 1; i < n; i++)
        out[i] = b[i-1];
    
    return out;
    
}//end-Isometry()

//ZZX Ring_Klein(const mat_RR& B, const mat_RR& BTilde, const vec_RR& Sigma, const vec_RR& Center) {
//    /**
//     * @param B - Basis B - (2*N0, 2*N0)
//     * @param BTilde - Gram-Schmidt Orthogonalization of B - (2*N0, 2*N0)
//     * @param Sigma - Standard deviation - (N0)
//     * @param Center - Center of the distribution - (m*N0 = 2*N0)
//     * @return A vector sampling according to D_{Span_R(B), \sigma, \center} - (m*N0 = 2*N0)
//     */
//    vec_RR c_i, d_i, inner_p, v_i;
//    vec_RR r, sigma, sigma_i;
//    CC mult[N0];
//    ZZX sample;
//    c_i = Center;
//    sigma = Sigma;
//    v_i.SetLength(Center.length());    
//    
//    for(int i = (2*N0-1); i >= 0; i++) {
//        InnerProduct(d_i, c, Block_BTilde[i]);
//        InnerProduct(inner_p, Block_BTilde[i], Block_BTilde[i]);
//        EuclideanDiv(d_i, r, d_i, inner_p);
//        FFTMulMod(mult, sigma, sigma);
//        for(j = 0; j < N0; j++)
//            sigma[i] = mult[i].real();
//        EuclideanDiv(sigma_i, r, sigma, inner_p);
////        Ring_Peikert();
//        
//    }//end-for
//    
//    return sample;
//     
//}//end-Ring_Klein

/* DZCreatePartition defines the x- and y-axes of the rectangles in the Gaussian distribution */
RR ZCreatePartition(vec_RR& X, int m, RR sigma, long n, RR tail) {
    
    RR::SetPrecision(n);
    
    /* The Ziggurat algorithm was designed for centered Gaussian distributions; 
     i.e., c = 0. */
    
    /*
     * Parameters description:
     * m: number of rectangles
     * sigma: Gaussian distribution parameter
     * n: bit precision
     */
    
    Vec<RR> bestX, auxX;
    RR bestdiff, r, v, z;
    
    RR overM = to_RR(1)/to_RR(m);
    RR minusOne = to_RR(-1);
    RR zero = to_RR(0);
    
    int i;
    int first = 1;
    
    auxX.SetLength(m);
    bestX.SetLength(m);            
    
    bestX[m-1] = minusOne;    
    z = minusOne;    
    r = (tail*sigma)/2; // r = x_m, the m-th x-value
    v = 0;
    bestdiff = r; // Arbitrary initial value
    
    while(z != zero && r > zero) { // If r = 0 then the area v is also 0 (what doesn't make sense...)

        z = ZRecursion(auxX, m, r, sigma, v);        

        if(z == minusOne && first) { // Error in "inv" or square root computation
            
            first = 0;
            add(r, r, 2*overM);
            div(overM, overM, to_RR(m));
            
            while(z != zero && r > zero) { // If r = 0 then the area v is also 0 (what doesn't make sense...)
                
                z = ZRecursion(auxX, m, r, sigma, v);

                if(z == minusOne)            
                    break;

                if(abs(z) < abs(bestdiff)) { // If the actual z is closest to zero, then that's the better partitioning until now
                    for(int i = 0; i < m; i++)
                        bestX[i] = auxX[i];                
                    bestdiff = z;
                }//end-if

                sub(r, r, overM);
                
            }//end-while            
            
        }//end-if
        
        if(z == minusOne && !first)
            break;
        
        if(abs(z) < abs(bestdiff)) { // If the actual z is closest to zero, then that's the better partitioning until now
            for(int i = 0; i < m; i++)
                bestX[i] = auxX[i];                
            bestdiff = z;
        }//end-if
        
        sub(r, r, overM);
        
    }//end-while
    
    if(z == zero)
        for(int i = 0; i < m; i++)
            bestX[i] = auxX[i];                    
               
    if(bestX[m-1] != -1) { // Some partitioning was found
        
        X.SetLength(m);
        for(i = 0; i < X.length(); i++)
            X[i] = bestX[i];
        
    }//end-if
    
    return v;
    
}//end-DZCreatePartition()

/* It is used in DZCreatePartition to define the distance y0 */
RR ZRecursion(Vec<RR>& X, int m, RR r, RR sigma, RR& v) {
    
    RR zero = to_RR(0);
    RR c = zero; // Center of distribution
    RR interm, overPi;
    RR minusOne = to_RR(-1);
    div(overPi, minusOne, ComputePi_RR());
    
    X[m-1] = r;
    
    v = r*Probability(r, sigma, c) /*+ integrate r to infinity Probability(x) */;
    
    for(int i = (m-2); i >= 0; i--) {
        
        if(X[i+1] == zero)
            return minusOne;
        
        // TODO: general formula for rho^{-1}(x)
        // This inversion of rho(x) works only when variance is equal to 1/2*pi
        interm = overPi * log(v/X[i+1] + Probability(X[i+1], sigma, c));
        
        if(interm < zero)
            return minusOne;
            
        X[i] = sqrt(interm);
       
    }//end-for
    
    return (v - X[0] + X[0]*Probability(X[0], sigma, c));
        
}//end-DZCreatePartition()

/* Given a threshold r, this function returns a value in the tail region (Thomas et al., 2007) */
RR NewMarsagliaTailMethod(RR r) {
    
    /* r is the right most x-coordinate in Ziggurat partitioning */
    RR a, b;
    RR x, y;
    double ulong_size = (double)(sizeof(unsigned long)*8);
    
    do {                        
        // Two variables with values in (0, 1)
        a = to_RR((RandomBnd(ulong_size-1) + 1)/ulong_size);
        b = to_RR((RandomBnd(ulong_size-1) + 1)/ulong_size);
        div(x, -log(a), r);
        y = -log(b);
    } while(y + y > x*x);
    
    return (r+x);
    
}//end-NewMarsagliaTailMethod()

RR Ziggurat(const vec_RR& X, int m, RR sigma, long precision, RR v) {

/*
 * Important: run ZCreatePartition() procedure before calling the Ziggurat sampler.
 */
    
#ifdef DEBUG
    cout << "[*] Ziggurat status: ";
#endif
    
    RR::SetPrecision(precision);
    
    RR c, probX, probX1, probZ, U, z;
    int i;
    double ulong_size = (double)(sizeof(unsigned long)*8);
    
    c = to_RR(0);
    
    for(;;) {
    
        i = NTL::RandomBnd((long)(m-1)) + 1;
        U = to_RR((NTL::RandomBnd((long)(ulong_size-1)) + 1)/ulong_size); // Uniform sample in (0, 1)
        
        if(i < (m-1))
            NTL::mul(z, U, X[i]);
        else // Base strip
            NTL::div(z, U*v, Probability(X[i], sigma, c));
            
        if(z < X[i-1]) // The sample is in the left-most part of the i-th rectangle
            return z;

        if(i == (m-1))
            return NewMarsagliaTailMethod(X[m-1]);
        
        probX = Probability(X[i], sigma, c);
        probX1 = Probability(X[i+1], sigma, c);
        probZ = Probability(z, sigma, c);

        U = to_RR((RandomBnd(ulong_size-1) + 1)/ulong_size);

        // The sample can be in the right-most part of the i-th rectangle        
        if(i > 0 && U*(probX - probX1) < (probZ - probX1))
            return z;
            
    }//end-for
    
#ifdef DEBUG
    cout << "Pass!" << endl;
#endif
    
}//end-Ziggurat()

int KnuthYao(const Vec< Vec<int> >& P, const Vec<int>& begin, int tailcut, RR sigma, RR c) {

    int bound, center, col, d, invalidSample, pNumRows, pNumCols, S, signal;
    unsigned enable, hit;
    unsigned long r;
    
    bound = tailcut*to_int(sigma);
    center = to_int(c);
    d = 0;
    hit = 0;
    signal = 1 - 2*RandomBits_long(1);
    invalidSample = bound+1;
    pNumRows = P.length();
    pNumCols = P[0].length();    
    
    Vec<int> randomBits;
    randomBits.SetLength(pNumRows);
    
    int i, index, j, length;
    length = sizeof(unsigned long)*8; 
    
    index = 0;
    for(i = 0; i < (pNumRows/length+1); i++) {
        r = RandomWord();
        for(j = 0; j < length && index < pNumRows; j++, r >>= 1)
            randomBits[index++] = (r & 1);
    }//end-for

    S = 0;    
    for(int row = 0; row < pNumRows; row++) {
        
        d = 2*d + randomBits[row]; // Distance calculus
        
        for(col = begin[row]; col < pNumCols; col++) {
            
            d = d - P[row][col];
            
            enable = (unsigned)(d + 1); // "enable" turns 0 iff d = -1
            enable = (1 ^ ((enable | -enable) >> 31)) & 1; // "enable" turns 1 iff "enable" was 0
             
            /* When enable&!hit becomes 1, "col" is added to "S";
             * e.g. enable = 1 and hit = 0 */
            S += Select(invalidSample, col, (enable & !hit));
            hit += (enable & !hit);
                            
        }//end-for
        
    }//end-for
    
    /* Note: the "col" value is in [0, bound]. So, the invalid sample must be 
     * greater than bound. */
    S %= invalidSample;
    S = S - bound + center;
    S *= signal;
    
    return S;
    
}//end-Knuth-Yao()

/* This method build the probability matrix for samples in the range 
 * [-tailcut*\floor(sigma), +tailcut*\floor(sigma)] */
void BuildProbabilityMatrix(Vec< Vec<int> >& P, Vec<int>& begin, int precision, int tailcut, RR sigma, RR c) {
    
    RR::SetPrecision(to_long(precision));

    Vec< Vec<int> > auxP;
    Vec<int> auxBegin;
    
    // The random variable consists of elements in [c-tailcut*sigma, c+tailcut*sigma]
    int i, j, bound, pNumCols, pNumRows, x;
    vec_RR probOfX;
    RR pow;
    
    bound = tailcut*to_int(sigma);
    
    probOfX.SetLength(bound+1);
       
    auxP.SetLength(precision);
    for(i = 0; i < auxP.length(); i++)
        auxP[i].SetLength(bound+1);

    for(x = bound; x > 0; x--)
        probOfX[bound-x] = Probability(to_RR(x) + c, sigma, c);
    div(probOfX[bound], Probability(to_RR(0) + c, sigma, c), to_RR(2));
    
    i = -1;
    for(j = 0; j < precision; j++) {
        pow = power2_RR(i--); // 2^{i}
        for(x = bound; x >= 0; x--) {
            auxP[j][bound-x] = 0;                
            if(probOfX[bound-x] >= pow) {
                auxP[j][bound-x] = 1;
                probOfX[bound-x] -= pow;
            }//end-if
        }//end-for
    }//end-while
    
    P = auxP;
    
    pNumCols = P[0].length();
    pNumRows = P.length();
    
    auxBegin.SetLength(pNumRows);
    
    // Computing in which position the non-zero values in P start and end 
    for(i = 0; i < pNumRows; i++) {
        
        auxBegin[i] = pNumCols-1;
        
        for(j = 0; j < pNumCols; j++)
            if(P[i][j] == 1) {
                auxBegin[i] = j;
                break;
            }//end-if
        
    }//end-for
    
    begin = auxBegin;
                
}//end-BuildProbabilityMatrix()

// If the bit is zero, the output becomes "a" 
int Select(int a, int b, unsigned bit) {
    
    unsigned mask;
    int output;
    
    mask = -bit;
    output = mask & (a ^ b);
    output = output ^ a;
    
    return output;
    
}//end-Select()

RR Probability(RR x, RR sigma, RR c) {
    
    RR S = sigma*sqrt(2*ComputePi_RR());
    RR overS = 1/S;
    
    if(x == to_RR(0))
        return overS;
    
    return overS*exp(-(power((x-c)/sigma, 2))/2.0);
    
}//end-Probability()
