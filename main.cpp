/* 
 * File:   main.cpp
 * Author: jnortiz
 *
 * Created on November 26, 2015, 9:42 AM
 */

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>

#include "params.h"
#include "io.h"
#include "FFT.h"
#include "Sampling.h"
#include "Random.h"
#include "Algebra.h"
#include "Scheme.h"

#define PRECISION 128
//#define PRINT_ALL 1

using namespace std;
using namespace NTL;

int main() {

    srand(rdtsc());
    
    ZZX MSK[4];
    ZZ_pX MPK;
    
    MSK_Data *MSKD = new MSK_Data;
    
    ZZ_p::init(q1);
    zz_p::init(q0);

    cout << "[!] Parameters (q, N0) = (" << q0 << ", " << N0 << ")" << endl;
    
    Keygen(MPK, MSK);
    CompleteMSK(MSKD, MSK);

    MPK.kill();     
    
    RR::SetPrecision(PRECISION);

    mat_RR B;
    B.SetDims(4*N0, N0);

    int col, i, j, row;    
    row = 0;

    for(i = 0; i < N0; i++) {
        for(j = 0, col = 0; j < N0; j++, col++)
            B[row][col] = MSKD->B[i][j];
        row++;
    }//end-for

    for(i = 0; i < N0; i++) {
        for(j = N0, col = 0; j < 2*N0; j++, col++)
            B[row][col] = MSKD->B[i][j];
        row++;
    }//end-for

    for(i = N0; i < 2*N0; i++) {
        for(j = 0, col = 0; j < N0; j++, col++)
            B[row][col] = MSKD->B[i][j];
        row++;
    }//end-for

    for(i = N0; i < 2*N0; i++) {
        for(j = N0, col = 0; j < 2*N0; j++, col++)
            B[row][col] = MSKD->B[i][j];
        row++;
    }//end-for

    mat_RR BTilde; // Gram-Schmidt orthogonalization of basis B
    BlockGSO(BTilde, B, N0, PRECISION);

    mat_RR Block_BTilde;
    Block_BTilde.SetDims(2*N0, 2*N0);

    for(i = 0; i < N0; i++) {
        for(j = 0; j < N0; j++)
            Block_BTilde[i][j] = BTilde[i][j];
    }

    for(i = N0; i < 2*N0; i++) {
        for(j = 0; j < N0; j++)
            Block_BTilde[j][i] = BTilde[i][j];
    }

    for(i = N0; i < 2*N0; i++) {
        for(j = 0; j < N0; j++)
            Block_BTilde[i][j] = BTilde[N0+i][j];
    }

    for(i = N0; i < 2*N0; i++) {
        for(j = N0; j < 2*N0; j++)
            Block_BTilde[i][j] = BTilde[2*N0+i][j-N0];
    }

    NTL::clear(B);
    B.SetDims(2*N0, 2*N0);

    for(i = 0; i < 2*N0; i++)
        for(j = 0; j < 2*N0; j++)
            B[i][j] = MSKD->B[i][j];

    BTilde.kill();

    mat_RR innerp;
    vec_RR sigma, sigma_squared;

    PrecomputationRingKlein(sigma_squared, innerp, Block_BTilde, MSK[0], MSK[1]);

    mat_RR b;
    vec_RR X;
    RR eta, v, sigma0, tailcut;
    sigma0 = to_RR(3.195);
    tailcut = to_RR(13);

    PrecomputationRingPeikert(b, X, v, eta, innerp, sigma_squared, sigma0, 
                PRECISION, tailcut);

    vec_RR c, sample;

    /* Random center */
    long int sqr = (long int)sqrt(q0);
    c.SetLength(2*N0);
    for(i = 0; i < 2*N0; i++)
        c[i] = RandomBnd(sqr);

    sample = RingKlein(innerp, B, Block_BTilde, b, X, sigma_squared, c, 
                sigma0, eta, v, PRECISION);

    cout << "[>] Sample from lattice: " << sample << endl;

    B.kill();
    innerp.kill();
    sigma.kill(); sigma_squared.kill();
    b.kill(); X.kill();
    Block_BTilde.kill();
    sample.kill();

    delete(MSKD);
    
    return 0;
    
}//end-main()