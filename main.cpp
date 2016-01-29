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

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp() {
    struct timespec now;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &now);
    return now.tv_nsec + (timestamp_t)now.tv_sec * 1000000000.0;
}//end-get_timestamp()

int main() {

    srand(rdtsc());
    
    ZZX MSK[4];
    ZZ_pX MPK;
    
    MSK_Data *MSKD = new MSK_Data;
    timestamp_t average, ts_start, ts_end;
    
    ZZ_p::init(q1);
    zz_p::init(q0);

    cout << "[!] Parameters (q, N0) = (" << q0 << ", " << N0 << ")" << endl;
    int it, n_iterations;
    n_iterations = 1;
    
    average = 0.0;
    for(it = 0; it < n_iterations; it++) {
        
        ts_start = get_timestamp();
        Keygen(MPK, MSK);
        CompleteMSK(MSKD, MSK);
        ts_end = get_timestamp();

        average += (ts_end - ts_start);
        
    }//end-for

    MPK.kill();     
    
    cout << "[!] KeyGen average running time for " << n_iterations << " iterations: " << (float)(average/((float)(n_iterations)*1000000000.0)) << " s." << endl;           

    int function_id = 4;
    
    switch(function_id) {
        case 1: {
            /*
             * Ring_Klein Gaussian sampler 
             */            
            
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

            average = 0.0;
            n_iterations = 1;
            
            for(it = 0; it < n_iterations; it++) {

                ts_start = get_timestamp();
                BlockGSO(BTilde, B, N0, PRECISION);
                ts_end = get_timestamp();

                average += (ts_end - ts_start);

            }//end-for

            cout << "[!] Block-GSO average running time for " << n_iterations << " iterations: " << (float)(average/((float)(n_iterations)*1000000000.0)) << " s." << endl;      

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
            
        #ifdef PRINT_ALL
            for(i = 0; i < 4*N0; i++) {
                if(i % N0 == 0)
                    cout << endl;
                for(j = 0; j < N0; j++)
                    cout << BTilde[i][j] << " ";
                cout << endl;
            }//end-for
        #endif            
            
//            vec_RR sigma;
//            sigma.SetLength(N0);
//            
//            for(i = 0; i < N0; i++)
//                sigma[i] = RandomBnd(sqrt(q0));
//            
//            RR sigma0 = to_RR(3.145);
//            
//            OfflineRingPeikert(sigma0, sigma, Block_BTilde);
//
//            sigma.kill();            
            B.kill();
            Block_BTilde.kill();

            break;
        }//end-case-1
        
        case 2: {
            
            cout << endl;    
            /*
             * Polynomial multiplication
             */

            CC out[2*N0];
            vec_RR a_, b_;
            a_.SetLength(N0);
            b_.SetLength(N0);

            a_[0] = -1;
            a_[1] = 5;

            b_[2] = -2;

            FFTMul(out, a_, b_);

            cout << "/** FFTMul **/" << endl;
            for(int i = 0; i < 2*N0; i++)
                cout << out[i] << " ";
            cout << endl;

            a_.kill();
            b_.kill();
            
            break;
        }//end-case-2
        case 3: {
            
            cout << endl;    
            /*
             * Modular polynomial multiplication
             */

            vec_RR a_x, b_x;
            CC mult[N0];
            
            cout << "[>] Multiplication of a(x)*b(x) % phi:" << endl;
            FFTMulMod(mult, a_x, b_x);

            for(int i = 0; i < N0; i++)
                cout << mult[i] << " ";
            cout << endl;
            
            a_x.kill(); b_x.kill();
            
            break;
        }
        case 4: {
            
            RR::SetOutputPrecision(4);
            RR::SetPrecision(150);
            
            vec_RR a, b, innerp;
            a.SetLength(2*N0); b.SetLength(2*N0);
            innerp.SetLength(N0);
            
            /* Defining a(x) and b(x) as random binary polynomials */
            for(int i = 0; i < (2*N0); i++) {
                a[i] = NTL::random_RR();
                b[i] = NTL::random_RR();
            }//end-for                
            
            cout << "\na in H: " << a << endl;
            cout << "b in H: " << b << endl;

            /* Computing the inner product of a(x) and b(x) \in H */
            InnerProduct(innerp, a, b);
            cout << endl << "<a, b>_K: " << innerp << endl;
            innerp.kill();

            /* Reducing a(x) and b(x) into K */
            a.SetLength(N0);
            b.SetLength(N0);
            
            cout << "\na in K: " << a << endl;
            cout << "b in K: " << b << endl;
            
            /* Trying to invert a(x) in K */
            vec_RR inverse;
            CC one[N0];
            inverse.SetLength(N0);            
            Inverse(inverse, a);            
            cout << endl << "a^{-1}: " << inverse << endl;
            
            /* Verifying if inversion was successful */
            FFTMulMod(one, a, inverse);
            inverse.kill();            
            cout << "a times a^{-1}: ";
            for(int i = 0; i < N0; i++)
                cout << one[i] << " ";
            cout << endl;

            /* Square-root in K */
            vec_RR sqr;
            SquareRoot(sqr, a);
            a.kill();
            FFTMulMod(one, sqr, sqr);
            cout << endl << "sqrt(a): " << sqr << endl;
            
            cout << "a in K: ";
            for(int i = 0; i < N0; i++)
                cout << one[i].real() << " ";
            cout << endl;
            
            sqr.kill();
            
            break;
        }
        default:
            break;      
    }//end-switch
    
    delete(MSKD);
    
    return 0;
    
}//end-main()