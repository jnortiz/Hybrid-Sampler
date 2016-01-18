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

int main(int argc, char** argv) {

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

    int precision_id = 1;
    
    switch(precision_id) {
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
            
            mat_RR b;
            vec_RR Sigma, X;
            RR v;

            Sigma.SetLength(N0);
            
            for(i = 0; i < N0; i++)
                Sigma[i] = RandomBnd(q0*10);
            
            cout << "Sigma: " << Sigma << endl;
            
            long precision = 107;
            int m = 64;
            int tailcut = 13;
            
            clear(b);
            
            OfflineRing_Peikert(X, b, v, Block_BTilde, Sigma, precision, m, tailcut);

            X.kill();
            b.kill();
            
            B.kill();
            Block_BTilde.kill();
            Sigma.kill();

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

            cout << "\n/** FFTMul **/" << endl;
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
             * Euclidean division of polynomials
             */             

            vec_RR q, r, a_x, b_x;

        //    int k = NextPowerOfTwo(N0+1);
        //    int N = (int)pow(2.0, (float)k);
            a_x.SetLength(N0);
            b_x.SetLength(N0);

            int sr_q = sqrt(q0);
            
            for(int i = 0; i < N0; i++) {
                a_x[i] = RandomBnd(sr_q);
//                b_x[i] = RandomBnd(sr_q);
            }
            
            b_x = a_x;
//            a_x[0] = 1;
//            a_x[1] = 1;
//            a_x[2] = 1;
//            a_x[3] = 1;
//            a_x[4] = 1;
//            a_x[5] = 0;
//            a_x[6] = 0;
//            a_x[7] = 0;
//
//            b_x[0] = 1;
//            b_x[1] = 0;
//            b_x[2] = 0;
//            b_x[3] = 0;
//            b_x[4] = 1;
//            b_x[5] = 0;
//            b_x[6] = 0;
//            b_x[7] = 0;

            EuclideanDiv(q, r, a_x, b_x);

            cout << "a(x) = " << a_x << endl;
            cout << "b(x) = " << b_x << endl;
            cout << "q(x) = " << q << endl;
            cout << "r(x) = " << r << endl;

            q.kill(); r.kill();
            a_x.kill(); b_x.kill();
            
            ZZX test;
            test.SetLength(N0);
            test[0] = conv<ZZ>(10);
            
            cout << endl << test << endl;

            NTL::sub(test, conv<ZZ>(-4), test);
            
            cout << endl << test << endl;
            
            break;            
        }//end-case-3
        case 4: {
            
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
        case 5: {
            
            cout << endl;
            /*
             * XGCD
             */
            
            vec_RR a_x, b_x;
            a_x.SetLength(N0+1);
            b_x.SetLength(N0); // Cyclotomic polynomial
            
            a_x[0] = 1;
            a_x[N0] = 1;

            b_x[0] = 29;
            b_x[1] = -119;
            b_x[2] = -209;
            b_x[3] = -50;
            b_x[4] = 159;
            b_x[5] = -609;
            b_x[6] = -381;
            b_x[7] = 211;            
            
            vec_RR d, s, t;

            XGCD(d, s, t, a_x, b_x);

            cout << "\na(x) = " << a_x << endl;
            cout << "b(x) = " << b_x << endl;
            cout << "gcd(a,b) = " << d << endl;
            cout << "inv(a) = " << t << endl;
            cout << "inv(b) = " << s << endl;

            int q = 23;
            
            for(int i = 0; i < t.length(); i++)
                cout << (conv<int>(conv<float>(t[i])) % q) << " ";
            cout << endl;

            for(int i = 0; i < s.length(); i++)
                cout << (conv<int>(conv<float>(s[i])) % q) << " ";
            cout << endl;
            
            a_x.kill();
            b_x.kill();
            d.kill();
            s.kill();
            t.kill();
            
            break;
        }//end-case-5
        case 6: {
            
            vec_RR sigma, sqr;
            sigma.SetLength(N0);
            
            sigma[0] = 1;
            sigma[1] = 2;
            sigma[2] = 3;
            sigma[3] = 0;
            sigma[4] = -1;
            sigma[5] = -2;
            sigma[6] = 1;
            sigma[7] = 0;
            
            SquareRoot(sqr, sigma);
            
            sigma.kill(); sqr.kill();

            break;
        }
        default:
            break;      
    }//end-switch
    
    delete(MSKD);
    
    return 0;
    
}//end-main()