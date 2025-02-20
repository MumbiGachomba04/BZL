#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <chrono>
#include <limits>
#ifdef _OMP
#include <omp.h>
#endif

#ifdef _RAVETRACE
#include "sdv_tracing.h"
#endif

#define PI 3.14159265358979323846


inline double fast_sin(double x) {
    const double B = 4.0 / PI;
    const double C = -4.0 / (PI * PI);
    return B * x + C * x * std::abs(x);
}

void solveTridiagonal(double &gse, std::vector<double> &gsv, const std::vector<double> &alpha, const std::vector<double> &beta, int nIters) {
    std::vector<double> T(nIters * nIters, 0.0);
    T[0]=alpha[0];
    for (int i = 1; i < nIters; i++) {
            T[i * nIters + i] = alpha[i];  
            T[i * nIters + (i - 1)] = beta[i - 1]; // Subdiagonal
            T[(i - 1) * nIters + i] = beta[i - 1]; // Superdiagonal
    }
   

    // Power iteration to estimate the eigenvalue and eigenvector
    std::vector<double> eigenvec(nIters, 1.0); // Random initial vector
    double eigenval = 0.0, prevEigenval = 0.0;
    const double tol = 1e-10;

    for (int iter = 0; iter < 1000; iter++) {
        // T * eigenvec
        std::vector<double> temp(nIters, 0.0);
    double cumulate;
#ifdef _OPENMPWORKSHARING
#pragma omp parallel for reduction(+:cumulate)
#endif
#ifdef _OPENMPSIMD
#pragma omp simd reduction(+:cumulate)
#endif  
        for (int i = 0; i < nIters; ++i) {
            cumulate=0.0;
            for (int j = 0; j < nIters; ++j) {
                cumulate += T[i * nIters + j] * eigenvec[j];
            }
            temp[i] = cumulate;
        }

        // Normalize
        double norm = std::sqrt(std::accumulate(temp.begin(), temp.end(), 0.0, [](double sum, double val) {
            return sum + val * val;
        }));
#ifdef _OPENMPWORKSHARING
#pragma omp parallel for  
#endif      
#ifdef _OPENMPSIMD
#pragma omp simd 
#endif 
        for (int i = 0; i < nIters; ++i) {
            temp[i] /= norm;
        }

        // Estimate eigenvalue
        eigenval = std::inner_product(temp.begin(), temp.end(), eigenvec.begin(), 0.0);

        // Check for convergence
        if (std::fabs(eigenval - prevEigenval) < tol) {
            break;
        }

        prevEigenval = eigenval;
        eigenvec = temp;
    }

    gse = eigenval; // The estimated eigenvalue
    gsv = eigenvec; // The corresponding eigenvector
}

void lanczos(int maxiter,int size_basetwo) {
    
    double convCrit = 1e-10;
    int nIters = 0;
    uint64_t nstates = 1 << size_basetwo; // Number of states (typically a power of 2)
    std::vector<double> alpha(maxiter, 0.0);
    std::vector<double> beta(maxiter, 0.0);
    std::vector<double> gsv(maxiter, 0.0);
    std::vector<double> r(nstates, 1.0 / std::sqrt(nstates)); // Initial guess for eigenvector
    std::vector<double> q(nstates, 0.0);
    std::vector<double> matrix(nstates * nstates, 0.0);
    double dold = std::numeric_limits<double>::max();
    double gse;
    int threads=1;

#ifdef _RAVETRACE
    int values[] = {0,1,2,3,4,5,6,7,8};
    const char * v_names[] = {"Other","fill_matrix","init_q","mat_vec","update_alpha","orthoganalise","update_beta","normalise","tridiagonal"};
    trace_name_event_and_values(1000, "code_region", 9, values, v_names);
    trace_init();
#endif

    // Fill matrix
#ifdef _OPENMPWORKSHARING
#pragma omp parallel 
{
    #pragma omp single
    threads = omp_get_num_threads();
    
}
#pragma omp parallel for
#endif
#ifdef _RAVETRACE
trace_event_and_value(1000,1);
#endif
    for (uint64_t i = 0; i < nstates; i++) {
#ifdef _OPENMPSIMD
#pragma omp simd 
#endif
        for (uint64_t k = 0; k < nstates; k++) {
            matrix[i * nstates + k] += fast_sin((double)i + (double)k); 
        }
    }
#ifdef _RAVETRACE
trace_event_and_value(1000,0);
#endif
     // Initialize q with the first residual vector
#ifdef _OPENMPWORKSHARING
#pragma omp parallel for 
#endif 
#ifdef _OPENMPSIMD
#pragma omp simd 
#endif
#ifdef _RAVETRACE
trace_event_and_value(1000,2);
#endif
    for (uint64_t i = 0; i < nstates; i++) {
        q[i] = r[i];
    }
#ifdef _RAVETRACE
trace_event_and_value(1000,0);
#endif
    double time=0.0;
auto start = std::chrono::high_resolution_clock::now();
    // Lanczos iteration
    int n;
    for (n = 0; n < maxiter - 1; n++) {
  
    // Matvec multiplication
        std::vector<double> temp(nstates, 0.0);
        double cumulate ;
#ifdef _OPENMPWORKSHARING
#pragma omp parallel for schedule(runtime) reduction(+:cumulate)
#endif 
#ifdef _OPENMPSIMD
#pragma omp simd reduction(+:cumulate)
#endif        
#ifdef _RAVETRACE
trace_event_and_value(1000,3);
#endif
        for (uint64_t i = 0; i < nstates; i++) {
            cumulate =0.0;
            for (uint64_t k = 0; k < nstates; k++) {
                cumulate += matrix[i * nstates + k] * q[k];
            }
            temp[i]=cumulate;
        }
#ifdef _RAVETRACE
trace_event_and_value(1000,0);
#endif

        // Update alpha
        double tempalpha = 0.0;
#ifdef _OPENMPWORKSHARING
#pragma omp parallel for reduction(+ : tempalpha) 
#endif 
#ifdef _OPENMPSIMD
#pragma omp simd reduction(+:tempalpha)
#endif  
#ifdef _RAVETRACE
trace_event_and_value(1000,4);
#endif
        for (uint64_t i = 0; i < nstates; i++) {
            tempalpha += temp[i] * q[i];
        }
        alpha[n] = tempalpha;
#ifdef _RAVETRACE
trace_event_and_value(1000,0);
#endif
#ifdef _OPENMPWORKSHARING
#pragma omp parallel for 
#endif 
#ifdef _OPENMPSIMD
#pragma omp simd 
#endif 
#ifdef _RAVETRACE
trace_event_and_value(1000,5);
#endif	
        for (uint64_t i = 0; i < nstates; i++) {
            r[i] = temp[i] - alpha[n] * q[i];
            if (n > 0) {
                r[i] -= beta[n - 1] * r[i]; // Orthogonalize against previous r
            }
        }
#ifdef _RAVETRACE
trace_event_and_value(1000,0);
#endif

        // Update beta
        double tempbeta = 0.0;
#ifdef _OPENMPWORKSHARING
#pragma omp parallel for reduction(+ : tempbeta)
#endif  
#ifdef _OPENMPSIMD
#pragma omp simd reduction(+ : tempbeta) 
#endif 
#ifdef _RAVETRACE
trace_event_and_value(1000,6);
#endif
        for (uint64_t i = 0; i < nstates; i++) {
            tempbeta += r[i] * r[i];
        }
        beta[n] = std::sqrt(tempbeta);
#ifdef _RAVETRACE
trace_event_and_value(1000,0);
#endif
        // Normalize r
#ifdef _OPENMPWORKSHARING
#pragma omp parallel for 
#endif 
#ifdef _OPENMPSIMD
#pragma omp simd 
#endif
#ifdef _RAVETRACE
trace_event_and_value(1000,7);
#endif	
        for (uint64_t i = 0; i < nstates; i++) {
            r[i] /= beta[n];
        }
#ifdef _RAVETRACE
trace_event_and_value(1000,0);
#endif

        // Swap q and r for next iteration
        std::swap(q, r);

        // Solve the tridiagonal matrix for eigenvalues and eigenvectors
        nIters = n + 1; // Adjust for the 0-based indexing
#ifdef _RAVETRACE
trace_event_and_value(1000,8);
#endif        
        solveTridiagonal(gse, gsv, alpha, beta, nIters);
#ifdef _RAVETRACE
trace_event_and_value(1000,0);
#endif
        std::cout << "Eigenvalue after iteration " << nIters << ": " << gse << std::endl;

        // Convergence check
        if (n > 5 && std::fabs(gse - dold) < convCrit ) {
            std::cout << "Converged after " << nIters<< " iterations!" << std::endl;
            break;
        }
        dold = gse;
    }
     auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    time= elapsed.count();
     std::cout << "Number of threads: " << threads << std::endl;  
     std::cout << "Elapsed time: " << time << " seconds" << std::endl;
     std::cout << "Throughput: " << n/time << " iters/second" << std::endl;
}


int main(int argc, char *argv[]) {
    int maxiter, size_basetwo;
    if (argc != 3) {
      std::cout<< "Error: parsing command line arguments"<< std::endl;
      std::cout<< "./executable <MaximumIterations> <Matrix Size (base2)>"<< std::endl;
      exit(EXIT_FAILURE);
   }
    
    maxiter = atoi(argv[1]);
    size_basetwo = atoi(argv[2]);

    lanczos(maxiter,size_basetwo);
    return 0;
}
