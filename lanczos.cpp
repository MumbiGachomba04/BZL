#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>

void solveTridiagonal(double &gse, std::vector<double> &gsv, const std::vector<double> &alpha, const std::vector<double> &beta, int nIters) {
    std::vector<double> T(nIters * nIters, 0.0);
    for (int i = 0; i < nIters; i++) {
        T[i * nIters + i] = alpha[i]; 
        if (i > 0) {
            T[i * nIters + (i - 1)] = beta[i - 1]; // Subdiagonal
            T[(i - 1) * nIters + i] = beta[i - 1]; // Superdiagonal
        }
    }

    // Power iteration to estimate the eigenvalue and eigenvector
    std::vector<double> eigenvec(nIters, 1.0); // Random initial vector
    double eigenval = 0.0, prevEigenval = 0.0;
    const double tol = 1e-10;

    for (int iter = 0; iter < 1000; iter++) {
        // T * eigenvec
        std::vector<double> temp(nIters, 0.0);
        for (int i = 0; i < nIters; ++i) {
            for (int j = 0; j < nIters; ++j) {
                temp[i] += T[i * nIters + j] * eigenvec[j];
            }
        }

        // Normalize
        double norm = std::sqrt(std::accumulate(temp.begin(), temp.end(), 0.0, [](double sum, double val) {
            return sum + val * val;
        }));
        
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

void lanczos() {
    int maxiter = 80;
    double convCrit = 1e-10;
    int nIters = 0;
    uint64_t nstates = 1 << 10; // Number of states (typically a power of 2)
    std::vector<double> alpha(maxiter, 0.0);
    std::vector<double> beta(maxiter, 0.0);
    std::vector<double> gsv(maxiter, 0.0);
    std::vector<double> r(nstates, 1.0 / std::sqrt(nstates)); // Initial guess for eigenvector
    std::vector<double> q(nstates, 0.0);
    std::vector<double> matrix(nstates * nstates, 0.0);
    double dold = std::numeric_limits<double>::max();
    double gse;

    // Fill matrix
    for (uint64_t i = 0; i < nstates; i++) {
        for (uint64_t k = 0; k < nstates; k++) {
            matrix[i * nstates + k] += std::sin((double)i + (double)k); 
        }
    }

     // Initialize q with the first residual vector
    for (uint64_t i = 0; i < nstates; i++) {
        q[i] = r[i];
    }

    // Lanczos iteration
    for (int n = 0; n < maxiter - 1; n++) {
  

        std::vector<double> temp(nstates, 0.0);
        for (uint64_t i = 0; i < nstates; i++) {
            for (uint64_t k = 0; k < nstates; k++) {
                temp[i] += matrix[i * nstates + k] * q[k];
            }
        }

        // Update alpha
        double tempalpha = 0.0;
        for (uint64_t i = 0; i < nstates; i++) {
            tempalpha += temp[i] * q[i];
        }
        alpha[n] = tempalpha;

   
        for (uint64_t i = 0; i < nstates; i++) {
            r[i] = temp[i] - alpha[n] * q[i];
            if (n > 0) {
                r[i] -= beta[n - 1] * r[i]; // Orthogonalize against previous r
            }
        }


        // Update beta
        double tempbeta = 0.0;
        for (uint64_t i = 0; i < nstates; i++) {
            tempbeta += r[i] * r[i];
        }
        beta[n] = std::sqrt(tempbeta);

        // Normalize r
        for (uint64_t i = 0; i < nstates; i++) {
            r[i] /= beta[n];
        }

        // Swap q and r for next iteration
        std::swap(q, r);

        // Solve the tridiagonal matrix for eigenvalues and eigenvectors
        nIters = n + 1; // Adjust for the 0-based indexing
        solveTridiagonal(gse, gsv, alpha, beta, nIters);

        std::cout << "Eigenvalue after iteration " << nIters << ": " << gse << std::endl;

        // Convergence check
        if (n > 5 && std::fabs(gse - dold) < convCrit ) {
            std::cout << "Converged after " << nIters<< " iterations!" << std::endl;
            break;
        }
        dold = gse;
    }
}


int main() {
    lanczos();
    return 0;
}