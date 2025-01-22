#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>

void solveTridiagonal(double &gse, std::vector<double> &gsv, const std::vector<double> &alpha, const std::vector<double> &beta, int nIters) {
    // Construct the tridiagonal matrix
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
    int maxiter = 500;
    int j = -1; // Start with -1 so that the first increment brings it to 0
    double convCrit = 1e-6;
    int nIters = 0;
    uint64_t nstates = 1 << 8; // Number of states (typically a power of 2)
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

    // Lanczos iteration
    for (int n = 0; n < maxiter - 1; n++) {
        // Increment j at the start of each iteration
        j++;

        // Orthogonalize q and r
        if (j != 0) {
            for (uint64_t i = 0; i < nstates; i++) {
                double t = r[i];
                r[i] = q[i] / beta[j - 1];
                q[i] = -beta[j - 1] * t;
            }
        }

        // Matrix-vector multiplication
        std::fill(q.begin(), q.end(), 0.0);
        for (uint64_t i = 0; i < nstates; i++) {
            for (uint64_t k = 0; k < nstates; k++) {
                q[i] += matrix[i * nstates + k] * r[k];
            }
        }

        // Update alpha
        double tempalpha = 0.0;
        for (uint64_t i = 0; i < nstates; i++) {
            tempalpha += r[i] * q[i];
        }
        alpha[j] = tempalpha;

        // Adjust q
        for (uint64_t i = 0; i < nstates; i++) {
            q[i] -= alpha[j] * r[i];
        }

        // Update beta
        double tempbeta = 0.0;
        for (uint64_t i = 0; i < nstates; i++) {
            tempbeta += q[i] * q[i];
        }
        beta[j] = std::sqrt(tempbeta);

        // Solve the tridiagonal matrix for eigenvalues and eigenvectors
        nIters = j + 1; // Adjust for the 0-based indexing
        solveTridiagonal(gse, gsv, alpha, beta, nIters);

        std::cout << "Eigenvalue after iteration " << n + 1 << ": " << gse << std::endl;

        // Compute the residual
        std::vector<double> residual(nstates, 0.0);
        for (uint64_t i = 0; i < nstates; i++) {
            for (uint64_t k = 0; k < nstates; k++) {
                residual[i] += matrix[i * nstates + k] * gsv[k];
            }
            residual[i] -= gse * gsv[i]; // Subtract eigenvalue * eigenvector
        }

        double residual_norm = std::sqrt(std::inner_product(residual.begin(), residual.end(), residual.begin(), 0.0));
        std::cout << "Residual norm after iteration " << n + 1 << ": " << residual_norm << std::endl;

        // Convergence check
        if (n > 5 && std::fabs(gse - dold) < convCrit && residual_norm < convCrit) {
            std::cout << "Converged after " << n + 1 << " iterations!" << std::endl;
            break;
        }
        dold = gse;
    }
}


int main() {
    lanczos();
    return 0;
}
