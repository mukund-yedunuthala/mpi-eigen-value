#include <cstdlib>
#include <iomanip>
#include <random>
#include <iostream>
#include <algorithm>
#include <chrono>
#include "linalg.h"

Matrix genLinkMatrix(unsigned int& size, int& seed, double& linkProbability) {
    std::mt19937 generator(seed);
    std::bernoulli_distribution distribution(linkProbability);
    Matrix linkMatrix(size, size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            // Generate a random 0 or 1 based on the link probability
            linkMatrix.setValue(i, j, distribution(generator));
        }
    }
    return linkMatrix;
}

Matrix modifyLinkMatrix(unsigned int size, Matrix& linkMatrix) {
    Matrix modifiedLinkMatrix(size, size); //P
    Vectors nj (size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            double nji = nj.getDataAtIndex(i);
            nji += linkMatrix.getValue(j,i);
            nj.setValueAtIndex(i,nji);
        }
    }
    // Construct Q, d and e
    Matrix Q = linkMatrix;
    Vectors d (size);
    Vectors e (size); e.setToValue(1.);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (nj.getDataAtIndex(j)>0) {
                Q.setValue(i,j,(Q.getValue(i,j)/nj.getDataAtIndex(j)));
            }
            else if (nj.getDataAtIndex(j)==0) {
                d.setValueAtIndex(j, 1.);
            }
        }
    }
    // edT
    Matrix edT (size, size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            edT.setValue(i,j,
                (d.getDataAtIndex(j)*e.getDataAtIndex(i))/size
            );
        }
    }
    // adding up
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            modifiedLinkMatrix.setValue(i,j,
                Q.getValue(i,j) + edT.getValue(i,j)
            );
        }
    }
    return modifiedLinkMatrix;
}

Vectors powerIterations(Matrix& P, unsigned int size, int max_it, double tol) {
    Vectors rk (size), rk1 (size), qk (size);
    rk.setToValue(1.0/size);
    for (int it = 0; it < max_it; it++) {
        qk = P.dot(rk);
        for (int i = 0; i<size; i++) {
            rk1.setValueAtIndex(i, qk.getDataAtIndex(i)/qk.abssum());
        }
        double max_diff = 0;
        for (int i = 0; i<size; i++) {
            double diff = rk1.getDataAtIndex(i)- rk.getDataAtIndex(i);
            max_diff = std::max(max_diff,std::abs(diff));
        }
        std::cout << "Iteration: " << it+1 << "\t Diff: " << max_diff << "\n";
        if (max_diff < tol) {
            break;
        }
        rk = rk1;
    }
    return rk;
}

double computeEigenValue(Matrix& P, Vectors& r) {
    Vectors temp = P.dot(r);
    double lbda = (r.dot(temp))/(r.norm()*r.norm());
    return lbda;
}

int main(int argc, char **argv) {
    unsigned int size {};
    size = std::atoi(argv[1]);
    auto start = std::chrono::high_resolution_clock::now();
    int seed = 135356;
    double linkProbability = 0.7;
    Matrix linkMatrix = genLinkMatrix(size, seed, linkProbability);
    Matrix P = modifyLinkMatrix(size, linkMatrix);
    Vectors r(size);
    r = powerIterations(P, size, 100, 1e-10);
    double theta_k = computeEigenValue(P,r);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout<< "Eigen value:\t" << theta_k << "\n";
    std::cout << "Time elapsed: " << std::setprecision(5);
    std::cout << duration.count() << "ms\n";

    return  0;
}
