#include <cstdlib>
#include <iomanip>
#include <random>
#include <iostream>
#include <algorithm>
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
    Matrix modifiedLinkMatrix(size, size);
    for (int i = 0; i < size; ++i) {
        int sum = 0;
        for (int j = 0; j < size; ++j) {
            sum = sum + linkMatrix.getValue(j,i);
        }
        for (int j = 0; j < size; ++j) {
            if (sum != 0) {
                modifiedLinkMatrix.setValue(j, i, linkMatrix.getValue(j, i)/sum);
            }
            else {
                modifiedLinkMatrix.setValue(j,i,1.0/size);
            }
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
            rk1.setValueAtIndex(i, qk.getDataAtIndex(i)/qk.sum());
        }
        double max_diff = 0;
        for (int i = 0; i<size; i++) {
            double diff = rk1.getDataAtIndex(i)- rk.getDataAtIndex(i);
            max_diff = std::max(max_diff,std::abs(diff));
        }
        std::cout << std::setprecision(6);
        std::cout << "Iteration: " << it+1 << "\t Diff: " << max_diff << "\n";
        if (max_diff < tol) {
            break;
        }
        rk1.print();
        rk = rk1;
    }
    return rk;
}

double computeEigenValue(Matrix& P, Vectors& r) {
    Vectors temp = P.dot(r);
    double lbda = (r.dot(temp))/(r.norm()*r.norm());
    return lbda;
}

int main() {
    unsigned int size {6};
    int seed = 322;
    double linkProbability = 0.5;
    Matrix linkMatrix = genLinkMatrix(size, seed, linkProbability);
    std::cout<<"Link Matrix:\n";
    linkMatrix.print();
    Matrix P = modifyLinkMatrix(size, linkMatrix);
    std::cout<<"Modified Matrix:\n";
    P.print();
    Vectors r(size);
    r = powerIterations(P, size, 30, 1e-7);
    double theta_k = computeEigenValue(P,r);
    std::cout<< "Eigen value:\t" << theta_k << "\n";
    return  0;
}
