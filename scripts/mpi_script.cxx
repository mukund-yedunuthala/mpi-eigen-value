#include <cstdlib>
#include <iomanip>
#include <random>
#include <iostream>
#include <algorithm>
#include <mpi.h>
#include "linalg.h"

Matrix genLinkMatrix(unsigned int& size, int& seed, double& linkProbability) {
    std::mt19937 generator(seed);
    std::bernoulli_distribution distribution(linkProbability);
    // srand(seed);
    Matrix linkMatrix(size, size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            // Generate a random 0 or 1 based on the link probability
            linkMatrix.setValue(i, j, distribution(generator));
            // linkMatrix.setValue(i, j, rand() % 2);

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
    // std::cout<<"Normalized link Matrix:\n"; Q.print();
    // std::cout<<"Dangling nodes:\n"; d.print();
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
    // std::cout<<"Left stochastic matrix P:\n";modifiedLinkMatrix.print();
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
    MPI_Init(&argc, &argv);
    unsigned int size {};
    int rank{}, commSize {};
    double start_time, end_time;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    if (rank == 0) {start_time = MPI_Wtime();}
    
    size = std::atoi(argv[1]);
    int seed = 1;
    double linkProbability = 0.4;
    Matrix linkMatrix = genLinkMatrix(size, seed, linkProbability);
    // std::cout<<"Link Matrix:\n"; linkMatrix.print();
    Matrix P = modifyLinkMatrix(size, linkMatrix);
    // // std::cout<<"Modified Matrix:\n"; P.print();
    Vectors r(size);
    r = powerIterations(P, size, 100, 1e-10);
    // r.print();
    double theta_k = computeEigenValue(P,r);
    std::cout<< "Eigen value:\t" << theta_k << "\n";
    if (rank == 0) {
        end_time = MPI_Wtime();
        std::cout << "Time elapsed: " << std::setprecision(5);
        std::cout << end_time - start_time << "s\n";
    }

    MPI_Finalize();
    return  0;
}
