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
Vectors parallelPowerIterations(Matrix& P, int &rank, int& commSize, unsigned int& size, int max_it, double tol) {
    int ROOT {};
    unsigned int rowsToSend {size/commSize};
    unsigned long int bufferSize {rowsToSend * size};
    Vectors rk (size), rk1 (size);
    rk.setToValue(1.0/size);

    double* sendBuffer = new double[size*size]{};
    double* recvBuffer = new double[rowsToSend*size]{};
    for (unsigned long int i = 0; i < size; i++) {
        for (unsigned long int j = 0; j < size; j++) {
            sendBuffer[i*size+j] = P.getValue(i,j);
        }
    }
    MPI_Scatter(
        sendBuffer, bufferSize, MPI_DOUBLE, 
        recvBuffer, bufferSize, MPI_DOUBLE, 
        ROOT, MPI_COMM_WORLD
    );
    double* rk_i = new double[rowsToSend]{};
    double* qk_i = new double[rowsToSend]{};
    double* resBuffer = new double[size]; 

    
    for (int it = 0; it < max_it; it++) {
        
        MPI_Scatter(
            rk.getData(), rowsToSend, MPI_DOUBLE, 
            rk_i, rowsToSend, MPI_DOUBLE, 
            ROOT, MPI_COMM_WORLD
        );
        
        for (unsigned int i = 0; i < rowsToSend; i++) {
            qk_i[i] = 0;
            for (unsigned int j = 0; j < size; j++) {
                qk_i[i] += (recvBuffer[i*size+j]*rk_i[j]);
            }
        }
        // Local sum
        double qk_i_sum{};
        for (int i = 0; i<rowsToSend; i++) {
            qk_i_sum += std::abs(qk_i[i]) ;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        double L1_qk{};
        MPI_Reduce(
            &qk_i_sum, &L1_qk, 1, 
            MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD
        );
        MPI_Bcast(&L1_qk, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
        for (int i = 0; i<rowsToSend; i++) {
            rk_i[i] =  qk_i[i]/L1_qk ;
        }
        
        MPI_Allgather(
            rk_i, rowsToSend, MPI_DOUBLE,
            resBuffer, rowsToSend, MPI_DOUBLE,
            MPI_COMM_WORLD
        );

        rk1.setData(resBuffer);
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
    delete [] resBuffer;
    delete [] qk_i;
    delete [] rk_i;
    delete [] recvBuffer;
    delete [] sendBuffer;
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
    Matrix P = modifyLinkMatrix(size, linkMatrix);


    Vectors r(size);
    r = parallelPowerIterations(P, rank, commSize, size, 50, 1e-10);

    if (rank == 0){ 
        double theta_k = computeEigenValue(P,r);
        std::cout<< "Eigen value:\t" << theta_k << "\n";
        end_time = MPI_Wtime();
        std::cout << "Time elapsed: " << std::setprecision(5);
        std::cout << end_time - start_time << "s\n";
    }

    MPI_Finalize();
    return  0;
}
