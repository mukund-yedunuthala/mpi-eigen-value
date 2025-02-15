/**
 * \file 
 * \brief Source file for the parallel implementation using MPI. 
 * \author Mukund Yedunuthala
 */

#include <cstdlib>
#include <iomanip>
#include <random>
#include <iostream>
#include <algorithm>
#include <mpi.h>
#include "linalg.h"

/**
 * \brief Generates a link matrix with random connections based on link probability.
 * 
 * This function initializes a square matrix of the given size, where each element
 * is randomly assigned 0 or 1 based on the given link probability. The random
 * values are generated using a Bernoulli distribution with a given seed.
 * 
 * \param size Reference to the size of the matrix (number of rows and columns).
 * \param seed Reference to the seed for the random number generator.
 * \param linkProbability Reference to the probability of a link (1) occurring.
 * \return A Matrix object representing the generated link matrix.
 */
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

/**
 * \brief Modifies the link matrix to create a transition probability matrix.
 * 
 * This function takes an existing link matrix and modifies it to ensure that it
 * represents a valid stochastic matrix. It constructs the matrices Q, d, e, and edT
 * to handle cases where nodes have no outgoing links.
 * 
 * \param size The size of the matrix (number of rows and columns).
 * \param linkMatrix Reference to the input link matrix.
 * \return A Matrix object representing the modified transition probability matrix.
 */
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

/**
 * \brief Computes the eigenvector using power iteration.
 * 
 * This function applies the power iteration method to approximate the dominant 
 * eigenvector of the given matrix P. The process iterates until the difference
 * between successive approximations is below the given tolerance or the maximum
 * number of iterations is reached.
 * 
 * \param P Reference to the transition probability matrix.
 * \param rank The rank of the process.
 * \param commSize The number of processes present in the swarm.
 * \param size The size of the matrix (number of rows and columns).
 * \param max_it The maximum number of iterations.
 * \param tol The convergence tolerance.
 * \return A Vectors object representing the eigenvector.
 */
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
                qk_i[i] += (recvBuffer[i*size+j]*rk.getDataAtIndex(j));
            }
        }
        // Local sum
        double qk_i_sum{};
        for (int i = 0; i<rowsToSend; i++) {
            qk_i_sum += std::abs(qk_i[i]) ;
        }
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
        double max_diff {};
        for (int i = 0; i<size; i++) {
            double diff = rk1.getDataAtIndex(i)- rk.getDataAtIndex(i);
            max_diff = std::max(max_diff,std::abs(diff));
        }
        std::cout << std::setprecision(6);
        if (rank==ROOT) std::cout << "Iteration: " << it+1 << "\t Diff: " << max_diff << "\n";
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

/**
 * \brief Computes the eigenvalue of the matrix P using Rayleigh quotient.
 * 
 * \param P Reference to the transition probability matrix.
 * \param r Reference to the dominant eigenvector obtained from power iteration.
 * \return The computed eigenvalue.
 */
double computeEigenValue(Matrix& P, Vectors& r) {
    Vectors temp = P.dot(r);
    double lbda = (r.dot(temp))/(r.norm()*r.norm());
    return lbda;
}

/**
 * \brief The main function to compute the eigenvalue of a generated matrix.
 * 
 * The program takes the size of the matrix as a command-line argument, generates a
 * link matrix with random connections, modifies it into a transition probability matrix,
 * and then applies the power iteration method to compute the eigenvector.
 * Finally, the eigenvalue is calculated and displayed along with the elapsed time.
 * 
 * \param argc The number of command-line arguments.
 * \param argv The command-line arguments (expects the matrix size as an argument).
 * \return An integer indicating the exit status (0 for success).
 */
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
