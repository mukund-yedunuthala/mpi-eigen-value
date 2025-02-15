/**
 * \file 
 * \brief Source file for the sequential implementation. 
 * \author Mukund Yedunuthala
 */

#include <cstdlib>
#include <iomanip>
#include <random>
#include <iostream>
#include <algorithm>
#include <chrono>
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
    Matrix linkMatrix(size, size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            // Generate a random 0 or 1 based on the link probability
            linkMatrix.setValue(i, j, distribution(generator));
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
 * This function applies the power iteration method to approximate the  
 * eigenvector of the given matrix P. The process iterates until the difference
 * between successive approximations is below the given tolerance or the maximum
 * number of iterations is reached.
 * 
 * \param P Reference to the transition probability matrix.
 * \param size The size of the matrix (number of rows and columns).
 * \param max_it The maximum number of iterations.
 * \param tol The convergence tolerance.
 * \return A Vectors object representing the eigenvector.
 */
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
