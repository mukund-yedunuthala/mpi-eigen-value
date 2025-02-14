#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <mpi.h>

using namespace std;

// Generate random matrix
vector<vector<double>> rand_matrix(int n) {
    vector<vector<double>> random_matrix(n, vector<double>(n));

    srand(1); // Seed for reproducibility

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            random_matrix[i][j] = rand() % 2; // Random values of 0 or 1
        }
    }

    return random_matrix;
}

// Vector and vector multiplication
void localP_r_multiply(const vector<double>& local_matrix, const vector<double>& r,
                       vector<double>& local_r_new, int process_rows, int n) {
    for (int i = 0; i < process_rows; ++i) {
        local_r_new[i] = 0;
        for (int j = 0; j < n; ++j) {
            local_r_new[i] += local_matrix[i * n + j] * r[j]; // Multiply each local matrix row with the global vector r
        }
    }
}

// Matrix and vector multiplication using a flattened local matrix
void matrix_vec_multiplication(const vector<double>& local_matrix, vector<double>& r, vector<double>& result_vec, int process_rows, int n) {
    for (int i = 0; i < process_rows; ++i) {
        result_vec[i] = 0;
        for (int j = 0; j < n; ++j) {
            result_vec[i] += local_matrix[i * n + j] * r[j]; // Matrix-vector multiplication
        }
    }
}

// Normalizing - L1 norm
void normalize_vec(vector<double>& vec, int n) {
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        sum += vec[i];
    }
    for (int i = 0; i < n; ++i) {
        vec[i] /= sum;
    }
}

double dot_product(const vector<double>& rk, const vector<double>& P_rk, int n) {
    double result = 0.0;
    for (int i = 0; i < n; ++i) {
        result += rk[i] * P_rk[i];
    }
    return result;
}

// L2 norm
double norm_2(const vector<double>& rk, int n) {
    double sum_of_squares = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_of_squares += rk[i] * rk[i];
    }
    return sqrt(sum_of_squares);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int max_iters = 30;
    int n;
    double start_time, end_time, time_taken;

    if (rank == 0) {
        start_time = MPI_Wtime();
        // cout << "Enter the size of the matrix (n): ";
        // cin >> n;
        n = 5;
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<double> d(n, 0);
    vector<double> e(n, 1.0);
    vector<vector<double>> P(n, vector<double>(n, 0));
    vector<double> flat_P(n * n, 0);

    // Generate the random matrix
    vector<vector<double>> L = rand_matrix(n);
   
    // Check for zero rows and construct the vector d
    for (int i = 0; i < n; ++i) {
        bool isZeroRow = true;
        for (int j = 0; j < n; ++j) {
            if (L[i][j] != 0) {
                isZeroRow = false;
                break;
            }
        }
        d[i] = isZeroRow ? 1 : 0;
    }

    // Compute d^T * e (outer product)
    vector<vector<double>> deT(n, vector<double>(n, 0)); // Matrix to store d^T * e
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            deT[i][j] = d[i] * e[j]; // Outer product d^T * e
        }
    }

    // Construct P matrix: P = Q + d^T * e
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            P[i][j] = L[i][j] + deT[i][j]; // Add Q and d^T * e
        }
    }


    // Calculate the outgoing link count for each page (n_j)
    vector<int> outgoing_links(n, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            outgoing_links[i] += P[i][j];  // Count outgoing links from page i
        }
    }

    // Construct matrix Q and handle dangling nodes (no outgoing links)
    for (int i = 0; i < n; ++i) {
        if (outgoing_links[i] > 0) {
            for (int j = 0; j < n; ++j) {
                P[i][j] /= outgoing_links[i];
            }
        }
    }
    
    // Flatten the matrix P for scattering
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            flat_P[i * n + j] = P[i][j];
        }
    }
        /*
        if (rank == 0) {
        cout << "P matrix is:\n";
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                cout << P[i][j] << " ";
            }
            cout << endl;
        }
    }*/

    // Calculate the number of rows each process will handle
    int rows_per_process = n / size; // Base number of rows
    int extra_rows = n % size; // Extra rows for some processes

    // Prepare to scatter the matrix
    vector<int> sendcounts(size);
    vector<int> displs(size);

    for (int i = 0; i < size; ++i) {
        sendcounts[i] = rows_per_process * n; // Each row has n elements
        if (i < extra_rows) {
            sendcounts[i] += n; // Distribute extra rows
        }
        displs[i] = (i == 0) ? 0 : displs[i - 1] + sendcounts[i - 1];
    }

    // Each process will receive its corresponding subset of flat_P
    int recv_size = sendcounts[rank];
    vector<double> local_matrix(recv_size);

    // Scatter the matrix
    MPI_Scatterv(flat_P.data(), sendcounts.data(), displs.data(), MPI_DOUBLE,
                 local_matrix.data(), recv_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /*
    // Now each process has its own portion of the matrix in local_matrix
    cout << "Process " << rank << " received elements: ";
    for (int i = 0; i < recv_size; ++i) {
        cout << local_matrix[i] << " ";
    }
    cout << endl;*/

    vector<double> r(n, 1.0 / n);
    vector<double> local_r_new(rows_per_process + (rank < extra_rows ? 1 : 0), 0.0); // Local result for each process
    vector<double> r_new(n, 0.0);  // Gathered result on rank 0

    MPI_Bcast(r.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int iter = 0; iter < max_iters; ++iter) {
        // Step 1: Multiply the matrix P by the current rank vector r
        localP_r_multiply(local_matrix, r, local_r_new, rows_per_process + (rank < extra_rows ? 1 : 0), n);

        /*
    // Print the local_r_new for each process
    cout << "Process " << rank << " local_r_new: ";
    for (double value : local_r_new) {
        cout << value << " ";
    }
    cout << endl;
        */
        // Gather the results from all processes
        vector<int> send_counts(size);
        vector<int> displ_r(size);

        // Calculate send_counts and displ_r
        for (int i = 0; i < size; ++i) {
                send_counts[i] = rows_per_process + (i < extra_rows ? 1 : 0); // Size of local_r_new for each process
                displ_r[i] = (i == 0) ? 0 : displ_r[i - 1] + send_counts[i - 1]; // Displacement
        }

        MPI_Gatherv(local_r_new.data(), local_r_new.size(), MPI_DOUBLE, // Local results from each process
            r_new.data(), send_counts.data(), displ_r.data(), MPI_DOUBLE,
            0, MPI_COMM_WORLD);

        if (rank == 0) {
                        normalize_vec(r_new, n);
                        }
                        // Broadcast the normalized vector to all processes
                MPI_Bcast(r_new.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        r = r_new;  // Update r for the next iteration

        }
        /*
        if (rank == 0){
                cout << "Root received r_new: ";
        for (double value : r_new) {
             cout << value << " ";
        }
        cout << endl;
        }*/

        // Step 4: Compute P*r_k locally
    vector<double> local_P_rk(local_r_new.size(), 0.0);  // Local result for P*r_k

    localP_r_multiply(local_matrix, r, local_P_rk, rows_per_process + (rank < extra_rows ? 1 : 0), n);
        /*
    // Print the local_P_rk for each process
    cout << "Process " << rank << " local_P_rk: ";
    for (double value : local_P_rk) {
        cout << value << " ";
    }
    cout << endl;
        */
        // Step 5: Gather the partial results of P*r_k from all processes
    vector<double> P_rk(n, 0.0);

        // Step 5: Gather the partial results of P*r_k from all processes
        vector<int> send_counts(size);
        vector<int> displ_rk(size);
        // Calculate send_counts and displ_rk for gathering the results of P*r_k
        for (int i = 0; i < size; ++i) {
                send_counts[i] = rows_per_process + (i < extra_rows ? 1 : 0);  // Size of local_P_rk for each process
                displ_rk[i] = (i == 0) ? 0 : displ_rk[i - 1] + send_counts[i - 1];  // Displacement for gathering
        }

        // Gather the local_P_rk vectors from all processes into P_rk at the root process
        MPI_Gatherv(local_P_rk.data(), local_P_rk.size(), MPI_DOUBLE,
                                P_rk.data(), send_counts.data(), displ_rk.data(), MPI_DOUBLE,
                                0, MPI_COMM_WORLD);

        // Step 6: Calculate the dot product of r and P_rk on the root process
    double numerator = 0.0;
    for (int i = 0; i < n; ++i) {
        cout << r[i] << " ";
    }
    cout << "\n";
    if (rank == 0) {
            numerator = dot_product(r, P_rk, n);
    }

    // Broadcast the numerator to all processes
    MPI_Bcast(&numerator, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Step 7: Calculate the norm-2 of r on the root process
    double denominator = norm_2(r, n);

    // Step 8: Calculate the Rayleigh quotient (eigenvalue approximation)
    double eigenvalue = numerator / (denominator * denominator);

        if (rank == 0) {
        end_time = MPI_Wtime();
        time_taken = end_time - start_time;
        // Print the calculated eigenvalue
        cout << "Approximated eigenvalue (Rayleigh quotient): " << eigenvalue << endl;
        cout << "Total computation time with " << size << " processors: " << time_taken << " seconds." << endl;
        }

        MPI_Finalize();
    return 0;
}