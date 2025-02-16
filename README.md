## **Overview**  
This project implements a **parallel dense matrix and vector format** in **C++** using **MPI** and the **C++ Standard Template Library (STL)**. The implementation includes:  
- **Custom C++ matrix and vector classes**, utilizing dynamically allocated arrays for efficient data handling.  
- **Link matrix generation** using a **Bernoulli distribution** with a user-defined **random seed**.  
- **Parallel power iteration** to approximate the **largest eigenvalue** and **left stochastic matrix** \( P \) derived from the generated **link matrix** \( L \).  

## **Features**  
- **C++ Matrix and Vector Classes**:  
  - Custom implementations with dynamically allocated arrays.  
  - Efficient memory management and operator overloading.  
- **Link Matrix Generation**:  
  - Constructed using a **Bernoulli distribution** (0/1 values).  
  - Accepts a **user-defined random seed** for reproducibility.  
  - Transformed into a **left stochastic matrix** \( P \) by normalizing row sums.  
- **Parallel Matrix-Vector Computation**:  
  - Utilizes **MPI** for distributed Rayleigh method iterations.  
  - Supports **scatter/gather** for efficient communication.  
- **Power Iteration Algorithm**:  
  - Iteratively approximates the **dominant eigenvalue**.
  - Implements convergence checks for numerical stability.  

## **Dependencies**  
- **MPI Library** (e.g., OpenMPI or MPICH)  
- **CMake**
- **C++ Standard Template Library (STL)**  

## **Installation**  
1. Install an MPI implementation and CMake:  
   ```sh
   sudo apt install mpich  # Debian/Ubuntu  
   sudo dnf install openmpi  # Fedora 
   sudo pacman -S cmake # Arch 
   ```  
2. Clone the repository:  
   ```sh
   git clone <repository_url>
   cd <repository>
   ```  
3. Compile the program using CMake:
    
    - For just a sequential implementation (which does not depend on MPI): 
        ```sh
        cmake -S . -B build -DWITH_MPI=0 
        ```  

    - For a parallel implementation (which depends on MPI): 
        ```sh
        cmake -S . -B build -DWITH_MPI=1 
        ```  

    - For just a sequential implementation and to build tests: 
        ```sh
        cmake -S . -B build -DWITH_MPI=0 -DBUILD_TESTING=1
        ```  
## **Usage**  
Run the program with the desired matrix size:  
```sh
mpirun -np <num_processes> ./EVP <matrix_size>
```
where:  
- `<num_processes>` is the number of MPI processes.  
- `<matrix_size>` defines the dimension \( N \times N \) of the link matrix.  

## **Example**  
```sh
mpirun -np 4 ./EVP 100
```  

## **Implementation Details**  

### **Matrix and Vector Classes**  
- **Matrix Class (`Matrix`)**:  
  - Uses a **dynamically allocated 2D array** (or 1D row-major storage).  
  - Implements matrix operations (e.g., matrix-vector multiplication).  
  - Includes **MPI distribution functions** for parallel computation.  
- **Vector Class (`Vector`)**:  
  - Uses a **dynamically allocated array**.  
  - Implements vector operations (e.g., dot product, normalization).  
  - Overloads **operators** for clean syntax.  

### **Parallelization Strategy**  
- **Row-wise distribution** of the matrix among MPI processes.  
- Uses **MPI Scatter/Gather** for communication efficiency.  

## **Performance Considerations**  
- **Minimized MPI communication overhead** via optimized data distribution.  
- **Memory-efficient storage** using dynamically allocated arrays.  
- **Scalability**: Supports large matrices distributed across multiple processes.  

## **References**  
- **MPI Documentation**: [https://www.mpi-forum.org/](https://www.mpi-forum.org/)  
- **Power Iteration Method**: [https://en.wikipedia.org/wiki/Power_iteration](https://en.wikipedia.org/wiki/Power_iteration)  
