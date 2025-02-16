/**
 * \file 
 * \brief Source file consisting a collection of basic linear algebra classes. 
 * \author Mukund Yedunuthala
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "../include/linalg.h"

// C'tor
/**
 * \brief Constructs a Vectors object with the given size.
 * \param s Reference to the size of the vector.
 */
Vectors::Vectors(unsigned int& s) {
    size = s;
    data = new double[s]{};
}

/**
 * \brief Destructor to free allocated memory.
 */
Vectors::~Vectors() {
    delete [] data;
}

// Copy
/**
 * \brief Copy constructor for deep copying another Vectors object.
 * \param other The vector to be copied.
 */
Vectors::Vectors(const Vectors& other) {
    size = other.size;
    for (unsigned int i = 0; i < size; i++) {
        data[i] = other.data[i];
    }
}

// Copy assignment
/**
 * \brief Copy assignment operator for deep copying another Vectors object.
 * \param other The vector to be assigned.
 * \return Reference to the assigned vector.
 */
Vectors& Vectors::operator=(const Vectors& other){
    if (this == &other) {
        return  *this;
    }
    delete [] data;
    size = other.size;
    data = new double[size];
    for (unsigned int i = 0; i < size; i++) {
        data[i] = other.data[i];
    }
    return *this;
}

// Setters
/**
 * \brief Sets the vector data using an external array.
 * \param inputData Reference to the pointer of the input data.
 */
void Vectors::setData(double*& inputData) {
    for (unsigned int i = 0; i < size; i++) {
        data[i] = inputData[i];
    }
}

/**
 * \brief Sets the value at a specific index.
 * \param i The index where the value should be set.
 * \param value The value to set.
 */
void Vectors::setValueAtIndex(unsigned int i, double value) {
    data[i] = value;
}

/**
 * \brief Sets all elements of the vector to a given value.
 * \param value The value to set all elements to.
 */
void Vectors::setToValue(double value) {
    for (unsigned int i = 0; i < size; i++) {
        data[i] = value;
    }
}

// Getters
/**
 * \brief Retrieves the value at a specific index in the vector.
 * \param i The index to access.
 * \return The value at the specified index.
 */
double Vectors::getDataAtIndex(unsigned int i) {
    return data[i];
}

/**
 * \brief Gets the size of the vector.
 * \return The size of the vector.
 */
int Vectors::getSize() {
    return size;
}

/**
 * \brief Retrieves the pointer to the data array.
 * \return Pointer to the vector data.
 */
double* Vectors::getData() {
    return data;
}

// Methods
/**
 * \brief Computes the sum of all elements in the vector.
 * \return The sum of vector elements.
 */
double Vectors::sum() {
    double sum = 0.0;
    for (unsigned int i = 0; i < size; i++) {
        sum = sum + data[i] ;
    }
    return sum;
}

/**
 * \brief Computes the absolute sum of all elements.
 * \return The sum of absolute values of all elements.
 */
double Vectors::abssum() {
    double sum = 0.0;
    for (unsigned int i = 0; i < size; i++) {
        sum = sum + std::abs(data[i]);
    }
    return sum;
}

/**
 * \brief Computes the Euclidean norm (magnitude) of the vector.
 * \return The norm of the vector.
 */
double Vectors::norm() {
    double norm = 0.0;
    for (unsigned int i = 0; i < size; i++) {
        norm = norm + (data[i] * data[i]) ;
    }
    return std::sqrt(norm);
}

/**
 * \brief Computes the dot product of this vector with another vector.
 * \param other The other vector to compute the dot product with.
 * \return The dot product result.
 */
double Vectors::dot(Vectors& other) {
    double value = 0.0;
    if (other.getSize() == size) {
        for (unsigned int i = 0; i < size; i++) {
            value = value + (data[i]*other.getDataAtIndex(i));
        }
    }
    return  value;
}

/**
 * \brief Prints the vector elements to the console.
 */
void Vectors::print() {
    for (unsigned int i = 0; i < size; i++) {
        std::cout << data[i] << "\t";
    }
    std::cout << "\n";
}

/**
 * \brief Constructs a Matrix object with the given dimensions.
 * \param r Reference to the number of rows.
 * \param c Reference to the number of columns.
 */
Matrix::Matrix(unsigned int& r, unsigned int& c) {
    rows = r;
    cols = c;
    data = new double [rows*cols]{};
}

/**
 * \brief Copy constructor for deep copying another Matrix object.
 * \param other The matrix to be copied.
 */
Matrix::Matrix(const Matrix& other) {
    rows = other.rows;
    cols = other.cols;
    data = new double [rows*cols]{};
    for (unsigned int i = 0; i<rows; i++) {
        for (unsigned int j = 0; j<cols; j++) {
            data[i*cols+j] = other.data[i*cols+j];
        }
    }
}

/**
 * \brief Copy assignment operator for deep copying another Matrix object.
 * \param other The matrix to be assigned.
 * \return Reference to the assigned matrix.
 */
Matrix& Matrix::operator=(const Matrix& other) {
    if (this == &other) {
        return *this;
    }
    delete [] data;
    rows = other.rows;
    cols = other.cols;
    data = new double [rows*cols]{};
    for (unsigned int i = 0; i<rows; i++) {
        for (unsigned int j = 0; j<cols; j++) {
            data[i*cols+j] = other.data[i*cols+j];
        }
    }
    return *this;
}

/**
 * \brief Destructor to free allocated memory.
 */
Matrix::~Matrix() {
    delete [] data;
}

// Setters
/**
 * \brief Sets the value at a specific row and column.
 * \param i The row index.
 * \param j The column index.
 * \param value The value to be set.
 */
void Matrix::setValue(unsigned int i, unsigned int j, double value) {
    data[i*cols+j] = value;
}

// Getters
/**
 * \brief Retrieves the value at a specific row and column.
 * \param i The row index.
 * \param j The column index.
 * \return The value at the specified position.
 */
double Matrix::getValue(unsigned int i, unsigned int j) {
    return data[(i*cols)+j];
}

/**
 * \brief Gets the total number of elements in the matrix.
 * \return The total number of elements (rows * cols).
 */
int Matrix::getSize() {
    return rows*cols;
}

/**
 * \brief Prints the matrix elements to the console.
 */
void Matrix::print() {
    for (unsigned int i = 0; i<rows; i++) {
        for (unsigned int j = 0; j<cols; j++) {
            std::cout << std::setprecision(3) << data[i*cols+j] << "\t";
        }
        std::cout << "\n";
    }
}

// Methods
/**
 * \brief Computes the matrix-vector multiplication.
 * \param vec The vector to multiply with the matrix.
 * \return A new Vectors object containing the result of the multiplication.
 */
Vectors Matrix::dot(Vectors& vec) {
    unsigned int dim = vec.getSize();
    Vectors result(dim);
    if (rows == dim) {
        for (unsigned int i = 0; i < rows; i++) {
            double sum{};
            for (unsigned int j = 0; j < rows; j++) {
                sum = sum + (vec.getDataAtIndex(j) * data[i*cols+j]);
                result.setValueAtIndex(i, sum);
            }
        }
    } else {

    }
    return result;
}
