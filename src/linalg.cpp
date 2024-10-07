#include <iostream>
#include <iomanip>
#include <cmath>

#include "../include/linalg.h"

// C'tor
Vectors::Vectors(unsigned int& s) {
    size = s;
    data = new double[s]{};
}

// D'tor
Vectors::~Vectors() {
    delete [] data;
}

// Copy
Vectors::Vectors(const Vectors& other) {
    size = other.size;
    for (unsigned int i = 0; i < size; i++) {
        data[i] = other.data[i];
    }
}

// Copy assignment
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
void Vectors::setData(double*& inputData) {
    for (unsigned int i = 0; i < size; i++) {
        data[i] = inputData[i];
    }
}

void Vectors::setValueAtIndex(unsigned int i, double value) {
    data[i] = value;
}

void Vectors::setToValue(double value) {
    for (unsigned int i = 0; i < size; i++) {
        data[i] = value;
    }
}
// Getters
double Vectors::getDataAtIndex(unsigned int i) {
    return data[i];
}

int Vectors::getSize() {
    return size;
}

// Methods
double Vectors::sum() {
    double sum = 0.0;
    for (unsigned int i = 0; i < size; i++) {
        sum = sum + data[i] ;
    }
    return sum;
}
double Vectors::norm() {
    double norm = 0.0;
    for (unsigned int i = 0; i < size; i++) {
        norm = norm + (data[i] * data[i]) ;
    }
    return std::sqrt(norm);
}

double Vectors::dot(Vectors& other) {
    double value = 0.0;
    if (other.getSize() == size) {
        for (unsigned int i = 0; i < size; i++) {
            value = value + (data[i]*other.getDataAtIndex(i));
        }
    }
    return  value;
}

void Vectors::print() {
    for (unsigned int i = 0; i < size; i++) {
        std::cout << data[i] << " ";
    }
    std::cout << "\n";
}

Matrix::Matrix(unsigned int& r, unsigned int& c) {
    rows = r;
    cols = c;
    data = new double [rows*cols]{};
}

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

Matrix::~Matrix() {
    delete [] data;
}

// Setters
void Matrix::setValue(unsigned int i, unsigned int j, double value) {
    data[i*cols+j] = value;
}
// Getters
double Matrix::getValue(unsigned int i, unsigned int j) {
    return data[(i*cols)+j];
}

int Matrix::getSize() {
    return rows*cols;
}

void Matrix::print() {
    for (unsigned int i = 0; i<rows; i++) {
        for (unsigned int j = 0; j<cols; j++) {
            std::cout << std::setprecision(3) << data[i*cols+j] << "\t";
        }
        std::cout << "\n";
    }
}

// Methods
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
