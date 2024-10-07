#include "../include/matrix.h"
#include <iostream>

Matrix::Matrix(unsigned int& r, unsigned int& c) {
    rows = r;
    cols = c;
    data = new double [rows*cols]{};
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
            std::cout << data[i*cols+j] << " ";
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
    return  result;
}
