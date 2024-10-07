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
