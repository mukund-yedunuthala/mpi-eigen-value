#include "../include/vectors.h"
#include <iostream>
#include <cmath>
Vectors::Vectors(unsigned int& s) {
    size = s;
    data = new double[s];
    for (unsigned int i = 0; i < size; i++) {
        data[i] = 0.0;
    }
}

Vectors::~Vectors() {
    delete [] data;
}

void Vectors::setData(double*& inputData) {
    for (unsigned int i = 0; i < size; i++) {
        data[i] = inputData[i];
    }
}

double Vectors::getDataAtIndex(unsigned int i) {
    return data[i];
}

double Vectors::norm() {
    double norm = 0.0;
    for (unsigned int i = 0; i < size; i++) {
        norm = norm + (data[i] * data[i]) ;
    }
    return std::sqrt(norm);
}

void Vectors::printVec() {
    for (unsigned int i = 0; i < size; i++) {
        std::cout << data[i] << " ";
    }
    std::cout << "\n";
}
