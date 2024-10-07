#include "../include/vectors.h"
#include <iostream>
#include <cmath>
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

// Getters
double Vectors::getDataAtIndex(unsigned int i) {
    return data[i];
}

int Vectors::getSize() {
    return size;
}

// Methods
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
