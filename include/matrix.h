#include "../include/vectors.h"

class Matrix {
private:
    int rows, cols;
    double* data;

public:
    Matrix(unsigned int& r, unsigned int& c);
    Matrix(const Matrix& other);
    Matrix& operator=(const Matrix& other);
    ~Matrix();

    double getValue(unsigned int i, unsigned int j);
    int getSize();

    void setValue(unsigned int i, unsigned int j, double value);

    Vectors dot(Vectors& vec);
    void print();
};
