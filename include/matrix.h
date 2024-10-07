class Matrix {
private:
    int rows, cols;
    double* data;

public:
    Matrix(unsigned int& r, unsigned int& c);
    void setValue(unsigned int i, unsigned int j, double value);
    double getValue(unsigned int i, unsigned int j);
    int getSize();
    void print();
    ~Matrix();
};
