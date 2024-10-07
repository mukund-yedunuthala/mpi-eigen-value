class Vectors
{
private:
    unsigned int size;
    double* data = nullptr;
public:
    Vectors(unsigned int& s);
    ~Vectors();
    Vectors(const Vectors& other);
    Vectors& operator=(const Vectors& other);
    double getDataAtIndex(unsigned int i);
    int getSize();

    void setData(double*& inputData);
    void setValueAtIndex(unsigned int i, double value);
    void setToValue(double value);

    double sum();
    double norm();
    double dot(Vectors& other);
    void print();
};
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
