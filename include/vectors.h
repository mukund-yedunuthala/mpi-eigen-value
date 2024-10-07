class Vectors
{
private:
    unsigned int size;
    double* data = nullptr;
public:
    Vectors(unsigned int& s);
    ~Vectors();

    double getDataAtIndex(unsigned int i);
    int getSize();

    void setData(double*& inputData);
    void setValueAtIndex(unsigned int i, double value);

    double norm();
    double dot(Vectors& other);
    void print();
};
