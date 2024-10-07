class Vectors
{
private:
    unsigned int size;
    double* data = nullptr;
public:
    Vectors(unsigned int& s);
    ~Vectors();
    void setData(double*& inputData);
    double getDataAtIndex(unsigned int i);
    void setValueAtIndex(unsigned int i, double value);
    double norm();
    void printVec();
};
