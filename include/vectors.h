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
    double norm();
    void printVec();
};
