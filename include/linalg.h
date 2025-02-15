/**
 * \file 
 * \brief Header file consisting a collection of basic linear algebra classes. 
 * \author Mukund Yedunuthala
 */

/**
 * \class Vectors
 * \brief A class representing a dynamic array of double values with vector operations.
 * \author Mukund Yedunuthala
 */
class Vectors
{
private:
    unsigned int size;  ///< The size of the vector.
    double* data = nullptr;  ///< Pointer to the dynamically allocated array of vector elements.

public:
    /**
     * \brief Constructs a Vectors object with the given size.
     * \param s Reference to the size of the vector.
     */
    Vectors(unsigned int& s);

    /**
     * \brief Destructor to free allocated memory.
     */
    ~Vectors();

    /**
     * \brief Copy constructor for deep copying another Vectors object.
     * \param other The vector to be copied.
     */
    Vectors(const Vectors& other);

    /**
     * \brief Copy assignment operator for deep copying another Vectors object.
     * \param other The vector to be assigned.
     * \return Reference to the assigned vector.
     */
    Vectors& operator=(const Vectors& other);

    /**
     * \brief Retrieves the value at a specific index in the vector.
     * \param i The index to access.
     * \return The value at the specified index.
     */
    double getDataAtIndex(unsigned int i);

    /**
     * \brief Gets the size of the vector.
     * \return The size of the vector.
     */
    int getSize();

    /**
     * \brief Retrieves the pointer to the data array.
     * \return Pointer to the vector data.
     */
    double* getData();

    /**
     * \brief Sets the vector data using an external array.
     * \param inputData Reference to the pointer of the input data.
     */
    void setData(double*& inputData);

    /**
     * \brief Sets the value at a specific index.
     * \param i The index where the value should be set.
     * \param value The value to set.
     */
    void setValueAtIndex(unsigned int i, double value);

    /**
     * \brief Sets all elements of the vector to a given value.
     * \param value The value to set all elements to.
     */
    void setToValue(double value);

    /**
     * \brief Computes the sum of all elements in the vector.
     * \return The sum of vector elements.
     */
    double sum();

    /**
     * \brief Computes the absolute sum of all elements.
     * \return The sum of absolute values of all elements.
     */
    double abssum();

    /**
     * \brief Computes the Euclidean norm (magnitude) of the vector.
     * \return The norm of the vector.
     */
    double norm();

    /**
     * \brief Computes the dot product of this vector with another vector.
     * \param other The other vector to compute the dot product with.
     * \return The dot product result.
     */
    double dot(Vectors& other);

    /**
     * \brief Prints the vector elements to the console.
     */
    void print();
};

/**
 * \class Matrix
 * \brief A class representing a 2D matrix with basic operations.
 */
class Matrix
{
private:
    int rows, cols;  ///< Number of rows and columns in the matrix.
    double* data;  ///< Pointer to the dynamically allocated array of matrix elements.

public:
    /**
     * \brief Constructs a Matrix object with the given dimensions.
     * \param r Reference to the number of rows.
     * \param c Reference to the number of columns.
     */
    Matrix(unsigned int& r, unsigned int& c);

    /**
     * \brief Copy constructor for deep copying another Matrix object.
     * \param other The matrix to be copied.
     */
    Matrix(const Matrix& other);

    /**
     * \brief Copy assignment operator for deep copying another Matrix object.
     * \param other The matrix to be assigned.
     * \return Reference to the assigned matrix.
     */
    Matrix& operator=(const Matrix& other);

    /**
     * \brief Destructor to free allocated memory.
     */
    ~Matrix();

    /**
     * \brief Retrieves the value at a specific row and column.
     * \param i The row index.
     * \param j The column index.
     * \return The value at the specified position.
     */
    double getValue(unsigned int i, unsigned int j);

    /**
     * \brief Gets the total number of elements in the matrix.
     * \return The total number of elements (rows * cols).
     */
    int getSize();

    /**
     * \brief Sets the value at a specific row and column.
     * \param i The row index.
     * \param j The column index.
     * \param value The value to be set.
     */
    void setValue(unsigned int i, unsigned int j, double value);

    /**
     * \brief Computes the matrix-vector multiplication.
     * \param vec The vector to multiply with the matrix.
     * \return A new Vectors object containing the result of the multiplication.
     */
    Vectors dot(Vectors& vec);

    /**
     * \brief Prints the matrix elements to the console.
     */
    void print();
};
