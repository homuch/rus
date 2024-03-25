#include <vector>
#include <iostream>
#include <cassert>
#include <complex>
#include "matrix.h"
#include "numeric_definition.h"

template <typename num_type>
MyMatrix<num_type>::MyMatrix(size_t rows, size_t cols, num_type val) : rows(rows), cols(cols), data(rows * cols, val) {}

template <typename num_type>
MyMatrix<num_type>::MyMatrix(size_t rows, size_t cols) : rows(rows), cols(cols), data(rows * cols) {}

template <typename num_type>
void MyMatrix<num_type>::to_identity(const num_type &one, const num_type &zero)
{
    if (rows != cols)
    {
        std::cerr << "Matrix is not square" << std::endl;
        return;
    }
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            (*this)(i, j) = (i == j) ? one : zero;
        }
    }
}

template <typename num_type>
void MyMatrix<num_type>::set_all(num_type val)
{
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            (*this)(i, j) = val;
        }
    }
}

template <typename num_type>
void MyMatrix<num_type>::print() const
{
    for (size_t i = 0; i < rows; i++)
    {
        std::cout << "| ";
        for (size_t j = 0; j < cols; j++)
        {
            std::cout << (*this)(i, j) << " ";
        }
        std::cout << " |\n";
    }
    std::cout << std::flush;
}

template <typename num_type>
std::vector<num_type> MyMatrix<num_type>::diagonal() const
{
    using std::min;
    const size_t length = min(rows, cols);

    std::vector<num_type> result;
    for (size_t i = 0; i < length; i++)
    {
        result.push_back((*this)(i, i));
    }
    return result;
}

template <typename num_type>
std::vector<num_type> MyMatrix<num_type>::row(size_t i) const
{
    std::vector<num_type> result;
    for (size_t j = 0; j < cols; j++)
    {
        result.push_back((*this)(i, j));
    }
    return result;
}

template <typename num_type>
std::vector<num_type> MyMatrix<num_type>::col(size_t j) const
{
    std::vector<num_type> result;
    for (size_t i = 0; i < rows; i++)
    {
        result.push_back((*this)(i, j));
    }
    return result;
}

template <typename num_type>
void MyMatrix<num_type>::swap_row(size_t i, size_t j)
{
    using std::swap;
    if (i == j)
        return;
    for (size_t k = 0; k < cols; k++)
    {
        swap((*this)(i, k), (*this)(j, k));
    }
}

template <typename num_type>
void MyMatrix<num_type>::swap_col(size_t i, size_t j)
{
    using std::swap;
    if (i == j)
        return;
    for (size_t k = 0; k < rows; k++)
    {
        swap((*this)(k, i), (*this)(k, j));
    }
}

template class MyMatrix<std::complex<real_t>>;