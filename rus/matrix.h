#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <vector>
#include <iostream>
#include <cassert>

template <typename num_type>
class MyMatrix
{
public:
    typedef typename std::vector<num_type>::size_type size_t;
    MyMatrix(size_t rows, size_t cols, num_type val);
    MyMatrix(size_t rows, size_t cols);
    inline num_type &operator()(size_t i, size_t j)
    {
        assert(i < rows && j < cols);
        return data[i * cols + j];
    }
    inline const num_type &operator()(size_t i, size_t j) const
    {
        assert(i < rows && j < cols);
        return data[i * cols + j];
    }
    inline size_t get_rows() const { return rows; }
    inline size_t get_cols() const { return cols; }
    void to_identity(const num_type &one = 1, const num_type &zero = 0);
    void set_all(num_type val);
    void print() const;
    std::vector<num_type> diagonal() const;
    std::vector<num_type> row(size_t i) const;

    std::vector<num_type> col(size_t j) const;

    void swap_row(size_t i, size_t j);

    void swap_col(size_t i, size_t j);

private:
    size_t rows, cols;
    std::vector<num_type> data;
};

#endif // __MATRIX_H__