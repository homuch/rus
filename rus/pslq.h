#ifndef __PSLQ_COMPLEX_H__
#define __PSLQ_COMPLEX_H__

#include <complex>
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
#include <functional>
#include <gmpxx.h>
#include "matrix.h"
#include "numeric_definition.h"

mpf_class round(const mpf_class &x);

// todo: correct epsilon(that reflect the difference of phase)
class PslqComplex
{
public:
    typedef std::complex<real_t> Complex;
    typedef std::vector<Complex> ComplexVector;
    typedef MyMatrix<Complex> ComplexMatrix;
    typedef ComplexVector::size_type size_t;
    enum STATUS
    {
        ITERATION,
        SUCCESS,
        FAILURE
    };
    struct Id_Val
    {
        size_t id;
        Complex val;
    };
    PslqComplex(
        const ComplexVector &input_vector,
        const real_t &epsilon,
        const unsigned &max_iterations,
        int debug_level = 0);
    void run();

    ComplexVector get_result() const;

private:
    void pslq_iteration();
    void pslq_init();

    /**
     * @brief Compute the partial sum of the input vector
     * @param v input vector
     * @return partial sum of the input: sqrt(sum(v[i:n] . v[i:n))
     */
    ComplexVector partial_sum(const ComplexVector &v) const;
    void init_h(const ComplexVector &s);
    void init_reduction();

    void iter_reduction(const size_t &im);

    void update_h(const size_t &im);

    Id_Val min(const ComplexVector &v) const;
    Id_Val max(const ComplexVector &v) const;

    Id_Val min_for(const ComplexVector &v, std::function<real_t(const Complex &, const size_t &)> trans) const;
    Id_Val max_for(const ComplexVector &v, std::function<real_t(const Complex &, const size_t &)> trans) const;
    void print_vector(const ComplexVector &v) const;
    const size_t n;
    int idb; // debug level
    int nwds;
    unsigned itm; // maximum number of iterations
    real_t eps;
    const real_t gamma;
    ComplexVector x, y, r;
    ComplexMatrix b, h;
    STATUS status;
};

#endif // __PSLQ_COMPLEX_H__