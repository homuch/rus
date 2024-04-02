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
#include <set>
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
        int debug_level = 0,
        int effort_level = 1,
        std::function<real_t(const Complex &min_y_val, const ComplexVector &b_col)> calculate_error = [](const Complex &min_y_val, const ComplexVector &b_col) -> real_t
        {
            using std::abs;
            return abs(min_y_val) * abs(max(b_col).val);
        }

    );
    void run();

    std::vector<ComplexVector> get_results() const;

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

    static Id_Val min(const ComplexVector &v);
    static Id_Val max(const ComplexVector &v);

    static Id_Val min_for(const ComplexVector &v, std::function<real_t(const Complex &, const size_t &)> trans);
    static Id_Val max_for(const ComplexVector &v, std::function<real_t(const Complex &, const size_t &)> trans);
    void print_vector(const ComplexVector &v) const;
    static bool compareComplexVectorRealPart(const PslqComplex::ComplexVector &a, const PslqComplex::ComplexVector &b);
    const size_t n;
    int idb; // debug level
    int nwds;
    int effort;   // effort level (collect more results if effort > 1)
    unsigned itm; // maximum number of iterations
    real_t eps;
    const real_t gamma;
    ComplexVector x, y;
    ComplexMatrix b, h;
    STATUS status;
    std::set<ComplexVector, decltype(compareComplexVectorRealPart) *> r_sets;
    // error function (e.g. For the phase error)
    std::function<real_t(const Complex &min_y_val, const ComplexVector &b_col)> get_error;
};

#endif // __PSLQ_COMPLEX_H__