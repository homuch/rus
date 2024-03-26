#include <complex>
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
#include <functional>
#include <gmpxx.h>
#include "pslq.h"
#include "matrix.h"
#include "numeric_definition.h"

mpf_class round(const mpf_class &x)
{
    return floor(x + 0.5);
}

PslqComplex::PslqComplex(
    const ComplexVector &input_vector,
    const real_t &epsilon,
    const unsigned &max_iterations,
    int debug_level,
    std::function<real_t(const Complex &min_y_val, const ComplexVector &b_col)> calculate_error)
    : n(input_vector.size()), eps(epsilon), itm(max_iterations), idb(debug_level),
      gamma(sqrt(real_t(4.0) / real_t(3.0))), status(STATUS::ITERATION),
      x(input_vector), y(n), b(n, n), h(n, n - 1), get_error(calculate_error)
{
}

void PslqComplex::run()
{
    unsigned it_count = itm;
    unsigned db_period = 500;

    if (idb >= 2)
    {
        std::cout << "input x: \n";
        print_vector(x);
    }

    pslq_init();

    while (status == STATUS::ITERATION && it_count-- > 0)
    {
        if (it_count % db_period == 0 && idb > 1)
            std::cout << it_count << std::endl;

        pslq_iteration();
    }
    size_t min_y_id = min(y).id;
    r = b.col(min_y_id);
    switch (status)
    {
    case STATUS::SUCCESS:
    {
        std::cout << "Success" << std::endl;
        break;
    }
    case STATUS::ITERATION:
        std::cout << "Iteration limit exceeded " << itm << std::endl;
        break;

    default:
        std::cout << "Failed" << std::endl;
        break;
    }
}
PslqComplex::ComplexVector PslqComplex::get_result() const
{
    return r;
}

void PslqComplex::pslq_iteration()
{
    using std::abs;
    using std::swap;

    const size_t im = max_for(h.diagonal(), [this](const Complex &c, const size_t &si) -> real_t
                              { return pow(this->gamma.get_d(), si) * abs(c); })
                          .id; // assume si won't be too large

    const size_t im1 = im + 1;
    swap(y[im], y[im1]);

    b.swap_col(im, im1);
    h.swap_row(im, im1);

    if (im < n - 2)
        update_h(im);

    iter_reduction(im);

    Id_Val min_y = min(y);

    if (get_error(min_y.val, b.col(min_y.id)) < eps)
        status = STATUS::SUCCESS;
}
void PslqComplex::pslq_init()
{
    using std::abs;

    b.to_identity(Complex(real_t(1)), Complex(real_t(0)));
    ComplexVector s = partial_sum(x);
    const Complex sum = s[0];
    for (auto &ele : s)
        ele /= sum;

    y = x;
    for (auto &ele : y)
        ele /= sum;

    init_h(s);
    init_reduction();
    if (abs(min(y).val) < eps)
    {
        status = STATUS::SUCCESS;
    }
    if (idb >= 3)
    {
        std::cout << __func__ << ": y: \n";
        print_vector(y);
        std::cout << __func__ << ": s: \n";
        print_vector(s);
        std::cout << __func__ << ": b: \n";
        b.print();
        std::cout << __func__ << ": h: \n";
        h.print();
    }
}

/**
 * @brief Compute the partial sum of the input vector
 * @param v input vector
 * @return partial sum of the input: sqrt(sum(v[i:n] . v[i:n))
 */
PslqComplex::ComplexVector PslqComplex::partial_sum(const ComplexVector &v) const
{
    using std::sqrt;

    ComplexVector result(v.size());
    Complex sum = real_t(0.0);
    for (int i = v.size() - 1; i >= 0; --i)
    {
        sum += v[i] * v[i];
        result[i] = sqrt(sum);
    }
    return result;
}
void PslqComplex::init_h(const ComplexVector &s)
{
    for (int j = 0; j < n - 1; ++j)
    {
        for (int i = 0; i <= j - 1; ++i)
            h(i, j) = real_t(0.0);

        h(j, j) = s[j + 1] / s[j];
        const Complex temp = y[j] / (s[j] * s[j + 1]);

        for (int i = j + 1; i < n; ++i)
            h(i, j) = -y[i] * temp;
    }
}
void PslqComplex::init_reduction()
{
    using std::abs;
    using std::round;
    for (int i = 1; i < n; ++i)
    {
        for (int j = i - 1; j >= 0; --j)
        {
            Complex ratio = h(i, j) / h(j, j);
            Complex temp = Complex(round(ratio.real()), round(ratio.imag()));
            if (temp != Complex(0.0))
            {
                y[j] = y[j] + temp * y[i];

                for (int k = i; k < n; ++k)
                {
                    b(k, j) += temp * b(k, i);
                }

                for (int k = 0; k <= j; ++k)
                {
                    h(i, k) -= temp * h(j, k);
                }
            }
        }
    }
}

void PslqComplex::iter_reduction(const size_t &im)
{
    using std::min;
    using std::round;

    const size_t im1 = im + 1;
    for (size_t i = im1; i < n; i++)
    {
        size_t j1 = min(i - 1, im1);

        for (int j = j1; j >= 0; j--)
        {
            Complex ratio = h(i, j) / h(j, j);
            Complex temp = Complex(round(ratio.real()), round(ratio.imag()));
            if (temp != Complex(0.0))
            {
                y[j] = y[j] + temp * y[i];

                for (int k = 0; k < n; k++)
                {
                    b(k, j) += temp * b(k, i);
                }

                for (int k = 0; k <= j; k++)
                {
                    h(i, k) -= temp * h(j, k);
                }
            }
        }
    }
}

void PslqComplex::update_h(const size_t &im)
{
    using std::sqrt;

    const size_t im1 = im + 1;
    Complex t1, t2, t3, t4;
    if (im <= n - 2)
    {
        t1 = h(im, im);
        t2 = h(im, im1);
        t3 = sqrt(t1 * t1 + t2 * t2);
        t1 = t1 / t3;
        t2 = t2 / t3;

        for (size_t i = im; i < n; i++)
        {
            t3 = h(i, im);
            t4 = h(i, im1);
            h(i, im) = t1 * t3 + t2 * t4;
            h(i, im1) = -t2 * t3 + t1 * t4;
        }
    }
}

PslqComplex::Id_Val PslqComplex::min(const ComplexVector &v)
{
    using std::abs;
    Id_Val result;
    result.val = v[0];
    result.id = 0;
    real_t minimum = abs(result.val);

    for (size_t i = 1; i < v.size(); i++)
    {

        if (abs(v[i]) < minimum)
        {
            result.val = v[i];
            result.id = i;
            minimum = abs(v[i]);
        }
    }
    return result;
}
PslqComplex::Id_Val PslqComplex::max(const ComplexVector &v)
{
    using std::abs;
    Id_Val result;
    result.val = v[0];
    result.id = 0;
    real_t maximum = abs(result.val);

    for (size_t i = 1; i < v.size(); i++)
    {
        if (maximum < abs(v[i]))
        {
            result.val = v[i];
            result.id = i;
            maximum = abs(v[i]);
        }
    }
    return result;
}

PslqComplex::Id_Val PslqComplex::min_for(const ComplexVector &v, std::function<real_t(const Complex &, const size_t &)> trans)
{
    using std::abs;
    Id_Val result;
    result.val = trans(v[0], 0);
    result.id = 0;
    real_t minimum = abs(result.val);

    for (size_t i = 1; i < v.size(); i++)
    {
        if (trans(v[i], i) < minimum)
        {
            result.val = v[i];
            result.id = i;
            minimum = abs(v[i]);
        }
    }
    return result;
}
PslqComplex::Id_Val PslqComplex::max_for(const ComplexVector &v, std::function<real_t(const Complex &, const size_t &)> trans)
{
    using std::abs;
    Id_Val result;
    result.val = trans(v[0], 0);
    result.id = 0;
    real_t maximum = abs(result.val);

    for (size_t i = 1; i < v.size(); i++)
    {
        if (maximum < trans(v[i], i))
        {
            result.val = v[i];
            result.id = i;
            maximum = abs(v[i]);
        }
    }
    return result;
}
void PslqComplex::print_vector(const ComplexVector &v) const
{
    for (size_t i = 0; i < v.size(); i++)
    {
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
}
