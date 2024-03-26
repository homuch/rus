#include <iostream>
#include <vector>
#include <complex>
#include <gmpxx.h>
#include "rus/pslq.h"
#include "rus/normalization.h"
#include "es/exactdecomposer.h"

using std::complex;
using std::vector;
vector<complex<mpf_class>> create_vector(double theta)
{
    using std::sin;
    vector<complex<mpf_class>> v;
    for (size_t i = 0; i < 4; i++)
    {
        v.push_back(mpf_class(sin(theta / 2 + i * M_PI / 4)));
    }
    return v;
}

Normalization::OmegaRing create_omega_ring(const vector<complex<mpf_class>> &v)
{
    vector<mpz_class> vz;
    for (auto &z : v)
    {
        vz.push_back(mpz_class(z.real()));
    }
    return Normalization::OmegaRing(vz[0], vz[1], vz[2], vz[3]);
}

void print_OmegaRing(const Normalization::OmegaRing &z)
{
    std::cout << "OmegaRing: " << z[0] << "+" << z[1] << "*w +" << z[2] << "*w^2 +" << z[3] << "*w^3" << std::endl;
}

int main()
{
    std::cout << "Stage 1: PSLQ\n";
    const double THETA = M_PI / 64;
    // Define the input vector x
    vector<complex<mpf_class>> x = create_vector(THETA);

    // ! custom error function often leads to a real epsilon which is too small
    // ! this may cause a longer circuit/T depth
    auto get_error = [THETA](const complex<mpf_class> &min_y_val, const vector<complex<mpf_class>> &b_col) -> mpf_class
    {
        Normalization::OmegaRing z = create_omega_ring(b_col);
        auto cz = z.toComplex(0);
        return abs(2 * ((decltype(cz)(cos(THETA / 2), sin(THETA / 2))) * cz).imag() / abs(cz));
    };

    PslqComplex pslq(x, 1e-5, 10000, 1, get_error);
    pslq.run();
    PslqComplex::ComplexVector r = pslq.get_result();

    // std::cout << "x . r = \n";
    // for (size_t i = 0; i < x.size(); i++)
    // {
    //     std::cout << x[i] << " * " << r[i];
    //     std::cout << ((i == x.size() - 1) ? " = \n" : " + ");
    // }
    // complex<mpf_class> result = mpf_class(0.0);
    // for (size_t i = 0; i < x.size(); i++)
    // {
    //     result += x[i] * r[i];
    //     std::cout << x[i] * r[i] << " ";
    //     std::cout << ((i == x.size() - 1) ? " = \n" : " + ");
    // }
    // std::cout << std::endl;
    // std::cout << result << std::endl;

    Normalization::OmegaRing z = create_omega_ring(r);
    print_OmegaRing(z);

    // todo: check the real epsilon
    auto cz = z.toComplex(0);
    std::cout << "real epsilon = " << 2 * ((decltype(cz)(cos(THETA / 2), sin(THETA / 2))) * cz).imag() / abs(cz) << std::endl;

    std::cout << "Stage 2: Normalization\n";

    Normalization normSolver(z);
    normSolver.solve();
    matrix2x2<mpz_class> result2 = normSolver.get_result();

    std::cout << "Stage 3: decomposition\n";
    circuit cir = exactDecomposer::decompose(result2);
    std::ostringstream oss;
    cir.toStream(oss);
    std::string str = oss.str();
    int count = std::count(str.begin(), str.end(), 'T');
    std::cout << "Number of 'T' in the stream: " << count << std::endl;
    std::cout << "Circuit depth: " << std::count(str.begin(), str.end(), '\n') << std::endl;
    return 0;
}