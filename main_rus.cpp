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
    // Define the input vector x
    vector<complex<mpf_class>> x = create_vector(M_PI / 8);

    PslqComplex pslq(x, 1e-6, 10000, 1);
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

    std::cout << "Stage 2: Normalization\n";

    Normalization normSolver(z);
    normSolver.solve();
    matrix2x2<mpz_class> result2 = normSolver.get_result();

    std::cout << "Stage 3: decomposition\n";
    circuit cir = exactDecomposer::decompose(result2);
    cir.toStream(std::cout);
    return 0;
}