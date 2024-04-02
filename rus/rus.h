#ifndef __RUS__H__
#define __RUS__H__
#include <iostream>
#include <vector>
#include <complex>
#include <gmpxx.h>
#include <mpfr.h>
#include <functional>
#include <algorithm>
#include <string>
#include <memory>
#include <stdexcept>
#include "pslq.h"
#include "normalization.h"
#include "es/exactdecomposer.h"
#include "gatelibrary.h"
class RUS
{
public:
    // todo: add depth as criterion
    enum CRITERION
    {
        T_COUNT, // T count
        G_COUNT, // gate count
    };
    RUS(int debug_level, const mpf_class &epsilon, int effort_level, CRITERION criterion);

    RUS(int debug_level, const mpf_class &epsilon, int effort_level, int pslq_iters, CRITERION criterion);

    void run(const mpf_class &theta, circuit &best_cir);

    static mpf_class parse_theta(const std::string &theta_str);
    static mpf_class sin(const mpf_class &theta);
    static mpf_class cos(const mpf_class &theta);

private:
    static std::vector<std::complex<mpf_class>> create_vector(const double &theta);
    static std::vector<std::complex<mpf_class>> create_vector(const mpf_class &theta);
    static Normalization::OmegaRing create_omega_ring(const std::vector<std::complex<mpf_class>> &v);
    static inline void print_OmegaRing(const Normalization::OmegaRing &z)
    {
        std::cout << "OmegaRing: " << z[0] << "+" << z[1] << "*w +" << z[2] << "*w^2 +" << z[3] << "*w^3" << std::endl;
    }
    static std::function<real_t(const PslqComplex::Complex &min_y_val, const PslqComplex::ComplexVector &b_col)>
    gen_error_function(const mpf_class &theta);
    static int count_t(const circuit &cir);
    static int count_gate(const circuit &cir);
    int idb; // debug level
    mpf_class eps;
    int effort; // effort level, >=1
    int pslq_max_iter;
    CRITERION crit;
};

class QasmGenerator
{
public:
    static void to_qasm(circuit &cir, std::ostream &os, std::string qbit_def = "", std::string ancil_def = "", bool without_header = true);

private:
    static std::vector<std::string> gen_name_qasm();

    static void fill_in_gatename(std::vector<std::string> &name_qasm, const std::string &gate_name);
    static std::string split_by_semicolon(const std::string &str);
};

#endif