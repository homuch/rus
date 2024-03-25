#ifndef __NUMERIC_DEFINITION_H__
#define __NUMERIC_DEFINITION_H__
#include <complex>
#include <gmpxx.h>
typedef mpf_class real_t;
typedef mpz_class int_t;
typedef std::complex<real_t> complex_t;

// add _mpf to the end if change real_t to mpf_class
static const auto SQRT2 = 1.41421356237309504880168872420969807856967187537694807317667973799_mpf;

static const auto OMEGA = std::complex<real_t>(SQRT2 / 2, SQRT2 / 2);

static const auto RHO = SQRT2;
static const auto NU = decltype(SQRT2)(SQRT2 + 1);

#endif // __NUMERIC_DEFINITION_H__