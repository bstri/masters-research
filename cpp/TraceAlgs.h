#ifndef TRACE_ALGS_H
#define TRACE_ALGS_H

#include <gmp.h>
#include "Matrix.h"
#include "MPList.h"

using namespace std;

MPList<mpq_t> TraceAlgorithm(mpz_t out, const Matrix<mpz_t>& A);

MPList<mpq_t> TraceAlgorithm1Optimized(mpz_t out, const Matrix<mpz_t>& A);

MPList<mpz_t> TraceAlgorithm2(mpz_t out, const Matrix<mpz_t>& A);

MPList<mpz_t> TraceAlgorithm3(const Matrix<mpz_t>& A);

MPList<mpz_t> CompiledTraceAlgorithm3(const Matrix<mpz_t>& A);

MPList<mpz_t> CompiledTraceAlgorithm(const Matrix<mpz_t>& M);

MPList<mpz_t> CompiledTraceAlgorithmShorthand(const Matrix<mpz_t>& M);

// Counts number of 3-cycles in complete graphs (M is incidence matrix)
MPList<mpz_t> CompleteGraph3Cycles(const Matrix<mpz_t>& M);

#endif