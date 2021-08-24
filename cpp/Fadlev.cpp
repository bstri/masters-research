#include "Fadlev.h"
#include "Matrix.h"
#include "MPList.h"
#include <cstdlib>
#include <gmp.h>

using namespace std;

void Fadlev(MPList<mpz_t>& c, Matrix<mpq_t>& Ainv, const Matrix<mpz_t>& A){
	const int n = A.Rows;
	
	mpz_set_ui(c[0], 1);
	Matrix<mpz_t> B(A);
	Matrix<mpz_t> M(n);

	for(unsigned int i = 1; i < (unsigned int) n; i++){
		// c[i] = -trace(B)/i
		B.Trace(c[i]);
		mpz_divexact_ui(c[i], c[i], i);
		mpz_neg(c[i], c[i]);

		// M = B + c[i]*I
		for(int j = 0; j < n; j++){
			for(int k = 0; k < j; k++)
				mpz_set(M.Index(j,k), B.Index(j,k));
			for(int k = j; k < n; k++)
				mpz_set(M.Index(j,k), B.Index(j,k));
			mpz_add(M.Index(j,j), B.Index(j,j), c[i]);
		}

		// B = A*M
		A.RightMultiply(B, M);
	}
	B.Trace(c[n]);
	mpz_divexact_ui(c[n], c[n], n);
	mpz_neg(c[n], c[n]);

	// Ainv = -M/c[n]
	// Make sure it's invertible
	if(mpz_sgn(c[n]) != 0){
		mpq_t d;
		mpq_init(d);
		mpq_set_z(d, c[n]);
		mpq_inv(d, d);
		mpq_neg(d, d);
		Ainv.Copy(M);
		Ainv.Scale(d);
	}
}

void Fadlev(MPList<mpq_t>& c, Matrix<mpq_t>& Ainv, const Matrix<mpq_t>& A){
	const int n = A.Rows;
	
	mpq_set_ui(c[0], 1, 1);
	Matrix<mpq_t> B(A);
	Matrix<mpq_t> M(n);

	mpq_t div;
	mpq_init(div);
	for(int i = 1; i < n; i++){
		// c[i] = -trace(B)/i
		B.Trace(c[i]);
		mpq_set_si(div, -1, i);
		mpq_mul(c[i], c[i], div);

		// M = B + c[i]*I
		for(int j = 0; j < n; j++){
			for(int k = 0; k < j; k++)
				mpq_set(M.Index(j,k), B.Index(j,k));
			for(int k = j; k < n; k++)
				mpq_set(M.Index(j,k), B.Index(j,k));
			mpq_add(M.Index(j,j), B.Index(j,j), c[i]);
		}

		// B = A*M
		A.RightMultiply(B, M);
	}
	B.Trace(c[n]);
	mpq_set_si(div, -1, n);
	mpq_mul(c[n], c[n], div);

	// Ainv = -M/c[n]
	// Make sure it's invertible
	if(mpq_sgn(c[n]) != 0){
		mpq_set(div, c[n]);
		mpq_inv(div, div);
		mpq_neg(div, div);
		Ainv.Copy(M);
		Ainv.Scale(div);
	}
}