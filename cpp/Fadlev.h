#ifndef FADLEV_H
#define FADLEV_H

#include <gmp.h>
#include "Matrix.h"
#include "MPList.h"

using namespace std;

/* 
The Faddeev-Leverrier algorithm for computing coefficients of the characteristic equation
The first parameter is the MPList of coefficients to be calculated, the last coefficient being the constant.
This coefficient is equal to the determinant of the input matrix, up to a minus sign. 
The second parameter will be where the inverse of the input matrix is stored, 
since it is computed with little extra cost.
Note that this algorithm supports only square matrices
*/

void Fadlev(MPList<mpz_t>& coef, Matrix<mpq_t>& Ainv, const Matrix<mpz_t>& A);

void Fadlev(MPList<mpq_t>& coef, Matrix<mpq_t>& Ainv, const Matrix<mpq_t>& A);

#endif