#include <gmp.h>
#include "Matrix.h"
#include "MPList.h"

using namespace std;

/* Usage: "<prog> <input hollow symmetric matrix file>" */

MPList<mpz_t> CompleteGraph3Cycles(const Matrix<mpz_t>& M);

int main(int argc, char* argv[]){
	if(argc != 2){
		fprintf(stderr, "Usage: %s <input matrix>\n", argv[0]);
		return 1;
	}
	Matrix<mpz_t> A = Matrix<mpz_t>::FromFileHollowSymmetric(argv[1]);
	printf("A =\n");
	A.Print();
	printf("\n");

	MPList<mpz_t> terms = CompleteGraph3Cycles(A);
	terms.Print("\n");

	return 0;
}

MPList<mpz_t> CompleteGraph3Cycles(const Matrix<mpz_t>& A){
	Matrix<mpz_t> _T1(A.Rows, A.Columns);
	Matrix<mpz_t> _T2(A.Rows, A.Columns);
	MPList<mpz_t> terms(1);
	
	A.RightMultiply(_T1, A);
	Matrix<mpz_t> A3(A.Rows, A.Columns);
	A.RightMultiply(A3, _T1);
	
	// A3
	A3.Trace(terms[0]);
	
	return terms;
}