#include <gmp.h>
#include <cstdio>
#include "Matrix.h"
#include "MPList.h"
#include <vector>

using namespace std;

void basis4(const Matrix<mpz_t>& M, MPList<mpq_t>& terms);

int main(int argc, char* argv[]){
	int subGraphSize = 4;
	int size = subGraphSize;
	int dimension = 3; // how many connected even subgraphs (matrix expressions) of size 'subGraphSize'
	int totalTerms = 3;
	Matrix<mpq_t> system = Matrix<mpq_t>(totalTerms, totalTerms+1);
	MPList<mpq_t> terms = MPList<mpq_t>(totalTerms); 
	mpz_t edges, temp, val;
	mpz_inits(edges, temp, val, NULL);
	mpq_t valQ;
	mpq_init(valQ);

	for (int basisIndex = 0; basisIndex < totalTerms; basisIndex++) {
		Matrix<mpz_t> A = Matrix<mpz_t>::FromCompleteGraph(size);
		basis4(A,terms);
		terms.Print();
		for (int i = 0; i < totalTerms; i++){
			mpq_set(system.Index(basisIndex, i), terms[i]);
		}

		Matrix<mpq_t> A_q = Matrix<mpq_t>(A.Rows, A.Columns);
		A.ToRationalMatrix(A_q);
		Matrix<mpq_t>::Determinant(A_q, system.Index(basisIndex, totalTerms));
		mpq_out_str(stdout, 10, system.Index(basisIndex, totalTerms));
		printf("\n\n");

		size++;

	}
	system.Print();
	printf("\n");
	vector<int> pivots;
	Matrix<mpq_t>::RowEchelonForm(system, pivots);
	system.Print();
	printf("\n");
	system.BackSubstitution(pivots);
	system.Print();

	return 0;
}

/*  these are basis vectors for even graphs of 4 edges
 */
void basis4(const Matrix<mpz_t>& A, MPList<mpq_t>& terms){
	Matrix<mpz_t> _T1(A.Rows, A.Columns);
	Matrix<mpz_t> _T2(A.Rows, A.Columns);
	
	Matrix<mpz_t> A2(A.Rows, A.Columns);
	A.RightMultiply(A2, A);
	A.RightMultiply(_T2, A2);
	Matrix<mpz_t> A4(A.Rows, A.Columns);
	A.RightMultiply(A4, _T2);
	
	Matrix<mpz_t> A_2(A.Rows, A.Columns);
	A.MultiplyEntrywise(A_2, A);
	
	// A4
	A4.Trace(terms[0]);
	
	// d(A2)*d(A2)
	Matrix<mpz_t> _T3(A.Rows, A.Columns);
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	_T1.MultiplyEntrywise(_T3, _T2);
	_T3.Trace(terms[1]);
	
	// A_2@A_2
	A_2.RightMultiply(_T1, A_2);
	_T1.Trace(terms[2]);

	// // 2 + 2 disconnected
	// A2.Trace(terms[3]);
	// mpq_mul(terms[3], terms[3], terms[3]);

	// A2.RightMultiply(_T1, A);
	// _T1.Trace(terms[4]);

	// A2.Trace(terms[5]);

	// A.Trace(terms[6]);
}