#include <gmp.h>
#include <cstdio>
#include "Matrix.h"
#include "MPList.h"
#include <vector>

using namespace std;

void basis4(const Matrix<mpz_t>& M, MPList<mpq_t>& terms);

int main(int argc, char* argv[]){
	int subGraphSize = 4;
	int dimension = 3; // how many connected even subgraphs (matrix expressions) of size 'subGraphSize'
	Matrix<mpq_t> system = Matrix<mpq_t>(dimension, dimension + 1);
	MPList<mpq_t> terms = MPList<mpq_t>(dimension); 
	mpz_t edges, temp, val;
	mpz_inits(edges, temp, val, NULL);

	for (int i = 0; i < dimension; i++) {
		int size = subGraphSize + i;
		Matrix<mpz_t> A = Matrix<mpz_t>::FromCompleteGraph(size);
		printf("-----iteration %d-----\n", i+1);
		basis4(A,terms);
		terms.Print();
		for (int j = 0; j < dimension; j++){
			mpq_set(system.Index(i, j), terms[j]);
		}

		// # m-cycles in K_n = (n P m)/(2m)
		// m-permutations of n vertices. Rotations and reversal double counted, so divide by 2m
		mpz_set_ui(temp, size);
		mpz_set_ui(val, 1);
		for (int j = 0; j < subGraphSize; j++) {
			mpz_mul(val, val, temp);
			mpz_sub_ui(temp, temp, 1);
		}
		mpz_divexact_ui(val, val, 2*subGraphSize);
		mpq_set_z(system.Index(i, dimension), val);
	}
	printf("\n\n");
	system.Print();
	printf("\n\n");
	vector<int> pivots;
	Matrix<mpq_t>::RowEchelonForm(system, pivots);
	system.BackSubstitution(pivots);
	system.Print();

	return 0;
}

/*  these are basis vectors for even graphs of 4 edges
 binoculars
 paper #2
 book #2, equivalent to above
A2@d(A2)
 4-edge
 paper/book #3
 equivalent to above
A@A_3
 acorn (odd)
A@(A_2)@A */
void basis4(const Matrix<mpz_t>& A, MPList<mpq_t>& terms){
	Matrix<mpz_t> _T1(A.Rows, A.Columns);
	Matrix<mpz_t> _T2(A.Rows, A.Columns);
	mpz_t term;
	mpz_init(term);
	
	Matrix<mpz_t> A2(A.Rows, A.Columns);
	A.RightMultiply(A2, A);
	A.RightMultiply(_T2, A2);
	Matrix<mpz_t> A4(A.Rows, A.Columns);
	A.RightMultiply(A4, _T2);
	
	Matrix<mpz_t> A_2(A.Rows, A.Columns);
	A.MultiplyEntrywise(A_2, A);
	
	// A4
	A4.Trace(term);
	mpq_set_z(terms[0], term);
	
	// d(A2)*d(A2)
	Matrix<mpz_t> _T3(A.Rows, A.Columns);
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	_T1.MultiplyEntrywise(_T3, _T2);
	_T3.Trace(term);
	mpq_set_z(terms[1], term);
	
	// A_2@A_2 
	A_2.RightMultiply(_T1, A_2);
	_T1.Trace(term);
	mpq_set_z(terms[2], term);
	
}