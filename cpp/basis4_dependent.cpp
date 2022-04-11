#include <gmp.h>
#include <cstdio>
#include "Matrix.h"
#include "MPList.h"

using namespace std;

void basis4(const Matrix<mpz_t>& M, MPList<mpz_t>& terms);

int main(int argc, char* argv[]){
	int subGraphSize = 4;
	int size = subGraphSize;
	int dimension = 3; // how many connected even subgraphs (matrix expressions) of size 'subGraphSize'
	Matrix<mpq_t> system = Matrix<mpq_t>(dimension);
	MPList<mpq_t> values = MPList<mpq_t>(dimension);
	MPList<mpz_t> terms = MPList<mpz_t>(dimension); 
	MPList<mpq_t> termsQ = MPList<mpq_t>(dimension);
	mpz_t edges, temp, val;
	mpz_inits(edges, temp, val, NULL);
	mpq_t valQ;
	mpq_init(valQ);

	while (true) {
		Matrix<mpz_t> A = Matrix<mpz_t>::FromCompleteGraph(size);
		basis4(A,terms);
		terms.Print();
		terms.ToRationalList(termsQ);
		for (int i = 0; i < dimension; i++){
			mpq_set(system.Index(dimension - 1, i), termsQ[i]);
		}

		// # m-cycles in K_n = (n P m)/(2m)
		// m-permutations of n vertices. Rotations and reversal double counted, so divide by 2m
		mpz_set_ui(temp, size);
		mpz_set_ui(val, 1);
		for (int i = 0; i < subGraphSize; i++) {
			mpz_mul(val, val, temp);
			mpz_sub_ui(temp, temp, 1);
		}
		mpz_divexact_ui(val, val, 2*subGraphSize);
		mpq_set_z(valQ, val);
		mpq_set(values[dimension-1], valQ);
		
		printf("num %d-cycles in K_%d: ", subGraphSize, size);
		mpImpl::out_str(stdout, 10, val);
		printf("\n");

		int singular = Matrix<mpq_t>::RowEchelonForm(system, values, true);

		// A.~Matrix<mpz_t>();

		if (!singular) {
			printf("\n");
			break;
		}
		if (size - subGraphSize + 1 >= dimension)
			printf("dependent\n\n");
		else
			printf("\n");

		size++;

		if (size > 2*dimension + subGraphSize){
			printf("INCONSISTENT?\n");
			return -1; // system might be inconsistent
		}
	}
	MPList<mpq_t> solution = system.BackSubstitution(values);
	solution.Print();

	return 0;
}

/*  these are basis vectors for even graphs of 4 edges
d(A2)*d(A2)
A2@d(A2)
A@(A_2)@A */
void basis4(const Matrix<mpz_t>& A, MPList<mpz_t>& terms){
	Matrix<mpz_t> _T1(A.Rows, A.Columns);
	Matrix<mpz_t> _T2(A.Rows, A.Columns);
	
	A.RightMultiply(_T1, A);
	A.RightMultiply(_T2, _T1);
	Matrix<mpz_t> A4(A.Rows, A.Columns);
	A.RightMultiply(A4, _T2);
	
	Matrix<mpz_t> A_2(A.Rows, A.Columns);
	A.MultiplyEntrywise(A_2, A);
	Matrix<mpz_t> A_3(A.Rows, A.Columns);
	A.MultiplyEntrywise(A_3, A_2);
	
	// A4
	A4.Trace(terms[0]);
	
	// A_2@A_2
	A_2.RightMultiply(_T1, A_2);
	_T1.Trace(terms[1]);
	
	// A@A_3
	A.RightMultiply(_T1, A_3);
	_T1.Trace(terms[2]);
	
}