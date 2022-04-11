#include <gmp.h>
#include <cstdio>
#include "Matrix.h"
#include "MPList.h"
#include <vector>
#include <fstream>

using namespace std;

void basis6(const Matrix<mpz_t>& M, MPList<mpq_t>& terms);

int main(int argc, char* argv[]){
	int subGraphSize = 6;
	int subGraphIndex = 1; // if the subgraph is a basis graph, what index it is. -1 for none
	int dimension = 10; // how many connected even subgraphs (matrix expressions) of size 'subGraphSize'
	Matrix<mpq_t> system = Matrix<mpq_t>(dimension, dimension + 1);
	MPList<mpq_t> terms = MPList<mpq_t>(dimension); 
	Matrix<mpz_t> cycle = Matrix<mpz_t>::FromCycle(subGraphSize);

	ifstream inFile("evenGraphs6Book.txt");
	string line;
	for (int i = 0; i < dimension; i++) {
		do {
			getline(inFile, line);
		} while (line.find_first_not_of(" ") == string::npos);
		Matrix<mpz_t> basisGraph = Matrix<mpz_t>::FromLineHollowSymmetric(line);
		Matrix<mpz_t> combined = cycle.Combine(basisGraph, false);
		printf("-----iteration %d-----\n", i+1);

		basis6(combined, terms);
		terms.Print();
		for (int j = 0; j < dimension; j++){
			mpq_set(system.Index(i, j), terms[j]);
		}

		mpq_set_ui(system.Index(i, dimension), i == subGraphIndex ? 2 : 1, 1);
	}
	inFile.close();

	printf("\n\n");
	system.Print();
	printf("\n\n");
	vector<int> pivots;
	Matrix<mpq_t>::RowEchelonForm(system, pivots);
	system.BackSubstitution(pivots);
	system.Print();

	return 0;
}

void basis6(const Matrix<mpz_t>& A, MPList<mpq_t>& terms){
	Matrix<mpz_t> _T1(A.Rows, A.Columns);
	Matrix<mpz_t> _T2(A.Rows, A.Columns);
	mpz_t term;
	mpz_init(term);
	
	Matrix<mpz_t> A2(A.Rows, A.Columns);
	A.RightMultiply(A2, A);
	Matrix<mpz_t> A3(A.Rows, A.Columns);
	A.RightMultiply(A3, A2);
	Matrix<mpz_t> A4(A.Rows, A.Columns);
	A.RightMultiply(A4, A3);
	A.RightMultiply(_T2, A4);
	Matrix<mpz_t> A6(A.Rows, A.Columns);
	A.RightMultiply(A6, _T2);
	
	Matrix<mpz_t> A_2(A.Rows, A.Columns);
	A.MultiplyEntrywise(A_2, A);
	Matrix<mpz_t> A_3(A.Rows, A.Columns);
	A.MultiplyEntrywise(A_3, A_2);
	A.MultiplyEntrywise(_T1, A_3);
	Matrix<mpz_t> A_5(A.Rows, A.Columns);
	A.MultiplyEntrywise(A_5, _T1);
	
	// A3*A3
	A3.MultiplyEntrywise(_T1, A3);
	_T1.Trace(term);
	mpq_set_z(terms[0], term);
	
	// A6
	A6.Trace(term);
	mpq_set_z(terms[1], term);
	
	// A_2@A_2@A_2
	A_2.RightMultiply(_T1, A_2);
	_T1.RightMultiply(_T2, A_2);
	_T2.Trace(term);
	mpq_set_z(terms[2], term);
	
	// A@d(A@d(A2)@A)@A
	A2.GetDiagonal(_T1);
	A.RightMultiply(_T2, _T1);
	_T2.RightMultiply(_T1, A);
	_T1.GetDiagonal(_T2);
	A.RightMultiply(_T1, _T2);
	_T1.RightMultiply(_T2, A);
	_T2.Trace(term);
	mpq_set_z(terms[3], term);
	
	// d(A2)@d(A2)@A2
	Matrix<mpz_t> _T3(A.Rows, A.Columns);
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	_T1.RightMultiply(_T3, _T2);
	_T3.RightMultiply(_T1, A2);
	_T1.Trace(term);
	mpq_set_z(terms[4], term);
	
	// A@A_5
	A.RightMultiply(_T1, A_5);
	_T1.Trace(term);
	mpq_set_z(terms[5], term);
	
	// d(A2)@A_2@A_2
	A2.GetDiagonal(_T1);
	_T1.RightMultiply(_T2, A_2);
	_T2.RightMultiply(_T1, A_2);
	_T1.Trace(term);
	mpq_set_z(terms[6], term);
	
	// A_3@A3
	A_3.RightMultiply(_T1, A3);
	_T1.Trace(term);
	mpq_set_z(terms[7], term);
	
	// (A2*A2)@A_2
	A2.MultiplyEntrywise(_T1, A2);
	_T1.RightMultiply(_T2, A_2);
	_T2.Trace(term);
	mpq_set_z(terms[8], term);
	
	// d(A2)@A4
	A2.GetDiagonal(_T1);
	_T1.RightMultiply(_T2, A4);
	_T2.Trace(term);
	mpq_set_z(terms[9], term);
	
}