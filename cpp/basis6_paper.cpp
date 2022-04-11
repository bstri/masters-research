#include <gmp.h>
#include <cstdio>
#include "Matrix.h"
#include "MPList.h"
#include <fstream>
#include <string>

using namespace std;

void basis6(const Matrix<mpz_t>& M, MPList<mpz_t>& terms);

int main(int argc, char* argv[]){
	int subGraphSize = 6;
	int size = subGraphSize;
	int cycleIndex = 9;
	int dimension = 10; // how many connected even subgraphs (matrix expressions) of size 'subGraphSize'
	Matrix<mpq_t> system = Matrix<mpq_t>(dimension);
	MPList<mpq_t> values = MPList<mpq_t>(dimension);
	MPList<mpz_t> terms = MPList<mpz_t>(dimension); 
	MPList<mpq_t> termsQ = MPList<mpq_t>(dimension);
	mpz_t edges, temp, val;
	mpz_inits(edges, temp, val, NULL);
	mpq_t valQ;
	mpq_init(valQ);
	Matrix<mpz_t> cycle = Matrix<mpz_t>::FromCycle(size);

	ifstream inFile("evenGraphs6.txt");
	string line;
	int singular;

	for (int basisIndex = 0; basisIndex < dimension; basisIndex++) {
		do {
			getline(inFile, line);
		} while (line.find_first_not_of(" ") == string::npos);
		Matrix<mpz_t> basisGraph = Matrix<mpz_t>::FromLineHollowSymmetric(line, true);
		Matrix<mpz_t> combined = cycle.Combine(basisGraph, false);
		printf("-----iteration %d-----\n", basisIndex+1);
		// basisGraph.Print();
		// combined.Print();

		basis6(combined,terms);

		terms.Print();
		for (int i = 0; i < dimension; i++){
			mpq_set_z(system.Index(dimension - 1, i), terms[i]);
		}

		mpq_set_ui(values[dimension-1], basisIndex == cycleIndex ? 2 : 1, 1);

		singular = Matrix<mpq_t>::RowEchelonForm(system, values, true);

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
	inFile.close();
	MPList<mpq_t> solution = system.BackSubstitution(values);
	solution.Print();

	return 0;
}

void basis6(const Matrix<mpz_t>& A, MPList<mpz_t>& terms){
	Matrix<mpz_t> _T1(A.Rows, A.Columns);
	Matrix<mpz_t> _T2(A.Rows, A.Columns);
	
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
	
	// d(A2)*d(A2)*d(A2)
	Matrix<mpz_t> _T3(A.Rows, A.Columns);
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	_T1.MultiplyEntrywise(_T3, _T2);
	A2.GetDiagonal(_T1);
	_T3.MultiplyEntrywise(_T2, _T1);
	_T2.Trace(terms[0]);
	
	// A_3@A3
	A_3.RightMultiply(_T1, A3);
	_T1.Trace(terms[1]);
	
	// d(A2)*d(A_2@A_2)
	A2.GetDiagonal(_T1);
	A_2.RightMultiply(_T2, A_2);
	_T2.GetDiagonal(_T3);
	_T1.MultiplyEntrywise(_T2, _T3);
	_T2.Trace(terms[2]);
	
	// A_3@A_3
	A_3.RightMultiply(_T1, A_3);
	_T1.Trace(terms[3]);
	
	// d(A2)*d(A@d(A2)@A)
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	A.RightMultiply(_T3, _T2);
	_T3.RightMultiply(_T2, A);
	_T2.GetDiagonal(_T3);
	_T1.MultiplyEntrywise(_T2, _T3);
	_T2.Trace(terms[4]);
	
	// A_2@A_2@A_2
	A_2.RightMultiply(_T1, A_2);
	_T1.RightMultiply(_T2, A_2);
	_T2.Trace(terms[5]);
	
	// (A2*A2)@A_2
	A2.MultiplyEntrywise(_T1, A2);
	_T1.RightMultiply(_T2, A_2);
	_T2.Trace(terms[6]);
	
	// d(A3)*d(A3)
	A3.GetDiagonal(_T1);
	A3.GetDiagonal(_T2);
	_T1.MultiplyEntrywise(_T3, _T2);
	_T3.Trace(terms[7]);
	
	// d(A2)*d(A4)
	A2.GetDiagonal(_T1);
	A4.GetDiagonal(_T2);
	_T1.MultiplyEntrywise(_T3, _T2);
	_T3.Trace(terms[8]);
	
	// A6
	A6.Trace(terms[9]);
	
}