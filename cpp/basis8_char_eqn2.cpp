#include <gmp.h>
#include <cstdio>
#include "Matrix.h"
#include "MPList.h"
#include <fstream>
#include <string>
#include "Fadlev.h"
#include <vector>

using namespace std;

void basis8(const Matrix<mpz_t>& M, MPList<mpq_t>& terms);

int main(int argc, char* argv[]){
	int subGraphSize = 8;
	int cycleIndex = 6;
	int size = subGraphSize;
	int dimension = 43; // how many connected even subgraphs (matrix expressions) of size 'subGraphSize'
	int totalGraphs = 67;
	Matrix<mpq_t> system = Matrix<mpq_t>(totalGraphs, totalGraphs + 1);
	MPList<mpq_t> terms = MPList<mpq_t>(totalGraphs); 
	mpz_t edges, temp, val;
	mpz_inits(edges, temp, val, NULL);
	Matrix<mpz_t> k8 = Matrix<mpz_t>::FromCompleteGraph(size);

	ifstream inFile("evenGraphs8.txt");
	string line;

	for (int basisIndex = 0; basisIndex < dimension; basisIndex++) {
		do {
			getline(inFile, line);
		} while (line.find_first_not_of(" ") == string::npos);
		Matrix<mpz_t> basisGraph = Matrix<mpz_t>::FromLineHollowSymmetric(line);
		Matrix<mpz_t> combined = k8.Combine(basisGraph, false);
		printf("-----iteration %d-----\n", basisIndex+1);
		// basisGraph.Print();
		// combined.Print();

		basis8(combined,terms);

		terms.Print();
		for (int i = 0; i < totalGraphs; i++){
			mpq_set(system.Index(basisIndex, i), terms[i]);
		}

		Matrix<mpq_t> Ainv = Matrix<mpq_t>(combined.Rows, combined.Columns);
		MPList<mpz_t> coefs = MPList<mpz_t>(combined.Rows + 1);
		Fadlev(coefs, Ainv, combined);
		mpq_set_z(system.Index(basisIndex, totalGraphs), coefs[8]);
	}
	inFile.close();

	ifstream inFile2("disconnectedGraphs8.txt");
	for (int i = dimension; i < totalGraphs - 1; i++) {
		do {
			getline(inFile2, line);
		} while (line.find_first_not_of(" ") == string::npos);
		Matrix<mpz_t> basisGraph = Matrix<mpz_t>::FromLineHollowSymmetric(line);
		Matrix<mpz_t> combined = k8.Combine(basisGraph, false);
		printf("-----iteration %d-----\n", i+1);

		basis8(combined,terms);

		terms.Print();
		for (int j = 0; j < totalGraphs; j++){
			mpq_set(system.Index(i, j), terms[j]);
		}

		Matrix<mpq_t> Ainv = Matrix<mpq_t>(combined.Rows, combined.Columns);
		MPList<mpz_t> coefs = MPList<mpz_t>(combined.Rows + 1);
		Fadlev(coefs, Ainv, combined);
		mpq_set_z(system.Index(i, totalGraphs), coefs[8]);
	}
	inFile2.close();

	printf("\n\n");
	system.Print();
	printf("\n\n");
	vector<int> pivots;
	Matrix<mpq_t>::RowEchelonForm(system, pivots);
	// system.Print();
	// printf("\n\n");
	system.BackSubstitution(pivots);
	system.Print();

	// MPList<mpq_t> solution = system.BackSubstitution(values);
	// printf("here2\n");
	// solution.Print();

	return 0;
}

// regex I used to handle converting the trace stuff
// (Trace\(term)s\[(\d+)\]\);
// $1);\n\tmpq_set_z(terms[$2], term);
void basis8(const Matrix<mpz_t>& A, MPList<mpq_t>& terms){
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
	Matrix<mpz_t> A5(A.Rows, A.Columns);
	A.RightMultiply(A5, A4);
	Matrix<mpz_t> A6(A.Rows, A.Columns);
	A.RightMultiply(A6, A5);
	A.RightMultiply(_T2, A6);
	Matrix<mpz_t> A8(A.Rows, A.Columns);
	A.RightMultiply(A8, _T2);
	
	Matrix<mpz_t> A_2(A.Rows, A.Columns);
	A.MultiplyEntrywise(A_2, A);
	Matrix<mpz_t> A_3(A.Rows, A.Columns);
	A.MultiplyEntrywise(A_3, A_2);
	Matrix<mpz_t> A_4(A.Rows, A.Columns);
	A.MultiplyEntrywise(A_4, A_3);
	Matrix<mpz_t> A_5(A.Rows, A.Columns);
	A.MultiplyEntrywise(A_5, A_4);
	A.MultiplyEntrywise(_T1, A_5);
	Matrix<mpz_t> A_7(A.Rows, A.Columns);
	A.MultiplyEntrywise(A_7, _T1);
	
	// (A@d(A2)@A)*(A@d(A2)@A)
	Matrix<mpz_t> _T3(A.Rows, A.Columns);
	A2.GetDiagonal(_T1);
	A.RightMultiply(_T2, _T1);
	_T2.RightMultiply(_T1, A);
	A2.GetDiagonal(_T2);
	A.RightMultiply(_T3, _T2);
	_T3.RightMultiply(_T2, A);
	_T1.MultiplyEntrywise(_T3, _T2);
	_T3.Trace(term);
	mpq_set_z(terms[0], term);
	
	// A3@d(A@d(A3)@A)
	A3.GetDiagonal(_T1);
	A.RightMultiply(_T2, _T1);
	_T2.RightMultiply(_T1, A);
	_T1.GetDiagonal(_T2);
	A3.RightMultiply(_T1, _T2);
	_T1.Trace(term);
	mpq_set_z(terms[1], term);
	
	// (A2*A)@(A2*A3)
	A2.MultiplyEntrywise(_T1, A);
	A2.MultiplyEntrywise(_T2, A3);
	_T1.RightMultiply(_T3, _T2);
	_T3.Trace(term);
	mpq_set_z(terms[2], term);
	
	// (A2*A2*A2)@A2
	A2.MultiplyEntrywise(_T1, A2);
	_T1.MultiplyEntrywise(_T2, A2);
	_T2.RightMultiply(_T1, A2);
	_T1.Trace(term);
	mpq_set_z(terms[3], term);
	
	// A4@d(A4)
	A4.GetDiagonal(_T1);
	A4.RightMultiply(_T2, _T1);
	_T2.Trace(term);
	mpq_set_z(terms[4], term);
	
	// A5@d(A3)
	A3.GetDiagonal(_T1);
	A5.RightMultiply(_T2, _T1);
	_T2.Trace(term);
	mpq_set_z(terms[5], term);
	
	// A8
	A8.Trace(term);
	mpq_set_z(terms[6], term);
	
	// A_2@(A2*A2)@A_2
	A2.MultiplyEntrywise(_T1, A2);
	A_2.RightMultiply(_T2, _T1);
	_T2.RightMultiply(_T1, A_2);
	_T1.Trace(term);
	mpq_set_z(terms[7], term);
	
	// A_2@(A3*A)@A_2
	A3.MultiplyEntrywise(_T1, A);
	A_2.RightMultiply(_T2, _T1);
	_T2.RightMultiply(_T1, A_2);
	_T1.Trace(term);
	mpq_set_z(terms[8], term);
	
	// A_2@A_2@d(A4)
	A_2.RightMultiply(_T1, A_2);
	A4.GetDiagonal(_T2);
	_T1.RightMultiply(_T3, _T2);
	_T3.Trace(term);
	mpq_set_z(terms[9], term);
	
	// A2@d(A2@d(A2)@A2)
	A2.GetDiagonal(_T1);
	A2.RightMultiply(_T2, _T1);
	_T2.RightMultiply(_T1, A2);
	_T1.GetDiagonal(_T2);
	A2.RightMultiply(_T1, _T2);
	_T1.Trace(term);
	mpq_set_z(terms[10], term);
	
	// A@d(A2)@A@d(A2)@A2
	A2.GetDiagonal(_T1);
	A.RightMultiply(_T2, _T1);
	_T2.RightMultiply(_T1, A);
	A2.GetDiagonal(_T2);
	_T1.RightMultiply(_T3, _T2);
	_T3.RightMultiply(_T1, A2);
	_T1.Trace(term);
	mpq_set_z(terms[11], term);
	
	// (A@d(A2)@A)*A4
	A2.GetDiagonal(_T1);
	A.RightMultiply(_T2, _T1);
	_T2.RightMultiply(_T1, A);
	_T1.MultiplyEntrywise(_T2, A4);
	_T2.Trace(term);
	mpq_set_z(terms[12], term);
	
	// d(A2)@d(A2)@A4
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	_T1.RightMultiply(_T3, _T2);
	_T3.RightMultiply(_T1, A4);
	_T1.Trace(term);
	mpq_set_z(terms[13], term);
	
	// A_4@(A2*A2)
	A2.MultiplyEntrywise(_T1, A2);
	A_4.RightMultiply(_T2, _T1);
	_T2.Trace(term);
	mpq_set_z(terms[14], term);
	
	// A_2@(A2*(A@d(A2)@A))
	A2.GetDiagonal(_T1);
	A.RightMultiply(_T2, _T1);
	_T2.RightMultiply(_T1, A);
	A2.MultiplyEntrywise(_T2, _T1);
	A_2.RightMultiply(_T1, _T2);
	_T1.Trace(term);
	mpq_set_z(terms[15], term);
	
	// A_2@d(A2)@(A2*A2)
	A2.GetDiagonal(_T1);
	A_2.RightMultiply(_T2, _T1);
	A2.MultiplyEntrywise(_T1, A2);
	_T2.RightMultiply(_T3, _T1);
	_T3.Trace(term);
	mpq_set_z(terms[16], term);
	
	// A_3@A@A_3@A
	A_3.RightMultiply(_T1, A);
	_T1.RightMultiply(_T2, A_3);
	_T2.RightMultiply(_T1, A);
	_T1.Trace(term);
	mpq_set_z(terms[17], term);
	
	// A_3@A2@A_3
	A_3.RightMultiply(_T1, A2);
	_T1.RightMultiply(_T2, A_3);
	_T2.Trace(term);
	mpq_set_z(terms[18], term);
	
	// ((A_3@A)*A_2)@A2
	A_3.RightMultiply(_T1, A);
	_T1.MultiplyEntrywise(_T2, A_2);
	_T2.RightMultiply(_T1, A2);
	_T1.Trace(term);
	mpq_set_z(terms[19], term);
	
	// A_5@A3
	A_5.RightMultiply(_T1, A3);
	_T1.Trace(term);
	mpq_set_z(terms[20], term);
	
	// A_3@A@d(A2)@A2
	A_3.RightMultiply(_T1, A);
	A2.GetDiagonal(_T2);
	_T1.RightMultiply(_T3, _T2);
	_T3.RightMultiply(_T1, A2);
	_T1.Trace(term);
	mpq_set_z(terms[21], term);
	
	// d(A2)@A_3@A3
	A2.GetDiagonal(_T1);
	_T1.RightMultiply(_T2, A_3);
	_T2.RightMultiply(_T1, A3);
	_T1.Trace(term);
	mpq_set_z(terms[22], term);
	
	// A_3@A@d(A_2@A_2)
	A_3.RightMultiply(_T1, A);
	A_2.RightMultiply(_T2, A_2);
	_T2.GetDiagonal(_T3);
	_T1.RightMultiply(_T2, _T3);
	_T2.Trace(term);
	mpq_set_z(terms[23], term);
	
	// A_7@A
	A_7.RightMultiply(_T1, A);
	_T1.Trace(term);
	mpq_set_z(terms[24], term);
	
	// d(A2)@A_5@A
	A2.GetDiagonal(_T1);
	_T1.RightMultiply(_T2, A_5);
	_T2.RightMultiply(_T1, A);
	_T1.Trace(term);
	mpq_set_z(terms[25], term);
	
	// d(A2)@d(A2)@A_3@A
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	_T1.RightMultiply(_T3, _T2);
	_T3.RightMultiply(_T1, A_3);
	_T1.RightMultiply(_T3, A);
	_T3.Trace(term);
	mpq_set_z(terms[26], term);
	
	// d(A2)@d(A2)@A@d(A2)@A
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	_T1.RightMultiply(_T3, _T2);
	_T3.RightMultiply(_T1, A);
	A2.GetDiagonal(_T3);
	_T1.RightMultiply(_T2, _T3);
	_T2.RightMultiply(_T1, A);
	_T1.Trace(term);
	mpq_set_z(terms[27], term);
	
	// d(A2)@d(A2)@d(A2)@A2
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	_T1.RightMultiply(_T3, _T2);
	A2.GetDiagonal(_T1);
	_T3.RightMultiply(_T2, _T1);
	_T2.RightMultiply(_T3, A2);
	_T3.Trace(term);
	mpq_set_z(terms[28], term);
	
	// A_2@A_2@A_2@A_2
	A_2.RightMultiply(_T1, A_2);
	_T1.RightMultiply(_T2, A_2);
	_T2.RightMultiply(_T1, A_2);
	_T1.Trace(term);
	mpq_set_z(terms[29], term);
	
	// d(A_2@A_2)@A@d(A2)@A
	A_2.RightMultiply(_T1, A_2);
	_T1.GetDiagonal(_T2);
	_T2.RightMultiply(_T1, A);
	A2.GetDiagonal(_T2);
	_T1.RightMultiply(_T3, _T2);
	_T3.RightMultiply(_T1, A);
	_T1.Trace(term);
	mpq_set_z(terms[30], term);
	
	// d(A2)@A_3@d(A2)@A
	A2.GetDiagonal(_T1);
	_T1.RightMultiply(_T2, A_3);
	A2.GetDiagonal(_T1);
	_T2.RightMultiply(_T3, _T1);
	_T3.RightMultiply(_T2, A);
	_T2.Trace(term);
	mpq_set_z(terms[31], term);
	
	// A_4@A_2@A_2
	A_4.RightMultiply(_T1, A_2);
	_T1.RightMultiply(_T2, A_2);
	_T2.Trace(term);
	mpq_set_z(terms[32], term);
	
	// d(A_2@A_2@A_2)@A2
	A_2.RightMultiply(_T1, A_2);
	_T1.RightMultiply(_T2, A_2);
	_T2.GetDiagonal(_T1);
	_T1.RightMultiply(_T2, A2);
	_T2.Trace(term);
	mpq_set_z(terms[33], term);
	
	// A_3@A5
	A_3.RightMultiply(_T1, A5);
	_T1.Trace(term);
	mpq_set_z(terms[34], term);
	
	// (A2*A_2)@A4
	A2.MultiplyEntrywise(_T1, A_2);
	_T1.RightMultiply(_T2, A4);
	_T2.Trace(term);
	mpq_set_z(terms[35], term);
	
	// A_2@(A3*A3)
	A3.MultiplyEntrywise(_T1, A3);
	A_2.RightMultiply(_T2, _T1);
	_T2.Trace(term);
	mpq_set_z(terms[36], term);
	
	// d(A2)@A6
	A2.GetDiagonal(_T1);
	_T1.RightMultiply(_T2, A6);
	_T2.Trace(term);
	mpq_set_z(terms[37], term);
	
	// A_3@A@d(A3)@A
	A_3.RightMultiply(_T1, A);
	A3.GetDiagonal(_T2);
	_T1.RightMultiply(_T3, _T2);
	_T3.RightMultiply(_T1, A);
	_T1.Trace(term);
	mpq_set_z(terms[38], term);
	
	// ((A_2@(A2*A))*A)@A2
	A2.MultiplyEntrywise(_T1, A);
	A_2.RightMultiply(_T2, _T1);
	_T2.MultiplyEntrywise(_T1, A);
	_T1.RightMultiply(_T2, A2);
	_T2.Trace(term);
	mpq_set_z(terms[39], term);
	
	// (A_3@A2)*A3
	A_3.RightMultiply(_T1, A2);
	_T1.MultiplyEntrywise(_T2, A3);
	_T2.Trace(term);
	mpq_set_z(terms[40], term);
	
	// d(A3)@A@d(A2)@A2
	A3.GetDiagonal(_T1);
	_T1.RightMultiply(_T2, A);
	A2.GetDiagonal(_T1);
	_T2.RightMultiply(_T3, _T1);
	_T3.RightMultiply(_T2, A2);
	_T2.Trace(term);
	mpq_set_z(terms[41], term);
	
	// d(A2)@d(A3)@A3
	A2.GetDiagonal(_T1);
	A3.GetDiagonal(_T2);
	_T1.RightMultiply(_T3, _T2);
	_T3.RightMultiply(_T1, A3);
	_T1.Trace(term);
	mpq_set_z(terms[42], term);

	// c_8(A), coef on lambda^8 term
	// Matrix<mpq_t> Ainv = Matrix<mpq_t>(A.Rows, A.Columns);
	// MPList<mpz_t> coefs = MPList<mpz_t>(A.Rows + 1);
	// Fadlev(coefs, Ainv, A);
	// mpq_set_z(terms[44], coefs[A.Rows - 8]);

	// 6 + 2 disconnected

	mpz_t traceA2;
	mpz_init(traceA2);
	A2.Trace(traceA2);

	// d(A2)*d(A2)*d(A2)
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	_T1.MultiplyEntrywise(_T3, _T2);
	A2.GetDiagonal(_T1);
	_T3.MultiplyEntrywise(_T2, _T1);
	_T2.Trace(term);
	mpz_mul(term, term, traceA2);
	mpq_set_z(terms[43], term);
	
	// A_3@A3
	A_3.RightMultiply(_T1, A3);
	_T1.Trace(term);
	mpz_mul(term, term, traceA2);
	mpq_set_z(terms[44], term);
	
	// d(A2)*d(A_2@A_2)
	A2.GetDiagonal(_T1);
	A_2.RightMultiply(_T2, A_2);
	_T2.GetDiagonal(_T3);
	_T1.MultiplyEntrywise(_T2, _T3);
	_T2.Trace(term);
	mpz_mul(term, term, traceA2);
	mpq_set_z(terms[45], term);
	
	// A_3@A_3
	A_3.RightMultiply(_T1, A_3);
	_T1.Trace(term);
	mpz_mul(term, term, traceA2);
	mpq_set_z(terms[46], term);
	
	// d(A2)*d(A@d(A2)@A)
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	A.RightMultiply(_T3, _T2);
	_T3.RightMultiply(_T2, A);
	_T2.GetDiagonal(_T3);
	_T1.MultiplyEntrywise(_T2, _T3);
	_T2.Trace(term);
	mpz_mul(term, term, traceA2);
	mpq_set_z(terms[47], term);
	
	// A_2@A_2@A_2
	A_2.RightMultiply(_T1, A_2);
	_T1.RightMultiply(_T2, A_2);
	_T2.Trace(term);
	mpz_mul(term, term, traceA2);
	mpq_set_z(terms[48], term);
	
	// (A2*A2)@A_2
	A2.MultiplyEntrywise(_T1, A2);
	_T1.RightMultiply(_T2, A_2);
	_T2.Trace(term);
	mpz_mul(term, term, traceA2);
	mpq_set_z(terms[49], term);
	
	// d(A3)*d(A3)
	A3.GetDiagonal(_T1);
	A3.GetDiagonal(_T2);
	_T1.MultiplyEntrywise(_T3, _T2);
	_T3.Trace(term);
	mpz_mul(term, term, traceA2);
	mpq_set_z(terms[50], term);
	
	// d(A2)*d(A4)
	A2.GetDiagonal(_T1);
	A4.GetDiagonal(_T2);
	_T1.MultiplyEntrywise(_T3, _T2);
	_T3.Trace(term);
	mpz_mul(term, term, traceA2);
	mpq_set_z(terms[51], term);
	
	// A6
	A6.Trace(term);
	mpz_mul(term, term, traceA2);
	mpq_set_z(terms[52], term);

	// 5 + 3 disconnected

	mpz_t traceA3;
	mpz_init(traceA3);
	A3.Trace(traceA3);

	// A5
	A5.Trace(term);
	mpz_mul(term, term, traceA3);
	mpq_set_z(terms[53], term);

	// A2*A3
	A2.MultiplyEntrywise(_T1, A3);
	_T1.Trace(term);
	mpz_mul(term, term, traceA3);
	mpq_set_z(terms[54], term);
    
	// A2@A_3
	A2.RightMultiply(_T1, A_3);
	_T1.Trace(term);
	mpz_mul(term, term, traceA3);
	mpq_set_z(terms[55], term);

	// 4 + 4 disconnected, and 4 + 2 + 2 disconnected

	// A4
	A4.Trace(term);
	mpq_set_z(terms[56], term);
	mpz_mul(term, term, traceA2);
	mpz_mul(term, term, traceA2);
	mpq_set_z(terms[62], term);
	
	// d(A2)*d(A2)
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	_T1.MultiplyEntrywise(_T3, _T2);
	_T3.Trace(term);
	mpq_set_z(terms[57], term);
	mpz_mul(term, term, traceA2);
	mpz_mul(term, term, traceA2);
	mpq_set_z(terms[63], term);
	
	// A_2@A_2
	A_2.RightMultiply(_T1, A_2);
	_T1.Trace(term);
	mpq_set_z(terms[58], term);
	mpz_mul(term, term, traceA2);
	mpz_mul(term, term, traceA2);
	mpq_set_z(terms[64], term);

	mpq_mul(terms[59], terms[56], terms[57]);
	mpq_mul(terms[60], terms[56], terms[58]);
	mpq_mul(terms[61], terms[57], terms[58]);
	mpq_mul(terms[56], terms[56], terms[56]);
	mpq_mul(terms[57], terms[57], terms[57]);
	mpq_mul(terms[58], terms[58], terms[58]);

	// 3 + 3 + 2 disconnected

	mpz_mul(term, traceA3, traceA3);
	mpz_mul(term, term, traceA2);
	mpq_set_z(terms[65], term);

	// 2 + 2 + 2 + 2 disconnected

	mpz_mul(term, traceA2, traceA2);
	mpz_mul(term, term, traceA2);
	mpz_mul(term, term, traceA2);
	mpq_set_z(terms[66], term);
}