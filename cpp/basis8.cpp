#include <gmp.h>
#include <cstdio>
#include "Matrix.h"
#include "MPList.h"
#include <vector>
#include <fstream>

using namespace std;

void basis8(const Matrix<mpz_t>& M, MPList<mpq_t>& terms);

int main(int argc, char* argv[]){
	int subGraphSize = 8;
	int subGraphIndex = 6; // if the subgraph is a basis graph, what index it is. -1 for none
	int dimension = 67; // how many connected even subgraphs (matrix expressions) of size 'subGraphSize'
	Matrix<mpq_t> system = Matrix<mpq_t>(dimension, dimension + 1);
	MPList<mpq_t> terms = MPList<mpq_t>(dimension); 
	Matrix<mpz_t> cycle = Matrix<mpz_t>::FromCycle(subGraphSize);

	ifstream inFile("evenGraphs8.txt");
	string line;
	for (int i = 0; i < dimension; i++) {
		do {
			getline(inFile, line);
		} while (line.find_first_not_of(" ") == string::npos);
		Matrix<mpz_t> basisGraph = Matrix<mpz_t>::FromLineHollowSymmetric(line);
		Matrix<mpz_t> combined = cycle.Combine(basisGraph, false);
		printf("-----iteration %d-----\n", i+1);

		basis8(combined, terms);
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

void basis8(const Matrix<mpz_t>& A, MPList<mpq_t>& terms){
	Matrix<mpz_t> _T1(A.Rows, A.Columns);
	Matrix<mpz_t> _T2(A.Rows, A.Columns);
	mpz_t term, _t1, _t2;
	mpz_inits(term, _t1, _t2, NULL);
	
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
	
	// (A_2@A_2)*(A_2@A_2)
	A_2.RightMultiply(_T1, A_2);
	A_2.RightMultiply(_T2, A_2);
	_T1.MultiplyEntrywise(_T3, _T2);
	_T3.Trace(term);
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
	
	// tr(d(A2)*d(A2)*d(A2))*tr(A2)
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	_T1.MultiplyEntrywise(_T3, _T2);
	A2.GetDiagonal(_T1);
	_T3.MultiplyEntrywise(_T2, _T1);
	_T2.Trace(_t1);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[43], _t1);
	
	// tr(A_3@A3)*tr(A2)
	A_3.RightMultiply(_T1, A3);
	_T1.Trace(_t1);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[44], _t1);
	
	// tr(d(A2)*d(A_2@A_2))*tr(A2)
	A2.GetDiagonal(_T1);
	A_2.RightMultiply(_T2, A_2);
	_T2.GetDiagonal(_T3);
	_T1.MultiplyEntrywise(_T2, _T3);
	_T2.Trace(_t1);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[45], _t1);
	
	// tr(A_3@A_3)*tr(A2)
	A_3.RightMultiply(_T1, A_3);
	_T1.Trace(_t1);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[46], _t1);
	
	// tr(d(A2)*d(A@d(A2)@A))*tr(A2)
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	A.RightMultiply(_T3, _T2);
	_T3.RightMultiply(_T2, A);
	_T2.GetDiagonal(_T3);
	_T1.MultiplyEntrywise(_T2, _T3);
	_T2.Trace(_t1);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[47], _t1);
	
	// tr(A_2@A_2@A_2)*tr(A2)
	A_2.RightMultiply(_T1, A_2);
	_T1.RightMultiply(_T2, A_2);
	_T2.Trace(_t1);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[48], _t1);
	
	// tr((A2*A2)@A_2)*tr(A2)
	A2.MultiplyEntrywise(_T1, A2);
	_T1.RightMultiply(_T2, A_2);
	_T2.Trace(_t1);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[49], _t1);
	
	// tr(d(A3)*d(A3))*tr(A2)
	A3.GetDiagonal(_T1);
	A3.GetDiagonal(_T2);
	_T1.MultiplyEntrywise(_T3, _T2);
	_T3.Trace(_t1);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[50], _t1);
	
	// tr(d(A2)*d(A4))*tr(A2)
	A2.GetDiagonal(_T1);
	A4.GetDiagonal(_T2);
	_T1.MultiplyEntrywise(_T3, _T2);
	_T3.Trace(_t1);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[51], _t1);
	
	// tr(A6)*tr(A2)
	A6.Trace(_t1);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[52], _t1);
	
	// tr(A5)*tr(A3)
	A5.Trace(_t1);
	A3.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[53], _t1);
	
	// tr(A2*A3)*tr(A3)
	A2.MultiplyEntrywise(_T1, A3);
	_T1.Trace(_t1);
	A3.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[54], _t1);
	
	// tr(A2@A_3)*tr(A3)
	A2.RightMultiply(_T1, A_3);
	_T1.Trace(_t1);
	A3.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[55], _t1);
	
	// tr(A4)*tr(A4)
	A4.Trace(_t1);
	A4.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[56], _t1);
	
	// tr(d(A2)*d(A2))*tr(d(A2)*d(A2))
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	_T1.MultiplyEntrywise(_T3, _T2);
	_T3.Trace(_t1);
	A2.GetDiagonal(_T3);
	A2.GetDiagonal(_T1);
	_T3.MultiplyEntrywise(_T2, _T1);
	_T2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[57], _t1);
	
	// tr(A_2@A_2)*tr(A_2@A_2)
	A_2.RightMultiply(_T1, A_2);
	_T1.Trace(_t1);
	A_2.RightMultiply(_T1, A_2);
	_T1.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[58], _t1);
	
	// tr(A4)*tr(d(A2)*d(A2))
	A4.Trace(_t1);
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	_T1.MultiplyEntrywise(_T3, _T2);
	_T3.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[59], _t1);
	
	// tr(A4)*tr(A_2@A_2)
	A4.Trace(_t1);
	A_2.RightMultiply(_T1, A_2);
	_T1.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[60], _t1);
	
	// tr(d(A2)*d(A2))*tr(A_2@A_2)
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	_T1.MultiplyEntrywise(_T3, _T2);
	_T3.Trace(_t1);
	A_2.RightMultiply(_T3, A_2);
	_T3.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[61], _t1);
	
	// tr(A4)*tr(A2)*tr(A2)
	A4.Trace(_t1);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[62], _t1);
	
	// tr(d(A2)*d(A2))*tr(A2)*tr(A2)
	A2.GetDiagonal(_T1);
	A2.GetDiagonal(_T2);
	_T1.MultiplyEntrywise(_T3, _T2);
	_T3.Trace(_t1);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[63], _t1);
	
	// tr(A_2@A_2)*tr(A2)*tr(A2)
	A_2.RightMultiply(_T1, A_2);
	_T1.Trace(_t1);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[64], _t1);
	
	// tr(A3)*tr(A3)*tr(A2)
	A3.Trace(_t1);
	A3.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[65], _t1);
	
	// tr(A2)*tr(A2)*tr(A2)*tr(A2)
	A2.Trace(_t1);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	A2.Trace(_t2);
	mpz_mul(_t1, _t1, _t2);
	mpq_set_z(terms[66], _t1);
	
}