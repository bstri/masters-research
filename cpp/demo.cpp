#include <cstdio> // for printf
#include <gmp.h>
#include <ctime> // for time, clock
#include <cstdlib> // for atoi
#include <cctype>
#include "Matrix.h"
#include "Fadlev.h"
#include "MPList.h"

using namespace std;

void Demo(const Matrix<mpz_t>& A, const Matrix<mpz_t>& B);

int main(int argc, char* argv[]){
	if(argc != 2){
		fprintf(stderr, "Error. Usage: %s <int | fileName>", argv[0]);
		return -1;
	}

	mpz_t maxRnd, rnd;
	mpz_inits(maxRnd, rnd, NULL);
	mpz_set_ui(maxRnd, 3);
	gmp_randstate_t r;
	gmp_randinit_default(r);
	gmp_randseed_ui(r, time(NULL));

	if(isdigit(argv[1][0])){
		int n = atoi(argv[1]);
		Matrix<mpz_t> A = Matrix<mpz_t>::Random(3, 5, r, maxRnd);
		Matrix<mpz_t> B = Matrix<mpz_t>::Random(n, r, maxRnd);
		Demo(A, B);
		
	} else {
		Matrix<mpz_t> A = Matrix<mpz_t>::FromFile(argv[1]);
		Matrix<mpz_t> B = Matrix<mpz_t>::Random(A.Rows, r, maxRnd);
		Demo(A, B);
	}

	return 0;
}

void Demo(const Matrix<mpz_t>& A, const Matrix<mpz_t>& B){
	int n = A.Rows;

	printf("A = \n");
	// A.Print();

	// printf("\nB = \n");
	// // B.Print();
	// printf("\n");

	// Matrix<mpz_t> M(n);
	// Matrix<mpq_t> Ainv(n);

	// printf("A*B = \n");
	// clock_t t1 = clock();
	// A.RightMultiply(M, B);
	// t1 = clock() - t1;
	// // M.Print();

	// A.MultiplyEntrywise(M, B);
	// printf("\nHadamard (entrywise) product of A and B = \n");
	// M.Print();

	// printf("\nFadlev(A) = \n");
	// MPList<mpz_t> c(n + 1);
	// clock_t t2 = clock();
	// Fadlev(c, Ainv, A);
	// t2 = clock() - t2;
	// c.Print();
	// printf("Constant term from Fadlev(A) = \n");
	// mpz_out_str(stdout, 10, c[n]);
	// printf("\nDet(A) from Gaussian elimination = ");
	mpq_t det;
	mpq_init(det);
	Matrix<float> q(A.Rows, A.Columns);
	// Matrix<mpq_t> q(A.Rows, A.Columns);
	A.ToRationalMatrix(q);
	q.Print();
	// Matrix<float>::RowEchelonForm(q, false);
	// printf("\nAfter row echelon form...\n\n");
	// q.Print();
	// float _det;
	// clock_t tDet = clock();
	// // Matrix<mpq_t>::Determinant(q, det);
	// Matrix<float>::Determinant(q, _det);
	// tDet = clock() - tDet;
	// printf("%.0f\n", _det);
	// mpq_out_str(stdout, 10, det);

	// printf("\nInv(A) = \n");
	// Ainv.Print();

	// mpz_t trace;
	// mpz_init(trace);
	// A.Trace(trace);
	// printf("\nTrace(A) = ");
	// mpz_out_str(stdout, 10, trace);

	// Matrix<mpz_t> D = Matrix<mpz_t>::Diagonal(A);
	// printf("\n\nDiag(A) = \n");
	// // D.Print();

	// Matrix<mpz_t> L = Matrix<mpz_t>::Lower(A);
	// printf("\nL = lower triangular of A = \n");
	// L.Print();
	// printf("\nFadlev(L) = \n");
	// clock_t t3 = clock();
	// Fadlev(c, Ainv, L);
	// t3 = clock() - t3;
	// c.Print();
	// printf("Constant term from Fadlev(L) = \n");
	// mpz_out_str(stdout, 10, c[n]);
	// printf("\nDet(L) from diagonal product = \n");
	// mpz_t d;
	// mpz_init(d);
	// mpz_set(d, L.Index(0,0));
	// for(int i = 1; i < L.Rows; i++)
	// 	mpz_mul(d, d, L.Index(i,i));
	// mpz_out_str(stdout, 10, d);
	// printf("\n");

	// printf("\nA*B took %f seconds\n", ((float) t1)/CLOCKS_PER_SEC);
	// printf("Fadlev(A) took %f seconds\n", ((float) t2)/CLOCKS_PER_SEC);
	// printf("Det(A) took %f seconds\n", ((float) tDet)/CLOCKS_PER_SEC);
	// printf("Fadlev(L) took %f seconds\n", ((float) t3)/CLOCKS_PER_SEC);
}