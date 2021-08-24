#include <cstdio> 
#include <cctype>
#include "Matrix.h"
#include "Fadlev.h"
#include "MPList.h"
#include "TraceAlgs.h"
#include <string>
#include <gmp.h>

using namespace std;

/*
Forms a nXn matrix of traces, n being the number of matrices input,
then solves for x in Ax=y, with y the column vector of Fadlev's constant coefficients for each input matrix

Usage
./traceMatrixSolve [-f=precBits] <file1> <file2> <file3> <file4>
*/

void UsageError();

int main(int argc, char* argv[]){
	bool f = false;
	bool d = false;
	int fn = 0;
	string arg;
	Matrix<mpz_t> matrices[4](5);
	int matricesRead = 0;
	for(int i = 1; i < argc; i++){
		arg = argv[i];
		if(arg[0] == '-'){
			if(arg[1] == 'f' && arg[2] == '='){
				if(isdigit(arg[3])){
					f = true;
					fn = stoi(arg.substr(3));
				} else {
					UsageError();
					return -1;
				}
			} else if(arg[1] == 'd')
				d = true;
			else {
				UsageError();
				return -1;
			}
		} else {
			Matrix<mpz_t> A = Matrix<mpz_t>::FromFileHollowSymmetric(arg.c_str());
			matrices[matricesRead] = A;
			matricesRead++;
		}
	}
	if(matricesRead != 4){
		UsageError();
		return -1;
	}

	Matrix<mpq_t> y(4,1);
	Matrix<mpq_t> x(4,1);
	int n = matrices[0].Rows;
	MPList<mpz_t> coefficients(n + 1);
	Matrix<mpq_t> Inv(n);
	Matrix<mpz_t> A(matricesRead);
	for(int i = 0; i < matricesRead; i++){
		Fadlev(coefficients, Inv, matrices[i]);
		mpq_set_z(y.Index(i,0), coefficients[n]);
		MPList<mpz_t> terms = TraceAlgorithm3(matrices[i]);
		for(int j = 0; j < 4; j++){
			mpz_set(A.Index(i,j), terms[j]);
		}
		if(d){
			printf("Trace basis for Matrix %d:\n",i);
			terms.Print();
		}
	}
	if(d){
		printf("\nA = \n");
		A.Print();
		printf("\ny = \n");
		y.Print();
		printf("\n");
	}
	Matrix<mpq_t> Ainv(A.Rows, A.Columns);
	MPList<mpz_t> coef(5);
	Fadlev(coef, Ainv, A);
	Ainv.RightMultiply(x, y);
	printf("x = \n");
	if(f)
		x.PrintSciNot(fn);
	else
		x.Print();

	return 0;
}

void UsageError(){
	printf("Usage:\n"
		"./traceMatrixSolve [-f=precBits] [-d] <file1> <file2> <file3> <file4>\n"
		"-d: debugging mode (print intermediate info\n"
		"f=n: print results in scientific notation with n bits of mantissa precision\n");
}