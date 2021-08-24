#include <cstdio> 
#include <cctype>
#include "Matrix.h"
#include "Matrix_new.h"
#include "Fadlev.h"
#include "MPList.h"
#include "TraceAlgs.h"
#include <string>
#include <gmp.h>
#include <gmpxx.h>
#include <ctime> // for time, clock

using namespace std;

/*
Run the Faddeev-Leverrier algorithm on a square matrix

Program Options
-i = Print the inverse
-q = Use a rational matrix, as opposed to an integer matrix
-h = The input matrix shall be hollow (zeros on main diagonal)
-s = The input matrix shall be symmetric
-f=<n> = Display output as a number with a decimal, using n bits of mantissa precision
-y <file> = Parse file for a column vector y and solve for x in y = Ax
-d = debugging mode, print out each term when running trace algorithm
-r=<n> = Randomly generate an nXn matrix to use as input
<file> = file to be used as matrix input

Note: Whether the -h and -s flags are used will change the way the input file, if supplied, is interpreted
For -s, you need only provide the upper triangular portion of the matrix. 
If -h is also included, then the main diagonal may be omitted
-h does nothing if neither -s nor -r are used (i.e. the whole matrix must be specified in the input file)

Note: If both -s and -h are used (but not -q), the code will run an algorithm unique to hollow symmetric matrices in the integers. 
To choose which algorithm is run, edit the function definition of traceAlgorithm() below
*/

int GetArgAfterEquals(string& arg, string original);

int ParseNumericalOption(int& val, string arg);

void UsageError();

template <typename T>
void fadlevAlgorithm(const Matrix<T>& A, bool printInverse, bool printAsFloat, mp_bitcnt_t precBits, bool solve, string vectorFile);

void traceAlgorithm(const Matrix<mpz_t>& A, bool debugging);

int main(int argc, char* argv[]){
	bool r, inv, q, h, s, f, fd, y, yd, d, k;
	r = inv = q = h = s = f = fd = y = yd = d = k = false;
	int rn = 0;
	int fn = 0;
	int kn = 0;
	char* arg = argv[1];
	int n = 1;
	char flag;
	string fileName, vectorFile, strArg, numArg;
	while(arg != NULL){
		strArg = arg;
		if(arg[0] == '-'){
			flag = arg[1];
			switch(flag){
				case 'i': inv = true; break;
				case 'q': q = true; break;
				case 'h': h = true; break;
				case 's': s = true; break;
				case 'd': d = true; break;
				case 'f': 
					f = true; 
					GetArgAfterEquals(numArg, strArg);
					if(ParseNumericalOption(fn, numArg) < 0){
						UsageError();
						return -1;
					} 
					break;
				case 'r': 
					r = true; 
					GetArgAfterEquals(numArg, strArg);
					if(ParseNumericalOption(rn, numArg) < 0){
						UsageError();
						return -1;
					} 
					break;
				case 'k':
					k = true; 
					GetArgAfterEquals(numArg, strArg);
					if(ParseNumericalOption(kn, numArg) < 0){
						UsageError();
						return -1;
					} 
					break;
				case 'y': y = true;	break;
				default: UsageError(); return -1;
			}
		} else {
			if(y){
				vectorFile = strArg;
				yd = true;
			}
			else if(!fd){
				fileName = strArg;
				fd = true;
			} else { // two file names have been input but the y flag has not been set
				UsageError();
				return -1;
			}
		}
		n++;
		arg = argv[n];
	}
	if(y && !yd){
		fprintf(stderr, "Usage error: -y was specified but not column vector was input\n");
		return -1;
	}
	if(r + k + fd != 1){
		fprintf(stderr, "Usage error: Exactly one of -r or <input file> must be given\n");
		return -1;
	}
	if(fd){
		if(q){ // rational matrix
			if(s){
				if(h){
					Matrix<mpq_t> A = Matrix<mpq_t>::FromFileHollowSymmetric(fileName.c_str());
					fadlevAlgorithm(A, inv, f, fn, y, vectorFile);
				} else {
					Matrix<mpq_t> A = Matrix<mpq_t>::FromFileSymmetric(fileName.c_str());
					fadlevAlgorithm(A, inv, f, fn, y, vectorFile);
				}
			} else {
				Matrix<mpq_t> A = Matrix<mpq_t>::FromFile(fileName.c_str());
				fadlevAlgorithm(A, inv, f, fn, y, vectorFile);
			}
		} else { // integer matrix
			if(s){
				if(h){
					Matrix<mpz_t> A = Matrix<mpz_t>::FromFileHollowSymmetric(fileName.c_str());
					fadlevAlgorithm(A, inv, f, fn, y, vectorFile);
					traceAlgorithm(A, d);
				} else {
					Matrix<mpz_t> A = Matrix<mpz_t>::FromFileSymmetric(fileName.c_str());
					fadlevAlgorithm(A, inv, f, fn, y, vectorFile);
				}
			} else {
				Matrix<mpz_t> A = Matrix<mpz_t>::FromFile(fileName.c_str());
				fadlevAlgorithm(A, inv, f, fn, y, vectorFile);
			}
		}
	} else {
		if (r) { // Make random matrix
			mpz_t maxRnd, rnd;
			mpz_inits(maxRnd, rnd, NULL);
			mpz_set_ui(maxRnd, 10);
			gmp_randstate_t r;
			gmp_randinit_default(r);
			gmp_randseed_ui(r, time(NULL));

			if(q){
				Matrix<mpq_t> A = Matrix<mpq_t>::Random((int) rn, r, maxRnd);
				if(h){
					for(int r = 0; r < A.Rows; r++){
						mpq_set_ui(A.Index(r,r), 0, 1);
					}
				}
				if(s){
					for(int r = 0; r < A.Rows; r++){
						for(int c = r + 1; c < A.Columns; c++){
							mpq_set(A.Index(c,r), A.Index(r,c));
						}
					}
				}
				printf("A = \n");
				A.Print();
				printf("\n");
				fadlevAlgorithm(A, inv, f, fn, y, vectorFile);
			} else {
				Matrix<mpz_t> A = Matrix<mpz_t>::Random(rn, r, maxRnd);
				if(h){
					for(int r = 0; r < A.Rows; r++){
						mpz_set_ui(A.Index(r,r), 0);
					}
				}
				if(s){
					for(int r = 0; r < A.Rows; r++){
						for(int c = r + 1; c < A.Columns; c++){
							mpz_set(A.Index(c,r), A.Index(r,c));
						}
					}
				}
				fadlevAlgorithm(A, inv, f, fn, y, vectorFile);
				if(s && h)
					traceAlgorithm(A, d);
			}
		} else if (k) {
			// Matrix_new<mpz_class> B = Matrix_new<mpz_class>::FromCompleteGraph(kn);
			// Matrix_new<mpz_class> B2(kn);
			// Matrix_new<mpz_class> B3(kn);
			// clock_t t1 = clock();
			// B.RightMultiply(B2,B);
			// B2.RightMultiply(B3,B);
			// mpz_class output = B3.Trace()/6;
			// t1 = clock() - t1;
			// cout << "Number of 3-cycles in K" << kn << " - " << output;
			// printf("\nIt took %f seconds\n", ((float) t1)/CLOCKS_PER_SEC);

			Matrix<mpz_t> A = Matrix<mpz_t>::FromCompleteGraph(kn);
			Matrix<mpz_t> A2(kn);
			Matrix<mpz_t> A3(kn);
			mpz_t t;
			mpz_init(t);
			mpq_t output, sixth, tq;
			mpq_inits(output, sixth, tq, NULL); 
			mpq_set_ui(sixth, 1, 6);
			
			clock_t t1 = clock();
			A.RightMultiply(A2, A);
			A2.RightMultiply(A3, A);
			A3.Trace(t);
			mpq_set_z(tq, t);
			mpq_mul(output, tq, sixth);
			t1 = clock() - t1;

			printf("Number of 3-cycles in K%d - ", kn);
			mpq_out_str(stdout, 10, output);
			printf("\nIt took %f seconds\n", ((float) t1)/CLOCKS_PER_SEC);
		}
	}

	return 0;
}

int GetArgAfterEquals(string& arg, string original){
	unsigned int start;
	if((start = original.find_first_of('=')) == string::npos)
		return -1;
	start++;
	arg = original.substr(start);
	return 0;
}

int ParseNumericalOption(int& val, string arg){
	if(!isdigit(arg[0]))
		return -1;
	val = atoi(arg.c_str());
	return 0;
}

void UsageError(){
	fprintf(stderr, "Usage\n"
		"-i: show the inverse of the matrix at no extra cost\n"
		"-q: interpret the matrix in the rational numbers\n"
		"-h: use a hollow matrix (zeros on main diagonal)\n"
		"-s: use a symmetric matrix\n"
		"-d: debugging mode, print out each term when running trace algorithm\n"
		"-f=<n>: output numbers in exponential notation with n bits of mantissa precision\n"
		"-y <file>: Parse file for a column vector y and solve for x in y = Ax\n"
		"(-r=<n>: use random nXn matrix] | <input matrix file>)\n");
}

template <typename T>
void fadlevAlgorithm(const Matrix<T>& A, bool printInverse, bool printAsFloat, mp_bitcnt_t precBits, bool solve, string vectorFile){
	printf("A = \n");
	A.Print();
	Matrix<mpq_t> Ainv(A.Rows, A.Columns);
	MPList<T> c(A.Rows + 1);

	Fadlev(c, Ainv, A);

	if(printInverse){
		printf("\nInv(A) = \n");
		if(printAsFloat)
			Ainv.PrintSciNot(precBits);
		else
			Ainv.Print();
	}
	
	printf("\nFadlev(A) = \n");
	if(printAsFloat)
		c.PrintSciNot(precBits);
	else
		c.Print();

	if(solve){
		Matrix<mpq_t> y = Matrix<mpq_t>::FromFile(vectorFile.c_str());
		Matrix<mpq_t> x(A.Columns, y.Columns);
		printf("\nSolve for x in Ax = y: \n");
		Ainv.RightMultiply(x, y);
		if(printAsFloat)
			x.PrintSciNot(precBits);
		else
			x.Print();
	}
}

void traceAlgorithm(const Matrix<mpz_t>& A, bool debugging){
	mpz_t p;
	mpz_init(p);
	printf("\nTrace Algorithm\n");
	MPList<mpz_t> terms = TraceAlgorithm2(p, A);
	if(debugging)
		terms.Print("\n");
	printf("= ");
	mpz_out_str(stdout, 10, p);
	printf("\n");
}
