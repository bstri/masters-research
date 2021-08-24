#ifndef MATRIX_NEW_H
#define MATRIX_NEW_H

#include <cstdio>
#include <cstdlib>
#include <gmp.h>
#include <gmpxx.h>
#include <fstream>
#include <string>
#include "mpImpl.h"
#include <iostream>

using namespace std;



/* NOTE:
A matrix such as
|1 2|
|3 4|
will be stored as an array like
{1, 2, 3, 4}
*/

// return successive space-delimited and newline-delimited entries in the file stream
// string getNextTokenFromFile(ifstream& inFile, size_t& start, string& line);

template <typename T>
class Matrix_new{
public:
	int Rows;
	int Columns;
	T* Data;

	// Constructor
	Matrix_new(int rows, int columns){
		init(rows, columns);
	}
	// Square (sizeXsize) matrix constructor
	Matrix_new(int size){
		init(size, size);
	}
	// Copy constructor
	// Matrix_new(const Matrix_new<T>& M){
	// 	init(M.Rows, M.Columns);
	// 	Copy(M);
	// }

	// Construct the sizeXsize Identity matrix
	static Matrix_new<T> Identity(int size){
		Matrix_new<T> I(size);
		for(int r = 0; r < size; r++){
			for(int c = 0; c < size; c++)
				I.Index(r,c) = 1;
		}
		return I;
	}

	// Constructs a lower triangular matrix from M
	static Matrix_new<T> Lower(const Matrix_new<T>& M){
		Matrix_new<T> L(M.Rows, M.Columns);
		M.pasteLower(L);
		return L;
	}

	// Make L a lower triangular matrix from this
	void GetLower(Matrix_new<T>& L) const {
		pasteLower(L);
		for(int r = 0; r < L.Rows; r++){
			for(int c = r + 1; c < L.Columns; c++)
                L.Index(r,c) = 0;
		}
	}

	// Constructs an upper triangular matrix from M
	static Matrix_new<T> Upper(const Matrix_new<T>& M){
		Matrix_new<T> U(M.Rows, M.Columns);
		M.pasteUpper(U);
		return U;
	}

	// Make U an upper triangular matrix from this
	void GetUpper(Matrix_new<T>& U) const {
		pasteUpper(U);
		for(int r = 0; r < U.Rows; r++){
			for(int c = 0; c < r; c++)
                U.Index(r,c) = 0;
		}
	}

	// Constructs a diagonal matrix from M
	static Matrix_new<T> Diagonal(const Matrix_new<T>& M){
		Matrix_new<T> D(M.Rows, M.Columns);
		M.pasteDiagonal(D);
		return D;
	}

	// Make D a diagonal matrix from this
	void GetDiagonal(Matrix_new<T>& D) const {
		pasteDiagonal(D);
		for(int r = 0; r < D.Rows; r++){
			for(int c = 0; c < r; c++)
                D.Index(r,c) = 0;
			for(int c = r + 1; c < D.Columns; c++)
                D.Index(r,c) = 0;
		}
	}

	// Constructs a matrix with random integer entries between 1 and maxRnd
	// static Matrix_new<T> Random(int rows, int columns, gmp_randstate_t rndState, const mpz_t maxRnd){
	// 	Matrix_new<T> R(rows, columns);
	// 	R.randomize(rndState, maxRnd);
	// 	return R;
	// }
	// // Square random matrix
	// static Matrix_new<T> Random(int size, gmp_randstate_t rndState, const mpz_t maxRnd){
	// 	Matrix_new<T> R(size);
	// 	R.randomize(rndState, maxRnd);
	// 	return R;
	// }

	// Construct matrix from text file
	// First line should have the format: <num rows> <num columns>
	// Spaces delimit entries, newlines delimit rows
	// Note: rational matrix input does not support decimal points, only fractions of the form "a/b"
	// static Matrix_new<T> FromFile(const char* fileName){
	// 	ifstream inFile(fileName);
	// 	string line;
	// 	getline(inFile, line);

	// 	size_t start = 0;
	// 	string r = getNextTokenFromFile(inFile, start, line);
	// 	string c = getNextTokenFromFile(inFile, start, line);
	// 	Matrix_new<T> F(stoi(r), stoi(c));

	// 	int rowsRead = 0;
	// 	int column = 0;
	// 	string token;
	// 	while(inFile){
	// 		token = getNextTokenFromFile(inFile, start, line);
	// 		mpImpl::set_str(F.Index(rowsRead, column), token.c_str(), 10);
	// 		mpImpl::canonicalize(F.Index(rowsRead, column));
	// 		column++;
	// 		if(start == 0){
	// 			rowsRead++;
	// 			column = 0;
	// 		}
	// 	}

	// 	inFile.close();
	// 	return F;
	// }

	// Constructs symmetric matrix from file, with only the upper triangular portion of the matrix specified
	// static Matrix_new<T> FromFileSymmetric(const char* fileName){
	// 	ifstream inFile(fileName);
	// 	string line;
	// 	getline(inFile, line);

	// 	size_t start = 0;
	// 	string r = getNextTokenFromFile(inFile, start, line);
	// 	string c = getNextTokenFromFile(inFile, start, line);
	// 	Matrix_new<T> F(stoi(r), stoi(c));

	// 	int rowsRead = 0;
	// 	int column = 0;
	// 	string token;
	// 	while(inFile){
	// 		token = getNextTokenFromFile(inFile, start, line);
	// 		mpImpl::set_str(F.Index(rowsRead, column), token.c_str(), 10);
	// 		mpImpl::canonicalize(F.Index(rowsRead, column));
	// 		if(rowsRead != column){
	// 			mpImpl::set_str(F.Index(column, rowsRead), token.c_str(), 10);
	// 			mpImpl::canonicalize(F.Index(column, rowsRead));
	// 		}
	// 		column++;
	// 		if(start == 0){
	// 			rowsRead++;
	// 			column = rowsRead;
	// 		}
	// 	}

	// 	inFile.close();
	// 	return F;
	// }

	// Constructs hollow symmetric matrix from file, with only the upper triangular portion specified, not including the main diagonal
	// static Matrix_new<T> FromFileHollowSymmetric(const char* fileName){
	// 	ifstream inFile(fileName);
	// 	string line;
	// 	getline(inFile, line);

	// 	size_t start = 0;
	// 	string r = getNextTokenFromFile(inFile, start, line);
	// 	string c = getNextTokenFromFile(inFile, start, line);
	// 	Matrix_new<T> F(stoi(r), stoi(c));

	// 	int rowsRead = 0;
	// 	int column = 1;
	// 	string token;
	// 	while(inFile){
	// 		token = getNextTokenFromFile(inFile, start, line);
	// 		mpImpl::set_str(F.Index(rowsRead, column), token.c_str(), 10);
	// 		mpImpl::canonicalize(F.Index(rowsRead, column));
	// 		if(rowsRead != column){
	// 			mpImpl::set_str(F.Index(column, rowsRead), token.c_str(), 10);
	// 			mpImpl::canonicalize(F.Index(column, rowsRead));
	// 		}
	// 		column++;
	// 		if(start == 0){
	// 			rowsRead++;
	// 			column = rowsRead + 1;
	// 		}
	// 	}

	// 	inFile.close();
	// 	return F;
	// }

	static Matrix_new<T> FromCompleteGraph(int n){
		Matrix_new<T> K(n);
		for (int c = 0; c < n; c++){
			for (int r = 0; r < n; r++){
				if (r == c) {
                    K.Index(r,c) = 0;
				} else {
                    K.Index(r,c) = 1;
				}
			}
		}
		return K;
	}

	// Overload assignment operator
	Matrix_new<T>& operator=(const Matrix_new<T>& M) {
		if(this == &M)
			return *this;
		if(Rows != M.Rows || Columns != M.Columns){
			Deallocate();
			init(M.Rows, M.Columns);
		}
		Copy(M);
		return *this;
	}

	// Destructor
	~Matrix_new(){
		Deallocate();
	}

	// Free up memory used by this matrix
	void Deallocate(){
		// for(int i = 0; i < Rows*Columns; i++)
		// 	mpz_clear(Data[i]);
		free(Data);
		Data = NULL;
	}

	// Returns the row,column element
	T& Index(int row, int column) const {
		return Data[row*Columns + column];
	}

	// Multiplies all entries by the scalar
	void Scale(const T k){
		for(int r = 0; r < Rows; r++){
			for(int c = 0; c < Columns; c++)
                Index(r,c) *= k;
		}
	}

	// Multiplies this and B and stores the result in M
	void RightMultiply(Matrix_new<T>& M, const Matrix_new<T>& B) const {
        mpz_class sum;
		for(int c = 0; c < B.Columns; c++){
			for(int r = 0; r < Rows; r++){
                sum = 0;
				for(int k = 0; k < Columns; k++){
                    sum += B.Index(k,c)*Index(r,k);
				}
                M.Index(r,c) = sum;
			}
		}
	}

	// AKA Hadamard product
	// Multiplies this and B entrywise and stores the result in M
	void MultiplyEntrywise(Matrix_new<T>& M, const Matrix_new<T>& B) const {
		// T prod;
		// mpImpl::init(prod);
		for(int r = 0; r < Rows; r++){
			for(int c = 0; c < Columns; c++){
                M.Index(r,c) = Index(r,c) * B.Index(r,c);
				// mpImpl::mul(prod, Index(r,c), B.Index(r,c));
				// mpImpl::set(M.Index(r,c), prod);
			}
		}
	}

	// Prints the matrix to stdout
	void Print() const {
		for(int r = 0; r < Rows; r++){
			for(int c = 0; c < Columns; c++){
                cout << Index(r,c) << " ";
			}
			printf("\n");
		}
	}

	// Prints the matrix in scientific notation
	// void PrintSciNot(mp_bitcnt_t precBits) const {
	// 	mpf_t f;
	// 	mpf_init2(f, precBits);
	// 	for(int r = 0; r < Rows; r++){
	// 		for(int c = 0; c < Columns; c++){
	// 			mpImpl::to_float(f, Index(r,c));
	// 			mpf_out_str(stdout, 10, 0, f);
	// 			printf(" ");
	// 		}
	// 		printf("\n");
	// 	}
	// }

	// Prints the r,c element to stdout
	// Note: compilation with '-O0' is necessary to turn off optimization
	// void PrintElement(int r, int c) const {
	// 	cout << Index(r,c) << "\n";
	// }

	// Computes the trace of this and stores it in trace
	T Trace() const {
		mpz_class trace = Index(0,0);
		for(int r = 1; r < Rows; r++)
			trace += Index(r,r);
        return trace;
	}
	
	// void Copy(const Matrix_new<mpz_t>& M){
	// 	T converted;
	// 	mpImpl::init(converted);
	// 	for(int r = 0; r < Rows; r++){
	// 		for(int c = 0; c < Columns; c++){
	// 			mpImpl::from_mpz(converted, M.Index(r,c));
	// 			mpImpl::set(Index(r,c), converted);
	// 		}
	// 	}
	// }
	// void Copy(const Matrix_new<mpq_t>& M){
	// 	T converted;
	// 	mpImpl::init(converted);
	// 	for(int r = 0; r < Rows; r++){
	// 		for(int c = 0; c < Columns; c++){
	// 			mpImpl::from_mpq(converted, M.Index(r,c));
	// 			mpImpl::set(Index(r,c), converted);
	// 		}
	// 	}
	// }
	
private:
	// Allocate memory and initialize mp vars
	void init(int rows, int columns){
		Rows = rows;
		Columns = columns;
		Data = (T*) malloc(rows*columns*sizeof(T));
	}

	// Copy the lower triangular of this to L, not overwriting other entries
	void pasteLower(Matrix_new<T>& L) const {
		for(int r = 0; r < L.Rows; r++){
			for(int c = 0; c <= r; c++)
				L.Index(r,c) = Index(r,c);
		}
	}

	// Copy the upper triangular of this to U, not overwriting other entries
	void pasteUpper(Matrix_new<T>& U) const {
		for(int r = 0; r < U.Rows; r++){
			for(int c = r; c < U.Columns; c++)
				U.Index(r,c) = Index(r,c);
		}
	}

	// Copy the main diagonal of this to D, not overwriting other entries
	void pasteDiagonal(Matrix_new<T>& D) const {
		for(int r = 0; r < D.Rows; r++){
			D.Index(r,r) = Index(r,r);
		}
	}

	// Randomize entries to between 1 and maxRnd
	// void randomize(gmp_randstate_t rndState, const mpz_t maxRnd){
	// 	mpz_t rnd;
	// 	mpImpl::init(rnd);
	// 	T conv;
	// 	mpImpl::init(conv);
	// 	for(int r = 0; r < Rows; r++){
	// 		for(int c = 0; c < Columns; c++){
	// 			mpz_urandomm(rnd, rndState, maxRnd);
	// 			mpz_add_ui(rnd, rnd, 1);
	// 			mpImpl::from_mpz(conv, rnd);
	// 			mpImpl::set(Index(r,c), conv);
	// 		}
	// 	}
	// }

	// return successive space-delimited and newline-delimited entries in the file stream
	static string getNextTokenFromFile(ifstream& inFile, size_t& start, string& line){
		string sub;
		string ws = " ";
		if(start == 0){
			while(line.length() == 0 || (start = line.find_first_not_of(ws)) == string::npos)
				getline(inFile, line);
			line = line.substr(0, line.find_last_not_of(ws) + 1);
		}
		size_t end = line.find_first_of(ws, start);
		sub = line.substr(start, end - start);
		if(end == string::npos){
			getline(inFile, line);
			start = 0;
		} else
			start = line.find_first_not_of(ws, end);
		return sub;
	}
};


#endif