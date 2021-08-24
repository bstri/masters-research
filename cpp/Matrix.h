#ifndef MATRIX_H
#define MATRIX_H

#include <cstdio>
#include <cstdlib>
#include <gmp.h>
#include <fstream>
#include <string>
#include "mpImpl.h"

using namespace std;

/* gmp functions tend not to return things, but instead use the first parameter as a return value,
 so I've adopted that same style in the Matrix class */

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
class Matrix{
public:
	int Rows;
	int Columns;
	T* Data;

	// Constructor
	Matrix(int rows, int columns){
		init(rows, columns);
	}
	// Square (sizeXsize) matrix constructor
	Matrix(int size){
		init(size, size);
	}
	// Copy constructor
	Matrix(const Matrix<T>& M){
		init(M.Rows, M.Columns);
		Copy(M);
	}

	// Construct the sizeXsize Identity matrix
	static Matrix<T> Identity(int size){
		Matrix<T> I(size, size);
		for(int r = 0; r < size; r++){
			for(int c = 0; c < size; c++)
				mpImpl::set_ui(I.Index(r,c), 1);
		}
		return I;
	}

	// Constructs a lower triangular matrix from M
	static Matrix<T> Lower(const Matrix<T>& M){
		Matrix<T> L(M.Rows, M.Columns);
		M.pasteLower(L);
		return L;
	}

	// Make L a lower triangular matrix from this
	void GetLower(Matrix<T>& L) const {
		pasteLower(L);
		for(int r = 0; r < L.Rows; r++){
			for(int c = r + 1; c < L.Columns; c++)
				mpImpl::set_ui(L.Index(r,c), 0);
		}
	}

	// Constructs an upper triangular matrix from M
	static Matrix<T> Upper(const Matrix<T>& M){
		Matrix<T> U(M.Rows, M.Columns);
		M.pasteUpper(U);
		return U;
	}

	// Make U an upper triangular matrix from this
	void GetUpper(Matrix<T>& U) const {
		pasteUpper(U);
		for(int r = 0; r < U.Rows; r++){
			for(int c = 0; c < r; c++)
				mpImpl::set_ui(U.Index(r,c), 0);
		}
	}

	// Constructs a diagonal matrix from M
	static Matrix<T> Diagonal(const Matrix<T>& M){
		Matrix<T> D(M.Rows, M.Columns);
		M.pasteDiagonal(D);
		return D;
	}

	// Make D a diagonal matrix from this
	void GetDiagonal(Matrix<T>& D) const {
		pasteDiagonal(D);
		for(int r = 0; r < D.Rows; r++){
			for(int c = 0; c < r; c++)
				mpImpl::set_ui(D.Index(r,c), 0);
			for(int c = r + 1; c < D.Columns; c++)
				mpImpl::set_ui(D.Index(r,c), 0);
		}
	}

	// Constructs a matrix with random integer entries between 1 and maxRnd
	static Matrix<T> Random(int rows, int columns, gmp_randstate_t rndState, const mpz_t maxRnd){
		Matrix<T> R(rows, columns);
		R.randomize(rndState, maxRnd);
		return R;
	}
	// Square random matrix
	static Matrix<T> Random(int size, gmp_randstate_t rndState, const mpz_t maxRnd){
		Matrix<T> R(size);
		R.randomize(rndState, maxRnd);
		return R;
	}

	// Construct matrix from text file
	// First line should have the format: <num rows> <num columns>
	// Spaces delimit entries, newlines delimit rows
	// Note: rational matrix input does not support decimal points, only fractions of the form "a/b"
	static Matrix<T> FromFile(const char* fileName){
		ifstream inFile(fileName);
		string line;
		getline(inFile, line);

		size_t start = 0;
		string r = getNextTokenFromFile(inFile, start, line);
		string c = getNextTokenFromFile(inFile, start, line);
		Matrix<T> F(stoi(r), stoi(c));

		int rowsRead = 0;
		int column = 0;
		string token;
		while(inFile){
			token = getNextTokenFromFile(inFile, start, line);
			mpImpl::set_str(F.Index(rowsRead, column), token.c_str(), 10);
			mpImpl::canonicalize(F.Index(rowsRead, column));
			column++;
			if(start == 0){
				rowsRead++;
				column = 0;
			}
		}

		inFile.close();
		return F;
	}

	// Constructs symmetric matrix from file, with only the upper triangular portion of the matrix specified
	static Matrix<T> FromFileSymmetric(const char* fileName){
		ifstream inFile(fileName);
		string line;
		getline(inFile, line);

		size_t start = 0;
		string r = getNextTokenFromFile(inFile, start, line);
		string c = getNextTokenFromFile(inFile, start, line);
		Matrix<T> F(stoi(r), stoi(c));

		int rowsRead = 0;
		int column = 0;
		string token;
		while(inFile){
			token = getNextTokenFromFile(inFile, start, line);
			mpImpl::set_str(F.Index(rowsRead, column), token.c_str(), 10);
			mpImpl::canonicalize(F.Index(rowsRead, column));
			if(rowsRead != column){
				mpImpl::set_str(F.Index(column, rowsRead), token.c_str(), 10);
				mpImpl::canonicalize(F.Index(column, rowsRead));
			}
			column++;
			if(start == 0){
				rowsRead++;
				column = rowsRead;
			}
		}

		inFile.close();
		return F;
	}

	// Constructs hollow symmetric matrix from file, with only the upper triangular portion specified, not including the main diagonal
	static Matrix<T> FromFileHollowSymmetric(const char* fileName){
		ifstream inFile(fileName);
		string line;
		getline(inFile, line);

		size_t start = 0;
		string r = getNextTokenFromFile(inFile, start, line);
		string c = getNextTokenFromFile(inFile, start, line);
		Matrix<T> F(stoi(r), stoi(c));

		int rowsRead = 0;
		int column = 1;
		string token;
		while(inFile){
			token = getNextTokenFromFile(inFile, start, line);
			mpImpl::set_str(F.Index(rowsRead, column), token.c_str(), 10);
			mpImpl::canonicalize(F.Index(rowsRead, column));
			if(rowsRead != column){
				mpImpl::set_str(F.Index(column, rowsRead), token.c_str(), 10);
				mpImpl::canonicalize(F.Index(column, rowsRead));
			}
			column++;
			if(start == 0){
				rowsRead++;
				column = rowsRead + 1;
			}
		}

		inFile.close();
		return F;
	}

	static Matrix<T> FromCompleteGraph(int n){
		Matrix<T> K(n);
		for (int c = 0; c < n; c++){
			for (int r = 0; r < n; r++){
				if (r == c) {
					mpz_set_ui(K.Index(r,c), 0);
				} else {
					mpz_set_ui(K.Index(r,c), 1);
				}
			}
		}
		return K;
	}

	static Matrix<mpq_t> ToRationalMatrix(const Matrix<mpz_t>& src) {
		Matrix<mpq_t> Q(src.Rows, src.Columns);
		for (int r = 0; r < Q.Rows; r++) {
			for (int c = 0; c < Q.Columns; c++) {
				mpImpl::to_mpq(Q.Index(r,c), src.Index(r,c));
			}
		}
		return Q;
	}
	static Matrix<mpq_t> ToRationalMatrix(const Matrix<mpq_t>& src) {
		Matrix<mpq_t> Q(src.Rows, src.Columns);
		for (int r = 0; r < Q.Rows; r++) {
			for (int c = 0; c < Q.Columns; c++) {
				mpImpl::to_mpq(Q.Index(r,c), src.Index(r,c));
			}
		}
		return Q;
	}

	// Overload assignment operator
	Matrix<T>& operator=(const Matrix<T>& M) {
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
	~Matrix(){
		Deallocate();
	}

	// Free up memory used by this matrix
	void Deallocate(){
		for(int i = 0; i < Rows*Columns; i++)
			mpImpl::clear(Data[i]);
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
				mpImpl::mul(Index(r,c), Index(r,c), k);
		}
	}

	// Multiplies this and B and stores the result in M
	void RightMultiply(Matrix<T>& M, const Matrix<T>& B) const {
		T sum, prod;
		mpImpl::init(sum);
		mpImpl::init(prod);
		for(int c = 0; c < B.Columns; c++){
			for(int r = 0; r < Rows; r++){
				mpImpl::set_ui(sum, 0);
				for(int k = 0; k < Columns; k++){
					mpImpl::mul(prod, B.Index(k,c), Index(r,k));
					mpImpl::add(sum, sum, prod);
				}
				mpImpl::set(M.Index(r,c), sum);
			}
		}
	}

	// AKA Hadamard product
	// Multiplies this and B entrywise and stores the result in M
	void MultiplyEntrywise(Matrix<T>& M, const Matrix<T>& B) const {
		T prod;
		mpImpl::init(prod);
		for(int r = 0; r < Rows; r++){
			for(int c = 0; c < Columns; c++){
				mpImpl::mul(prod, Index(r,c), B.Index(r,c));
				mpImpl::set(M.Index(r,c), prod);
			}
		}
	}

	void Determinant(mpq_t& det) const {
		Matrix<mpq_t> Q = this->ToRationalMatrix(*this);
		mpq_set_ui(det, 1, 1);
		mpq_t zero, coefficient, sub;
		mpq_inits(zero, coefficient, sub, NULL);
		for (int r = 0; r < Q.Rows - 1; r++) {
			// find first non-zero pivot and swap row with current row
			// for (int r2 = r; r2 < Q.Rows; r2++) {
			// 	if (mpq_cmp(Q.Index(r2, r), zero) != 0) {
			// 		if (r2 == r)
			// 			break;
			// 		// swap row r2 with row r
			// 		for (int c = r; c < Q.Columns; c++) {
			// 			mpq_swap(Q.Index(r2, c), Q.Index(r, c));
			// 		}
			// 	}
			// }

			mpq_mul(det, det, Q.Index(r, r));
			if (mpq_cmp(det, zero) != 0)
				return; // determinant is 0

			// perform row reduction
			for (int r2 = r + 1; r2 < Q.Rows; r2++) {
				mpq_div(coefficient, Q.Index(r2, r), Q.Index(r, r));
				for (int c = r; c < Q.Columns; c++) {
					mpq_mul(sub, coefficient, Q.Index(r, c));
					mpq_sub(Q.Index(r2, c), Q.Index(r2, c), sub);
				}
			}
		}
		mpq_mul(det, det, Q.Index(Q.Rows - 1, Q.Columns - 1));
	}

	// Prints the matrix to stdout
	void Print() const {
		for(int r = 0; r < Rows; r++){
			for(int c = 0; c < Columns; c++){
				mpImpl::out_str(stdout, 10, Index(r,c));
				printf(" ");
			}
			printf("\n");
		}
	}

	// Prints the matrix in scientific notation
	void PrintSciNot(mp_bitcnt_t precBits) const {
		mpf_t f;
		mpf_init2(f, precBits);
		for(int r = 0; r < Rows; r++){
			for(int c = 0; c < Columns; c++){
				mpImpl::to_float(f, Index(r,c));
				mpf_out_str(stdout, 10, 0, f);
				printf(" ");
			}
			printf("\n");
		}
	}

	// Prints the r,c element to stdout
	// Note: compilation with '-O0' is necessary to turn off optimization
	void PrintElement(int r, int c) const {
		mpImpl::out_str(stdout, 10, Index(r,c));
		printf("\n");
	}

	// Computes the trace of this and stores it in trace
	void Trace(T trace) const {
		mpImpl::set(trace, Index(0,0));
		for(int r = 1; r < Rows; r++)
			mpImpl::add(trace, trace, Index(r,r));
	}
	
	void Copy(const Matrix<mpz_t>& M){
		T converted;
		mpImpl::init(converted);
		for(int r = 0; r < Rows; r++){
			for(int c = 0; c < Columns; c++){
				mpImpl::from_mpz(converted, M.Index(r,c));
				mpImpl::set(Index(r,c), converted);
			}
		}
	}
	void Copy(const Matrix<mpq_t>& M){
		T converted;
		mpImpl::init(converted);
		for(int r = 0; r < Rows; r++){
			for(int c = 0; c < Columns; c++){
				mpImpl::from_mpq(converted, M.Index(r,c));
				mpImpl::set(Index(r,c), converted);
			}
		}
	}
	
private:
	// Allocate memory and initialize mp vars
	void init(int rows, int columns){
		Rows = rows;
		Columns = columns;
		Data = (T*) malloc(rows*columns*sizeof(T));
		for(int r = 0; r < rows; r++){
			for(int c = 0; c < columns; c++)
				mpImpl::init(Index(r,c));
		}
	}

	// Copy the lower triangular of this to L, not overwriting other entries
	void pasteLower(Matrix<T>& L) const {
		for(int r = 0; r < L.Rows; r++){
			for(int c = 0; c <= r; c++)
				mpImpl::set(L.Index(r,c), Index(r,c));
		}
	}

	// Copy the upper triangular of this to U, not overwriting other entries
	void pasteUpper(Matrix<T>& U) const {
		for(int r = 0; r < U.Rows; r++){
			for(int c = r; c < U.Columns; c++)
				mpImpl::set(U.Index(r,c), Index(r,c));
		}
	}

	// Copy the main diagonal of this to D, not overwriting other entries
	void pasteDiagonal(Matrix<T>& D) const {
		for(int r = 0; r < D.Rows; r++){
			mpImpl::set(D.Index(r,r), Index(r,r));
		}
	}

	// Randomize entries to between 1 and maxRnd
	void randomize(gmp_randstate_t rndState, const mpz_t maxRnd){
		mpz_t rnd;
		mpImpl::init(rnd);
		T conv;
		mpImpl::init(conv);
		for(int r = 0; r < Rows; r++){
			for(int c = 0; c < Columns; c++){
				mpz_urandomm(rnd, rndState, maxRnd);
				mpz_add_ui(rnd, rnd, 1);
				mpImpl::from_mpz(conv, rnd);
				mpImpl::set(Index(r,c), conv);
			}
		}
	}

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