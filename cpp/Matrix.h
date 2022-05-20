#ifndef MATRIX_H
#define MATRIX_H

#include <cstdio>
#include <cstdlib>
#include <gmp.h>
#include <fstream>
#include <string>
#include "mpImpl.h"
#include "MPList.h"
#include <vector>

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

	// Weights the non-zero entries of the matrix, keeping symmetry and hollow-ness
	Matrix<T> RandomlyWeight(gmp_randstate_t rndState, unsigned int maxRnd) const {
		Matrix<T> weightedCopy = *this;
		mpz_t rand, maxRand_z, sign, two;
		mpz_inits(rand, maxRand_z, sign, two, NULL);
		mpz_set_ui(maxRand_z, maxRnd);
		mpz_set_ui(two, 2);
		for (int i = 1; i < Rows; i++) {
			for (int j = 0; j < i; j++) {
				mpz_urandomm(sign, rndState, two);
				mpz_mul_ui(sign, sign, 2);
				mpz_sub_ui(sign, sign, 1);
				mpz_urandomm(rand, rndState, maxRand_z);
				mpz_add_ui(rand, rand, 1);
				mpz_mul(rand, rand, sign);
				mpImpl::mul(weightedCopy.Index(i,j), weightedCopy.Index(i,j), rand);
				mpImpl::set(weightedCopy.Index(j,i), weightedCopy.Index(i,j));
			}
		}
		return weightedCopy;
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

	// flatten: if true, makes all nonzero elements 1
	static Matrix<T> FromLineHollowSymmetric(string line, bool flatten=false){
		vector<string> elements;
		string token;
		size_t start = 0;
		int row = 1;
		int col = 0;
		T element;
		mpImpl::init(element);
		do {
			token = getNextTokenFromLine(start, line);
			if (flatten) {
				if (token != "0") {
					token = "1";
				}
			}
			elements.push_back(token);
			col++;
			if (col == row) {
				col = 0;
				row++;
			}
		} while (start != string::npos);
		Matrix<T> M(row);
		row = 1;
		col = 0;
		mpImpl::set_ui(M.Index(0,0), 0);
		for (size_t i = 0; i < elements.size(); i++) {
			mpImpl::set_str(element, elements[i].c_str(), 10);
			mpImpl::canonicalize(element);
			mpImpl::set(M.Index(row,col), element);
			mpImpl::set(M.Index(col,row), element);
			col++;
			if (row == col) {
				col = 0;
				row++;
			}
		}
		return M;
	}

	static Matrix<T> FromCompleteGraph(int n){
		Matrix<T> K(n);
		for (int c = 0; c < n; c++){
			for (int r = 0; r < n; r++){
				if (r == c) {
					mpImpl::set_ui(K.Index(r,c), 0);
				} else {
					mpImpl::set_ui(K.Index(r,c), 1);
				}
			}
		}
		return K;
	}

	static Matrix<T> FromCycle(int n){
		Matrix<T> C(n);
		for (int i = 1; i < n; i++) {
			mpImpl::set_ui(C.Index(i, i-1), 1);
			mpImpl::set_ui(C.Index(i-1, i), 1);
		}
		mpImpl::set_ui(C.Index(0, n-1), 1);
		mpImpl::set_ui(C.Index(n-1, 0), 1);
		return C;
	}

	Matrix<T> Combine(const Matrix<T> other, bool disconnected) const {
		int offset = disconnected ? 0 : 1;
		Matrix<T> combined(Rows + other.Rows - offset, Columns + other.Columns - offset);
		for (int r = 0; r < Rows; r++){
			for (int c = 0; c < Columns; c++) {
				mpImpl::set(combined.Index(r,c), Index(r,c));
			}
		}
		for (int r = Rows - offset; r < Rows + other.Rows - offset; r++) {
			for (int c = Columns - offset; c < Columns + other.Columns - offset; c++) {
				mpImpl::set(combined.Index(r,c), other.Index(r - Rows + offset, c - Columns + offset));
			}
		}
		return combined;
	}

	Matrix<T> AppendRow(MPList<T> row) const {
		Matrix<T> newMatrix = Matrix<T>(Rows + 1, Columns);
		for (int r = 0; r < Rows; r++){
			for (int c = 0; c < Columns; c++){
				mpImpl::set(newMatrix.Index(r,c), Index(r,c));
			}
		}
		for (int i = 0; i < Columns; i++){
			mpImpl::set(newMatrix.Index(Rows,i), row[i]);
		}
		return newMatrix;
	}

	void ToRationalMatrix(Matrix<mpz_t>& out) const {
		for (int r = 0; r < out.Rows; r++) {
			for (int c = 0; c < out.Columns; c++) {
				mpImpl::to_mpq(out.Index(r,c), this->Index(r,c));
			}
		}
	}
	void ToRationalMatrix(Matrix<mpq_t>& out) const {
		for (int r = 0; r < out.Rows; r++) {
			for (int c = 0; c < out.Columns; c++) {
				mpImpl::to_mpq(out.Index(r,c), this->Index(r,c));
			}
		}
	}
	void ToRationalMatrix(Matrix<float>& out) const {
		for (int r = 0; r < out.Rows; r++) {
			for (int c = 0; c < out.Columns; c++) {
				mpImpl::to_float(out.Index(r,c), this->Index(r,c));
			}
		}
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

	// T is expected to be float or mpq_t type
	// If matrix is singular, returnEarly will stop this method early, returning 1
	// otherwise, 0 will be returned
	static int RowEchelonForm(Matrix<T>& Q, MPList<T>& augment, bool returnEarly){
		T zero, coefficient, sub;
		mpImpl::init(zero);
		mpImpl::init(coefficient);
		mpImpl::init(sub);
		int row = 0;
		for (int c = 0; c < Q.Columns; c++) {
			// find first non-zero pivot and swap row with current row
			bool noNonzero = true;
			for (int r = row; r < Q.Rows; r++) {
				if (mpImpl::cmp(Q.Index(r, c), zero) != 0) {
					noNonzero = false;
					if (r == row)
						break;
					// swap row r with row row
					for (int c2 = c; c2 < Q.Columns; c2++) {
						mpImpl::swap(Q.Index(row, c2), Q.Index(r, c2));
					}
					mpImpl::swap(augment[r], augment[row]);
					break;
				}
			}
			if (noNonzero) {
				if (returnEarly)
					return 1;
				continue;
			}
			// perform row reduction
			for (int r = row + 1; r < Q.Rows; r++) {
				mpImpl::div(coefficient, Q.Index(r, c), Q.Index(row, c));
				for (int c2 = c; c2 < Q.Columns; c2++) {
					mpImpl::mul(sub, coefficient, Q.Index(row, c2));
					mpImpl::sub(Q.Index(r, c2), Q.Index(r, c2), sub);
				}
				mpImpl::mul(sub, coefficient, augment[row]);
				mpImpl::sub(augment[r], augment[r], sub);
			}
			row++;
		}
		return 0;
	}

	// T is expected to be float or mpq_t type
	// reduce specifies whether the first entry in each pivot column should be set to 1
	// If a column contains no pivot and returnEarly is set, this method returns 1
	// otherwise, 0 will be returned
	static int RowEchelonForm(Matrix<T>& Q, vector<int>& pivots, bool reduce = true, bool returnEarly = false){
		T zero, coefficient, sub;
		mpImpl::init(zero);
		mpImpl::init(coefficient);
		mpImpl::init(sub);
		int row = 0;
		for (int c = 0; c < Q.Columns; c++) {
			// find first non-zero pivot and swap row with current row
			bool noNonzero = true;
			for (int r = row; r < Q.Rows; r++) {
				if (mpImpl::cmp(Q.Index(r, c), zero) != 0) {
					noNonzero = false;
					if (r == row)
						break;
					// swap row r with row row	
					for (int c2 = c; c2 < Q.Columns; c2++) {
						mpImpl::swap(Q.Index(row, c2), Q.Index(r, c2));
					}
					break;
				}
			}
			if (noNonzero) {
				printf("no pivot in column %d\n", c+1);
				if (returnEarly)
					return 1;
				continue;
			}
			pivots.push_back(c);
			// printf("pivot in column %d, term is ", c);
			// mpImpl::out_str(stdout, 10, Q.Index(row, c));
			// printf("\n");
			// perform row reduction
			if (reduce) {
				for (int c2 = c + 1; c2 < Q.Columns; c2++) {
					mpImpl::div(coefficient, Q.Index(row, c2), Q.Index(row, c));
					mpImpl::set(Q.Index(row, c2), coefficient);
				}
				mpImpl::set_ui(Q.Index(row, c), 1);
			}
			for (int r = row + 1; r < Q.Rows; r++) {
				mpImpl::div(coefficient, Q.Index(r, c), Q.Index(row, c));
				for (int c2 = c; c2 < Q.Columns; c2++) {
					mpImpl::mul(sub, coefficient, Q.Index(row, c2));
					mpImpl::sub(Q.Index(r, c2), Q.Index(r, c2), sub);
				}
			}
			row++;
			if (row == Q.Rows)
				break;
		}
		return 0;
	}

	MPList<mpq_t> BackSubstitution(MPList<mpq_t> values){
		MPList<mpq_t> solution = MPList<mpq_t>(values.Len);
		mpq_t temp1, temp2;
		mpq_inits(temp1, temp2, NULL);
		
		for (int r = values.Len-1; r >= 0; r--){
			mpq_set(temp1, values[r]);
			for (int c = values.Len-1; c > r; c--) {
				mpq_mul(temp2, Index(r,c), solution[c]);
				mpq_sub(temp1, temp1, temp2);
			}
			mpq_div(solution[r], temp1, Index(r,r));
		}

		return solution;
	}

	// expects self to have pivots set to 1
	void BackSubstitution(vector<int>& pivots){
		T coefficient, sub;
		mpImpl::init(coefficient);
		mpImpl::init(sub);
		for (int row = pivots.size() - 1; row > 0; row--) {
			int c = pivots[row];
			for (int r = row - 1; r >= 0; r--) {
				mpImpl::set(coefficient, Index(r, c));
				for (int c2 = c; c2 < Columns; c2++) {
					mpImpl::mul(sub, Index(row, c2), coefficient);
					mpImpl::sub(Index(r, c2), Index(r, c2), sub);
				}
			}
		}
	}

	// T is expected to be float or mpq_t type
	static void Determinant(Matrix<T>& Q, T& det) {
		T zero;
		mpImpl::init(zero);

		vector<int> _;
		int result = RowEchelonForm(Q, _, false, true);
		if (result == 1) {
			mpImpl::set(det, zero);
			return;
		}

		mpImpl::set_ui(det, 1, 1);
		for (int i = 0; i < Q.Rows; i++){
			mpImpl::mul(det, det, Q.Index(i, i));
		}
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
				mpImpl::to_mpf(f, Index(r,c));
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
	void Trace(mpz_t trace) const {
		T sum;
		mpImpl::init(sum);
		mpImpl::set(sum, Index(0,0));
		for(int r = 1; r < Rows; r++)
			mpImpl::add(sum, sum, Index(r,r));
		mpImpl::to_mpz(trace, sum);
	}
	void Trace(mpq_t trace) const {
		T sum;
		mpImpl::init(sum);
		mpImpl::set(sum, Index(0,0));
		for(int r = 1; r < Rows; r++)
			mpImpl::add(sum, sum, Index(r,r));
		mpImpl::to_mpq(trace, sum);
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
				// mpz_add_ui(rnd, rnd, 1);
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

	static string getNextTokenFromLine(size_t& start, string& line){
		// if (start < 0 || start >= line.size())
		// 	return "";
		string ws = " ";
		start = line.find_first_not_of(ws, start);
		size_t end = line.find_first_of(ws, start);
		string token = line.substr(start, end - start);
		start = end;
		return token;
	}
};


#endif