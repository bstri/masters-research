#include <gmp.h>
#include "Matrix.h"
#include "TraceAlgs.h"
#include "MPList.h"

using namespace std;

/*
Pointers for writing your own trace algorithms:
Write the algorithm in Matlab first. It will be easy to read, and it serves as a comparison for the c++ versions
Start with a list of desired matrix definitions. I like the following convention:
A2 = A^2
A3 = A2*A
...
A_2 = A.*A
A_3 = A_2.*A
...
D2 = diag(A2)
D3 = diag(A3)
...
Note that diag() is not the same as Matlab's "diag", but the result of A2.GetDiagonal (Matrix.h)

Then, create two temporary matrices, T1 and T2.
The strategy now is easiest to see in TraceAlgorithm2. Pick an operation, then store it alternatingly in T1 or T2:
A_2.RightMultiply(T1, A) --> T1 = A_2*A
T1.RightMultiply(T2, D2) --> T2 = T1*D2
T2.RightMultiply(T1, A) --> T1 = T2*A
Result --> T1 = A_2*A*D2*A

Algorithms written with this strategy may be candidates for (memory) optimization. See TraceAlgorithm1Optimized
*/

/*
Note: If it is known that all the terms must be integers for hollow symmetric matrices, 
then the program can be simplified by removing the mention of rationals and sped up since integer operations are quicker
*/

/* 
Note about algorithm1:
these algorithms assume the final answer is an integer for an integer input matrix
However, it assumes nothing about the individual terms that comprise the answer, 
and therefore does arithmetic on rationals. 
*/
MPList<mpq_t> TraceAlgorithm1(mpz_t out, const Matrix<mpz_t>& A){
	Matrix<mpz_t> A2(A.Rows, A.Columns); // A^2
	A.RightMultiply(A2, A);
	Matrix<mpz_t> A3(A.Rows, A.Columns); // A^3
	A.RightMultiply(A3, A2);
	Matrix<mpz_t> A4(A.Rows, A.Columns); // A^4
	A.RightMultiply(A4, A3);
	Matrix<mpz_t> A5(A.Rows, A.Columns); // A^5
	A.RightMultiply(A5, A4);

	Matrix<mpz_t> A_2(A.Rows, A.Columns); // A.A (Hadamard product)
	A.MultiplyEntrywise(A_2, A);
	Matrix<mpz_t> A_3(A.Rows, A.Columns); // A.A.A
	A.MultiplyEntrywise(A_3, A_2);
	Matrix<mpz_t> A_4(A.Rows, A.Columns); // A.A.A.A
	A.MultiplyEntrywise(A_4, A_3);
	Matrix<mpz_t> A_5(A.Rows, A.Columns); // A.A.A.A.A
	A.MultiplyEntrywise(A_5, A_4);

	Matrix<mpz_t> D2 = Matrix<mpz_t>::Diagonal(A2); // diag(A^2)
	Matrix<mpz_t> D3 = Matrix<mpz_t>::Diagonal(A3); // ...
	Matrix<mpz_t> D4 = Matrix<mpz_t>::Diagonal(A4);
	Matrix<mpz_t> D5 = Matrix<mpz_t>::Diagonal(A5);

	Matrix<mpz_t> T1(A.Rows, A.Columns);
	Matrix<mpz_t> T2(A.Rows, A.Columns);
	mpz_t t;
	mpz_init(t);
	mpq_t output, half;
	mpq_inits(output, half, NULL);
	mpq_set_ui(half, 1, 2);
	MPList<mpq_t> terms(14);

	// 2tr(A*A_3*A_3)
	A.RightMultiply(T1, A_3);
	T1.RightMultiply(T2, A_3);
	T2.Trace(t);
	mpz_mul_ui(t, t, 2);
	mpq_set_z(terms[0], t);
	mpq_add(output, output, terms[0]);

	// .5tr(A*(A^2.A^2.A^2))
	A2.MultiplyEntrywise(T1, A2);
	A2.MultiplyEntrywise(T2, T1);
	A.RightMultiply(T1, T2);
	T1.Trace(t);
	mpq_set_z(terms[1], t);
	mpq_mul(terms[1], terms[1], half);
	mpq_add(output, output, terms[1]);

	// 2tr(A^2*A_5)
	A2.RightMultiply(T1, A_5);
	T1.Trace(t);
	mpz_mul_ui(t, t, 2);
	mpq_set_z(terms[2], t);
	mpq_add(output, output, terms[2]);

	// -.5tr(D3.D4)
	D3.MultiplyEntrywise(T1, D4);
	T1.Trace(t);
	mpz_neg(t, t);
	mpq_set_z(terms[3], t);
	mpq_mul(terms[3], terms[3], half);
	mpq_add(output, output, terms[3]);

	// -tr(D3.diag((A_2)^2))
	A_2.RightMultiply(T1, A_2);
	Matrix<mpz_t> D = Matrix<mpz_t>::Diagonal(T1);
	D3.MultiplyEntrywise(T2, D);
	T2.Trace(t);
	mpz_neg(t, t);
	mpq_set_z(terms[4], t);
	mpq_add(output, output, terms[4]);

	// -4tr(A^2*D2*A_3)
	D2.RightMultiply(T1, A_3);
	A2.RightMultiply(T2, T1);
	T2.Trace(t);
	mpz_mul_si(t, t, -4);
	mpq_set_z(terms[5], t);
	mpq_add(output, output, terms[5]);

	// -.5tr(A*D2*A*A_3)
	A.RightMultiply(T1, A_3);
	D2.RightMultiply(T2, T1);
	A.RightMultiply(T1, T2);
	T1.Trace(t);
	mpz_neg(t, t);
	mpq_set_z(terms[6], t);
	mpq_mul(terms[6], terms[6], half);
	mpq_add(output, output, terms[6]);

	// tr(D3.D2.D2)
	D2.MultiplyEntrywise(T1, D2);
	D3.MultiplyEntrywise(T2, T1);
	T2.Trace(t);
	mpq_set_z(terms[7], t);
	mpq_add(output, output, terms[7]);

	// -2tr((A^2.A_2*A_2)*A)
	A2.MultiplyEntrywise(T1, A_2);
	T1.RightMultiply(T2, A_2);
	T2.RightMultiply(T1, A);
	T1.Trace(t);
	mpz_mul_si(t, t, -2);
	mpq_set_z(terms[8], t);
	mpq_add(output, output, terms[8]);

	// .5tr(A2*D2*A*D2)
	A2.RightMultiply(T1, D2);
	T1.RightMultiply(T2, A);
	T2.RightMultiply(T1, D2);
	T1.Trace(t);
	mpq_set_z(terms[9], t);
	mpq_mul(terms[9], terms[9], half);
	mpq_add(output, output, terms[9]);

	// 1.5tr((A^3.A^2)*A_2)
	A3.MultiplyEntrywise(T1, A2);
	T1.RightMultiply(T2, A_2);
	T2.Trace(t);
	mpz_mul_ui(t, t, 3);
	mpq_set_z(terms[10], t);
	mpq_mul(terms[10], terms[10], half);
	mpq_add(output, output, terms[10]);

	// .5tr(A*D3*A*D2)
	A.RightMultiply(T1, D3);
	T1.RightMultiply(T2, A);
	T2.RightMultiply(T1, D2);
	T1.Trace(t);
	mpq_set_z(terms[11], t);
	mpq_mul(terms[11], terms[11], half);
	mpq_add(output, output, terms[11]);

	// .5tr(A^4*A_3)
	A4.RightMultiply(T1, A_3);
	T1.Trace(t);
	mpq_set_z(terms[12], t);
	mpq_mul(terms[12], terms[12], half);
	mpq_add(output, output, terms[12]);

	// -.5tr(D5*A^2)
	D5.RightMultiply(T1, A2);
	T1.Trace(t);
	mpz_neg(t, t);
	mpq_set_z(terms[13], t);
	mpq_mul(terms[13], terms[13], half);
	mpq_add(output, output, terms[13]);

	mpz_set_q(out, output);
	return terms;
}

// This function is marginally faster than the original, but it also only uses 9 matrices (the original uses 14)
MPList<mpq_t> TraceAlgorithm1Optimized(mpz_t out, const Matrix<mpz_t>& A){
	Matrix<mpz_t> a(A.Rows, A.Columns); // A squared
	A.RightMultiply(a, A);
	Matrix<mpz_t> b(A.Rows, A.Columns); // A^3
	a.RightMultiply(b, A);
	Matrix<mpz_t> c(A.Rows, A.Columns); // A^4
	b.RightMultiply(c, A);
	Matrix<mpz_t> h(A.Rows, A.Columns); // A.A (Hadamard product)
	A.MultiplyEntrywise(h, A);
	Matrix<mpz_t> i(A.Rows, A.Columns); // A.A.A
	h.MultiplyEntrywise(i, A);
	Matrix<mpz_t> w = Matrix<mpz_t>::Diagonal(a); // diag(A^2)
	Matrix<mpz_t> x = Matrix<mpz_t>::Diagonal(b); // diag(A^3)
	Matrix<mpz_t> M1(A.Rows, A.Columns); // temporary storage
	Matrix<mpz_t> M2(A.Rows, A.Columns); // temporary storage
	
	mpz_t t;
	mpz_init(t);
	mpq_t total, half;
	mpq_inits(total, half, NULL);
	mpq_set_ui(half, 1, 2);
	MPList<mpq_t> terms(14);

	// 1.5tr((A^3.A^2)*A_2)
	b.MultiplyEntrywise(M1, a);
	M1.RightMultiply(M2, h);
	M2.Trace(t);
	mpz_mul_ui(t, t, 3);
	mpq_set_z(terms[0], t);
	mpq_mul(terms[0], terms[0], half);
	mpq_add(total, total, terms[0]);

	// .5tr(A^4*A_3)
	c.RightMultiply(M1, i);
	M1.Trace(t);
	mpq_set_z(terms[1], t);
	mpq_mul(terms[1], terms[1], half);
	mpq_add(total, total, terms[1]);

	c.GetDiagonal(M1); // M1 = D4

	// -.5tr(D3.D4)
	x.MultiplyEntrywise(M2, M1);
	M2.Trace(t);
	mpz_neg(t, t);
	mpq_set_z(terms[2], t);
	mpq_mul(terms[2], terms[2], half);
	mpq_add(total, total, terms[2]);

	c.RightMultiply(M1,A);
	M1.GetDiagonal(c); // c = D5

	// -.5tr(D5*A^2)
	c.RightMultiply(M1, a);
	M1.Trace(t);
	mpz_neg(t, t);
	mpq_set_z(terms[3], t);
	mpq_mul(terms[3], terms[3], half);
	mpq_add(total, total, terms[3]);

	i.MultiplyEntrywise(M2,A);
	M2.MultiplyEntrywise(b,A); // b = A_5

	// 2tr(A^2*A_5)
	a.RightMultiply(M1, b);
	M1.Trace(t);
	mpz_mul_ui(t, t, 2);
	mpq_set_z(terms[4], t);
	mpq_add(total, total, terms[4]);

	// -4tr(A^2*D2*A_3)
	w.RightMultiply(M1, i);
	a.RightMultiply(M2, M1);
	M2.Trace(t);
	mpz_mul_si(t, t, -4);
	mpq_set_z(terms[5], t);
	mpq_add(total, total, terms[5]);

	// -.5tr(A*D2*A*A_3)
	w.RightMultiply(b,A); // b = D2*A
	b.RightMultiply(M1, i);
	A.RightMultiply(M2, M1);
	M2.Trace(t);
	mpz_neg(t, t);
	mpq_set_z(terms[6], t);
	mpq_mul(terms[6], terms[6], half);
	mpq_add(total, total, terms[6]);

	// 2tr(A*A_3*A_3)
	A.RightMultiply(M1, i);
	M1.RightMultiply(M2, i);
	M2.Trace(t);
	mpz_mul_ui(t, t, 2);
	mpq_set_z(terms[7], t);
	mpq_add(total, total, terms[7]);

	// .5tr(A2*D2*A*D2)
	a.RightMultiply(M1, b);
	M1.RightMultiply(M2, w);
	M2.Trace(t);
	mpq_set_z(terms[8], t);
	mpq_mul(terms[8], terms[8], half);
	mpq_add(total, total, terms[8]);

	// .5tr(A*(A^2.A^2.A^2))
	a.MultiplyEntrywise(M1, a);
	a.MultiplyEntrywise(M2, M1);
	A.RightMultiply(M1, M2);
	M1.Trace(t);
	mpq_set_z(terms[9], t);
	mpq_mul(terms[9], terms[9], half);
	mpq_add(total, total, terms[9]);

	// -2tr((A^2.A_2*A_2)*A)
	a.MultiplyEntrywise(M1, h);
	M1.RightMultiply(M2, h);
	M2.RightMultiply(M1, A);
	M1.Trace(t);
	mpz_mul_si(t, t, -2);
	mpq_set_z(terms[10], t);
	mpq_add(total, total, terms[10]);

	// -tr(D3.diag((A_2)^2))
	h.RightMultiply(M1, h);
	M1.GetDiagonal(M2);
	x.MultiplyEntrywise(M1, M2);
	M1.Trace(t);
	mpz_neg(t, t);
	mpq_set_z(terms[11], t);
	mpq_add(total, total, terms[11]);

	// tr(D3.D2.D2)
	w.MultiplyEntrywise(M1, w);
	x.MultiplyEntrywise(M2, M1);
	M2.Trace(t);
	mpq_set_z(terms[12], t);
	mpq_add(total, total, terms[12]);

	// .5tr(A*D3*A*D2)
	A.RightMultiply(M1, x);
	M1.RightMultiply(M2, A);
	M2.RightMultiply(M1, w);
	M1.Trace(t);
	mpq_set_z(terms[13], t);
	mpq_mul(terms[13], terms[13], half);
	mpq_add(total, total, terms[13]);
	
	mpz_set_q(out, total);
	return terms;
}

// This function assumes that all division is exact with no remainder
// Non-exact division is undefined, so results will be unusable if this happens
MPList<mpz_t> TraceAlgorithm2(mpz_t out, const Matrix<mpz_t>& A){
	Matrix<mpz_t> A2(A.Rows, A.Columns); // A^2
	A.RightMultiply(A2, A);
	Matrix<mpz_t> A3(A.Rows, A.Columns); // A^3
	A.RightMultiply(A3, A2);
	Matrix<mpz_t> A4(A.Rows, A.Columns); // A^4
	A.RightMultiply(A4, A3);
	Matrix<mpz_t> A6(A.Rows, A.Columns); // A^6
	A3.RightMultiply(A6, A3);

	Matrix<mpz_t> A_2(A.Rows, A.Columns); // A.A (Hadamard product)
	A.MultiplyEntrywise(A_2, A);
	Matrix<mpz_t> A_3(A.Rows, A.Columns); // A.A.A
	A.MultiplyEntrywise(A_3, A_2);

	Matrix<mpz_t> D2 = Matrix<mpz_t>::Diagonal(A2); // diag(A^2)
	Matrix<mpz_t> D3 = Matrix<mpz_t>::Diagonal(A3); // ...
	Matrix<mpz_t> D4 = Matrix<mpz_t>::Diagonal(A4);

	Matrix<mpz_t> T1(A.Rows, A.Columns);
	Matrix<mpz_t> T2(A.Rows, A.Columns);
	MPList<mpz_t> terms(10);

	// 1/3*tr(d2.*d2.*d2)
	D2.MultiplyEntrywise(T1, D2);
	T1.MultiplyEntrywise(T2, D2);
	T2.Trace(terms[0]);
	mpz_divexact_ui(terms[0], terms[0], 3);
	mpz_add(out, out, terms[0]);

	// 1/2*tr(A_3*A3)
	A_3.RightMultiply(T1, A3);
	T1.Trace(terms[1]);
	mpz_divexact_ui(terms[1], terms[1], 2);
	mpz_add(out, out, terms[1]);

	// -tr(d2.*diag(A_2*A_2))
	A_2.RightMultiply(T1, A_2);
	T1.GetDiagonal(T2);
	D2.MultiplyEntrywise(T1, T2);
	T1.Trace(terms[2]);
	mpz_neg(terms[2], terms[2]);
	mpz_add(out, out, terms[2]);

	// 1/3*tr(A_3*A_3)
	A_3.RightMultiply(T1, A_3);
	T1.Trace(terms[3]);
	mpz_divexact_ui(terms[3], terms[3], 3);
	mpz_add(out, out, terms[3]);

	// 1/4*tr(d2.*diag(A*d2*A))
	A.RightMultiply(T1, D2);
	T1.RightMultiply(T2, A);
	T2.GetDiagonal(T1);
	D2.MultiplyEntrywise(T2, T1);
	T2.Trace(terms[4]);
	mpz_divexact_ui(terms[4], terms[4], 4);
	mpz_add(out, out, terms[4]);

	// -1/3*tr(A_2*A_2*A_2)
	A_2.RightMultiply(T1, A_2);
	T1.RightMultiply(T2, A_2);
	T2.Trace(terms[5]);
	mpz_divexact_ui(terms[5], terms[5], 3);
	mpz_neg(terms[5], terms[5]);
	mpz_add(out, out, terms[5]);

	// 3/4*tr((A2.*A2)*A_2)
	A2.MultiplyEntrywise(T1, A2);
	T1.RightMultiply(T2, A_2);
	T2.Trace(terms[6]);
	mpz_mul_ui(terms[6], terms[6], 3);
	mpz_divexact_ui(terms[6], terms[6], 4);
	mpz_add(out, out, terms[6]);

	// -1/4*tr(d3.*d3)
	D3.MultiplyEntrywise(T1, D3);
	T1.Trace(terms[7]);
	mpz_divexact_ui(terms[7], terms[7], 4);
	mpz_neg(terms[7], terms[7]);
	mpz_add(out, out, terms[7]);

	// -1/2*tr(d2.*d4);
	D2.MultiplyEntrywise(T1, D4);
	T1.Trace(terms[8]);
	mpz_divexact_ui(terms[8], terms[8], 2);
	mpz_neg(terms[8], terms[8]);
	mpz_add(out, out, terms[8]);

	// 1/12*tr(A6);
	A6.Trace(terms[9]);
	mpz_divexact_ui(terms[9], terms[9], 12);
	mpz_add(out, out, terms[9]);

	return terms;
}

MPList<mpz_t> TraceAlgorithm3(const Matrix<mpz_t>& A){
	Matrix<mpz_t> A2(A.Rows, A.Columns);
	A.RightMultiply(A2, A);
	Matrix<mpz_t> T1(A.Rows, A.Columns);
	Matrix<mpz_t> T2(A.Rows, A.Columns);
	MPList<mpz_t> terms(4);

	// tr(A2*A_3)
	A.MultiplyEntrywise(T1, A);
	T1.MultiplyEntrywise(T2, A); // T2 = A_3
	A2.RightMultiply(T1, T2); // T1 = A*A_3
	T1.Trace(terms[2]);

	// tr(A3)*tr(A2)
	A2.Trace(terms[3]);
	A2.RightMultiply(T1, A); // T1 = A3
	T1.Trace(terms[0]);
	mpz_mul(terms[3], terms[3], terms[0]);

	// tr(A5)
	A2.RightMultiply(T2, T1); // T2 = A5
	T2.Trace(terms[0]);

	T1.GetDiagonal(T2); // T2 = d3
	A2.GetDiagonal(T1); // T1 = d2
	T2.MultiplyEntrywise(A2, T1); // A2 = d3.*d2 since we don't need it anymore
	A2.Trace(terms[1]);

	return terms;
}

MPList<mpz_t> CompiledTraceAlgorithm3(const Matrix<mpz_t>& A){
    Matrix<mpz_t> _T1(A.Rows, A.Columns);
    Matrix<mpz_t> _T2(A.Rows, A.Columns);
    MPList<mpz_t> terms(3);

    // A@A@A@A@A
    A.RightMultiply(_T1, A);
    _T1.RightMultiply(_T2, A);
    _T2.RightMultiply(_T1, A);
    _T1.RightMultiply(_T2, A);
    _T2.Trace(terms[0]);

    // d(A@A@A)*d(A@A)
    Matrix<mpz_t> _T3(A.Rows, A.Columns);
    A.RightMultiply(_T1, A);
    _T1.RightMultiply(_T2, A);
    _T2.GetDiagonal(_T1);
    A.RightMultiply(_T2, A);
    _T2.GetDiagonal(_T3);
    _T1.MultiplyEntrywise(_T2, _T3);
    _T2.Trace(terms[1]);

    // (A@A)@(A*A*A)
    A.RightMultiply(_T1, A);
    A.MultiplyEntrywise(_T2, A);
    _T2.MultiplyEntrywise(_T3, A);
    _T1.RightMultiply(_T2, _T3);
    _T2.Trace(terms[2]);

    return terms;
}

MPList<mpz_t> CompiledTraceAlgorithm(const Matrix<mpz_t>& M){
	Matrix<mpz_t> _T1(M.Rows, M.Columns);
	Matrix<mpz_t> _T2(M.Rows, M.Columns);
	MPList<mpz_t> terms(5);
	
	// M*d(M*M)*M*(M@M@M)
	Matrix<mpz_t> _T3(M.Rows, M.Columns);
	M.MultiplyEntrywise(_T1, M);
	_T1.GetDiagonal(_T2);
	M.MultiplyEntrywise(_T1, _T2);
	_T1.MultiplyEntrywise(_T2, M);
	M.RightMultiply(_T1, M);
	_T1.RightMultiply(_T3, M);
	_T2.MultiplyEntrywise(_T1, _T3);
	_T1.Trace(terms[0]);
	
	// M*d(M*M*(M@M@M))*M
	M.MultiplyEntrywise(_T1, M);
	M.RightMultiply(_T2, M);
	_T2.RightMultiply(_T3, M);
	_T1.MultiplyEntrywise(_T2, _T3);
	_T2.GetDiagonal(_T1);
	M.MultiplyEntrywise(_T2, _T1);
	_T2.MultiplyEntrywise(_T1, M);
	_T1.Trace(terms[1]);
	
	// M*(M@M@M)*M*M*M
	M.RightMultiply(_T1, M);
	_T1.RightMultiply(_T2, M);
	M.MultiplyEntrywise(_T1, _T2);
	_T1.MultiplyEntrywise(_T2, M);
	_T2.MultiplyEntrywise(_T1, M);
	_T1.MultiplyEntrywise(_T2, M);
	_T2.Trace(terms[2]);
	
	// M*((M*M*M)@(M@M))*M
	M.MultiplyEntrywise(_T1, M);
	_T1.MultiplyEntrywise(_T2, M);
	M.RightMultiply(_T1, M);
	_T2.RightMultiply(_T3, _T1);
	M.MultiplyEntrywise(_T2, _T3);
	_T2.MultiplyEntrywise(_T3, M);
	_T3.Trace(terms[3]);
	
	// M*M*d(M*M)*M*M*M
	M.MultiplyEntrywise(_T1, M);
	M.MultiplyEntrywise(_T2, M);
	_T2.GetDiagonal(_T3);
	_T1.MultiplyEntrywise(_T2, _T3);
	_T2.MultiplyEntrywise(_T1, M);
	_T1.MultiplyEntrywise(_T2, M);
	_T2.MultiplyEntrywise(_T1, M);
	_T1.Trace(terms[4]);
	
	return terms;
}

MPList<mpz_t> CompiledTraceAlgorithmShorthand(const Matrix<mpz_t>& M){
	Matrix<mpz_t> _T1(M.Rows, M.Columns);
	Matrix<mpz_t> _T2(M.Rows, M.Columns);
	MPList<mpz_t> terms(5);
	
	Matrix<mpz_t> M2(M.Rows, M.Columns);
	M.RightMultiply(M2, M);
	Matrix<mpz_t> M3(M.Rows, M.Columns);
	M.RightMultiply(M3, M2);
	
	Matrix<mpz_t> M_2(M.Rows, M.Columns);
	M.MultiplyEntrywise(M_2, M);
	Matrix<mpz_t> M_3(M.Rows, M.Columns);
	M.MultiplyEntrywise(M_3, M_2);
	
	// M*d(M_2)*M*M3
	M_2.GetDiagonal(_T1);
	M.MultiplyEntrywise(_T2, _T1);
	_T2.MultiplyEntrywise(_T1, M);
	_T1.MultiplyEntrywise(_T2, M3);
	_T2.Trace(terms[0]);
	
	// M*d(M_2*M3)*M
	M_2.MultiplyEntrywise(_T1, M3);
	_T1.GetDiagonal(_T2);
	M.MultiplyEntrywise(_T1, _T2);
	_T1.MultiplyEntrywise(_T2, M);
	_T2.Trace(terms[1]);
	
	// M*M3*M_3
	M.MultiplyEntrywise(_T1, M3);
	_T1.MultiplyEntrywise(_T2, M_3);
	_T2.Trace(terms[2]);
	
	// M*(M_3@M2)*M
	M_3.RightMultiply(_T1, M2);
	M.MultiplyEntrywise(_T2, _T1);
	_T2.MultiplyEntrywise(_T1, M);
	_T1.Trace(terms[3]);
	
	// M_2*d(M_2)*M_3
	M_2.GetDiagonal(_T1);
	M_2.MultiplyEntrywise(_T2, _T1);
	_T2.MultiplyEntrywise(_T1, M_3);
	_T1.Trace(terms[4]);
	
	return terms;
}

MPList<mpz_t> CompleteGraph3Cycles(const Matrix<mpz_t>& A){
	Matrix<mpz_t> _T1(A.Rows, A.Columns);
	Matrix<mpz_t> _T2(A.Rows, A.Columns);
	MPList<mpz_t> terms(1);
	
	A.RightMultiply(_T1, A);
	Matrix<mpz_t> A3(A.Rows, A.Columns);
	A.RightMultiply(A3, _T1);
	
	// A3
	A3.Trace(terms[0]);
	
	return terms;
}