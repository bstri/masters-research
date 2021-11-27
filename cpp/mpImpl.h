#ifndef MPIMPL_H
#define MPIMPL_H

#include <gmp.h>

namespace mpImpl {
	void init(mpz_t op);
	void init(mpq_t op);
	void init(float& op);

	// These are not straightforward to implement
	// See: "cpp variadic functions"
	// void inits(mpz_t op, ...)
	// void inits(mpq_t op, ...)

	void clear(mpz_t op);
	void clear(mpq_t op);
	void clear(float op);

	void from_mpz(mpz_t dest, const mpz_t src);
	void from_mpz(mpq_t dest, const mpz_t src);
	void from_mpz(float& dest, const mpz_t src);

	void from_mpq(mpq_t dest, const mpq_t src);
	void from_mpq(mpz_t dest, const mpq_t src);

	void to_mpq(mpq_t dest, const mpz_t src);
	void to_mpq(mpq_t dest, const mpq_t src);

	void to_mpf(mpf_t dest, const mpz_t src);
	void to_mpf(mpf_t dest, const mpq_t src);

	void to_float(float& dest, const mpz_t src);
	void to_float(float& dest, const mpq_t src);

	void set(mpz_t dest, const mpz_t op);
	void set(mpq_t dest, const mpq_t op);

	void set_ui(mpz_t dest, unsigned long op);
	void set_ui(mpq_t dest, signed long n, unsigned long d = 1);
	void set_ui(float& dest, signed long n, unsigned long d = 1);

	int set_str(mpz_t rop, const char *str, int base);
	int set_str(mpq_t rop, const char *str, int base);

	void add(mpz_t sum, const mpz_t op1, const mpz_t op2);
	void add(mpq_t sum, const mpq_t op1, const mpq_t op2);

	void sub(mpq_t dif, const mpq_t op1, const mpq_t op2);
	void sub(mpz_t dif, const mpz_t op1, const mpz_t op2);
	void sub(float& dif, const float op1, const float op2);

	void mul(mpz_t product, const mpz_t op1, const mpz_t op2);
	void mul(mpq_t product, const mpq_t op1, const mpq_t op2);
	void mul(float& product, const float op1, const float op2);

	void div(mpq_t quotient, const mpq_t op1, const mpq_t op2);
	void div(float& quotient, const float op1, const float op2);

	int cmp(const mpz_t op1, const mpz_t op2);
	int cmp(const mpq_t op1, const mpq_t op2);
	int cmp(const float op1, const float op2);

	void swap(mpq_t op1, mpq_t op2);
	void swap(float& op1, float& op2);

	void out_str(FILE* stream, int base, const mpz_t op);
	void out_str(FILE* stream, int base, const mpq_t op);
	void out_str(FILE* stream, int base, float op);

	void canonicalize(mpz_t op);
	void canonicalize(mpq_t op);
}

#endif