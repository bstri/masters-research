#include "mpImpl.h"
#include <gmp.h>

using namespace std;

void mpImpl::init(mpz_t op){
	mpz_init(op);
}
void mpImpl::init(mpq_t op){
	mpq_init(op);
}
void mpImpl::init(float& op){
	op = 0;
}

// These are not straightforward to implement
// See: "cpp variadic functions"
// void mpImpl::inits(mpz_t op, ...){
// 	mpz_init(op);
// }
// void mpImpl::inits(mpq_t op, ...){
// 	mpq_init(op);
// }

void mpImpl::clear(mpz_t op){
	mpz_clear(op);
}
void mpImpl::clear(mpq_t op){
	mpq_clear(op);
}
void mpImpl::clear(float op){

}

void mpImpl::from_mpz(mpz_t dest, const mpz_t src){
	mpz_set(dest, src);
}
void mpImpl::from_mpz(mpq_t dest, const mpz_t src){
	mpq_set_z(dest, src);
}
void mpImpl::from_mpz(float& dest, const mpz_t src){
	dest = (float) mpz_get_d(src);
}

void mpImpl::from_mpq(mpq_t dest, const mpq_t src){
	mpq_set(dest, src);
}
void mpImpl::from_mpq(mpz_t dest, const mpq_t src){
	mpz_set_q(dest, src);
}

void mpImpl::to_mpz(mpz_t dest, const mpz_t src){
	mpz_set(dest, src);
}
void mpImpl::to_mpz(mpz_t dest, const mpq_t src){
	mpz_set_q(dest, src);
}

void mpImpl::to_mpq(mpq_t dest, const mpz_t src){
	mpq_set_z(dest, src);
}
void mpImpl::to_mpq(mpq_t dest, const mpq_t src){
	mpq_set(dest, src);
}

void mpImpl::to_mpf(mpf_t dest, const mpz_t src){
	mpf_set_z(dest, src);
}
void mpImpl::to_mpf(mpf_t dest, const mpq_t src){
	mpf_set_q(dest, src);
}

void mpImpl::to_float(float& dest, const mpz_t src){
	dest = mpz_get_d(src);
}
void mpImpl::to_float(float& dest, const mpq_t src){
	dest = mpq_get_d(src);
}

void mpImpl::set(mpz_t dest, const mpz_t op){
	mpz_set(dest, op);
}
void mpImpl::set(mpq_t dest, const mpq_t op){
	mpq_set(dest, op);
}
void mpImpl::set(float& dest, const float op){
	dest = op;
}

void mpImpl::set_ui(mpz_t dest, unsigned long int op){
	mpz_set_ui(dest, op);
}
void mpImpl::set_ui(mpq_t dest, signed long n, unsigned long d){
	mpq_set_ui(dest, n, d);
}
void mpImpl::set_ui(float& dest, signed long n, unsigned long d){
	dest = n/d;
}

int mpImpl::set_str(mpz_t rop, const char *str, int base){
	return mpz_set_str(rop, str, base);
}
int mpImpl::set_str(mpq_t rop, const char *str, int base){
	return mpq_set_str(rop, str, base);
}

void mpImpl::add(mpz_t sum, const mpz_t op1, const mpz_t op2){
	mpz_add(sum, op1, op2);
}
void mpImpl::add(mpq_t sum, const mpq_t op1, const mpq_t op2){
	mpq_add(sum, op1, op2);
}

void mpImpl::sub(mpq_t dif, const mpq_t op1, const mpq_t op2){
	mpq_sub(dif, op1, op2);
}
void mpImpl::sub(mpz_t dif, const mpz_t op1, const mpz_t op2){
	mpz_sub(dif, op1, op2);
}
void mpImpl::sub(float& dif, const float op1, const float op2){
	dif = op1 - op2;
}

void mpImpl::mul(mpz_t product, const mpz_t op1, const mpz_t op2){
	mpz_mul(product, op1, op2);
}
void mpImpl::mul(mpq_t product, const mpq_t op1, const mpq_t op2){
	mpq_mul(product, op1, op2);
}
void mpImpl::mul(float& product, const float op1, const float op2){
	product = op1*op2;
}

void mpImpl::div(mpq_t quotient, const mpq_t op1, const mpq_t op2){
	mpq_div(quotient, op1, op2);
}
void mpImpl::div(float& quotient, const float op1, const float op2){
	quotient = op1/op2;
}

int mpImpl::cmp(const mpz_t op1, const mpz_t op2){
	return mpz_cmp(op1, op2);
}
int mpImpl::cmp(const mpq_t op1, const mpq_t op2){
	return mpq_cmp(op1, op2);
}
int mpImpl::cmp(const float op1, const float op2){
	return op1 - op2;
}

void mpImpl::swap(mpq_t op1, mpq_t op2){
	mpq_swap(op1, op2);
}
void mpImpl::swap(float& op1, float& op2){
	float temp = op1;
	op1 = op2;
	op2 = temp;
}

void mpImpl::out_str(FILE* stream, int base, const mpz_t op){
	mpz_out_str(stream, base, op);
}
void mpImpl::out_str(FILE* stream, int base, const mpq_t op){
	mpq_out_str(stream, base, op);
}
void mpImpl::out_str(FILE* stream, int base, float op){
	printf("%.1f", op);
}

void mpImpl::canonicalize(mpz_t op){
	// do nothing
}
void mpImpl::canonicalize(mpq_t op){
	mpq_canonicalize(op);
}
