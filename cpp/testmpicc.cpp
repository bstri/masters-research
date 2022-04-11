// #include <mpi.h>
#include <gmp.h>

using namespace std;

int main(int argc, char* argv[]) {
    mpz_t a, b;
    mpz_inits(a,b,NULL);
    mpz_add(a,a,b);
    return 0;
}