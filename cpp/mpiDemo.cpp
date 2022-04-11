// mpicxx mpiDemo.cpp mpImpl.cpp -std=c++11 -I /opt/ohpc/pub/libs/spack/opt/spack/linux-centos7-sandybridge/gcc-4.8.5/gmp-6.1.2-hc7piarcvdeoqjtu3ikpsrapvkbivdbm/include/ -L /opt/ohpc/pub/libs/spack/opt/spack/linux-centos7-sandybridge/gcc-4.8.5/gmp-6.1.2-hc7piarcvdeoqjtu3ikpsrapvkbivdbm/lib/ -lgmpxx -lgmp -o mpiDemo 
// mpiexec -np 2 ./mpiDemo 4
#include <mpi.h>
#include <gmp.h>
#include "Matrix.h"
#include <ctime>
#include <cstdlib>
#include <cstdio>

using namespace std;

int main(int argc, char* argv[]){
    mpz_t maxRnd, rnd;
	mpz_inits(maxRnd, rnd, NULL);
	mpz_set_ui(maxRnd, 10);
	gmp_randstate_t r;
	gmp_randinit_default(r);
	gmp_randseed_ui(r, time(NULL));

    int n = atoi(argv[1]);
    Matrix<mpz_t> A = Matrix<mpz_t>::Random(n, r, maxRnd);
    Matrix<mpz_t> B = Matrix<mpz_t>::Random(n, r, maxRnd);
	Matrix<mpz_t> result(n);

	clock_t t = clock();

	MPI_Init(&argc, &argv);

	int np;
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	if (id == 0) {
		printf("there are %d processes\n", np);
	}
	printf("%p", (void *)&A);

	int a = 0;
	mpz_t sum, prod;
	mpz_inits(sum, prod, NULL);
	int i;
	while ((i = a*np + id) < A.Rows*B.Columns) {
		printf("RANK %d \t i %d \t a %d\n", id, i, a);
		div_t entry = div(i, B.Columns);
		mpz_set_ui(sum, 0);
		for (int j = 0; j < A.Columns; j++) {
			mpz_mul(prod, A.Index(entry.quot, j), B.Index(j, entry.rem));
			mpz_add(sum, sum, prod);
		}
		mpz_set(result.Data[i], sum);
		a++;
	}

	// if (id == 0) {
	// 	t = clock() - t;
	// }

	MPI_Finalize();

	if (id == 0) {
		printf("elapsed time - %d\n", clock() - t);

		if (n <= 10) {
			A.Print();
			printf("\n");
			B.Print();
			printf("\n");
			result.Print();
		}
	} else if (id == np - 1) {
		// printf("\n");
		// result.Print();
	}

	return 0;
}