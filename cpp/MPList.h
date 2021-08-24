#ifndef MPLIST_H
#define MPLIST_H

#include <cstdio>
#include <cstdlib>
#include <gmp.h>
#include "mpImpl.h"

using namespace std;

template <typename T>
class MPList{
public:
	int Len;
	T* Data;

	// Constructor
	MPList(int numEntries){
		Data = (T*) malloc(numEntries*sizeof(T));
		for(int i = 0; i < numEntries; i++)
			mpImpl::init(Data[i]);
		Len = numEntries;
	}

	// Destructor
	~MPList(){
		Deallocate();
	}

	// Allow indexing this list
	T& operator[](int index){
		return Data[index];
	}

	// Free up memory used by this list
	void Deallocate(){
		for(int i = 0; i < Len; i++)
			mpImpl::clear(Data[i]);
		free(Data);
		Data = NULL;
	}

	// Print this list to stdout
	void Print(const char* delim = "  "){
		for(int i = 0; i < Len; i++){
			mpImpl::out_str(stdout, 10, Data[i]);
			printf("%s", delim);
		}
		printf("\n");
	}

	// Print list in scientific notation
	void PrintSciNot(mp_bitcnt_t precBits, const char* delim = "  "){
		mpf_t f;
		mpf_init2(f, precBits);
		for(int i = 0; i < Len; i++){
			mpImpl::to_float(f, Data[i]);
			mpf_out_str(stdout, 10, 0, f);
			printf("%s", delim);
		}
		printf("\n");
	}
};

#endif