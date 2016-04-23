// BATCH 5
// 2012B5A7589P : RAJATH S
// 2012B3A7463P : KOKANDAKAR AJINKYA HEMANT

#include "bitvector.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

bitvector* create_bitvector(int len) {
	assert(len >= 0);

	int sz;
	if (len % 8 == 0) {
		sz = len/8;
	} else {
		sz = len/8 + 1;
	}

	bitvector* bvec = (bitvector*) malloc(sizeof(bitvector));
	bvec->bitarray = (char*) calloc(sz, sizeof(char));
	bvec->len = sz;
	bvec->sz = len;

	return bvec;
}

void destroy_bitvector(bitvector* bvec) {
	free(bvec->bitarray);
	free(bvec);
}

int getIndex(bitvector* bvec, int index) {
	assert(index >= 0);
	assert(index < bvec->sz);

	int div = index / 8;
	int rem = index % 8;

	char bucket = bvec->bitarray[div];
	char mask = (0x80 >> rem);

	if ((mask & bucket) == 0)
		return 0;
	else
		return 1;
}

void setIndex(bitvector* bvec, int index) {
	assert(index >= 0);
	assert(index < bvec->sz);

	int div = index / 8;
	int rem = index % 8;

	char bucket = bvec->bitarray[div];
	char mask = (0x80 >> rem);

	bucket = bucket | mask;

	bvec->bitarray[div] = bucket;
}

void unsetIndex(bitvector* bvec, int index) {
	assert(index >= 0);
	assert(index < bvec->sz);

	int div = index / 8;
	int rem = index % 8;

	char bucket = bvec->bitarray[div];
	char mask = ~(0x80 >> rem);

	bucket = bucket & mask;

	bvec->bitarray[div] = bucket;
}

bitvector* create_copy(bitvector* bvec) {
	bitvector* bvec1 = (bitvector*) malloc(sizeof(bitvector));
	bvec1->len = bvec->len;
	bvec1->sz = bvec->sz;
	bvec1->bitarray = (char*) malloc(bvec1->len * sizeof(char));
	memcpy(bvec1->bitarray, bvec->bitarray, bvec->len);

	return bvec1;
}
