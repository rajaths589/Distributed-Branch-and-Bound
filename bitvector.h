#ifndef BIT_VEC_DEFN
#define BIT_VEC_DEFN

typedef struct bitvector {
	char* bitarray;
	int len;
	int sz;
} bitvector;

bitvector* create_bitvector(int len);
void destroy_bitvector(bitvector* bvec);
int getIndex(bitvector* bvec, int index);
void setIndex(bitvector* bvec, int index);
void unsetIndex(bitvector* bvec, int index);
bitvector* create_copy(bitvector* bvec);

#endif
