#ifndef SYNTAX_COMFORTS
#define SYNTAX_COMFORTS

#include <mpi.h>

#define NEW(TYPE, X) TYPE* X = (TYPE*) malloc(sizeof(TYPE))
//#define ABSVAL(X, Y) (X>Y)?(X-Y):(Y-X)

int MIN(int a, int b);
int MAX(int a, int b);
int ABSVAL(int a, int b);
// assumes that temp is declared
#define SWAP(X, Y) temp = X; X = Y; Y = temp

extern int err;

#define REPORT_ERROR(X, ERR_MSG, RET_VAL) \
					if ((err = X) != MPI_SUCCESS) { \
						printf(ERR_MSG); \
						MPI_Abort(MPI_COMM_WORLD, err); \
						return RET_VAL; \
					}

#endif
