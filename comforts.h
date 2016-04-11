#ifndef SYNTAX_COMFORTS
#define SYNTAX_COMFORTS

#define NEW(TYPE, X) TYPE* X = (TYPE*) malloc(sizeof(TYPE))
#define MIN(X, Y) X < Y ? X:Y

// assumes that temp is declared
#define SWAP(X, Y) temp = X; X = Y; Y = temp

#endif