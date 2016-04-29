all : tsp  vertex_color

clean :
	rm -rf *.o
	rm -rf tsp
	rm -rf vertex_color

tsp : tsp.c solver.c leftist_heap.c graph.c bitvector.c list.c
	mpicc -g -Wall -fopenmp -o tsp tsp.c solver.c leftist_heap.c graph.c bitvector.c list.c -lm

vertex_color: graph_color.c solver.c leftist_heap.c graph.c bitvector.c
	mpicc -g -Wall -fopenmp -o vertex_color graph_color.c solver.c leftist_heap.c graph.c bitvector.c list.c -lm
