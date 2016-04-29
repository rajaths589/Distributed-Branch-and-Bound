all : tsp  vertex_color
debug : tsp-debug vertex_color-debug

clean :
	rm -rf *.o
	rm -rf tsp
	rm -rf vertex_color

tsp : tsp.c solver.c leftist_heap.c graph.c bitvector.c list.c
	mpicc -g -Wall -fopenmp -o tsp tsp.c solver.c leftist_heap.c graph.c bitvector.c list.c -lm

vertex_color: graph_color.c solver.c leftist_heap.c graph.c bitvector.c list.c
	mpicc -g -Wall -fopenmp -o vertex_color graph_color.c solver.c leftist_heap.c graph.c bitvector.c list.c -lm

tsp-debug: tsp.c solver.c leftist_heap.c graph.c bitvector.c list.c
	vtcc -vt:cc mpicc -openmp -Wall -fopenmp -o tsp tsp.c solver.c leftist_heap.c graph.c bitvector.c list.c -lm

vertex_color-debug: graph_color.c solver.c leftist_heap.c graph.c bitvector.c list.c
	vtcc -vt:cc mpicc -openmp -Wall -fopenmp -o vertex_color graph_color.c solver.c leftist_heap.c graph.c bitvector.c list.c -lm
