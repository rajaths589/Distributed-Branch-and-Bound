all : tsp  vertex_color

clean :
	rm -rf *.o
	rm -rf tsp
	rm -rf vertex_color

tsp : tsp.c
	gcc -g -c tsp.c

vertex_color: vertex_color.c
	gcc -g -c vertex_color.c

leftist_heap.o : leftist_heap.c
	gcc -g -c -fopenmp leftist_heap.c

solver.o : solver.c
	gcc -g -c -fopenmp solver.c
