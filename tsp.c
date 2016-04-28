#include "generics.h"
#include "comforts.h"
#include "graph.h"
#include "bitvector.h"
#include "priority_queue.h"

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

int max_length;
int solution_vector_size = 0;

typedef struct tsp_data {
	graph* g;
	int start_point;
	int max_len;	// equal to num vertices
} tsp_data;

struct bound_array {
	int* path;
	int curr_length;
	int max_length;
	bitvector* used_vertices;
};

typedef struct bound_array tsp_path;

tsp_path* create_soln_copy(tsp_path* p) {
	NEW(tsp_path, p1);
	p1->curr_length = p->curr_length;
	p1->max_length = p->max_length;
	p1->path = calloc(p1->max_length, sizeof(int));
	memcpy(p1->path, p->path, p1->max_length*sizeof(int));
	p1->used_vertices = create_copy(p->used_vertices);

	return p1;
}

void* populate_domain_data(int argc, char** argv) {
	assert(argc == 2);
	FILE* fp = fopen(argv[1], "r");

	assert(fp != NULL);

	int num_vertices, num_edges;
	NEW(tsp_data, data);

	fscanf(fp, "%d %d", &num_vertices, &num_edges);
	graph* g = create_graph(num_vertices);

	int from, to;
	float weight;

	for (int i = 0; i < num_edges; i++) {
		fscanf(fp, "%d %d %f", &from, &to, &weight);
		add_edge(g, from, to, weight);
	}

	data->g = g;
	data->max_len = g->num_vertices;

	fscanf(fp, "%d", &data->start_point);

	assert(data->start_point >= 0);
	assert(data->start_point < num_vertices);

	solution_vector_size = sizeof(float) + sizeof(int) + (data->max_len+1)*sizeof(int);

	return (void*) data;
}

solution_vector get_root_partial_solution(void* domain_specific_data) {
	NEW(tsp_path, empty_solution);

	tsp_data* data = (tsp_data*) domain_specific_data;

	empty_solution->max_length = data->max_len + 1;
	empty_solution->curr_length = 1;
	empty_solution->path = (int*) calloc(data->max_len + 1, sizeof(int));
	empty_solution->path[0] = data->start_point;
	empty_solution->used_vertices = create_bitvector(data->max_len);
	setIndex(empty_solution->used_vertices, data->start_point);

	return (void*) empty_solution;
}

int construct_candidates(solution_vector partial_solution, float partial_soln_score,
						 queue* private_queue, void* domain_specific_data) {

	tsp_path* partial = (tsp_path*) partial_solution;
	tsp_data* data = (tsp_data*) domain_specific_data;

	if (partial->curr_length == data->max_len+1)
		return 1;

	list_node* l;
	tsp_path* path1;
	l = data->g->adjacency_list[partial->path[partial->curr_length-1]];

	if (partial->curr_length == data->max_len) {
		while (l != NULL) {
			if (l->to == data->start_point) {
				partial->path[partial->curr_length] = l->to;
				partial->curr_length++;
				pq_insert_nc(private_queue, partial_soln_score + l->weight, partial);
				break;
			}
			l = l->next;
		}
	} else {
		while (l != NULL) {
			if (!getIndex(partial->used_vertices, l->to)) {
				path1 = create_soln_copy(partial);
				path1->path[path1->curr_length] = l->to;
				path1->curr_length++;
				setIndex(path1->used_vertices, l->to);
				pq_insert_nc(private_queue, partial_soln_score + l->weight, path1);
			}
			l = l->next;
		}

		free(partial);
	}

	return 0;
}

void print_solution(solution_vector sol, float score) {
	tsp_path* solution = (tsp_path*) sol;

	printf("Curr length : %d\n", solution->curr_length);
	for (int i = 0; i < solution->curr_length; i++) {
		printf("%d\t", solution->path[i]);
	}
	printf("\n");
	printf("Cost : %f\n", score);
}

int safe_solution_buffer_size(tsp_data* data) {
	return (data->max_len+1)*sizeof(int) + sizeof(float) + sizeof(int);
}

//buff size will be decided in consultation with user
void pack_solution(void* buff, int buff_size, int* pos, solution_vector vec, float score,
				   MPI_Comm comm, void* problem_data) {

	tsp_path* path = (tsp_path*) vec;
//	print_solution(path, score);
	MPI_Pack(&score, 1, MPI_FLOAT, buff, buff_size, pos, comm);
	MPI_Pack(&path->curr_length, 1, MPI_INT, buff, buff_size, pos, comm);
	MPI_Pack(path->path, path->max_length, MPI_INT, buff, buff_size, pos, comm);
}

solution_vector unpack_solution(void* buff, int buff_size, MPI_Comm comm, int* pos, float*
score,
						void* problem_data) {

	tsp_path* p = (tsp_path*) get_root_partial_solution(problem_data);
	MPI_Unpack(buff, buff_size, pos, score, 1, MPI_FLOAT, comm);
	MPI_Unpack(buff, buff_size, pos, &p->curr_length, 1, MPI_INT, comm);
	MPI_Unpack(buff, buff_size, pos, p->path, p->max_length, MPI_INT, comm);
//	print_solution(p, *score);

	//populate bitvector
//	printf("%d %d\n", p->curr_length, p->max_length);
//	assert(p->curr_length < p->max_length);
	for (int i = 0; i < p->curr_length; i++)
		setIndex(p->used_vertices, p->path[i]);

	return (solution_vector) p;
}
