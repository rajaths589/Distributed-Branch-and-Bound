#include "generics.h"
#include "comforts.h"
#include "graph.h"
#include "bitvector.h"
#include "priority_queue.h"

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>

int solution_vector_size = 0;

typedef struct graph_color_data {
	graph* g;
	int max_colors;
} graph_color_data;

struct bound_array {
	int* vertex_colors;
	int curr_length;
	int max_length;
};

typedef struct bound_array color_assignment;

color_assignment* create_soln_copy(color_assignment* c) {
	NEW(color_assignment, c1);
	c1->max_length = c->max_length;
	c1->vertex_colors = (int*) calloc(c->max_length, sizeof(int));
	memcpy(c1->vertex_colors, c->vertex_colors, c->max_length*sizeof(int));
	c1->curr_length = c->curr_length;

	return c1;
}

void* populate_domain_data(int argc, char** argv) {
	assert(argc == 2);
	FILE* fp = fopen(argv[1], "r");

	assert(fp != NULL);

	int num_vertices, num_edges;
	NEW(graph_color_data, data);

	fscanf(fp, "%d %d", &num_vertices, &num_edges);
	graph* g = create_graph(num_vertices);

	int from, to, discard;

	for (int i = 0; i < num_edges; i++) {
		fscanf(fp, "%d %d %d", &from, &to, &discard);
		assert(from != to);
		add_edge(g, from, to, 0);
//		add_edge(g, to, from, 0);
	}

	data->g = g;
	data->max_colors = g->num_vertices;

	solution_vector_size = sizeof(float) + sizeof(int) + data->max_colors * sizeof(int);

	return (void*) data;
}

float get_root_partial_soln_score(void* domain_specific_data) {
	return 1.0;
}

solution_vector get_root_partial_solution(void* domain_specific_data) {
	NEW(color_assignment, empty_solution);

	graph_color_data* data = (graph_color_data*) domain_specific_data;

	empty_solution->max_length = data->max_colors;
	empty_solution->curr_length = 1;
	empty_solution->vertex_colors = (int*) calloc(data->max_colors, sizeof(int));
	empty_solution->vertex_colors[0] = 0;

	return (void*) empty_solution;
}

int construct_candidates(solution_vector partial_solution, float partial_soln_score,
						 queue* private_queue, void* domain_specific_data) {

	color_assignment* partial = (color_assignment*) partial_solution;
	graph_color_data* data = (graph_color_data*) domain_specific_data;

	if (partial->curr_length == data->max_colors)
		return 1;

	list_node* l;
	color_assignment* extension;

	l = data->g->adjacency_list[partial->curr_length];
	bitvector* clashing_colors = create_bitvector(data->max_colors);

	while (l != NULL) {
		if (l->to < partial->curr_length) {
			setIndex(clashing_colors, partial->vertex_colors[l->to]);
		}
		l = l->next;
	}
	bitvector* used_colors = create_bitvector(data->max_colors);
	for (int i = 0; i < partial->curr_length; i++) {
		setIndex(used_colors, partial->vertex_colors[i]);
	}

	for (int i = 0; i < data->max_colors; i++) {
		if (!getIndex(clashing_colors, i)) {
			extension = create_soln_copy(partial);
			extension->vertex_colors[extension->curr_length] = i;
			extension->curr_length++;
			if (getIndex(used_colors, i)) {
				pq_insert_nc(private_queue, partial_soln_score, extension);
			} else {
				pq_insert_nc(private_queue, partial_soln_score + 1, extension);
			}
		}
	}

	return 0;
}

void print_solution(solution_vector solution, float score) {
	color_assignment* cassign = (color_assignment*) solution;

	printf("Solution:\n");
	for (int i = 0; i < cassign->max_length; i++) {
		printf("%d\t", cassign->vertex_colors[i]);
	}
	printf("\n");
	printf("Colors: %d\n", (int) score);
}

void pack_solution(void* buff, int buff_size, int* pos, solution_vector vec, float score,
				   MPI_Comm comm, void* problem_data) {

	color_assignment* cassign = (color_assignment*) vec;
	MPI_Pack(&score, 1, MPI_FLOAT, buff, buff_size, pos, comm);
	MPI_Pack(&cassign->curr_length, 1, MPI_INT, buff, buff_size, pos, comm);
	MPI_Pack(cassign->vertex_colors, cassign->max_length, MPI_INT, buff, buff_size, pos, comm);
}

solution_vector unpack_solution(void* buff, int buff_size, MPI_Comm comm, int* pos,
								float* score, void* problem_data) {
	color_assignment* cassign = (color_assignment*) get_root_partial_solution(problem_data);

	MPI_Unpack(buff, buff_size, pos, score, 1, MPI_FLOAT, comm);
	MPI_Unpack(buff, buff_size, pos, &cassign->curr_length, 1, MPI_INT, comm);
	MPI_Unpack(buff, buff_size, pos, cassign->vertex_colors, cassign->max_length, MPI_INT, comm);
	return (solution_vector) cassign;
}
