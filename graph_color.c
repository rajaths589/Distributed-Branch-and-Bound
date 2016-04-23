#include "generics.h"
#include "comforts.h"
#include "graph.h"
#include "bitvector.h"
#include "priority_queue.h"

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>


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
	memcpy(c1->vertex_colors, c->vertex_colors, c->max_length);
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

	int from, to;

	for (int i = 0; i < num_edges; i++) {
		fscanf(fp, "%d %d %f", from, to, weight);
		assert(from != to);
		add_edge(g, from, to);
		add_edge(g, to, from);
	}

	data->g = g;
	data->max_colors = g->num_vertices;

	return (void*) data;
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

	l = data->g->adjacency_list[partial->curr_length-1];
	bitvector* clashing_colors = create_bitvector(data->max_colors);

	while (l != NULL) {
		if (l->to < partial->curr_length) {
			setIndex(clashing_colors, partial->vertex_colors[l->to]);
		}
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
	printf("Solution:\n");
	for (int i = 0; i < solution->max_length; i++) {
		printf("%d\t", solution->vertex_colors[i]);
	}
	printf("\n");
	printf("Colors: %d\n", (int) score);
}
