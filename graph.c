#include "graph.h"

#include <stdio.h>
#include <assert.h>

graph* create_graph(int vertices) {
	NEW(graph, g);
	g->num_vertices = vertices;
	g->adjacency_list = (list_node**) calloc(g->num_vertices, sizeof(list_node*));
	g->num_edges = 0;

	return g;
}

void destroy_graph(graph* g) {
	if (g == NULL)
		return;

	list_node* l, l1;

	for (int i = 0; i < g->num_vertices; i++) {
		l = g->adjacency_list[i];

		while (l != NULL) {
			l1 = l;
			l = l->next;

			free(l1);
		}
	}

	free(g->adjacency_list);
	free(g);
}

void add_edge(graph* g, int from, int to, float weight) {
	assert(g != NULL);
	assert(from >= 0);
	assert(to >= 0);
	assert(from < g->num_vertices);
	assert(to < g->num_vertices);
	assert(weight >= 0);

	NEW(list_node, new_edge);
	new_edge->from = from;
	new_edge->to = to;
	new_edge->weight = weight;

	list_node* l;
	l = g->adjacency_list[from];

	while (l != NULL && l->next != NULL) {
		assert(l->to != to);
		l = l->next;
	}

	if (l != NULL) {
		assert(l->to != to);
		l->next = new_edge;
	} else {
		g->adjacency_list[from] = new_edge;
	}

	g->num_edges ++;
}
