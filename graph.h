#ifndef GRAPH_ADT_DEFN
#define GRAPH_ADT_DEFN

struct list_node {
	struct list_node* next;
	int to;
	int from;
	float weight;
};

struct graph {
	int num_vertices;
	int num_edges;
	struct list_node** adjacency_list;
};

typedef struct graph graph;

graph* create_graph(int vertices);
void destroy_graph(graph* g);
void add_edge(graph* g, int from, int to, float weight);

#endif
