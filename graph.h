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
typedef struct list_node list_node;

graph* create_graph(int vertices);
void destroy_graph(graph* g);
void add_edge(graph* g, int from, int to, float weight);

void print_graph(graph* g);

#endif
