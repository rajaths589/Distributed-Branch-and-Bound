#ifndef PRIORITY_QUEUE_DEFN
#define PRIORITY_QUEUE_DEFN

#include <omp.h>
#include "generics.h"

struct queue {
	int length;
	omp_lock_t global_lock;	// lock for coarse-grained parallelism
	struct queue_head* root_node;
};

struct queue_head {

	float priority;			// float for generality or int sufficient for our purpose ?
	solution_vector* partial_solution;

	int distance;
	struct queue_head* left_subtree;
	struct queue_head* right_subtree;
};

solution_vector* pq_min_extract(struct queue* q);
void pq_insert(struct queue* q, float priority, solution_vector* partial_solution);
void pq_merge(struct queue* q1, struct queue_head* q2, int length);
void pq_prune(struct queue* q, float min_bound);

//add function for bulk extract for load balancing

#endif