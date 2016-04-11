#ifndef PRIORITY_QUEUE_DEFN
#define PRIORITY_QUEUE_DEFN

#include <omp.h>
#include "generics.h"

#define LOCK(Q) omp_set_lock(&Q->global_lock)
#define UNLOCK(Q) omp_unset_lock(&Q->global_lock)

struct queue {
	int length;
	omp_lock_t global_lock;	// lock for coarse-grained parallelism
	struct queue_head* root_node;
};

struct queue_head {

	float priority;			// float for generality or int sufficient for our purpose ?
	solution_vector* partial_solution;

	int distance;
	int length;
	struct queue_head* left_subtree;
	struct queue_head* right_subtree;
};

struct queue* create_queue();
void destroy_queue(struct queue* q);
int pq_length(struct queue* q);

solution_vector* pq_min_extract(struct queue* q);
void pq_insert(struct queue* q, float priority, solution_vector* partial_solution);
void pq_merge(struct queue* q1, struct queue_head* q2, int length);
void pq_prune(struct queue* q, float min_bound);

//add function for bulk extract for load balancing


typedef struct queue queue;
typedef struct queue_head queue_head;

#endif