#include <stdlib.h>
#include <omp.h>
#include <assert.h>

#include "priority_queue.h"
#include "comforts.h"

queue_head* merge_helper(queue_head* q1, queue_head* q2);
queue_head* prune_helper(queue_head* qh, float min_bound);
void destroy_node(queue_head* qh);

int node_distance(queue_head* q) {
	if (q == NULL)
		return -1;

	return q->distance;
}

int node_length(queue_head* q) {
	if (q == NULL)
		return 0;

	return q->length;
}

queue* create_queue() {
	NEW(queue, q);
	q->length = 0;
	q->root_node = NULL;
	omp_init_lock(&q->global_lock);

	return q;
}

void destroy_queue(queue* q) {
	destroy_node(q->root_node);
	omp_destroy_lock(&q->global_lock);
	free(q);
}

void destroy_node(queue_head* qh) {
	if (qh != NULL) {
		destroy_node(qh->left_subtree);
		destroy_node(qh->right_subtree);

		free(qh);
	}
}

int pq_length(queue* q) {
	int l;

	LOCK(q);
	l = q->length;
	UNLOCK(q);

	return l;
}

solution_vector pq_min_extract(queue* q, float* pr) {
	solution_vector* min_vector;

	LOCK(q);

	if (q->length == 0) {
		min_vector = NULL;
	} else {
		min_vector = q->root_node->partial_solution;
		*pr = q->root_node->priority;

		queue_head *left, *right;
		left = q->root_node->left_subtree;
		right = q->root_node->right_subtree;

		q->root_node = merge_helper(left, right);
		q->length --;
	}

	UNLOCK(q);

	return min_vector;
}

//not concurrent
void pq_insert_nc(queue* q, float priority, solution_vector partial_solution) {
	NEW(queue_head, qh);

	qh->priority = priority;
	qh->partial_solution = partial_solution;
	qh->left_subtree = NULL;
	qh->right_subtree = NULL;
	qh->distance = 0;

	q->root_node = merge_helper(q->root_node, qh);
	q->length ++;
}


void pq_merge(queue* q1, queue* q2) {
	LOCK(q1);

	q1->root_node = merge_helper(q1->root_node, q2->root_node);
	q1->length += q2->length;

	UNLOCK(q1);

	q2->length = 0;
	q2->root_node = NULL;
}

void pq_prune(queue* q, float min_bound) {
	LOCK(q);

	q->root_node = prune_helper(q->root_node, min_bound);
	if (q->root_node == NULL)
		q->length = 0;
	else
		q->length = q->root_node->length;

	UNLOCK(q);
}

queue_head* merge_helper(queue_head* q1, queue_head* q2) {
	queue_head* temp;

	if (q2 == NULL)
		return q1;

	if (q1 == NULL)
		return q2;

	if (q2->priority < q1->priority) {
		SWAP(q1, q2);
	}

	q1->right_subtree = merge_helper(q1->right_subtree, q2);
	if (node_distance(q1->right_subtree) > node_distance(q1->left_subtree)) {
		SWAP(q1->left_subtree, q1->right_subtree);
	}

	if (q1->right_subtree == NULL)
		q1->distance = 0;
	else
		q1->distance = 1 + node_distance(q1->right_subtree);

	return q1;
}

queue_head* prune_helper(queue_head* q, float min_bound) {
	queue_head* temp;

	if (q == NULL)
		return NULL;

	if (q->priority > min_bound)
		return NULL;

	q->left_subtree = prune_helper(q->left_subtree, min_bound);
	q->right_subtree = prune_helper(q->right_subtree, min_bound);

	if (node_distance(q->right_subtree) > node_distance(q->left_subtree)) {
		SWAP(q->left_subtree, q->right_subtree);
	}

	if (q->right_subtree == NULL) {
		q->distance = 0;
		if (q->left_subtree != NULL)
			q->length = 1 + q->left_subtree->length;
		else
			q->length = 1;
	}
	else {
		q->distance = 1 + node_distance(q->right_subtree);
		q->length = 1 + q->left_subtree->length + q->right_subtree->length;
	}

	return q;
}

queue_head* extract_helper(queue_head* qh, int num) {
	queue_head *temp_r, *temp;

	if (MAX(node_length(qh->left_subtree), node_length(qh->right_subtree)) ==
		node_length(qh->left_subtree)) {
		if (node_length(qh->left_subtree) == num) {
			temp_r = qh->left_subtree;
			qh->left_subtree = NULL;
		} else {
			temp_r = extract_helper(qh->left_subtree, num);
		}
		} else {
			if (node_length(qh->right_subtree) == num) {
				temp_r = qh->right_subtree;
				qh->right_subtree = NULL;
			} else {
				temp_r = extract_helper(qh->right_subtree, num);
			}
		}

		if (node_distance(qh->right_subtree) > node_distance(qh->left_subtree)) {
			SWAP(qh->left_subtree, qh->right_subtree);
		}
		qh->distance = MIN(node_distance(qh->left_subtree), node_distance(qh->right_subtree)) +
		1;
		qh->length = node_length(qh->left_subtree) + node_distance(qh->right_subtree);

	return temp_r;
}

queue_head* pq_extract(struct queue* q, int num) {
	assert(num > 0);

	LOCK(q);

	assert(num < q->length);

	queue_head* qh = extract_helper(q->root_node, num);
	q->length = q->root_node->length;

	UNLOCK(q);

	return qh;
}
