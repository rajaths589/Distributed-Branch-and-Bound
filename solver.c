#include <omp.h>
#include <mpi.h>

#include "priority_queue.h"

void expand_partial_solution(queue* shared_queue, queue* private_queue, data domain_data) {
	float score;
	solution_vector* partial_solution = pq_min_extract(shared_queue, &score);

	if (partial_solution != NULL) {
		construct_candidates(partial_solution, score, private_queue, domain_data);

		pq_merge(shared_queue, private_queue);
	}
}

void solution_main(queue* shared_queue, data domain_data) {

}